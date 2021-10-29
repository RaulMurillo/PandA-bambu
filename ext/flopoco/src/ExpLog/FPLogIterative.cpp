/*
  An FP logarithm for FloPoCo

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

*/

// TODO List:
//  * test cases for boundary cases pfinal etc
//  * finetune pipeline
//  * Port back the Arith paper
#include <fstream>
#include <sstream>
#include <math.h>	// for NaN
#include "FPLogIterative.hpp"
#include "FPLog.hpp"
//#include "TestBenches/FPNumber.hpp"
#include "ShiftersEtc/LZOC.hpp"
#include "ShiftersEtc/Normalizer.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "Table.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "utils.hpp"
//#include "IntMult/IntSquarer.hpp"
//#include "ConstMult/IntIntKCM.hpp"
#include "UserInterface.hpp"


using namespace std;


namespace flopoco{

#define OLD 0


#if 0 // I keep this because it looks like it used to be important
	
	int FPLogIterative::FirstInvTable::check_accuracy(int wF) {
		int i;
		mpz_class j;
		double x1,x2,y,e1,e2;
		double maxerror=0.0;
		double prod=0.0;

		maxMulOut=0;
		minMulOut=2;

		for (i=minIn; i<=maxIn; i++) {
			// x1 and x2 are respectively the smallest and largest FP possible
			// values leading to input i
			x1=input2double(i);
			if(i>>(wIn-1)) //MSB of input
				x2= - negPowOf2(wF)          // <--wF --->
					+ ((double)(i+1+(1<<wIn))) //   11 11...11 (binary)
					/ ((double)(1<<(wIn+1))); // 0.11 11...11 (binary)
			else
				x2= - negPowOf2(wF-1)
					+ ((double)(i+1+(1<<wIn))) //  10 11...11 (binary)
					/ ((double)(1<<(wIn))); // 1.0 11...11 (binary)
			j=function(i);
			y=output2double(j);
			if(UserInterface::verbose)
				cout << "i="<<i<< " ("<<input2double(i)<<") j="<<j
					  <<" min="<< x1*y <<" max="<< x2*y<< endl;
			prod=x1*y; if (prod<minMulOut) minMulOut=prod;
			prod=x2*y; if (prod>maxMulOut) maxMulOut=prod;
			e1= fabs(x1*y-1); if (e1>maxerror) maxerror=e1;
			e2= fabs(x2*y-1); if (e2>maxerror) maxerror=e2;
		}
		cout << "FirstInvTable: Max error=" <<maxerror << "  log2=" << log2(maxerror) <<endl;
		cout << "               minMulOut=" <<minMulOut << " maxMulOut=" <<maxMulOut  <<endl;

		printf("%1.30e\n", log2(maxerror));

		return (int) (ceil(log2(maxerror)));
	}

#endif





	FPLogIterative::FPLogIterative(OperatorPtr parentOp, Target* target,
																 int wE, int wF,
																 int inTableSize)
		: FPLog(parentOp, target, wE, wF)
	{

		setCopyrightString("F. de Dinechin, C. Klein  (2008-2011)");

		ostringstream o;
		srcFileName = "FPLogIterative";

		o << "FPLogIterative_" << wE << "_" << wF << "_" << inTableSize << "_";
		if(getTarget()->isPipelined())
			o << getTarget()->frequencyMHz() ;
		else
			o << "comb";
		setNameWithFreqAndUID(o.str());

		addFPInput("X", wE, wF);
		addFPOutput("R", wE, wF, 2); // 2 because faithfully rounded

		int i, bitsPerStage;

		if (inTableSize==0) {
			if (wF>=20)
				bitsPerStage=12;
			else
				bitsPerStage=8; // somehow arbitrary
		}
		else
			bitsPerStage=inTableSize;

		if(bitsPerStage <6)
			throw string("FPLogIterative error: tables need at least 6 input bits");



		// First compute the parameters of each iteration

		// stage 0, needs a specific inverter table
		p[0] = 0;
		a[0] = bitsPerStage;
		// How many bits have been zeroed?  a0=5 is a lucky case but it messes up the following
		p[1] = a[0]-2; // if a0==5, would be a0-1

		// The number of guard bit
		// For each stage: 0.5 ulp for rounding the log table
		// For each stage after the first one:
		//                 1 ulp for truncating Z
		//                 1 ulp for truncating P
		//                 1 ulp for truncating EiY
		//        total 3.5 bit/stage
		// Plus for the rest of the computation:
		//                 1 ulp for the approximation error of neglecting the 3rd-order term
		//                 1 ulp for truncating Z before input to the square
		//                 1 ulp for truncating Z^2

		// the iterations i=0 and i=1 lead to no truncation

		i=1;
		gLog=4;
		while(3*p[i]+1 <= p[i]+wF+gLog){ // while the third-order term of log(1+Zi) is not negligible
			if(i==1) { 	// stage 1 cannot be more accurate than 2p1-1
				a[1] = p[1];
				p[2] = 2*p[1]-1;
			}
			else {
				a[i] = bitsPerStage;
				p[i+1] = p[i] + a[i] - 1; // and we zero out a[i]-1 bits
			}
			i++;
			gLog=max(4, intlog2(3+0.5+0.5+3*i-1));
		}

		// The number of stages, not counting stage 0
		stages = i-1;
		gLog=max(4, intlog2(3+0.5+0.5+3*stages));

		if(UserInterface::verbose>=2) {
			cerr << "> FPLogIterative\t Initial parameters:" << endl;
			for(i=0; i<=stages; i++) {
				cerr << "> FPLogIterative\t";
				cerr<<"\tp"<<i<<"=" << p[i];
				cerr<<"\ta"<<i<<"=" << a[i];
				cerr <<endl;
			}
		}
		// Now we probably have reduced too far
		pfinal = p[stages+1];
		int extraBits = pfinal - ((wF+gLog-2)>>1);
		int extraBitsperstage =  floor(((double)extraBits) / ((double)(stages+1)));
		if(UserInterface::verbose)
			cerr << "> FPLogIterative\t before optimization:  pfinal=" << pfinal << "--- extraBits=" << extraBits << "---extraBitsperstage=" << extraBitsperstage << endl;
		if(bitsPerStage>6) {
			for (i=0; i<= stages; i++)
				a[i] -= extraBitsperstage;
		}
		else  {
			for (i=2; i<= stages; i++)
				a[i] -= extraBitsperstage;
		}
		// Recompute the pi
		p[1] = a[0]-2;

		for(i=1; i<=stages; i++){ // for faithful rounding
			if(i==1)   	// stage 1 cannot be more accurate than 2p1-1
				p[2] = p[1] + a[1] -1; // 2p[1]-1
			else
				p[i+1] = p[i] + a[i] - 1;
		}
		pfinal = p[stages+1];
		extraBits = pfinal -  ((wF+gLog-2)>>1);
		if(stages>=2) { // remove extra bits from stages 2 to stages
			extraBitsperstage =  floor(((double)extraBits) / ((double)(stages-1)));
			int extraBitsStage2 = extraBits - (stages-1)*extraBitsperstage;
			a[2] -= extraBitsperstage + extraBitsStage2;
			p[2+1] = p[2] + a[2] - 1;
			for(i=3; i<= stages; i++) {
				a[i] -= extraBitsperstage;
				p[i+1] = p[i] + a[i] - 1;
			}
			pfinal = p[stages+1];
		}

		extraBits = pfinal -  ((wF+gLog-2)>>1);

		if(UserInterface::verbose>=2)
			cerr << "> FPLogIterative\t after optimization:   pfinal=" << pfinal << "--- extraBits=" << extraBits << endl;


		if(UserInterface::verbose)
			cerr << "> FPLogIterative"<<tab<<"Guard bits: " << gLog << endl;


		target_prec = wF + pfinal +gLog;
		if(UserInterface::verbose==2)
			cerr << "> FPLogIterative"<<tab<<"Target precision: " << target_prec << endl;

		s[0] = wF+2;
		psize[0] = s[0] + a[0]+1;
		//	sfullZ[0] = wF+2;
		sbt[1] = psize[0] - p[1] -2; // -2 because P0 begins with 01   was: wF+2 ;
		s[1] = sbt[1];
		t[1] = 0;
		//	sfullZ[1] = sfullZ[0]      +    a[0] + 1;
		//         size of Y0          size of approx inv of A0

		for(i=1; i<=stages; i++) {
			// size of Z before truncation; Zi = 0Bi - 0AiZi + EiYi ;
			// AiZi has weight -2*p[i], and size a[i]+s[i]
			// EiYi has weight -2*p[i]-1 and size 1+p[i]+s[i]
			// except i=1: same weight but size may be 1+p[i]+s[i]+1
			if(i==1)
				sbt[i+1] =    max( a[i]+s[i]+1, 1+p[i]+sbt[i]+1);
			else
				sbt[i+1] =    max( a[i]+s[i]+1, 1+p[i]+sbt[i]);

			if(p[i+1]+sbt[i+1] <= target_prec)
				{ // No truncation at all
					psize[i] = s[i];
					s[i+1] = sbt[i+1];
					t[i+1] = 0;
				}
			else
				{ // Truncate everybody to targetprec :
					// we compute Z[i+1] = B[i] - A[i]Z[i] + (1+Z[i])>>Ei[i]
					// Product  A[i]Z[i] has MSB 2*p[i], LSB target_prec, therefore size target_prec - 2*p[i]
					// We need only this many bits of Z[i] to compute it.
					psize[i] = target_prec - 2*p[i];
					if (psize[i]>s[i]) // in the first iterations
						psize[i]=s[i];
					s[i+1] = target_prec - p[i+1];
					t[i+1] = sbt[i+1] - s[i+1];
				}
		}

		sfinal =  s[stages+1];

		// MSB of squarer input is p[stages+1];
		// LSB will be target_prec
		// size will be target_prec - p[stages+1]


		if(UserInterface::verbose)
			cerr<<"> FPLogIterative\t needs 1+"<<stages<<" range reduction stages"<<endl;
		if(UserInterface::verbose>=2) {
			for(i=0; i<=stages; i++) {
				cerr << "> FPLogIterative\t";
				cerr<<"\tp"<<i<<"=" << p[i];
				cerr<<"\ta"<<i<<"=" << a[i];
				cerr<<"\ts"<<i<<"=" << s[i];
				cerr<<"\tpsize"<<i<<"=" << psize[i];
				cerr <<endl;
			}
			cerr << "> FPLogIterative\t\tsfinal=" << sfinal << "\tpfinal=" << pfinal << endl;

		}


		// On we go with the vhdl
		addConstant("g",  "positive",          gLog);
		addConstant("wE", "positive",          wE);
		addConstant("wF", "positive",          wF);
		addConstant("log2wF", "positive",     intlog2(wF));
		addConstant("targetprec", "positive", target_prec);
		addConstant("sfinal", "positive",     s[stages+1]);
		addConstant("pfinal", "positive",     p[stages+1]);

		vhdl << tab << declare("XExnSgn", 3) << " <=  X(wE+wF+2 downto wE+wF);" << endl;
		vhdl << tab << declare("FirstBit") << " <=  X(wF-1);" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(), // this is a MUX
													 "Y0", wF+2) << " <= \"1\" & X(wF-1 downto 0) & \"0\" when FirstBit = '0' else \"01\" & X(wF-1 downto 0);" << endl;
		vhdl << tab << declare("Y0h", wF) << " <= Y0(wF downto 1);" << endl;
		vhdl << tab << "-- Sign of the result;" << endl;
		vhdl << tab << 	declare(getTarget()->logicDelay(),"sR") << " <= '0'   when  (X(wE+wF-1 downto wF) = ('0' & (wE-2 downto 0 => '1')))  -- binade [1..2)" << endl
		                              << "     else not X(wE+wF-1);                -- MSB of exponent" << endl;
		
		vhdl << tab << declare(getTarget()->adderDelay(wF-pfinal+2),
													 "absZ0", wF-pfinal+2) << " <=   Y0(wF-pfinal+1 downto 0)          when (sR='0') else" << endl
			  << "             ((wF-pfinal+1 downto 0 => '0') - Y0(wF-pfinal+1 downto 0));" << endl;
		vhdl << tab << declare(getTarget()->adderDelay(wE),
													 "E", wE) << " <= (X(wE+wF-1 downto wF)) - (\"0\" & (wE-2 downto 1 => '1') & (not FirstBit));" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(),
													 "absE", wE) << " <= ((wE-1 downto 0 => '0') - E)   when sR = '1' else E;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(),
													 "EeqZero") << " <= '1' when E=(wE-1 downto 0 => '0') else '0';" << endl;
		
		newInstance("LZOC",
								"lzoc1",
								"wIn=" + to_string(wF) + " countType=-1",
								"I=>Y0h, OZB=>FirstBit",
								"O=>lzo");

		vhdl << tab << declare("pfinal_s", intlog2(wF)) << " <= \"" << unsignedBinary(mpz_class(pfinal), intlog2(wF)) << "\";"<<endl;
		vhdl << tab << declare(getTarget()->adderDelay(intlog2(wF)+1),
													 "shiftval", intlog2(wF)+1) << " <= ('0' & lzo) - ('0' & pfinal_s); " << endl;
		vhdl << tab << declare("shiftvalinL", intlog2(wF-pfinal+2))     << " <= shiftval(" << intlog2(wF-pfinal+2)-1 << " downto 0);" << endl;
		vhdl << tab << declare("shiftvalinR", intlog2(sfinal-pfinal+1)) << " <= shiftval(" << intlog2(sfinal-pfinal+1)-1 << " downto 0);" << endl;
		vhdl << tab << declare("doRR") << " <= shiftval(log2wF); -- sign of the result" << endl;

		//done in parallel with the shifter
		vhdl << tab << declare(getTarget()->logicDelay(),"small") << " <= EeqZero and not(doRR);" << endl;

		// ao stands for "almost one"
		vhdl << tab << "-- The left shifter for the 'small' case" <<endl;


		newInstance("Shifter",
								"small_lshift",
								"wX=" + to_string(wF-pfinal+2) + " maxShift=" + to_string(wF-pfinal+2) + " dir=0",
								"X=>absZ0,S=>shiftvalinL",
								"R=>small_absZ0_normd_full");



		int small_absZ0_normd_size = getSignalByName("small_absZ0_normd_full")->width() - (wF-pfinal+2);
		vhdl << tab << declare("small_absZ0_normd", small_absZ0_normd_size) << " <= small_absZ0_normd_full" << range(small_absZ0_normd_size -1, 0) << "; -- get rid of leading zeroes" << endl;

		vhdl << tab << "---------------- The range reduction box ---------------" << endl;

		vhdl << tab << declare("A0", a[0]) << " <= X" << range(wF-1,  wF-a[0]) << ";" << endl;

		vhdl << tab << "-- First inv table" << endl;

		vector<mpz_class> it0Content;
		for(int x=0; x<(1<<a[0]); x++) {
			// convert x to a double
			double xd;
			if(x>>(a[0]-1)) //MSB of input
				xd= ((double)(x+(1<<a[0]))) //      11xxxx (binary)
				/  ((double)(1<<(a[0]+1))); // 0.11xxxxx (binary)
			else
				xd= ((double)(x+(1<<a[0]))) //   10xxxx (binary)
					/  ((double)(1<<a[0])); // 1.0xxxxx (binary)
			mpz_class r =  ceil( ((double)(1<<a[0])) / xd); // double rounding, but who cares really
		// The following line allows us to prove that the case a0=5 is lucky enough to need one bit less than the general case
		//cout << "FirstInvTable> x=" << x<< "  r=" <<r << "  error=" << ceil( ((double)(1<<(wOut-1))) / d)  - ( ((double)(1<<(wOut-1))) / d) << endl;
			it0Content.push_back(r);
		}
		Table::newUniqueInstance(this, "A0", "InvA0",
														 it0Content,
														 "InvA0Table",
														 a[0], a[0]+1);
		





		if(false && getTarget()->isPipelined() && a[0]+1>=9) {
			// TODO currently disabled, check it works
			// This is a full (untruncated multiplier)
			newInstance("IntMultiplier",
								"p0_mult",
								"wX=" + to_string(a[0]+1) + " wY=" + to_string(wF+2),
								"X=>InvA0, Y=>Y0",
								"R=>P0");
		}
		else {
			REPORT(DETAILED, "unpipelined multiplier for P0, implemented as * in VHDL");
			vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(psize[0]), // TODO very approximate timing for a very flat mult
														 "P0",  psize[0]) << " <= InvA0 * Y0;" <<endl <<endl;
		}


		vhdl << tab << declare("Z1", s[1]) << " <= P0" << range (psize[0] - p[1]-3,  0) << ";" << endl;

		for (i=1; i<= stages; i++) {

			vhdl <<endl;
			//computation
			vhdl << tab << declare(join("A",i), a[i]) <<      " <= " << join("Z",i) << range(s[i]-1,      s[i]-a[i]) << ";"<<endl;
			vhdl << tab << declare(join("B",i), s[i]-a[i]) << " <= " << join("Z",i) << range(s[i]-a[i]-1, 0        ) << ";"<<endl;

			vhdl << tab << declare(join("ZM",i), psize[i]) << " <= " << join("Z",i) ;
			if(psize[i] == s[i])
				vhdl << ";"<<endl;
			else
				vhdl << range(s[i]-1, s[i]-psize[i])  << ";" << endl;

			if(false){ // TODO currently unplugged, resurrect when the rest works 

			newInstance("IntMultiplier",
								join("p",i)+"_mult",
								"wX=" + to_string(a[0]+1) + " wY=" + to_string(a[0]+1),
								"X=>InvA0, Y=>Y0",
								"R=>P0");
			}
			else {
				REPORT(DETAILED,"Unpipelined multiplier for P"<<i<<", implemented as * in VHDL");
				vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(psize[i]+a[i]), // TODO very approximate timing for a very flat mult
															 join("P",i),  psize[i] + a[i]) << " <= " << join("A",i) << "*" << join("ZM",i) << ";" << endl;
			}

			int yisize = s[i]+p[i]+1; // +1 because implicit 1

			// While the multiplication takes place, we may prepare the other term

			string Yi= join("Y",i);
			string Zi= join("Z",i);
			vhdl << tab << declare(Yi, yisize) << " <= \"1\" & " << rangeAssign(p[i]-1, 0, "'0'") << " & " << Zi <<";"<<endl;

			// We truncate EiY to a position well above target_prec
			if(i==1) {
				//manageCriticalPath( getTarget()->lutDelay());
				vhdl << tab << declare(getTarget()->logicDelay(),
															 join( "EiY",i), s[i+1])
						 << " <= " << Yi ;
				// now we may need to truncate or to pad Yi
				if(yisize > s[i+1]) // truncate Yi
					vhdl << range(yisize-1, yisize-s[i+1]);
				else if (yisize < s[i+1]) // pad Yi
					vhdl << " & " << rangeAssign(s[i+1] - yisize -1, 0, "'0'");
				vhdl  << "  when " << join("A",i) << of(a[i]-1) << " = '1'" << endl
						<< tab << "  else  \"0\" & " << Yi;
				if(yisize > s[i+1]-1) // truncate Yi
					vhdl << range(yisize-1, yisize-s[i+1]+1);
				vhdl 	  << ";" << endl ;
			}
			else { // i>1, general case
				vhdl << tab << declare(join("EiY",i), s[i+1]) << " <= " ;
				vhdl << rangeAssign(2*p[i] -p[i+1] -2,  0,  "'0'")  << " & " << Yi << range(yisize-1,  yisize - (s[i+1] - (2*p[i]-p[i+1]) )-1) << ";" << endl ;
			}

			vhdl << tab << declare( join("addXIter",i), s[i+1] ) << " <= \"0\" & " << join("B",i);
			if (s[i+1] > 1+(s[i]-a[i]))  // need to padd Bi
				vhdl << " & " << rangeAssign(s[i+1] - 1-(s[i]-a[i]) -1,  0 , "'0'");
			vhdl <<";"<<endl;

			newInstance("IntAdder",
									join("addIter1_",i),
									"wIn="+to_string(s[i+1]),
									"X=>"+join("addXIter",i) + ", Y=>"+join("EiY",i),
									"R=>"+join("EiYPB",i),
									"Cin=>'0'");

			vhdl << tab << declare(getTarget()->logicDelay(),
														 join("Pp", i), s[i+1])  << " <= " << rangeAssign(p[i]-a[i],  0,  "'1'") << " & not(" << join("P", i);
			// either pad, or truncate P
			if(p[i]-a[i]+1  + psize[i]+a[i]  < s[i+1]) // size of leading 0s + size of p
				vhdl << " & "<< rangeAssign(s[i+1] - (p[i]-a[i]+1  + psize[i]+a[i]) - 1,    0,  "'0'");  // Pad
			if(p[i]-a[i]+1  + psize[i]+a[i]  > s[i+1])
				//truncate
				vhdl << range(psize[i]+a[i]-1,    p[i]-a[i]+1  + psize[i]+a[i] - s[i+1]);
			vhdl << ");"<<endl;

			newInstance("IntAdder",
									join("addIter2_",i),
									"wIn="+to_string(s[i+1]),
									"X=>"+join("EiYPB",i) + ", Y=>"+join("Pp",i),
									"R=>"+join("Z",i+1),
									"Cin=>'1'");

		}

		vhdl << tab << declare("Zfinal", s[stages+1]) << " <= " << join("Z", stages+1) << ";" << endl;


		// In the small path we need Z2O2 accurate to  (wF+gLog+2) - pfinal
		// In the RR path we need Z2O2 accurate to sfinal-pfinal
		// Take the max. This is useful  for small precs only
		int squarerInSize;
		if((wF+gLog+2) - pfinal > sfinal-pfinal)
			squarerInSize = (wF+gLog+2) - pfinal;
		else
			squarerInSize = sfinal-pfinal;

		vhdl << tab << declare(getTarget()->logicDelay(), "squarerIn", squarerInSize) << " <= "
			 << "Zfinal(sfinal-1 downto sfinal-"<< squarerInSize << ") when doRR='1'" << endl;
		if(squarerInSize>small_absZ0_normd_size)
			vhdl << tab << "                 else (small_absZ0_normd & " << rangeAssign(squarerInSize-small_absZ0_normd_size-1, 0, "'0'") << ");  " << endl;
		else  // sfinal-pfinal <= small_absZ0_normd_size
			vhdl << tab << "                 else small_absZ0_normd" << range(small_absZ0_normd_size-1, small_absZ0_normd_size - squarerInSize) << ";  " << endl<< endl;

#if OLD
		IntSquarer* sq = new IntSquarer(target, squarerInSize);

		inPortMap  (sq, "X", "squarerIn");
		outPortMap (sq, "R", "Z2o2_full");
		vhdl << instance(sq, "squarer");
		
#else // OLD
		// TODO replace with a squarer
		vhdl << tab << declare("Z2o2_full", 2*squarerInSize) << " <= squarerIn*squarerIn;" << endl; 
#endif
		
		vhdl << tab << declare("Z2o2_full_dummy", 2*squarerInSize) << " <= Z2o2_full;" << endl;

		vhdl << tab << declare("Z2o2_normal", sfinal-pfinal-1) << " <= Z2o2_full_dummy ("<< 2*squarerInSize-1 << "  downto " << 2*squarerInSize - (sfinal-pfinal-1) << ");" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(), "addFinalLog1pY", sfinal) << " <= (pfinal downto 0  => '1') & not(Z2o2_normal);" <<endl;


		newInstance("IntAdder",
								"addFinalLog1p_normalAdder",
								"wIn="+to_string(sfinal),
								"X=>Zfinal, Y=>addFinalLog1pY",
								"R=>Log1p_normal",
								"Cin=>'1'");

		vhdl << endl << tab << "-- Now the log tables, as late as possible" << endl;

		// the log tables have small input and large outputs, so we had rather register the inputs.
		// We have to compute the sum of the outputs of the log tables, and we want this sum (AlmostLog) to be synchronized to Log1pNormal, to which it will be added.
		// To get this synchro right, we have to do a first dummy evaluation of the pipeline to profile its depth
		// We first do it in a dummy way, starting from cycle 0, to measure the depth of this sub-pipeline

		// First log table: this used to be better structured code 
		vector<mpz_class> lt0Content;
		for(int x=0; x<(1<<a[0]); x++) {
			mpz_class result;
			mpfr_t i,l;
			mpz_t r;
			mpfr_init(i);
			mpfr_init2(l,target_prec);
			mpz_init2(r,400);

			mpz_class fitout = it0Content[x];
			double apprinv = ((double)fitout.get_si()) /  ((double)(1<<a[0])); 
			mpfr_set_d(i, apprinv, GMP_RNDN);
			mpfr_log(l, i, GMP_RNDN);
			mpfr_neg(l, l, GMP_RNDN);
			
			// Remove the sum of small offsets that are added to the other log tables
			for(int j=1; j<=stages; j++){
				mpfr_set_d(i, 1.0, GMP_RNDN);
				int pi=p[j];
				mpfr_mul_2si(i, i, -2*pi, GMP_RNDN);
				mpfr_sub(l, l, i,  GMP_RNDN);
			}
			// code the log in 2's compliment
			mpfr_mul_2si(l, l, target_prec, GMP_RNDN);
			mpfr_get_z(r, l, GMP_RNDN);
			result = mpz_class(r); // signed

			// This is a very inefficient way of converting
			mpz_class t = mpz_class(1) << target_prec;
			result = t+result;
			if(result>t) result-=t;
			//  cout << "x="<<x<<" apprinv="<<apprinv<<" logapprinv="<<log(apprinv)<<" result="<<result<<endl;
			mpfr_clear(i);
			mpfr_clear(l);
			mpz_clear(r);
			lt0Content.push_back(result);
		}
		Table::newUniqueInstance(this, "A0", "L0",
														 lt0Content,
														 "LogTable0",
														 a[0], target_prec);

		vhdl << tab << declare("S1", target_prec) << " <= L0;"<<endl;

		for (i=1; i<= stages; i++) {

			vector<mpz_class> ltiContent;
			for(int x=0; x<(1<<a[i]); x++) {
				mpz_class result;
				mpfr_t xmp,l;
				mpz_t r;

				int wIn=a[i];
				int wOut=target_prec - p[i];
				mpfr_init(xmp);
				mpfr_init2(l,wOut);
				mpz_init2(r,400);

				// input x to double 
				double  xd;
				// computation in double is exact as long as we don't want a quad operator...
				double Ei;
				if((i>1) || (i==1 && (1==(x>>(wIn-1)))) )
					Ei = 1.0 / ((double) (((uint64_t) 1)<<(2*p[i])));
				else
					Ei = 1.0 / ((double) (((uint64_t) 1)<<(2*p[i]+1)));

				xd = ((double) (-x))   /   ((double) (((uint64_t) 1)<<(p[i]+wIn)));

				//cout << endl << d << " " << Ei << "   " ;
				xd += Ei;
				mpfr_set_d(xmp, xd, GMP_RNDN);
				mpfr_log1p(l, xmp, GMP_RNDN);
				mpfr_neg(l, l, GMP_RNDN);
				// cout << "i=" << i <<  " div" << (p[i]+wIn+1) << "  x=" << x << "  apprinv=" << apprinv << "  l=" << mpfr_get_d(l, GMP_RNDN) << endl;
				// Add the small offset that ensures that it never gets negative (the sum of these offsets will be subtracted to L0)
				mpfr_set_d(xmp, 1.0, GMP_RNDN);
				mpfr_mul_2si(xmp, xmp, -2*p[i], GMP_RNDN);
				mpfr_add(l, l, xmp,  GMP_RNDN);

				mpfr_mul_2si(l, l, p[i]+wOut, GMP_RNDN);
				mpfr_get_z(r, l, GMP_RNDN);
				result=mpz_class(r);
				mpfr_clear(xmp);
				mpfr_clear(l);
				mpz_clear(r);

				ltiContent.push_back(result);
			}
		Table::newUniqueInstance(this, join("A",i), join("L",i),
														 ltiContent,
														 join("LogTable",i),
														 a[i], target_prec - p[i]);



		vhdl << tab << declare(join("sopX",i), target_prec) << " <= (" << rangeAssign(target_prec-1, target_prec-p[i],  "'0'") << " & " << join("L",i) <<");"<<endl;

		newInstance("IntAdder",
								join("adderS",i),
								"wIn="+to_string(target_prec),
								"X=>"+join("S",i)+", Y=>"+join("sopX",i),
								"R=>"+join("S",i+1),
								"Cin=>'0'");

		}

		vhdl << tab << declare("almostLog", target_prec) << " <= " << join("S",stages+1) << ";" << endl;

		vhdl << tab << declare( "adderLogF_normalY", target_prec ) << " <= ((targetprec-1 downto sfinal => '0') & Log1p_normal);" << endl;

		newInstance("IntAdder",
								"adderLogF_normal",
								"wIn="+to_string(target_prec),
								"X=>almostLog, Y=>adderLogF_normalY",
								"R=>LogF_normal",
								"Cin=>'0'");

		newInstance("FixRealKCM",
								"MulLog2",
								// unsigned here, the conversion to signed comes later
								 "signedIn=0 msbIn=" + to_string(wE-1)
								+ " lsbIn=0"
								+ " lsbOut=" + to_string(-wF-gLog) 
								+ " constant=log(2)"
								+ " targetUlpError=1.0",
								"X=>absE",
								"R=>absELog2");
				

		vhdl << tab << declare("absELog2_pad", wE+target_prec) << " <=   "
			 << "absELog2 & (targetprec-wF-g-1 downto 0 => '0');       " << endl;

		vhdl << tab << declare("LogF_normal_pad", wE+target_prec) << " <= (wE-1  downto 0 => LogF_normal(targetprec-1))  & LogF_normal;" << endl;

		vhdl << tab << declare("lnaddX", wE+target_prec) << " <= absELog2_pad;"<<endl;


		vhdl << tab << declare("lnaddY", wE+target_prec) << " <= LogF_normal_pad when sR='0' else not(LogF_normal_pad); "<<endl;

		newInstance("IntAdder",
								"lnadder",
								"wIn="+to_string(wE+target_prec),
								"X=>lnaddX, Y=>lnaddY, Cin=>sR",
								"R=>Log_normal"
								);

		newInstance("Normalizer",
								"final_norm",
								"wX=" + to_string(wE+target_prec) + " wR=" + to_string(target_prec) + " maxShift=" + to_string(wE+(wF>>1)) + " countType=0",
								"X=>Log_normal",
								"R=>Log_normal_normd, Count=>E_normal");


		int Z2o2_small_size=(wF+gLog+2) - pfinal; // we need   (wF+gLog+2) - pfinal bits of Z2O2

		// if (syncCycleFromSignal("Z2o2_full_dummy"))
		// 	setCriticalPath(cpZ2o2_full_dummy);

		vhdl << tab << declare("Z2o2_small_bs", Z2o2_small_size)  << " <= Z2o2_full_dummy" << range(2*squarerInSize -1, 2*squarerInSize -Z2o2_small_size) << ";" << endl;

		newInstance("Shifter",
								"ao_rshift",
								"wX=" + to_string(Z2o2_small_size) + " maxShift=" + to_string(sfinal-pfinal+1) + " dir=1", // dir=right
								"X=>Z2o2_small_bs, S=>shiftvalinR",
								"R=>Z2o2_small_s");

		vhdl << tab << "  -- send the MSB to position pfinal" << endl;
		int Z2o2_small_sSize = getSignalByName("Z2o2_small_s")->width();
		vhdl << tab << declare("Z2o2_small", wF+gLog+2) << " <=  (pfinal-1 downto 0  => '0') & Z2o2_small_s"
			  << range(Z2o2_small_sSize-1,  Z2o2_small_sSize - (wF+gLog+2) + pfinal) << ";" << endl;

		vhdl << tab << "-- mantissa will be either Y0-z^2/2  or  -Y0+z^2/2,  depending on sR  " << endl;
		vhdl << tab << declare("Z_small", wF+gLog+2) << " <= small_absZ0_normd & " << rangeAssign((wF+gLog+2)-small_absZ0_normd_size-1, 0, "'0'") << ";" << endl;

		// manageCriticalPath( getTarget()->localWireDelay() + getTarget()->lutDelay() );
		vhdl << tab << declare("Log_smallY", wF+gLog+2) << " <= Z2o2_small when sR='1' else not(Z2o2_small);" << endl;
		vhdl << tab << declare("nsRCin",1, false) << " <= not ( sR );" << endl;

		newInstance("IntAdder",
								"log_small_adder",
								"wIn="+to_string(wF+gLog+2),
								"X=>Z_small, Y=>Log_smallY, Cin=>nsRCin",
								"R=>Log_small"
								);

		vhdl << tab << "-- Possibly subtract 1 or 2 to the exponent, depending on the LZC of Log_small" << endl;
		vhdl << tab << declare("E0_sub", 2) << " <=   \"11\" when Log_small(wF+g+1) = '1'" << endl
			  << "          else \"10\" when Log_small(wF+g+1 downto wF+g) = \"01\"" << endl
			  << "          else \"01\" ;" << endl;

		// Is underflow possible?
		vhdl << tab <<	"-- The smallest log will be log(1+2^{-wF}) \\approx 2^{-wF}  = 2^" << -wF <<  endl
			  << tab << "-- The smallest representable number is 2^{1-2^(wE-1)} = 2^" << 1 -(1<<(wE-1))<< endl;
		double cpE_small;
		if(1 -(1<<(wE-1)) < -wF) {
			vhdl << tab <<	"-- No underflow possible" <<  endl;
			vhdl << tab << declare("ufl") << " <= '0';" << endl;
			// manageCriticalPath( getTarget()->adderDelay(wE) );
			// cpE_small = getCriticalPath();
			vhdl << tab << declare("E_small", wE) << " <=  (\"0\" & (wE-2 downto 2 => '1') & E0_sub)  -  ";
			if(wE>getSignalByName("lzo")->width())
				vhdl << "((wE-1 downto " << getSignalByName("lzo")->width() << " => '0') & lzo) ;" << endl;
			else
				vhdl << "lzo;" << endl;
		}
		else{
			vhdl << tab <<	"-- Underflow may happen" <<  endl;
			// manageCriticalPath( getTarget()->adderDelay(wE) );
			// cpE_small = getCriticalPath();
			vhdl << tab << declare("E_small", wE+1) << " <=  (\"00\" & (wE-2 downto 2 => '1') & E0_sub)  -  (";
			vhdl << "'0' & lzo);" << endl;
			vhdl << tab << declare("ufl") << " <= E_small(wE);" << endl;
		}

		// double cpLog_small_normd;
		// manageCriticalPath( getTarget()->lutDelay() );
		// cpLog_small_normd = getCriticalPath();
		vhdl << tab << declare("Log_small_normd", wF+gLog) << " <= Log_small(wF+g+1 downto 2) when Log_small(wF+g+1)='1'" << endl
			  << "           else Log_small(wF+g downto 1)  when Log_small(wF+g)='1'  -- remove the first zero" << endl
			  << "           else Log_small(wF+g-1 downto 0)  ; -- remove two zeroes (extremely rare, 001000000 only)" << endl ;


		int E_normalSize = getSignalByName("E_normal")->width();

		// manageCriticalPath(getTarget()->lutDelay() + getTarget()->adderDelay(wE));
		// double cpER = getCriticalPath();

		vhdl << tab << declare("E0offset", wE) << " <= \"" << unsignedBinary((mpz_class(1)<<(wE-1)) -2 + wE , wE) << "\"; -- E0 + wE "<<endl;
		vhdl << tab << declare("ER", wE) << " <= E_small" << range(wE-1,0) << " when small='1'" << endl;
		if(wE>E_normalSize)
			vhdl << "      else E0offset - (" << rangeAssign(wE-1,  E_normalSize, "'0'") << " & E_normal);" << endl;
		else
			vhdl << "      else E0offset - E_normal;" << endl;

		// //////////////////////////////////////////// now back to log
		// setCycleFromSignal("E_normal");
		// setCriticalPath( cpE_normal );
		// if ( syncCycleFromSignal("Log_small_normd"))
		// 	setCriticalPath( cpLog_small_normd );

		// manageCriticalPath( getTarget()->lutDelay() );
		vhdl << tab << declare("Log_g", wF+gLog) << " <=  Log_small_normd(wF+g-2 downto 0) & \"0\" when small='1'           -- remove implicit 1" << endl
			  << "      else Log_normal_normd(targetprec-2 downto targetprec-wF-g-1 );  -- remove implicit 1" << endl ;

		// Sticky is always 1 for a transcendental function !
		// vhdl << tab << declare("sticky") << " <= '0' when Log_g(g-2 downto 0) = (g-2 downto 0 => '0')    else '1';" << endl;
		vhdl << tab << declare("round") << " <= Log_g(g-1) ; -- sticky is always 1 for a transcendental function " << endl;


		vhdl << tab << "-- if round leads to a change of binade, the carry propagation magically updates both mantissa and exponent" << endl;

		// if (syncCycleFromSignal("ER") )
		// 	setCriticalPath( cpER );

		vhdl << tab << declare("fraX", wE+wF) << " <= (ER & Log_g(wF+g-1 downto g)) ; " << endl;
		vhdl << tab << declare("fraY", wE+wF) << " <= ((wE+wF-1 downto 1 => '0') & round); " << endl;

		newInstance("IntAdder",
								"finalRoundAdder",
								"wIn="+to_string(wE+wF),
								"X=>fraX, Y=>fraY",
								"R=>EFR",
								"Cin=>'0'");

		vhdl << tab << declare(getTarget()->logicDelay(),"Rexn", 3)
				 << " <= \"110\" when ((XExnSgn(2) and (XExnSgn(1) or XExnSgn(0))) or (XExnSgn(1) and XExnSgn(0))) = '1' else" << endl
				 << "                              \"101\" when XExnSgn(2 downto 1) = \"00\"  else" << endl
				 << "                              \"100\" when XExnSgn(2 downto 1) = \"10\"  else" << endl
				 << "                              \"00\" & sR when (((Log_normal_normd(targetprec-1)='0') and (small='0')) or ( (Log_small_normd (wF+g-1)='0') and (small='1'))) or (ufl = '1') else" << endl
				 << "                               \"01\" & sR;" << endl;
		vhdl << tab << "R<=  Rexn & EFR;" << endl;

		REPORT(3, "Leaving constructor");
	}

	FPLogIterative::~FPLogIterative()
	{
	}

}
