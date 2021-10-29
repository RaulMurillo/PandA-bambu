/*
  Floating Point Divider for FloPoCo


  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Authors: Maxime Christ, Jeremie Detrey, Florent de Dinechin

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2015.
  All rights reserved.

 */




/*
Benchmarks for ./flopoco frequency=2 fpdiv we=8 wf=23 

current best is 42 : 543 LUTs  32.263ns   

(best FPAdder is 329 LUTs, 11,7ns )

 
REMARK: 
extraBit could possibly be reduced by 1: 
it is possible to replace the late  normalization step by an early comparison of 1.FX and 1.FY, 
 and 1-bit shift left 1.FX if 1.FX < 1.FY (it turns out this won't break the convergence condition on first step) 

Florent tested it and it is not worth it: we get a small saving in latency, no saving in area for odd wF, e.g.23: 544 LUT,  31.487ns
BUT large  overhead in latency and area for even wF (when nbDigits is not reduced). 
In short when nDigits is not reduced we have a  mantissa comparison that is a strict overhead.

 */

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "FPDiv.hpp"
#include "Table.hpp"

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0

	// BEWARE, only works for dMax<1
	vector<mpz_class> FPDiv::selFunctionTable(double dMin, double dMax, int nbBitD, int nbBitW, int alpha, int radix) {
		int wIn = nbBitD+nbBitW;
		int wOut = 1+intlog2(alpha); //Nb of bit for the integer part of w
		double rho = ((double)alpha)/(radix-1);

		if (dMax>1) {
			REPORT(0, endl << "FPDiv::selFunctionTable: Code untested for dMax>1" << endl);
		}
		// what is the weight of the MSB of d, the bits of D passed to selFunTable?
		// dMin between 0.5 and 0.75 -> msbWeight=0.5
		// dMin between 0.75 and 0.111_2 - > msbwWeight=1/4
		double dmsbweight=.5; // 
		while(dMin >= (1-dmsbweight/2))	{   // e.g when dMin=0.5 or 0.6, no loop
	 		dmsbweight *= 0.5;
		}
		vector<mpz_class> t;
		for (int x=0; x<(1<<wIn); x++) {					
			int w = (x>>nbBitD);	 // splitting x into  w and d
			int d = x - (w<<nbBitD); // when nbitD=0, this ensures d=0
			int wSign = (w >> (nbBitW-1)) << (nbBitW-1); //getting the sign bit
			w -= 2*wSign; // recovering the signed value

			// Now compute the intervals that these values represent: any number in [realwMin, realwMax) may have been truncated to w.
			double realwMin = 2*radix * ((double)w) / ((double)(1<<nbBitW)); // exact, in [-beta, beta). This bound is reached
			double realwMax = 2*radix * ((double)w+1) / ((double)(1<<nbBitW)); // exact, in [-beta, beta). This bound is never reached
			double epsilon = ( (double) 1)  / (1ULL<<40) ; 
			realwMax -= epsilon; // dirty hack 

			double realdMin = (1-dmsbweight)  +  dmsbweight * ((double)d) / ((double)(1<<nbBitD)); // exact, in [0, 1). This bound is reached
			double realdMax = (1-dmsbweight)  +  dmsbweight *  ((double)d+1) / ((double)(1<<nbBitD)); // exact, in [0, 1). This bound is reached
			realdMax -= epsilon; // dirty hack
			
			realdMin = max(realdMin,dMin); // ignore whatever is on the left of dMin
			realdMax = min(realdMax,dMax); // ignore whatever is on the right of dMax
			int result;
			double uSlope, lSlope, uBound, lBound;
					
			result=1<<30;

			// was: 			for(int k = -alpha; k <= alpha; k++)		{
			// but we want to favor negative numbers when Ri<0
			int kstart, kstop, delta;
			if(wSign>=0) {
				kstart = -alpha;
				kstop = alpha+1;
				delta = +1;
			}
			else {
				kstart = alpha;
				kstop = -alpha-1;
				delta = -1;
			}
			for(int k = kstart; k!= kstop; k+=delta) {
				uSlope = k+rho;
				lSlope = k-rho;

				uBound = uSlope * (uSlope>=0? realdMin: realdMax);	
				lBound = lSlope * (lSlope>=0? realdMax: realdMin);	

				/* the actual test is the first line. 
					 The second and third are there to pad with -alpha and alpha the unreachable parts 
					 of the domain, becaause we don't have don't cares in FloPoCo yet. 
					 Incidentally, it protects us from the Pentium bug */
				//cerr << "k=" << k<< " realwMax=" << realwMax << " <=  uBound="<< uBound << "? :"<< (realwMax <= uBound) << "   realwMin=" << realwMin << ">= lBound="<<lBound << " ?" <<endl;
				if( (realwMax <= uBound && realwMin >= lBound)
						|| (k == alpha && realwMin >= lBound)  
						|| (k == -alpha && realwMax <= uBound) 
						)
					{
						result = k;
						break;
					}
		}
			if(result==1<<30) {
				THROWERROR("digit selection failed for w=" << w << " and d=" << d);
			}
		int result2c = result;
			if(result < 0)	
				result2c+=(1<<wOut); //switch to two's complement

			mpz_class mpzresult;
						
			mpzresult = mpz_class(result2c);

			//cerr << " w="<< w  << "  "<<  realwMin  << "  "<<  realwMax  <<" d=" << d << "  "<<  realdMin  << "  "<<  realdMax << " qi=" << result << "  " << result << endl;
			t.push_back(mpzresult);
		}
		return t;
	}
	

	
	FPDiv::FPDiv(OperatorPtr parentOp, Target* target, int wE, int wF, int srt) :
		Operator(parentOp, target), wE(wE), wF(wF) {

		int i;
		ostringstream name;
		setCopyrightString("Maxime Christ, Florent de Dinechin (2015)");

		srcFileName="FPDiv";
		name<<"FPDiv_"<<wE<<"_"<<wF;
		setNameWithFreqAndUID(name.str());

		if(srt!=42 && srt!=43 && srt!=87){
		THROWERROR("Invalid value for srt: " << srt  );
		}

		if(srt==42) {
			radix=4;
			alpha=2;
		}
		if(srt==43) {
			radix=4;
			alpha=3;
		}
		if(srt==87) {
			radix=8;
			alpha=7;
		}


		addFPInput ("X", wE, wF);
		addFPInput ("Y", wE, wF);
		addFPOutput("R", wE, wF);
			

		vhdl << tab << declare("fX",wF+1) << " <= \"1\" & X(" << wF-1 << " downto 0);" << endl;
		vhdl << tab << declare("fY",wF+1) << " <= \"1\" & Y(" << wF-1 << " downto 0);" << endl;

		vhdl << tab << "-- exponent difference, sign and exception combination computed early, to have fewer bits to pipeline" << endl;

		vhdl << tab << declare("expR0", wE+2) << " <= (\"00\" & X(" << wE+wF-1 << " downto " << wF << ")) - (\"00\" & Y(" << wE+wF-1 << " downto " << wF<< "));" << endl;
		vhdl << tab << declare(getTarget()->lutDelay(), "sR") << " <= X(" << wE+wF << ") xor Y(" << wE+wF<< ");" << endl;
		vhdl << tab << "-- early exception handling " <<endl;
		vhdl << tab << declare("exnXY",4) << " <= X(" << wE+wF+2 << " downto " << wE+wF+1  << ") & Y(" << wE+wF+2 << " downto " << wE+wF+1 << ");" <<endl;
		
		vhdl << tab << "with exnXY select" <<endl
				 << tab << tab << declare(getTarget()->lutDelay(),"exnR0", 2) << " <= " << endl
				 << tab << tab << tab << "\"01\"	 when \"0101\",										-- normal" <<endl
				 << tab << tab << tab << "\"00\"	 when \"0001\" | \"0010\" | \"0110\", -- zero" <<endl
				 << tab << tab << tab << "\"10\"	 when \"0100\" | \"1000\" | \"1001\", -- overflow" <<endl
				 << tab << tab << tab << "\"11\"	 when others;										-- NaN" <<endl;

		int extraBit = 0;
		int nbBitsD, nbBitsW;
		int qSize; // size in bits of Q
		int qMsbToDiscard; // how many known MSB zeroes in Q due to initial alignment to respect R0<rho D
	
		/*****************************RADIX 8 SRT **********************************************/
		if(srt==87)
		{
			// parameter setup
			radix=8;
			alpha=7;
			nbBitsD=2;
			nbBitsW=5;

			extraBit+=2; //Here we'll prescale by 5/4 => 2 right extra bits
			extraBit+=1; //The sticky bit
			extraBit+=1; //The result will be in [1/2, 2[ => 1 more bit (2^0)
			extraBit+=2; //To round correctly the result
			extraBit+=3; //floor() and the bits cut to get a result depending on wF instead of nDigit (cf. last step before normalization)

			nDigit = floor(((double)(wF + extraBit))/3);

			/////////////////////////////////////////////////////////////////////////Prescaling
			//TODO : maybe we can reduce fX and fY
			vhdl << tab << " -- Prescaling" << endl;
			vhdl << tab << "with fY " << range(wF-1, wF-2) << " select" << endl;
			vhdl << tab << tab << declare(getTarget()->adderDelay(wF+3), "prescaledfY", wF+3)
					 << " <= " << endl //potentially *5/4 => we need 2 more bits
					 << tab << tab << tab << "(\"0\" & fY & \"0\") + (fY & \"00\") when \"00\","<<endl /////////[1/2, 5/8[ * 3/2 => [3/4, 15/16[
					 << tab << tab << tab << "(\"00\" & fY) + (fY & \"00\") when \"01\","<<endl ////////////////[5/8, 3/4[ * 5/4 => [25/32, 15/16[
					 << tab << tab << tab << "fY &\"00\" when others;"<<endl; /////////no prescaling

			vhdl << tab << "with fY " << range(wF-1, wF-2) << " select" << endl
					 << tab << tab << declare(getTarget()->adderDelay(wF+4), "prescaledfX", wF+4) << " <= " << endl
					 << tab << tab << tab << "(\"00\" & fX & \"0\") + (\"0\" & fX & \"00\") when \"00\","<<endl  /////////[1/2, 5/[*3/2 => [3/4, 15/16[
					 << tab << tab << tab << "(\"000\" & fX) + (\"0\" & fX & \"00\") when \"01\","<<endl  ////////////////[5/8, 3/4[*5/4 => [25/32, 15/16[
			     << tab << tab << tab << "\"0\" & fX &\"00\" when others;"<<endl; /////////no prescaling

			vhdl << tab << declare(join("w", nDigit-1), wF+6) << " <=  \"00\" & prescaledfX;" << endl; //TODO : review that, maybe MSB 0 to save

			vector<mpz_class> tableContent = selFunctionTable(0.75, 1.0, nbBitsD, nbBitsW, alpha, radix);
			Table* selfunctiontable = new Table(this, target, tableContent,"selFunction7_4", nbBitsD+nbBitsW, 4);

			for(i=nDigit-1; i>=1; i--) {

				REPORT(DEBUG, "Entering iteration " << i);
				// TODO: get rid of all the ostringstream on the model of qi
				string qi =join("q", i);						//actual quotient digit, LUT's output

				ostringstream wi, wim1, seli, wipad, wim1full, wim1fulla;
				wi << "w" << i;						//actual partial remainder
				wim1 << "w" << i-1;					//partial remainder for the next iteration, = left shifted wim1full
				seli << "sel" << i;					//constructed as the wi's first 4 digits and D's first, LUT's input
				wipad << "w" << i << "pad";			//1-left-shifted wi
				wim1full << "w" << i-1 << "full";	//partial remainder after this iteration, = wi+qi*D
				wim1fulla << "w" << i-1 << "fulla";	//partial remainder after this iteration, = wi+qi*D
				string tInstance = "SelFunctionTable" + to_string(i);

				vhdl << tab << declare(seli.str(),7) << " <= " << wi.str() << range( wF+5, wF+1)<<" & prescaledfY"<< range(wF, wF-1) <<";" << endl;

				newSharedInstance(selfunctiontable , tInstance, "X=>"+seli.str(), "Y=>"+ qi);
				// REPORT(DEBUG, "After table instance " << i);
				
				vhdl << tab << declare(wipad.str(), wF+7) << " <= " << wi.str() << " & '0';" << endl;

				// qui has a fanout of(wF+7), which we add to both its uses 
				vhdl << tab << "with " << qi << range(1,0) << " select " << endl
						 << tab << declare(getTarget()->adderDelay(wF+7)+getTarget()->fanoutDelay(2*(wF+7)), wim1fulla.str(), wF+7) << " <= " << endl
						 << tab << tab << wipad.str() << " - (\"0000\" & prescaledfY)			when \"01\"," << endl
						 << tab << tab << wipad.str() << " + (\"0000\" & prescaledfY)			when \"11\"," << endl
						 << tab << tab << wipad.str() << " + (\"000\" & prescaledfY & \"0\")		when \"10\"," << endl
						 << tab << tab << wipad.str() << "							when others;" << endl;

				REPORT(DEBUG, "After with 1 ");

#if 0 // Splitting the following logic into two levels gives better synthesis results...
				vhdl << tab << "with " << qi << range(3,1) << " select " << endl;
				vhdl << tab << declare(getTarget()->adderDelay(wF+7)+getTarget()->fanoutDelay(2*(wF+7)), wim1full.str(), wF+7) << " <= " << endl;
				vhdl << tab << tab << wim1fulla.str() << " - (\"00\" & prescaledfY & \"00\")			when \"001\" | \"010\"," << endl;
				vhdl << tab << tab << wim1fulla.str() << " - (\"0\" & prescaledfY & \"000\")			when \"011\"," << endl;
				vhdl << tab << tab << wim1fulla.str() << " + (\"00\" & prescaledfY & \"00\")			when \"110\" | \"101\"," << endl;
				vhdl << tab << tab << wim1fulla.str() << " + (\"0\" & prescaledfY & \"000\")			when \"100\"," << endl;
				vhdl << tab << tab << wim1fulla.str() << " 			   		  when others;" << endl;
#else
				string fYdec =join("fYdec", i-1);	//
				vhdl << tab << "with " << qi << range(3,1) << " select " << endl
						 << tab << declare(getTarget()->lutDelay()+getTarget()->fanoutDelay(2*(wF+7)), fYdec, wF+7) << " <= " << endl
						 << tab << tab << "(\"00\" & prescaledfY & \"00\")			when \"001\" | \"010\" | \"110\"| \"101\"," << endl
						 << tab << tab << "(\"0\" & prescaledfY & \"000\")			when \"011\"| \"100\"," << endl
						 << tab << tab << rangeAssign(wF+6,0,"'0'") << "when others;" << endl;
				
				vhdl << tab << "with " << qi << of(3) << " select" << endl; // Remark here: seli(6)==qi(3) but it we get better results using the latter.
				vhdl << tab << declare(getTarget()->adderDelay(wF+7), wim1full.str(), wF+7) << " <= " << endl;
				vhdl << tab << tab << wim1fulla.str() << " - " << fYdec << "			when '0'," << endl;
				vhdl << tab << tab << wim1fulla.str() << " + " << fYdec << "			when others;" << endl;

#endif
				vhdl << tab << declare(wim1.str(),wF+6) << " <= " << wim1full.str()<<range(wF+3,0)<<" & \"00\";" << endl;
			}


			vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wF+3), "q0",4) << "(3 downto 0) <= \"0000\" when  w0 = (" << wF+5 << " downto 0 => '0')" << endl;
			vhdl << tab << "             else w0(" << wF+5 << ") & \"010\";" << endl;

			for(i=nDigit-1; i>=1; i--) {
				string qi =join("q", i);
				ostringstream qPi, qMi;
				qPi << "qP" << i;
				qMi << "qM" << i;
				vhdl << tab << declare(qPi.str(), 3) <<" <=      " << qi << "(2 downto 0);" << endl;
				vhdl << tab << declare(qMi.str(), 3)<<" <=      " << qi << "(3) & \"00\";" << endl;
			}

			vhdl << tab << declare("qP0", 3) << " <= q0(2 downto 0);" << endl;
			vhdl << tab << declare("qM0", 3) << " <= q0(3)  & \"00\";" << endl;

			vhdl << tab << declare("qP", 3*nDigit) << " <= qP" << nDigit-1;
			for (i=nDigit-2; i>=0; i--)
				vhdl << " & qP" << i;
			vhdl << ";" << endl;

			vhdl << tab << declare("qM", 3*nDigit) << " <= qM" << nDigit-1 << range(1,0);
			for (i=nDigit-2; i>=0; i--)
				vhdl << " & qM" << i;
			vhdl << " & \"0\";" << endl;

			vhdl << tab << declare(getTarget()->adderDelay(3*nDigit), "quotient", 3*nDigit) << " <= qP - qM;" << endl;

			//The last +3 in computing nDigit is for this part
			// Here we get rid of the leading bit, which is a known zero, we keep
			// 1 bit for the norm, 1+wF fraction bits, 1 round bit, and one sticky bit out of the LSBs			
			qSize=3*nDigit; // where nDigit is  floor((wF + extraBit)/3);
			qMsbToDiscard=1;
			
			}



		

		else 
			/*****************************RADIX 4 SRT **********************************************/
			// For alpha=3, we need 2 LUTs/bit in the recurrence but the table is small 
			// For alpha=2 with no prescaling we have a large sel table (10 bits in) but 1 LUT/bit recurrence.
			//             with prescaling the tables get much smaller (8 input bits)
		{
			vector<mpz_class> tableContent;
			Table* selfunctiontable;
			// -------- Parameter set up -----------------

			prescaling=0 ; // failed to work so far

			if(alpha==3) {
				nbBitsD=1;
				nbBitsW=4;
				extraBit=6; //
				tableContent = selFunctionTable(0.5, 1.0, nbBitsD, nbBitsW, alpha, radix);
			}
			
			else if(alpha==2){
				if(prescaling==1) {
					nbBitsD=2;
					nbBitsW=6;
					tableContent = selFunctionTable(0.75, 9./8., nbBitsD, nbBitsW, alpha, radix);
				}
				else if(prescaling==0) { // no prescaling
					nbBitsD=3;
					nbBitsW=6;
					tableContent = selFunctionTable(0.5, 1.0, nbBitsD, nbBitsW, alpha, radix);
				}
				else THROWERROR("prescaling="<< prescaling << " is not an option");

				extraBit=7; // one more than for alpha=3 due to the initial alignment to respect rhoD
			}
			else THROWERROR("alpha="<< alpha << " is not an option");
			
			nDigit = (wF+extraBit) >> 1;

			int dSize, subSize;
			if(prescaling==0) {
				dSize=wF+1;
				vhdl << tab << declare("D", dSize) << " <= fY ;"<< endl;
				vhdl << tab << declare("psX", dSize+1) << " <= \"0\" & fX ;"<< endl;				
			}
			else if(prescaling==1) {
				vhdl << tab << " -- prescaling " << endl;
				dSize=wF+2;
				vhdl << tab << declare("D", dSize) << " <=  (fY & \"0\") when fY" << of(wF-1) << " = '1' -- D when D was in [1.5,2)"<< endl
						 << tab << tab << " else (\"0\" & fY) + (fY & \"0\") ; -- 3/2*D when D was in [1, 1.5)"  << endl;
				vhdl << tab << " -- Now D is in [3/2, 9/4) and one bit wider " << endl;
				
				vhdl << tab << declare("psX", dSize+1) << " <=  (\"0\" & fX & \"0\") when fY" << of(wF-1) << " = '1'"<< endl
						 << tab << tab << " else (\"00\" & fX) + (\"0\" & fX & \"0\") ; "<< endl;
			}
			else THROWERROR("prescaling="<< prescaling << " is not an option");
			
			if(alpha==3) {
			  vhdl << tab << " -- compute 3D" << endl;
			  vhdl << tab << declare(getTarget()->adderDelay(dSize+1), "Dx3",dSize+1)
					 << " <= (\"0\" & D) + (D & \"0\");" << endl;
      }

			subSize=dSize+3; // Size of the subtraction in the main iteration
			

#if 0 // experiments with prescaling
			nbBitsD=0;
			nbBitsW=5;
			tableContent = selFunctionTable(1.-1./64., 1.+1./8., nbBitsD, nbBitsW, alpha, radix);
#endif
#if 0 // experiments with prescaling
			// 42 with prescaling by 3/2 when D<42/64 allows to save one table input bit. After synthesis of just this experiment the
			// LUTs are reduced from 564 to 533. Since the prescaling itself will cost 2 x 27 LUTs, it is a loss. 
			nbBitsD=3;
			nbBitsW=5;
			tableContent = selFunctionTable(42./64, 1, nbBitsD, nbBitsW, alpha, radix);
#endif


			selfunctiontable = new Table(this, target, tableContent,"selFunction", nbBitsD+nbBitsW, 3);
			selfunctiontable->setShared();

			////////////////////// Main SRT loop, unrolled ///////////////////////

			for(i=nDigit-1; i>=1; i--) {

				string qi =join( "q", i);						//current quotient digit, LUT's output
				string wi = join("betaw", i);						// current partial remainder
				string wifull =join("w", i);						// current partial remainder
				string seli = join("sel", i);					//constructed as the wi's first 4 digits and D's first, LUT's input
				string qiTimesD = join("absq", i)+"D";		//qi*Y
				string wim1full = join("w", i-1);	//partial remainder after this iteration, = wi+qi*D

				/*
						Detailed algorithm for alpha=3 :
					*	building seli for the selection function
						seli = wi (25 downto 22) & fY(22), i.e. the first 4 digits of the remainder and the first useful digit of the divisor
					*	deducing the value of qi out of seli
					*	deducing the value of qi*D out of qi
						qi is bounded to [-3,3], so you can compute qi*D in 1 addition
					*	D : 24 bits, qi : 3 bits => qi*D : potentially 27 bits
rox P						or wi is 26 bits long
						=>computing wipad = wi left shifted (27 bits)
					*	computing the remainder of this step
						wi-1full = wi-qi*D
					*	left shifting wi-1full to obtain wi-1, next partial remainder to work on
				*/
				if(i==nDigit-1){
					if(alpha==2) {
						vhdl << tab << declare(wi, subSize) << " <=  \"00\" & psX;" << endl;
					} 
					else { // alpha=3
						vhdl << tab << declare(wi, subSize) << " <=  \"0\" & psX & \"0\";" << endl;
					}
				}
				else {
					//					
					vhdl << tab << declare(wi,subSize) << " <= " << wifull<<range(subSize-3,0)<<" & \"00\"; -- multiplication by the radix" << endl;
				}

				if(prescaling==1) { // now that D may exceed 1, we need to consider its top bit as well
					vhdl << tab << declare(seli, nbBitsW+nbBitsD) << " <= " << wi << range(subSize-1, subSize-nbBitsW) << " & D" << range(dSize-1,dSize-nbBitsD)  << ";" << endl;
				}
				else {
					vhdl << tab << declare(seli, nbBitsW+nbBitsD) << " <= " << wi << range(subSize-1, subSize-nbBitsW) << " & D" << range(dSize-2,dSize-1-nbBitsD)  << ";" << endl;				 
				}
				
				newSharedInstance(selfunctiontable , "SelFunctionTable" + to_string(i), "X=>"+seli, "Y=>"+ qi);
				vhdl << endl;

				if(alpha==3) {
					// Two options for here. More experiments are needed, the best is probably target-dependent 
#if 1  // The following leads to higher frequency and higher resource usage: 
					// For (8,23) on Virtex6 with ISE this gives 466Mhz, 1083 regs+ 1029 LUTs 
					vhdl << tab << "with " << qi << " select" << endl;
					vhdl << tab << tab << declare(getTarget()->fanoutDelay(subSize) + getTarget()->adderDelay(subSize), qiTimesD,subSize)
							 << " <= "<< endl 
							 << tab << tab << tab << "\"000\" & D  		   when \"001\" | \"111\"," << endl
							 << tab << tab << tab << "\"00\" & D & \"0\"	 when \"010\" | \"110\"," << endl
							 << tab << tab << tab << "\"00\" & Dx3    	   when \"011\" | \"101\"," << endl
							 << tab << tab << tab << "(" << subSize-1 << " downto 0 => '0')	when others;" << endl<< endl;
#else // Recompute 3Y locally to save the registers: the LUT is used anyway (wrong! on Virtex6 ISE finds a MUX) 
					// For (8,23) on Virtex6 with ISE this gives 345Mhz, 856 regs+ 1051 LUTs 
					// Note that this option probably scales worse if we pipeline this addition  
					vhdl << tab << "with " << qi << " select" << endl
							 << tab << tab << declare(getTarget()->fanoutDelay(subSize) + getTarget()->lutDelay(), join("addendA",i),subSize)
							 << " <= " << endl 
							 << tab << tab << tab << "\"000\" & D              when \"001\" | \"111\" | \"011\" | \"101\"," << endl
							 << tab << tab << tab << "(" << subSize-1 << " downto 0 => '0')  when others;" << endl;

					vhdl << tab << "with " << qi << " select" << endl
							 << tab << tab << declare(join("addendB",i),subSize) << " <= "<< endl 
							 << tab << tab << tab << "\"00\" & D & \"0\"       when \"010\" | \"110\"| \"011\" | \"101\"," << endl
							 << tab << tab << tab << "(" << subSize-1 << " downto 0 => '0')  when others;" << endl;
					
					vhdl << tab << tab << declare(getTarget()->adderDelay(subSize), qiTimesD,subSize)
							 << " <= " << join("addendA",i) << " + " << join("addendB",i) << ";"<< endl << endl;
#endif				

					vhdl << tab << "with " << qi << "(2) select" << endl;
					vhdl << tab << declare(getTarget()->adderDelay(subSize), wim1full, subSize)
							 << "<= " << wi << " - " << qiTimesD << " when '0'," << endl
							 << tab << "      " << wi << " + " << qiTimesD << " when others;" << endl << endl;

				} // end if alpha=3

				
				else { // alpha=2
					vhdl << tab << "with " << qi << " select" << endl;
					// no delay for qiTimesD, it should be merged in the following addition
					vhdl << tab << tab << declare(qiTimesD,subSize) << " <= "<< endl 
							 << tab << tab << tab << "\"000\" & D						 when \"001\" | \"111\", -- mult by 1" << endl
							 << tab << tab << tab << "\"00\" & D & \"0\"			   when \"010\" | \"110\", -- mult by 2" << endl
							 << tab << tab << tab << "(" << subSize-1 << " downto 0 => '0')	 when others;        -- mult by 0" << endl << endl;
					
					//				vhdl << tab << declare(wi, subSize) << " <= " << wi << " & \"0\";" << endl;

					vhdl << tab << "with " << qi << "(2) select" << endl
							 << tab << declare(getTarget()->adderDelay(subSize), wim1full, subSize)
							 << "<= " << wi << " - " << qiTimesD << " when '0'," << endl
							 << tab << "      " << wi << " + " << qiTimesD << " when others;" << endl << endl;


				}
			} // end loop

			vhdl << tab << declare("wfinal", wF+2) << " <= w0" <<range(wF+1,0) << ";" << endl;
			vhdl << tab << declare("qM0") << " <= wfinal" << of(wF+1) << "; -- rounding bit is the sign of the remainder" << endl;

			for(i=nDigit-1; i>=1; i--) {
				ostringstream qPi, qMi;
				string qi = join("q",i);
				qPi << "qP" << i;
				qMi << "qM" << i;
				vhdl << tab << declare(qPi.str(), 2) <<" <=      " << qi << "(1 downto 0);" << endl;
				vhdl << tab << declare(qMi.str(), 2)<<" <=      " << qi << "(2) & \"0\";" << endl;
			}

			vhdl << tab << declare("qP", 2*nDigit-2) << " <= qP" << nDigit-1;
			for (i=nDigit-2; i>=1; i--)
				vhdl << " & qP" << i;
			vhdl << ";" << endl;

			vhdl << tab << declare("qM", 2*nDigit-2) << " <= qM" << nDigit-1 << "(0)";
			for (i=nDigit-2; i>=1; i--)
				vhdl << " & qM" << i;
			vhdl << " & qM0;" << endl;


			// TODO an IntAdder here
			vhdl << tab << declare(getTarget()->adderDelay(2*nDigit), "quotient", 2*nDigit-2) << " <= qP - qM;" << endl;

			// preparing the extraction of a mantissa from q
			qSize=2*nDigit-2; // where nDigit is  floor((wF + extraBit)/2)
			if(alpha==2)
				qMsbToDiscard=1;// due to initial alignment to respect R0<rho D with rho=2/3
			else
				qMsbToDiscard=0; 
		}
		
		vhdl << tab << "-- We need a mR in (0, -wf-2) format: 1+wF fraction bits, 1 round bit, and 1 guard bit for the normalisation," << endl;
		vhdl << tab << "-- quotient is the truncation of the exact quotient to at least 2^(-wF-2) bits" << endl;
		vhdl << tab << "-- now discarding its possible known MSB zeroes, and dropping the possible extra LSB bit (due to radix 4) " << endl;
		int lsbSize = qSize-qMsbToDiscard-(wF+3);
		vhdl << tab << declare(getTarget()->lutDelay(), "mR", wF+3) << " <= quotient(" << qSize-1-qMsbToDiscard << " downto "<< lsbSize <<"); " << endl;

		vhdl << tab << "-- normalisation" << endl;		
		vhdl << tab << declare(getTarget()->lutDelay(), "fRnorm", wF+1) << " <=    mR(" << wF+1 << " downto 1)  when mR" << of(wF+2) << "= '1'" << endl;
		vhdl << tab << "        else mR(" << wF << " downto 0);  -- now fRnorm is a (-1, -wF-1) fraction" << endl;
		
		vhdl << tab << declare(getTarget()->lutDelay(), "round") << " <= fRnorm(0); " << endl;

		vhdl << tab << declare(getTarget()->adderDelay(wE+2), "expR1", wE+2) << " <= expR0"
						 << " + (\"000\" & (" << wE-2 << " downto 1 => '1') & mR(" << wF+2 << ")); -- add back bias" << endl;

		vhdl << tab << "-- final rounding" <<endl;
		vhdl << tab <<  declare("expfrac", wE+wF+2) << " <= "
				 << "expR1 & fRnorm(" <<wF << " downto " << 1 << ") ;" << endl;
		
		vhdl << tab << declare("expfracR", wE+wF+2) << " <= "
				 << "expfrac + ((" << wE+wF+1 << " downto 1 => '0') & round);" << endl;

		vhdl << tab <<  declare(getTarget()->lutDelay(), "exnR", 2)
				 << " <=      \"00\"  when expfracR(" << wE+wF+1 << ") = '1'   -- underflow" <<endl;
		vhdl << tab << "        else \"10\"  when  expfracR(" << wE+wF+1 << " downto " << wE+wF << ") =  \"01\" -- overflow" <<endl;
		vhdl << tab << "        else \"01\";      -- 00, normal case" <<endl;

		// Normalisation is common to both radix 4 and radix 8
		vhdl << tab << "with exnR0 select" <<endl;
		vhdl << tab << tab << declare(getTarget()->lutDelay(), "exnRfinal", 2) << " <= " <<endl;
		//vhdl << tab << tab << tab << "exnR   when \"01\", -- normal" <<endl;
		vhdl << tab << tab << tab << "exnR   when \"01\", -- normal" <<endl;
		vhdl << tab << tab << tab << "exnR0  when others;" <<endl;
		vhdl << tab << "R <= exnRfinal & sR & "
				 << "expfracR(" << wE+wF-1 << " downto 0);" <<endl;
	}

	FPDiv::~FPDiv() {
	}



	// Various functions that used to be in NbBitsMin
	void SRTDivNbBitsMin::computeNbBit (int radix, int digitSet)
	{

		double ro = (double)digitSet/(radix-1);
		cout<<"Rendundancy coefficiant is rho="<<ro<<endl<<endl;

		// eq. 5.82 p. 293 of Digital Arithmetic
		double exactDeltaMin = 1-log2((2*ro - 1)/(2*(digitSet-1)));

		int deltaMinMinus = exactDeltaMin;
		int deltaMinPlus = deltaMinMinus + 1;

		double nbBitMinus = -log2((2*ro-1)/2 - (digitSet-ro)*pow(2, 0-deltaMinMinus))+deltaMinMinus+log2(radix);
		int nbBitM = ceil(nbBitMinus);

		double nbBitPlus = -log2((2*ro-1)/2 - (digitSet-ro)*pow(2, 0-deltaMinPlus))+deltaMinPlus+log2(radix);
		int nbBitP = ceil(nbBitPlus);

		cout<<"There are 2 possibilities :"<<endl;
		cout<<"-Delta = "<<deltaMinPlus<<", nbBits = "<<nbBitP<<endl;
		cout<<"-Delta = "<<deltaMinMinus<<", nbBits = "<<nbBitM<<endl<<endl;

		int delta = (nbBitP-nbBitM < 0 ? deltaMinPlus : deltaMinMinus);
		int nbBit = (nbBitP-nbBitM < 0 ? nbBitP : nbBitM);

		cout<<"Therefore, the min is for Delta = "<<delta<<endl;
		cout<<"You'll need "<<nbBit<<" bits : ";
		cout<<delta-1<<" bits for D and ";
		cout<<nbBit-delta+1;
		cout<<" bits for the partial remainder"<<endl;

		cout<<"Checking for better configurations..."<<endl<<endl;

		//Optimized parameters for a digit set [-7,7] with a radix 8
		//delta = 3;
		//nbBit = 7;
		//radix = 8;
		//digitSet = 7;

		if(checkDistrib(delta, nbBit-delta+1, radix, digitSet))
		{
			int gain = 0;
			cout<<"Optimization found!"<<endl;
			while(checkDistrib(delta, nbBit-delta, radix, digitSet))
			{
				nbBit--;
				gain++;
			}
			while(checkDistrib(delta-1, nbBit-delta+1, radix, digitSet))
			{
				delta--;
				nbBit--;
				gain++;
			}

			cout<<"You lost "<<gain<<" extra bits"<<endl;
		}
		else
		{
			cout<<"No better configuration found"<<endl;
		}


		cout<<"Final values are : ";
		cout<<delta-1<<" bits for D and ";
		cout<<nbBit-delta+1;
		cout<<" bits for the partial remainder"<<endl;

		plotPDDiagram(delta-1, nbBit-delta+1, radix, digitSet);
		//plotPDDiagram(3, 5, 8, 7);

		cout<<"An implementation of this configuration is approximately "<<estimateCost(nbBit, radix, digitSet)<<" times larger than the actual implementation"<<endl;

	}

	//Produces the P-D Diagram corresponding to the previous analysis in svg format
	void SRTDivNbBitsMin::plotPDDiagram(int delta, int t, int radix, int digitSet)
	{
		double ro = (double)digitSet/((double)radix-1);

		ofstream svg("PDDiagram.svg", ios::out|ios::trunc);
		const int width = 1024;
		const int height = 750;
		int maxH = (-digitSet-ro)*(height/2)/(digitSet+1); //The height doesn't always match with the extrema of Uk and Lk,
		int minH = (digitSet+ro)*(height/2)/(digitSet+1);  //so those value are used to draw a better graph

		//Doc info
		svg<<"<?xml version=\"1.0\"?>"<<endl;
		svg<<"<!DOCTYPE  svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"<<endl;
		svg<<"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" >"<<endl;
		svg<<"<text y=\"15\">Les cases incluent leurs bords bas et gauche</text>"<<endl;
		svg<<"<text y=\"40\">		radix"<<radix<<", [-"<<digitSet<<","<<digitSet<<"]"<<"</text>"<<endl;

		//Lk-Uk
		for(int i = -digitSet; i <= digitSet; i++)
		{
			//Lk
			svg<<"	<line x1=\""<<50<<"\" y1=\""<<50+height/2<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2-(i-ro)*(height/2)/(digitSet+1);
			svg<<"\" style=\"stroke:rgb(0,200,200);stroke-width:2\"/>"<<endl;
			//Uk
			svg<<"	<line x1=\""<<50<<"\" y1=\""<<50+height/2<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2-(i+ro)*(height/2)/(digitSet+1);
			svg<<"\" style=\"stroke:rgb(0,200,200);stroke-width:2\"/>"<<endl;
			//Overlapping
			svg<<"	<line x1=\""<<50+width+(i+digitSet+1)*10<<"\" y1=\""<<50+height/2-(i-ro)*(height/2)/(digitSet+1);
			svg<<"\" x2=\""<<50+width+(i+digitSet+1)*10<<"\" y2=\""<<50+height/2-(i+ro)*(height/2)/(digitSet+1);
			svg<<"\" style=\"stroke:rgb(0,0,200);stroke-width:2\"/>"<<endl;
		}


		//Horizontal axis subdiv
		int nbSubdiv = pow(2, delta);
		int pitch = width/nbSubdiv;

		for(int i = 0; i <= nbSubdiv; i++)
		{
			svg<<"	<line x1=\""<<50+width/2+i*pitch/2<<"\" y1=\""<<50+height/2+maxH;//50+3*height/4-((delta%2==0?i:i-1)*(-digitSet-ro)*(height/2)/(digitSet+1)/2/nbSubdiv);
			svg<<"\" x2=\""<<50+width/2+i*pitch/2<<"\" y2=\""<<50+height/2+minH;// /4-((delta%2==0?i:i-1)*(digitSet+ro)*(height/2)/(digitSet+1)/2/nbSubdiv);
			svg<<"\" style=\"stroke:rgb(100,100,100);stroke-width:2\"/>"<<endl;
		}


		//Vertical axis subdiv
		nbSubdiv = pow(2, t);

		for(int i = -nbSubdiv/2; i <= nbSubdiv/2; i++)
		{
			svg<<"	<line x1=\""<<50+width/2<<"\" y1=\""<<50+height/2+i*2*minH/nbSubdiv;
			svg<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2+i*2*minH/nbSubdiv;
			svg<<"\" style=\"stroke:rgb(100,100,100);stroke-width:1\"/>"<<endl;
		}


		//Axis
		svg<<"	<line x1=\""<<50<<"\" y1=\""<<50<<"\" x2=\""<<50<<"\" y2=\""<<50+height<<"\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>"<<endl;//vertical 0
		svg<<"	<line x1=\""<<50<<"\" y1=\""<<50+height/2<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2<<"\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>"<<endl;//horizontal

		svg<<"	<line x1=\""<<50+width<<"\" y1=\""<<50+height/2+minH;
		svg<<"\" x2=\""<<50+width<<"\" y2=\""<<50+height/2+maxH<<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>"<<endl;//vertical 1

		svg<<"	<line x1=\""<<50+width/2<<"\" y1=\""<<50+height/2+minH/2;
		svg<<"\" x2=\""<<50+width/2<<"\" y2=\""<<50+height/2+maxH/2<<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>"<<endl;//vertical 0.5

		//Prescaling (cas particulier base 8, digitSet 7)
	//	svg<<"	<line x1=\""<<50+width/2<<"\" y1=\""<<50+height/2+minH/2;
	//	svg<<"\" x2=\""<<50+width/2<<"\" y2=\""<<50+height/2+maxH/2<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 0.5
	//
	//	svg<<"	<line x1=\""<<50+5*width/8-3<<"\" y1=\""<<50+height/2+5*minH/8;
	//	svg<<"\" x2=\""<<50+5*width/8-3<<"\" y2=\""<<50+height/2+5*maxH/8<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 0.625-
	//
	//	svg<<"	<line x1=\""<<50+5*width/8+3<<"\" y1=\""<<50+height/2+5*minH/8;
	//	svg<<"\" x2=\""<<50+5*width/8+3<<"\" y2=\""<<50+height/2+5*maxH/8<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 0.625+
	//
	//	svg<<"	<line x1=\""<<50+3*width/4+3<<"\" y1=\""<<50+height/2+3*minH/4;
	//	svg<<"\" x2=\""<<50+3*width/4+3<<"\" y2=\""<<50+height/2+3*maxH/4<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 0.75+
	//
	//	svg<<"	<line x1=\""<<50+3*width/4-3<<"\" y1=\""<<50+height/2+3*minH/4;
	//	svg<<"\" x2=\""<<50+3*width/4-3<<"\" y2=\""<<50+height/2+3*maxH/4<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 0.75-
	//
	//	svg<<"	<line x1=\""<<50+25*width/32<<"\" y1=\""<<50+height/2+25*minH/32;
	//	svg<<"\" x2=\""<<50+25*width/32<<"\" y2=\""<<50+height/2+25*maxH/32<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 25/32
	//
	//	svg<<"	<line x1=\""<<50+15*width/16-3<<"\" y1=\""<<50+height/2+15*minH/16;
	//	svg<<"\" x2=\""<<50+15*width/16-3<<"\" y2=\""<<50+height/2+15*maxH/16<<"\" style=\"stroke:rgb(0,100,0);stroke-width:5\"/>"<<endl;//vertical 15/16-
	//
	//	svg<<"	<line x1=\""<<50+15*width/16+3<<"\" y1=\""<<50+height/2+15*minH/16;
	//	svg<<"\" x2=\""<<50+15*width/16+3<<"\" y2=\""<<50+height/2+15*maxH/16<<"\" style=\"stroke:rgb(0,200,0);stroke-width:5\"/>"<<endl;//vertical 15/16+


		//Case pb, cas particulier base 8 digitSet 7
	//	svg<<"	<rect x=\""<<50+5*width/8<<"\" y=\""<<50+3*minH/2<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+5*width/8<<"\" y=\""<<50+minH/2-minH/16<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+minH/2<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+3*minH/2-minH/16<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+3*minH/2-minH/8<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;
	//	svg<<"	<rect x=\""<<50+9*width/16<<"\" y=\""<<50+9*minH/16+1<<"\" width=\""<<width/16<<"\" height=\""<<minH/16<<"\" style=\"fill:rgb(255,0,0);fill-opacity:0.3;\" />"<<endl;

		//Fin
		svg<<"</svg>"<<endl;
	}


	bool SRTDivNbBitsMin::checkDistrib(int delta, int t, int radix, int digitSet)
	{
		double ro = digitSet/(radix-1);
		double wMax = U(digitSet, ro, 1);

		for(int k = -digitSet+1; k <= digitSet; k++)
		{

			//the position of the left intersection with Lk
			double leftL = L(k, ro, 0.5);
			double leftCeiledL = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*leftL);
			int leftLineL = (wMax-leftCeiledL)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram
			int leftCol = 0;
			bool leftCornerL = (leftL-leftCeiledL == 0);
			vector<int> crossedBoxes;

			//the position of the left intersection with Uk-1
			double leftU = U(k, ro, 0.5);
			double leftCeiledU = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*leftU);
			int leftLineU = (wMax-leftCeiledU)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram
			bool leftCornerU = (leftU-leftCeiledU == 0);

			int index;

			for(double i = 0.5; i < 1 ; i += pow(2, -delta))
			{
				//the position of the right intersection with Lk
				double rightL = L(k, ro, i+pow(2, -delta));//value of Lk for d=i
				double rightCeiledL = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*rightL);
				int rightLineL = (wMax-rightCeiledL)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram
				int rightCol = (i-0.5+pow(2, -delta))*pow(2, delta-1)/0.5;

				//the position of the right intersection with Uk-1
				double rightU = U(k-1, ro, i+pow(2, -delta));//value of Uk-1 for d=i
				double rightCeiledU = pow(2,-(t-log2(radix)-1))*ceil(pow(2, (t-log2(radix)-1))*rightU);
				int rightLineU = (wMax-rightCeiledU)*pow(2,t)/(2*wMax);//line number in the grid from P-D diagram

				bool rightCornerL = (rightL-rightCeiledL == 0);
				bool rightCornerU = (rightU-rightCeiledU == 0);

				if(leftCornerL && k-ro>0)
				{						//positive coeff =>the bottom-left corner, the box crossed is above the computed-one
										//       _______________
										//	    | /				|   <====== crossed box
					leftLineL--;		//		|/______________|
										//		/               |   <====== computed crossed box
										//	   /|_______________|
				}
				if(rightCornerL && k-ro<0)
				{
					rightLineL--;
				}


				if(leftCornerU && k-ro>0)
				{
					leftLineU--;
				}
				if(rightCornerU && k-ro<0)
				{
					rightLineU--;
				}




				for(int c = (k-ro<0 ? leftLineL : rightLineL); c <= (k-ro<0 ? rightLineL : leftLineL); c++)
				{
					index = leftCol + pow(2,delta-1)*c;
					crossedBoxes.push_back(index);
				}

				for(int c = (k-ro<0 ? leftLineU : rightLineU); c <= (k-ro<0 ? rightLineU : leftLineU); c++)
				{
					index = leftCol + pow(2,delta-1)*c;
					for(vector<int>::iterator it = crossedBoxes.begin(); it != crossedBoxes.end(); ++it)
					{
						if (index == *it)
						{
							return false;
						}
					}
				}



				rightLineL+=(rightCornerL && k-ro<0 ? 1 : 0);
				rightLineU+=(rightCornerU && k-ro<0 ? 1 : 0);

				leftL = rightL;
				leftCeiledL = rightCeiledL;
				leftLineL = rightLineL;
				leftU = rightU;
				leftCeiledU = rightCeiledU;
				leftLineU = rightLineU;
				leftCol = rightCol;

				crossedBoxes.clear();
			}

		}

		return true;
	}

	double SRTDivNbBitsMin::L(int k, double ro, double d)
	{
		return (k-ro)*d;
	}

	double SRTDivNbBitsMin::U(int k, double ro, double d)
	{
		return (k+ro)*d;
	}

	double SRTDivNbBitsMin::estimateCost(int nbBit, int radix, int digitSet)
	{
		int nbIter = (23+6)>>1;
		double flopocoCost = nbIter*(2*27 + 3);
		cout<<"FlopocoCost = "<<flopocoCost<<endl;
		double res = (nbIter*2/log2(radix))*((ceil(log2(digitSet))+1)*ceil(pow(2, nbBit-6))+27*(1+ceil(pow(2, log2(digitSet)-7))));
		cout<<"Cost of this configuration = "<<res<<endl;
		return res/flopocoCost;
	}
	
	void FPDiv::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		/* Compute correct value */
		FPNumber fpx(wE, wF), fpy(wE, wF);
		fpx = svX;
		fpy = svY;
		mpfr_t x, y, r;
		mpfr_init2(x, 1+wF);
		mpfr_init2(y, 1+wF);
		mpfr_init2(r, 1+wF);
		fpx.getMPFR(x); // fake 0
		fpy.getMPFR(y);
		mpfr_div(r, x, y, GMP_RNDN);
		FPNumber  fpr(wE, wF, r);

		/* Set outputs */
		mpz_class svR = fpr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		mpfr_clears(x, y, r, NULL);
	}



	void FPDiv::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		// Regression tests
		tc = new TestCase(this);
		tc->addFPInput("X", 1.0);
		tc->addFPInput("Y", FPNumber::plusDirtyZero);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);
		tc->addFPInput("X", FPNumber::minusDirtyZero);
		tc->addFPInput("Y", FPNumber::plusInfty);
		emulate(tc);
		tcl->add(tc);


	}

	OperatorPtr FPDiv::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int wE;
		UserInterface::parseStrictlyPositiveInt(args, "wE", &wE);
		int wF;
		UserInterface::parseStrictlyPositiveInt(args, "wF", &wF);
		int srt;
		UserInterface::parsePositiveInt(args, "srt", &srt);
		return new FPDiv(parentOp, target, wE, wF, srt);
	}

	TestList FPDiv::unitTest(int index)
	{
		// the static list of mandatory tests
		TestList testStateList;
		vector<pair<string,string>> paramList;
		
		if(index==-1) 
		{ // The unit tests

			for(int wF=5; wF<53; wF+=1) // test various input widths
			{
					int wE = 6+(wF/10);
					while(wE>wF)
					{
						wE -= 2;
					}
					paramList.push_back(make_pair("wF",to_string(wF)));
					paramList.push_back(make_pair("wE",to_string(wE)));
					paramList.push_back(make_pair("srt","42"));
					testStateList.push_back(paramList);
					paramList.clear();
					paramList.push_back(make_pair("wF",to_string(wF)));
					paramList.push_back(make_pair("wE",to_string(wE)));
					paramList.push_back(make_pair("srt","43"));
					testStateList.push_back(paramList);
					paramList.clear();
					paramList.push_back(make_pair("wF",to_string(wF)));
					paramList.push_back(make_pair("wE",to_string(wE)));
					paramList.push_back(make_pair("srt","87"));
					testStateList.push_back(paramList);
					paramList.clear();
			}
		}
		else     
		{
				// finite number of random test computed out of index
		}	
		return testStateList;
	}


	
	void FPDiv::registerFactory(){
		UserInterface::add("FPDiv", // name
											 "A correctly rounded floating-point division.",
											 "BasicFloatingPoint", // categories
											 "http://www.cs.ucla.edu/digital_arithmetic/files/ch5.pdf",
											 "wE(int): exponent size in bits; \
wF(int): mantissa size in bits;\
srt(int)=42: Can be 42, 43 or 87 so far. Default 42 means radix 4 with digits between -2 and 2. Other choices may have a better area/speed trade-offs",
											"The algorithm used here is the division by digit recurrence (SRT). In radix 4, we use a maximally redundant digit set. In radix 8, we use split-digits in [-7,7], and a bit of prescaling.",
											 FPDiv::parseArguments,
											 FPDiv::unitTest

											 ) ;

	}



	
	OperatorPtr SRTDivNbBitsMin::NbBitsMinParseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int radix, digitSet;
		UserInterface::parseStrictlyPositiveInt(args, "radix", &radix);
		UserInterface::parseStrictlyPositiveInt(args, "alpha", &digitSet);
		computeNbBit(radix, digitSet);
		return NULL;
	}

	void SRTDivNbBitsMin::registerFactory(){
		UserInterface::add("SRTDivNbBitsMin", // name
											 "A SRT design tool",
											 "Miscellaneous", // categories
											 "",
											 "radix(int): It has to be 2^n; \
alpha(int): digit set is [-alpha, alpha]",
											 "",
											 SRTDivNbBitsMin::NbBitsMinParseArguments
											 ) ;

	}


}
