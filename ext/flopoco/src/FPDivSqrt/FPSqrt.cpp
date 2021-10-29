/*
  Floating Point Square Root for FloPoCo

  Authors :
  Jeremie Detrey, Florent de Dinechin (digit-recurrence version)

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

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
#include "FPSqrt.hpp"

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0
	//#define LESS_DSPS

	FPSqrt::FPSqrt(OperatorPtr parentOp, Target* target, int wE_, int wF_, int method_) :
		Operator(parentOp,target), wE(wE_), wF(wF_), method(method_) {

		ostringstream name;

		name<<"FPSqrt_"<<wE<<"_"<<wF;

		uniqueName_ = name.str();

		if(method==0) {
			REPORT(0, "Floating-point square root using plain restoring algorithm");
		}
		else if(method==1) {
			REPORT(0, "Floating-point square root using binary nonrestoring algorithm");
		}
		else  {
			THROWERROR("Wrong value passed to parameter: method=" << method);
		}

		// -------- Parameter set up -----------------

		addFPInput ("X", wE, wF);
		addFPOutput("R", wE, wF);


		vhdl << tab << declare("fracX", wF) << " <= X" << range(wF-1, 0) << "; -- fraction"  << endl;
		vhdl << tab << declare("eRn0", wE) << " <= \"0\" & X" << range(wE+wF-1, wF+1) << "; -- exponent" << endl;
		vhdl << tab << declare("xsX", 3) << " <= X"<< range(wE+wF+2, wE+wF) << "; -- exception and sign" << endl;

		vhdl << tab << declare("eRn1", wE) << " <= eRn0 + (\"00\" & " << rangeAssign(wE-3, 0, "'1'") << ") + X(" << wF << ");" << endl;


		
		if(method==0) {
			vhdl << tab << "-- now implementing the recurrence: d_i has position -i " << endl;
			// R1 = 2R0 - 2S0 -2^-1 d_0
			vhdl << tab << "--  this is a binary restoring algorithm, see e.g. Parhami book 2nd ed. p. 438" << endl;
			vhdl << tab << "--  recurrence is R_i = 2R_{i-1} -2d_iS_{i-1} - 2^{-i}di^2  for i in {1..n}" << endl;
			vhdl << tab << "--  so try  T_i = 2R_{i-1} -2S_{i-1} - 2^{-i} " << endl;
			vhdl << tab << "--  If T_i>=0 then d_{i}=1 and R_i=T_i " << endl;
			vhdl << tab << "--  If T_i<0 then d_{i}=0 and R_i=2R_{i-1} " << endl;
			vhdl << tab << "--  this is a binary restoring algorithm, see e.g. Parhami book 2nd ed. p. 438" << endl;
			// R_i has a constant size: wF+1+1+1 (implicit bit, + norm bit + round bit), format ufix(1,-wF-1)
			// S_i has size i+1  e.g. S0 = 1, format ufix(0,-i);
			// d_i \in |{0,1\}  has position -i in S_i
			// 2S_{i-1} has format ufix(1,-i+2)
			// d_i^2 == d_i is added at position -i
			// so tentative subtrahand is  2S_{i-1} + 2^{-i}    is   2S_{i-1} 01, format ufix(1,-i)  
			// but subtraction has result on sfix(2,-wF-1) to keep the sign,
			// so effective subtraction on format sfix(2, -i) : size i+3 
			vhdl << tab << " -- initialization " << endl;
			vhdl << tab << declare("d0") << " <= '1';" << endl;
			vhdl << tab << declare("S0", 1) << " <= \"1\";" << endl;
			vhdl << tab << declare(getTarget()->lutDelay(), "R0", wF+3)
					 << " <= \"00\" & fracX & \"0\" when X(" << wF << ") = '1' else   -- parity of EX-E0 is opposite to that of EX" << endl
					 << tab << "       fracX"<<of(wF-1) <<" & (not fracX"<<of(wF-1) <<") & fracX" << range(wF-2,0) << " & \"00\"; -- pre-normalization" << endl;
		//		vhdl << tab << declare(join("d",wF+3)) << " <= '0';" << endl;
		//		vhdl << tab << declare(join("s",wF+3)) << " <= '1';" << endl;
			int maxstep=wF+1;
			for(int i=1; i<=maxstep; i++) {
				vhdl << tab << "-- Step " << i << endl;
				string TwoRim1 = "R" + to_string(i-1) + "s";
				string Ti = join("T", i);
				string di = join("d", i);
				string Ri = join("R", i);
				string Rim1 = join("R", i-1);
				string TwoRim1H = TwoRim1 + "_h";
				string TwoRim1L = TwoRim1 + "_l";
				string Si = join("S", i);
				string Sim1 = join("S", i-1);

				int subsize=i+3; // this is a sfix(2,-i)
				vhdl << tab << declare(TwoRim1,wF+4) << " <= " << Rim1 << " & \"0\";" << endl;
				vhdl << tab << declare(TwoRim1H,subsize) << " <= " << TwoRim1 << range(wF+4-1, wF+4-subsize) << ";" << endl;
				if(i <= wF) {
					vhdl << tab << declare(TwoRim1L,wF+4-subsize) << " <= "  << TwoRim1 << range(wF+4-subsize-1, 0) << ";" << endl;
				}
				vhdl << tab << declare(getTarget()->adderDelay(subsize), Ti, subsize)
						 << " <= " << TwoRim1H << " - (\"0\" & "<< Sim1 << " & \"01\"); -- tentative subtraction "<< endl;
				vhdl << tab << declare(di) << " <= not " << Ti << of(subsize-1) << "; -- next digit"<< endl;
				vhdl << tab << declare(Si, i+1) << " <= " << Sim1 << " & " << di << "; "<< endl;
				string Rih = "R" + to_string(i) + "_h";
				vhdl << tab << declare(getTarget()->logicDelay(), Rih, subsize-1)
						 << " <= " << Ti << range(subsize-2,0)
						 << "   when " << di << "= '1' else" << endl
						 << tab << "       " << TwoRim1H << range(subsize-2,0) <<"; " << endl;
				vhdl << tab << declare(Ri, wF+3)
						 << " <= " << Rih  << (i <= wF? " & " + TwoRim1L: "") <<"; " << endl;
				
			} // end  digit iteration

			vhdl << tab << declare(target->lutDelay(), "fR", wF) << " <= "<< join("S", maxstep) <<range(wF, 1) << ";-- removing leading 1" << endl;
			vhdl << tab << declare("round") << " <= " << join("d", maxstep) << "; -- round bit" << endl;

		}



		else if(method==1) {
			// Digit-recurrence implementation recycled from FPLibrary: works better!
			// Sorry for the completely inconsistent signal names in the C++,
			// this code was incrementally modified to match the mames and indices in the ASA book, and history shows.
			// TODO refactor the change of variable i -> i-1 in the C++, and rename R into T etc
			vhdl << tab << declare(getTarget()->lutDelay(),
														 "fracXnorm", wF+4) << " <= \"1\" & fracX & \"000\" when X(" << wF << ") = '0' else" << endl
					 << tab << "      \"01\" & fracX&\"00\"; -- pre-normalization" << endl;
			vhdl << tab << declare("S0", 2) << " <= \"01\";" << endl;

			vhdl << tab << declare(getTarget()->adderDelay(4),
														 "T1", wF+4) << " <= (\"0111\" + fracXnorm" << range(wF+3, wF)<< ") & fracXnorm" << range(wF-1, 0)<< ";"<< endl;

			
		//		vhdl << tab << declare(join("d",wF+3)) << " <= '0';" << endl;
		//		vhdl << tab << declare(join("s",wF+3)) << " <= '1';" << endl;
			vhdl << tab << "-- now implementing the recurrence " << endl;
			//			vhdl << tab << "--  w_{i} = 2w_{i-1} -2s_{i}S_{i-1} - 2^{-i-1}s_{i}^2  for i in {1..n}" << endl;
			vhdl << tab << "--  this is a binary non-restoring algorithm, see ASA book" << endl;
			int maxstep=wF+2;
			for(int i=3; i<=maxstep; i++) {
				double stageDelay= getTarget()->adderDelay(i);
				REPORT(2, "estimated delay for stage "<< i << " is " << stageDelay << "s");
				// was: int i = wF+3-step; // to have the same indices as FPLibrary
				vhdl << tab << "-- Step " << i-1 << endl;
				string di = join("d", i-2);
				string TwoRim1 = "T" + to_string(i-2) + "s";
				string Ri = join("T", i-1);
				string Rim1 = join("T", i-2);
				string Si = join("S", i-2);
				string Sim1 = join("S", i-3);
				//			string zs = join("zs", i);
				string ds = join("U", i-2);
				string TwoRim1H = TwoRim1 + "_h";
				string TwoRim1L = TwoRim1 + "_l";
				string wh = "T" + to_string(i) + "_h";
				vhdl << tab << declare(di) << " <= not "<< Rim1 << "("<< wF+3<<"); --  bit of weight "<< -(i-2) << endl;
				vhdl << tab << declare(TwoRim1,wF+5) << " <= " << Rim1 << " & \"0\";" << endl;
				vhdl << tab << declare(TwoRim1H,i+3) << " <= " << TwoRim1 << range(wF+4, wF+2-i) << ";" << endl;
				if(i <= wF+1) {
					vhdl << tab << declare(TwoRim1L,wF+2-i) << " <= "  << TwoRim1 << range(wF+1-i, 0) << ";" << endl;
				}
				vhdl << tab << declare(ds,i+3) << " <=  \"0\" & ";
				if (i>1)
					vhdl 	<< Sim1 << " & ";
				vhdl << di  << " & (not " << di << ")" << " & \"1\"; " << endl;
				vhdl << tab <<  declare(stageDelay, wh, i+3) << " <=   " << TwoRim1H << " - " << ds << " when " << di << "='1'" << endl
						 << tab << tab << "  else " << TwoRim1H << " + " << ds << ";" << endl;
				vhdl << tab << declare(Ri, wF+4) << " <= " << wh << range(i+1,0);
				if(i <= wF+1)
					vhdl << " & " << TwoRim1L << ";" << endl;
				else
					vhdl << ";" << endl;
				vhdl << tab << declare(Si, i) << " <= ";
				if(i==1)
					vhdl << "\"\" & " << di << ";"<< endl;
				else
					vhdl << Sim1 /*<< range(i-1,1)*/ << " & " << di << "; -- here -1 becomes 0 and 1 becomes 1"<< endl;
			}
			string dfinal=join("d", maxstep);
			vhdl << tab << declare(dfinal) << " <= not "<< join("T", maxstep-1) << of(wF+3)<<" ; -- the sign of the remainder will become the round bit" << endl;
			vhdl << tab << declare("mR", wF+3) << " <= "<< join("S", maxstep-2)<<" & "<<dfinal<<"; -- result significand" << endl;

			// end of component FPSqrt_Sqrt in fplibrary
			vhdl << tab << declare("fR", wF) << " <= mR" <<range(wF, 1) << ";-- removing leading 1" << endl;
			vhdl << tab << declare("round") << " <= mR(0); -- round bit" << endl;


		} // End method select


			vhdl << tab << declare(target->adderDelay(wF),
														 "fRrnd", wF) << " <= fR" <<" + (" << rangeAssign(wF-1, 1, "'0'") << " & round); -- rounding sqrt never changes exponents " << endl;
		vhdl << tab << declare("Rn2", wE+wF) << " <= eRn1 & fRrnd;" << endl;

		vhdl << tab << "-- sign and exception processing" << endl;
		vhdl << tab <<  "with xsX select" << endl
		     << tab << tab << declare(target->lutDelay(), "xsR", 3) << " <= \"010\"  when \"010\",  -- normal case" << endl
		     << tab << tab <<  "       \"100\"  when \"100\",  -- +infty" << endl
		     << tab << tab <<  "       \"000\"  when \"000\",  -- +0" << endl
		     << tab << tab <<  "       \"001\"  when \"001\",  -- the infamous sqrt(-0)=-0" << endl
		     << tab << tab <<  "       \"110\"  when others; -- return NaN" << endl;

		vhdl << tab << "R <= xsR & Rn2; " << endl;
	}

  FPSqrt::~FPSqrt() {
  }






		void FPSqrt::emulate(TestCase * tc)
		{
			/* Get I/O values */
			mpz_class svX = tc->getInputValue("X");

			/* Compute correct value */
			FPNumber fpx(wE, wF);
			fpx = svX;
			mpfr_t x, r;
			mpfr_init2(x, 1+wF);
			mpfr_init2(r, 1+wF);
			fpx.getMPFR(x);

			if(correctRounding) {
				mpfr_sqrt(r, x, GMP_RNDN);
				FPNumber  fpr(wE, wF, r);
				/* Set outputs */
				mpz_class svr= fpr.getSignalValue();
				tc->addExpectedOutput("R", svr);
			}
			else { // faithful rounding
				mpfr_sqrt(r, x, GMP_RNDU);
				FPNumber  fpru(wE, wF, r);
				mpz_class svru = fpru.getSignalValue();
				tc->addExpectedOutput("R", svru);

				mpfr_sqrt(r, x, GMP_RNDD);
				FPNumber  fprd(wE, wF, r);
				mpz_class svrd = fprd.getSignalValue();
				/* Set outputs */
				tc->addExpectedOutput("R", svrd);
			}

			mpfr_clears(x, r, NULL);
		}



		void FPSqrt::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		// Regression tests
		tc = new TestCase(this);
		tc->addFPInput("X", 1.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);
		tc->addFPInput("X", 25.0);
		emulate(tc);
		tcl->add(tc);

		tc = new TestCase(this);
		tc->addFPInput("X", 2.0);
		emulate(tc);
		tcl->add(tc);

	}




		// One test out of 4 fully random (tests NaNs etc)
		// All the remaining ones test positive numbers.
		TestCase* FPSqrt::buildRandomTestCase(int i){

			TestCase *tc;
			mpz_class a;

			tc = new TestCase(this);
			/* Fill inputs */
			if ((i & 3) == 0)
				a = getLargeRandom(wE+wF+3);
			else
				a  = getLargeRandom(wE+wF) + (mpz_class(1)<<(wE+wF+1)); // 010xxxxxx
			tc->addInput("X", a);

			/* Get correct outputs */
			emulate(tc);

			return tc;
		}

		OperatorPtr FPSqrt::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
			int wE;
			UserInterface::parseStrictlyPositiveInt(args, "wE", &wE);
			int wF;
			UserInterface::parseStrictlyPositiveInt(args, "wF", &wF);
			int method;
			UserInterface::parsePositiveInt(args, "method", &method);
			return new FPSqrt(parentOp, target, wE, wF, method);
		}


	TestList FPSqrt::unitTest(int index)
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

	void FPSqrt::registerFactory(){
			UserInterface::add("FPSqrt", // name
												 "A correctly rounded floating-point square root function.",
												 "BasicFloatingPoint", // categories
												 "",
												 "wE(int): exponent size in bits; \
wF(int): mantissa size in bits; \
method(int)=1: 0 for plain restoring, 1 for nonrestoring method",
												 "",
												 FPSqrt::parseArguments,
												 FPSqrt::unitTest
												 ) ;

		}
	}
