/*
  A Posit Divider for FloPoCo
  
  Author:  Raul Murillo

  This file is part of the FloPoCo project
  
  Initial software.
  Copyright © UCM, 
  2021.
  All rights reserved.

*/

#include <iostream>
#include <sstream>

#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "PositSqrt.hpp"
#include "Posit/Decoder/PositFastDecoder.hpp"
#include "Posit/Encoder/PositFastEncoder.hpp"
// #include "IntAddSubCmp/IntAdder.hpp"
#include "FixFunctions/FixDiv.hpp"
#include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	// PositSqrt::PositSqrt(OperatorPtr parentOp, Target *target, int width, int wES, int iters, bool useGoldschmidt, float dspOccupationThreshold, int LUT_out) : Operator(parentOp, target), width_(width), wES_(wES), iters_(iters), useGoldschmidt_(useGoldschmidt), dspOccupationThreshold(dspOccupationThreshold), LUT_out_(LUT_out)
	PositSqrt::PositSqrt(OperatorPtr parentOp, Target *target, int width, int wES, float dspOccupationThreshold) : Operator(parentOp, target), width_(width), wES_(wES), dspOccupationThreshold(dspOccupationThreshold)
	{
		setCopyrightString("Raul Murillo (2022)");
		ostringstream name;
		srcFileName = "PositSqrt";
		// // Signed multiplication needs this library
		// useNumericStd();

		if (width_ < 3)
		{
			throw std::string("PositSqrt Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositSqrt Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositSqrt_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositSqrt \n");
		// REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES << ", " << iters << ", " << useGoldschmidt);
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositSqrt");

		//====================================================================|
		addFullComment("Decode X operand");
		//====================================================================|
		newInstance("PositFastDecoder",
					"X_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X",
					"Sign=>X_sgn,SF=>X_sf,Frac=>X_f,NZN=>X_nzn");

		addComment("Sign and Special Cases Computation");
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_nzn", 1, false) << " <= NOT(X_sgn) AND X_nzn;" << endl;
		vhdl << tab << declare("XY_finalSgn", 1, false) << " <= X_sgn;" << endl;

		//====================================================================|
		addFullComment("Exponent computation");
		//====================================================================|
		// To obtain an integer, exponent E should be even.
		// Consequently, if E is odd, we utilize frac/2 and E-1.
		// vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "X_sf_1", wE_ + 1) << " <= "
		// 	<< "(X_sf(X_sf'high) & X_sf) when odd_exp='1' else (X_sf(X_sf'high) & X_sf) + 1;" << endl;
		vhdl << tab << declare("odd_exp", 1, false) << " <= X_sf(0);" << endl;
		addComment("Divide exponent by 2");
		vhdl << tab << declare("X_sf_3", wE_ + 1) << " <= "
			<< "X_sf(X_sf'high) & X_sf(X_sf'high) & X_sf" << range(wE_-1, 1) << ";" << endl;


		//====================================================================|
		addFullComment("Sqrt of the fraction");
		//====================================================================|
		// TODO: https://ieeexplore.ieee.org/document/624623?arnumber=624623

#if 1
		//====================================================================|
		addFullComment("Non-Restoring algorithm");
		//====================================================================|

		const int w = 3+wF_;
		vhdl << tab << declare("one_bit", 1, false) << " <= '1';" << endl;

		// Ensure that X < 1. Need to adjust exponent.
		// We divide frac by 4, so sqrt(1+f) = sqrt((1+f)/4)*2
		// If exponent E is odd, we consider 2*(1+f) as initial fraction, so just divide by 2 to get sqrt(2*(1+f))=sqrt((1+f)/2)*2
		// TODO: append more bits (?)
		// // vhdl << tab << declare("append_0", 1, false) << " <= '0' when odd_exp='1' else X_f" << of(0) << ";" << endl;
		vhdl << tab << declare("r_0", w+1) << " <= (\"001\" & X_f & '0') when odd_exp='1' else (\"0001\" & X_f);" << endl;
		vhdl << tab << declare("q_0", w) << " <= (others => '0');" << endl;
		vhdl << tab << declare("real_q_0", w) << " <= (others => '0');" << endl;
		vhdl << tab << declare("pow_2_0", w+1) << " <= \"01" << zg(w-2+1, 1) << ";" << endl;

		string r_i, two_r_i, r_i1, q_i, q_i1, two_q_i, s_i, pow_2i, pow_2i1, n_i, rem_z, z_i, z_i1, real_q_i, real_q_i1;
		for (int i = 0; i < w; ++i){
			addComment("Iteration "+std::to_string(i+1));
			r_i = join("r_", i);
			two_r_i = join("two_r_", i);
			q_i = join("q_", i);
			two_q_i = join("two_q_", i);
			r_i1 = join("r_", i+1);
			q_i1 = join("q_", i+1);
			s_i = join("s_", i);
			pow_2i = join("pow_2_", i);
			pow_2i1 = join("pow_2_", i+1);
			n_i = join("n_", i);
			rem_z = join("rem_z_", i);
			z_i = join("z_", i);
			z_i1 = join("z_", i+1);
			real_q_i = join("real_q_", i);
			real_q_i1 = join("real_q_", i+1);

			// Check sign of r_i
			vhdl << tab << declare(s_i, 1, false) << " <= "<< r_i << of(w-1+1) << ";" << endl;

			// q_{i+1} = NOT(s_i)
			// r_{i+1} = 2*r_i -+ (2*Q_i +- 2^{-i})
			// 2*r_i
			// vhdl << tab << declare(two_r_i, w) << " <= " << r_i << range(w-2, 0) << " & '0';" << endl;
			if (i==0) {
				vhdl << tab << declare(getTarget()->logicDelay(1), q_i1, w) << " <= NOT(" << s_i << ") & " << zg(w-1) << ";" << endl;
				vhdl << tab << declare(real_q_i1, w) << " <= (" << s_i << ") & " << zg(w-1) << ";" << endl;
				
				vhdl << tab << declare(two_r_i, w+1) << " <= " << r_i << range(w-2+1, 0) << " & '0';" << endl;
				vhdl << tab << declare(two_q_i, w+1) << " <= (others => '0');" << endl;
			}
			else {
				//  TODO: This is not reusable. Use an OR operation instead.
				vhdl << tab << declare(getTarget()->logicDelay(2), q_i1, w) << " <= " << q_i << range(w-1,w-i) << " & NOT(" << s_i << " OR " << z_i << ") & " << zg(w-i-1) << ";" << endl;
				vhdl << tab << declare(getTarget()->logicDelay(2), real_q_i1, w) << " <= " << q_i1 << range(w-2, w-(i+1)) << " & '1' & " << zg(w-(i+1)) << " when " << z_i << "='0' else " << real_q_i << ";" << endl;
				// vhdl << tab << declare(getTarget()->logicDelay(2), q_i1, wF_) << " <= " << q_i << " when " << z_i1 << "='1' else (" << q_i << " OR " << pow_2i1 << range(wF_-1, 0) << ");" << endl;

				vhdl << tab << declare(two_r_i, w+1) << " <= " << r_i << range(w-2+1, 0) << " & '0';" << endl;
				if (i==1) {
					vhdl << tab << declare(two_q_i, w+1) << " <= '0' & '1' & " << zg(w-(1+i)+1) << ";" << endl;
				}
				else {
					vhdl << tab << declare(two_q_i, w+1) << " <= '0' & " << q_i << range(w-2, w-i) << " & '1' & " << zg(w-(1+i)+1) << ";" << endl;
				}
			}
			// vhdl << tab << declare(two_q_i, w) << " <= '0' & " << q_i << range(w-1, 1) << ";" << endl;
			// vhdl << tab << declare(pow_2i, w) << " <= \"00\" & " << zg(i) << " & '1' & " << zg(w-i-3) << ";" << endl;
			vhdl << tab << declare(pow_2i1, w+1) << " <= '0' & " << pow_2i << range(w-1+1, 1) << ";" << endl;	// Right shift pow_2i

			vhdl << tab << declare(getTarget()->adderDelay(w+1)+getTarget()->logicDelay(1), n_i, w+1) << " <= (" << two_q_i << " + NOT(" << pow_2i1 << ")) when " << s_i << "='1' else NOT(" << two_q_i << " + " << pow_2i1 << ");" << endl;
			
			newInstance("IntAdder", 
						join("sub_", i+1),
						"wIn="+to_string(w+1),
						"X=>"+two_r_i+",Cin=>one_bit,Y=>" + n_i,
						"R=>"+r_i1
						);

			// Detect if remainder n_i1 is zero. In such a case, do not modify partial result q_i.
			// TODO: Revise. In such a case, final result is different
			vhdl << tab << declare(getTarget()->eqConstComparatorDelay(w+1), rem_z, 1, false) << " <= '1' when " << r_i1 << " = 0 else '0';" << endl;
			if (i==0){
				vhdl << tab << declare(z_i1, 1, false) << " <= " << rem_z << ";" << endl;
			}
			else{
				vhdl << tab << declare(getTarget()->logicDelay(2), z_i1, 1, false) << " <= " << rem_z << " OR " << z_i << ";" << endl;
			}
		}

		addComment("Convert the quotient to the digit set {0,1}");
		// // q_pos is the positive term (original Q, since -1 digits are stored as zeros)
		// // q_neg is the negative term (one's complement on the original Q)
		// vhdl << tab << declare("q_pos", w) << " <= " << q_i1 << ";" << endl;
		// vhdl << tab << declare("q_neg", w) << " <= " << real_q_i1 << ";" << endl;
		// // vhdl << tab << declare("zeros", w) << " <= '0'";
		// // for (int i = 1; i<w; ++i){
		// // 	vhdl << " & " << join("z_", i);
		// // }
		// // vhdl << " ;" << endl;
		// // vhdl << tab << declare("q_neg", w) << " <= " << q_i1 << " NOR zeros;" << endl;
		// // quotient = q_pos - q_neg
		// newInstance("IntAdder",
		// 			"sub_quotient",
		// 			"wIn="+to_string(w),
		// 			"X=>q_pos,Cin=>one_bit,Y=>q_neg",
		// 			"R=>sqrt_f"
		// 			);
		// Shortcut: In this case, we can assume the following: (see Israel Koren book)
		vhdl << tab << declare("sqrt_f", w) << " <= " << q_i1 << range(w-2, 0) << " & '1' when " << z_i << "='0' else " << real_q_i1 << "; -- get the double of sqrt: first bit (=0) shifted out" << endl;

#else
		int wF = wF_;
		vhdl << tab << declare("fracX", wF) << " <= X_f;" << endl;

		// FloPoCo FPSqrt ****
			// Digit-recurrence implementation recycled from FPLibrary: works better!
			// Sorry for the completely inconsistent signal names in the C++,
			// this code was incrementally modified to match the mames and indices in the ASA book, and history shows.
			// TODO refactor the change of variable i -> i-1 in the C++, and rename R into T etc
			vhdl << tab << declare(getTarget()->lutDelay(),
														 "fracXnorm", wF+4) << " <= \"1\" & fracX & \"000\" when odd_exp = '0' else" << endl
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

#endif

		//====================================================================|
		addFullComment("Generate final posit");
		//====================================================================|

		// vhdl << tab << declare("XY_frac", wF_) << " <= XY_normF" << range(divSize - 4, divSize - 4 - wF_ + 1) << ";" << endl;
		// vhdl << tab << declare("grd", 1, false) << " <= XY_normF" << of(divSize - 4 - wF_) << ";" << endl;
		// vhdl << tab << declare(getTarget()->eqConstComparatorDelay(divSize - 4 - wF_), "stk", 1, false) << " <= "
		// 	 << "'0' when (XY_normF" << range(divSize - 4 - wF_ - 1, 0) << " = " << zg(divSize - 4 - wF_) << ") else '1';" << endl;


		vhdl << tab << declare("XY_sf", wE_ + 1) << " <= X_sf_3;" << endl;
		vhdl << tab << declare("XY_frac", wF_) << " <= sqrt_f" << range(w-2,2) << ";" << endl;	// Multiply by 2
		vhdl << tab << declare("grd", 1, false) << " <= sqrt_f" << of(1) << ";" << endl;
		vhdl << tab << declare("stk", 1, false) << " <= sqrt_f" << of(0) << ";" << endl;
		// vhdl << tab << declare("XY_frac", wF_) << " <= fR;" << endl;	// Multiply by 2
		// vhdl << tab << declare("grd", 1, false) << " <= round;" << endl;
		// vhdl << tab << declare("stk", 1, false) << " <= '0';" << endl;

		newInstance("PositFastEncoder",
					"PositEncoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"Sign=>XY_finalSgn, SF=>XY_sf, Frac=>XY_frac, Guard=>grd, Sticky=>stk, NZN=>XY_nzn",
					"R=>R");
		REPORT(DETAILED, "Finished generating operator PositSqrt");

		addFullComment("End of vhdl generation");
	}

	PositSqrt::~PositSqrt() {}

	void PositSqrt::emulate(TestCase *tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		
		/* Compute correct value */
		PositNumber posx(width_, wES_, svX);
		mpfr_t x, r;
		mpfr_init2(x, 1000*width_ -2);
		mpfr_init2(r, 1000*width_ -2);
		posx.getMPFR(x);
		mpfr_sqrt(r, x, GMP_RNDN);
		
		// Set outputs
		PositNumber posr(width_, wES_, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		// clean up
		mpfr_clears(x, r, NULL);
	}

	void PositSqrt::buildStandardTestCases(TestCaseList* tcl)
	{
		TestCase *tc;
		mpz_class x, y;

		// 1/1
		x = mpz_class(1) << (width_ - 2);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// 0/1
		x = mpz_class(0);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// 1/0
		x = mpz_class(1) << (width_ - 2);
		y = mpz_class(0);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1/-1
		x = mpz_class(3) << (width_ - 2);
		y = mpz_class(3) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1/1
		x = mpz_class(3) << (width_ - 2);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// nan/1
		x = mpz_class(1) << (width_ - 1);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// nan/0
		x = mpz_class(1) << (width_ - 1);
		y = mpz_class(0);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue/maxvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = (mpz_class(1) << (width_ - 1))-1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue/-maxvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = (mpz_class(1) << (width_ - 1))+1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -minvalue/maxvalue
		x = (mpz_class(1) << (width_))-1;
		y = (mpz_class(1) << (width_ - 1))-1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue/minvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = mpz_class(1);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}

	OperatorPtr PositSqrt::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES, iters, LUT_out;
		double dspOccupationThreshold=0.0;
		bool useGoldschmidt;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		// UserInterface::parsePositiveInt(args, "iters", &iters);
		// UserInterface::parseBoolean(args, "useGoldschmidt", &useGoldschmidt);
		// UserInterface::parsePositiveInt(args, "LUT_out", &LUT_out);
		UserInterface::parseFloat(args, "dspThreshold", &dspOccupationThreshold);
		// return new PositSqrt(parentOp, target, width, wES, iters, useGoldschmidt, dspOccupationThreshold, LUT_out);
		return new PositSqrt(parentOp, target, width, wES, dspOccupationThreshold);
	}

	void PositSqrt::registerFactory()
	{
		UserInterface::add("PositSqrt", // name
						   "A correctly rounded posit multiplier.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
						   wES(int): posit exponent size in bits;\
						   iters(int)=0: if positive, implements an iterative algorithm;\
						   useGoldschmidt(bool)=false: if positive and iters>0, implements Goldschmidt algorithm;\
						   dspThreshold(real)=0.0: threshold of relative occupation ratio of a DSP multiplier to be used or not;\
                           LUT_out(int)=10: if using an iterative algorithm, indicates the output size of LUT for initial guess",
						   "", // htmldoc
						   PositSqrt::parseArguments);
	}

} // namespace flopoco
