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

#include "PositDiv.hpp"
#include "Posit/Decoder/PositFastDecoder.hpp"
#include "Posit/Encoder/PositFastEncoder.hpp"
// #include "IntAddSubCmp/IntAdder.hpp"
#include "FixFunctions/FixDiv.hpp"
#include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositDiv::PositDiv(OperatorPtr parentOp, Target *target, int width, int wES, int iters, bool useGoldschmidt, float dspOccupationThreshold, int LUT_out) : Operator(parentOp, target), width_(width), wES_(wES), iters_(iters), useGoldschmidt_(useGoldschmidt), dspOccupationThreshold(dspOccupationThreshold), LUT_out_(LUT_out)
	{
		setCopyrightString("Raul Murillo (2022)");
		ostringstream name;
		srcFileName = "PositDiv";
		// // Signed multiplication needs this library
		// useNumericStd();

		if (width_ < 3)
		{
			throw std::string("PositDiv Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositDiv Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositDiv_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositDiv \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES << ", " << iters << ", " << useGoldschmidt);
		REPORT(DEBUG, "debug of PositDiv");

		//====================================================================|
		addFullComment("Decode X & Y operands");
		//====================================================================|
		newInstance("PositFastDecoder",
					"X_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X",
					"Sign=>X_sgn,SF=>X_sf,Frac=>X_f,NZN=>X_nzn");

		newInstance("PositFastDecoder",
					"Y_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>Y",
					"Sign=>Y_sgn,SF=>Y_sf,Frac=>Y_f,NZN=>Y_nzn");

		//====================================================================|
		addFullComment("Divide X & Y");
		//====================================================================|
		addComment("Sign and Special Cases Computation");
		// vhdl << declare(getTarget()->logicDelay(2), "XY_sgn", 1, false) << " <= X_sgn XOR Y_sgn;" << endl;
		// vhdl << declare(getTarget()->logicDelay(2), "XY_z", 1, false) << " <= X_z OR Y_z;" << endl;
		// vhdl << declare(getTarget()->logicDelay(2), "XY_nar", 1, false) << " <= X_nar OR Y_nar;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_nzn", 1, false) << " <= X_nzn AND Y_nzn;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "X_nar", 1, false) << " <= X_sgn AND NOT(X_nzn);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "Y_nar", 1, false) << " <= Y_sgn AND NOT(Y_nzn);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "Y_zero", 1, false) << " <= NOT(Y_sgn OR Y_nzn);" << endl;


		// Using FloPoCo FixDivider
		addComment("Divide the fractions (using FloPoCo module FixDivider)");
		// getTarget()->setPlainVHDL(true);
		// getTarget()->setTilingMethod();
		vhdl << tab << declare(getTarget()->logicDelay(1), "XX_f", wF_ + 2) << " <= X_sgn & NOT(X_sgn) & X_f;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(1), "YY_f", wF_ + 2) << " <= Y_sgn & NOT(Y_sgn) & Y_f;" << endl;
		newInstance("FixDiv",
					"FracDivider",
					// "ints=1 frac=" + to_string(wF_) + " dspThreshold="+to_string(dspOccupationThreshold) + " truncate=false iters=0 useGoldschmidt=false",
					"ints=1 frac=" + to_string(wF_) + " dspThreshold="+to_string(dspOccupationThreshold) + " truncate=false iters=" + to_string(iters_) + " useGoldschmidt=" + to_string(useGoldschmidt_) + " LUT_out=" + to_string(LUT_out_),
					"X=>XX_f, Y=>YY_f",
					"R=>XY_f");

		// bool iterative = true;
		if (iters_ > 0){
			int divSize = 2 * (wF_ + 2);
			addComment("Normalize fraction");
			vhdl << tab << declare("XY_sgn", 1, false) << " <= XY_f" << of(divSize - 1) << ";" << endl;
			vhdl << tab << declare("XY_util_frac", 2*wF_+2) << " <= XY_f" << range(2*wF_+1, 0) << ";" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(2), "shift_1", 1, false) << " <= XY_f" << of(2*wF_+1) << " XNOR XY_f" << of(2*wF_) << ";" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(3), "shift_2", 1, false) << " <= XY_f" << of(2*wF_+1) << " AND XY_f" << of(2*wF_)  << " AND XY_f" << of(2*wF_-1) << ";" << endl; // Case XY_f = -0.5
			vhdl << tab << declare("shift_case", 2) << " <= shift_1 & shift_2;" << endl;

			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "XY_frac", wF_) << " <= " << endl
				<< tab << tab << "XY_util_frac" << range(2*wF_-3, wF_-2) << " when \"11\"," << endl	// shift 2 bits
				<< tab << tab << "XY_util_frac" << range(2*wF_-2, wF_-1) << " when \"10\"," << endl	// shift 1 bit
				<< tab << tab << "XY_util_frac" << range(2*wF_-1, wF_) << " when others;" << endl;	// no shift
			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "grd", 1, false) << " <= " << endl
				<< tab << tab << "XY_util_frac" << of(wF_-3) << " when \"11\"," << endl
				<< tab << tab << "XY_util_frac" << of(wF_-2) << " when \"10\"," << endl
				<< tab << tab << "XY_util_frac" << of(wF_-1) << " when others;" << endl;
			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "stk_tmp", wF_-1) << " <= " << endl;
			if (wF_>3) 
				vhdl << tab << tab << "XY_util_frac" << range(wF_-4, 0) << " & \"00\" when \"11\"," << endl;
			else
				vhdl << tab << tab << "\"00\" when \"11\"," << endl;
				vhdl << tab << tab << "XY_util_frac" << range(wF_-3, 0) << " & '0' when \"10\"," << endl
				<< tab << tab << "XY_util_frac" << range(wF_-2, 0) << " when others;" << endl;
			vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wF_-1), "stk", 1, false) << " <= "
				<< "'0' when (stk_tmp=0) else '1';" << endl;
			// vhdl << tab << "with shift_case select "  << declare(getTarget()->logicDelay(2)+getTarget()->eqConstComparatorDelay(wF_-1), "stk", 1, false) << " <= " << endl
			// 	 << tab << tab << "('0' when (XY_util_frac" << range(wF_-4, 0) << "=0) else '1') when \"11\"," << endl
			// 	 << tab << tab << "('0' when (XY_util_frac" << range(wF_-3, 0) << "=0) else '1') when \"10\"," << endl
			// 	 << tab << tab << "('0' when (XY_util_frac" << range(wF_-2, 0) << "=0) else '1') when others;" << endl;

		}
		else {
			// Sequential division
			// XY_f has same size as inputs X_f, Y_f (wF_ + 2)
#if 0
			// int divSize = (wF_ + 2);
			int divSize = (wF_ + 2+1);
			addComment("Normalize fraction");
			// vhdl << tab << declare("XY_sgn", 1, false) << " <= XY_f" << of(divSize - 1) << ";" << endl;
			// vhdl << tab << declare(getTarget()->logicDelay(2), "shift_1", 1, false) << " <= XY_f" << of(wF_+1) << " XNOR XY_f" << of(wF_) << ";" << endl;
			// vhdl << tab << declare(getTarget()->logicDelay(3), "shift_2", 1, false) << " <= XY_f" << of(wF_+1) << " AND XY_f" << of(wF_)  << " AND XY_f" << of(wF_-1) << ";" << endl; // Case XY_f = -0.5
			// vhdl << tab << declare("shift_case", 2) << " <= shift_1 & shift_2;" << endl;

			// vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "XY_frac", wF_) << " <= " << endl
			// 	<< tab << tab << "XY_f" << range(wF_-3, 0) << " & \"00\" when \"11\"," << endl	// shift 2 bits
			// 	<< tab << tab << "XY_f" << range(wF_-2, 0) << " & \"0\" when \"10\"," << endl	// shift 1 bit
			// 	<< tab << tab << "XY_f" << range(wF_-1, 0) << " when others;" << endl;	// no shift

			vhdl << tab << declare("XY_sgn", 1, false) << " <= XY_f" << of(divSize - 1) << ";" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(2), "shift_1", 1, false) << " <= XY_f" << of(wF_+1+1) << " XNOR XY_f" << of(wF_+1) << ";" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(3), "shift_2", 1, false) << " <= XY_f" << of(wF_+1+1) << " AND XY_f" << of(wF_+1)  << " AND XY_f" << of(wF_-1+1) << ";" << endl; // Case XY_f = -0.5
			vhdl << tab << declare("shift_case", 2) << " <= shift_1 & shift_2;" << endl;

			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "XY_frac", wF_) << " <= " << endl
				<< tab << tab << "XY_f" << range(wF_-3+1, 0) << " & \"0\" when \"11\"," << endl	// shift 2 bits
				<< tab << tab << "XY_f" << range(wF_-2+1, 0) << " when \"10\"," << endl	// shift 1 bit
				<< tab << tab << "XY_f" << range(wF_-1+1, 0+1) << " when others;" << endl;	// no shift

			// vhdl << tab << declare("grd", 1, false) << " <= '0';" << endl;
			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "grd", 1, false) << " <= " << endl
				<< tab << tab << "'0' when \"11\"," << endl	// shift 2 bits
				<< tab << tab << "'0' when \"10\"," << endl	// shift 1 bit
				<< tab << tab << "XY_f" << of(0) << " when others;" << endl;	// no shift
			vhdl << tab << declare("stk", 1, false) << " <= '0';" << endl;
#else		// Extra bits for the quotient
			int divSize = (wF_ + 2 + 4)-1;
			addComment("Normalize fraction");
			vhdl << tab << declare("XY_sgn", 1, false) << " <= XY_f" << of(divSize - 1) << ";" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(2), "shift_1", 1, false) << " <= XY_f" << of(divSize-1) << " XNOR XY_f" << of(divSize-2) << ";" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(3), "shift_2", 1, false) << " <= XY_f" << of(divSize-1) << " AND XY_f" << of(divSize-2)  << " AND XY_f" << of(divSize-3) << ";" << endl; // Case XY_f = -0.5
			vhdl << tab << declare("shift_case", 2) << " <= shift_1 & shift_2;" << endl;

			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "XY_frac", wF_) << " <= " << endl
				<< tab << tab << "XY_f" << range(divSize-5,divSize-wF_-4) << " when \"11\"," << endl	// shift 2 bits
				<< tab << tab << "XY_f" << range(divSize-4,divSize-wF_-3) << " when \"10\"," << endl	// shift 1 bit
				<< tab << tab << "XY_f" << range(divSize-3,divSize-wF_-2) << " when others;" << endl;	// no shift

			// Guard and Sticky bits
			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "grd", 1, false) << " <= " << endl
				<< tab << tab << "XY_f" << of(divSize-wF_-5) << " when \"11\"," << endl
				<< tab << tab << "XY_f" << of(divSize-wF_-4) << " when \"10\"," << endl
				<< tab << tab << "XY_f" << of(divSize-wF_-3) << " when others;" << endl;
			// vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "stk_tmp", 3) << " <= " << endl
			// 	<< tab << tab << "XY_f" << of(divSize-wF_-6) << " & \"00\" when \"11\"," << endl
			// 	<< tab << tab << "XY_f" << range(divSize-wF_-5, 0) << " & '0' when \"10\"," << endl
			// 	<< tab << tab << "XY_f" << range(divSize-wF_-4, 0) << " when others;" << endl;
			vhdl << tab << "with shift_case select " << declare(getTarget()->logicDelay(2), "stk_tmp", 3) << " <= " << endl
				<< tab << tab << "\"000\" when \"11\"," << endl
				<< tab << tab << "XY_f" << range(divSize-wF_-5, 0) << " & \"00\" when \"10\"," << endl
				<< tab << tab << "XY_f" << range(divSize-wF_-4, 0) << " & \"0\" when others;" << endl;
			vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wF_-1), "stk", 1, false) << " <= "
				<< "'0' when (stk_tmp=0) else '1';" << endl;
#endif
		}
	
		addComment("Subtract the exponent values");
		vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "XY_sf", wE_ + 1) << " <= "
			<< "(X_sf(X_sf'high) & X_sf) - (Y_sf(Y_sf'high) & Y_sf) - shift_1 - shift_2;" << endl;


		//====================================================================|
		addFullComment("Generate final posit");
		//====================================================================|
		vhdl << tab << declare(getTarget()->logicDelay(3), "XY_finalSgn", 1, false) << " <= "
			 << "XY_sgn when XY_nzn = '1' else (X_nar OR Y_nar OR Y_zero);" << endl;

		// vhdl << tab << declare("XY_frac", wF_) << " <= XY_normF" << range(divSize - 4, divSize - 4 - wF_ + 1) << ";" << endl;
		// vhdl << tab << declare("grd", 1, false) << " <= XY_normF" << of(divSize - 4 - wF_) << ";" << endl;
		// vhdl << tab << declare(getTarget()->eqConstComparatorDelay(divSize - 4 - wF_), "stk", 1, false) << " <= "
		// 	 << "'0' when (XY_normF" << range(divSize - 4 - wF_ - 1, 0) << " = " << zg(divSize - 4 - wF_) << ") else '1';" << endl;

		newInstance("PositFastEncoder",
					"PositEncoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"Sign=>XY_finalSgn, SF=>XY_sf, Frac=>XY_frac, Guard=>grd, Sticky=>stk, NZN=>XY_nzn",
					"R=>R");
		REPORT(DETAILED, "Finished generating operator PositDiv");

		addFullComment("End of vhdl generation");
	}

	PositDiv::~PositDiv() {}

	void PositDiv::emulate(TestCase *tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
		/* Compute correct value */
		PositNumber posx(width_, wES_, svX);
		PositNumber posy(width_, wES_, svY);
		mpfr_t x, y, r;
		mpfr_init2(x, 1000*width_ -2);
		mpfr_init2(y, 1000*width_ -2);
		mpfr_init2(r, 1000*width_ -2);
		posx.getMPFR(x);
		posy.getMPFR(y);
		mpfr_div(r, x, y, GMP_RNDN);
		
		// Set outputs
		PositNumber posr(width_, wES_, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		// clean up
		mpfr_clears(x, y, r, NULL);
	}

	void PositDiv::buildStandardTestCases(TestCaseList* tcl)
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

	OperatorPtr PositDiv::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES, iters, LUT_out;
		double dspOccupationThreshold=0.0;
		bool useGoldschmidt;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		UserInterface::parsePositiveInt(args, "iters", &iters);
		UserInterface::parseBoolean(args, "useGoldschmidt", &useGoldschmidt);
		UserInterface::parseFloat(args, "dspThreshold", &dspOccupationThreshold);
		UserInterface::parsePositiveInt(args, "LUT_out", &LUT_out);
		return new PositDiv(parentOp, target, width, wES, iters, useGoldschmidt, dspOccupationThreshold, LUT_out);
	}

	void PositDiv::registerFactory()
	{
		UserInterface::add("PositDiv", // name
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
						   PositDiv::parseArguments);
	}

} // namespace flopoco
