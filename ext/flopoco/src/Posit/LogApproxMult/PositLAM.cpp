/*
  A Posit logarithm-approximate multiplier (LAM) for FloPoCo
  
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

#include "PositLAM.hpp"
#include "Posit/Decoder/PositDecoder.hpp"
#include "Posit/Encoder/PositEncoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "IntAddSubCmp/IntAdder.hpp"
#include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositLAM::PositLAM(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositLAM";

		if (width_ < 3)
		{
			throw std::string("PositLAM Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositLAM Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositLAM_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositLAM \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositLAM");

		//=========================================================================|
		addFullComment("Decode X & Y operands");
		// ========================================================================|
		newInstance("PositDecoder",
					"X_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X",
					"Sign=>sign_X,SF=>sf_X,Frac=>f_X,Zero=>z_X,Inf=>inf_X,Abs_in=>open");

		newInstance("PositDecoder",
					"Y_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>Y",
					"Sign=>sign_Y,SF=>sf_Y,Frac=>f_Y,Zero=>z_Y,Inf=>inf_Y,Abs_in=>open");

		addComment("Gather operands");
		vhdl << tab << declare("op_X", wE_ + wF_) << " <= sf_X & f_X;" << endl;
		vhdl << tab << declare("op_Y", wE_ + wF_) << " <= sf_Y & f_Y;" << endl;

		//=========================================================================|
		addFullComment("Sign and Special cases computation");
		// ========================================================================|
		vhdl << tab << declare(getTarget()->logicDelay(2), "sign", 1, false) << " <= sign_X XOR sign_Y;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "zero", 1, false) << " <= z_X OR z_Y;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "inf", 1, false) << " <= inf_X OR inf_Y;" << endl;

		//=========================================================================|
		addFullComment("Add exponents & fractions all together");
		// ========================================================================|
		// TODO: Use IntAdder here(?)
		vhdl << tab << declare(getTarget()->adderDelay(wE_ + wF_ + 1), "add_r", wE_ + wF_ + 1) << " <= (op_X(op_X'high) & op_X) + (op_Y(op_Y'high) & op_Y);" << endl;

		vhdl << tab << declare("sf", wE_ + 1) << " <= add_r" << range(wE_ + wF_, wF_) << ";" << endl;
		vhdl << tab << declare("frac", wF_) << " <= add_r" << range(wF_ - 1, 0) << ";" << endl;

		//=========================================================================|
		addFullComment("Data Encoding");
		// ========================================================================|
		vhdl << tab << declare("rnd", 1, false) << " <= '0';" << endl;
		vhdl << tab << declare("stk", 1, false) << " <= '0';" << endl;

		newInstance("PositEncoder",
					"R_encoding",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"Sign=>sign, SF=>sf, Frac=>frac, Round=>rnd, Sticky=>stk, Zero=>zero, Inf=>inf",
					"R=>R");

		addFullComment("End of vhdl generation");
	}

	PositLAM::~PositLAM() {}

	void PositLAM::emulate(TestCase *tc)
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
		mpfr_mul(r, x, y, GMP_RNDN);
		
		// Set outputs
		PositNumber posr(width_, wES_, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		// clean up
		mpfr_clears(x, y, r, NULL);
	}

	void PositLAM::buildStandardTestCases(TestCaseList* tcl)
	{
		TestCase *tc;
		mpz_class x, y;

		// 1*1
		x = mpz_class(1) << (width_ - 2);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// 1*0
		x = mpz_class(1) << (width_ - 2);
		y = mpz_class(0);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1*-1
		x = mpz_class(3) << (width_ - 2);
		y = mpz_class(3) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1*1
		x = mpz_class(3) << (width_ - 2);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// nan*1
		x = mpz_class(1) << (width_ - 1);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// nan*0
		x = mpz_class(1) << (width_ - 1);
		y = mpz_class(0);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue*maxvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = (mpz_class(1) << (width_ - 1))-1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue*-maxvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = (mpz_class(1) << (width_ - 1))+1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -minvalue*maxvalue
		x = (mpz_class(1) << (width_))-1;
		y = (mpz_class(1) << (width_ - 1))-1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue*minvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = mpz_class(1);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}

	OperatorPtr PositLAM::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositLAM(parentOp, target, width, wES);
	}

	void PositLAM::registerFactory()
	{
		UserInterface::add("PositLAM", // name
						   "An approximate posit divider.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositLAM::parseArguments);
	}

} // namespace flopoco
