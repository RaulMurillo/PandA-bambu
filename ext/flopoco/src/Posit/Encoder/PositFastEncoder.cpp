/*
  Posit encoder
  
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

#include "PositFastEncoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositFastEncoder::PositFastEncoder(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositFastEncoder";

		if (width_ < 3)
		{
			throw std::string("PositFastEncoder Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositFastEncoder Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int regSize = intlog2(width_ - 1) + 1;
		wE_ = regSize + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositFastEncoder_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("Sign");
		addInput("SF", wE_ + 1);
		addInput("Frac", wF_);
		addInput("Guard");
		addInput("Sticky");
		addInput("NZN");
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositFastEncoder \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositFastEncoder");

		//====================================================================|
		addFullComment("Get value of regime");
		//====================================================================|
		vhdl << tab << declare("rc", 1, false) << " <= SF(SF'high);" << endl;
		vhdl << tab << declare("rcVect", regSize) << " <= (others => rc);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "k", regSize) << " <= "
			 << "SF" << range(wE_ - 1, wES_) << " XOR rcVect;" << endl;
		if (wES_ > 0)
		{
			vhdl << tab << declare("sgnVect", wES_) << " <= (others => Sign);" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(2), "exp", wES_) << " <= SF" << range(wES_ - 1, 0) << " XOR sgnVect;" << endl;
		}
		addComment("Check for regime overflow");
		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(regSize), "ovf", 1, false) << " <= "
			 << "'1' when (k > \"" << unsignedBinary(width_ - 3, regSize) << "\") else '0';" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(1), "regValue", regSize - 1) << " <= "
			 << "k" << range(regSize - 2, 0) << " when ovf = '0' else "
			 << "\"" << unsignedBinary(width_ - 2, regSize - 1) << "\";" << endl;

		//====================================================================|
		addFullComment("Generate regime - shift out exponent and fraction");
		//====================================================================|
		vhdl << tab << declare(getTarget()->logicDelay(2), "regNeg", 1, false) << " <= Sign XOR rc;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(1), "padBit", 1, false) << " <= NOT(regNeg);" << endl;

		vhdl << tab << declare("inputShifter", width_ - 1) << " <= regNeg ";
		if (wES_ > 0)
		{
			vhdl << "& exp ";
		}
		vhdl << "& Frac & Guard;" << endl;

		newInstance("Shifter",
					"RegimeGenerator",
					"wX=" + to_string(width_ - 1) + " wR=" + to_string(width_ - 1) + " maxShift=" + to_string(width_ - 1) + " dir=1 inputPadBit=true computeSticky=true",
					"X=>inputShifter, S=>regValue, padBit=>padBit",
					"R=>shiftedPosit, Sticky=>stkBit");

		vhdl << tab << declare("unroundedPosit", width_ - 1) << " <= padBit & shiftedPosit" << range(width_ - 2, 1) << ";" << endl;

		//====================================================================|
		addFullComment("Round to nearest even");
		//====================================================================|
		vhdl << tab << declare("lsb", 1, false) << " <= shiftedPosit" << of(1) << ";" << endl;
		vhdl << tab << declare("rnd", 1, false) << " <= shiftedPosit" << of(0) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "stk", 1, false) << " <= stkBit OR Sticky;" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(4), "round", 1, false) << " <= rnd AND (lsb OR stk OR ovf);" << endl;
		vhdl << tab << declare(getTarget()->adderDelay(width_ - 1), "roundedPosit", width_ - 1) << " <= "
			 << "unroundedPosit + round;" << endl;

		//====================================================================|
		addFullComment("Check sign & Special Cases");
		//====================================================================|
		vhdl << tab << declare("unsignedPosit", width_ - 1) << " <= "
			 << "roundedPosit when NZN = '1' else (others => '0');" << endl;

		// Prepare output
		vhdl << tab << "R <= Sign & unsignedPosit;" << endl;

		addFullComment("End of vhdl generation");
	}

	PositFastEncoder::~PositFastEncoder() {}

	void PositFastEncoder::emulate(TestCase *tc) {}

	OperatorPtr PositFastEncoder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositFastEncoder(parentOp, target, width, wES);
	}

	void PositFastEncoder::registerFactory()
	{
		UserInterface::add("PositFastEncoder", // name
						   "A posit encoder with a single architecture.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositFastEncoder::parseArguments);
	}

} // namespace flopoco
