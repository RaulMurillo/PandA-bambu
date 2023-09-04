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

#include "PositEncoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositEncoder::PositEncoder(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021-2022)");
		ostringstream name;
		srcFileName = "PositEncoder";

		if (width_ < 3)
		{
			throw std::string("PositEncoder Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositEncoder Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int regSize = intlog2(width_ - 1) + 1;
		wE_ = regSize + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositEncoder_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("Sign");
		addInput("SF", wE_ + 1); // +1 to handle overflow
		addInput("Frac", wF_);
		addInput("Round");
		addInput("Sticky");
		addInput("Zero");
		addInput("Inf");
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositEncoder \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositEncoder");

		//====================================================================|
		addFullComment("Get value of regime");
		//====================================================================|
		// Sanity check
		if (wE_ - wES_ != regSize)
			throw std::string("Sorry, but parameters for posit encoder do not match sizes.");

		if (wES_ > 0)
		{
			vhdl << tab << declare("e", wES_) << " <= SF" << range(wES_ - 1, 0) << ";" << endl;
		}
		vhdl << tab << declare("k", regSize) << " <= SF" << range(wE_ - 1, wES_) << ";" << endl;
		vhdl << tab << declare("rc", 1, false) << " <= SF" << of(wE_) << ";" << endl;
		vhdl << tab << declare("v_rc", regSize) << " <= (others => rc) ;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "offset_tmp", regSize) << " <= k XOR v_rc;" << endl;

		addComment("Check for regime overflow");
		uint32_t maxShift = width_ - 2;
		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(regSize), "reg_ovf", 1, false)
			 //  << " <= '1' when (offset_tmp > CONV_STD_LOGIC_VECTOR(" << maxShift - 1 << "," << regSize << ")) else '0';" << endl;
			 << " <= '1' when (offset_tmp >= " << maxShift << ") else '0';" << endl;

		//=========================================================================|
		addFullComment("Generate regime - shift out exponent and fraction");
		// ========================================================================|
		vhdl << tab << declare(getTarget()->logicDelay(1), "pad", 1, false) << " <= not rc;" << endl;
		vhdl << tab << declare("input_shifter", 2 + wES_ + wF_ + 1) << " <= pad & rc ";
		if (wES_ > 0)
		{
			vhdl << "& e ";
		}
		vhdl << "& Frac & Round;" << endl;

		// Generate exp-frac offset to append regime
		int shift_size = intlog2(width_ - 2);
		vhdl << tab << declare(getTarget()->logicDelay(1), "shift_offset", shift_size) << " <= "
			 << "\"" << unsignedBinary(maxShift, shift_size) << "\" when reg_ovf = '1' else "
			 << "offset_tmp" << range(shift_size - 1, 0) << ";" << endl;

		//Sanity check
		if (width_ != wES_ + wF_ + 1 + 2)
			throw std::string("Sorry, but parameters for posit encoder do not match sizes.");

		ostringstream param_shift, inmap_shift, outmap_shift;
		param_shift << "wX=" << width_;
		param_shift << " maxShift=" << maxShift;
		param_shift << " wR=" << width_; // if(computeSticky) wOut =  wIn;
		param_shift << " dir=" << Shifter::Right;
		param_shift << " computeSticky=true";
		param_shift << " inputPadBit=true";
		inmap_shift << "X=>input_shifter,S=>shift_offset,padBit=>pad";
		outmap_shift << "R=>shifted_posit,Sticky=>stkBit";

		newInstance("Shifter",
					"RegimeGenerator",
					param_shift.str(),
					inmap_shift.str(),
					outmap_shift.str());

		//=========================================================================|
		addFullComment("Round to nearest even");
		// ========================================================================|
		vhdl << tab << declare("lsb", 1, false) << " <= shifted_posit" << of(1) << ";" << endl;
		vhdl << tab << declare("rnd", 1, false) << " <= shifted_posit" << of(0) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "stk", 1, false) << " <= stkBit OR Sticky;" << endl;

		// Check for regime underflow & round
		vhdl << tab << declare(getTarget()->logicDelay(4), "round_r", 1, false) << " <= rnd AND (lsb OR stk OR reg_ovf);" << endl;
		vhdl << tab << declare(getTarget()->adderDelay(width_ - 1), "rounded_p", width_ - 1) << " <= shifted_posit" << range(width_ - 1, 1) << " + round_r;" << endl;

		//=========================================================================|
		addFullComment("Check sign & Special cases");
		// ========================================================================|
		addComment("Two's complement if posit is negative");
		vhdl << tab << declare("vSign", width_ - 1) << " <= (others => Sign);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2) + getTarget()->adderDelay(width_ - 1), "final_p", width_ - 1) << " <= "
			 << "(vSign XOR rounded_p) + Sign;" << endl;

		// Assign final result
		vhdl << tab << declare(getTarget()->logicDelay(2), "result", width_)
			 // << tab << tab << "'1' & " << zg(width_ - 1) << " when Inf = '1' else " << endl
			 // << tab << tab << zg(width_) << " when Zero = '1' else" << endl
			 << " <= (" << (width_ - 1) << " => Inf, others => '0') "
			 << " when (Zero OR Inf) = '1' else (Sign & final_p);" << endl;

		vhdl << tab << "R <= result;" << endl;

		addFullComment("End of vhdl generation");
	}

	PositEncoder::~PositEncoder() {}

	void PositEncoder::emulate(TestCase *tc) {}

	OperatorPtr PositEncoder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositEncoder(parentOp, target, width, wES);
	}

	void PositEncoder::registerFactory()
	{
		UserInterface::add("PositEncoder", // name
						   "A posit encoder with a single architecture.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositEncoder::parseArguments);
	}

} // namespace flopoco
