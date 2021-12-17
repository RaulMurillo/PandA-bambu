/*
  Posit decoder
  
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

#include "PositDecoder.hpp"
#include "ShiftersEtc/Normalizer.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositDecoder::PositDecoder(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositDecoder";

		if (width_ < 3)
		{
			throw std::string("PositDecoder Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositDecoder Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int regSize = intlog2(width_ - 1) + 1;
		wE_ = regSize + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositDecoder_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addOutput("Sign");
		addOutput("SF", wE_);
		addOutput("Frac", wF_);
		addOutput("Zero");
		addOutput("Inf");
		addOutput("Abs_in", width_ - 1);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositDecoder \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositDecoder");

		//====================================================================|
		addFullComment("Sign bit & special cases");
		//====================================================================|
		vhdl << tab << declare(0., "s", 1, false) << " <= X" << of(width_ - 1) << ";" << endl;
		vhdl << tab << declare(0., "remainP", width_ - 1) << " <= X" << range(width_ - 2, 0) << ";" << endl;

		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_ - 1), "special", 1, false) << " <= "
			 << "'1' when (remainP = " << zg(width_ - 1) << ") else '0';" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "is_zero", 1, false) << " <= not(s) AND special;" << endl; // 1 if X is zero
		vhdl << tab << declare(getTarget()->logicDelay(2), "is_NAR", 1, false) << "<= s AND special;" << endl;		  // 1 if X is infinity

		//====================================================================|
		addFullComment("Get absolute value of the Posit");
		//====================================================================|
		vhdl << tab << declare(0., "v_sign", width_ - 1) << " <= (others => s);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2) + getTarget()->adderDelay(width_ - 1), "p_abs", width_ - 1) << " <= (v_sign XOR remainP) + s;" << endl;

		//====================================================================|
		addFullComment("Count leading zeros/ones of regime & shift it out");
		//====================================================================|
		vhdl << tab << declare(0., "rc", 1, false) << " <= p_abs(p_abs'high);" << endl; // Regime check
		vhdl << tab << declare(0., "regPosit", width_ - 2) << " <= p_abs" << range(width_ - 3, 0) << ";" << endl;

		Normalizer *lzocs = (Normalizer *)
			newInstance("Normalizer",
						"RegimeCounter",
						"wX=" + to_string(width_ - 2) + " wR=" + to_string(width_ - 2) + " maxShift=" + to_string(width_ - 2) + " countType=-1",
						"X=>regPosit, OZb=>rc",
						"Count=>regLength, R=>shiftedPosit");

		//=========================================================================|
		addFullComment("Determine the scaling factor - regime & exp");
		// ========================================================================|
		int wCount = lzocs->getCountWidth();
		vhdl << tab << declare(getTarget()->logicDelay(1), "k", regSize) << " <= " << endl
			 << tab << tab << zg(regSize - wCount) << " & regLength when rc = '1' else" << endl
			 << tab << tab << og(regSize - wCount) << " & NOT(regLength);" << endl;

		vhdl << tab << declare(0., "pSF", wE_) << " <= k";
		if (wES_ > 0)
		{
			vhdl << " & shiftedPosit" << range(wF_ + wES_ - 1, wF_);
		}
		vhdl << ";" << endl;

		//=========================================================================|
		addFullComment("Extract fraction");
		// ========================================================================|
		vhdl << tab << declare(0., "pFrac", wF_) << " <= shiftedPosit" << range(wF_ - 1, 0) << ";" << endl;

		addComment("Prepare outputs");
		vhdl << tab << "Sign <= s;" << endl;
		vhdl << tab << "SF <= pSF;" << endl;
		vhdl << tab << "Frac <= pFrac;" << endl;
		vhdl << tab << "Zero <= is_zero;" << endl;
		vhdl << tab << "Inf <= is_NAR;" << endl;
		vhdl << tab << "Abs_in <= p_abs;" << endl;

		addFullComment("End of vhdl generation");
	}

	PositDecoder::~PositDecoder() {}

	void PositDecoder::emulate(TestCase *tc) {}

	OperatorPtr PositDecoder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositDecoder(parentOp, target, width, wES);
	}

	void PositDecoder::registerFactory()
	{
		UserInterface::add("PositDecoder", // name
						   "A posit decoder with a single architecture.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositDecoder::parseArguments);
	}

} // namespace flopoco
