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

#include "PositFastDecoder.hpp"
#include "ShiftersEtc/Normalizer.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositFastDecoder::PositFastDecoder(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021-2022)");
		ostringstream name;
		srcFileName = "PositFastDecoder";

		if (width_ < 3)
		{
			throw std::string("PositFastDecoder Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositFastDecoder Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int regSize = intlog2(width_ - 1) + 1;
		wE_ = regSize + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositFastDecoder_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addOutput("Sign");
		addOutput("SF", wE_);
		addOutput("Frac", wF_);
		addOutput("NZN");

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositFastDecoder \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositFastDecoder");

		//====================================================================|
		addFullComment("Sign bit & special cases");
		//====================================================================|
		vhdl << tab << declare("sgn", 1, false) << " <= X" << of(width_ - 1) << ";" << endl;
		// Neither Zero nor NaR - This is 1 when the posit encodes a normal value
		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_ - 1), "pNZN", 1, false) << " <= '0' when (X" << range(width_ - 2, 0) << " = " << zg(width_ - 1) << ") else '1';" << endl;

		//====================================================================|
		addFullComment("Count leading zeros/ones of regime & shift it out");
		//====================================================================|
		vhdl << tab << declare("rc", 1, false) << " <= X" << of(width_ - 2) << ";" << endl;
		vhdl << tab << declare("regPosit", width_ - 2) << " <= X" << range(width_ - 3, 0) << ";" << endl;

		Normalizer *lzocs = (Normalizer *)
			newInstance("Normalizer",
						"RegimeCounter",
						"wX=" + to_string(width_ - 2) + " wR=" + to_string(width_ - 2) + " maxShift=" + to_string(width_ - 2) + " countType=-1",
						"X=>regPosit, OZb=>rc",
						"Count=>regLength, R=>shiftedPosit");

		//=========================================================================|
		addFullComment("Determine the scaling factor - regime & exp");
		// ========================================================================|
		vhdl << tab << declare(getTarget()->logicDelay(1), "k", regSize) << " <= "
			 << zg(regSize - lzocs->getCountWidth()) << " & regLength when rc /= sgn else " // k = len(reg)-1
			 << og(regSize - lzocs->getCountWidth()) << " & NOT(regLength);" << endl;		// k = -len(reg)

		if (wES_ > 0)
		{
			// vhdl << tab << declare("exp", wES_) << " <=  shiftedPosit" << range(wF_ + wES_ - 1, wF_) << ";" << endl;
			vhdl << tab << declare("sgnVect", wES_) << " <= (others => sgn);" << endl;
			// vhdl << tab << declare(getTarget()->logicDelay(2), "actualExp", wES_) << " <=  (sgnVect XOR exp);" << endl;
			vhdl << tab << declare(getTarget()->logicDelay(2), "exp", wES_) << " <= shiftedPosit" << range(wF_ + wES_ - 1, wF_) << " XOR sgnVect;" << endl;
		}
		vhdl << tab << declare("pSF", wE_) << " <= k";
		if (wES_ > 0)
		{
			// vhdl << " & actualExp";
			vhdl << " & exp";
		}
		vhdl << ";" << endl;

		//=========================================================================|
		addFullComment("Extract fraction");
		// ========================================================================|
		// Discard negated regime bit and exponent bits (if any)
		vhdl << tab << declare("pFrac", wF_) << " <= shiftedPosit" << range(wF_ - 1, 0) << ";" << endl;

		// Prepare outputs
		vhdl << tab << "Sign <= sgn;" << endl;
		vhdl << tab << "SF <= pSF;" << endl;
		vhdl << tab << "Frac <= pFrac;" << endl;
		vhdl << tab << "NZN <= pNZN;" << endl;

		addFullComment("End of vhdl generation");
	}

	PositFastDecoder::~PositFastDecoder() {}

	void PositFastDecoder::emulate(TestCase *tc) {}

	OperatorPtr PositFastDecoder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositFastDecoder(parentOp, target, width, wES);
	}

	void PositFastDecoder::registerFactory()
	{
		UserInterface::add("PositFastDecoder", // name
						   "A posit decoder with a single architecture.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositFastDecoder::parseArguments);
	}

} // namespace flopoco
