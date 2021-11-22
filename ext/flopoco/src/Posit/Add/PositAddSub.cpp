/*
  A Posit adder/subtractor for FloPoCo
  
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

#include "PositAddSub.hpp"
#include "Posit/Decoder/PositDecoder.hpp"
#include "Posit/Encoder/PositEncoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "ShiftersEtc/Normalizer.hpp"
// #include "IntAddSubCmp/IntAdder.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositAddSub::PositAddSub(OperatorPtr parentOp, Target *target, int width, int wES, int Sub) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositAddSub";

		if (width_ < 3)
		{
			throw std::string("PositAddSub Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositAddSub Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------
		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositAddSub_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositAddSub \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositAddSub");

		//=========================================================================|
		addFullComment("Decode X & Y operands");
		// ========================================================================|
		newInstance("PositDecoder",
					"X_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X",
					"Sign=>sign_X,SF=>sf_X,Frac=>frac_X,Zero=>z_X,Inf=>inf_X,Abs_in=>abs_X");

		newInstance("PositDecoder",
					"Y_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>Y",
					"Sign=>sign_Y,SF=>sf_Y,Frac=>frac_Y,Zero=>z_Y,Inf=>inf_Y,Abs_in=>abs_Y");

		//=========================================================================|
		addFullComment("Sign and Special cases computation");
		// ========================================================================|
		if(Sub)
		{
			vhdl << tab << declare(getTarget()->logicDelay(3), "OP", 1, false) << " <= NOT(sign_X XOR sign_Y);" << endl;
		}
		else
		{
			vhdl << tab << declare(getTarget()->logicDelay(2), "OP", 1, false) << " <= (sign_X XOR sign_Y);" << endl;
		}

		vhdl << tab << declare(getTarget()->logicDelay(2), "inf", 1, false) << " <= inf_X OR inf_Y;" << endl;
		// No need to check if any zero operand

		//=========================================================================|
		addFullComment("Compare operands and adjust values");
		// ========================================================================|
		vhdl << tab << declare(getTarget()->eqComparatorDelay(width - 1), "is_larger", 1, false) << " <= '1' when abs_X > abs_Y else '0';" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(1), "larger_sign", 1, false) << " <= " << endl
			 << tab << tab << "sign_X when '1'," << endl
			 << tab << tab << "sign_Y when '0'," << endl
			 << tab << tab << "'-' when others;" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(1), "larger_sf", wE_) << " <= " << endl
			 << tab << tab << "sf_X when '1'," << endl
			 << tab << tab << "sf_Y when '0'," << endl
			 << tab << tab << "\"" << string(wE_, '-') << "\" when others;" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(1), "larger_frac", wF_ + 1) << " <= " << endl
			 << tab << tab << "((NOT z_X) & frac_X) when '1'," << endl
			 << tab << tab << "((NOT z_Y) & frac_Y) when '0'," << endl
			 << tab << tab << "\"" << string(wF_ + 1, '-') << "\" when others;" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(1), "smaller_frac", wF_ + 3) << " <= " << endl
			 << tab << tab << "((NOT z_Y) & frac_Y & \"00\")  when '1'," << endl
			 << tab << tab << "((NOT z_X) & frac_X & \"00\")  when '0'," << endl
			 << tab << tab << "\"" << string(wF_ + 3, '-') << "\" when others;" << endl;

		// Compute smaller fraction offset as the absolute value of Scaling Factors' difference, (and saturate it)
		addComment("Exponents difference");
		vhdl << tab << declare(getTarget()->adderDelay(wE_), "sf_diff", wE_) << " <= sf_X - sf_Y;" << endl;
		vhdl << tab << declare("sf_sign", 1, false) << " <= sf_diff(sf_diff'high);" << endl;
		vhdl << tab << declare(.0, "sf_signVect", wE_) << " <= (others => sf_sign);" << endl;
		vhdl << tab << declare(getTarget()->adderDelay(wE_) + getTarget()->logicDelay(), "offset", wE_) << " <= (sf_signVect XOR sf_diff) + sf_sign;" << endl;

		int maxshiftsize = intlog2(wF_ + 3);
		if ((wE_) > maxshiftsize)
		{
			vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wE_ - maxshiftsize + 1), "shift_saturate", 1, false) << " <= "
				 << "'0' when (offset" << range(wE_ - 1, maxshiftsize) << " = 0) else '1';" << endl;

			vhdl << tab << declare(getTarget()->logicDelay(1), "frac_offset", maxshiftsize) << " <= "
				 << "CONV_STD_LOGIC_VECTOR(" << wF_ + 3 << "," << maxshiftsize << ") when shift_saturate = '1' else "
				 << "offset" << range(maxshiftsize - 1, 0) << ";" << endl;
		}
		else
		{ // Only for wES==0. In this case (wE_) == maxshiftsize
			if ((wE_) != maxshiftsize)
			{ // Sanity check
				// Unreachable point
				throw std::string("Sorry, but parameters for posit encoder do not match sizes.");
				// cerr << "Wrong sizes: " << (wE_) << " and " << maxshiftsize << endl;
			}
			vhdl << tab << declare("frac_offset", maxshiftsize) << " <= offset;" << endl;
		}

		//=========================================================================|
		addFullComment("Align fractions");
		// ========================================================================|
		ostringstream param1, inmap1, outmap1;
		param1 << "wX=" << wF_ + 3;
		param1 << " maxShift=" << wF_ + 3; // Adjust this. Is there a lower maximum length to shift?
		param1 << " wR=" << wF_ + 3;
		param1 << " dir=" << Shifter::Right;
		param1 << " computeSticky=true";
		//param1 << " inputPadBit=true";

		//inmap1 << "X=>smaller_frac,S=>frac_offset,padBit=>OP";
		inmap1 << "X=>smaller_frac,S=>frac_offset";

		outmap1 << "R=>shifted_frac,Sticky=>sticky";

		newInstance("Shifter",
					"RightShifterFraction",
					param1.str(),
					inmap1.str(),
					outmap1.str());

		//=========================================================================|
		addFullComment("Add aligned fractions");
		//=========================================================================|
		vhdl << tab << "with OP select " << declare(getTarget()->adderDelay(wF_ + 5), "add_frac", wF_ + 5) << " <= " << endl
			 << tab << tab << "('0' & larger_frac & \"000\") + ('0' & shifted_frac & sticky) when '0'," << endl
			 << tab << tab << "('0' & larger_frac & \"000\") - ('0' & shifted_frac & sticky) when '1'," << endl
			 << tab << tab << "\"" << string(wF_ + 5, '-') << "\" when others;" << endl;

		addComment("Normalization of add_frac");
		ostringstream param2, inmap2, outmap2;
		int wCount = intlog2(wF_ + 5);
		param2 << "wX=" << wF_ + 5;   //??
		param2 << " wR=" << wF_ + 5; //??
		param2 << " maxShift=" << wF_ + 5;
		param2 << " countType=" << 0;

		inmap2 << "X=>add_frac";
		outmap2 << "Count=>nZerosNew,R=>normFrac";

		Normalizer *lzocs = (Normalizer *)
			newInstance("Normalizer",
						"FractionNormalizer",
						param2.str(),
						inmap2.str(),
						outmap2.str());

		addComment("Adjust exponent");
		vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "sf_tmp", wE_ + 1) << "<= (larger_sf(larger_sf'high) & larger_sf) + '1';" << endl;
		vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "sf_final", wE_ + 1) << "<= sf_tmp - (" << zg(wE_ + 1 - lzocs->getCountWidth()) << " & nZerosNew);" << endl;

		//=========================================================================|
		addFullComment("Data Rounding & Encoding");
		// ========================================================================|
		vhdl << tab << declare("frac", wF_) << " <= normFrac" << range(wF_ + 3, 4) << ";" << endl;
		vhdl << tab << declare("rnd", 1, false) << "<= normFrac" << of(3) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(3), "stk", 1, false) << "<= normFrac" << of(2) << " or normFrac" << of(1) << " or normFrac" << of(0) << ";" << endl;

		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(lzocs->getCountWidth()), "z", 1, false) << " <= '1' when nZerosNew = " << og(lzocs->getCountWidth(), 0) << " else '0';" << endl;

		newInstance("PositEncoder",
					"R_encoding",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"Sign=>larger_sign, SF=>sf_final, Frac=>frac, Round=>rnd, Sticky=>stk, Zero=>z, Inf=>inf",
					"R=>R");

		addFullComment("End of vhdl generation");
	};

	PositAddSub::~PositAddSub() {}

	void PositAddSub::emulate(TestCase *tc) {}

	OperatorPtr PositAddSub::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES, Sub;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		UserInterface::parsePositiveInt(args, "Sub", &Sub);
		return new PositAddSub(parentOp, target, width, wES, Sub);
	}

	void PositAddSub::registerFactory()
	{
		UserInterface::add("PositAddSub", // name
						   "A posit adder/subtractor with a single architecture.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
						   wES(int): posit exponent size in bits;\
						   Sub(int): indicates if perform addition (0) or subtraction (1)",
						   "", // htmldoc
						   PositAddSub::parseArguments);
	}

} // namespace flopoco
