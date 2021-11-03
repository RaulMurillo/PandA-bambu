/*
  Conversion from Floating-point to Posit format
  
  Author:  Raul Murillo

  This file is part of the FloPoCo project
  
  Initial software.
  Copyright © UCM, 
  2021.
  All rights reserved.

*/
#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FP2Posit.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "TestBenches/PositNumber.hpp"
#include "TestBenches/IEEENumber.hpp"

using namespace std;
namespace flopoco
{
	void FP2Posit::computeFP2PositWidths(int const widthP, int const wESP, int *wE, int *wF)
	{
		*wE = intlog2(widthP) + 1 + wESP;
		*wF = widthP - (wESP + 3);
	}

	FP2Posit::FP2Posit(Operator *parentOp, Target *target,
					   int wEF,
					   int wFF,
					   int widthP,
					   int wESP) : Operator(parentOp, target), wEF_(wEF), wFF_(wFF), widthP_(widthP), wESP_(wESP)
	{
		setCopyrightString("Raul Murillo, 2021");
		srcFileName = "FP2Posit";

		ostringstream name;
		name << "Float" << wEF_ << "_" << wFF_ << "_to_Posit" << widthP_ << "_" << wESP_;
		setNameWithFreqAndUID(name.str());

		// SET UP THE IO SIGNALS

		if (widthP < 3)
		{
			throw std::string("FP2Posit Constructor : widthP is too small, should be greater than two");
		}
		int freeWidth = widthP - 1;
		if (wESP >= freeWidth)
		{
			//Avoid posits without even one bit of precision
			throw std::string("FP2Posit Constructor : invalid value of wESP");
		}
		computeFP2PositWidths(widthP_, wESP_, &wE_, &wF_);
		if (wE_ < wEF_)
		{
			// Possible issues
			REPORT(DEBUG, "WARNING - FP2Posit Constructor : exponent field length of the output is too small, possible range issues");
		}
		if (wF_ < wFF_)
		{
			REPORT(DEBUG, "WARNING - FP2Posit Constructor : fraction field length of the output is too small, will round rightmost bits");
		}

		int widthI = 1 + wEF_ + wFF_;
		addInput("I", widthI);
		addOutput("O", widthP_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of FP2Posit \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << wEF << ", " << wFF << ", " << widthP << ", " << wESP);
		REPORT(DEBUG, "debug of FP2Posit");

		//====================================================================|
		addFullComment("Split sign, exponent & fraction");
		//====================================================================|
		vhdl << tab << declare(.0, "sign", 1, false) << " <= I" << of(widthI - 1) << ";" << endl;

		// vhdl << "with sign select " <<
		// 	declare(target->logicDelay(widthI), "absI", widthI - 1) << " <= " << endl <<
		// 	tab <<  "I" << range(widthI - 2, 0) << " when '0'," << endl <<
		// 	tab << "not(I" << range(widthI - 2, 0) << ") + 1 when '1'," << endl <<
		// 	tab << "\"" << string(widthI - 1, '-') << "\" when others;" << endl;

		vhdl << tab << declare(.0, "biased_exp", wEF_) << " <= I" << range(widthI - 2, wFF_) << ";" << endl;

		vhdl << tab << declare(.0, "fraction", wF_ + 1) << " <= "; // Include round bit
		if (wF_ + 1 < wFF)
		{ // Need to round mantissa
			int trunc_bits = wFF - (wF_ + 1);
			vhdl << "I" << range(wFF - 1, trunc_bits) << ";" << endl;
			vhdl << tab << declare(target->eqConstComparatorDelay(trunc_bits), "trunc_sticky", 1, false) << " <= "
				 << "'0' when I" << range(trunc_bits - 1, 0) << " = \"" << string(trunc_bits, '0') << "\" else '1';" << endl;
		}
		else
		{ // Simple case, pad with 0's
			vhdl << "I" << range(wFF - 1, 0) << " & \"" << string((wF_ + 1) - wFF, '0') << "\";" << endl;
		}

		addComment("Special cases");
		// // Consider subnormal numbers
		// vhdl << declare(target->eqConstComparatorDelay(wFF),"is_special", 1, false) << " <= '1' when I" << range(wFF-1, 0) << " = \"" << string(wFF,'0') << "\" else '0';" << endl;
		// vhdl << declare(target->eqConstComparatorDelay(wEF),"zero_sub", 1, false) << "<= '1' when biased_exp = \"" << string(wEF,'0') << "\" else '0';" << endl;
		// vhdl << declare(target->logicDelay(2),"is_zero", 1, false) << "<= is_special AND zero_sub;" << endl;
		// vhdl << declare(target->logicDelay(2),"is_sub", 1, false) << "<= NOT(is_special) AND zero_sub;" << endl;
		// vhdl << declare(target->eqConstComparatorDelay(wEF),"is_NaN", 1, false) << "<= '1' when biased_exp = \"" << string(wEF,'1') << "\" else '0';" << endl;

		// Do not consider subnormal numbers
		vhdl << tab << declare(target->eqConstComparatorDelay(wEF), "is_zero", 1, false) << "<= '1' when biased_exp = \"" << string(wEF, '0') << "\" else '0';" << endl;
		vhdl << tab << declare(target->eqConstComparatorDelay(wEF), "is_NaN", 1, false) << "<= '1' when biased_exp = \"" << string(wEF, '1') << "\" else '0';" << endl;

		addComment("Compute unbiased exponent");
		const uint64_t bias = (1 << (wEF - 1)) - 1;
		ostringstream concat;
		if (wE_ < wEF)
		{ // FP has larger dynamic range than Posit format
			vhdl << tab << declare(target->adderDelay(wEF), "unbiased_exp", wEF) << " <= biased_exp - " << bias << ";" << endl;
			// If a 32-bit number is representable with 16 bits, i.e., in the range [-32768,+32767] (32768 = 2^15), then the 32-15=17 msbs will all be the same.
			int esMSB = wEF - wE_;
			vhdl << tab << declare(target->eqConstComparatorDelay(esMSB + 1), "notTooBig", 1, false) << " <= '1' when unbiased_exp" << range(wEF - 1, wE_ - 1) << " = \"" << string(esMSB + 1, '0') << "\" else '0';" << endl;
			vhdl << tab << declare(target->eqConstComparatorDelay(esMSB + 1), "notTooSmall", 1, false) << " <= '1' when unbiased_exp" << range(wEF - 1, wE_ - 1) << " = \"" << string(esMSB + 1, '1') << "\" else '0';" << endl;
			vhdl << tab << declare(target->logicDelay(2), "expFit", 1, false) << " <= notTooBig or notTooSmall;" << endl;

			vhdl << tab << "with expFit select " << declare(target->logicDelay(wE_ + wF_ + 1), "exponent", wE_) << " <= " << endl
				 << tab << tab << "unbiased_exp" << range(wE_ - 1, 0) << " when '1'," << endl
				 << tab << tab << "unbiased_exp(unbiased_exp'high) & \"" << string(wE_ - 1, '1') << "\" when '0'," << endl
				 << tab << tab << "\"" << string(wE_, '-') << "\" when others;" << endl;
		}
		else
		{ // Pad exponent
			concat << "\"" << string(wE_ - wEF, '0') << "\" & ";
			vhdl << tab << declare(target->adderDelay(wE_), "exponent", wE_) << " <= (" << concat.str() << "biased_exp) - " << bias << ";" << endl;
		}

		if (wESP > 0)
		{
			vhdl << tab << declare(0., "partial_exponent_us", wESP) << " <= exponent" << range(wESP - 1, 0) << ";" << endl;
		}

		//====================================================================|
		addFullComment("Generate regime field");
		//====================================================================|
		// how much to shift
		int wCount = intlog2(widthP);
		vhdl << tab << declare(0., "bin_regime", wCount) << " <= exponent" << range(wE_ - 2, wESP) << ";" << endl;
		vhdl << tab << declare(0., "first_regime", 1, false) << " <= exponent" << of(wE_ - 1) << ";" << endl;

		vhdl << tab << "with first_regime select " << declare(target->logicDelay(wCount), "regime", wCount) << " <= " << endl
			 << tab << tab << "bin_regime when '0', " << endl
			 << tab << tab << "not bin_regime when '1', " << endl
			 << tab << tab << "\"" << string(wCount, '-') << "\" when others;" << endl;

		// what to start with (negative or positive exp)
		vhdl << tab << declare(target->logicDelay(1), "pad_bit", 1, false) << " <= not(first_regime);" << endl;
		vhdl << tab << declare(target->logicDelay(1), "start_regime", 2) << " <= not(first_regime) & first_regime;" << endl;

		// Sanity check
		if (widthP != 2 + wESP + wF_ + 1)
			throw std::string("Sorry, but parameters for FP2Posit do not match sizes.");
		vhdl << tab << declare(0., "input_shifter", widthP) << " <= start_regime ";
		if (wESP > 0)
		{
			vhdl << "& partial_exponent_us ";
		}
		vhdl << "& fraction;" << endl;

		ostringstream param, inmap, outmap;
		param << "wX=" << widthP;
		param << " wR=" << widthP;
		param << " maxShift=" << widthP;
		param << " dir=" << Shifter::Right;
		param << " computeSticky=true";
		param << " inputPadBit=true";
		inmap << "X=>input_shifter,S=>regime,padBit=>pad_bit";
		outmap << "R=>extended_posit,Sticky=>pre_sticky";

		newInstance("Shifter", "rshift", param.str(), inmap.str(), outmap.str());

		//====================================================================|
		addFullComment("Round shifted result");
		//====================================================================|
		vhdl << tab << declare(0., "truncated_posit", widthP - 1) << " <= extended_posit" << range(widthP - 1, 1) << ";" << endl;
		vhdl << tab << declare(0., "lsb", 1, false) << " <= extended_posit" << of(1) << ";" << endl;
		vhdl << tab << declare(0., "guard", 1, false) << " <= extended_posit" << of(0) << ";" << endl;
		if (wF_ + 1 < wFF)
		{
			vhdl << tab << declare(target->logicDelay(2), "sticky", 1, false) << " <= trunc_sticky OR pre_sticky;" << endl;
		}
		else
		{
			vhdl << tab << declare(0., "sticky", 1, false) << " <= pre_sticky;" << endl;
		}

		vhdl << tab << declare(target->logicDelay(3), "round_bit", 1, false) << " <= guard and (sticky or lsb);" << endl;

		vhdl << tab << declare(target->adderDelay(widthP - 1), "rounded_reg_exp_frac", widthP - 1) << " <= truncated_posit + round_bit;" << endl;

		vhdl << tab << "with sign select " << declare(target->logicDelay(widthP - 1), "rounded_posit", widthP) << " <= " << endl
			 << tab << tab << "sign & rounded_reg_exp_frac when '0'," << endl
			 << tab << tab << "sign & not(rounded_reg_exp_frac) + 1 when '1'," << endl
			 << tab << tab << "\"" << string(widthP, '-') << "\" when others;" << endl;

		vhdl << tab << declare(target->logicDelay(2), "rounded_posit_zero", widthP) << " <= rounded_posit when (is_zero or is_NaN) = '0' else (is_NaN & \"" << string(widthP - 1, '0') << "\");" << endl;
		vhdl << tab << "O <= rounded_posit_zero;" << endl;

		addFullComment("End of vhdl generation");
	};

	void FP2Posit::emulate(TestCase *tc)
	{
		// mpz_class si = tc->getInputValue("I");
		// IEEENumber num(wEF_, wFF_, si);
		// mpfr_t val;
		// mpfr_init2(val, wFF_ + 1);
		// num.getMPFR(val);

		// // This yields to error due to dynamic range overflow
		// PositNumber posit(widthP_, wESP_, val);
		// tc->addExpectedOutput("O", posit.getSignalValue());
		// mpfr_clear(val);
	}

	OperatorPtr FP2Posit::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int widthP, wESP;
		int wEF, wFF;
		UserInterface::parseStrictlyPositiveInt(args, "width", &widthP);
		UserInterface::parseInt(args, "es", &wESP);
		UserInterface::parseStrictlyPositiveInt(args, "floatES", &wEF);
		UserInterface::parseStrictlyPositiveInt(args, "floatF", &wFF);
		return new FP2Posit(parentOp, target, wEF, wFF, widthP, wESP);
	}

	void FP2Posit::registerFactory()
	{
		UserInterface::add("FP2Posit",								   // name
						   "Converts Floating Point Format to Posits", // description, string
						   "Conversions",							   // category, from the list defined in UserInterface.cpp
						   "",										   //seeAlso
						   "floatES(int): float exponent field length;\
						   floatF(int): float fraction field length;\
						   width(int): The size of the posit;\
						   es(int): The exponent size (for the posit);",
						   "",
						   FP2Posit::parseArguments);
	}

} //namespace
