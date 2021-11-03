/*
  Conversion from Posit to Floating-point format
  
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

#include "Posit2FP.hpp"
#include "ShiftersEtc/Normalizer.hpp"
#include "TestBenches/PositNumber.hpp"
#include "TestBenches/IEEENumber.hpp"

using namespace std;

namespace flopoco
{
	void Posit2FP::computePositWidths(int const widthP, int const esI, int *wE, int *wF)
	{
		*wE = intlog2(widthP - 1) + 1 + esI;
		*wF = widthP - (esI + 3);
	}

	Posit2FP::Posit2FP(Operator *parentOp, Target *target,
					   int widthP,
					   int wESP,
					   int wEF,
					   int wFF) : Operator(parentOp, target), widthP_(widthP), wESP_(wESP), wEF_(wEF), wFF_(wFF)
	{
		setCopyrightString("Raul Murillo, 2021");
		srcFileName = "Posit2FP";

		ostringstream name;
		name << "Posit" << widthP << "_" << wESP << "_to_Float" << wEF << "_" << wFF;
		setNameWithFreqAndUID(name.str());

		if (widthP < 3)
		{
			throw std::string("Posit2FP Constructor : widthP is too small, should be greater than two");
		}
		int freeWidth = widthP - 1;
		if (wESP >= freeWidth)
		{
			//Avoid posits without even one bit of precision
			throw std::string("Posit2FP Constructor : invalid value of wESP");
		}
		computePositWidths(widthP, wESP, &wE_, &wF_);
		if (wEF_ < wE_)
		{
			throw std::string("Posit2FP Constructor : floatES=%d is too small, should be greater", wEF_);
			// TODO
		}
		if (wFF_ < wF_)
		{
			REPORT(DEBUG, "WARNING - Posit2FP Constructor : fraction field size of the output is too small, will round rightmost bits");
			// TODO
		}

		addInput("I", widthP);
		addOutput("O", 1 + wEF_ + wFF_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of Posit2FP \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << widthP << ", " << wESP << ", " << wEF << ", " << wFF);
		REPORT(DEBUG, "debug of Posit2FP");

		//====================================================================|
		addFullComment("Sign and Special cases");
		//====================================================================|
		vhdl << tab << declare(.0, "sign", 1, false) << " <= I" << of(widthP - 1) << ";" << endl;
		vhdl << tab << declare(.0, "encoding", widthP - 1) << " <= I" << range(widthP - 2, 0) << ";" << endl;

		vhdl << tab << declare(target->eqConstComparatorDelay(widthP - 1), "isNZN", 1, false) << " <= "
			 << "'1' when encoding=\"" << string(widthP - 1, '0') << "\" else '0';" << endl;
		vhdl << tab << declare(target->logicDelay(2), "isNAR", 1, false) << " <= sign AND isNZN;" << endl;
		vhdl << tab << declare(target->logicDelay(2), "isZero", 1, false) << " <= NOT(sign) AND isNZN;" << endl;

		//====================================================================|
		addFullComment("Regime extraction");
		//====================================================================|
		vhdl << tab << "with sign select "
			 << declare(target->logicDelay(widthP), "absoluteEncoding", widthP - 1) << " <= " << endl
			 << tab << tab << "encoding when '0'," << endl
			 << tab << tab << "not(encoding) + 1 when '1'," << endl
			 << tab << tab << "\"" << string(widthP - 1, '-') << "\" when others;" << endl;

		vhdl << declare(0., "exponentSign", 1, false) << " <= absoluteEncoding" << of(widthP - 2) << ";" << endl;
		vhdl << declare(0., "encodingTail", widthP - 2, true) << " <= absoluteEncoding" << range(widthP - 3, 0) << ";" << endl;

		ostringstream param, inmap, outmap;
		int wCount = intlog2(widthP - 2);
		param << "wX=" << widthP - 2;
		param << " wR=" << widthP - 2;
		param << " maxShift=" << widthP - 2;
		inmap << "X=>encodingTail,OZb=>exponentSign";
		outmap << "Count=>lzCount,R=>shiftedResult";

		newInstance("Normalizer", "lzoc", param.str(), inmap.str(), outmap.str());

		vhdl << tab << "with exponentSign select "
			 << declare(target->logicDelay(wCount), "rangeExp", wCount + 1, true) << " <= " << endl
			 << tab << tab << "('1' & not (lzCount)) when '0'," << endl
			 << tab << tab << "('0' & lzCount) when '1'," << endl
			 << tab << tab << "\"" << string(wCount + 1, '-') << "\" when others;" << endl;

		//====================================================================|
		addFullComment("Generate biased exponent");
		//====================================================================|
		ostringstream concat, preconcat;

		if (wESP > 0)
		{
			addComment("Extract exponent");
			vhdl << tab << declare(0., "exponentVal", wESP, true)
				 << " <= shiftedResult" << range(wF_ + wESP - 1, wF_) << ";" << endl;

			concat << " & exponentVal";
		}
		if (wEF_ > wE_)
		{
			addComment("Pad exponent");
			vhdl << tab << "with exponentSign select " << declare(target->logicDelay(), "exponentPad", wEF_ - wE_, true) << " <= " << endl
				 << tab << tab << "\"" << string(wEF_ - wE_, '1') << "\" when '0'," << endl
				 << tab << tab << "\"" << string(wEF_ - wE_, '0') << "\" when '1'," << endl
				 << tab << tab << "\"" << string(wEF_ - wE_, '-') << "\" when others;" << endl;

			preconcat << "exponentPad & ";
		}
		uint64_t bias = (1 << (wEF_ - 1)) - 1;
		vhdl << tab << declare(target->adderDelay(wEF_), "exponent", wEF_, true)
			 << " <= (" << preconcat.str() << "rangeExp" << concat.str() << ") + " << bias << ";" << endl;

		//====================================================================|
		addFullComment("Generate mantissa");
		//====================================================================|
		if (wFF_ >= wF_)
		{ // Simple case, pad with 0's
			vhdl << tab << declare(0., "mantissa", wFF_, true) << " <= shiftedResult" << range(wF_ - 1, 0) << " & \"" << string(wFF_ - wF_, '0') << "\";" << endl;
			// Maybe change first bit of mantissa to isNaR to get IEEE NaN (?)
		}
		else
		{ // Need to round mantissa
			int truncBits = wF_ - wFF_;
			vhdl << tab << declare(0., "lsb", 1, false) << " <= shiftedResult" << of(truncBits) << ";" << endl;
			vhdl << tab << declare(0., "guard", 1, false) << " <= shiftedResult" << of(truncBits - 1) << ";" << endl;
			if (truncBits > 1)
			{
				vhdl << tab << declare(target->eqConstComparatorDelay(truncBits - 1), "sticky", 1, false) << " <= '0' when shiftedResult" << range(truncBits - 2, 0) << "=\"" << string(truncBits - 1, '0') << "\" else '1';" << endl;
				vhdl << tab << declare(target->logicDelay(3), "round_bit", 1, false) << " <= guard and (sticky or lsb);" << endl;
			}
			else
			{
				vhdl << tab << declare(target->logicDelay(2), "round_bit", 1, false) << " <= guard and lsb;" << endl;
			}
			vhdl << tab << declare(target->adderDelay(wFF_), "mantissa", wFF_, true) << " <= "
				 << "shiftedResult" << range(wF_ - 1, truncBits) << " + round_bit;" << endl;
		}

		addComment("Final result");
		vhdl << tab << declare(target->logicDelay(2), "finalSignExp", wEF_ + 1, true) << " <= "
			 << "(finalSignExp'high=>'0', others => isNAR) when (isZero OR isNAR) = '1' else sign & exponent;" << endl;
		vhdl << tab << "O <= (finalSignExp & mantissa);" << endl;
	}

	void Posit2FP::emulate(TestCase *tc)
	{
		mpz_class si = tc->getInputValue("I");
		PositNumber posit(widthP_, wESP_, si);
		mpfr_t val;
		mpfr_init2(val, widthP_ - 2);
		posit.getMPFR(val);
		IEEENumber num(wE_, wF_, val);
		tc->addExpectedOutput("O", num.getSignalValue());
		mpfr_clear(val);
	}

	OperatorPtr Posit2FP::parseArguments(Operator *parentOp, Target *target,
										 std::vector<std::string> &args)
	{
		int widthP, wESP;
		int wEF, wFF;
		UserInterface::parseStrictlyPositiveInt(args, "width", &widthP);
		UserInterface::parsePositiveInt(args, "es", &wESP);
		UserInterface::parseStrictlyPositiveInt(args, "floatES", &wEF);
		UserInterface::parseStrictlyPositiveInt(args, "floatF", &wFF);
		return new Posit2FP(parentOp, target, widthP, wESP, wEF, wFF);
	}

	void Posit2FP::registerFactory()
	{
		UserInterface::add("Posit2FP",
						   "Convert Posit to floating point",
						   "Conversions",
						   "",
						   "width(int): total size of the encoding;\
							es(int): posit exponent field length;\
							floatES(int): float exponent field length;\
							floatF(int): float fraction field length;",
						   "",
						   Posit2FP::parseArguments);
	}

} //namespace
