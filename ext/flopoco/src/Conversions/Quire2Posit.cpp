/*
  Conversion from Quire to Posit format
  
  Author:  Raul Murillo

  This file is part of the FloPoCo project
  
  Initial software.
  Copyright © UCM,  
  2021.
  All rights reserved.

*/

// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
// general c++ library
// #include <vector>
// #include <math.h>
// #include <string.h>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

// include the header of the Operator
#include "Quire2Posit.hpp"

#include "Posit/Encoder/PositFastEncoder.hpp"
#include "ShiftersEtc/Normalizer.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"


using namespace std;

namespace flopoco{

#define DEBUGVHDL 0

	Quire2Posit::Quire2Posit(OperatorPtr parentOp, Target* target, int width, int wES, int carry) :
        Operator(parentOp, target), width_(width), wES_(wES), C_(carry) {


		setCopyrightString("Raul Murillo (2021)");

        // useNumericStd();

		ostringstream name;

        srcFileName="Quire2Posit";

		if (width_ < 3)
		{
			throw std::string("Quire2Posit Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("Quire2Posit Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

        int sizeRegime = intlog2(width_ - 1) + 1;
        wE_ = sizeRegime + wES_;
        wF_ = width_ - 3 - wES_;
        int maxExp = (1 << wES_) * (width_ - 2);
        int minExp = -maxExp;
        wQuire_ = 1 + C_ + 1 + 2 * maxExp - 2 * minExp;

		name<<"Quire2Posit_" << width_ << "_" << wES_ << "_Quire_" << wQuire_;
        setNameWithFreqAndUID(name.str());

		addInput ("X", wQuire_);
		addOutput("R", width_);

        addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of Quire2Posit \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES << " and " << carry);
		REPORT(DEBUG, "debug of Quire2Posit");

        //====================================================================|
        addFullComment("Extract Sign bit");
        //====================================================================|
        vhdl << tab << declare("sgn", 1, false) << " <= X" << of(wQuire_ - 1) << ";" << endl;
        // Special cases
        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(maxExp), "stkTmp", 1, false) << " <= "
             << "'0' when (X" << range(maxExp - 1, 0) << " = " << zg(maxExp) << ") else '1';" << endl;
        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wQuire_-1), "nzn", 1, false) << " <= "
             << "'0' when (X" << range(wQuire_ - 2, 0) << " = " << zg(wQuire_-1) << ") else '1';" << endl;

        //====================================================================|
        addFullComment("Check for overflow");
        //====================================================================|
        vhdl << tab << declare("carryCheck", C_ + maxExp) << " <= (others => sgn);" << endl;
        vhdl << tab << declare(getTarget()->eqComparatorDelay(C_ + maxExp), "ovf", 1, false) << " <= "
             << "'0' when (X" << range(wQuire_ - 2, wQuire_ - 1 - (C_ + maxExp)) << " = carryCheck) else '1';" << endl;

        //====================================================================|
        addFullComment("Count leading zeros/ones & extract fraction");
        //====================================================================|
        int pRange = 2*maxExp + 1;
        vhdl << tab << declare("positRange", pRange) << " <= X" << range(wQuire_ - 1 - (C_ + maxExp) - 1, wQuire_ - 1 - (C_ + maxExp) - pRange) << ";" << endl;

        Normalizer *lzocs = (Normalizer *)
            newInstance("Normalizer",
                        "LZOCAndShifter",
                        "wX=" + to_string(pRange) + " wR=" + to_string(pRange) + " maxShift=" + to_string(pRange) + " countType=-1",
                        "X=>positRange, OZb=>sgn",
                        "Count=>intExp, R=>tmpFrac");

        // First bit of tmpFrac is not from final fraction
        vhdl << tab << declare("frac", wF_) << " <= tmpFrac" << range(pRange - 2, pRange - wF_ - 1) << ";" << endl;
        vhdl << tab << declare("rnd", 1, false) << " <= tmpFrac" << of(pRange - wF_ - 2) << ";" << endl;
        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(pRange - wF_ - 2), "stkBit", 1, false) << " <= "
             << "'0' when (tmpFrac" << range(pRange - wF_ - 3, 0) << " = " << zg(pRange - wF_ - 2) << ") else '1';" << endl;
        vhdl << tab << declare(getTarget()->logicDelay(2), "stk", 1, false) << " <= stkBit OR stkTmp;" << endl;

        //=========================================================================|
        addFullComment("Determine the scaling factor - regime & exp");
        //=========================================================================|
        vhdl << tab << declare("maxExp", wE_+1) << " <= \"" << unsignedBinary(maxExp, wE_+1) << "\";" << endl;
        vhdl << tab << declare(getTarget()->adderDelay(wE_+1), "sfTmp", wE_+1) << " <= "
             << "maxExp - (" << zg(wE_ + 1 - lzocs->getCountWidth()) << " & intExp);" << endl;

        vhdl << tab << declare(getTarget()->logicDelay(), "sf", wE_+1) << " <= "
             << "sfTmp when ovf='0' else maxExp;" << endl;

        //====================================================================|
        addFullComment("Generate final posit");
        //====================================================================|
        newInstance("PositFastEncoder",
                    "PositEncoder",
                    "width=" + to_string(width_) + " wES=" + to_string(wES_),
                    "Sign=>sgn, SF=>sf, Frac=>frac, Guard=>rnd, Sticky=>stk, NZN=>nzn",
                    "R=>R");

        addFullComment("End of vhdl generation");
    }

	Quire2Posit::~Quire2Posit() {
	}

    void Quire2Posit::emulate(TestCase * tc){
    }

	OperatorPtr Quire2Posit::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
        int width, wES, carry;
        UserInterface::parseStrictlyPositiveInt(args, "width", &width);
        UserInterface::parsePositiveInt(args, "wES", &wES);
        UserInterface::parseStrictlyPositiveInt(args, "carry", &carry);
        return new Quire2Posit(parentOp, target, width, wES, carry);
	}

	void Quire2Posit::registerFactory()
	{
		UserInterface::add("Quire2Posit", // name
                            "Conversion from posit to quire.",
                            "Conversions",
                            "", //seeAlso
                            "width(int): posit size in bits;\
                             wES(int): posit exponent size in bits;\
                             carry(int): extra bits to savely allow the sum up to 2^carry products",
                            "", // htmldoc
                            Quire2Posit::parseArguments);
	}

} // namespace flopoco
