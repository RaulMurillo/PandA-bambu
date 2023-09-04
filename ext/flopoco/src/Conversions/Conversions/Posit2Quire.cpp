/*
  Conversion from Posit to Quire format
  
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
#include "Posit2Quire.hpp"

#include "ShiftersEtc/Normalizer.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"


using namespace std;

namespace flopoco{

#define DEBUGVHDL 0

	Posit2Quire::Posit2Quire(OperatorPtr parentOp, Target* target, int width, int wES, int carry) :
        Operator(parentOp, target), width_(width), wES_(wES), C_(carry) {


		setCopyrightString("Raul Murillo (2021)");

        useNumericStd();

		ostringstream name;

        srcFileName="Posit2Quire";

		if (width_ < 3)
		{
			throw std::string("Posit2Quire Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("Posit2Quire Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

        int sizeRegime = intlog2(width_ - 1) + 1;
        wE_ = sizeRegime + wES_;
        wF_ = width_ - 3 - wES_;
        int maxExp = (1 << wES_) * (width_ - 2);
        int minExp = -maxExp;
        wQuire_ = 1 + C_ + 1 + 2 * maxExp - 2 * minExp;

		name<<"Posit2Quire_" << width_ << "_" << wES_ << "_Quire_" << wQuire_;
        setNameWithFreqAndUID(name.str());

		addInput ("X", width_);
		addOutput("R", wQuire_);

        addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of Posit2Quire \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES << " and " << carry);
		REPORT(DEBUG, "debug of Posit2Quire");

        //====================================================================|
        addFullComment("Extract Sign bit");
        //====================================================================|
        vhdl << tab << declare("sgn", 1, false) << " <= X" << of(width_ - 1) << ";" << endl;
        // Special cases
        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_-1), "nzn", 1, false) << " <= '0' when (X" << range(width_ - 2, 0) << " = " << zg(width_-1) << ") else '1';" << endl;

        //====================================================================|
        addFullComment("Count leading zeros/ones of regime & shift it out");
        //====================================================================|
        vhdl << tab << declare("rc", 1, false) << " <= X" << of(width_ - 2) << ";" << endl;
        vhdl << tab << declare("regPosit", width_ - 2) << " <= X" << range(width_ - 3, 0) << ";" << endl;

        Normalizer* lzocs = (Normalizer*)
			newInstance("Normalizer",
                        "LZOCAndShifter",
                        "wX=" + to_string(width_-2) + " wR=" + to_string(width_-2) + " maxShift=" + to_string(width_-2) + " countType=-1",
                        "X=>regPosit, OZb=>rc",
                        "Count=>regLength, R=>shiftedPosit");

        //=========================================================================|
		addFullComment("Determine the scaling factor - regime & exp");
		// ========================================================================|
        int unsignedRegSize = lzocs->getCountWidth();
        vhdl << tab << declare(getTarget()->logicDelay(1), "neg_sf", 1, false) << " <= '0' when (rc /= sgn) else '1';" << endl;
        vhdl << tab << declare(getTarget()->logicDelay(1), "k", unsignedRegSize) << " <= "
             << "regLength when neg_sf='0' else NOT(regLength);" << endl;

		if (wES_ > 0)
        {
			vhdl << tab << declare("exp", wES_) << " <= shiftedPosit" << range(wF_ + wES_ - 1, wF_) << ";" << endl;
            vhdl << tab << declare("sgnVect", wES_) << " <= (others => sgn);" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "actualExp", wES_) << " <= (sgnVect XOR exp);" << endl;
            vhdl << tab << declare("sf", unsignedRegSize+wES_) << " <= k & actualExp;" << endl;
        }
        else
            vhdl << tab << declare("sf", unsignedRegSize) << " <= k;" << endl;

        //====================================================================|
        addFullComment("Shift fraction into corresponding quire format");
        //====================================================================|
        // Sanity check
        if (unsignedRegSize+wES_ != intlog2(maxExp))
            throw std::string("Sorry, but parameters for int to posit converter do not match sizes.");

        // Use a single small Right shifter, adjust sfBiased according to sf sign.
        vhdl << tab << declare(getTarget()->adderDelay(intlog2(maxExp)) + getTarget()->logicDelay(), "sfBiased", intlog2(maxExp)) << " <= "
             << "std_logic_vector(\"" << unsignedBinary(maxExp, intlog2(maxExp)) << "\" - unsigned(sf)) when neg_sf='0' else NOT(sf);" << endl;

        vhdl << tab << declare(getTarget()->logicDelay(1), "paddedFrac", maxExp + 1 + wF_) << " <= NOT(sgn) & shiftedPosit" << range(wF_ - 1, 0) << " & " << zg(maxExp) << ";" << endl;
        newInstance("Shifter",
                    "Frac_RightShifter",
                    "wX=" + to_string(maxExp + 1 + wF_) + " wR=" + to_string(maxExp + 1 + wF_) + " maxShift=" + to_string(maxExp) + " dir=1 inputPadBit=true",
                    "X=>paddedFrac, S=>sfBiased, padBit=>sgn",
                    "R=>fixedPosit");

        vhdl << tab << declare(getTarget()->logicDelay(1), "quirePosit", 2*(maxExp + 1) + wF_) << " <= "
             << "(fixedPosit & " << zg(maxExp + 1) << ") when neg_sf='0' else "
             << "((" << 2*(maxExp + 1) + wF_ - 1 << " downto " << maxExp + 1 + wF_ << " => sgn) & fixedPosit);" << endl;

        vhdl << tab << declare(getTarget()->logicDelay(1), "quire", wQuire_-1) << " <= "
            << "((" << wQuire_-2 << " downto " << 3*maxExp + 1 << " => sgn) & " // generate carry bits & pad with 0's
            << "quirePosit & (" << maxExp-wF_-2 << " downto 0 => '0')) when nzn='1' else (others => '0');" << endl;

        vhdl << tab << "R <= sgn & quire;" << endl;

        addFullComment("End of vhdl generation");
    }

	Posit2Quire::~Posit2Quire() {
	}

    void Posit2Quire::emulate(TestCase * tc){
    }

	OperatorPtr Posit2Quire::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
        int width, wES, carry;
        UserInterface::parseStrictlyPositiveInt(args, "width", &width);
        UserInterface::parsePositiveInt(args, "wES", &wES);
        UserInterface::parseStrictlyPositiveInt(args, "carry", &carry);
        return new Posit2Quire(parentOp, target, width, wES, carry);
	}

	void Posit2Quire::registerFactory()
	{
		UserInterface::add("Posit2Quire", // name
                            "Conversion from posit to quire.",
                            "Conversions",
                            "", //seeAlso
                            "width(int): posit size in bits;\
                             wES(int): posit exponent size in bits;\
                             carry(int): extra bits to savely allow the sum up to 2^carry products",
                            "", // htmldoc
                            Posit2Quire::parseArguments);
	}

} // namespace flopoco
