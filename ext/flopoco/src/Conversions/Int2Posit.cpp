/*
  Conversion from Posit to integer format
  
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
#include "Int2Posit.hpp"

#include "ShiftersEtc/Normalizer.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

    Int2Posit::Int2Posit(OperatorPtr parentOp, Target *target, int widthI, int widthO, int wES, bool trunc_p) : Operator(parentOp, target), widthINT_(widthI), width_(widthO), wES_(wES), trunc_p(trunc_p)
    {

        setCopyrightString("Raul Murillo (2021)");

        ostringstream name;

        srcFileName = "Int2Posit";

        name << "Int2Posit_" << widthINT_ << "_to_" << width_ << "_" << wES_ << "_" << (trunc_p == true ? "T" : "NT");
        setNameWithFreqAndUID(name.str());

        if (width_ < 3)
        {
            throw std::string("Int2Posit Constructor: width is too small, should be greater than two");
        }
        if (wES_ >= width_ - 3)
        {
            //Avoid posits without even one bit of precision
            throw std::string("Int2Posit Constructor: invalid value of wES");
        }

        // -------- Parameter set up -----------------

        int sizeRegime = intlog2(width_ - 1) + 1;
        wE_ = sizeRegime + wES_;
        wF_ = width_ - 3 - wES_;

        addInput("X", widthINT_);
        addOutput("R", width_);

        addFullComment("Start of vhdl generation");

        REPORT(INFO, "Declaration of Int2Posit \n");
        REPORT(DETAILED, "this operator has received the following parameters: " << widthI << ", " << widthO << ", " << wES << " and " << trunc_p);
        REPORT(DEBUG, "debug of Int2Posit");

        //====================================================================|
        addFullComment("Extract Sign bit");
        //====================================================================|
        vhdl << tab << declare("sgn", 1, false) << " <= X" << of(widthINT_ - 1) << ";" << endl;
        vhdl << tab << declare("unsgnInt", widthINT_ - 1) << " <= X" << range(widthINT_ - 2, 0) << ";" << endl;

        // Special case - Neither Zero nor NaR
        // This is 1 when input is a normal int value
        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(widthINT_ - 1), "nzn", 1, false) << " <= "
             << "'0' when (unsgnInt = " << zg(widthINT_ - 1) << ") else '1';" << endl;

        //====================================================================|
        addFullComment("Count leading zeros/ones & extract fraction");
        //====================================================================|
        Normalizer *lzocs = (Normalizer *)
            newInstance("Normalizer",
                        "LZOCAndShifter",
                        "wX=" + to_string(widthINT_ - 1) + " wR=" + to_string(widthINT_ - 1) + " maxShift=" + to_string(widthINT_ - 1) + " countType=-1",
                        "X=>unsgnInt, OZb=>sgn",
                        "Count=>intExp, R=>tmpFrac");

        // First bit of tmpFrac is not from final fraction
        int shifterSize;
        if (wF_ < (widthINT_ - 2))
        {
            int fracSize = wF_;
            shifterSize = 2 + wES_ + fracSize;
            vhdl << tab << declare("frac", fracSize) << " <= tmpFrac" << range(widthINT_ - 3, widthINT_ - 2 - fracSize) << ";" << endl;
            if (!trunc_p)
            {
                shifterSize++; // Add round bit
                int stickyBits = widthINT_ - 2 - wF_;
                vhdl << tab << declare("rndBit", 1, false) << " <= tmpFrac" << of(stickyBits - 1) << ";" << endl;
                if (stickyBits > 1)
                {
                    vhdl << tab << declare(getTarget()->eqConstComparatorDelay(stickyBits - 1), "stkBit", 1, false) << " <= "
                         << "'0' when (tmpFrac" << range(stickyBits - 2, 0) << " = " << zg(stickyBits - 1) << ") else '1';" << endl;
                }
            }
        }
        else
        {
            int fracSize = widthINT_ - 2;
            shifterSize = 2 + wES_ + fracSize;
            vhdl << tab << declare("frac", fracSize) << " <= tmpFrac" << range(widthINT_ - 3, 0) << ";" << endl;
            if (!trunc_p)
            {
                shifterSize++; // Add round bit
                vhdl << tab << declare("rndBit", 1, false) << " <= '0';" << endl;
            }
        }
        // Sanity check
        if (shifterSize > width_)
            throw std::string("Sorry, but parameters for int to posit converter do not match sizes.");

        //=========================================================================|
        addFullComment("Determine the scaling factor - regime & exp");
        //=========================================================================|
        vhdl << tab << declare(getTarget()->adderDelay(wE_), "sf", wE_) << " <= "
             << "\"" << unsignedBinary(widthINT_ - 2, wE_) << "\" - (" << zg(wE_ - lzocs->getCountWidth()) << " & intExp);" << endl;
        // Special handling when input int is -1
        vhdl << tab << declare("minusOne", 1, false) << " <= sf(sf'high);" << endl;

        int paddedSize = width_ - int(trunc_p) - 1;
        int offsetSize = intlog2(paddedSize);
        vhdl << tab << declare("reg", offsetSize) << " <= ";
        if ((wE_ - wES_) > offsetSize)
        {
            vhdl << "sf" << range(offsetSize + wES_ - 1, wES_);
        }
        else
        {
            vhdl << zg(offsetSize - (wE_ - wES_)) << " & sf" << range(wE_ - 1, wES_);
        }
        vhdl << ";" << endl;

        if (wES_ > 0)
        {
            vhdl << tab << declare("sgnVect", wES_) << " <= (others => sgn);" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "expBits", wES_) << " <= "
                 << "sf" << range(wES_ - 1, 0) << " XOR sgnVect;" << endl;
        }

        //=========================================================================|
        addFullComment("Shift out fraction according to scaling factor");
        //=========================================================================|

        vhdl << tab << declare(getTarget()->logicDelay(1), "padBit", 1, false) << " <= NOT sgn;" << endl;

        // An extra padBit must be (left) appended after Shifter - this saves some area & power
        vhdl << tab << declare("paddedFrac", paddedSize) << " <= sgn";
        if (wES_ > 0)
        {
            vhdl << " & expBits";
        };
        vhdl << " & frac";
        if (!trunc_p) // Add round bit
        {
            vhdl << " & rndBit";
        }
        if (shifterSize < paddedSize)
        {
            vhdl << " & " << zg(paddedSize - shifterSize);
        }
        vhdl << ";" << endl;

        if (trunc_p)
        {
            newInstance("Shifter",
                        "RightShifterComponent",
                        "wX=" + to_string(paddedSize) + " wR=" + to_string(width_ - 2) + " maxShift=" + to_string(paddedSize) + " dir=1 inputPadBit=true",
                        "X=>paddedFrac, S=>reg, padBit=>padBit",
                        "R=>shiftedFrac");

            vhdl << tab << declare("positNumber", width_ - 1) << " <= padBit & shiftedFrac;" << endl;
        }
        else
        {
            newInstance("Shifter",
                        "RightShifterComponent",
                        "wX=" + to_string(paddedSize) + " wR=" + to_string(width_ - 1) + " maxShift=" + to_string(paddedSize) + " dir=1 inputPadBit=true computeSticky=true",
                        "X=>paddedFrac, S=>reg, padBit=>padBit",
                        "R=>shiftedFrac, Sticky=>sticky");

            vhdl << tab << declare("validFrac", width_ - 1) << " <= padBit & shiftedFrac" << range(width_ - 2, 1) << ";" << endl;
            vhdl << tab << declare("lsb", 1, false) << " <= shiftedFrac" << of(1) << ";" << endl;
            vhdl << tab << declare("rnd", 1, false) << " <= shiftedFrac" << of(0) << ";" << endl;

            if ((wF_ + 1) < (widthINT_ - 2))
            {
                vhdl << tab << declare(getTarget()->logicDelay(4), "round", 1, false) << " <= "
                     << "rnd AND (lsb OR sticky OR stkBit);" << endl;
            }
            else
            {
                vhdl << tab << declare(getTarget()->logicDelay(3), "round", 1, false) << " <= "
                     << "rnd AND (lsb OR sticky);" << endl;
            }

            vhdl << tab << declare(getTarget()->adderDelay(width_ - 1), "positNumber", width_ - 1) << " <= "
                 << "validFrac + round;" << endl;
        }

        vhdl << tab << declare(getTarget()->logicDelay(2), "result", width_ - 1) << " <= "
             << "positNumber when (nzn AND NOT(minusOne)) = '1' else (nzn & " << zg(width_ - 2) << ");" << endl;
        vhdl << tab << "R <= sgn & result;" << endl;

        addFullComment("End of vhdl generation");
    }

    Int2Posit::~Int2Posit() {}

    void Int2Posit::emulate(TestCase *tc) {}

    OperatorPtr Int2Posit::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
    {
        int widthI, widthO, wES;
        bool trunc;
        UserInterface::parseStrictlyPositiveInt(args, "widthI", &widthI);
        UserInterface::parseStrictlyPositiveInt(args, "widthO", &widthO);
        UserInterface::parsePositiveInt(args, "wES", &wES);
        UserInterface::parseBoolean(args, "trunc", &trunc);
        return new Int2Posit(parentOp, target, widthI, widthO, wES, trunc);
    }

    void Int2Posit::registerFactory()
    {
        UserInterface::add("Int2Posit", // name
                           "Conversion from posit to integer.",
                           "Conversions",
                           "", //seeAlso
                           "widthI(int): integer size in bits;\
                            widthO(int): output posit size in bits;\
                            wES(int): posit exponent size in bits;\
                            trunc(bool)=false: true means truncated (cheaper), false means rounded",
                           "", // htmldoc
                           Int2Posit::parseArguments);
    }

} // namespace flopoco
