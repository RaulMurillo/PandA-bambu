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
#include "Posit2Int.hpp"

#include "ShiftersEtc/Normalizer.hpp"
#include "ShiftersEtc/Shifters.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"


using namespace std;

namespace flopoco{

#define DEBUGVHDL 0

	Posit2Int::Posit2Int(OperatorPtr parentOp, Target* target, int widthI, int wESI, int widthO, bool trunc_p) :
        Operator(parentOp, target), width_(widthI), wES_(wESI), widthO_(widthO), trunc_p(trunc_p) {


		setCopyrightString("Raul Murillo (2021)");		

		ostringstream name;

        srcFileName="Posit2Int";

		name<<"Posit2Int_" << width_ << "_" << wES_ << "_to_" << widthO_ << "_" << (trunc_p == true ? "T" : "NT");
        setNameWithFreqAndUID(name.str());

		if (width_ < 3)
		{
			throw std::string("Posit2Int Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("Posit2Int Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

        int sizeRegime = intlog2(width_ - 1) + 1;
        wE_ = sizeRegime + wES_;
        wF_ = width_ - 3 - wES_;

        unsigned int positExp = (1 << wES_) * (width_ - 2);
        unsigned int intExp = (widthO_ - 1);
        bool possibleOverflow = (positExp >= intExp); // Max int value is 2**(widthO_ - 1) - 1
        if (possibleOverflow){
            REPORT(INFO, "Need to handle possible overflow. Max posit exponent is " <<  positExp << " >= " << intExp);
        }
        else{
            REPORT(INFO, "NO Need to handle overflow. Max posit exponent is " <<  positExp << " < " << intExp);
        }

		addInput ("X", width_);
		addOutput("R", widthO_);

        addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of Posit2Int \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << widthI << ", " << wESI << ", " << widthO << " and " << trunc_p);
		REPORT(DEBUG, "debug of Posit2Int");

        //====================================================================|
        addFullComment("Extract Sign bit");
        //====================================================================|
        vhdl << tab << declare("sgn", 1, false) << " <= X" << of(width_ - 1) << ";" << endl;

        // Special cases
        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_-1), "zn", 1, false) << " <= '1' when (X" << range(width_ - 2, 0) << " = " << zg(width_-1) << ") else '0';" << endl;

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
        vhdl << tab << declare(getTarget()->logicDelay(1), "k", sizeRegime) << " <= "
                << zg(sizeRegime - lzocs->getCountWidth()) << " & regLength when rc /= sgn else "
                << og(sizeRegime - lzocs->getCountWidth()) << " & NOT(regLength);" << endl;

		if (wES_ > 0)
		{
			vhdl << tab << declare("exp", wES_) << " <=  shiftedPosit" << range(wF_ + wES_ - 1, wF_) << ";" << endl;
            vhdl << tab << declare("sgnVect", wES_) << " <=  (others => sgn);" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "actualExp", wES_) << " <=  (sgnVect XOR exp);" << endl;
		}
        vhdl << tab << declare("sf", wE_) << " <=  k";
		if (wES_ > 0)
		{
			vhdl << " & actualExp";
		}
		vhdl << ";" << endl;

        //=========================================================================|
		addFullComment("Extract fraction");
		// ========================================================================|
        // Discard negated regime bit and exponent bits (if any)
		vhdl << tab << declare(getTarget()->logicDelay(1), "sgnFrac", wF_+1) << " <= NOT(sgn) & shiftedPosit" << range(wF_ - 1, 0) << ";" << endl;

        //=========================================================================|
		addFullComment("Shift out fraction according to scaling factor");
		// ========================================================================|
        int shifterSize = widthO_-1;
        if (!trunc_p){
            shifterSize ++; // Add round bit
        }
        int sizeRightShift = intlog2(shifterSize);

        vhdl << tab << declare("paddedFrac", shifterSize) << " <= sgnFrac & " << zg(shifterSize - (wF_+1)) << ";" << endl;

        // shiftVal = the number of positions that frac must be shifted to the right
        if (possibleOverflow){
            vhdl << tab << declare(getTarget()->adderDelay(wE_), "shiftVal_tmp", wE_) << " <= "
                << "\"" << unsignedBinary(widthO_-2, wE_) << "\" - sf;" << endl; 

            vhdl << tab << declare("ovf", 1, false) << " <= shiftVal_tmp(shiftVal_tmp'high);" << endl;    // Overflow
            vhdl << tab << declare(getTarget()->logicDelay(3) + getTarget()->eqConstComparatorDelay(lzocs->getCountWidth()), "nudf", 1, false) << " <= "
                << "'1' when ((regLength = " << zg(lzocs->getCountWidth()) << ") OR (rc /= sgn)) else '0';" << endl;  // 0 if Underflow, useful to get final sign bit

            vhdl << tab << declare("shiftVal", sizeRightShift) << " <= shiftVal_tmp" << range(sizeRightShift-1, 0) << ";" << endl;
        }
        else{
            if (sizeRightShift > wE_)
            {
                int padBits = sizeRightShift - wE_;
                vhdl << tab << declare("sfExtended", sizeRightShift) << " <= (" << padBits-1 << " downto 0 => sf'high) & sf;" << endl;
            }

            vhdl << tab << declare(getTarget()->adderDelay(sizeRightShift), "shiftVal", sizeRightShift) << " <= "
                << "\"" << unsignedBinary(widthO_-2, sizeRightShift) << "\" - sf";
            if (sizeRightShift > wE_)
            {
                vhdl << "Extended";
            }
            vhdl << ";" << endl;
        }
        
        if (trunc_p){
            newInstance("Shifter",
                        "RightShifterComponent",
                        "wX=" + to_string(shifterSize) + " wR=" + to_string(shifterSize) + " maxShift=" + to_string(widthO_) + " dir=1 inputPadBit=true",
                        "X=>paddedFrac, S=>shiftVal, padBit=>sgn",
                        "R=>shiftedFrac");
            
            vhdl << tab << declare("intNumber", widthO_) << " <= sgn & shiftedFrac;" << endl;
        }
        else{
            newInstance("Shifter",
                        "RightShifterComponent",
                        "wX=" + to_string(shifterSize) + " wR=" + to_string(shifterSize) + " maxShift=" + to_string(shifterSize) + " dir=1 inputPadBit=true computeSticky=true",
                        "X=>paddedFrac, S=>shiftVal, padBit=>sgn",
                        "R=>shiftedFrac, Sticky=>stk");

            vhdl << tab << declare("rnd", 1, false) << " <= shiftedFrac" << of(0) << ";" << endl;
            vhdl << tab << declare("lsb", 1, false) << " <= shiftedFrac" << of(1) << ";" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(3), "round", 1, false) << " <= rnd AND (lsb OR stk);" << endl;
            
            vhdl << tab << declare(getTarget()->adderDelay(widthO_), "intNumber", widthO_) << " <= (sgn & shiftedFrac" << range(shifterSize-1, 1) << ") + round;" << endl;
        }

        if (possibleOverflow){
            vhdl << tab << declare(getTarget()->logicDelay(4), "result", widthO_) << " <= intNumber when (zn OR ovf OR NOT(nudf)) = '0' else ( ((sgn AND zn) OR (NOT(zn) AND nudf)) & " << zg(widthO_-1) << ");" << endl;
        }
        else{
            vhdl << tab << declare(getTarget()->logicDelay(1), "result", widthO_) << " <= intNumber when zn = '0' else (sgn & " << zg(widthO_-1) << ");" << endl;
        }

        vhdl << tab << "R <= result;" << endl;

        addFullComment("End of vhdl generation");
    }

	Posit2Int::~Posit2Int() {
	}

    void Posit2Int::emulate(TestCase * tc){
    }

	OperatorPtr Posit2Int::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
        int widthI, wES, widthO;
        bool trunc;
		UserInterface::parseStrictlyPositiveInt(args, "widthI", &widthI);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		UserInterface::parseStrictlyPositiveInt(args, "widthO", &widthO);
		UserInterface::parseBoolean(args, "trunc", &trunc);
		return new Posit2Int(parentOp, target, widthI, wES, widthO, trunc);
	}

	void Posit2Int::registerFactory()
	{
		UserInterface::add("Posit2Int", // name
                            "Conversion from posit to integer.",
                            "Conversions",
                            "", //seeAlso
                            "widthI(int): posit size in bits;\
                            wES(int): posit exponent size in bits;\
                            widthO(int): output integer size in bits;\
                            trunc(bool)=false: true means truncated (cheaper), false means rounded",
                            "", // htmldoc
                            Posit2Int::parseArguments);
	}

} // namespace flopoco
