/*
  A Posit Multiplier for FloPoCo
  
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

#include "PositMult.hpp"
#include "Posit/Decoder/PositFastDecoder.hpp"
#include "Posit/Encoder/PositFastEncoder.hpp"
// #include "IntAddSubCmp/IntAdder.hpp"
#include "IntMult/IntMultiplier.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositMult::PositMult(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositMult";
		// Signed multiplication needs this library
		useNumericStd();

		if (width_ < 3)
		{
			throw std::string("PositMult Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositMult Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositMult_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositMult \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositMult");

		//====================================================================|
		addFullComment("Decode X & Y operands");
		//====================================================================|
		newInstance("PositFastDecoder",
					"X_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X",
					"Sign=>X_sgn,SF=>X_sf,Frac=>X_f,NZN=>X_nzn");

		newInstance("PositFastDecoder",
					"Y_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>Y",
					"Sign=>Y_sgn,SF=>Y_sf,Frac=>Y_f,NZN=>Y_nzn");

		//====================================================================|
		addFullComment("Multiply X & Y");
		//====================================================================|
		addComment("Sign and Special Cases Computation");
		// vhdl << declare(getTarget()->logicDelay(2), "XY_sgn", 1, false) << " <= X_sgn XOR Y_sgn;" << endl;
		// vhdl << declare(getTarget()->logicDelay(2), "XY_z", 1, false) << " <= X_z OR Y_z;" << endl;
		// vhdl << declare(getTarget()->logicDelay(2), "XY_nar", 1, false) << " <= X_nar OR Y_nar;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_nzn", 1, false) << " <= X_nzn AND Y_nzn;" << endl;

		int multSize = 2 * (wF_ + 2);
		// addComment("Multiply the fractions (using naive * operator, may improve in the future)");
		// TODO: Use IntAdder with carry??
		// TODO: Add delay

		// vhdl << declare("XY_f", multSize) << " <= "
		//      << "std_logic_vector(signed(X_sgn & NOT(X_sgn) & X_f) * signed(Y_sgn & NOT(Y_sgn) & Y_f));" << endl;

		// Another approach - using FloPoCo IntMultiplier
		addComment("Multiply the fractions (using FloPoCo IntMultiplier)");
		// getTarget()->setPlainVHDL(true);
		// getTarget()->setTilingMethod();
		vhdl << tab << declare(getTarget()->logicDelay(1), "XX_f", wF_ + 2) << " <= X_sgn & NOT(X_sgn) & X_f;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(1), "YY_f", wF_ + 2) << " <= Y_sgn & NOT(Y_sgn) & Y_f;" << endl;
		newInstance("IntMultiplier",
					"FracMultiplier",
					"wX=" + to_string(wF_ + 2) + " wY=" + to_string(wF_ + 2) + " wOut=" + to_string(multSize) + " signedIO=true",
					"X=>XX_f, Y=>YY_f",
					"R=>XY_f");

		// Adjust overflow
		vhdl << tab << declare("XY_sgn", 1, false) << " <= XY_f" << of(multSize - 1) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_ovfExtra", 1, false) << " <= "
			 << "NOT(XY_sgn) AND XY_f" << of(multSize - 2) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_ovf", 1, false) << " <= "
			 << " (XY_sgn XOR XY_f" << of(multSize - 3) << ");" << endl;

		// Normalize fraction
		// Without including hidden bits
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_normF", multSize - 3) << " <= "
			 << "XY_f" << range(multSize - 4, 0) << " when (XY_ovfExtra OR XY_ovf) = '1' else "
			 << "(XY_f" << range(multSize - 5, 0) << " & '0');" << endl;

		addComment("Add the exponent values");
		// // TODO: Use IntAdder with carry??
		// vhdl << declare(getTarget()->adderDelay(wE_+1), "XY_sf", wE_+1) << " <= "
		//      << "(X_sf(X_sf'high) & X_sf) + (Y_sf(Y_sf'high) & Y_sf) + XY_ovf + XY_ovfExtra;" << endl;

		// vhdl << declare(getTarget()->logicDelay(1), "XY_ovfBits", 2) << " <= "
		//  << "\"10\" when XY_ovfExtra = '1' else ('0' & XY_ovf);" << endl;
		vhdl << tab << declare("XY_ovfBits", 2) << " <= XY_ovfExtra & XY_ovf;" << endl;
		vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "XY_sf", wE_ + 1) << " <= "
			 << "std_logic_vector(unsigned(X_sf(X_sf'high) & X_sf) + unsigned(Y_sf(Y_sf'high) & Y_sf) + unsigned(XY_ovfBits));" << endl;

		//====================================================================|
		addFullComment("Generate final posit");
		//====================================================================|
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_finalSgn", 1, false) << " <= "
			 << "XY_sgn when XY_nzn = '1' else (X_sgn AND NOT(X_nzn)) OR (Y_sgn AND NOT(Y_nzn));" << endl;

		vhdl << tab << declare("XY_frac", wF_) << " <= XY_normF" << range(multSize - 4, multSize - 4 - wF_ + 1) << ";" << endl;
		vhdl << tab << declare("grd", 1, false) << " <= XY_normF" << of(multSize - 4 - wF_) << ";" << endl;
		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(multSize - 4 - wF_), "stk", 1, false) << " <= "
			 << "'0' when (XY_normF" << range(multSize - 4 - wF_ - 1, 0) << " = " << zg(multSize - 4 - wF_) << ") else '1';" << endl;

		newInstance("PositFastEncoder",
					"PositEncoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"Sign=>XY_finalSgn, SF=>XY_sf, Frac=>XY_frac, Guard=>grd, Sticky=>stk, NZN=>XY_nzn",
					"R=>R");

		addFullComment("End of vhdl generation");
	}

	PositMult::~PositMult() {}

	void PositMult::emulate(TestCase *tc) {}

	OperatorPtr PositMult::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositMult(parentOp, target, width, wES);
	}

	void PositMult::registerFactory()
	{
		UserInterface::add("PositMult", // name
						   "A correctly rounded posit multiplier.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositMult::parseArguments);
	}

} // namespace flopoco
