/*
  A Posit FMA for FloPoCo
  
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

#include "PositMAC.hpp"
#include "Posit/Decoder/PositFastDecoder.hpp"
#include "ShiftersEtc/Normalizer.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "IntMult/IntMultiplier.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositMAC::PositMAC(OperatorPtr parentOp, Target *target, int width, int wES, int carry) : Operator(parentOp, target), width_(width), wES_(wES), C_(carry)
	{

		setCopyrightString("Raul Murillo (2021)");
		addHeaderComment("Inputs: this FMA computes A*B+C");

		// Signed multiplication needs this library
		useNumericStd();

		ostringstream name;

		srcFileName = "PositMAC";

		if (width_ < 3)
		{
			throw std::string("PositMAC Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositMAC Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;
		wQuire_ = 1 + C_ + 1 + 2 * maxExp - 2 * minExp;

		name << "PositMAC_" << width_ << "_" << wES_ << "_Quire_" << wQuire_;
		setNameWithFreqAndUID(name.str());

		addInput("A", width_);
		addInput("B", width_);
		addInput("C", wQuire_);
		addOutput("R", wQuire_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositMAC \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES << " and " << carry);
		REPORT(DEBUG, "debug of PositMAC");

		//====================================================================|
		addFullComment("Decode A & B operands");
		//====================================================================|
		newInstance("PositFastDecoder",
					"A_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>A",
					"Sign=>A_sgn,SF=>A_sf,Frac=>A_f,NZN=>A_nzn");

		newInstance("PositFastDecoder",
					"B_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>B",
					"Sign=>B_sgn,SF=>B_sf,Frac=>B_f,NZN=>B_nzn");

		//====================================================================|
		addFullComment("Multiply A & B");
		//====================================================================|
		addComment("Sign and Special Cases Computation");
		vhdl << tab << declare(getTarget()->logicDelay(2), "AB_nzn", 1, false) << " <= A_nzn AND B_nzn;" << endl;
		// Really neccessary??
		vhdl << tab << declare(getTarget()->logicDelay(4), "AB_nar", 1, false) << " <= (A_sgn AND NOT(A_nzn)) OR (B_sgn AND NOT(B_nzn));" << endl;

		addComment("Multiply the fractions");
		int multSize = 2 * (wF_ + 2);
		// TODO: Add delay
		// vhdl << tab << declare("AB_f", multSize) << " <= "
		//      << "std_logic_vector(signed(A_sgn & NOT(A_sgn) & A_f) * signed(B_sgn & NOT(B_sgn) & B_f));" << endl;

		// Another approach - using FloPoCo IntMultiplier
		vhdl << tab << declare(getTarget()->logicDelay(1), "AA_f", wF_ + 2) << " <= A_sgn & NOT(A_sgn) & A_f;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(1), "BB_f", wF_ + 2) << " <= B_sgn & NOT(B_sgn) & B_f;" << endl;
		newInstance("IntMultiplier",
					"FracMultiplier",
					"wX=" + to_string(wF_ + 2) + " wY=" + to_string(wF_ + 2) + " wOut=" + to_string(multSize) + " signedIO=true",
					"X=>AA_f, Y=>BB_f",
					"R=>AB_f");
		// vhdl << tab << declare(getTarget()->adderDelay(multSize), "AB_f", multSize) << " <= std_logic_vector(signed(AA_f) * signed(BB_f));" << endl;

		// Adjust overflow
		vhdl << tab << declare("AB_sgn", 1, false) << " <= AB_f" << of(multSize - 1) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "AB_ovfExtra", 1, false) << " <= "
			 << "NOT(AB_sgn) AND AB_f" << of(multSize - 2) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "AB_ovf", 1, false) << " <= "
			 << "AB_ovfExtra OR (AB_sgn XOR AB_f" << of(multSize - 3) << ");" << endl;
		// Normalize fraction
		// Without sign bit
		vhdl << tab << declare(getTarget()->logicDelay(2), "AB_normF", multSize - 3) << " <= "
			 << "AB_f" << range(multSize - 4, 0) << " when AB_ovf = '1' else (AB_f" << range(multSize - 5, 0) << " & '0');" << endl;

		addComment("Add the exponent values");
		// TODO: Use IntAdder with carry??
		// vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "AB_sf", wE_ + 1) << " <= "
		//      << "std_logic_vector(unsigned(A_sf(A_sf'high) & A_sf) + unsigned(B_sf(B_sf'high) & B_sf)) when AB_ovfExtra='0' else "
		//      << "std_logic_vector(unsigned(A_sf(A_sf'high) & A_sf) + unsigned(B_sf(B_sf'high) & B_sf) + 1);" << endl;
		// vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "AB_sf_tmp", wE_ + 1) << " <= "
		//      << "std_logic_vector(unsigned(A_sf(A_sf'high) & A_sf) + unsigned(B_sf(B_sf'high) & B_sf) + (AB_ovfExtra &\"\"));" << endl;

		// vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "AB_sf", wE_ + 1) << " <= "
		//      << "std_logic_vector(unsigned(AB_sf_tmp) + (AB_ovf & \"\"));" << endl;

		// Another approach - using FloPoCo IntAdder
		vhdl << tab << declare("AA_sf", wE_ + 1) << " <= A_sf(A_sf'high) & A_sf;" << endl;
		vhdl << tab << declare("BB_sf", wE_ + 1) << " <= B_sf(B_sf'high) & B_sf;" << endl;
		newInstance("IntAdder",
					"SFAdder",
					"wIn=" + to_string(wE_ + 1),
					"X=>AA_sf,Y=>BB_sf,Cin=>AB_ovfExtra",
					"R=>AB_sf_tmp");

		newInstance("IntAdder",
					"RoundingAdder",
					"wIn=" + to_string(wE_ + 1),
					"X=>AB_sf_tmp,Cin=>AB_ovf",
					"R=>AB_sf",
					"Y=>" + zg(wE_ + 1));

		//====================================================================|
		addFullComment("Shift AB fraction into corresponding quire format");
		//====================================================================|
		// int bias = 2 * maxExp;
		// int offsetSize = intlog2(2 * bias + 1);
		// vhdl << tab << declare(getTarget()->adderDelay(offsetSize), "AB_sfBiased", offsetSize) << " <= "
		//      << "std_logic_vector(\"" << unsignedBinary(bias, offsetSize) << "\" - unsigned(AB_sf)) when AB_ovf='0' else "
		//      << "std_logic_vector(\"" << unsignedBinary(bias, offsetSize) << "\" - unsigned(AB_sf) - 1);" << endl;

		// // Add sign bit to fraction & pad with 0's
		// vhdl << tab << declare("paddedFrac", wQuire_ - (C_ + 1)) << " <= NOT(AB_sgn) & AB_normF & " << zg(wQuire_ - (C_ + 1) - (multSize - 2)) << ";" << endl;

		// // Approach: shift only the fixed point part - generate carry bits appart
		// // TODO: Another approach; use 2 smaller Right shifters, one for pos exp, another for neg exp.
		// newInstance("Shifter",
		//             "QuireRightShifter",
		//             "wX=" + to_string(wQuire_ - (C_ + 1)) + " wR=" + to_string(wQuire_ - (C_ + 1)) + " maxShift=" + to_string(2 * bias + 1) + " dir=1 inputPadBit=true",
		//             "X=>paddedFrac, S=>AB_sfBiased, padBit=>AB_sgn",
		//             "R=>fixedPosit");

		// vhdl << tab << declare("zeros", wQuire_ - 1) << " <= " << zg(wQuire_ - 1) << ";" << endl;
		// vhdl << tab << declare(getTarget()->logicDelay(), "AB_quire", wQuire_) << " <= "
		//      << "(" << wQuire_ - 1 << " downto " << wQuire_ - (C_ + 1) << " => AB_sgn) & fixedPosit when AB_nzn='1' else "
		//      << "(AB_nar & zeros);" << endl;

		/////////////////////////////

		// Another approach
		// Use a single small Right shifter, adjust sfBiased according to sf sign.
		// Sanity check // CAUTION: this does not work when width = 9, 17, 33, etc.
		if (wE_ != intlog2(2 * maxExp))
			throw std::string("Sorry, but parameters for posit MAC do not match sizes.");
		vhdl << tab << declare("neg_sf", 1, false) << " <= AB_sf" << of(wE_) << ";" << endl;
		vhdl << tab << declare("AB_effectiveSF", wE_) << " <= AB_sf" << range(wE_ - 1, 0) << ";" << endl;

		// vhdl << tab << declare(getTarget()->adderDelay(intlog2(2*maxExp)) + getTarget()->logicDelay(), "AB_sfBiased", intlog2(2*maxExp)) << " <= "
		//      << "std_logic_vector(\"" << unsignedBinary(2*maxExp, intlog2(2*maxExp)) << "\" - unsigned(AB_effectiveSF)) when neg_sf='0' else "
		//      << "std_logic_vector(unsigned(NOT(AB_effectiveSF)));" << endl;

		// Another aproach - shift only in the integer/frac part of the quire
		vhdl << tab << declare(getTarget()->logicDelay(1), "adderInput", intlog2(2 * maxExp)) << " <= ";
		if (wE_ < intlog2(2 * maxExp))
		{
			vhdl << "(" << intlog2(2 * maxExp) - 1 << " downto " << wE_ << " => NOT(neg_sf)) & ";
		}
		vhdl << "NOT(AB_effectiveSF);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(1), "adderBias", intlog2(2 * maxExp)) << " <= \"" << unsignedBinary(2 * maxExp, intlog2(2 * maxExp)) << "\" when neg_sf='0' else "
			 << og(wE_) << ";" << endl;
		vhdl << tab << declare("ob", 1, false) << " <= '1';" << endl;
		newInstance("IntAdder",
					"BiasedSFAdder",
					"wIn=" + to_string(intlog2(2 * maxExp)),
					"X=>adderInput,Y=>adderBias,Cin=>ob",
					"R=>AB_sfBiased");

		vhdl << tab << declare(getTarget()->logicDelay(1), "paddedFrac", 2 * (maxExp + 1 + wF_)) << " <= "
			 << "(NOT(AB_sgn) & AB_normF & " << zg(2 * maxExp) << ") when neg_sf='0' else "
			 << "((" << 2 * (maxExp + 1 + wF_) - 1 << " downto " << 2 * maxExp << " => AB_sgn) & NOT(AB_sgn) & AB_normF & " << zg(2 * maxExp - 2 * (1 + wF_)) << ");" << endl;

		newInstance("Shifter",
					"Frac_RightShifter",
					"wX=" + to_string(2 * (maxExp + 1 + wF_)) + " wR=" + to_string(2 * (maxExp + 1 + wF_)) + " maxShift=" + to_string(2 * maxExp) + " dir=1 inputPadBit=true",
					"X=>paddedFrac, S=>AB_sfBiased, padBit=>AB_sgn",
					"R=>fixedPosit");

		vhdl << tab << declare(getTarget()->logicDelay(1), "quirePosit", 4 * maxExp + 1) << " <= "
			 << "(fixedPosit & " << zg(2 * (maxExp - 1 - wF_) + 1) << ") when neg_sf='0' else "
			 << "((" << 4 * maxExp << " downto " << 2 * (maxExp + 1 + wF_) << " => AB_sgn) & fixedPosit);" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(1), "AB_quire", wQuire_) << " <= "
			 << "((" << wQuire_ - 1 << " downto " << 4 * maxExp + 1 << " => AB_sgn) & " // generate carry bits & pad with 0's
			 << "quirePosit) when AB_nzn='1' else AB_nar & (" << wQuire_ - 2 << " downto 0 => '0');" << endl;

		//====================================================================|
		addFullComment("Add quires");
		//====================================================================|
		// vhdl << tab << declare(getTarget()->adderDelay(wQuire_), "ABC_add", wQuire_) << " <= "
		//      << "std_logic_vector(unsigned(AB_quire) + unsigned(C));" << endl;

		vhdl << tab << declare("zb", 1, false) << " <= '0';" << endl;
		newInstance("IntAdder",
					"QuireAdder",
					"wIn=" + to_string(wQuire_),
					"X=>AB_quire,Y=>C,Cin=>zb",
					"R=>ABC_add");
		// vhdl << tab << declare(getTarget()->adderDelay(wQuire_), "ABC_add", wQuire_) << " <= std_logic_vector(unsigned(AB_quire) + unsigned(C));" << endl;

		// Special case
		vhdl << tab << declare("zeros", wQuire_ - 1) << " <= (others => '0');" << endl;
		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wQuire_ - 1), "C_nar", 1, false) << " <= "
			 << "C" << of(wQuire_ - 1) << " when (C" << range(wQuire_ - 2, 0) << " = zeros) else '0';" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "ABC_nar", 1, false) << " <= "
			 << "AB_nar OR C_nar;" << endl;

		vhdl << tab << declare(getTarget()->logicDelay(1), "result", wQuire_) << " <= "
			 << "ABC_add when ABC_nar='0' else ('1' & zeros);" << endl;

		vhdl << tab << "R <= result;" << endl;

		addFullComment("End of vhdl generation");
	}

	PositMAC::~PositMAC() {}

	void PositMAC::emulate(TestCase *tc) {}

	OperatorPtr PositMAC::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES, carry;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		UserInterface::parsePositiveInt(args, "carry", &carry);
		return new PositMAC(parentOp, target, width, wES, carry);
	}

	void PositMAC::registerFactory()
	{
		UserInterface::add("PositMAC", // name
						   "A posit FMA unit without rounding.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits;\
                            carry(int): extra bits to savely allow the sum up to 2^carry products",
						   "", // htmldoc
						   PositMAC::parseArguments);
	}

} // namespace flopoco
