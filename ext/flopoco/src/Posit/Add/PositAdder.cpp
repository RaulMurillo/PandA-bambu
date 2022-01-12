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

#include "PositAdder.hpp"
#include "Posit/Decoder/PositFastDecoder.hpp"
#include "Posit/Encoder/PositFastEncoder.hpp"
#include "ShiftersEtc/Shifters.hpp"
#include "ShiftersEtc/Normalizer.hpp"
// #include "IntAddSubCmp/IntAdder.hpp"
#include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositAdder::PositAdder(OperatorPtr parentOp, Target *target, int width, int wES, bool sub) : Operator(parentOp, target), width_(width), wES_(wES), sub_(sub)
	{
		setCopyrightString("Raul Murillo (2021)");
		srcFileName = "PositAdder";
		ostringstream name;
		if(sub_)
		{
			name << "PositSubtractor_";
		}
		else
		{
			name << "PositAdder_";
		}

		if (width_ < 3)
		{
			throw std::string("PositAdder Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositAdder Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------
		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositAdder \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositAdder");

		//=========================================================================|
		addFullComment("Decode X & Y operands");
		// ========================================================================|
		string pmY = "Y";
		if (sub_) {	// minus Y
			vhdl << tab << declare(getTarget()->logicDelay(1) + getTarget()->adderDelay(width_), "mY", width_)
					<< " <= NOT(Y) + 1;"<<endl;
			pmY = "mY";
		}
		newInstance("PositFastDecoder",
					"X_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X",
					"Sign=>X_sgn,SF=>X_sf,Frac=>X_f,NZN=>X_nzn");
		
		newInstance("PositFastDecoder",
					"Y_decoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>"+pmY,
					"Sign=>Y_sgn,SF=>Y_sf,Frac=>Y_f,NZN=>Y_nzn");

		addFullComment("Check for Zeros and NaRs");
		vhdl << tab << declare(getTarget()->logicDelay(2), "X_not_zero", 1, false) << " <= X_sgn OR X_nzn;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "X_nar", 1, false) << " <= X_sgn AND NOT(X_nzn);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "Y_not_zero", 1, false) << " <= Y_sgn OR Y_nzn;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "Y_nar", 1, false) << " <= Y_sgn AND NOT(Y_nzn);" << endl;

		//=========================================================================|
		addFullComment("Compare operands and adjust values");
		// ========================================================================|
		vhdl << tab << declare(getTarget()->eqComparatorDelay(wE_), "is_larger", 1, false) << " <= '1' when (signed(X_sf) > signed(Y_sf)) else '0';" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(1), "larger_sf", wE_) << " <= " << endl
			 << tab << tab << "X_sf when '1'," << endl
			 << tab << tab << "Y_sf when '0'," << endl
			 << tab << tab << "\"" << string(wE_, '-') << "\" when others;" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(1), "smaller_sf", wE_) << " <= " << endl
			 << tab << tab << "Y_sf when '1'," << endl
			 << tab << tab << "X_sf when '0'," << endl
			 << tab << tab << "\"" << string(wE_, '-') << "\" when others;" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(2), "larger_frac", 2 + wF_) << " <= " << endl
			 << tab << tab << "(X_sgn & (NOT(X_sgn) AND X_not_zero) & X_f) when '1'," << endl
			 << tab << tab << "(Y_sgn & (NOT(Y_sgn) AND Y_not_zero) & Y_f) when '0'," << endl
			 << tab << tab << "\"" << string(2 + wF_, '-') << "\" when others;" << endl;

		vhdl << tab << "with is_larger select " << declare(getTarget()->logicDelay(2), "smaller_frac", 2 + wF_) << " <= " << endl
			 << tab << tab << "(Y_sgn & (NOT(Y_sgn) AND Y_not_zero) & Y_f) when '1'," << endl
			 << tab << tab << "(X_sgn & (NOT(X_sgn) AND X_not_zero) & X_f) when '0'," << endl
			 << tab << tab << "\"" << string(2 + wF_, '-') << "\" when others;" << endl;

		//=========================================================================|
		addFullComment("Compute exponents difference & align fractions");
		// ========================================================================|
		/*
		// as larger_sf >= smaller_sf, the difference is always positive
		// and can be represented as an unsigned int with wE_ bits
		*/
		vhdl << tab << declare(getTarget()->adderDelay(wE_), "offset", wE_) << " <= larger_sf - smaller_sf;" << endl;
		
		int maxshiftsize = intlog2(wF_ + 4);
		if ((wE_) > maxshiftsize)
		{
			addComment("Saturate maximum offset");
			vhdl << tab << declare(getTarget()->eqConstComparatorDelay(wE_ - maxshiftsize + 1), "shift_saturate", 1, false) << " <= "
				 << "'0' when (offset" << range(wE_ - 1, maxshiftsize) << " = " << zg(wE_ - maxshiftsize + 1)<< ") else '1';" << endl;

			vhdl << tab << declare(getTarget()->logicDelay(1), "frac_offset", maxshiftsize) << " <= "
				 << "CONV_STD_LOGIC_VECTOR(" << wF_ + 4 << "," << maxshiftsize << ") when shift_saturate = '1' else "
				 << "offset" << range(maxshiftsize - 1, 0) << ";" << endl;
		}
		else
		{ // Only for wES==0. In this case (wE_) == maxshiftsize
			if ((wE_) != maxshiftsize)
			{ // Sanity check - Unreachable point
				throw std::string("Sorry, but parameters for posit encoder do not match sizes.");
			}
			vhdl << tab << declare("frac_offset", maxshiftsize) << " <= offset;" << endl;
		}

		addComment("Align fractions - right shift the smaller one");
		// append two extra bits, so rounding info is not missed during shifting
		vhdl << tab << declare("input_shifter", 2 + wF_ + 2) << " <= smaller_frac & \"00\";" << endl;
		vhdl << tab << declare("pad_bit", 1, false) << " <= smaller_frac(smaller_frac'high);" << endl;

		ostringstream param1, inmap1, outmap1;
		param1 << "wX=" << wF_ + 4;
		param1 << " maxShift=" << wF_ + 4;
		param1 << " wR=" << wF_ + 4;
		param1 << " dir=" << Shifter::Right;
		param1 << " computeSticky=true";
		param1 << " inputPadBit=true";

		inmap1 << "X=>input_shifter,S=>frac_offset,padBit=>pad_bit";

		outmap1 << "R=>shifted_frac,Sticky=>stk_tmp";

		newInstance("Shifter",
					"RightShifterFraction",
					param1.str(),
					inmap1.str(),
					outmap1.str());

		vhdl << tab << declare("smaller_frac_sh", wF_ + 2) << " <= shifted_frac" << range(wF_ + 3, 2) << ";" << endl;
		vhdl << tab << declare("grd_tmp", 1, false) << " <= shifted_frac(1);" << endl;
		vhdl << tab << declare("rnd_tmp", 1, false) << " <= shifted_frac(0);" << endl;

		//=========================================================================|
		addFullComment("Add fractions");
		//=========================================================================|
		vhdl << tab << declare(getTarget()->adderDelay(wF_ + 3), "add_frac", wF_ + 3) << " <= "
				<< "(larger_frac" << of(wF_ + 1) << " & larger_frac) + (smaller_frac_sh" << of(wF_ + 1) << " & smaller_frac_sh);" << endl;
		vhdl << tab << declare("grd_bit", 1, false) << " <= grd_tmp;" << endl;
		vhdl << tab << declare("rnd_bit", 1, false) << " <= rnd_tmp;" << endl;
		vhdl << tab << declare("stk_bit", 1, false) << " <= stk_tmp;" << endl;

		addComment("Normalization of fraction");
		vhdl << tab << declare("count_type", 1, false) << " <= add_frac" << of(wF_ + 2) << ";" << endl;
		vhdl << tab << declare("add_frac_shift", wF_ + 5) << " <= add_frac" << range(wF_ + 1, 0) << " & grd_bit & rnd_bit & stk_bit;" << endl;
		ostringstream param_norm, inmap_norm, outmap_norm;
		param_norm << "wX=" << wF_ + 5;
		param_norm << " wR=" << wF_ + 5;//wF_ + 5;
		param_norm << " maxShift=" << wF_ + 5 - 1;	// Not sure if count -1 is ok; exhaustive test (> 8-bits) is needed (for sure  maxShift=wF_ + 5 is correct)
		param_norm << " countType=" << -1;

		inmap_norm << "X=>add_frac_shift,OZb=>count_type";
		outmap_norm << "Count=>count,R=>norm_frac_tmp";
		Normalizer *lzocs = (Normalizer *)
			newInstance("Normalizer",
						"FractionNormalizer",
						param_norm.str(),
						inmap_norm.str(),
						outmap_norm.str());
		int wCount = lzocs->getCountWidth();

		addComment("Correct final exponent");
		vhdl << tab << declare(getTarget()->adderDelay(wE_ + 1), "add_sf", wE_ + 1) << " <= (larger_sf" << of(wE_-1) << " & larger_sf) - (" << zg(wE_ + 1 - wCount) << " & count) + 1;" << endl;

		//=========================================================================|
		addFullComment("Data Rounding & Encoding");
		// ========================================================================|
		vhdl << tab << declare(getTarget()->eqConstComparatorDelay(intlog2(wF_ + 5)), "is_not_zero", 1, false) << " <= count_type when (count = " << og(wCount) << ") else '1';" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "is_nar", 1, false) << " <= X_nar OR Y_nar;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "XY_nzn", 1, false) << " <= is_not_zero AND NOT(is_nar);" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(3), "sign", 1, false) << " <= is_nar OR (is_not_zero AND add_frac" << of(wF_ + 2) << ");" << endl;

		vhdl << tab << declare("norm_frac", wF_) << " <= norm_frac_tmp" << range(wF_+3,4) << ";" << endl;
		vhdl << tab << declare("grd", 1, false) << " <= norm_frac_tmp" << of(3) << ";" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(3), "stk", 1, false) << " <= norm_frac_tmp(2) OR norm_frac_tmp(1) OR norm_frac_tmp(0);" << endl;

		newInstance("PositFastEncoder",
					"PositEncoder",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"Sign=>sign, SF=>add_sf, Frac=>norm_frac, Guard=>grd, Sticky=>stk, NZN=>XY_nzn",
					"R=>R");

		addFullComment("End of vhdl generation");
	};

	PositAdder::~PositAdder() {}

	void PositAdder::emulate(TestCase *tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		
		/* Compute correct value */
		PositNumber posx(width_, wES_, svX);
		PositNumber posy(width_, wES_, svY);
		mpfr_t x, y, r;
		mpfr_init2(x, 1000*width_ -2);
		mpfr_init2(y, 1000*width_ -2);
		mpfr_init2(r, 1000*width_ -2);
		posx.getMPFR(x);
		posy.getMPFR(y);
		if(sub_)
		{
			mpfr_sub(r, x, y, GMP_RNDN);
		}
		else
		{
			mpfr_add(r, x, y, GMP_RNDN);
		}
		
		// Set outputs
		PositNumber posr(width_, wES_, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		// clean up
		mpfr_clears(x, y, r, NULL);
	}

	void PositAdder::buildStandardTestCases(TestCaseList* tcl)
	{
		TestCase *tc;
		mpz_class x, y;

		// 1*1
		x = mpz_class(1) << (width_ - 2);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// 1*0
		x = mpz_class(1) << (width_ - 2);
		y = mpz_class(0);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1*-1
		x = mpz_class(3) << (width_ - 2);
		y = mpz_class(3) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1*1
		x = mpz_class(3) << (width_ - 2);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// nan*1
		x = mpz_class(1) << (width_ - 1);
		y = mpz_class(1) << (width_ - 2);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// nan*0
		x = mpz_class(1) << (width_ - 1);
		y = mpz_class(0);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue*maxvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = (mpz_class(1) << (width_ - 1))-1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue*-maxvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = (mpz_class(1) << (width_ - 1))+1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -minvalue*maxvalue
		x = (mpz_class(1) << (width_))-1;
		y = (mpz_class(1) << (width_ - 1))-1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// maxvalue*minvalue
		x = (mpz_class(1) << (width_ - 1))-1;
		y = mpz_class(1);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}

	OperatorPtr PositAdder::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		bool sub;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		UserInterface::parseBoolean(args, "sub", &sub);
		return new PositAdder(parentOp, target, width, wES, sub);
	}

	void PositAdder::registerFactory()
	{
		UserInterface::add("PositAdder", // name
						   "A posit adder/subtractor with a single architecture.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits;\
                            sub(bool)=false: true means subtractor, false means adder",
						   "", // htmldoc
						   PositAdder::parseArguments);
	}

} // namespace flopoco
