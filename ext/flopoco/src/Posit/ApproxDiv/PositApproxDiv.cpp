/*
  An approximate Posit Divider for FloPoCo
  
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

#include "PositApproxDiv.hpp"
#include "Posit/Mult/PositMult.hpp"
#include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositApproxDiv::PositApproxDiv(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{

		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositApproxDiv";

		if (width_ < 3)
		{
			throw std::string("PositApproxDiv Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositApproxDiv Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------

		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		int maxExp = (1 << wES_) * (width_ - 2);
		int minExp = -maxExp;

		name << "PositApproxDiv_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", width_);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositApproxDiv \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositApproxDiv");

		//====================================================================|
		addFullComment("Compute approximate reciprocal of Y");
		//====================================================================|
		// Add '0' at begining to check if Y is 0/NaR, sso recipY should be NaRo recipY should be NaR
		vhdl << tab << declare(getTarget()->logicDelay(1) + getTarget()->adderDelay(width_), "compY", width_) << " <= ('0' & not(Y" << range(width_ - 2, 0) << ")) + 1;" << endl;
		vhdl << tab << declare(getTarget()->logicDelay(2), "recipY", width_) << " <= (Y(Y'high) OR compY(compY'high)) & compY" << range(width_-2, 0) << ";" << endl;

		//====================================================================|
		addFullComment("Perform multiplication with reciprocal of Y");
		//====================================================================|
		newInstance("PositMult",
					"ReciprocalMultiplier",
					"width=" + to_string(width_) + " wES=" + to_string(wES_),
					"X=>X, Y=>recipY",
					"R=>R");

		addFullComment("End of vhdl generation");
	}

	PositApproxDiv::~PositApproxDiv() {}

	void PositApproxDiv::emulate(TestCase *tc)
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
		mpfr_div(r, x, y, GMP_RNDN);
		
		// Set outputs
		PositNumber posr(width_, wES_, r);
		mpz_class svR = posr.getSignalValue();
		tc->addExpectedOutput("R", svR);
		
		// clean up
		mpfr_clears(x, y, r, NULL);
	}

	void PositApproxDiv::buildStandardTestCases(TestCaseList* tcl)
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


	OperatorPtr PositApproxDiv::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositApproxDiv(parentOp, target, width, wES);
	}

	void PositApproxDiv::registerFactory()
	{
		UserInterface::add("PositApproxDiv", // name
						   "An approximate posit divider.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositApproxDiv::parseArguments);
	}

} // namespace flopoco
