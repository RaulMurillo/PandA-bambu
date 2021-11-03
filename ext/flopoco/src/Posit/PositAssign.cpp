/*
  A Posit naive assign operator
  
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

#include "PositAssign.hpp"
// #include "TestBenches/PositNumber.hpp"
// #include "TestBenches/IEEENumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositAssign::PositAssign(OperatorPtr parentOp, Target *target, int width, int wES) : Operator(parentOp, target), width_(width), wES_(wES)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositAssign";

		if (width_ < 3)
		{
			throw std::string("PositAssign Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositAssign Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------
		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		// int maxExp = (1 << wES_) * (width_ - 2);
		// int minExp = -maxExp;

		name << "PositAssign_" << width_ << "_" << wES_;
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addOutput("R", width_);

		REPORT(INFO, "Declaration of PositAssign \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositAssign");

		/*	VHDL code description	*/
		vhdl << tab << "R <= X;" << endl;
	};

	PositAssign::~PositAssign() {}

	void PositAssign::emulate(TestCase *tc) {}

	OperatorPtr PositAssign::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		return new PositAssign(parentOp, target, width, wES);
	}

	void PositAssign::registerFactory()
	{
		UserInterface::add("PositAssign", // name
						   "A naive identity operator for posit numbers.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits",
						   "", // htmldoc
						   PositAssign::parseArguments);
	}

} // namespace flopoco
