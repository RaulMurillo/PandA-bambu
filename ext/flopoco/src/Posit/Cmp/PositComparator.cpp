/*
  A Posit comparator unit
  
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

#include "PositComparator.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

	PositComparator::PositComparator(OperatorPtr parentOp, Target *target, int width, int wES, int criteria) : Operator(parentOp, target), width_(width), wES_(wES), criteria_(criteria)
	{
		setCopyrightString("Raul Murillo (2021)");
		ostringstream name;
		srcFileName = "PositComparator";

		if (width_ < 3)
		{
			throw std::string("PositComparator Constructor: width is too small, should be greater than two");
		}
		if (wES_ >= width_ - 3)
		{
			//Avoid posits without even one bit of precision
			throw std::string("PositComparator Constructor: invalid value of wES");
		}

		// -------- Parameter set up -----------------
		int sizeRegime = intlog2(width_ - 1) + 1;
		wE_ = sizeRegime + wES_;
		wF_ = width_ - 3 - wES_;
		// int maxExp = (1 << wES_) * (width_ - 2);
		// int minExp = -maxExp;

		switch (criteria)
		{
		case -2:
			name << "PositComparator_" << width_ << "_" << wES_ << "_"
				 << "less";
			break;
		case -1:
			name << "PositComparator_" << width_ << "_" << wES_ << "_"
				 << "leq";
			break;
		case 0:
			// TODO: Optimize, addapt for constant comparison
			name << "PositComparator_" << width_ << "_" << wES_ << "_"
				 << "eq";
			break;
		case 1:
			name << "PositComparator_" << width_ << "_" << wES_ << "_"
				 << "geq";
			break;
		case 2:
			name << "PositComparator_" << width_ << "_" << wES_ << "_"
				 << "greater";
			break;
		default:
			name << "PositComparator_" << width_ << "_" << wES_ << "_"
				 << "eq";
		}
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		addOutput("R", 1);

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of PositComparator \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << width << ", " << wES);
		REPORT(DEBUG, "debug of PositComparator");

		switch (criteria)
		{
		case -2:
			vhdl << tab << declare(getTarget()->eqComparatorDelay(width), "cmp_r", 1, false) << " <= '1' when (X < Y) else '0';" << endl;
			break;
		case -1:
			vhdl << tab << declare(getTarget()->eqComparatorDelay(width), "cmp_r", 1, false) << " <= '1' when (X <= Y) else '0';" << endl;
			break;
		case 0:
			vhdl << tab << declare(getTarget()->eqComparatorDelay(width), "cmp_r", 1, false) << " <= '1' when (X = Y) else '0';" << endl;
			break;
		case 1:
			vhdl << tab << declare(getTarget()->eqComparatorDelay(width), "cmp_r", 1, false) << " <= '1' when (X >= Y) else '0';" << endl;
			break;
		case 2:
			vhdl << tab << declare(getTarget()->eqComparatorDelay(width), "cmp_r", 1, false) << " <= '1' when (X > Y) else '0';" << endl;
			break;
		default:
			throw std::string("PositComparator Constructor: invalid criteria. Expecting an int between -2 and 2");
		}
		vhdl << tab << "R(0) <= cmp_r;" << endl;
	};

	PositComparator::~PositComparator() {}

	void PositComparator::emulate(TestCase *tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");

		mpz_class svR;

		switch (criteria_)
		{
		case -2:
			if (svX < svY)
				svR = 1;
			else
				svR = 0;
			break;
		case -1:
			if (svX <= svY)
				svR = 1;
			else
				svR = 0;
			break;
		case 0:
			if (svX == svY)
				svR = 1;
			else
				svR = 0;
			break;
		case 1:
			if (svX >= svY)
				svR = 1;
			else
				svR = 0;
			break;
		case 2:
			if (svX > svY)
				svR = 1;
			else
				svR = 0;
			break;
		default:;
		}

		tc->addExpectedOutput("R", svR);
	}

	OperatorPtr PositComparator::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int width, wES;
		UserInterface::parseStrictlyPositiveInt(args, "width", &width);
		UserInterface::parsePositiveInt(args, "wES", &wES);
		int criteria;
		UserInterface::parseInt(args, "criteria", &criteria);
		return new PositComparator(parentOp, target, width, wES, criteria);
	}

	void PositComparator::registerFactory()
	{
		UserInterface::add("PositComparator", // name
						   "An approximate posit divider.",
						   "Posit",
						   "", //seeAlso
						   "width(int): posit size in bits;\
                            wES(int): posit exponent size in bits;\
							criteria(int): -2=lesser than, -1=lesser or equal, 0=equal, 1=greater or equal, 2=greater",
						   "", // htmldoc
						   PositComparator::parseArguments);
	}

} // namespace flopoco
