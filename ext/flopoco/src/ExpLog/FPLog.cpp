/*
  An FP logarithm for FloPoCo

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

*/

// TODO List:
//  * test cases for boundary cases pfinal etc
//  * finetune pipeline
//  * Port back the Arith paper
#include <fstream>
#include <sstream>
#include <math.h>	// for NaN
#include "FPLog.hpp"
#include "FPLogIterative.hpp"
#include "TestBenches/FPNumber.hpp"
#include "utils.hpp"
//#include "IntMult/IntSquarer.hpp"
//#include "ConstMult/IntIntKCM.hpp"
#include "UserInterface.hpp"


using namespace std;


namespace flopoco{
	
	FPLog::FPLog(OperatorPtr parentOp, Target* target, int wE, int wF):
		Operator(parentOp, target), wE(wE), wF(wF)
	{}

	FPLog::~FPLog() {}

	
	void FPLog::emulate(TestCase * tc)
	{
		/* Get I/O values */
		mpz_class svX = tc->getInputValue("X");

		/* Compute correct value */
		FPNumber fpx(wE, wF);
		fpx = svX;

		mpfr_t x, ru,rd;
		mpfr_init2(x,  1+wF);
		mpfr_init2(ru, 1+wF);
		mpfr_init2(rd, 1+wF);
		fpx.getMPFR(x);
		mpfr_log(rd, x, GMP_RNDD);
		mpfr_log(ru, x, GMP_RNDU);
#if 0
		mpfr_out_str (stderr, 10, 30, x, GMP_RNDN); cerr << " ";
		mpfr_out_str (stderr, 10, 30, rd, GMP_RNDN); cerr << " ";
		mpfr_out_str (stderr, 10, 30, ru, GMP_RNDN); cerr << " ";
		cerr << endl;
#endif
		FPNumber  fprd(wE, wF, rd);
		FPNumber  fpru(wE, wF, ru);
		mpz_class svRD = fprd.getSignalValue();
		mpz_class svRU = fpru.getSignalValue();
		tc->addExpectedOutput("R", svRD);
		tc->addExpectedOutput("R", svRU);
		mpfr_clears(x, ru, rd, NULL);
	}


	// TEST FUNCTIONS


	void FPLog::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;
		mpz_class x;


		tc = new TestCase(this);
		tc->addFPInput("X", 1.0);
		tc->addComment("1.0");
		emulate(tc);
		tcl->add(tc);

#if 0
		tc = new TestCase(this);
		tc->addComment("The worst case of the error analysis: max cancellation, and full range reduction");
		x = (mpz_class(1) << wF) - (mpz_class(1) << (wF-pfinal+2)) // mantissa
			+ (((mpz_class(1) << (wE-1)) -2) << wF)  // exponent
			+ (mpz_class(1) << (wE+wF+1))	; // exn=010
		tc->addInput("X", x);
		emulate(tc);
		tcl->add(tc);
#endif


	}



	// One test out of 8 fully random (tests NaNs etc)
	// All the remaining ones test positive numbers.
	// with special treatment for exponents 0 and -1,
	// and for the range reduction worst case.

	TestCase* FPLog::buildRandomTestCase(int i){

		TestCase *tc;
		mpz_class a;
		mpz_class normalExn = mpz_class(1)<<(wE+wF+1);

		tc = new TestCase(this);
		/* Fill inputs */
		if ((i & 7) == 0)
			a = getLargeRandom(wE+wF+3);
		else if ((i & 7) == 1) // exponent of 1
			a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-1)) << wF) + normalExn;
		else if ((i & 7) == 2) // exponent of 0.5
			a  = getLargeRandom(wF) + ((((mpz_class(1)<<(wE-1))-2)) << wF) + normalExn;
#if 0 // To resurrect somehow some day
		else if ((i & 7) == 3) { // worst case for range reduction
			tc->addComment("The worst case of the error analysis: max cancellation, and full range reduction");
			a = (mpz_class(1) << wF) - (mpz_class(1) << (wF-pfinal+2)) + getLargeRandom(wF-pfinal+2) // mantissa
				+ (((mpz_class(1) << (wE-1)) -2) << wF)  // exponent
				+ (mpz_class(1) << (wE+wF+1))	; // exn=010
		}
#endif
		else
			a  = getLargeRandom(wE+wF)  + normalExn; // 010xxxxxx

		tc->addInput("X", a);
		/* Get correct outputs */
		emulate(tc);
		// add to the test case list
		return tc;
	}




		TestList FPLog::unitTest(int index)
		{
		// the static list of mandatory tests
			TestList testStateList;
			vector<pair<string,string>> paramList;
			
			if(index==-1) 
		{ // The unit tests

			// First test with plainVHDL, then with cool multipliers
			for(int wF=5; wF<53; wF+=1) // test various input widths
			{ 
				int nbByteWE = 6+(wF/10);
				while(nbByteWE>wF){
					nbByteWE -= 2;
				}

				paramList.push_back(make_pair("wF",to_string(wF)));
				paramList.push_back(make_pair("wE",to_string(nbByteWE)));
				paramList.push_back(make_pair("plainVHDL","true")); 
				testStateList.push_back(paramList);
				paramList.clear();
			}
			for(int wF=5; wF<53; wF+=1) // test various input widths
			{ 
				int nbByteWE = 6+(wF/10);
				while(nbByteWE>wF){
					nbByteWE -= 2;
				}

				paramList.push_back(make_pair("wF",to_string(wF)));
				paramList.push_back(make_pair("wE",to_string(nbByteWE)));
				testStateList.push_back(paramList);
				paramList.clear();
			}
		}
		else     
		{
				// finite number of random test computed out of index
		}	

		return testStateList;
	}



	
	OperatorPtr FPLog::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int wE;
		UserInterface::parseStrictlyPositiveInt(args, "wE", &wE);
		int wF;
		UserInterface::parseStrictlyPositiveInt(args, "wF", &wF);
		int method;
		UserInterface::parseInt(args, "method", &method);
		int inTableSize;
		UserInterface::parseInt(args, "inTableSize", &inTableSize);
		if(method==0)
			return new FPLogIterative(parentOp, target, wE, wF, inTableSize);
		else if(method==1)
			{
				throw(string("TODO for Correntin"));
			}
		else 
			{
				throw(string("FPLog: the method parameter should currently be 0 or 1"));
			}
	}



	
	void FPLog::registerFactory(){
		UserInterface::add("FPLog", // name
											 "Floating-point logarithm",
											 "ElementaryFunctions", // categories
											 "",
											 "wE(int): exponent size in bits; \
                        wF(int): mantissa size in bits; \
                        method(int)=0: 0 for iterative, 1 for polynomial; \
                        inTableSize(int)=0: The input size to the tables of the iterative method, in bits, between 6 and 16. 0 choses a a sensible value",
											 "For details on the technique used, see <a href=\"bib/flopoco.html#DetDinPuj2007:Arith\">this article</a> and <a href=\"bib/flopoco.html#2010-RR-FPLog\">this research report</a>.",
											 FPLog::parseArguments,
											 FPLog::unitTest
											 ) ;

	}

}
