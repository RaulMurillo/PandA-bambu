/*
  Polynomial Function Evaluator for FloPoCo
	This version uses piecewise polynomial approximation for a trade-off between tables and multiplier hardware.
	
  Authors: Florent de Dinechin (rewrite from scratch of original code by Mioara Joldes and Bogdan Pasca, see the attic directory)

  This file is part of the FloPoCo project
	launched by the Arénaire/AriC team of Ecole Normale Superieure de Lyon
  currently developed by the Socrate team at CITILab/INSA de Lyon
 
  Initial software.
  Copyright © ENS-Lyon, INSA-Lyon, INRIA, CNRS, UCBL,  
  2008-2014.
  All rights reserved.

  */

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <limits.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include <gmpxx.h>
#include "../utils.hpp"


#include "FixFunctionByPiecewisePoly.hpp"
#include "Table.hpp"
#include "FixFunctionByTable.hpp"
#include "FixHornerEvaluator.hpp"

using namespace std;

namespace flopoco{

	/* TODO for the brave
		 In many cases, we have some of the sigma_i who have a constant sign (see the reports).
		 This info is stored in sigmaSign, and currently unused.
		 In such cases, we could get back to reducing x to unsigned, and get back to fully unsigned arithmetic, at least in some of the steps.
		 There is probably a market for that.
		 This needs to be done cleanly (i.e. without bloating the code/making it unreadable)
		 Maybe we should wait for power-of-two segmentation

	 */

	
	/* Error analysis:
		 Target error is exp2(lsbOut).
		 Final rounding may entail up to exp2(lsbOut-1).
		 So approximation+rounding error budget is  exp2(lsbOut-1).
		 We call PiecewisePolyApprox with approximation target error bound exp2(lsbOut-2).
		 It reports pwp->LSB: may be lsbOut-2, may be up to lsbOut-2- intlog2(degree) 
		 It also reports pwp->approxErrorBound
		 so rounding error budget is  exp2(lsbOut-1) - pwp->approxErrorBound
		 Staying within this budget is delegated to FixHornerEvaluator: see there for the details
			 
	 */

	
#define DEBUGVHDL 0


	FixFunctionByPiecewisePoly::FixFunctionByPiecewisePoly(OperatorPtr parentOp, Target* target, string func, int lsbIn_, int lsbOut_, int degree_, bool finalRounding_, double approxErrorBudget_):
		Operator(parentOp, target), degree(degree_), lsbIn(lsbIn_), lsbOut(lsbOut_), finalRounding(finalRounding_), approxErrorBudget(approxErrorBudget_){

		srcFileName="FixFunctionByPiecewisePoly";
		setNameWithFreqAndUID("FixFunctionByPiecewisePoly"); 
		setCopyrightString("Florent de Dinechin (2014-2020)");

		if(finalRounding==false){
			THROWERROR("FinalRounding=false not implemented yet" );
		}
		useNumericStd();

		f=new FixFunction(func, false, lsbIn, lsbOut); // this will provide emulate etc.
		msbOut = f->msbOut;
		signedOut = f->signedOut;
		
		addHeaderComment("Evaluator for " +  f-> getDescription() + "\n");
		REPORT(DETAILED, "Entering constructor, FixFunction description: " << f-> getDescription());
 		int wX=-lsbIn;
		addInput("X", wX); // TODO manage signedIn

		int outputSize = msbOut-lsbOut+1; // OK both for signed and unsigned out 
		addOutput("Y" ,outputSize , 2);

		if(degree==0){ // This is a simple table
			REPORT(DETAILED, "Degree 0: building a simple table");
			ostringstream params;
			params << "f="<<func
						 << " signedIn=false" << " lsbIn=" << lsbIn
						 << " lsbOut=" << lsbOut; 
			newInstance("FixFunctionByTable",
									"simpleTable",
									params.str(),
									"X=>X",
									"Y=>YY" );

			vhdl << tab << "Y <= YY;" << endl; 
		}
		else{
			if(degree==1){ // For degree 1, MultiPartite could work better
			REPORT(INFO, "Degree 1: You should consider using FixFunctionByMultipartiteTable");
			}

			
			// Build the polynomial approximation
			double targetAcc= approxErrorBudget*pow(2, lsbOut);
			REPORT(INFO, "Computing polynomial approximation of degree " << degree << " for target accuracy "<< targetAcc);
			pwp = new UniformPiecewisePolyApprox(func, targetAcc, degree);
			alpha =  pwp-> alpha; // coeff table input size 
						
			// Build the coefficient table out of the vector of polynomials. This is also where we add the final rounding bit
			buildCoeffTable();

			// What remains of the error budget for the evaluation phase ?
			double roundingErrorBudget=exp2(lsbOut-1)-pwp->approxErrorBound;
			// (this is actually useless because we will jointly optimize approx+round errors on each subinterval)
			REPORT(INFO, "Overall error budget = " << exp2(lsbOut) << "  of which approximation error = " << pwp->approxErrorBound
						 << " hence rounding error budget = "<< roundingErrorBudget );

			// The VHDL that splits the input into A and Z
			vhdl << tab << declare("A", alpha)  << " <= X" << range(wX-1, wX-alpha) << ";" << endl;
			vhdl << tab << declare("Z", wX-alpha)  << " <= X" << range(wX-alpha-1, 0) << ";" << endl; // called Y in the book
			vhdl << tab << declare("Zs", wX-alpha)  << " <= (not Z(" << wX-alpha-1 << ")) & Z" << range(wX-alpha-2, 0) << "; -- centering the interval" << endl;

			Table::newUniqueInstance(this, "A", "Coeffs",
															 coeffTableVector, "coeffTable", alpha, polyTableOutputSize );
			
			addComment(" Split the table output into each coefficient, adding back the constant signs if any");
			int currentShift=0;
			for(int i=pwp->degree; i>=0; i--) {
				int actualSize = pwp->MSB[i] - pwp->LSB + (pwp->coeffSigns[i]==0 ? 1 : 0);
				vhdl << tab << declare(join("A",i), pwp->MSB[i] - pwp->LSB +1)
						 << " <= ";
				if (pwp->coeffSigns[i]!=0) // add constant sign back
					vhdl << (pwp->coeffSigns[i]==1? "\"0\"" :  "\"1\"") << " & " ;				
				vhdl << "Coeffs" << range(currentShift + actualSize-1, currentShift) << ";" << endl;
				currentShift += actualSize;
			}

			
			// What follows is related to Horner evaluator
			// Here I wish I could plug other (more parallel) evaluators. 


			REPORT(INFO, "Now building the Horner evaluator for rounding error budget "<< roundingErrorBudget);

		// This is the same order as newInstance() would do, but does not require to write a factory for this Operator, which wouldn't make sense
		schedule();
		inPortMap("Y", "Zs");
		for(int i=0; i<=pwp->degree; i++) {
			inPortMap(join("A",i),  join("A",i));
		}
		outPortMap("R", "HornerOutput");
		OperatorPtr h = new  FixHornerEvaluator(this, target, 
																						lsbIn+alpha+1,
																						msbOut,
																						lsbOut,
																						pwp->poly // provides degree and coeff formats
																						// do we need to pass the constant signs of the coefficients? No, they have been added back
																						// and everybody is signed.
																						);
		vhdl << instance(h, "Horner", false);
			
		vhdl << tab << "Y <= " << "std_logic_vector(HornerOutput);" << endl;

		}
	}



	FixFunctionByPiecewisePoly::~FixFunctionByPiecewisePoly() {
		delete f;
	}





	
	void FixFunctionByPiecewisePoly::buildCoeffTable() {
		// First compute the table output size
		polyTableOutputSize=0;
		for (int i=0; i<=degree; i++) {
			polyTableOutputSize += pwp->MSB[i] - pwp->LSB + (pwp->coeffSigns[i]==0? 1 : 0);
		} 
		REPORT(DETAILED, "Poly table input size  = " << alpha);
		REPORT(DETAILED, "Poly table output size = " << polyTableOutputSize);

		int x;
		for(x=0; x<(1<<alpha); x++) {
			mpz_class z=0;
			int currentShift=0;
			for(int i=pwp->degree; i>=0; i--) {
				mpz_class coeff = pwp-> getCoeffAsPositiveMPZ(x, i); // coeff of degree i from poly number x
				if (pwp->coeffSigns[i] != 0) {// sign is constant among all the coefficients: remove it from here, it will be added back as a constant in the VHDL
					mpz_class mask = (mpz_class(1)<<(pwp->MSB[i] - pwp->LSB) ) - 1; // size is msb-lsb+1
					coeff = coeff & mask; 
				}
				z += coeff << currentShift; // coeff of degree i from poly number x
				// REPORT(DEBUG, "i=" << i << "   z=" << unsignedBinary(z, 64));
				if(i==0 && finalRounding){ // coeff of degree 0
					int finalRoundBitPos = lsbOut-1;
					z += mpz_class(1)<<(currentShift + finalRoundBitPos - pwp->LSB); // add the round bit
					//REPORT(DEBUG, "i=" << i << " + z=" << unsignedBinary(z, 64));
					// This may entail an overflow of z, in the case -tiny -> +tiny, e.g. first polynomial of atan
					// This is OK modulo 2^wOut (two's complement), but will break the vhdl output of Table: fix it here.
					z = z & ((mpz_class(1)<<polyTableOutputSize) -1);
					// REPORT(INFO, "Adding final round bit at position " << finalRoundBitPos-pwp->LSB);
				}
				currentShift +=  pwp->MSB[i] - pwp->LSB + (pwp->coeffSigns[i]==0? 1: 0);
			}
			coeffTableVector.push_back(z);
		}

	}



	
	

 


	void FixFunctionByPiecewisePoly::emulate(TestCase* tc){
		f->emulate(tc);
	}

	void FixFunctionByPiecewisePoly::buildStandardTestCases(TestCaseList* tcl){
		TestCase *tc;

		tc = new TestCase(this); 
		tc->addInput("X", 0);
		emulate(tc);
		tcl->add(tc);
		tc = new TestCase(this); 
		tc->addInput("X", (mpz_class(1) << f->wIn) -1);
		emulate(tc);
		tcl->add(tc);
	}


	TestList FixFunctionByPiecewisePoly::unitTest(int index)
	{
		// the static list of mandatory tests
		TestList testStateList;
		vector<string> functionList;
		functionList.push_back("sin(pi/2*x)");
		functionList.push_back("sin(x)");
		functionList.push_back("exp(x)");
		functionList.push_back("log(x+1)");
		functionList.push_back("log(x+0.5)"); // a function that gets negative
		functionList.push_back("tanh(4*x)");

		vector<pair<string,string>> paramList;
		if(index==-1) 
			{ // The unit tests
				for (size_t i=0; i<functionList.size(); i++) {
					// first deg 2 and 3, 15 bits, exhaustive test, then deg 5 for 25 bits 
					string f = functionList[i];
					paramList.push_back(make_pair("f","\"" + f + "\""));
					paramList.push_back(make_pair("plainVHDL","true"));
					paramList.push_back(make_pair("lsbOut","-15"));
					paramList.push_back(make_pair("lsbIn","-15"));
					paramList.push_back(make_pair("d","2"));
					paramList.push_back(make_pair("TestBench n=","-2"));
					testStateList.push_back(paramList);
					paramList.clear();

					paramList.push_back(make_pair("f","\"" + f + "\""));
					paramList.push_back(make_pair("plainVHDL","true"));
					paramList.push_back(make_pair("lsbOut","-24"));
					paramList.push_back(make_pair("lsbIn","-24"));
					paramList.push_back(make_pair("d","3"));
					testStateList.push_back(paramList);
					paramList.clear();

					paramList.push_back(make_pair("f","\"" + f + "\""));
					paramList.push_back(make_pair("plainVHDL","true"));
					paramList.push_back(make_pair("lsbOut","-30"));
					paramList.push_back(make_pair("lsbIn","-30"));
					paramList.push_back(make_pair("d","5"));
					testStateList.push_back(paramList);
					paramList.clear();

					paramList.push_back(make_pair("f","\"" + f + "\""));
					paramList.push_back(make_pair("plainVHDL","true"));
					paramList.push_back(make_pair("lsbOut","-32"));
					paramList.push_back(make_pair("lsbIn","-32"));
					paramList.push_back(make_pair("d","5"));
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

	
	OperatorPtr FixFunctionByPiecewisePoly::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
		int lsbIn, lsbOut, d;
		string f;
		double approxErrorBudget;
		UserInterface::parseString(args, "f", &f); 
		UserInterface::parseInt(args, "lsbIn", &lsbIn);
		UserInterface::parseInt(args, "lsbOut", &lsbOut);
		UserInterface::parsePositiveInt(args, "d", &d);
		UserInterface::parseFloat(args, "approxErrorBudget", &approxErrorBudget);
		return new FixFunctionByPiecewisePoly(parentOp, target, f, lsbIn, lsbOut, d, true, approxErrorBudget);
	}

	void FixFunctionByPiecewisePoly::registerFactory(){
		UserInterface::add("FixFunctionByPiecewisePoly", // name
											 "Evaluator of function f on [0,1), using a piecewise polynomial of degree d with Horner scheme.",
											 "FunctionApproximation",
											 "",
											 "f(string): function to be evaluated between double-quotes, for instance \"exp(x*x)\";\
                        lsbIn(int): weight of input LSB, for instance -8 for an 8-bit input;\
                        lsbOut(int): weight of output LSB;\
                        d(int): degree of the polynomial;\
                        approxErrorBudget(real)=0.25: error budget in ulp for the approximation, between 0 and 0.5",                        
											 "This operator uses a table for coefficients, and Horner evaluation with truncated multipliers sized just right.<br>For more details, see <a href=\"bib/flopoco.html#DinJolPas2010-poly\">this article</a>.",
											 FixFunctionByPiecewisePoly::parseArguments,
											 FixFunctionByPiecewisePoly::unitTest
											 ) ;
		
	}

}
	

