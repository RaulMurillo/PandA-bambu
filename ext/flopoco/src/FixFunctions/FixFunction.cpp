/*
  FixFunction object for FloPoCo

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Aric team at Ecole Normale Superieure de Lyon
	then by the Socrate team at INSA de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,

  All rights reserved.

*/

#include "FixFunction.hpp"
#include <sstream>

namespace flopoco{


	FixFunction::FixFunction(string sollyaString_, bool signedIn_, int lsbIn_, int lsbOut_):
		sollyaString(sollyaString_), lsbIn(lsbIn_),  lsbOut(lsbOut_),  signedIn(signedIn_)
	{
		if(signedIn)
			wIn=-lsbIn+1; // add the sign bit at position 0
		else
			wIn=-lsbIn;
		ostringstream completeDescription;
		completeDescription << sollyaString_;
		if(signedIn)
			completeDescription << " on [-1,1)";
		else
			completeDescription << " on [0,1)";


		// Now do the parsing in Sollya
		fS = sollya_lib_parse_string(sollyaString.c_str());
		/* If  parse error throw an exception */
		if (sollya_lib_obj_is_error(fS) || !sollya_lib_obj_is_function(fS))
			throw(string("FixFunction: Unable to parse input function: ")+sollyaString);

		initialize();

		wOut=msbOut-lsbOut+1;

		if(lsbIn!=0) {// we have an IO specification
			completeDescription << " for lsbIn=" << lsbIn << " (wIn=" << wIn << "), msbout=" << msbOut << ", lsbOut="
													<< lsbOut << " (wOut=" << wOut << "). ";
		}
		completeDescription << outputDescription;
		description = completeDescription.str();
}




	FixFunction::FixFunction(sollya_obj_t fS_, bool signedIn_):
		signedIn(signedIn_), fS(fS_)
	{
		initialize();
	}



void	FixFunction::initialize()
	{
#if 1 // this sometimes enlarges the interval. 
		if(signedIn)
			inputRangeS = sollya_lib_parse_string("[-1;1]");
		else			
			inputRangeS = sollya_lib_parse_string("[0;1]");
#else // This is tighter : interval is [O, 1-1b-l]  but it causes more problems than it solves 
		string maxvalIn="1-1b"+to_string(lsbIn);
		ostringstream uselessNoise;
		uselessNoise << "[" << (signedIn?"-1":"0") << ";" << maxvalIn << "]";
		inputRangeS = sollya_lib_parse_string(uselessNoise.str().c_str());
#endif
		
		sollya_obj_t outIntervalS = sollya_lib_evaluate(fS,inputRangeS);
		sollya_obj_t supS = sollya_lib_sup(outIntervalS);
		sollya_obj_t infS = sollya_lib_inf(outIntervalS);
		mpfr_t supMP, infMP, tmp;
		mpfr_init2(supMP, 1000); // no big deal if we are not accurate here 
		mpfr_init2(infMP, 1000); // no big deal if we are not accurate here 
		mpfr_init2(tmp, 1000); // no big deal if we are not accurate here 
		sollya_lib_get_constant(supMP, supS);
		sollya_lib_get_constant(infMP, infS);
		if(mpfr_sgn(infMP) >= 0)	{
			signedOut=false;
			}
		else {
			signedOut=true;
		}
		ostringstream t; // write it before we take the absolute value below
		t << "Out interval: [" << mpfr_get_d(infMP,MPFR_RNDD) << "; "<< mpfr_get_d(supMP,MPFR_RNDU) << "]";
		//		cerr << " " << t.str() << endl;
		// Now recompute the MSB explicitely.
		mpfr_abs(supMP, supMP, GMP_RNDU);
		mpfr_abs(infMP, infMP, GMP_RNDU);
		mpfr_max(tmp, infMP, supMP, GMP_RNDU);
		mpfr_log2(tmp, tmp, GMP_RNDU);
		mpfr_floor(tmp, tmp);
		msbOut = mpfr_get_si(tmp, GMP_RNDU);
		if(signedOut)
			msbOut++;

		t << ". Output is " << (signedOut?"signed":"unsigned");
		outputDescription=t.str();

		sollya_lib_clear_obj(outIntervalS);
		sollya_lib_clear_obj(supS);
		sollya_lib_clear_obj(infS);
 		mpfr_clears(supMP,infMP,tmp, NULL);
	}

	FixFunction::~FixFunction()
	{
	  sollya_lib_clear_obj(fS);
	  sollya_lib_clear_obj(inputRangeS);
	  sollya_lib_clear_obj(outputRangeS);
	}

	string FixFunction::getDescription() const
	{
		return description;
	}

	void FixFunction::eval(mpfr_t r, mpfr_t x) const
	{
		sollya_lib_evaluate_function_at_point(r, fS, x, NULL);
	}


	double FixFunction::eval(double x) const
	{
		mpfr_t mpX, mpR;
		double r;

		mpfr_inits(mpX, mpR, NULL);
		mpfr_set_d(mpX, x, GMP_RNDN);
		sollya_lib_evaluate_function_at_point(mpR, fS, mpX, NULL);
		r = mpfr_get_d(mpR, GMP_RNDN);

		mpfr_clears(mpX, mpR, NULL);
		return r;
	}


	void FixFunction::eval(mpz_class x, mpz_class &rNorD, mpz_class &ru, bool correctlyRounded) const
	{
		int precision=100*(wIn+wOut);
		sollya_lib_set_prec(sollya_lib_constant_from_int(precision));

		mpfr_t mpX, mpR;
		mpfr_init2(mpX,wIn+2);
		mpfr_init2(mpR,precision);

		if(signedIn) {
			mpz_class negateBit = mpz_class(1) << (wIn);
			if ((x >> (-lsbIn)) !=0)
				x -= negateBit;
		}
		/* Convert x to an mpfr_t in [0,1[ */
		mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
		mpfr_div_2si(mpX, mpX, -lsbIn, GMP_RNDN);

		/* Compute the function */
		eval(mpR, mpX);
		//		REPORT(FULL,"function() input is:"<<sPrintBinary(mpX));
		//cerr << 100*(wIn+wOut) <<" function("<<mpfr_get_d(mpX, GMP_RNDN)<<") output before rounding is:"<<mpfr_get_d(mpR, GMP_RNDN) << " " ;
		/* Compute the signal value */
		mpfr_mul_2si(mpR, mpR, -lsbOut, GMP_RNDN);

		/* So far we have a highly accurate evaluation. Rounding to target size happens only now
		 */
		if(correctlyRounded){
			mpfr_get_z(rNorD.get_mpz_t(), mpR, GMP_RNDN);
			// convert to two's complement
			if(rNorD<0) {
				rNorD += (mpz_class(1)<<wOut);
			}
			
		}
		else{
			mpfr_get_z(rNorD.get_mpz_t(), mpR, GMP_RNDD);
			if(rNorD<0) {
				rNorD += (mpz_class(1)<<(wOut));
			}
			mpfr_get_z(ru.get_mpz_t(), mpR, GMP_RNDU);
			if(ru<0) {
				ru += (mpz_class(1)<<wOut);
			}
		}

		//		REPORT(FULL,"function() output r = ["<<rd<<", " << ru << "]");
		mpfr_clear(mpX);
		mpfr_clear(mpR);
	}





	void FixFunction::emulate(TestCase * tc, bool correctlyRounded){
			mpz_class x = tc->getInputValue("X");
			mpz_class rNorD,ru;
			eval(x,rNorD,ru,correctlyRounded);
			//cerr << " x=" << x << " -> " << rNorD << " " << ru << endl; // for debugging
			tc->addExpectedOutput("Y", rNorD);
			if(!correctlyRounded)
				tc->addExpectedOutput("Y", ru);
	}
} //namespace
