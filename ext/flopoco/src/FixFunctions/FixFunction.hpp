#ifndef _FIXFUNCTION_HPP_
#define _FIXFUNCTION_HPP_

#include <string>
#include <iostream>

#include <sollya.h>
#include <gmpxx.h>
#include "../TestBenches/TestCase.hpp"

using namespace std;

/* Stylistic convention here: all the sollya_obj_t have names that end with a capital S */
namespace flopoco{
	
	/** The FixFunction objects holds a fixed-point function. 
			It provides an interface to Sollya services such as 
			parsing it,
			evaluating it at arbitrary precision,
			providing a polynomial approximation on an interval
	*/

	class FixFunction {
	public:


		/**
			 The FixFunctionByTable constructor
			 @param[string] func    a string representing the function
			 @param[bool]   signedIn: if true, input range is [0,1], else input range is [0,1]
			 @param[int]    lsbIn    input LSB weight (-lsbX is the input size)
			 @param[int]    lsbOut  output LSB weight

			 There are defaults for lsbOut and msbOut for situations when they are computed afterwards.
		 */
		FixFunction(string sollyaString, bool signedIn, int lsbIn=0, int lsbOut=0);
		FixFunction(sollya_obj_t fS, bool signedIn);

		virtual ~FixFunction();

		string getDescription() const;

		/** helper method: computes the value of the function on a double; no range check */
		double eval(double x) const;

		/** helper method: computes the value of the function on an MPFR; no range check */
		void eval(mpfr_t r, mpfr_t x) const;

		/** if correctlyRounded is true, compute the CR result in rNorD; otherwise computes the RD in rNorD and the RU in ru.
				x as well as the result(s) are given as bit vectors.
				These bit vectors are held in a mpz_class which is always positive, although it may hold a two's complement value.
				(this is the usual behaviour of emulate())
		*/
		void eval(mpz_class x, mpz_class &rNorD, mpz_class &ru, bool correctlyRounded=false) const;

		void emulate(TestCase * tc,	bool correctlyRounded=false /**< if true, correctly rounded RN; if false, faithful function */);

		// All the following public, not good practice I know, but life is complicated enough
		// All these public attributes are read at some point by Operator classes
		string sollyaString;
		int lsbIn;   
		int wIn;   
		int msbOut; /*< Now computed by the constructor; assuming signed out */
		int lsbOut;
		int wOut;
		bool signedIn;
		bool signedOut; /**< computed by the constructor */
		string description;
		sollya_obj_t fS;
		sollya_obj_t inputRangeS;  /**< computed by the constructor */
		sollya_obj_t outputRangeS; /**< computed by the constructor */
	private:
		void initialize();
		string outputDescription;
	};

}
#endif // _FIXFUNCTION_HH_
