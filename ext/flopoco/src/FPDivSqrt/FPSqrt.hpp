#ifndef FPSQRT_HPP
#define FPSQRT_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "TestBenches/FPNumber.hpp"

namespace flopoco{
	/** The FPSqrt class */
	class FPSqrt : public Operator
	{
	public:
		/**
		 * The FPSqrt constructor
		 * @param[in]		target		the target device
		 * @param[in]		wE			the the with of the exponent for the f-p number X
		 * @param[in]		wF			the the with of the fraction for the f-p number X
		 */
		FPSqrt(OperatorPtr parentOp, Target* target, int wE, int wF, int method=0);

		/**
		 * FPSqrt destructor
		 */
		~FPSqrt();



		/**
		 * Emulate a correctly rounded square root using MPFR.
		 * @param tc a TestCase partially filled with input values
		 */
		void emulate(TestCase * tc);

		/* Overloading the Operator method to limit testing of NaNs and negative numbers*/
		TestCase* buildRandomTestCase(int i);
		
		void buildStandardTestCases(TestCaseList* tcl);

		// User-interface stuff
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
		static TestList unitTest(int index);

		static void registerFactory();


	private:
		/** The width of the exponent for the input X */
		int wE;
		/** The width of the fraction for the input X */
		int wF;
		/** an int that selects the method */
		int method;
		/** A boolean selecting between IEEE-compliant correct rounding
			 or faithful (last-bit accurate) result  */
		bool correctRounding;
	}
		;
}
#endif //FPSQRT_HPP
