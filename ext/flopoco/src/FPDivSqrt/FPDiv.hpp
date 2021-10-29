#ifndef FPDIV_HPP
#define FPDIV_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "TestBenches/FPNumber.hpp"

namespace flopoco{

	/** The FPDiv class */
	class FPDiv : public Operator
	{
	public:
		/**
		 * The FPDiv constructor
		 * @param[in]		target		the target device
		 * @param[in]		wE			the width of the exponent for the f-p number X
		 * @param[in]		wF			the width of the fraction for the f-p number X
		 */
		FPDiv(OperatorPtr parentOp, Target* target, int wE, int wF, int radix=0);

		/**
		 * FPDiv destructor
		 */
		~FPDiv();


		/**
		 * compute the selection function table 
		 * See the Digital Arithmetic book by Ercegovac and Lang for details
     * Some of these parameters are computed by the utility NbBitsMin

		 * @param dMin 0.5  if there is no prescaling; In case of prescaling, can be larger
		 * @param dMax 1 if there is no prescaling; In case of prescaling, can be larger
		 * @param nbBitD  number of bits of d to input to the selection table
		 * @param number of bits of w to input to the selection table
		 * @param alpha digit set is {-alpha, .. alpha} 
		 * @param radix radix used
		 */

		vector<mpz_class> selFunctionTable(double dMin, double dMax, int nbBitD, int nbBitW, int alpha, int radix);
		
		/**
		 * Emulate a correctly rounded division using MPFR.
		 * @param tc a TestCase partially filled with input values
		 */
		void emulate(TestCase * tc);

		/* Overloading the Operator method */
		void buildStandardTestCases(TestCaseList* tcl);

		static TestList unitTest(int index);

		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);
		static void registerFactory();


		
	private:
		/** The width of the exponent for the input X */
		int wE;
		/** The width of the fraction for the input X */
		int wF;
		/** the radix */
		int radix;
		/** the digit set parameter: digits will be in {-alpha... alpha} */
		int alpha;
		/** The number of iterations */
		int nDigit;
		/** prescaling parameter: 0 means no prescaling; 1 means prescaling with 1 addition; 2 means prescaling with 2 additions */
		int prescaling;
		
	};


	class SRTDivNbBitsMin : public Operator //not an operator, should be disabled in Factories/ -- just to benefit from the lousy parser
	{
	private:
		static void plotPDDiagram(int delta, int t, int radix, int digitSet);
		static bool checkDistrib(int delta, int t, int radix, int digitSet);
		static double L(int k, double ro, double d);
		static double U(int k, double ro, double d);
		static double estimateCost(int nbBit, int radix, int digitSet);
		static void computeNbBit(int radix, int digitSet);
	public:
		static void registerFactory();
		static OperatorPtr NbBitsMinParseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);
	};

}
#endif //FPDIV_HPP
