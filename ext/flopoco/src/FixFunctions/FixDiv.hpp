#ifndef FIXDIV_HPP
#define FIXDIV_HPP

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco
{
	class FixDiv : public Operator
	{
	public:
		/**
		 * @brief 	The constructor
		 * @param[in]	parentOp		parent operator in the instance hierarchy
		 * @param[in]	target			the target device
		 * @param[in]	ints		    the total width of the inputs
		 * @param[in]	frac			the exponent field size of the inputs
		 * @param[in]   dspOccupationThreshold	the threshold of relative occupation ratio of a DSP multiplier
		 * @param[in]   truncate		the output will be truncated to have same length as inputs
		 */
		FixDiv(OperatorPtr parentOp, Target *target, int ints, int frac, float dspOccupationThreshold=0.0, bool truncate=false, int iters=0, bool useGoldschmidt=false, int LUT_out=10);

		/**
		 * FixDiv destructor
		 */
		~FixDiv();

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

		/** Factory register method */
		static void registerFactory();

		/** emulate() function to be shared by various implementations */
		void emulate(TestCase *tc);

		//void buildStandardTestCases(TestCaseList* tcl);

	private:
		int width_;
		int ints_;
		int frac_;
		int iters_;
		float dspOccupationThreshold; /**< threshold of relative occupation ratio of a DSP multiplier to be used or not */
        bool truncate_;
		// bool iterative;
		bool goldschmidt;
		int LUT_out_;
	};
}
#endif
