#ifndef PositApproxDiv_HPP
#define PositApproxDiv_HPP

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco
{
	class PositApproxDiv : public Operator
	{
	public:
		/**
		 * @brief 		The constructor
		 * @param[in]	parentOp	parent operator in the instance hierarchy
		 * @param[in]	target		the target device
		 * @param[in]	width		the total width of the input posits
		 * @param[in]	wES			the exponent field size of the input posits
		 */
		PositApproxDiv(OperatorPtr parentOp, Target *target, int width, int wES);

		/**
		 * PositApproxDiv destructor
		 */
		~PositApproxDiv();

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

		/** Factory register method */
		static void registerFactory();

		/** emulate() function to be shared by various implementations */
		void emulate(TestCase *tc);

		void buildStandardTestCases(TestCaseList* tcl);

	private:
		int width_;
		int wES_;
		int wE_;
		int wF_;
	};
}
#endif
