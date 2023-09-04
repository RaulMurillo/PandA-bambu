#ifndef Quire2Posit_HPP
#define Quire2Posit_HPP

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco
{
	class Quire2Posit : public Operator
	{
	public:
		/**
		 * @brief			The constructor
		 * @param[in]		parentOp	parent operator in the instance hierarchy
		 * @param[in]		target		the target device
		 * @param[in]		width		the total width of the output posit
		 * @param[in]		wES			the exponent field size of the output posit
		 * @param[in]		carry		extra carry guard bits to savely allow the sum up to 2^carry products
		 */
		Quire2Posit(OperatorPtr parentOp, Target *target, int width, int wES, int carry);

		/**
		 * Quire2Posit destructor
		 */
		~Quire2Posit();

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

		/** Factory register method */
		static void registerFactory();

		/** emulate() function to be shared by various implementations */
		void emulate(TestCase *tc);

	private:
		int width_;
		int wES_;
		int wE_;
		int wF_;
		int C_;
		int wQuire_;
	};
}
#endif
