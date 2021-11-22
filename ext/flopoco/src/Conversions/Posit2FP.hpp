#ifndef POSIT2FP_HPP
#define POSIT2FP_HPP

#include "Operator.hpp"

namespace flopoco
{
	class Posit2FP : public Operator
	{
	public:
		/**
		 * @brief Constructor
		 * @param[in]	target	the target device
		 * @param[in]	widthP	the total width of the input posit
		 * @param[in]	wESP	the exponent field size of the input posit
		 * @param[in]	wEF		the exponent field size of the convertion float result
		 * @param[in]	wFF		the fraction field size of the convertion float result
		 */
		Posit2FP(Operator* parentOp, Target* target, int widthP, int wESP, int wEF, int wFF);

		void emulate(TestCase *tc);

		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);
		static void registerFactory();

	private:
		int widthP_;
		int wESP_;
		int wE_;
		int wF_;
		int wEF_;
		int wFF_;

		static void computePositWidths(int const widthI, int const esI, int* wE, int* wF);
	};
}

#endif
