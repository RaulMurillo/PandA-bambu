#ifndef FP2POSIT_HPP
#define FP2POSIT_HPP

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco
{
	class FP2Posit : public Operator
	{
	public:
		/**
		 * @brief Constructor
		 * @param[in]	target	the target device
		 * @param[in]	wEF		the exponent field size of the input float
		 * @param[in]	wFF		the fraction field size of the input float
		 * @param[in]	widthP	the total width of the convertion posit result
		 * @param[in]	wESP	the exponent field size of the convertion posit result
		 */
		FP2Posit(Operator *parentOp, Target *target, int wEF, int wFF, int widthP, int wESP);

		// destructor
		// ~FP2Posit();

		void emulate(TestCase *tc);

		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);
		static void registerFactory();

	private:
		//int widthI_;
		int wEF_;
		int wFF_;
		int wE_;
		int wF_;
		int widthP_;
		int wESP_;

		static void computeFP2PositWidths(int const widthO, int const wES, int *wE, int *wF);
	};

} //namespace

#endif