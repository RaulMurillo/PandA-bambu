#ifndef FPMULTIPLIERS_HPP
#define FPMULTIPLIERS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "IntMult//IntMultiplier.hpp"
#include "IntAddSubCmp/IntAdder.hpp"

namespace flopoco{

	/** The FPMult class */
	class FPMult : public Operator
	{
	public:

		/**
		 * The FPMutliplier constructor
		 * @param[in]		target		the target device
		 * @param[in]		wEX			the width of the exponent for the f-p number X
		 * @param[in]		wFX			the width of the fraction for the f-p number X
		 * @param[in]		wEY			the width of the exponent for the f-p number Y
		 * @param[in]		wFY			the width of the fraction for the f-p number Y
		 * @param[in]		wER			the width of the exponent for the multiplication result
		 * @param[in]		wFR			the width of the fraction for the multiplication result
		 **/
		FPMult(OperatorPtr parentOp, Target* target, int wEX, int wFX, int wEY, int wFY, int wER, int wFR,
		             bool norm = true, bool correctlyRounded=true, float dspOccupationThreshold = 0.0);

		/**
		 * FPMult destructor
		 */
		~FPMult();

		/**
		 * Emulate the operator using MPFR.
		 * @param tc a TestCase partially filled with input values
		 */
		void emulate(TestCase * tc);

		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

		static void registerFactory();

	protected:

		int  wEX_;                  /**< The width of the exponent for the input X */
		int  wFX_;                  /**< The width of the fraction for the input X */
		int  wEY_;                  /**< The width of the exponent for the input Y */
		int  wFY_;                  /**< The width of the fraction for the input Y */
		int  wER_;                  /**< The width of the exponent for the output R */
		int  wFR_;                  /**< The width of the fraction for the output R */
		bool normalized_;	          /**< Signal if the output of the operator is to be or not normalized*/
		bool correctlyRounded_;	    /**< true: operator computes correct rounding; false: operator computes faithful rounding */
		float dspOccupationThreshold; /**< threshold of relative occupation ratio of a DSP multiplier to be used or not */

	};
}

#endif
