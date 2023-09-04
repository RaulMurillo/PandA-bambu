#ifndef INT2POSIT_HPP
#define INT2POSIT_HPP

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco{
  class Int2Posit : public Operator {
	public:
		/**
		 * @brief The  constructor
		 * @param[in]		parentOp	parent operator in the instance hierarchy
		 * @param[in]		target		the target device
		 * @param[in]	    widthI	    the total width of the input integer
		 * @param[in]		widthO		the total width of the output posit for the conversion result
		 * @param[in]	    wES		    the exponent field size of the output posit
         * @param[in]		trunc_p	    the output is not rounded when trunc_p is true
		 */
        Int2Posit(OperatorPtr parentOp, Target* target, int widthI, int widthO, int wES, bool trunc_p=false);

		/**
		 * Int2Posit destructor
		 */
		~Int2Posit();


		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

		/** Factory register method */ 
		static void registerFactory();

		/** emulate() function to be shared by various implementations */
		void emulate(TestCase * tc);

	private:
        int width_;
        int wES_;
        int wE_;
        int wF_;
        int widthINT_;
        bool trunc_p;

	};
}
#endif
