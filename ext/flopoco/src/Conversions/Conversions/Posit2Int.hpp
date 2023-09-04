#ifndef POSIT2INT_HPP
#define POSIT2INT_HPP

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "utils.hpp"
#include "Operator.hpp"

namespace flopoco{
  class Posit2Int : public Operator {
	public:
		/**
		 * @brief The  constructor
		 * @param[in]		parentOp	parent operator in the instance hierarchy
		 * @param[in]		target		the target device
		 * @param[in]	    widthI	    the total width of the input posit
		 * @param[in]	    esI		    the exponent field size of the input posit
		 * @param[in]		widthO		the total width of the output integer for the conversion result
         * @param[in]		trunc_p	    the output is not rounded when trunc_p is true
		 */
        Posit2Int(OperatorPtr parentOp, Target* target, int widthI, int wESI, int widthO, bool trunc_p=false);

		/**
		 * Posit2Int destructor
		 */
		~Posit2Int();


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
        int widthO_;
        bool trunc_p;

	};
}
#endif
