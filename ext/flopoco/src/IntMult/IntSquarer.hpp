#ifndef IntSquarerS_HPP
#define IntSquarerS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "Operator.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
#include "IntMult//IntMultiplier.hpp"

namespace flopoco{

	/** The IntSquarer class for experimenting with squarers. 
	 */
	class IntSquarer : public Operator
	{
	public:
		/**
		 * The IntSquarer constructor
		 * @param[in] target the target device
		 * @param[in] wIn    the with of the inputs and output
		 * @param[in] wOut the delays for each input; -1 means: exact squarer
		 **/
		IntSquarer(OperatorPtr parentOp, Target* target, int wIn, bool signedIn_=false, int wOut=0);

		/**
		 *  Destructor
		 */
		~IntSquarer();


		void emulate(TestCase* tc);

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

		/** Factory register method */ 
		static void registerFactory();

	protected:
		int wIn;                         /**< the width of X and Y */
		int signedIn;
		int wOut;                         /**< the actual width of R*/
		
	private:
		int guardBits;
		bool faithfulOnly;      /**< true for a truncated bit heap, but false for a tabulated square */
	};
}
#endif
