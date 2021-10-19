#ifndef Normalizer_HPP
#define Normalizer_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

/* 

Use cases:
  usually wR=wX (standard FP operators) and then computeSticky doesn't make sense.
  wR>wX does not make sense: this operator removes meaningless zeroes, it doesn't add information
	wR<wX makes sense at least in one circumstance: FPLargeAcc.
  	There we want to normalize the large accumulator into a (much smaller) FP format.
		Only in such case does computeSticky make sense: it is computed out of the bits to the right of what we keep.
		Optimizing the sticky computation should be low priority, though: who cares about ensuring the ties to even rule in such a non-standard operator? 

	
Refactoring needed to change the interface from wCount to maxShift
It will impact the following operators (plus FPLarge Acc)

src/Conversions/Posit2FP.o
src/Conversions/Fix2FP.o
src/FPAddSub/FPAddDualPath.o
src/FPAddSub/FPAddSinglePath.o
src/ExpLog/FPLogIterative.o
src/Conversions/Posit2PIF.o
src/Posit/Add/PIFAdd.o


For ref here is the perf on 7k70t-fbv484 with old version of 
./flopoco FPAdd frequency=1 we=8 wf=23
 329 LUTs 11.709ns 


*/

namespace flopoco{

	/** 
	 * A leading zero/one counter + shifter + sticky bit computer for FloPoCo
	 */ 
	class Normalizer : public Operator
	{
	public:
	
		/**  
		 *  Normalizer constructor, used in FPLog
		 * @param[in] target the target device for this operator
		 * @param[in] wX the width of the mantissa input
		 * @param[in] wR the width of the mantissa output
		 * @param[in] wCount the numbers of bits to count, often equal to wX but sometimes less (see FPLog)
		 * @param[in] computeSticky Should the operator compute a sticky bit out of the shifted-out bits?
		 * @param[in] countType 0: count zeroes, 1: count ones; -1: have a dynamic OZb input that tells what to count 
		 */
		Normalizer(OperatorPtr parentOp, Target* target, int wX, int wR, int maxShift, bool computeSticky=false, const int countType=-1);
	
		/** The Normalizer destructor */
		~Normalizer();

	
		/** Returns the number of bits of the count
		 * @return the number of bits of the count
		 */
		int getCountWidth() const;
	
	
	
	
		void emulate(TestCase* tc);


				/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

		/** Factory register method */ 
		static void registerFactory();

	private:
	
		int          wX;                   /**< The number of bits of the input */
		int          wR;                  /**< The number of bits of the shifted output */
		int          maxShift;              /**< max number of bits to count */ 
		int          wCount;                /**< The number of bits of the count output */
		int          countType;             /**< -1|0|1. If 0, count zeroes (LZC). If 1, count ones (LOC). If -1, generic LZOC is instantiated, with an input stating what to count */
		bool         computeSticky;         /**<  if true, a sticky bit will be computed and output */
	};
}
#endif
