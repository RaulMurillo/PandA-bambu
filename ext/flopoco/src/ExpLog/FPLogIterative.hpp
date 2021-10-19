#ifndef __FPLOGITERATIVE_HPP
#define __FPLOGITERATIVE_HPP
#include <vector>
#include <sstream>

#include "FPLog.hpp"


namespace flopoco{

#define MAXRRSTAGES 2000 // 4000 bits of accuracy should be enough for anybody

	class FPLogIterative : public FPLog
	{
	protected:


	public:
		FPLogIterative(OperatorPtr parentOp, Target* target, int wE, int wF, int inTableSize=0);
		~FPLogIterative();

	private:
		/** The input sizes to the successive tables*/
		int a[MAXRRSTAGES];

		/** The intermediate precision: at step i, the exp part is bounded by
			 1 <= m < 1+2^-p[i]*/
		int p[MAXRRSTAGES];

		/** The useful size of the product Ai*Zi, and hence of the subword of Zi input to the mult.  */
		int psize[MAXRRSTAGES];


		/** size of the virtual mantissa:
			 1.000(p[i] zeroes)00000zzzzzzz
			 the zzzz are the bits of Zi, which will be actually computed  */
		int sfullZ[MAXRRSTAGES];

		/** size of Zi before truncation, should be equal to s[i]+t[i] */
		int sbt[MAXRRSTAGES];

		/** The number of non-zero bits of Zi, will be limited to wF+g */
		int s[MAXRRSTAGES];

		/** The numbers of bits of the full Zi to truncate to limit the size to wF+g */
		int t[MAXRRSTAGES];


		// The max number of stages
		int stages;

		int sfinal; // equal to s[stages+1]

		int pfinal;  // equal to s[stages+1]

		// The guard bits
		int gLog;


		// A global boolean flag disabling the simulation of the fullsize mantissa
		// as soon as it would be more than 64 bits
		int fullSizeSim;

		// The target precision: numbers may be truncated so that their LSB has weight -target_prec
		int target_prec;

	};
}
#endif
