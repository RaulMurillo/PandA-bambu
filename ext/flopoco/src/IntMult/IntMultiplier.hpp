#ifndef IntMultiplierS_HPP
#define IntMultiplierS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Table.hpp"
#include "BitHeap/BitHeap.hpp"

#include "IntMult/MultiplierBlock.hpp"
#include "IntMult/BaseMultiplierCategory.hpp"
#include "TilingStrategy.hpp"
#include "MultiplierTileCollection.hpp"

namespace flopoco {
	class IntMultiplier : public Operator {

	public:

		/**
		 * The IntMultiplier constructor
		 * @param[in] target           the target device
		 * @param[in] wX             X multiplier size (including sign bit if any)
		 * @param[in] wY             Y multiplier size (including sign bit if any)
		 * @param[in] wOut           wOut size for a truncated multiplier (0 means full multiplier)
		 * @param[in] signedIO       false=unsigned, true=signed
		 * @param[in] texOutput      true=generate a tek file with the found tiling solution
		 **/
		IntMultiplier(Operator *parentOp, Target* target, int wX, int wY, int wOut=0, bool signedIO = false, float dspOccupationThreshold=0.0, int maxDSP=-1, bool superTiles=false, bool use2xk=false, bool useirregular=false, bool useLUT=true, bool useDSP=true, bool useKaratsuba=false, int beamRange=0, bool optiTrunc=true, bool minStages=true);

		/**
		 * The emulate function.
		 * @param[in] tc               a test-case
		 */
		void emulate ( TestCase* tc );

		void buildStandardTestCases(TestCaseList* tcl);

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

		/** Factory register method */
		static void registerFactory();

		/**
		 * @brief Compute the size required to store the untruncated product of inputs of a given width
		 * @param wX size of the first input
		 * @param wY size of the second input
		 * @return the number of bits needed to store a product of I<wX> * I<WY>
		 */
		static unsigned int prodsize(unsigned int wX, unsigned int wY, bool signedX, bool signedY);

		static TestList unitTest(int index);

		BitHeap* getBitHeap(void) {return bitHeap;}

	protected:

		unsigned int wX;                         /**< the width for X after possible swap such that wX>wY */
		unsigned int wY;                         /**< the width for Y after possible swap such that wX>wY */
		unsigned int wFullP;                     /**< size of the full product: wX+wY  */
		unsigned int wOut;                       /**< size of the output, to be used only in the standalone constructor and emulate.  */
		bool signedIO;                   /**< true if the IOs are two's complement */
		bool negate;                    /**< if true this multiplier computes -xy */
		float dspOccupationThreshold;   /**< threshold of relative occupation ratio of a DSP multiplier to be used or not */
		int maxDSP;            /**< limit the number of DSP-Blocks used in multiplier */
		BitHeap *bitHeap;

	private:
//		Operator* parentOp;  			/**< For a virtual multiplier, adding bits to some external BitHeap,
		/**
		 * Realise a tile by instantiating a multiplier, selecting the inputs and connecting the output to the bitheap
		 *
		 * @param tile: the tile to instentiate
		 * @param idx: the tile identifier for unique name
		 * @param output_name: the name to which the output of this tile should be mapped
		 */
		Operator* realiseTile(
				TilingStrategy::mult_tile_t & tile,
				size_t idx,
				string output_name
			);

		/** returns the amount of consecutive bits, which are not constantly zero
		 * @param bm:                          current BaseMultiplier
		 * @param xPos, yPos:                  position of lower right corner of the BaseMultiplier
		 * @param totalOffset:                 see placeSingleMultiplier()
		 * */
		unsigned int getOutputLengthNonZeros(
				BaseMultiplierParametrization const & parameter,
				unsigned int xPos,
				unsigned int yPos,
				unsigned int totalOffset
			);

		unsigned int getLSBZeros(
				BaseMultiplierParametrization const & parameter,
				unsigned int xPos,
				unsigned int yPos,
				unsigned int totalOffset,
				int mode
			);

		/**
		 * @brief Compute the number of bits below the output msb that we need to keep in the summation
		 * @param wX first input width
		 * @param wY second input width
		 * @param wOut number of bits kept in the output
		 * @return the the number of bits below the output msb that we need to keep in the summation to ensure faithful rounding
		 */
		unsigned int computeGuardBits(unsigned int wX, unsigned int wY, unsigned int wOut);

        /**
         * @brief Compute several parameters for a faithfully rounding truncated multiplier
         * @param wFull width of result of a non-truncated multiplier with the same input widths
         * @param wOut requested output width of the result vector of the truncated multiplier
         * @param g the number of bits below the output LSB that we need to keep in the summation
         * @param k number of bits to keep in in the column with weight w-g
         * @param errorBudget maximal permissible weight of the sum of the omitted partial products (as they would appear in an array multiplier)
         * @param constant to recenter the truncation error around 0 since it can otherwise only be negative, since there are only partial products left out. This allows a larger error, so more products can be omitted
         * @return none
         */
        void computeTruncMultParams(unsigned wFull, unsigned wOut, unsigned &g, unsigned &k, unsigned long long &errorBudget, unsigned long long &constant);

        void computeTruncMultParamsMPZ(unsigned wFull, unsigned wOut, unsigned &g, unsigned &k, mpz_class &errorBudget, mpz_class &constant);

		/**
		 * add a unique identifier for the multiplier, and possibly for the block inside the multiplier
		 */
		string addUID(string name, int blockUID=-1);

		int multiplierUid;

        /**
         * @brief Define and calculate the size of the output signals of each multiplier tile in the solution to the BitHeap
         * @param bh BitHeap instance, where the partial results form the multiplier tiles should compressed
         * @param solution list of the placed tiles with their parametrization and anchor point
         * @param bitheapLSBWeight weight (2^bitheapLSBWeight) of the LSB that should be compressed on BH. It is supposed to be 0 for regular multipliers, but can be higher for truncated multipliers.
         * @return none
         */
		void branchToBitheap(BitHeap* bh, list<TilingStrategy::mult_tile_t> &solution , unsigned int bitheapLSBWeight);

        /**
         * @brief Checks if a tiling for a truncated multiplier meets the error budget as required for faithfulness
         * @param solution list of the placed tiles with their parametrization and anchor point
         * @param guardBits the number of bits below the output LSB that we need to keep in the summation
         * @param errorBudget maximal permissible weight of the sum of the omitted partial products (as they would appear in an array multiplier)
         * @param constant to recenter the truncation error around 0 since it can otherwise only be negative, since there are only partial products left out. This allows a larger error, so more products can be omitted
         * @return none
         */
        mpz_class checkTruncationError(list<TilingStrategy::mult_tile_t> &solution, unsigned int guardBits, const mpz_class& errorBudget, const mpz_class& constant) const;

        /**
         * @brief Calculate the LSB of the BitHeap required to maintain faithfulness, so that unnecessary LSBs to meet the error budget of multiplier tiles can be omitted from compression
         * @param solution list of the placed tiles with their parametrization and anchor point
         * @param guardBits the number of bits below the output LSB that we need to keep in the summation
         * @param errorBudget maximal permissible weight of the sum of the omitted partial products (as they would appear in an array multiplier)
         * @param constant to recenter the truncation error around 0 since it can otherwise only be negative, since there are only partial products left out. This allows a larger error, so more products can be omitted
         * @return none
         */
        static int calcBitHeapLSB(list<TilingStrategy::mult_tile_t> &solution, unsigned guardBits, const mpz_class& errorBudget, const mpz_class& constant, const mpz_class& actualTruncError);

        /**
         * @brief calculate the width of the diagonal of a rectangle, that is equivalent to the number of partial product bits in the hypothetical resulting bitheap of a Bough-Wooley-Multiplier
         * @param wX width of the rectangle
         * @param wY height of the rectangle
         * @param col sum of the coordinates (x+y) that describe the diagonal
         * @param wFull full product width of the rectangular multiplier
         * @return number of partial products in the diagonal
         */
        static unsigned int widthOfDiagonalOfRect(unsigned int wX, unsigned int wY, unsigned int col, unsigned wFull);
    };

}
#endif
