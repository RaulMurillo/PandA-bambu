#pragma once

#include "TilingStrategy.hpp"
#include "./../BitHeap/CompressionStrategy.hpp"

#ifdef HAVE_SCALP
#include <ScaLP/Solver.h>
#include <ScaLP/Exception.h>    // ScaLP::Exception
#include <ScaLP/SolverDynamic.h> // ScaLP::newSolverDynamic
#endif //HAVE_SCALP
#include <iomanip>
#include "BaseMultiplier.hpp"
#include "BaseMultiplierDSPSuperTilesXilinx.hpp"
#include "MultiplierTileCollection.hpp"

namespace flopoco {

/*!
 * The TilingAndCompressionOptILP class
 */
class TilingAndCompressionOptILP : public TilingStrategy, public CompressionStrategy
{

public:
    using TilingStrategy::TilingStrategy;
    TilingAndCompressionOptILP(
        unsigned int wX,
        unsigned int wY,
        unsigned int wOut,
        bool signedIO,
        BaseMultiplierCollection* bmc,
		base_multiplier_id_t prefered_multiplier,
		float occupation_threshold,
		int maxPrefMult,
        MultiplierTileCollection tiles_,
        BitHeap* bitheap,
        unsigned guardBits,
        unsigned keepBits,
        mpz_class errorBudget,
        mpz_class &centerErrConstant,
        bool performOptimalTruncation,
        bool minStages);

    void solve() override;
	void compressionAlgorithm() override;

private:
    float occupation_threshold_;
    int dpX, dpY, dpS, wS, dpK, dpC, dpERC, dpSt, s_max, max_pref_mult_;
    unsigned prodWidth, guardBits, keepBits;
    mpz_class &centerErrConstant, eBudget;
    unsigned long long  errorBudget;
	vector<BaseMultiplierCategory*> tiles;
	bool performOptimalTruncation, minLutsFocus;
#ifdef HAVE_SCALP
    BasicCompressor* flipflop;
    void constructProblem(int s_max);
    bool addFlipFlop();
    ScaLP::Result bestResult;
    ScaLP::Solver *solver;

    void resizeBitAmount(unsigned int stages);

     /**
     * @brief Writes a file "result.txt" containing the ILP-solution
     * @param solver reference to solver
     * @return none
     */
    void writeSolutionFile(ScaLP::Result result);

    /**
     * @brief Checks if a tiling for a truncated multiplier meets the error budget as required for faithfulness
     * @param solution list of the placed tiles with their parametrization and anchor point
     * @param guardBits the number of bits below the output LSB that we need to keep in the summation
     * @param errorBudget maximal permissible weight of the sum of the omitted partial products (as they would appear in an array multiplier)
     * @param constant to recenter the truncation error around 0 since it can otherwise only be negative, since there are only partial products left out. This allows a larger error, so more products can be omitted
     * @return true when the error bound is met, otherwise false
     */
    bool checkTruncationError(list<TilingStrategy::mult_tile_t> &solution, unsigned int guardBits, mpz_class errorBudget, mpz_class constant);


    void parseTilingSolution(long long &optTruncNumericErr);

#endif


    };

}
