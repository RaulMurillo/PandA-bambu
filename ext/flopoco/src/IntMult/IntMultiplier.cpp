/*
  An integer multiplier for FloPoCo

  Authors:  Martin Kumm with parts of previous IntMultiplier of Bogdan Pasca, F de Dinechin, Matei Istoan, Kinga Illyes and Bogdan Popa spent

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2017.
  All rights reserved.
*/

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntMultiplier.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
#include "BaseMultiplierCollection.hpp"

#include "TilingStrategyOptimalILP.hpp"
#include "TilingStrategyBasicTiling.hpp"
#include "TilingStrategyGreedy.hpp"
#include "TilingStrategyXGreedy.hpp"
#include "TilingStrategyBeamSearch.hpp"
#include "TilingAndCompressionOptILP.hpp"

using namespace std;

namespace flopoco {


    IntMultiplier::IntMultiplier (Operator *parentOp, Target* target_, int wX_, int wY_, int wOut_, bool signedIO_, float dspOccupationThreshold, int maxDSP, bool superTiles, bool use2xk, bool useirregular, bool useLUT, bool useDSP, bool useKaratsuba, int beamRange, bool optiTrunc, bool minStages):
		Operator ( parentOp, target_ ),wX(wX_), wY(wY_), wOut(wOut_),signedIO(signedIO_), dspOccupationThreshold(dspOccupationThreshold) {
        srcFileName = "IntMultiplier";
        setCopyrightString("Martin Kumm, Florent de Dinechin, Kinga Illyes, Bogdan Popa, Bogdan Pasca, 2012");

        ostringstream name;
        name << "IntMultiplier";
        setNameWithFreqAndUID(name.str());

        // the addition operators need the ieee_std_signed/unsigned libraries
        useNumericStd();

        multiplierUid = parentOp->getNewUId();
        wFullP = prodsize(wX, wY, signedIO_, signedIO_);

        if (wOut == 0)
            wOut = prodsize(wX, wY, signedIO_, signedIO_);

        unsigned int guardBits = 0, keepBits = 0;
        mpz_class errorBudget = 0, centerErrConstant = 0;
        if(wFullP - wOut)
            computeTruncMultParamsMPZ(wFullP, wOut, guardBits, keepBits, errorBudget, centerErrConstant);
        //cout << " guardBits=" << guardBits << " keepBits=" << keepBits << " errorBudget=" << errorBudget << " centerErrConstant=" << centerErrConstant << endl;

        REPORT(INFO, "IntMultiplier(): Constructing a multiplier of size " <<
					wX << "x" << wY << " faithfully rounded to bit " << wOut <<
					". Will use " << guardBits << " guard bits and DSP threshhold of " << dspOccupationThreshold)

		string xname="X";
		string yname="Y";

		// Set up the IO signals
		addInput(xname, wX, true);
		addInput(yname, wY, true);
		addOutput("R", wOut, 2, true);

		// The larger of the two inputs
		vhdl << tab << declare(addUID("XX"), wX, true) << " <= " << xname << " ;" << endl;
		vhdl << tab << declare(addUID("YY"), wY, true) << " <= " << yname << " ;" << endl;

		if(target_->plainVHDL()) {
			vhdl << tab << declareFixPoint("XX",signedIO,-1, -wX) << " <= " << (signedIO?"signed":"unsigned") << "(" << xname <<");" << endl;
			vhdl << tab << declareFixPoint("YY",signedIO,-1, -wY) << " <= " << (signedIO?"signed":"unsigned") << "(" << yname <<");" << endl;
			vhdl << tab << declareFixPoint("RR",signedIO,-1, -wX-wY) << " <= XX*YY;" << endl;
			vhdl << tab << "R <= std_logic_vector(RR" << range(wX+wY-1, wX+wY-wOut) << ");" << endl;

			return;
		}

        BitHeap bitHeap(this, wFullP);
		//BitHeap bitHeap(this, wOut + guardBits);
        unsigned bitHeapLSBWeight = 0;

		REPORT(INFO, "Creating BaseMultiplierCollection")
//		BaseMultiplierCollection baseMultiplierCollection(parentOp->getTarget(), wX, wY);
		BaseMultiplierCollection baseMultiplierCollection(getTarget());

//		baseMultiplierCollection.print();


		MultiplierTileCollection multiplierTileCollection(getTarget(), &baseMultiplierCollection, wX, wY, superTiles, use2xk, useirregular, useLUT, useDSP, useKaratsuba);


		string tilingMethod = getTarget()->getTilingMethod();

		REPORT(INFO, "Creating TilingStrategy using tiling method " << tilingMethod)

		TilingStrategy* tilingStrategy;
		if(tilingMethod.compare("heuristicBasicTiling") == 0) {
			tilingStrategy = new TilingStrategyBasicTiling(
					wX,
					wY,
					wOut + guardBits,
					signedIO,
					&baseMultiplierCollection,
					baseMultiplierCollection.getPreferedMultiplier(),
					dspOccupationThreshold,
                    ((maxDSP<0)?(unsigned)INT_MAX:(unsigned)maxDSP)
			);
		}
		else if(tilingMethod.compare("heuristicGreedyTiling") == 0) {
			tilingStrategy = new TilingStrategyGreedy(
					wX,
					wY,
					wOut,
					signedIO,
					&baseMultiplierCollection,
					baseMultiplierCollection.getPreferedMultiplier(),
					dspOccupationThreshold,
					maxDSP,
					useirregular,
					use2xk,
					superTiles,
					useKaratsuba,
					multiplierTileCollection,
                    guardBits,
                    keepBits
			);
		}
		else if(tilingMethod.compare("heuristicXGreedyTiling") == 0) {
			tilingStrategy = new TilingStrategyXGreedy(
					wX,
					wY,
					wOut,
					signedIO,
					&baseMultiplierCollection,
					baseMultiplierCollection.getPreferedMultiplier(),
					dspOccupationThreshold,
					maxDSP,
					useirregular,
					use2xk,
					superTiles,
					useKaratsuba,
					multiplierTileCollection,
                    guardBits,
                    keepBits
			);
		}
		else if(tilingMethod.compare("heuristicBeamSearchTiling") == 0) {
			tilingStrategy = new TilingStrategyBeamSearch(
					wX,
					wY,
					wOut,
					signedIO,
					&baseMultiplierCollection,
					baseMultiplierCollection.getPreferedMultiplier(),
					dspOccupationThreshold,
					maxDSP,
					useirregular,
					use2xk,
					superTiles,
					useKaratsuba,
					multiplierTileCollection,
					beamRange,
                    guardBits,
                    keepBits
			);
		} else if(tilingMethod.compare("optimal") == 0){
			tilingStrategy = new TilingStrategyOptimalILP(
					wX,
					wY,
					wOut,
					signedIO,
					&baseMultiplierCollection,
					baseMultiplierCollection.getPreferedMultiplier(),
					dspOccupationThreshold,
					maxDSP,
					multiplierTileCollection,
                    guardBits,
                    keepBits,
                    errorBudget,
                    centerErrConstant,
                    optiTrunc
			);

		}  else if(tilingMethod.compare("optimalTilingAndCompression") == 0){
			tilingStrategy = new TilingAndCompressionOptILP(
					wX,
					wY,
					wOut,
					signedIO,
					&baseMultiplierCollection,
					baseMultiplierCollection.getPreferedMultiplier(),
					dspOccupationThreshold,
					maxDSP,
					multiplierTileCollection,
					&bitHeap,
                    guardBits,
                    keepBits,
                    errorBudget,
                    centerErrConstant,
                    optiTrunc,
                    minStages
			);

		} else {
			THROWERROR("Tiling strategy " << tilingMethod << " unknown");
		}

		REPORT(DEBUG, "Solving tiling problem")
		tilingStrategy->solve();

		list<TilingStrategy::mult_tile_t> &solution = tilingStrategy->getSolution();
		if(signedIO)
            for(auto & tile : solution) {
                //resize DSPs to be aligned with left and bottom border of the tiled area to allow the correct handling of the sign
                tile.first = tile.first.shrinkFitDSP(tile.second.first,tile.second.second, wX, wY);
                //Set signedness of individual tiles according to their position
                tile.first = tile.first.setSignStatus(tile.second.first,tile.second.second, wX, wY, signedIO);
            }
        tilingStrategy->printSolution();
		auto solLen = solution.size();
		REPORT(DETAILED, "Found solution has " << solLen << " tiles")
		if (target_->generateFigures()) {
			ofstream texfile;
			texfile.open("multiplier_tiling.tex");
			if ((texfile.rdstate() & ofstream::failbit) != 0) {
				cerr << "Error when opening multiplier_tiling.tex file for output. Will not print tiling configuration."
					 << endl;
			} else {
				tilingStrategy->printSolutionTeX(texfile, wOut, false);
				texfile.close();
			}

			texfile.open("multiplier_tiling.svg");
			if ((texfile.rdstate() & ofstream::failbit) != 0) {
				cerr << "Error when opening multiplier_tiling.svg file for output. Will not print tiling configuration."
					 << endl;
			} else {
				tilingStrategy->printSolutionSVG(texfile, wOut, false);
				texfile.close();
			}

			texfile.open("multiplier_shape.tex");
			if ((texfile.rdstate() & ofstream::failbit) != 0) {
				cerr << "Error when opening multiplier_shape.tex file for output. Will not print tiling configuration."
					 << endl;
			} else {
				tilingStrategy->printSolutionTeX(texfile, wOut, true);
				texfile.close();
			}
		}

		schedule();

       	if (guardBits > 0) {
            //Check truncated solution
            mpz_class actualTruncError = checkTruncationError(solution, guardBits, errorBudget, centerErrConstant);
            cout << "calc min req weight is=" << prodsize(wX, wY, signedIO, signedIO) - (wOut + guardBits) << endl;
            bitHeapLSBWeight = (dynamic_cast<CompressionStrategy*>(tilingStrategy))?0:calcBitHeapLSB(solution, guardBits, errorBudget, centerErrConstant, actualTruncError);
            guardBits = wFullP - wOut - bitHeapLSBWeight; //To select result bits, because the dynamic ilp does not consider guardBits

		    //this is the rounding bit for a faithfully rounded truncated multiplier
			bitHeap.addConstantOneBit(static_cast<int>(guardBits) - 1);
			//these are the constant bits to recenter the average error around 0 and allow for more truncation error
            mpz_class colweight, bitstate;
            int i = wFullP - wOut - guardBits ;
            do{
                mpz_pow_ui(colweight.get_mpz_t(), mpz_class(2).get_mpz_t(), i);
                mpz_and(bitstate.get_mpz_t(), colweight.get_mpz_t(), centerErrConstant.get_mpz_t());
                if (bitstate) {
                    bitHeap.addConstantOneBit(i - (wFullP - wOut - guardBits));
                    REPORT(DEBUG,  "Adding constant bit with weight=" << i << " BitHeap col=" << i - (wFullP - wOut - guardBits) << "to recenter the truncation error at 0");
                    cout << "height at pos " << i - (wFullP - wOut - guardBits) << ": " << bitHeap.getColumnHeight( i - (wFullP - wOut - guardBits)) << endl;
                }
                i++;
            } while(colweight <= centerErrConstant);
		}

        branchToBitheap(&bitHeap, solution, bitHeapLSBWeight);

		if (dynamic_cast<CompressionStrategy*>(tilingStrategy)) {
			std::cout << "Class is derived from CompressionStrategy, passing result for compressor tree.\n";
			bitHeap.startCompression(dynamic_cast<CompressionStrategy*>(tilingStrategy));
		} else {
			bitHeap.startCompression();
		}

		vhdl << tab << "R" << " <= " << bitHeap.getSumName() << range(wOut-1 + guardBits, guardBits) << ";" << endl;

		delete tilingStrategy;
	}

	unsigned int IntMultiplier::computeGuardBits(unsigned int wX, unsigned int wY, unsigned int wOut)
	{
		unsigned int minW = (wX < wY) ? wX : wY;
		unsigned int maxW = (wX < wY) ? wY : wX;

		auto ps = prodsize(wX, wY, signedIO, signedIO);
		auto nbDontCare = ps - wOut;

		mpz_class errorBudget{1};
		errorBudget <<= (nbDontCare >= 1) ? nbDontCare - 1 : 0;
		errorBudget -= 1;

		REPORT(DEBUG, "computeGuardBits: error budget is " << errorBudget.get_str())

		unsigned int nbUnneeded;
		for (nbUnneeded = 1 ; nbUnneeded <= nbDontCare ; nbUnneeded += 1) {
			unsigned int bitAtCurCol;
			if (nbUnneeded < minW) {
				bitAtCurCol = nbUnneeded;
			} else if (nbUnneeded <= maxW) {
				bitAtCurCol = minW;
			} else {
				bitAtCurCol = ps -nbUnneeded;
			}
			REPORT(DETAILED, "computeGuardBits: Nb bit in column " << nbUnneeded << " : " << bitAtCurCol)
			mpz_class currbitErrorAmount{bitAtCurCol};
			currbitErrorAmount <<= (nbUnneeded - 1);

			REPORT(DETAILED, "computeGuardBits: Local error for column " << nbUnneeded << " : " << currbitErrorAmount.get_str())

			errorBudget -= currbitErrorAmount;
			REPORT(DETAILED, "computeGuardBits: New error budget: " << errorBudget.get_str())
			if(errorBudget < 0) {
				break;
			}
		}
		nbUnneeded -= 1;

		return nbDontCare - nbUnneeded;
	}

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
    void IntMultiplier::computeTruncMultParams(unsigned wFull, unsigned wOut, unsigned &g, unsigned &k, unsigned long long &errorBudget, unsigned long long &constant){
        // first loop iterates over the columns, right to left
        unsigned w = wFull - wOut; //weight of the LSB of the result, relative to the LSB of a non-truncated multiplier
        errorBudget = (1ULL<<(w-1)); //tiling error budget
        if(w == 0) return;
        unsigned col = 0, height;
        unsigned long long weightedSumOfTruncatedBits = 0;   //actual error
        bool loop=true;
        while(loop){
            col++;                                                          //bitheap column
            height = (col > (wFull/2))?wFull-col:col;                       //number of partial products in column
            weightedSumOfTruncatedBits += height * (1ULL<<(col-1));
            constant = (1ULL<<(w-1)) - (1ULL<<(col-1));
            //cout << "col=" << col << " height=" << height << " wstb=" << weightedSumOfTruncatedBits << " errorBudget=" << errorBudget << " C=" << constant << endl;
            loop = (weightedSumOfTruncatedBits < errorBudget + constant);
        } // when we exit the loop, we have found g
        g = w-(col-1);
        // Now add back bits in rigthtmost column, one by one
        k = 0;
        while(weightedSumOfTruncatedBits >= errorBudget + constant) {
            weightedSumOfTruncatedBits -= (1ULL<<(col-1));
            k++;
        }
    }

    void IntMultiplier::computeTruncMultParamsMPZ(unsigned wFull, unsigned wOut, unsigned &g, unsigned &k, mpz_class &errorBudget, mpz_class &constant){
        // first loop iterates over the columns, right to left
        unsigned w = wFull - wOut; //weight of the LSB of the result, relative to the LSB of a non-truncated multiplier
        mpz_class colweight;//, max64;
        mpz_pow_ui(errorBudget.get_mpz_t(), mpz_class(2).get_mpz_t(), w-1);
        if(w == 0) return;
        unsigned col = 0, height;
        mpz_class weightedSumOfTruncatedBits = 0;   //actual error
        bool loop=true;
        while(loop){
            col++;                                                          //bitheap column
            //height = (col > (wFull/2))?wFull-col:col;                       //number of partial products in column (square multiplier)
            height = widthOfDiagonalOfRect(wX, wY, col, wFull);             //number of partial products in column
            mpz_pow_ui(colweight.get_mpz_t(), mpz_class(2).get_mpz_t(), col-1);
            weightedSumOfTruncatedBits += mpz_class(height) * colweight;
            constant = errorBudget - colweight;
            //cout << "col=" << col << " height=" << height << " wstb=" << weightedSumOfTruncatedBits << " errorBudget=" << errorBudget << " C=" << C << endl;
            loop = (weightedSumOfTruncatedBits < errorBudget + constant);
        } // when we exit the loop, we have found g
        g = w-(col-1);
        // Now add back bits in rigthtmost column, one by one
        k = 0;
        while(weightedSumOfTruncatedBits >= errorBudget + constant) {
            weightedSumOfTruncatedBits -= colweight;
            k++;
        }

/*        unsigned long long max64u = UINT64_MAX;
        mpz_import(max64.get_mpz_t(), 1, -1, sizeof max64u, 0, 0, &max64u);
        if(errorBudget <= max64 && constant <= max64){
            mpz_export(&errorBudget, 0, -1, sizeof errorBudget, 0, 0, errorBudget.get_mpz_t());
            mpz_export(&constant, 0, -1, sizeof constant, 0, 0, errorBudget.get_mpz_t());
        } else {
            cout << "WARNING: errorBudget or constant exceeds the number range of uint64, use optiTrunc=0 or results will be faulty!" << endl;
        }*/

    }

    unsigned int IntMultiplier::widthOfDiagonalOfRect(unsigned wX, unsigned wY, unsigned col, unsigned wFull){
        if(0 <= col && col <= min((int)wX, (int)wY)){
            return col;
        } else if(min((int)wX,(int)wY) < col && col <= max((int)wX, (int)wY)){
            return min((int)wX,(int)wY);
        } else {
            return wFull-col;
        }
    }

    /**
     * @brief Compute the size required to store the untruncated product of inputs of a given width
     * @param wX size of the first input
     * @param wY size of the second input
     * @return the number of bits needed to store a product of I<wX> * I<WY>
     */
    unsigned int IntMultiplier::prodsize(unsigned int wX, unsigned int wY, bool signedX, bool signedY)
    {
        if(wX == 0 || wY == 0)
            return 0;

        if (wX == 1 && wY == 1) {
            return 1;
        }

        if(wX == 1 && !signedX)
            return wY;

        if(wY == 1 && !signedY)
            return wX;

        return wX + wY;
    }

	/**
	 * @brief getZeroMSB if the tile goes out of the multiplier,
	 *		  compute the number of extra zeros that should be inputed to the submultiplier as MSB
	 * @param anchor lsb of the tile
	 * @param tileWidth width of the tile
	 * @param multiplierLimit msb of the overall multiplier
	 * @return
	 */
	inline unsigned int getZeroMSB(int anchor, unsigned int tileWidth, int multiplierLimit)
	{
		unsigned int result = 0;
		if (static_cast<int>(tileWidth) + anchor > multiplierLimit) {
			result = static_cast<unsigned int>(
						static_cast<int>(tileWidth) + anchor - multiplierLimit
						);
		}
		return result;
	}

	/**
	 * @brief getInputSignal Perform the selection of the useful bit and pad if necessary the input signal to fit a certain format
	 * @param msbZeros number of zeros to add left of the selection
	 * @param selectSize number of bits to get from the inputSignal
	 * @param offset where to start taking bits from the input signal
	 * @param lsbZeros number of zeros to append to the selection for right padding
	 * @param signalName input signal name from which to select the bits
	 * @return the string corresponding to the signal selection and padding
	 */
	inline string getInputSignal(unsigned int msbZeros, unsigned int selectSize, unsigned int offset, unsigned int lsbZeros, string const & signalName)
	{
		stringstream s;
		if (msbZeros > 0)
			s << zg(static_cast<int>(msbZeros)) << " & ";

		s << signalName << range(
				 static_cast<int>(selectSize - 1 + offset),
				 static_cast<int>(offset)
			);

		if (lsbZeros > 0)
			s << " & " << zg(static_cast<int>(lsbZeros));

		return s.str();
	}

	void IntMultiplier::branchToBitheap(BitHeap* bitheap, list<TilingStrategy::mult_tile_t> &solution, unsigned int bitheapLSBWeight)
	{
		size_t i = 0;
		stringstream oname, ofname;
		for(auto & tile : solution) {
			auto& parameters = tile.first;
			auto& anchor = tile.second;

			int xPos = anchor.first;
			int yPos = anchor.second;

			int LSBWeight = xPos + yPos + parameters.getRelativeResultLSBWeight();
			unsigned int outLSBWeight = (LSBWeight < 0)?0:static_cast<unsigned int>(LSBWeight);                         // calc result LSB weight corresponding to tile position
			unsigned int truncated = (outLSBWeight < bitheapLSBWeight) ? bitheapLSBWeight - outLSBWeight : 0;           // calc result LSBs to be ignored
			unsigned int bitHeapOffset = (outLSBWeight < bitheapLSBWeight) ? 0 : outLSBWeight - bitheapLSBWeight;       // calc bits between the tiles output LSB and the bitheaps LSB

			//unsigned int toSkip = lsbZerosXIn + lsbZerosYIn + truncated;                                                // calc LSB bits to be ignored in the tiles output
			unsigned int toSkip = ((LSBWeight < 0) ? static_cast<unsigned int>(-LSBWeight) : 0) + truncated;            // calc LSB bits to be ignored in the tiles output
			unsigned int tokeep = parameters.getRelativeResultMSBWeight() - toSkip - parameters.getRelativeResultLSBWeight()+1;                                     // the tiles MSBs that are actually used
			assert(tokeep > 0); //A tiling should not give a useless tile

			oname.str("");
			oname << "tile_" << i << "_output";
			realiseTile(tile, i, oname.str());

			ofname.str("");
			ofname << "tile_" << i << "_filtered_output";

			if(parameters.getOutputWeights().size()){
			    //The multiplier tile has more then one output signal
				for(unsigned i = 0; i < parameters.getOutputWeights().size(); i++){
					if(i){
						vhdl << declare(.0, ofname.str() + to_string(i), 41) << " <= " << "" << oname.str() + to_string(i) << "(40 downto 0)" << ";" << endl;
						getSignalByName(ofname.str() + to_string(i))->setIsSigned();
						bitheap->addSignal(ofname.str() + to_string(i), bitHeapOffset+parameters.getOutputWeights()[i]);
					} else {
						vhdl << declare(.0, ofname.str(), 41) << " <= " << "" << oname.str() << "(40 downto 0)" << ";" << endl;
						bitheap->addSignal(ofname.str(), bitHeapOffset+parameters.getOutputWeights()[i]);
					}
					cout << "output (" << i+1 << "/" << parameters.getOutputWeights().size() << "): " << ofname.str() + to_string(i) << " shift " << bitHeapOffset+parameters.getOutputWeights()[i] << endl;
				}
			} else {
                //The multiplier tile has only a single output signal
				unsigned int xInputLength = parameters.getTileXWordSize();
				unsigned int yInputLength = parameters.getTileYWordSize();
				bool xIsSigned = parameters.isSignedMultX();
				bool yIsSigned = parameters.isSignedMultY();

				bool bothOne = (xInputLength == 1) && (yInputLength == 1);
				bool signedCase = (bothOne and (xIsSigned != yIsSigned)) or ((not bothOne) and (xIsSigned or yIsSigned));

				vhdl << tab << declareFixPoint(.0, ofname.str(), signedCase, tokeep-1, 0) << " <= " << ((signedCase) ? "" : "un") << "signed(" << oname.str() <<
						range(toSkip + tokeep - 1, toSkip) << ");" << endl;


				bitheap->addSignal(ofname.str(), bitHeapOffset);

			}
			i += 1;
		}
	}

	Operator* IntMultiplier::realiseTile(TilingStrategy::mult_tile_t & tile, size_t idx, string output_name)
	{
		auto& parameters = tile.first;
		auto& anchor = tile.second;
		int xPos = anchor.first;
		int yPos = anchor.second;

		unsigned int xInputLength = parameters.getTileXWordSize();
		unsigned int yInputLength = parameters.getTileYWordSize();
		//unsigned int outputLength = parameters.getOutWordSize();

		unsigned int lsbZerosXIn = (xPos < 0) ? static_cast<unsigned int>(-xPos) : 0;
		unsigned int lsbZerosYIn = (yPos < 0) ? static_cast<unsigned int>(-yPos) : 0;

		unsigned int msbZerosXIN = getZeroMSB(xPos, xInputLength, wX);
		unsigned int msbZerosYIn = getZeroMSB(yPos, yInputLength, wY);

		unsigned int selectSizeX = xInputLength - (msbZerosXIN + lsbZerosXIn);
		unsigned int selectSizeY = yInputLength - (msbZerosYIn + lsbZerosYIn);

		unsigned int effectiveAnchorX = (xPos < 0) ? 0 : static_cast<unsigned int>(xPos);
		unsigned int effectiveAnchorY = (yPos < 0) ? 0 : static_cast<unsigned int>(yPos);

		auto tileXSig = getInputSignal(msbZerosXIN, selectSizeX, effectiveAnchorX, lsbZerosXIn, "X");
		auto tileYSig = getInputSignal(msbZerosYIn, selectSizeY, effectiveAnchorY, lsbZerosYIn, "Y");

		stringstream nameX, nameY, nameOutput;
		nameX << "tile_" << idx << "_X";
		nameY << "tile_" << idx << "_Y";

		nameOutput << "tile_" << idx << "_Out";

		vhdl << tab << declare(0., nameX.str(), xInputLength) << " <= " << tileXSig << ";" << endl;
		vhdl << tab << declare(0., nameY.str(), yInputLength) << " <= " << tileYSig << ";" << endl;

		string multIn1SigName = nameX.str();
		string multIn2SigName = nameY.str();

		if (parameters.isFlippedXY())
			std::swap(multIn1SigName, multIn2SigName);

		nameOutput.str("");
		nameOutput << "tile_" << idx << "_mult";

		inPortMap("X", multIn1SigName);
		inPortMap("Y", multIn2SigName);
		outPortMap("R", output_name);
		for(unsigned i = 1; i < parameters.getOutputWeights().size(); i++){
			outPortMap("R" + to_string(i), output_name + to_string(i));
		}

		auto mult = parameters.generateOperator(this, getTarget());

		vhdl << instance(mult, nameOutput.str(), false) <<endl;
		return mult;
	}

	string IntMultiplier::addUID(string name, int blockUID)
	{
		ostringstream s;

		s << name << "_m" << multiplierUid;
		if(blockUID != -1)
			s << "b" << blockUID;
		return s.str() ;
	}

	unsigned int IntMultiplier::getOutputLengthNonZeros(
			BaseMultiplierParametrization const & parameter,
			unsigned int xPos,
			unsigned int yPos,
			unsigned int totalOffset
		)
	{
//        unsigned long long int sum = 0;
		mpfr_t sum;
		mpfr_t two;

		mpfr_inits(sum, two, NULL);

		mpfr_set_si(sum, 1, GMP_RNDN); //init sum to 1 as at the end the log(sum+1) has to be computed
		mpfr_set_si(two, 2, GMP_RNDN);

/*
	cout << "sum=";
	mpfr_out_str(stdout, 10, 0, sum, MPFR_RNDD);
	cout << endl;
	flush(cout);
*/

		for(unsigned int i = 0; i < parameter.getTileXWordSize(); i++){
			for(unsigned int j = 0; j < parameter.getTileYWordSize(); j++){
				if(i + xPos >= totalOffset && j + yPos >= totalOffset){
					if(i + xPos < (wX + totalOffset) && j + yPos < (wY + totalOffset)){
						if(parameter.shapeValid(i,j)){
							//sum += one << (i + j);
							mpfr_t r;
							mpfr_t e;
							mpfr_inits(r, e, NULL);
							mpfr_set_si(e, i+j, GMP_RNDD);
							mpfr_pow(r, two, e, GMP_RNDD);
							mpfr_add(sum, sum, r, GMP_RNDD);
//                            cout << " added 2^" << (i+j) << endl;
						}
						else{
							//cout << "not true " << i << " " << j << endl;
						}
					}
					else{
						//cout << "upper borders " << endl;
					}
				}
				else{
					//cout << "lower borders" << endl;
				}
			}
		}

		mpfr_t length;
		mpfr_inits(length, NULL);
		mpfr_log2(length, sum, GMP_RNDU);
		mpfr_ceil(length,length);

		unsigned long length_ul = mpfr_get_ui(length, GMP_RNDD);
//        cout << " outputLength has a length of " << length_ul << endl;
		return (int) length_ul;
	}

	/*this function computes the amount of zeros at the lsb position (mode 0). This is done by checking for every output bit position, if there is a) an input of the current BaseMultiplier, and b) an input of the IntMultiplier. if a or b or both are false, then there is a zero at this position and we check the next position. otherwise we break. If mode==1, only condition a) is checked */
	/*if mode is 2, only the offset for the bitHeap is computed. this is done by counting how many diagonals have a position, where is an IntMultiplierInput but not a BaseMultiplier input. */
	unsigned int IntMultiplier::getLSBZeros(
			BaseMultiplierParametrization const & param,
			unsigned int xPos,
			unsigned int yPos,
			unsigned int totalOffset,
			int mode
		){
		//cout << "mode is " << mode << endl;
		unsigned int max = min(int(param.getMultXWordSize()), int(param.getMultYWordSize()));
		unsigned int zeros = 0;
		unsigned int mode2Zeros = 0;
		bool startCountingMode2 = false;

		for(unsigned int i = 0; i < max; i++){
			unsigned int steps = i + 1; //for the lowest bit make one step (pos 0,0), second lowest 2 steps (pos 0,1 and 1,0) ...
			//the relevant positions for every outputbit lie on a diogonal line.


			for(unsigned int j = 0; j < steps; j++){
				bool bmHasInput = false;
				bool intMultiplierHasBit = false;
				if(param.shapeValid(j, steps - (j + 1))){
					bmHasInput = true;
				}
				if(bmHasInput && mode == 1){
					return zeros;   //only checking condition a)
				}
				//x in range?
				if((xPos + j >= totalOffset && xPos + j < (wX + totalOffset))){
					//y in range?
					if((yPos + (steps - (j + 1))) >= totalOffset && (yPos + (steps - (j + 1))) < (wY + totalOffset)){
						intMultiplierHasBit = true;
					}
				}
				if(mode == 2 && yPos + xPos + (steps - 1) >= 2 * totalOffset){
					startCountingMode2 = true;
				}
				//cout << "position " << j << "," << (steps - (j + 1)) << " has " << bmHasInput << " and " << intMultiplierHasBit << endl;
				if(bmHasInput && intMultiplierHasBit){
					//this output gets some values. So finished computation and return
					if(mode == 0){
						return zeros;
					}
					else{	//doesnt matter if startCountingMode2 is true or not. mode2Zeros are being incremented at the end of the diagonal
						return mode2Zeros;
					}

				}
			}

			zeros++;
			if(startCountingMode2 == true){
				mode2Zeros++;
			}
		}
		if(mode != 2){
			return zeros;       //if reached, that would mean the the bm shares no area with IntMultiplier
		}
		else{
			return mode2Zeros;
		}
	}

	void IntMultiplier::emulate (TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		mpz_class svY = tc->getInputValue("Y");
		mpz_class svR;

//		cerr << "X : " << svX.get_str(10) << " (" << svX.get_str(2) << ")" << endl;
//		cerr << "Y : " << svY.get_str(10) << " (" << svY.get_str(2) << ")" << endl;


		if(signedIO)
		{
			// Manage signed digits
			mpz_class big1 = (mpz_class{1} << (wX));
			mpz_class big1P = (mpz_class{1} << (wX-1));
			mpz_class big2 = (mpz_class{1} << (wY));
			mpz_class big2P = (mpz_class{1} << (wY-1));

			if(svX >= big1P) {
				svX -= big1;
				//cerr << "X is neg. Interpreted value : " << svX.get_str(10) << endl;
			}

			if(svY >= big2P) {
				svY -= big2;
				//cerr << "Y is neg. Interpreted value : " << svY.get_str(10) << endl;
			}
		}
		svR = svX * svY;
		//cerr << "Computed product : " << svR.get_str() << endl;

		if(negate)
			svR = -svR;

		// manage two's complement at output
		if(svR < 0) {
			svR += (mpz_class(1) << wFullP);
			//cerr << "Prod is neg : encoded value : " << svR.get_str(2) << endl;
		}

		if(wOut>=wFullP) {
			tc->addExpectedOutput("R", svR);
			//cerr << "Exact result is required" << endl;
		}
		else {
			// there is truncation, so this mult should be faithful
			// TODO : this is not faithful but almost faithful :
			// currently only test if error <= ulp, faithful is < ulp
			mpz_class svRtrunc = (svR >> (wFullP-wOut));
			tc->addExpectedOutput("R", svRtrunc);
			if ((svRtrunc << (wFullP - wOut)) != svR) {
				svRtrunc++;
				svRtrunc &= (mpz_class(1) << (wOut)) -1;
				tc->addExpectedOutput("R", svRtrunc);
			}
		}
	}


	void IntMultiplier::buildStandardTestCases(TestCaseList* tcl)
	{
		TestCase *tc;

		mpz_class x, y;

		// 1*1
		x = mpz_class(1);
		y = mpz_class(1);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// -1 * -1
		x = (mpz_class(1) << wX) -1;
		y = (mpz_class(1) << wY) -1;
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);

		// The product of the two max negative values overflows the signed multiplier
		x = mpz_class(1) << (wX -1);
		y = mpz_class(1) << (wY -1);
		tc = new TestCase(this);
		tc->addInput("X", x);
		tc->addInput("Y", y);
		emulate(tc);
		tcl->add(tc);
	}





	OperatorPtr IntMultiplier::parseArguments(OperatorPtr parentOp, Target *target, std::vector<std::string> &args) {
		int wX,wY, wOut, maxDSP;
		bool signedIO,superTile, use2xk, useirregular, useLUT, useDSP, useKaratsuba, optiTrunc, minStages;
		double dspOccupationThreshold=0.0;
		int beamRange = 0;

		UserInterface::parseStrictlyPositiveInt(args, "wX", &wX);
		UserInterface::parseStrictlyPositiveInt(args, "wY", &wY);
		UserInterface::parsePositiveInt(args, "wOut", &wOut);
		UserInterface::parseBoolean(args, "signedIO", &signedIO);
		UserInterface::parseBoolean(args, "superTile", &superTile);
		UserInterface::parseBoolean(args, "use2xk", &use2xk);
		UserInterface::parseBoolean(args, "useirregular", &useirregular);
		UserInterface::parseBoolean(args, "useLUT", &useLUT);
		UserInterface::parseBoolean(args, "useDSP", &useDSP);
		UserInterface::parseBoolean(args, "useKaratsuba", &useKaratsuba);
		UserInterface::parseFloat(args, "dspThreshold", &dspOccupationThreshold);
		UserInterface::parseInt(args, "maxDSP", &maxDSP);
        UserInterface::parseBoolean(args, "optiTrunc", &optiTrunc);
        UserInterface::parseBoolean(args, "minStages", &minStages);
		UserInterface::parsePositiveInt(args, "beamRange", &beamRange);

		return new IntMultiplier(parentOp, target, wX, wY, wOut, signedIO, dspOccupationThreshold, maxDSP, superTile, use2xk, useirregular, useLUT, useDSP, useKaratsuba, beamRange, optiTrunc, minStages);
	}



	void IntMultiplier::registerFactory(){
		UserInterface::add("IntMultiplier", // name
											 "A pipelined integer multiplier.",
											 "BasicInteger", // category
											 "", // see also
											 "wX(int): size of input X; wY(int): size of input Y;\
						wOut(int)=0: size of the output if you want a truncated multiplier. 0 for full multiplier;\
						signedIO(bool)=false: inputs and outputs can be signed or unsigned;\
						maxDSP(int)=-1: limits the number of DSP-Tiles used in Multiplier;\
						use2xk(bool)=false: if true, attempts to use the 2xk-LUT-Multiplier with relatively high efficiency;\
						useirregular(bool)=false: if true, attempts to use the irregular-LUT-Multipliers with higher area/lut efficiency than the rectangular versions;\
						useLUT(bool)=true: if true, attempts to use the LUT-Multipliers for tiling;\
						useDSP(bool)=true: if true, attempts to use the DSP-Multipliers for tiling;\
						useKaratsuba(bool)=false: if true, attempts to use rectangular Karatsuba for tiling;\
						superTile(bool)=false: if true, attempts to use the DSP adders to chain sub-multipliers. This may entail lower logic consumption, but higher latency.;\
						dspThreshold(real)=0.0: threshold of relative occupation ratio of a DSP multiplier to be used or not;\
                        optiTrunc(bool)=true: if true, considers the Truncation error dynamicly, instead of defining a hard border for tiling, like in th ARITH paper;\
                        minStages(bool)=true: if true, minimizes stages in combined opt. of tiling an comp., otherwise try to find a sol. with less LUTs and more stages;\
						beamRange(int)=3: range for beam search", // This string will be parsed
											 "", // no particular extra doc needed
											IntMultiplier::parseArguments,
											IntMultiplier::unitTest
											 ) ;
	}


	TestList IntMultiplier::unitTest(int )
	{
		// the static list of mandatory tests
		TestList testStateList;
		vector<pair<string,string>> paramList;

		list<pair<int,int>> wordSizes = {{1,1},{2,2},{3,3},{4,4},{8,8},{16,8},{8,16},{13,17},{24,24},{27,41},{35,35},{53,53},{64,64},{10,99},{99,10},{100,100}};

		for (auto wordSizePair : wordSizes)
			{
				int wX = wordSizePair.first;
				int wY = wordSizePair.second;
				for(int sign=0; sign < 2; sign++)
					{
						paramList.push_back(make_pair("wX", to_string(wX)));
						paramList.push_back(make_pair("wY", to_string(wY)));
						paramList.push_back(make_pair("signedIO", sign ? "true" : "false"));
#if 0		// commented by Florent: better leave the default if you don't actually test various values

						paramList.push_back(make_pair("dspThreshold", to_string(1.0)));
						paramList.push_back(make_pair("maxDSP", to_string(0)));
#endif

						testStateList.push_back(paramList);

						// same parameters, but truncate (to the size of the smallest input for no particular reason)
						paramList.push_back(make_pair("wOut", to_string(wX)));
						testStateList.push_back(paramList);

						paramList.clear();
					}
			}

		return testStateList;
	}

    mpz_class IntMultiplier::checkTruncationError(list<TilingStrategy::mult_tile_t> &solution, unsigned guardBits, const mpz_class& errorBudget, const mpz_class& constant) const{
        std::vector<std::vector<bool>> mulArea(wX, std::vector<bool>(wY,false));

        for(auto & tile : solution) {
            auto &parameters = tile.first;
            auto &anchor = tile.second;
            int xPos = anchor.first;
            int yPos = anchor.second;

            for(int x = 0; x < parameters.getMultXWordSize(); x++){
                for(int y = 0; y < parameters.getMultYWordSize(); y++){
                    if(0 <= xPos+x && 0 <= yPos+y && xPos+x < wX && yPos+y < wY)
                        mulArea[xPos+x][yPos+y] = mulArea[xPos+x][yPos+y] || parameters.shapeValid(x,y);
                }
            }
        }

        mpz_class truncError, maxErr = errorBudget+constant;
        truncError = mpz_class(0);
        for(int y = 0; y < wY; y++){
            for(int x = wX-1; 0 <= x; x--){
                if(!mulArea[x][y])
                    truncError += (mpz_class(1)<<(x+y));
                cout << (mulArea[x][y] ? 1 : 0);
            }
            cout << endl;
        }

        if(truncError <= maxErr){
            cout << "OK: actual truncation error=" << truncError << " is smaller than the max. permissible error=" << maxErr << " by " << maxErr-truncError << "." << endl;
        } else {
            cout << "ERROR: actual truncation error=" << truncError << " is larger than the max. permissible error=" << maxErr << " by " << truncError-maxErr << "." << endl;
        }
        return truncError;
	}

    int IntMultiplier::calcBitHeapLSB(list<TilingStrategy::mult_tile_t> &solution, unsigned guardBits, const mpz_class& errorBudget, const mpz_class& constant, const mpz_class& actualTruncError){
        int lsbs[solution.size()], msbs[solution.size()], prunableBits[solution.size()], weight=0, nBits = 0, i=0;;
        for(auto & tile : solution) {
            auto &parameters = tile.first;
            auto &anchor = tile.second;
            int xPos = anchor.first;
            int yPos = anchor.second;
            lsbs[i] = xPos + yPos + parameters.getRelativeResultLSBWeight();
            msbs[i] = lsbs[i] + parameters.getOutWordSize();
            prunableBits[i]=0;
            cout << "Tile " << i << " " << parameters.getMultType() << " at (" << xPos << "," << yPos << ") has an LSB of" << lsbs[i] << endl;
            i++;
        }
        double error = 0;
        do{
            cout << " min weight=" << weight << endl;
            for(i = 0; i < solution.size(); i++){
                if(lsbs[i] <= weight && weight < msbs[i]) {
                    prunableBits[i]++;
                    if(weight < lsbs[i]+(msbs[i]-lsbs[i])/2 ){
                        error += prunableBits[i]*pow(2,weight);
                    } else {
                        error += ((msbs[i]-lsbs[i])-prunableBits[i])*pow(2,weight);
                    }
                    cout << "trying to prune " << prunableBits[i] << " bits in tile " << i << " with weight " << weight << " error is "  << error << " additional permissible error is "  << errorBudget + constant - actualTruncError << endl;
                }
            }
            weight++;
        } while(actualTruncError + error < errorBudget + constant || 0 == error);
        weight--;
        cout << "min req weight is=" << weight << endl;
        return weight;
    }

}
