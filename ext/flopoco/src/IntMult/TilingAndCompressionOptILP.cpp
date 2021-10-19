#include "IntMult/TilingAndCompressionOptILP.hpp"
#include "IntMultiplier.hpp"
#include "BaseMultiplierLUT.hpp"
#include "MultiplierTileCollection.hpp"

using namespace std;
namespace flopoco {

TilingAndCompressionOptILP::TilingAndCompressionOptILP(
		unsigned int wX_,
		unsigned int wY_,
		unsigned int wOut_,
		bool signedIO_,
		BaseMultiplierCollection* bmc,
		base_multiplier_id_t prefered_multiplier,
		float occupation_threshold,
		int maxPrefMult,
        MultiplierTileCollection mtc_,
        BitHeap* bitheap,
        unsigned guardBits,
        unsigned keepBits,
        mpz_class errorBudget,
        mpz_class &centerErrConstant,
        bool performOptimalTruncation,
        bool minStages):TilingStrategy(
			wX_,
			wY_,
			wOut_,
			signedIO_,
			bmc),
        CompressionStrategy(bitheap),
		max_pref_mult_ {maxPrefMult},
		occupation_threshold_{occupation_threshold},
		tiles{mtc_.MultTileCollection},
        guardBits{guardBits},
        keepBits{keepBits},
        eBudget{errorBudget},
        centerErrConstant{centerErrConstant},
        performOptimalTruncation{performOptimalTruncation},
        minLutsFocus{!minStages}
	{
	    //Check if error bound fits into UINT64 and otherwise use optiTrunc=0
        cout << errorBudget << endl;
        mpz_class max64;
        unsigned long long max64u = (1ULL << 52)-1ULL;//UINT64_MAX; //Limit to dynamic of double type
        mpz_import(max64.get_mpz_t(), 1, -1, sizeof max64u, 0, 0, &max64u);
        if(errorBudget <= max64){
            mpz_export(&this->errorBudget, 0, -1, sizeof this->errorBudget, 0, 0, errorBudget.get_mpz_t());
        } else {
            if(performOptimalTruncation)
                cout << "WARNING: errorBudget or constant exceeds the number range of uint64, switching to optiTrunc=0" << endl;
            this->errorBudget = 0;
            this->performOptimalTruncation = false;
        }
        cout << this->errorBudget << endl;

        cout << "guardBits " << guardBits << " keepBits " << keepBits << endl;

        for(auto &p:tiles)
        {
                cout << p->getLUTCost(0, 0, wX, wY, false) << " " << p->getType() << endl;
        }
 	}

void TilingAndCompressionOptILP::solve()
{

#ifndef HAVE_SCALP
    throw "Error, TilingAndCompressionOptILP::solve() was called but FloPoCo was not built with ScaLP library";
#else
    cout << "using ILP solver " << target->getILPSolver() << endl;
    solver = new ScaLP::Solver(ScaLP::newSolverDynamic({target->getILPSolver(),"Gurobi","CPLEX","SCIP","LPSolve"}));
    solver->timeout = target->getILPTimeout();

    //for debugging it might be better to order the compressors by efficiency
    orderCompressorsByCompressionEfficiency();

    cout << "available compressors" << endl;
    for(unsigned int i = 0; i < possibleCompressors.size(); i++){
        cout << possibleCompressors[i]->getStringOfIO() << " efficiency " << possibleCompressors[i]->getEfficiency() << endl;

    }

    ScaLP::status stat;
    solver->reset();
    bool tryOneMoreStage = minLutsFocus;
    s_max = 0;
    do{ //increase stages till a feasible solution is found to get the minimal possible number of compressor stages
        //Try to find a cheaper result with one more stage, check if previous solution was feasible
        if((stat == ScaLP::status::OPTIMAL || stat == ScaLP::status::FEASIBLE || stat == ScaLP::status::TIMEOUT_FEASIBLE) && tryOneMoreStage){
            cerr << "The result with " << s_max << " stages was " << stat << " but attempting to get a cheaper result with one stage more." << endl;
            bestResult = solver->getResult();
            tryOneMoreStage = false;
        }

        s_max++;
        solver->reset();
        TilingStrategy::solution.clear();
        constructProblem(s_max);

        // Try to solve
        cout << "starting solver, this might take a while..." << endl;
        solver->quiet = false;
        stat = solver->solve();

        // print results
        cerr << "The result is " << stat << endl;
        //cerr << solver->getResult() << endl;
        writeSolutionFile(solver->getResult());

    } while(stat == ScaLP::status::INFEASIBLE || tryOneMoreStage);

    //Select the least expensive solution, when a solution with an additional compression stage was tested
    ScaLP::Result res = solver->getResult();
    if(!bestResult.empty()){
        if(bestResult.objectiveValue <= res.objectiveValue){
            cout << "The solution with " << s_max-1 << " stages is less or equally expensive (" << bestResult.objectiveValue << ") then with " << s_max << " stages (" << res.objectiveValue << ")" << endl;
            s_max--;
            writeSolutionFile(bestResult);
        } else {
            cout << "The solution with " << s_max-1 << " stages is more expensive (" << bestResult.objectiveValue << ") then with " << s_max << " stages (" << res.objectiveValue << ")" << endl;
            bestResult = res;
        }
    } else bestResult = res;

    bool truncErrorBoundsNotMet = false;
    long long optTruncNumericErr = 0;
    do {                                //Loop as a workaround for numeric solver problems during optimal truncation to ensure that the truncation error margin is met
        if(truncErrorBoundsNotMet){
            solver->reset();
            TilingStrategy::solution.clear();
            constructProblem(s_max);

            // Try to solve
            cout << "starting solver, this might take a while..." << endl;
            solver->quiet = false;
            stat = solver->solve();

            // print results
            cerr << "The result is " << stat << endl;
            //cerr << solver->getResult() << endl;
            writeSolutionFile(solver->getResult());
        }
        truncErrorBoundsNotMet = true;  //iterate, if error truncation bounds are not met

        parseTilingSolution(optTruncNumericErr);

    } while(performOptimalTruncation && (wOut < (int)prodWidth) && !checkTruncationError(TilingStrategy::solution, guardBits, eBudget, centerErrConstant) && 0 < optTruncNumericErr);

#endif
}

    void TilingAndCompressionOptILP::compressionAlgorithm() {
#ifndef HAVE_SCALP
		throw "Error, TilingAndCompressionOptILP::compressionAlgorithm() was called but FloPoCo was not built with ScaLP library";
#else
		//adds the Bits to stages and columns
		orderBitsByColumnAndStage();

		//populates bitAmount
		fillBitAmounts();

		//prints out how the inputBits of the bitheap looks like
		printBitAmounts();

		//new solution
		CompressionStrategy::solution = BitHeapSolution();
		CompressionStrategy::solution.setSolutionStatus(BitheapSolutionStatus::OPTIMAL_PARTIAL);

		float compressor_cost = 0;
		vector<vector<int>> zeroInputsVector(s_max, vector<int>((int)prodWidth,0));
		resizeBitAmount(s_max-1);

		for(auto &p:bestResult.values)
		{
			if(p.second > 0.5){     //parametrize all multipliers at a certain position, for which the solver returned 1 as solution, to flopoco solution structure
				std::string var_name = p.first->getName();
				//cout << var_name << "\t " << p.second << endl;
				//if(var_name.substr(0,1).compare("k") != 0) continue;
				switch(var_name.substr(0,1).at(0)) {
					case 'k':{      //decision variables 'k' for placed compressors
						int sta_id = stoi(var_name.substr(2, dpSt));
						int com_id = stoi(var_name.substr(2 + dpSt + 1, dpK));
						int col_id = stoi(var_name.substr(2 + dpSt + 1 + dpK + 1, dpC));
						//cout << "is compressor" << com_id << " stage " << sta_id << " column" << col_id << endl;
						compressor_cost += ((com_id < (int)possibleCompressors.size() && sta_id < s_max - 1)?possibleCompressors[com_id]->area:0) * p.second;
						if (possibleCompressors[com_id] != flipflop && sta_id < s_max - 1) {
							float instances = p.second;
							while (instances-- > 0.5) {
								CompressionStrategy::solution.addCompressor(sta_id, col_id, possibleCompressors[com_id]);
								cout << "is compressor" << com_id << " stage " << sta_id << " column" << col_id
									 << " compressor type " << possibleCompressors[com_id]->getStringOfIO() << endl;
							}
						}
						break;
					}
					case 'Z':{          //compressor input bits that are set zero have to be recorded for VHDL-Generation
						int sta_id = stoi(var_name.substr(2, dpSt));
						int col_id = stoi(var_name.substr(2 + dpSt + 1, dpC));
						zeroInputsVector[sta_id][col_id] = (int)(-1);
						break;
					}
					default:            //Other decision variables are not needed for VHDL-Generation
						break;
				}
				for(int stage = 0; stage < s_max; stage++){
					CompressionStrategy::solution.setEmptyInputsByRemainingBits(stage, zeroInputsVector[stage]);
				}
			}
		}

		cout << "Total compressor LUT-cost: " << compressor_cost << endl;

		//reports the area in LUT-equivalents
		printSolutionStatistics();

		//here the VHDL-Code for the compressors as well as the bits->compressors->bits are being written.
		applyAllCompressorsFromSolution();
#endif
}


#ifdef HAVE_SCALP
void TilingAndCompressionOptILP::constructProblem(int s_max)
{
    cout << "constructing problem formulation with " << s_max << " stages..." << endl;
    wS = tiles.size();


    //Assemble cost function, declare problem variables
    cout << "   assembling cost function, declaring problem variables..." << endl;
    ScaLP::Term obj;
    int x_neg = 0, y_neg = 0, keepBits_ = keepBits;
    for(int s = 0; s < wS; s++){
        x_neg = (x_neg < (int)tiles[s]->wX())?tiles[s]->wX() - 1:x_neg;
        y_neg = (y_neg < (int)tiles[s]->wY())?tiles[s]->wY() - 1:y_neg;
    }
    prodWidth = IntMultiplier::prodsize(wX, wY, signedIO, signedIO);
    int nx = wX-1, ny = wY-1, ns = wS-1; dpX = 1; dpY = 1; dpS = 1;     //calc number of decimal places, for var names
    int nk = possibleCompressors.size()+4, nc = prodWidth + 1, nerc = prodWidth-wOut, nst = s_max; dpK = 1; dpC = 1; dpERC = 1; dpSt = 1;
    nx = (x_neg > nx)?x_neg:nx;                                         //in case the extend in negative direction is larger
    ny = (y_neg > ny)?y_neg:ny;
    while (nx /= 10)
        dpX++;
    while (ny /= 10)
        dpY++;
    while (ns /= 10)
        dpS++;
    while (nk /= 10)
        dpK++;
    while (nc /= 10)
        dpC++;
    while (nc /= 10)
        dpERC++;
    while (nst /= 10)
        dpSt++;

    vector<ScaLP::Term> bitsinColumn(prodWidth + 1);
    vector<vector<vector<ScaLP::Variable>>> solve_Vars(wS, vector<vector<ScaLP::Variable>>(wX+x_neg, vector<ScaLP::Variable>(wY+y_neg)));
    ScaLP::Term maxEpsTerm, constVecTerm;
    __uint64_t sumOfPosEps = 0;
    // add the Constraints
    cout << "   adding the constraints to problem formulation..." << endl;
    for(int y = 0; y < wY; y++){
        for(int x = 0; x < wX; x++){
            stringstream consName;
            consName << "p" << setfill('0') << setw(dpX) << x << setfill('0') << setw(dpY) << y;            //one constraint for every position in the area to be tiled
            ScaLP::Term pxyTerm;
            for(int s = 0; s < wS; s++){					//for every available tile...
                for(int ys = 0 - tiles[s]->wY() + 1; ys <= y; ys++){					//...check if the position x,y gets covered by tile s located at position (xs, ys) = (x-wtile..x, y-htile..y)
                    for(int xs = 0 - tiles[s]->wX() + 1; xs <= x; xs++){
                        if(occupation_threshold_ == 1.0 && ((wX - xs) < (int)tiles[s]->wX() || (wY - ys) < (int)tiles[s]->wY())) break;
                        if(signedIO && (wX < xs+tiles[s]->wX() || wY < ys+tiles[s]->wY()) ) break;                      //Avoid considering protruding cases in signed case, tiling is valid without
                        if(signedIO && (wX == xs+tiles[s]->wX() && !tiles[s]->signSupX() || wY == ys+tiles[s]->wY() && !tiles[s]->signSupY() )) break; //Avoid placing tiles without signed support at the bottom and left |_ edge of the area to be tiled.
                        if(tiles[s]->shape_contribution(x, y, xs, ys, wX, wY, signedIO)){
//                            if((wOut < (int)prodWidth) && ((xs+tiles[s]->wX()+ys+tiles[s]->wY()-2) < ((int)prodWidth-wOut-guardBits))) break;
                            if(tiles[s]->shape_utilisation(xs, ys, wX, wY, signedIO) >=  occupation_threshold_ ){
                                if(solve_Vars[s][xs+x_neg][ys+y_neg] == nullptr){
                                    stringstream nvarName;
                                    nvarName << "d" << setfill('0') << setw(dpS) << s << ((xs < 0)?"m":"") << setfill('0') << setw(dpX) << ((xs<0)?-xs:xs) << ((ys < 0)?"m":"")<< setfill('0') << setw(dpY) << ((ys<0)?-ys:ys) ;
                                    //std::cout << nvarName.str() << endl;
                                    ScaLP::Variable tempV = ScaLP::newBinaryVariable(nvarName.str());
                                    solve_Vars[s][xs+x_neg][ys+y_neg] = tempV;
                                    obj.add(tempV, (double)tiles[s]->ownLUTCost(xs, ys, wX, wY, signedIO));    //append variable to cost function

                                    int col_min = xs+ys+tiles[s]->getRelativeResultLSBWeight(tiles[s]->getParametrisation());
                                    int col_max = xs+ys+tiles[s]->getRelativeResultMSBWeight(tiles[s]->getParametrisation(),signedIO && xs+tiles[s]->wX() == wX, signedIO && ys+tiles[s]->wY() == wY);
                                    for(col_min = ((col_min < 0)?0:col_min); col_min <= (((int)prodWidth<col_max)?(int)prodWidth:col_max); col_min++){
                                        bitsinColumn[col_min].add(tempV, 1);                                      //Output bits from tiles, that have to be compressed on the compressor tree, when a tile is used
                                    }

                                    if(signedIO && (xs+tiles[s]->wX() == wX || ys+tiles[s]->wY() == wY) ){              //Handling of the sign extension bits in dynamically calculated constant bit vector
                                        if(tiles[s]->wX() == 1 && tiles[s]->wY() == 1 && xs == wX-1 && ys == wY-1) break; //the 1x1 tile with (1,1) signedness should not get a sign extension vector.
                                        constVecTerm.add(solve_Vars[s][xs+x_neg][ys+y_neg],(double)(((1ULL<<(wX+wY))-(1ULL<<(xs+ys+tiles[s]->getRelativeResultMSBWeight(tiles[s]->getParametrisation(),signedIO && xs+tiles[s]->wX() == wX, signedIO && ys+tiles[s]->wY() == wY))))));
                                        //printf("s=%02d, x=%02d, y=%02d, wx=%01d, wy=%01d, wOut=%02d, msb=%02d, 0x%08llx\n", s, xs, ys, tiles[s]->wX(), tiles[s]->wY(), tiles[s]->getRelativeResultMSBWeight(tiles[s]->getParametrisation(),signedIO && xs+tiles[s]->wX() == wX, signedIO && ys+tiles[s]->wY() == wY)-tiles[s]->getRelativeResultLSBWeight(tiles[s]->getParametrisation())+1, xs+ys+tiles[s]->getRelativeResultMSBWeight(tiles[s]->getParametrisation(),signedIO && xs+tiles[s]->wX() == wX, signedIO && ys+tiles[s]->wY() == wY), ((1ULL<<(wX+wY))-(1ULL<<(xs+ys+tiles[s]->getRelativeResultMSBWeight(tiles[s]->getParametrisation(),signedIO && xs+tiles[s]->wX() == wX, signedIO && ys+tiles[s]->wY() == wY)))));
                                    }
                                }
                                pxyTerm.add(solve_Vars[s][xs+x_neg][ys+y_neg], 1);
                            }
                        }
                    }
                }
            }

            ScaLP::Constraint c1Constraint;
            if(performOptimalTruncation && (wOut < (int)prodWidth) && ((x + y) < ((int)prodWidth - wOut))) {
                stringstream nvarName;
                nvarName << " b" << setfill('0') << setw(dpX) << x << setfill('0') << setw(dpY) << y;
                ScaLP::Variable tempV = ScaLP::newBinaryVariable(nvarName.str());
                maxEpsTerm.add(tempV, (-1LL) * (1LL << (x + y)));
                maxEpsTerm.add(1ULL << (x + y));

                c1Constraint = pxyTerm - tempV == 0;
            } else if(!performOptimalTruncation && (wOut < (int)prodWidth) && ((x + y) < ((int)prodWidth - wOut - guardBits))){
                //c1Constraint = pxyTerm <= (bool)1;
            } else if(!performOptimalTruncation && (wOut < (int)prodWidth) && ((x + y) == ((int)prodWidth - wOut - guardBits))){
                if((keepBits_)?keepBits_--:0){
                    c1Constraint = pxyTerm == (bool)1;
                    cout << "keepBit at" << x << "," << y << endl;
                } else {
                    cout << "NO keepBit at" << x << "," << y << endl;
                }
            } else {
                c1Constraint = pxyTerm == (bool)1;
            }

            c1Constraint.name = consName.str();
            solver->addConstraint(c1Constraint);
        }
    }

    //limit use of DSPs
    if(0 <= max_pref_mult_) {
        //check if DSP tiles are available
        int nDSPTiles = 0;
        for (int s = 0; s < wS; s++)
            if (tiles[s]->getDSPCost())
                nDSPTiles++;
        if (nDSPTiles) {
            cout << "   adding the constraint to limit the use of DSP-Blocks to " << max_pref_mult_ << " instances..."
                 << endl;
            stringstream consName;
            consName << "limDSP";
            ScaLP::Term pxyTerm;
            for (int y = 0 - 24 + 1; y < wY; y++) {
                for (int x = 0 - 24 + 1; x < wX; x++) {
                    for (int s = 0; s < wS; s++)
                        if (solve_Vars[s][x + x_neg][y + y_neg] != nullptr)
                            for (int c = 0; c < tiles[s]->getDSPCost(); c++)
                                pxyTerm.add(solve_Vars[s][x + x_neg][y + y_neg], 1);
                }
            }
            ScaLP::Constraint c1Constraint = pxyTerm <= max_pref_mult_;     //set max usage equ.
            c1Constraint.name = consName.str();
            solver->addConstraint(c1Constraint);
        }
    }

    //make shure the available precision is present in case of truncation
    if(performOptimalTruncation && (wOut < (int)prodWidth))
    {
        cout << "   multiplier is truncated by " << (int)prodWidth-wOut << " bits (err=" << (unsigned long)wX*(((unsigned long)1<<((int)wOut-guardBits))) << "), ensure sufficient precision..." << endl;
        //cout << "   guardBits=" << guardBits << endl;
        //cout << "   g=" << guardBits << " k=" << keepBits << " errorBudget=" << errorBudget << " difference to conservative est: " << errorBudget-(long long)(((unsigned long)1)<<(prodWidth-(int)wOut-1)-1) << endl;

        stringstream nvarName;
        nvarName << "C";
        ScaLP::Variable Cvar = ScaLP::newIntegerVariable(nvarName.str());
        constVecTerm.add(Cvar, 1);

        // C = 2^i*c_i + 2^(i+1)*c_(i+1)...
        ScaLP::Term cTerm;
        vector<ScaLP::Variable> cVars(prodWidth-wOut-1);
        for(int i = 0; i < prodWidth-wOut-1; i++){
            stringstream cvarName;
            cvarName << "c" << i;
            //cout << cvarName.str() << endl;
            cVars[i] = ScaLP::newBinaryVariable(cvarName.str());
            cTerm.add(cVars[i],  (double)( 1ULL << i));
            //constVecTerm.add(cVars[i],  (double)( 1ULL << i));
            //obj.add(cVars[i], (double)0.65);    //append variable to cost function (not needed for combined optimization)
        }
        ScaLP::Constraint cConstraint = cTerm - Cvar == 0;
        stringstream cName;
        cName << "cCons";
        cConstraint.name = cName.str();
        solver->addConstraint(cConstraint);

        //Limit value of constant to center the error around 0 to the error budget
        ScaLP::Constraint cLimConstraint = Cvar < errorBudget;
        stringstream cLimName;
        cLimName << "cLimCons";
        cLimConstraint.name = cLimName.str();
        solver->addConstraint(cLimConstraint);

        //Limit the error budget
        //cout << "  maxErr=" << errorBudget << endl;
        ScaLP::Constraint truncConstraint = maxEpsTerm - Cvar < errorBudget;
        stringstream consName;
        consName << "maxEps";
        truncConstraint.name = consName.str();
        solver->addConstraint(truncConstraint);
    }
    //Add rounding bit to constant vector
    if(wOut < (int)prodWidth) constVecTerm.add((double)( 1ULL << (prodWidth-wOut-1)));
    //When the optimal truncation is not used, the truncation error recentering constant is added to the constant vector
    if((wOut < (int)prodWidth) && !performOptimalTruncation) constVecTerm.add(centerErrConstant.get_d());

    //constant vector bits for sign extension and truncation constants
    vector<ScaLP::Variable> cvBits(prodWidth+5);
    if(signedIO || wOut < (int)prodWidth){
        for(int i = 0; i < prodWidth+5; i++){
            stringstream cvarName;
            cvarName << "v" << setfill('0') << setw(dpC) << i;
            //cout << cvarName.str() << " weight " << (double)(1ULL << i) << endl;
            cvBits[i] = ScaLP::newBinaryVariable(cvarName.str());
            constVecTerm.add(cvBits[i], -(double)(1ULL << i));
        }
        //Calculate individual constant vector bits
        ScaLP::Constraint cvConstraint = constVecTerm == 0;
        stringstream consName;
        consName << "cVec";
        cvConstraint.name = consName.str();
        solver->addConstraint(cvConstraint);
    }


    // BitHeap compression part of ILP formulation
    bitheap->final_add_height = 3;   //Select if dual or ternary adder is used in final stage of bitheap compression
    addFlipFlop();      //Add FF to list of compressors
    vector<vector<ScaLP::Variable>> bitsInColAndStage(s_max+1, vector<ScaLP::Variable>(bitsinColumn.size()+1));
    //ScaLP::Term selectLastStage;
    for(int s = 0; s <= s_max; s++){
        vector<ScaLP::Term> bitsinNextColumn(prodWidth + 1);
        vector<ScaLP::Term> bitsinCurrentColumn(prodWidth + 1);
        vector<ScaLP::Term> bitsinLastStageColumn(prodWidth + 1);
        vector<ScaLP::Term> rcdDependencies(prodWidth + 1);
/*        stringstream stage;                 //
        stage << "D_" << s;
        ScaLP::Variable stageV = ScaLP::newBinaryVariable(stage.str());
        selectLastStage.add(stageV, 1);*/
        for(unsigned c = 0; c <= prodWidth; c++){        //one bit more for carry of ripple carry adder
            if(s < s_max - 1){
                for(unsigned e = 0; e < possibleCompressors.size(); e++){
                    stringstream nvarName;
                    nvarName << "k_" << setfill('0') << setw(dpSt) << s << "_" << setfill('0') << setw(dpK) << e << "_" << setfill('0') << setw(dpC) << c;
                    //std::cout << nvarName.str() << endl;
                    ScaLP::Variable tempV = ScaLP::newIntegerVariable(nvarName.str(), 0, ScaLP::INF());
                    obj.add(tempV,  possibleCompressors[e]->area);    //append variable to cost function
                    for(int ce = 0; ce < (int) possibleCompressors[e]->getHeights() && ce < (int)bitsinCurrentColumn.size() - (int)c; ce++){   //Bits that can be removed by compressor e in stage s in column c for constraint C1
                        //cout << possibleCompressors[e]->getHeightsAtColumn((unsigned) ce, false) << " c: " << c+ce << endl;
                        bitsinCurrentColumn[c+ce].add(tempV, possibleCompressors[e]->getHeightsAtColumn((unsigned) ce, false));
                    }
                    for(int ce = 0; ce < (int) possibleCompressors[e]->getOutHeights() && ce < (int)bitsinNextColumn.size() - (int)c; ce++){   //Bits that are added by compressor e in stage s+1 in column c for constraint C2
                        //cout << possibleCompressors[e]->getOutHeightsAtColumn((unsigned) ce, false) << " c: " << c+ce << endl;
                        bitsinNextColumn[c+ce].add(tempV, possibleCompressors[e]->getOutHeightsAtColumn((unsigned) ce, false));
                    }
                }
            }
            if(s == s_max-1){                   //consideration of the ripple carry adder in final stage
                for(unsigned e = (c == prodWidth) ? possibleCompressors.size() + 2 : possibleCompressors.size() - 1; e <= ((c == prodWidth) ? possibleCompressors.size() + 2 : possibleCompressors.size() + 1); e++) {   //front element of ripple carry adder should only appear in the MSB column
                    stringstream nvarName;
                    nvarName << "k_" << setfill('0') << setw(dpSt) << s << "_" << setfill('0') << setw(dpK) << e << "_" << setfill('0') << setw(dpC) << c;      //for final FFs or placeholder for 2-input ripple carry adder
                    //std::cout << nvarName.str() << endl;
                    ScaLP::Variable tempV = ScaLP::newBinaryVariable(nvarName.str());
                    obj.add(tempV, (e == possibleCompressors.size()-1)?0.0:1);    //append variable to cost function, r.c.a.-area (cost) is 1 6LUT, FFs cost 0.0LUT
                    int takes_carry = (e == possibleCompressors.size() && bitheap->final_add_height == 2)?1:0;      //The dual input adder can process an additional carry input bit in its LSB column
                    bitsinCurrentColumn[c].add(tempV, ((e == possibleCompressors.size()-1)?1: (c == prodWidth) ? 0 : (bitheap->final_add_height + takes_carry) ));   //FFs remove only one bit from current stage, but the ripple carry adder two, the front element of ripple carry adder does not remove any bit but just provides the carry
                    bitsinNextColumn[c + 0].add(tempV, 1);
                    //if(c == bitsinColumn.size()-1)
                    //    bitsinNextColumn[c + 1].add(tempV, 1);
                    if(0 < c && possibleCompressors.size()+1 <= e)
                        rcdDependencies[c-1].add(tempV, -1);
                    if(c < bitsinColumn.size()-1 && possibleCompressors.size() <= e && e < possibleCompressors.size()+2)
                        rcdDependencies[c].add(tempV, 1);
                    if(c == bitsinColumn.size() && e == possibleCompressors.size()+2)
                        rcdDependencies[c].add(tempV, 1);
                }
            }
            stringstream curBits;
            if(bitsInColAndStage[s][c] == nullptr){                                                                 //N_s_c: Bits that enter current compressor stage
                curBits << "N_" << s << "_" << c;
                //cout << curBits.str() << endl;
                bitsInColAndStage[s][c] = ScaLP::newIntegerVariable(curBits.str(), 0, ScaLP::INF());
            }
            if(s == 0 && bitsInColAndStage[s][c] != nullptr){
                stringstream consName0;
                consName0 << "C0_" << s << "_" << c;
                bitsinColumn[c].add(bitsInColAndStage[s][c], -1);      //Add variable to account for the output bits from sub-multipliers
                if(cvBits[c] != nullptr) bitsinColumn[c].add(cvBits[c], 1);      //Output bits from constant vector
                ScaLP::Constraint c0Constraint = bitsinColumn[c] == 0;     //C0_s_c
                c0Constraint.name = consName0.str();
                solver->addConstraint(c0Constraint);
            }

            if(s < s_max){
                stringstream consName1, consName2, zeroBits, nextBits;
                consName1 << "C1_" << s << "_" << c;
                zeroBits << "Z_" << setfill('0') << setw(dpSt) << s << "_" << setfill('0') << setw(dpC) << c;
                //cout << zeroBits.str() << endl;
                bitsinCurrentColumn[c].add(ScaLP::newIntegerVariable(zeroBits.str(), 0, ScaLP::INF()), -1);      //Unused compressor input bits, that will be set zero
                bitsinCurrentColumn[c].add(bitsInColAndStage[s][c], -1);      //Bits arriving in current stage of the compressor tree
                ScaLP::Constraint c1Constraint = bitsinCurrentColumn[c] == 0;     //C1_s_c
                c1Constraint.name = consName1.str();
                solver->addConstraint(c1Constraint);
                consName2 << "C2_" << s << "_" << c;
                //cout << consName2.str() << endl;
                if(bitsInColAndStage[s+1][c] == nullptr){
                    nextBits << "N_" << s+1 << "_" << c;
                    //cout << nextBits.str() << endl;
                    bitsInColAndStage[s+1][c] = ScaLP::newIntegerVariable(nextBits.str(), 0, ScaLP::INF());
                }
                bitsinNextColumn[c].add(bitsInColAndStage[s+1][c], -1); //Output Bits of compressors to next stage
                ScaLP::Constraint c2Constraint = bitsinNextColumn[c] == 0;     //C2_s_c
                c2Constraint.name = consName2.str();
                solver->addConstraint(c2Constraint);
            }
            if(s_max-1 <= s){
                if(s_max == s) {
                    stringstream consName3;
                    consName3 << "C3_" << s << "_" << c;
                    //bitsinLastStageColumn[c].add(stageV, 10000);
                    bitsinLastStageColumn[c].add(bitsInColAndStage[s][c], 1);
                    //if(s < s_max - 1)
                    //    bitsinLastStageColumn[c].add(bitsInColAndStage[s+1][c], 4);
                    ScaLP::Constraint c3Constraint =
                            bitsinLastStageColumn[c] <= ((s == s_max - 1) ? bitheap->final_add_height : 1);     //C3_s_c
                    c3Constraint.name = consName3.str();
                    solver->addConstraint(c3Constraint);
                    //cout << consName3.str() << " " << stage.str() << endl;
                }
                if(0 < c && s < s_max){
                    stringstream consName5;
                    consName5 << "C5_" << s << "_" << c;
                    ScaLP::Constraint c5Constraint = rcdDependencies[c-1] == 0;     //C5_s_c
                    c5Constraint.name = consName5.str();
                    solver->addConstraint(c5Constraint);
                }
            }

        }
    }
 /*   ScaLP::Constraint c4Constraint = selectLastStage == 1;     //C4
    c4Constraint.name = "C4";
    solver->addConstraint(c4Constraint);
*/
    // Set the Objective
    cout << "   setting objective (minimize cost function)..." << endl;
    solver->setObjective(ScaLP::minimize(obj));

    // Write Linear Program to file for debug purposes
    cout << "   writing LP-file for debuging..." << endl;
    solver->writeLP("tile.lp");
}


    bool TilingAndCompressionOptILP::addFlipFlop(){
        //BasicCompressor* flipflop;
        bool foundFlipflop = false;
        for(unsigned int e = 0; e < possibleCompressors.size(); e++){
            //inputs
            if(possibleCompressors[e]->getHeights() == 1 && possibleCompressors[e]->getHeightsAtColumn(0) == 1){
                if(possibleCompressors[e]->getOutHeights() == 1 && possibleCompressors[e]->getOutHeightsAtColumn(0) == 1){
                    foundFlipflop = true;
                    flipflop = possibleCompressors[e];
                    break;
                }
            }
        }

        if(!foundFlipflop){
            //add flipflop at back of possibleCompressor
            //vector<int> newVect;
            //BasicCompressor *newCompressor;
            //int col0=1;
            //newVect.push_back(col0);
            BasicCompressor *newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), vector<int> {1}, 0.0, "combinatorial", true);
            possibleCompressors.push_back(newCompressor);

            flipflop = newCompressor;
        }

        return !foundFlipflop;
    }

    void TilingAndCompressionOptILP::resizeBitAmount(unsigned int stages){

        stages++;	//we need also one stage for the outputbits

        unsigned int columns = bitAmount[bitAmount.size() - 1].size();
        //we need one stage more for the
        while(bitAmount.size() < stages){
            bitAmount.resize(bitAmount.size() + 1);
            bitAmount[bitAmount.size() - 1].resize(columns, 0);
        }
    }

    bool TilingAndCompressionOptILP::checkTruncationError(list<TilingStrategy::mult_tile_t> &solution, unsigned guardBits, mpz_class errorBudget, mpz_class constant){
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
                if(mulArea[x][y] == false)
                    truncError += (mpz_class(1)<<(x+y));
                //cout << ((mulArea[x][y] == true)?1:0);
            }
            //cout << endl;
        }

        if(truncError <= maxErr){
            cout << "OK: actual truncation error=" << truncError << " is smaller than the max. permissible error=" << maxErr << " by " << maxErr-truncError << "." << endl;
            return true;
        } else {
            cout << "WARNING: actual truncation error=" << truncError << " is larger than the max. permissible error=" << maxErr << " by " << truncError-maxErr << "." << endl;
            return false;
        }
    }

    void TilingAndCompressionOptILP::parseTilingSolution(long long &optTruncNumericErr) {
        //parse tiling solution
        unsigned long long int new_constant = 0;
        double sum[4] = {0, 0, 0, 0};                   //variable to sum the rounded and double decision variables to check for numeric problems
        if(wOut < (int) prodWidth) cout << "centerErrConstant was: " << centerErrConstant << endl;//new error re-centering constant for truncation from solution
        double total_cost = 0 ;
        int dsp_cost = 0, own_lut_cost=0;
        for(auto &p:bestResult.values)
        {
            if(p.second > 0.5){     //parametrize all multipliers at a certain position, for which the solver returned 1 as solution, to flopoco solution structure
                string var_name = p.first->getName();
                //cout << var_name << "\t " << p.second << endl;
                switch(var_name.substr(0,1).at(0)){
                    case 'd':{
                        int mult_id = stoi(var_name.substr(1, dpS));
                        int x_negative = (var_name.substr(1 + dpS, 1) == "m") ? 1 : 0;
                        int m_x_pos = stoi(var_name.substr(1 + dpS + x_negative, dpX)) * ((x_negative) ? (-1) : 1);
                        int y_negative = (var_name.substr(1 + dpS + x_negative + dpX, 1) == "m") ? 1 : 0;
                        int m_y_pos = stoi(var_name.substr(1 + dpS + dpX + x_negative + y_negative, dpY)) * ((y_negative) ? (-1) : 1);
                        cout << "is true:  " << setfill(' ') << setw(dpY) << mult_id << " " << setfill(' ') << setw(dpY) << m_x_pos << " " << setfill(' ') << setw(
                                dpY) << m_y_pos << " cost: " << setfill(' ') << setw(5) << tiles[mult_id]->getLUTCost(m_x_pos, m_y_pos, wX, wY, signedIO) << endl;

                        total_cost += (double) tiles[mult_id]->getLUTCost(m_x_pos, m_y_pos, wX, wY, signedIO);
                        own_lut_cost += tiles[mult_id]->ownLUTCost(m_x_pos, m_y_pos, wX, wY, signedIO);
                        dsp_cost += tiles[mult_id]->getDSPCost();
                        auto coord = make_pair(m_x_pos, m_y_pos);
                        TilingStrategy::solution.push_back(make_pair(
                                tiles[mult_id]->getParametrisation().tryDSPExpand(m_x_pos, m_y_pos, wX, wY, signedIO), coord));
                        break;
                    }
                    case 'k':{
                        //the compressors are processed in compressionAlgorithm()
                        break;
                    }
                    case 'c':{
                        int c_id = stoi(var_name.substr(1, dpC));
                        new_constant |= (1ULL << c_id);
                        cout << var_name << " pos " << c_id << " dpC " << dpC << endl;
                        break;
                    }
                    default:
                        continue;
                }
            }
            //check variables for numeric derivations due to rounding for optimal truncation.
            if (p.first->getName().substr(1, 1) == "b") {
                int m_x_pos = stoi(p.first->getName().substr(2, dpX)) ;
                int m_y_pos = stoi(p.first->getName().substr(2 + dpX, dpY));
                sum[0] += p.second * (double)(1LL << (m_x_pos + m_y_pos));
                if (p.second >= 0.5) sum[1] += (double)(1LL << (m_x_pos + m_y_pos));
            }
            if (p.first->getName().substr(0, 1) == "c") {
                int shift = stoi(p.first->getName().substr(1, dpC));
                sum[2] += p.second * (1 << shift);
                if (p.second >= 0.5) sum[3] += (1 << shift);
            }
        }
        cout << "Total LUT cost:" << total_cost << endl;
        cout << "Own LUT cost:" << own_lut_cost << endl;
        cout << "Total DSP cost:" << dsp_cost << endl;

        if(performOptimalTruncation){
            //refresh centerErrConstant according to ILP solution
            mpz_import(centerErrConstant.get_mpz_t(), 1, -1, sizeof new_constant, 0, 0, &new_constant);
            cout << "centerErrConstant now is: " << centerErrConstant << endl;

            //check for numeric errors in solution
            optTruncNumericErr = (long long)ceil(fabs((sum[0] + sum[2]) - (sum[1] + sum[3])));
            if(0 < optTruncNumericErr){
                cout << "Warning: Numeric problems in solution, repeating ILP with Offset for Error of " << optTruncNumericErr << endl;
                errorBudget -= optTruncNumericErr;      //Reduce errorBudget by scope of numeric derivation for next iteration
            }
        }
    }

    void TilingAndCompressionOptILP::writeSolutionFile(ScaLP::Result result) {
        ofstream result_file;
        result_file.open("result.txt");
        result_file << result;
        result_file.close();
    }

#endif


}   //end namespace flopoco
