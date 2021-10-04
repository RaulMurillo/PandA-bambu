/*

  This file is part of the FloPoCo project
  initiated by the Aric team at Ecole Normale Superieure de Lyon
  and developed by the Socrate team at Institut National des Sciences Appliquées de Lyon

  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2008-2014.
  All rights reserved.

*/
#include <iostream>
#include <iomanip>

#include "FixHornerEvaluator.hpp"
#include "IntMult/FixMultAdd.hpp"

using namespace std;


	/*



		 This is a simplified version of the computation in the ASAP 2010 paper, simplified because x is in [-1,1)

		 S_d  =  a_d
		 S_i  =  a_i + x^C_i * S_{i+1}     \forall i \in \{ 0...d-1\}
		 p(y)      =  \sigma_0


		 2/  each Horner step may entail up to two errors: one when truncating X to Xtrunc, one when rounding/truncating the product. 
		 step i rounds to lsbMult[i], such that lsbMult[i] <=  polyApprox->LSB
		 Heuristic to determine lsbMult[i] is as follows:
		 - set them all to polyApprox->LSB.
		 - Compute the cumulated rounding error (as per the procedure below)
		 - while it exceeds the rounding error budget, decrease lsbMult[i], starting with the higher degrees (which will have the lowest area/perf impact)
		 - when we arrive to degree 1, increase again.

		 The error of one step is the sum of two terms: the truncation of X, and the truncation of the product.
		 The first is sometimes only present in lower-degree steps (for later steps, X may be used untruncated
		 For the second, two cases:
		 * If plainVHDL is true, we 
		   - get the full product, 
			 - truncate it to lsbMult[i]-1, 
			 - add it to the coefficient, appended with a rounding bit in position lsbMult[i]-1
			 - truncate the result to lsbMult[i], so we have effectively performed a rounding to nearest to lsbMult[i]:
			 epsilonMult = exp2(lsbMult[i]-1)
		 * If plainVHDL is false, we use a FixMultAdd faithful to lsbMult[i], so the mult rounding error is twice as high as in the plainVHDL case:
			 epsilonMult = exp2(lsbMult[i])


		 We still should consider the DSP granularity to truncate x to the smallest DSP-friendly size larger than |sigma_{j-1}|
		 TODO

a 		 So all that remains is to compute the parameters of the FixMultAdd.
		 We need
		 * size of \sigma_d: MSB is that of  a_{d-j}, plus 1 for overflows (sign extended).
							 LSB is lsbMult[i]

		 * size of x truncated x^T_i : 
		 * size of the truncated result:
		 * weight of the LSB of a_{d-j}
		 * weight of the MSB of a_{d-j}
		 */


#define LARGE_PREC 1000 // 1000 bits should be enough for everybody


namespace flopoco{




	FixHornerEvaluator::FixHornerEvaluator(OperatorPtr parentOp, Target* target, 
																				 int lsbIn_,
																				 int msbOut_,
																				 int lsbOut_,
																				 vector<BasicPolyApprox*> poly_,
																				 bool finalRounding):
		FixPolyEval(parentOp, target, lsbIn_, msbOut_, lsbOut_, poly_, finalRounding)
	{
		initialize();
		computeArchitecturalParameters();
		generateVHDL();
	}





	void FixHornerEvaluator::initialize(){
		setNameWithFreqAndUID("FixHornerEvaluator");		
		setCopyrightString("F. de Dinechin (2014-2020)");
		srcFileName="FixHornerEvaluator";
	}



	void FixHornerEvaluator::computeArchitecturalParameters(){
		// initialize the worst case parameters so that we can use array notation. All dummy values
		wcSumSign = vector<int>(degree+1, 17);
		wcMaxAbsSum = vector<double>(degree+1, -1);	
		wcSumMSB = vector<int>(degree+1,INT_MIN);
		wcSumLSB = vector<int>(degree+1,INT_MAX);
		wcYLSB = vector<int>(degree,INT_MAX);
		isZero = vector<bool>(degree+1, true) ;
		//wcProductMSB = vector<int>(degree,0);
		//wcProductLSB = vector<int>(degree,0);

		sollya_obj_t yS = sollya_lib_build_function_free_variable();		
		sollya_obj_t rangeS = sollya_lib_parse_string("[-1;1]");		

		REPORT(DEBUG, "Entering computeArchitecturalParameters, for " << poly.size() << " intervals" );
	
		// iterate over all the polynomials to implement the error analysis on each interval
		for (size_t k=0; k<poly.size(); k++) {
			REPORT(DETAILED, "Error analysis on interval " << k  << " of " << poly.size()-1 );
			// detecting zero coefficients
			for(int i=0; i<=degree; i++) {
				if(!poly[k] ->  getCoeff(i) -> isZero())
					isZero[i]=false;
			}
			//for(int i=0; i<=degree; i++) {	REPORT(0, "isZero["<< i << "] = " << isZero[i]);  }

			// First, compute the max abs value of the d intermediate sums
			vector<double> maxAbsSum(degree+1, -1);
			vector<int> sumMSB(degree+1, INT_MIN);
			// S_d=C_d in the ASA book
			sollya_obj_t sS =	sollya_lib_constant(poly[k] -> getCoeff(degree) -> fpValue);
			maxAbsSum[degree] = abs(mpfr_get_d(poly[k] -> getCoeff(degree) -> fpValue, MPFR_RNDN)); // should be round away from 0 but nobody will notice
			sumMSB[degree] = poly[k] ->  getCoeff(degree) -> MSB;
			wcSumMSB[degree]     = max(wcSumMSB[degree],     sumMSB[degree]); // this one is probably useless

			for(int i=degree-1; i>=0; i--) {
				int sumSign;
				sollya_obj_t cS = sollya_lib_constant(poly[k] -> getCoeff(i) -> fpValue);
				sollya_obj_t pS = sollya_lib_mul(yS, sS);
				sollya_lib_clear_obj(sS); // it has been used
				sS = sollya_lib_add(cS, pS);
				sollya_lib_clear_obj(pS); // it has been used
				sollya_lib_clear_obj(cS); // it has been used
				// REPORT(0, "interval " << k << ": expression of S_"<<i);
				// sollya_lib_printf("%b\n", sS);

				sollya_obj_t sIntervalS = sollya_lib_evaluate(sS,rangeS);
				sollya_obj_t supS = sollya_lib_sup(sIntervalS);
				sollya_obj_t infS = sollya_lib_inf(sIntervalS);
				mpfr_t supMP, infMP, tmp;
				mpfr_init2(supMP, 1000); // no big deal if we are not accurate here 
				mpfr_init2(infMP, 1000); // no big deal if we are not accurate here 
				mpfr_init2(tmp, 1000); // no big deal if we are not accurate here 
				sollya_lib_get_constant(supMP, supS);
				sollya_lib_get_constant(infMP, infS);

				if(mpfr_sgn(infMP) >=0 )
					sumSign = 1;
				else if(mpfr_sgn(supMP) <0 )
					sumSign = -1;
				else 
					sumSign = 0;
				
				// Now recompute the MSB explicitely.
				mpfr_abs(supMP, supMP, GMP_RNDU);
				mpfr_abs(infMP, infMP, GMP_RNDU);
				mpfr_max(supMP, infMP, supMP, GMP_RNDU); // now we have the supnorm
				maxAbsSum[i] = mpfr_get_d(supMP, GMP_RNDU);
				mpfr_log2(tmp, supMP, GMP_RNDU);
				mpfr_floor(tmp, tmp);
				sumMSB[i] = 1+ mpfr_get_si(tmp, GMP_RNDU); // 1+ because we assume signed arithmetic for s

				REPORT(DETAILED, "interval " << k <<":  maxAbsSum[" << i << "] = " << maxAbsSum[i] << "  sumMSB[" << i << "] = " << sumMSB[i]);
				
				sollya_lib_clear_obj(sIntervalS);
				sollya_lib_clear_obj(supS);
				sollya_lib_clear_obj(infS);
				mpfr_clears(supMP,infMP,tmp, NULL);

				if(wcSumSign[i]==17)
					wcSumSign[i]=sumSign;
				else {
					if(sumSign!=wcSumSign[i])
						wcSumSign[i]=0; // otherwise leave it as it is
				} 
		
				// Finally update the worst-case values 
				wcSumMSB[i]     = max(wcSumMSB[i],     sumMSB[i]);
				//				wcSumLSB[i]     = min(wcSumLSB[i],     sumLSB);
				//wcProductMSB[i] = max(wcProductMSB[i], pMSB[i]);
				//wcProductLSB[i] = min(wcProductLSB[i], pLSB[i]);
			} // end loop on degree 
			// and free remaining memory
			sollya_lib_clear_obj(sS);

			
			REPORT(DEBUG, "OK, now we have the max Si, we may implement the error analysis, for approxErrorBound="<< poly[k]->getApproxErrorBound());
			// initialization
			double evalErrorBudget = exp2(lsbOut-1) - poly[k]->getApproxErrorBound();
			int lsb = lsbOut;
			vector<int> lsbY(degree,0);

			bool evalErrorNotOK=true;
			while(evalErrorNotOK) {
				// we have delta_Y=2^lsbY, and we want to balance the error term maxS[i]deltaY
				// with delta_multadd = 2^(lsb-1) if plainVHDL, 2^lsb otherwise    
				for(int i=degree-1; i>=0; i--) {
					lsbY[i] = max(lsb - sumMSB[i+1], lsbIn); 
				}
				double evalerror = 0;
				for(int i=degree-1; i>=0; i--) {
					double multaddError, yTruncationError;
					if(getTarget()->plainVHDL()) {
						multaddError=exp2(lsb-1); // full multipliers => correct rounding
					}
					else {
						multaddError=exp2(lsb); // faithful truncated multipliers
					}
					if(lsbY[i]==lsbIn) {
						yTruncationError = 0.0; // no truncation, no error
					}
					else {
						yTruncationError = exp2(lsbY[i])*maxAbsSum[i+1];
					}
					evalerror += multaddError + yTruncationError;
				}
				if(evalerror < evalErrorBudget)
					evalErrorNotOK=false;
					
				REPORT(DETAILED, "Interval " << k  << "  evalErrorBudget=" << evalErrorBudget  << "   lsb=" << lsb
							 << " => evalError=" << evalerror << (evalErrorNotOK?":  increasing lsb... " : ":  OK!"));
				if(evalErrorNotOK) {
					lsb--;
				}
			}

#if 0 // this is just to check that the error computation is tight
			// If I plug this code the autotest fails 12% of the cases
			REPORT(0,"******************** Sabotage ! ********************************");
			lsb++;
			for(int i=degree-1; i>=0; i--) {
				lsbY[i] =  max(lsb - sumMSB[i+1], lsbIn);
			}	
#endif
			
			// OK, we have the lsbY and lsb for this polynomial, now update the worst-case
			for(int i=degree-1; i>=0; i--) {
				wcYLSB[i] = min(wcYLSB[i], lsbY[i]);
			}	
			for(int i=degree; i>=0; i--) {
				wcSumLSB[i] = min(wcSumLSB[i], lsb);
			}	
			
		} // closes the for loop on k (the intervals)
 
		// Final reporting
		REPORT(INFO, "Architecture parameters:")
		for(int i=degree-1; i>=0; i--) {
			REPORT(INFO,"  Horner step " << i << ": YLSB=" << setw(3) << wcYLSB[i]
						 << "   SSgn=" << setw(2) << wcSumSign[i] << " SMSB=" << setw(3) << wcSumMSB[i] << " SLSB=" << setw(3) << wcSumLSB[i]
						 << "\t Mult size " << 1-wcYLSB[i] << "x" << wcSumMSB[i+1]-wcSumLSB[i+1]+1)
			}	
		

		sollya_lib_clear_obj(yS); 
		sollya_lib_clear_obj(rangeS);
	} 



	

	
	void FixHornerEvaluator::generateVHDL(){
		addInput("Y", -lsbIn+1);
		vhdl << tab << declareFixPoint("Ys", true, 0, lsbIn) << " <= signed(Y);" << endl;
		for (int j=0; j<=degree; j++) {
			addInput(join("A",j), coeffMSB[j]-coeffLSB[j] +1);
		}

		// declaring outputs
		addOutput("R", msbOut-lsbOut+1);


		// convert the coefficients to signed. Remark: constant signs have been inserted by the caller
		for(int i=0; i<=degree; i++) {
			vhdl << tab << declareFixPoint(join("As", i), true, coeffMSB[i], coeffLSB[i])
					 << " <= " << "signed(" << join("A",i) << ");" <<endl;
		}

		// Initialize the Horner recurrence
		resizeFixPoint(join("S", degree), join("As", degree), wcSumMSB[degree], wcSumLSB[degree]);

		// Main loop of the Horner recurrence
		for(int i=degree-1; i>=0; i--) {
			resizeFixPoint(join("YsTrunc", i), "Ys", 0, wcYLSB[i]);

			//  assemble faithful operators (either FixMultAdd, or truncated mult)

			if(getTarget()->plainVHDL()) {	// no pipelining here
				int pMSB = 0 + wcSumMSB[i+1] + 1; // not attempting to save the MSB bit that could be saved.
				int pLSB = wcYLSB[i] + wcSumLSB[i+1];	
				vhdl << tab << declareFixPoint(join("P", i), true, pMSB,  pLSB)
						 <<  " <= "<< join("YsTrunc", i) <<" * S" << i+1 << ";" << endl;
				// Align before addition
				resizeFixPoint(join("Ptrunc", i), join("P", i), wcSumMSB[i], wcSumLSB[i]-1);
				if(isZero[i]) {
					vhdl << tab <<	declareFixPoint(join("Aext", i), true, wcSumMSB[i], wcSumLSB[i]) << " <= " <<  zg( wcSumMSB[i] -  wcSumLSB[i] +1) << ";" << endl; 
					}
				else {
					resizeFixPoint(join("Aext", i), join("As", i), wcSumMSB[i], wcSumLSB[i]-1); // -1 to make space for the round bit
				}

				vhdl << tab << declareFixPoint(join("SBeforeRound", i), true, wcSumMSB[i], wcSumLSB[i]-1)
						 << " <= " << join("Aext", i) << " + " << join("Ptrunc", i) << "+'1';" << endl;
				resizeFixPoint(join("S", i), join("SBeforeRound", i), wcSumMSB[i], wcSumLSB[i]);
			}

			else { // using FixMultAdd
				THROWERROR("Sorry, use the plainVHDL option until we revive FixMultAdd");
#if 0
				//REPORT(DEBUG, "*** iteration " << i << );
				FixMultAdd::newComponentAndInstance(this,
																						join("Step",i),     // instance name
																						join("XsTrunc",i),  // x
																						join("S", i+1), // y
																						join("As", i),       // a
																						join("S", i),   // result
																						wcSumMSB[i], wcSumLSB[i]
																						);
#endif
			}
		}
		if(finalRounding)
			resizeFixPoint("Rs", "S0",  msbOut, lsbOut);

		vhdl << tab << "R <= " << "std_logic_vector(Rs);" << endl;

	}



	
	FixHornerEvaluator::~FixHornerEvaluator(){}


}
