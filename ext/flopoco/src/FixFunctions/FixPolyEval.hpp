/*
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2020 - .
  All rights reserved.

*/
#ifndef __FIXPOLYEVAL_HPP
#define __FIXPOLYEVAL_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "BasicPolyApprox.hpp"


namespace flopoco{

	/** An abstract polynomial evaluator.
	 It assumes the input X is an signed number in [-1, 1[ so lsbX=-wX.
	*/

  class FixPolyEval : public Operator
  {
  public:
    /** The constructor with manual control of all options.
     * @param    lsbIn input lsb weight, 
			 @param    msbOut  output MSB weight, used to determine wOut
			 @param    lsbOut  output LSB weight
			 @param    poly vector of BasicPolyApprox
			 @param    finalRounding If true, the result is rounded on msbout-lsbout+1 bits. If false, it is unrounded and the guard bits are output as well
     */

																				
    FixPolyEval(OperatorPtr parentOp, Target* target, 
								int lsbIn,
								int msbOut,
								int lsbOut,
								vector<BasicPolyApprox*> p, // these should all be the same degree and all have the same format  
								bool finalRounding=true);

		~FixPolyEval();
		
  protected:
		vector<BasicPolyApprox*> poly;
    int degree;                       /**< degree of the polynomial, extracted by the constructor */
		int lsbIn;                        /** LSB of input. Input is assumed in [0,1], so unsigned and MSB=-1 */
		int msbOut;                        /** MSB of output  */
		int lsbOut;                        /** LSB of output */
    vector<int> coeffMSB;             /**< vector of size degree: MSB weight for the coefficients of degree i; extracted by te constructor */
    vector<int> coeffLSB;             /**< vector of MSB weights for each coefficient, extracted by te constructor */
		bool finalRounding;               /** If true, the operator returns a rounded result (i.e. add the half-ulp then truncate)
																					If false, the operator returns the full, unrounded results including guard bits */

  };

}
#endif
