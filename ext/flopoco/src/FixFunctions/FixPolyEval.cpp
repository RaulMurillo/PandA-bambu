/*

  This file is part of the FloPoCo project
  initiated by the Aric team at Ecole Normale Superieure de Lyon
  and developed by the Socrate team at Institut National des Sciences Appliquées de Lyon

  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2020-
  All rights reserved.

*/
#include <iostream>

#include "FixPolyEval.hpp"

using namespace std;


	/*
		This is an abstract class from which actual polynomial evaluators should inherit.
		 */



namespace flopoco{

    FixPolyEval::FixPolyEval(OperatorPtr parentOp_, Target* target_, 
														 int lsbIn_,
														 int msbOut_,
														 int lsbOut_,
														 vector<BasicPolyApprox*> poly_, // these should all be the same degree and all have the same format  
														 bool finalRounding_) :
			Operator(parentOp_, target_),
			poly(poly_),
			lsbIn(lsbIn_),
			msbOut(msbOut_),
			lsbOut(lsbOut_),
			finalRounding(finalRounding_)
		{
			// Extraction of shared polynomial data, + sanity checks
			// The degree
			degree = poly[0]->getDegree();
			for (auto p: poly) {
				if(p->getDegree() !=degree) {
					THROWERROR("Something wrong in FixPolyEval, it seems the polynomials do not all have the same degree");
				}
			}

			// The coefficient formats
			for(int i=0; i<=degree; i++) {
				coeffMSB.push_back( poly[0]->getCoeff(i)->MSB );
				coeffLSB.push_back( poly[0]->getCoeff(i)->LSB );
			}
			for (auto p: poly) {
				for(int i=0; i<=degree; i++) {
					if(coeffMSB[i] != p->getCoeff(i)->MSB) {
						THROWERROR("Something wrong in FixPolyEval, it seems the polynomials do not all have the same MSB for coeff "<< i);
					}
					if(coeffLSB[i] != p->getCoeff(i)->LSB) {
						THROWERROR("Something wrong in FixPolyEval, it seems the polynomials do not all have the same LSB for coeff "<< i);
					}
				}
			}
			
		}


	FixPolyEval::~FixPolyEval(){}


}
