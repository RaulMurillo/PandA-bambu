/*

  A class that manages fixed-point  polynomial approximation for FloPoCo (and possibly later for metalibm).

	At the moment, only on [0,1], no possibility to have it on [-1,1].
	Rationale: no point really as soon as we are not in the degenerate case alpha=0.

  Authors: Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL

  All rights reserved.

*/


/*
	 The function is assumed to have inputs in [0,1]

	 Stylistic remark: use index i for the subintervals, and j for the degree

*/
#include "UniformPiecewisePolyApprox.hpp"
#include <sstream>
#include <iomanip>
#include <limits.h>
#include <float.h>

namespace flopoco{

	UniformPiecewisePolyApprox::UniformPiecewisePolyApprox(FixFunction *f_, double targetAccuracy_, int degree_):
		degree(degree_), f(f_), targetAccuracy(targetAccuracy_)
	{
		needToFreeF = false;
		srcFileName="UniformPiecewisePolyApprox"; // should be somehow static but this is too much to ask me
		build();
	}


	UniformPiecewisePolyApprox::UniformPiecewisePolyApprox(string sollyaString_, double targetAccuracy_, int degree_):
		degree(degree_), targetAccuracy(targetAccuracy_)
	{
		//  parsing delegated to FixFunction
		f = new FixFunction(sollyaString_, false /* on [0,1]*/);
		needToFreeF = true;
		srcFileName="UniformPiecewisePolyApprox"; // should be somehow static but this is too much to ask me
		build();
	}



	UniformPiecewisePolyApprox::~UniformPiecewisePolyApprox()
	{
		if(needToFreeF)	free(f);
	}



#if 0
	/** a local function to build g_i(x) = f(2^(-alpha)*x + i*2^-alpha) defined on [0,1] */

	sollya_obj_t buildSubIntervalFunction(sollya_obj_t fS, int alpha, int i){
		stringstream s;
		s << "(1b-" << alpha << ")*x + ("<< i << "b-" << alpha << ")";
		string ss = s.str(); // do not use c_str directly on the stringstream, it is too transient (?)
		sollya_obj_t newxS = sollya_lib_parse_string(ss.c_str());
		sollya_obj_t giS = sollya_lib_substitute(fS,newxS);
		sollya_lib_clear_obj(newxS);
		return giS;
	}
#else
	/** a local function to build g_i(x) = f(2^(-alpha-1)*x + i*2^(-alpha) + 2^(-alpha-1)) defined on [-1,1] */
	sollya_obj_t UniformPiecewisePolyApprox::buildSubIntervalFunction(sollya_obj_t fS, int alpha, int i){
		stringstream s;

		s << "(1b-" << alpha+1 << ")*x + ("<< i << "b-" << alpha << " + 1b-" << alpha+1 << ")";
		string ss = s.str(); // do not use c_str directly on the stringstream, it is too transient (?)
		sollya_obj_t newxS = sollya_lib_parse_string(ss.c_str());
		sollya_obj_t giS = sollya_lib_substitute(fS,newxS);

		sollya_lib_clear_obj(newxS);

		return giS;
	}
#endif


	
	// split into smaller and smaller intervals until the function can be approximated by a polynomial of degree given by degree.
	void UniformPiecewisePolyApprox::build()
	{
		ostringstream cacheFileName;
		cacheFileName << "PiecewisePoly_"<<vhdlize(f->description) << "_" << degree << "_" << targetAccuracy << ".cache";

		// Test existence of cache file, create it if necessary
		openCacheFile(cacheFileName.str());

		if(!cacheFile.is_open())
		{
			//********************** Do the work, then write the cache *********************
			sollya_obj_t fS = f->fS; // no need to free this one
			sollya_obj_t rangeS;

			rangeS  = sollya_lib_parse_string("[-1;1]");
			// TODO test with [-1,1] which is the whole point of current refactoring.
			// There needs to be a bit of logic here because rangeS should be [-1,1] by default to exploit signed arithmetic,
			// except in the case when alpha=0 because then rangeS should be f->rangeS (and should not be freed)

			// Limit alpha to 24, because alpha will be the number of bits input to a table
			// it will take too long before that anyway
			bool alphaOK;
			for (alpha=0; alpha<24; alpha++)
			{
				nbIntervals = 1<<alpha;
				alphaOK = true;
				REPORT(DETAILED, " Testing alpha=" << alpha );
				for (int i=0; i<nbIntervals; i++) {
					// The worst case is typically on the left (i==0) or on the right (i==nbIntervals-1).
					// To test these two first, we do this small rotation of i
					int ii=(i+nbIntervals-1) & ((1<<alpha)-1);

					// First build g_i(x) = f(2^(-alpha)*x + i*2^(-alpha))
					sollya_obj_t giS = buildSubIntervalFunction(fS, alpha, ii);

					if(DEBUG <= UserInterface::verbose)
						sollya_lib_printf("> UniformPiecewisePolyApprox: alpha=%d, ii=%d, testing  %b \n", alpha, ii, giS);
					// Now what degree do we need to approximate gi?
					int degreeInf, degreeSup;
					BasicPolyApprox::guessDegree(giS, rangeS, targetAccuracy, &degreeInf, &degreeSup);
					// REPORT(DEBUG, " guessDegree returned (" << degreeInf <<  ", " << degreeSup<<")" ); // no need to report, it is done by guessDegree()
					sollya_lib_clear_obj(giS);
					// For now we only consider degreeSup. Is this a TODO?
					if(degreeSup>degree) {
						REPORT(DEBUG, "   alpha=" << alpha << " failed." );
						alphaOK=false;
						break;
					}
				} // end for loop on i

				// Did we succeed?
				if (alphaOK)
					break;
			} // end for loop on alpha

			if (alphaOK)
				REPORT(INFO, "Found alpha=" << alpha);

			// Compute the LSB of each coefficient. Minimum value is:
			LSB = floor(log2(targetAccuracy*degree));
			REPORT(DEBUG, "To obtain target accuracy " << targetAccuracy << " with a degree-"<<degree
					<<" polynomial, we compute coefficients accurate to LSB="<<LSB);
			// It is pretty sure that adding intlog2(degree) bits is enough for FPMinimax.

			// The main loop starts with the almost hopeless LSB defined above, then tries to push it down, a0 first, then all the others.
			// If guessdegree returned an interval, it tries lsbAttemptsMax times, then gives up and tries to increase the degree.
			// Otherwise lsbAttemptsMax is ignored, but in practice success happens earlier anyway

			int lsbAttemptsMax = intlog2(degree)+1;
			int lsbAttempts=0; // a counter of attempts to move the LSB down, caped by lsbAttemptsMax
			// bool a0Increased=false; // before adding LSB bits to everybody we try to add them to a0 only.

			bool success=false;
			while(!success) {
				// Now fill the vector of polynomials, computing the coefficient parameters along.
				//		LSB=INT_MAX; // very large
				// MSB=INT_MIN; // very small
				approxErrorBound = 0.0;
				BasicPolyApprox *p;

				REPORT(DETAILED, "Computing the actual polynomials ");
				// initialize the vector of MSB weights
				for (int j=0; j<=degree; j++) {
					MSB.push_back(INT_MIN);
				}

				for (int i=0; i<nbIntervals; i++) {
					REPORT(DETAILED, " ... computing polynomial approx for interval " << i << " / "<< nbIntervals);
					// Recompute the substitution. No big deal.
					sollya_obj_t giS = buildSubIntervalFunction(fS, alpha, i);

					p = new BasicPolyApprox(giS, degree, LSB, true);
					poly.push_back(p);
					if (approxErrorBound < p->getApproxErrorBound()){
						REPORT(DEBUG, "   new approxErrorBound=" << p->getApproxErrorBound() );
						approxErrorBound = p->getApproxErrorBound();
					}
					if (approxErrorBound>targetAccuracy){
						break; // fail, exit for loop
					}

					// Now compute the englobing MSB for each coefficient
					for (int j=0; j<=degree; j++) {
						// if the coeff is zero, we can set its MSB to anything, so we exclude this case
						if (  (!p->getCoeff(j)->isZero())  &&  (p->getCoeff(j)->MSB > MSB[j])  )
							MSB[j] = p->getCoeff(j)->MSB;
					}
					
					// Now set the MSB and LSB of the zero coefficients, for consistency
					// This prevents a future crash in this case
					for (int j=0; j<=degree; j++) {
						// if the coeff is zero, we can set its MSB to anything, so we exclude this case
						if (p->getCoeff(j)->isZero()) {
							p->getCoeff(j)->MSB = MSB[j]; // before that it was set arbitrary to 1
							p->getCoeff(j)->LSB = LSB; // before that it was set arbitrary to 0
							p->getCoeff(j)->width = MSB[j]-LSB+1;
						}
					}
					
				} // end for loop on i

				if (approxErrorBound < targetAccuracy) {
					REPORT(INFO, " *** Success! Final approxErrorBound=" << approxErrorBound << "  is smaller than target accuracy: " << targetAccuracy  );
					success=true;
				}
				else {
					REPORT(INFO, "With LSB="<<LSB<<", approx error:" << approxErrorBound << " is larger than target accuracy: " << targetAccuracy
							<< ". Decreasing LSB and starting over. Thank you for your patience");
					//empty poly
					for (auto i:poly)
						free(i);
					while(!poly.empty())
						poly.pop_back();

					if(lsbAttempts<=lsbAttemptsMax) {
						lsbAttempts++;
						LSB--;
					}
					else {
						LSB+=lsbAttempts;
						lsbAttempts=0;
						alpha++;
						nbIntervals=1<<alpha;
						REPORT(INFO, "guessDegree mislead us, increasing alpha to " << alpha << " and starting over");
					}
				}
			} // end while(!success)

			// Now we need to resize all the coefficients of degree i to the largest one
			// TODO? we could also check if one of the coeffs is always positive or negative, and optimize generated code accordingly
			for (int i=0; i<nbIntervals; i++) {
				for (int j=0; j<=degree; j++) {
					// REPORT(DEBUG "Resizing MSB of coeff " << j << " of poly " << i << " : from " << poly[i] -> getCoeff(j) -> MSB << " to " <<  MSB[j]);
					poly[i] -> getCoeff(j) -> changeMSB(MSB[j]);
					// REPORT(DEBUG, "   Now  " << poly[i] -> getCoeff(j) -> MSB);
				}
			}

			// Write the cache file
			writeToCacheFile(cacheFileName.str());

			//cleanup the sollya objects
			sollya_lib_clear_obj(rangeS);
		}
		else
		{
			REPORT(INFO, "Polynomial data cache found: " << cacheFileName.str());
			//********************** Just read the cache *********************
			readFromCacheFile(cacheFileName.str());
			REPORT(INFO, "Polynomial data read: " << cacheFileName.str());
		} // end if cache

		// Check if all the coefficients of a given degree are of the same sign
		checkCoefficientsSign();


		// A bit of reporting
		createPolynomialsReport();

#if 0 // display the coefficients in pgfplot format			cout << "\addplot  coordinates { ";
		cout << "xmax=" << nbIntervals-1 << endl;
		for (int k=0; k<=degree; k++) {
			cout << "\\addplot[only marks, mark=+]  coordinates { ";
			for(int i=0; i<nbIntervals; i++) {
				cout << "(" << i << ", "
						 << poly[i]->getCoeff(k)->getConstantAsMPZ()
						 <<  ") ";
			}
			cout << "} ;" << endl
					 <<  "\\legend{$C_" << k << "$} ;" << endl;
		}
#endif		
	}

	
	mpz_class UniformPiecewisePolyApprox::getCoeffAsPositiveMPZ(int i, int d){
		BasicPolyApprox* p = poly[i];
		FixConstant* c = p->getCoeff(d);
		return c->getBitVectorAsMPZ();
	}


	void UniformPiecewisePolyApprox::openCacheFile(string cacheFileName)
	{
		cacheFile.open(cacheFileName.c_str(), ios::in);

		// check for bogus .cache file
		if(cacheFile.is_open() && cacheFile.peek() == std::ifstream::traits_type::eof())
		{
			cacheFile.close();
			std::remove(cacheFileName.c_str());
			cacheFile.open(cacheFileName.c_str(), ios::in); // of course this one will fail
		}
	}


	void UniformPiecewisePolyApprox::writeToCacheFile(string cacheFileName)
	{
		cacheFile.open(cacheFileName.c_str(), ios::out);

		cacheFile << "Polynomial data cache for " << cacheFileName << endl;
		cacheFile << "Erasing this file is harmless, but do not try to edit it." <<endl;

		cacheFile << degree <<endl;
		cacheFile << alpha <<endl;
		cacheFile << LSB <<endl;

		for (int j=0; j<=degree; j++) {
			cacheFile << MSB[j] << endl;
		}

		cacheFile << approxErrorBound << endl;

		// now write the coefficients themselves
		for(int i=0; i<(1<<alpha); i++)
		{
			for (int j=0; j<=degree; j++)
			{
				cacheFile <<  poly[i] -> getCoeff(j) -> getConstantAsMPZ() << endl;
			}
			cacheFile << poly[i] -> getApproxErrorBound() << endl;
		}
	}


	void UniformPiecewisePolyApprox::readFromCacheFile(string cacheFileName)
	{
		string line;

		getline(cacheFile, line); // ignore the first line which is a comment
		getline(cacheFile, line); // ignore the second line which is a comment

		cacheFile >> degree;
		cacheFile >> alpha;
		nbIntervals = 1<<alpha;
		cacheFile >> LSB;

		for (int j=0; j<=degree; j++) {
			int msb;
			cacheFile >> msb;
			MSB.push_back(msb);
		}

		cacheFile >> approxErrorBound;

		for (int i=0; i<(1<<alpha); i++) {
			vector<mpz_class> coeff;
			for (int j=0; j<=degree; j++) {
				mpz_class c;
				cacheFile >> c;
				coeff.push_back(c);
			}
			BasicPolyApprox* p = new BasicPolyApprox(degree,MSB,LSB,coeff);
			poly.push_back(p);
			double aeb;
			cacheFile >> aeb;
			p->setApproxErrorBound(aeb);
		}
	}

	void UniformPiecewisePolyApprox::checkCoefficientsSign()
	{
		for (int j=0; j<=degree; j++) {
			mpz_class mpzsign = (poly[0]->getCoeff(j)->getBitVectorAsMPZ()) >> (MSB[j]-LSB);
			coeffSigns.push_back((mpzsign==0?+1:-1));
			for (int i=1; i<(1<<alpha); i++) {
				mpzsign = (poly[i]->getCoeff(j)->getBitVectorAsMPZ()) >> (MSB[j]-LSB);
				int sign = (mpzsign==0 ? 1 : -1);
				if (sign != coeffSigns[j])
					coeffSigns[j] = 0;
				// the following line is probably never useful: if one of the coeff is the unique
				// negative number that has no positive opposite in two's complement,
				// then there is nothing to win to convert it to unsigned			 
				if (poly[i]->getCoeff(j)->getBitVectorAsMPZ() == (mpz_class(1) << (MSB[j]-LSB)))
					coeffSigns[j] = 0;
			}
		}
		// Currently we just report the constant signs in this class.
		// Their actual exploitation is delegated to the hardware-generating classes
	}
	

	void UniformPiecewisePolyApprox::createPolynomialsReport()
	{
		int totalOutputSize;

		REPORT(INFO,"Parameters of the approximation polynomials: ");
		REPORT(INFO,"  Degree=" << degree	<< "  alpha=" << alpha
				<< "    maxApproxErrorBound=" << approxErrorBound  << "    common coeff LSB="  << LSB);

		totalOutputSize=0;
		for (int j=0; j<=degree; j++) {
			int size = MSB[j]-LSB + (coeffSigns[j] == 0);
			totalOutputSize += size ;
			REPORT(INFO,"  Coeff"<<setw(2) << j<<":  signedMSB =" <<setw(3)<< MSB[j]
						 << (coeffSigns[j]==0? ",  variable sign " : ", constant sign "+string(coeffSigns[j]==1?"+":"-") ) 
						 << "   => stored size ="<<setw(3) << size << " bits"
						 );
		}

		REPORT(INFO, "  Total size of the table is " << nbIntervals << " x " << totalOutputSize << " = " << nbIntervals*totalOutputSize << " bits");
	}


	OperatorPtr UniformPiecewisePolyApprox::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		string f;
		double ta;
		int d;

		UserInterface::parseString(args, "f", &f);
		UserInterface::parseFloat(args, "targetAcc", &ta);
		UserInterface::parseInt(args, "d", &d);

		UniformPiecewisePolyApprox *ppa = new UniformPiecewisePolyApprox(f, ta, d);
		cout << "Accuracy is " << ppa->approxErrorBound << " ("<< log2(ppa->approxErrorBound) << " bits)";

		return NULL;
	}



	void UniformPiecewisePolyApprox::registerFactory()
	{
		UserInterface::add("UniformPiecewisePolyApprox", // name
											 "Helper/Debug feature, does not generate VHDL. Uniformly segmented piecewise polynomial approximation of function f, accurate to targetAcc on [0,1)",
											 "FunctionApproximation",
											 "",
											 "\
f(string): function to be evaluated between double-quotes, for instance \"exp(x*x)\";\
targetAcc(real): the target approximation errror of the polynomial WRT the function;\
d(int): the degree to use",
											 "",
											 UniformPiecewisePolyApprox::parseArguments
											 ) ;
	}



} //namespace
