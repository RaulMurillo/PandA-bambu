/*
  A normalizer (leading zero/one counter + shifter + optional sticky bit) for FloPoCo

  Author: Florent de Dinechin

   This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2011.
  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Normalizer.hpp"

using namespace std;



namespace flopoco{


	Normalizer::Normalizer(OperatorPtr parentOp_, Target* target_, int wX_, int wR_, int maxShift_, bool computeSticky_, const int countType_) :
		Operator(parentOp_, target_), wX(wX_), wR(wR_), maxShift(maxShift_),  computeSticky(computeSticky_), countType(countType_) {

		wCount = intlog2(maxShift);  
		// -------- Parameter set up -----------------
		srcFileName = "Normalizer";
		setCopyrightString("Florent de Dinechin, (2007-2020)");

		REPORT(DETAILED, "wX="<<wX << " wR="<<wR << " maxShift="<<maxShift << " computeSticky=" << computeSticky  << " countType=" << countType);

		if(wR>wX) {
			THROWERROR( "As far as we know, wR>wX doesn't make sense. If you believe otherwise, get in touch!");
		}
		if(wR==wX && computeSticky) {
			THROWERROR( "As far as we know, computeSticky only makes sense if wR<wX doesn't make sense. If you believe otherwise, get in touch!");
		}
		 
		ostringstream name;
		name << "Normalizer_" << (countType<0?"ZO":((countType>0)?"O":"Z"))
				 << (computeSticky?"Stk":"") << "_" << wX << "_"<<wR<<"_"<<maxShift;
		setNameWithFreqAndUID(name.str());

		addInput ("X", wX);
		if (countType==-1) addInput ("OZb"); /* if we generate a generic LZOC */
		addOutput("Count", wCount);
		addOutput("R", wR);
		if (computeSticky)   addOutput("Sticky"); /* if we require a sticky bit computation */



		vhdl << tab << declare(join("level",wCount), wX) << " <= X ;"   <<endl;
		if (countType==-1) vhdl << tab << declare("sozb") << "<= OZb;"<<endl;
		if ((computeSticky)&&(wR<wX))   vhdl << tab << declare(join("sticky",wCount)) << " <= '0' ;"<<endl; //init sticky


		// Now comes the main loop.
		// i is the level index. Level i counts 2^i bits, and shifts by 2^i
		int currLevSize=wX, prevLevSize=0;
		for (int i=wCount-1; i>=0; i--){
			prevLevSize = currLevSize;

			// level(k) = max ( max (2^k, wR) + 2^k -1) , wX)
			currLevSize = (wR>intpow2(i)?wR:intpow2(i));
			currLevSize += (intpow2(i)-1);
			currLevSize = (currLevSize > wX? wX: currLevSize);

			// Delay evaluation.
			// As we output the count bits, their computation will not be merged inside the shift
			//REPORT( DEBUG, "currSize="<<currLevSize);

			double countBitDelay;
			if (countType>=0)
				countBitDelay = getTarget()-> eqConstComparatorDelay( intpow2(i) )  ;
			else
				countBitDelay = getTarget()->eqComparatorDelay( intpow2(i) ) ;
			
			vhdl << tab << declare(countBitDelay, join("count",i))
					 << "<= '1' when " <<join("level",i+1)<<range(prevLevSize-1,prevLevSize - intpow2(i))<<" = "
					 <<"("<<prevLevSize-1<<" downto "<<prevLevSize - intpow2(i)<<"=>"<< (countType==-1? "sozb": countType==0?"'0'":"'1'")<<") else '0';"<<endl;

			// The shift will take at most one LUT delay per level. We don't take into account that shift level can be merged: TODO ? It seems non-trivial.
			double shiftDelay = getTarget()->logicDelay(3);
			vhdl << tab << declare(shiftDelay,join("level",i),currLevSize)
					 << "<= " << join("level",i+1)<<"("<<prevLevSize-1<<" downto "<< prevLevSize-currLevSize << ")"
					 << " when " << join("count",i) << "='0' else ";
			int l,r;
			l = prevLevSize - intpow2(i) - 1;
			r = (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 );
			if (l>=r)
				vhdl << join("level",i+1) << "("<<prevLevSize - intpow2(i) - 1 <<" downto "<< (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) <<")";
			
			if (prevLevSize - intpow2(i) < currLevSize )
				vhdl << (l>=r?" & ":"") << rangeAssign(currLevSize -(prevLevSize - intpow2(i))-1,0,"'0'");
			vhdl << ";"<<endl;

			if ((computeSticky)&&(wR<wX)) {

				// Delay computation. Here we try to compute as much of the sticky in each level.
				double levelStickyDelay;
				// n is the size on which we compute the sticky bit
				int n = max( prevLevSize-currLevSize,
										 (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - int(intpow2(i)) ) - currLevSize : 0 ))  ;
				if ( countType == -1 )
					levelStickyDelay= getTarget()->eqComparatorDelay(n);
				else
					levelStickyDelay= getTarget()->eqConstComparatorDelay(n);
				

				vhdl << tab << declare(join("sticky_high_",i)) << "<= '0'";
				if (prevLevSize-currLevSize > 0)
					vhdl << "when " << join("level",i+1)<<"("<<prevLevSize-currLevSize -1 <<" downto "<< 0 <<") = CONV_STD_LOGIC_VECTOR(0,"<< prevLevSize-currLevSize <<") else '1'";
				vhdl << ";"<<endl;

   			vhdl << tab << declare(join("sticky_low_",i)) << "<= '0'";
				if ((currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) > 0)
					vhdl << "when " <<join("level",i+1)<<"("<<(currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) -1
						  <<" downto "<< 0 <<") = CONV_STD_LOGIC_VECTOR(0,"<< (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) <<") else '1'";
				vhdl << ";"<<endl;

				vhdl << tab << declare(levelStickyDelay, join("sticky",i))
						 << "<= " << join("sticky",i+1) << " or " << join("sticky_high_",i)
						 << " when " << join("count",i) << "='0' else " << join("sticky",i+1) << " or " << join("sticky_low_",i)<<";"<<endl;
			}

			vhdl <<endl;
		}

		//assign back the value to wR
		vhdl << tab << "R <= "<< join("level",0)
			  << (wR<=wX?"":join("&",rangeAssign(wR-wX-1,0,"'0'")))<<";"<<endl;


		vhdl << tab << declare("sCount",wCount) <<(wCount==1?"(0)":"")<<" <= ";
		for (int i=wCount-1; i>=0; i--){
			vhdl <<join("count",i);
			vhdl << (i>0?" & ":join(";","\n"));
		}

		if((1<<wCount)-1 > wX) {
			vhdl << tab << "Count <= sCount;"<<endl;

			/*			vhdl << tab << "Count <= CONV_STD_LOGIC_VECTOR("<<wX_<<","<<wCount_<<") when sCount=CONV_STD_LOGIC_VECTOR("<<intpow2(wCount_)-1<<","<<wCount_<<")"<<endl
				  << tab << tab << "else sCount;"<<endl;*/
		}
		else {
			vhdl << tab << "Count <= sCount;"<<endl;
		}


		if (computeSticky){
			if (wR>=wX)
				vhdl << tab << "Sticky <= '0';"<<endl;
			else
				vhdl << tab << "Sticky <= sticky0;"<<endl;
		}
		REPORT( DEBUG, "Leaving Normalizer");

	}

	Normalizer::~Normalizer() {
	}


	int Normalizer::getCountWidth() const{
		return wCount;
	}



	void Normalizer::emulate(TestCase* tc)
	{
		mpz_class inputValue  = tc->getInputValue("I");

		mpz_class sozb = 42; //dummy value
		if (countType == -1)
			sozb = tc->getInputValue("OZb");

		int sticky=0;
		int count =0;
		mpz_class shiftOutputValue = inputValue;


		mpz_class bit = (countType == -1) ? sozb : (countType == 0 ? 0 : 1); /* what are we counting in the specific case */

		int j=wX-1;
		while ((count < (1<<wCount)-1) &&  (j>=0)  && mpz_tstbit(inputValue.get_mpz_t(), j) == bit)   {
			count ++;
			j--;
			shiftOutputValue = shiftOutputValue <<1;
			}

		// Now reformat the output value to its size, and compute the sticky of the remaining bits.
		// The max size of shiftOutputValue is ((1<<wCount)-1) + wX
		mpz_class outputMask = (mpz_class(1) << wR) -1;
		int numBitsForSticky = wX - wR;
		if(numBitsForSticky >= 0) {// should be the typical use case where we need to compute a sticky bit
			mpz_class stickyMask = (mpz_class(1) << numBitsForSticky) -1;
			mpz_class bitsForSticky = shiftOutputValue & stickyMask;
			sticky = (bitsForSticky==0? 0 :1 );
			shiftOutputValue = (shiftOutputValue >> numBitsForSticky) & outputMask;
		}
		else {
			shiftOutputValue = (shiftOutputValue << numBitsForSticky) & outputMask;
			sticky=0;
		}
			
		
		tc->addExpectedOutput("O", shiftOutputValue);
		tc->addExpectedOutput("Count", count);

		if (computeSticky)
			tc->addExpectedOutput("Sticky",sticky);
	}


	

	OperatorPtr Normalizer::parseArguments(OperatorPtr parentOp, Target *target, std::vector<std::string> &args) {
		int wX, wR, maxShift, countType;
		bool computeSticky;
		UserInterface::parseStrictlyPositiveInt(args, "wX", &wX);
		UserInterface::parseStrictlyPositiveInt(args, "wR", &wR);
		UserInterface::parseStrictlyPositiveInt(args, "maxShift", &maxShift);
		UserInterface::parseBoolean(args, "computeSticky", &computeSticky);
		UserInterface::parseInt(args, "countType", &countType);
		return new Normalizer(parentOp, target, wX, wR, maxShift, computeSticky, countType);
	}


	
	void Normalizer::registerFactory(){
		UserInterface::add("Normalizer", // name
											 "A combined leading zero/one counter and left shifter, useful for floating-point normalization.",
											 "ShiftersLZOCs",  // category
											 "", // see also
											 "wX(int): input size in bits;\
                        wR(int): output size in bits, with wR <= wX;\
                        maxShift(int): how many bits to count, with maxShift<= wX ;\
                        computeSticky(bool)=false: if true and wR<wX, a sticky bit is computed out of the discarded bits;\
                        countType(int)=-1:  0 to count zeroes, 1 to count ones, -1 to have a dynamic OZb input that tells what to count", // This string will be parsed
											 "", // no particular extra doc needed
											 Normalizer::parseArguments
											 ) ;
		
	}


}
