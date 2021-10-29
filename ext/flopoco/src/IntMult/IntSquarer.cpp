/*
  An integer squarer for FloPoCo

  Author: Bogdan Pasca

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "IntSquarer.hpp"
#include "BitHeap/BitHeap.hpp"

using namespace std;

namespace flopoco{

	/** compute the guard bits and the size of the last column for the faithful truncation of an exact bit heap.
			input: a vector bhc of bit heap columns, for an integer bit heap (lsb at position 0)
	 a positive integer l defining the position of the lsb of the result after truncation (overall error should be strictly smaller than 2^l 
returns: an integer pext that defines the position of the last column
	an integer t that defines the number of bits to keep in the last column */
	void compute_truncation_params(vector<int> bhc, int l, int& lext, int& t) {
		mpz_class bound = mpz_class(1) << l;
		mpz_class deltaMinor, newDeltaMinor;
		// First loop: determine lext.
		lext=-1;
		newDeltaMinor=0;
		do {
			lext+=1;
			deltaMinor=newDeltaMinor;
			newDeltaMinor += (mpz_class(bhc[lext])<<lext);
			//cerr << "newdeltaminor= " << newDeltaMinor << " lext=" << lext << " bound=" << bound << endl;
		}
    while(newDeltaMinor + (mpz_class(1)<<lext)  < bound);
		// cerr << "deltaminor= " << deltaMinor << " lext=" << lext << endl;
		// on exit deltaMinor is the weighted sum of all the bits of position < lext

		//  Second loop: lext and Deltaminor no longer move, increase t
		t=0;
		while (deltaMinor + (mpz_class(t+1)<<lext) < bound) {
			t+=1;
		}
	}


	vector<int> buildExactSquareColumnHeights(int wIn) {
		vector<int> bhc;

		// mostly for The Book, naive bit-level version
		for (int i=0; i<2*wIn;i++) { bhc.push_back(0); }
		// simpler version of the loops in the constructor.
		// There is some annoying redundancy here but I don't know how to avoid it.
		// The diagonal bits
		for (int i=0; i<wIn; i++) {
			bhc[2*i] += 1 ;
		}
		// The triangle bits 
		for (int i=0; i<wIn-1; i++) {
			for (int j=i+1; j<wIn; j++) {
				bhc[i+j+1] += 1;
			}
		}
		return bhc;
	}

	
	// a small aux function that tells us if this bit must be stored in the bit heap or if we can discard it.
	bool needThisBit(int pos, int lext, int&k) {
		if(pos>lext)
			return true;
		if(pos==lext && k>0) {
			k -=1;
			return true;
		}
		return false;
	}


	
	IntSquarer::IntSquarer(OperatorPtr parentOp_, Target* target_,  int wIn_, bool signedIn_, int wOut_):
		Operator(parentOp_, target_), wIn(wIn_), signedIn(signedIn_), wOut(wOut_)
	{
		srcFileName = "IntSquarer";
		if(signedIn){
			THROWERROR("Sorry, signedIn not implemented yet");
		}
		if(wOut==0){
			wOut=2*wIn; // valid whatever signedIn
		}
		ostringstream name;
		name << "IntSquarer_" << wIn << "_" << wOut;
		setNameWithFreqAndUID(name.str());
		setCopyrightString("Florent de Dinechin (2021)");

		vector<int> bhc = buildExactSquareColumnHeights(wIn);
		int l=0;    // initialization that works for the exact case
		int lext=0;
		int t=0; // trucated bits in the rightmost column
		guardBits=0;
		if(wOut!=2*wIn){
			l=2*wIn-wOut;
			compute_truncation_params(bhc,  l, lext,  t);
			guardBits=l-lext;
			REPORT(INFO, "Truncation parameters: lext=" << lext << " (g=" << guardBits << ")  t=" << t);
			faithfulOnly=true;
		} else {
			lext=0;
			faithfulOnly=false; // this is an exact squarer	
		}

		// If you add t++; here simulation fails so we seem to touch the optimal.
		int k  = bhc[lext]-t; // kept bits in the rightmost column
		


		// Set up the IO signals
		addInput ("X"  , wIn);
		addOutput("R"  , wOut);

		// use the fix-point constructor
		BitHeap bh(this, 2*wIn-1, lext);

		// Naive bit-level version, to get drawings for The Book
		// The diagonal bits
		for (int i=0; i<wIn; i++) {
			if(needThisBit(2*i, lext, k)) {
				string bit="x"+to_string(i)+"sq";
				vhdl << tab << declare(getTarget()->logicDelay(), bit) << " <= X" << of(i) << ";" <<endl;
				bh.addBit(bit, 2*i);
			}
		}
		// The triangle bits 
		for (int i=0; i<wIn-1; i++) {
			for (int j=i+1; j<wIn; j++) {
				if(needThisBit(i+j+1, lext, k)) {
					string bit="x"+to_string(i)+"x"+to_string(j);
					vhdl << tab << declare(getTarget()->logicDelay(), bit) << " <= X" << of(i) << " and X" << of(j) << ";" <<endl;
					bh.addBit(bit, i+j+1);
				}	
			}
		}
		// for (int i=lext; i<2*wIn; i++) {	cerr << i << " "  << bh.getColumnHeight(i) << endl ; }
		
		// the correction constant + round bit
		if(wOut!=2*wIn) {
			for (int i=lext; i<l; i++) {
				REPORT(0, "constant bit at pos " << i);
				bh.addConstantOneBit(i);
			}
		}


		if(wOut==2*wIn){ // sanity check
			for (int i=0; i<2*wIn; i++) {
				if(bhc[i] != bh.getColumnHeight(i)) {
					THROWERROR("failed sanity check on bit heap heights for i=" << i << " : " << bhc[i] << " vs " <<  bh.getColumnHeight(i));
				}
			}
		}

		
		bh.startCompression();
		string bhr=bh.getSumName();


		if(wOut==2*wIn) {
			vhdl << tab << "R <= " << bhr << ";" <<endl;
		}
		else {
			vhdl << tab << "R <= " << bhr<< range(wOut+guardBits-1, guardBits) << ";" <<endl;
		}			
#if 0 // Old code
		
		if (wIn <= 17 ) {
			vhdl << tab << declare( "sX", wIn) << " <= X;" << endl;
			vhdl << tab << declare( "sY", wIn) << " <= X;" << endl;

			manageCriticalPath( getTarget()->LogicToDSPWireDelay() + getTarget()->DSPMultiplierDelay() );
			vhdl << tab << "R <= sX * sY;" << endl;
			getOutDelayMap()["R"] = getCriticalPath();
		}

		else if ((wIn>17) && (wIn<=33)) {
			vhdl << tab << declare("x0_16",18) << " <= \"0\" & X" << range(16,0) << ";" << endl;
			vhdl << tab << declare("x17_32",18) << " <= \"00\" & " << zg(33-wIn,0) << " & X" << range(wIn-1,17) << ";" << endl;
			vhdl << tab << declare("x17_32_shr",18) << " <= \"0\" & " << zg(33-wIn,0) << " & X" << range(wIn-1,17) << " & \"0\"" << ";" << endl;


			manageCriticalPath( getTarget()->LogicToDSPWireDelay() + getTarget()->DSPMultiplierDelay() );
			vhdl << tab << declare("p0",36) << " <= x0_16 * x0_16;" <<endl;
			vhdl << tab << declare("p1_x2",36) << " <= x17_32_shr * x0_16;" <<endl;

			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("s1",36) << " <= p1_x2 + ( \"00000000000000000\" & p0" << range(35,17) <<");" <<endl;
			vhdl << tab << declare("p2",36) << " <= x17_32 * x17_32;" <<endl;

			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("s2",36) << " <= p2 + ( \"00000000000000000\" & s1" << range(35,17) <<");" <<endl;

			manageCriticalPath( getTarget()->DSPToLogicWireDelay() );
			vhdl << tab << "R <= s2" << range(2*wIn-34-1,0) << " & s1"<<range(16,0) << " & p0"<<range(16,0)<<";" << endl;
			getOutDelayMap()["R"] = getCriticalPath();

		}else if ((wIn>17) && (wIn<=34)) {
			//FIXME: bottleneck
			vhdl << tab << declare("x0_16",17) << " <= X" << range(16,0) << ";" << endl;
			if (wIn<34)
				vhdl << tab << declare("x17_33",17) << " <= "<<zg(34-wIn,0) << " & X" << range(wIn-1,17) << ";" << endl;
			else
				vhdl << tab << declare("x17_33",17) << " <= X" << range(33,17) << ";" << endl;

			manageCriticalPath( getTarget()->LogicToDSPWireDelay() + getTarget()->DSPMultiplierDelay() );
			vhdl << tab << declare("p0",34) << " <= x0_16 * x0_16;" <<endl;
			vhdl << tab << declare("f0",18) << " <= p0" << range(17,0) << ";" << endl;
			vhdl << tab << declare("p1",34) << " <= x17_33 * x0_16;" <<endl;

			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("s1",34) << " <= p1 + ( \"0\" & p0" << range(33,18) <<");" <<endl;
			vhdl << tab << declare("f1",16) << " <= s1" << range(15,0) << ";" << endl;
			vhdl << tab << declare("p2",34) << " <= x17_33 * x17_33;" <<endl;

			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("s2",34) << " <= p2 + ( \"0\" & s1" << range(33,16) <<");" <<endl;

			if (wIn<34)
				vhdl << tab << "R <= s2" << range(2*wIn-34-1,0) << " & f1 & f0;" << endl;
			else
				vhdl << tab << "R <= s2 & f1 & f0;" << endl;

			getOutDelayMap()["R"] = getCriticalPath();
		}
		else if ((wIn>34) && (wIn<=51)) {
			// --------- Sub-components ------------------

			if (wIn<51)
				vhdl << tab << declare("sigX",51) << "<= "<<zg(51-wIn,0) <<" & X;"<<endl;
			else
				vhdl << tab << declare("sigX",51) << "<= X;"<<endl;

			manageCriticalPath( getTarget()->LogicToDSPWireDelay() + getTarget()->DSPMultiplierDelay() );
			vhdl << tab << declare("x0_16_sqr",34) << "<= sigX" << range(16,0) << " * sigX" << range(16,0)<<";"<<endl;
			vhdl << tab << declare("x17_33_sqr",34) << "<= sigX" << range(33,17) << " * sigX" << range(33,17)<<";"<<endl;
			vhdl << tab << declare("x34_50_sqr",34) << "<= sigX" << range(50,34) << " * sigX" << range(50,34)<<";"<<endl;
			vhdl << tab << declare("x0_16_x17_33",34) << "<= sigX" << range(16,0) << " * sigX" << range(33,17)<<";"<<endl;
			vhdl << tab << declare("x0_16_x34_50_prod",34) << "<= sigX" << range(16,0) << " * sigX" << range(50,34)<<";"<<endl;
			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("x0_16_x34_50",34) << "<= x0_16_x17_33" << range(33,17) << " + x0_16_x34_50_prod;"<<endl;
			vhdl << tab << declare("x17_33_x34_50_prod",34) << "<= sigX" << range(33,17) << " * sigX" << range(50,34)<<";"<<endl;
			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("x17_33_x34_50",34) << "<= x0_16_x34_50" << range(33,17) << " + x17_33_x34_50_prod;"<<endl;

//			nextCycle(); ////////////////////////////////////////////////
			vhdl << tab << declare("op1",84) << "<= x34_50_sqr & x17_33_sqr & x0_16_sqr" << range(33,18)<<";"<<endl;
			vhdl << tab << declare("op2",84) << "<= \"0000000000000000\" & x17_33_x34_50 & x0_16_x34_50" << range(16,0)<<" & x0_16_x17_33" << range(16,0)<<";"<<endl;

			// WAS:			intadder = new IntAdder(target, 84, inDelayMap("X", getTarget()->DSPToLogicWireDelay() + getCriticalPath()  ) );
			intadder = new IntAdder(target, 84 );

			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << tab << instance(intadder, "ADDER1");

			syncCycleFromSignal("adderOutput", false);
			setCriticalPath( intadder->getOutputDelay("R"));
			getOutDelayMap()["R"] = getCriticalPath();
			vhdl << tab << "R <= adderOutput" << range(2*wIn-19,0) << " & x0_16_sqr" << range(17,0)<<";"<<endl;
		}
		else if ((wIn>51) && (wIn<=68) && (wIn!=53) ) {
			// --------- Sub-components ------------------

			if (wIn<68)
				vhdl << tab << declare("sigX",68) << "<= "<<zg(68-wIn,0) <<" & X;"<<endl;
			else
				vhdl << tab << declare("sigX",68) << "<= X;"<<endl;

			setCriticalPath( getMaxInputDelays(inputDelays) );

			manageCriticalPath( getTarget()->LogicToDSPWireDelay() + getTarget()->DSPMultiplierDelay() );
			vhdl << tab << declare("x0_16_x17_33",34) << "<= sigX" << range(16,0) << " * sigX" << range(33,17)<<";"<<endl;
			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("x0_16_sqr",36) << "<= (\"00\" & x0_16_x17_33" << range(15,0)<<" & \"000000000000000000\") + "
				  << "( \"0\" & sigX" << range(16,0) << ") * (\"0\" & sigX" << range(16,0)<<");"<<endl;
			vhdl << tab << declare("x17_33_x34_50",34) << "<= sigX" << range(33,17) << " * sigX" << range(50,34)<<";"<<endl;
			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("x17_33_sqr",36) << "<= (\"00\" & x17_33_x34_50" << range(15,0)<<" & x0_16_x17_33" << range(33,16) <<") + x0_16_sqr(34) + "
				  << "( \"0\" & sigX" << range(33,17) << ") * (\"0\" & sigX" << range(33,17)<<");"<<endl;
			vhdl << tab << declare("x51_67_x34_50",34) << "<= sigX" << range(67,51) << " * sigX" << range(50,34)<<";"<<endl;
			vhdl << tab << declare("x_0_16_34_50",34) << " <= ( sigX" << range(16,0) << ") * (sigX" << range(50,34)<<");"<<endl;
			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("x34_50_sqr",36) << "<= (\"00\" & x51_67_x34_50" << range(15,0)<<" & x17_33_x34_50" << range(33,16) <<") + x17_33_sqr(34) + "
				  << "( \"0\" & sigX" << range(50,34) << ") * (\"0\" & sigX" << range(50,34)<<");"<<endl;
			vhdl << tab << declare("x_0_16_51_67_pshift", 34) << " <= x_0_16_34_50" << range(33,17) << " + "
				  << "( sigX" << range(16,0) << ") * (sigX" << range(67,51)<<");"<<endl;
			manageCriticalPath( getTarget()->DSPCascadingWireDelay() + getTarget()->DSPAdderDelay() );
			vhdl << tab << declare("x51_67_sqr",34) << "<= ( \"00000000000000\" & x51_67_x34_50" << range(33,16) <<") + x34_50_sqr(34) + "
				  << "( sigX" << range(67,51) << ") * (sigX" << range(67,51)<<");"<<endl;
			vhdl << tab << declare("x_17_33_51_67_pshift", 34) << " <= x_0_16_51_67_pshift" << range(33,17) << " + "
				  << "( sigX" << range(33,17) << ") * (sigX" << range(67,51)<<");"<<endl;

//			nextCycle(); ////////////////////////////////////////////////
			vhdl << tab << declare("op1",101) << "<= x51_67_sqr & x34_50_sqr" << range(33,0) << " & x17_33_sqr" << range(33,1) <<  ";"<<endl;
			vhdl << tab << declare("op2",101) << "<="<< zg(101-68,0)<<" & x_17_33_51_67_pshift & x_0_16_51_67_pshift" << range(16,0)<<" & x_0_16_34_50" << range(16,0)<<";"<<endl;
			//WAS			intadder = new IntAdder(target, 101, inDelayMap("X", getTarget()->DSPToLogicWireDelay() + getCriticalPath() ));
			intadder = new IntAdder(target, 101);

			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << tab << instance(intadder, "ADDER1");

			syncCycleFromSignal("adderOutput", false);
			setCriticalPath( intadder->getOutputDelay("R"));

			vhdl << tab << "R <= adderOutput" << range(2*wIn-36,0) << " & x17_33_sqr" << range(0,0) << " & x0_16_sqr" << range(33,0)<<";"<<endl;
			getOutDelayMap()["R"] = getCriticalPath();
		}
		else if (wIn==53){
			//TODO -> port to new pipeline framework
			//instantiate a 51bit squarer
			intsquarer = new IntSquarer(target, 51);

			bool tempPipelineStatus = getTarget()->isPipelined();
			bool tempDSPStatus = getTarget()->hasHardMultipliers();
			getTarget()->setPipelined(false);
			getTarget()->setUseHardMultipliers(false);

			if (tempPipelineStatus)
				getTarget()->setPipelined();
			if (tempDSPStatus)
				getTarget()->setUseHardMultipliers(true);

			intadder = new IntAdder(target, 54);
			intadd2 = new IntAdder(target, 53);


			vhdl << tab << declare("sigX",53) << "<= X;"<<endl;

			inPortMapCst(intsquarer, "X", "sigX(50 downto 0)");
			outPortMap(intsquarer, "R", "out_Squarer_51");
			vhdl << tab << instance(intsquarer, "SQUARER51");


			vhdl << tab << declare("op1mul2",53) << "<= (\"00\" & sigX" << range(50,0) <<") when sigX(51)='1' else "<<zg(53,0) << ";"<<endl;
			vhdl << tab << declare("op2mul2",53) << "<= (\"0\" & sigX" << range(50,0) <<" & \"0\") when sigX(52)='1' else "<<zg(53,0) << ";"<<endl;

			nextCycle(); ////////////////////////////////////////////////

			inPortMap(intadd2, "X", "op1mul2");
			inPortMap(intadd2, "Y", "op2mul2");
			inPortMapCst(intadd2, "Cin", "'0'");
			outPortMap(intadd2, "R", "x51_52_times_x_0_50");
			vhdl << tab << instance(intadd2, "MULT2");

			syncCycleFromSignal("out_Squarer_51",true);

			vhdl << tab << declare("x51_52_sqr",4) << " <= sigX" << range(52,51) << " * sigX" << range(52,51) <<";"<<endl;

			nextCycle(); ////////////////////////////////////////////////
			vhdl << tab << declare("op1",54) << "<= x51_52_sqr & out_Squarer_51" << range(101,52) <<  ";"<<endl;
			vhdl << tab << declare("op2",54) << "<="<< zg(1,0)<<" & x51_52_times_x_0_50;"<<endl;

			inPortMap(intadder, "X", "op1");
			inPortMap(intadder, "Y", "op2");
			inPortMapCst(intadder, "Cin", "'0'");
			outPortMap(intadder, "R", "adderOutput");
			vhdl << tab << instance(intadder, "ADDER54");

			syncCycleFromSignal("adderOutput", false);

			vhdl << tab << "R <= adderOutput & out_Squarer_51" << range(51,0)<<";"<<endl;
		} else {
			cerr << " For the moment IntSquarer does not support inputs larger than 68 bits. " << endl;
			exit (EXIT_FAILURE);
		}
#endif
	}

	IntSquarer::~IntSquarer() {
	}




	void IntSquarer::emulate(TestCase* tc)
	{
		mpz_class svX = tc->getInputValue("X");
		if(signedIn) {
				svX = bitVectorToSigned(svX, wIn);
		}
		mpz_class svR = svX * svX ;
		if(faithfulOnly) {//faithful rounding
			mpz_class svRRD, svRRU;
			svRRD = svR >> (2*wIn-wOut); // truncation
			if (svR-(svRRD << (2*wIn-wOut))==0) { // if all shifted out bits were 0 
				svRRU = svRRD;
			}
			else {
				svRRU = svRRD+1;
			}
			tc->addExpectedOutput("R", svRRD);
			tc->addExpectedOutput("R", svRRU);
		}
		else { 
			if(wOut==2*wIn) {
				tc->addExpectedOutput("R", svR);
			}
			else {
				THROWERROR("TODO: manage correctly-rounded truncated squarers");
			}
		}
	}


	OperatorPtr IntSquarer::parseArguments(OperatorPtr parentOp, Target *target, std::vector<std::string> &args) {
		int wIn, wOut;
		bool signedIn;
		UserInterface::parseStrictlyPositiveInt(args, "wIn", &wIn);
		UserInterface::parseInt(args, "wOut", &wOut);
		UserInterface::parseBoolean(args, "signedIn", &signedIn);
		return new IntSquarer(parentOp,target, wIn, signedIn, wOut);
	}


	
	void IntSquarer::registerFactory(){
		UserInterface::add("IntSquarer", // name
											 "An integer squarer.",
											 "BasicInteger", // category
											 "", // see also
											 "wIn(int): size of input in bits;\
						            wOut(int)=0: size of the output if you want a truncated squarer. 0 for exact (full) squarer; \
						            signedIn(bool)=false: inputs can be signed or unsigned (output always unsigned);", // This string will be parsed
											 "", // no particular extra doc needed
											 IntSquarer::parseArguments
											 ) ;
	}


}
