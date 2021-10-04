/*
  A barrel shifter for FloPoCo

  Authors: Florent de Dinechin, Bogdan Pasca

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
#include <gmp.h>
#include <mpfr.h>
#include <stdio.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "Shifters.hpp"

using namespace std;

// TODO there is a small inefficiency here, as most bits of s don't need to be copied all the way down

namespace flopoco{


	Shifter::Shifter(OperatorPtr parentOp, Target* target, int wX_, int maxShift_, ShiftDirection direction_, int wR_, bool computeSticky_, bool inputPadBit_) :
		Operator(parentOp, target), wX(wX_), maxShift(maxShift_), direction(direction_), wR(wR_), computeSticky(computeSticky_), inputPadBit(inputPadBit_)
	{
		setCopyrightString ( "Bogdan Pasca (2008-2011), Florent de Dinechin (2008-2019)" );
		srcFileName = "Shifters";

		//cout << endl << "! Shifter::Shifter, wout=" << wR << endl << endl;

		// -------- Parameter set up -----------------
		if(wR==-1){ // The user asked a sensible default for wR, here it is
			if(computeSticky)
				wR =  wX; // no information loss since we compute the sticky
			else
				wR = wX + maxShift;  // no information loss
		}
		//cout << endl << "!! Shifter::Shifter, wout=" << wR << endl << endl;

		//Sanity check -- there should probably be more
		if(wR>wX+maxShift)
			THROWERROR("Sorry, but a shifter of wX=" << wX << " input bits by maxShift=" << maxShift << " bits should have wR<=wX+maxShift=" << wX+maxShift << ". \n Somebody asked for wR=" << wR);
		if(wR<wX)
			THROWERROR("Sorry, but I can't see the point in a shifter with wR<wX. Please truncate your input first.  wX=" << wX << " wR=" << wR);
		
		ostringstream name;
		if(direction==Left) name <<"Left";
		else                 name <<"Right";
		name << "Shifter";
		if(computeSticky)
			name << "Sticky";
		name<<wX<<"_by_max_"<<maxShift;
		setNameWithFreqAndUID(name.str());

		
		REPORT(DETAILED, " wX="<<wX<<" maxShift="<<maxShift<<" direction="<< (direction == Right?  "RightShifter": "LeftShifter") );


		wShiftIn     = intlog2(maxShift);

		addInput ("X", wX);
		addInput ("S", wShiftIn);
		if(inputPadBit)
					addInput("padBit");
		addOutput("R", wR);
		if(computeSticky)
					addOutput("Sticky");

		string padbit = (inputPadBit? "padBit" : "'0'");
		//vhdl << tab << declare(getTarget()->localWireDelay(), "level0", wX) << "<= X;" << endl;
		//vhdl << tab << declare(getTarget()->localWireDelay(), "ps", wShiftIn) << "<= S;" << endl;
		vhdl << tab << declare("ps", wShiftIn) << "<= S;" << endl;

		// Pipelining
		// The theory is that the number of shift levels that the tools should be able to pack in a row of LUT-k is ceil((k-1)/2)
		// Actually there are multiplexers in CLBs etc that can also be used, so let us keep it simple at k/2. It works on Virtex6, TODO test on Stratix.
		int levelPerLut = (getTarget()->lutInputs()-1)/2; // this is a floor so it is the same as the formula above
		REPORT(DETAILED, "Trying to pack " << levelPerLut << " levels in a row of LUT" << getTarget()->lutInputs())
		int levelInLut=0;
		double levelDelay;
		double totalDelay=0; // for reporting

		if (computeSticky) {
			// TODO if wX<wR, begin with shifts without sticky, then start sticky computation
			// In the mean time we have padded the input to wR, hopefully with few bits only...
			if(wX<wR) {
				string pad="("+to_string(wR-wX-1) + " downto 0 => " + padbit + ")";
				vhdl << tab << declare("Xpadded", wR) << " <= "<< (direction == Left? pad+"&":"") << "X"<< (direction == Right? "&"+pad : "") <<";"  << endl;
			}
			else // with the earlier sanity checks we should have wX=wR
			 {
				 vhdl << tab << declare("Xpadded", wR) << " <= X;" << endl;
			}	
			// It is better to start the loop with larger bits so as to give time to sticky computation.
			// The datapath size is constant anyway
			vhdl << tab << declare(join("level", wShiftIn), wR) << "<= Xpadded;" << endl;
			for(int i=wShiftIn-1; i>=0; i--){
				levelInLut ++;
				if (levelInLut >= levelPerLut) {
					levelDelay=getTarget()->logicDelay();
					totalDelay += levelDelay;
					REPORT(DETAILED, "level delay is " << levelDelay << "   total delay (if no pipeline occurs) is " << totalDelay);
					levelInLut=0;
				}
				else {// this level incurs no delay
					levelDelay=0;
				}
				//first compute the sticky
				vhdl << tab << declare(getTarget()->logicDelay() + getTarget()->eqConstComparatorDelay(intpow2(i)),     join("stk", i)) 
						 << " <= '1' when (" << "level" << i+1 << range(intpow2(i)-1,0)  << "/=" << zg(intpow2(i)) << " and ps" << of(i) << "='1')";
				if(i<wShiftIn-1)
					vhdl << " or stk" << i+1 << " ='1'";
				vhdl	 << "   else '0';" <<  endl;
				// then the shift
				vhdl << tab << declare(levelDelay, 
															 join("level", i), wR) << " <= "
						 << " level" << i+1 << " when  ps" << of(i) << "='0'"
						 << "    else (" << intpow2(i)-1 << " downto 0 => " << padbit << ") ";
				if(wR-1>= intpow2(i)) // because when wR is a power of two, that's all folks 
					vhdl << "& level" << i+1 << range(wR-1, intpow2(i));
				vhdl << ";" << endl;
			}
		  vhdl << tab << "R <= level0;"<<endl;
			vhdl << tab << "Sticky <= stk0;"<<endl;
		}

		
		else{//  no sticky computation, better to start with small shifts to minimize the overall datapath
			
				vhdl << tab << declare("level0", wX) << "<= X;" << endl;
				for(int currentLevel=0; currentLevel<wShiftIn; currentLevel++){
					levelInLut ++;
					if (levelInLut >= levelPerLut) {
						levelDelay=getTarget()->logicDelay() + getTarget()->fanoutDelay(wX+intpow2(currentLevel+1)-1);
						totalDelay += levelDelay;
						REPORT(DETAILED, "level delay is " << levelDelay << "   total delay (if no pipeline occurs) is " << totalDelay);
						levelInLut=0;
					}
					else // this level incurs no delay
						levelDelay=0;
			
					REPORT (DEBUG, "delay of level " << currentLevel <<" is " << levelDelay); 
					ostringstream currentLevelName, nextLevelName;
					currentLevelName << "level"<<currentLevel;
					nextLevelName << "level"<<currentLevel+1;
					if (direction==Right){
						vhdl << tab << declare(levelDelay,
																	 nextLevelName.str(),wX+intpow2(currentLevel+1)-1 )
								 <<" <=  ("<<intpow2(currentLevel)-1 <<" downto 0 => " << padbit << ") & "<<currentLevelName.str()<<" when ps";
						if (wShiftIn > 1)
							vhdl << "(" << currentLevel << ")";
						vhdl << " = '1' else "
								 << tab << currentLevelName.str() <<" & ("<<intpow2(currentLevel)-1<<" downto 0 => '0');"<<endl;
					}
					else {//(direction==Left)
						if(computeSticky) {
							THROWERROR("Nobody ever asked for a left shifter that also computes a sticky bit. Please implement it.")
								}
						else {
							vhdl << tab << declare(levelDelay, nextLevelName.str(),wX+intpow2(currentLevel+1)-1 )
									 << "<= " << currentLevelName.str() << " & ("<<intpow2(currentLevel)-1 <<" downto 0 => '0') when ps";
							if (wShiftIn>1)
								vhdl << "(" << currentLevel<< ")";
							vhdl << "= '1' else "
									 << tab <<" ("<<intpow2(currentLevel)-1<<" downto 0 => "  << padbit << ") & "<< currentLevelName.str() <<";"<<endl;
						}
					}
					ostringstream lastLevelName;
					lastLevelName << "level"<<wShiftIn;
					if (direction==Right)
						vhdl << tab << "R <= "<<lastLevelName.str()<<"("<< wX + intpow2(wShiftIn)-1-1 << " downto " << wX + intpow2(wShiftIn)-1 - wR <<");"<<endl;
					else
						vhdl << tab << "R <= "<<lastLevelName.str()<<"("<< wR-1 << " downto 0);"<<endl;
				}
			}
		}

	Shifter::~Shifter() {
	}


	void Shifter::emulate(TestCase* tc)
	{
		mpz_class sx = tc->getInputValue("X");
		mpz_class ss = tc->getInputValue("S");
		mpz_class sr ;

		mpz_class shiftedInput = sx;
		int i;

		if (direction==Left){
			mpz_class shiftAmount = ss;
			for (i=0;i<shiftAmount;i++)
				shiftedInput=shiftedInput*2;

			for (i= wX+intpow2(wShiftIn)-1-1; i>=wR;i--)
				if ( mpzpow2(i) <= shiftedInput )
					shiftedInput-=mpzpow2(i);
		}else{
			mpz_class shiftAmount = maxShift-ss;

			if (shiftAmount > 0){
				for (i=0;i<shiftAmount;i++)
					shiftedInput=shiftedInput*2;
			}else{
				for (i=0;i>shiftAmount;i--)
					shiftedInput=shiftedInput/2;
			}
		}

		sr=shiftedInput;
		tc->addExpectedOutput("R", sr);
	}



	OperatorPtr Shifter::parseArguments(OperatorPtr parentOp, Target *target, std::vector<std::string> &args) {
		int wX, wR, maxShift;
		bool dirArg, computeSticky, inputPadBit;
		UserInterface::parseStrictlyPositiveInt(args, "wX", &wX);
		UserInterface::parseInt(args, "wR", &wR);
		UserInterface::parseStrictlyPositiveInt(args, "maxShift", &maxShift);
		UserInterface::parseBoolean(args, "dir", &dirArg);
		UserInterface::parseBoolean(args, "computeSticky", &computeSticky);
		UserInterface::parseBoolean(args, "inputPadBit", &inputPadBit);
		ShiftDirection dir = (dirArg?Shifter::Right:Shifter::Left);
		return new Shifter(parentOp, target, wX, maxShift, dir, wR, computeSticky, inputPadBit);
	}



	void Shifter::registerFactory(){
		UserInterface::add("Shifter", // name
											 "A flexible shifter.",
											 "ShiftersLZOCs",
											 "",
											 "wX(int): input size in bits;\
											  maxShift(int): maximum shift distance in bits;\
											  dir(bool): 0=left, 1=right;	\
											  wR(int)=-1: size of the shifted output , -1 means computed, will be equal to wX+maxShift;\
											  computeSticky(bool)=false: if true and wR<wX+maxShift, shifted-out bits are ORed into a sticky bit;\
											  inputPadBit(bool)=false: if true, add an input bit used for left-padding, as in sign extension",
											 "", // no particular extra doc needed
											 Shifter::parseArguments
											 ) ;

	}


}
