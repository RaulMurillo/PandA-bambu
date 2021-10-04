
#include "CompressionStrategy.hpp"
#include <climits>
#include "PrimitiveComponents/Xilinx/XilinxGPC.hpp"
#include "assert.h"

using namespace std;

namespace flopoco{


	CompressionStrategy::CompressionStrategy(BitHeap *bitheap_) :
		bitheap(bitheap_)
	{
		stringstream s;

		 //for REPORT to work
		srcFileName = bitheap->getOp()->getSrcFileName() + ":BitHeap:CompressionStrategy";
		guid = Operator::getNewUId();
		s << bitheap->name << "_CompressionStrategy_" << guid;
		uniqueName_ = s.str();

		//no chunks already compressed
		compressionDoneIndex = 0;

		if(UserInterface::pipelineActive_) {
			//the initial delay between the soonest bit and the rest of the bits
			//  to be compressed is equal to the delay of a basic compressor
			compressionDelay = bitheap->op->getTarget()->tableDelay(bitheap->op->getTarget()->lutInputs() + 2, 5, true);
			stagesPerCycle = (1.0/bitheap->op->getTarget()->frequency()) / compressionDelay;
		}
		else {
			compressionDelay = 0;
			stagesPerCycle = INT_MAX;
		}

		REPORT(DEBUG, "compressionDelay is " << compressionDelay << " and stagesPerCycle is " << stagesPerCycle);

		//generate all the compressors that can be used
		generatePossibleCompressors();

		//create the plotter
		bitheapPlotter = new BitheapPlotter(bitheap);
	}

	CompressionStrategy::~CompressionStrategy(){
	}

	void CompressionStrategy::orderBitsByColumnAndStage()
	{
		REPORT(DEBUG, "in orderBitsByColumnAndStage")
		double objectiveDelay = 1.0/bitheap->getOp()->getTarget()->frequency();


		bitheap->op->schedule();
		bitheap->sortBitsInColumns();

		REPORT(DEBUG, "after scheduling of inputs the bits in the bitheap are the following ")
		printBitsInBitheap();
		//first get the minimal cycle and maximum cycle
		int minCycle = 1E6;
		int maxCycle = -1;
		Bit* minBit = nullptr;
		Bit* maxBit = nullptr;

		for(unsigned int i = 0; i < bitheap->bits.size(); i++)
		{
			for(unsigned int j = 0; j < bitheap->bits[i].size(); j++)
			{
				int currentCycle = getStageOfArrivalForBit(bitheap->bits[i][j]);

				REPORT(DEBUG, "bit at column " << i << " and j " << j << " is called " << bitheap->bits[i][j]->signal->getName() << ". cycle is " << bitheap->bits[i][j]->signal->getCycle() << " and critical path is " << bitheap->bits[i][j]->signal->getCriticalPath() << " currentCycle=" << currentCycle);


				if(currentCycle > maxCycle){
					maxCycle = getStageOfArrivalForBit(bitheap->bits[i][j]);
					maxBit = bitheap->bits[i][j];
				}

				if(currentCycle < minCycle){
					minCycle = getStageOfArrivalForBit(bitheap->bits[i][j]);
					minBit = bitheap->bits[i][j];
				}

			}
		}
		REPORT(DEBUG, "max is cycle " << maxCycle);
		REPORT(DEBUG, "min is cycle " << minCycle);
		assert(minBit != nullptr);
		assert(maxBit != nullptr);
		REPORT(DEBUG, "minBit is " << minBit->signal->getName() << " and maxBit is " << maxBit->signal->getName());

		REPORT(DEBUG, "objectiveDelay = " << objectiveDelay);
		REPORT(DEBUG, "compressionDelay = " << compressionDelay);
		//get the first and last stage, where bits arrive.

		unsigned int maxStage;
		unsigned int minStage;
		maxStage = getStageOfArrivalForBit(maxBit);
		minStage = getStageOfArrivalForBit(minBit);

		maxStage -= minStage;
		//in [0 ... maxStage] can bits arrive
		REPORT(DEBUG, "maxStage is " << maxStage << " with a offset(minStage) of " << minStage);

		orderedBits.resize(maxStage + 1);
		for(unsigned s = 0; s < orderedBits.size(); s++){
			orderedBits[s].resize(bitheap->width);
		}

		//add the bits to orderedBits
		for(unsigned int c = 0; c < bitheap->bits.size(); c++){
			for(unsigned int i = 0; i < bitheap->bits[c].size(); i++){
				unsigned int tempStage = getStageOfArrivalForBit(bitheap->bits[c][i]);
				tempStage -= minStage;  //subtract the offset
				orderedBits[tempStage][c].push_back(bitheap->bits[c][i]); //TODO: check if we set the new criticalPath accordingly
				REPORT(DEBUG, "added Bit " << bitheap->bits[c][i]->signal->getName() << " with cycle " << bitheap->bits[c][i]->signal->getCycle() << " and criticalPath " << bitheap->bits[c][i]->signal->getCriticalPath() << " to stage " << tempStage << " and column " << c);
			}
		}
	}

	void CompressionStrategy::fillBitAmounts(){
		bitAmount.resize(orderedBits.size());
		for(unsigned int s = 0; s < orderedBits.size(); s++){
			bitAmount[s].resize(orderedBits[s].size() + 4, 0);  //add four because we allow compressors to have outputbits in higher columns than the MSB of the bitheap. Those outputbits will be removed.
		}
		for(unsigned int s = 0; s < orderedBits.size(); s++){
			for(unsigned int c = 0; c < orderedBits[s].size(); c++){
				bitAmount[s][c] = orderedBits[s][c].size();
			}
		}
	}

	double CompressionStrategy::getCompressionEfficiency(unsigned int stage, unsigned int column, BasicCompressor* compressor, unsigned int middleLength){
		int inputBits = 0;
		int outputBits = 0;
		for(unsigned int c = 0; c < compressor->getHeights(middleLength); c++){
			if(bitAmount.size() > stage && bitAmount[stage].size() > column + c){
				//position exists, now check how many bits there are to compress
				if(bitAmount[stage][column + c] >= (int)compressor->getHeightsAtColumn(c, middleLength)){
					inputBits += compressor->getHeightsAtColumn(c, middleLength);
				}
				else{
					if(bitAmount[stage][column + c] > 0){ //preventing negative values (holes)
						inputBits += bitAmount[stage][column + c];
					}
				}
			}
		}
		for(unsigned int c = 0; c < compressor->getOutHeights(middleLength); c++){
			if(compressor->getOutHeightsAtColumn(c, middleLength) > 0){
				outputBits += compressor->getOutHeightsAtColumn(c, middleLength);
			}
		}
		double ratio = 0.0;
		if(compressor->getArea(middleLength) != 0.0){
			ratio = ((double)(inputBits - outputBits)) / compressor->getArea(middleLength);
		}   //else return 0.0
		return ratio;
	}

	void CompressionStrategy::placeCompressor(unsigned int stage, unsigned int column, BasicCompressor* compressor, unsigned int middleLength){
		for(unsigned int c = 0; c < compressor->getHeights(middleLength); c++){
			if(bitAmount.size() > stage && bitAmount[stage].size() > column + c){
				bitAmount[stage][column + c] -= compressor->getHeightsAtColumn(c, middleLength);
			}
		}
		for(unsigned int c = 0; c < compressor->getOutHeights(middleLength); c++){
			if(bitAmount.size() <= stage + 1){
				THROWERROR("tried to add bits into a stage, which does not exist");
			}
			else{
				if(bitAmount[stage + 1].size() > column + c){   //drop outputbits which go beyond the border
					bitAmount[stage + 1][column + c] += compressor->getOutHeightsAtColumn(c, middleLength);
				}
			}
		}
		solution.addCompressor(stage, column, compressor, middleLength);
	}

	void CompressionStrategy::printBitsInBitheap(){
		if(UserInterface::verbose < DEBUG){
			return;
		}
		for(unsigned int c = 0; c < bitheap->bits.size(); c++){
			for(unsigned int i = 0; i < bitheap->bits[c].size(); i++){
				REPORT(DEBUG, "bit at column " << c << " and i " << i << " is called " << bitheap->bits[c][i]->signal->getName() << ". cycle is " << bitheap->bits[c][i]->signal->getCycle() << " and critical path is " << bitheap->bits[c][i]->signal->getCriticalPath());
			}
		}
	}

	void CompressionStrategy::printBitAmounts(){
		if(UserInterface::verbose < DETAILED){
			return;
		}
		//first find max value length and maxcolumn
		unsigned int maxValueDigits = 1;
		unsigned int maxColumn = 1;
		for(unsigned int s = 0; s < bitAmount.size(); s++){
			if(bitAmount[s].size() > maxColumn){
				maxColumn = bitAmount[s].size();
			}
			for(unsigned int c = 0; c < bitAmount[s].size(); c++){
				unsigned int currentMaxValueDigits = 0;
				int value = bitAmount[s][c];
				if(value < 0){
					value *= (-1);
				}
				while(value > 0){
					currentMaxValueDigits++;
					value /= 10;
				}

				if(bitAmount[s][c] < 0){
					currentMaxValueDigits++;
				}

				if(currentMaxValueDigits > maxValueDigits){
					maxValueDigits = currentMaxValueDigits;
				}
			}

		}

		for(unsigned int s = 0; s < bitAmount.size(); s++){
			ostringstream line;
			line << "bitheap column height in stage " + to_string(s) << ": ";
			unsigned int tempColumn = maxColumn;
			int totalBits = 0;
			while(tempColumn > 0){
				tempColumn--;
				unsigned int addZeros = maxValueDigits + 1;
				if(bitAmount[s].size() > tempColumn){
					int value = bitAmount[s][tempColumn];
					while(value != 0){
						addZeros--;
						value /= 10;
					}
					if(bitAmount[s][tempColumn] <= 0){
						addZeros--; //- sign or zero
					}
					for(unsigned int k = 0; k < addZeros; k++){
						line << " ";
					}
					line << bitAmount[s][tempColumn] << " ";
					totalBits += bitAmount[s][tempColumn];
				}
			}
			cerr << line.str() << ", totalBits=" << totalBits << endl;
		}
	}

	void CompressionStrategy::printSolutionStatistics(){
		REPORT(DEBUG, "calculating compression area");
		unsigned int totalArea = 0;  //area in Compressor.hpp
		for(unsigned int s = 0; s < bitAmount.size() - 1; s++) {

			for (unsigned int c = 0; c < bitheap->width; c++) {

				vector<pair<BasicCompressor *, unsigned int> > tempVector;
				tempVector = solution.getCompressorsAtPosition(s, c);

				for (unsigned int j = 0; j < tempVector.size(); j++) {
					/*
					if (tempVector[j].first==tempVector[j+1].first){
						tempVector[j+1].second=+tempVector[j].second+1;
						continue;
					}
					*/
					REPORT(DEBUG, "applying " << tempVector[j].second+1 << " compressors of type " << tempVector[j].first->getStringOfIO() <<
												 " at stage " << s << " and column " << c << " with a combined LUT-area cost of " << (tempVector[j].second+1)*tempVector[j].first->getArea());
					totalArea+=(tempVector[j].second+1)*tempVector[j].first->getArea();


				}
			}
		}
		REPORT(DETAILED, "total area of the compression is equivalent to " << totalArea << " LUTs");
		REPORT(DETAILED, "total number of stages is " << bitAmount.size()-1);
	}



	void CompressionStrategy::applyAllCompressorsFromSolution(){
		REPORT(DEBUG, "applying all compressors");

		unsigned int colorCount = 0; //for bitheapPlotter
		for(unsigned int s = 0; s < bitAmount.size() - 1; s++){
			//get the zeros of this column
			vector<unsigned int> emptyInputs = solution.getEmptyInputsByStage(s);
			//get bits for plotter
			Bit *soonestBit, *soonestCompressibleBit;
			soonestBit = getSoonestBit(0, bitheap->width-1);
			soonestCompressibleBit = getSoonestCompressibleBit(0, bitheap->width - 1, compressionDelay);
			for(unsigned int c = 0; c < bitheap->width; c++){   //drop all compressors which start > MSB
				vector<pair<BasicCompressor*, unsigned int> > tempVector;
				tempVector = solution.getCompressorsAtPosition(s, c);
				REPORT(DEBUG, "at stage " << s << " and column " << c << " there are " << tempVector.size() << " compressors");
				for(unsigned int j = 0; j < tempVector.size(); j++){
					REPORT(DEBUG, "applying compressor " << tempVector[j].first->getStringOfIO());
					//applyCompressor
					Compressor* realCompressor = tempVector[j].first->getCompressor();  //TODO: consider middleLength
					//unsigned int middleLength = tempVector[j].second;
					vector<vector<Bit*> > tempBitVector;
					tempBitVector.resize(realCompressor->heights.size());
					//fill Bitvectors
					//note c is the starting column of the compressor in the bitheap. cTemp is the column within the compressor
					for(unsigned int cTemp = 0; cTemp < realCompressor->heights.size(); cTemp++){

						int maxBitsToCompress = realCompressor->heights[cTemp];
						//handle empty inputs
						if(maxBitsToCompress > 0 && emptyInputs[c + cTemp] > 0){
							if(emptyInputs[c + cTemp] >= (unsigned int) maxBitsToCompress){
								//there are more or the same amount of remaining empty inputs in this column than inputs of this compressor at this column
								//therefore set all the inputs of the compressor to zero
								emptyInputs[c + cTemp] -= (unsigned int) maxBitsToCompress;
								maxBitsToCompress = 0;
							}
							else{
								//some of the inputs of the compressor should be set to zero.
								maxBitsToCompress -= (int) emptyInputs[c + cTemp];
								emptyInputs[c + cTemp] = 0;
							}
						}
						REPORT(DEBUG, "after emptyInputs maxBitsToCompress is " << maxBitsToCompress);
						unsigned int bitsFound = 0;
						//check if the position of bits (c + cTemp) is inside the range of the bitheap
						if(c + cTemp < bitheap->bits.size()){
							for(unsigned int k = 0; k < bitheap->bits[c + cTemp].size(); k++){
								if(bitsFound < (unsigned) maxBitsToCompress){
									Bit* tempBit = bitheap->bits[c + cTemp][k];
									if(tempBit->type == BitType::free){
										bitsFound++;
										tempBitVector[cTemp].push_back(tempBit);
									}
								}
							}
						}
					}
					applyCompressor(tempBitVector, realCompressor, c);
					bitheap->op->schedule();
					bitheap->sortBitsInColumns();
				}
			}

			//assume that in every stage we compress bits
			colorCount++;
			bitheap->colorBits(BitType::justAdded, colorCount);
			bitheapPlotter->takeSnapshot(soonestBit, soonestCompressibleBit);

			bitheap->removeCompressedBits();
			bitheap->markBitsForCompression();

			bitheapPlotter->takeSnapshot(soonestBit, soonestCompressibleBit);

			REPORT(DEBUG, "finished stage " << s << " of " << bitAmount.size() - 1);
			REPORT(DEBUG, "the bits in the bitheap are ordered as following");
			printBitsInBitheap();
		}
		concatenateLSBColumns();

		for(unsigned int c = 0; c < bitheap->bits.size(); c++){
			if(bitheap->bits[c].size() > 3){
				THROWERROR("in column " << c << " are more than two bits (" << bitheap->bits[c].size() << ") for the final adder ");
			}
		}
	}



	bool CompressionStrategy::checkAlgorithmReachedAdder(unsigned int adderHeight, unsigned int stage){
		if(stage >= bitAmount.size()){
			THROWERROR("Doing the check if the algorithm for generating the compressortree is finished. Tried to access stage " << stage << " but there aren't that many stages");
		}
		for(unsigned int c = 0; c < bitAmount[stage].size(); c++){
			if(bitAmount[stage][c] > (int) adderHeight){
				return false;   //column exists where there are more bits left
			}
		}

		//now check if in all other stages there are no bits left
		for(unsigned int s = 0; s < bitAmount.size(); s++){
			if(s != stage){
				for(unsigned int c = 0; c < bitAmount[s].size(); c++){
					if(bitAmount[s][c] > 0){
						return false;
					}
				}
			}
		}

		//mark solution as done
		solution.markSolutionAsComplete();
		//both tests were successful
		return true;
	}


	unsigned CompressionStrategy::getStageOfArrivalForBit(Bit* bit){
		return 0; //!!! This is a workaround. BitHeap compression currently does not work when bits are distributed to different stages. TODO: fix
		if(!UserInterface::pipelineActive_) {
			return 0;
		}
		//for a bit to be in stage i, its criticalPath delay has to be smaller or equal to boundaries[i].
		//(for the same cycle)
		double objectiveDelay = 1.0/bitheap->getOp()->getTarget()->frequency();
		vector<double> boundaries;    //TODO: Instead of creating this array, a simple division of compressionDelay / objectiveDelay (+/- 1) should do!
		double remainingTime = objectiveDelay;
		while(remainingTime >= compressionDelay){
			remainingTime -= compressionDelay;
			boundaries.insert(boundaries.begin(), remainingTime);
		}

		unsigned int stage = bit->signal->getCycle() * stagesPerCycle;  //maxStage = number of stages - 1
		double tempCriticalPath = bit->signal->getCriticalPath();
		if(tempCriticalPath > boundaries[stagesPerCycle - 1]){
			//bit is in next cycle
			stage += stagesPerCycle;
		}
		else{
			for(unsigned i = 0; i < boundaries.size(); i++){
				if(tempCriticalPath > boundaries[i]){
					stage++;
				}
			}
		}
//      cerr << "objectiveDelay=" << objectiveDelay << ", compressionDelay=" << compressionDelay << ", bit->signal->getCycle()=" << bit->signal->getCycle() << ", stagesPerCycle=" << stagesPerCycle << ", stage=" << stage << endl;
		return stage;
	}



	void CompressionStrategy::orderCompressorsByCompressionEfficiency(){
		vector<double> compressionRatio;
		vector<bool> alreadyChosen;

		for(unsigned int i = 0; i < possibleCompressors.size(); i++){
			double ratio = possibleCompressors[i]->getEfficiency();
			compressionRatio.push_back(ratio);
			alreadyChosen.push_back(false);
		}

		vector<BasicCompressor*> rightOrder;
		for(unsigned int z = 0; z < possibleCompressors.size(); z++){
			unsigned int pos = 0;
			double currentMaxRatio = -1.0;
			for(unsigned int i = 0; i < possibleCompressors.size(); i++){
				if((alreadyChosen[i] == false) && (compressionRatio[i] > currentMaxRatio)){
					currentMaxRatio = compressionRatio[i];
					pos = i;
				}
			}
			rightOrder.push_back(possibleCompressors[pos]);
			alreadyChosen[pos] = true;
			REPORT(DETAILED, "found compressor " << possibleCompressors[pos]->getStringOfIO() << " with compressionRatio of " << currentMaxRatio << ". Original compressor position is " << pos);
		}
		possibleCompressors = rightOrder;
	}



	void CompressionStrategy::startCompression()
	{
		double delay;

		// Add the constant bits
		REPORT(DEBUG, "Adding the constant bits");
		bitheap-> op->vhdl << endl << tab << "-- Adding the constant bits " << endl;
		bool isConstantNonzero = false;
		for (int w = bitheap->lsb; w <= bitheap->msb; w++){
			if (1 == ((bitheap->constantBits>>(w-bitheap->lsb)) & 1) ){
				Bit* bit = bitheap->addBit("'1'", w);
				//set the signal to constant type, with declaration
				bit->signal->setType(Signal::constantWithDeclaration);
				//initialize the signals predecessors and successors
				bit->signal->resetPredecessors();
				bit->signal->resetSuccessors();
				//set the timing for the constant signal, at cycle 0, criticalPath 0, criticalPathContribution 0
				bit->signal->setCycle(0);
				bit->signal->setCriticalPath(0.0);
				bit->signal->setCriticalPathContribution(0.0);
				bit->signal->setHasBeenScheduled(true);

				isConstantNonzero = true;
			}
		}
		//when the constant bits are all zero, report it
		if(!isConstantNonzero){
			REPORT(DEBUG, "All the constant bits are zero, nothing to add");
			bitheap-> op->vhdl << tab << tab << "-- All the constant bits are zero, nothing to add" << endl;
		}
		bitheap-> op->vhdl << endl;


		//take a snapshot of the bitheap, before the start of the compression
		bitheapPlotter->takeSnapshot(getSoonestBit(0, bitheap->width-1), getSoonestBit(0, bitheap->width-1));

		//first, set the delay to the delay of a compressor
		delay = compressionDelay;

		REPORT(DEBUG, "current delay is: delay=" << delay);

		//compress the bitheap to height 2
		//  (the last column can be of height 3)
		if(!bitheap->isCompressed || bitheap->compressionRequired())
		{
			REPORT(DEBUG, "starting compressionAlgorithm");
			compressionAlgorithm();
			REPORT(DEBUG, "compressionAlgorithm ended");

		}else{
			REPORT(DEBUG, "Bitheap already compressed, so startCompression has nothing else to do.");
		}
		//generate the final addition
		generateFinalAddVHDL(bitheap->getOp()->getTarget()->getVendor() == "Xilinx");
		//concatenate all the chunks and create the final result
		concatenateChunks();
		//mark the bitheap as compressed
		bitheap->isCompressed = true;
		//plot the bitheap
		if(bitheap->getOp()->getTarget()->generateFigures()) {
		    bitheapPlotter->plotBitHeap();
            printSolutionTeX();
		}
	}

	bool CompressionStrategy::compress(double delay, Bit *soonestCompressibleBit)
	{
		bool compressionPerformed = false;

		//if there is no soonest bit, the there's nothing else to do
		//  the bit wit the smallest (cycle, critical path), which is also compressible
		if(soonestCompressibleBit == nullptr)
			return false;

		//try to apply the compressors in the decreasing
		//  order of their compression ratio
		for(unsigned i=0; i<possibleCompressors.size(); i++)
		{
			//go through all the lines and try to apply the compressor
			//  to bits that are within the given delay
			for(unsigned j=compressionDoneIndex; j<bitheap->bits.size(); j++)
			{
				vector<Bit*> compressorBitVector = canApplyCompressor(j, i, soonestCompressibleBit, delay);

				while(compressorBitVector.size() > 0)
				{
					Compressor* realCompressor = possibleCompressors[i]->getCompressor();
					applyCompressor(compressorBitVector, realCompressor, bitheap->lsb+j);
					compressorBitVector.clear();
					compressorBitVector = canApplyCompressor(j, i, soonestCompressibleBit, delay);
					compressionPerformed = true;
				}
			}
		}

		return compressionPerformed;
	}



	void CompressionStrategy::applyCompressor(vector<Bit*> bitVector, Compressor* compressor, int weight)
	{
		vector<string> compressorInputs;
		ostringstream inputName, compressorIONames;
		int instanceUID = Operator::getNewUId();
		unsigned count = 0;

		if(bitVector.size() == 0)
			THROWERROR("Bit vector to compress is empty in applyCompressor");
		if(compressor == nullptr)
			THROWERROR("Compressor empty in applyCompressor");

		//build the inputs to the compressor
		for(unsigned i=0; i<compressor->heights.size(); i++)
		{
			inputName.str("");

			if(compressor->heights[i] > 0)
			{
				if(count >= bitVector.size())
					THROWERROR("Bit vector does not containing sufficient bits, "
							<< "as requested by the compressor is applyCompressor.");

				inputName << bitVector[count]->getName();
				count++;
			}
			for(int j=1; j<compressor->heights[i]; j++)
			{
				if(count >= bitVector.size())
					THROWERROR("Bit vector does not containing sufficient bits, "
							<< "as requested by the compressor is applyCompressor.");

				inputName << " & " << bitVector[count]->getName();
				count++;
			}

			compressorInputs.push_back(string(inputName.str()));
		}

		//create the signals for the compressor inputs
		//  and create the port mappings
		bitheap->getOp()->vhdl << endl;
		for(unsigned i=0; i<compressor->heights.size(); i++)
		{
			if(compressor->heights[i] > 0)
			{
				compressorIONames.str("");
				compressorIONames << compressor->getName() << "_bh"
						<< bitheap->guid << "_uid" << instanceUID << "_In" << i;
				bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(compressorIONames.str(), compressor->heights[i])
						<< " <= \"\" & " << compressorInputs[i] << ";" << endl;
				bitheap->getOp()->inPortMap(join("X",i), compressorIONames.str());
			}
		}

		//create the signals for the compressor output
		//  and create the port mapping
		compressorIONames.str("");
		compressorIONames << compressor->getName() << "_bh"
				<< bitheap->guid << "_uid" << instanceUID << "_Out";
		// Following line commented by F2D: outPortMap does the declaration
		// bitheap->getOp()->declare(compressorIONames.str(), compressor->wOut, (compressor->wOut > 1));
		bitheap->getOp()->outPortMap("R", compressorIONames.str());

		//create the compressor instance
		bitheap->getOp()->vhdl << bitheap->getOp()->instance(compressor, join(compressor->getName(), "_uid", instanceUID)) << endl;

		//mark the bits that were at the input of the compressor as having been compressed
		//  so that they can be eliminated from the bitheap
		for(unsigned i=0; i<bitVector.size(); i++)
			bitheap->markBit(bitVector[i], BitType::compressed);

		//mark the compressor as having been used
		compressor->markUsed();

		//add the result of the compression to the bitheap
		bitheap->addSignal(compressorIONames.str(), weight);

		//mark the bits of the result of the compression as just added
		bitheap->markBits(bitheap->getOp()->getSignalByName(compressorIONames.str()), BitType::justAdded, weight);
	}

	void CompressionStrategy::applyCompressor(vector<vector<Bit*> > inputBits, Compressor* compressor, int weight)
	{
		vector<string> compressorInputs;
		ostringstream inputName, vectorName;
		int instanceUID = Operator::getNewUId();


		if(compressor == nullptr){
			THROWERROR("Compressor empty in applyCompressor");
		}

		//inputs for compressor
		for(unsigned int c = 0; c < compressor->heights.size(); c++){
			inputName.str("");
			for(unsigned int i = 0; i < (unsigned) compressor->heights[c]; i++){
				if(i == 0){
					if(inputBits[c].size() > i){
						inputName << inputBits[c][i]->getName();
					}
					else{
						inputName << "\"0\"";
					}
				}
				else{
					if(inputBits[c].size() > i){
						inputName << " & " << inputBits[c][i]->getName();
					}
					else{
						inputName << " & \"0\"";
					}
				}
			}
			compressorInputs.push_back(string(inputName.str()));
		}

		//input signals
		bitheap->getOp()->vhdl << endl;
		for(unsigned i=0; i<compressor->heights.size(); i++)
		{
			if(compressor->heights[i] > 0)
			{
				vectorName.str("");
				vectorName << compressor->getName() << "_bh"
						<< bitheap->guid << "_uid" << instanceUID << "_In" << i;
				bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(vectorName.str(), compressor->heights[i])
						<< " <= \"\" & " << compressorInputs[i] << ";" << endl;
				bitheap->getOp()->inPortMap(join("X",i), vectorName.str());
			}
		}

		//create the signals for the compressor outputs. Can handle compressors who have
		//more than one outputbit per column. If that is the case, its assumed that the outputvector
		//for e.g. the second bits spans from lsb of the compressor to the closest occurrence of the in this
		//example second bit. Only the bits which are not constantly zero will be added to the bitheap.
		unsigned int tempHeight = 0;
		while(true){
			tempHeight++;
			//does this height occur in the compressor
			int lastOccurence = -1;
			for(unsigned int c = 0; c < compressor->outHeights.size(); c++){
				if((unsigned) compressor->outHeights[c] >= tempHeight){
					lastOccurence = c;
				}
			}
			if(lastOccurence == -1){
				break;  //not found
			}


			vectorName.str("");
			vectorName << compressor->getName() << "_bh" << bitheap->guid << "_uid"
					<< instanceUID << "_Out" << (tempHeight - 1);
			//          bitheap->getOp()->declare(vectorName.str(), lastOccurence + 1, (lastOccurence + 1 > 1));

			//check if there are more than one outputVectors in total
			bool singleOutputVector = true;
			for(unsigned int c = 0; c < compressor->outHeights.size(); c++){
				if(compressor->outHeights[c] > 1){
					singleOutputVector = false;
					break;
				}
			}
			if(singleOutputVector == true){
				bitheap->getOp()->outPortMap("R", vectorName.str());
			}
			else{
				ostringstream tempR;
				tempR.str("");
				tempR << "R" << (tempHeight - 1);
				bitheap->getOp()->outPortMap(tempR.str(), vectorName.str());
			}

			//create the compressor
			bitheap->getOp()->vhdl << bitheap->getOp()->instance(compressor, join(compressor->getName(), "_uid", instanceUID)) << endl;

			//add the outputBits to the bitheap
			for(unsigned int c = 0; c <= (unsigned) lastOccurence; c++){
				if(c + weight < bitheap->width && (unsigned)compressor->outHeights[c] >= tempHeight){
					bitheap->addBit(vectorName.str() + of(c), c + weight + bitheap->lsb);
				}
				bitheap->markBits(bitheap->getOp()->getSignalByName(vectorName.str()), BitType::justAdded, weight + bitheap->lsb);
			}

		}

		//set the inputbits to compressed
		for(unsigned int c = 0; c < inputBits.size(); c++){
			for(unsigned int i = 0; i < inputBits[c].size(); i++){
				bitheap->markBit(inputBits[c][i], BitType::compressed);
			}
		}
		compressor->markUsed();
	}


	void CompressionStrategy::generateFinalAddVHDL(bool isXilinx)
	{
		//check if the last two lines are already compressed
		if(bitheap->getMaxHeight() < 2)
		{
			//an adder isn't necessary; concatenateLSBColumns should do the rest
			//concatenateLSBColumns(); //it is already called in
		}else if(bitheap->getMaxHeight() == 2  || heightIs2exceptCin()){
			ostringstream adderIn0, adderIn0Name, adderIn1, adderIn1Name, adderOutName, adderCin, adderCinName;
			int adderStartIndex = compressionDoneIndex + 1;

			//create the names for the inputs/output of the adder
			adderIn0Name << "bitheapFinalAdd_bh" << bitheap->guid << "_In0";
			adderIn1Name << "bitheapFinalAdd_bh" << bitheap->guid << "_In1";
			adderCinName << "bitheapFinalAdd_bh" << bitheap->guid << "_Cin";
			adderOutName << "bitheapFinalAdd_bh" << bitheap->guid << "_Out";

			//add the first bits to the adder inputs
			if(compressionDoneIndex <= bitheap->width)
			{
				if(bitheap->bits[bitheap->width-1].size() > 1)
				{
					adderIn0 << bitheap->bits[bitheap->width-1][0]->getName();
					adderIn1 << bitheap->bits[bitheap->width-1][1]->getName();
				}else if(bitheap->bits[bitheap->width-1].size() > 0)
				{
					adderIn0 << bitheap->bits[bitheap->width-1][0]->getName();
					adderIn1 << "\"0\"";
				}else{
					adderIn0 << "\"0\"";
					adderIn1 << "\"0\"";
				}
			}
			//add the rest of the terms to the adder inputs
			for(int i=bitheap->width-2; i>(int)compressionDoneIndex; i--)
			{
				if(bitheap->bits[i].size() > 1)
				{
					adderIn0 << " & " << bitheap->bits[i][0]->getName();
					adderIn1 << " & " << bitheap->bits[i][1]->getName();
				}else if(bitheap->bits[i].size() > 0)
				{
					adderIn0 << " & " << bitheap->bits[i][0]->getName();
					adderIn1 << " & " << "\"0\"";
				}else{
					adderIn0 << " & " << "\"0\"";
					adderIn1 << " & " << "\"0\"";
				}
			}
			//add the last column to the adder, if compression is still required on that column
			if((compressionDoneIndex == 0) && (bitheap->bits[compressionDoneIndex].size() > 1))
			{
				adderIn0 << " & " << bitheap->bits[compressionDoneIndex][0]->getName();
				adderIn1 << " & " << bitheap->bits[compressionDoneIndex][1]->getName();

				if(bitheap->bits[compressionDoneIndex].size() > 2)
					adderCin << bitheap->bits[compressionDoneIndex][2]->getName();
				else
					adderCin << "\'0\'";

				adderStartIndex--;
			}else{
				if(bitheap->bits[compressionDoneIndex+1].size() > 2)
					adderCin << bitheap->bits[compressionDoneIndex+1][2]->getName();
				else
					adderCin << "\'0\'";
			}

			//create the signals for the inputs/output of the adder
			bitheap->getOp()->vhdl << endl;
			bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderIn0Name.str(), bitheap->msb-adderStartIndex+1+1-bitheap->lsb)
					<< " <= \"0\" & " << adderIn0.str() << ";" << endl;
			bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderIn1Name.str(), bitheap->msb-adderStartIndex+1+1-bitheap->lsb)
					<< " <= \"0\" & " << adderIn1.str() << ";" << endl;
			bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderCinName.str())
					<< " <= " << adderCin.str() << ";" << endl;
			bitheap->getOp()->vhdl << endl;

			//declare the adder
#if 0
			IntAdder* adder;

			//create the adder
			bitheap->getOp()->declare(adderOutName.str(), bitheap->msb-adderStartIndex+1+1);
			adder = new IntAdder(bitheap->op, bitheap->getOp()->getTarget(), bitheap->msb-adderStartIndex+1+1);

			//create the port maps for the adder
			bitheap->getOp()->inPortMap("X", adderIn0Name.str());
			bitheap->getOp()->inPortMap("Y", adderIn1Name.str());
			bitheap->getOp()->inPortMap("Cin", adderCinName.str());
			bitheap->getOp()->outPortMap("R", adderOutName.str());


			//create the instance of the adder
			bitheap->getOp()->vhdl << bitheap->getOp()->instance(adder, join("bitheapFinalAdd_bh", bitheap->guid)) << endl;
#else
			bitheap->getOp()->newInstance("IntAdder",
																					"bitheapFinalAdd_bh"+to_string(bitheap->guid),
																					"wIn=" + to_string(bitheap->msb-adderStartIndex+1+1-bitheap->lsb),
																					"X=>"+ adderIn0Name.str()
																					+ ",Y=>"+adderIn1Name.str()
																					+ ",Cin=>" + adderCinName.str(),
																					"R=>"+ adderOutName.str()   );
#endif
			//add the result of the final add as the last chunk
			chunksDone.push_back(join(adderOutName.str(), range(bitheap->msb-adderStartIndex-bitheap->lsb, 0)));

		} else if(bitheap->getMaxHeight() == 3){
            ostringstream adderIn0, adderIn0Name, adderIn1, adderIn1Name, adderIn2, adderIn2Name, adderOutName, adderCin, adderCinName;
            int adderStartIndex = compressionDoneIndex + 1;

            //create the names for the inputs/output of the adder
            adderIn0Name << "bitheapFinalAdd_bh" << bitheap->guid << "_In0";
            adderIn1Name << "bitheapFinalAdd_bh" << bitheap->guid << "_In1";
            adderIn2Name << "bitheapFinalAdd_bh" << bitheap->guid << "_In2";
            adderCinName << "bitheapFinalAdd_bh" << bitheap->guid << "_Cin";
            adderOutName << "bitheapFinalAdd_bh" << bitheap->guid << "_Out";

            //add the first bits to the adder inputs
            if(compressionDoneIndex <= bitheap->width)
            {
                if(bitheap->bits[bitheap->width-1].size() > 2)
                {
                    adderIn0 << bitheap->bits[bitheap->width-1][0]->getName();
                    adderIn1 << bitheap->bits[bitheap->width-1][1]->getName();
                    adderIn2 << bitheap->bits[bitheap->width-1][2]->getName();
                }else if(bitheap->bits[bitheap->width-1].size() == 2)
                {
                    adderIn0 << bitheap->bits[bitheap->width-1][0]->getName();
                    adderIn1 << bitheap->bits[bitheap->width-1][1]->getName();
                    adderIn2 << "\"0\"";
                }else if(bitheap->bits[bitheap->width-1].size() == 1)
                {
                    adderIn0 << bitheap->bits[bitheap->width-1][0]->getName();
                    adderIn1 << "\"0\"";
                    adderIn2 << "\"0\"";
                }else{
                    adderIn0 << "\"0\"";
                    adderIn1 << "\"0\"";
                    adderIn2 << "\"0\"";
                }
            }
            //add the rest of the terms to the adder inputs
            for(int i=bitheap->width-2; i>(int)compressionDoneIndex; i--)
            {
                if(bitheap->bits[i].size() > 2)
                {
                    adderIn0 << " & " << bitheap->bits[i][0]->getName();
                    adderIn1 << " & " << bitheap->bits[i][1]->getName();
                    adderIn2 << " & " << bitheap->bits[i][2]->getName();
                }else if(bitheap->bits[i].size() == 2)
                {
                    adderIn0 << " & " << bitheap->bits[i][0]->getName();
                    adderIn1 << " & " << bitheap->bits[i][1]->getName();
                    adderIn2 << " & " << "\"0\"";
                }else if(bitheap->bits[i].size() == 1)
                {
                    adderIn0 << " & " << bitheap->bits[i][0]->getName();
                    adderIn1 << " & " << "\"0\"";
                    adderIn2 << " & " << "\"0\"";
                }else{
                    adderIn0 << " & " << "\"0\"";
                    adderIn1 << " & " << "\"0\"";
                    adderIn2 << " & " << "\"0\"";
                }
            }
            //add the last column to the adder, if compression is still required on that column
            if((compressionDoneIndex == 0) && (bitheap->bits[compressionDoneIndex].size() > 1))
            {
                adderIn0 << " & " << bitheap->bits[compressionDoneIndex][0]->getName();
                adderIn1 << " & " << bitheap->bits[compressionDoneIndex][1]->getName();
                if(bitheap->bits[compressionDoneIndex].size() > 2)
                    adderIn2 << " & " << bitheap->bits[compressionDoneIndex][2]->getName();
                else
                    adderIn2 << " & \'0\'";

                if(bitheap->bits[compressionDoneIndex].size() > 3)
                    adderCin << bitheap->bits[compressionDoneIndex][3]->getName();
                else
                    adderCin << "\'0\'";

                adderStartIndex--;
            }else{
                if(bitheap->bits[compressionDoneIndex+1].size() > 2)
                    adderCin << bitheap->bits[compressionDoneIndex+1][2]->getName();
                else
                    adderCin << "\'0\'";
            }

            //create the signals for the inputs/output of the adder
            bitheap->getOp()->vhdl << endl;
            bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderIn0Name.str(), bitheap->msb-adderStartIndex+1+1)
                                   << " <= \"0\" & " << adderIn0.str() << ";" << endl;
            bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderIn1Name.str(), bitheap->msb-adderStartIndex+1+1)
                                   << " <= \"0\" & " << adderIn1.str() << ";" << endl;
            bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderIn2Name.str(), bitheap->msb-adderStartIndex+1+1)
                                   << " <= \"0\" & " << adderIn2.str() << ";" << endl;
            bitheap->getOp()->vhdl << tab << bitheap->getOp()->declare(adderCinName.str())
                                   << " <= " << adderCin.str() << ";" << endl;
            bitheap->getOp()->vhdl << endl;

            //declare the adder
            bitheap->getOp()->newInstance("XilinxTernaryAddSub",
                                          "bitheapFinalAdd_bh"+to_string(bitheap->guid),
                                          "wIn=" + to_string(bitheap->msb-adderStartIndex+1+1),
                                          "x_i=>"+ adderIn0Name.str()
                                          + ",y_i=>"+adderIn1Name.str()
                                          + ",z_i=>"+adderIn2Name.str(),
                                          //+ ",carry_in=>" + adderCinName.str(),
                                          "sum_o=>"+ adderOutName.str()   );

            //add the result of the final add as the last chunk
			chunksDone.push_back(join(adderOutName.str(), range(bitheap->msb-adderStartIndex-bitheap->lsb, 0)));
		}
	}


	Bit* CompressionStrategy::getSoonestBit(unsigned lsbColumn, unsigned msbColumn)
	{
		Bit* soonestBit = nullptr;
		unsigned count = lsbColumn;

		if(msbColumn > bitheap->width - 1)
			THROWERROR("Invalid arguments for getSoonest bit: lsbColumn="
					<< lsbColumn << " msbColumn=" << msbColumn);
		if(bitheap->getMaxHeight() == 0)
			REPORT(DEBUG, "Warning: trying to obtain the soonest bit from an empty bitheap!");

		//determine the first non-empty bit column
		while((count < msbColumn) && (bitheap->bits[count].size() == 0))
			count++;
		//select the first bit from the column
		soonestBit = bitheap->bits[count][0];

		//determine the soonest bit
		for(unsigned i=count; i<=msbColumn; i++)
			for(unsigned j=0; j<bitheap->bits[i].size(); j++)
			{
				Bit* currentBit = bitheap->bits[i][j];

				//only consider bits that are available for compression
				if(currentBit->type != BitType::free)
					continue;

				if((soonestBit->signal->getCycle() > currentBit->signal->getCycle()) ||
						((soonestBit->signal->getCycle() == currentBit->signal->getCycle())
								&& (soonestBit->signal->getCriticalPath() > currentBit->signal->getCriticalPath())))
				{
					soonestBit = currentBit;
				}
			}

		return soonestBit;
	}


	Bit* CompressionStrategy::getSoonestCompressibleBit(unsigned lsbColumn, unsigned msbColumn, double delay)
	{
		Bit *soonestBit = nullptr, *soonestCompressibleBit = nullptr;
		unsigned count = lsbColumn;
		vector<Bit*> appliedCompressor;

		if(msbColumn > bitheap->width - 1)
			THROWERROR("Invalid arguments for getSoonestCompressibleBit bit: lsbColumn="
					<< lsbColumn << " msbColumn=" << msbColumn);
		if(bitheap->getMaxHeight() == 0)
			REPORT(DEBUG, "Warning: trying to obtain the soonest bit from an empty bitheap!");

		//determine the first non-empty bit column
		while((count < msbColumn) && (bitheap->bits[count].size() == 0))
			count++;

		//search in each column for the soonest bit
		for(unsigned i=count; i<=msbColumn; i++)
		{
			//clear the possible content of the compressor
			appliedCompressor.clear();

			//initialize the soonest bit
			if(bitheap->bits[i].size())
			{
				soonestBit = bitheap->bits[i][0];
			}else
			{
				//empty column, so no compressor can be applied
				continue;
			}

			//try to apply a compressor to the current column
			for(unsigned j=0; j<possibleCompressors.size(); j++)
			{
				appliedCompressor = canApplyCompressor(i, j, soonestBit, delay);
				if(appliedCompressor.size() > 0)
					break;
			}

			//if a compressor can be applied on this column, then check
			//  if the current soonest bit is earlier
			if(soonestCompressibleBit == nullptr)
			{
				soonestCompressibleBit = soonestBit;
				continue;
			}
			if(appliedCompressor.size() > 0)
				if((soonestBit->signal->getCycle() > soonestCompressibleBit->signal->getCycle()) ||
						((soonestBit->signal->getCycle() == soonestCompressibleBit->signal->getCycle())
								&& (soonestBit->signal->getCriticalPath() > soonestCompressibleBit->signal->getCriticalPath())))
				{
					soonestCompressibleBit = soonestBit;
				}
		}

		return soonestCompressibleBit;
	}


	vector<Bit*> CompressionStrategy::canApplyCompressor(unsigned columnNumber, unsigned compressorNumber, Bit* soonestBit, double delay)
	{
		BasicCompressor *compressor;
		vector<Bit*> returnValue;
		double soonestBitTotalTime = 0.0;

		if(compressorNumber > possibleCompressors.size())
		{
			THROWERROR("Compressor index out of range: compressorNumber="
					<< compressorNumber << " maximum index=" << possibleCompressors.size()-1);
		}else{
			compressor = possibleCompressors[compressorNumber];
			soonestBitTotalTime = soonestBit->signal->getCycle() * (1.0/bitheap->getOp()->getTarget()->frequency()) +
									soonestBit->signal->getCriticalPath();
		}

		//check there are sufficient columns of bits
		returnValue.clear();
		if(compressor->heights.size()+columnNumber-1 >= bitheap->bits.size())
			return returnValue;

		for(unsigned i=0; i<compressor->heights.size(); i++)
		{
			int bitCount = 0;
			unsigned columnIndex = 0;

			while((bitCount < compressor->heights[i])
					&& (columnIndex < bitheap->bits[columnNumber+i].size()))
			{
				Bit *currentBit = bitheap->bits[columnNumber+i][columnIndex];
				double currentBitTotalTime = 0.0;

				//only consider bits that are free
				if(currentBit->type != BitType::free)
				{
					//pass to the next bit
					columnIndex++;
					continue;
				}

				//the total time of the current bit
				currentBitTotalTime = currentBit->signal->getCycle() * 1.0/bitheap->getOp()->getTarget()->frequency() +
										currentBit->signal->getCriticalPath();
				//check if the bit can be compressed
				if(currentBitTotalTime-soonestBitTotalTime <= delay)
				{
					returnValue.push_back(currentBit);
					bitCount++;
				}
				//pass to the next bit in the column
				columnIndex++;
			}
			//if the current column does not contain sufficient bits,
			//  then the compressor cannot be applied
			if(bitCount < compressor->heights[i])
			{
				returnValue.clear();
				return returnValue;
			}
		}

		return returnValue;
	}


	void CompressionStrategy::generatePossibleCompressors()
	{
		int col0, col1;

#if 0 // The generic alternative
		//Generate all "workable" compressors for 2 columns
		for(col0=maxCompressibleBits; col0>=3; col0--)
			for(col1=col0; col1>=0; col1--)
				if((col0 + col1<=maxCompressibleBits) && (intlog2(col0 + 2*col1)<=3))
					{
						vector<int> newVect;

						REPORT(DEBUG, "Generating compressor for col1=" << col1 <<", col0=" << col0);

						newVect.push_back(col0);
						newVect.push_back(col1);
						Compressor* bc = new Compressor(op->getTarget(), newVect);
						possibleCompressors.push_back(bc);
					}

#endif

		//generate the compressors, in decreasing order
		//  of their compression ratio
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=6;

			newVect.push_back(col0);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 3.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=4;
			col1=1;

			newVect.push_back(col0);
			newVect.push_back(col1);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 2.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=5;

			newVect.push_back(col0);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 2.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=3;
			col1=1;

			newVect.push_back(col0);
			newVect.push_back(col1);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 2.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=4;

			newVect.push_back(col0);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 2.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=3;
			col1=2;

			newVect.push_back(col0);
			newVect.push_back(col1);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 2.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		{
			vector<int> newVect;
			BasicCompressor *newCompressor;

			col0=3;

			newVect.push_back(col0);
			newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), newVect, 1.0, "combinatorial", true);
			//newCompressor->setShared();
			possibleCompressors.push_back(newCompressor);
		}
		if(bitheap->getOp()->getTarget()->useTargetOptimizations() && (bitheap->getOp()->getTarget()->getVendor() == "Xilinx"))
		{
			REPORT(DEBUG,"Adding target optimized GPCs for Xilinx FPGAs");

//          possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {6,0,6})); //outdated since (6,0,7;5) GPC
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {7,0,6}));
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {5,1,6}));
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {3,2,6}));
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {5,2,3,1}));
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {5,1,4,1}));
//            possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {6,0,4,1})); //outdated since (1,4,0,7;5) GPC
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {7,0,4,1}));
			possibleCompressors.push_back(new BasicXilinxGPC(bitheap->getOp(), bitheap->getOp()->getTarget(), {7,1,1,2}));
		}
	}


	void CompressionStrategy::concatenateLSBColumns()
	{
		unsigned count = compressionDoneIndex;

		//check up until which index has the bitheap already been compressed
		while((count < bitheap->bits.size()) && (bitheap->bits[count].size() <= 1))
			count++;
		//all the columns might have been already compressed
		if(count == bitheap->bits.size())
			count--;
		//if the current column needs compression, then do not consider it
		if((bitheap->bits[count].size() > 1) && (count > 0))
			count--;

		if(count == compressionDoneIndex)
		{
			if((!bitheap->compressionRequired()) && (count == 0) && (bitheap->bits[count].size() <= 1))
			{
				//the bitheap is compressed, all that is left
				//  is to add the lsb bit to the result
				bitheap->getOp()->vhdl << tab
						<< bitheap->getOp()->declare(join("tmp_bitheapResult_bh", bitheap->guid, "_", count),1,true);
				if(bitheap->bits[count].size() > 0)
					bitheap->getOp()->vhdl << "(0) <= " << bitheap->bits[count][0]->getName() << ";" << endl;
				else
					bitheap->getOp()->vhdl << " <= \'0\';" << endl;

				chunksDone.insert(chunksDone.begin(), join("tmp_bitheapResult_bh", bitheap->guid, "_", count));
			}else{
				//no new columns that have been compressed
				return;
			}
		}else if(count > compressionDoneIndex){
			//new columns have been compressed

			//create the signals for a new chunk
			if(count-compressionDoneIndex+(compressionDoneIndex == 0 ? 1 : 0) > 1)
			{
				bitheap->getOp()->vhdl << tab
					<< bitheap->getOp()->declare(join("tmp_bitheapResult_bh", bitheap->guid, "_", count),
							count-compressionDoneIndex+(compressionDoneIndex == 0 ? 1 : 0), true) << " <= " ;
			}
			else
			{
				bitheap->getOp()->vhdl << tab
					<< bitheap->getOp()->declare(join("tmp_bitheapResult_bh", bitheap->guid, "_", count)) << " <= " ;
			}

			//add the bits to the chunk
			//  add the first bit
			if(bitheap->bits[count].size() > 0)
				bitheap->getOp()->vhdl << bitheap->bits[count][0]->getName();
			else
				bitheap->getOp()->vhdl << "\"0\"";
			//  add the rest of the bits
			for(unsigned i=count-1; i>compressionDoneIndex; i--)
			{
				if(bitheap->bits[i].size() > 0)
					bitheap->getOp()->vhdl << " & " << bitheap->bits[i][0]->getName();
				else
					bitheap->getOp()->vhdl << " & \"0\"";
			}
			//add one last bit, if necessary
			if(compressionDoneIndex == 0)
			{
				if(bitheap->bits[0].size() > 0)
					bitheap->getOp()->vhdl << " & " << bitheap->bits[0][0]->getName();
				else
					bitheap->getOp()->vhdl << " & \"0\"";
			}
			bitheap->getOp()->vhdl << ";" << endl;

			//add the new chunk to the vector of compressed chunks
			chunksDone.push_back(join("tmp_bitheapResult_bh", bitheap->guid, "_", count));

			//advance the compression done index
			compressionDoneIndex = count;
		}
	}


	void CompressionStrategy::concatenateChunks()
	{
		if(chunksDone.size() == 0)
			THROWERROR("No chunks to concatenate in concatenateChunksDone");

		//create the signal containing the result of the bitheap compression
		bitheap->getOp()->vhdl << tab
				<< bitheap->getOp()->declare(join("bitheapResult_bh", bitheap->guid), bitheap->width) << " <= ";
		// the chunks in reverse order, from msb to lsb
		bitheap->getOp()->vhdl << chunksDone[chunksDone.size()-1];
		for(int i=(int)chunksDone.size()-2; i>=0; i--)
			bitheap->getOp()->vhdl << " & " << chunksDone[i];
		bitheap->getOp()->vhdl << ";" << endl;
	}

    int CompressionStrategy::heightIs2exceptCin() {
        for(int i=bitheap->width-1; i>(int)compressionDoneIndex; i--)
        {
            if(bitheap->bits[i].size() > 2)
                return 0;                       //Dual input adder can only handle 2 inputs
        }
        if(bitheap->bits[compressionDoneIndex].size() > 3)
            return 0;                           //But the adder LSB can process 3 bits due to the Cin, but not more
        return 1;
    }

    void CompressionStrategy::printSolutionTeX(){
        stringstream compr;
        compr <<
              tab << "\\documentclass[margin=0pt]{standalone}" << endl <<
              tab << "\\usepackage{graphicx}" << endl <<
              tab << "\\usepackage{xcolor,tikz,tikzscale}" << endl <<
              tab << "\\begin{document}" << endl <<
              tab << "\\begin{tikzpicture}" << endl;

        #define N_TIKZ_COLORS 14
        #define MAX_STRING_SIZE 20
        char tikz_colors[N_TIKZ_COLORS][MAX_STRING_SIZE] = { "red", "green", "blue", "cyan", "magenta", "yellow", "brown", "lime", "olive", "orange", "pink", "purple", "teal", "violet"};
        int color = 0;

        vector<vector<int> > bh_bit(solution.getNumberOfStages()+1, vector<int>(bitheap->width+4, 0));
        vector<vector<int> > coveredBits(solution.getNumberOfStages()+1, vector<int>(bitheap->width+4, 0));
        vector<int> bh_constants(bitheap->width+4, 0);

        int lastMaxHeight = 0, previousMaxHeight;
        for(unsigned int c = 0; c < orderedBits[0].size(); c++) {
            bh_bit[0][c] = orderedBits[0][c].size();
            lastMaxHeight = (lastMaxHeight < orderedBits[0][c].size())?orderedBits[0][c].size():lastMaxHeight;
            //cout << "col=" << c << " height=" << orderedBits[0][c].size() << endl;
            for(unsigned int b = 0; b < orderedBits[0][c].size(); b++) {
                if(orderedBits[0][c][b]->signal->type() == Signal::SignalType::constantWithDeclaration || orderedBits[0][c][b]->signal->type() == Signal::SignalType::constant) bh_constants[c]++;
                //cout << (orderedBits[0][c][b]->signal->getName()) << " is constant " << (orderedBits[0][c][b]->signal->type() == Signal::SignalType::constantWithDeclaration) << endl;
            }
        }

        int shift = 0;
        int finalAdderStart = 0, afound = 0;
        for(unsigned int s = 0; s <= solution.getNumberOfStages(); s++) {
            int maxHeight = 0;
            for (unsigned int c = 0; c < bitheap->width; c++) {

                if(s < solution.getNumberOfStages()){
                    vector<pair<BasicCompressor *, unsigned int> > tempVector;
                    tempVector = solution.getCompressorsAtPosition(s, c);
                    for (unsigned int j = 0; j < tempVector.size(); j++) {      //for every compressor in current column and stage
                        int cc;
                        //cout << tempVector[j].first->heights.size() << " " << tempVector[j].first->heights[0] << " " << tempVector[j].first->getStringOfIO() << endl;
                        compr << tab <<  "\\draw[fill="<< tikz_colors[color++%N_TIKZ_COLORS] << "!50,opacity=0.75] ("<<-(int)c*0.25L+0.125L<<","<<(coveredBits[s][c]*0.25L)-0.125L-shift*0.25L<<")  {[rounded corners=2.5pt]";
                        for(cc = 0; cc < tempVector[j].first->heights.size(); cc++){
                            compr << "--("<<-(int)(c+cc)*0.25+0.125<<","<<(coveredBits[s][c+cc]+tempVector[j].first->heights[cc])*0.25L-0.125L-shift*0.25L<<")";
                            if(cc+1 == tempVector[j].first->heights.size() || cc+1 < tempVector[j].first->heights.size() && (coveredBits[s][c+cc]+tempVector[j].first->heights[cc]) != coveredBits[s][c+cc+1]+tempVector[j].first->heights[cc+1])
                                compr << "--("<<-(int)(c+cc+1)*0.25+0.125<<","<<(coveredBits[s][c+cc]+tempVector[j].first->heights[cc])*0.25L-0.125L-shift*0.25L<<")";
                            if(cc+1 == tempVector[j].first->heights.size()){
                                compr << "--("<<-(int)(c+cc+1)*0.25+0.125<<","<<coveredBits[s][c+cc]*0.25L-0.125L-shift*0.25L<<")";
                            }
                            coveredBits[s][c+cc] += tempVector[j].first->heights[cc];
                        }
                        for(cc--; 0 <= cc; cc--){
                            compr << "--("<<-(int)(c+cc)*0.25+0.125<<","<<(((int)coveredBits[s][c+cc]-(int)tempVector[j].first->heights[cc]))*0.25L-0.125L-shift*0.25L<<")";
                            if(0 <= cc-1 && ((int)coveredBits[s][c+cc]-(int)tempVector[j].first->heights[cc]) != ((int)coveredBits[s][c+cc-1]-(int)tempVector[j].first->heights[cc-1])){
                                compr << "--("<<-(int)(c+cc)*0.25+0.125<<","<<(((int)coveredBits[s][c+cc-1]-(int)tempVector[j].first->heights[cc-1]))*0.25L-0.125L-shift*0.25L<<")";
                            }
                        }
                        compr << "--cycle};" << endl;
                        int prevBitPos;
                        for(cc = 0; cc < tempVector[j].first->outHeights.size(); cc++){
                            if(0 < cc){ //draw connection for output bits of compressors
                                compr << tab << "\\draw[black,thick] (-0.25*" << c+cc << "," << bh_bit[s+1][c+cc]*0.25L-shift*0.25L-(lastMaxHeight)*0.25L << ")"
                                << "--(-0.25*" << c+cc-1 << "," << prevBitPos*0.25L-shift*0.25L-(lastMaxHeight)*0.25L << ") ;" << endl;
                            }
                            prevBitPos = bh_bit[s+1][c+cc];
                            bh_bit[s+1][c+cc] += tempVector[j].first->outHeights[cc];
                            maxHeight = (maxHeight<bh_bit[s+1][c+cc])?bh_bit[s+1][c+cc]:maxHeight;
                        }
                    }
                    //copy bits that are not compressed in current stage to subsequent stage
                    bh_bit[s + 1][c] += (0 < bh_bit[s][c] - coveredBits[s][c]) ? bh_bit[s][c] - coveredBits[s][c] : 0;
                    maxHeight = (maxHeight<bh_bit[s+1][c])?bh_bit[s+1][c]:maxHeight;
                } else {    //draw final adder
                    if(1 < bh_bit[s][c] && !afound){
                        finalAdderStart = c; afound = 1;
                        compr << tab <<  "\\draw[fill=blue!50,opacity=0.75] ("<<-(int)finalAdderStart*0.25L+0.125L<<","<<-0.125L-(shift)*0.25L<<")  {[rounded corners=2.5pt]";
                        compr << "--("<<-(int)bitheap->width*0.25L+0.125L<<","<<-0.125L-(shift)*0.25L<<")";
                        compr << "--("<<-(int)bitheap->width*0.25L+0.125L<<","<<+0.25*bitheap->final_add_height-0.125L-(shift)*0.25L<<")";
                        compr << "--("<<-(int)finalAdderStart*0.25L+0.125L<<","<<+0.25*bitheap->final_add_height-0.125L-(shift)*0.25L<<")--cycle};";
                    }
                }

                if(bh_bit[s][c]){       //bits in column of current stage
                    compr << tab <<  tab <<"\\foreach \\i in {0,...," << bh_bit[s][c]-1 << "}" << endl <<
                          tab << tab << tab <<"\\draw[fill=blue!20] (-0.25*" << c << ",0.25*\\i+" << -shift*0.25L << ") circle (2pt);" << endl;
                    if(s == 0 && bh_constants[c] != 0) for(int b=0; b < bh_constants[c]; b++){
                        compr << tab <<  tab <<"\\node[] at (-0.25*" << c << ",0+" << 0.25*b-shift*0.25L << ") {\\textsf{\\tiny c}};" << endl;
                    }
                }
            }
            previousMaxHeight = lastMaxHeight;
            shift += lastMaxHeight;
            lastMaxHeight = maxHeight;

            if(s == solution.getNumberOfStages()) { //Output of final stage of the compressor tree
                compr << tab << "\\draw[black,thick] (-0.25*" << finalAdderStart << "," << -shift*0.25L << ")"
                      << "--(-0.25*" << bitheap->width-1 << "," << -shift*0.25L << ") ;" << endl;
                compr << tab <<  tab <<"\\foreach \\i in {0,...," << bitheap->width-1 << "}" << endl <<
                         tab <<  tab << tab <<"\\draw[fill=blue!20] (-0.25*\\i," << -shift*0.25L << ") circle (2pt);" << endl;
            }

        }

        compr << tab << "\\end{tikzpicture}" << endl <<
        tab << "\\end{document}" << endl << endl;

        ofstream result_file;
        result_file.open("compression.tex");
        result_file << compr.str();
        result_file.close();
    }


}
