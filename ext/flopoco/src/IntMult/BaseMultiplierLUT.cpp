#include "BaseMultiplierLUT.hpp"
#include "IntMultiplierLUT.hpp"

namespace flopoco {

Operator* BaseMultiplierLUT::generateOperator(
		Operator *parentOp,
		Target* target,
		Parametrization const & parameters) const
{
	return new IntMultiplierLUT(
			parentOp,
			target,
			parameters.getMultXWordSize(),
			parameters.getMultYWordSize(),
			parameters.isSignedMultX(),
			parameters.isSignedMultY(),
			parameters.isFlippedXY()
		);
}

	double BaseMultiplierLUT::getLUTCost(int x_anchor, int y_anchor, int wMultX, int wMultY, bool signedIO){
        bool signedX = signedIO && (wMultX == x_anchor+wX());
        bool signedY = signedIO && (wMultY == y_anchor+wY());
		//int ws = (wX()==1)?wY():((wY()==1)?wX():wX()+wY());
		int ws = (wX()==1 && !signedX)?wY():((wY()==1 && !signedY)?wX():wX()+wY());     //The wOut of 2x1 or 2x1 is 1 bit larger, when the short or both sides are signed.
		int luts = ((ws <= 5)?ws/2+ws%2:ws);

		int x_min = ((x_anchor < 0)?0: x_anchor);
		int y_min = ((y_anchor < 0)?0: y_anchor);
		int lsb = x_min + y_min;

		int x_max = ((wMultX < x_anchor + (int)wX())?wMultX: x_anchor + (int)wX());
		int y_max = ((wMultY < y_anchor + (int)wY())?wMultY: y_anchor + (int)wY());
		int msb = x_max+y_max;
		ws = (x_max-x_min==1)?y_max-y_min:((y_max-y_min==1)?x_max-x_min:msb - lsb);

		return luts + ws*getBitHeapCompressionCostperBit();
	}

	int BaseMultiplierLUT::ownLUTCost(int x_anchor, int y_anchor, int wMultX, int wMultY, bool signedIO) {
        bool signedX = signedIO && (wMultX == x_anchor+wX());
        bool signedY = signedIO && (wMultY == y_anchor+wY());
        //int ws = (wX()==1)?wY():((wY()==1)?wX():wX()+wY());
        int ws = (wX()==1 && !signedX)?wY():((wY()==1 && !signedY)?wX():wX()+wY());
		int luts = ((ws <= 5)?ws/2+ws%2:ws);
		return luts;
	}


}   //end namespace flopoco

