#include "BaseMultiplierXilinx2xk.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CARRY4.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h"

namespace flopoco {

Operator* BaseMultiplierXilinx2xk::generateOperator(
		Operator *parentOp,
		Target* target,
		Parametrization const & parameters) const
{
	return new BaseMultiplierXilinx2xkOp(
			parentOp,
			target,
            parameters.isSignedMultX(),
			parameters.isSignedMultY(),
            parameters.getMultXWordSize(),
            parameters.getMultYWordSize()
		);
}

double BaseMultiplierXilinx2xk::getLUTCost(int x_anchor, int y_anchor, int wMultX, int wMultY, bool signedIO){
    int luts = ((wX() < wY())?wY():wX()) + 1;

    int x_min = ((x_anchor < 0)?0: x_anchor);
    int y_min = ((y_anchor < 0)?0: y_anchor);
    int lsb = x_min + y_min;

    int x_max = ((wMultX < x_anchor + (int)wX())?wMultX: x_anchor + wX());
    int y_max = ((wMultY < y_anchor + (int)wY())?wMultY: y_anchor + wY());
    int msb = (x_max==1)?y_max:((y_max==1)?x_max:x_max+y_max);

    if(signedIO && ((wMultX-x_anchor-(int)wX())== 0 || (wMultY-y_anchor-(int)wY())== 0)){
        luts++;    //The 2xk-multiplier needs one additional LUT in the signed case
    }

    return luts + (msb - lsb)*getBitHeapCompressionCostperBit();
}

int BaseMultiplierXilinx2xk::ownLUTCost(int x_anchor, int y_anchor, int wMultX, int wMultY, bool signedIO) {
    int luts = ((wX() < wY())?wY():wX()) + 1;
    if(signedIO && ((wMultX-x_anchor-(int)wX())== 0 || (wMultY-y_anchor-(int)wY())== 0)){
        luts++;    //The 2xk-multiplier needs one additional LUT in the signed case
    }
    return luts;
}

OperatorPtr BaseMultiplierXilinx2xk::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
{
    int wX, wY;
	bool xIsSigned,yIsSigned;
    UserInterface::parseStrictlyPositiveInt(args, "wX", &wX);
    UserInterface::parseStrictlyPositiveInt(args, "wY", &wY);
	UserInterface::parseBoolean(args,"xIsSigned",&xIsSigned);
	UserInterface::parseBoolean(args,"yIsSigned",&yIsSigned);

	return new BaseMultiplierXilinx2xkOp(parentOp,target,xIsSigned,yIsSigned, wX, wY);
}

void BaseMultiplierXilinx2xk::registerFactory()
{
    UserInterface::add("BaseMultiplierXilinx2xk", // name
                        "Implements a 2xY-LUT-Multiplier that can be realized efficiently on some Xilinx-FPGAs",
                       "BasicInteger", // categories
                        "",
                       "wX(int): size of input X;\
                        wY(int): size of input Y;\
						xIsSigned(bool)=0: input X is signed;\
						yIsSigned(bool)=0: input Y is signed;",
                       "",
                       BaseMultiplierXilinx2xk::parseArguments
//                       ,BaseMultiplierXilinx2xk::unitTest
    ) ;
}

void BaseMultiplierXilinx2xkOp::emulate(TestCase* tc)
{
    mpz_class svX = tc->getInputValue("X");
    mpz_class svY = tc->getInputValue("Y");
    mpz_class svR = svX * svY;
    tc->addExpectedOutput("R", svR);
}

TestList BaseMultiplierXilinx2xk::unitTest(int index)
{
    // the static list of mandatory tests
    TestList testStateList;
    vector<pair<string,string>> paramList;

    //test square multiplications:
    for(int w=1; w <= 6; w++)
    {
        paramList.push_back(make_pair("wY", to_string(w)));
        paramList.push_back(make_pair("wX", to_string(2)));
        testStateList.push_back(paramList);
        paramList.clear();
    }
    for(int w=1; w <= 6; w++)
    {
        paramList.push_back(make_pair("wX", to_string(w)));
        paramList.push_back(make_pair("wY", to_string(2)));
        testStateList.push_back(paramList);
        paramList.clear();
    }

    return testStateList;
}

BaseMultiplierXilinx2xkOp::BaseMultiplierXilinx2xkOp(Operator *parentOp, Target* target, bool isSignedX, bool isSignedY, int wX, int wY) : Operator(parentOp,target)
{
    ostringstream name;
    string in1,in2;
    int width;
    bool in2_signed, in1_signed;

    if(wX == 2)
    {
        in1 = "Y";
        in2 = "X";
        name << "BaseMultiplierXilinx2x" << wY;
        width = wY;
        in2_signed = isSignedX;         //X is the short side
        in1_signed = isSignedY;         //Y is the long side
    }
    else
    {
        in1 = "X";
        in2 = "Y";
        name << "BaseMultiplier" << wX << "x2";
        width = wX;
        in2_signed = isSignedY;         //Y is the short side
        in1_signed = isSignedX;         //X is the long side
    }
    setNameWithFreqAndUID(name.str());

    addInput("X", wX, true);
    addInput("Y", wY, true);

    addOutput("R", width+2, 1, true);

    if((wX != 2) && (wY != 2)) throw string("One of the input widths of the BaseMultiplierXilinx2xk has to be 2!");

    int add_lut = (in2_signed || in1_signed) ? 1 : 0;               //The signed cases require a additional LUT
    int needed_luts = width+1;//no. of required LUTs
    int needed_cc = ( (needed_luts+add_lut) / 4 ) + ( (needed_luts+add_lut)  % 4 > 0 ? 1 : 0 ); //no. of required carry chains

    declare( target->logicDelay(5),"cc_s", needed_cc * 4 );                         //TODO Check if delays are actually correct
    declare( target->logicDelay(5),"cc_di", needed_cc * 4 );
    declare( width * target->carryPropagateDelay(), "cc_co", needed_cc * 4 );
    declare( width * target->carryPropagateDelay(),"cc_o", needed_cc * 4 );

    //create the LUTs:
    for(int i=0; i < needed_luts + add_lut; i++)
    {
        //LUT content of the LUTs:
        lut_op lutop_o5, lutop_o6;
        if(in1_signed && !in2_signed) {         //The long side is signed
            if (i == needed_luts) {
                lutop_o6 = lut_in(5);
                lutop_o5 = lut_in(4);
            } else if (i == needed_luts - 1) {
                lutop_o6 = ~(lut_in(0) & lut_in(1)) ^ ~(lut_in(2) & lut_in(1)); //xor of two partial products
                lutop_o5 = ~(lut_in(0) & lut_in(1)); //and of first partial product
            } else {
                lutop_o6 = (lut_in(0) & lut_in(1)) ^ (lut_in(2) & lut_in(3)); //xor of two partial products
                lutop_o5 = lut_in(0) & lut_in(1); //and of first partial product
            }
        } else if(!in1_signed && in2_signed) {       //The short side is signed
            if(i == 0){
                lutop_o6 = (lut_in(3) & ~lut_in(0) & lut_in(2)) | (lut_in(3) & lut_in(0) & lut_in(2));
                lutop_o5 = (~lut_in(3) & lut_in(0) & lut_in(2));
            } else if(i == 1){
                lutop_o6 = (lut_in(3) & ~lut_in(0) & lut_in(2)) | (lut_in(1) & lut_in(0) & ~lut_in(2)) | (~lut_in(3) & lut_in(0) & lut_in(2));
                lutop_o5 = (~lut_in(1) & lut_in(0) & ~lut_in(2));
            } else if(i==needed_luts-1){
                lutop_o6 = (~lut_in(3) & lut_in(0) & ~lut_in(2)) | (lut_in(0) & lut_in(2));
                lutop_o5 = lut_in(4);
            } else if(i==needed_luts){
                lutop_o6 = (lut_in(0) & ~lut_in(2)) | (lut_in(0) & lut_in(2));
                lutop_o5 = lut_in(4);
            } else {
                lutop_o6 = (lut_in(3) & ~lut_in(0) & lut_in(2)) | (~lut_in(1) & lut_in(0) & ~lut_in(2)) | (~lut_in(3) & lut_in(0) & lut_in(2));
                lutop_o5 = lut_in(4);
            }
        } else if(in1_signed && in2_signed) {       //Both sides are signed
            if(i==needed_luts || i==needed_luts-1) {
                lutop_o6 = lut_in(0) & ~lut_in(1) | lut_in(2) & lut_in(1);
                lutop_o5 = lut_in(4);
            } else if(i==1) {
                lutop_o6 = (lut_in(0) & lut_in(1) &  ~(lut_in(2) & lut_in(3))) | ( ~lut_in(0) & lut_in(2) & lut_in(3));
                lutop_o5 = lut_in(0) & ( ~lut_in(1) | (lut_in(2) & lut_in(3)));
            } else if(i==0) {
                lutop_o6 = lut_in(2) & lut_in(3);
                lutop_o5 = lut_in(1) & ~lut_in(3) & lut_in(0) & lut_in(2);
            } else {
                lutop_o6 = lut_in(0) ^ lut_in(0) & lut_in(1) ^ lut_in(2) & lut_in(3);
                lutop_o5 = lut_in(3) & ~lut_in(1) & lut_in(0) & lut_in(2);
            }
        } else {
            lutop_o6 = (lut_in(0) & lut_in(1)) ^ (lut_in(2) & lut_in(3));   //partial product
            lutop_o5 = lut_in(0) & lut_in(1);                                       //carry
        }
        lut_init lutop( lutop_o5, lutop_o6 );

		Xilinx_LUT6_2 *cur_lut = new Xilinx_LUT6_2( this,target );
        cur_lut->setGeneric( "init", lutop.get_hex(), 64 );

        inPortMap("i0",in2 + of(1));
        inPortMap("i2",in2 + of(0));

        if((!in1_signed && !in2_signed) || (in1_signed && !in2_signed)) {
            if(i==0 || i==needed_luts)
                inPortMapCst("i1","'0'"); //connect 0 at LSB position
            else
                inPortMap("i1",in1 + of(i-1));

            if (i == needed_luts - 1 || i == needed_luts)
                inPortMapCst("i3", "'0'"); //connect 0 at MSB position
            else
                inPortMap("i3", in1 + of(i));
        } else {
            if(i==0)
                inPortMap("i1",in1 + of(i+1));
            else if(i==needed_luts)
                inPortMap("i1",in1 + of(i-2));
            else
                inPortMap("i1",in1 + of(i-1));

            if(i==needed_luts)
                inPortMap("i3", in1 + of(i-2));
            else if(i==needed_luts-1)
                inPortMap("i3", in1 + of(i-1));
            else
                inPortMap("i3", in1 + of(i));;
        }

        inPortMapCst("i4","'0'");
        inPortMapCst("i5","'1'");

        outPortMap("o5","cc_di" + of(i));
        outPortMap("o6","cc_s" + of(i));
        vhdl << cur_lut->primitiveInstance( join("lut",i)) << endl;
    }

    //create the carry chain:
    for( int i = 0; i < needed_cc; i++ ) {
		Xilinx_CARRY4 *cur_cc = new Xilinx_CARRY4( this,target );

        inPortMapCst("cyinit", "'0'" );
        if( i == 0 ) {
            inPortMapCst("ci", "'0'" ); //carry-in can not be used as AX input is blocked!!
        } else {
            inPortMap("ci", "cc_co" + of( i * 4 - 1 ) );
        }
        inPortMap("di", "cc_di" + range( i * 4 + 3, i * 4 ) );
        inPortMap("s", "cc_s" + range( i * 4 + 3, i * 4 ) );
        outPortMap("co", "cc_co" + range( i * 4 + 3, i * 4 ));
        outPortMap("o", "cc_o" + range( i * 4 + 3, i * 4 ));

        stringstream cc_name;
        cc_name << "cc_" << i;
        vhdl << cur_cc->primitiveInstance( cc_name.str());
    }
    vhdl << endl;

    declare( width * target->carryPropagateDelay()+target->logicDelay(5) + 9e-10,"result", width+2);
    if(add_lut){
        vhdl << tab << "result <= cc_o(" << width+1 << " downto 0);" << endl;
    } else {
        vhdl << tab << "result <= cc_co(" << width << ") & cc_o(" << width << " downto 0);" << endl;
    }
    vhdl << tab << "R <= result;" << endl;
}

}   //end namespace flopoco

