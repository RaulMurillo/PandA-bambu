#include "BaseMultiplierIrregularLUTXilinx.hpp"
#include "Table.hpp"
#include "IntMultiplier.hpp"

namespace flopoco
{

    Operator* BaseMultiplierIrregularLUTXilinx::generateOperator(
            Operator *parentOp,
            Target* target,
            Parametrization const & parameters) const
    {
        return new BaseMultiplierIrregularLUTXilinxOp(
                parentOp,
                target,
                parameters.isSignedMultX(),
                parameters.isSignedMultY(),
                (TILE_SHAPE)parameters.getShapePara()
        );
    }

    const int BaseMultiplierIrregularLUTXilinx::tile_properties[8][12] = {  {3,3,5,5,5,5,5,5,5,5,1,8},  // A (x,y,r_yuxu,r_yuxs,r_ysxu,r_ysxs,MSB_yuxu,MSB_yuxs,MSB_ysxu,MSB_ysxs,LSB,area)
                                                                            {3,3,3,4,4,4,4,5,5,5,2,5},
                                                                            {3,2,3,4,3,3,3,4,3,3,1,4},
                                                                            {2,3,3,3,4,3,3,3,4,3,1,4},
                                                                            {3,2,4,5,4,4,3,4,3,3,0,5},
                                                                            {3,2,4,4,4,4,4,4,4,4,1,5},
                                                                            {2,3,4,4,5,4,3,3,4,3,0,5},
                                                                            {2,3,4,4,4,4,4,4,4,4,1,5}}; // H

    int BaseMultiplierIrregularLUTXilinx::getRelativeResultLSBWeight(Parametrization const& param) const
    {
        if(!param.getShapePara() || param.getShapePara() > 8)
            throw string("Error in ") + string("srcFileName") + string(": shape unknown");
        return getRelativeResultLSBWeight((TILE_SHAPE)param.getShapePara());
    }

    int BaseMultiplierIrregularLUTXilinx::getRelativeResultMSBWeight(Parametrization const& param) const
    {
        if(!param.getShapePara() || param.getShapePara() > 8)
            throw string("Error in ") + string("srcFileName") + string(": shape unknown");
        return getRelativeResultMSBWeight((TILE_SHAPE)param.getShapePara(), param.isSignedMultX(), param.isSignedMultY());
    }

    double BaseMultiplierIrregularLUTXilinx::getLUTCost(int x_anchor, int y_anchor, int wMultX, int wMultY, bool signedIO){
        bool signedX = signedIO && (wMultX == x_anchor+wX);
        bool signedY = signedIO && (wMultY == y_anchor+wY);
        int word_size = getRelativeResultMSBWeight(this->shape, signedX, signedY) - getRelativeResultLSBWeight(this->shape) + 1;
        int lut_required = (this->wX+this->wY <= 5)?word_size/2+word_size%2:word_size;

        return lut_required + word_size*getBitHeapCompressionCostperBit();

        //TODO Implement position dependent method, although LUT-Multiplier dont seem to be placed so that they are protruding
        /*
        int x_min = ((x_anchor < 0)?0: x_anchor);
        int y_min = ((y_anchor < 0)?0: y_anchor);
        int lsb = x_min + y_min;

        int x_max = ((wMultX < x_anchor + this->wX)?wMultX: x_anchor + this->wX);
        int y_max = ((wMultY < y_anchor + this->wY)?wMultY: y_anchor + this->wY);
        int msb = x_max+y_max;

        int ws = (x_max-x_min==1)?y_max-y_min:((y_max-y_min==1)?x_max-x_min:msb - lsb);
        cout << lut_required + ws*0.65 << " " << shape << " " << ws << " " << lut_required << endl;
        return lut_required + ws*0.65;*/
    }

    bool BaseMultiplierIrregularLUTXilinx::shapeValid(Parametrization const& param, unsigned x, unsigned y) const
    {
        int pattern = BaseMultiplierIrregularLUTXilinxOp::get_pattern((TILE_SHAPE)param.getShapePara());
        if(0 < param.getTileXWordSize() && x < param.getTileXWordSize() && 0 < param.getTileYWordSize() && y < param.getTileYWordSize())
            return (pattern & (1 << (y * param.getTileXWordSize() + x)));
        return false;
    }

    bool BaseMultiplierIrregularLUTXilinx::shapeValid(int x, int y)
    {
        int pattern = BaseMultiplierIrregularLUTXilinxOp::get_pattern(this->shape);
        //cout << "pattern " << pattern << " valid " << (pattern & (1 << (y * wX + x))) << " mask " << ( (1 << (y * wX + x)))  << " wX " << wX   << " wY " << wY << endl;
        if(0 <= x && x < wX && 0 <= y && y < wY)
            return (pattern & (1 << (y * wX + x)));
        return false;
    }

    bool BaseMultiplierIrregularLUTXilinx::shapeValid(TILE_SHAPE shape, int x, int y)
    {
        unsigned short bit_pattern[8] = {0x1fe, 0xf4, 0x1e, 0x1e, 0x1f, 0x3e, 0x1f, 0x3e};
        int wX = get_wX((TILE_SHAPE)shape);
        int wY = get_wY((TILE_SHAPE)shape);
        if(0 <= x && x < wX && 0 <= y && y < wY)
            return (bit_pattern[(shape)-1] & (1 << (y * wX + x)));
        return false;
    }

    BaseMultiplierIrregularLUTXilinxOp::BaseMultiplierIrregularLUTXilinxOp(Operator *parentOp, Target* target, bool isSignedX, bool isSignedY, BaseMultiplierIrregularLUTXilinx::TILE_SHAPE shape) : Operator(parentOp,target), shape(shape), wX(BaseMultiplierIrregularLUTXilinx::get_wX(shape)), wY(BaseMultiplierIrregularLUTXilinx::get_wY(shape)), isSignedX(isSignedX), isSignedY(isSignedY)
    {
        int wR = BaseMultiplierIrregularLUTXilinx::get_wR(shape, isSignedX, isSignedY);

        ostringstream name;
        name << "BaseMultiplierIrregularLUTXilinx_" << wX << (isSignedX==1 ? "_signed" : "") << "x" << wY  << (isSignedY==1 ? "_signed" : "");

        setNameWithFreqAndUID(name.str());

        addInput("X", wX);
        addInput("Y", wY);
        addOutput("R", wR);

        vector<mpz_class> val;
        REPORT(DEBUG, "Filling table for a non-rectangular LUT multiplier of max size " << wX << "x" << wY << " (out put size is " << wR << ")")
        for (int yx=0; yx < 1<<(wX+wY); yx++)
        {
            val.push_back(function(yx));
        }
        Operator *op = new Table(this, target, val, "MultTable", wX+wY, wR);
        op->setShared();
        UserInterface::addToGlobalOpList(op);

        vhdl << declare(0.0,"Xtable",wX+wY) << " <= Y & X;" << endl;

        inPortMap("X", "Xtable");
        outPortMap("Y", "Y1");

        vhdl << tab << "R <= Y1;" << endl;
        vhdl << instance(op, "TableMult");
    }

    const unsigned short BaseMultiplierIrregularLUTXilinxOp::bit_pattern[8] = {0x1fe, 0xf4, 0x1e, 0x1e, 0x1f, 0x3e, 0x1f, 0x3e};

    mpz_class BaseMultiplierIrregularLUTXilinxOp::function(int yx)
    {
        int temp = 0;
        int y = yx>>wX;
        int x = yx -(y<<wX);
        x = (isSignedX && x & (1<<(wX-1)))?x-2*(1<<(wX-1)):x;
        y = (isSignedY && y & (1<<(wY-1)))?y-2*(1<<(wY-1)):y;

        int raw_result = x * y;
        //printf("0x%03x: x=%+03d, y=%+03d, r=%+03d", yx, x, y, raw_result);

        for( int yp = 0; yp < wY; yp++){
            for( int xp = 0; xp < wX; xp++){
                if(!(bit_pattern[(this->shape)-1]&(1<<(yp*wX+xp)))){
                    if((isSignedX && xp==(wX-1)) && (isSignedY && yp==(wY-1)) || !(isSignedX && xp==(wX-1)) && !(isSignedY && yp==(wY-1)) ){
                        raw_result -= (x&(1<<xp))*(y&(1<<yp));
                    } else {
                        raw_result += (x&(1<<xp))*(y&(1<<yp));
                    }
                }
            }
        }
        //printf(", 0x%03x, 0x%03x, msb=%d", raw_result, (1<<wX+wY)-1, BaseMultiplierIrregularLUTXilinx::getRelativeResultMSBWeight(this->shape, this->isSignedX, this->isSignedY));
        raw_result &= (1<<(BaseMultiplierIrregularLUTXilinx::getRelativeResultMSBWeight(this->shape, this->isSignedX, this->isSignedY)+1))-1;
        //printf(", 0x%03x, shift=%01d, 0x%03x\n", raw_result, (unsigned)BaseMultiplierIrregularLUTXilinx::getRelativeResultLSBWeight(this->shape), (unsigned)(raw_result >> BaseMultiplierIrregularLUTXilinx::getRelativeResultLSBWeight(this->shape)));
        mpz_class r  = (raw_result >> BaseMultiplierIrregularLUTXilinx::getRelativeResultLSBWeight(this->shape));

        REPORT(DEBUG, "Value for x=" << x << ", y=" << y << " : " << r.get_str(2))
        return r;

    }

    OperatorPtr BaseMultiplierIrregularLUTXilinx::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
    {
        int wS;
        bool isSignedX, isSignedY;
        UserInterface::parseStrictlyPositiveInt(args, "wS", &wS);
        UserInterface::parseBoolean(args,"xSigned", &isSignedX);
        UserInterface::parseBoolean(args,"ySigned", &isSignedY);

        return new BaseMultiplierIrregularLUTXilinxOp(parentOp,target, isSignedX, isSignedY, (BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)wS);
    }

    void BaseMultiplierIrregularLUTXilinx::registerFactory(){
        UserInterface::add("BaseMultiplierIrregularLUTXilinx", // name
                           "Implements a non rectangular LUT multiplier from a set that yields a relatively high efficiency compared to recangular LUT multipliers \n"
            "  _ _     _        _ _       _     _ _ _    _ _      _ _     _\n"
            " |_|_|_  |_|_     |_|_|_    |_|_  |_|_|_|  |_|_|_   |_|_|   |_|_\n"
            " |_|_|_| |_|_|_     |_|_|   |_|_|   |_|_|  |_|_|_|  |_|_|   |_|_|\n"
            " |_|_|_|   |_|_|              |_|                     |_|   |_|_|\n"
            " shape 1  shape 2 shape 3 shape 4 shape 5  shape 6 shape 7  shape 8\n"

                           ,
                           "BasicInteger", // categories
                           "",
                           "wS(int): shape ID;\
                            xSigned(bool)=false: input X can be signed or unsigned;\
                            ySigned(bool)=false: input Y can be signed or unsigned;",
                           "",
                           BaseMultiplierIrregularLUTXilinx::parseArguments,
                           BaseMultiplierIrregularLUTXilinxOp::unitTest
        ) ;

    }

    int BaseMultiplierIrregularLUTXilinx::ownLUTCost(int x_anchor, int y_anchor, int wMultX, int wMultY, bool signedIO) {
        bool signedX = signedIO && (wMultX == x_anchor+wX);
        bool signedY = signedIO && (wMultY == y_anchor+wY);
        int word_size = getRelativeResultMSBWeight(this->shape, signedX, signedY) - getRelativeResultLSBWeight(this->shape) + 1;
        int lut_required = (this->wX+this->wY <= 5)?word_size/2+word_size%2:word_size;

        return lut_required;
    }

    void BaseMultiplierIrregularLUTXilinxOp::emulate(TestCase* tc)
    {
        mpz_class svX = tc->getInputValue("X");
        mpz_class svY = tc->getInputValue("Y");
        mpz_class svR = 0;

        svR = svX * svY;

        for( int yp = 0; yp < wY; yp++){
            for( int xp = 0; xp < wX; xp++){
                if(!(bit_pattern[(this->shape)-1]&(1<<(yp*wX+xp)))){
                    if((isSignedX && xp==(wX-1)) && (isSignedY && yp==(wY-1)) || !(isSignedX && xp==(wX-1)) && !(isSignedY && yp==(wY-1)) ){
                        svR -= (svX&(1<<xp))*(svY&(1<<yp));
                    } else {
                        svR += (svX&(1<<xp))*(svY&(1<<yp));
                    }
                }
            }
        }

        svR &= (1<<(BaseMultiplierIrregularLUTXilinx::getRelativeResultMSBWeight(this->shape, this->isSignedX, this->isSignedY)+1))-1;
        svR >>= BaseMultiplierIrregularLUTXilinx::getRelativeResultLSBWeight(this->shape);
        tc->addExpectedOutput("R", svR);
    }

    TestList BaseMultiplierIrregularLUTXilinxOp::unitTest(int index)
    {
        // the static list of mandatory tests
        TestList testStateList;
        vector<pair<string,string>> paramList;

        //test square multiplications:
        for(int w=1; w <= 6; w++)
        {
            paramList.push_back(make_pair("wS", to_string(w)));
            testStateList.push_back(paramList);
            paramList.clear();
        }

        return testStateList;
    }

    /*    int BaseMultiplierIrregularLUTXilinx::getRelativeResultMSBWeight(TILE_SHAPE shape, bool isSignedX, bool isSignedY){
        unsigned short bit_pattern[8] = {0x1fe, 0xf4, 0x1e, 0x1e, 0x1f, 0x3e, 0x1f, 0x3e};
        int wX = get_wX((TILE_SHAPE)shape);
        int wY = get_wY((TILE_SHAPE)shape);
        int min=INT16_MAX, max=INT16_MIN, min_bits=INT16_MAX, n_max_bits=0, p_max_bits=0, max_bits=0;
        for (int yx=0; yx < 1<<(wX+wY); yx++)
        {
            int y = yx>>wX;
            int x = yx -(y<<wX);
            x = (isSignedX && x & (1<<(wX-1)))?x-2*(1<<(wX-1)):x;
            y = (isSignedY && y & (1<<(wY-1)))?y-2*(1<<(wY-1)):y;

            int result, raw_result = x * y;
            //printf("0x%03x: x=%+03d, y=%+03d, r=%+03d", yx, x, y, raw_result);
            for( int yp = 0; yp < wY; yp++){
                for( int xp = 0; xp < wX; xp++){
                    if(!(bit_pattern[(shape)-1]&(1<<(yp*wX+xp)))){
                        if((isSignedX && xp==(wX-1)) && (isSignedY && yp==(wY-1)) || !(isSignedX && xp==(wX-1)) && !(isSignedY && yp==(wY-1)) ){
                            raw_result -= (x&(1<<xp))*(y&(1<<yp));
                        } else {
                            raw_result += (x&(1<<xp))*(y&(1<<yp));
                        }
                    }
                }
                //raw_result -= (!(bit_pattern[(this->shape)-1]&(1<<(yp*wX+xp))))?(x&(1<<xp))*(y&(1<<yp)):0;
            }
            //printf(", 0x%03x, 0x%03x, msb=%d", raw_result, (1<<wX+wY)-1, getRelativeResultMSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape, isSignedX || isSignedY));
            //raw_result &= (1<<(getRelativeResultMSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape, isSignedX || isSignedY)+1))-1;
            //printf(", 0x%03x, shift=%01d, 0x%03x\n", raw_result, (unsigned)getRelativeResultLSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape), (unsigned)(raw_result >> getRelativeResultLSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape)));
            //result  = (raw_result >> getRelativeResultLSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape));
            result = raw_result;
            min = (result < min)?result:min;
            max = (max < result)?result:max;
            for(int i = 0; i < 9; i++){
                n_max_bits = (result < 0)?((n_max_bits<=i && ((result & (1<<i+1)) != 0) && ((result & (1<<(i))) == 0))?i+1:n_max_bits):n_max_bits;
                p_max_bits = (result < 0)?p_max_bits:(p_max_bits<=i && ((result & (1<<i)) != 0))?i:p_max_bits;
                min_bits = (i < min_bits && (result & (1<<i)))?i:min_bits;
                //msbp = (result < 0)?((max_bits<=i && ((result & (1<<i+1)) != 0) && ((result & (1<<(i))) == 0))?0:msbp):(max_bits<=i && ((result & (1<<i)) != 0))?1:msbp;
            }

        }
        max_bits = (p_max_bits < n_max_bits)?n_max_bits:(n_max_bits == p_max_bits)?p_max_bits+1:(min<0)?p_max_bits+1:p_max_bits;
        //max_bits = (msbp && min < 0)?max_bits+1:max_bits; //If the MSB was detected a a positive value one additional bit is required in signed case
        //printf("shape=%d, signedX=%d, singedY=%d, r_min=%+03d, r_max=%+03d, min_bits=%d, max_bits=%d, bits_r=%d\n", shape, isSignedX, isSignedY, min, max, min_bits, max_bits, max_bits-min_bits+1);
        return max_bits;
    }

    void BaseMultiplierIrregularLUTXilinx::draw_tile(TILE_SHAPE shape, bool isSignedX, bool isSignedY, stringstream *stream_handle ){
        stringstream outline;
        ofstream result_file;
        if(stream_handle == nullptr){
            result_file.open("irregular_tiles.tex");
            outline <<
            tab << "\\documentclass[margin=0pt]{standalone}" << endl <<
            tab << "\\usepackage{graphicx}" << endl <<
            tab << "\\usepackage{xcolor,tikz,tikzscale}" << endl <<
            tab << "\\begin{document}" << endl <<
            tab << "\\begin{tikzpicture}[yscale=-1,xscale=-1]" << endl;
        }

        unsigned short bit_pattern[8] = {0x1fe, 0xf4, 0x1e, 0x1e, 0x1f, 0x3e, 0x1f, 0x3e};
        int wX = get_wX((TILE_SHAPE)shape);
        int wY = get_wY((TILE_SHAPE)shape);

        outline << tab <<"\\foreach \\x in {0,...," << wX << "}{" << endl <<
        tab << tab << "\\draw[very thin, dotted, gray] (0+\\x, 0) -- (0+\\x, 0+" << wY << ");" << endl <<
        tab << tab <<"}" << endl <<
        tab << tab <<"\\foreach \\y in {0,...," << wY << "}{" << endl <<
        tab << tab <<"\\draw[very thin, dotted, gray] (0+0, 0+\\y) -- (0+" << wX << ", 0+\\y);" << endl <<
        tab << tab <<"}" << endl;

        int sx = 0, sy = 0, cx = 0, cy = 0, lx = 0, ly = 0, xdir = 1, ydir=-1;
        while((bit_pattern[(shape)-1]&(1<<(sy*wX+sx))) == 0)
            sy++;
        while((bit_pattern[(shape)-1]&(1<<((1+sy)*wX+sx))) != 0)
            sy++;
        cy = sy;
        cout << "x=" << cx << " y=" << cy << " xdir=" << xdir << " ydir=" << ydir << endl;
        outline << tab << tab << "\\draw[black, thick, fill=lightgray, opacity=0.3] (" << cx << "," << cy << ")--";
        do{
            lx = cx; ly = cy;
            if       (!shapeValid(shape, cx,cy+ydir) && !shapeValid(shape, cx+xdir,cy) && !shapeValid(shape, cx+xdir,cy+ydir)){
                xdir *= (-1); ydir *= (-1);
                outline << "(" << (cx+((-1==xdir)?1:0)) << "," << (cy+((-1==ydir)?1:0)) << ")--";
            } else if(!shapeValid(shape, cx,cy+ydir) && !shapeValid(shape, cx+xdir,cy) &&  shapeValid(shape, cx+xdir,cy+ydir)){
                cx += xdir; cy += ydir;
            } else if(!shapeValid(shape, cx,cy+ydir) &&  shapeValid(shape, cx+xdir,cy) && !shapeValid(shape, cx+xdir,cy+ydir)){
                cx += xdir;
            } else if(!shapeValid(shape, cx,cy+ydir) &&  shapeValid(shape, cx+xdir,cy) &&  shapeValid(shape, cx+xdir,cy+ydir)){
                cx += xdir;
            } else if( shapeValid(shape, cx,cy+ydir) && !shapeValid(shape, cx+xdir,cy) && !shapeValid(shape, cx+xdir,cy+ydir)){
                cy += ydir;
            } else if( shapeValid(shape, cx,cy+ydir) && !shapeValid(shape, cx+xdir,cy) &&  shapeValid(shape, cx+xdir,cy+ydir)){
                cy += ydir;
            } else if( shapeValid(shape, cx,cy+ydir) &&  shapeValid(shape, cx+xdir,cy) && !shapeValid(shape, cx+xdir,cy+ydir)){
                cy += ydir;
            } else if( shapeValid(shape, cx,cy+ydir) &&  shapeValid(shape, cx+xdir,cy) &&  shapeValid(shape, cx+xdir,cy+ydir)){
                cy += ydir;
            }
            cout << "x=" << cx << " y=" << cy << " xdir=" << xdir << " ydir=" << ydir << endl;
            outline << "(" << (cx+((-1==xdir)?1:0)) << "," << (cy+((1==ydir)?1:0)) << ")--";

        } while(sx != cx || sy != cy);
        xdir *= (-1); ydir *= (-1);
        outline << "(" << (cx+((-1==xdir)?1:0)) << "," << (cy+((-1==ydir)?1:0)) << ")--";
        outline << "cycle;\n";

        for( int yp = 0; yp < wY; yp++){
            for( int xp = 0; xp < wX; xp++){
                if(bit_pattern[(shape)-1]&(1<<(yp*wX+xp))){
                    int xval = (isSignedX && xp==(wX-1))?-(1<<xp):(1<<xp);
                    int yval = (isSignedY && yp==(wY-1))?-(1<<yp):(1<<yp);
                    outline << tab << tab << "\\node[] at ("<< xp+0.85 << "," << yp+0.15 <<") {\\tiny"<< xval << "};\n";
                    outline << tab << tab << "\\node[] at ("<< xp+0.15 << "," << yp+0.85 <<") {\\tiny"<< yval << "};\n";
                    outline << tab << tab << "\\node[] at ("<< xp+0.5 << "," << yp+0.5 <<") {"<< xval*yval << "};\n";
                }
            }
        }

        //outline << tab << tab << "\\draw (-0.25," << wY/2 << "+0.6)  node {\\footnotesize $\\downarrow$};\n";
        outline << tab << tab << "\\draw (-0.25," << wY/2.0 << ")  node[rotate=90] {\\tiny $\\leftarrow y=\\left\\{" << ((isSignedY)?-(1<<(wY-1)):0) << ".." << ((isSignedY)?(1<<(wY-1))-1:(1<<wY)-1)  << "\\right\\}$};\n";
        outline << tab << tab << "\\draw (" << wX/2.0 << ",-0.2)  node {\\footnotesize $\\leftarrow$\\tiny $x=\\left\\{" << ((isSignedX)?-(1<<(wX-1)):0) << ".." << ((isSignedX)?(1<<(wX-1))-1:(1<<wX)-1)  << "\\right\\}$};\n";

        if(stream_handle == nullptr){
            outline << tab << "\\end{tikzpicture}" << endl <<
            "\\end{document}" << endl;

            result_file << outline.str();
            result_file.close();
        } else {
            *stream_handle  << outline.str();
        }
    }

    void BaseMultiplierIrregularLUTXilinx::calc_tile_properties(TILE_SHAPE shape, bool isSignedX, bool isSignedY, int &min, int &max, int &min_bits, int &max_bits, int &bits){
        unsigned short bit_pattern[8] = {0x1fe, 0xf4, 0x1e, 0x1e, 0x1f, 0x3e, 0x1f, 0x3e};
        int wX = BaseMultiplierIrregularLUTXilinx::get_wX(shape);
        int wY = BaseMultiplierIrregularLUTXilinx::get_wY(shape);
        int n_max_bits=0, p_max_bits=0; min=INT16_MAX; max=INT16_MIN; min_bits=INT16_MAX; max_bits=0;
        for (int yx=0; yx < 1<<(wX+wY); yx++)
        {
            int y = yx>>wX;
            int x = yx -(y<<wX);
            x = (isSignedX && x & (1<<(wX-1)))?x-2*(1<<(wX-1)):x;
            y = (isSignedY && y & (1<<(wY-1)))?y-2*(1<<(wY-1)):y;

            int result, raw_result = x * y;
            //printf("0x%03x: x=%+03d, y=%+03d, r=%+03d", yx, x, y, raw_result);
            for( int yp = 0; yp < wY; yp++){
                for( int xp = 0; xp < wX; xp++){
                    if(!(bit_pattern[(shape)-1]&(1<<(yp*wX+xp)))){
                        if((isSignedX && xp==(wX-1)) && (isSignedY && yp==(wY-1)) || !(isSignedX && xp==(wX-1)) && !(isSignedY && yp==(wY-1)) ){
                            raw_result -= (x&(1<<xp))*(y&(1<<yp));
                        } else {
                            raw_result += (x&(1<<xp))*(y&(1<<yp));
                        }
                    }
                }
                //raw_result -= (!(bit_pattern[(this->shape)-1]&(1<<(yp*wX+xp))))?(x&(1<<xp))*(y&(1<<yp)):0;
            }
            //printf(", 0x%03x, 0x%03x, msb=%d", raw_result, (1<<wX+wY)-1, getRelativeResultMSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape, isSignedX || isSignedY));
            //raw_result &= (1<<(getRelativeResultMSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape, isSignedX || isSignedY)+1))-1;
            //printf(", 0x%03x, shift=%01d, 0x%03x\n", raw_result, (unsigned)getRelativeResultLSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape), (unsigned)(raw_result >> getRelativeResultLSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape)));
            //result  = (raw_result >> getRelativeResultLSBWeight((BaseMultiplierIrregularLUTXilinx::TILE_SHAPE)shape));
            result = raw_result;
            min = (result < min)?result:min;
            max = (max < result)?result:max;
            for(int i = 0; i < 9; i++){
                n_max_bits = (result < 0)?((n_max_bits<=i && ((result & (1<<i+1)) != 0) && ((result & (1<<(i))) == 0))?i+1:n_max_bits):n_max_bits;
                p_max_bits = (result < 0)?p_max_bits:(p_max_bits<=i && ((result & (1<<i)) != 0))?i:p_max_bits;
                min_bits = (i < min_bits && (result & (1<<i)))?i:min_bits;
                //msbp = (result < 0)?((max_bits<=i && ((result & (1<<i+1)) != 0) && ((result & (1<<(i))) == 0))?0:msbp):(max_bits<=i && ((result & (1<<i)) != 0))?1:msbp;
            }

        }
        max_bits = (p_max_bits < n_max_bits)?n_max_bits:(n_max_bits == p_max_bits)?p_max_bits+1:(min<0)?p_max_bits+1:p_max_bits;
        //max_bits = (msbp && min < 0)?max_bits+1:max_bits; //If the MSB was detected a a positive value one additional bit is required in signed case
        //printf("shape=%d, signedX=%d, singedY=%d, r_min=%+03d, r_max=%+03d, min_bits=%d, max_bits=%d, bits_r=%d\n", shape, isSignedX, isSignedY, min, max, min_bits, max_bits, max_bits-min_bits+1);
        bits = max_bits-min_bits+1;
    }

    void BaseMultiplierIrregularLUTXilinx::draw_property_sheet(void){
        stringstream outline;
        ofstream result_file;
        result_file.open("irregular_tile_properties.tex");
        outline <<
        tab << "\\documentclass[margin=0pt]{standalone}" << endl <<
        tab << "\\usepackage{graphicx}" << endl <<
        tab << "\\usepackage{xcolor,tikz,tikzscale}" << endl <<
        tab << "\\begin{document}" << endl <<
        tab << "\\begin{tikzpicture}[yscale=-1,xscale=-1]" << endl;

        for(int shape = 1; shape <=8; shape++){
            for(int isSignedY = 0; isSignedY<2; isSignedY++){
                for(int isSignedX = 0; isSignedX<2; isSignedX++){
                    outline << tab << "\\begin{scope}[xshift=" << -4*(2*isSignedY+1*isSignedX) << "cm, yshift=" << 5*shape << "cm]"  << endl;
                    int min, max, min_bits, max_bits, bits;
                    calc_tile_properties((TILE_SHAPE)shape, isSignedX, isSignedY, min, max, min_bits, max_bits, bits);
                    draw_tile((TILE_SHAPE)shape ,isSignedX, isSignedY, &outline);
                    outline << tab << tab << "\\node[] at ("<< 1.5 << "," << 3.4 <<") {$r=\\left\\{" << min << ".." << max << "\\right\\}$};\n";
                    outline << tab << tab << "\\node[] at ("<< 1.5 << "," << 3.8 <<") {MSB=" << max_bits << ", LSB=" << min_bits << "};\n";
                    outline << tab << tab << "\\node[] at ("<< 1.5 << "," << 4.2 <<") {$w_{Out}=" << bits << "$ (LUT)};\n";
                    outline << tab << "\\end{scope}"  << endl;
                }
            }
        }

        outline << tab <<"\\foreach \\g in {-1,...," << 3 << "}{" << endl <<
        tab << tab << "\\draw[black] (-0.5-4*\\g, +3.5) -- (-0.5-4*\\g, -0.5+" << 9*5 << ");" << endl <<
        tab << tab <<"}" << endl <<
        tab << tab <<"\\foreach \\s in {1,...," << 9 << "}{" << endl <<
        tab << tab <<"\\draw[black] (4.5+0, -0.5+5*\\s) -- (-0.5+" << -4*3 << ", -0.5+5*\\s);" << endl <<
        tab << tab <<"}" << endl <<
        tab << tab <<"\\foreach \\s in {1,...," << 8 << "}{" << endl <<
        tab << tab << "\\node[] at ("<< 4 << ", 5*\\s+2.0) {\\Huge $\\s$};" << endl <<
        tab << tab <<"}" << endl <<
        tab << tab <<"\\foreach \\ys in {0,...," << 1 << "}{" << endl <<
        tab << tab <<"\\foreach \\xs in {0,...," << 1 << "}{" << endl <<
        tab << tab << "\\node[align=center] at (1.5-8*\\ys-4*\\xs , 4) {isSignedX=\\xs \\\\ isSignedY=\\ys};" << endl <<
        tab << tab <<"}}" << endl;

        outline << tab << "\\end{tikzpicture}" << endl <<
        "\\end{document}" << endl;

        result_file << outline.str();
        result_file.close();
    }*/


}
