//
// Created by jk on 9/21/20.
//

/*
 *  Description:
 *  IIR filter with multiplierless option (graph in rpag format mandatory for multiplierless)
 *  If no coeffa is given the filter turns into a FIR filter (feedback path is omitted)
 *  IIR and FIR variant use faithful rounding
 *  Caution: shifta/b denotes the number of (right)shitfs, not scaling. Always give positive number
 *
 *
 *  Optimizations:
 *  Word sizes for adder and register have conservative computation. Using wcpg can optimze the wordsizes for each
 *  adder. Conservative computation does not apply for feedbackAdd0, this one is good
 *
 */

#include "FixIIRShiftAdd.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>

#if HAVE_WCPG
extern "C"
{
#include "wcpg.h"
}
#endif

// this sign-macro is introduced because the initial computation ignored the sign bit resp. only works with positive values
// so e.g 576 is okay with 10 bits, but since it is converted into SIGNED integer, we need one bit more
// Later implementation may not use singed values, so it can be easily set to zero using this macro
#define TEMP_SIGN_FOR_RESIZE 1

using namespace std;

// test with
// flopoco FixIIRShiftAdd msbIn=16 lsbIn=-8 lsbOut=-8 shifta=3 shiftb=0 coeffb="576:512:128:7:8" coeffa="2:4:1:7" guardBits=2 verbose=2 testbench n=0
// make -j4 && ./flopoco FixIIRShiftAdd msbIn=4 lsbIn=-4 lsbOut=-4 shifta=3 shiftb=2 coeffb="7:15:7:7:15" coeffa="2:4:1:7" guardBits=0 verbose=0 grapha="{{'R',[1],1,[1],0},{\'A\',[7],1,[1],0,3,[-1],0,0},{'O',[1],1,[1],1,0},{'O',[2],1,[1],1,1},{'O',[4],1,[1],1,2},{'O',[7],1,[7],1,0}}" graphb="{{\'A\',[7],1,[1],0,3,[-1],0,0},{\'A\',[15],1,[1],0,4,[-1],0,0},{'O',[7],1,[7],1,0},{'O',[15],1,[15],1,0}}" method="multiplierless" testbench n=0
// is coeffa 1/2:1/4 bc shifta = 3 (scalea = 1/(2^3) = 0.125)

namespace flopoco {
    FixIIRShiftAdd::FixIIRShiftAdd(OperatorPtr parentOp, Target* target, int lsbIn_, int lsbOut_, int msbIn_, int guardBits_, vector<string> coeffb_, vector<string> coeffa_, int64_t shifta_, int64_t shiftb_, string method_, string grapha_, string graphb_, double H_, double Heps_) :
            Operator(parentOp, target), lsbIn(lsbIn_), lsbOut(lsbOut_), msbIn(msbIn_), guardBits(guardBits_), coeffb(coeffb_), coeffa(coeffa_), shifta(shifta_), shiftb(shiftb_), method(method_), grapha(grapha_), graphb(graphb_), H(H_), Heps(Heps_) {
        srcFileName = "FixIIRShiftAdd";
        setCopyrightString("Jonas Kuehle (2020)");

        useNumericStd();
        useNumericStd_Signed();
        useNumericStd_Unsigned();

        ostringstream name;
        name << "FixIIRShiftAdd";
        setNameWithFreqAndUID(name.str());

        m = coeffa.size();
        n = coeffb.size();

        coeffa_d = (double *) malloc(m * sizeof(double));
        coeffa_si = (int64_t *) malloc(m * sizeof(int64_t));
        coeffb_d = (double *) malloc(n * sizeof(double));
        coeffb_si = (int64_t *) malloc(n * sizeof(int64_t));
        coeffa_mp_f = (mpfr_t *) malloc(m * sizeof(mpfr_t));
        coeffb_mp_f = (mpfr_t *) malloc(n * sizeof(mpfr_t));
        coeffa_mp_f_scaled = (mpfr_t *) malloc(m * sizeof(mpfr_t));
        coeffb_mp_f_scaled = (mpfr_t *) malloc(n * sizeof(mpfr_t));
        registerForwardTB_mpz = (mpz_t*)malloc(sizeof(mpz_t) * n);
        registerFeedbackTB_mpz = (mpz_t*)malloc(sizeof(mpz_t) * (m+1));

        for (uint32_t i = 0; i < n; i++) {
            // parse the coeffs from the string, with Sollya parsing
            sollya_obj_t node;

            node = sollya_lib_parse_string(coeffb[i].c_str());
            if (node == 0)                    // If conversion did not succeed (i.e. parse error)
            THROWERROR("Unable to parse string " << coeffb[i] << " as a numeric constant");

            mpfr_init2(coeffb_mp_f[i], 10000);
            mpfr_init2(coeffb_mp_f_scaled[i], 10000);
            sollya_lib_get_constant(coeffb_mp_f[i], node);
            //mpz_t test;
            //sollya_lib_get_constant_as_mpz(test, node);
            coeffb_d[i] = mpfr_get_d(coeffb_mp_f[i], GMP_RNDN);
            coeffb_si[i] = mpfr_get_si(coeffb_mp_f[i], GMP_RNDN);

            REPORT(DETAILED, "b[" << i << "]=" << setprecision(15) << scientific << coeffb_d[i]);
        }

        for (uint32_t i = 0; i < m; i++) {
            // parse the coeffs from the string, with Sollya parsing
            sollya_obj_t node;

            node = sollya_lib_parse_string(coeffa[i].c_str());
            if (node == 0)                    // If conversion did not succeed (i.e. parse error)
            THROWERROR("Unable to parse string " << coeffb[i] << " as a numeric constant");

            mpfr_init2(coeffa_mp_f[i], 10000);
            mpfr_init2(coeffa_mp_f_scaled[i], 10000);
            sollya_lib_get_constant(coeffa_mp_f[i], node);
            coeffa_d[i] = mpfr_get_d(coeffa_mp_f[i], GMP_RNDN);
            coeffa_si[i] = mpfr_get_si(coeffa_mp_f[i], GMP_RNDN);
            REPORT(DETAILED, "a[" << i << "]=" << scientific << setprecision(15) << coeffa_d[i]);
        }

        // check if FIR or IIR, if coeffa is -1 or shifta is -1 its detected as FIR
        if ((m == 1 && coeffa_si[0] == -1) || shifta == -1)
            isFIR = true;

        wcpgBitGain = 0;

        if(!isFIR) // wcpg and guard bits only necessary for IIR
        {
#if HAVE_WCPG
            coeffa_d_wcpg = (double *) malloc(m * sizeof(double));
            coeffb_d_wcpg = (double *) malloc(n * sizeof(double));
            if(H == 0) // H not set by argument
            {
                REPORT(INFO, "Computing worst-case peak gain");
                for (uint32_t i = 0; i < m; i++) {
                    coeffa_d_wcpg[i] = coeffa_d[i] * (1 / (pow(2, shifta)));
                    REPORT(DETAILED,
                           "coeffa_d_wcpg: " << coeffa_d_wcpg[i] << ", coeffa: " << coeffa_d[i] << ", shifta: "
                                             << shifta)
                }
                for (uint32_t i = 0; i < n; i++)
                    coeffb_d_wcpg[i] = coeffb_d[i] * (1 / (pow(2, shiftb)));

                if (!WCPG_tf(&H, coeffb_d_wcpg, coeffa_d_wcpg, n, m, (int) 0)) THROWERROR("Could not compute H");
                REPORT(INFO, "Computed filter worst-case peak gain: H=" << H)
            }
            else // H set by argument -> overwrite
            {
                REPORT(LIST, "H=" << H);
            }

            wcpgBitGain = intlog2(H);
            REPORT(INFO, "wcpg bit gain: " << wcpgBitGain);

            // ################# COMPUTE GUARD BITS AND HEPS #########################################
            if(guardBits < 0)
            {
                if (Heps == 0) // if guard bits not set by argument then overwrite
                {
                    REPORT(DETAILED, "computing guard bits");
                    double one_d[1] = {1.0};
                    if (!WCPG_tf(&Heps, one_d, coeffa_d_wcpg, 1, m, (int) 0)) THROWERROR("Could not compute Heps");

                    REPORT(LIST, "Computed error amplification worst-case peak gain: Heps=" << Heps);
                }
                else
                {
                    REPORT(LIST, "Heps=" << Heps);
                }
                guardBits = intlog2(Heps) + 1; // +1 for last bit accuracy
            }

    #else
                THROWERROR("WCPG was not found (see cmake output), cannot compute worst-case peak gain H. Compile FloPoCo with WCPG");
    #endif

            REPORT(INFO, "No of guard bits=" << guardBits);
        }
        wIn = msbIn - lsbIn + 1;
        msbInt = ceil(msbIn + wcpgBitGain);
        msbOutIIR = msbInt;
        lsbInt = lsbOut - guardBits;
        //int wInt = wOut + guardBits; // not needed yet
        REPORT(INFO, "msbIn: \t\t\t\t" << msbIn << ", lsbIn: \t\t\t\t " << lsbIn)
        REPORT(INFO, "msbInt: \t\t\t\t" << msbInt << ", lsbInt: \t\t\t\t" << lsbInt)

        // ################################### WORD SIZES FORWARD PATH ############################################

        msbRegisterForwardOut = (int *) malloc(n * sizeof(int));
        msbAdditionsForwardOut = (int *) malloc(n * sizeof(int));
        msbMultiplicationsForwardOut = (int *) malloc(n * sizeof(int));
        lsbForwardPath = lsbIn; // lsb is constant in forward path until scaling

        /* for n coefficients there are n multiplication, but n-1 registers and additions
         * There is a special case for sign computation for plain
         * If input it -128 and coeff is -2, there resulting number will be 256 which need two more bits, but log2(abs(-2)) = 1
         * This applies to all negative powers of two. It affects word sizes of multiplications, registers and additions
         * This cannot appear in multiplierless since only positive coeffs are used
         */
        // multiplications
        for(uint32_t i = 0; i < n; i++)
        {
            if(coeffb_si[i] == 0)
                msbMultiplicationsForwardOut[i] = msbIn;
            else
                if(method == "plain" && isPowerOfTwo(coeffb_si[i], true))
                    msbMultiplicationsForwardOut[i] = msbIn + ceil(log2(abs(coeffb_si[i])+1));
                else
                    msbMultiplicationsForwardOut[i] = msbIn + ceil(log2(abs(coeffb_si[i])));
            REPORT(DEBUG, "msbMultiplicationsForwardOut["<<i<<"]: " << msbMultiplicationsForwardOut[i])
        }

        // registers (only word size is needed for vhdl code generation, but it can be derived from msb and lsb)
        for(int i = n-2; i >= 0; i--)
        {
            int64_t bitGainForRegister = 0;
            for (uint32_t j = i + 1; j < n; j++)
            {
                if(method == "plain" && isPowerOfTwo(coeffb_si[j], true)) // see comment multiplications forward path
                    bitGainForRegister += abs(coeffb_si[j])+1;
                else
                    bitGainForRegister += abs(coeffb_si[j]);
            }
            if(bitGainForRegister == 0) // avoid log2(0)
                msbRegisterForwardOut[i] = msbIn;
            else
                msbRegisterForwardOut[i] = msbIn + ceil(log2((abs(bitGainForRegister))));
            REPORT(DEBUG, "msbRegisterForwardOut["<<i<<"]: " << msbRegisterForwardOut[i])
        }

        // additions, last addition (addN-1) is special case since it has two inputs
        // assert necessary for implementation reasons
        // It would be possible to copy msb from preceding register
        if (n < 2)
            THROWERROR("Number of coefficients in forward path (coeffb) must be larger or equal two");
        if(abs(coeffb_si[n - 1]) + abs(coeffb_si[n - 2]) == 0)
            msbAdditionsForwardOut[n - 1] = msbIn + 0;
        else
            msbAdditionsForwardOut[n - 1] = msbIn + ceil(log2((abs(coeffb_si[n - 1]) + abs(coeffb_si[n - 2]))));
        REPORT(DEBUG, "msbAdditionsForwardOut["<<n - 1<<"]: " << msbAdditionsForwardOut[n-1])
        bitGainForAdditons = 0;
        for (int i = n - 2; i >= 0; i--) // for additions except addN-1
        {
            bitGainForAdditons = 0;
            for (uint32_t j = i; j < n; j++)
                if(method == "plain" && isPowerOfTwo(coeffb_si[j], true)) // see comment multiplications forward path
                    bitGainForAdditons += abs(coeffb_si[j])+1;
                else
                    bitGainForAdditons += abs(coeffb_si[j]);
            REPORT(INFO, "bitGainForAdditons: " << ceil(log2((abs(bitGainForAdditons)))));
            if(bitGainForAdditons == 0)
                msbAdditionsForwardOut[i] = msbIn;
            else
                msbAdditionsForwardOut[i] = msbIn + ceil(log2((abs(bitGainForAdditons))));
            REPORT(DEBUG, "msbAdditionsForwardOut["<<i<<"]: " << msbAdditionsForwardOut[i])
        }

        // scaleb
        msbScaleBOut = msbAdditionsForwardOut[0] - abs(shiftb);
        lsbScaleBOut = lsbForwardPath - abs(shiftb);
        REPORT(INFO, "msbScaleBOut: \t\t\t" << msbScaleBOut << ", lsbScaleBOut: \t\t\t" << lsbScaleBOut)

        // msbOut is computed different for FIR and IIR since the output is at different position
        msbOutFIR = msbScaleBOut;
        if(!isFIR)
            wOut = msbOutIIR - lsbOut + 1;
        else
            wOut = msbOutFIR - lsbOut + 1;
        cout << "isFir: " << isFIR << endl;

        // ################################### WORD SIZES BACKPATH ###############################################
        if(!isFIR) {
            // check if replaceble by wcpg
            int64_t wordsizeGainFeedback = 0; // growth of wordsize in feedback path. Used multiple times for wird size / msb / lsb calculation
            for (uint32_t i = 0; i < m; i++)
                wordsizeGainFeedback += (abs(coeffa_si[i]));
            wordsizeGainFeedback = ceil(log2(wordsizeGainFeedback));
            REPORT(INFO, "wordsizeGainFeedback: " << wordsizeGainFeedback);

            msbRegisterFeedbackOut = (int *) malloc(m * sizeof(int));
            lsbRegisterFeedbackOut = (int *) malloc(m * sizeof(int));
            msbMultiplicationsFeedbackOut = (int *) malloc(m * sizeof(int));
            lsbMultiplicationsFeedbackOut = (int *) malloc(m * sizeof(int));
            msbAdditionsFeedbackOut = (int *) malloc(m * sizeof(int));
            lsbAdditionsFeedbackOut = (int *) malloc(m * sizeof(int));

            // Multiplications
            for (uint32_t i = 0; i < m; i++) {
                if (coeffa_si[i] == 0) // avoid log(0)
                    msbMultiplicationsFeedbackOut[i] =
                            msbInt - shifta +
                            1; // use minus shifta since scaling is after truncation and must be kept in mind
                else if (method == "plain" &&
                         isPowerOfTwo(coeffa_si[i], false)) // see comment multiplications forward path
                    msbMultiplicationsFeedbackOut[i] = msbInt - shifta + ceil(log2(abs(coeffa_si[i]) + 1));
                else
                    msbMultiplicationsFeedbackOut[i] = msbInt - shifta + ceil(log2(abs(coeffa_si[i])));
                lsbMultiplicationsFeedbackOut[i] = lsbInt - shifta;
                REPORT(INFO, "msbMultiplicationsFeedbackOut[" << i << "]: \t" << msbMultiplicationsFeedbackOut[i]
                                                              << ", lsbMultiplicationsFeedbackOut[" << i << "] \t"
                                                              << lsbMultiplicationsFeedbackOut[i])
            }

            // registers
            // values of registers are computed first -> values for additions can be derived
            for (int i = m - 1; i >= 0; i--) {
                int64_t tmp = 0;
                for (uint32_t j = i; j < m; j++) {
                    if (method == "plain" &&
                        isPowerOfTwo(coeffa_si[j], false)) // see comment multiplications forward path
                        tmp += abs(coeffa_si[j]) + 1;
                    else
                        tmp += abs(coeffa_si[j]);
                }

                if(tmp == 0) // avoid log2(0)
                    msbRegisterFeedbackOut[i] = msbInt - shifta + 1;
                else
                    msbRegisterFeedbackOut[i] = msbInt - shifta + ceil(log2((abs(tmp))));

                lsbRegisterFeedbackOut[i] = lsbInt - shifta;
                REPORT(INFO,
                       "msbRegisterFeedbackOut[" << i << "]: \t" << msbRegisterFeedbackOut[i]
                                                 << ", lsbRegisterFeedbackOut["
                                                 << i << "]: \t\t" << lsbRegisterFeedbackOut[i])
            }

            // Additions, first addition is special case, computation depends on register size
            // Output from an addition is equal to the input of the succeeding register and registers values are already computed
            // exception for add0, see comment for addition 0 below
            for (int i = m - 1; i > 0; i--) {
                msbAdditionsFeedbackOut[i] = msbRegisterFeedbackOut[i -
                                                                    1]; // no need for special case handling power of two since register msb is used
                lsbAdditionsFeedbackOut[i] = lsbInt - shifta;
                REPORT(INFO,
                       "msbAdditionsFeedbackOut[" << i << "]: \t" << msbAdditionsFeedbackOut[i]
                                                  << ", lsbAdditionsFeedbackOut["
                                                  << i << "]: \t" << lsbAdditionsFeedbackOut[i])
            }

            /*
             * Addition 0
             * The computation for the msb and word sizes in the forward and feedback path is very conservative.
             * This may lead to the fact, that one input word size of the adder is larger than the result
             * The maximum word size is given by msbIn + wcpgBitGain, so the maximum output of the adder is
             * msbIn + wcpgBitGain + shifta since it is the latest operation (msbInt = msbIn + wcpgBitGain)
             * wcpg determines maximum word size after scalea
             */

            msbAdditionsFeedbackOut[0] = msbIn + wcpgBitGain;
            lsbAdditionsFeedbackOut[0] = min(lsbInt - shifta, lsbScaleBOut);
            REPORT(INFO,
                   "msbAdditionsFeedbackOut[0]: \t" << msbAdditionsFeedbackOut[0] << ", lsbAdditionsFeedbackOut[0]: \t"
                                                    << lsbAdditionsFeedbackOut[0])

            // truncationAfterScaleA
            msbTruncationIntOut = msbInt;
            lsbTruncationIntOut = lsbInt;
            REPORT(INFO,
                   "msbTruncationIntOut: \t\t" << msbTruncationIntOut << ", lsbTruncationIntOut: \t\t"
                                               << lsbTruncationIntOut)

            // ScaleA, only shift fixed point, no code generation
            msbScaleAOut = msbTruncationIntOut - (int) shifta;
            lsbScaleAOut = lsbTruncationIntOut - (int) shifta;
            REPORT(INFO, "msbScaleAOut: \t\t\t" << msbScaleAOut << ", lsbScaleAOut: \t\t\t" << lsbScaleAOut)
        }
            addInput("X", wIn, true);
            addOutput("Result", wOut);
            REPORT(INFO, "input size (wIn): " << wIn)
            REPORT(INFO, "wordsize size scaleb: " << msbScaleBOut-lsbScaleBOut+1)
            REPORT(INFO, "output size (wOut): " << wOut)

        // ################################# Sign computation for adders #########################################
        vector<adderSigns> adderFeedbackpPath;
        vector<adderSigns> adderForwardPath;

        for(uint32_t i = 0; i < n; i++) // initialize adder default with + for use in plain
        {
            adderSigns tmp;
            adderForwardPath.push_back(tmp);
        }

        for(uint32_t i = 0; i < m; i++) // initialize adder default with + for use in plain
        {
            adderSigns tmp;
            adderFeedbackpPath.push_back(tmp);
        }

        if(method=="multiplierless")
        {
            /*
             * Fedback path
             * If a coefficient in feedbackpath is positive, the result of the multiplication is a negative number (x*(-an))
             * In this case the inputs of the succeeding adder (default both "+") must changed so that the right input is "-"
             */

            for (uint32_t i = 1; i < m; i++)
                if (coeffa_si[i - 1] >= 0)
                    adderFeedbackpPath.at(i).signRight = "-";

            if (coeffa_si[m - 1] >= 0)
                adderFeedbackpPath.at(m - 1).signLeft = "-";

            for(uint32_t i = 0; i < m; i++)
                REPORT(LIST, "feedback add" << i << " left: " << adderFeedbackpPath.at(i).signLeft << " right: " << adderFeedbackpPath.at(i).signRight)

            /*
             * Forward path
             * If coefficient in forward path is negativ, then the result of the multiplication is negative
             * so change input sign of adder to "-"
             */

            for (uint32_t i = 0; i < n-1; i++)
                if (coeffb_si[i] <= 0)
                    adderForwardPath.at(i).signRight = "-";

            if (coeffb_si[n - 1] <= 0)
                adderForwardPath.at(n - 1).signLeft = "-";

            /*
             * It can happen, that both inputs are "-" which shouldnt
             * If that happens, change both "-" inputs to "+" and set the sign of the next higher adder to "-"
             */

            preventDoubleNegativeInputsAdder(m, &adderFeedbackpPath);
            preventDoubleNegativeInputsAdder(n, &adderForwardPath);

            for(uint32_t i = 0; i < m; i++)
                REPORT(LIST, "feedback add" << i << " left: " << adderFeedbackpPath.at(i).signLeft << " right: " << adderFeedbackpPath.at(i).signRight)

        }

        // ################################# TESTBENCH INIT #########################################
        // In VHDL there are for n coefficients n-1 registers, but in testbench register0 stores the input x
        // so there are n registers for n coefficients. So register0 has word length of msbIn msbIn - lsbIn +1

        //mpz_init2(registerForwardTB_mpz[0] , msbIn - lsbIn + 1);
        mpz_init(registerForwardTB_mpz[0]);
        mpz_set_si(registerForwardTB_mpz[0], 0);

        for(uint32_t i = 1; i < n; i++)
        {
            //mpz_init2(registerForwardTB_mpz[i], msbRegisterForwardOut[i-1] - lsbForwardPath + 1); // mpz_reg[i] has msb of reg[i-1] since reg0 in testbench is emulated as input
            mpz_init(registerForwardTB_mpz[i]);
            mpz_set_si(registerForwardTB_mpz[i], 0); // necessary of 0 by default?
        }

        for(uint32_t i = 0; i < m+1; i++)
        {
            //mpz_init2(registerFeedbackTB_mpz[i], msbRegisterFeedbackOut[i] - lsbRegisterFeedbackOut[i] + 1 + guardBits);
            mpz_init(registerFeedbackTB_mpz[i]);
            mpz_set_si(registerFeedbackTB_mpz[0], 0); // necessary of 0 by default?
        }

        // Initialisations for the emulate (faithfully rounded)
        xHistory  = (mpfr_t*) malloc(n * sizeof(mpfr_t));
        yHistory  = (mpfr_t*) malloc((m+2) * sizeof(mpfr_t)); // We need to memorize the m previous y, and the current output. Plus one because it helps debugging

        hugePrec = 10*(1+msbOutIIR+-lsbOut+guardBits);
        currentIndex=0x0FFFFFFFFFFFFFFFUL; // it will be decremented, let's start from high

        for (uint32_t i = 0; i<m+2; i++)
        {
            mpfr_init2 (yHistory[i], hugePrec);
            mpfr_set_d(yHistory[i], 0.0, GMP_RNDN);
        }
        for (uint32_t i=0; i<n; i++) {
            mpfr_init2 (xHistory[i], hugePrec);
            mpfr_set_d(xHistory[i], 0.0, GMP_RNDN);
        }

        // faithful testbench
        for(uint32_t i = 0; i < n; i++)
            mpfr_div_2ui(coeffb_mp_f_scaled[i], coeffb_mp_f[i], shiftb, GMP_RNDN);

        for(uint32_t i = 0; i < m; i++)
            mpfr_div_2ui(coeffa_mp_f_scaled[i], coeffa_mp_f[i], shifta, GMP_RNDN);

        // ################################# CODE GENERATION FORWARD PATH #########################################
        /*
         *
         *  Multiplications
         *  TEMP_SIGN_FOR_RESIZE adds one dedicated bit for the resize operation. It is necessary, since the word size
         *  computation is made for positive integers. Additional space for the sign bit can be obtained by
         *  increasing the msbIn by 1. Since VHDL computed the resulting word size of a multiplication by adding the
         *  word sizes of both factors, additional scaling at the end of the line is necessary to resize the word to the
         *  originated word size (-1 there is mandatory)
         */

        if(method == "plain")
        {
            // multiplications
            for (uint32_t i = 0; i < n; i++) {
                vhdl << declare(join("forMul", i), msbMultiplicationsForwardOut[i] - lsbForwardPath +1)
                     << " <= std_logic_vector(resize(signed(X) * (signed(std_logic_vector(to_signed(" << coeffb_si[i] << ","
                     << ceil(log2(abs(coeffb_si[i]) + TEMP_SIGN_FOR_RESIZE)) + TEMP_SIGN_FOR_RESIZE << ")))),"<<msbMultiplicationsForwardOut[i] - lsbForwardPath +1<<"));" << endl ;
            }
        }

        else if (method=="multiplierless")
        {
            vector<string> outportmapsBuilder;
            stringstream outPortMaps;
            for (uint32_t i = 0; i < n; i++)
            {
                if (FixIIRShiftAdd::getIndexCoeff(coeffb_si, n, coeffb_si[i]) == i) // first time occurence of element, if the result if getIndexCoeff() is not i it means that the element was found earlier in the array, hence duplicate
                {
                    declare(join("forMul", i), msbMultiplicationsForwardOut[i] - lsbForwardPath +1);
                    if(coeffb_si[i] != 0)
                        outportmapsBuilder.push_back(join("R_c", abs(coeffb_si[i]), "=>forMul", i));
                }
                else // duplicate detected, set the duplicate to the signal where it was first found, so if coeff_b[0] and coeff_b[2] are equal, forMul2 is set to forMul0
                {
                    REPORT(LIST, "duplicate detected (" << coeffb_si[i] << ")")
                    vhdl << declare(join("forMul", i), msbMultiplicationsForwardOut[i] - lsbForwardPath + 1) << "<= std_logic_vector(forMul" << FixIIRShiftAdd::getIndexCoeff(coeffb_si, n, coeffb_si[i]) << ");" <<  endl;
                }
            }

            // add comma between signals in port map. Must be done separatly since leading coefficients can be 0 what doesnt result in a signal
            for(uint32_t i = 0; i < outportmapsBuilder.size(); i++)
            {
                outPortMaps << outportmapsBuilder.at(i);
                if(i != outportmapsBuilder.size()-1)
                    outPortMaps << ",";
            }

            for (uint32_t i = 0; i < n; i++)
                if(coeffb_si[i] == 0)
                    vhdl << "forMul" << i << " <= (others => '0');" << endl;

            stringstream parameters;
            parameters << "wIn=" << wIn << " graph=" << graphb;
            string inPortMaps = "X0=>X";
            REPORT(LIST, "outPortMaps (forward): " << outPortMaps.str())
            cout << "parameters: " << parameters.str() << endl;
            newInstance("IntConstMultShiftAdd", "IntConstMultShiftAddComponentForward", parameters.str(), inPortMaps,
                        outPortMaps.str());

            REPORT(LIST, "IntConstMultShiftAdd for forward path created")
        }

        vhdl << endl;

        // additions

        string forwardpathOperandLeft = "forAdd";
        string forwardpathOperandRight = "forRegOut";

        for (uint32_t i = 0; i < n - 1; i++)
        {
            if(adderForwardPath.at(i).signLeft == "-" && adderFeedbackpPath.at(i).signRight == "+") // treat special case, see comments addition feedbackpath
            {
                vhdl << declare(join(forwardpathOperandRight, i), msbAdditionsForwardOut[i] - lsbForwardPath + 1)
                     << " <= std_logic_vector(resize(signed("
                     << join(forwardpathOperandLeft, i) << "), " << msbAdditionsForwardOut[i] - lsbForwardPath + 1
                     << ") " << adderForwardPath.at(i).signLeft << " signed(" << join("forMul", i) << "));"
                     << endl;
            }
            else {
                vhdl << declare(join(forwardpathOperandLeft, i), msbAdditionsForwardOut[i] - lsbForwardPath + 1)
                     << " <= std_logic_vector(resize(signed("
                     << join(forwardpathOperandRight, i) << "), " << msbAdditionsForwardOut[i] - lsbForwardPath + 1
                     << ") " << adderForwardPath.at(i).signRight << " signed(" << join("forMul", i) << "));"
                     << endl;
            }
        }

        vhdl << endl;
        // registers
        for (uint32_t i = 0; i < n - 2; i++) {
            newInstance("ShiftReg", join("forReg", i), join("n=1 reset=2 w=", msbRegisterForwardOut[i] - lsbForwardPath + 1),
                        join("X=>forAdd", i + 1), join("Xd1=>forRegOut", i));
        }
        // treat special case for regN-1 since it has different input (from mult not from add)
        newInstance("ShiftReg", join("forReg", n - 2), join("n=1 reset=2 w=", msbRegisterForwardOut[n-2] - lsbForwardPath + 1),
                    join("X=>forMul", n - 1), join("Xd1=>forRegOut", n - 2));
        vhdl << endl;

        if(!isFIR) {
            // ################################# CODE GENERATION BACKPATH #########################################
            // truncation, truncate to lsbInt
            vhdl << declare("truncationAfterScaleA", msbInt - lsbInt + 1) << " <= feedbackAdd0("
                 << msbInt - lsbInt + 1 + (lsbInt - lsbScaleAOut) - 1 << " downto " << lsbInt - lsbScaleAOut << ");"
                 << endl;

            // multiplications, see comment multiplications forward path for explanation of TEMP_SIGN_FOR_RESIZE
            if (method == "plain") {
                for (uint32_t i = 0; i < m; i++) {
                    vhdl << declare(join("feedbackMul", i),
                                    msbMultiplicationsFeedbackOut[i] - lsbMultiplicationsFeedbackOut[i] + 1)
                         << " <= std_logic_vector(resize(signed(truncationAfterScaleA) * (signed(std_logic_vector(to_signed("
                         << -coeffa_si[i] << ","
                         << ceil(abs(log2((abs(coeffa_si[i]) + TEMP_SIGN_FOR_RESIZE)))) + TEMP_SIGN_FOR_RESIZE
                         << "))))," << msbMultiplicationsFeedbackOut[i] - lsbMultiplicationsFeedbackOut[i] + 1 << "));"
                         << endl;
                }
            } else if (method == "multiplierless") {
                stringstream outPortMaps;
                vector<string> outportmapsBuilder;
                for (uint32_t i = 0; i < m; i++) {
                    if (FixIIRShiftAdd::getIndexCoeff(coeffa_si, m, coeffa_si[i]) ==
                        i) // first time occurence of element, if the result if getIndexCoeff() is not i it means that the element was found earlier in the array, hence duplicate
                    {
                        declare(join("feedbackMul", i),
                                msbMultiplicationsFeedbackOut[i] - lsbMultiplicationsFeedbackOut[i] + 1);
                        if (coeffa_si[i] != 0)
                            outportmapsBuilder.push_back(join("R_c", abs(coeffa_si[i]), "=>feedbackMul", i));
                    } else // duplicate detected, set the duplicate to the signal where it was first found, so if coeff_b[0] and coeff_b[2] are equal, feedbackMul is set to feedbackMul
                    {
                        REPORT(LIST, "duplicate detected (" << coeffa_si[i] << ")")
                        vhdl << declare(join("feedbackMul", i),
                                        msbMultiplicationsFeedbackOut[i] - lsbMultiplicationsFeedbackOut[i] + 1)
                             << "<= std_logic_vector(feedbackMul"
                             << FixIIRShiftAdd::getIndexCoeff(coeffa_si, m, coeffa_si[i]) << ");" << endl;
                    }
                }

                for (uint32_t i = 0; i < m; i++)
                    if (coeffa_si[i] == 0)
                        vhdl << "feedbackMul" << i << " <= (others => '0');" << endl;

                // add comma between signals in port map. Must be done separatly since leading coefficients can be 0 what doesnt result in a signal
                for (uint32_t i = 0; i < outportmapsBuilder.size(); i++) {
                    REPORT(LIST, "outputmapsbuilder: " << outportmapsBuilder.at(i))
                    outPortMaps << outportmapsBuilder.at(i);
                    if (i != outportmapsBuilder.size() - 1)
                        outPortMaps << ",";
                }

                stringstream parameters;
                parameters << "wIn=" << msbScaleAOut - lsbScaleAOut + 1 << " graph=" << grapha;
                string inPortMaps = "X0=>truncationAfterScaleA";
                REPORT(LIST, "outPortMaps (feedback): " << outPortMaps.str())
                newInstance("IntConstMultShiftAdd", "IntConstMultShiftAddComponentFeedback", parameters.str(),
                            inPortMaps,
                            outPortMaps.str());

                REPORT(DETAILED, "IntConstMultShiftAdd for feedback path created")
            }

            /*
             * Additions
             * addition0 ist special case since it has scaleB as input and the other input is possible too large bc of the
             * conservative computation in the feedback path (see comment msb computation for this addition)
             * --> feedbackRegOut0 has to be rescaled to max output of adder (later it can be optimized to mInt + wcpgGain - (msbScaleBOut - msbIn) + shifta
             * both inputs of the additions must have same amount of lsb-bits
             *
             * There are 3 cases for the signs. left "+" and right "+", left "+" and right "-" and left "-" and right "+". The possible fourth case is handled
             * in preventDoubleNegativeInputsAdder. Left "-" and right "+" is special case and must be treated differently since -a+b ist not feasible using
             * IntConstMultShift, so use b-a.
             */

            string operandNameFeedbackLeft = "feedbackRegOut"; // except feedbackAdd0
            string operandNameFeedbackRight = "feedbackMul"; // except feedbackAdd0

            int64_t numLsbBitsDifference = abs(
                    lsbRegisterFeedbackOut[0] - lsbScaleBOut); // determine how often to shift
            REPORT(DEBUG, "numLsbBitsDifference: " << numLsbBitsDifference)
            if (lsbScaleBOut < lsbRegisterFeedbackOut[0]) {
                /* lsbScaleBOut has more lsb bits, so add lsb-bits to regOut0
                 * resize and shifting (to the left), scale this value to the size of the output of feedbackAdd0, since this value cannot be larger
                 * the word size cannot be larger than the word size of the adder (see comment msb computation for this addition)
                 */
                wordsizefeedbackAdd0ResizedReg0 = min(msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1,
                                                      msbRegisterFeedbackOut[0] - lsbRegisterFeedbackOut[0] + 1 +
                                                      numLsbBitsDifference);
                vhdl << declare("feedbackAdd0ResizedReg0", wordsizefeedbackAdd0ResizedReg0)
                     << " <= std_logic_vector(shift_left(resize(signed(feedbackRegOut0), "
                     << wordsizefeedbackAdd0ResizedReg0 << ")," << numLsbBitsDifference << "));" << endl;
                if (adderFeedbackpPath.at(0).signLeft == "-" &&
                    adderFeedbackpPath.at(0).signRight == "+") // treat special case and switch operands and change sign
                {
                    vhdl << declare("feedbackAdd0", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1)
                         << " <= std_logic_vector(resize(signed(feedbackAdd0ResizedReg0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << ") "
                         << adderFeedbackpPath.at(0).signLeft << " signed(forAdd0));" << endl;
                } else {
                    vhdl << declare("feedbackAdd0", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1)
                         << " <= std_logic_vector(resize(signed(forAdd0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << ") "
                         << adderFeedbackpPath.at(0).signRight << " signed(feedbackAdd0ResizedReg0));" << endl;
                }
            } else if (lsbScaleBOut > lsbRegisterFeedbackOut[0]) {
                //vhdl << declare("feedbackAdd0ResizedScaleB", msbScaleBOut - lsbScaleBOut + numLsbBitsDifference + 1) << " <= std_logic_vector(shift_left(signed(scaleBOut(" << msbScaleBOut - lsbScaleBOut +1  << " downto 0))), to_integer(to_signed(1, 2)));" << endl;
                wordsizefeedbackAdd0ResizedScaleB = min(msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1,
                                                        msbScaleBOut - lsbScaleBOut + 1 + numLsbBitsDifference);
                vhdl << declare("feedbackAdd0ResizedScaleB", wordsizefeedbackAdd0ResizedScaleB)
                     << " <= std_logic_vector(shift_left(resize(signed(forAdd0), " << wordsizefeedbackAdd0ResizedScaleB
                     << ")," << numLsbBitsDifference << "));" << endl;
                if (adderFeedbackpPath.at(0).signLeft == "-" &&
                    adderFeedbackpPath.at(0).signRight == "+") // treat special case and switch operands and change sign
                {
                    vhdl << declare("feedbackAdd0", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1)
                         << " <= std_logic_vector(resize(signed(feedbackRegOut0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << ") "
                         << adderFeedbackpPath.at(0).signLeft << " resize(signed(feedbackAdd0ResizedScaleB), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << "));" << endl;
                } else {
                    vhdl << declare("feedbackAdd0", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1)
                         << " <= std_logic_vector(resize(signed(feedbackAdd0ResizedScaleB), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << ") "
                         << adderFeedbackpPath.at(0).signRight << " resize(signed(feedbackRegOut0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << "));" << endl;
                }
            } else // if both inputs are well aligned, use regular inputs
            {
                if (adderFeedbackpPath.at(0).signLeft == "-" &&
                    adderFeedbackpPath.at(0).signRight == "+") // treat special case and switch operands and change sign
                {
                    vhdl << declare("feedbackAdd0", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1)
                         << " <= std_logic_vector(resize(signed(feedbackRegOut0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << ") "
                         << adderFeedbackpPath.at(0).signLeft << " resize(signed(forAdd0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << "));" << endl;
                } else {
                    vhdl << declare("feedbackAdd0", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1)
                         << " <= std_logic_vector(resize(signed(forAdd0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << ") "
                         << adderFeedbackpPath.at(0).signRight << " resize(signed(feedbackRegOut0), "
                         << msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 << "));" << endl;
                }
            }

            for (uint32_t i = 1; i < m; i++) {
                if (adderFeedbackpPath.at(i).signLeft == "-" &&
                    adderFeedbackpPath.at(i).signRight == "+") // treat special case and switch operands and change sign
                {
                    vhdl << declare(join("feedbackAdd", i), msbAdditionsFeedbackOut[i] - lsbAdditionsFeedbackOut[i] + 1)
                         << " <= std_logic_vector(resize(signed("
                         << join(operandNameFeedbackRight, i - 1) << "), "
                         << msbAdditionsFeedbackOut[i] - lsbAdditionsFeedbackOut[i] + 1 << ") "
                         << adderFeedbackpPath.at(i).signLeft
                         << " signed(" << join(operandNameFeedbackLeft, i) << "));"
                         << endl;
                } else {
                    vhdl << declare(join("feedbackAdd", i), msbAdditionsFeedbackOut[i] - lsbAdditionsFeedbackOut[i] + 1)
                         << " <= std_logic_vector(resize(signed("
                         << join(operandNameFeedbackLeft, i) << "), "
                         << msbAdditionsFeedbackOut[i] - lsbAdditionsFeedbackOut[i] + 1 << ") "
                         << adderFeedbackpPath.at(i).signRight
                         << " signed(" << join(operandNameFeedbackRight, i - 1) << "));"
                         << endl;
                }
            }

            newInstance("ShiftReg", join("feedbackReg", m - 1),
                        join("n=1 reset=1 w=", msbRegisterFeedbackOut[m - 1] - lsbRegisterFeedbackOut[m - 1] + 1),
                        join("X=>feedbackMul", m - 1), join("Xd1=>feedbackRegOut", m - 1));

            for (uint32_t i = 0; i < m - 1; i++) {
                newInstance("ShiftReg", join("feedbackReg", i),
                            join("n=1 reset=1 w=", msbRegisterFeedbackOut[i] - lsbRegisterFeedbackOut[i] + 1),
                            join("X=>feedbackAdd", i + 1), join("Xd1=>feedbackRegOut", i));
            }

            // determine bit position for the bit that has to be added
            uint16_t bitPositionRoundingBit = abs(lsbAdditionsFeedbackOut[0] - lsbInt) + guardBits - 1;
            // cutting bit so it survises truncation
            vhdl << declare("faithfulRoundingBitValue", 1) << " <= feedbackAdd0(" << bitPositionRoundingBit
                 << " downto " << bitPositionRoundingBit << " ); " << endl;

            // shift feedback0 to the right, resp. cut two more bits, guard bits also considered (truncation)
            vhdl << declare("endresult", msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0] + 1 -
                                         abs(lsbAdditionsFeedbackOut[0] - lsbInt - guardBits)) << " <= feedbackAdd0("
                 << (msbAdditionsFeedbackOut[0] - lsbAdditionsFeedbackOut[0]) << " downto "
                 << abs(lsbAdditionsFeedbackOut[0] - lsbInt) + guardBits << ");" << endl;
            //vhdl << "Result <= std_logic_vector(resize(signed(endresult)," << wOut << "));" << endl; // not rounded result

            vhdl << "Result <= std_logic_vector(resize(signed(endresult + faithfulRoundingBitValue)," << wOut << "));"
                 << endl; // faithfully rounded result
        }

        // resulting signal for FIR filter
        if(isFIR)
        {
            if(lsbOut - lsbScaleBOut > 0) // if there are more lsbScaleBOut bits than lsbOut,
            {
                REPORT(DEBUG, "cut " << lsbOut - lsbScaleBOut << " bit. Range from " << msbScaleBOut-lsbScaleBOut+1-1 << " to: " << lsbOut - lsbScaleBOut-1 << endl);
                uint16_t bitPositionRoundingBit = msbScaleBOut-lsbScaleBOut-1;
                vhdl << declare("faithfulRoundingBitValue", 1) << " <= forAdd0(" << bitPositionRoundingBit << " downto " << bitPositionRoundingBit << " ); " << endl;
                // no +1 in lsbOut - lsbScaleBOut and the  other since for downto 1 must be suctracted again
                vhdl << declare("resultNotRounded", wOut) << " <= std_logic_vector(forAdd0)(" << msbScaleBOut-lsbScaleBOut << "  downto " << lsbOut - lsbScaleBOut << "); " << endl;
                vhdl << "Result <= std_logic_vector(resultNotRounded+faithfulRoundingBitValue);" << endl;
            }
            else if(lsbOut - lsbScaleBOut < 0)
            {
                REPORT(DEBUG, "pad " << abs(lsbOut - lsbScaleBOut) << "bit on right side")
                vhdl << "Result <= std_logic_vector(shift_left(resize(signed(forAdd0), " << wOut << "), " << abs(lsbOut - lsbScaleBOut) << ")); " << endl;
            }
            else
            {
                REPORT(DEBUG, "proper alignment");
                vhdl << "Result <= std_logic_vector(forAdd0);" << endl;
            }
        }
    }

    FixIIRShiftAdd::~FixIIRShiftAdd()
    {
        delete(coeffb_d);
        delete(coeffa_d);
        delete(coeffa_si);
        delete(coeffb_si);
        delete(coeffa_d_wcpg);
        delete(coeffb_d_wcpg);
        delete(msbRegisterForwardOut);
        delete(msbAdditionsForwardOut);
        delete(msbMultiplicationsForwardOut);
        delete(msbRegisterFeedbackOut);
        delete(lsbRegisterFeedbackOut);
        delete(msbMultiplicationsFeedbackOut);
        delete(lsbMultiplicationsFeedbackOut);
        delete(msbAdditionsFeedbackOut);
        delete(lsbAdditionsFeedbackOut);

        mpz_clear(resultForwardPathTB_mpz);
        mpz_clear(feedbackAdd0);

        for (uint32_t i = 0; i < n; i++)
        {
            mpfr_clear(coeffb_mp_f[i]);
            mpz_clear(registerForwardTB_mpz[i]);
        }

        for (uint32_t i = 0; i < m; i++)
        {
            mpfr_clear(coeffa_mp_f[i]);
            mpz_clear(registerForwardTB_mpz[i]);
        }
    };

    void FixIIRShiftAdd::emulate(TestCase * tc) {
        // my old testbench (for bit level testing)
#if 0
        mpz_class sx_ = tc->getInputValue("X");            // get the input bit vector as an integer
        sx_ = bitVectorToSigned(sx_, msbIn - lsbIn + 1);    // convert it to a signed mpz_class
        mpz_t x_;
        //mpz_init2(x_, msbIn-lsbIn+1); // restrict size
        mpz_init(x_);
        mpz_set_si(x_, sx_.get_si()); // convert input to mpz_t
        //mpz_init2(resultForwardPathTB_mpz, msbScaleBOut - lsbScaleBOut + 1);
        mpz_init(resultForwardPathTB_mpz);
        mpz_t coeffResult;
        //mpz_init2(coeffResult, msbScaleBOut - lsbInt +1); // use this as word width since it cannot be larger than that, can be optimzed but not neccessary
        mpz_init(coeffResult);

        mpz_set_si(resultForwardPathTB_mpz, 0);
        mpz_set_si(registerForwardTB_mpz[0], mpz_get_si(x_));

        for (int i = 0; i < n; i++) {
            // no add_mul_si() for mpz_t
            mpz_mul_si(coeffResult, registerForwardTB_mpz[i], coeffb_si[i]); // coeffb_d[i] * registerforwardTB[i];
            mpz_add(resultForwardPathTB_mpz, resultForwardPathTB_mpz, coeffResult);
        }

        for (int i = n - 1; i >= 1; i--) {
            mpz_set_si(registerForwardTB_mpz[i], mpz_get_si(registerForwardTB_mpz[i - 1]));
        }

        // ######################################## TESTBENCH BACK PATH ###############################################

        mpz_init(feedbackAdd0);
        mpz_set_si(feedbackAdd0, 0);

        // align both inputs for addition feedbackadd0
        // input from forward path is affected by scaleb, input from feedbackpath is affected by scalea and guard bits

        int numLsbBitsDifferenceTB = abs(lsbRegisterFeedbackOut[0] - lsbScaleBOut);
        REPORT(DEBUG, "numLsbBitsDifferenceTB: " << numLsbBitsDifferenceTB)
        if (lsbScaleBOut < lsbRegisterFeedbackOut[0]) {
            // lsbScaleBOut has more lsb bits, so add lsb-bits to regOut0
            mpz_mul_2exp(registerFeedbackTB_mpz[0], registerFeedbackTB_mpz[0], numLsbBitsDifferenceTB);
        } else if (lsbScaleBOut > lsbRegisterFeedbackOut[0]) {
            // output of forwardpath has more lsb, so add lsb bits to output of forward path (result_forwardpath)
            mpz_mul_2exp(resultForwardPathTB_mpz, resultForwardPathTB_mpz, numLsbBitsDifferenceTB);
        }

        // add 0
        mpz_add(feedbackAdd0, registerFeedbackTB_mpz[0], resultForwardPathTB_mpz);
        REPORT(DEBUG,
               "feedbackAdd0: " << feedbackAdd0 << ", registerFeedbackTB_mpz[0]: " << registerFeedbackTB_mpz[0]
                                << ", resultForwardPathTB_mpz: " << resultForwardPathTB_mpz)

        // set result and shift to align to lsbOut since lsbAdditionsFeedbackOut[0] can different lsb value (only more lsb)
        mpz_t result_mpz;
        mpz_init_set(result_mpz, feedbackAdd0);

        /*
         * faithful rounding
         * determine the bit which must be added (do it before truncating)
         * Shift result (resop. temp result value) to right so that this bit is the lsb
         * mod2 to determine if 0 or 1
         * if 1, add to result
         */
        // determine most-left truncated bit position and value
        //REPORT(LIST, "msb: " << msbAdditionsFeedbackOut[0] << ", lsb: " << lsbAdditionsFeedbackOut[0])
        //REPORT(LIST, "most-left truncated bit position: " << abs(lsbAdditionsFeedbackOut[0]-lsbInt) + guardBits-1) // if there are 9 truncated bits, bit 8..0 are truncated, so add bit 8

        mpz_t resultTempForFaithfulRounding;
        mpz_init_set(resultTempForFaithfulRounding, feedbackAdd0);
        mpz_t faithfulRoundingBitValue;
        mpz_init_set_ui(faithfulRoundingBitValue, 2); // set to two for debugging purposes
        // shift right so that most-left truncated bit is bot 0
        mpz_fdiv_q_2exp(resultTempForFaithfulRounding, resultTempForFaithfulRounding,
                        abs(lsbAdditionsFeedbackOut[0] - lsbInt) + guardBits - 1);
        mpz_mod_ui(faithfulRoundingBitValue, resultTempForFaithfulRounding, 2); // use mod2 to determine if 0 or 1
        //REPORT(LIST, "this value must be added: " << mpz_get_si(faithfulRoundingBitValue))

        mpz_fdiv_q_2exp(result_mpz, result_mpz,
                        abs(lsbAdditionsFeedbackOut[0] - lsbInt)); // realign result, truncation to wout
        mpz_fdiv_q_2exp(result_mpz, result_mpz, guardBits); // shift back bc of guard bits
        REPORT(DEBUG,
               "RESULT: " << result_mpz << ", shifted: " << (mpz_get_si(result_mpz)) / (int) pow(2, abs(lsbOut)))
        mpz_class result_class(result_mpz);
        //tc->addExpectedOutput("Result", result_class);

        // add this bit value to the truncated value to get the faithfully rounded result
        mpz_add(result_mpz, result_mpz, faithfulRoundingBitValue);
        mpz_class result_class_fathfully(result_mpz);
        tc->addExpectedOutput("Result", result_class_fathfully); // use this for faithful rounded result

        REPORT(DEBUG, "")

        // scaling (only in code, msb and lsb is altered)
        // truncation
        mpz_t truncated;
        mpz_init(truncated);
        // current lsb is lsb from feedbackAdd0

        REPORT(DEBUG, "lsbAdditionsFeedbackOut: " << lsbAdditionsFeedbackOut[0] << ", scaleAOut: " << lsbScaleAOut
                                                  << ", lsbInt: " << lsbInt);
        // truncate to lsbInt
        REPORT(DEBUG, "truncate " << feedbackAdd0 << " by " << abs(lsbScaleAOut - lsbInt) << " bits")
        mpz_fdiv_q_2exp(truncated, feedbackAdd0, abs(lsbScaleAOut - lsbInt)); // dont use tdiv, use fdiv!

        REPORT(DEBUG, "truncated: " << truncated)
        REPORT(DEBUG, "add0: " << feedbackAdd0)
        for (int i = 0; i < m; i++) REPORT(DEBUG, "register" << i << ": " << registerFeedbackTB_mpz[i])

        for (int i = 0; i < m; i++) {
            mpz_t coeffR;
            mpz_init(coeffR);
            mpz_mul_si(coeffR, truncated, coeffa_si[i]);
            // mpz_mul_si(coeffR, coeffR, -1);
            REPORT(DEBUG, "coeffr mul" << i << ": " << mpz_get_si(coeffR))
            mpz_sub(registerFeedbackTB_mpz[i], registerFeedbackTB_mpz[i + 1], coeffR);
        }
#endif
        if (1) {
            mpz_class sx;
            mpfr_t x, s, t;

            mpfr_init2(s, hugePrec);
            mpfr_init2(t, hugePrec);
            mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0

            mpfr_init2(x, msbIn - lsbIn +1);

            sx = tc->getInputValue("X");        // get the input bit vector as an integer
            //sx = bitVectorToSigned(sx, 1 - lsbIn);                        // convert it to a signed mpz_class
            sx = bitVectorToSigned(sx, msbIn - lsbIn +1);                        // convert it to a signed mpz_class
            //REPORT(LIST, "input: " << sx)
            mpfr_set_z(x, sx.get_mpz_t(), GMP_RNDD);                // convert this integer to an MPFR; this rounding is exact
            mpfr_div_2si(x, x, -lsbIn, GMP_RNDD);                        // multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact
            //REPORT(LIST, "x incl. lsb bits: " << mpfr_get_d(x, GMP_RNDN))
            mpfr_set(xHistory[currentIndex % n], x, GMP_RNDN); // exact

            //REPORT(LIST, "output variable: " << mpfr_get_d(s, GMP_RNDN))
            // TODO CHECK HERE
            for (uint32_t i = 0; i < n; i++) {
                mpfr_mul(t, xHistory[(currentIndex + i) % n], coeffb_mp_f_scaled[i],
                         GMP_RNDZ);                    // Here rounding possible, but precision used is ridiculously high so it won't matter
                //REPORT(LIST, "forward: " << mpfr_get_d(t, GMP_RNDN)<< " = " << mpfr_get_d(xHistory[(currentIndex+i)%n], GMP_RNDN) << " * " << mpfr_get_d(coeffb_mp_f_scaled[i], GMP_RNDN))
                mpfr_add(s, s, t, GMP_RNDN);                            // same comment as above
            }
            if(!isFIR)
            {
                for (uint32_t i = 0; i < m; i++) {
                    mpfr_mul(t, yHistory[(currentIndex + i + 1) % (m + 2)], coeffa_mp_f_scaled[i],
                             GMP_RNDZ);                    // Here rounding possible, but precision used is ridiculously high so it won't matter
                    //REPORT(LIST, "feedback: " << mpfr_get_d(t, GMP_RNDN)<< " = " << mpfr_get_d(yHistory[(currentIndex +i+1)%(m+2)], GMP_RNDN) << " * " << mpfr_get_d(coeffa_mp_f_scaled[i], GMP_RNDN))
                    mpfr_sub(s, s, t, GMP_RNDZ);                            // same comment as above
                }
            //REPORT(LIST, "output variable: " << mpfr_get_d(s, GMP_RNDZ))
            mpfr_set(yHistory[(currentIndex + 0) % (m + 2)], s, GMP_RNDN);
            }

#if 0// debugging the emulate
            cout << "x=" << 	mpfr_get_d(xHistory[currentIndex % n], GMP_RNDN);
            cout << " //// y=" << 	mpfr_get_d(s,GMP_RNDN) << "  ////// ";
            for (uint32_t i=0; i< n; i++)		{
                cout << "  x" << i<< "c" << i<<  "=" <<
                    mpfr_get_d(xHistory[(currentIndex+i)%n],GMP_RNDN) << "*" << mpfr_get_d(coeffb_mp_f[i],GMP_RNDN);
            }
            cout << "  // ";
            for (uint32_t i=0; i<m; i++) {
                cout <<"  ya" << i+1 << "=" <<
                    mpfr_get_d(yHistory[(currentIndex +i+1)%(m+2)],GMP_RNDN) << "*" << mpfr_get_d(coeffa_mp_f[i],GMP_RNDN);
            }
            cout << endl;

#endif

            currentIndex--;

            //	coeff		  1 2 3
            //      yh      y 0 0 0
            // now we should have in s the (exact in most cases) sum
            // round it up and down

            // debug: with this we observe if the simulation diverges
            double d = mpfr_get_d(s, GMP_RNDD);
            miny = min(d, miny);
            maxy = max(d, maxy);
            //		cout << "y=" << d <<  "\t  log2(|y|)=" << (ceil(log2(abs(d)))) << endl;

            // make s an integer -- no rounding here
            mpfr_mul_2si(s, s, -lsbOut, GMP_RNDN);
            //REPORT(LIST, "output variable scaled with lsbOut: " << mpfr_get_d(s, GMP_RNDZ))

            // We are waiting until the first meaningful value comes out of the IIR
            int bitVectorSize; // bit vector size is different for FIR since the position of the output is different
            if(!isFIR)
                bitVectorSize = msbOutIIR - lsbOut + 1;
            else
                bitVectorSize = msbScaleBOut - lsbOut + 1;

        if(!isFIR)
        {
            mpz_class rdz, ruz;
            mpfr_get_z(rdz.get_mpz_t(), s, GMP_RNDD);                    // there can be a real rounding here

#if 1 // to unplug the conversion that fails to see if it diverges further
            rdz = signedToBitVector(rdz, bitVectorSize);
            tc->addExpectedOutput("Result", rdz);

            mpfr_get_z(ruz.get_mpz_t(), s, GMP_RNDU);                    // there can be a real rounding here
            ruz = signedToBitVector(ruz, bitVectorSize);
            tc->addExpectedOutput("Result", ruz);
#endif
        }
            mpfr_clears(x, t, s, NULL);
        }
    }

    void FixIIRShiftAdd::buildStandardTestCases(TestCaseList* tcl)
    {
    if(0) {
        TestCase *tc;
        // compute the impulse response
        computeImpulseResponse();
        // Now fill with a signal that follows the sign alternance of the impulse response: this builds a worst-case signal
        miny = 0;
        maxy = 0;
        int storageOffset = n + m;
        uint32_t kmax = vanishingK - storageOffset;
        //REPORT(LIST, "kmax: " << kmax)
        for (uint32_t i = 0; i < kmax; i++) {
            mpz_class val;
#if 0
            if(yi[kmax-i]<0) {
                val = ((mpz_class(1)<<(-lsbIn)) -1) ; // 011111
            }
            else {
                val = ((mpz_class(1)<<(-lsbIn)) +1); // 100001
            }
#else // multiplying by 1 and -1 ensures no rounding error in the FIR part
            // hence the following code
            // But no observable difference...
            //val = ((mpz_class(1)<<(-lsbIn)) -1) * 9 / 10; // 011111;  *9/10 to trigger rounding errors
            val = ((mpz_class(1)<<(-lsbIn)+msbIn) -1) * 9 / 10; // 011111;  *9/10 to trigger rounding errors
            //REPORT(LIST, "val not complement: " << val )
            if (yi[kmax - i] >= 0) {
                // two's complement
                //val = ((mpz_class(1)<<(-lsbIn+1)) -1) -val +1 ; // 111111 - val + 1
                val = ((mpz_class(1)<<(-lsbIn+1)+msbIn) -1) -val +1 ; // 111111 - val + 1
                //REPORT(LIST, "twos complement: " << val )
            }

#endif
            tc = new TestCase(this);
            tc->addInput("X", val);
            emulate(tc);
            tcl->add(tc);
        }

        REPORT(0, "Filter output remains in [" << miny << ", " << maxy << "]");
    }
        if(0) {
            int largestPossibleValue = pow(2, msbIn - lsbIn + 1) - 1;
            int numTestcases = 10;

            int *testcases = (int *) malloc(numTestcases * sizeof(int));
            for (int i = 0; i < numTestcases; i++)
            {
                testcases[i] = largestPossibleValue;
            }

            TestCase *tc;
            for (int i = 0; i < numTestcases; i++) {
                tc = new TestCase(this);
                tc->addInput("X", mpz_class(testcases[i]));
                emulate(tc);
                tcl->add(tc);
            }
        }
    }

    OperatorPtr FixIIRShiftAdd::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int msbIn;
        UserInterface::parseInt(args, "msbIn", &msbIn);
        int lsbIn;
        UserInterface::parseInt(args, "lsbIn", &lsbIn);
        int lsbOut;
        UserInterface::parseInt(args, "lsbOut", &lsbOut);
        int shifta;
        UserInterface::parseInt(args, "shifta", &shifta);
        int shiftb;
        UserInterface::parseInt(args, "shiftb", &shiftb);
        double H;
        UserInterface::parseFloat(args, "H", &H);
        double Heps;
        UserInterface::parseFloat(args, "Heps", &Heps);
        int guardBits;
        UserInterface::parseInt(args, "guardbits", &guardBits);
        vector<string> coeffb;
        vector<string> coeffa;
        string in;
        UserInterface::parseString(args, "coeffb", &in);
        stringstream ssa(in);
        while( ssa.good() )	{
            string substr;
            getline( ssa, substr, ':' );
            coeffb.push_back( substr );
        }

        UserInterface::parseString(args, "coeffa", &in);
        stringstream ssb(in);
        while( ssb.good() )	{
            string substr;
            getline( ssb, substr, ':' );
            coeffa.push_back( substr );
        }

        string method;
        UserInterface::parseString(args, "method", &method);

        string grapha;
        UserInterface::parseString(args, "grapha", &grapha);

        string graphb;
        UserInterface::parseString(args, "graphb", &graphb);

        return new FixIIRShiftAdd(parentOp, target, lsbIn, lsbOut, msbIn, guardBits, coeffb, coeffa, shifta, shiftb, method, grapha, graphb, H, Heps);
    }

    TestList FixIIRShiftAdd::unitTest(int index) {
        TestList testStateList;
        vector<pair<string, string>> paramList;

        // FIR testbench
        paramList.push_back(make_pair("msbIn", "5"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:2:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "5"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-6"));
        paramList.push_back(make_pair("coeffb", "1:2:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "5"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-10"));
        paramList.push_back(make_pair("coeffb", "1:2:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "5"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:2:0"));
        paramList.push_back(make_pair("shiftb", "3"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "5"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-6"));
        paramList.push_back(make_pair("coeffb", "1:2:0"));
        paramList.push_back(make_pair("shiftb", "3"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "5"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-10"));
        paramList.push_back(make_pair("coeffb", "1:2:0"));
        paramList.push_back(make_pair("shiftb", "3"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("shiftb", "0"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "32:-32:16"));
        paramList.push_back(make_pair("shiftb", "9"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "32:-32:16"));
        paramList.push_back(make_pair("shiftb", "9"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "2:-2:1"));
        paramList.push_back(make_pair("shiftb", "3"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "2:2:1"));
        paramList.push_back(make_pair("shiftb", "3"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "2:6"));
        paramList.push_back(make_pair("shiftb", "16"));
        testStateList.push_back(paramList);
        paramList.clear();

        // IIR testbench

        if (1) {// old testbench, use again later
            paramList.push_back(make_pair("msbIn", "21"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "576:512:128:7:8"));
            paramList.push_back(make_pair("coeffa", "2097152:4194304:1048576:7340032"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "23"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "21"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "576:512:128:7:8"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "21"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "21"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "21"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            //paramList.push_back(make_pair("guardbits",  "1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "21"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            //paramList.push_back(make_pair("guardbits",  "5"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "17"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            //paramList.push_back(make_pair("guardbits",  "1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "17"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-10"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            //paramList.push_back(make_pair("guardbits",  "1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "17"));
            paramList.push_back(make_pair("lsbIn", "-6"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            //paramList.push_back(make_pair("guardbits",  "1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "17"));
            paramList.push_back(make_pair("lsbIn", "-6"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
            paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
            paramList.push_back(make_pair("shiftb", "0"));
            paramList.push_back(make_pair("shifta", "3"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            // multiplierless

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-16"));
            paramList.push_back(make_pair("lsbOut", "-16"));
            paramList.push_back(make_pair("coeffb", "32:-32:16"));
            paramList.push_back(make_pair("coeffa", "96:48"));
            paramList.push_back(make_pair("shiftb", "9"));
            paramList.push_back(make_pair("shifta", "9"));
            paramList.push_back(make_pair("graphb", "\"{{\'O\',[32],1,[1],0,5},{\'O\',[16],1,[1],0,4}}\""));
            paramList.push_back(make_pair("grapha",
                                          "\"{{\'O\',[96],2,[3],1,5},{\'O\',[48],2,[3],1,4},{\'A\',[3],1,[1],0,1,[1],0,0}}\""));
            paramList.push_back(make_pair("guardbits", "-1"));
            paramList.push_back(make_pair("method", "multiplierless"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-16"));
            paramList.push_back(make_pair("lsbOut", "-16"));
            paramList.push_back(make_pair("coeffb", "128:144:32"));
            paramList.push_back(make_pair("coeffa", "0:128"));
            paramList.push_back(make_pair("shiftb", "12"));
            paramList.push_back(make_pair("shifta", "11"));
            paramList.push_back(make_pair("graphb",
                                          "\"{{\'O\',[128],2,[1],0,7},{\'O\',[144],2,[9],1,4},{\'O\',[32],2,[1],0,5},{\'A\',[9],1,[1],0,3,[1],0,0}}\""));
            paramList.push_back(make_pair("grapha", "\"{{\'O\',[128],1,[1],0,7}}\""));
            paramList.push_back(make_pair("guardbits", "-1"));
            paramList.push_back(make_pair("method", "multiplierless"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-16"));
            paramList.push_back(make_pair("lsbOut", "-16"));
            paramList.push_back(make_pair("coeffb", "32:-32:16"));
            paramList.push_back(make_pair("coeffa", "96:48"));
            paramList.push_back(make_pair("shiftb", "9"));
            paramList.push_back(make_pair("shifta", "9"));
            paramList.push_back(make_pair("graphb", "\"{{\'O\',[32],1,[1],0,5},{\'O\',[16],1,[1],0,4}}\""));
            paramList.push_back(make_pair("grapha",
                                          "\"{{\'O\',[96],2,[3],1,5},{\'O\',[48],2,[3],1,4},{\'A\',[3],1,[1],0,1,[1],0,0}}\""));
            paramList.push_back(make_pair("method", "multiplierless"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-16"));
            paramList.push_back(make_pair("lsbOut", "-16"));
            paramList.push_back(make_pair("coeffb", "32:32:16"));
            paramList.push_back(make_pair("coeffa", "-96:48"));
            paramList.push_back(make_pair("shiftb", "9"));
            paramList.push_back(make_pair("shifta", "9"));
            paramList.push_back(make_pair("graphb", "\"{{\'O\',[32],1,[1],0,5},{\'O\',[16],1,[1],0,4}}\""));
            paramList.push_back(make_pair("grapha",
                                          "\"{{\'O\',[96],2,[3],1,5},{\'O\',[48],2,[3],1,4},{\'A\',[3],1,[1],0,1,[1],0,0}}\""));
            paramList.push_back(make_pair("method", "multiplierless"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "32:32:16"));
            paramList.push_back(make_pair("coeffa", "-96:48"));
            paramList.push_back(make_pair("shiftb", "9"));
            paramList.push_back(make_pair("shifta", "9"));
            paramList.push_back(make_pair("graphb", "\"{{\'O\',[32],1,[1],0,5},{\'O\',[16],1,[1],0,4}}\""));
            paramList.push_back(make_pair("grapha",
                                          "\"{{\'O\',[96],2,[3],1,5},{\'O\',[48],2,[3],1,4},{\'A\',[3],1,[1],0,1,[1],0,0}}\""));
            paramList.push_back(make_pair("method", "multiplierless"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-16"));
            paramList.push_back(make_pair("lsbOut", "-16"));
            paramList.push_back(make_pair("coeffb", "128:144:32"));
            paramList.push_back(make_pair("coeffa", "0:128"));
            paramList.push_back(make_pair("shiftb", "12"));
            paramList.push_back(make_pair("shifta", "11"));
            paramList.push_back(make_pair("graphb",
                                          "\"{{\'O\',[128],2,[1],0,7},{\'O\',[144],2,[9],1,4},{\'O\',[32],2,[1],0,5},{\'A\',[9],1,[1],0,3,[1],0,0}}\""));
            paramList.push_back(make_pair("grapha", "\"{{\'O\',[128],1,[1],0,7}}\""));
            paramList.push_back(make_pair("method", "multiplierless"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "128:144:32"));
            paramList.push_back(make_pair("coeffa", "0:128"));
            paramList.push_back(make_pair("shiftb", "12"));
            paramList.push_back(make_pair("shifta", "11"));
            paramList.push_back(make_pair("graphb",
                                          "\"{{\'O\',[128],2,[1],0,7},{\'O\',[144],2,[9],1,4},{\'O\',[32],2,[1],0,5},{\'A\',[9],1,[1],0,3,[1],0,0}}\""));
            paramList.push_back(make_pair("grapha", "\"{{\'O\',[128],1,[1],0,7}}\""));
            paramList.push_back(make_pair("method", "multiplierless"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            // test negative power of two thing in plain
            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-2:-2:2"));
            paramList.push_back(make_pair("coeffa", "6:3"));
            paramList.push_back(make_pair("shiftb", "3"));
            paramList.push_back(make_pair("shifta", "3"));
            paramList.push_back(make_pair("method", "plain"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "-2:-2:-2"));
            paramList.push_back(make_pair("coeffa", "6:3"));
            paramList.push_back(make_pair("shiftb", "3"));
            paramList.push_back(make_pair("shifta", "3"));
            paramList.push_back(make_pair("method", "plain"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();

            paramList.push_back(make_pair("msbIn", "-1"));
            paramList.push_back(make_pair("lsbIn", "-8"));
            paramList.push_back(make_pair("lsbOut", "-8"));
            paramList.push_back(make_pair("coeffb", "2:2:2"));
            paramList.push_back(make_pair("coeffa", "6:3"));
            paramList.push_back(make_pair("shiftb", "3"));
            paramList.push_back(make_pair("shifta", "3"));
            paramList.push_back(make_pair("method", "plain"));
            paramList.push_back(make_pair("guardbits", "-1"));
            testStateList.push_back(paramList);
            paramList.clear();
        }
        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
        paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "-576:512:128:7:8:5:6"));
        paramList.push_back(make_pair("coeffa", "2:4:1:7:2:2"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("coeffa", "4:2"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("coeffa", "4:2"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("coeffa", "256:128"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "9"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("coeffa", "256:128"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "9"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "17"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("coeffa", "4:2"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "1:0:0"));
        paramList.push_back(make_pair("coeffa", "4:2:2"));
        paramList.push_back(make_pair("shiftb", "0"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "32:-32:16"));
        paramList.push_back(make_pair("coeffa", "96:48"));
        paramList.push_back(make_pair("shiftb", "9"));
        paramList.push_back(make_pair("shifta", "9"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "32:-32:16"));
        paramList.push_back(make_pair("coeffa", "96:48"));
        paramList.push_back(make_pair("shiftb", "9"));
        paramList.push_back(make_pair("shifta", "9"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "2:-2:1"));
        paramList.push_back(make_pair("coeffa", "6:3"));
        paramList.push_back(make_pair("shiftb", "3"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-8"));
        paramList.push_back(make_pair("lsbOut", "-8"));
        paramList.push_back(make_pair("coeffb", "2:2:1"));
        paramList.push_back(make_pair("coeffa", "6:3"));
        paramList.push_back(make_pair("shiftb", "3"));
        paramList.push_back(make_pair("shifta", "3"));
        paramList.push_back(make_pair("guardbits", "-1"));
        testStateList.push_back(paramList);
        paramList.clear();

        paramList.push_back(make_pair("msbIn", "-1"));
        paramList.push_back(make_pair("lsbIn", "-16"));
        paramList.push_back(make_pair("lsbOut", "-16"));
        paramList.push_back(make_pair("coeffb", "2:6"));
        paramList.push_back(make_pair("coeffa", "6:10:6"));
        paramList.push_back(make_pair("shiftb", "16"));
        paramList.push_back(make_pair("shifta", "16"));
        paramList.push_back(make_pair("guardbits", "-1"));
        paramList.push_back(make_pair("grapha",
                                      "\"{{\'O\',[6],3,[3],2,1},{\'O\',[10],3,[5],1,1},{\'A\',[3],2,[1],0,-1,[5],1,-1},{\'A\',[5],1,[1],0,2,[1],0,0}}\""));
        paramList.push_back(make_pair("graphb", "\"{{\'O\',[2],2,[1],0,1},{\'O\',[6],2,[3],1,1},{\'A\',[3],1,[1],0,2,[-1],0,0}}\""));
        paramList.push_back(make_pair("method", "multiplierless"));

        testStateList.push_back(paramList);
        paramList.clear();

        return testStateList;
    }

    void FixIIRShiftAdd::registerFactory(){
        UserInterface::add("FixIIRShiftAdd",
                           "An Infinite Impulse Response filter generator using IntConstMultShiftAdd (optional).",
                           "FiltersEtc", // categories
                           "",
                           "msbIn(int): input most significant bit;\
                        lsbIn(int): input least significant bit;\
                        lsbOut(int): output least significant bit;\
                        H(real)=0: worst-case peak gain. if 0, it will be computed by the WCPG library;\
                        Heps(real)=0: worst-case peak gain of the feedback loop. if 0, it will be computed by the WCPG library;\
                        guardbits(int)=-1: number of guard bits for computation in recursive feedback path (-1: automatic);\
                        coeffa(string)=-1: colon-separated list of real coefficients using Sollya syntax. Example: coeffa=\"1.234567890123:sin(3*pi/8)\";\
                        coeffb(string): colon-separated list of real coefficients using Sollya syntax. Example: coeffb=\"1.234567890123:sin(3*pi/8)\";\
                        shifta(int)=-1: Num of rightshifts for coeffa (must be positive);\
                        shiftb(int): Num of rightshifts for coeffb (must be positive);\
                        method(string)=plain: plain or multiplierless;\
                        grapha(string)=emptya: graph in rpag format for coeffa;\
                        graphb(string)=emptyb: graph in rpag format for coeffb;\
                        ",
                           "",
                           FixIIRShiftAdd::parseArguments,
                           FixIIRShiftAdd::unitTest
        ) ;
    }

    int64_t FixIIRShiftAdd::getIndexCoeff(int64_t* coeff, int64_t arrayLength, int64_t val)
    {
        for(int i = 0; i < arrayLength; i++)
        {
            if(abs(coeff[i]) == abs(val)) // use abs() since output for x and -x is the same (there will be only one signal fot x and -x)
                return i;
        }
        return -1;
    }

    void FixIIRShiftAdd::preventDoubleNegativeInputsAdder(uint16_t numAdder, vector<adderSigns> *adder)
    {
        for (int i = numAdder - 1; i >= 0; i--) {
            if (adder->at(i).signRight == "-" && adder->at(i).signLeft == "-") {
                adder->at(i).signRight = "+";
                adder->at(i).signLeft = "+";
                // if signRight is already -, change signLeft
                if (adder->at(i - 1).signRight == "+")
                    adder->at(i - 1).signRight = "-";
                else
                    adder->at(i - 1).signLeft = "-";
            }
        }
    }

    void FixIIRShiftAdd::computeImpulseResponse() {
        // simulate the filter on an impulsion for long enough, until
        // double threshold = 0.5/(1<<-lsbOut);
        double threshold = 0; //soyons fous
        double epsilon=1e15; // initialize with a large value
        uint64_t k;
        //initialize the ui and yi
        for(uint32_t i=0; i<n+m; i++) {
            ui.push_back(0);
            yi.push_back(0);
        }
        ui.push_back(1); // input impulse
        yi.push_back(0);

        k=0;
        int storageOffset=n+m;

        while (epsilon>threshold) {
            // make room
            ui.push_back(0);
            yi.push_back(0);
            // compute the new y
            double y=0;
            for(uint32_t i=0; i<n; i++) {
                y += ui[storageOffset+ k-i]*coeffb_d[i] ;
            }
            for(uint32_t i=0; i<m; i++) {
                //		cout << "    k=" << k <<" i=" << i <<  "  yi[storageOffset+ k-i] =" << yi[storageOffset+ k-i] << endl;
                y -= yi[storageOffset+ k-i]*coeffa_d[i] ;
            }
            k++;
            yi[storageOffset+k] = y;

            epsilon = abs(y);
            //cout << "k=" << k << " yi=" << y << endl;
            if(k>=300000){
                REPORT(0, "computeImpulseResponse: giving up for k=" <<k << " with epsilon still at " << epsilon << ", it seems hopeless");
                epsilon=0;
            }
        }
        vanishingK=k;
        REPORT(0, "Impulse response vanishes for k=" << k);
    }

    bool FixIIRShiftAdd::isPowerOfTwo(int64_t number, bool checkNegativeNumber)
    {
        // coefficients are limited to 2^64
        // if number if power of two, it is 1 followed by zeros and (number -1) has inverted bits
        // so use bitwise and
        if(checkNegativeNumber) // negate number so it works for negative numbers only (forwardpath), feedpack path must test for positive coefficient
            number *= (-1);
        // only works for positive numbers
        if(number == 0 || number == 1)
            return false;
        return ((number & (number-1))) == 0? true : false;
    }

}