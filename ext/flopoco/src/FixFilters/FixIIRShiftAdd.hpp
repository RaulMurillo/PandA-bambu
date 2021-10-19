//
// Created by jk on 9/21/20.
//

#ifndef FLOPOCO_FIXIIRSHIFTADD_H
#define FLOPOCO_FIXIIRSHIFTADD_H

#endif //FLOPOCO_FIXIIRSHIFTADD_H

#include "Operator.hpp"
#include "utils.hpp"
#define TB_INPUT(X) (int)(X*pow(2, abs(lsbIn)))

namespace flopoco{

    class FixIIRShiftAdd : public Operator {

    private:
        struct adderSigns /**< adder signs for the use with IntConstMultShiftAdd component */
        {
            string signLeft = "+";
            string signRight = "+";
        };
        void computeImpulseResponse(); // evaluates the filter on an impulsion
        /** @brief checks if an integer number is a power of two.
         * checkNegativeNumber indicates whether ONLY positive or neg. numbers are checked*/
        bool isPowerOfTwo(int64_t, bool);
    public:
        /** @brief Constructor ; you must use bitheap in case of negative coefficient*/
        FixIIRShiftAdd(OperatorPtr parentOp, Target* target, int lsbIn, int lsbOut, int msbIn, int guardBits, vector<string> coeffb, vector<string> coeffa, int64_t shifta, int64_t shiftb, string method, string grapha, string graphb, double H, double Heps);

        /** @brief Destructor */
        ~FixIIRShiftAdd();

        // Below all the functions needed to test the operator
        /**
         * @brief the emulate function is used to simulate in software the operator
         * in order to compare this result with those outputed by the vhdl opertator
         */
        void emulate(TestCase * tc);

        /** @brief function used to create Standard testCase defined by the developper */
        void buildStandardTestCases(TestCaseList* tcl);

        /** @brief returns the index for a given value in a given array, -1 if not found */
        int64_t getIndexCoeff(int64_t*, int64_t, int64_t);

        void preventDoubleNegativeInputsAdder(uint16_t, vector<adderSigns>*);

        // User-interface stuff
        /** Factory method */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
        static TestList unitTest(int index);
        static void registerFactory();

    private:
        int lsbIn;					/**< weight of the LSB in the input, considered as a signed number in (-1,1) */
        int msbOutIIR;					/**< weight of the MSB in the result */
        int msbOutFIR;
        int lsbOut;					/**< weight of the LSB in the result */
        int msbIn;
        int lsbInt;
        int msbInt;
        int wOut;
        int64_t bitGainForAdditons;
        uint64_t wIn;                     /**<  input word size */
        int guardBits;                  /**< number of guard bits for feedback path */
        vector<string> coeffb;			/**< the b_i coefficients as strings */
        vector<string> coeffa;			/**< the a_i coefficients as strings */
        int64_t shifta;             /**< number of shifts for coeffa (not scaling factor) */
        int64_t shiftb;             /**< number of shifts for coeffb (not scaling factor) */
        string method;
        string grapha;
        string graphb;
        bool buildWorstCaseTestBench; /**< if true, build the worst-case test bench */
        uint32_t n;							/**< number of taps on the numerator */
        uint32_t m;							/**< number of taps on the denominator */
        double H;                /**< worst case peak gain of filter */
        double Heps;
        int wcpgBitGain;            /**< required amount of bits for wcpg */
        mpfr_t* coeffb_mp_f;			/**< the coefficients as MPFR numbers */
        mpfr_t* coeffa_mp_f;			/**< the coefficients as MPFR numbers */
        double* coeffb_d;           /**< version of coeffb as C-style arrays of double, because WCPG needs it this way */
        double* coeffa_d;           /**< version of coeffa as C-style arrays of double, because WCPG needs it this way */
        int64_t* coeffa_si;         /**< coefficients as signed int */
        int64_t* coeffb_si;         /**< coefficients as signed int */
        double* coeffa_d_wcpg;      /**< coefficients as double for wcpg */
        double* coeffb_d_wcpg;      /**< coefficients as double for wcpg */
        int* msbRegisterForwardOut; /**< output msb of registers in forward path  */
        int* msbAdditionsForwardOut; /**< output msb of additions in forward path  */
        int* msbMultiplicationsForwardOut; /**< output msb of multiplications in forward path  */
        int msbScaleBOut;             /**< output msb of scaleB  */
        int lsbScaleBOut;             /**< output lsb of scaleB  */
        int lsbForwardPath;             /**< lsb in forward path (constant im forward path)  */
        int* msbRegisterFeedbackOut;  /**< output msb of registers in feedback path  */
        int* lsbRegisterFeedbackOut;  /**< output lsb of registers in feedback path  */
        int* msbMultiplicationsFeedbackOut; /**< output msb of multiplications in feedback path  */
        int* lsbMultiplicationsFeedbackOut; /**< output lsb of multiplications in feedback path  */
        int* msbAdditionsFeedbackOut; /**< output msb of additions in feedback path  */
        int* lsbAdditionsFeedbackOut; /**< output lsb of additions in feedback path  */
        int msbScaleAOut;               /**< output msb of scaleA  */
        int lsbScaleAOut;               /**< output msb of scaleA  */
        int msbTruncationIntOut;           /**< msb after the truncation to wInt*/
        int lsbTruncationIntOut;           /**< lsb after the truncation to wInt */
        int wordsizefeedbackAdd0ResizedReg0; /**< word size of the aligned word of feedbackReg0 for input for feedbackAdd0 */
        int wordsizefeedbackAdd0ResizedScaleB; /**< word size of the aligned word of scaleB for input for feedbackAdd0 */
        mpz_t* registerForwardTB_mpz;
        mpz_t* registerFeedbackTB_mpz;
        mpz_t resultForwardPathTB_mpz; /**< result of the forward path in testbench */
        mpz_t feedbackAdd0;
        bool isFIR = false; /**< if not all parameters for feedback path given FIR is assumed instead of IIR */

        /* faithfully rounded testbench */
        vector<double> yi;  // outputs in the trace of simulation in double precision
        uint64_t vanishingK; /**< faithul testbench */
        double miny, maxy; /**< faithul testbench */
        vector<double> ui;  // inputs in the trace of simulation in double precision
        uint64_t currentIndex;       // used for round-robin access to the history
        mpfr_t* xHistory; // history of x used by emulate
        mpfr_t* yHistory; // history of y (result) used by emulate
        int hugePrec;
        mpfr_t* coeffb_mp_f_scaled;			/**< the scaled coefficients as MPFR numbers */
        mpfr_t* coeffa_mp_f_scaled;			/**< the scaled coefficients as MPFR numbers */
    };

}