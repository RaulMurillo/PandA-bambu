/*
  A Fixed-Point divider for FloPoCo

  Author:  Raul Murillo

  This file is part of the FloPoCo project

  Initial software.
  Copyright © UCM,
  2022.
  All rights reserved.

*/

#include <iostream>
#include <sstream>

#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "FixDiv.hpp"
// #include "Posit/Decoder/PositFastDecoder.hpp"
// #include "Posit/Encoder/PositFastEncoder.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
#include "IntMult/IntMultiplier.hpp"
// #include "TestBenches/PositNumber.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

    FixDiv::FixDiv(OperatorPtr parentOp, Target *target, int ints, int frac, float dspOccupationThreshold, bool truncate, int iters, bool useGoldschmidt, int LUT_out) : Operator(parentOp, target), ints_(ints), frac_(frac), dspOccupationThreshold(dspOccupationThreshold), truncate_(truncate), iters_(iters), goldschmidt(useGoldschmidt), LUT_out_(LUT_out)
	{
		setCopyrightString("Raul Murillo (2022)");
		ostringstream name;
		srcFileName = "FixDiv";
		// Signed multiplication needs this library
		// useNumericStd();

        width_ = 1+ints_+frac_; // Add sign bit

		if (width_ < 3)
		{
			throw std::string("PositDiv Constructor: width is too small, should be greater than two");
		}

		// -------- Parameter set up -----------------

		// int sizeRegime = intlog2(width_ - 1) + 1;
		// wE_ = sizeRegime + wES_;
		// wF_ = width_ - 3 - wES_;
		// int maxExp = (1 << wES_) * (width_ - 2);
		// int minExp = -maxExp;
        mpfr_t t;
        mpfr_init(t);
        mpfr_set_d (t, 2.0, MPFR_RNDD);
        // const string two = signedFixPointNumber(t, ints_+frac_, 0);   // Doesn't work as expected!!
        const string two = zg(ints_-1, -1)+"10"+zg(frac_, 1);   // 2.0 in fixed point
        const string one = zg(ints_-1, -1)+"01"+zg(frac_, 1);   // 1.0 in fixed point

		name << "FixDiv_" << ints_ << "_" << frac_;
        if (truncate_)
            name << "_trunc";
        if (iters_ > 0){
            if (goldschmidt)
                name << "_Goldschmidt";
            else
                name << "_Newton_Raphson";
        }
		setNameWithFreqAndUID(name.str());

		addInput("X", width_);
		addInput("Y", width_);
		if (iters <= 0)
            addOutput("R", width_+4-1);
		else if (truncate_ )
            addOutput("R", width_);
        else
            addOutput("R", 2*width_);

		// if (computeSticky_)
        //     addOutput("Sticky"); /* if we require a sticky bit computation */

		addFullComment("Start of vhdl generation");

		REPORT(INFO, "Declaration of FixDiv \n");
		REPORT(DETAILED, "this operator has received the following parameters: " << ints << ", " << frac);
		REPORT(DEBUG, "debug of FixDiv");

        if (iters_ > 0){
            // const int LUT_out = 10;
            if (goldschmidt){   // 2n bits precision with rectangular multipliers
                //====================================================================|
                addFullComment("Goldschmidt's division algorithm");
                //====================================================================|
                /*
                / q_0 = X; d_0 = Y; q_k+1 = quotient
                / d_{i+1} = d_{i}*p_{i}
                / q_{i+1} = q_{i}*p_{i}
                / p_{i+1} = 2 - d_{i+1}
                */
                // We can use rectangular multipliers, i.e. multiply signals of different size
                // Thus, we remove trailing zeros from initial guess p_0

                // const string guess_pos = zg(ints_+1, -1)+"11"+zg(frac_-2, 1);   // 0.75 in fixed point
                // const string guess_neg = og(ints_+1, -1)+"01"+zg(frac_-2, 1);   // -0.75 in fixed point
                const string guess_pos = zg(ints_+1, -1)+"11"+zg(0, 1);   // 0.75 in fixed point
                const string guess_neg = og(ints_+1, -1)+"01"+zg(0, 1);   // -0.75 in fixed point

                // p_0 = 1/d
                const int EK = LUT_out_-4; // Total = EK+4
                // Just need to check 2 bits (excluding the leftmost)
                vhdl << tab << "with Y" << range(width_-1, width_-(4+EK)) << " select " << declare(getTarget()->logicDelay(2), "p_0", 4+EK) << " <= " << endl
                    << tab << tab << og(ints_+1, -1)<<"1"<<zg(1, 1) << " when \"00\"," << endl // -0.5 when Y=-2.0
                    << tab << tab << guess_neg << " when \"01\"," << endl // -0.75 when Y=-1.0-0.x
                    // << tab << tab << og(ints_, -1)<<"1"<<zg(frac_, 1) << " when \"110\"," << endl // -1.0 when Y=-1.0    // Not possible
                    << tab << tab << zg(ints_, -1)<<"1"<<zg(2, 1) << " when \"10\"," << endl // 1.0 when Y=1.0+0.x
                    << tab << tab << guess_pos << " when \"11\"," << endl // 0.75 when Y=1.5+0.x
                    << tab << tab << "\"" << string(4+EK, '-') << "\" when others;" << endl;

                string p_i, d_i, q_i, p_i1, d_i1, q_i1, d_i1_tmp, q_i1_tmp, round_p_i;
                p_i = "p_0";
                d_i = "d";
                q_i = "q";
                d_i1 = "d_0";
                d_i1_tmp = "tmp_d_0";
                q_i1 = "q_0";
                q_i1_tmp = "tmp_q_0";
                // d_0 = Y
                vhdl << tab << declare("d", width_) << " <= Y;" << endl;
                // q_0 = X
                vhdl << tab << declare("q", width_) << " <= X;" << endl;

                int w_p, w_d, w_p_next, w_d_next, rem_bits;
                w_p = 2+EK;	// fraction bits of p_i
                w_d = frac_;// fraction bits of q_i and d_i

    /*
                for (int i = 0; i<iters_; ++i){
                    addComment("Iteration "+std::to_string(i));
                    w_d_next = min(w_d+w_p, frac_);   // Limit maximum size of d_{i}, q_{i}, p_{i}
                    w_p_next = w_d_next;//w_p + 8;
                    p_i = join("p_", i);
                    p_i1 = join("p_", i+1);
                    round_p_i = join("round_p_", i+1);

                    // q_{i+1} = q_{i}*p_{i}
                    newInstance("IntMultiplier", join("mult_q_", i), "wX="+to_string(1+ints_+w_d)+" wY="+to_string(1+ints_+w_p)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i+",Y=>"+p_i, "R=>"+q_i1_tmp);
                    // d_{i+1} = d_{i}*p_{i}
                    newInstance("IntMultiplier", join("mult_d_", i), "wX="+to_string(1+ints_+w_d)+" wY="+to_string(1+ints_+w_p)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),"X=>"+d_i+",Y=>"+p_i, "R=>"+d_i1_tmp);

                    rem_bits = w_d+w_p-w_d_next;
                    if (rem_bits>0){    // Round d_{i+1} and q_{i+1}
                        vhdl << tab << declare(getTarget()->adderDelay(1+ints_+w_p_next), d_i1, 1+ints_+w_d_next) << " <= " << d_i1_tmp << range(ints_+w_d_next+rem_bits, rem_bits) << " + " << d_i1_tmp << of(rem_bits-1) << ";" << endl;
                        vhdl << tab << declare(getTarget()->adderDelay(1+ints_+w_p_next), q_i1, 1+ints_+w_d_next) << " <= " << q_i1_tmp << range(ints_+w_d_next+rem_bits, rem_bits) << " + " << q_i1_tmp << of(rem_bits-1) << ";" << endl;
                        
                        // p_{i+1} = 2 - d_{i+1}
                        // The computation of p_{i+1} is very simple (by two’s complementing d_{i+1}, keeping the sign).
                        int jkl = w_d_next+rem_bits-w_p_next;
                        vhdl << tab << declare(getTarget()->eqConstComparatorDelay(rem_bits), round_p_i, 1, false) << " <= '0' when " << d_i1_tmp << range(jkl-1, 0) << " = 0 else '1';" << endl;
                        vhdl << tab << declare(getTarget()->adderDelay(ints_+w_p_next), p_i1, 1+ints_+w_p_next) << " <= '0' & NOT(" << d_i1_tmp << range(ints_+w_d_next+rem_bits-1, jkl) << ")+NOT(" << round_p_i << ");" << endl;

                    }
                    else{   // Now unreachable
                        vhdl << tab << declare(d_i1, 1+ints_+w_d_next) << " <= " << d_i1_tmp << range(ints_+w_d_next, 0) <<";" << endl;
                        vhdl << tab << declare(q_i1, 1+ints_+w_d_next) << " <= " << q_i1_tmp << range(ints_+w_d_next, 0) <<";" << endl;

                        // p_{i+1} = 2 - d_{i+1}
                        // The computation of p_{i+1} is very simple (by two’s complementing d_{i+1}, keeping the sign).
                        vhdl << tab << declare(getTarget()->adderDelay(ints_+w_p_next), p_i1, 1+ints_+w_p_next) << " <= '0' & NOT(" << d_i1 <<range(ints_+w_p_next-1,0)<< ")+1;" << endl;
                    }

                    w_d = w_d_next;
                    w_p = w_p_next;

                    d_i = join("d_", i);
                    q_i = join("q_", i);
                    d_i1 = join("d_", i+1);
                    q_i1 = join("q_", i+1);
                    d_i1_tmp = join("tmp_d_", i+1);
                    q_i1_tmp = join("tmp_q_", i+1);

                }
                p_i = join("p_", iters_);
                addComment("Iteration "+std::to_string(iters_));
                // q_{i+1} = q_{i}*p_{i}
                newInstance("IntMultiplier", join("mult_q_", iters_), "wX="+to_string(1+ints_+w_d)+" wY="+to_string(1+ints_+w_p)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i+",Y=>"+p_i, "R=>"+q_i1_tmp);

    */
                

                    // addComment("Iteration "+std::to_string(i+1));
                    // // r_{i+1} = r_i*(2-Y*r_i)
                    // ri = join("r_", i);
                    // mi = join("m_", i);
                    // neg_mi = join("neg_m_", i);
                    // si = join("s_", i);
                    // tmp_next_ri = join("tmp_r_", i+1);
                    // next_ri = join("r_", i+1);
                    // round_i = join("round_", i);
                    // // Y*r_i
                    // newInstance("IntMultiplier", join("mult_0_", i), "wX="+to_string(width_)+" wY="+to_string(r_size)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>"+ri, "R=>"+mi);
                    // // (2-Y*r_i)
                    // // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_-1), si, width_) << " <= " << mi << of(2*frac_+ints_+frac_) << " & (NOT(" << mi << range(2*frac_+ints_-1+frac_, frac_) << ") + 1) ;" << endl;
                    // // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->eqConstComparatorDelay(width_+2), round_i, 1, false) << " <= '0' when " << mi << range(width_+1, 0) << "=0 else '1';" << endl;
                    // // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_), si, width_) << " <= ('1' & " << zg(width_-1) << ") - ('0' & " << mi << range(2*width_-1, width_+2) << " & " << round_i << " ) ;" << endl;
                    // // vhdl << tab << declare(getTarget()->logicDelay(), si, width_) << " <= '0' & NOT(" << mi << range(2*width_-4, width_-2) << ");" << endl;

                    // vhdl << tab << declare(getTarget()->logicDelay(), si, width_) << " <= '0' & NOT(" << mi << range(width_+r_size-4, r_size-2) << ");" << endl;

                    // // r_{i+1} = r_i*(2-Y*r_i)
                    // newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(width_)+" wY="+to_string(r_size)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+si+",Y=>"+ri, "R=>"+tmp_next_ri);
                    // // Rounding r_{i+1}

                    // int next_r_s = i==iters_-1 ? width_ : min(r_size+r_0_size-2, width_);
                    
                    // vhdl << tab << declare(next_ri, next_r_s) << " <= " << tmp_next_ri << range(width_+r_size-3, width_+r_size-2-next_r_s) << " + " << tmp_next_ri << of(width_+r_size-2-next_r_s-1) << ";" << endl;
                    // r_size = next_r_s;

                const int w_p_0 = w_p;

                for (int i = 0; i<iters_; ++i){
                    addComment("Iteration "+std::to_string(i));
                    p_i = join("p_", i);
                    p_i1 = join("p_", i+1);

                    int mult_p_size = 2*(1+ints_)+w_d+w_p;
                    // q_{i+1} = q_{i}*p_{i}
                    newInstance("IntMultiplier", join("mult_q_", i), "wX="+to_string(1+ints_+w_d)+" wY="+to_string(1+ints_+w_p)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i+",Y=>"+p_i, "R=>"+q_i1_tmp);
                    vhdl << tab << declare(getTarget()->adderDelay(1+ints_+w_d), q_i1, 1+ints_+w_d) << " <= " << q_i1_tmp << range(mult_p_size-3,mult_p_size-3-ints_-w_d) << " + " << q_i1_tmp << of(mult_p_size-3-ints_-w_d-1) << ";" << endl;
                    // d_{i+1} = d_{i}*p_{i}
                    newInstance("IntMultiplier", join("mult_d_", i), "wX="+to_string(1+ints_+w_d)+" wY="+to_string(1+ints_+w_p)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),"X=>"+d_i+",Y=>"+p_i, "R=>"+d_i1_tmp);
                    vhdl << tab << declare(getTarget()->adderDelay(1+ints_+w_d), d_i1, 1+ints_+w_d) << " <= " << d_i1_tmp << range(mult_p_size-3,mult_p_size-3-ints_-w_d) << " + " << d_i1_tmp << of(mult_p_size-3-ints_-w_d-1) << ";" << endl;

                    // p_{i+1} = 2 - d_{i+1}
                    // The computation of p_{i+1} is very simple (by two’s complementing d_{i+1}, keeping the sign).
                    int w_p_next = i==iters_-1 ? width_-(1+ints_) : w_p+w_p_0;
                    vhdl << tab << declare(getTarget()->logicDelay(), p_i1, 1+ints_+w_p_next) << " <= '0' & NOT(" << d_i1_tmp << range(mult_p_size-4, mult_p_size-3-ints-w_p_next) << ");" << endl;
                    w_p = w_p_next; // update

                    // Adjust names
                    d_i = join("d_", i);
                    q_i = join("q_", i);
                    d_i1 = join("d_", i+1);
                    q_i1 = join("q_", i+1);
                    d_i1_tmp = join("tmp_d_", i+1);
                    q_i1_tmp = join("tmp_q_", i+1);
                }
                addComment("Iteration "+std::to_string(iters_));
                // q_{i+1} = q_{i}*p_{i}
                newInstance("IntMultiplier",
                            join("mult_q_", iters_),
                            "wX="+to_string(1+ints_+w_d)+" wY="+to_string(1+ints_+w_p)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),
                            "X=>"+q_i+",Y=>"+p_i1,
                            "R=>"+q_i1_tmp
                            );

                
                /*
            addComment("Iteration "+std::to_string(iters_)+" - generate final quotient");
                // quotient = q_{k+1} = q_{k}*p_{k}
                newInstance("IntMultiplier", join("mult_q_", iters_), "wX="+to_string(width_+frac_)+" wY="+to_string(width_+frac_)+ " signedIO=true dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i1+",Y=>"+p_i1, "R=>quotient_long");
            vhdl << tab << declare("quotient", 2*width_) << " <= quotient_long" << range(2*(width_+frac_)-1, 2*frac_) << " + quotient_long" << of(2*frac_-1) << ";" << endl;
                */
            // Rounding final quotient, if necessary
            int w_q = (1+ints_+w_d) + (1+ints_+w_p);
            if (w_q <= 2*width_) {
                int extrabits = 2*width_ - w_q;
                vhdl << tab << declare("quotient", 2*width_) << " <= " << q_i1_tmp << " & " << zg(extrabits) << ";" << endl;
            }
            else {
                int extrabits = w_q - 2*width_;
                vhdl << tab << declare("quotient", 2*width_) << " <= " << q_i1_tmp << range(w_q-1, extrabits) << " + " << q_i1_tmp << of(extrabits-1) << ";" << endl;
            }
            }
            else if (goldschmidt && false){   // 2n bits precision
                //====================================================================|
                addFullComment("Goldschmidt's division algorithm");
                //====================================================================|
                /*
                / q_0 = X; d_0 = Y; q_k = quotient
                / d_{i+1} = d_i*p_{i}
                / q_{i+1} = q_i*p_{i}
                / p_{i+1} = 2 - d_{i+1}
                */
                // vhdl << tab << declare("Y_sign", 1, false) << " <= Y" << of(width_-1) << ";" << endl;
                // vhdl << tab << "with Y_sign select " << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_), "d_0", width_) << " <= " << endl
                //      << tab << tab << "'0' & Y" << range(width_-1, 1) << " when '0'," << endl
                //      << tab << tab << "'0' & (NOT(Y" << range(width_-1, 1) << ")+1) when others;" << endl;
                // vhdl << tab << "with Y_sign select " << declare(getTarget()->adderDelay(width_), "q_0", width_) << " <= " << endl
                //      << tab << tab << "X" << of(width_-1) << " & X" << range(width_-1, 1) << " when '0'," << endl
                //      << tab << tab << "X" << of(width_-1) << " & (NOT(X" << range(width_-1, 1) << ")+1) when others;" << endl;

                // vhdl << tab << declare("d_0", width_) << " <= Y" << of(width_-1) << " & Y" << range(width_-1, 1) << ";" << endl;
                // vhdl << tab << declare("q_0", width_) << " <= X" << of(width_-1) << " & X" << range(width_-1, 1) << ";" << endl;

                const string guess_pos = zg(ints_+1, -1)+"11"+zg(frac_-2+frac_, 1);   // 0.75 in fixed point
                const string guess_neg = og(ints_+1, -1)+"01"+zg(frac_-2+frac_, 1);   // -0.75 in fixed point
                // vhdl << tab << declare("two_fix", width_) << " <= " << two <<" ;" << endl;   // 2.0 in fixed point
                // vhdl << tab << declare("one_bit") << " <= '1' ;" << endl;

                // p_0 = 1/d
                vhdl << tab << "with Y" << range(width_-1, width_-3) << " select " << declare(getTarget()->logicDelay(3), "p_0", width_+frac_) << " <= " << endl
                    << tab << tab << og(ints_+1, -1)<<"1"<<zg(frac_-1+frac_, 1) << " when \"100\"," << endl // -0.5 when Y=-2.0
                    << tab << tab << guess_neg << " when \"101\"," << endl // -0.75 when Y=-1.0-0.x
                    // << tab << tab << og(ints_, -1)<<"1"<<zg(frac_, 1) << " when \"110\"," << endl // -1.0 when Y=-1.0    // Not possible
                    << tab << tab << zg(ints_, -1)<<"1"<<zg(frac_+frac_, 1) << " when \"010\"," << endl // 1.0 when Y=1.0+0.x
                    << tab << tab << guess_pos << " when \"011\"," << endl // 0.75 when Y=1.5+0.x
                    << tab << tab << "\"" << string(width_+frac_, '-') << "\" when others;" << endl;
                // newInstance("IntMultiplier", "mult_q_0", "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>X,Y=>p_0", "R=>q_0_tmp");
                // newInstance("IntMultiplier", "mult_d_0", "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>p_0", "R=>d_0_tmp");
                // vhdl << tab << declare("q_0", width_) << " <= q_0_tmp"<<range(2*frac_+ints_, frac_)<<";" << endl;
                // d_0 = Y
                vhdl << tab << declare("d_0", width_+frac_) << " <= Y &" << zg(frac_) << ";" << endl;
                // vhdl << tab << declare("d_0", width_) << " <= d_0_tmp"<<range(2*frac_+ints_, frac_)<<";" << endl;
                // q_0 = X
                vhdl << tab << declare("q_0", width_+frac_) << " <= X &" << zg(frac_) << ";" << endl;

                string p_i, d_i, q_i, p_i1, d_i1, q_i1, d_i1_tmp, q_i1_tmp, neg_d_i1;
                for (int i = 0; i<iters_; ++i){
                    p_i = join("p_", i);
                    d_i = join("d_", i);
                    q_i = join("q_", i);
                    p_i1 = join("p_", i+1);
                    d_i1 = join("d_", i+1);
                    neg_d_i1 = join("neg_d_", i+1);
                    q_i1 = join("q_", i+1);
                    d_i1_tmp = join("tmp_d_", i+1);
                    q_i1_tmp = join("tmp_q_", i+1);
                    // d_{i+1} = d_i*p_{i}
                    newInstance("IntMultiplier", join("mult_d_", i+1), "wX="+to_string(width_+frac_)+" wY="+to_string(width_+frac_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+d_i+",Y=>"+p_i, "R=>"+d_i1_tmp);
                    // q_{i+1} = q_i*p_{i}
                    newInstance("IntMultiplier", join("mult_q_", i+1), "wX="+to_string(width_+frac_)+" wY="+to_string(width_+frac_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i+",Y=>"+p_i, "R=>"+q_i1_tmp);
                    vhdl << tab << declare(d_i1, width_+frac_) << " <= " << d_i1_tmp << range(4*frac_+ints_, 2*frac_) <<";" << endl;
                    vhdl << tab << declare(q_i1, width_+frac_) << " <= " << q_i1_tmp << range(4*frac_+ints_, 2*frac_) <<";" << endl;
                    // the computation of p_{i+1} is very simple (by two’s complementing d_{i+1}).
                    // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_-1), p_i1, width_) << " <= '0' & (NOT("<<d_i<<range(width_-2,0)<<")+1);" << endl;
                    vhdl << tab << declare(getTarget()->adderDelay(width_-1), p_i1, width_+frac_) << " <= '0' & NOT(" << d_i1 <<range(width_-2+frac_,0)<< ")+1;" << endl;
                    // newInstance("IntAdder", join("sub_", i), "wIn="+to_string(width_),"X=>two_fix,Cin=>one_bit,Y=>" + neg_d_i1, "R=>"+p_i1);
                    // if(i < iters_-1){
                    //     vhdl << tab << declare(d_i1, width_) << " <= " << d_i1_tmp << range(2*frac_+ints_, frac_) <<";" << endl;
                    //     vhdl << tab << declare(q_i1, width_) << " <= " << q_i1_tmp << range(2*frac_+ints_, frac_) <<";" << endl;
                    // }
                }
                // vhdl << tab << declare("quotient", 2*width_) << " <= " << q_i1_tmp << ";" << endl;
                newInstance("IntMultiplier", join("mult_q_", iters_+1), "wX="+to_string(width_+frac_)+" wY="+to_string(width_+frac_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i1+",Y=>"+p_i1, "R=>quotient_long");
                vhdl << tab << declare("quotient", 2*width_) << " <= quotient_long" << range(2*(width_+frac_)-1, 2*frac_) << ";" << endl;


                // vhdl << tab << declare("sign_ext", ints_+2) << " <= (others=>NOT(X"<<of(width_-1)<<"));" << endl;
                // vhdl << tab << "with Y" << range(width_-1, width_-3) << " select " << declare(getTarget()->adderDelay(width_), "quotient", 2*width_) << " <= " << endl
                //     << tab << tab << "sign_ext & (NOT(X)+1) & " << zg(frac_-1)<< " when \"100\"," << endl // -2.0
                //     << tab << tab << n_i1_tmp << " when others;" << endl;
            }
            else if (goldschmidt && false){    // n bits precision - not enough
                //====================================================================|
                addFullComment("Goldschmidt's division algorithm");
                //====================================================================|
                /*
                / q_0 = X; d_0 = Y; q_k = quotient
                / d_{i+1} = d_i*p_{i}
                / q_{i+1} = q_i*p_{i}
                / p_{i+1} = 2 - d_{i+1}
                */
                // vhdl << tab << declare("Y_sign", 1, false) << " <= Y" << of(width_-1) << ";" << endl;
                // vhdl << tab << "with Y_sign select " << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_), "d_0", width_) << " <= " << endl
                //      << tab << tab << "'0' & Y" << range(width_-1, 1) << " when '0'," << endl
                //      << tab << tab << "'0' & (NOT(Y" << range(width_-1, 1) << ")+1) when others;" << endl;
                // vhdl << tab << "with Y_sign select " << declare(getTarget()->adderDelay(width_), "q_0", width_) << " <= " << endl
                //      << tab << tab << "X" << of(width_-1) << " & X" << range(width_-1, 1) << " when '0'," << endl
                //      << tab << tab << "X" << of(width_-1) << " & (NOT(X" << range(width_-1, 1) << ")+1) when others;" << endl;

                // vhdl << tab << declare("d_0", width_) << " <= Y" << of(width_-1) << " & Y" << range(width_-1, 1) << ";" << endl;
                // vhdl << tab << declare("q_0", width_) << " <= X" << of(width_-1) << " & X" << range(width_-1, 1) << ";" << endl;

                const string guess_pos = zg(ints_+1, -1)+"11"+zg(frac_-2, 1);   // 0.75 in fixed point
                const string guess_neg = og(ints_+1, -1)+"01"+zg(frac_-2, 1);   // -0.75 in fixed point
                vhdl << tab << declare("two_fix", width_) << " <= " << two <<" ;" << endl;   // 2.0 in fixed point
                vhdl << tab << declare("one_bit") << " <= '1' ;" << endl;

                // p_0 = 1/d
                vhdl << tab << "with Y" << range(width_-1, width_-3) << " select " << declare(getTarget()->logicDelay(3), "p_0", width_) << " <= " << endl
                    << tab << tab << og(ints_+1, -1)<<"1"<<zg(frac_-1, 1) << " when \"100\"," << endl // -0.5 when Y=-2.0
                    << tab << tab << guess_neg << " when \"101\"," << endl // -0.75 when Y=-1.0-0.x
                    // << tab << tab << og(ints_, -1)<<"1"<<zg(frac_, 1) << " when \"110\"," << endl // -1.0 when Y=-1.0    // Not possible
                    << tab << tab << zg(ints_, -1)<<"1"<<zg(frac_, 1) << " when \"010\"," << endl // 1.0 when Y=1.0+0.x
                    << tab << tab << guess_pos << " when \"011\"," << endl // 0.75 when Y=1.5+0.x
                    << tab << tab << "\"" << string(width_, '-') << "\" when others;" << endl;
                // newInstance("IntMultiplier", "mult_q_0", "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>X,Y=>p_0", "R=>q_0_tmp");
                // newInstance("IntMultiplier", "mult_d_0", "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>p_0", "R=>d_0_tmp");
                // vhdl << tab << declare("q_0", width_) << " <= q_0_tmp"<<range(2*frac_+ints_, frac_)<<";" << endl;
                // d_0 = Y
                vhdl << tab << declare("d_0", width_) << " <= Y;" << endl;
                // vhdl << tab << declare("d_0", width_) << " <= d_0_tmp"<<range(2*frac_+ints_, frac_)<<";" << endl;
                // q_0 = X
                vhdl << tab << declare("q_0", width_) << " <= X;" << endl;

                string p_i, d_i, q_i, p_i1, d_i1, q_i1, d_i1_tmp, q_i1_tmp, neg_d_i1;
                for (int i = 0; i<iters_; ++i){
                    p_i = join("p_", i);
                    d_i = join("d_", i);
                    q_i = join("q_", i);
                    p_i1 = join("p_", i+1);
                    d_i1 = join("d_", i+1);
                    neg_d_i1 = join("neg_d_", i+1);
                    q_i1 = join("q_", i+1);
                    d_i1_tmp = join("tmp_d_", i+1);
                    q_i1_tmp = join("tmp_q_", i+1);
                    // d_{i+1} = d_i*p_{i}
                    newInstance("IntMultiplier", join("mult_d_", i+1), "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+d_i+",Y=>"+p_i, "R=>"+d_i1_tmp);
                    // q_{i+1} = q_i*p_{i}
                    newInstance("IntMultiplier", join("mult_q_", i+1), "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i+",Y=>"+p_i, "R=>"+q_i1_tmp);
                    vhdl << tab << declare(d_i1, width_) << " <= " << d_i1_tmp << range(2*frac_+ints_, frac_) <<";" << endl;
                    vhdl << tab << declare(q_i1, width_) << " <= " << q_i1_tmp << range(2*frac_+ints_, frac_) <<";" << endl;
                    // the computation of p_{i+1} is very simple (by two’s complementing d_{i+1}).
                    // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_-1), p_i1, width_) << " <= '0' & (NOT("<<d_i<<range(width_-2,0)<<")+1);" << endl;
                    vhdl << tab << declare(getTarget()->adderDelay(width_-1), p_i1, width_) << " <= '0' & NOT(" << d_i1 <<range(width_-2,0)<< ")+1;" << endl;
                    // newInstance("IntAdder", join("sub_", i), "wIn="+to_string(width_),"X=>two_fix,Cin=>one_bit,Y=>" + neg_d_i1, "R=>"+p_i1);
                    // if(i < iters_-1){
                    //     vhdl << tab << declare(d_i1, width_) << " <= " << d_i1_tmp << range(2*frac_+ints_, frac_) <<";" << endl;
                    //     vhdl << tab << declare(q_i1, width_) << " <= " << q_i1_tmp << range(2*frac_+ints_, frac_) <<";" << endl;
                    // }
                }
                // vhdl << tab << declare("quotient", 2*width_) << " <= " << q_i1_tmp << ";" << endl;
                newInstance("IntMultiplier", join("mult_q_", iters_+1), "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+q_i1+",Y=>"+p_i1, "R=>quotient");


                // vhdl << tab << declare("sign_ext", ints_+2) << " <= (others=>NOT(X"<<of(width_-1)<<"));" << endl;
                // vhdl << tab << "with Y" << range(width_-1, width_-3) << " select " << declare(getTarget()->adderDelay(width_), "quotient", 2*width_) << " <= " << endl
                //     << tab << tab << "sign_ext & (NOT(X)+1) & " << zg(frac_-1)<< " when \"100\"," << endl // -2.0
                //     << tab << tab << n_i1_tmp << " when others;" << endl;
            }
            else{
            //====================================================================|
            addFullComment("Newton-Raphson algorithm");
            //====================================================================|
            addComment("Get initial guess for reciprocal");
            // Constant value

            // We can use rectangular multipliers, i.e. multiply signals of different size
            // Thus, we remove trailing zeros from initial guess r_0

            // /* Errors: 26095 */
            // const string guess_pos = zg(ints_+1, -1)+"11\"";//+zg(frac_-2, 1);   // 0.75 in fixed point
            // const string guess_neg = og(ints_+1, -1)+"01\"";//+zg(frac_-2, 1);   // -0.75 in fixed point
            // vhdl << tab << "with Y" << of(width_-1) << " select " << declare(getTarget()->logicDelay(1), "r_0", width_-frac_+2) << " <= " << endl
            // 	 << tab << tab << guess_neg << " when '1'," << endl
            // 	 << tab << tab << guess_pos << " when '0'," << endl
            // 	 << tab << tab << "\"" << string(width_-frac_+2, '-') << "\" when others;" << endl;

            /* Errors: 24075 */
            /**
            
            int r_0_size = width_-frac_+2;
            vhdl << tab << "with Y" << range(width_-2, width_-3) << " select " << declare(getTarget()->logicDelay(2), "r_0", width_-frac_+2) << " <= " << endl
                << tab << tab << og(ints_+1, -1)<<"10"<<zg(0, 1) << " when \"00\"," << endl // -0.5 when Y = -2.0 + 0.x
                << tab << tab << og(ints_+1, -1)<<"01"<<zg(0, 1) << " when \"01\"," << endl // -0.75 when Y = -1.x
                // << tab << tab << og(ints_, -1)<<"1"<<zg(frac_, 1) << " when \"110\"," << endl // -1.0 when Y = -1.0    // Not possible
                << tab << tab << zg(ints_, -1)<<"10"<<zg(1, 1) << " when \"10\"," << endl // 1.0 when Y = 1.0+0.x
                << tab << tab << zg(ints_+1, -1)<<"11"<<zg(0, 1) << " when \"11\"," << endl // 0.75 when Y = 1.5+0.x
                << tab << tab << "\"" << string(width_-frac_+2, '-') << "\" when others;" << endl;

            */

            // /* Errors: 27758 */
            // vhdl << tab << "with Y" << range(width_-1, width_-3) << " select " << declare(getTarget()->logicDelay(3), "r_0", width_-frac_+2) << " <= " << endl
            //     << tab << tab << og(ints_+1, -1)<<"10\"" << " when \"100\"," << endl // -2.0
            //     << tab << tab << "(Y"<<of(width_-1)<<" & (NOT(Y"<<range(width_-2, frac_-2)<<") + 1)) when others;" << endl;
            //     // << tab << tab << og(ints_+1, -1)<<"1"<<zg(frac_-1, 1) << " when \"100\"," << endl // -2.0
            //     // << tab << tab << "(Y"<<of(width_-1)<<" & (NOT(Y"<<range(width_-2,0)<<") + 1)) when others;" << endl;
            
            
            /* ANOTHER APPROACH FOR LUTs*/
            /**/
            // Just need to check 2 bits (excluding the leftmost)
            int r_0_size = LUT_out_;//width_-frac_+2;
            // vhdl << tab << "with Y" << range(width_-2, width_-3) << " select " << declare(getTarget()->logicDelay(2), "r_0", r_0_size) << " <= " << endl
            // vhdl << tab << "with Y" << range(width_-1, width_-r_0_size) << " select " << declare(getTarget()->logicDelay(2), "r_0", r_0_size) << " <= " << endl
            vhdl << tab << "with Y" << range(width_-1, width_-r_0_size) << " select " << declare(getTarget()->logicDelay(2), "r_0", r_0_size) << " <= " << endl
                << tab << tab << og(ints_+1, -1)<<"10"<<zg(0, 1) << " when \"00\"," << endl // -0.5 when Y = -2.0 + 0.x
                << tab << tab << og(ints_+1, -1)<<"01"<<zg(0, 1) << " when \"01\"," << endl // -0.75 when Y = -1.x
                // << tab << tab << og(ints_, -1)<<"1"<<zg(frac_, 1) << " when \"110\"," << endl // -1.0 when Y = -1.0    // Not possible
                << tab << tab << zg(ints_, -1)<<"10"<<zg(1, 1) << " when \"10\"," << endl // 1.0 when Y = 1.0+0.x
                << tab << tab << zg(ints_+1, -1)<<"11"<<zg(0, 1) << " when \"11\"," << endl // 0.75 when Y = 1.5+0.x
                << tab << tab << "\"" << string(r_0_size, '-') << "\" when others;" << endl;
            /**/

            //====================================================================|
            addComment("Main loop of NR algorithm");
            //====================================================================|
#if 1
            /*
            /  r_0: Initial guess for 1/Y
            /  r_{i+1} = r_i*(2-Y*r_i)
            */
            string ri, mi, neg_mi, si, tmp_next_ri, next_ri, round_i;
            int r_size = r_0_size;
            for (int i = 0; i<iters_; ++i){
                addComment("Iteration "+std::to_string(i+1));
                // r_{i+1} = r_i*(2-Y*r_i)
                ri = join("r_", i);
                mi = join("m_", i);
                neg_mi = join("neg_m_", i);
                si = join("s_", i);
                tmp_next_ri = join("tmp_r_", i+1);
                next_ri = join("r_", i+1);
                round_i = join("round_", i);

                // if (i==0){  // First iteration uses smaller multipliers
                //     // Y*r_i
                //     newInstance("IntMultiplier", join("mult_0_", i), "wX="+to_string(width_)+" wY="+to_string(r_0_size)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>"+ri, "R=>"+mi);
                //     // (2-Y*r_i)
                //     vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_+2-1), si, frac_+(r_0_size-1-ints_)+1+1) << " <= " << mi << of(frac_+(r_0_size-1-ints_)+1) << " & (NOT(" << mi << range(frac_+(r_0_size-1-ints_), 0) << ") + 1) ;" << endl;
                //     // r_{i+1} = r_i*(2-Y*r_i)
                //     newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(r_0_size)+" wY="+to_string(frac_+(r_0_size-1-ints_)+1+1)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+ri+",Y=>"+si, "R=>"+tmp_next_ri);
                //     // Rounding r_{i+1}
                //     int prodsize = IntMultiplier::prodsize(r_0_size, frac_+(r_0_size-1-ints_)+1+1, true, true);
                //     int zz = prodsize - 1 - ints_;
                //     if (4 <= frac_ && zz <= width_+frac_) {  // Normalize r_{i+1} size to use same mult size in following iterations
                //         vhdl << tab << declare(next_ri, width_+frac_) << " <= " << tmp_next_ri << range(zz-1, 0) << " & " << zg(width_+frac_-zz) << ";" << endl;
                //     }
                //     else {
                //         int discard = zz - (width_+frac_);
                //         vhdl << tab << declare(getTarget()->adderDelay(width_+frac_), next_ri, width_+frac_) << " <= " << tmp_next_ri << range(zz-1, discard) << " + " << tmp_next_ri << of(discard-1) << ";" << endl;
                //     }
                //     // vhdl << tab << declare(next_ri, width_+frac_) << " <= " << tmp_next_ri << range(2*frac_+ints_+frac_, frac_) << " + " << tmp_next_ri << of(frac_-1) <<";" << endl;
                // }
                // else{
                //     // Y*r_i
                //     newInstance("IntMultiplier", join("mult_0_", i), "wX="+to_string(width_)+" wY="+to_string(r_size)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>"+ri, "R=>"+mi);
                //     // (2-Y*r_i)
                //     vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_-1), si, width_+frac_) << " <= " << mi << of(2*frac_+ints_+frac_) << " & (NOT(" << mi << range(2*frac_+ints_-1+frac_, frac_) << ") + 1) ;" << endl;
                //     // r_{i+1} = r_i*(2-Y*r_i)
                //     newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(width_+frac_)+" wY="+to_string(width_+frac_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+ri+",Y=>"+si, "R=>"+tmp_next_ri);
                //     // Rounding r_{i+1}
                //     vhdl << tab << declare(getTarget()->adderDelay(width_+frac_), next_ri, width_+frac_) << " <= " << tmp_next_ri << range(2*frac_+ints_+2*frac_, 2*frac_) << " + " << tmp_next_ri << of(2*frac_-1) << ";" << endl;
                //     r_size = width_+frac_;
                // }

                // Y*r_i
                newInstance("IntMultiplier", join("mult_0_", i), "wX="+to_string(width_)+" wY="+to_string(r_size)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>"+ri, "R=>"+mi);
                // (2-Y*r_i)
                // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_-1), si, width_) << " <= " << mi << of(2*frac_+ints_+frac_) << " & (NOT(" << mi << range(2*frac_+ints_-1+frac_, frac_) << ") + 1) ;" << endl;
                // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->eqConstComparatorDelay(width_+2), round_i, 1, false) << " <= '0' when " << mi << range(width_+1, 0) << "=0 else '1';" << endl;
                // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_), si, width_) << " <= ('1' & " << zg(width_-1) << ") - ('0' & " << mi << range(2*width_-1, width_+2) << " & " << round_i << " ) ;" << endl;
                // vhdl << tab << declare(getTarget()->logicDelay(), si, width_) << " <= '0' & NOT(" << mi << range(2*width_-4, width_-2) << ");" << endl;

                vhdl << tab << declare(getTarget()->logicDelay(), si, width_) << " <= '0' & NOT(" << mi << range(width_+r_size-4, r_size-2) << ");" << endl;

                // r_{i+1} = r_i*(2-Y*r_i)
                newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(width_)+" wY="+to_string(r_size)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+si+",Y=>"+ri, "R=>"+tmp_next_ri);
                // Rounding r_{i+1}

                int next_r_s = i==iters_-1 ? width_ : min(r_size+r_0_size-2, width_);
                
                vhdl << tab << declare(next_ri, next_r_s) << " <= " << tmp_next_ri << range(width_+r_size-3, width_+r_size-2-next_r_s) << " + " << tmp_next_ri << of(width_+r_size-2-next_r_s-1) << ";" << endl;
                r_size = next_r_s;
            }

            // string next_ri_round = join(next_ri, "_round");
            // vhdl << tab << declare(next_ri_round, width_) << " <= " << next_ri << range(width_+frac_-1, frac_) << " + " << next_ri << of(frac_-1) << ";" << endl;
            //====================================================================|
            addFullComment("Multiply dividend by reciprocal of divisor");
            //====================================================================|
            newInstance("IntMultiplier",
                        "SignificandMultiplication",
                        // "wX="+to_string(width_)+" wY="+to_string(width_+frac_)+" signedIO=true"+" dspThreshold="+to_string(dspOccupationThreshold),
                        "wX="+to_string(width_)+" wY="+to_string(width_)+" signedIO=true"+" dspThreshold="+to_string(dspOccupationThreshold),
                        "X=>X,Y=>"+next_ri,//_round,
                        "R=>quotient_tmp"
                        );
            // vhdl << tab << declare("quotient", 2*width_) << " <= quotient_tmp" << range(2*width_+frac_-1, frac_) << " + quotient_tmp" << of(frac_-1) << ";" << endl;
            vhdl << tab << declare("quotient", 2*width_) << " <= quotient_tmp;" << endl;

#else
            /*
            /  r_0: Initial guess for 1/Y
            /  r_{i+1} = r_i+r_i*(1-Y*r_i)
            */
            string ri, mi, neg_mi, si, tmp_next_ri, next_ri, ext_ri, tmp_ki, ki, c0_i, c1_i;
            vhdl << tab << declare("one_val", width_) << " <= \"01\" & "<<zg(frac_)<<";" << endl;  // revise
            vhdl << tab << declare("one_bit", 1, false) << " <= '1';" << endl;
            vhdl << tab << declare("zero_bit", 1, false) << " <= '0';" << endl;

            for (int i = 0; i<iters_; ++i){
                addComment("Iteration "+std::to_string(i+1));
                // r_{i+1} = r_i+r_i*(1-Y*r_i)
                ri = join("r_", i);
                ext_ri = join("ext_r_", i);
                mi = join("m_", i);
                neg_mi = join("neg_m_", i);
                si = join("s_", i);
                tmp_ki = join("tmp_k_", i);
                ki = join("k_", i);
                c0_i = join("carry0_", i);
                c1_i = join("carry1_", i);
                // tmp_next_ri = join("tmp_r_", i+1);
                next_ri = join("r_", i+1);

                if (i==0){  // First iteration uses smaller multipliers
                    // Y*r_i
                    newInstance("IntMultiplier", join("mult_0_", i), "wX="+to_string(width_)+" wY="+to_string(width_-frac_+2)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>"+ri, "R=>"+mi);
                    // (1-Y*r_i)
                    vhdl << tab << declare(getTarget()->logicDelay(), neg_mi, width_)  << " <= NOT " << mi << range(1+ints_+frac_+1, 2) << ";" << endl;
                    // vhdl << tab << declare(getTarget()->logicDelay(2), c0_i, 1, false)  << " <= NOT(" << mi << of(1) << " OR " << mi << of(0) << ");" << endl;
                    // newInstance("IntAdder", join("add_0_", i), "wIn="+to_string(width_),"X=>one_val,Cin=>"+c0_i+",Y=>" + neg_mi, "R=>"+si);
                    newInstance("IntAdder", join("add_0_", i), "wIn="+to_string(width_),"X=>one_val,Cin=>one_bit,Y=>" + neg_mi, "R=>"+si);
                    // vhdl << tab << declare(getTarget()->logicDelay()+getTarget()->adderDelay(width_+2-1), si, width_+2) << " <= " << mi << of(2+frac_+ints_) << " & (NOT(" << mi << range(2+frac_+ints_-1, 0) << ") + 1) ;" << endl;
                    // r_i*(1-Y*r_i)
                    newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(width_-frac_+2)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+ri+",Y=>"+si, "R=>"+tmp_ki);
                    // Truncate k_{i}
                    vhdl << tab << declare(ki, width_) << " <= " << tmp_ki << range(1+ints_+frac_+1, 2) << ";" << endl;
                    vhdl << tab << declare(ext_ri, width_) << " <= " << ri << " & " << zg(frac_-2) << ";" << endl;
                    // r_{i+1} = r_i*(2-Y*r_i)
                    // vhdl << tab << declare(getTarget()->eqConstComparatorDelay(2), c1_i, 1, false)  << " <= '0' when " << tmp_ki << range(1,0) << " = 0 else '1';" << endl;
                    // newInstance("IntAdder", join("add_1_", i), "wIn="+to_string(width_),"X=>"+ext_ri+",Cin=>"+c1_i+",Y=>" + ki, "R=>"+next_ri);
                    newInstance("IntAdder", join("add_1_", i), "wIn="+to_string(width_),"X=>"+ext_ri+",Cin=>zero_bit,Y=>" + ki, "R=>"+next_ri);
                    // newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(width_-frac_+2)+" wY="+to_string(width_+2)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+ri+",Y=>"+si, "R=>"+tmp_next_ri);
                }
                else{
                    // Y*r_i
                    newInstance("IntMultiplier", join("mult_0_", i), "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>Y,Y=>"+ri, "R=>"+mi);
                    // (1-Y*r_i)
                    vhdl << tab << declare(getTarget()->logicDelay(), neg_mi, width_)  << " <= NOT " << mi << range(2*frac_+ints_, frac_) << ";" << endl;
                    vhdl << tab << declare(getTarget()->eqConstComparatorDelay(frac_), c0_i, 1, false)  << " <= '1' when " << mi << range(frac_-1, 0) << " = 0 else '0';" << endl;
                    newInstance("IntAdder", join("add_0_", i), "wIn="+to_string(width_),"X=>one_val,Cin=>"+c0_i+",Y=>" + neg_mi, "R=>"+si);
                    // r_i*(1-Y*r_i)
                    newInstance("IntMultiplier", join("mult_1_", i), "wX="+to_string(width_)+" wY="+to_string(width_)+ " signedIO=true" +" dspThreshold="+to_string(dspOccupationThreshold),"X=>"+ri+",Y=>"+si, "R=>"+tmp_ki);
                    // Truncate k_{i}
                    vhdl << tab << declare(ki, width_) << " <= " << tmp_ki << range(2*frac_+ints_, frac_) << ";" << endl;
                    // r_{i+1} = r_i*(2-Y*r_i)
                    vhdl << tab << declare(getTarget()->eqConstComparatorDelay(frac_), c1_i, 1, false)  << " <= '0' when " << tmp_ki << range(frac_-1,0) << " = 0 else '1';" << endl;
                    newInstance("IntAdder", join("add_1_", i), "wIn="+to_string(width_),"X=>"+ri+",Cin=>"+c1_i+",Y=>" + ki, "R=>"+next_ri);
                }
            }
            //====================================================================|
            addFullComment("Multiply dividend by reciprocal of divisor");
            //====================================================================|
            newInstance("IntMultiplier",
                        "SignificandMultiplication",
                        "wX="+to_string(width_)+" wY="+to_string(width_)+" signedIO=true"+" dspThreshold="+to_string(dspOccupationThreshold),
                        "X=>X,Y=>"+next_ri,
                        "R=>quotient_tmp"
                        );
            vhdl << tab << declare("quotient", 2*width_) << " <= quotient_tmp;" << endl;

#endif
        }
        } // if (iters_ > 0)
        else {  // Sequential division
            //====================================================================|
            addFullComment("Non-Restoring Division algorithm");
            //====================================================================|
            // vhdl << tab << declare("n_0", width_) << " <= X;" << endl;

            // Ensure that |Dividend| < |Divisor|
            vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_), "X_minus_two", 1, false) << " <= '1' when X = " << og(1, -1) << zg(width_-1, 1) << " else '0';" << endl;
            vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_), "Y_plus_one", 1, false) << " <= '1' when Y = " << zg(1, -1)  << og(1, -2) << zg(width_-2, 1) << " else '0';" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "corner_case", 1, false) << " <= X_minus_two AND Y_plus_one;" << endl;

            vhdl << tab << declare("X_sign", 1, false) << " <= X"<<of(width_-1)<<";" << endl;
            vhdl << tab << declare("div", width_) << " <= Y;" << endl;
            vhdl << tab << declare("div_sign", 1, false) << " <= Y"<<of(width_-1)<<";" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "diff_signs", 1, false) << " <= X_sign XOR div_sign;" << endl;
            vhdl << tab << declare(getTarget()->eqComparatorDelay(width_-1), "div_ge_tmp", 1, false) << " <= '1' when Y"<<range(width_-2,0)<<" > X"<<range(width_-2,0)<<" else '0';" << endl;
            // vhdl << tab << declare(getTarget()->logicDelay(2), "div_gr", 1, false) << " <= div_ge_tmp when diff_signs='0' else (div_ge_tmp XOR X_sign);" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "div_gr", 1, false) << " <= '0';" << endl;   // Always divide by 2

            vhdl << tab << declare(getTarget()->logicDelay(1), "n_0", width_) << " <= X when div_gr='1' else (X" << of(width_-1)<<" & X" << range(width_-1,1)<<");" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(1), "append_0", 1, false) << " <= '0' when div_gr='1' else X" << of(0)<<";" << endl;

            vhdl << tab << declare(getTarget()->logicDelay(), "neg_div", width_) << " <= NOT(Y);" << endl;
            vhdl << tab << declare("one_bit", 1, false) << " <= '1' ;" << endl;
            string r_i, n_i, n_i1, d_i, s_i, z_i, z_i1, rem_z;
#if 1
            for (int i = 0; i<frac_+1+4-1; ++i){
                addComment("Iteration "+std::to_string(i+1));
                r_i = join("r_", i);
                n_i = join("n_", i);
                n_i1 = join("n_", i+1);
                s_i = join("s_", i+1);
                d_i = join("pm_div_", i+1);
                z_i = join("z_", i);
                z_i1 = join("z_", i+1);
                rem_z = join("rem_z_", i);
                // Need an extra bit (Or 2) for dealing wit n_i sign after addition
                // Check if r_i and div have same sign
                vhdl << tab << declare(getTarget()->logicDelay(2), s_i, 1, false) << " <= NOT("<< n_i << of(width_-1) << " XOR div_sign);" << endl;
                // 2*r_i
                if (i==0) {
                    vhdl << tab << declare(r_i, width_) << " <= " << n_i << range(width_-2, 0) << " & append_0;" << endl;
                }
                else {
                    vhdl << tab << declare(r_i, width_) << " <= " << n_i << range(width_-2, 0) << " & '0';" << endl;
                }
                // +div or -div
                vhdl << tab << declare(getTarget()->logicDelay(), d_i, width_) << " <= div when " << s_i << "='0' else neg_div;" << endl;
                // r_{i+1} = 2*r_i +/- div
                newInstance("IntAdder", join("sub_", i), "wIn="+to_string(width_),"X=>"+r_i+",Cin=>"+s_i+",Y=>" + d_i, "R=>"+n_i1);
                vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_), rem_z, 1, false) << " <= '1' when "<<n_i1<<" = 0 else '0';" << endl;
                // Detect if remainder n_i1 is zero
                if (i==0){
                    vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_), z_i1, 1, false) << " <= "<<rem_z<<";" << endl;
                }
                else{
                    vhdl << tab << declare(getTarget()->eqConstComparatorDelay(width_)+getTarget()->logicDelay(2), z_i1, 1, false) << " <= "<<rem_z<<" OR " << z_i << ";" << endl;
                }

                // vhdl << tab << "q_1" << of(width_-i-ints_-2) << " <= " << s_i << ";" << endl;
                // vhdl << tab << "q_2" << of(width_-i-ints_-2) << " <= " << s_i << ";" << endl;
            }
            addComment("Convert the quotient to the digit set {0,1}");
            // q_1 is the positive term (original Q, since -1 digits are stored as zeros)
            // q_2 is the negative term (one's complement on the original Q)
            // quotient = q_1 - q_2
            vhdl << tab << declare("q_1", width_+1+4-1) << " <= " << zg(1+ints_);
            for (int i = 1; i<=frac_+1+4-1; ++i){
                vhdl << " & " << join("s_", i);
            }
            vhdl << " ;" << endl;
            vhdl << tab << declare("q_2", width_+1+4-1) << " <= " << og(1+ints_);
            for (int i = 1; i<=frac_+1+4-1; ++i){
                vhdl << " & " << join("s_", i);
            }
            vhdl << " ;" << endl;
            newInstance("IntAdder",
                        "sub_quotient",
                        "wIn="+to_string(width_+1+4-1),
                        "X=>q_1,Cin=>one_bit,Y=>q_2",
                        "R=>quotient_tmp"
                        );

            addComment("Correction step");
            // In general, if the final remainder and the dividend have opposite signs,
            // a correction step is needed.
            // If the dividend and divisor have the same sign, quotient is corrected by subtracting ulp.
            // If the dividend and divisor have opposite signs, quotient is corrected by adding ulp.
            // Another case that needs correction is a zero remainder in an intermediate step.
            vhdl << tab << declare("rem_sign", 1, false) << " <= " << n_i1 << of(width_-1) << ";" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(3), "rem_div_sign", 1, false) << " <= NOT(" << rem_z << ") AND (rem_sign XOR X_sign);" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "rem_dvr_sign", 1, false) << " <= rem_sign XOR div_sign;" << endl;
            vhdl << tab << declare(getTarget()->logicDelay(2), "div_div_sign", 1, false) << " <= X_sign XOR div_sign;" << endl;
            // vhdl << tab << declare("last_rem_z", 1, false) << " <= " << rem_z << ";" << endl;
            vhdl << tab << declare("interm_zero_rem", 1, false) << " <= NOT(" << rem_z << ") AND " << z_i1 << ";" << endl;

            // vhdl << tab << declare(getTarget()->logicDelay(3), "correct_q", 1, false) << " <= NOT(last_rem_z) AND (rem_div_sign OR any_rem_z);" << endl;



            // vhdl << tab << declare("q_config", 4) << " <= rem_div_sign & div_div_sign & interm_zero_rem & rem_sign;" << endl;

            // vhdl << tab << declare("zz", 1, false) << " <= '0';" << endl;
            // vhdl << tab << "with q_config select " << declare(getTarget()->logicDelay(4), "sub_add_ulp", width_) << " <= " << endl
            //     << tab << tab << og(width_) << " when \"1000\" | \"1001\" | \"1010\" | \"1011\"," << endl // "10XX" subtract ulp when dividend and divisor have the same sign
            //     << tab << tab << zg(width_-1) << " & '1' when \"1100\" | \"1101\" | \"1110\" | \"1111\"," << endl // "11xx" add ulp when dividend and divisor have opposite signs
            //     << tab << tab << og(width_) << " when \"0011\" | \"0111\"," << endl // "0X11" subtract ulp when remainder is negative
            //     << tab << tab << zg(width_-1) << " & '1' when \"0010\" | \"0110\"," << endl // "0X10" add ulp when remainder is positive (!= 0)
            //     << tab << tab << zg(width_) << " when others;" << endl; // Do nothing if not need to correct quotient

            // vhdl << tab << declare("q_config", 4) << " <= 'rem_div_sign' & div_div_sign & interm_zero_rem & rem_dvr_sign;" << endl;
            vhdl << tab << declare("q_config", 5) << " <= '0' & div_div_sign & interm_zero_rem & rem_dvr_sign & corner_case;" << endl;

            vhdl << tab << declare("zz", 1, false) << " <= '0';" << endl;
            vhdl << tab << "with q_config select " << declare(getTarget()->logicDelay(5), "sub_add_ulp", width_+1+4-1) << " <= " << endl
                << tab << tab << og(width_+1+4-1) << " when \"10000\" | \"10010\" | \"10100\" | \"10110\"," << endl // "10XX" subtract ulp when dividend and divisor have the same sign
                << tab << tab << zg(width_+1+3-1) << " & '1' when \"11000\" | \"11010\" | \"11100\" | \"11110\"," << endl // "11xx" add ulp when dividend and divisor have opposite signs
                << tab << tab << og(width_+1+4-1) << " when \"00110\" | \"01110\"," << endl // "0X11" subtract ulp when interm_zero_rem and rem_dvr_sign differ
                << tab << tab << zg(width_+1+3-1) << " & '1' when \"00100\" | \"01100\"," << endl // "0X10" add ulp when remainder rem_dvr_sign are the same
                << tab << tab << og(width_+1+4-1) << " when \"00001\" | \"00011\" | \"00101\" | \"00111\" | \"01001\" | \"01011\" | \"01101\" | \"01111\" | \"10001\" | \"10011\" | \"10101\" | \"10111\" | \"11001\" | \"11011\" | \"11101\" | \"11111\"," << endl // "xxxx1" subtract ulp when X=-2, Y=1
                << tab << tab << zg(width_+1+4-1) << " when others;" << endl; // Do nothing if not need to correct quotient

            newInstance("IntAdder",
                        "correct_quotient",
                        "wIn="+to_string(width_+1+4-1),
                        "X=>quotient_tmp,Cin=>zz,Y=>sub_add_ulp",
                        "R=>quotient_aux"
                        );
            
            vhdl << tab << declare(getTarget()->logicDelay(1), "quotient", width_+4-1) << " <= quotient_aux" << range(width_+4-1,1)<<" when div_gr='1' else (quotient_aux" << range(width_-1+4-1,0)<<");" << endl;


#else   // Extra bits for the quotient
            for (int i = 0; i<frac_+4; ++i){
                addComment("Iteration "+std::to_string(i+1));
                r_i = join("r_", i);
                n_i = join("n_", i);
                n_i1 = join("n_", i+1);
                s_i = join("s_", i+1);
                d_i = join("pm_div_", i+1);
                // Check if r_i and div have same sign
                vhdl << tab << declare(getTarget()->logicDelay(2), s_i, 1, false) << " <= NOT("<< n_i << of(width_-1) << " XOR div" << of(width_-1) << ");" << endl;
                // 2*r_i
                vhdl << tab << declare(r_i, width_) << " <= "<<n_i<<range(width_-2,0)<<" & '0';" << endl;
                // +div or -div
                vhdl << tab << declare(getTarget()->logicDelay(), d_i, width_) << " <= div when " << s_i << "='0' else neg_div;" << endl;
                newInstance("IntAdder", join("sub_", i), "wIn="+to_string(width_),"X=>"+r_i+",Cin=>"+s_i+",Y=>" + d_i, "R=>"+n_i1);

                // vhdl << tab << "q_1" << of(width_-i-ints_-2) << " <= " << s_i << ";" << endl;
                // vhdl << tab << "q_2" << of(width_-i-ints_-2) << " <= " << s_i << ";" << endl;
            }
            addComment("Convert the quotient to the digit set {0,1}");
            // q_1 is the positive term (original Q, since -1 digits are stored as zeros)
            // q_2 is the negative term (one's complement on the original Q)
            // quotient = q_1 - q_2
            vhdl << tab << declare("q_1", width_+4) << " <= " << zg(ints_+1);
            for (int i = 1; i<=frac_+4; ++i){
                vhdl << " & " << join("s_", i);
            }
            vhdl << " ;" << endl;
            vhdl << tab << declare("q_2", width_+4) << " <= " << og(ints_+1);
            for (int i = 1; i<=frac_+4; ++i){
                vhdl << " & " << join("s_", i);
            }
            vhdl << " ;" << endl;
            newInstance("IntAdder",
                        "sub_quotient",
                        "wIn="+to_string(width_+4),
                        "X=>q_1,Cin=>one_bit,Y=>q_2",
                        "R=>quotient"
                        );
#endif

        }
        if (truncate_)
            vhdl << tab << "R <= quotient" << range(2*frac_+ints_, frac_) << ";" << endl;
        else
            vhdl << tab << "R <= quotient;" << endl;

        addFullComment("End of vhdl generation");

    }

	FixDiv::~FixDiv() {}

    void FixDiv::emulate(TestCase *tc) {}

	OperatorPtr FixDiv::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int ints, frac, iters, LUT_out;
		double dspOccupationThreshold=0.0;
        bool truncate, useGoldschmidt;
		UserInterface::parseStrictlyPositiveInt(args, "ints", &ints);
		UserInterface::parseStrictlyPositiveInt(args, "frac", &frac);
		UserInterface::parseFloat(args, "dspThreshold", &dspOccupationThreshold);
		UserInterface::parseBoolean(args, "truncate", &truncate);
		UserInterface::parsePositiveInt(args, "iters", &iters);
		UserInterface::parseBoolean(args, "useGoldschmidt", &useGoldschmidt);
		UserInterface::parsePositiveInt(args, "LUT_out", &LUT_out);
		return new FixDiv(parentOp, target, ints, frac, dspOccupationThreshold, truncate, iters, useGoldschmidt, LUT_out);
	}

	void FixDiv::registerFactory()
	{
		UserInterface::add("FixDiv", // name
						   "A correctly rounded fixed-point divider.",
						   "FixPoint",
						   "", //seeAlso
						   "ints(int): size of integer part in bits;\
						   frac(int): size of fraction part in bits;\
						   dspThreshold(real)=0.0: threshold of relative occupation ratio of a DSP multiplier to be used or not;\
                           truncate(bool)=false: if true, the output has the same length as inputs, and other bits are discarded;\
                           iters(int)=0: if greater than zero, implements an iterative division algorithm;\
						   useGoldschmidt(bool)=false: if true, the Goldschmidt division algorithm is used;\
                           LUT_out(int)=10: if using an iterative algorithm, indicates the output size of LUT for initial guess",
						   "", // htmldoc
						   FixDiv::parseArguments);
	}


} // namespace flopoco
