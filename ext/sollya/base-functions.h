/*

  Copyright 2014-2018 by

  Centre de recherche INRIA Sophia Antipolis Mediterranee,
  Equipes APICS, FACTAS,
  Sophia Antipolis, France

  and by

  Sorbonne Université
  CNRS, Laboratoire d'Informatique de Paris 6, LIP6
  F - 75005 Paris
  France.

  Contributors S. Chevillard, Ch. Lauter

  sylvain.chevillard@ens-lyon.org
  christoph.lauter@ens-lyon.org

  This software is a computer program whose purpose is to provide an
  environment for safe floating-point code development. It is
  particularily targeted to the automatized implementation of
  mathematical floating-point libraries (libm). Amongst other features,
  it offers a certified infinity norm, an automatic polynomial
  implementer and a fast Remez algorithm.

  This software is governed by the CeCILL-C license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-C
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-C license and that you accept its terms.

  This program is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*/

#ifndef BASE_FUNCTIONS_H
#define BASE_FUNCTIONS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gmp.h>
#include <mpfr.h>
#include "mpfi-compat.h"

#define SQRT 0
#define EXP 1
#define LOG 2
#define LOG_2 3
#define LOG_10 4
#define SIN 5
#define COS 6
#define TAN 7
#define ASIN 8
#define ACOS 9
#define ATAN 10
#define SINH 11
#define COSH 12
#define TANH 13
#define ASINH 14
#define ACOSH 15
#define ATANH 16
#define ABS 17
#define DOUBLE 18
#define DOUBLEDOUBLE 19
#define TRIPLEDOUBLE 20
#define ERF 21
#define ERFC 22
#define LOG_1P 23
#define EXP_M1 24
#define DOUBLEEXTENDED 25
#define CEIL 26
#define FLOOR 27
#define SINGLE 28
#define NEARESTINT 39
#define HALFPRECISION 30
#define QUAD 31

#define DECREASING 0         /* Indicates that for any x<y in the domain of f, f(x)>f(y) */
#define NONINCREASING 1      /* Indicates that for any x<y in the domain of f, f(x)>=f(y) */
#define MONOTONICITY_NONE 2  /* Indicates that f has no particular monotonic behavior */
#define NONDECREASING 3      /* Indicates that for any x<y in the domain of f, f(x)<=f(y) */
#define INCREASING 4         /* Indicates that for any x<y in the domain of f, f(x)<f(y) */

typedef struct nodeStruct * basefun_node_ptr_t;

typedef struct baseFunctionStruct baseFunction;
struct baseFunctionStruct
{
  int baseFunctionCode; /* The unique code defined above */
  char *functionName;   /* The name of the function as it is used inside Sollya */
  char *xmlString;      /* The xml code for the function */
  char *mpfrName;       /* The name of the function in mpfr. Used to generate code, e.g. as in implementconstant */
  int handledByImplementconst; /* A boolean. Functions that must not be handled by implementconstant are those functions having a discontinuity at a representable point */
  int isDefinedEverywhere; /* A boolean. True if the function is defined and finite on the whole real line */
  int isDifferentiableEverywhere; /* A boolean. True if the function is differentiable at every interior point of its domain */
  int onlyZeroIsZero; /* A boolean. True if the only real zero of the function is zero */
  int doesNotVanish; /* A boolean. True if the function never takes the value zero on the real line */
  int monotonicity; /* One of the monotonicity code, indicating the behavior of the function on its domain */
  int faithEvaluationOptimizedSupported; /* A boolean. When set to true, the functions getRecurseCutoff and getRecursePrec (see below) are used to practically propagate precision and cutoff across the function call. When set to false, the functions getRecurseCutoff and getRecursePrec are not used and nothing smart is done in order to make the evaluation faster. */
  mp_exp_t (*getRecurseCutoff)(mp_exp_t cutoff, mp_prec_t prec); /* This function implements the behavior of __tryFaithEvaluationOptimizedUnivariateGetRecurseCutoff for the particular case of the present base function f. See the full description of the semantics of this function in infnorm.c. This function is normally called only when faithEvaluationOptimizedSupported is set to true, but must be defined in all cases. Good practice is to use getRecurseCutoff_default as a default in case faithEvaluationOptimizedSupported is set to false. */
  mp_prec_t (*getRecursePrec)(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff); /* This function implements the behavior of __tryFaithEvaluationOptimizedUnivariateGetRecursePrec for the particular case of the present base function f. See the full description of the semantics of this function in infnorm.c. This function is normally called only when faithEvaluationOptimizedSupported is set to true, but must be defined in all cases. Good practice is to use getRecursePrec_default as a default in case faithEvaluationOptimizedSupported is set to false. */
  void (*baseAutodiff)(sollya_mpfi_t *, sollya_mpfi_t, int, int *); /* Computes the vector of the f^(k)(x0)/k!, k=0..n. The last parameter is a silent parameter */
  int (*interval_eval)(sollya_mpfi_t, sollya_mpfi_t); /* Performs an interval evaluation ``à la'' mpfi */
  int (*point_eval)(mpfr_ptr y, mpfr_srcptr x); /* Evaluates f(x) with infinite precision, and rounds the value to the nearest representable number that fits in variable y. Returns a ternary value indicating what has been the direction of the rounding (i.e., this function behaves like an evaluation ``à la'' mpfr, with MPFR_RNDN) */
  int (*try_exact_rational_eval)(mpq_t, mpq_t); /* Returns 1 and sets y=f(x) if it is provably an exact rational. Returns 0 and do nothing otherwise. */
  basefun_node_ptr_t (*diff_expr)(basefun_node_ptr_t ); /* If g is the argument, returns a tree representing diff(f o g) */
  basefun_node_ptr_t (*simplify)(basefun_node_ptr_t ); /* If g is the argument (supposed already simplified as much as possible), returns a tree representing a simplification
                                (without introducing errors) of f(g).
                                Notice that g must either be used (eaten up) to construct the new tree, or be freed by the function  */
  int (*evalsign)(int *, basefun_node_ptr_t ); /* If s and g are the arguments, tries to determine the sign of f(g) assuming that g is a constant expression. In case of success, the sign is assigned to variable s and 1 is returned. Otherwise, s is left unchanged and 0 is returned. Notice that returning 0 without even trying to do anything is a legitimate option. Especially, it is reasonable to assume here that the expression f(g) has already been simplified when possible. For some functions f (e.g., rounding functions, abs) still having f as the head of the tree after simplification usually means that there is no hope of easily getting the sign, in which case, doing nothing is a reasonable option. */
};

#ifndef BASE_FUNCTION_IMPLEMENTATION 

const baseFunction *get_basefun_sqrt();
const baseFunction *get_basefun_exp();
const baseFunction *get_basefun_log();
const baseFunction *get_basefun_log2();
const baseFunction *get_basefun_log10();
const baseFunction *get_basefun_sin();
const baseFunction *get_basefun_cos();
const baseFunction *get_basefun_tan();
const baseFunction *get_basefun_asin();
const baseFunction *get_basefun_acos();
const baseFunction *get_basefun_atan();
const baseFunction *get_basefun_sinh();
const baseFunction *get_basefun_cosh();
const baseFunction *get_basefun_tanh();
const baseFunction *get_basefun_asinh();
const baseFunction *get_basefun_acosh();
const baseFunction *get_basefun_atanh();
const baseFunction *get_basefun_abs();
const baseFunction *get_basefun_double();
const baseFunction *get_basefun_single();
const baseFunction *get_basefun_halfprecision();
const baseFunction *get_basefun_quad();
const baseFunction *get_basefun_doubledouble();
const baseFunction *get_basefun_tripledouble();
const baseFunction *get_basefun_erf();
const baseFunction *get_basefun_erfc();
const baseFunction *get_basefun_log1p();
const baseFunction *get_basefun_expm1();
const baseFunction *get_basefun_doubleextended();
const baseFunction *get_basefun_ceil();
const baseFunction *get_basefun_floor();
const baseFunction *get_basefun_nearestint();

#define basefun_sqrt (get_basefun_sqrt())
#define basefun_exp (get_basefun_exp())
#define basefun_log (get_basefun_log())
#define basefun_log2 (get_basefun_log2())
#define basefun_log10 (get_basefun_log10())
#define basefun_sin (get_basefun_sin())
#define basefun_cos (get_basefun_cos())
#define basefun_tan (get_basefun_tan())
#define basefun_asin (get_basefun_asin())
#define basefun_acos (get_basefun_acos())
#define basefun_atan (get_basefun_atan())
#define basefun_sinh (get_basefun_sinh())
#define basefun_cosh (get_basefun_cosh())
#define basefun_tanh (get_basefun_tanh())
#define basefun_asinh (get_basefun_asinh())
#define basefun_acosh (get_basefun_acosh())
#define basefun_atanh (get_basefun_atanh())
#define basefun_abs (get_basefun_abs())
#define basefun_double (get_basefun_double())
#define basefun_single (get_basefun_single())
#define basefun_halfprecision (get_basefun_halfprecision())
#define basefun_quad (get_basefun_quad())
#define basefun_doubledouble (get_basefun_doubledouble())
#define basefun_tripledouble (get_basefun_tripledouble())
#define basefun_erf (get_basefun_erf())
#define basefun_erfc (get_basefun_erfc())
#define basefun_log1p (get_basefun_log1p())
#define basefun_expm1 (get_basefun_expm1())
#define basefun_doubleextended (get_basefun_doubleextended())
#define basefun_ceil (get_basefun_ceil())
#define basefun_floor (get_basefun_floor())
#define basefun_nearestint (get_basefun_nearestint())

#endif

basefun_node_ptr_t makeSqrt(basefun_node_ptr_t);
basefun_node_ptr_t makeExp(basefun_node_ptr_t);
basefun_node_ptr_t makeLog(basefun_node_ptr_t);
basefun_node_ptr_t makeLog2(basefun_node_ptr_t);
basefun_node_ptr_t makeLog10(basefun_node_ptr_t);
basefun_node_ptr_t makeSin(basefun_node_ptr_t);
basefun_node_ptr_t makeCos(basefun_node_ptr_t);
basefun_node_ptr_t makeTan(basefun_node_ptr_t);
basefun_node_ptr_t makeAsin(basefun_node_ptr_t);
basefun_node_ptr_t makeAcos(basefun_node_ptr_t);
basefun_node_ptr_t makeAtan(basefun_node_ptr_t);
basefun_node_ptr_t makeAbs(basefun_node_ptr_t);
basefun_node_ptr_t makeDouble(basefun_node_ptr_t);
basefun_node_ptr_t makeSingle(basefun_node_ptr_t);
basefun_node_ptr_t makeQuad(basefun_node_ptr_t);
basefun_node_ptr_t makeHalfPrecision(basefun_node_ptr_t);
basefun_node_ptr_t makeDoubledouble(basefun_node_ptr_t);
basefun_node_ptr_t makeTripledouble(basefun_node_ptr_t);
basefun_node_ptr_t makeErf(basefun_node_ptr_t);
basefun_node_ptr_t makeErfc(basefun_node_ptr_t);
basefun_node_ptr_t makeLog1p(basefun_node_ptr_t);
basefun_node_ptr_t makeExpm1(basefun_node_ptr_t);
basefun_node_ptr_t makeDoubleextended(basefun_node_ptr_t);
basefun_node_ptr_t makeCeil(basefun_node_ptr_t);
basefun_node_ptr_t makeFloor(basefun_node_ptr_t);
basefun_node_ptr_t makeNearestInt(basefun_node_ptr_t);
basefun_node_ptr_t makeSinh(basefun_node_ptr_t);
basefun_node_ptr_t makeCosh(basefun_node_ptr_t);
basefun_node_ptr_t makeTanh(basefun_node_ptr_t);
basefun_node_ptr_t makeAsinh(basefun_node_ptr_t);
basefun_node_ptr_t makeAcosh(basefun_node_ptr_t);
basefun_node_ptr_t makeAtanh(basefun_node_ptr_t);

#endif /* ifdef BASE_FUNCTIONS */
