/*

  Copyright 2006-2018 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2),

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

  Sorbonne Université
  CNRS, Laboratoire d'Informatique de Paris 6, LIP6
  F - 75005 Paris
  France

  and by

  Centre de recherche INRIA Sophia Antipolis Mediterranee,
  Equipes APICS, FACTAS,
  Sophia Antipolis, France.

  Contributors Ch. Lauter, M. Joldes, S. Chevillard

  christoph.lauter@ens-lyon.org
  joldes@laas.fr
  sylvain.chevillard@ens-lyon.org

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

#define BASE_FUNCTION_IMPLEMENTATION (1)

#include <stdlib.h>
#include "general.h"
#include "base-functions.h"
#include "autodiff.h"
#include "mpfi-compat.h"
#include "infnorm.h"
#include "double.h"

/******************************************************************************/
/*                                                                            */
/*                     DECLARATIONS OF THE INTERNAL FUNCTIONS                 */
/*                                                                            */
/******************************************************************************/

static void sqrt_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void exp_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void expm1_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void log_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void log1p_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void log2_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void log10_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void sin_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void cos_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void sinh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void cosh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void tan_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void tanh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void atan_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void atanh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void asin_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void acos_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void asinh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void acosh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void erf_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void erfc_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void abs_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void single_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void quad_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void halfprecision_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void double_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void doubledouble_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void tripledouble_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void doubleextended_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void ceil_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void floor_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static void nearestint_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent);
static basefun_node_ptr_t sqrt_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t exp_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t log_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t log2_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t log10_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t sin_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t cos_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t tan_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t asin_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t acos_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t atan_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t sinh_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t cosh_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t tanh_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t asinh_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t acosh_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t atanh_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t erf_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t erfc_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t log1p_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t expm1_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t abs_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t double_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t single_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t quad_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t halfprecision_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t doubledouble_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t tripledouble_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t doubleextended_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t ceil_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t floor_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t nearestint_diff_expr(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_generic(baseFunction const *f, basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_integer_rounding_functions(baseFunction const *f, basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_exp(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_log2(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_log10(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_sin(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_cos(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_tan(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_sinh(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_cosh(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_tanh(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_asinh(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_acosh(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_atanh(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_erf(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_erfc(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_log1p(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_expm1(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_sqrt(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_log(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_asin(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_acos(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_atan(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_abs(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_rounding_function(baseFunction const *rounding_func, basefun_node_ptr_t g, mp_prec_t min_precx, mp_prec_t precy);
static basefun_node_ptr_t simplify_double(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_single(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_quad(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_halfprecision(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_doubledouble(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_tripledouble(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_doubleextended(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_floor(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_ceil(basefun_node_ptr_t g);
static basefun_node_ptr_t simplify_nearestint(basefun_node_ptr_t g);
static int evaluateSignTrigoUnsafe(int *s, basefun_node_ptr_t child, int baseFunctionCode);
static int sqrt_evalsign(int *sign, basefun_node_ptr_t c);
static int log_evalsign(int *sign, basefun_node_ptr_t c);
static int log2_evalsign(int *sign, basefun_node_ptr_t c);
static int log10_evalsign(int *sign, basefun_node_ptr_t c);
static int sin_evalsign(int *sign, basefun_node_ptr_t c);
static int cos_evalsign(int *sign, basefun_node_ptr_t c);
static int tan_evalsign(int *sign, basefun_node_ptr_t c);
static int asin_evalsign(int *sign, basefun_node_ptr_t c);
static int acos_evalsign(int *sign, basefun_node_ptr_t c);
static int odd_increasing_function_evalsign(int *sign, basefun_node_ptr_t c);
static int atan_evalsign(int *sign, basefun_node_ptr_t c);
static int sinh_evalsign(int *sign, basefun_node_ptr_t c);
static int tanh_evalsign(int *sign, basefun_node_ptr_t c);
static int asinh_evalsign(int *sign, basefun_node_ptr_t c);
static int erf_evalsign(int *sign, basefun_node_ptr_t c);
static int expm1_evalsign(int *sign, basefun_node_ptr_t c);
static int positive_function_evalsign(int *sign, basefun_node_ptr_t c);
static int cosh_evalsign(int *sign, basefun_node_ptr_t c);
static int erfc_evalsign(int *sign, basefun_node_ptr_t c);
static int exp_evalsign(int *sign, basefun_node_ptr_t c);
static int acosh_evalsign(int *sign, basefun_node_ptr_t c);
static int atanh_evalsign(int *sign, basefun_node_ptr_t c);
static int void_evalsign(int *sign, basefun_node_ptr_t c);
static int abs_evalsign(int *sign, basefun_node_ptr_t c)  ;
static int double_evalsign(int *sign, basefun_node_ptr_t c);
static int single_evalsign(int *sign, basefun_node_ptr_t c);
static int quad_evalsign(int *sign, basefun_node_ptr_t c);
static int halfprecision_evalsign(int *sign, basefun_node_ptr_t c);
static int doubledouble_evalsign(int *sign, basefun_node_ptr_t c);
static int tripledouble_evalsign(int *sign, basefun_node_ptr_t c);
static int doubleextended_evalsign(int *sign, basefun_node_ptr_t c);
static int log1p_evalsign(int *sign, basefun_node_ptr_t c);
static int ceil_evalsign(int *sign, basefun_node_ptr_t c);
static int floor_evalsign(int *sign, basefun_node_ptr_t c);
static int nearestint_evalsign(int *sign, basefun_node_ptr_t c);
static int sqrt_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int exp_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int log_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int log2_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int log10_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int sin_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int cos_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int tan_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int asin_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int acos_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int atan_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int sinh_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int cosh_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int tanh_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int asinh_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int acosh_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int atanh_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int abs_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int erf_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int erfc_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int log1p_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int expm1_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int ceil_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int floor_point_eval(mpfr_ptr res, mpfr_srcptr x);
static int nearestint_point_eval(mpfr_ptr res, mpfr_srcptr x);
static mp_exp_t getRecurseCutoff_default(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_linearContractive(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_linearNotContractive(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_nonZeroWellConditioned(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_sqrt(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_cos(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_cosh(mp_exp_t cutoff, mp_prec_t prec);
static mp_exp_t getRecurseCutoff_acosh(mp_exp_t cutoff, mp_prec_t prec);
static mp_prec_t getRecursePrec_default(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_sqrt(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_exp(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_log(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_log2(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_log10(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_sin(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_cos(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_tan(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_asin(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_acos(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_atan(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_sinh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_cosh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_tanh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_asinh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_acosh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_atanh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_erf(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_erfc(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_log1p(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static mp_prec_t getRecursePrec_expm1(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff);
static int try_exact_rational_eval_special_cases(baseFunction const *f, mpq_t res, mpq_t a);
static int try_exact_rational_eval_piecewise_constant_functions(baseFunction const *f, mpq_t res, mpq_t a);
static int try_exact_rational_eval_f0_eq_0(mpq_t res, mpq_t a);
static int try_exact_rational_eval_f0_eq_1(mpq_t res, mpq_t a);
static int try_exact_rational_eval_f1_eq_0(mpq_t res, mpq_t a);
static int try_exact_rational_eval_abs(mpq_t res, mpq_t x);
static int try_exact_rational_eval_sqrt(mpq_t res, mpq_t x);
static int try_exact_rational_eval_log2(mpq_t res, mpq_t x);
static int try_exact_rational_eval_log10(mpq_t res, mpq_t x);
static int try_exact_rational_eval_double(mpq_t res, mpq_t x);
static int try_exact_rational_eval_doubledouble(mpq_t res, mpq_t x);
static int try_exact_rational_eval_tripledouble(mpq_t res, mpq_t x);
static int try_exact_rational_eval_doubleextended(mpq_t res, mpq_t x);
static int try_exact_rational_eval_quad(mpq_t res, mpq_t x);
static int try_exact_rational_eval_single(mpq_t res, mpq_t x);
static int try_exact_rational_eval_halfprecision(mpq_t res, mpq_t x);
static int try_exact_rational_eval_floor(mpq_t res, mpq_t x);
static int try_exact_rational_eval_ceil(mpq_t res, mpq_t x);
static int try_exact_rational_eval_nearestint(mpq_t res, mpq_t x);


/******************************************************************************/
/*                                                                            */
/*                     DEFINITIONS OF THE BASEFUN OBJECTS                     */
/*                                                                            */
/******************************************************************************/

static const baseFunction basefun_sqrt[1] = { {
    .baseFunctionCode = SQRT,
    .functionName = "sqrt",
    .xmlString = "<root/>\n",
    .mpfrName = "mpfr_sqrt",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_sqrt,
    .getRecursePrec = getRecursePrec_sqrt,
    .baseAutodiff = sqrt_diff,
    .interval_eval = sollya_mpfi_sqrt,
    .point_eval = sqrt_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_sqrt,
    .diff_expr = sqrt_diff_expr,
    .simplify = simplify_sqrt,
    .evalsign = sqrt_evalsign
  } };

static const baseFunction basefun_exp[1] = { {
    .baseFunctionCode = EXP,
    .functionName = "exp",
    .xmlString = "<exp/>\n",
    .mpfrName = "mpfr_exp",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 1,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_nonZeroWellConditioned,
    .getRecursePrec = getRecursePrec_exp,
    .baseAutodiff = exp_diff,
    .interval_eval = sollya_mpfi_exp,
    .point_eval = exp_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_1,
    .diff_expr = exp_diff_expr,
    .simplify = simplify_exp,
    .evalsign = exp_evalsign
  } };

static const baseFunction basefun_log[1] = { {
    .baseFunctionCode = LOG,
    .functionName = "log",
    .xmlString = "<ln/>\n",
    .mpfrName = "mpfr_log",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_default,
    .getRecursePrec = getRecursePrec_log,
    .baseAutodiff = log_diff,
    .interval_eval = sollya_mpfi_log,
    .point_eval = log_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f1_eq_0,
    .diff_expr = log_diff_expr,
    .simplify = simplify_log,
    .evalsign = log_evalsign
  } };

static const baseFunction basefun_log2[1] = { {
    .baseFunctionCode = LOG_2,
    .functionName = "log2",
    .xmlString = "<log/><logbase><cn>2</cn></logbase>\n",
    .mpfrName = "mpfr_log2",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_default,
    .getRecursePrec = getRecursePrec_log2,
    .baseAutodiff = log2_diff,
    .interval_eval = sollya_mpfi_log2,
    .point_eval = log2_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_log2,
    .diff_expr = log2_diff_expr,
    .simplify = simplify_log2,
    .evalsign = log2_evalsign
  } };

static const baseFunction basefun_log10[1] = { {
    .baseFunctionCode = LOG_10,
    .functionName = "log10",
    .xmlString = "<log/>\n",
    .mpfrName = "mpfr_log10",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_default,
    .getRecursePrec = getRecursePrec_log10,
    .baseAutodiff = log10_diff,
    .interval_eval = sollya_mpfi_log10,
    .point_eval = log10_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_log10,
    .diff_expr = log10_diff_expr,
    .simplify = simplify_log10,
    .evalsign = log10_evalsign
  } };

static const baseFunction basefun_sin[1] = { {
    .baseFunctionCode = SIN,
    .functionName = "sin",
    .xmlString = "<sin/>\n",
    .mpfrName = "mpfr_sin",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = MONOTONICITY_NONE,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearContractive,
    .getRecursePrec = getRecursePrec_sin,
    .baseAutodiff = sin_diff,
    .interval_eval = sollya_mpfi_sin,
    .point_eval = sin_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = sin_diff_expr,
    .simplify = simplify_sin,
    .evalsign = sin_evalsign
  } };

static const baseFunction basefun_cos[1] = { {
    .baseFunctionCode = COS,
    .functionName = "cos",
    .xmlString = "<cos/>\n",
    .mpfrName = "mpfr_cos",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = MONOTONICITY_NONE,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_cos,
    .getRecursePrec = getRecursePrec_cos,
    .baseAutodiff = cos_diff,
    .interval_eval = sollya_mpfi_cos,
    .point_eval = cos_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_1,
    .diff_expr = cos_diff_expr,
    .simplify = simplify_cos,
    .evalsign = cos_evalsign
  } };

static const baseFunction basefun_tan[1] = { {
    .baseFunctionCode = TAN,
    .functionName = "tan",
    .xmlString = "<tan/>\n",
    .mpfrName = "mpfr_tan",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = MONOTONICITY_NONE,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearNotContractive,
    .getRecursePrec = getRecursePrec_tan,
    .baseAutodiff = tan_diff,
    .interval_eval = sollya_mpfi_tan,
    .point_eval = tan_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = tan_diff_expr,
    .simplify = simplify_tan,
    .evalsign = tan_evalsign
  } };

static const baseFunction basefun_asin[1] = { {
    .baseFunctionCode = ASIN,
    .functionName = "asin",
    .xmlString = "<arcsin/>\n",
    .mpfrName = "mpfr_asin",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearNotContractive,
    .getRecursePrec = getRecursePrec_asin,
    .baseAutodiff = asin_diff,
    .interval_eval = sollya_mpfi_asin,
    .point_eval = asin_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = asin_diff_expr,
    .simplify = simplify_asin,
    .evalsign = asin_evalsign
  } };

static const baseFunction basefun_acos[1] = { {
    .baseFunctionCode = ACOS,
    .functionName = "acos",
    .xmlString = "<arccos/>\n",
    .mpfrName = "mpfr_acos",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = DECREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_nonZeroWellConditioned,
    .getRecursePrec = getRecursePrec_acos,
    .baseAutodiff = acos_diff,
    .interval_eval = sollya_mpfi_acos,
    .point_eval = acos_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f1_eq_0,
    .diff_expr = acos_diff_expr,
    .simplify = simplify_acos,
    .evalsign = acos_evalsign
  } };

static const baseFunction basefun_atan[1] = { {
    .baseFunctionCode = ATAN,
    .functionName = "atan",
    .xmlString = "<arctan/>\n",
    .mpfrName = "mpfr_atan",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearContractive,
    .getRecursePrec = getRecursePrec_atan,
    .baseAutodiff = atan_diff,
    .interval_eval = sollya_mpfi_atan,
    .point_eval = atan_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = atan_diff_expr,
    .simplify = simplify_atan,
    .evalsign = atan_evalsign
  } };

static const baseFunction basefun_sinh[1] = { {
    .baseFunctionCode = SINH,
    .functionName = "sinh",
    .xmlString = "<sinh/>\n",
    .mpfrName = "mpfr_sinh",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearNotContractive,
    .getRecursePrec = getRecursePrec_sinh,
    .baseAutodiff = sinh_diff,
    .interval_eval = sollya_mpfi_sinh,
    .point_eval = sinh_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = sinh_diff_expr,
    .simplify = simplify_sinh,
    .evalsign = sinh_evalsign
  } };

static const baseFunction basefun_cosh[1] = { {
    .baseFunctionCode = COSH,
    .functionName = "cosh",
    .xmlString = "<cosh/>\n",
    .mpfrName = "mpfr_cosh",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 1,
    .monotonicity = MONOTONICITY_NONE,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_cosh,
    .getRecursePrec = getRecursePrec_cosh,
    .baseAutodiff = cosh_diff,
    .interval_eval = sollya_mpfi_cosh,
    .point_eval = cosh_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_1,
    .diff_expr = cosh_diff_expr,
    .simplify = simplify_cosh,
    .evalsign = cosh_evalsign
  } };

static const baseFunction basefun_tanh[1] = { {
    .baseFunctionCode = TANH,
    .functionName = "tanh",
    .xmlString = "<tanh/>\n",
    .mpfrName = "mpfr_tanh",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearContractive,
    .getRecursePrec = getRecursePrec_tanh,
    .baseAutodiff = tanh_diff,
    .interval_eval = sollya_mpfi_tanh,
    .point_eval = tanh_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = tanh_diff_expr,
    .simplify = simplify_tanh,
    .evalsign = tanh_evalsign
  } };

static const baseFunction basefun_asinh[1] = { {
    .baseFunctionCode = ASINH,
    .functionName = "asinh",
    .xmlString = "<arcsinh/>\n",
    .mpfrName = "mpfr_asinh",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearContractive,
    .getRecursePrec = getRecursePrec_asinh,
    .baseAutodiff = asinh_diff,
    .interval_eval = sollya_mpfi_asinh,
    .point_eval = asinh_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = asinh_diff_expr,
    .simplify = simplify_asinh,
    .evalsign = asinh_evalsign
  } };

static const baseFunction basefun_acosh[1] = { {
    .baseFunctionCode = ACOSH,
    .functionName = "acosh",
    .xmlString = "<arccosh/>\n",
    .mpfrName = "mpfr_acosh",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_acosh,
    .getRecursePrec = getRecursePrec_acosh,
    .baseAutodiff = acosh_diff,
    .interval_eval = sollya_mpfi_acosh,
    .point_eval = acosh_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f1_eq_0,
    .diff_expr = acosh_diff_expr,
    .simplify = simplify_acosh,
    .evalsign = acosh_evalsign
  } };

static const baseFunction basefun_atanh[1] = { {
    .baseFunctionCode = ATANH,
    .functionName = "atanh",
    .xmlString = "<arctanh/>\n",
    .mpfrName = "mpfr_atanh",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearNotContractive,
    .getRecursePrec = getRecursePrec_atanh,
    .baseAutodiff = atanh_diff,
    .interval_eval = sollya_mpfi_atanh,
    .point_eval = atanh_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = atanh_diff_expr,
    .simplify = simplify_atanh,
    .evalsign = atanh_evalsign
  } };

static const baseFunction basefun_abs[1] = { {
    .baseFunctionCode = ABS,
    .functionName = "abs",
    .xmlString = "<abs/>\n",
    .mpfrName = "mpfr_abs",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = MONOTONICITY_NONE,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_linearContractive, /* Simply transmit the cutoff unchanged */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = abs_diff,
    .interval_eval = sollya_mpfi_abs,
    .point_eval = abs_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_abs,
    .diff_expr = abs_diff_expr,
    .simplify = simplify_abs,
    .evalsign = abs_evalsign
  } };

static const baseFunction basefun_erf[1] = { {
    .baseFunctionCode = ERF,
    .functionName = "erf",
    .xmlString = "<csymbol definitionURL=\"http://www.openmath.org/CDs/errorFresnelInts.ocd\" encoding=\"OpenMath\">erf</csymbol>\n",
    .mpfrName = "mpfr_erf",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearNotContractive, /* Zero at zero, first Taylor coefficient less than 2, series has alternating sign */
    .getRecursePrec = getRecursePrec_erf,
    .baseAutodiff = erf_diff,
    .interval_eval = sollya_mpfi_erf,
    .point_eval = erf_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = erf_diff_expr,
    .simplify = simplify_erf,
    .evalsign = erf_evalsign
  } };

static const baseFunction basefun_erfc[1] = { {
    .baseFunctionCode = ERFC,
    .functionName = "erfc",
    .xmlString = "<csymbol definitionURL=\"http://www.openmath.org/CDs/errorFresnelInts.ocd\" encoding=\"OpenMath\">erfc</csymbol>\n",
    .mpfrName = "mpfr_erfc",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 1,
    .monotonicity = DECREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_nonZeroWellConditioned,
    .getRecursePrec = getRecursePrec_erfc,
    .baseAutodiff = erfc_diff,
    .interval_eval = sollya_mpfi_erfc,
    .point_eval = erfc_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_1,
    .diff_expr = erfc_diff_expr,
    .simplify = simplify_erfc,
    .evalsign = erfc_evalsign
  } };

static const baseFunction basefun_log1p[1] = { {
    .baseFunctionCode = LOG_1P,
    .functionName = "log1p",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">log1p</csymbol>\n",
    .mpfrName = "mpfr_log1p",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 0,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearContractive,
    .getRecursePrec = getRecursePrec_log1p,
    .baseAutodiff = log1p_diff,
    .interval_eval = sollya_mpfi_log1p,
    .point_eval = log1p_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = log1p_diff_expr,
    .simplify = simplify_log1p,
    .evalsign = log1p_evalsign
  } };

static const baseFunction basefun_expm1[1] = { {
    .baseFunctionCode = EXP_M1,
    .functionName = "expm1",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">expm1</csymbol>\n",
    .mpfrName = "mpfr_expm1",
    .handledByImplementconst = 1,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 1,
    .onlyZeroIsZero = 1,
    .doesNotVanish = 0,
    .monotonicity = INCREASING,
    .faithEvaluationOptimizedSupported = 1,
    .getRecurseCutoff = getRecurseCutoff_linearNotContractive,
    .getRecursePrec = getRecursePrec_expm1,
    .baseAutodiff = expm1_diff,
    .interval_eval = sollya_mpfi_expm1,
    .point_eval = expm1_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_f0_eq_0,
    .diff_expr = expm1_diff_expr,
    .simplify = simplify_expm1,
    .evalsign = expm1_evalsign
  } };

static const baseFunction basefun_double[1] = { {
    .baseFunctionCode = DOUBLE,
    .functionName = "double",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">double</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = double_diff,
    .interval_eval = sollya_mpfi_round_to_double,
    .point_eval = sollya_mpfr_round_to_double,
    .try_exact_rational_eval = try_exact_rational_eval_double,
    .diff_expr = double_diff_expr,
    .simplify = simplify_double,
    .evalsign = double_evalsign
  } };

static const baseFunction basefun_single[1] = { {
    .baseFunctionCode = SINGLE,
    .functionName = "single",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">single</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = single_diff,
    .interval_eval = sollya_mpfi_round_to_single,
    .point_eval = sollya_mpfr_round_to_single,
    .try_exact_rational_eval = try_exact_rational_eval_single,
    .diff_expr = single_diff_expr,
    .simplify = simplify_single,
    .evalsign = single_evalsign
  } };


static const baseFunction basefun_halfprecision[1] = { {
    .baseFunctionCode = HALFPRECISION,
    .functionName = "halfprecision",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">halfprecision</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = halfprecision_diff,
    .interval_eval = sollya_mpfi_round_to_halfprecision,
    .point_eval = sollya_mpfr_round_to_halfprecision,
    .try_exact_rational_eval = try_exact_rational_eval_halfprecision,
    .diff_expr = halfprecision_diff_expr,
    .simplify = simplify_halfprecision,
    .evalsign = halfprecision_evalsign
  } };

static const baseFunction basefun_quad[1] = { {
    .baseFunctionCode = QUAD,
    .functionName = "quad",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">quad</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = quad_diff,
    .interval_eval = sollya_mpfi_round_to_quad,
    .point_eval = sollya_mpfr_round_to_quad,
    .try_exact_rational_eval = try_exact_rational_eval_quad,
    .diff_expr = quad_diff_expr,
    .simplify = simplify_quad,
    .evalsign = quad_evalsign
  } };

static const baseFunction basefun_doubledouble[1] = { {
    .baseFunctionCode = DOUBLEDOUBLE,
    .functionName = "doubledouble",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">doubledouble</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = doubledouble_diff,
    .interval_eval = sollya_mpfi_round_to_doubledouble,
    .point_eval = sollya_mpfr_round_to_doubledouble,
    .try_exact_rational_eval = try_exact_rational_eval_doubledouble,
    .diff_expr = doubledouble_diff_expr,
    .simplify = simplify_doubledouble,
    .evalsign = doubledouble_evalsign
  } };


static const baseFunction basefun_tripledouble[1] = { {
    .baseFunctionCode = TRIPLEDOUBLE,
    .functionName = "tripledouble",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">tripledouble</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = tripledouble_diff,
    .interval_eval = sollya_mpfi_round_to_tripledouble,
    .point_eval = sollya_mpfr_round_to_tripledouble,
    .try_exact_rational_eval = try_exact_rational_eval_tripledouble,
    .diff_expr = tripledouble_diff_expr,
    .simplify = simplify_tripledouble,
    .evalsign = tripledouble_evalsign
  } };

static const baseFunction basefun_doubleextended[1] = { {
    .baseFunctionCode = DOUBLEEXTENDED,
    .functionName = "doubleextended",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">doubleextended</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = doubleextended_diff,
    .interval_eval = sollya_mpfi_round_to_doubleextended,
    .point_eval = sollya_mpfr_round_to_doubleextended,
    .try_exact_rational_eval = try_exact_rational_eval_doubleextended,
    .diff_expr = doubleextended_diff_expr,
    .simplify = simplify_doubleextended,
    .evalsign = doubleextended_evalsign
  } };

static const baseFunction basefun_ceil[1] = { {
    .baseFunctionCode = CEIL,
    .functionName = "ceil",
    .xmlString = "<ceiling/>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = ceil_diff,
    .interval_eval = sollya_mpfi_ceil,
    .point_eval = ceil_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_ceil,
    .diff_expr = ceil_diff_expr,
    .simplify = simplify_ceil,
    .evalsign = ceil_evalsign
  } };

static const baseFunction basefun_floor[1] = { {
    .baseFunctionCode = FLOOR,
    .functionName = "floor",
    .xmlString = "<floor/>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = floor_diff,
    .onlyZeroIsZero = 0,
    .interval_eval = sollya_mpfi_floor,
    .point_eval = floor_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_floor,
    .diff_expr = floor_diff_expr,
    .simplify = simplify_floor,
    .evalsign = floor_evalsign
  } };

static const baseFunction basefun_nearestint[1] = { {
    .baseFunctionCode = NEARESTINT,
    .functionName = "nearestint",
    .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">nearestint</csymbol>\n",
    .mpfrName = "",
    .handledByImplementconst = 0,
    .isDefinedEverywhere = 1,
    .isDifferentiableEverywhere = 0,
    .onlyZeroIsZero = 0,
    .doesNotVanish = 0,
    .monotonicity = NONDECREASING,
    .faithEvaluationOptimizedSupported = 0,
    .getRecurseCutoff = getRecurseCutoff_default, /* could do something smarter */
    .getRecursePrec = getRecursePrec_default,
    .baseAutodiff = nearestint_diff,
    .interval_eval = sollya_mpfi_nearestint,
    .point_eval = nearestint_point_eval,
    .try_exact_rational_eval = try_exact_rational_eval_nearestint,
    .diff_expr = nearestint_diff_expr,
    .simplify = simplify_nearestint,
    .evalsign = nearestint_evalsign
  } };



/******************************************************************************/
/*                                                                            */
/*                    AUTOMATIC DIFFERENTIATION FUNCTIONS                     */
/*                                                                            */
/******************************************************************************/

static void sqrt_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t oneHalf;
  mpfr_init2(oneHalf, 12);
  mpfr_set_d(oneHalf, 0.5, GMP_RNDN);
  constantPower_diff(res, x, oneHalf, n, silent);
  mpfr_clear(oneHalf);
  return;
}

static void exp_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_exp(temp, x);
  for(i=0;i<=n;i++) {
    sollya_mpfi_set(res[i], temp);
    sollya_mpfi_div_ui(temp, temp, i+1);
  }

  sollya_mpfi_clear(temp);
  return;
}

static void expm1_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  exp_diff(res, x, n, silent);
  sollya_mpfi_sub_ui(res[0], res[0], 1);
  return;
}


static void log_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t minusOne;
  mp_prec_t prec;
  int i;

  prec = getToolPrecision();

  sollya_mpfi_log(res[0], x);

  if(n>=1) {
    mpfr_init2(minusOne, prec);
    mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(res+1, x, minusOne, n-1, silent);
    mpfr_clear(minusOne);
  }
  for(i=1;i<=n;i++) sollya_mpfi_div_ui(res[i], res[i], i);
  return;
}

static void log1p_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t minusOne;
  sollya_mpfi_t u;
  int i;
  mp_prec_t prec;

  prec = getToolPrecision();

  sollya_mpfi_log1p(res[0], x);

  if(n>=1) {
    sollya_mpfi_init2(u, prec);
    sollya_mpfi_add_ui(u, x, 1);
    mpfr_init2(minusOne, prec);
    mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(res+1, u, minusOne, n-1, silent);
    mpfr_clear(minusOne);
    sollya_mpfi_clear(u);
  }

  for(i=1;i<=n;i++) sollya_mpfi_div_ui(res[i], res[i], i);

  return;
}

/* log2(x) = log(x) * (1/log(2)) */
static void log2_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t log2;
  mp_prec_t prec;

  prec = getToolPrecision();
  sollya_mpfi_init2(log2, prec);

  sollya_mpfi_set_ui(log2, 2); sollya_mpfi_log(log2, log2);
  log_diff(res,x,n,silent);
  for(i=0;i<=n;i++) sollya_mpfi_div(res[i], res[i], log2);

  sollya_mpfi_clear(log2);
  return;
}

/* idem */
static void log10_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t log10;
  mp_prec_t prec;

  prec = getToolPrecision();
  sollya_mpfi_init2(log10, prec);

  sollya_mpfi_set_ui(log10, 10); sollya_mpfi_log(log10, log10);
  log_diff(res,x,n,silent);
  for(i=0;i<=n;i++) sollya_mpfi_div(res[i], res[i], log10);

  sollya_mpfi_clear(log10);
  return;
}

static void sin_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_sin(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  for(i=2; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);

  if(n>=1) {
    sollya_mpfi_cos(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
    for(i=3; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);
  }

  return;
}

static void cos_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_cos(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  for(i=2; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);

  if(n>=1) {
    sollya_mpfi_sin(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
    for(i=1; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);
  }

  return;
}

static void sinh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_sinh(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);

  if(n>=1) {
    sollya_mpfi_cosh(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  }

  return;
}

static void cosh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_cosh(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);

  if(n>=1) {
    sollya_mpfi_sinh(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  }

  return;
}

/*  u=tan(x), tan^(n) / n! = p_(n)(u) with
    p_0 = u;

    recurrence formula: p_(n+1)(u) = (p_n(u))' / (n+1) = p'_n(u) * (1+u^2) / (n+1)
    -> p_n of degree n+1
*/

static void tan_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array;
  sollya_mpfi_t u;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( (n+2),sizeof(sollya_mpfi_t));

  for (index=0; index<=n+1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
  }
  sollya_mpfi_init2(u, prec);

  sollya_mpfi_tan(u, x);
  sollya_mpfi_set_ui(coeffs_array[0], 0);
  sollya_mpfi_set_ui(coeffs_array[1], 1);

  sollya_mpfi_set(res[0], u);

  for(index=1; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array, coeffs_array, index);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
    /* now it represents p_(index-1)' */

    for(i=index+1; i>=2; i--) {
      sollya_mpfi_add(coeffs_array[i], coeffs_array[i], coeffs_array[i-2]);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, u, index+1);
  }

  for (index=0; index<=n+1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
  }
  sollya_mpfi_clear(u);
  safeFree(coeffs_array);

  return;
}

/*  u=tanh(x), tanh^(n) / n! = p_(n)(u) with
    p_0 = u;

    recurrence formula: p_(n+1)(u) = (p_n(u))' / (n+1) = p'_n(u) * (1-u^2) / (n+1)
    -> p_n of degree n+1
*/

static void tanh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array;
  sollya_mpfi_t u;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( (n+2),sizeof(sollya_mpfi_t));

  for (index=0; index<=n+1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
  }
  sollya_mpfi_init2(u, prec);

  sollya_mpfi_tanh(u, x);
  sollya_mpfi_set_ui(coeffs_array[0], 0);
  sollya_mpfi_set_ui(coeffs_array[1], 1);

  sollya_mpfi_set(res[0], u);

  for(index=1; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array, coeffs_array, index);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
    /* now it represents p_(index-1)' */

    for(i=index+1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], coeffs_array[i-2]);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, u, index+1);
  }

  for (index=0; index<=n+1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
  }
  sollya_mpfi_clear(u);
  safeFree(coeffs_array);

  return;
}

/*
  Some thoughts:
  For all holonomic functions, but especially for those below, it would certainly be more interesting to use the recurrence relation that exists between the successive values of the derivative at point x, instead of using a recurrence on the polynomials giving the closed form of the n-th derivative.

  The fact is that the recurrence on the polynomials is more general, as it gives an expression for the derivative *at any point*, but since we are only interested at evaluating the successive derivatives at one single point, it is probably not a good idea as each polynomial requires O(n) operations to be constructed and evaluated, hence a total O(degree^2) operations. This should be compared with O(degree) operations with the recurrence on the coefficients, but it might be that it is not suitable for interval computations, as it may involve recurrence with two or more terms.

  Moreover, the code could be automatically generated (this of few interest here, but could turn very interesting if it were to be recoded inside COQ for instance) for any holonomic function. Here is something I began to write on this topic, on the example of asin:
  deq := holexprtodiffeq(arcsin(x+x0),y(x));
  rec := diffeqtorec(deq,y(x),a(n));
  recc := op(1, rec); cond := [seq(op(i, rec), i=2..nops(rec))];
  recc := collect(recc,a(n)); # a(n) is asin^(n)(x0)/n! and recc is the recurrence it satisfies
  recc := add(factor(op(i,recc)),i=1..nops(recc));
*/



/* atan_diff : reccurence formula: p_(n+1) = (p'_n * (1+x^2) - 2nx * p_n) / (n+1)
   atan^(0) = atan(x)
   atan^(n) / n! = p_(n)(x)/((1+x^2)^n)
   p_1=1;

   --> degree of p_n is (n-1)
*/
static void atan_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_atan(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_add_ui(u, u, 1);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_add(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1));
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1));
    sollya_mpfi_sub(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, index);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}


/* atanh_diff : reccurence formula: p_(n+1) = (p'_n * (1-x^2) + 2nx * p_n)/ (n+1)
   atanh^(0) = atanh(x)
   atanh^(n)/n! = p_(n)(x)/((1-x^2)^n)
   p_1=1;

   --> degree of p_n is (n-1)
*/
static void atanh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_atanh(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_sub_ui(u, u, 1);
    sollya_mpfi_neg(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1));
      sollya_mpfi_add(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1));
    sollya_mpfi_add(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, index);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}



/* asin_diff : recurrence formula: p_(n+1)= (p'_n * (1-x^2) + (2n-1)x * p_n)/(n+1)
   asin^(0) = asin(x)
   asin^(n) / n! = p_(n)(x) / (1-x^2)^((2n-1)/2)
   p_1=1;

   --> degree of p_n is (n-1)
*/
static void asin_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_asin(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_sub_ui(u, u, 1);
    sollya_mpfi_neg(u, u);
    sollya_mpfi_sqrt(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1)-1);
      sollya_mpfi_add(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1)-1);
    sollya_mpfi_add(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, 2*index-1);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}


/* acos_diff : except for the res[0], all the terms are equal to -asin^(n)(x)/n! */
static void acos_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  asin_diff(res,x,n,silent);

  sollya_mpfi_acos(res[0], x);

  for (i=1; i<=n;i++)  sollya_mpfi_neg(res[i], res[i]);

  return;
}


/* asinh_diff : recurrence formula: p_(n+1) = (p'_n * (1+x^2) - (2n-1)x * p_n) / (n+1)
   asinh^(0) = asinh(x)
   asinh^(n)/n! = p_(n)(x) / (1+x^2)^((2n-1)/2)
   p_1=1;

   --> degree of p_n is (n-1)
*/
static void asinh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_asinh(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_add_ui(u, u, 1);
    sollya_mpfi_sqrt(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_add(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1)-1);
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1)-1);
    sollya_mpfi_sub(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, 2*index-1);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}


/* acosh_diff : recurrence formula: p_(n+1) = (p'_n * (x^2-1) - (2n-1)x * p_n) / (n+1)
   acosh^(0) = acosh(x)
   acosh^(n)/n! = p_(n)(x) / (x^2-1)^((2n-1)/2)
   p_1=1;

   --> degree of p_n is (n-1)
*/
static void acosh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_acosh(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_sub_ui(u, u, 1);
    sollya_mpfi_sqrt(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i-2], coeffs_array_diff[i]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1)-1);
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1)-1);
    sollya_mpfi_add(coeffs_array[1], temp, coeffs_array_diff[1]);
    sollya_mpfi_neg(coeffs_array[1], coeffs_array[1]);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_neg(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, 2*index-1);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}

/* erf^(n)(x)/n! = p_n(x)*e^(-x^2)             */
/* with p_1(x) = 2/sqrt(pi)                    */
/* and p_(n+1)(x) = (p_n'(x) - 2xp_n(x))/(n+1) */
/*  -> degree of p_n is n-1                    */
static void erf_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_erf(res[0], x);

  if(n>=1) {
    sollya_mpfi_const_pi(temp);
    sollya_mpfi_sqrt(temp, temp);
    sollya_mpfi_ui_div(temp, 2, temp);

    sollya_mpfi_sqr(u, x);
    sollya_mpfi_neg(u, u);
    sollya_mpfi_exp(u, u);

    sollya_mpfi_mul(res[1], temp, u);

    sollya_mpfi_set(coeffs_array[0], temp);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=1; i--) {
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2);
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_mul(res[index], res[index], u);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}

static void erfc_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  erf_diff(res, x, n, silent);

  sollya_mpfi_erfc(res[0], x);

  for (i=1; i<=n;i++)  sollya_mpfi_neg(res[i],res[i]);

  return;
}

static void abs_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  mpfr_t temp2;
  mp_prec_t prec;

  prec = getToolPrecision();

  sollya_mpfi_abs(res[0], x);
  if(n >= 1) {
    if (sollya_mpfi_has_zero(x))  sollya_mpfi_interv_si(res[1], -1, 1);
    else sollya_mpfi_set_si(res[1], sollya_mpfi_is_nonneg(x) ? 1 : (-1));
  }

  if(n >= 2) {
    mpfr_init2(temp2, prec);
    mpfr_set_nan(temp2);

    if (!(*silent)) {
      *silent = 1;
      printMessage(1, SOLLYA_MSG_ABS_NOT_TWICE_DIFFERENTIABLE, "Warning: the absolute value is not twice differentiable.\n");
      printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
    }
    for(i=2;i<=n;i++) sollya_mpfi_set_fr(res[i], temp2);
    mpfr_clear(temp2);
  }

  return;
}

static void single_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_SINGLE_NOT_DIFFERENTIABLE, "Warning: the single rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void quad_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_QUAD_NOT_DIFFERENTIABLE, "Warning: the quad rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void halfprecision_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_HALF_NOT_DIFFERENTIABLE, "Warning: the half-precision rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void double_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_DOUBLE_NOT_DIFFERENTIABLE, "Warning: the double rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void doubledouble_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_DOUBLE_DOUBLE_NOT_DIFFERENTIABLE, "Warning: the doubledouble rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void tripledouble_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_TRIPLE_DOUBLE_NOT_DIFFERENTIABLE, "Warning: the tripledouble rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void doubleextended_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_DOUBLEEXTENDED_NOT_DIFFERENTIABLE, "Warning: the doubleextended rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void ceil_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_CEIL_NOT_DIFFERENTIABLE, "Warning: the ceil operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void floor_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_FLOOR_NOT_DIFFERENTIABLE, "Warning: the floor operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

static void nearestint_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_NEARESTINT_NOT_DIFFERENTIABLE, "Warning: the nearestint operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}






/******************************************************************************/
/*                                                                            */
/*                    EXPRESSION DIFFERENTIATION FUNCTIONS                    */
/*                                                                            */
/******************************************************************************/

static basefun_node_ptr_t sqrt_diff_expr(basefun_node_ptr_t g) {
  return makeDiv( differentiateUnsimplified(g),
                  makeMul(makeConstantDouble(2.0),
                          makeSqrt(copyTree(g))
                          )
                  );
}

static basefun_node_ptr_t exp_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeExp(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t log_diff_expr(basefun_node_ptr_t g) {
  return makeDiv( differentiateUnsimplified(g),
                  copyTree(g)
                  );
}

static basefun_node_ptr_t log2_diff_expr(basefun_node_ptr_t g) {
  return makeDiv( differentiateUnsimplified(g),
                  makeMul(copyTree(g), makeLog(makeConstantDouble(2.0)))
                  );
}

static basefun_node_ptr_t log10_diff_expr(basefun_node_ptr_t g) {
  return makeDiv( differentiateUnsimplified(g),
                  makeMul(copyTree(g), makeLog(makeConstantDouble(10.0)))
                  );
}

static basefun_node_ptr_t sin_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeCos(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t cos_diff_expr(basefun_node_ptr_t g) {
  return makeNeg( makeMul( makeSin(copyTree(g)),
                           differentiateUnsimplified(g)
                           )
                  );
}

static basefun_node_ptr_t tan_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeAdd(makeConstantDouble(1.0),
                          makePow(makeTan(copyTree(g)), makeConstantDouble(2.0))
                          ),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t asin_diff_expr(basefun_node_ptr_t g) {
  return makeDiv( differentiateUnsimplified(g),
                  makeSqrt(makeSub(makeConstantDouble(1.0),
                                   makePow(copyTree(g), makeConstantDouble(2.0))
                                   )
                           )
                  );
}

static basefun_node_ptr_t acos_diff_expr(basefun_node_ptr_t g) {
  return makeNeg( makeDiv( differentiateUnsimplified(g),
                           makeSqrt(makeSub(makeConstantDouble(1.0),
                                            makePow(copyTree(g), makeConstantDouble(2.0))
                                            )
                                    )
                           )
                  );
}

static basefun_node_ptr_t atan_diff_expr(basefun_node_ptr_t g) {
  return makeDiv( differentiateUnsimplified(g),
                  makeAdd(makeConstantDouble(1.0),
                          makePow(copyTree(g), makeConstantDouble(2.0))
                          )
                  );
}

static basefun_node_ptr_t sinh_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeCosh(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t cosh_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeSinh(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t tanh_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeSub(makeConstantDouble(1.0),
                          makePow(makeTanh(copyTree(g)), makeConstantDouble(2.0))
                          ),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t asinh_diff_expr(basefun_node_ptr_t g) {
  return makeDiv(differentiateUnsimplified(g),
                 makeSqrt(makeAdd(makeConstantDouble(1.0),
                                  makePow(copyTree(g), makeConstantDouble(2.0))
                                  )
                          )
                 );
}

static basefun_node_ptr_t acosh_diff_expr(basefun_node_ptr_t g) {
  return makeDiv(differentiateUnsimplified(g),
                 makeSqrt(makeSub(makePow(copyTree(g), makeConstantDouble(2.0)),
                                  makeConstantDouble(1.0)
                                  )
                          )
                 );
}

static basefun_node_ptr_t atanh_diff_expr(basefun_node_ptr_t g) {
  return makeDiv(differentiateUnsimplified(g),
                 makeSub(makeConstantDouble(1.0),
                         makePow(copyTree(g), makeConstantDouble(2.0))
                         )
                 );
}

static basefun_node_ptr_t erf_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeDiv(makeConstantDouble(2.0), makeSqrt(makePi())),
                  makeMul(makeExp(makeNeg(makePow(copyTree(g), makeConstantDouble(2.0)))),
                          differentiateUnsimplified(g)
                          )
                  );
}

static basefun_node_ptr_t erfc_diff_expr(basefun_node_ptr_t g) {
  return makeNeg(makeMul( makeDiv(makeConstantDouble(2.0), makeSqrt(makePi())),
                          makeMul(makeExp(makeNeg(makePow(copyTree(g), makeConstantDouble(2.0)))),
                                  differentiateUnsimplified(g)
                                  )
                          )
                 );
}

static basefun_node_ptr_t log1p_diff_expr(basefun_node_ptr_t g) {
  return makeDiv(differentiateUnsimplified(g),
                 makeAdd(makeConstantDouble(1.0),
                         copyTree(g)
                         )
                 );
}

static basefun_node_ptr_t expm1_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeExp(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t abs_diff_expr(basefun_node_ptr_t g) {
  return makeMul( makeDiv( copyTree(g),
                           makeAbs(copyTree(g))
                           ),
                  differentiateUnsimplified(g)
                  );
}

static basefun_node_ptr_t double_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_DOUBLE_NOT_DIFFERENTIABLE,
               "Warning: the double rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t single_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_SINGLE_NOT_DIFFERENTIABLE,
               "Warning: the single rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t quad_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_QUAD_NOT_DIFFERENTIABLE,
               "Warning: the quad rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t halfprecision_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_HALF_NOT_DIFFERENTIABLE,
               "Warning: the half-precision rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t doubledouble_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_DOUBLE_DOUBLE_NOT_DIFFERENTIABLE,
               "Warning: the double-double rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t tripledouble_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_TRIPLE_DOUBLE_NOT_DIFFERENTIABLE,
               "Warning: the triple-double rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t doubleextended_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_DOUBLEEXTENDED_NOT_DIFFERENTIABLE,
               "Warning: the double-extended rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t ceil_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);	printMessage(1,SOLLYA_MSG_CEIL_NOT_DIFFERENTIABLE,
                                     "Warning: the ceil operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t floor_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_FLOOR_NOT_DIFFERENTIABLE,
               "Warning: the floor operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

static basefun_node_ptr_t nearestint_diff_expr(basefun_node_ptr_t g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_NEARESTINT_NOT_DIFFERENTIABLE,
               "Warning: the nearestint operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}






/******************************************************************************/
/*                                                                            */
/*                   FUNCTIONS FOR SIMPLIFYING EXPRESSIONS                    */
/*                                                                            */
/******************************************************************************/

/* Generic algorithm for simplifying an expression of the form f(g) where g is already simplified */
static basefun_node_ptr_t simplify_generic(baseFunction const *f, basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;

  simplified = allocateNode();

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value, tools_precision);
    if ( (f->point_eval(*value, *(accessThruMemRef(g)->value)) == 0) &&
         mpfr_number_p(*value) )  {
      simplified->nodeType = CONSTANT;
      simplified->value = value;
      free_memory(g);
      return simplified;
    }
    else {
      mpfr_clear(*value);
      safeFree(value);
    }
  }

  /* Default case: we do nothing but constructing f(g) */
  simplified->nodeType = UNARY_BASE_FUNC;
  simplified->baseFun = f;
  simplified->child1 = g;
  return simplified;
}

/* Algorithm for simplifying an expression of the form f(g) where f is nearestint|floor|ceil and g is already simplified */
static basefun_node_ptr_t simplify_integer_rounding_functions(baseFunction const *f, basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value = NULL;
  mpfi_t xI, yI, zI;
  rangetype xrange, yrange;
  mp_prec_t prec;
  int ok = 0; /* Boolean: stores whether we manage to exactly evaluate the value of f(g) in variable *value */

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    prec = mpfr_get_prec(*(accessThruMemRef(g)->value));
    prec = (tools_precision<prec)?prec:tools_precision;
    mpfr_init2(*value, prec);
    ok = ( (f->point_eval(*value, *(accessThruMemRef(g)->value)) == 0) &&
           mpfr_number_p(*value) );
  }
  else if (isConstant(g)) {
    mpfi_init2(xI, 12);
    mpfi_init2(yI, 128);
    mpfi_init2(zI, ((tools_precision > 128)? tools_precision : 128));
    mpfi_set_ui(xI, 1);
    /* HACK ALERT: For performance reasons, we will access the internals of an mpfi_t !!! */
    /* Also, since mpfr_t is an array of 1 __mpfr_struct, the address of a __mpfr_struct can
       be both interpreted as a mpfr_t or a (mpfr_t *) */
    xrange.a = (mpfr_t *)(&(xI->left)); xrange.b = (mpfr_t *)(&(xI->right));
    yrange.a = (mpfr_t *)(&(yI->left)); yrange.b = (mpfr_t *)(&(yI->right));
    evaluateRangeFunction(yrange, g, xrange, 128);
    f->interval_eval(zI, yI);
    ok = ( mpfr_number_p(&(zI->left)) &&
           mpfr_number_p(&(zI->right)) &&
           (mpfr_cmp(&(zI->left), &(zI->right)) == 0) );
    if (ok) {
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value, ((tools_precision > 128)? tools_precision : 128));
      mpfr_set(*value, &(zI->left),GMP_RNDN); /* Exact */
    }
    mpfi_clear(xI);
    mpfi_clear(yI);
    mpfi_clear(zI);
  }

  if (ok) {
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    simplified->value = value;
    free_memory(g);
    return simplified;
  }

  /* Else... (*value) cannot be used. Maybe g can be structurally recognized as an integer, even though it is not representable? */
  if (value != NULL) {
    mpfr_clear(*value);
    safeFree(value);
  }
  if (isIntegerConstant(g)) return g;

  /* Default case: we do nothing but constructing f(g) */
  simplified = allocateNode();
  simplified->nodeType = UNARY_BASE_FUNC;
  simplified->baseFun = f;
  simplified->child1 = g;
  return simplified;
}

static basefun_node_ptr_t simplify_exp(basefun_node_ptr_t g) { return simplify_generic(basefun_exp, g); }
static basefun_node_ptr_t simplify_log2(basefun_node_ptr_t g) { return simplify_generic(basefun_log2, g); }
static basefun_node_ptr_t simplify_log10(basefun_node_ptr_t g) { return simplify_generic(basefun_log10, g); }
static basefun_node_ptr_t simplify_sin(basefun_node_ptr_t g) { return simplify_generic(basefun_sin, g); }
static basefun_node_ptr_t simplify_cos(basefun_node_ptr_t g) { return simplify_generic(basefun_cos, g); }
static basefun_node_ptr_t simplify_tan(basefun_node_ptr_t g) { return simplify_generic(basefun_tan, g); }
static basefun_node_ptr_t simplify_sinh(basefun_node_ptr_t g) { return simplify_generic(basefun_sinh, g); }
static basefun_node_ptr_t simplify_cosh(basefun_node_ptr_t g) { return simplify_generic(basefun_cosh, g); }
static basefun_node_ptr_t simplify_tanh(basefun_node_ptr_t g) { return simplify_generic(basefun_tanh, g); }
static basefun_node_ptr_t simplify_asinh(basefun_node_ptr_t g) { return simplify_generic(basefun_asinh, g); }
static basefun_node_ptr_t simplify_acosh(basefun_node_ptr_t g) { return simplify_generic(basefun_acosh, g); }
static basefun_node_ptr_t simplify_atanh(basefun_node_ptr_t g) { return simplify_generic(basefun_atanh, g); }
static basefun_node_ptr_t simplify_erf(basefun_node_ptr_t g) { return simplify_generic(basefun_erf, g); }
static basefun_node_ptr_t simplify_erfc(basefun_node_ptr_t g) { return simplify_generic(basefun_erfc, g); }
static basefun_node_ptr_t simplify_log1p(basefun_node_ptr_t g) { return simplify_generic(basefun_log1p, g); }
static basefun_node_ptr_t simplify_expm1(basefun_node_ptr_t g) { return simplify_generic(basefun_expm1, g); }

static basefun_node_ptr_t simplify_sqrt(basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;
  mp_prec_t prec, p;

  simplified = allocateNode();
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    prec = tools_precision;
    p = mpfr_get_prec(*(accessThruMemRef(g)->value));
    if (p > prec) prec = p;
    prec += 10;
    if (prec > 256 * tools_precision) prec = 256 * tools_precision;
    mpfr_init2(*value,prec);
    simplified->value = value;
    if ((mpfr_sqrt(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
        (!mpfr_number_p(*value))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_sqrt;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if ((accessThruMemRef(g)->nodeType == POW) &&
        (accessThruMemRef(accessThruMemRef(g)->child2)->nodeType == CONSTANT) &&
        (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(g)->child2)->value),2.0) == 0.0) &&
        (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(g)->child2)->value)))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_abs;
      simplified->child1 = copyTree(accessThruMemRef(g)->child1);
      free_memory(g);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_sqrt;
      simplified->child1 = g;
    }
  }
  return simplified;
}


static basefun_node_ptr_t simplify_log(basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;

  simplified = allocateNode();
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    simplified->value = value;
    if ((mpfr_log(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
        (!mpfr_number_p(*value))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_log;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if ((accessThruMemRef(g)->nodeType == ADD) &&
        (accessThruMemRef(accessThruMemRef(g)->child1)->nodeType == CONSTANT) &&
        (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(g)->child1)->value),1.0) == 0.0) &&
        (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(g)->child1)->value)))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_log1p;
      simplified->child1 = copyTree(accessThruMemRef(g)->child2);
      free_memory(g);
    } else {
      if ((accessThruMemRef(g)->nodeType == ADD) &&
          (accessThruMemRef(accessThruMemRef(g)->child2)->nodeType == CONSTANT) &&
          (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(g)->child2)->value),1.0) == 0.0) &&
          (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(g)->child2)->value)))) {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_log1p;
        simplified->child1 = copyTree(accessThruMemRef(g)->child1);
        free_memory(g);
      } else {
        if ( (accessThruMemRef(g)->nodeType == UNARY_BASE_FUNC) &&
             (accessThruMemRef(g)->baseFun->baseFunctionCode == EXP) ) {
          freeNode(simplified);
          simplified = copyTree(accessThruMemRef(g)->child1);
          free_memory(g);
        } else {
          simplified->nodeType = UNARY_BASE_FUNC;
          simplified->baseFun = basefun_log;
          simplified->child1 = g;
        }
      }
    }
  }
  return simplified;
}




static basefun_node_ptr_t simplify_asin(basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;
  mpfr_t temp;

  simplified = allocateNode();
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    mpfr_init2(temp,53);
    mpfr_set_ui(temp,1,GMP_RNDN);
    if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
      free_memory(g);
      simplified->child2 = allocateNode();
      simplified->child2->nodeType = CONSTANT;
      simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(simplified->child2->value),53);
      mpfr_set_si(*(simplified->child2->value),2,GMP_RNDN);
      simplified->nodeType = DIV;
      simplified->child1 = allocateNode();
      simplified->child1->nodeType = PI_CONST;
    } else {
      mpfr_set_si(temp,-1,GMP_RNDN);
      if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
        free_memory(g);
        simplified->child2 = allocateNode();;
        simplified->child2->nodeType = CONSTANT;
        simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(simplified->child2->value),53);
        mpfr_set_si(*(simplified->child2->value),-2,GMP_RNDN);
        simplified->nodeType = DIV;
        simplified->child1 = allocateNode();
        simplified->child1->nodeType = PI_CONST;
      } else {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,tools_precision);
        simplified->value = value;
        if ((mpfr_asin(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
            (!mpfr_number_p(*value))) {
          simplified->nodeType = UNARY_BASE_FUNC;
          simplified->baseFun = basefun_asin;
          simplified->child1 = g;
          mpfr_clear(*value);
          safeFree(value);
        } else {
          free_memory(g);
        }
      }
    }
    mpfr_clear(temp);
  } else {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_asin;
    simplified->child1 = g;
  }
  return simplified;
}


static basefun_node_ptr_t simplify_acos(basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;
  mpfr_t temp;

  simplified = allocateNode();
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    mpfr_init2(temp,53);
    mpfr_set_si(temp,-1,GMP_RNDN);
    if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
      free_memory(g);
      simplified->nodeType = PI_CONST;
    } else {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      if ((mpfr_acos(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
          (!mpfr_number_p(*value))) {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_acos;
        simplified->child1 = g;
        mpfr_clear(*value);
        safeFree(value);
      } else {
        free_memory(g);
      }
    }
    mpfr_clear(temp);
  } else {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_acos;
    simplified->child1 = g;
  }
  return simplified;
}


static basefun_node_ptr_t simplify_atan(basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;
  mpfr_t temp;

  simplified = allocateNode();
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    mpfr_init2(temp,53);
    mpfr_set_ui(temp,1,GMP_RNDN);
    if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
      free_memory(g);
      simplified->child2 = allocateNode();
      simplified->child2->nodeType = CONSTANT;
      simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(simplified->child2->value),53);
      mpfr_set_si(*(simplified->child2->value),4,GMP_RNDN);
      simplified->nodeType = DIV;
      simplified->child1 = allocateNode();
      simplified->child1->nodeType = PI_CONST;
    } else {
      mpfr_set_si(temp,-1,GMP_RNDN);
      if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
        free_memory(g);
        simplified->child2 = allocateNode();
        simplified->child2->nodeType = CONSTANT;
        simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(simplified->child2->value),53);
        mpfr_set_si(*(simplified->child2->value),-4,GMP_RNDN);
        simplified->nodeType = DIV;
        simplified->child1 = allocateNode();
        simplified->child1->nodeType = PI_CONST;
      } else {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,tools_precision);
        simplified->value = value;
        if ((mpfr_atan(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
            (!mpfr_number_p(*value))) {
          simplified->nodeType = UNARY_BASE_FUNC;
          simplified->baseFun = basefun_atan;
          simplified->child1 = g;
          mpfr_clear(*value);
          safeFree(value);
        } else {
          free_memory(g);
        }
      }
    }
    mpfr_clear(temp);
  } else {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_atan;
    simplified->child1 = g;
  }
  return simplified;
}



static basefun_node_ptr_t simplify_abs(basefun_node_ptr_t g) {
  basefun_node_ptr_t simplified;
  mpfr_t *value;
  int s = 17;

  if ( (accessThruMemRef(g)->nodeType == CONSTANT) &&
       mpfr_number_p(*(accessThruMemRef(g)->value)) ) {
    s = mpfr_sgn(*(accessThruMemRef(g)->value));
    if (s >= 0) return g;
    else {
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value, mpfr_get_prec(*(accessThruMemRef(g)->value)));
      mpfr_abs(*value, *(accessThruMemRef(g)->value), GMP_RNDN); /* exact */
      simplified = allocateNode();
      simplified->nodeType = CONSTANT;
      simplified->value = value;
      free_memory(g);
      return simplified;
    }
  }

  /* else... */
  if (evaluateSign(&s, g)) { /* If g is not a constant expression, returns 0 */
    /* Here, g is a constant expression and we managed to find its sign */
    if (s == 0) {
      free_memory(g);
      return makeConstantDouble(0.0);
    }
    else {
      if (s > 0) return g;
      else return makeNeg(g);
    }
  }

  /* abs(abs(something)) --> abs(something) */
  if ( (accessThruMemRef(g)->nodeType == UNARY_BASE_FUNC) &&
       (accessThruMemRef(g)->baseFun->baseFunctionCode == ABS) ) {
    simplified = makeAbs(copyTree(accessThruMemRef(g)->child1));
    free_memory(g);
    return simplified;
  }

  /* Default case: we do nothing but constructing f(g) */
  return makeAbs(g);
}


static basefun_node_ptr_t simplify_rounding_function(baseFunction const *rounding_func, basefun_node_ptr_t g, mp_prec_t min_precx, mp_prec_t precy) {
  basefun_node_ptr_t simplified;
  mpfr_t *value = NULL;
  int ok = 0; /* Boolean: stores whether we manage to exactly evaluate the value of f(g) in variable *value */
  mpfi_t xI, yI, zI;
  rangetype xrange, yrange;
  mp_prec_t precx = (tools_precision > min_precx) ? tools_precision : min_precx;

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value, precx);
    rounding_func->point_eval(*value, *(accessThruMemRef(g)->value));
    ok = mpfr_number_p(*value);
  }
  else {
    if (isConstant(g)) {
      mpfi_init2(xI, 12);
      mpfi_init2(yI, precy);
      mpfi_init2(zI, precx);
      mpfi_set_ui(xI, 1);
      /* HACK ALERT: For performance reasons, we will access the internals of an mpfi_t !!! */
      /* Also, since mpfr_t is an array of 1 __mpfr_struct, the address of a __mpf
	 be both interpreted as a mpfr_t or a (mpfr_t *) */
      xrange.a = (mpfr_t *)(&(xI->left)); xrange.b = (mpfr_t *)(&(xI->right));
      yrange.a = (mpfr_t *)(&(yI->left)); yrange.b = (mpfr_t *)(&(yI->right));
      evaluateRangeFunction(yrange, g, xrange, precy);
      rounding_func->interval_eval(zI, yI);
      ok = ( mpfr_number_p(&(zI->left)) &&
             mpfr_number_p(&(zI->right)) &&
             (mpfr_cmp(&(zI->left), &(zI->right)) == 0) );
      if (ok) {
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value, precx);
        mpfr_set(*value, &(zI->left), GMP_RNDN); /* Exact */
      }
      mpfi_clear(xI);
      mpfi_clear(yI);
      mpfi_clear(zI);
    }
  }

  simplified = allocateNode();
  if (ok) {
    simplified->nodeType = CONSTANT;
    simplified->value = value;
    free_memory(g);
    return simplified;
  }

  /* Default case. If value has been allocated we clean it, and we return the
     standard simplified expression */
  if (value != NULL) {
    mpfr_clear(*value);
    safeFree(value);
  }
  simplified->nodeType = UNARY_BASE_FUNC;
  simplified->baseFun = rounding_func;
  simplified->child1 = g;
  return simplified;
}

static basefun_node_ptr_t simplify_double(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_double, g, 64, 96); }
static basefun_node_ptr_t simplify_single(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_single, g, 64, 96); }
static basefun_node_ptr_t simplify_quad(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_quad, g, 128, 160); }
static basefun_node_ptr_t simplify_halfprecision(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_halfprecision, g, 64, 32); }
static basefun_node_ptr_t simplify_doubledouble(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_doubledouble, g, 129, 161); }
static basefun_node_ptr_t simplify_tripledouble(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_tripledouble, g, 200, 224); }
static basefun_node_ptr_t simplify_doubleextended(basefun_node_ptr_t g) { return simplify_rounding_function(basefun_doubleextended, g, 128, 96); }
static basefun_node_ptr_t simplify_floor(basefun_node_ptr_t g) { return simplify_integer_rounding_functions(basefun_floor, g); }
static basefun_node_ptr_t simplify_ceil(basefun_node_ptr_t g) { return simplify_integer_rounding_functions(basefun_ceil, g); }
static basefun_node_ptr_t simplify_nearestint(basefun_node_ptr_t g) { return simplify_integer_rounding_functions(basefun_nearestint, g); }




/******************************************************************************/
/*                                                                            */
/*                     FUNCTIONS FOR EVALUATING THE SIGN                      */
/*                                                                            */
/******************************************************************************/

static int evaluateSignTrigoUnsafe(int *s, basefun_node_ptr_t child, int baseFunctionCode) {
  mpfr_t value, value2;
  mpfr_t piHalf;
  mpfr_t dummyX;
  int okay, res;
  basefun_node_ptr_t tempNode;
  int signA;

  okay = 0;

  mpfr_init2(value,defaultprecision);
  mpfr_init2(piHalf,defaultprecision);
  mpfr_init2(dummyX,12);
  mpfr_set_ui(dummyX,1,GMP_RNDN);
  if (evaluateFaithful(value, child, dummyX, defaultprecision) &&
      mpfr_number_p(value)) {
    mpfr_const_pi(piHalf,GMP_RNDN);
    mpfr_div_2ui(piHalf,piHalf,1,GMP_RNDN);
    mpfr_div(value,value,piHalf,GMP_RNDN);
    mpfr_rint(value,value,GMP_RNDN);
    mpfr_div_2ui(value,value,1,GMP_RNDN);
    /* Here, diff is approximately value * pi
       and value * 2 is an integer
    */
    tempNode = makeMul(makeConstant(value),makePi());
    if (compareConstant(&signA, child, tempNode, NULL, 0)) {
      if (signA == 0) {
	/* Here, we have proven that child is equal to value * pi
	 */
	mpfr_init2(value2,defaultprecision);
	mpfr_rint(value2,value,GMP_RNDN);      /* exact, same precision */
	mpfr_sub(value,value,value2,GMP_RNDN); /* exact, Sterbenz */
	/* Here, we know that child is equal to (n + value) * pi for
	   some integer n. We know that value can only be 0 or +/- 0.5
	*/
	switch (baseFunctionCode) {
	case SIN:
	  /* sin is zero for all n * pi, n in Z */
	  if (mpfr_zero_p(value)) {
	    okay = 1;
	    res = 0;
	  }
	  break;
	case COS:
	  /* cos is zero for all (n + 1/2) * pi, n in Z */
	  if (!mpfr_zero_p(value)) {
	    okay = 1;
	    res = 0;
	  }
	  break;
	case TAN:
	  /* tan is zero for all n * pi, n in Z */
	  if (mpfr_zero_p(value)) {
	    okay = 1;
	    res = 0;
	  }
	  break;
	default:
	  sollyaFprintf(stderr,"Error: evaluateSignTrigoUnsafe: unknown identifier (%d) in the tree\n",baseFunctionCode);
	  exit(1);
	}
	mpfr_clear(value2);
      }
    }
    free_memory(tempNode);
  }
  mpfr_clear(dummyX);
  mpfr_clear(piHalf);
  mpfr_clear(value);

  if (okay) *s = res;
  return okay;
}

static int sqrt_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA >= 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

static int log_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA, okayB, signB;
  basefun_node_ptr_t  tempNode;
  tempNode = makeConstantDouble(1.0);
  okayA = compareConstant(&signA, c, tempNode, NULL, 0);
  free_memory(tempNode);
  okayB = evaluateSign(&signB, c);
  if (okayA && okayB && (signB > 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

static int log2_evalsign(int *sign, basefun_node_ptr_t c) { return log_evalsign(sign, c); }
static int log10_evalsign(int *sign, basefun_node_ptr_t c) { return log_evalsign(sign, c); }

static int sin_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA == 0)) {
    *sign = 0;
    return 1;
  }
  /* else */
  return evaluateSignTrigoUnsafe(sign, c, SIN);
}

static int cos_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA == 0)) {
    *sign = 0;
    return 1;
  }
  /* else */
  return evaluateSignTrigoUnsafe(sign, c, COS);
}

static int tan_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA == 0)) {
    *sign = 0;
    return 1;
  }
  /* else */
  return evaluateSignTrigoUnsafe(sign, c, TAN);
}


static int asin_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA, okayB, signB;
  basefun_node_ptr_t tempNode, tempNode2;
  okayA = evaluateSign(&signA, c);
  tempNode = makeAbs(copyTree(c));
  tempNode2 = makeConstantDouble(1.0);
  okayB = compareConstant(&signB, tempNode, tempNode2, NULL, 0);
  free_memory(tempNode);
  free_memory(tempNode2);

  if (okayA && okayB && (signB <= 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

static int acos_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA, okayB, signB, okayC, signC;
  basefun_node_ptr_t tempNode, tempNode2;
  okayA = evaluateSign(&signA, c);
  tempNode = makeAbs(copyTree(c));
  tempNode2 = makeConstantDouble(1.0);
  okayB = compareConstant(&signB, tempNode, tempNode2, NULL, 0);
  okayC = compareConstant(&signC, c, tempNode2, NULL, 0);
  free_memory(tempNode);
  free_memory(tempNode2);

  if (okayA && okayB && okayC && (signB <= 0)) {
    if (signC == 0) *sign = 0; else *sign = 1;
    return 1;
  }
  return 0;
}

static int odd_increasing_function_evalsign(int *sign, basefun_node_ptr_t c) {
  return evaluateSign(sign, c);
}

static int atan_evalsign(int *sign, basefun_node_ptr_t c) { return odd_increasing_function_evalsign(sign, c); }
static int sinh_evalsign(int *sign, basefun_node_ptr_t c) { return odd_increasing_function_evalsign(sign, c); }
static int tanh_evalsign(int *sign, basefun_node_ptr_t c) { return odd_increasing_function_evalsign(sign, c); }
static int asinh_evalsign(int *sign, basefun_node_ptr_t c) { return odd_increasing_function_evalsign(sign, c); }
static int erf_evalsign(int *sign, basefun_node_ptr_t c) { return odd_increasing_function_evalsign(sign, c); }
static int expm1_evalsign(int *sign, basefun_node_ptr_t c) { return odd_increasing_function_evalsign(sign, c); }

static int positive_function_evalsign(int *sign, basefun_node_ptr_t c) {
  int signA;
  if (evaluateSign(&signA, c)) {
    *sign = 1;
    return 1;
  }
  return 0;
}

static int cosh_evalsign(int *sign, basefun_node_ptr_t c) { return positive_function_evalsign(sign, c); }
static int erfc_evalsign(int *sign, basefun_node_ptr_t c) { return positive_function_evalsign(sign, c); }
static int exp_evalsign(int *sign, basefun_node_ptr_t c) {  return positive_function_evalsign(sign, c); }

static int acosh_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA;
  basefun_node_ptr_t tempNode;
  tempNode = makeConstantDouble(1.0);
  okayA = compareConstant(&signA, c, tempNode, NULL, 0);
  free_memory(tempNode);
  if (okayA && (signA >= 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

static int atanh_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA, okayB, signB;
  basefun_node_ptr_t tempNode, tempNode2;
  okayA = evaluateSign(&signA, c);
  tempNode = makeAbs(copyTree(c));
  tempNode2 = makeConstantDouble(1.0);
  okayB = compareConstant(&signB, tempNode, tempNode2, NULL, 0);
  free_memory(tempNode);
  free_memory(tempNode2);
  if (okayA && okayB && (signB < 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

static int void_evalsign(int *sign, basefun_node_ptr_t c) {
  UNUSED_PARAM(sign);
  UNUSED_PARAM(c);
  return 0;
}


static int abs_evalsign(int *sign, basefun_node_ptr_t c)  { return void_evalsign(sign, c); }
static int double_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }
static int single_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }
static int quad_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }
static int halfprecision_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }
static int doubledouble_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }
static int tripledouble_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }
static int doubleextended_evalsign(int *sign, basefun_node_ptr_t c) { return void_evalsign(sign, c); }

static int log1p_evalsign(int *sign, basefun_node_ptr_t c) {
  int okayA, signA, okayB, signB;
  basefun_node_ptr_t tempNode;
  tempNode = makeConstantDouble(-1.0);
  okayA = compareConstant(&signA, c, tempNode, NULL, 0);
  okayB = evaluateSign(&signB, c);
  free_memory(tempNode);
  if (okayA && okayB && (signA > 0)) {
    *sign = signB;
    return 1;
  }
  return 0;
}

static int ceil_evalsign(int *sign, basefun_node_ptr_t c) {
  int okay, s, okayA, signA, okayB, signB;
  basefun_node_ptr_t tempNode;
  okay = 0;
  okayA = evaluateSign(&signA, c);
  tempNode = makeConstantDouble(-1.0);
  if (okayA)
    okayB = compareConstant(&signB, c, tempNode, NULL, 0);
  else
    okayB = 0;
  if (okayA && okayB) {
    okay = 1;
    if (signB <= 0) {
      s = -1;
    } else {
      if (signA <= 0) {
        s = 0;
      } else {
        s = 1;
      }
    }
  }
  free_memory(tempNode);
  if (okay) *sign = s;
  return okay;
}

static int floor_evalsign(int *sign, basefun_node_ptr_t c) {
  int okay, s, okayA, signA, okayB, signB;
  basefun_node_ptr_t tempNode;
  okay = 0;
  okayA = evaluateSign(&signA, c);
  tempNode = makeConstantDouble(1.0);
  if (okayA)
    okayB = compareConstant(&signB, c, tempNode, NULL, 0);
  else
    okayB = 0;
  if (okayA && okayB) {
    okay = 1;
    if (signA < 0) {
      s = -1;
    } else {
      if (signB < 0) {
        s = 0;
      } else {
        s = 1;
      }
    }
  }
  free_memory(tempNode);
  if (okay) *sign = s;
  return okay;
}

static int nearestint_evalsign(int *sign, basefun_node_ptr_t c) {
  int okay, s, okayA, signA, okayB, signB;
  basefun_node_ptr_t tempNode;
  tempNode = NULL;
  okay = 0;
  okayA = evaluateSign(&signA, c);
  if (okayA) {
    if (signA == 0) {
      okayB = 1;
      signB = 0;
    }
    else if (signA > 0) {
      tempNode = makeConstantDouble(0.5);
      okayB = compareConstant(&signB, c, tempNode, NULL, 0);
    }
    else {
      tempNode = makeConstantDouble(-0.5);
      okayB = compareConstant(&signB, c, tempNode, NULL, 0);
    }
  }
  else okayB = 0;

  if (okayA && okayB) {
    okay = 1;
    if ( (signA < 0) && (signB < 0) ) s = -1;
    else if ( (signA > 0) && (signB > 0) ) s = 1;
    else s = 0;
  }
  if (tempNode != NULL) free_memory(tempNode);
  if (okay) *sign = s;
  return okay;
}




/******************************************************************************/
/*                                                                            */
/*                       POINTWISE EVALUATION FUNCTIONS                       */
/*                                                                            */
/******************************************************************************/
static int sqrt_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_sqrt(res, x, MPFR_RNDN); }
static int exp_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_exp(res, x, MPFR_RNDN); }
static int log_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_log(res, x, MPFR_RNDN); }
static int log2_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_log2(res, x, MPFR_RNDN); }
static int log10_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_log10(res, x, MPFR_RNDN); }
static int sin_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_sin(res, x, MPFR_RNDN); }
static int cos_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_cos(res, x, MPFR_RNDN); }
static int tan_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_tan(res, x, MPFR_RNDN); }
static int asin_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_asin(res, x, MPFR_RNDN); }
static int acos_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_acos(res, x, MPFR_RNDN); }
static int atan_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_atan(res, x, MPFR_RNDN); }
static int sinh_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_sinh(res, x, MPFR_RNDN); }
static int cosh_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_cosh(res, x, MPFR_RNDN); }
static int tanh_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_tanh(res, x, MPFR_RNDN); }
static int asinh_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_asinh(res, x, MPFR_RNDN); }
static int acosh_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_acosh(res, x, MPFR_RNDN); }
static int atanh_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_atanh(res, x, MPFR_RNDN); }
static int abs_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_abs(res, x, MPFR_RNDN); }
static int erf_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_erf(res, x, MPFR_RNDN); }
static int erfc_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_erfc(res, x, MPFR_RNDN); }
static int log1p_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_log1p(res, x, MPFR_RNDN); }
static int expm1_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_expm1(res, x, MPFR_RNDN); }
static int ceil_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_rint_ceil(res, x, MPFR_RNDN); }
static int floor_point_eval(mpfr_ptr res, mpfr_srcptr x) { return mpfr_rint_floor(res, x, MPFR_RNDN); }
static int nearestint_point_eval(mpfr_ptr res, mpfr_srcptr x) { return sollya_mpfr_rint_nearestint(res, x, MPFR_RNDN); }



/******************************************************************************/
/*                                                                            */
/*                          HOW TO TRANSMIT CUTOFF                            */
/*                                                                            */
/******************************************************************************/

/* For the functions that are infinitely bad conditioned at zero (e.g. log) or
   those for which one does not know what to return
   --> reset cutoff to smallest possible value */
static mp_exp_t getRecurseCutoff_default(mp_exp_t cutoff, mp_prec_t prec) {
  UNUSED_PARAM(cutoff);
  UNUSED_PARAM(prec);
  return mpfr_get_emin_min();
}

/* For functions f such that f(x) ~ x at zero and |f(x)|<=|x| on [-1/2, 1/2]
   --> typically f(x) = sum a_i*x^i with a_0=0, a_1=1 and (a_i) has alternating signs */
static mp_exp_t getRecurseCutoff_linearContractive(mp_exp_t cutoff, mp_prec_t prec) {
  UNUSED_PARAM(prec);
  return cutoff;
}

/* For functions f such that f(x) ~ x at zero and |f(x)|<=2|x| on [-1/2, 1/2]
   --> typically f(x) = sum a_i*x^i with a_0=0, a_1=1 and a_i >= 0 decreasing fast enough */
static mp_exp_t getRecurseCutoff_linearNotContractive(mp_exp_t cutoff, mp_prec_t prec) {
  mp_exp_t res;
  UNUSED_PARAM(prec);
  res = cutoff - 1;
  if ((res >= 0) || (res < mpfr_get_emin_min())) res = mpfr_get_emin_min();
  return res;
}

/* For functions that behaves like a+bx around zero with |a/b| not smaller than 2 */
static mp_exp_t getRecurseCutoff_nonZeroWellConditioned(mp_exp_t cutoff, mp_prec_t prec) {
  UNUSED_PARAM(cutoff);
  return -(prec + 3);
}

/* Zero at zero, x < 2^(c * 2) implies sqrt(x) < 2^c */
static mp_exp_t getRecurseCutoff_sqrt(mp_exp_t cutoff, mp_prec_t prec) {
  mp_exp_t res, temp;
  UNUSED_PARAM(prec);
  temp = mpfr_get_emin_min();
  temp >>= 1;
  temp++;
  if (cutoff <= temp) return mpfr_get_emin_min();
  res = cutoff << 2;
  return res;
}

/* No zero at zero, take Taylor development around zero and precision */
static mp_exp_t getRecurseCutoff_cos(mp_exp_t cutoff, mp_prec_t prec) {
  UNUSED_PARAM(cutoff);
  return -((prec >> 1) + 4);
}

/* No zero at zero, take Taylor development around zero and precision */
static mp_exp_t getRecurseCutoff_cosh(mp_exp_t cutoff, mp_prec_t prec) {
  UNUSED_PARAM(cutoff);
  return -(prec + 3); /* TODO: why not the same as cos? */
}

/* Not defined around zero, cutoff = -1 suffices to sees that */
static mp_exp_t getRecurseCutoff_acosh(mp_exp_t cutoff, mp_prec_t prec) {
  UNUSED_PARAM(cutoff);
  UNUSED_PARAM(prec);
  return -1;
}



/******************************************************************************/
/*                                                                            */
/*                  HOW TO TRANSMIT PRECISION ACCROSS FUNCTIONS               */
/*                                                                            */
/******************************************************************************/

/* To use when one does not know what to do */
static mp_prec_t getRecursePrec_default(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  return prec+10;
}

static mp_prec_t getRecursePrec_sqrt(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere */
  return prec;
}
static mp_prec_t getRecursePrec_exp(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Range: double precision range */
  return prec + 12;
}
static mp_prec_t getRecursePrec_log(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around 1 */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_log2(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around 1 */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_log10(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around 1 */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_sin(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around zeros of sin */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_cos(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around zeros of cos */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_tan(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around zeros and poles of tan */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_asin(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere but around +/- 1 */
  return prec + 10;
}
static mp_prec_t getRecursePrec_acos(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  /* Fine everywhere but around + 1 */
  if (considerCutoff) return cutoffPrec;
  return prec + 10;
}
static mp_prec_t getRecursePrec_atan(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere */
  return prec;
}
static mp_prec_t getRecursePrec_sinh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Range: double precision range */
  return prec + 12;
}
static mp_prec_t getRecursePrec_cosh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Range: double precision range */
  return prec + 12;
}
static mp_prec_t getRecursePrec_tanh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* When x is large, very difficult function */
  return prec + 50;
}
static mp_prec_t getRecursePrec_asinh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere */
  return prec;
}
static mp_prec_t getRecursePrec_acosh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere but around + 1 */
  return prec + 10;
}
static mp_prec_t getRecursePrec_atanh(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere but around +/- 1 */
  return prec + 10;
}
static mp_prec_t getRecursePrec_erf(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere */
  return prec;
}
static mp_prec_t getRecursePrec_erfc(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere in negative range, okay up to x = 100 */
  return prec + 15;
}
static mp_prec_t getRecursePrec_log1p(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Fine everywhere but around - 1 */
  return prec + 10;
}
static mp_prec_t getRecursePrec_expm1(mp_prec_t cutoffPrec, mp_prec_t prec, int considerCutoff) {
  UNUSED_PARAM(considerCutoff);
  UNUSED_PARAM(cutoffPrec);
  /* Range: double precision range */
  return prec + 12;
}



/******************************************************************************/
/*                                                                            */
/*                     EXACT RATIONAL EVALUATION ATTEMPT                      */
/*                                                                            */
/******************************************************************************/


static int try_exact_rational_eval_special_cases(baseFunction const *f, mpq_t res, mpq_t a) {
  int success = 1;
  mpfr_t num, den, numres, denres;
  int numtern, dentern;
  mpq_t numresq, denresq;

  mpfr_init2(num, 12); mpz_to_mpfr(num, mpq_numref(a));
  mpfr_init2(numres, mpfr_get_prec(num));

  mpfr_init2(den, 12); mpz_to_mpfr(den, mpq_denref(a));
  mpfr_init2(denres, mpfr_get_prec(den));

  numtern = f->point_eval(numres, num);
  dentern = f->point_eval(denres, den);

  if ((numtern != 0) || (dentern != 0) ||
      (!mpfr_number_p(numres)) ||
      (!mpfr_number_p(denres))) {
    success = 0;
    mpfr_clear(denres); mpfr_clear(numres);
  }
  mpfr_clear(den); mpfr_clear(num);
  if (!success) return 0;

  mpq_init(numresq); sollya_mpfr_to_mpq(numresq, numres);
  mpq_init(denresq); sollya_mpfr_to_mpq(denresq, denres);
  mpfr_clear(denres); mpfr_clear(numres);
  if (f->baseFunctionCode == SQRT) {
    if (mpq_sgn(denresq) == 0)
      success = 0;
    else
      mpq_div(res, numresq, denresq);
  }
  else
    mpq_sub(res, numresq, denresq); /* This corresponds to the case when f->baseFunctionCode == LOG_{2|10} */

  mpq_clear(denresq); mpq_clear(numresq);
  return success;
}


static int try_exact_rational_eval_piecewise_constant_functions(baseFunction const *f, mpq_t res, mpq_t a) {
  sollya_mpfi_t aI, rI;
  mpfr_t aF, rF;
  mp_prec_t prec;
  mp_exp_t E;

  /* Creates a (possibly almost) point interval aI containing a */
  mpfr_init2(aF, 12);
  if (mpqHoldsOnMpfr(aF, a)) {
    sollya_mpfi_init2(aI, mpfr_get_prec(aF));
    sollya_mpfi_set_fr(aI, aF);
  } else {
    sollya_mpfi_init2(aI, 12);
    sollya_mpfi_set_q(aI, a);
    if (sollya_mpfi_has_nan(aI) ||
        sollya_mpfi_has_infinity(aI)) {
      sollya_mpfi_set_prec(aI, ((mp_prec_t) 66000ll));
      sollya_mpfi_set_q(aI, a);
    } else {
      E = sollya_mpfi_max_exp(aI) + ((mp_exp_t) 10);
      if (E < ((mp_exp_t) 0)) {
        E = (mp_exp_t) 0;
      }
      prec = ((mp_prec_t) E) + ((mp_prec_t) 400ll);
      if (prec < ((mp_prec_t) 400)) {
        prec = 400;
      }
      sollya_mpfi_set_prec(aI, prec);
      sollya_mpfi_set_q(aI, a);
    }
  }
  mpfr_clear(aF);

  if (sollya_mpfi_has_nan(aI) || sollya_mpfi_has_infinity(aI)) {
    sollya_mpfi_clear(aI);
    return 0;
  }

  /* Evaluates rI = f(aI) */
  sollya_mpfi_init2(rI, sollya_mpfi_get_prec(aI) + 10);
  f->interval_eval(rI, aI);
  sollya_mpfi_clear(aI);
  if (!sollya_mpfi_is_point_and_real(rI)) {
    sollya_mpfi_clear(rI);
    return 0;
  }

  /* Convert rI to a mpfr and then back to a rational */
  mpfr_init2(rF, sollya_mpfi_get_prec(rI));
  sollya_mpfi_get_left(rF, rI);
  sollya_mpfi_clear(rI);
  if (!mpfr_number_p(rF)) {
    mpfr_clear(rF);
    return 0;
  }
  sollya_mpfr_to_mpq(res, rF);
  mpfr_clear(rF);
  return 1;
}

static int try_exact_rational_eval_f0_eq_0(mpq_t res, mpq_t a) {
  if (mpq_cmp_ui(a, 0, 1) == 0) {
    mpq_set_ui(res, 0, 1);
    return 1;
  }
  return 0;
}

static int try_exact_rational_eval_f0_eq_1(mpq_t res, mpq_t a) {
  if (mpq_cmp_ui(a, 0, 1) == 0) {
    mpq_set_ui(res, 1, 1);
    return 1;
  }
  return 0;
}

static int try_exact_rational_eval_f1_eq_0(mpq_t res, mpq_t a) {
  if (mpq_cmp_ui(a, 1, 1) == 0) {
    mpq_set_ui(res, 0, 1);
    return 1;
  }
  return 0;
}

static int try_exact_rational_eval_abs(mpq_t res, mpq_t x) {
  mpq_abs(res, x); return 1;
}

static int try_exact_rational_eval_sqrt(mpq_t res, mpq_t x) { return try_exact_rational_eval_special_cases(basefun_sqrt, res, x); }
static int try_exact_rational_eval_log2(mpq_t res, mpq_t x) { return try_exact_rational_eval_special_cases(basefun_log2, res, x); }
static int try_exact_rational_eval_log10(mpq_t res, mpq_t x) { return try_exact_rational_eval_special_cases(basefun_log10, res, x); }

static int try_exact_rational_eval_double(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_double, res, x); }
static int try_exact_rational_eval_doubledouble(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_doubledouble, res, x); }
static int try_exact_rational_eval_tripledouble(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_tripledouble, res, x); }
static int try_exact_rational_eval_doubleextended(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_doubleextended, res, x); }
static int try_exact_rational_eval_quad(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_quad, res, x); }
static int try_exact_rational_eval_single(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_single, res, x); }
static int try_exact_rational_eval_halfprecision(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_halfprecision, res, x); }
static int try_exact_rational_eval_floor(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_floor, res, x); }
static int try_exact_rational_eval_ceil(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_ceil, res, x); }
static int try_exact_rational_eval_nearestint(mpq_t res, mpq_t x) { return try_exact_rational_eval_piecewise_constant_functions(basefun_nearestint, res, x); }



/******************************************************************************/
/*                                                                            */
/*                     EXTERNAL ACCESSORS FOR INTERNAL FUNCTIONS              */
/*                                                                            */
/******************************************************************************/

const baseFunction *get_basefun_sqrt() { return basefun_sqrt; }
const baseFunction *get_basefun_exp() { return basefun_exp; }
const baseFunction *get_basefun_log() { return basefun_log; }
const baseFunction *get_basefun_log2() { return basefun_log2; }
const baseFunction *get_basefun_log10() { return basefun_log10; }
const baseFunction *get_basefun_sin() { return basefun_sin; }
const baseFunction *get_basefun_cos() { return basefun_cos; }
const baseFunction *get_basefun_tan() { return basefun_tan; }
const baseFunction *get_basefun_asin() { return basefun_asin; }
const baseFunction *get_basefun_acos() { return basefun_acos; }
const baseFunction *get_basefun_atan() { return basefun_atan; }
const baseFunction *get_basefun_sinh() { return basefun_sinh; }
const baseFunction *get_basefun_cosh() { return basefun_cosh; }
const baseFunction *get_basefun_tanh() { return basefun_tanh; }
const baseFunction *get_basefun_asinh() { return basefun_asinh; }
const baseFunction *get_basefun_acosh() { return basefun_acosh; }
const baseFunction *get_basefun_atanh() { return basefun_atanh; }
const baseFunction *get_basefun_abs() { return basefun_abs; }
const baseFunction *get_basefun_double() { return basefun_double; }
const baseFunction *get_basefun_single() { return basefun_single; }
const baseFunction *get_basefun_halfprecision() { return basefun_halfprecision; }
const baseFunction *get_basefun_quad() { return basefun_quad; }
const baseFunction *get_basefun_doubledouble() { return basefun_doubledouble; }
const baseFunction *get_basefun_tripledouble() { return basefun_tripledouble; }
const baseFunction *get_basefun_erf() { return basefun_erf; }
const baseFunction *get_basefun_erfc() { return basefun_erfc; }
const baseFunction *get_basefun_log1p() { return basefun_log1p; }
const baseFunction *get_basefun_expm1() { return basefun_expm1; }
const baseFunction *get_basefun_doubleextended() { return basefun_doubleextended; }
const baseFunction *get_basefun_ceil() { return basefun_ceil; }
const baseFunction *get_basefun_floor() { return basefun_floor; }
const baseFunction *get_basefun_nearestint() { return basefun_nearestint; }

/******************************************************************************/
/*                                                                            */
/*                     SHORTCUTS FOR CONSTRUCTING FUNCTIONS                   */
/*                                                                            */
/******************************************************************************/

node *makeSqrt(node *op1) {
  return makeUnary(op1,basefun_sqrt);
}

node *makeExp(node *op1) {
  return makeUnary(op1,basefun_exp);
}

node *makeLog(node *op1) {
  return makeUnary(op1,basefun_log);
}

node *makeLog2(node *op1) {
  return makeUnary(op1,basefun_log2);
}

node *makeLog10(node *op1) {
  return makeUnary(op1,basefun_log10);
}

node *makeSin(node *op1) {
  return makeUnary(op1,basefun_sin);
}

node *makeCos(node *op1) {
  return makeUnary(op1,basefun_cos);
}

node *makeTan(node *op1) {
  return makeUnary(op1,basefun_tan);
}

node *makeAsin(node *op1) {
  return makeUnary(op1,basefun_asin);
}

node *makeAcos(node *op1) {
  return makeUnary(op1,basefun_acos);
}

node *makeAtan(node *op1) {
  return makeUnary(op1,basefun_atan);
}

node *makeAbs(node *op1) {
  return makeUnary(op1,basefun_abs);
}

node *makeDouble(node *op1) {
  return makeUnary(op1,basefun_double);
}

node *makeSingle(node *op1) {
  return makeUnary(op1,basefun_single);
}

node *makeQuad(node *op1) {
  return makeUnary(op1,basefun_quad);
}

node *makeHalfPrecision(node *op1) {
  return makeUnary(op1,basefun_halfprecision);
}

node *makeDoubledouble(node *op1) {
  return makeUnary(op1,basefun_doubledouble);
}

node *makeTripledouble(node *op1) {
  return makeUnary(op1,basefun_tripledouble);
}

node *makeErf(node *op1 ) {
  return makeUnary(op1,basefun_erf);
}

node *makeErfc(node *op1) {
  return makeUnary(op1,basefun_erfc);
}

node *makeLog1p(node *op1) {
  return makeUnary(op1,basefun_log1p);
}

node *makeExpm1(node *op1) {
  return makeUnary(op1,basefun_expm1);
}

node *makeDoubleextended(node *op1) {
  return makeUnary(op1,basefun_doubleextended);
}

node *makeCeil(node *op1) {
  return makeUnary(op1,basefun_ceil);
}

node *makeFloor(node *op1) {
  return makeUnary(op1,basefun_floor);
}

node *makeNearestInt(node *op1) {
  return makeUnary(op1,basefun_nearestint);
}

node *makeSinh(node *op1) {
  return makeUnary(op1,basefun_sinh);
}

node *makeCosh(node *op1) {
  return makeUnary(op1,basefun_cosh);
}

node *makeTanh(node *op1) {
  return makeUnary(op1,basefun_tanh);
}

node *makeAsinh(node *op1) {
  return makeUnary(op1,basefun_asinh);
}

node *makeAcosh(node *op1) {
  return makeUnary(op1,basefun_acosh);
}

node *makeAtanh(node *op1) {
  return makeUnary(op1,basefun_atanh);
}
