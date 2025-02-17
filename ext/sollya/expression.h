/*

  Copyright 2006-2018 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  Boîte Courrier 169
  4, place Jussieu
  F-75252 Paris Cedex 05
  France,

  Sorbonne Université
  CNRS, Laboratoire d'Informatique de Paris 6, LIP6
  F - 75005 Paris
  France,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2)

  and by

  Centre de recherche INRIA Sophia Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors Ch. Lauter, S. Chevillard

  christoph.lauter@ens-lyon.org
  sylvain.chevillard@ens-lyon.org

  This software is a computer program whose purpose is to provide an
  environment for safe floating-point code development. It is
  particularly targeted to the automated implementation of
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

#ifndef EXPRESSION_H
#define EXPRESSION_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <mpfr.h>
#include <stdio.h>
#include "base-functions.h"
#include "chain.h"
#include "library.h"
#include "hooks.h"
#include "polynomials.h"

/* A macro to handle compiler warnings for intended fall-thrus in switch statements */

#ifdef HAVE_FALLTHRU_ATTRIBUTE
#      if HAVE_FALLTHRU_ATTRIBUTE
#         define SOLLYA_INTENDED_FALLTHRU    __attribute__((fallthrough));
#      else
#         define SOLLYA_INTENDED_FALLTHRU
#      endif
#else
#      define SOLLYA_INTENDED_FALLTHRU
#endif

/* End of the fall-thru macro */


enum __point_eval_enum_t {
  POINT_EVAL_FAILURE = 0,
  POINT_EVAL_EXACT,
  POINT_EVAL_CORRECTLY_ROUNDED,
  POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT,
  POINT_EVAL_FAITHFULLY_ROUNDED,
  POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT,
  POINT_EVAL_BELOW_CUTOFF
};
typedef enum __point_eval_enum_t point_eval_t;

/* Possible values for nodeType */
#define VARIABLE 0
#define CONSTANT 1
#define ADD 2
#define SUB 3
#define MUL 4
#define DIV 5
#define NEG 6
#define UNARY_BASE_FUNC 7
#define POW 8
#define PI_CONST 9
#define POLYNOMIAL 10
#define LIBRARYFUNCTION 11
#define LIBRARYCONSTANT 12
#define PROCEDUREFUNCTION 13
/* Attention: other values for nodeType defined in execute.h, and also #define MEMREF 278 below */

#define SOLLYA_MAX_ARG_ARRAY_ALLOC_SIZE 33554432ull

#ifndef NODE_TYPEDEF
#define NODE_TYPEDEF
typedef struct nodeStruct node;
#endif

typedef struct memRefCacheStruct *memRefCache;
struct memRefCacheStruct {
  sollya_mpfi_t *evalCacheX, *evalCacheY;
  mp_prec_t evalCachePrec;
  node *derivCache;
  node *derivUnsimplCache;
  int simplifyCacheDoesNotSimplify;
  int simplifyCacheRationalMode;
  node *simplifyCache;
  node *hornerCache;
  int isCorrectlyTyped;
  eval_hook_t *evaluationHook;
  polynomial_t polynomialRepresentation;
  int memRefChildFromPolynomial;
  mpfr_t *pointEvalCacheX, *pointEvalCacheY;
  mp_exp_t pointEvalCacheCutoff;
  point_eval_t pointEvalCacheResultType;
  uint64_t hash;
  int hashComputed;
  int treeSizeCache;
  int treeSizeCacheFilled;
  node *substituteCacheY;
  node *substituteCacheX;
  int substituteCacheMaySimplify;
  int containsNotANumbersIsCached;
  int containsNotANumbersCacheResult;
  int containsHooksIsCached;
  int containsHooksCacheResult;
  unsigned int containsHooksCall;
  int isPureTreeIsCached;
  int isPureTreeCacheResult;
  int isConstantIsCached;
  int isConstantCacheResult;
  int faithEvalOptFunSupportCached;
  int faithEvalOptFunSupportCacheResult;
};

struct nodeStruct
{
  int nodeType;
  mpfr_t *value;
  node *child1;
  node *child2;
  libraryFunction *libFun;
  baseFunction const *baseFun;
  int libFunDeriv;
  char *string;
  chain *arguments;
  libraryProcedure *libProc;
  node **argArray;
  int argArraySize;
  size_t argArrayAllocSize;
  memRefCache cache;
};

/* HELPER TYPE FOR THE PARSER */
typedef struct doubleNodeStruct doubleNode;
struct doubleNodeStruct
{
  node *a;
  node *b;
};

/* END HELPER TYPE */


typedef struct rangetypeStruct rangetype;

struct rangetypeStruct
{
  mpfr_t *a;
  mpfr_t *b;
};

#define MEMREF 278

void *safeMalloc(size_t);
void safeFree(void *);

/* Highly important and often used functions */

static inline node* allocateNode() {
  node *res;
  res = (node *) safeMalloc(sizeof(node));
  return res;
}

static inline void freeNode(node *ptr) {
  safeFree(ptr);
}

static inline void freeMemRef(node *ptr) {
  safeFree(ptr);
}

static inline node* addMemRefEvenOnNull(node *tree) {
  void *ptr, *ptrCache;
  node *res;

  ptr = safeMalloc(sizeof(node) + sizeof(struct memRefCacheStruct));
  ptrCache = ptr + sizeof(node);
  res = (node *) ptr;

  res->nodeType = MEMREF;
  res->child1 = tree;
  res->libFunDeriv = 1;
  res->child2 = NULL;
  res->arguments = NULL;
  res->value = NULL;
  res->cache = (memRefCache) ptrCache;

  res->cache->evalCacheX = NULL;
  res->cache->evalCacheY = NULL;
  res->cache->evalCachePrec = 0;
  res->cache->derivCache = NULL;
  res->cache->derivUnsimplCache = NULL;
  res->cache->simplifyCacheDoesNotSimplify = -1;
  res->cache->simplifyCacheRationalMode = -1;
  res->cache->simplifyCache = NULL;
  res->cache->hornerCache = NULL;
  res->cache->isCorrectlyTyped = 0;
  res->cache->evaluationHook = NULL;
  res->cache->polynomialRepresentation = NULL;
  res->cache->memRefChildFromPolynomial = 0;
  res->cache->pointEvalCacheX = NULL;
  res->cache->pointEvalCacheY = NULL;
  res->cache->hashComputed = 0;
  res->cache->treeSizeCache = -1;
  res->cache->treeSizeCacheFilled = 0;
  res->cache->substituteCacheY = NULL;
  res->cache->substituteCacheX = NULL;
  res->cache->substituteCacheMaySimplify = 0;
  res->cache->containsNotANumbersIsCached = 0;
  res->cache->containsNotANumbersCacheResult = 1;
  res->cache->containsHooksIsCached = 0;
  res->cache->containsHooksCacheResult = 0;
  res->cache->containsHooksCall = 0u;
  res->cache->isPureTreeIsCached = 0;
  res->cache->isPureTreeCacheResult = 0;
  res->cache->isConstantIsCached = 0;
  res->cache->isConstantCacheResult = 0;
  res->cache->faithEvalOptFunSupportCached = 0;
  res->cache->faithEvalOptFunSupportCacheResult = 0;

  return res;
}

static inline node* addMemRef(node *tree) {
  if (tree == NULL) return NULL;
  if (tree->nodeType == MEMREF) return tree;
  return addMemRefEvenOnNull(tree);
}

static inline node* getMemRefChild(node *tree) {
  if (tree->nodeType != MEMREF) return tree->child1;
  if (tree->child1 != NULL) return tree->child1;
  if (tree->cache->polynomialRepresentation == NULL) return NULL;
  tree->child1 = polynomialGetExpressionExplicit(tree->cache->polynomialRepresentation);
  tree->cache->memRefChildFromPolynomial = 1;
  return tree->child1;
}

static inline node* accessThruMemRef(node *tree) {
  node *res;
  if (tree == NULL) return NULL;
  res = tree;
  while (res->nodeType == MEMREF) {
    res = getMemRefChild(res);
  }
  return res;
}

/* End of highly used functions */

void printTree(node *tree);
node* differentiate(node *tree);
node* differentiateUnsimplified(node *tree);
node* simplifyTree(node *tree);
void free_memory(node *tree);
int evaluateConstantExpression(mpfr_t result, node *tree, mp_prec_t prec);
void evaluate(mpfr_t result, node *tree, mpfr_t x, mp_prec_t prec);
void printValue(mpfr_t *value);
node* copyTree(node *tree);
int treeContainsHooks(node *tree);
node* horner(node *tree);
int getDegree(node *tree);
int getDegreeSilent(node *tree);
int getDegreeMpz(mpz_t res, node *tree);
int getDegreeMpzVerified(mpz_t res, node *tree);
int isPolynomialExtraSafe(node *tree);
node* expand(node *tree);
node* simplifyTreeErrorfree(node *tree);
int isNotUniformlyZero(node *tree);
int isNotUniformlyInfinite(node *tree);
int isIntegerConstant(node *tree);
int getNumeratorDenominator(node **numerator, node **denominator, node *tree);
node *substitute(node* tree, node *t);
node *substituteEnhanced(node* tree, node *t, int doNotEvaluate, int maySimplify);
void composePolynomials(node **, chain **, node *, node *, mp_prec_t);
int readDyadic(mpfr_t res, char *c);
int readHexadecimal(mpfr_t res, char *c);
int isPolynomial(node *tree);
int isRationalFunction(node *tree);
int isAffine(node *tree);
int arity(node *tree);
void fprintValue(FILE *fd, mpfr_t value);
void fprintTree(FILE *fd, node *tree);
int isSyntacticallyEqual(node *tree1, node *tree2);
void fprintHeadFunction(FILE *fd,node *tree, char *x, char *y);
int isConstant(node *tree);
void getCoefficients(int *degree, node ***coefficients, node *poly);
node *makePolynomial(mpfr_t *coefficients, int degree);
node *makePolynomialConstantExpressions(node **coeffs, int deg);
int treeSize(node *tree);
void printMpfr(mpfr_t x);
int highestDegreeOfPolynomialSubexpression(node *tree);
node *getIthCoefficient(node *poly, int i);
node *getSubpolynomial(node *poly, chain *monomials, int fillDegrees, mp_prec_t prec);
node *makeCanonical(node *func, mp_prec_t prec);
char *sprintTree(node *tree);
char *sprintValue(mpfr_t *value);
char *sPrintHexadecimal(mpfr_t x);
void printBinary(mpfr_t x);
int isHorner(node *tree);
int isCanonical(node *tree);
int isPowerOfVariable(node *);
char *sprintMidpointMode(mpfr_t a, mpfr_t b);
void fprintValueForXml(FILE *, mpfr_t );
void fprintValueWithPrintMode(FILE *, mpfr_t );
void fprintTreeWithPrintMode(FILE *, node *);
int readDecimalConstant(mpfr_t ,char *);
int getMaxPowerDivider(node *tree);
void simplifyMpfrPrec(mpfr_t rop, mpfr_t op);
node *simplifyRationalErrorfree(node *tree);
int tryEvaluateConstantTermToMpq(mpq_t res, node *tree);
node *simplifyAllButDivision(node *tree);
int sollya_mpfr_to_mpq( mpq_t y, mpfr_t x);
int sollya_mpfr_to_mpz( mpz_t y, mpfr_t x);
mp_prec_t getMpzPrecision(mpz_t x);
void mpz_to_mpfr(mpfr_t res, mpz_t a);
int mpqHoldsOnMpfr(mpfr_t res, mpq_t a);
int containsOnlyRealNumbers(node * tree);
void tryCopyTreeAnnotations(node *newTree, node *oldTree);
int containsNonDifferentiableSubfunctions(node *tree);
node *gcd(node *a, node *b);
node *eucldiv(node *a, node *b);
node *euclmod(node *a, node *b);

node *makeVariable();
node *makeConstant(mpfr_t x);
node *makeConstantDouble(double x);
node *makeConstantInt(int x);
node *makeConstantUnsignedInt(unsigned int x);
node *makeConstantMpz(mpz_t x);
node *makeAdd(node *op1, node *op2);
node *makeSub(node *op1, node *op2);
node *makeMul(node *op1, node *op2);
node *makeDiv(node *op1, node *op2);
node *makeNeg(node *op1);
node *makePow(node *op1, node *op2);
node *makePi();
node *makeUnary(node *op1, baseFunction const *baseFun);
node *makeBinary(node *op1, node *op2, int opType);

node* addMemRef(node *);
node* accessThruMemRef(node *);

void freeVariableCache();
void resetVariableCacheHandling();


#endif /* ifdef EXPRESSION_H*/
