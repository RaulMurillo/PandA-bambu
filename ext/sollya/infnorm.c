/*

  Copyright 2006-2018 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2),

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
  France

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

#include <mpfr.h>
#include "mpfi-compat.h"
#include "expression.h"
#include "infnorm.h"
#include "chain.h"
#include "double.h"
#include "general.h"
#include "proof.h"
#include "remez.h"
#include "execute.h"
#include "hooks.h"
#include <string.h>


#include <stdio.h> /* fprintf, fopen, fclose, */
#include <stdlib.h> /* exit, free, mktemp */
#include <errno.h>

#define DIFFSIZE 5000000

void printInterval(sollya_mpfi_t interval);

/* Evaluate a library constant function into an interval */
void libraryConstantToInterval(sollya_mpfi_t res, node *tree) {
  mpfr_t approx, lbound, rbound;
  mp_prec_t prec = sollya_mpfi_get_prec(res);

  mpfr_init2(approx, prec + 20); /* some guard bits may avoid reinit in tree->libFun */
  if (accessThruMemRef(tree)->libFun->hasData) {
    enterExternalCode();
    ((void (*)(mpfr_t, mp_prec_t, void *)) (accessThruMemRef(tree)->libFun->code))(approx, prec, accessThruMemRef(tree)->libFun->data);
    leaveExternalCode();
  } else {
    enterExternalCode();
    ((void (*)(mpfr_t, mp_prec_t)) (accessThruMemRef(tree)->libFun->code))(approx, prec);
    leaveExternalCode();
  }
  mpfr_init2(lbound, prec-2);
  mpfr_init2(rbound, prec-2);
  mpfr_set(lbound, approx, GMP_RNDD);
  mpfr_set(rbound, approx, GMP_RNDU);
  mpfr_nextbelow(lbound);
  mpfr_nextabove(rbound);

  sollya_mpfi_interv_fr(res, lbound, rbound);
  mpfr_clear(approx);
  mpfr_clear(lbound);
  mpfr_clear(rbound);
  return;
}

void evaluateNewtonMPFRWithStartPoint(mpfr_t result, node *tree, mpfr_t x, mp_prec_t prec, mpfr_t a, mpfr_t b) {
  mp_prec_t p;
  mpfr_t myX;

  if ((mpfr_number_p(x) &&
       mpfr_number_p(a) &&
       mpfr_number_p(b)) &&
      ((mpfr_cmp(a, x) <= 0) &&
       (mpfr_cmp(x, b) <= 0))) {
    evaluate(result, tree, x, prec);
  }

  p = mpfr_get_prec(a);
  if (mpfr_get_prec(b) > p) p = mpfr_get_prec(b);
  if (mpfr_get_prec(x) > p) p = mpfr_get_prec(x);
  mpfr_init2(myX, p);
  if (mpfr_cmp(x, a) < 0) {
    mpfr_set(myX, a, GMP_RNDN);
  } else {
    if (mpfr_cmp(b, x) < 0) {
      mpfr_set(myX, b, GMP_RNDN);
    } else {
      mpfr_set(myX, x, GMP_RNDN);
    }
  }

  evaluate(result, tree, myX, prec);
  mpfr_clear(myX);
}

int newtonMPFRWithStartPoint(mpfr_t res, node *tree, node *diff_tree, mpfr_t a, mpfr_t b, mpfr_t start, mp_prec_t prec) {
  mpfr_t x, x2, temp1, temp2, am, bm, tempX, tempY;
  unsigned long int n=1;
  int okay, lucky, hasZero, i, freeTrees;
  node *myTree, *myDiffTree;

  if (mpfr_sgn(a) * mpfr_sgn(b) <= 0) {
    mpfr_init2(tempX, 12);
    mpfr_init2(tempY, 12);
    mpfr_set_si(tempX, 0, GMP_RNDN);
    evaluateNewtonMPFRWithStartPoint(tempY, tree, tempX, prec, a, b);
    if (mpfr_zero_p(tempY)) {
      mpfr_set(res,tempY,GMP_RNDN);
      mpfr_clear(tempY);
      mpfr_clear(tempX);
      return 1;
    }
    mpfr_clear(tempY);
    mpfr_clear(tempX);
  }

  freeTrees = 0;
  if (accessThruMemRef(tree)->nodeType == DIV) {
    freeTrees = 1;
    myTree = copyTree(accessThruMemRef(tree)->child1);
    myDiffTree = differentiate(myTree);
  } else {
    myTree = tree;
    myDiffTree = diff_tree;
  }

  mpfr_init2(x,prec);
  mpfr_init2(x2,prec);
  mpfr_init2(temp1,prec);
  mpfr_init2(temp2,prec);
  mpfr_init2(am,prec/2);
  mpfr_init2(bm,prec/2);
  mpfr_set(am,a,GMP_RNDN);
  mpfr_nextbelow(am);
  mpfr_nextbelow(am);
  mpfr_set(bm,b,GMP_RNDN);
  mpfr_nextabove(am);
  mpfr_nextabove(bm);

  okay = 0;

  if (mpfr_sgn(a) != mpfr_sgn(b)) {
    mpfr_set_d(x,0.0,GMP_RNDN);
    evaluateNewtonMPFRWithStartPoint(temp1, myTree, x, prec, a, b);
    if (mpfr_zero_p(temp1)) {
      mpfr_set(res,x,GMP_RNDN);
      okay = 1;
    }
  }

  if (!okay) {
    evaluateNewtonMPFRWithStartPoint(temp1, myTree, a, prec, a, b);
    if (mpfr_zero_p(temp1)) {
      mpfr_set(res,a,GMP_RNDN);
      okay = 1;
    } else {
      evaluateNewtonMPFRWithStartPoint(temp2, myTree, b, prec, a, b);
      if (mpfr_zero_p(temp2)) {
	mpfr_set(res,b,GMP_RNDN);
	okay = 1;
      } else {

	mpfr_mul(temp1,temp1,temp2,GMP_RNDN);
	hasZero = (mpfr_sgn(temp1) <= 0);

	mpfr_set(x,start,GMP_RNDN);
	lucky = 0;

	i = 5000;
	while((n<=(unsigned long int) prec+25) && (mpfr_cmp(am,x) <= 0) && (mpfr_cmp(x,bm) <= 0) && (i > 0)) {
	  evaluateNewtonMPFRWithStartPoint(temp1, myTree, x, prec, a, b);
	  if (mpfr_zero_p(temp1)) {
	    lucky = 1;
	    break;
	  }
	  evaluateNewtonMPFRWithStartPoint(temp2, myDiffTree, x, prec, a, b);
	  mpfr_div(temp1, temp1, temp2, GMP_RNDN);
	  mpfr_sub(x2, x, temp1, GMP_RNDN);
	  if (mpfr_cmp(x2,x) == 0) break;
	  if (mpfr_zero_p(x) || mpfr_zero_p(x2)) {
	    n *= 2;
	  } else {
	    if (mpfr_get_exp(x) == mpfr_get_exp(x2)) {
	      n *= 2;
	    } else {
	      i--;
	    }
	  }
	  mpfr_set(x,x2,GMP_RNDN);
	}

	if (mpfr_cmp(x,a) < 0) {
	  mpfr_set(res,a,GMP_RNDN);
	  if (hasZero) {
	    okay = 1;
	  } else {
	    evaluateNewtonMPFRWithStartPoint(temp1, myTree, x, prec, a, b);
	    evaluateNewtonMPFRWithStartPoint(temp2, myDiffTree, x, prec, a, b);
	    mpfr_div(temp1, temp1, temp2, GMP_RNDN);
	    mpfr_sub(x, x, temp1, GMP_RNDN);
	    if (mpfr_cmp(x,a) >= 0) {
	      okay = 1;
	    } else {
	      okay = 0;
	    }
	  }
	} else {
	  if (mpfr_cmp(b,x) < 0) {
	    mpfr_set(res,b,GMP_RNDN);
	    if (hasZero) {
	      okay = 1;
	    } else {
	      evaluateNewtonMPFRWithStartPoint(temp1, myTree, x, prec, a, b);
	      evaluateNewtonMPFRWithStartPoint(temp2, myDiffTree, x, prec, a, b);
	      mpfr_div(temp1, temp1, temp2, GMP_RNDN);
	      mpfr_sub(x, x, temp1, GMP_RNDN);
	      if (mpfr_cmp(b,x) >= 0) {
		okay = 1;
	      } else {
		okay = 0;
	      }
	    }
	  } else {
	    mpfr_set(res,x,GMP_RNDN);
	    if (!lucky) {
	      evaluateNewtonMPFRWithStartPoint(temp1, myTree, x, prec, a, b);
	      evaluateNewtonMPFRWithStartPoint(temp2, myDiffTree, x, prec, a, b);
	      mpfr_div(temp1, temp1, temp2, GMP_RNDN);
	      mpfr_abs(temp1,temp1,GMP_RNDN);
	      mpfr_abs(x,x,GMP_RNDN);
	      mpfr_div_ui(x,x,1,GMP_RNDN);
	      okay = (mpfr_cmp(temp1,x) <= 0);
	    } else {
	      okay = 1;
	    }
	  }
	}
      }
    }
  }

  if (freeTrees) {
    free_memory(myTree);
    free_memory(myDiffTree);
  }

  mpfr_clear(x); mpfr_clear(temp1); mpfr_clear(temp2); mpfr_clear(x2); mpfr_clear(am); mpfr_clear(bm);
  return okay;
}

int newtonMPFR(mpfr_t res, node *tree, node *diff_tree, mpfr_t a, mpfr_t b, mp_prec_t prec) {
  mpfr_t start;
  int result;

  mpfr_init2(start,prec);
  mpfr_add(start,a,b,GMP_RNDN);
  mpfr_div_2ui(start,start,1,GMP_RNDN);

  result = newtonMPFRWithStartPoint(res, tree, diff_tree, a, b, start, prec);

  mpfr_clear(start);

  return result;
}


void makeMpfiAroundMpfr(sollya_mpfi_t res, mpfr_t x, unsigned int thousandUlps) {
  mpfr_t xp, xs;
  mp_prec_t prec;
  sollya_mpfi_t xI;


  prec = mpfr_get_prec(x);
  mpfr_init2(xp,prec);
  mpfr_init2(xs,prec);
  sollya_mpfi_init2(xI,prec);

  mpfr_set(xp,x,GMP_RNDD);
  mpfr_set(xs,x,GMP_RNDU);

  mpfr_nextbelow(xp);
  mpfr_nextabove(xs);

  sollya_mpfi_interv_fr(xI,xp,xs);

  sollya_mpfi_blow(xI,xI,(((double) thousandUlps) * 250.0)); /* THIS CAST IS DANGEROUS ON SYSTEMS WHERE INTS ARE 64 BITS */

  sollya_mpfi_set(res,xI);

  sollya_mpfi_clear(xI);
  mpfr_clear(xp);
  mpfr_clear(xs);
}

static inline sollya_mpfi_t *chooseAndInitMpfiPtr(sollya_mpfi_t *localPtr, mp_prec_t prec) {
  sollya_mpfi_t *ptr;

  ptr = getReusedGlobalMPFIVars(1, prec);
  if (ptr == NULL) {
    sollya_mpfi_init2(*localPtr, prec);
    ptr = localPtr;
  }
  return ptr;
}

static inline void clearChosenMpfiPtr(sollya_mpfi_t *ptr, sollya_mpfi_t *localPtr) {
  if (ptr == localPtr) {
    sollya_mpfi_clear(*localPtr);
    return;
  }
  returnReusedGlobalMPIVars(1);
}

static inline chain* evaluateIRec(sollya_mpfi_t result, node *tree, sollya_mpfi_t x, mp_prec_t prec, int simplifiesA, int simplifiesB, mpfr_t *hopitalPoint, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThinArg, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  sollya_mpfi_t stack1, stack2, tempI, tempI2;
  sollya_mpfi_t stack3, zI, numeratorInZI, denominatorInZI, newExcludeTemp, xMXZ, temp1, temp2, tempA, tempB;
  sollya_mpfi_t *newExclude;
  sollya_mpfi_t leftConstantTerm, rightConstantTerm;
  sollya_mpfi_t leftLinearTerm, rightLinearTerm;
  mpfr_t al, ar, bl, br, xl, xr, z, z2;
  mpfr_t *newHopitalPoint;
  node *derivNumerator, *derivDenominator, *tempNode;
  node *derivLeft, *derivRight;
  chain *leftExcludes, *rightExcludes, *excludes, *t1, *t2;
  chain *leftExcludesConstant, *rightExcludesConstant;
  chain *leftExcludesLinear, *rightExcludesLinear;
  exprBoundTheo *leftTheo, *rightTheo, *internalTheo;
  exprBoundTheo *leftTheoConstant, *rightTheoConstant, *leftTheoLinear, *rightTheoLinear;
  int isPolynom;
  int xIsPoint;
  sollya_mpfi_t tempInterval;
  mp_prec_t *precPtr;
  sollya_mpfi_t *intervalPtr;
  sollya_mpfi_t *reusedVars;
  node *tC1, *tC2, *temp;
  int workForThin;
  mp_prec_t pxppp, ppppp;
  int leftRecursionDone, rightRecursionDone;
  sollya_mpfi_t leftRecursionResult, rightRecursionResult;
  int evaluationLibraryfunctionDone;
  mp_prec_t effectivePrecInterval, ePA, ePB;
  sollya_mpfi_t stack3AltA, stack3AltB;

  leftRecursionDone = 0;
  rightRecursionDone = 0;

  workForThin = workForThinArg;

  if (tree->nodeType == MEMREF) {
    if ((theo != NULL) || (!noExcludes)) return evaluateIRec(result, getMemRefChild(tree), x, prec, simplifiesA, simplifiesB, hopitalPoint, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

    if ((tree->arguments != NULL) &&
	(*((mp_prec_t *) tree->arguments->value) >= prec)) {
      sollya_mpfi_set(result, *((sollya_mpfi_t *) tree->arguments->next->value));
      if (!(sollya_mpfi_has_nan(result) || sollya_mpfi_has_infinity(result))) return NULL;
    }

    if ((tree->cache->evalCacheX != NULL) &&
	(tree->cache->evalCacheY != NULL) &&
	(tree->cache->evalCachePrec >= prec) &&
	(sollya_mpfi_get_prec(*(tree->cache->evalCacheY)) >= sollya_mpfi_get_prec(result)) &&
	(sollya_mpfi_equal_p(*(tree->cache->evalCacheX), x))) {
      sollya_mpfi_set(result, *(tree->cache->evalCacheY));
      if (!(sollya_mpfi_has_nan(result) || sollya_mpfi_has_infinity(result))) return NULL;
    }

    if (evaluateWithEvaluationHook(result, x, prec, 1, tree->cache->evaluationHook)) {
      excludes = NULL;
    } else {
      if ((!noLazyHooks) && (evaluateWithEvaluationHook(result, x, prec, 0, tree->cache->evaluationHook))) {
	excludes = NULL;
	if (lazyHookUsed != NULL) {
	  *lazyHookUsed |= 1;
	  *lazyHookUsed = !!*lazyHookUsed;
	}
      } else {
	if (tree->cache->polynomialRepresentation != NULL) {
	  polynomialEvalMpfi(result, tree->cache->polynomialRepresentation, x);
	  excludes = NULL;
	} else {
	  excludes = evaluateIRec(result, getMemRefChild(tree), x, prec, simplifiesA, simplifiesB, hopitalPoint, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	}
      }
    }

    if ((excludes == NULL) && (!(sollya_mpfi_has_nan(result) || sollya_mpfi_has_infinity(result)))) {
      if (tree->cache->evalCacheX == NULL) {
	tree->cache->evalCacheX = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	tree->cache->evalCacheY = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(tree->cache->evalCacheX),sollya_mpfi_get_prec(x));
	sollya_mpfi_init2(*(tree->cache->evalCacheY),sollya_mpfi_get_prec(result));
	sollya_mpfi_set(*(tree->cache->evalCacheX), x);
	sollya_mpfi_set(*(tree->cache->evalCacheY), result);
	tree->cache->evalCachePrec = prec;
      } else {
	sollya_mpfi_set_prec(*(tree->cache->evalCacheX),sollya_mpfi_get_prec(x));
	sollya_mpfi_set_prec(*(tree->cache->evalCacheY),sollya_mpfi_get_prec(result));
	sollya_mpfi_set(*(tree->cache->evalCacheX), x);
	sollya_mpfi_set(*(tree->cache->evalCacheY), result);
	tree->cache->evalCachePrec = prec;
      }

      if (tree->arguments != NULL) {
	if (prec > *((mp_prec_t *) tree->arguments->value)) {
	  *((mp_prec_t *) tree->arguments->value) = prec;
	  sollya_mpfi_init2(tempInterval, sollya_mpfi_get_prec(*((sollya_mpfi_t *) tree->arguments->next->value)));
	  sollya_mpfi_set(tempInterval, *((sollya_mpfi_t *) tree->arguments->next->value));
	  sollya_mpfi_set_prec(*((sollya_mpfi_t *) tree->arguments->next->value), sollya_mpfi_get_prec(result));
	  sollya_mpfi_intersect(*((sollya_mpfi_t *) tree->arguments->next->value), tempInterval, result);
	  sollya_mpfi_clear(tempInterval);
	}
      } else {
	if (isConstant(tree)) {
	  precPtr = (mp_prec_t *) safeMalloc(sizeof(mp_prec_t));
	  *precPtr = prec;
	  intervalPtr = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*intervalPtr, sollya_mpfi_get_prec(result));
	  sollya_mpfi_set(*intervalPtr, result);
	  tree->arguments = addElement(addElement(NULL, intervalPtr), precPtr);
	}
      }
    }

    return excludes;
  }

  excludes = NULL;

  if (theo != NULL) nullifyExprBoundTheo(theo);

  isPolynom = 0;
  internalTheo = NULL;
  if (theo != NULL) {
    isPolynom = isPolynomial(tree);

    if (isPolynom) {
      theo->functionType = POLYNOMIAL;
      leftTheo = NULL;
      rightTheo = NULL;
    } else {
      internalTheo = theo;
      internalTheo->functionType = tree->nodeType;
      switch (arity(tree)) {
      case 1:
	leftTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	rightTheo = NULL;
	internalTheo->boundLeft = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(internalTheo->boundLeft),prec);
	break;
      case 2:
	leftTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	rightTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	internalTheo->boundLeft = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(internalTheo->boundLeft),prec);
	internalTheo->boundRight = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(internalTheo->boundRight),prec);
	break;
      default:
	leftTheo = NULL;
	rightTheo = NULL;
	break;
      }
      internalTheo->theoLeft = leftTheo;
      internalTheo->theoRight = rightTheo;
    }

    theo->x = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*(theo->x),sollya_mpfi_get_prec(x));
    sollya_mpfi_set(*(theo->x),x);
    theo->function = copyTree(tree);
  } else {
    leftTheo = NULL;
    rightTheo = NULL;
  }

  if ((theo == NULL) &&
      noExcludes) {
    switch (tree->nodeType) {
    case MUL:
    case DIV:
      if ((tree->child1->nodeType == MEMREF) &&
	  (tree->child2->nodeType == MEMREF)) {
	if ((tree->child1->cache->evalCacheX != NULL) &&
	    (tree->child1->cache->evalCacheY != NULL) &&
	    (tree->child2->cache->evalCacheX != NULL) &&
	    (tree->child2->cache->evalCacheY != NULL)) {
	  if (sollya_mpfi_equal_p(*(tree->child1->cache->evalCacheX), x) &&
	      sollya_mpfi_equal_p(*(tree->child2->cache->evalCacheX), x)) {
	    switch (tree->nodeType) {
	    case MUL:
	      if ((sollya_mpfi_is_zero(*(tree->child1->cache->evalCacheY)) ||
		   sollya_mpfi_is_zero(*(tree->child2->cache->evalCacheY))) &&
		  (!(sollya_mpfi_has_infinity(*(tree->child1->cache->evalCacheY)) ||
		     sollya_mpfi_has_infinity(*(tree->child2->cache->evalCacheY)) ||
		     sollya_mpfi_has_nan(*(tree->child1->cache->evalCacheY)) ||
		     sollya_mpfi_has_nan(*(tree->child2->cache->evalCacheY))))) {
		sollya_mpfi_set_si(result, 0);
		return NULL;
	      }
	      break;
	    case DIV:
	      if ((sollya_mpfi_is_zero(*(tree->child1->cache->evalCacheY)) &&
		   (!sollya_mpfi_has_zero(*(tree->child2->cache->evalCacheY)))) &&
		  (!(sollya_mpfi_has_infinity(*(tree->child1->cache->evalCacheY)) ||
		     sollya_mpfi_has_infinity(*(tree->child2->cache->evalCacheY)) ||
		     sollya_mpfi_has_nan(*(tree->child1->cache->evalCacheY)) ||
		     sollya_mpfi_has_nan(*(tree->child2->cache->evalCacheY))))) {
		sollya_mpfi_set_si(result, 0);
		return NULL;
	      }
	      break;
	    }
	  }
	}
      }
    }
  }

  if ((theo == NULL) &&
      noExcludes &&
      (sollya_mpfi_get_prec(result) >= prec)) {
    switch (tree->nodeType) {
    case VARIABLE:
      sollya_mpfi_set(result, x);
      return NULL;
      break;
    case CONSTANT:
      sollya_mpfi_set_fr(result,*(tree->value));
      return NULL;
      break;
    case PI_CONST:
      sollya_mpfi_const_pi(result);
      return NULL;
      break;
    case NEG:
      evaluateIRec(result, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, NULL, 1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      sollya_mpfi_neg(result, result);
      return NULL;
      break;
    case UNARY_BASE_FUNC:
      evaluateIRec(result, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, NULL, 1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      tree->baseFun->interval_eval(result, result);
      return NULL;
      break;
    default:
      break;
    }
  }

  if ((theo == NULL) && noExcludes && (!(sollya_mpfi_has_nan(x) || sollya_mpfi_has_infinity(x)))) {
    switch (tree->nodeType) {
    case VARIABLE:
      sollya_mpfi_set(result, x);
      return NULL;
      break;
    case CONSTANT:
      sollya_mpfi_set_fr(result,*(tree->value));
      return NULL;
      break;
    case PI_CONST:
      sollya_mpfi_const_pi(result);
      return NULL;
      break;
    case NEG:
      reusedVars = getReusedGlobalMPFIVars(1, prec);
      if (reusedVars == NULL) break;
      evaluateIRec(reusedVars[0], tree->child1, x, prec, simplifiesA, simplifiesB, NULL, NULL,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      if (sollya_mpfi_has_nan(reusedVars[0]) ||
	  sollya_mpfi_has_infinity(reusedVars[0])) {
	leftRecursionDone = 1;
	sollya_mpfi_init2(leftRecursionResult, sollya_mpfi_get_prec(reusedVars[0]));
	sollya_mpfi_set(leftRecursionResult, reusedVars[0]);
	returnReusedGlobalMPIVars(1);
	break;
      }
      sollya_mpfi_neg(result, reusedVars[0]);
      returnReusedGlobalMPIVars(1);
      return NULL;
      break;
    case UNARY_BASE_FUNC:
      reusedVars = getReusedGlobalMPFIVars(1, prec);
      if (reusedVars == NULL) break;
      evaluateIRec(reusedVars[0], tree->child1, x, prec, simplifiesA, simplifiesB, NULL, NULL,1, fastAddSub, workForThin,  noLazyHooks, cutoff, lazyHookUsed);
      if (sollya_mpfi_has_nan(reusedVars[0]) ||
	  sollya_mpfi_has_infinity(reusedVars[0])) {
	leftRecursionDone = 1;
	sollya_mpfi_init2(leftRecursionResult, sollya_mpfi_get_prec(reusedVars[0]));
	sollya_mpfi_set(leftRecursionResult, reusedVars[0]);
	returnReusedGlobalMPIVars(1);
	break;
      }
      tree->baseFun->interval_eval(result, reusedVars[0]);
      returnReusedGlobalMPIVars(1);
      return NULL;
      break;
    case ADD:
    case SUB:
      if (!fastAddSub) break;
      /* The fall-through is intended */
      /* fallthrough */
      SOLLYA_INTENDED_FALLTHRU
    case MUL:
    case POW:
      reusedVars = getReusedGlobalMPFIVars(2, prec);
      if (reusedVars != NULL) {
	evaluateIRec(reusedVars[0], tree->child1, x, prec, simplifiesA, simplifiesB, NULL, NULL,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	evaluateIRec(reusedVars[1], tree->child2, x, prec, simplifiesA, simplifiesB, NULL, NULL,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	if (sollya_mpfi_has_nan(reusedVars[0]) ||
	    sollya_mpfi_has_infinity(reusedVars[0]) ||
	    sollya_mpfi_has_nan(reusedVars[1]) ||
	    sollya_mpfi_has_infinity(reusedVars[1])) {
	  leftRecursionDone = 1;
	  sollya_mpfi_init2(leftRecursionResult, sollya_mpfi_get_prec(reusedVars[0]));
	  sollya_mpfi_set(leftRecursionResult, reusedVars[0]);
	  rightRecursionDone = 1;
	  sollya_mpfi_init2(rightRecursionResult, sollya_mpfi_get_prec(reusedVars[1]));
	  sollya_mpfi_set(rightRecursionResult, reusedVars[1]);
	  returnReusedGlobalMPIVars(2);
	} else {
	  switch (tree->nodeType) {
	  case ADD:
	    sollya_mpfi_add(result, reusedVars[0], reusedVars[1]);
	    break;
	  case SUB:
	    sollya_mpfi_sub(result, reusedVars[0], reusedVars[1]);
	    break;
	  case MUL:
	    sollya_mpfi_mul(result, reusedVars[0], reusedVars[1]);
	    break;
	  case POW:
	    sollya_mpfi_pow(result, reusedVars[0], reusedVars[1]);
	    break;
	  }
	  returnReusedGlobalMPIVars(2);
	  return NULL;
	}
      } else {
	sollya_mpfi_init2(stack1, prec);
	sollya_mpfi_init2(stack2, prec);
	evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, NULL,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	evaluateIRec(stack2, tree->child2, x, prec, simplifiesA, simplifiesB, NULL, NULL,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	if (sollya_mpfi_has_nan(stack1) ||
	    sollya_mpfi_has_infinity(stack1) ||
	    sollya_mpfi_has_nan(stack2) ||
	    sollya_mpfi_has_infinity(stack2)) {
	  leftRecursionDone = 1;
	  sollya_mpfi_init2(leftRecursionResult, sollya_mpfi_get_prec(stack1));
	  sollya_mpfi_set(leftRecursionResult, stack1);
	  rightRecursionDone = 1;
	  sollya_mpfi_init2(rightRecursionResult, sollya_mpfi_get_prec(stack2));
	  sollya_mpfi_set(rightRecursionResult, stack2);
	  sollya_mpfi_clear(stack2);
	  sollya_mpfi_clear(stack1);
	} else {
	  switch (tree->nodeType) {
	  case ADD:
	    sollya_mpfi_add(result, stack1, stack2);
	    break;
	  case SUB:
	    sollya_mpfi_sub(result, stack1, stack2);
	    break;
	  case MUL:
	    sollya_mpfi_mul(result, stack1, stack2);
	    break;
	  case POW:
	    sollya_mpfi_pow(result, stack1, stack2);
	    break;
	  }
	  sollya_mpfi_clear(stack2);
	  sollya_mpfi_clear(stack1);
	  return NULL;
	}
      }
      break;
    default:
      break;
    }
  }

  sollya_mpfi_init2(stack1, prec);
  sollya_mpfi_init2(stack2, prec);
  sollya_mpfi_init2(stack3, prec);
  mpfr_init2(al,prec);
  mpfr_init2(ar,prec);
  mpfr_init2(bl,prec);
  mpfr_init2(br,prec);

  sollya_mpfi_diam_abs(al,x);

  if (mpfr_zero_p(al)) xIsPoint = 1; else xIsPoint = 0;

  switch (tree->nodeType) {
  case VARIABLE:
    sollya_mpfi_set(stack3,x);
    excludes = NULL;
    break;
  case CONSTANT:
    sollya_mpfi_set_fr(stack3,*(tree->value));
    excludes = NULL;
    break;
  case ADD:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      leftExcludes = NULL;
    } else {
      leftExcludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    if (noExcludes && rightRecursionDone) {
      sollya_mpfi_set(stack2, rightRecursionResult);
      rightExcludes = NULL;
    } else {
      rightExcludes = evaluateIRec(stack2, tree->child2, x, prec, simplifiesA, simplifiesB, NULL, rightTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    sollya_mpfi_add(stack3, stack1, stack2);
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
      sollya_mpfi_set(*(internalTheo->boundRight),stack2);
    }
    if ((simplifiesA > 0) && (sollya_mpfi_has_zero(stack3)) && (!sollya_mpfi_has_zero(stack1)) && (!sollya_mpfi_has_zero(stack2)) && !xIsPoint) {

      if (internalTheo != NULL) {
	internalTheo->simplificationUsed = DECORRELATE;
	leftTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(leftTheoConstant);
	rightTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(rightTheoConstant);
	leftTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(leftTheoLinear);
	rightTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(rightTheoLinear);
	internalTheo->theoLeftConstant = leftTheoConstant;
	internalTheo->theoRightConstant = rightTheoConstant;
	internalTheo->theoLeftLinear = leftTheoLinear;
	internalTheo->theoRightLinear = rightTheoLinear;
      } else {
	leftTheoConstant = NULL;
	rightTheoConstant = NULL;
	leftTheoLinear = NULL;
	rightTheoLinear = NULL;
      }

      mpfr_init2(z,prec);
      sollya_mpfi_init2(zI,prec);
      sollya_mpfi_init2(leftConstantTerm,prec);
      sollya_mpfi_init2(rightConstantTerm,prec);
      sollya_mpfi_init2(leftLinearTerm,prec);
      sollya_mpfi_init2(rightLinearTerm,prec);
      sollya_mpfi_init2(xMXZ,prec);
      sollya_mpfi_init2(temp1,prec);
      sollya_mpfi_init2(temp2,prec);
      sollya_mpfi_init2(tempA,prec);
      sollya_mpfi_init2(tempB,prec);

      sollya_mpfi_mid(z,x);
      sollya_mpfi_set_fr(zI,z);

      leftExcludesConstant = evaluateIRec(leftConstantTerm, tree->child1, zI, prec, simplifiesA-1, simplifiesB, NULL, leftTheoConstant,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      rightExcludesConstant = evaluateIRec(rightConstantTerm, tree->child2, zI, prec, simplifiesA-1, simplifiesB, NULL, rightTheoConstant,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

      printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_DECORRELATION,"Information: Differentiating while evaluating for decorrelation.\n");

      derivLeft = differentiate(tree->child1);
      derivRight = differentiate(tree->child2);

      leftExcludesLinear = evaluateIRec(leftLinearTerm, derivLeft, x, prec, simplifiesA-1, simplifiesB, NULL, leftTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      rightExcludesLinear = evaluateIRec(rightLinearTerm, derivRight, x, prec, simplifiesA-1, simplifiesB, NULL, rightTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

      sollya_mpfi_add(tempA,leftConstantTerm,rightConstantTerm);
      sollya_mpfi_add(tempB,leftLinearTerm,rightLinearTerm);

      sollya_mpfi_sub(xMXZ,x,zI);

      sollya_mpfi_mul(temp2,xMXZ,tempB);
      sollya_mpfi_add(temp1,tempA,temp2);

      sollya_mpfi_get_left(al,temp1);
      sollya_mpfi_get_right(ar,temp1);

      if (mpfr_number_p(al) && mpfr_number_p(ar)) {

	printMessage(8,SOLLYA_MSG_DECORRELATION_INTERVAL_ADDITION_OR_SUBTRACTION,"Information: decorrelating an interval addition.\n");

	sollya_mpfi_intersect(stack3,stack3,temp1);
	if (!noExcludes) {
	  excludes = concatChains(leftExcludes,rightExcludes);
	  excludes = concatChains(excludes,leftExcludesConstant);
	  excludes = concatChains(excludes,rightExcludesConstant);
	  excludes = concatChains(excludes,leftExcludesLinear);
	  excludes = concatChains(excludes,rightExcludesLinear);
	}

	if (internalTheo != NULL) {
	  internalTheo->boundLeftConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundLeftConstant),prec);
	  internalTheo->boundRightConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundRightConstant),prec);
	  internalTheo->boundLeftLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundLeftLinear),prec);
	  internalTheo->boundRightLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundRightLinear),prec);
	  internalTheo->xZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->xZ),prec);
	  internalTheo->xMXZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->xMXZ),prec);
	  internalTheo->leftDerivative = copyTree(derivLeft);
	  internalTheo->rightDerivative = copyTree(derivRight);
	  sollya_mpfi_set(*(internalTheo->boundLeftConstant),leftConstantTerm);
	  sollya_mpfi_set(*(internalTheo->boundRightConstant),rightConstantTerm);
	  sollya_mpfi_set(*(internalTheo->boundLeftLinear),leftLinearTerm);
	  sollya_mpfi_set(*(internalTheo->boundRightLinear),rightLinearTerm);
	  sollya_mpfi_set(*(internalTheo->xZ),zI);
	  sollya_mpfi_set(*(internalTheo->xMXZ),xMXZ);
	}
      } else {
	if (internalTheo != NULL) {
	  internalTheo->simplificationUsed = 0;
	  freeExprBoundTheo(internalTheo->theoLeftConstant);
	  internalTheo->theoLeftConstant = NULL;
	  freeExprBoundTheo(internalTheo->theoRightConstant);
	  internalTheo->theoRightConstant = NULL;
	  freeExprBoundTheo(internalTheo->theoLeftLinear);
	  internalTheo->theoLeftLinear = NULL;
	  freeExprBoundTheo(internalTheo->theoRightLinear);
	  internalTheo->theoRightLinear = NULL;
	}
	excludes = concatChains(leftExcludes,rightExcludes);
	freeChain(leftExcludesLinear,freeMpfiPtr);
	freeChain(rightExcludesLinear,freeMpfiPtr);
	freeChain(leftExcludesConstant,freeMpfiPtr);
	freeChain(rightExcludesConstant,freeMpfiPtr);
      }

      free_memory(derivLeft);
      free_memory(derivRight);
      mpfr_clear(z);
      sollya_mpfi_clear(zI);
      sollya_mpfi_clear(xMXZ);
      sollya_mpfi_clear(temp1);
      sollya_mpfi_clear(temp2);
      sollya_mpfi_clear(tempA);
      sollya_mpfi_clear(tempB);
      sollya_mpfi_clear(leftConstantTerm);
      sollya_mpfi_clear(rightConstantTerm);
      sollya_mpfi_clear(leftLinearTerm);
      sollya_mpfi_clear(rightLinearTerm);
    } else {
      excludes = concatChains(leftExcludes,rightExcludes);
    }
    break;
  case SUB:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      leftExcludes = NULL;
    } else {
      leftExcludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    if (noExcludes && rightRecursionDone) {
      sollya_mpfi_set(stack2, rightRecursionResult);
      rightExcludes = NULL;
    } else {
      rightExcludes = evaluateIRec(stack2, tree->child2, x, prec, simplifiesA, simplifiesB, NULL, rightTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    sollya_mpfi_sub(stack3, stack1, stack2);
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
      sollya_mpfi_set(*(internalTheo->boundRight),stack2);
    }
    if ((simplifiesA > 0) && (sollya_mpfi_has_zero(stack3)) && (!sollya_mpfi_has_zero(stack1)) && (!sollya_mpfi_has_zero(stack2)) && !xIsPoint) {

      if (internalTheo != NULL) {
	internalTheo->simplificationUsed = DECORRELATE;
	leftTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(leftTheoConstant);
	rightTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(rightTheoConstant);
	leftTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(leftTheoLinear);
	rightTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	nullifyExprBoundTheo(rightTheoLinear);
	internalTheo->theoLeftConstant = leftTheoConstant;
	internalTheo->theoRightConstant = rightTheoConstant;
	internalTheo->theoLeftLinear = leftTheoLinear;
	internalTheo->theoRightLinear = rightTheoLinear;
      } else {
	leftTheoConstant = NULL;
	rightTheoConstant = NULL;
	leftTheoLinear = NULL;
	rightTheoLinear = NULL;
      }

      mpfr_init2(z,prec);
      sollya_mpfi_init2(zI,prec);
      sollya_mpfi_init2(leftConstantTerm,prec);
      sollya_mpfi_init2(rightConstantTerm,prec);
      sollya_mpfi_init2(leftLinearTerm,prec);
      sollya_mpfi_init2(rightLinearTerm,prec);
      sollya_mpfi_init2(xMXZ,prec);
      sollya_mpfi_init2(temp1,prec);
      sollya_mpfi_init2(temp2,prec);
      sollya_mpfi_init2(tempA,prec);
      sollya_mpfi_init2(tempB,prec);


      sollya_mpfi_mid(z,x);
      sollya_mpfi_set_fr(zI,z);

      leftExcludesConstant = evaluateIRec(leftConstantTerm, tree->child1, zI, prec, simplifiesA-1, simplifiesB, NULL, leftTheoConstant,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      rightExcludesConstant = evaluateIRec(rightConstantTerm, tree->child2, zI, prec, simplifiesA-1, simplifiesB, NULL, rightTheoConstant,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

      printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_DECORRELATION,"Information: Differentiating while evaluating for decorrelation.\n");

      derivLeft = differentiate(tree->child1);
      derivRight = differentiate(tree->child2);

      leftExcludesLinear = evaluateIRec(leftLinearTerm, derivLeft, x, prec, simplifiesA-1, simplifiesB, NULL, leftTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      rightExcludesLinear = evaluateIRec(rightLinearTerm, derivRight, x, prec, simplifiesA-1, simplifiesB, NULL, rightTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

      sollya_mpfi_sub(tempA,leftConstantTerm,rightConstantTerm);
      sollya_mpfi_sub(tempB,leftLinearTerm,rightLinearTerm);

      sollya_mpfi_sub(xMXZ,x,zI);

      sollya_mpfi_mul(temp2,xMXZ,tempB);
      sollya_mpfi_add(temp1,tempA,temp2);

      sollya_mpfi_get_left(al,temp1);
      sollya_mpfi_get_right(ar,temp1);

      if (mpfr_number_p(al) && mpfr_number_p(ar)) {

	printMessage(8,SOLLYA_MSG_DECORRELATION_INTERVAL_ADDITION_OR_SUBTRACTION,"Information: decorrelating an interval subtraction.\n");

	sollya_mpfi_intersect(stack3,stack3,temp1);
	if (!noExcludes) {
	  excludes = concatChains(leftExcludes,rightExcludes);
	  excludes = concatChains(excludes,leftExcludesConstant);
	  excludes = concatChains(excludes,rightExcludesConstant);
	  excludes = concatChains(excludes,leftExcludesLinear);
	  excludes = concatChains(excludes,rightExcludesLinear);
	}

	if (internalTheo != NULL) {
	  internalTheo->boundLeftConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundLeftConstant),prec);
	  internalTheo->boundRightConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundRightConstant),prec);
	  internalTheo->boundLeftLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundLeftLinear),prec);
	  internalTheo->boundRightLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->boundRightLinear),prec);
	  internalTheo->xZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->xZ),prec);
	  internalTheo->xMXZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(internalTheo->xMXZ),prec);
	  internalTheo->leftDerivative = copyTree(derivLeft);
	  internalTheo->rightDerivative = copyTree(derivRight);
	  sollya_mpfi_set(*(internalTheo->boundLeftConstant),leftConstantTerm);
	  sollya_mpfi_set(*(internalTheo->boundRightConstant),rightConstantTerm);
	  sollya_mpfi_set(*(internalTheo->boundLeftLinear),leftLinearTerm);
	  sollya_mpfi_set(*(internalTheo->boundRightLinear),rightLinearTerm);
	  sollya_mpfi_set(*(internalTheo->xZ),zI);
	  sollya_mpfi_set(*(internalTheo->xMXZ),xMXZ);
	}
      } else {
	if (internalTheo != NULL) {
	  internalTheo->simplificationUsed = 0;
	  freeExprBoundTheo(internalTheo->theoLeftConstant);
	  internalTheo->theoLeftConstant = NULL;
	  freeExprBoundTheo(internalTheo->theoRightConstant);
	  internalTheo->theoRightConstant = NULL;
	  freeExprBoundTheo(internalTheo->theoLeftLinear);
	  internalTheo->theoLeftLinear = NULL;
	  freeExprBoundTheo(internalTheo->theoRightLinear);
	  internalTheo->theoRightLinear = NULL;
	}
	excludes = concatChains(leftExcludes,rightExcludes);
	freeChain(leftExcludesLinear,freeMpfiPtr);
	freeChain(rightExcludesLinear,freeMpfiPtr);
	freeChain(leftExcludesConstant,freeMpfiPtr);
	freeChain(rightExcludesConstant,freeMpfiPtr);
      }

      free_memory(derivLeft);
      free_memory(derivRight);
      mpfr_clear(z);
      sollya_mpfi_clear(zI);
      sollya_mpfi_clear(xMXZ);
      sollya_mpfi_clear(temp1);
      sollya_mpfi_clear(temp2);
      sollya_mpfi_clear(tempA);
      sollya_mpfi_clear(tempB);
      sollya_mpfi_clear(leftConstantTerm);
      sollya_mpfi_clear(rightConstantTerm);
      sollya_mpfi_clear(leftLinearTerm);
      sollya_mpfi_clear(rightLinearTerm);
    } else {
      excludes = concatChains(leftExcludes,rightExcludes);
    }
    break;
  case MUL:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      leftExcludes = NULL;
    } else {
      leftExcludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    if (noExcludes && rightRecursionDone) {
      sollya_mpfi_set(stack2, rightRecursionResult);
      rightExcludes = NULL;
    } else {
      rightExcludes = evaluateIRec(stack2, tree->child2, x, prec, simplifiesA, simplifiesB, NULL, rightTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    if (xIsPoint &&
	(theo == NULL) &&
	((sollya_mpfi_has_infinity(stack1) && sollya_mpfi_is_zero(stack2)) ||
	 (sollya_mpfi_has_infinity(stack2) && sollya_mpfi_is_zero(stack1)))) {
      /* We have [0] times something that contains infinity */
      if (sollya_mpfi_has_infinity(stack1)) {
	tC1 = tree->child1;
	tC2 = tree->child2;
      } else {
	tC2 = tree->child1;
	tC1 = tree->child2;
      }
      /* Here, the problem is always with tC1 */
      if (accessThruMemRef(tC1)->nodeType == DIV) {
	temp = addMemRef(makeDiv(makeMul(copyTree(tC2), copyTree(accessThruMemRef(tC1)->child1)),
				 copyTree(accessThruMemRef(tC1)->child2)));
	freeChain(leftExcludes,freeMpfiPtr);
	freeChain(rightExcludes,freeMpfiPtr);
	excludes = evaluateIRec(stack3, temp, x, prec, simplifiesA, simplifiesB, NULL, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	free_memory(temp);
      } else {
	/* There's nothing we can do */
	sollya_mpfi_mul(stack3, stack1, stack2);
	if (internalTheo != NULL) {
	  sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
	  sollya_mpfi_set(*(internalTheo->boundRight),stack2);
	}
	excludes = concatChains(leftExcludes,rightExcludes);
      }
    } else {
      sollya_mpfi_mul(stack3, stack1, stack2);
      if (internalTheo != NULL) {
	sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
	sollya_mpfi_set(*(internalTheo->boundRight),stack2);
      }
      excludes = concatChains(leftExcludes,rightExcludes);
    }
    break;
  case DIV:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      leftExcludes = NULL;
    } else {
      leftExcludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    if (noExcludes && rightRecursionDone) {
      sollya_mpfi_set(stack2, rightRecursionResult);
      rightExcludes = NULL;
    } else {
      rightExcludes = evaluateIRec(stack2, tree->child2, x, prec, simplifiesA, simplifiesB, NULL, rightTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }

    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
      sollya_mpfi_set(*(internalTheo->boundRight),stack2);
    }

    sollya_mpfi_get_left(al,stack1);
    sollya_mpfi_get_right(ar,stack1);
    sollya_mpfi_get_left(bl,stack2);
    sollya_mpfi_get_right(br,stack2);
    if (mpfr_zero_p(al) &&
	mpfr_zero_p(ar)) {
      if (mpfr_zero_p(bl) &&
	  mpfr_zero_p(br) &&
	  (simplifiesB > 0)) {
	/* [0;0] / [0;0] */

	printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_HOPITALS_RULE,"Information: Differentiating while evaluating for Hopital's rule.\n");
	derivNumerator = differentiate(tree->child1);
	derivDenominator = differentiate(tree->child2);

	freeChain(leftExcludes,freeMpfiPtr);
	freeChain(rightExcludes,freeMpfiPtr);

	printMessage(8,SOLLYA_MSG_USING_HOPITALS_RULE_ON_POINT_DIVISION,"Information: using Hopital's rule on point division.\n");

	if (internalTheo != NULL) {
	  internalTheo->simplificationUsed = HOPITAL_ON_POINT;
	  internalTheo->leftDerivative = copyTree(derivNumerator);
	  internalTheo->rightDerivative = copyTree(derivDenominator);
	  leftTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	  rightTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	  internalTheo->theoLeftLinear = leftTheoLinear;
	  internalTheo->theoRightLinear = rightTheoLinear;
	} else {
	  leftTheoLinear = NULL;
	  rightTheoLinear = NULL;
	}

	leftExcludes = evaluateIRec(stack1, derivNumerator, x, prec, simplifiesA, simplifiesB-1, NULL, leftTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	rightExcludes = evaluateIRec(stack2, derivDenominator, x, prec, simplifiesA, simplifiesB-1, NULL, rightTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

	free_memory(derivNumerator);
	free_memory(derivDenominator);
	sollya_mpfi_div(stack3, stack1, stack2);
	excludes = concatChains(leftExcludes,rightExcludes);
      } else {
	/* [0;0] / [bl;br], bl,br != 0 */
	freeChain(rightExcludes,freeMpfiPtr);

	printMessage(8,SOLLYA_MSG_SIMPLIFYING_INTERVAL_DIV_WITH_ZERO_POINT_NUMERA, "Information: simplifying an interval division with 0 point numerator.\n");
	printMessage(12,SOLLYA_MSG_CONTINUATION,"Simplification on function\n%b\n",tree);

	if (sollya_mpfi_has_nan(stack2)) {
	  sollya_mpfi_set(stack3, stack2);
	} else {
	  sollya_mpfi_interv_d(stack3,0.0,0.0);
	}
	excludes = leftExcludes;
	if (internalTheo != NULL) {
	  internalTheo->simplificationUsed = NUMERATOR_IS_ZERO;
	  freeExprBoundTheo(internalTheo->theoRight);
	  internalTheo->theoRight = NULL;
	  freeMpfiPtr(internalTheo->boundRight);
	  internalTheo->boundRight = NULL;
	}
      }
    } else {
      if (sollya_mpfi_has_zero(stack2) &&
	  (simplifiesB > 0)) {
	ppppp = prec;
	pxppp = sollya_mpfi_get_prec(x);
	if (pxppp > ppppp) ppppp = pxppp;
	mpfr_init2(xl,ppppp);
	mpfr_init2(xr,ppppp);

	sollya_mpfi_get_left(xl,x);
	sollya_mpfi_get_right(xr,x);

	if ((mpfr_cmp(xl,xr) != 0) || ((mpfr_sgn(xl) == 0) && (mpfr_sgn(xr) == 0))) {

	  printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_HOPITALS_RULE,"Information: Differentiating while evaluating for Hopital's rule.\n");
	  derivDenominator = differentiate(tree->child2);

	  if ((simplifiesB == (hopitalrecursions + 1)) || (hopitalPoint == NULL)){
	    mpfr_init2(z,prec);
	    newtonMPFR(z,tree->child2,derivDenominator,xl,xr,prec);
	    newHopitalPoint = &z;
	  } else {
	    mpfr_init2(z,mpfr_get_prec(*hopitalPoint));
	    mpfr_set(z,*hopitalPoint,GMP_RNDN);
	    newHopitalPoint = hopitalPoint;
	  }

	  if (mpfr_number_p(z)) {
	    sollya_mpfi_init2(zI,prec);
	    sollya_mpfi_set_fr(zI,z);
	    sollya_mpfi_init2(numeratorInZI,prec);
	    sollya_mpfi_init2(denominatorInZI,prec);

	    if (internalTheo != NULL) {
              leftTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
              rightTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
	    } else {
	      leftTheoConstant = NULL;
	      rightTheoConstant = NULL;
	    }

	    t1 = evaluateIRec(numeratorInZI, tree->child1, zI, prec, simplifiesA, simplifiesB-1, newHopitalPoint, leftTheoConstant,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
	    t2 = evaluateIRec(denominatorInZI, tree->child2, zI, prec, simplifiesA, simplifiesB-1, newHopitalPoint, rightTheoConstant,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

	    freeChain(t1,freeMpfiPtr);
	    freeChain(t2,freeMpfiPtr);

	    sollya_mpfi_get_left(al,numeratorInZI);
	    sollya_mpfi_get_right(ar,numeratorInZI);
	    sollya_mpfi_get_left(bl,denominatorInZI);
	    sollya_mpfi_get_right(br,denominatorInZI);

	    if (mpfr_zero_p(al) && mpfr_zero_p(ar) && mpfr_zero_p(bl) && mpfr_zero_p(br)) {
	      /* Hopital's rule can be applied */

	      printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_HOPITALS_RULE,"Information: Differentiating while evaluating for Hopital's rule.\n");
	      derivNumerator = differentiate(tree->child1);

	      tempNode = allocateNode();
	      tempNode->nodeType = DIV;
	      tempNode->child1 = derivNumerator;
	      tempNode->child2 = copyTree(derivDenominator);

	      freeChain(leftExcludes,freeMpfiPtr);
	      freeChain(rightExcludes,freeMpfiPtr);

	      if (internalTheo != NULL) {
		leftTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
		nullifyExprBoundTheo(leftTheoLinear);
		internalTheo->simplificationUsed = HOPITAL;
		internalTheo->leftDerivative = copyTree(tempNode->child1);
		internalTheo->rightDerivative = copyTree(tempNode->child2);
		internalTheo->theoLeftConstant = leftTheoConstant;
		internalTheo->theoRightConstant = rightTheoConstant;
		internalTheo->theoLeftLinear = leftTheoLinear;
		internalTheo->boundLeftLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		sollya_mpfi_init2(*(internalTheo->boundLeftLinear),prec);
		internalTheo->boundLeftConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		sollya_mpfi_init2(*(internalTheo->boundLeftConstant),prec);
		internalTheo->boundRightConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		sollya_mpfi_init2(*(internalTheo->boundRightConstant),prec);
		sollya_mpfi_set(*(internalTheo->boundLeftConstant),numeratorInZI);
		sollya_mpfi_set(*(internalTheo->boundRightConstant),denominatorInZI);
		internalTheo->xZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		sollya_mpfi_init2(*(internalTheo->xZ),prec);
		sollya_mpfi_set(*(internalTheo->xZ),zI);
	      } else {
		leftTheoLinear = NULL;
	      }

	      if (simplifiesB == (hopitalrecursions + 1)) {
		printMessage(8,SOLLYA_MSG_USING_HOPITALS_RULE_IN_GENERAL_CASE,"Information: using Hopital's rule (general case) on denominator zero.\n");
	      }

	      if (simplifiesB != (hopitalrecursions + 1)) {
		printMessage(8,SOLLYA_MSG_RECURSION_ON_USE_OF_HOPITALS_RULE,"Information: recursion on use of Hopital's rule\n");
	      }

	      excludes = evaluateIRec(stack3, tempNode, x, prec, simplifiesA, simplifiesB-1, newHopitalPoint, leftTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

	      if (internalTheo != NULL) sollya_mpfi_set(*(internalTheo->boundLeftLinear),stack3);

	      free_memory(tempNode);
	    } else {
	      if (internalTheo != NULL) {
		freeExprBoundTheo(leftTheoConstant);
		freeExprBoundTheo(rightTheoConstant);
	      }

	      printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_HOPITALS_RULE,"Information: Differentiating while evaluating for Hopital's rule.\n");
	      derivNumerator = differentiate(tree->child1);

	      if ((simplifiesB == (hopitalrecursions + 1)) || (hopitalPoint == NULL)) {
		mpfr_init2(z2,prec);
		newtonMPFR(z2,tree->child1,derivNumerator,xl,xr,prec);
		newHopitalPoint = &z2;
	      } else {
		mpfr_init2(z2,mpfr_get_prec(*hopitalPoint));
		mpfr_set(z2,*hopitalPoint,GMP_RNDN);
		newHopitalPoint = hopitalPoint;
	      }

	      if (mpfr_number_p(z2)) {
		sollya_mpfi_set_fr(zI,z2);

		if (internalTheo != NULL) {
		  leftTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
		  rightTheoConstant = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
		} else {
		  leftTheoConstant = NULL;
		  rightTheoConstant = NULL;
		}

		t1 = evaluateIRec(numeratorInZI, tree->child1, zI, prec, simplifiesA, simplifiesB-1, newHopitalPoint, leftTheoConstant,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
		t2 = evaluateIRec(denominatorInZI, tree->child2, zI, prec, simplifiesA, simplifiesB-1, newHopitalPoint, rightTheoConstant,1, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

		freeChain(t1,freeMpfiPtr);
		freeChain(t2,freeMpfiPtr);

		sollya_mpfi_get_left(al,numeratorInZI);
		sollya_mpfi_get_right(ar,numeratorInZI);
		sollya_mpfi_get_left(bl,denominatorInZI);
		sollya_mpfi_get_right(br,denominatorInZI);

		if (mpfr_zero_p(al) && mpfr_zero_p(ar) && mpfr_zero_p(bl) && mpfr_zero_p(br)) {

		  tempNode = allocateNode();
		  tempNode->nodeType = DIV;
		  tempNode->child1 = derivNumerator;
		  tempNode->child2 = copyTree(derivDenominator);

		  freeChain(leftExcludes,freeMpfiPtr);
		  freeChain(rightExcludes,freeMpfiPtr);

		  if (internalTheo != NULL) {
		    leftTheoLinear = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
		    nullifyExprBoundTheo(leftTheoLinear);
		    internalTheo->simplificationUsed = HOPITAL;
		    internalTheo->leftDerivative = copyTree(tempNode->child1);
		    internalTheo->rightDerivative = copyTree(tempNode->child2);
		    internalTheo->theoLeftConstant = leftTheoConstant;
		    internalTheo->theoRightConstant = rightTheoConstant;
		    internalTheo->theoLeftLinear = leftTheoLinear;
		    internalTheo->boundLeftLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		    sollya_mpfi_init2(*(internalTheo->boundLeftLinear),prec);
		    internalTheo->boundLeftConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		    sollya_mpfi_init2(*(internalTheo->boundLeftConstant),prec);
		    internalTheo->boundRightConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		    sollya_mpfi_init2(*(internalTheo->boundRightConstant),prec);
		    sollya_mpfi_set(*(internalTheo->boundLeftConstant),numeratorInZI);
		    sollya_mpfi_set(*(internalTheo->boundRightConstant),denominatorInZI);
		    internalTheo->xZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		    sollya_mpfi_init2(*(internalTheo->xZ),prec);
		    sollya_mpfi_set(*(internalTheo->xZ),zI);
		  } else {
		    leftTheoLinear = NULL;
		  }

		  if (simplifiesB == (hopitalrecursions + 1)) {
		    printMessage(8,SOLLYA_MSG_USING_HOPITALS_RULE_IN_GENERAL_CASE,"Information: using Hopital's rule (general case) on numerator zero.\n");
		  }

		  if (simplifiesB != (hopitalrecursions + 1)) {
		    printMessage(8,SOLLYA_MSG_RECURSION_ON_USE_OF_HOPITALS_RULE,"Information: recursion on use of Hopital's rule\n");
		  }

		  excludes = evaluateIRec(stack3, tempNode, x, prec, simplifiesA, simplifiesB-1, newHopitalPoint, leftTheoLinear,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

		  if (internalTheo != NULL) sollya_mpfi_set(*(internalTheo->boundLeftLinear),stack3);

		  free_memory(tempNode);
		} else {

		  if (internalTheo != NULL) {
		    freeExprBoundTheo(leftTheoConstant);
		    freeExprBoundTheo(rightTheoConstant);
		  }

		  if (!noExcludes) {
		    newExclude = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		    sollya_mpfi_init2(*newExclude,prec);
		    makeMpfiAroundMpfr(*newExclude,z,16777216);
		    sollya_mpfi_init2(newExcludeTemp,prec);
		    makeMpfiAroundMpfr(newExcludeTemp,z2,16777216);
		    sollya_mpfi_union(*newExclude,*newExclude,newExcludeTemp);
		    sollya_mpfi_clear(newExcludeTemp);

		    excludes = concatChains(leftExcludes,rightExcludes);
		    excludes = addElement(excludes,newExclude);
		  }

		  free_memory(derivNumerator);
		  mpfr_clear(z2);
		  sollya_mpfi_div(stack3, stack1, stack2);
		}
	      } else {

		free_memory(derivNumerator);
		mpfr_clear(z2);

		if (!noExcludes) {
		  newExclude = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
		  sollya_mpfi_init2(*newExclude,prec);
		  makeMpfiAroundMpfr(*newExclude,z,16777216);
		  excludes = concatChains(leftExcludes,rightExcludes);
		  excludes = addElement(excludes,newExclude);
		}

		sollya_mpfi_div(stack3, stack1, stack2);
	      }
	    }

	    sollya_mpfi_clear(numeratorInZI);
	    sollya_mpfi_clear(denominatorInZI);
	    sollya_mpfi_clear(zI);
	  } else {
	    sollya_mpfi_div(stack3, stack1, stack2);
	    excludes = concatChains(leftExcludes,rightExcludes);
	  }
	  free_memory(derivDenominator);
	  mpfr_clear(z);
	} else {
	  sollya_mpfi_div(stack3, stack1, stack2);
	  excludes = concatChains(leftExcludes,rightExcludes);
	}
	mpfr_clear(xl);
	mpfr_clear(xr);
      } else {
	sollya_mpfi_div(stack3, stack1, stack2);
	excludes = concatChains(leftExcludes,rightExcludes);
      }
    }
    break;
  case NEG:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      excludes = NULL;
    } else {
      excludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    sollya_mpfi_neg(stack3, stack1);
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
    }
    break;
  case UNARY_BASE_FUNC:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      excludes = NULL;
    } else {
      excludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    tree->baseFun->interval_eval(stack3, stack1);
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
    }
    break;
  case POW:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      leftExcludes = NULL;
    } else {
      leftExcludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    if (noExcludes && rightRecursionDone) {
      sollya_mpfi_set(stack2, rightRecursionResult);
      rightExcludes = NULL;
    } else {
      rightExcludes = evaluateIRec(stack2, tree->child2, x, prec, simplifiesA, simplifiesB, NULL, rightTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    sollya_mpfi_pow(stack3, stack1, stack2);
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
      sollya_mpfi_set(*(internalTheo->boundRight),stack2);
    }
    excludes = concatChains(leftExcludes,rightExcludes);
    break;
  case LIBRARYFUNCTION:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      excludes = NULL;
    } else {
      excludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    evaluationLibraryfunctionDone = 0;
    if (sollya_mpfi_is_point_and_real(x) &&
	(!sollya_mpfi_is_point_and_real(stack1)) &&
	(!sollya_mpfi_has_zero(stack1)) &&
	(!sollya_mpfi_has_infinity(stack1)) &&
	(!sollya_mpfi_has_nan(stack1)) &&
	(!sollya_mpfi_is_empty(stack1))) {
      if (!sollya_mpfi_enclosure_accurate_enough(stack1, (sollya_mpfi_get_prec(stack1) > 22 ? (sollya_mpfi_get_prec(stack1) - 10) : 12))) {
	if (sollya_mpfi_get_effective_precision(&effectivePrecInterval, stack1)) {
	  if (effectivePrecInterval > 2) {
	    if (effectivePrecInterval < 12) effectivePrecInterval = 12;
	    if ((effectivePrecInterval + 20) < sollya_mpfi_get_prec(stack3)) {
	      sollya_mpfi_init2(stack3AltA, effectivePrecInterval + 10);
	      sollya_mpfi_init2(stack3AltB, effectivePrecInterval + 20);
	      mpfi_init2(tempI, sollya_mpfi_get_prec(stack3AltA));
	      if (tree->libFun->hasData) {
		enterExternalCode();
		((int (*)(mpfi_t, mpfi_t, int, void *)) (tree->libFun->code))(tempI, stack1, tree->libFunDeriv, tree->libFun->data);
		leaveExternalCode();
	      } else {
		enterExternalCode();
		((int (*)(mpfi_t, mpfi_t, int)) (tree->libFun->code))(tempI, stack1, tree->libFunDeriv);
		leaveExternalCode();
	      }
	      sollya_init_and_convert_interval(tempI2, tempI);
	      sollya_mpfi_set(stack3AltA, tempI2);
	      sollya_mpfi_clear(tempI2);
	      mpfi_clear(tempI);
	      mpfi_init2(tempI, sollya_mpfi_get_prec(stack3AltB));
	      if (tree->libFun->hasData) {
		enterExternalCode();
		((int (*)(mpfi_t, mpfi_t, int, void *)) (tree->libFun->code))(tempI, stack1, tree->libFunDeriv, tree->libFun->data);
		leaveExternalCode();
	      } else {
		enterExternalCode();
		((int (*)(mpfi_t, mpfi_t, int)) (tree->libFun->code))(tempI, stack1, tree->libFunDeriv);
		leaveExternalCode();
	      }
	      sollya_init_and_convert_interval(tempI2, tempI);
	      sollya_mpfi_set(stack3AltB, tempI2);
	      sollya_mpfi_clear(tempI2);
	      mpfi_clear(tempI);
	      if (!(sollya_mpfi_has_nan(stack3AltA) ||
		    sollya_mpfi_has_infinity(stack3AltA) ||
		    sollya_mpfi_has_nan(stack3AltB) ||
		    sollya_mpfi_has_infinity(stack3AltB))) {
		if (sollya_mpfi_get_effective_precision(&ePA, stack3AltA)) {
		  if (sollya_mpfi_get_effective_precision(&ePB, stack3AltB)) {
		    if (((ePA > ePB ? ePA : ePB) - (ePA <= ePB ? ePA : ePB)) < 25) {
		      sollya_mpfi_set(stack3, stack3AltB);
		      evaluationLibraryfunctionDone = 1;
		    }
		  }
		}
	      }
	      sollya_mpfi_clear(stack3AltB);
	      sollya_mpfi_clear(stack3AltA);
	    }
	  }
	}
      }
    }
    if (!evaluationLibraryfunctionDone) {
      mpfi_init2(tempI, sollya_mpfi_get_prec(stack3));
      if (tree->libFun->hasData) {
	enterExternalCode();
	((int (*)(mpfi_t, mpfi_t, int, void *)) (tree->libFun->code))(tempI, stack1, tree->libFunDeriv, tree->libFun->data);
	leaveExternalCode();
      } else {
	enterExternalCode();
	((int (*)(mpfi_t, mpfi_t, int)) (tree->libFun->code))(tempI, stack1, tree->libFunDeriv);
	leaveExternalCode();
      }
      sollya_init_and_convert_interval(tempI2, tempI);
      sollya_mpfi_set(stack3, tempI2);
      sollya_mpfi_clear(tempI2);
      mpfi_clear(tempI);
    }
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
    }
    break;
  case PROCEDUREFUNCTION:
    if (noExcludes && leftRecursionDone) {
      sollya_mpfi_set(stack1, leftRecursionResult);
      excludes = NULL;
    } else {
      excludes = evaluateIRec(stack1, tree->child1, x, prec, simplifiesA, simplifiesB, NULL, leftTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
    computeFunctionWithProcedure(stack3, tree->child2, stack1, (unsigned int) tree->libFunDeriv, tree);
    if (internalTheo != NULL) {
      sollya_mpfi_set(*(internalTheo->boundLeft),stack1);
    }
    break;

  case PI_CONST:
    sollya_mpfi_const_pi(stack3);
    excludes = NULL;
    break;
  case LIBRARYCONSTANT:
    libraryConstantToInterval(stack3, tree);
    excludes = NULL;
    break;
  default:
    sollyaFprintf(stderr,"Error: evaluateIRec: unknown identifier in the tree\n");
    exit(1);
  }

  sollya_mpfi_set(result,stack3);

  if (theo != NULL) {
    theo->y = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*(theo->y),sollya_mpfi_get_prec(result));
    sollya_mpfi_set(*(theo->y),result);
  }
  sollya_mpfi_clear(stack1);
  sollya_mpfi_clear(stack2);
  sollya_mpfi_clear(stack3);
  mpfr_clear(al);
  mpfr_clear(ar);
  mpfr_clear(bl);
  mpfr_clear(br);
  if (leftRecursionDone) {
    sollya_mpfi_clear(leftRecursionResult);
  }
  if (rightRecursionDone) {
    sollya_mpfi_clear(rightRecursionResult);
  }

  return excludes;
}

static inline chain* evaluateI(sollya_mpfi_t result, node *tree, sollya_mpfi_t x, mp_prec_t prec, int simplifiesA, int simplifiesB, mpfr_t *hopitalPoint, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThinArg, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  chain *res;

  res = evaluateIRec(result, tree, x, prec, simplifiesA, simplifiesB, hopitalPoint, theo, noExcludes, fastAddSub, workForThinArg, noLazyHooks, cutoff, lazyHookUsed);

  return res;
}

static inline chain* evaluateITaylor(sollya_mpfi_t result, node *func, node *deriv, sollya_mpfi_t x, mp_prec_t prec, int recurse, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThin, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed);

static inline chain* evaluateITaylorOnDiv(sollya_mpfi_t result, node *func, sollya_mpfi_t x, mp_prec_t prec, int recurse, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThin, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  node *numerator, *denominator, *derivNumerator, *derivDenominator;
  chain *excludes, *numeratorExcludes, *denominatorExcludes;
  exprBoundTheo *numeratorTheo, *denominatorTheo;
  sollya_mpfi_t resultNumerator, resultDenominator, resultIndirect;
  mpfr_t tempNaN;

  mpfr_init2(tempNaN, prec);

  if (accessThruMemRef(func)->nodeType == DIV) {
    numerator = accessThruMemRef(func)->child1;
    denominator = accessThruMemRef(func)->child2;
    derivNumerator = differentiate(numerator);
    derivDenominator = differentiate(denominator);
    sollya_mpfi_init2(resultNumerator, prec);
    sollya_mpfi_init2(resultDenominator, prec);
    sollya_mpfi_init2(resultIndirect, prec);

    if (theo != NULL) {
      numeratorTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
      denominatorTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
    } else {
      numeratorTheo = NULL;
      denominatorTheo = NULL;
    }

    numeratorExcludes = evaluateITaylor(resultNumerator, numerator, derivNumerator, x, prec, recurse, numeratorTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    denominatorExcludes = evaluateITaylor(resultDenominator, denominator, derivDenominator, x, prec, recurse, denominatorTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    excludes = concatChains(numeratorExcludes,denominatorExcludes);
    sollya_mpfi_div(resultIndirect, resultNumerator, resultDenominator);
    if (sollya_mpfi_bounded_p(resultIndirect)) {
      sollya_mpfi_set(result, resultIndirect);
      sollya_mpfi_nan_normalize(result);

      if (theo != NULL) {
	theo->functionType = func->nodeType;
	theo->boundLeft = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->boundLeft),prec);
	theo->boundRight = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->boundRight),prec);
	theo->theoLeft = numeratorTheo;
	theo->theoRight = denominatorTheo;
	theo->x = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->x),sollya_mpfi_get_prec(x));
	sollya_mpfi_set(*(theo->x),x);
	theo->function = copyTree(func);
	theo->y = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->y),sollya_mpfi_get_prec(result));
	sollya_mpfi_set(*(theo->y),result);
      }
    } else {
      freeChain(excludes,freeMpfiPtr);
      if (theo != NULL) {
	freeExprBoundTheo(numeratorTheo);
	freeExprBoundTheo(denominatorTheo);
      }
      excludes = evaluateI(result, func, x, prec, 0, hopitalrecursions+1, NULL, theo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
      sollya_mpfi_nan_normalize(result);
    }

    sollya_mpfi_clear(resultNumerator);
    sollya_mpfi_clear(resultDenominator);
    sollya_mpfi_clear(resultIndirect);
    free_memory(derivNumerator);
    free_memory(derivDenominator);
    mpfr_clear(tempNaN);
    return excludes;
  }
  else {
    excludes = evaluateI(result, func, x, prec, 0, hopitalrecursions+1, NULL, theo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    sollya_mpfi_nan_normalize(result);
    mpfr_clear(tempNaN);
    return excludes;
  }
}

static inline chain* evaluateITaylorInner(sollya_mpfi_t, node *, node *, sollya_mpfi_t, mp_prec_t, int, exprBoundTheo *, int, int, int, int, mp_exp_t *, int *);

static inline chain* evaluateITaylorWrapped(sollya_mpfi_t result, node *func, node *deriv, sollya_mpfi_t x, mp_prec_t prec, int recurse, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThin, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  chain *excludes;
  mp_prec_t *precPtr;
  sollya_mpfi_t *intervalPtr;
  sollya_mpfi_t tempInterval;

  if ((theo != NULL) || (!noExcludes)) return evaluateITaylorInner(result, func, deriv, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

  if ((func->nodeType == MEMREF) &&
      (func->arguments != NULL) &&
      (*((mp_prec_t *) func->arguments->value) >= prec)) {
    sollya_mpfi_set(result, *((sollya_mpfi_t *) func->arguments->next->value));
    if (!(sollya_mpfi_has_nan(result) || sollya_mpfi_has_infinity(result))) return NULL;
  }

  if ((func->nodeType == MEMREF) &&
      (func->cache->evaluationHook != NULL) &&
      evaluateWithEvaluationHook(result, x, prec, 1, func->cache->evaluationHook)) {
    excludes = NULL;
  } else {
    if ((func->nodeType == MEMREF) &&
	(func->cache->evaluationHook != NULL) &&
	(!noLazyHooks) &&
	evaluateWithEvaluationHook(result, x, prec, 0, func->cache->evaluationHook)) {
      excludes = NULL;
      if (lazyHookUsed != NULL) {
	*lazyHookUsed |= 1;
	*lazyHookUsed = !!*lazyHookUsed;
      }
    } else {
      excludes = evaluateITaylorInner(result, func, deriv, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
  }

  if ((excludes == NULL) && (func->nodeType == MEMREF) && (!(sollya_mpfi_has_nan(result) || sollya_mpfi_has_infinity(result)))) {
    if (func->arguments != NULL) {
      if (prec > *((mp_prec_t *) func->arguments->value)) {
	*((mp_prec_t *) func->arguments->value) = prec;
	sollya_mpfi_init2(tempInterval, sollya_mpfi_get_prec(*((sollya_mpfi_t *) func->arguments->next->value)));
	sollya_mpfi_set(tempInterval, *((sollya_mpfi_t *) func->arguments->next->value));
	sollya_mpfi_set_prec(*((sollya_mpfi_t *) func->arguments->next->value), sollya_mpfi_get_prec(result));
	sollya_mpfi_intersect(*((sollya_mpfi_t *) func->arguments->next->value), tempInterval, result);
	sollya_mpfi_clear(tempInterval);
      }
    } else {
      if (isConstant(func)) {
	precPtr = (mp_prec_t *) safeMalloc(sizeof(mp_prec_t));
	*precPtr = prec;
	intervalPtr = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*intervalPtr, sollya_mpfi_get_prec(result));
	sollya_mpfi_set(*intervalPtr, result);
	func->arguments = addElement(addElement(NULL, intervalPtr), precPtr);
      }
    }
  }

  return excludes;
}

static inline chain* evaluateITaylor(sollya_mpfi_t result, node *func, node *deriv, sollya_mpfi_t x, mp_prec_t prec, int recurse, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThin, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  if (deriv == NULL) {
    return evaluateITaylorWrapped(result, func, deriv, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
  }
  if (containsNonDifferentiableSubfunctions(func)) {
    return evaluateITaylorWrapped(result, func, NULL, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
  }
  return evaluateITaylorWrapped(result, func, deriv, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
}

static inline chain* evaluateITaylorStart(sollya_mpfi_t result, node *func, node *deriv, sollya_mpfi_t x, mp_prec_t prec, int recurse, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThin, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  if (deriv == NULL) {
    return evaluateITaylor(result, func, deriv, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
  }
  if (accessThruMemRef(func)->nodeType == DIV) {
    if (!sollya_mpfi_is_point_and_real(x)) {
      return evaluateITaylorOnDiv(result, func, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }
  }
  return evaluateITaylor(result, func, deriv, x, prec, recurse, theo, noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
}

static inline chain* evaluateITaylorInner(sollya_mpfi_t result, node *func, node *deriv, sollya_mpfi_t x, mp_prec_t prec, int recurse, exprBoundTheo *theo, int noExcludes, int fastAddSub, int workForThin, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  mpfr_t xZ, rTl, rTr, leftX, rightX;
  sollya_mpfi_t xZI, xZI2, constantTerm, linearTerm, resultTaylor, resultDirect, temp, temp2;
  chain *excludes, *directExcludes, *taylorExcludes, *taylorExcludesLinear, *taylorExcludesConstant;
  exprBoundTheo *constantTheo, *linearTheo, *directTheo;
  node *nextderiv;
  int size;

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ((deriv == NULL) || (mpfr_cmp(&(x->left),&(x->right)) == 0)) {
    if (deriv != NULL)
      printMessage(25,SOLLYA_MSG_AVOIDING_TAYLOR_EVALUATION_ON_POINT_INTERVAL,"Information: avoiding using Taylor's formula on a point interval.\n");
    else
      printMessage(25,SOLLYA_MSG_NO_TAYLOR_EVALUATION_AS_NO_DERIVATIVE_GIVEN,"Warning: no Taylor evaluation is possible because no derivative has been given.\n");

    excludes = evaluateI(result, func, x, prec, 1, hopitalrecursions+1, NULL, theo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    sollya_mpfi_nan_normalize(result);

    return excludes;
  }

  mpfr_init2(leftX,sollya_mpfi_get_prec(x));
  mpfr_init2(rightX,sollya_mpfi_get_prec(x));

  sollya_mpfi_get_left(leftX,x);
  sollya_mpfi_get_right(rightX,x);

  printMessage(13,SOLLYA_MSG_USING_TAYLOR_EVALUATION,"Information: evaluating a function in interval arithmetic using Taylor's formula.\n");

  if (theo != NULL) {
    nullifyExprBoundTheo(theo);

    if (!isPolynomial(func)) {
      constantTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
      nullifyExprBoundTheo(constantTheo);
      linearTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
      nullifyExprBoundTheo(linearTheo);
      directTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
      nullifyExprBoundTheo(directTheo);
      theo->functionType = func->nodeType;
    } else {
      constantTheo = NULL;
      linearTheo = NULL;
      directTheo = NULL;
      theo->functionType = POLYNOMIAL;
    }
    theo->function = copyTree(func);
    theo->x = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*(theo->x),prec);
    sollya_mpfi_set(*(theo->x),x);
    theo->y = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*(theo->y),prec);
  } else {
    constantTheo = NULL;
    linearTheo = NULL;
    directTheo = NULL;
  }

  mpfr_init2(xZ,prec);
  sollya_mpfi_init2(xZI,prec);
  sollya_mpfi_init2(temp,prec);
  sollya_mpfi_init2(temp2,prec);
  sollya_mpfi_init2(constantTerm,prec);
  sollya_mpfi_init2(linearTerm,prec);
  sollya_mpfi_init2(resultTaylor,prec);
  sollya_mpfi_init2(resultDirect,prec);

  mpfr_init2(rTl,prec);
  mpfr_init2(rTr,prec);

  sollya_mpfi_mid(xZ,x);
  sollya_mpfi_set_fr(xZI,xZ);


  if ((recurse > 0) && (accessThruMemRef(func)->nodeType != DIV)) {
    nextderiv = differentiate(deriv);
    size = treeSize(nextderiv);

    if (size > DIFFSIZE) {
      printMessage(1,SOLLYA_MSG_NO_TAYLOR_EVALUATION_AS_NO_DERIVATIVE_GETS_HUGE,"Waring: during recursive Taylor evaluation the expression of a derivative has become\n");
      printMessage(1,SOLLYA_MSG_CONTINUATION,"as great that it contains more than %d nodes.\n",DIFFSIZE);
      printMessage(1,SOLLYA_MSG_CONTINUATION,"Will now stop recursive Taylor evaluation on this expression.\n");
      printMessage(2,SOLLYA_MSG_CONTINUATION,"Information: the size of the derivative is %d, we had %d recursion(s) left.\n",size,recurse-1);
      taylorExcludesLinear = evaluateI(linearTerm, deriv, x, prec, 1, hopitalrecursions+1, NULL, linearTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    } else {
      taylorExcludesLinear = evaluateITaylor(linearTerm, deriv, nextderiv, x, prec, recurse - 1, linearTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    }

    free_memory(nextderiv);
  } else {
    taylorExcludesLinear = evaluateI(linearTerm, deriv, x, prec, 1, hopitalrecursions+1, NULL, linearTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
  }

  if ((sollya_mpfi_is_nonneg(linearTerm) || sollya_mpfi_is_nonpos(linearTerm)) && sollya_mpfi_bounded_p(linearTerm)) {

    printMessage(12,SOLLYA_MSG_DERIVATIVE_DOES_NOT_CHANGE_SIGN_ON_TAYLOR_EVAL,"Information: the linear term during Taylor evaluation does not change its sign.\n");
    printMessage(12,SOLLYA_MSG_CONTINUATION,"Simplifying by taking the convex hull of the evaluations on the endpoints.\n");


    sollya_mpfi_init2(xZI2,prec);

    sollya_mpfi_set_fr(xZI,leftX);
    sollya_mpfi_set_fr(xZI2,rightX);

    directExcludes = evaluateI(resultDirect, func, xZI, prec, 0, hopitalrecursions+1, NULL, directTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    taylorExcludesConstant = evaluateI(constantTerm, func, xZI2, prec, 1, hopitalrecursions+1, NULL, constantTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

    sollya_mpfi_union(result,resultDirect,constantTerm);

    if (theo != NULL) {
      if (theo->functionType != POLYNOMIAL) {
	theo->simplificationUsed = MONOTONOCITY;
	theo->theoLeft = directTheo;
	theo->theoRight = constantTheo;
	theo->theoLeftLinear = linearTheo;
	theo->leftDerivative = copyTree(deriv);
	theo->boundLeft = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->boundLeft),prec);
	sollya_mpfi_set(*(theo->boundLeft),resultDirect);
	theo->theoRight = constantTheo;
	theo->boundRight = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->boundRight),prec);
	sollya_mpfi_set(*(theo->boundRight),constantTerm);
	theo->boundLeftLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*(theo->boundLeftLinear),prec);
	sollya_mpfi_set(*(theo->boundLeftLinear),linearTerm);
      }
    }

    excludes = concatChains(directExcludes,taylorExcludesConstant);
    excludes = concatChains(taylorExcludesLinear,excludes);

    sollya_mpfi_clear(xZI2);

  } else {

    taylorExcludesConstant = evaluateI(constantTerm, func, xZI, prec, 1, hopitalrecursions+1, NULL, constantTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

    sollya_mpfi_sub(temp, x, xZI);
    sollya_mpfi_mul(temp2, temp, linearTerm);
    sollya_mpfi_add(resultTaylor, constantTerm, temp2);
    taylorExcludes = concatChains(taylorExcludesConstant, taylorExcludesLinear);

    if (deriv != NULL)
      directExcludes = evaluateITaylorOnDiv(resultDirect, func, x, prec, recurse, directTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);
    else
      directExcludes = evaluateI(resultDirect, func, x, prec, 0, hopitalrecursions+1, NULL, directTheo,noExcludes, fastAddSub, workForThin, noLazyHooks, cutoff, lazyHookUsed);

    sollya_mpfi_get_left(rTl,resultTaylor);
    sollya_mpfi_get_right(rTr,resultTaylor);

    if (mpfr_number_p(rTl) && mpfr_number_p(rTr)) {
      sollya_mpfi_intersect(result,resultTaylor,resultDirect);
      excludes = concatChains(directExcludes,taylorExcludes);
      if (theo != NULL) {
	if (theo->functionType != POLYNOMIAL) {
	  theo->simplificationUsed = TAYLORPROOF;
	  theo->theoLeft = directTheo;
	  theo->theoLeftConstant = constantTheo;
	  theo->theoLeftLinear = linearTheo;
	  theo->leftDerivative = copyTree(deriv);
	  theo->boundLeft = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->boundLeft),prec);
	  sollya_mpfi_set(*(theo->boundLeft),resultDirect);
	  theo->boundRight = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->boundRight),prec);
	  sollya_mpfi_set(*(theo->boundRight),resultTaylor);
	  theo->boundLeftConstant = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->boundLeftConstant),prec);
	  sollya_mpfi_set(*(theo->boundLeftConstant),constantTerm);
	  theo->boundLeftLinear = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->boundLeftLinear),prec);
	  sollya_mpfi_set(*(theo->boundLeftLinear),linearTerm);
	  theo->xZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->xZ),prec);
	  sollya_mpfi_set(*(theo->xZ),xZI);
	  theo->xMXZ = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->xMXZ),prec);
	  sollya_mpfi_set(*(theo->xMXZ),temp);
	}
      }
    } else {
      sollya_mpfi_set(result,resultDirect);
      freeChain(taylorExcludes,freeMpfiPtr);
      excludes = directExcludes;
      if (theo != NULL) {
	if (theo->functionType != POLYNOMIAL) {
	  theo->simplificationUsed = IMPLICATION;
	  theo->theoLeft = directTheo;
	  theo->boundLeft = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*(theo->boundLeft),prec);
	  sollya_mpfi_set(*(theo->boundLeft),resultDirect);
	  freeExprBoundTheo(constantTheo);
	  freeExprBoundTheo(linearTheo);
	}
      }
    }
  }

  sollya_mpfi_nan_normalize(result);

  if (theo != NULL) sollya_mpfi_set(*(theo->y),result);

  mpfr_clear(xZ);
  mpfr_clear(rTl);
  mpfr_clear(rTr);
  sollya_mpfi_clear(xZI);
  sollya_mpfi_clear(temp);
  sollya_mpfi_clear(temp2);
  sollya_mpfi_clear(constantTerm);
  sollya_mpfi_clear(linearTerm);
  sollya_mpfi_clear(resultTaylor);
  sollya_mpfi_clear(resultDirect);
  mpfr_clear(leftX);
  mpfr_clear(rightX);

  return excludes;
}



chain *findZerosUnsimplified(node *func, node *deriv, sollya_mpfi_t range, mp_prec_t prec, mpfr_t diam, chain **noZeroProofs) {
  mpfr_t rangeDiam, l,m,r;
  chain *res, *leftchain, *rightchain;
  sollya_mpfi_t *temp;
  sollya_mpfi_t lI, rI, y;
  chain *excludes;
  chain *leftProofs, *rightProofs;
  chain **leftProofsPtr, **rightProofsPtr;
  exprBoundTheo *theo;

  leftProofs = NULL;
  rightProofs = NULL;
  if (noZeroProofs != NULL) {
    leftProofsPtr = &leftProofs;
    rightProofsPtr = &rightProofs;
    theo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
  } else {
    leftProofsPtr = NULL;
    rightProofsPtr = NULL;
    theo = NULL;
  }

  mpfr_init2(rangeDiam,prec);
  sollya_mpfi_diam_abs(rangeDiam,range);

  mpfr_init2(l,prec);
  mpfr_init2(r,prec);
  sollya_mpfi_get_left(l,range);
  sollya_mpfi_get_right(r,range);
  mpfr_nextabove(l); mpfr_nextabove(l); mpfr_nextabove(l); mpfr_nextabove(l);

  if ((mpfr_cmp(rangeDiam,diam) <= 0) || (mpfr_cmp(l,r) >= 0)) {
    res = (chain *) safeMalloc(sizeof(chain));
    res->next = NULL;
    temp = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*temp,prec);
    sollya_mpfi_set(*temp,range);
    res->value = temp;
    if (theo != NULL) freeExprBoundTheo(theo);
  } else {
    sollya_mpfi_init2(y,prec);
    excludes = evaluateITaylor(y, func, deriv, range, prec, taylorrecursions, theo,1,0,0,0,NULL,NULL);
    freeChain(excludes,freeMpfiPtr);
    if (!sollya_mpfi_bounded_p(y)) {
      printMessage(1,SOLLYA_MSG_NAN_OR_INF_ON_DERIVATIVE,"Warning: during zero-search the derivative of the function evaluated to NaN or Inf in the interval %w.\nThe function might not be continuously differentiable in this interval.\n",range);
    }
    if ((!sollya_mpfi_bounded_p(y)) || sollya_mpfi_has_zero(y)) {
      mpfr_init2(m,prec);
      sollya_mpfi_get_left(l,range);

      sollya_mpfi_mid(m,range);
      sollya_mpfi_init2(lI,prec);
      sollya_mpfi_init2(rI,prec);

      if (mpfr_cmp(l, m)<=0) sollya_mpfi_interv_fr(lI,l,m);
      else sollya_mpfi_interv_fr(lI,m,l);

      if (mpfr_cmp(m, r)<=0) sollya_mpfi_interv_fr(rI,m,r);
      else sollya_mpfi_interv_fr(rI,r,m);

      if (theo != NULL) freeExprBoundTheo(theo);

      leftchain = findZerosUnsimplified(func,deriv,lI,prec,diam,leftProofsPtr);
      rightchain = findZerosUnsimplified(func,deriv,rI,prec,diam,rightProofsPtr);

      res = concatChains(leftchain,rightchain);

      if (noZeroProofs != NULL) {
	*noZeroProofs = concatChains(leftProofs,rightProofs);
      }


      mpfr_clear(m);
      sollya_mpfi_clear(lI);
      sollya_mpfi_clear(rI);
    } else {
      res = NULL;
      if (noZeroProofs != NULL) *noZeroProofs = addElement(*noZeroProofs,theo);
    }
    sollya_mpfi_clear(y);
  }
  mpfr_clear(l);
  mpfr_clear(r);
  mpfr_clear(rangeDiam);
  return res;
}

chain *findZeros(node *func, node *deriv, sollya_mpfi_t range, mp_prec_t prec, mpfr_t diam, noZeroTheo *theo) {
  node *funcSimplified, *derivSimplified;
  chain *temp;
  chain **noZeroProofs;

  funcSimplified = horner(func);
  derivSimplified = horner(deriv);

  if (theo != NULL) {
    theo->function = copyTree(func);
    theo->derivative = copyTree(deriv);
    theo->funcEqual = (equalityTheo *) safeMalloc(sizeof(equalityTheo));
    theo->funcEqual->expr1 = copyTree(func);
    theo->funcEqual->expr2 = copyTree(funcSimplified);
    theo->derivEqual = (equalityTheo *) safeMalloc(sizeof(equalityTheo));
    theo->derivEqual->expr1 = copyTree(deriv);
    theo->derivEqual->expr2 = copyTree(derivSimplified);
    noZeroProofs = &(theo->exprBoundTheos);
  } else {
    noZeroProofs = NULL;
  }

  printMessage(3,SOLLYA_MSG_INVOKING_RECURSIVE_INTERVAL_ZERO_SEARCH,"Information: invoking the recursive interval zero search.\n");
  temp = findZerosUnsimplified(funcSimplified,derivSimplified,range,prec,diam,noZeroProofs);
  printMessage(3,SOLLYA_MSG_RECURSIVE_INTERVAL_ZERO_SEARCH_HAS_FINISHED,"Information: the recursive interval zero search has finished.\n");

  free_memory(funcSimplified);
  free_memory(derivSimplified);

  return temp;
}



void printInterval(sollya_mpfi_t interval) {
  mpfr_t l,r;
  mp_prec_t prec;
  char *temp_string;

  prec = sollya_mpfi_get_prec(interval);
  mpfr_init2(l,prec);
  mpfr_init2(r,prec);
  sollya_mpfi_get_left(l,interval);
  sollya_mpfi_get_right(r,interval);

  if ((dyadic == 0) && (midpointMode == 1)) {
    temp_string = sprintMidpointMode(l, r);
    if (temp_string != NULL) {
      sollyaPrintf("%s ",temp_string);
      safeFree(temp_string);
    } else {
      sollyaPrintf("[");
      printValue(&l);
      sollyaPrintf(";");
      printValue(&r);
      sollyaPrintf("]");
    }
  } else {
    sollyaPrintf("[");
    printValue(&l);
    sollyaPrintf(";");
    printValue(&r);
    sollyaPrintf("]");
  }

  mpfr_clear(l);
  mpfr_clear(r);
}

char *sprintInterval(sollya_mpfi_t interval) {
  mpfr_t l,r;
  mp_prec_t prec;
  char *temp_string, *temp_string2;
  char *res;

  prec = sollya_mpfi_get_prec(interval);
  mpfr_init2(l,prec);
  mpfr_init2(r,prec);
  sollya_mpfi_get_left(l,interval);
  sollya_mpfi_get_right(r,interval);

  if ((dyadic == 0) && (midpointMode == 1)) {
    temp_string = sprintMidpointMode(l, r);
    if (temp_string != NULL) {
      res = temp_string;
    } else {
      temp_string = sprintValue(&l);
      temp_string2 = sprintValue(&r);
      res = safeCalloc(strlen(temp_string) + strlen(temp_string2) + 3 + 1,
		       sizeof(char));
      sprintf(res,"[%s;%s]",temp_string,temp_string2);
      safeFree(temp_string);
      safeFree(temp_string2);
    }
  } else {
    temp_string = sprintValue(&l);
    temp_string2 = sprintValue(&r);
    res = safeCalloc(strlen(temp_string) + strlen(temp_string2) + 3 + 1,
		     sizeof(char));
    sprintf(res,"[%s;%s]",temp_string,temp_string2);
    safeFree(temp_string);
    safeFree(temp_string2);
  }

  mpfr_clear(l);
  mpfr_clear(r);

  return res;
}




void fprintInterval(FILE *fd, sollya_mpfi_t interval) {
  mpfr_t l,r;
  mp_prec_t prec;


  prec = sollya_mpfi_get_prec(interval);
  mpfr_init2(l,prec);
  mpfr_init2(r,prec);
  sollya_mpfi_get_left(l,interval);
  sollya_mpfi_get_right(r,interval);
  sollyaFprintf(fd,"[");
  fprintValue(fd,l);
  sollyaFprintf(fd,";");
  fprintValue(fd,r);
  sollyaFprintf(fd,"]");

  mpfr_clear(l);
  mpfr_clear(r);
}


chain *joinAdjacentIntervals(chain *intervals, mpfr_t diam) {
  chain *newChain, *curr;
  sollya_mpfi_t *tempI;
  mp_prec_t prec, p;
  mpfr_t newLeft, newRight, l,r, mpfr_temp;

  if (intervals == NULL) return NULL;
  if (intervals->next == NULL) {
    tempI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*tempI,sollya_mpfi_get_prec(*((sollya_mpfi_t *) (intervals->value))));
    sollya_mpfi_set(*tempI,*((sollya_mpfi_t *) (intervals->value)));
    newChain = addElement(NULL,tempI);
    return newChain;
  }

  prec = sollya_mpfi_get_prec(*((sollya_mpfi_t *) (intervals->value)));
  curr = intervals->next;
  while (curr != NULL) {
    p = sollya_mpfi_get_prec(*((sollya_mpfi_t *) (curr->value)));
    if (p > prec) prec = p;
    curr = curr->next;
  }


  mpfr_init2(newLeft,prec);
  mpfr_init2(newRight,prec);
  sollya_mpfi_get_left(newLeft,*((sollya_mpfi_t *) (intervals->value)));
  sollya_mpfi_get_right(newRight,*((sollya_mpfi_t *) (intervals->value)));
  mpfr_init2(l,prec);
  mpfr_init2(r,prec);
  mpfr_init2(mpfr_temp,prec);

  newChain = NULL;
  curr = intervals->next;
  while (curr != NULL) {
    sollya_mpfi_get_left(l,*((sollya_mpfi_t *) (curr->value)));
    sollya_mpfi_get_right(r,*((sollya_mpfi_t *) (curr->value)));
    mpfr_sub(mpfr_temp,r,newLeft,GMP_RNDN);
    if ((mpfr_cmp(l,newRight) == 0) && (mpfr_cmp(mpfr_temp,diam) < 0)) {
      mpfr_set(newRight,r,GMP_RNDN);
    } else {
      tempI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
      sollya_mpfi_init2(*tempI,prec);
      sollya_mpfi_interv_fr(*tempI,newLeft,newRight);
      newChain = addElement(newChain,tempI);
      mpfr_set(newLeft,l,GMP_RNDN);
      mpfr_set(newRight,r,GMP_RNDN);
    }
    curr = curr->next;
  }
  tempI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
  sollya_mpfi_init2(*tempI,prec);
  sollya_mpfi_interv_fr(*tempI,newLeft,newRight);
  newChain = addElement(newChain,tempI);

  mpfr_clear(l);
  mpfr_clear(r);
  mpfr_clear(newLeft);
  mpfr_clear(newRight);
  mpfr_clear(mpfr_temp);
  return newChain;
}

chain *joinAdjacentIntervalsMaximally(chain *intervals) {
  chain *newChain, *curr;
  sollya_mpfi_t *tempI;
  mp_prec_t prec, p;
  mpfr_t newLeft, newRight, l,r, mpfr_temp;

  if (intervals == NULL) return NULL;
  if (intervals->next == NULL) {
    tempI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*tempI,sollya_mpfi_get_prec(*((sollya_mpfi_t *) (intervals->value))));
    sollya_mpfi_set(*tempI,*((sollya_mpfi_t *) (intervals->value)));
    newChain = addElement(NULL,tempI);
    return newChain;
  }

  prec = sollya_mpfi_get_prec(*((sollya_mpfi_t *) (intervals->value)));
  curr = intervals->next;
  while (curr != NULL) {
    p = sollya_mpfi_get_prec(*((sollya_mpfi_t *) (curr->value)));
    if (p > prec) prec = p;
    curr = curr->next;
  }


  mpfr_init2(newLeft,prec);
  mpfr_init2(newRight,prec);
  sollya_mpfi_get_left(newLeft,*((sollya_mpfi_t *) (intervals->value)));
  sollya_mpfi_get_right(newRight,*((sollya_mpfi_t *) (intervals->value)));
  mpfr_init2(l,prec);
  mpfr_init2(r,prec);
  mpfr_init2(mpfr_temp,prec);

  newChain = NULL;
  curr = intervals->next;
  while (curr != NULL) {
    sollya_mpfi_get_left(l,*((sollya_mpfi_t *) (curr->value)));
    sollya_mpfi_get_right(r,*((sollya_mpfi_t *) (curr->value)));
    mpfr_sub(mpfr_temp,r,newLeft,GMP_RNDN);
    if (mpfr_cmp(l,newRight) == 0) {
      mpfr_set(newRight,r,GMP_RNDN);
    } else {
      tempI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
      sollya_mpfi_init2(*tempI,prec);
      sollya_mpfi_interv_fr(*tempI,newLeft,newRight);
      newChain = addElement(newChain,tempI);
      mpfr_set(newLeft,l,GMP_RNDN);
      mpfr_set(newRight,r,GMP_RNDN);
    }
    curr = curr->next;
  }
  tempI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
  sollya_mpfi_init2(*tempI,prec);
  sollya_mpfi_interv_fr(*tempI,newLeft,newRight);
  newChain = addElement(newChain,tempI);

  mpfr_clear(l);
  mpfr_clear(r);
  mpfr_clear(newLeft);
  mpfr_clear(newRight);
  mpfr_clear(mpfr_temp);
  return newChain;
}



chain *excludeIntervals(chain *mainIntervals, chain *excludeIntervals) {
  chain *curr, *previous, *curr2, *temp;
  sollya_mpfi_t *interval, *exclude;
  mp_prec_t prec, p;
  mpfr_t il, ir, el, er;

  if (mainIntervals == NULL) return NULL;
  if (excludeIntervals == NULL) return mainIntervals;

  prec = 1;
  curr = mainIntervals;
  while (curr != NULL) {
    p = sollya_mpfi_get_prec(*((sollya_mpfi_t *) curr->value));
    if (p > prec) prec = p;
    curr = curr->next;
  }
  curr = excludeIntervals;
  while (curr != NULL) {
    p = sollya_mpfi_get_prec(*((sollya_mpfi_t *) curr->value));
    if (p > prec) prec = p;
    curr = curr->next;
  }

  prec += 5;

  mpfr_init2(il,prec);
  mpfr_init2(ir,prec);
  mpfr_init2(el,prec);
  mpfr_init2(er,prec);

  curr2 = excludeIntervals;
  while (curr2 != NULL) {
    exclude = (sollya_mpfi_t *) (curr2->value);
    sollya_mpfi_get_left(el,*exclude);
    sollya_mpfi_get_right(er,*exclude);
    curr = mainIntervals;
    previous = NULL;
    while (curr != NULL) {
      interval = (sollya_mpfi_t *) (curr->value);
      sollya_mpfi_get_left(il,*interval);
      sollya_mpfi_get_right(ir,*interval);
      if ((mpfr_cmp(el,ir) < 0) && (mpfr_cmp(il,er) < 0)) { /* [il;ir] inter [el;er] != empty */
	if ((mpfr_cmp(il,el) < 0) && (mpfr_cmp(er,ir) < 0)) {
	  /* We must produce two intervals [il;el] and [er;ir] */
	  sollya_mpfi_interv_fr(*interval,il,el);
	  interval = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	  sollya_mpfi_init2(*interval,prec);
	  sollya_mpfi_interv_fr(*interval,er,ir);
	  temp = (chain *) safeMalloc(sizeof(chain));
	  temp->value = interval;
	  temp->next = curr->next;
	  curr->next = temp;
	} else {
	  if (mpfr_cmp(il,el) < 0) {
	    /* We must produce one interval [il;el] */
	    sollya_mpfi_interv_fr(*interval,il,el);
	  } else {
	    if (mpfr_cmp(er,ir) < 0) {
	      /* We must produce one interval [er;ir] */
	      sollya_mpfi_interv_fr(*interval,er,ir);
	    } else {
	      /* We must remove the interval completely */
	      if (previous != NULL) {
		/* We are not the first interval in the chain */
		previous->next = curr->next;
		sollya_mpfi_clear(*interval);
		safeFree(interval);
		safeFree(curr);
		curr = previous;
	      } else {
		/* We are the first interval in the chain */
		if (curr->next != NULL) {
		  /* We have a successor that will become the head of the chain */
		  mainIntervals = curr->next;
		  sollya_mpfi_clear(*interval);
		  safeFree(interval);
		  safeFree(curr);
		  curr = mainIntervals;
		} else {
		  /* We are the first and the last element in the chain, which will be empty */

		  sollya_mpfi_clear(*interval);
		  safeFree(interval);
		  safeFree(curr);
		  mpfr_clear(il);
		  mpfr_clear(ir);
		  mpfr_clear(el);
		  mpfr_clear(er);
		  return NULL;
		}
	      }
	    }
	  }
	}
      }
      previous = curr;
      curr = curr->next;
    }
    curr2 = curr2->next;
  }

  mpfr_clear(il);
  mpfr_clear(ir);
  mpfr_clear(el);
  mpfr_clear(er);

  return mainIntervals;
}


void infnormI(sollya_mpfi_t infnormval, node *func, node *deriv,
	      node *numeratorDeriv, node *derivNumeratorDeriv,
	      sollya_mpfi_t range, mp_prec_t prec, mpfr_t diam,
	      chain *intervalsToExclude,
	      chain **mightExcludes,
	      infnormTheo *theo) {
  chain *curr, *zeros, *tempChain, *tempChain2, *tempChain3;
  sollya_mpfi_t *currInterval;
  sollya_mpfi_t evalFuncOnInterval, lInterv, rInterv;
  mpfr_t innerLeft, innerRight, outerLeft, outerRight, l, r, tl, tr;
  mpfr_t diamJoin;
  mp_prec_t rangePrec;
  chain *excludes, *excludesTemp;
  int i;
  noZeroTheo *noZeros;
  exprBoundTheo *evalLeftBound, *evalRightBound, *currZeroTheo;

  currZeroTheo = NULL;
  if (theo != NULL) {
    theo->function = copyTree(func);
    theo->derivative = copyTree(deriv);
    theo->numeratorOfDerivative = copyTree(numeratorDeriv);
    theo->derivativeOfNumeratorOfDerivative = copyTree(derivNumeratorDeriv);
    theo->excludedIntervals = copyChain(intervalsToExclude,copyMpfiPtr);
    noZeros = (noZeroTheo *) safeCalloc(1,sizeof(noZeroTheo));
    theo->noZeros = noZeros;
    evalLeftBound = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
    nullifyExprBoundTheo(evalLeftBound);
    evalRightBound = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
    nullifyExprBoundTheo(evalRightBound);
    theo->evalLeftBound = evalLeftBound;
    theo->evalRightBound = evalRightBound;
    theo->domain = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    theo->infnorm = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    theo->evalOnZeros = NULL;
    sollya_mpfi_init2(*(theo->domain),sollya_mpfi_get_prec(range));
    sollya_mpfi_init2(*(theo->infnorm),sollya_mpfi_get_prec(infnormval));
    sollya_mpfi_set(*(theo->domain),range);
  } else {
    noZeros = NULL;
    evalLeftBound = NULL;
    evalRightBound = NULL;
  }

  mpfr_init2(innerLeft, prec);
  mpfr_init2(innerRight, prec);
  mpfr_init2(outerLeft, prec);
  mpfr_init2(outerRight, prec);
  mpfr_init2(tl, prec);
  mpfr_init2(tr, prec);
  sollya_mpfi_init2(evalFuncOnInterval,prec);

  rangePrec = sollya_mpfi_get_prec(range);
  mpfr_init2(l,rangePrec);
  mpfr_init2(r,rangePrec);
  sollya_mpfi_get_left(l,range);
  sollya_mpfi_get_right(r,range);
  sollya_mpfi_init2(rInterv, rangePrec);
  sollya_mpfi_init2(lInterv, rangePrec);
  sollya_mpfi_set_fr(rInterv,r);
  sollya_mpfi_set_fr(lInterv,l);

  excludes = evaluateITaylor(evalFuncOnInterval, func, deriv, lInterv, prec, taylorrecursions, evalLeftBound,0,0,0,0,NULL,NULL);
  sollya_mpfi_get_left(outerLeft,evalFuncOnInterval);
  sollya_mpfi_get_right(outerRight,evalFuncOnInterval);
  mpfr_set(innerLeft,outerRight,GMP_RNDU);
  mpfr_set(innerRight,outerLeft,GMP_RNDD);
  excludesTemp = evaluateITaylor(evalFuncOnInterval, func, deriv, rInterv, prec, taylorrecursions, evalRightBound,0,0,0,0,NULL,NULL);
  excludes = concatChains(excludes,excludesTemp);
  sollya_mpfi_get_left(tl,evalFuncOnInterval);
  sollya_mpfi_get_right(tr,evalFuncOnInterval);
  sollya_mpfr_min(outerLeft,outerLeft,tl,GMP_RNDD);
  sollya_mpfr_max(outerRight,outerRight,tr,GMP_RNDU);
  sollya_mpfr_min(innerLeft,innerLeft,tr,GMP_RNDU);
  sollya_mpfr_max(innerRight,innerRight,tl,GMP_RNDD);

  printMessage(3,SOLLYA_MSG_INVOKING_RECURSIVE_INTERVAL_ZERO_SEARCH,"Information: invoking interval zero search.\n");
  tempChain = findZeros(numeratorDeriv,derivNumeratorDeriv,range,prec,diam,noZeros);

  printMessage(3,SOLLYA_MSG_RECURSIVE_INTERVAL_ZERO_SEARCH_HAS_FINISHED,"Information: interval zero search is done.\n");
  mpfr_init2(diamJoin,prec);
  mpfr_mul_2ui(diamJoin,diam,3,GMP_RNDN);
  tempChain2 = joinAdjacentIntervals(tempChain,diamJoin);
  tempChain3 = copyChain(tempChain2,copyMpfiPtr);
  mpfr_mul_2ui(diamJoin,diamJoin,2,GMP_RNDN);
  zeros = joinAdjacentIntervals(tempChain3,diamJoin);
  mpfr_clear(diamJoin);

  zeros = excludeIntervals(zeros,intervalsToExclude);

  i = 0;
  for (curr=zeros;curr!=NULL;curr=curr->next) i++;
  printMessage(2,SOLLYA_MSG_CERTAIN_NUM_OF_INTVALS_ENCLOSING_ZEROS_OF_DERIV,
               "Information: %d interval(s) have (has) been found that possibly contain(s) the zeros of the derivative.\n",i);

  curr = zeros;
  while (curr != NULL) {

    if (theo != NULL) {
      currZeroTheo = (exprBoundTheo *) safeCalloc(1,sizeof(exprBoundTheo));
      nullifyExprBoundTheo(currZeroTheo);
    } else {
      currZeroTheo = NULL;
    }
    currInterval = ((sollya_mpfi_t *) (curr->value));
    excludesTemp = evaluateITaylor(evalFuncOnInterval, func, deriv, *currInterval, prec, taylorrecursions, currZeroTheo,0,0,0,0,NULL,NULL);

    excludes = concatChains(excludes,excludesTemp);
    sollya_mpfi_get_left(tl,evalFuncOnInterval);
    sollya_mpfi_get_right(tr,evalFuncOnInterval);

    if (theo != NULL) {
      theo->evalOnZeros = addElement(theo->evalOnZeros,currZeroTheo);
    }

    if (mpfr_nan_p(tl) || mpfr_nan_p(tr)) {
      printMessage(1,SOLLYA_MSG_NAN_OR_INF_ON_DERIVATIVE,"Warning: NaNs occurred during the interval evaluation of the zeros of the derivative.\n");
    }

    sollya_mpfr_min(outerLeft,outerLeft,tl,GMP_RNDD);
    sollya_mpfr_max(outerRight,outerRight,tr,GMP_RNDU);
    sollya_mpfr_min(innerLeft,innerLeft,tr,GMP_RNDU);
    sollya_mpfr_max(innerRight,innerRight,tl,GMP_RNDD);
    curr = curr->next;
  }

  freeChain(zeros,freeMpfiPtr);
  freeChain(tempChain,freeMpfiPtr);
  freeChain(tempChain2,freeMpfiPtr);
  freeChain(tempChain3,freeMpfiPtr);

  if (mpfr_cmp(innerLeft,innerRight) >= 0) {
    mpfr_neg(outerLeft,outerLeft,GMP_RNDN);
    sollya_mpfr_max(tr,outerLeft,outerRight,GMP_RNDU);
    mpfr_set_d(tl,0.0,GMP_RNDD);

    if (mpfr_cmp(tl, tr)<=0) sollya_mpfi_interv_fr(infnormval,tl,tr);
    else sollya_mpfi_interv_fr(infnormval,tr,tl);

  } else {
    mpfr_neg(innerLeft,innerLeft,GMP_RNDN);
    mpfr_neg(outerLeft,outerLeft,GMP_RNDN);
    sollya_mpfr_max(tl,innerLeft,innerRight,GMP_RNDD);
    sollya_mpfr_max(tr,outerLeft,outerRight,GMP_RNDU);
    if (mpfr_cmp(tl, tr)<=0) sollya_mpfi_interv_fr(infnormval,tl,tr);
    else sollya_mpfi_interv_fr(infnormval,tr,tl);
  }

  if (mightExcludes == NULL) {
    freeChain(excludes,freeMpfiPtr);
  } else {
    *mightExcludes = excludes;
  }

  if (theo != NULL) {
    sollya_mpfi_set(*(theo->infnorm),infnormval);
  }

  mpfr_clear(tl);
  mpfr_clear(tr);
  sollya_mpfi_clear(lInterv);
  sollya_mpfi_clear(rInterv);
  mpfr_clear(l);
  mpfr_clear(r);
  sollya_mpfi_clear(evalFuncOnInterval);
  mpfr_clear(innerLeft);
  mpfr_clear(innerRight);
  mpfr_clear(outerLeft);
  mpfr_clear(outerRight);
}


int isTrivialInfnormCase(rangetype result, node *func) {
  int isTrivial;
  node *simplifiedFunc, *numerator, *denominator;

  isTrivial = 0;
  simplifiedFunc = horner(func);

  if (accessThruMemRef(simplifiedFunc)->nodeType == CONSTANT) {
    mpfr_set(*(result.a),*(accessThruMemRef(simplifiedFunc)->value),GMP_RNDD);
    mpfr_set(*(result.b),*(accessThruMemRef(simplifiedFunc)->value),GMP_RNDU);
    isTrivial = 1;
  } else {
    if (accessThruMemRef(simplifiedFunc)->nodeType == SUB) {
      if (isSyntacticallyEqual(accessThruMemRef(simplifiedFunc)->child1,accessThruMemRef(simplifiedFunc)->child2)) {
	mpfr_set_d(*(result.a),0.0,GMP_RNDN);
	mpfr_set_d(*(result.b),0.0,GMP_RNDN);
	isTrivial = 1;
      }
    } else {
      if (getNumeratorDenominator(&numerator, &denominator, simplifiedFunc)) {
	if ((accessThruMemRef(numerator)->nodeType == CONSTANT) &&
	    mpfr_zero_p(*(accessThruMemRef(numerator)->value))) {
	  mpfr_set_d(*(result.a),0.0,GMP_RNDN);
	  mpfr_set_d(*(result.b),0.0,GMP_RNDN);
	  isTrivial = 1;
	} else {
	  if (isSyntacticallyEqual(numerator, denominator)) {
	    mpfr_set_d(*(result.a),1.0,GMP_RNDN);
	    mpfr_set_d(*(result.b),1.0,GMP_RNDN);
	    isTrivial = 1;
	  }
	}
	free_memory(denominator);
      }
      free_memory(numerator);
    }
  }

  free_memory(simplifiedFunc);
  return isTrivial;
}

/* If x_result != NULL, the mpfr_t it points to is used to store the value of a point where the infnorm is reached.
   It is supposed to be already initialized.
*/
void uncertifiedInfnorm(mpfr_t result, mpfr_t *x_result, node *f, mpfr_t a, mpfr_t b, unsigned long int points, mp_prec_t prec) {
  mpfr_t current_x, argmax, x1, x2, x3, step, y1, y2, y3, max, cutoff;
  mpfr_t ystar, y1diff, y3diff, xstar;
  mpfr_t zero_mpfr;
  mpfr_t perturb;

  mp_prec_t prec_bound = prec;
  node *f_diff, *f_diff2;
  int r;
  int count=0;
  int stop_algo = 0;
  const mp_prec_t INIT_PREC = 20;
  gmp_randstate_t random_state;


  if (mpfr_get_prec(a) > prec_bound) prec_bound = mpfr_get_prec(a);
  if (mpfr_get_prec(b) > prec_bound) prec_bound = mpfr_get_prec(b);


  /**************************** Dealing with special cases ****************************/
  if ( (!mpfr_number_p(a)) || (!mpfr_number_p(b)) ) {
    printMessage(1,SOLLYA_MSG_DOMAIN_IS_NO_CLOSED_INTERVAL_ON_THE_REALS,"Warning: a bound of the interval is infinite or NaN.\n");
    printMessage(1,SOLLYA_MSG_CONTINUATION,"This command cannot handle such intervals.\n");
    mpfr_set_nan(result);
    return;
  }

  if (mpfr_equal_p(a,b)) {
    printMessage(1,SOLLYA_MSG_DOMAIN_IS_REDUCED_TO_A_POINT_WILL_SIMPLY_EVAL,"Warning: the given interval is reduced to one point.\n");
    evaluateFaithful(result,f,a,prec);
    mpfr_abs(result,result,GMP_RNDU);
    return;
  }

  if (mpfr_greater_p(a,b)) {
    printMessage(1,SOLLYA_MSG_DOMAIN_IS_EMPTY,"Warning: the interval is empty.\n");
    mpfr_set_d(result,0.,GMP_RNDN);
    return;
  }

  if (isConstant(f)) {
    printMessage(1,SOLLYA_MSG_EXPRESSION_IS_CONSTANT,"Warning: the expression is constant.\n");
    evaluateFaithful(result,f,a,prec);
    mpfr_abs(result,result,GMP_RNDU);
    return;
  }
  /************************************************************************************/


  gmp_randinit_default(random_state);
  gmp_randseed_ui(random_state, 65845285);
  mpfr_init2(perturb, prec);

  mpfr_init2(zero_mpfr, 53);
  mpfr_set_d(zero_mpfr, 0., GMP_RNDN);

  mpfr_init2(step, prec);
  mpfr_sub(step, b, a, GMP_RNDU); /* since a<b and step is computed with rounding upwards, step>O */
  mpfr_div_ui(step, step, points, GMP_RNDN);

  mpfr_init2(current_x, prec_bound);
  mpfr_init2(x1, prec_bound);
  mpfr_init2(x2, prec_bound);
  mpfr_init2(x3, prec_bound);
  mpfr_init2(argmax, prec_bound);
  mpfr_init2(y1, prec);
  mpfr_init2(y2, prec);
  mpfr_init2(y3, prec);

  mpfr_init2(ystar, prec);
  mpfr_init2(y1diff, 5);  /* Only the sign of this variable is required */
  mpfr_init2(y3diff, 5);  /* Only the sign of this variable is required */
  mpfr_init2(xstar, prec_bound);

  mpfr_init2(max, prec);
  mpfr_set_d(max, 0., GMP_RNDN);
  mpfr_init2(cutoff, prec);
  mpfr_set_d(cutoff, 0., GMP_RNDN);

  f_diff = differentiate(f);
  f_diff2 = NULL;


  /* Initial value of x1 */
  mpfr_set(x1, a, GMP_RNDN); /* exact */
  mpfr_set(current_x, a, GMP_RNDN);
  do {
    count++;
    if (mpfr_greaterequal_p(x1,b)) {
      mpfr_set(x1, b, GMP_RNDN); /* exact */
      stop_algo = 1;
    }
    r = evaluateFaithfulWithCutOffFast(y1, f, f_diff, x1, cutoff, INIT_PREC);
    if (r==0) mpfr_set_d(y1, 0. , GMP_RNDN);

    if (!mpfr_number_p(y1)) {
      printMessage(1,SOLLYA_MSG_EVALUATION_AT_POINT_GIVES_NAN_EXCLUDING_POINT,"Warning: the evaluation of the given function %b in %v gives NaN.\nThis (possibly maximum) point will be excluded from the infnorm result.\n",f,x1);
      mpfr_add(current_x, current_x, step, GMP_RNDU); /* rounding up ensures that x1(new) > x1(old) */
      mpfr_urandomb(perturb, random_state); mpfr_mul_2ui(perturb, perturb, 1, GMP_RNDN);
      mpfr_sub_ui(perturb, perturb, 1, GMP_RNDN); mpfr_div_2ui(perturb, perturb, 2, GMP_RNDN);
      mpfr_mul(perturb, perturb, step, GMP_RNDN); /* perturb \in [-step/4; step/4] */
      mpfr_add(x1, current_x, perturb, GMP_RNDU);
    }
  } while ( (!mpfr_number_p(y1)) && (!stop_algo) );

  mpfr_abs(max, y1, GMP_RNDU);
  mpfr_set(argmax, x1, GMP_RNDN);
  printMessage(3,SOLLYA_MSG_THE_CURRENT_MAXIMUM_IS_A_CERTAIN_VALUE,"Information: current max is %v and is reached at %v\n",max,argmax);

  mpfr_div_2ui(cutoff, max, 1, GMP_RNDU);


  /* Initial value of x2 */
  do {
    count++;
    mpfr_add(current_x, current_x, step, GMP_RNDU); /* rounding up ensures that x2 > x1 */
    mpfr_urandomb(perturb, random_state); mpfr_mul_2ui(perturb, perturb, 1, GMP_RNDN);
    mpfr_sub_ui(perturb, perturb, 1, GMP_RNDN); mpfr_div_2ui(perturb, perturb, 2, GMP_RNDN);
    mpfr_mul(perturb, perturb, step, GMP_RNDN); /* perturb \in [-step/4; step/4] */
    mpfr_add(x2, current_x, perturb, GMP_RNDU);

    if (mpfr_greaterequal_p(x2,b)) {
      mpfr_set(x2, b, GMP_RNDN); /* exact */
      stop_algo = 1;
    }

    r = evaluateFaithfulWithCutOffFast(y2, f, f_diff, x2, cutoff, INIT_PREC);
    if (r==2) mpfr_set_d(y2, 0. , GMP_RNDN); /* under the cutoff */

    if (!mpfr_number_p(y2)) {
      printMessage(1,SOLLYA_MSG_EVALUATION_AT_POINT_GIVES_NAN_EXCLUDING_POINT,"Warning: the evaluation of the given function %b in %v gives NaN.\nThis (possibly maximum) point will be excluded from the infnorm result.\n",f,x2);
    }
  } while ( ( (!mpfr_number_p(y2)) || mpfr_equal_p(y1,y2) )
	    &&
	    (!stop_algo)
	    );

  if (mpfr_cmpabs(y2, max) > 0) { /* evaluates to false when y2=NaN */
    mpfr_abs(max, y2, GMP_RNDU);
    mpfr_set(argmax, x2, GMP_RNDN);
    printMessage(3,SOLLYA_MSG_THE_CURRENT_MAXIMUM_IS_A_CERTAIN_VALUE,"Information: current max is %v and is reached at %v\n",max,argmax);
    mpfr_div_2ui(cutoff, max, 1, GMP_RNDU);
  }

  /* Remark: at this point, it is possible that y1=y2 and stop_algo is true. */
  /* In this case, it means that f is constant over the interval and it is   */
  /* hence not necessary to run Newton's algorithm anyway.                   */

  /* Main loop */
  while(!stop_algo) {
    do {
      count++;
      if (verbosity >= 2) {
	if( count % 100 == 0) printMessage(2,SOLLYA_MSG_CERTAIN_AMOUNT_OF_CASES_HANDLED,"Information: %d out of %d points have been handled.\n",count,points);
      }

      mpfr_add(current_x, current_x, step, GMP_RNDU); /* rounding up ensures that x3 > x2 */
      mpfr_urandomb(perturb, random_state); mpfr_mul_2ui(perturb, perturb, 1, GMP_RNDN);
      mpfr_sub_ui(perturb, perturb, 1, GMP_RNDN); mpfr_div_2ui(perturb, perturb, 2, GMP_RNDN);
      mpfr_mul(perturb, perturb, step, GMP_RNDN); /* perturb \in [-step/4; step/4] */
      mpfr_add(x3, current_x, perturb, GMP_RNDU);

      if (mpfr_greaterequal_p(x3,b)) {
	mpfr_set(x3, b, GMP_RNDN); /* exact */
	stop_algo = 1;
      }

      r = evaluateFaithfulWithCutOffFast(y3, f, f_diff, x3, cutoff, INIT_PREC);
      if (r==2) mpfr_set_d(y3, 0. , GMP_RNDN); /* under the cutoff */

      if (!mpfr_number_p(y3)) {
	printMessage(1,SOLLYA_MSG_EVALUATION_AT_POINT_GIVES_NAN_EXCLUDING_POINT,"Warning: the evaluation of the given function %b in %v gives NaN.\nThis (possibly maximum) point will be excluded from the infnorm result.\n",f,x3);
      }
    } while ( ( (!mpfr_number_p(y3))|| mpfr_equal_p(y2,y3) )
	      &&
	      (!stop_algo)
	      );

    if (mpfr_cmpabs(y3, max) > 0) { /* evaluates to false when y3=NaN */
      mpfr_abs(max, y3, GMP_RNDU);
      mpfr_set(argmax, x3, GMP_RNDN);
      printMessage(3,SOLLYA_MSG_THE_CURRENT_MAXIMUM_IS_A_CERTAIN_VALUE,"Information: current max is %v and is reached at %v\n",max,argmax);
      mpfr_div_2ui(cutoff, max, 1, GMP_RNDU);
    }

    /* Call to Newton's algorithm if necessary */
    if ( (mpfr_cmpabs(y2,y1)>=0) && (mpfr_cmpabs(y2,y3)>=0) && (mpfr_cmp_d(y2,0.)!=0) ) {

      if (f_diff2 == NULL) f_diff2 = differentiate(f_diff);

      r = evaluateFaithfulWithCutOffFast(y1diff, f_diff, f_diff2, x1, zero_mpfr, prec+10);
      if (r==0) mpfr_set_d(y1diff, 0. , GMP_RNDN);

      r = evaluateFaithfulWithCutOffFast(y3diff, f_diff, f_diff2, x3, zero_mpfr, prec+10);
      if (r==0) mpfr_set_d(y3diff, 0. , GMP_RNDN);

      if ( (!mpfr_number_p(y1diff)) || (!mpfr_number_p(y3diff)) ) {
	printMessage(1,SOLLYA_MSG_EVALUATION_OF_DERIVATIVE_GIVES_NAN_NO_NEWTON,"Warning: the evaluation of the derivative %b of the given function in %v or %v gives NaN.\nNewton's algorithm will not be used on this interval.\n",f_diff,x1,x3);
      }
      else if(mpfr_sgn(y1diff)*mpfr_sgn(y3diff)<0) { /* If y1diff=0 or y3diff=0, there is no need to */
                                                     /* use Newton's algorithm since we already have */
                                                     /* the zero. Moreover, note that y1 and y3 have */
                                                     /* already been taken into account in max.      */

	findZero(xstar, f_diff, f_diff2, x1, x3, mpfr_sgn(y1diff), NULL, 0, (prec_bound >> 1) + 10);

	/* If xstar = NaN, a warning has already been produced by Newton's algorithm. */
	/* There is no need to print a warning again here.                            */
	if(mpfr_number_p(xstar)) {
	  r = evaluateFaithfulWithCutOffFast(ystar, f, f_diff, xstar, cutoff, prec+10);
	  if (r==2) mpfr_set_d(ystar, 0. , GMP_RNDN); /* under the cutoff */

	  if (!mpfr_number_p(ystar)) {
	    printMessage(1,SOLLYA_MSG_EVALUATION_AT_POINT_GIVES_NAN_EXCLUDING_POINT,"Warning: the evaluation of the given function %b in %v gives NaN.\nThis (possibly maximum) point will be excluded from the infnorm result.\n",f,xstar);
	  }
	  if (mpfr_cmpabs(ystar, max) > 0) { /* evaluates to false when ystar=NaN */
	    mpfr_abs(max, ystar, GMP_RNDU);
            mpfr_set(argmax, xstar, GMP_RNDN);
	    printMessage(3,SOLLYA_MSG_THE_CURRENT_MAXIMUM_IS_A_CERTAIN_VALUE,"Information: current max is %v and is reached at %v\n",max,argmax);
	    mpfr_div_2ui(cutoff, max, 1, GMP_RNDU);
	  }
	}
      }
    }

    mpfr_set(x1,x2,GMP_RNDN); /* exact */
    mpfr_set(y1,y2,GMP_RNDN); /* exact */
    mpfr_set(x2,x3,GMP_RNDN); /* exact */
    mpfr_set(y2,y3,GMP_RNDN); /* exact */
  }

  mpfr_set(result, max, GMP_RNDU);
  if (x_result != NULL) mpfr_set(*x_result, argmax, GMP_RNDN);

  free_memory(f_diff);
  free_memory(f_diff2);

  gmp_randclear(random_state);
  mpfr_clear(perturb);
  mpfr_clear(current_x);
  mpfr_clear(x1);
  mpfr_clear(x2);
  mpfr_clear(x3);
  mpfr_clear(argmax);
  mpfr_clear(step);
  mpfr_clear(y1);
  mpfr_clear(y2);
  mpfr_clear(y3);
  mpfr_clear(ystar);
  mpfr_clear(y1diff);
  mpfr_clear(y3diff);
  mpfr_clear(xstar);
  mpfr_clear(zero_mpfr);
  mpfr_clear(max);
  mpfr_clear(cutoff);
}



rangetype infnorm(node *func, rangetype range, chain *excludes,
		  mp_prec_t prec, mpfr_t diam, FILE *proof) {
  rangetype res;
  sollya_mpfi_t rangeI, resI;
  sollya_mpfi_t *excludeI;
  node *deriv, *numeratorDeriv, *derivNumeratorDeriv, *denominatorDeriv, *derivDenominatorDeriv;
  mpfr_t rangeDiameter, z, ya,yb;
  chain *mightExcludes, *curr, *secondMightExcludes, *initialExcludes;
  int newtonWorked;
  mp_prec_t p, p2;
  infnormTheo *theo;
  int freeInitialExcludes;



  freeInitialExcludes = 1;
  res.a = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
  res.b = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res.a),prec);
  mpfr_init2(*(res.b),prec);

  if ((!mpfr_number_p(*(range.a))) || (!mpfr_number_p(*(range.b)))) {
    printMessage(1,SOLLYA_MSG_DOMAIN_IS_NO_CLOSED_INTERVAL_ON_THE_REALS,"Warning: the bounds of the range an infinity norm is to be computed on are not numbers.\n");
    if (proof != NULL) {
      printMessage(1,SOLLYA_MSG_NO_PROOF_WILL_BE_GENERATED,"Warning: no proof will be generated.\n");
    }
    mpfr_set_d(*(res.a),0.0,GMP_RNDN);
    mpfr_set_inf(*(res.b),1);
    return res;
  }

  if ((mpfr_cmp(*(range.a),*(range.b)) == 0) && (proof == NULL)) {
    evaluateRangeFunctionFast(res, func, NULL, range, prec);
    mpfr_abs(*(res.a),*(res.a),GMP_RNDN);
    mpfr_abs(*(res.b),*(res.b),GMP_RNDN);
    if (mpfr_cmp(*(res.a),*(res.b)) > 0) {
      mpfr_init2(z,prec);
      mpfr_set(z,*(res.b),GMP_RNDN);
      mpfr_set(*(res.b),*(res.a),GMP_RNDN);
      mpfr_set(*(res.a),z,GMP_RNDN);
      mpfr_clear(z);
    }
    return res;
  }

  if (isTrivialInfnormCase(res, func)) {
    if (proof != NULL) {
      printMessage(1,SOLLYA_MSG_INFNORM_RESULT_IS_TRIVIAL,"Warning: the infnorm on the given function is trivially calculable.\n");
      printMessage(1,SOLLYA_MSG_CONTINUATION,"No proof will be generated.\n");
    }
    return res;
  }

  curr = excludes;
  initialExcludes = NULL;
  while (curr != NULL) {
    excludeI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    p = mpfr_get_prec(*(((rangetype *) curr->value)->a));
    p2 = mpfr_get_prec(*(((rangetype *) curr->value)->b));
    if (p2 > p) p = p2;
    if (prec > p) p = prec;
    sollya_mpfi_init2(*excludeI,p);
    sollya_mpfi_interv_fr(*excludeI,*(((rangetype *) curr->value)->a),*(((rangetype *) curr->value)->b));
    initialExcludes = addElement(initialExcludes,(void *) excludeI);
    curr = curr->next;
  }

  sollya_mpfi_init2(rangeI,prec);
  sollya_mpfi_init2(resI,prec);
  mpfr_init2(rangeDiameter,prec);
  mpfr_sub(rangeDiameter,*(range.b),*(range.a),GMP_RNDD);
  mpfr_mul(rangeDiameter,rangeDiameter,diam,GMP_RNDD);
  sollya_mpfi_interv_fr(rangeI,*(range.a),*(range.b));
  deriv = differentiate(func);

  if (getNumeratorDenominator(&numeratorDeriv,&denominatorDeriv,deriv)) {
    printMessage(1,SOLLYA_MSG_DERIVATIVE_IS_QUOTIENT,"Warning: the derivative of the function is a quotient, thus possibly not continuous in the interval.\n");
    printMessage(1,SOLLYA_MSG_CONTINUATION,"Only the zeros of the numerator will be searched and pole detection may fail.\n");
    printMessage(1,SOLLYA_MSG_CONTINUATION,"Be sure that the function is twice continuously differentiable if trusting the infnorm result.\n");

    mpfr_init2(z,prec);
    mpfr_init2(ya,prec);
    mpfr_init2(yb,prec);

    derivDenominatorDeriv = differentiate(denominatorDeriv);

    newtonWorked = newtonMPFR(z, denominatorDeriv, derivDenominatorDeriv, *(range.a), *(range.b), prec);

    if (newtonWorked && mpfr_number_p(z)) {
      evaluate(ya,numeratorDeriv,z,prec);
      evaluate(yb,denominatorDeriv,z,prec);

      mpfr_abs(ya,ya,GMP_RNDN);
      mpfr_abs(yb,yb,GMP_RNDN);

      mpfr_mul_ui(yb,yb,2,GMP_RNDN);

      if (mpfr_cmp(ya,yb) <= 0) {
	printMessage(1,SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_EXTENSIBLE_SINGULARITY,"Warning: the derivative of the function seems to have a extensible singularity in %v.\nThe infnorm result might not be trustful if the derivative cannot actually\nbe extended in this point.\n",z);
      } else {
	printMessage(1,SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_SINGULARITY,"Warning: the derivative of the function seems to have a singularity in %v.\nThe infnorm result is likely to be wrong.\n",z);
      }
    } else {
      evaluate(ya,denominatorDeriv,*(range.a),prec);
      evaluate(yb,denominatorDeriv,*(range.b),prec);

      if (mpfr_sgn(ya) != mpfr_sgn(yb)) {
	printMessage(1,SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_EXTENSIBLE_SINGULARITY,"Warning: the derivative of the function seems to have a (extensible) singularity in the considered interval.\n");
	printMessage(1,SOLLYA_MSG_CONTINUATION,"The infnorm result might be not trustful if the function is not continuously differentiable.\n");
      } else {
	printMessage(2,SOLLYA_MSG_DERIVATIVE_SEEMS_NOT_TO_HAVE_ANY_POLE,"Information: the derivative seems to have no (false) pole in the considered interval.\n");
      }
    }

    mpfr_clear(z);
    mpfr_clear(ya);
    mpfr_clear(yb);
    free_memory(derivDenominatorDeriv);
    free_memory(denominatorDeriv);
  }
  derivNumeratorDeriv = differentiate(numeratorDeriv);
  mightExcludes = NULL;

  if (proof != NULL) {
    theo = (infnormTheo *) safeCalloc(1,sizeof(infnormTheo));
  } else {
    theo = NULL;
  }

  printMessage(3,SOLLYA_MSG_INVOKING_INFNORM_SUBFUNCTION,"Information: invoking the interval infnorm subfunction.\n");

  infnormI(resI,func,deriv,numeratorDeriv,derivNumeratorDeriv,rangeI,
	   prec,rangeDiameter,initialExcludes,&mightExcludes,theo);

  printMessage(3,SOLLYA_MSG_INFNORM_SUBFUNCTION_HAS_FINISHED,"Information: interval infnorm subfunction has finished.\n");

  secondMightExcludes = NULL;

  if (mightExcludes != NULL) {
    mightExcludes = concatChains(mightExcludes,initialExcludes);
    freeInitialExcludes = 0;

    if (theo != NULL) freeInfnormTheo(theo);
    if (proof != NULL) {
      theo = (infnormTheo *) safeCalloc(1,sizeof(infnormTheo));
    } else {
      theo = NULL;
    }

    printMessage(3,SOLLYA_MSG_INVOKING_INFNORM_SUBFUNCTION,"Information: invoking the interval infnorm subfunction on additional excludes.\n");

    infnormI(resI,func,deriv,numeratorDeriv,derivNumeratorDeriv,rangeI,
	     2*prec,rangeDiameter,mightExcludes,&secondMightExcludes,theo);

    printMessage(3,SOLLYA_MSG_INFNORM_SUBFUNCTION_HAS_FINISHED,"Information: interval infnorm subfunction on additional excludes has finished.\n");

  }

  if (proof != NULL) {
    printMessage(2,SOLLYA_MSG_STARTING_TO_WRITE_THE_PROOF,"Information: started writing the proof.\n");
    fprintInfnormTheo(proof,theo,1);
    printMessage(2,SOLLYA_MSG_THE_PROOF_HAS_BEEN_WRITTEN,"Information: proof written.\n");
  }

  if (theo != NULL) freeInfnormTheo(theo);
  freeChain(mightExcludes,freeMpfiPtr);
  freeChain(secondMightExcludes,freeMpfiPtr);
  if (freeInitialExcludes) freeChain(initialExcludes,freeMpfiPtr);
  sollya_mpfi_get_left(*(res.a),resI);
  sollya_mpfi_get_right(*(res.b),resI);
  free_memory(deriv);
  free_memory(numeratorDeriv);
  free_memory(derivNumeratorDeriv);
  sollya_mpfi_clear(rangeI);
  sollya_mpfi_clear(resI);
  mpfr_clear(rangeDiameter);
  return res;
}


void evaluateRangeFunctionFast(rangetype yrange, node *func, node *deriv, rangetype xrange, mp_prec_t prec) {
  sollya_mpfi_t x, y;
  chain *tempChain;
  mp_prec_t p, p2;

  p = prec;
  p2 = mpfr_get_prec(*(xrange.a));
  if (p2 > p) p = p2;
  p2 = mpfr_get_prec(*(xrange.b));
  if (p2 > p) p = p2;

  sollya_mpfi_init2(x,p);
  sollya_mpfi_init2(y,prec);
  sollya_mpfi_interv_fr(x,*(xrange.a),*(xrange.b));

  tempChain = evaluateITaylorStart(y, func, deriv, x, prec, taylorrecursions, NULL, 1,0,0,0,NULL,NULL);

  sollya_mpfi_get_left(*(yrange.a),y);
  sollya_mpfi_get_right(*(yrange.b),y);

  freeChain(tempChain,freeMpfiPtr);
  sollya_mpfi_clear(x);
  sollya_mpfi_clear(y);
}

void evaluateInterval(sollya_mpfi_t y, node *func, node *deriv, sollya_mpfi_t x) {
  mp_prec_t prec;

  prec = sollya_mpfi_get_prec(y);

  /* We need more precision in the first steps to get the precision in the end. */
  prec += 10;

  evaluateITaylorStart(y, func, deriv, x, prec, taylorrecursions, NULL, 1,0,0,0,NULL,NULL);
}

void evaluateIntervalPlain(sollya_mpfi_t y, node *func, sollya_mpfi_t x) {
  evaluateITaylorStart(y, func, NULL, x, sollya_mpfi_get_prec(y), taylorrecursions, NULL, 1,0,0,0,NULL,NULL);
}

static inline void evaluateIntervalInternalFast(sollya_mpfi_t y, node *func, node *deriv, sollya_mpfi_t x, int adaptPrecision, int noLazyHooks, mp_exp_t *cutoff, int *lazyHookUsed) {
  mp_prec_t prec;

  prec = sollya_mpfi_get_prec(y);

  /* We need more precision in the first steps to get the precision in the end. */
  prec += 10;

  evaluateITaylorStart(y, func, deriv, x, prec, taylorrecursions, NULL, 1,1,adaptPrecision,noLazyHooks,cutoff,lazyHookUsed);
}

void evaluateConstantExpressionToInterval(sollya_mpfi_t y, node *func) {
  sollya_mpfi_t x;

  if (!isConstant(func)) {
    printMessage(1,SOLLYA_MSG_THE_EXPRESSION_IS_NOT_CONSTANT,"Warning: the given expression is not constant. Evaluating it at 1.\n");
  }

  sollya_mpfi_init2(x,12);
  sollya_mpfi_set_si(x,1);

  evaluateInterval(y, func, NULL, x);

  sollya_mpfi_clear(x);
}

void evaluateConstantExpressionToSharpInterval(sollya_mpfi_t y, node *func) {
  sollya_mpfi_t x, myY, myY2;
  mpfr_t X, Y, Y2;
  mp_prec_t prec;

  if (!isConstant(func)) {
    printMessage(1,SOLLYA_MSG_THE_EXPRESSION_IS_NOT_CONSTANT,"Warning: the given expression is not constant. Evaluating it at 1.\n");
  }

  sollya_mpfi_init2(x,12);
  sollya_mpfi_set_si(x,1);
  mpfr_init2(X, 12);
  mpfr_set_si(X, 1,GMP_RNDN);

  prec = sollya_mpfi_get_prec(y);
  mpfr_init2(Y, prec + 5);

  if (evaluateFaithful(Y, func, X, prec + 10)) {
    mpfr_init2(Y2, mpfr_get_prec(Y));
    sollya_mpfi_init2(myY, prec + 1);
    sollya_mpfi_init2(myY2, prec + 1);
    mpfr_set(Y2, Y, GMP_RNDN); /* exact */
    mpfr_nextabove(Y);
    mpfr_nextbelow(Y2);
    sollya_mpfi_interv_fr(myY, Y2, Y);
    evaluateInterval(myY2, func, NULL, x);
    sollya_mpfi_intersect(myY, myY, myY2);
    sollya_mpfi_set(y,myY);
    sollya_mpfi_clear(myY);
    sollya_mpfi_clear(myY2);
    mpfr_clear(Y2);
  } else {
    sollya_mpfi_init2(myY, prec * 256);
    evaluateInterval(myY, func, NULL, x);
    sollya_mpfi_set(y,myY);
    sollya_mpfi_clear(myY);
  }

  sollya_mpfi_clear(x);
  mpfr_clear(X);
  mpfr_clear(Y);
}

void evaluateRangeFunction(rangetype yrange, node *func, rangetype xrange, mp_prec_t prec) {
  node *deriv, *temp2, *numerator, *denominator, *f;
  rangetype myrange;
  node *myderiv;

  if (containsNonDifferentiableSubfunctions(func)) {
    evaluateRangeFunctionFast(yrange,func,NULL,xrange,prec);
    return;
  }

  if (mpfr_equal_p(*(xrange.a),*(xrange.b))) {
    evaluateRangeFunctionFast(yrange,func,NULL,xrange,prec);
    return;
  }

  if ((!(mpfr_nan_p(*(xrange.a)) ||
	 mpfr_inf_p(*(xrange.a)) ||
	 mpfr_zero_p(*(xrange.a)) ||
	 mpfr_nan_p(*(xrange.b)) ||
	 mpfr_inf_p(*(xrange.b)) ||
	 mpfr_zero_p(*(xrange.b)))) &&
      ((mpfr_get_exp(*(xrange.a)) == mpfr_get_exp(*(xrange.b))) &&
       (mpfr_sgn(*(xrange.a)) == mpfr_sgn(*(xrange.b)))
       )) {
    evaluateRangeFunctionFast(yrange,func,NULL,xrange,prec);
    return;
  }

  deriv = differentiate(func);
  temp2 = copyTree(func);

  f = NULL;

  if (getNumeratorDenominator(&numerator,&denominator,temp2)) {
    if (isSyntacticallyEqual(numerator, denominator)) {
      if (!isConstant(numerator)) {
        myderiv = differentiate(numerator);
        myrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        myrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(myrange.a),mpfr_get_prec(*(yrange.a)));
        mpfr_init2(*(myrange.b),mpfr_get_prec(*(yrange.b)));
        evaluateRangeFunctionFast(myrange,numerator,myderiv,xrange,prec);
        free_memory(myderiv);
        if (mpfr_sgn(*(myrange.a)) * mpfr_sgn(*(myrange.b)) == 1) {
          mpfr_clear(*(myrange.a));
          mpfr_clear(*(myrange.b));
          safeFree(myrange.a);
          safeFree(myrange.b);
          mpfr_set_d(*(yrange.a),1.0,GMP_RNDD);
          mpfr_set_d(*(yrange.b),1.0,GMP_RNDU);
          free_memory(numerator);
          free_memory(denominator);
          free_memory(deriv);
          free_memory(temp2);
          return;
        } else {
          f = copyTree(temp2);
        }
      } else {
        myderiv = differentiate(numerator);
        myrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        myrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(myrange.a),mpfr_get_prec(*(yrange.a)));
        mpfr_init2(*(myrange.b),mpfr_get_prec(*(yrange.b)));
        evaluateRangeFunctionFast(myrange,numerator,myderiv,xrange,prec);
        free_memory(myderiv);
        if (mpfr_sgn(*(myrange.a)) * mpfr_sgn(*(myrange.b)) == 1) {
          mpfr_clear(*(myrange.a));
          mpfr_clear(*(myrange.b));
          safeFree(myrange.a);
          safeFree(myrange.b);
          mpfr_set_d(*(yrange.a),1.0,GMP_RNDD);
          mpfr_set_d(*(yrange.b),1.0,GMP_RNDU);
          free_memory(numerator);
          free_memory(denominator);
          free_memory(deriv);
          free_memory(temp2);
          return;
        } else {
          if (mpfr_zero_p(*(myrange.a)) && mpfr_zero_p(*(myrange.b))) {
            mpfr_clear(*(myrange.a));
            mpfr_clear(*(myrange.b));
            safeFree(myrange.a);
            safeFree(myrange.b);
            mpfr_set_nan(*(yrange.a));
            mpfr_set_nan(*(yrange.b));
            free_memory(numerator);
            free_memory(denominator);
            free_memory(deriv);
            free_memory(temp2);
            return;
          } else {
            mpfr_clear(*(myrange.a));
            mpfr_clear(*(myrange.b));
            safeFree(myrange.a);
            safeFree(myrange.b);
            f = copyTree(temp2);
          }
        }
      }
    } else {
      f = copyTree(temp2);
    }
    free_memory(numerator);
    free_memory(denominator);
  } else {
    free_memory(numerator);
    f = copyTree(temp2);
  }

  evaluateRangeFunctionFast(yrange,f,deriv,xrange,prec);
  free_memory(deriv);
  free_memory(temp2);
  free_memory(f);
}



chain* findZerosFunction(node *func, rangetype range, mp_prec_t prec, mpfr_t diam) {
  sollya_mpfi_t rangeI;
  node *deriv, *numerator, *denominator;
  mpfr_t rangeDiameter, diamJoin;
  chain *zerosI, *zeros, *curr;
  rangetype *tempRange;
  chain *tempChain, *tempChain2, *tempChain3;

  sollya_mpfi_init2(rangeI,prec);
  mpfr_init2(rangeDiameter,prec);
  mpfr_sub(rangeDiameter,*(range.b),*(range.a),GMP_RNDD);
  mpfr_mul(rangeDiameter,rangeDiameter,diam,GMP_RNDD);
  sollya_mpfi_interv_fr(rangeI,*(range.a),*(range.b));

  if (getNumeratorDenominator(&numerator,&denominator,func)) free_memory(denominator);

  deriv = differentiate(numerator);

  tempChain = findZeros(numerator,deriv,rangeI,prec,rangeDiameter,NULL);
  mpfr_init2(diamJoin,prec);
  mpfr_mul_2ui(diamJoin,diam,3,GMP_RNDN);
  tempChain2 = joinAdjacentIntervals(tempChain,diamJoin);
  tempChain3 = copyChain(tempChain2,copyMpfiPtr);
  mpfr_mul_2ui(diamJoin,diamJoin,2,GMP_RNDN);
  zerosI = joinAdjacentIntervals(tempChain3,diamJoin);
  mpfr_clear(diamJoin);

  zeros = NULL;
  curr = zerosI;
  while (curr != NULL) {
    tempRange = (rangetype *) safeMalloc(sizeof(rangetype));
    tempRange->a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    tempRange->b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*(tempRange->a),prec);
    mpfr_init2(*(tempRange->b),prec);
    sollya_mpfi_get_left(*(tempRange->a),*((sollya_mpfi_t *) (curr->value)));
    sollya_mpfi_get_right(*(tempRange->b),*((sollya_mpfi_t *) (curr->value)));
    zeros = addElement(zeros,tempRange);
    curr = curr->next;
  }

  freeChain(zerosI,freeMpfiPtr);
  freeChain(tempChain,freeMpfiPtr);
  freeChain(tempChain2,freeMpfiPtr);
  freeChain(tempChain3,freeMpfiPtr);
  free_memory(numerator);
  free_memory(deriv);
  sollya_mpfi_clear(rangeI);
  mpfr_clear(rangeDiameter);
  return zeros;
}


int checkInfnormI(node *func, node *deriv, sollya_mpfi_t infnormval, sollya_mpfi_t range, mpfr_t diam, mp_prec_t prec) {
  sollya_mpfi_t evaluateOnRange, rangeLeft, rangeRight;
  chain *tempChain;
  mpfr_t l,m,r, diamRange;
  int resultLeft, resultRight;

  sollya_mpfi_init2(evaluateOnRange,prec);

  tempChain = evaluateITaylorStart(evaluateOnRange, func, deriv, range, prec, taylorrecursions, NULL, 1,0,0,0,NULL,NULL);

  freeChain(tempChain,freeMpfiPtr);

  if (sollya_mpfi_is_inside(evaluateOnRange, infnormval)) {
    /* Simple end case: the interval evaluation is contained in the given interval for infnorm */
    sollya_mpfi_clear(evaluateOnRange);
    return 1;
  }

  mpfr_init2(diamRange,prec);
  sollya_mpfi_diam_abs(diamRange,range);

  if (mpfr_cmp(diamRange,diam) <= 0) {
    /* Simple end case: the range to test is already smaller than diam but we could not check */
    printMessage(2,SOLLYA_MSG_COULD_NOT_CHECK_INFNORM_ON_A_CERTAIN_INTERVAL,"Information: could not check the infinity norm on the domain\n%w\nThe function evaluates here to\n%w.\n",range,evaluateOnRange);
    sollya_mpfi_clear(evaluateOnRange);
    mpfr_clear(diamRange);
    return 0;
  }

  mpfr_init2(l,prec);
  mpfr_init2(m,prec);
  mpfr_init2(r,prec);
  sollya_mpfi_init2(rangeLeft,prec);
  sollya_mpfi_init2(rangeRight,prec);

  sollya_mpfi_get_left(l,range);
  sollya_mpfi_mid(m,range);
  sollya_mpfi_get_right(r,range);

  sollya_mpfi_interv_fr(rangeLeft,l,m);
  sollya_mpfi_interv_fr(rangeRight,m,r);

  /* Recurse on half the range */

  resultLeft = 0;
  resultRight = 0;

  resultLeft = checkInfnormI(func, deriv, infnormval, rangeLeft, diam, prec);
  if (resultLeft) resultRight = checkInfnormI(func, deriv, infnormval, rangeRight, diam, prec);

  sollya_mpfi_clear(rangeRight);
  sollya_mpfi_clear(rangeLeft);
  mpfr_clear(r);
  mpfr_clear(m);
  mpfr_clear(l);
  sollya_mpfi_clear(evaluateOnRange);
  mpfr_clear(diamRange);

  return (resultLeft && resultRight);
}


int checkInfnorm(node *func, rangetype range, mpfr_t infnormval, mpfr_t diam, mp_prec_t prec) {
  node *deriv;
  sollya_mpfi_t rangeI, infnormvalI;
  mpfr_t rangeDiameter, tempLeft, tempRight;
  int result;

  sollya_mpfi_init2(rangeI,prec);
  sollya_mpfi_init2(infnormvalI,prec);
  mpfr_init2(rangeDiameter,prec);
  mpfr_init2(tempLeft,prec);
  mpfr_init2(tempRight,prec);

  mpfr_sub(rangeDiameter,*(range.b),*(range.a),GMP_RNDD);
  mpfr_mul(rangeDiameter,rangeDiameter,diam,GMP_RNDD);
  sollya_mpfi_interv_fr(rangeI,*(range.a),*(range.b));
  mpfr_abs(tempRight,infnormval,GMP_RNDU);
  mpfr_neg(tempLeft,tempRight,GMP_RNDD);
  sollya_mpfi_interv_fr(infnormvalI,tempLeft,tempRight);
  deriv = differentiate(func);

  result = checkInfnormI(func, deriv, infnormvalI, rangeI, rangeDiameter, prec);

  free_memory(deriv);
  mpfr_clear(tempLeft);
  mpfr_clear(tempRight);
  mpfr_clear(rangeDiameter);
  sollya_mpfi_clear(infnormvalI);
  sollya_mpfi_clear(rangeI);

  return result;
}


void evaluateConstantWithErrorEstimate(mpfr_t res, mpfr_t err, node *func, mpfr_t x, mp_prec_t prec) {
  rangetype xrange, yrange;
  mpfr_t temp;

  xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));

  mpfr_init2(*(xrange.a), prec);
  mpfr_init2(*(xrange.b), prec);
  mpfr_init2(*(yrange.a), prec);
  mpfr_init2(*(yrange.b), prec);
  mpfr_init2(temp,prec + 10);

  mpfr_set(*(xrange.a),x,GMP_RNDD);
  mpfr_set(*(xrange.b),x,GMP_RNDU);

  evaluateRangeFunction(yrange, func, xrange, prec);

  mpfr_add(temp,*(yrange.a),*(yrange.b),GMP_RNDN);
  mpfr_div_2ui(temp,temp,1,GMP_RNDN);
  mpfr_set(res,temp,GMP_RNDN);

  if (mpfr_zero_p(res)) {
    if (mpfr_zero_p(*(yrange.a)) && mpfr_zero_p(*(yrange.b))) {
      mpfr_set_d(err,0.0,GMP_RNDN);
    } else {
      mpfr_set_d(temp,1.0,GMP_RNDN);
      mpfr_div(temp,temp,res,GMP_RNDN);
      mpfr_set(err,temp,GMP_RNDU);
    }
  } else {
    mpfr_abs(*(yrange.a),*(yrange.a),GMP_RNDN);
    mpfr_abs(*(yrange.b),*(yrange.b),GMP_RNDN);
    if (mpfr_cmp(*(yrange.b),*(yrange.a)) > 0) {
      mpfr_set(*(yrange.a),*(yrange.b),GMP_RNDN);
    }
    mpfr_abs(temp,temp,GMP_RNDN);
    mpfr_sub(*(yrange.a),*(yrange.a),temp,GMP_RNDU);
    mpfr_div(err,*(yrange.a),temp,GMP_RNDU);
  }


  mpfr_clear(*(xrange.a));
  mpfr_clear(*(xrange.b));
  mpfr_clear(*(yrange.a));
  mpfr_clear(*(yrange.b));
  mpfr_init2(temp,prec);
  safeFree(xrange.a);
  safeFree(xrange.b);
  safeFree(yrange.a);
  safeFree(yrange.b);
}

chain* findZerosByNewton(node *func, mpfr_t a, mpfr_t b, mp_prec_t prec) {
  node *deriv;
  mpfr_t resNewtonStep, ap, bp, step, yAp, yBp;
  int newtonOkay;
  mpfr_t *newZero;
  chain *fpZeros;

  fpZeros = NULL;
  deriv = differentiate(func);
  mpfr_init2(resNewtonStep,prec);
  mpfr_init2(ap,prec);
  mpfr_init2(bp,prec);
  mpfr_init2(step,prec);
  mpfr_init2(yAp,prec);
  mpfr_init2(yBp,prec);

  mpfr_sub(step,b,a,GMP_RNDU);
  if (mpfr_zero_p(step)) {
    evaluate(resNewtonStep,func,a,prec);
    newZero = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*newZero,prec);
    mpfr_set(*newZero,resNewtonStep,GMP_RNDN);
    fpZeros = addElement(fpZeros,newZero);
  } else {
    mpfr_div_ui(step,step,defaultpoints,GMP_RNDU);
    mpfr_set(ap,a,GMP_RNDD);
    while (mpfr_cmp(ap,b) < 0) {
      mpfr_add(bp,ap,step,GMP_RNDN);
      sollya_mpfr_min(bp,bp,b,GMP_RNDU);

      newtonOkay = newtonMPFR(resNewtonStep, func, deriv, ap, bp, prec);

      if (newtonOkay) {
	newZero = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*newZero,prec);
	mpfr_set(*newZero,resNewtonStep,GMP_RNDN);
	fpZeros = addElement(fpZeros,newZero);
      } else {
	evaluateFaithful(yAp, func, ap, prec);
	evaluateFaithful(yBp, func, bp, prec);
	if (mpfr_number_p(yAp) && mpfr_number_p(yBp) &&
	    (mpfr_sgn(yAp) != mpfr_sgn(yBp))) {
	  newZero = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*newZero,prec);
	  mpfr_set(*newZero,ap,GMP_RNDN);
	  mpfr_add(*newZero,*newZero,bp,GMP_RNDN);
	  mpfr_div_2ui(*newZero,*newZero,1,GMP_RNDN);
	  fpZeros = addElement(fpZeros,newZero);
	}
      }

      mpfr_set(ap,bp,GMP_RNDN);
    }
  }
  mpfr_clear(step);
  mpfr_clear(bp);
  mpfr_clear(ap);
  mpfr_clear(yAp);
  mpfr_clear(yBp);
  mpfr_clear(resNewtonStep);
  free_memory(deriv);
  return fpZeros;
}



chain* fpFindZerosFunction(node *func, rangetype range, mp_prec_t prec) {
  mpfr_t diam;
  chain *intervalZeros, *fpZeros, *temp, *fpZerosOnInterval, *fpZeros2, *curr;
  mpfr_t *newZero;
  mpfr_t before, after, yBefore, yAfter, y, compare;
  int addToList, removedFromList;
  unsigned int oldDefaultPoints;


  oldDefaultPoints = defaultpoints;
  defaultpoints = defaultpoints >> 4;

  mpfr_init2(diam,prec+50);
  mpfr_set_d(diam,DEFAULTDIAM2,GMP_RNDN);

  intervalZeros = findZerosFunction(func, range, prec, diam);

  mpfr_clear(diam);

  fpZeros = NULL;

  while (intervalZeros != NULL) {
    fpZerosOnInterval = findZerosByNewton(func,
					  *(((rangetype *) (intervalZeros->value))->a),
					  *(((rangetype *) (intervalZeros->value))->b),
					  4*prec);
    fpZeros = concatChains(fpZeros, fpZerosOnInterval);
    mpfr_clear(*(((rangetype *) (intervalZeros->value))->a));
    mpfr_clear(*(((rangetype *) (intervalZeros->value))->b));
    safeFree(((rangetype *) (intervalZeros->value))->a);
    safeFree(((rangetype *) (intervalZeros->value))->b);
    safeFree(intervalZeros->value);
    temp = intervalZeros->next;
    safeFree(intervalZeros);
    intervalZeros = temp;
  }


  mpfr_init2(compare,prec);
  mpfr_set_d(compare,1.0,GMP_RNDN);
  mpfr_div_ui(compare,compare,prec,GMP_RNDN);
  fpZeros2 = NULL;
  curr = fpZeros;
  while (curr != NULL) {
    while ((curr->next != NULL) && (mpfr_cmp(*((mpfr_t *) (curr->value)),*((mpfr_t *) (curr->next->value))) == 0))
      curr = curr->next;
    newZero = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*newZero,prec);
    mpfr_set(*newZero,*((mpfr_t *) (curr->value)),GMP_RNDN);
    fpZeros2 = addElement(fpZeros2,newZero);
    if ((!mpfr_zero_p(*newZero)) && (mpfr_cmpabs(*newZero,compare) <= 0)) {
      newZero = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*newZero,prec);
      mpfr_set_d(*newZero,0.0,GMP_RNDN);
      fpZeros2 = addElement(fpZeros2,newZero);
    }
    curr = curr->next;
  }
  mpfr_clear(compare);

  while (fpZeros != NULL) {
    mpfr_clear(*((mpfr_t*) (fpZeros->value)));
    safeFree((fpZeros->value));
    temp = fpZeros->next;
    safeFree(fpZeros);
    fpZeros = temp;
  }

  mpfr_init2(before,prec);
  mpfr_init2(after,prec);
  mpfr_init2(yAfter,prec);
  mpfr_init2(yBefore,prec);
  mpfr_init2(y,prec);

  removedFromList = 0;
  fpZeros = NULL;
  while (fpZeros2 != NULL) {

    addToList = 0;

    evaluateFaithful(y, func, *((mpfr_t *) (fpZeros2->value)), prec);

    if (mpfr_zero_p(y) || (!mpfr_number_p(y)) || mpfr_zero_p(*((mpfr_t *) (fpZeros2->value)))) {
      addToList = 1;
    } else {

      mpfr_set(before,*((mpfr_t *) (fpZeros2->value)),GMP_RNDN);
      mpfr_set(after,*((mpfr_t *) (fpZeros2->value)),GMP_RNDN);
      mpfr_nextabove(after);
      mpfr_nextbelow(before);

      evaluateFaithful(yAfter, func, after, prec);
      evaluateFaithful(yBefore, func, before, prec);

      if ((!mpfr_number_p(yAfter)) || (!mpfr_number_p(yBefore))) {
	addToList = 1;
      } else {
	if (mpfr_sgn(yAfter) != mpfr_sgn(yBefore)) {
	  addToList = 1;
	} else {
	  if (mpfr_number_p(y)) {
	    if (mpfr_sgn(y) != mpfr_sgn(yAfter)) {
	      addToList = 1;
	    } else {
	      removedFromList = 1;
	      printMessage(2,SOLLYA_MSG_REMOVING_A_POSSIBLE_ZERO_AT_SOME_POINT,"Information: removing possible zero in %v.\nInformation: removing because all signs are equal.\n",*((mpfr_t *) (fpZeros2->value)));
	    }
	  } else {
	    removedFromList = 1;
	    printMessage(2,SOLLYA_MSG_REMOVING_A_POSSIBLE_ZERO_AT_SOME_POINT,"Information: removing possible zero in %v.\nInformation: removing because predecessor and successor signs are equal.\n",*((mpfr_t *) (fpZeros2->value)));
	  }
	}
      }
    }

    if (addToList) {
      fpZeros = addElement(fpZeros,fpZeros2->value);
    }
    temp = fpZeros2->next;
    safeFree(fpZeros2);
    fpZeros2 = temp;
  }

  if (removedFromList) {
    printMessage(1,SOLLYA_MSG_ZERO_FILTER_HAS_REMOVED_AT_LEAST_ONE_ZERO,"Warning: actual zero filter has removed at least one possible zero of higher order.\n");
  }

  mpfr_clear(before);
  mpfr_clear(after);
  mpfr_clear(yBefore);
  mpfr_clear(yAfter);
  mpfr_clear(y);

  sortChain(fpZeros,  cmpMpfrPtr);

  fpZeros2 = NULL;
  curr = fpZeros;
  while (curr != NULL) {
    while ((curr->next != NULL) && (mpfr_cmp(*((mpfr_t *) (curr->value)),*((mpfr_t *) (curr->next->value))) == 0))
      curr = curr->next;
    newZero = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*newZero,mpfr_get_prec(*((mpfr_t *) (curr->value))));
    mpfr_set(*newZero,*((mpfr_t *) (curr->value)),GMP_RNDN);
    fpZeros2 = addElement(fpZeros2,newZero);
    curr = curr->next;
  }

  while (fpZeros != NULL) {
    mpfr_clear(*((mpfr_t*) (fpZeros->value)));
    safeFree((fpZeros->value));
    temp = fpZeros->next;
    safeFree(fpZeros);
    fpZeros = temp;
  }

  defaultpoints = oldDefaultPoints;

  return fpZeros2;
}


chain *uncertifiedZeroDenominators(node *tree, mpfr_t a, mpfr_t b, mp_prec_t prec) {
  chain *leftPoles, *rightPoles, *newZeros;
  rangetype range;

  if (tree == NULL) return NULL;
  if (tree->nodeType == MEMREF) return uncertifiedZeroDenominators(getMemRefChild(tree), a, b, prec);
  switch (tree->nodeType) {
  case VARIABLE:
    return NULL;
    break;
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    return NULL;
    break;
  case ADD:
    leftPoles = uncertifiedZeroDenominators(tree->child1,a,b,prec);
    rightPoles = uncertifiedZeroDenominators(tree->child2,a,b,prec);
    return concatChains(leftPoles,rightPoles);
    break;
  case SUB:
    leftPoles = uncertifiedZeroDenominators(tree->child1,a,b,prec);
    rightPoles = uncertifiedZeroDenominators(tree->child2,a,b,prec);
    return concatChains(leftPoles,rightPoles);
    break;
  case MUL:
    leftPoles = uncertifiedZeroDenominators(tree->child1,a,b,prec);
    rightPoles = uncertifiedZeroDenominators(tree->child2,a,b,prec);
    return concatChains(leftPoles,rightPoles);
    break;
  case DIV:
    leftPoles = uncertifiedZeroDenominators(tree->child1,a,b,prec);
    rightPoles = uncertifiedZeroDenominators(tree->child2,a,b,prec);
    range.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    range.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*(range.a),prec);
    mpfr_init2(*(range.b),prec);
    mpfr_set(*(range.a),a,GMP_RNDD);
    mpfr_set(*(range.b),b,GMP_RNDU);
    newZeros = fpFindZerosFunction(tree->child2, range, prec);
    mpfr_clear(*(range.a));
    mpfr_clear(*(range.b));
    safeFree(range.a);
    safeFree(range.b);
    leftPoles = concatChains(leftPoles,rightPoles);
    return concatChains(leftPoles,newZeros);
    break;
  case NEG:
    return uncertifiedZeroDenominators(tree->child1,a,b,prec);
    break;
  case UNARY_BASE_FUNC:
    return uncertifiedZeroDenominators(tree->child1,a,b,prec);
    break;
  case POW:
    leftPoles = uncertifiedZeroDenominators(tree->child1,a,b,prec);
    rightPoles = uncertifiedZeroDenominators(tree->child2,a,b,prec);
    return concatChains(leftPoles,rightPoles);
    break;
  case LIBRARYFUNCTION:
    return uncertifiedZeroDenominators(tree->child1,a,b,prec);
    break;
  case PROCEDUREFUNCTION:
    return uncertifiedZeroDenominators(tree->child1,a,b,prec);
    break;
  default:
    sollyaFprintf(stderr,"Error: uncertifiedZeroDenominators: unknown identifier (%d) in the tree\n",tree->nodeType);
    exit(1);
  }
  return NULL;
}


int isEvaluable(node *func, mpfr_t x, mpfr_t *y, mp_prec_t prec) {
  mpfr_t val;
  rangetype xrange, yrange;

  mpfr_init2(val,prec);
  evaluate(val,func,x,prec);
  if (mpfr_number_p(val)) {
    if (y != NULL) {
      mpfr_set(*y,val,GMP_RNDN);
    }
    mpfr_clear(val);
    return ISFLOATINGPOINTEVALUABLE;
  }

  xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));

  mpfr_init2(*(xrange.a),prec);
  mpfr_init2(*(xrange.b),prec);
  mpfr_init2(*(yrange.a),prec);
  mpfr_init2(*(yrange.b),prec);

  mpfr_set(*(xrange.a),x,GMP_RNDD);
  mpfr_set(*(xrange.b),x,GMP_RNDU);

  evaluateRangeFunction(yrange, func, xrange, prec);

  if (mpfr_number_p(*(yrange.a)) && mpfr_number_p(*(yrange.b))) {
    mpfr_add(val,*(yrange.a),*(yrange.b),GMP_RNDN);
    mpfr_div_2ui(val,val,1,GMP_RNDN);
    if (!mpfr_number_p(val)) {
      mpfr_clear(val);
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
      return ISNOTEVALUABLE;
    }
    if (y != NULL) {
      mpfr_set(*y,val,GMP_RNDN);
    }
    mpfr_clear(val);
    mpfr_clear(*(xrange.a));
    mpfr_clear(*(xrange.b));
    mpfr_clear(*(yrange.a));
    mpfr_clear(*(yrange.b));
    safeFree(xrange.a);
    safeFree(xrange.b);
    safeFree(yrange.a);
    safeFree(yrange.b);
    return ISHOPITALEVALUABLE;
  }

  mpfr_clear(val);
  mpfr_clear(*(xrange.a));
  mpfr_clear(*(xrange.b));
  mpfr_clear(*(yrange.a));
  mpfr_clear(*(yrange.b));
  safeFree(xrange.a);
  safeFree(xrange.b);
  safeFree(yrange.a);
  safeFree(yrange.b);
  return ISNOTEVALUABLE;
}

int evaluateFaithful(mpfr_t result, node *tree, mpfr_t x, mp_prec_t prec) {
  mp_prec_t startPrec, p;
  mpfr_t cutoff;
  int res;

  p = mpfr_get_prec(result);
  startPrec = p + 10;
  if (prec > startPrec) startPrec = prec;

  mpfr_init2(cutoff,startPrec);
  mpfr_set_si(cutoff,0,GMP_RNDN);

  res = evaluateFaithfulWithCutOffFast(result, tree, NULL, x, cutoff, startPrec);
  if (res==3) res=0;

  mpfr_clear(cutoff);

  if (!res) {
    printMessage(4,SOLLYA_MSG_FAITHFUL_EVALUATION_RETURNS_NAN,"Warning: evaluateFaithful returned NaN.\n");
    mpfr_set_nan(result);
  }

  return res;

}

int determineHeuristicTaylorRecursions(node *func) {
  int highestDegree, sizeOfFunc, sizeOfCurrDeriv, i;
  node *temp, *temp2;

  highestDegree = highestDegreeOfPolynomialSubexpression(func);

  sizeOfFunc = treeSize(func);

  temp = differentiate(func);
  sizeOfCurrDeriv = treeSize(temp);
  i = -1;

  while ((highestDegree >= 0) && (((double) sizeOfCurrDeriv) <= ((double) 4) * ((double) sizeOfFunc))) {
    temp2 = differentiate(temp);
    free_memory(temp);
    temp = temp2;
    sizeOfCurrDeriv = treeSize(temp);
    i++;
    highestDegree--;
  }

  free_memory(temp);

  return i < 0 ? 0 : i;
}



int accurateInfnorm(mpfr_t result, node *func, rangetype range, chain *excludes, mp_prec_t startPrec) {
  rangetype res;
  sollya_mpfi_t rangeI, resI;
  sollya_mpfi_t *excludeI;
  node *deriv, *numeratorDeriv, *derivNumeratorDeriv, *denominatorDeriv, *derivDenominatorDeriv;
  mpfr_t rangeDiameter, z, ya,yb;
  chain *curr, *initialExcludes;
  int newtonWorked;
  mp_prec_t p, p2, prec;
  mpfr_t startDiam, currDiameter, resultUp, resultDown, stopDiameter;
  int okay, oldtaylorrecursions, t;


  prec = startPrec;
  p = mpfr_get_prec(result);

  mpfr_init2(resultUp,p);
  mpfr_init2(resultDown,p);

  if (p > prec) {
    prec = p;
    printMessage(1,SOLLYA_MSG_INTERMEDIATE_PRECISION_HAS_BEEN_INCREASED,"Warning: starting intermediate precision increased to %d bits.\n",prec);
  }

  res.a = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
  res.b = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res.a),prec);
  mpfr_init2(*(res.b),prec);

  if (isTrivialInfnormCase(res, func)) {
    printMessage(2,SOLLYA_MSG_INFNORM_RESULT_IS_TRIVIAL,"Information: the infnorm on the given function is trivially calculable.\n");
    mpfr_set(result,*(res.a),GMP_RNDU);
    mpfr_clear(*(res.a));
    mpfr_clear(*(res.b));
    safeFree(res.a);
    safeFree(res.b);
    return 1;
  }

  oldtaylorrecursions = taylorrecursions;
  t = determineHeuristicTaylorRecursions(func);
  if ((t > oldtaylorrecursions) && (t < ((oldtaylorrecursions + 1) * 2))) {
    taylorrecursions = t;
    printMessage(3,SOLLYA_MSG_TAYLOR_RECURSION_TEMPORARILY_SET_TO_A_VALUE,"Information: the number of Taylor recursions has temporarily been set to %d.\n",taylorrecursions);
  }

  curr = excludes;
  initialExcludes = NULL;
  while (curr != NULL) {
    excludeI = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    p = mpfr_get_prec(*(((rangetype *) curr->value)->a));
    p2 = mpfr_get_prec(*(((rangetype *) curr->value)->b));
    if (p2 > p) p = p2;
    if (prec > p) p = prec;
    sollya_mpfi_init2(*excludeI,p);
    sollya_mpfi_interv_fr(*excludeI,*(((rangetype *) curr->value)->a),*(((rangetype *) curr->value)->b));
    initialExcludes = addElement(initialExcludes,(void *) excludeI);
    curr = curr->next;
  }

  sollya_mpfi_init2(rangeI,prec);
  sollya_mpfi_init2(resI,prec);
  mpfr_init2(rangeDiameter,prec);
  mpfr_sub(rangeDiameter,*(range.b),*(range.a),GMP_RNDD);

  mpfr_init2(startDiam,prec);
  mpfr_set_d(startDiam,DEFAULTDIAM,GMP_RNDD);

  mpfr_mul(rangeDiameter,rangeDiameter,startDiam,GMP_RNDD);

  mpfr_clear(startDiam);

  sollya_mpfi_interv_fr(rangeI,*(range.a),*(range.b));
  deriv = differentiate(func);

  if (getNumeratorDenominator(&numeratorDeriv,&denominatorDeriv,deriv)) {
    printMessage(1,SOLLYA_MSG_DERIVATIVE_IS_QUOTIENT,"Warning: the derivative of the function is a quotient, thus possibly not continuous in the interval.\n");
    printMessage(1,SOLLYA_MSG_CONTINUATION,"Only the zeros of the numerator will be searched and pole detection may fail.\n");
    printMessage(1,SOLLYA_MSG_CONTINUATION,"Be sure that the function is twice continuously differentiable if trusting the infnorm result.\n");

    mpfr_init2(z,prec);
    mpfr_init2(ya,prec);
    mpfr_init2(yb,prec);

    derivDenominatorDeriv = differentiate(denominatorDeriv);

    newtonWorked = newtonMPFR(z, denominatorDeriv, derivDenominatorDeriv, *(range.a), *(range.b), prec);

    if (newtonWorked && mpfr_number_p(z)) {
      evaluate(ya,numeratorDeriv,z,prec);
      evaluate(yb,denominatorDeriv,z,prec);

      mpfr_abs(ya,ya,GMP_RNDN);
      mpfr_abs(yb,yb,GMP_RNDN);

      mpfr_mul_ui(yb,yb,2,GMP_RNDN);

      if (mpfr_cmp(ya,yb) <= 0) {
	printMessage(1,SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_EXTENSIBLE_SINGULARITY,"Warning: the derivative of the function seems to have a extensible singularity in %v.\nThe infnorm result might not be trustful if the derivative cannot actually\nbe extended in this point.\n",z);
      } else {
	printMessage(1,SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_SINGULARITY,"Warning: the derivative of the function seems to have a singularity in %v.\nThe infnorm result is likely to be wrong.\n",z);
      }
    } else {
      evaluate(ya,denominatorDeriv,*(range.a),prec);
      evaluate(yb,denominatorDeriv,*(range.b),prec);

      if (mpfr_sgn(ya) != mpfr_sgn(yb)) {
	printMessage(1,SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_EXTENSIBLE_SINGULARITY,"Warning: the derivative of the function seems to have a (extensible) singularity in the considered interval.\n");
	printMessage(1,SOLLYA_MSG_CONTINUATION,"The infnorm result might be not trustful if the function is not continuously differentiable.\n");
      } else {
	printMessage(2,SOLLYA_MSG_DERIVATIVE_SEEMS_NOT_TO_HAVE_ANY_POLE,"Information: the derivative seems to have no (false) pole in the considered interval.\n");
      }
    }

    mpfr_clear(z);
    mpfr_clear(ya);
    mpfr_clear(yb);
    free_memory(derivDenominatorDeriv);
    free_memory(denominatorDeriv);
  }
  derivNumeratorDeriv = differentiate(numeratorDeriv);

  mpfr_init2(currDiameter, prec);
  mpfr_init2(stopDiameter, prec);

  mpfr_div_2ui(stopDiameter,rangeDiameter,20,GMP_RNDD);

  okay = 0;

  while (prec <= 512 * startPrec) {

    mpfr_set(currDiameter,rangeDiameter,GMP_RNDD);

    while (mpfr_cmp(currDiameter,stopDiameter) >= 0) {

      infnormI(resI,func,deriv,numeratorDeriv,derivNumeratorDeriv,rangeI,
	       prec,currDiameter,initialExcludes,NULL,NULL);

      sollya_mpfi_get_left(resultDown,resI);
      sollya_mpfi_get_right(resultUp,resI);

      if (mpfr_cmp(resultDown,resultUp) == 0) {
	okay = 1;
	break;
      }

      mpfr_nextabove(resultDown);

      if (mpfr_cmp(resultDown,resultUp) == 0) {
	okay = 1;
	break;
      }

      mpfr_div_2ui(currDiameter,currDiameter,2,GMP_RNDD);

      printMessage(4,SOLLYA_MSG_ABS_DIAM_AND_PREC_SET_TO_CERTAIN_VALUES,"Information: the absolute diameter is now %v.\nThe current intermediate precision is %d bits.\n",currDiameter,(int) prec);

    }

    if (okay) break;

    prec *= 2;

    printMessage(4,SOLLYA_MSG_INTERMEDIATE_PRECISION_HAS_BEEN_INCREASED,"Information: the intermediate precision is now %d bits.\n",(int) prec);

  }

  if (okay) mpfr_set(result,resultUp,GMP_RNDU);

  mpfr_clear(*(res.a));
  mpfr_clear(*(res.b));
  safeFree(res.a);
  safeFree(res.b);

  mpfr_clear(stopDiameter);
  mpfr_clear(currDiameter);
  free_memory(deriv);
  free_memory(numeratorDeriv);
  free_memory(derivNumeratorDeriv);
  sollya_mpfi_clear(rangeI);
  sollya_mpfi_clear(resI);
  mpfr_clear(rangeDiameter);
  mpfr_clear(resultUp);
  mpfr_clear(resultDown);

  taylorrecursions = oldtaylorrecursions;

  return okay;
}

int sollya_mpfi_have_common_real_point(sollya_mpfi_t a, sollya_mpfi_t b) {
  mp_prec_t ap, bp;
  mpfr_t al, ar, bl, br;
  int res;

  ap = sollya_mpfi_get_prec(a);
  bp = sollya_mpfi_get_prec(b);
  mpfr_init2(al, ap);
  mpfr_init2(ar, ap);
  mpfr_init2(bl, bp);
  mpfr_init2(br, bp);

  sollya_mpfi_get_left(al, a);
  sollya_mpfi_get_right(ar, a);
  sollya_mpfi_get_left(bl, b);
  sollya_mpfi_get_right(br, b);

  if (mpfr_number_p(al) &&
      mpfr_number_p(al) &&
      mpfr_number_p(al) &&
      mpfr_number_p(al)) {
    res = (mpfr_cmp(al, br) <= 0) && (mpfr_cmp(bl, ar) <= 0);
  } else {
    res = 0;
  }

  mpfr_clear(al);
  mpfr_clear(ar);
  mpfr_clear(bl);
  mpfr_clear(br);

  return res;
}

static inline mpfr_t *chooseAndInitMpfrPtr(mpfr_t *localPtr, mp_prec_t prec) {
  mpfr_t *ptr;

  ptr = getReusedGlobalMPFRVars(1, prec);
  if (ptr == NULL) {
    mpfr_init2(*localPtr, prec);
    ptr = localPtr;
  }
  return ptr;
}

static inline void clearChosenMpfrPtr(mpfr_t *ptr, mpfr_t *localPtr) {
  if (ptr == localPtr) {
    mpfr_clear(*localPtr);
    return;
  }
  returnReusedGlobalMPFRVars(1);
}

static inline point_eval_t __tryFaithEvaluationOptimizedDoIt(mpfr_t, node *, mpfr_t, mp_exp_t, mp_prec_t, mp_prec_t *, node *);

static inline void __tryFaithEvaluationOptimizedUpdateMaxPrec(mp_prec_t *maxPrec, mp_prec_t prec) {
  if (maxPrec == NULL) return;
  if (prec > *maxPrec) *maxPrec = prec;
}

static inline point_eval_t __tryFaithEvaluationOptimizedAddSubInner(int *retry, int *newPrecSet, mp_prec_t *newPrec,
								    mpfr_t y, int subtract, node *g, node *h, mpfr_t x, mp_exp_t cutoff,
								    mp_prec_t prec, mp_prec_t minPrec,
								    mp_prec_t *maxPrecUsed) {
  mpfr_t v_gy, v_hy, v_t;
  mpfr_t *gy, *hy, *t;
  point_eval_t resG, resH, res;
  int ternary, tern1, tern2;
  mp_exp_t recCutoff, recCutoffG, recCutoffH, recCutoffGP, recCutoffHP;
  mp_prec_t recMaxPrecUsedG, recMaxPrecUsedH;
  sollya_mpfi_t v_X, v_Y, v_Z;
  sollya_mpfi_t *X, *Y, *Z;
  int zeroG, zeroH;
  mp_exp_t expBeforeCancel, expAfterCancel;
  mp_prec_t lostPrec, newPrecCutoff, precG, precH;

  /* Make compiler happy */
  X = NULL;
  Y = NULL;
  Z = NULL;
  /* End of compiler happiness */


  *retry = 0;
  *newPrecSet = 0;

  if (cutoff <= mpfr_get_emin_min()) {
    recCutoff = mpfr_get_emin_min();
  } else {
    recCutoff = cutoff - 3;
    if ((recCutoff >= 0) || (recCutoff < mpfr_get_emin_min())) {
      recCutoff = mpfr_get_emin_min();
    }
  }
  recCutoffG = recCutoff;
  recCutoffH = recCutoff;
  precG = prec;
  precH = prec;

  gy = chooseAndInitMpfrPtr(&v_gy, precG);
  recMaxPrecUsedG = 0;
  resG = __tryFaithEvaluationOptimizedDoIt(*gy, g, x, recCutoffG, minPrec, &recMaxPrecUsedG, g);
  switch (resG) {
  case POINT_EVAL_EXACT:
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    if (mpfr_number_p(*gy) && (!mpfr_zero_p(*gy))) {
      recCutoffHP = mpfr_get_exp(*gy) - mpfr_get_prec(y) - 5;
      if (recCutoffHP >= 0) recCutoffHP = -1;
      if (recCutoffHP > recCutoffH) recCutoffH = recCutoffHP;
    }
    break;
  default:
    break;
  }
  hy = chooseAndInitMpfrPtr(&v_hy, precH);
  recMaxPrecUsedH = 0;
  resH = __tryFaithEvaluationOptimizedDoIt(*hy, h, x, recCutoffH, minPrec, &recMaxPrecUsedH, h);
  if (resG == POINT_EVAL_FAILURE) {
    switch (resH) {
    case POINT_EVAL_EXACT:
    case POINT_EVAL_CORRECTLY_ROUNDED:
    case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    case POINT_EVAL_FAITHFULLY_ROUNDED:
    case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
      if (mpfr_number_p(*hy) && (!mpfr_zero_p(*hy))) {
	recCutoffGP = mpfr_get_exp(*hy) - mpfr_get_prec(y) - 5;
	if (recCutoffGP >= 0) recCutoffGP = -1;
	if (recCutoffGP > recCutoffG) {
	  recCutoffG = recCutoffGP;
	  recMaxPrecUsedH = 0;
	  resG = __tryFaithEvaluationOptimizedDoIt(*gy, g, x, recCutoffG, minPrec, &recMaxPrecUsedG, g);
	}
      }
      break;
    default:
      break;
    }
  }

  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsedG);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsedH);

  if ((resG == POINT_EVAL_FAILURE) ||
      (resH == POINT_EVAL_FAILURE)) {
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    *retry = 0;
    return POINT_EVAL_FAILURE;
  }

  if ((resG == POINT_EVAL_EXACT) &&
      (resH == POINT_EVAL_EXACT)) {
    if (subtract) {
      ternary = mpfr_sub(y, *gy, *hy, GMP_RNDN);
    } else {
      ternary = mpfr_add(y, *gy, *hy, GMP_RNDN);
    }
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }

  if ((!mpfr_number_p(*gy)) ||
      (!mpfr_number_p(*hy))) {
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    *retry = 0;
    return POINT_EVAL_FAILURE;
  }

  switch (resG) {
  case POINT_EVAL_FAILURE:
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(*gy));
    sollya_mpfi_set_fr(*X, *gy);
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(*gy) + 1);
    sollya_mpfi_set_fr(*X, *gy);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*gy));
    sollya_mpfi_set_fr(*X, *gy);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_BELOW_CUTOFF:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(*gy));
    sollya_mpfi_interv_si_2exp(*X, -1, recCutoffG, 1, recCutoffG);
    break;
  }

  switch (resH) {
  case POINT_EVAL_FAILURE:
    clearChosenMpfiPtr(X, &v_X);
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy));
    sollya_mpfi_set_fr(*Y, *hy);
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy) + 1);
    sollya_mpfi_set_fr(*Y, *hy);
    sollya_mpfi_blow_1ulp(*Y);
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy));
    sollya_mpfi_set_fr(*Y, *hy);
    sollya_mpfi_blow_1ulp(*Y);
    break;
  case POINT_EVAL_BELOW_CUTOFF:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy));
    sollya_mpfi_interv_si_2exp(*Y, -1, recCutoffH, 1, recCutoffH);
    break;
  }
  Z = chooseAndInitMpfiPtr(&v_Z, (mpfr_get_prec(*gy) > mpfr_get_prec(*hy)? mpfr_get_prec(*gy) : mpfr_get_prec(*hy)) + 10);

  if (subtract) {
    sollya_mpfi_sub(*Z, *X, *Y);
  } else {
    sollya_mpfi_add(*Z, *X, *Y);
  }

  t = chooseAndInitMpfrPtr(&v_t, mpfr_get_prec(y));

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tern1 = mpfr_set(y, &((*Z)->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(*t, &((*Z)->right), GMP_RNDN); /* rounds to final precision */

  if (mpfr_number_p(*t) && mpfr_number_p(y)) {
    if (mpfr_equal_p(*t, y)) {
      if ((tern1 == 0) && (tern2 == 0)) {
	res = POINT_EVAL_EXACT;
      } else {
	if ((tern1 != 0) && (tern1 == tern2)) {
	  res = POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
	} else {
	  res = POINT_EVAL_CORRECTLY_ROUNDED;
	}
      }
    } else {
      if (mpfr_cmp(y, *t) < 0) {
	mpfr_nextbelow(*t);
	if (mpfr_equal_p(*t, y)) {
	  if (!sollya_mpfi_fr_in_interval(y, *Z)) {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT;
	  } else {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED;
	  }
	} else {
	  if ((!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z))) && (sollya_mpfi_max_exp(*Z) < cutoff)) {
	    res = POINT_EVAL_BELOW_CUTOFF;
	  } else {
	    res = POINT_EVAL_FAILURE;
	  }
	}
      } else {
	if ((!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z))) && (sollya_mpfi_max_exp(*Z) < cutoff)) {
	  res = POINT_EVAL_BELOW_CUTOFF;
	} else {
	  res = POINT_EVAL_FAILURE;
	}
      }
    }
  } else {
    if ((!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z))) && (sollya_mpfi_max_exp(*Z) < cutoff)) {
      res = POINT_EVAL_BELOW_CUTOFF;
    } else {
      res = POINT_EVAL_FAILURE;
    }
  }

  /* If we couldn't compute an answer, try find a precision at which
     we might be able to compute an answer
  */
  if ((res == POINT_EVAL_FAILURE) &&
      (resG != POINT_EVAL_FAILURE) &&
      (resH != POINT_EVAL_FAILURE) &&
      mpfr_number_p(*gy) &&
      mpfr_number_p(*hy) &&
      (!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z)))) {
    zeroG = mpfr_zero_p(*gy);
    zeroH = mpfr_zero_p(*hy);
    if (zeroG && zeroH) {
      /* f(x) = g(x) + h(x), g(x) and h(x) both cancel, cutoff not yet
	 met.
      */
      *retry = 1;
      *newPrecSet = 0;
    } else {
      if (!(zeroG || zeroH)) {
	/* Got exponent of both terms g(x) and h(x) */
	if (sollya_mpfi_is_zero(*Z)) {
	  /* Cancels out completely, need at least all cancelled bits
	     plus the bits of y to get the bits of y
	  */
	  *newPrec = mpfr_get_prec(y) + sollya_mpfi_get_prec(*Z) + 5;
	  *newPrecSet = 1;
	  *retry = 1;
	} else {
	  /* Does not cancel out completely; estimate the number of
	     bits that cancelled.
	  */
	  if (zeroG) {
	    expBeforeCancel = mpfr_get_exp(*hy);
	  } else {
	    if (zeroH) {
	      expBeforeCancel = mpfr_get_exp(*gy);
	    } else {
	      expBeforeCancel = mpfr_get_exp(*gy);
	      if (mpfr_get_exp(*hy) > expBeforeCancel) expBeforeCancel = mpfr_get_exp(*hy);
	    }
	  }
	  expAfterCancel = sollya_mpfi_max_exp(*Z);
	  lostPrec = expBeforeCancel - expAfterCancel + 4;
	  *newPrec = prec + lostPrec + 10;
	  if (*newPrec < 12) *newPrec = 12;
	  if (*newPrec < minPrec) *newPrec = minPrec;
	  if (*newPrec < (prec + 10)) *newPrec = prec + 10;
	  if ((lostPrec < 1) || (*newPrec > 2 * prec + 10)) *newPrec = 2 * prec + 10;
	  if ((expBeforeCancel > mpfr_get_emin_min() + 10 + 2 * prec) &&
	      (expBeforeCancel - (2 * prec + 10) <= cutoff) &&
	      (expBeforeCancel > cutoff)) {
	    newPrecCutoff = expBeforeCancel - cutoff;
	    if (newPrecCutoff > *newPrec) {
	      *newPrec = newPrecCutoff;
	    }
	  }
	  *newPrecSet = 1;
	  *retry = 1;
	}
      }
    }
  }


  clearChosenMpfrPtr(t, &v_t);
  clearChosenMpfiPtr(Z, &v_Z);
  clearChosenMpfiPtr(Y, &v_Y);
  clearChosenMpfiPtr(X, &v_X);
  clearChosenMpfrPtr(hy, &v_hy);
  clearChosenMpfrPtr(gy, &v_gy);
  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedAddSubMain(mpfr_t y, int subtract, node *g, node *h, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  mp_prec_t prec, newPrec;
  point_eval_t res;
  int retry, newPrecSet;
  mp_prec_t recMaxPrecUsed;

  prec = mpfr_get_prec(y) + 20;
  if (prec < minPrec) prec = minPrec;

  /* 1st try */
  retry = 0;
  newPrecSet = 0;
  newPrec = prec;
  recMaxPrecUsed = 0;
  res = __tryFaithEvaluationOptimizedAddSubInner(&retry, &newPrecSet, &newPrec, y, subtract, g, h, x, cutoff, prec, minPrec, &recMaxPrecUsed);
  if (res != POINT_EVAL_FAILURE) {
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return res;
  }
  if (!retry) {
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return POINT_EVAL_FAILURE;
  }
  if (!newPrecSet) {
    newPrec = prec * 2;
  }
  prec = newPrec;
  if (prec < minPrec) prec = minPrec;

  /* 2nd try */
  retry = 0;
  newPrecSet = 0;
  newPrec = prec;
  recMaxPrecUsed = 0;
  res = __tryFaithEvaluationOptimizedAddSubInner(&retry, &newPrecSet, &newPrec, y, subtract, g, h, x, cutoff, prec, minPrec, &recMaxPrecUsed);
  if (res != POINT_EVAL_FAILURE) {
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return res;
  }
  if (!retry) {
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return POINT_EVAL_FAILURE;
  }
  if (!newPrecSet) {
    newPrec = prec * 2;
  }
  prec = newPrec;
  if (prec < minPrec) prec = minPrec;

  /* 3rd and last try */
  retry = 0;
  newPrecSet = 0;
  newPrec = prec;
  recMaxPrecUsed = 0;
  res = __tryFaithEvaluationOptimizedAddSubInner(&retry, &newPrecSet, &newPrec, y, subtract, g, h, x, cutoff, prec, minPrec, &recMaxPrecUsed);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);

  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedAddSub(mpfr_t y, int subtract, node *g, node *h, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  node *myG, *myH;
  int swapped;
  point_eval_t res;
  int ternary;

  if ((accessThruMemRef(g)->nodeType == CONSTANT) &&
      (accessThruMemRef(h)->nodeType == CONSTANT)) {
    if (subtract) {
      ternary = mpfr_sub(y, *(accessThruMemRef(g)->value), *(accessThruMemRef(h)->value), GMP_RNDN);
    } else {
      ternary = mpfr_add(y, *(accessThruMemRef(g)->value), *(accessThruMemRef(h)->value), GMP_RNDN);
    }
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((accessThruMemRef(g)->nodeType == CONSTANT) &&
      (accessThruMemRef(h)->nodeType == VARIABLE)) {
    if (subtract) {
      ternary = mpfr_sub(y, *(accessThruMemRef(g)->value), x, GMP_RNDN);
    } else {
      ternary = mpfr_add(y, *(accessThruMemRef(g)->value), x, GMP_RNDN);
    }
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((accessThruMemRef(g)->nodeType == VARIABLE) &&
      (accessThruMemRef(h)->nodeType == CONSTANT)) {
    if (subtract) {
      ternary = mpfr_sub(y, x, *(accessThruMemRef(h)->value), GMP_RNDN);
    } else {
      ternary = mpfr_add(y, x, *(accessThruMemRef(h)->value), GMP_RNDN);
    }
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((accessThruMemRef(g)->nodeType == VARIABLE) &&
      (accessThruMemRef(h)->nodeType == VARIABLE)) {
    if (subtract) {
      ternary = mpfr_sub(y, x, x, GMP_RNDN);
    } else {
      ternary = mpfr_add(y, x, x, GMP_RNDN);
    }
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }

  myG = g; myH = h;
  swapped = 0;
  switch (accessThruMemRef(h)->nodeType) {
  case CONSTANT:
  case VARIABLE:
    myG = h;
    myH = g;
    swapped = 1;
  default:
    break;
  }

  res = __tryFaithEvaluationOptimizedAddSubMain(y, subtract, myG, myH, x, cutoff, minPrec, maxPrecUsed);

  if (swapped && subtract && (res != POINT_EVAL_FAILURE)) {
    mpfr_neg(y, y, GMP_RNDN); /* exact */
  }

  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedMulDivInner(int *retry, mpfr_t y, int divide, node *g, node *h, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  mp_prec_t precG, precH;
  mpfr_t v_gy, v_hy;
  mpfr_t *gy, *hy;
  point_eval_t resG, resH;
  int ternary;
  mp_exp_t recCutoff;
  mp_prec_t recMaxPrecUsed;

  if (cutoff >= mpfr_get_emin_min() + 64) {
    recCutoff = cutoff - 64;
    if ((recCutoff >= 0) || (recCutoff < mpfr_get_emin_min())) recCutoff = mpfr_get_emin_min();
  } else {
    recCutoff = mpfr_get_emin_min();
  }
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    precG = mpfr_get_prec(*(accessThruMemRef(g)->value));
  } else {
    if (accessThruMemRef(g)->nodeType == VARIABLE) {
      precG = mpfr_get_prec(x);
    } else {
      precG = mpfr_get_prec(y) + 4;
    }
  }
  if (accessThruMemRef(h)->nodeType == CONSTANT) {
    precH = mpfr_get_prec(*(accessThruMemRef(h)->value));
  } else {
    if (accessThruMemRef(h)->nodeType == VARIABLE) {
      precH = mpfr_get_prec(x);
    } else {
      precH = mpfr_get_prec(y) + 4;
    }
  }
  if (precG < minPrec) precG = minPrec;
  if (precH < minPrec) precH = minPrec;
  gy = chooseAndInitMpfrPtr(&v_gy, precG);
  hy = chooseAndInitMpfrPtr(&v_hy, precH);
  recMaxPrecUsed = 0;
  resG = __tryFaithEvaluationOptimizedDoIt(*gy, g, x, recCutoff, minPrec, &recMaxPrecUsed, g);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
  if (resG == POINT_EVAL_FAILURE) {
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    *retry = 0;
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return POINT_EVAL_FAILURE;
  }
  recMaxPrecUsed = 0;
  resH = __tryFaithEvaluationOptimizedDoIt(*hy, h, x, (divide ? mpfr_get_emin_min() : recCutoff), minPrec, &recMaxPrecUsed, h);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
  if (resH == POINT_EVAL_FAILURE) {
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    *retry = 0;
    return POINT_EVAL_FAILURE;
  }
  if ((resG == POINT_EVAL_EXACT) &&
      (resH == POINT_EVAL_EXACT)) {
    if (divide) {
      if (mpfr_zero_p(*hy)) {
	clearChosenMpfrPtr(hy, &v_hy);
	clearChosenMpfrPtr(gy, &v_gy);
	*retry = 0;
	return POINT_EVAL_FAILURE;
      }
      ternary = mpfr_div(y, *gy, *hy, GMP_RNDN);
    } else {
      ternary = mpfr_mul(y, *gy, *hy, GMP_RNDN);
    }
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((resG == POINT_EVAL_EXACT) &&
      mpfr_zero_p(*gy)) {
    mpfr_set_si(y, 0, GMP_RNDN);
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_EXACT;
  }
  if ((resH == POINT_EVAL_EXACT) &&
      mpfr_zero_p(*hy) &&
      (!divide)) {
    mpfr_set_si(y, 0, GMP_RNDN);
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_EXACT;
  }
  if ((resG == POINT_EVAL_BELOW_CUTOFF) ||
      (resH == POINT_EVAL_BELOW_CUTOFF)) {
    if ((resH == POINT_EVAL_BELOW_CUTOFF) && divide) {
      /* This case should never happen */
      clearChosenMpfrPtr(hy, &v_hy);
      clearChosenMpfrPtr(gy, &v_gy);
      *retry = 0;
      return POINT_EVAL_FAILURE;
    }
    if ((resG == POINT_EVAL_BELOW_CUTOFF) &&
	(resH == POINT_EVAL_BELOW_CUTOFF)) {
      clearChosenMpfrPtr(hy, &v_hy);
      clearChosenMpfrPtr(gy, &v_gy);
      mpfr_set_si(y, 0, GMP_RNDN); /* exact */
      return POINT_EVAL_BELOW_CUTOFF;
    }
    if (cutoff <= mpfr_get_emin_min() + 32) {
      clearChosenMpfrPtr(hy, &v_hy);
      clearChosenMpfrPtr(gy, &v_gy);
      *retry = 1;
      return POINT_EVAL_FAILURE;
    }
    if (resG == POINT_EVAL_BELOW_CUTOFF) {
      if (divide) {
	if (mpfr_number_p(*hy) &&
	    (!mpfr_zero_p(*hy)) &&
	    (((mp_exp_t) 2) - mpfr_get_exp(*hy) < cutoff - recCutoff)) {
	  clearChosenMpfrPtr(hy, &v_hy);
	  clearChosenMpfrPtr(gy, &v_gy);
	  mpfr_set_si(y, 0, GMP_RNDN); /* exact */
	  return POINT_EVAL_BELOW_CUTOFF;
	}
      } else {
	if (mpfr_number_p(*hy) &&
	    (mpfr_zero_p(*hy) ||
	     (mpfr_get_exp(*hy) + ((mp_exp_t) 2) < cutoff - recCutoff))) {
	  clearChosenMpfrPtr(hy, &v_hy);
	  clearChosenMpfrPtr(gy, &v_gy);
	  mpfr_set_si(y, 0, GMP_RNDN); /* exact */
	  return POINT_EVAL_BELOW_CUTOFF;
	}
      }
    }
    if (resH == POINT_EVAL_BELOW_CUTOFF) {
      if (!divide) {
	if (mpfr_number_p(*gy) &&
	    (mpfr_zero_p(*gy) ||
	     (mpfr_get_exp(*gy) + ((mp_exp_t) 2) < cutoff - recCutoff))) {
	  clearChosenMpfrPtr(hy, &v_hy);
	  clearChosenMpfrPtr(gy, &v_gy);
	  mpfr_set_si(y, 0, GMP_RNDN); /* exact */
	  return POINT_EVAL_BELOW_CUTOFF;
	}
      }
    }
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    *retry = 1;
    return POINT_EVAL_FAILURE;
  }
  /* TODO: Proof required */
  if (divide) {
    ternary = mpfr_div(y, *gy, *hy, GMP_RNDN);
  } else {
    ternary = mpfr_mul(y, *gy, *hy, GMP_RNDN);
  }
  clearChosenMpfrPtr(hy, &v_hy);
  clearChosenMpfrPtr(gy, &v_gy);
  /* TODO: work for exactness predicates */
  if (ternary == 0) return POINT_EVAL_CORRECTLY_ROUNDED;
  return POINT_EVAL_FAITHFULLY_ROUNDED;
}

static inline point_eval_t __tryFaithEvaluationOptimizedMulDiv(mpfr_t y, int divide, node *g, node *h, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  int retry;
  point_eval_t res;
  mp_prec_t recMaxPrecUsed;
  int ternary;

  if ((accessThruMemRef(g)->nodeType == CONSTANT) &&
      (accessThruMemRef(h)->nodeType == CONSTANT)) {
    if (divide) {
      ternary = mpfr_div(y, *(accessThruMemRef(g)->value), *(accessThruMemRef(h)->value), GMP_RNDN);
    } else {
      ternary = mpfr_mul(y, *(accessThruMemRef(g)->value), *(accessThruMemRef(h)->value), GMP_RNDN);
    }
    if (mpfr_number_p(y)) {
      if (ternary == 0) return POINT_EVAL_EXACT;
      return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    }
  }
  if ((accessThruMemRef(g)->nodeType == CONSTANT) &&
      (accessThruMemRef(h)->nodeType == VARIABLE)) {
    if (divide) {
      ternary = mpfr_div(y, *(accessThruMemRef(g)->value), x, GMP_RNDN);
    } else {
      ternary = mpfr_mul(y, *(accessThruMemRef(g)->value), x, GMP_RNDN);
    }
    if (mpfr_number_p(y)) {
      if (ternary == 0) return POINT_EVAL_EXACT;
      return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    }
  }
  if ((accessThruMemRef(g)->nodeType == VARIABLE) &&
      (accessThruMemRef(h)->nodeType == CONSTANT)) {
    if (divide) {
      ternary = mpfr_div(y, x, *(accessThruMemRef(h)->value), GMP_RNDN);
    } else {
      ternary = mpfr_mul(y, x, *(accessThruMemRef(h)->value), GMP_RNDN);
    }
    if (mpfr_number_p(y)) {
      if (ternary == 0) return POINT_EVAL_EXACT;
      return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    }
  }
  if ((accessThruMemRef(g)->nodeType == VARIABLE) &&
      (accessThruMemRef(h)->nodeType == VARIABLE)) {
    if (divide) {
      ternary = mpfr_div(y, x, x, GMP_RNDN);
    } else {
      ternary = mpfr_mul(y, x, x, GMP_RNDN);
    }
    if (mpfr_number_p(y)) {
      if (ternary == 0) return POINT_EVAL_EXACT;
      return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    }
  }

  retry = 0;
  recMaxPrecUsed = 0;
  res = __tryFaithEvaluationOptimizedMulDivInner(&retry, y, divide, g, h, x, cutoff, minPrec, &recMaxPrecUsed);
  if (res != POINT_EVAL_FAILURE) {
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return res;
  }
  if (!retry) {
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
    return POINT_EVAL_FAILURE;
  }

  retry = 0;
  recMaxPrecUsed = 0;
  res = __tryFaithEvaluationOptimizedMulDivInner(&retry, y, divide, g, h, x, mpfr_get_emin_min(), minPrec, &recMaxPrecUsed);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedPow(mpfr_t y, node *g, node *h, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  mp_prec_t precG, precH;
  mpfr_t v_gy, v_hy;
  mpfr_t *gy, *hy;
  point_eval_t resG, resH;
  int ternary, tern1, tern2;
  sollya_mpfi_t v_X, v_Y, v_Z;
  sollya_mpfi_t *X, *Y, *Z;
  point_eval_t res;
  mp_exp_t cutoffH;
  mp_prec_t recMaxPrecUsed;

  /* Make compiler happy */
  X = NULL;
  Y = NULL;
  /* End of compiler happiness */

  if ((accessThruMemRef(g)->nodeType == CONSTANT) &&
      (accessThruMemRef(h)->nodeType == CONSTANT)) {
    ternary = mpfr_pow(y, *(accessThruMemRef(g)->value), *(accessThruMemRef(h)->value), GMP_RNDN);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((accessThruMemRef(g)->nodeType == CONSTANT) &&
      (accessThruMemRef(h)->nodeType == VARIABLE)) {
    ternary = mpfr_pow(y, *(accessThruMemRef(g)->value), x, GMP_RNDN);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((accessThruMemRef(g)->nodeType == VARIABLE) &&
      (accessThruMemRef(h)->nodeType == CONSTANT)) {
    ternary = mpfr_pow(y, x, *(accessThruMemRef(h)->value), GMP_RNDN);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  if ((accessThruMemRef(g)->nodeType == VARIABLE) &&
      (accessThruMemRef(h)->nodeType == VARIABLE)) {
    ternary = mpfr_pow(y, x, x, GMP_RNDN);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    precG = mpfr_get_prec(*(accessThruMemRef(g)->value));
  } else {
    if (accessThruMemRef(g)->nodeType == VARIABLE) {
      precG = mpfr_get_prec(x);
    } else {
      precG = mpfr_get_prec(y) + 25;
    }
  }
  if (accessThruMemRef(h)->nodeType == CONSTANT) {
    precH = mpfr_get_prec(*(accessThruMemRef(h)->value));
  } else {
    if (accessThruMemRef(h)->nodeType == VARIABLE) {
      precH = mpfr_get_prec(x);
    } else {
      precH = mpfr_get_prec(y) + 25;
    }
  }
  if (precG < minPrec) precG = minPrec;
  if (precH < minPrec) precH = minPrec;
  gy = chooseAndInitMpfrPtr(&v_gy, precG);
  hy = chooseAndInitMpfrPtr(&v_hy, precH);
  recMaxPrecUsed = 0;
  resG = __tryFaithEvaluationOptimizedDoIt(*gy, g, x, mpfr_get_emin_min(), minPrec, &recMaxPrecUsed, g);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
  if ((resG == POINT_EVAL_FAILURE) ||
      (resG == POINT_EVAL_BELOW_CUTOFF)) {
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_FAILURE;
  }
  cutoffH = -(2 * mpfr_get_prec(y) + 25);
  if ((cutoffH >= 0) || (cutoffH < mpfr_get_emin_min())) cutoffH = mpfr_get_emin_min();
  recMaxPrecUsed = 0;
  resH = __tryFaithEvaluationOptimizedDoIt(*hy, h, x, cutoffH, minPrec, &recMaxPrecUsed, h);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
  if (resH == POINT_EVAL_FAILURE) {
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_FAILURE;
  }
  if ((resG == POINT_EVAL_EXACT) &&
      (resH == POINT_EVAL_EXACT)) {
    ternary = mpfr_pow(y, *gy, *hy, GMP_RNDN);
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }
  switch (resG) {
  case POINT_EVAL_FAILURE:
  case POINT_EVAL_BELOW_CUTOFF:
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(*gy));
    sollya_mpfi_set_fr(*X, *gy);
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(*gy) + 1);
    sollya_mpfi_set_fr(*X, *gy);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(*gy));
    sollya_mpfi_set_fr(*X, *gy);
    sollya_mpfi_blow_1ulp(*X);
    break;
  }
  switch (resH) {
  case POINT_EVAL_FAILURE:
    clearChosenMpfiPtr(X, &v_X);
    clearChosenMpfrPtr(hy, &v_hy);
    clearChosenMpfrPtr(gy, &v_gy);
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy));
    sollya_mpfi_set_fr(*Y, *hy);
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy) + 1);
    sollya_mpfi_set_fr(*Y, *hy);
    sollya_mpfi_blow_1ulp(*Y);
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy));
    sollya_mpfi_set_fr(*Y, *hy);
    sollya_mpfi_blow_1ulp(*Y);
    break;
  case POINT_EVAL_BELOW_CUTOFF:
    Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(*hy));
    sollya_mpfi_interv_si_2exp(*Y, -1, cutoffH, 1, cutoffH);
    break;
  }
  Z = chooseAndInitMpfiPtr(&v_Z, mpfr_get_prec(y) + 15);
  sollya_mpfi_pow(*Z, *X, *Y);

  /* Reuse hy */
  mpfr_set_prec(*hy, mpfr_get_prec(y));

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tern1 = mpfr_set(y, &((*Z)->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(*hy, &((*Z)->right), GMP_RNDN); /* rounds to final precision */

  if (mpfr_number_p(*hy) && mpfr_number_p(y)) {
    if (mpfr_equal_p(*hy, y)) {
      if ((tern1 == 0) && (tern2 == 0)) {
	res = POINT_EVAL_EXACT;
      } else {
	if ((tern1 != 0) && (tern1 == tern2)) {
	  res = POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
	} else {
	  res = POINT_EVAL_CORRECTLY_ROUNDED;
	}
      }
    } else {
      if (mpfr_cmp(y, *hy) < 0) {
	mpfr_nextbelow(*hy);
	if (mpfr_equal_p(*hy, y)) {
	  if (!sollya_mpfi_fr_in_interval(y, *Z)) {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT;
	  } else {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED;
	  }
	} else {
	  if ((!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z))) && (sollya_mpfi_max_exp(*Z) < cutoff)) {
	    res = POINT_EVAL_BELOW_CUTOFF;
	  } else {
	    res = POINT_EVAL_FAILURE;
	  }
	}
      } else {
	if ((!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z))) && (sollya_mpfi_max_exp(*Z) < cutoff)) {
	  res = POINT_EVAL_BELOW_CUTOFF;
	} else {
	  res = POINT_EVAL_FAILURE;
	}
      }
    }
  } else {
    if ((!(sollya_mpfi_has_nan(*Z) || sollya_mpfi_has_infinity(*Z))) && (sollya_mpfi_max_exp(*Z) < cutoff)) {
      res = POINT_EVAL_BELOW_CUTOFF;
    } else {
      res = POINT_EVAL_FAILURE;
    }
  }

  clearChosenMpfiPtr(Z, &v_Z);
  clearChosenMpfiPtr(Y, &v_Y);
  clearChosenMpfiPtr(X, &v_X);
  clearChosenMpfrPtr(hy, &v_hy);
  clearChosenMpfrPtr(gy, &v_gy);

  return res;
}

/* Try to determine a cutoff for the recursive evaluation of a
   univariate function.

   For an evaluation of a function f(g(t)), x := g(t), we
   have to find (and return) a cutoff k such that

   |g(t)| = |x| <= 2^k

   implies either that

   *   |f(x)| <= 2^cutoff       where cutoff is the resp. argument

   or that for X = g(t) + delta = x + delta satisfying also |X| <= 2^k,

   *   f(x + delta) =  f(x) + delta * f'(xi)
   =  f(x) * (1 + delta * f'(xi)/f(x))   with xi between x and x + delta
   =: f(x) * (1 + eps)

   we got eps = delta * f'(xi)/f(x) smaller than 2^-prec, i.e.

   |delta| * |f'(xi)/f(x)| <= 2^-prec    where prec is the resp. argument.

   In this second case, typically we have |delta| <= 2^(k + 1).

   This (C) function is based on heuristics that do not endanger
   safeness, though. This means that if we return a (dummy) value that
   does not satisfy the above conditions, no mathematically wrong
   result will be computed, only evaluation performance will be
   affected.

*/
static inline mp_exp_t __tryFaithEvaluationOptimizedUnivariateGetRecurseCutoff(node *f, mp_exp_t cutoff, mp_prec_t prec) {
  mp_exp_t res;

  /* Check that incoming cutoff and precision are within resonable
     range
  */
  if ((cutoff > -1) ||
      (prec < 1)) {
    return mpfr_get_emin_min();
  }
  switch (f->nodeType) {
  case UNARY_BASE_FUNC:
    return f->baseFun->getRecurseCutoff(cutoff, prec);
    break;
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    /* Actually, we have no idea how the function behaves */
    res = cutoff - 1;
    if ((res >= 0) || (res < mpfr_get_emin_min())) res = mpfr_get_emin_min();
    if (-(prec + 5) < res) res = -(prec + 5);
    if ((res >= 0) || (res < mpfr_get_emin_min())) res = mpfr_get_emin_min();
    return res;
    break;
  default:
    break;
  }
  return mpfr_get_emin_min();
}


/* Try to determine a precision for the recursive evaluation of a
   univariate function.

   For an evaluation of a function f(g(t)), where g gets
   evaluated with a precision such that

   X = g(t) * (1 + eps) = x * (1 + eps)

   we have

   f(x * (1 + eps)) = f(x + x * eps)
   = f(x) + x * eps * f'(xi)   with xi between x and x + x * eps
   = f(x) * (1 + eps * (x * f'(xi))/f(x))

   We have to find and return a precision k such that for

   |eps| <= 2^-k

   we got

   |eps * (x * f'(xi))/f(x)| = |eps| * 2^(log2(|x * f'(xi)/f(x)|)) <= 2^-prec,

   where prec is an argument to the (C) function.

   We are going to solve the problem approximatively, which does not endanger
   safeness, as a check on the appropriate accuracy will be performed anyway.

   This means that we are going to consider the term above at first
   order only and take xi = x.

   Additionnally, for certain functions

   |x * f'(x)/f(x)|

   does not stay bounded for some x (as x is not bounded whatsoever,
   or just be the possible domain of MPFR numbers), we only consider
   cases when x is bounded in a reasonable, heuristically determined
   domain. An example of such a function is exponential, for
   which we have |x * f'(x)/f(x)| = |x| and where we bound x with
   |x| <= 2^12.

   In cases when f has a zero at x, |x * f'(x)/f(x)| may not stay
   bounded. In this case, this (C) function may return, in lieu
   of a precision k satisfying the above defined conditions, a precision
   k such that for |eps| <= 2^-k

   f(x * (1 + eps)) = f(x) + x * eps * f'(xi)

   is smaller than 2^cutoff, i.e.

   |f(x * (1 + eps))| = |f(x) + x * eps * f'(xi)|
   = |x * eps * f'(xi)|             (as f(x) = 0 in this case)
   = |eps| * 2^(log2(|x * f'(xi)|)) <= 2^cutoff.

*/
static inline mp_prec_t __tryFaithEvaluationOptimizedUnivariateGetRecursePrec(node *f, mp_prec_t prec, mp_exp_t cutoff) {
  int considerCutoff;
  mp_prec_t cutoffPrec;

  /* Make compiler happy */
  cutoffPrec = 12;
  /* End of compiler happiness */

  /* Consider log2(abs(x * diff(f)/f)) over some reasonable range

     Consider cutoff if it implies precision in some reasonable range.

  */
  considerCutoff = 0;
  if (cutoff >= -(2 * prec + 10)) {
    cutoffPrec = -cutoff + 10;
    if (cutoffPrec < 12) cutoffPrec = 12;
    if (cutoffPrec > 2 * prec + 10) cutoffPrec = 2 * prec + 10;
    considerCutoff = 1;
  }
  switch (f->nodeType) {
  case UNARY_BASE_FUNC:
    return f->baseFun->getRecursePrec(cutoffPrec, prec, considerCutoff);
    break;
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    /* Actually, we have no idea how the function behaves */
    return prec + 15;
    break;
  default:
    break;
  }
  /* Some default answer */
  return prec + 10;
}

static inline point_eval_t __tryFaithEvaluationOptimizedUnivariateImpreciseArg(mpfr_t y, node *f, mpfr_t x, point_eval_t err, mp_exp_t cutoffY, mp_exp_t cutoffX, node *mr) {
  sollya_mpfi_t v_X, v_Y;
  sollya_mpfi_t *X, *Y;
  point_eval_t res;
  int tern1, tern2;
  mpfi_t temp; /* The type is not a mistake, we use a library function */
  X = NULL;

  if (err == POINT_EVAL_FAILURE) return POINT_EVAL_FAILURE;
  if (!mpfr_number_p(x)) return POINT_EVAL_FAILURE;

  switch (err) {
  case POINT_EVAL_FAILURE:
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));
    sollya_mpfi_set_fr(*X, x);
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x) + 1);
    sollya_mpfi_set_fr(*X, x);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));
    sollya_mpfi_set_fr(*X, x);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_BELOW_CUTOFF:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));
    sollya_mpfi_interv_si_2exp(*X, -1, cutoffX, 1, cutoffX);
    break;
  }
  Y = chooseAndInitMpfiPtr(&v_Y, mpfr_get_prec(y) + 8);
  switch (f->nodeType) {
  case UNARY_BASE_FUNC:
    if (f->baseFun->faithEvaluationOptimizedSupported) f->baseFun->interval_eval(*Y, *X);
    break;
  case LIBRARYFUNCTION:
    mpfi_init2(temp, sollya_mpfi_get_prec(*Y));
    if (f->libFun->hasData) {
      enterExternalCode();
      ((int (*)(mpfi_t, mpfi_t, int, void *)) (f->libFun->code))(temp, *X, f->libFunDeriv, f->libFun->data);
      leaveExternalCode();
    } else {
      enterExternalCode();
      ((int (*)(mpfi_t, mpfi_t, int)) (f->libFun->code))(temp, *X, f->libFunDeriv);
      leaveExternalCode();
    }
    mpfi_to_sollya_mpfi(*Y, temp);
    mpfi_clear(temp);
    break;
  case PROCEDUREFUNCTION:
    computeFunctionWithProcedure(*Y, f->child2, *X, (unsigned int) f->libFunDeriv, mr);
    break;
  default:
    clearChosenMpfiPtr(Y, &v_Y);
    clearChosenMpfiPtr(X, &v_X);
    return POINT_EVAL_FAILURE;
    break;
  }

  /* This destroys x */
  mpfr_set_prec(x, mpfr_get_prec(y));

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tern1 = mpfr_set(y, &((*Y)->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(x, &((*Y)->right), GMP_RNDN); /* rounds to final precision */

  if (mpfr_number_p(x) && mpfr_number_p(y)) {
    if (mpfr_equal_p(x, y)) {
      if ((tern1 == 0) && (tern2 == 0)) {
	res = POINT_EVAL_EXACT;
      } else {
	if ((tern1 != 0) && (tern1 == tern2)) {
	  res = POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
	} else {
	  res = POINT_EVAL_CORRECTLY_ROUNDED;
	}
      }
    } else {
      if (mpfr_cmp(y, x) < 0) {
	mpfr_nextbelow(x);
	if (mpfr_equal_p(x, y)) {
	  if (!sollya_mpfi_fr_in_interval(y, *Y)) {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT;
	  } else {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED;
	  }
	} else {
	  if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoffY)) {
	    res = POINT_EVAL_BELOW_CUTOFF;
	  } else {
	    res = POINT_EVAL_FAILURE;
	  }
	}
      } else {
	if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoffY)) {
	  res = POINT_EVAL_BELOW_CUTOFF;
	} else {
	  res = POINT_EVAL_FAILURE;
	}
      }
    }
  } else {
    if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoffY)) {
      res = POINT_EVAL_BELOW_CUTOFF;
    } else {
      res = POINT_EVAL_FAILURE;
    }
  }
  clearChosenMpfiPtr(Y, &v_Y);
  clearChosenMpfiPtr(X, &v_X);
  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedUnivariate(mpfr_t y, node *f, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_exp_t *maxPrecUsed, node *mr) {
  mpfr_srcptr gyptr;
  int ternary;
  mp_prec_t prec;
  mpfr_t v_t;
  mpfr_t *t;
  point_eval_t resG, res;
  mp_exp_t cutoffX;
  mp_prec_t recMaxPrecUsed;

  int nodeType = f->nodeType;
  node *g = f->child1;

  /* Make compiler happy */
  gyptr = NULL;
  /* End of compiler happiness */

  /* Handle the case when g(x) = c or g(x) = x */
  if ((nodeType != LIBRARYFUNCTION) &&
      (nodeType != PROCEDUREFUNCTION)) {
    switch (accessThruMemRef(g)->nodeType) {
    case CONSTANT:
    case VARIABLE:
      switch (accessThruMemRef(g)->nodeType) {
      case CONSTANT:
	gyptr = *(accessThruMemRef(g)->value);
	break;
      case VARIABLE:
	gyptr = x;
      }
      if ( (nodeType == UNARY_BASE_FUNC) &&
           (f->baseFun->faithEvaluationOptimizedSupported) )
        ternary = f->baseFun->point_eval(y, gyptr);
      else
	return POINT_EVAL_FAILURE;

      __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, mpfr_get_prec(y));
      if (ternary == 0) return POINT_EVAL_EXACT;
      return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
      break;
    default:
      break;
    }
  }

  /* Handle the general case */
  prec = __tryFaithEvaluationOptimizedUnivariateGetRecursePrec(f, mpfr_get_prec(y), cutoff);
  if (prec < minPrec) prec = minPrec;
  cutoffX = __tryFaithEvaluationOptimizedUnivariateGetRecurseCutoff(f, cutoff, mpfr_get_prec(y));
  if ((cutoffX >= 0) || (cutoffX < mpfr_get_emin_min())) cutoffX = mpfr_get_emin_min();
  t = chooseAndInitMpfrPtr(&v_t, prec);
  recMaxPrecUsed = 0;
  resG = __tryFaithEvaluationOptimizedDoIt(*t, g, x, cutoffX, minPrec, &recMaxPrecUsed, g);
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, recMaxPrecUsed);
  switch (resG) {
  case POINT_EVAL_FAILURE:
    clearChosenMpfrPtr(t, &v_t);
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    ternary = -1;
    switch (nodeType) {
    case  UNARY_BASE_FUNC:
      if (f->baseFun->faithEvaluationOptimizedSupported)
        ternary = f->baseFun->point_eval(y, *t);
      break;
    case LIBRARYFUNCTION:
    case PROCEDUREFUNCTION:
      res = __tryFaithEvaluationOptimizedUnivariateImpreciseArg(y, f, *t, resG, cutoff, cutoffX, mr);
      clearChosenMpfrPtr(t, &v_t);
      return res;
      break;
    default:
      clearChosenMpfrPtr(t, &v_t);
      return POINT_EVAL_FAILURE;
      break;
    }
    clearChosenMpfrPtr(t, &v_t);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
  case POINT_EVAL_BELOW_CUTOFF:
    res = __tryFaithEvaluationOptimizedUnivariateImpreciseArg(y, f, *t, resG, cutoff, cutoffX, mr);
    clearChosenMpfrPtr(t, &v_t);
    return res;
    break;
  default:
    clearChosenMpfrPtr(t, &v_t);
    return POINT_EVAL_FAILURE;
    break;
  }
  return POINT_EVAL_FAILURE;
}

static inline int __tryFaithEvaluationOptimizedFuncSupported(node *f) {
  int res;

  if (f == NULL) return 0;
  switch (f->nodeType) {
  case MEMREF:
    if (f->cache->faithEvalOptFunSupportCached) {
      return f->cache->faithEvalOptFunSupportCacheResult;
    } else {
      if (f->cache->polynomialRepresentation != NULL) {
	res = 1;
      } else {
	res = __tryFaithEvaluationOptimizedFuncSupported(getMemRefChild(f));
      }
      f->cache->faithEvalOptFunSupportCacheResult = res;
      f->cache->faithEvalOptFunSupportCached = 1;
      return f->cache->faithEvalOptFunSupportCacheResult;
    }
    break;
  case VARIABLE:
  case CONSTANT:
  case PI_CONST:
    return 1;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    return (__tryFaithEvaluationOptimizedFuncSupported(f->child1) &&
	    __tryFaithEvaluationOptimizedFuncSupported(f->child2));
    break;
  case NEG:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return __tryFaithEvaluationOptimizedFuncSupported(f->child1);
    break;
  case UNARY_BASE_FUNC:
    return (f->baseFun->faithEvaluationOptimizedSupported &&
            __tryFaithEvaluationOptimizedFuncSupported(f->child1));
    break;
  default:
    return 0;
  }
  return 0;
}

static inline point_eval_t __tryFaithEvaluationOptimizedPolynomialRepresentation(mpfr_t y, polynomial_t p, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  sollya_mpfi_t v_Y, v_X;
  sollya_mpfi_t *Y, *X;
  mp_prec_t prec;
  point_eval_t res;
  mpfr_t v_t;
  mpfr_t *t;
  int tern1, tern2;

  prec = mpfr_get_prec(y) + 12;
  prec = prec + (prec >> 1);
  if (prec < minPrec) prec = minPrec;
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, prec);

  Y = chooseAndInitMpfiPtr(&v_Y, prec);
  X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));

  sollya_mpfi_set_fr(*X, x);
  polynomialEvalMpfi(*Y, p, *X);

  t = chooseAndInitMpfrPtr(&v_t, mpfr_get_prec(y));
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tern1 = mpfr_set(y, &((*Y)->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(*t, &((*Y)->right), GMP_RNDN); /* rounds to final precision */

  if (mpfr_number_p(*t) && mpfr_number_p(y)) {
    if (mpfr_equal_p(*t, y)) {
      if ((tern1 == 0) && (tern2 == 0)) {
	res = POINT_EVAL_EXACT;
      } else {
	if ((tern1 != 0) && (tern1 == tern2)) {
	  res = POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
	} else {
	  res = POINT_EVAL_CORRECTLY_ROUNDED;
	}
      }
    } else {
      if (mpfr_cmp(y, *t) < 0) {
	mpfr_nextbelow(*t);
	if (mpfr_equal_p(*t, y)) {
	  if (!sollya_mpfi_fr_in_interval(y, *Y)) {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT;
	  } else {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED;
	  }
	} else {
	  if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoff)) {
	    res = POINT_EVAL_BELOW_CUTOFF;
	  } else {
	    res = POINT_EVAL_FAILURE;
	  }
	}
      } else {
	if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoff)) {
	  res = POINT_EVAL_BELOW_CUTOFF;
	} else {
	  res = POINT_EVAL_FAILURE;
	}
      }
    }
  } else {
    if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoff)) {
      res = POINT_EVAL_BELOW_CUTOFF;
    } else {
      res = POINT_EVAL_FAILURE;
    }
  }

  clearChosenMpfrPtr(t, &v_t);
  clearChosenMpfiPtr(X, &v_X);
  clearChosenMpfiPtr(Y, &v_Y);

  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedHooksInner(mpfr_t y, eval_hook_t *hook, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed, int tight, mp_prec_t extraPrec, int *hookWorked, mp_prec_t *estimatePrec, int *knowEstimatePrec) {
  sollya_mpfi_t v_Y, v_X;
  sollya_mpfi_t *Y, *X;
  mp_prec_t prec;
  point_eval_t res;
  mpfr_t v_t;
  mpfr_t *t;
  int tern1, tern2;
  int hookRes;
  mp_exp_t mER, mEL, mE, tE;

  *hookWorked = 0;

  if (hook == NULL) return POINT_EVAL_FAILURE;

  prec = mpfr_get_prec(y) + 20;
  prec = prec + (prec >> 1);
  if (prec < minPrec) prec = minPrec;
  prec += extraPrec;
  __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, prec);

  Y = chooseAndInitMpfiPtr(&v_Y, (tight ? (mpfr_get_prec(y) + 10 + extraPrec) : ((prec > (mpfr_get_prec(y) + 10 + extraPrec)) ? prec : (mpfr_get_prec(y) + 10 + extraPrec))));
  X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));

  sollya_mpfi_set_fr(*X, x);
  hookRes = evaluateWithEvaluationHook(*Y, *X, prec, tight, hook);
  if (!hookRes) {
    clearChosenMpfiPtr(X, &v_X);
    clearChosenMpfiPtr(Y, &v_Y);
    return POINT_EVAL_FAILURE;
  }
  *hookWorked = 1;

  t = chooseAndInitMpfrPtr(&v_t, 12);
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ((estimatePrec != NULL) &&
      (knowEstimatePrec != NULL)) {
    *knowEstimatePrec = 0;
    if (!(sollya_mpfi_has_zero(*Y) ||
	  sollya_mpfi_has_infinity(*Y) ||
	  sollya_mpfi_has_nan(*Y) ||
	  sollya_mpfi_is_empty(*Y))) {
      /* Try to get an idea of the accuracy of the result */
      mpfr_sub(*t, &((*Y)->right), &((*Y)->left), GMP_RNDU);
      if (mpfr_number_p(*t) && (!mpfr_zero_p(*t))) {
	mEL = mpfr_get_exp(&((*Y)->left));
	mER = mpfr_get_exp(&((*Y)->right));
	mE = mEL; if (mER < mE) mE = mER;
	tE = mpfr_get_exp(*t);
	if (mE >= tE) {
	  *estimatePrec = (mp_prec_t) (mE - tE);
	  *knowEstimatePrec = 1;
	}
      }
    }
  }

  mpfr_set_prec(*t, mpfr_get_prec(y));
  tern1 = mpfr_set(y, &((*Y)->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(*t, &((*Y)->right), GMP_RNDN); /* rounds to final precision */

  if (mpfr_number_p(*t) && mpfr_number_p(y)) {
    if (mpfr_equal_p(*t, y)) {
      if ((tern1 == 0) && (tern2 == 0)) {
	res = POINT_EVAL_EXACT;
      } else {
	if ((tern1 != 0) && (tern1 == tern2)) {
	  res = POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
	} else {
	  res = POINT_EVAL_CORRECTLY_ROUNDED;
	}
      }
    } else {
      if (mpfr_cmp(y, *t) < 0) {
	mpfr_nextbelow(*t);
	if (mpfr_equal_p(*t, y)) {
	  if (!sollya_mpfi_fr_in_interval(y, *Y)) {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT;
	  } else {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED;
	  }
	} else {
	  if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoff)) {
	    res = POINT_EVAL_BELOW_CUTOFF;
	  } else {
	    res = POINT_EVAL_FAILURE;
	  }
	}
      } else {
	if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoff)) {
	  res = POINT_EVAL_BELOW_CUTOFF;
	} else {
	  res = POINT_EVAL_FAILURE;
	}
      }
    }
  } else {
    if ((!(sollya_mpfi_has_nan(*Y) || sollya_mpfi_has_infinity(*Y))) && (sollya_mpfi_max_exp(*Y) < cutoff)) {
      res = POINT_EVAL_BELOW_CUTOFF;
    } else {
      res = POINT_EVAL_FAILURE;
    }
  }

  clearChosenMpfrPtr(t, &v_t);
  clearChosenMpfiPtr(X, &v_X);
  clearChosenMpfiPtr(Y, &v_Y);

  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedHooks(mpfr_t y, eval_hook_t *hook, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed, int *hookGivesSomeAccuracy) {
  int hookWorked;
  point_eval_t res;
  mp_prec_t extraPrec, estimatePrec, firstPrec;
  int knowEstimatePrec, knowFirstPrec;

  *hookGivesSomeAccuracy = 0;

  if (hook == NULL) return POINT_EVAL_FAILURE;

  extraPrec = 0;
  firstPrec = 0;
  hookWorked = 0;
  res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 1, extraPrec, &hookWorked, NULL, NULL);
  if (res != POINT_EVAL_FAILURE) return res;

  extraPrec = 0;
  hookWorked = 0;
  knowEstimatePrec = 0;
  res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 0, extraPrec, &hookWorked, &estimatePrec, &knowEstimatePrec);
  if (res != POINT_EVAL_FAILURE) return res;
  if (!hookWorked) return res;
  *hookGivesSomeAccuracy = 1;
  if (knowEstimatePrec) {
    knowFirstPrec = 1;
    firstPrec = estimatePrec;
  }


  extraPrec = 10;
  hookWorked = 0;
  knowEstimatePrec = 0;
  res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 0, extraPrec, &hookWorked, &estimatePrec, &knowEstimatePrec);
  if (res != POINT_EVAL_FAILURE) return res;
  if (!hookWorked) return res;
  if (knowEstimatePrec) {
    if (knowFirstPrec) {
      if (estimatePrec <= firstPrec) return res;
      firstPrec = estimatePrec;
    } else {
      knowFirstPrec = 1;
      firstPrec = estimatePrec;
    }
  }


  extraPrec = 20;
  hookWorked = 0;
  knowEstimatePrec = 0;
  res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 0, extraPrec, &hookWorked, &estimatePrec, &knowEstimatePrec);
  if (res != POINT_EVAL_FAILURE) return res;
  if (!hookWorked) return res;
  if (knowEstimatePrec) {
    if (knowFirstPrec) {
      if (estimatePrec <= firstPrec) return res;
      firstPrec = estimatePrec;
    } else {
      knowFirstPrec = 1;
      firstPrec = estimatePrec;
    }
  }


  extraPrec = 40;
  hookWorked = 0;
  knowEstimatePrec = 0;
  res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 0, extraPrec, &hookWorked, &estimatePrec, &knowEstimatePrec);
  if (res != POINT_EVAL_FAILURE) return res;
  if (!hookWorked) return res;
  if (knowEstimatePrec) {
    if (knowFirstPrec) {
      if (estimatePrec <= firstPrec) return res;
      firstPrec = estimatePrec;
    } else {
      knowFirstPrec = 1;
      firstPrec = estimatePrec;
    }
  }


  extraPrec = 60;
  hookWorked = 0;
  knowEstimatePrec = 0;
  res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 0, extraPrec, &hookWorked, &estimatePrec, &knowEstimatePrec);
  if (res != POINT_EVAL_FAILURE) return res;
  if (!hookWorked) return res;
  if (knowEstimatePrec) {
    if (knowFirstPrec) {
      if (estimatePrec <= firstPrec) return res;
      firstPrec = estimatePrec;
    } else {
      knowFirstPrec = 1;
      firstPrec = estimatePrec;
    }
  }


  extraPrec = mpfr_get_prec(y) + 10;
  if (extraPrec < 80) extraPrec = 80;
  while (extraPrec <= 5 * mpfr_get_prec(y)) {
    hookWorked = 0;
    knowEstimatePrec = 0;
    res = __tryFaithEvaluationOptimizedHooksInner(y, hook, x, cutoff, minPrec, maxPrecUsed, 0, extraPrec, &hookWorked, &estimatePrec, &knowEstimatePrec);
    if (res != POINT_EVAL_FAILURE) return res;
    if (!hookWorked) return res;
    if (knowEstimatePrec) {
      if (knowFirstPrec) {
	if (estimatePrec <= firstPrec) return res;
	firstPrec = estimatePrec;
      } else {
	knowFirstPrec = 1;
	firstPrec = estimatePrec;
      }
    }
    extraPrec += ((extraPrec >> 2) > 10 ? (extraPrec >> 2) : 10);
  }

  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedDeducedLowerPrecResult(mpfr_t y, mpfr_t x, point_eval_t approx, mp_exp_t approxCutoff, mp_exp_t cutoff) {
  int ternary, tern1, tern2;
  sollya_mpfi_t v_X;
  sollya_mpfi_t *X;
  mpfr_t v_t;
  mpfr_t *t;
  point_eval_t res;

  /* Make compiler happy */
  X = NULL;
  /* End of compiler happiness */

  if (approx == POINT_EVAL_FAILURE) return POINT_EVAL_FAILURE;
  if (!mpfr_number_p(x)) return POINT_EVAL_FAILURE;
  if (approx == POINT_EVAL_EXACT) {
    ternary = mpfr_set(y, x, GMP_RNDN);
    if (ternary == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
  }

  switch (approx) {
  case POINT_EVAL_FAILURE:
    return POINT_EVAL_FAILURE;
    break;
  case POINT_EVAL_EXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));
    sollya_mpfi_set_fr(*X, x);
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x) + 1);
    sollya_mpfi_set_fr(*X, x);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));
    sollya_mpfi_set_fr(*X, x);
    sollya_mpfi_blow_1ulp(*X);
    break;
  case POINT_EVAL_BELOW_CUTOFF:
    X = chooseAndInitMpfiPtr(&v_X, mpfr_get_prec(x));
    sollya_mpfi_interv_si_2exp(*X, -1, approxCutoff, 1, approxCutoff);
    break;
  }

  t = chooseAndInitMpfrPtr(&v_t, mpfr_get_prec(y));
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tern1 = mpfr_set(y, &((*X)->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(*t, &((*X)->right), GMP_RNDN); /* rounds to final precision */

  if (mpfr_number_p(*t) && mpfr_number_p(y)) {
    if (mpfr_equal_p(*t, y)) {
      if ((tern1 == 0) && (tern2 == 0)) {
	res = POINT_EVAL_EXACT;
      } else {
	if ((tern1 != 0) && (tern1 == tern2)) {
	  res = POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
	} else {
	  res = POINT_EVAL_CORRECTLY_ROUNDED;
	}
      }
    } else {
      if (mpfr_cmp(y, *t) < 0) {
	mpfr_nextbelow(*t);
	if (mpfr_equal_p(*t, y)) {
	  if (!sollya_mpfi_fr_in_interval(y, *X)) {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT;
	  } else {
	    res = POINT_EVAL_FAITHFULLY_ROUNDED;
	  }
	} else {
	  if ((!(sollya_mpfi_has_nan(*X) || sollya_mpfi_has_infinity(*X))) && (sollya_mpfi_max_exp(*X) < cutoff)) {
	    res = POINT_EVAL_BELOW_CUTOFF;
	  } else {
	    res = POINT_EVAL_FAILURE;
	  }
	}
      } else {
	if ((!(sollya_mpfi_has_nan(*X) || sollya_mpfi_has_infinity(*X))) && (sollya_mpfi_max_exp(*X) < cutoff)) {
	  res = POINT_EVAL_BELOW_CUTOFF;
	} else {
	  res = POINT_EVAL_FAILURE;
	}
      }
    }
  } else {
    if ((!(sollya_mpfi_has_nan(*X) || sollya_mpfi_has_infinity(*X))) && (sollya_mpfi_max_exp(*X) < cutoff)) {
      res = POINT_EVAL_BELOW_CUTOFF;
    } else {
      res = POINT_EVAL_FAILURE;
    }
  }

  clearChosenMpfrPtr(t, &v_t);
  clearChosenMpfiPtr(X, &v_X);
  return res;
}

static inline point_eval_t __tryFaithEvaluationOptimizedDoIt(mpfr_t y, node *f, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed, node *mr) {
  point_eval_t res;
  int hookGivesSomeAccuracy;

  switch (f->nodeType) {
  case MEMREF:
    if ((f->cache->pointEvalCacheX != NULL) &&
	(f->cache->pointEvalCacheY != NULL) &&
	(f->cache->pointEvalCacheResultType != POINT_EVAL_FAILURE) &&
	(f->cache->pointEvalCacheCutoff <= cutoff) &&
	(mpfr_get_prec(*(f->cache->pointEvalCacheY)) >= mpfr_get_prec(y)) &&
	mpfr_number_p(x) &&
	mpfr_number_p(*(f->cache->pointEvalCacheX)) &&
	mpfr_number_p(*(f->cache->pointEvalCacheY)) &&
	mpfr_equal_p(*(f->cache->pointEvalCacheX), x)) {
      if (mpfr_get_prec(*(f->cache->pointEvalCacheY)) == mpfr_get_prec(y)) {
	mpfr_set(y, *(f->cache->pointEvalCacheY), GMP_RNDN); /* exact */
	return f->cache->pointEvalCacheResultType;
      } else {
	res = __tryFaithEvaluationOptimizedDeducedLowerPrecResult(y, *(f->cache->pointEvalCacheY), f->cache->pointEvalCacheResultType, f->cache->pointEvalCacheCutoff, cutoff);
	if (res != POINT_EVAL_FAILURE) return res;
      }
    }
    if ((f->cache->polynomialRepresentation != NULL) && (f->child1 == NULL)) {
      res = __tryFaithEvaluationOptimizedPolynomialRepresentation(y, f->cache->polynomialRepresentation, x, cutoff, minPrec, maxPrecUsed);
    } else {
      hookGivesSomeAccuracy = 0;
      res = __tryFaithEvaluationOptimizedHooks(y, f->cache->evaluationHook, x, cutoff, minPrec, maxPrecUsed, &hookGivesSomeAccuracy);
      if (hookGivesSomeAccuracy) return res;
      if (res == POINT_EVAL_FAILURE) {
	res = __tryFaithEvaluationOptimizedDoIt(y, getMemRefChild(f), x, cutoff, minPrec, maxPrecUsed, f);
      }
    }
    if ((res != POINT_EVAL_FAILURE) && (f->libFunDeriv >= 2) && mpfr_number_p(y)) {
      if (f->cache->pointEvalCacheX == NULL) {
	f->cache->pointEvalCacheX = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(f->cache->pointEvalCacheX), mpfr_get_prec(x));
      }
      if (f->cache->pointEvalCacheY == NULL) {
	f->cache->pointEvalCacheY = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(f->cache->pointEvalCacheY), mpfr_get_prec(y));
      }
      mpfr_set_prec(*(f->cache->pointEvalCacheX), mpfr_get_prec(x));
      mpfr_set_prec(*(f->cache->pointEvalCacheY), mpfr_get_prec(y));
      mpfr_set(*(f->cache->pointEvalCacheX), x, GMP_RNDN); /* exact */
      mpfr_set(*(f->cache->pointEvalCacheY), y, GMP_RNDN); /* exact */
      f->cache->pointEvalCacheCutoff = cutoff;
      f->cache->pointEvalCacheResultType = res;
    }
    return res;
    break;
  case VARIABLE:
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, mpfr_get_prec(y));
    if (mpfr_set(y, x, GMP_RNDN) == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    break;
  case CONSTANT:
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, mpfr_get_prec(y));
    if (mpfr_set(y, *(f->value), GMP_RNDN) == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    break;
  case PI_CONST:
    __tryFaithEvaluationOptimizedUpdateMaxPrec(maxPrecUsed, mpfr_get_prec(y));
    if (mpfr_const_pi(y, GMP_RNDN) == 0) return POINT_EVAL_EXACT;
    return POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT;
    break;
  case ADD:
  case SUB:
    return __tryFaithEvaluationOptimizedAddSub(y, (f->nodeType == SUB), f->child1, f->child2, x, cutoff, minPrec, maxPrecUsed);
    break;
  case MUL:
  case DIV:
    return __tryFaithEvaluationOptimizedMulDiv(y, (f->nodeType == DIV), f->child1, f->child2, x, cutoff, minPrec, maxPrecUsed);
    break;
  case POW:
    return __tryFaithEvaluationOptimizedPow(y, f->child1, f->child2, x, cutoff, minPrec, maxPrecUsed);
    break;
  case NEG:
    /* So terribly trivial that we do it here */
    res = __tryFaithEvaluationOptimizedDoIt(y, f->child1, x, cutoff, minPrec, maxPrecUsed, f->child1);
    if (res != POINT_EVAL_FAILURE) {
      mpfr_neg(y, y, GMP_RNDN); /* exact */
    }
    return res;
    break;
  case UNARY_BASE_FUNC:
    if (f->baseFun->faithEvaluationOptimizedSupported)
      return __tryFaithEvaluationOptimizedUnivariate(y, f, x, cutoff, minPrec, maxPrecUsed, mr);
    break;
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return __tryFaithEvaluationOptimizedUnivariate(y, f, x, cutoff, minPrec, maxPrecUsed, mr);
    break;
  default:
    return POINT_EVAL_FAILURE;
  }
  return POINT_EVAL_FAILURE;
}

static inline int __tryFaithEvaluationOptimized(int *retVal, mpfr_t y, node *func, mpfr_t x, mp_exp_t cutoff, mp_prec_t minPrec, mp_prec_t *maxPrecUsed) {
  point_eval_t res;

  /* Refuse work if x and y are the same MPFR variable */
  if (((void *) x) == ((void *) y)) return 0;

  /* Refuse work for func = NULL */
  if (func == NULL) return 0;

  /* Refuse work if x is not a real */
  if (!mpfr_number_p(x)) return 0;

  /* Call inner function and translate the success information */
  res = __tryFaithEvaluationOptimizedDoIt(y, func, x, cutoff, minPrec, maxPrecUsed, func);
  switch (res) {
  case POINT_EVAL_FAILURE:
    return 0;
    break;
  case POINT_EVAL_EXACT:
    if (!mpfr_number_p(y)) return 0;
    *retVal = 4;
    return 1;
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED:
    if (!mpfr_number_p(y)) return 0;
    *retVal = 6;
    return 1;
    break;
  case POINT_EVAL_CORRECTLY_ROUNDED_PROVEN_INEXACT:
    if (!mpfr_number_p(y)) return 0;
    *retVal = 7;
    return 1;
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED:
    if (!mpfr_number_p(y)) return 0;
    *retVal = 1;
    return 1;
    break;
  case POINT_EVAL_FAITHFULLY_ROUNDED_PROVEN_INEXACT:
    if (!mpfr_number_p(y)) return 0;
    *retVal = 5;
    return 1;
    break;
  case POINT_EVAL_BELOW_CUTOFF:
    if (!mpfr_number_p(y)) return 0;
    *retVal = 2;
    return 1;
  default:
    return 0;
    break;
  }

  return 0;
}


/* A performance optimized short-circuit variant to the function
   evaluateFaithfulWithCutOffFastInternalImplementation defined below.

   This function is not required to work in all cases. If it works, it
   has to return a non-zero value and assign the evaluation result to
   y and to set retVal to the "return values" defined (in comment)
   below.

   This function is supposed not to allocate memory when allocation
   can be avoided or amortized.

*/
static inline int firstTryEvaluateFaithfulWithCutOffFastInternalImplementation(int *retVal, mpfr_t y, node *func, mpfr_t x, mp_prec_t startprec, mp_exp_t *cutoff) {
  mp_prec_t pX, pY;
  int tern1, tern2;
  mp_exp_t myCutoff;
  mp_prec_t maxPrecUsed;

  /* Refuse work if x and y are the same MPFR variable */
  if (((void *) x) == ((void *) y)) return 0;

  /* Try to use a still more optimized function */
  if (__tryFaithEvaluationOptimizedFuncSupported(func)) {
    if (cutoff == NULL) {
      myCutoff = mpfr_get_emin_min();
    } else {
      myCutoff = *cutoff;
    }
    maxPrecUsed = 0;
    if (__tryFaithEvaluationOptimized(retVal, y, func, x, myCutoff, startprec, &maxPrecUsed)) {

      return 1;
    }
  }

  /* Get the precisions of the x and y arguments */
  pX = mpfr_get_prec(x);
  pY = mpfr_get_prec(y);

  /* Check if we have any chance to get a faithful rounding in y with
     evaluations at startprec bits
  */
  if (pY > startprec) return 0;

  /* Initialize interval variables for x and for the evaluation of f
     over this interval x

     If we are in a recursive call, simply fail and let the usual
     evaluation function do the work. Recursive calls are only
     possible in funny cases when a procedure function or a library
     function uses pointwise faithfully rounded evaluations to
     evaluate itself. If we are not in a recursive call, get the token
     to prevent other recursive calls.
  */
  if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used) return 0;
  __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 1;
  if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized) {
    sollya_mpfi_set_prec(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x, pX);
  } else {
    sollya_mpfi_init2(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x, pX);
    __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized = 1;
  }
  if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized) {
    sollya_mpfi_set_prec(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y, startprec);
  } else {
    sollya_mpfi_init2(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y, startprec);
    __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized = 1;
  }

  /* Set the interval x to the point x */
  sollya_mpfi_set_fr(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x, x); /* exact, same precision */

  /* Perform an interval evaluation of f over the interval x, yielding an interval y */
  evaluateIntervalInternalFast(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y, func, NULL, __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x, 1, 0, cutoff, NULL);

  /* If the interval we got is not bounded, we simply indicate failure
     and let the usual evaluation function do the work.
  */
  if (!sollya_mpfi_bounded_p(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y)) {
    __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 0;
    return 0;
  }

  /* Here, the interval y bounding f(x) contains nothing but real
     numbers.

     We now check if correct rounding is already possible.

     We DO NOT check if faithful rounding is already possible as this
     might require reevaluation (see the strategy in the function
     evaluateFaithfulWithCutOffFastInternalImplementation below).

     The test whether or not correct rounding is possible
     requires a temporary MPFR variable with *exactly* the
     same precision as the MPFR y variable.

  */
  if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized) {
    mpfr_set_prec(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp, pY);
  } else {
    mpfr_init2(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp, pY);
    __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized = 1;
  }

  /* Round both endpoints of the interval y to the nearest. If both roundings
     yield the same value and are real numbers, correct rounding is possible.
  */
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tern1 = mpfr_set(y, &(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y->left), GMP_RNDN); /* rounds to final precision */
  tern2 = mpfr_set(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp, &(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y->right), GMP_RNDN); /* rounds to final precision */

  /* Check if both roundings are real and equal */
  if (mpfr_number_p(y) &&
      mpfr_number_p(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp) &&
      mpfr_equal_p(y, __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp)) {
    /* Here, we could determine a correct rounding.
       The numerical result variable y has already been set.
       Set the "return value" to "correctly rounded proven inexact" (7),
       "correctly rounded might be exact" (6) or "exact" (4).
       Return "success".

       * The result is (proven) exact if none of the roundings was
       inexact.
       * The result is proven inexact if none of the roundings was
       exact and both roundings went in the same direction.
       * Otherwise, we don't know.
       */
    if ((tern1 == 0) && (tern2 == 0)) {
      /* the result was "exact" */
      *retVal = 4;
    } else {
      if (tern1 * tern2 > 0) {
	/* the result was proven "inexact" */
	*retVal = 7;
      } else {
	/* We do not know if the result is exact or not. */
	*retVal = 6;
      }
    }
    __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 0;
    return 1;
  }

  /* We could not easily determine a faithful (correct) rounding */
  __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 0;
  return 0;
}

/* Return values:

   0 -> evaluation did not allow for a faithful evaluation nor to get
   below the threshold but the proof interval was a closed subset
   of the reals or an exact infinity

   1 -> evaluation gave a faithful result, we do not know if it is
   exact or inexact

   2 -> evaluation got the proof interval below the cutoff

   3 -> evaluation did not allow to conclude on a function as it
   evaluated to NaN or Inf

   4 -> evaluation gave a faithfully resp. correctly rounded result
   and we know that it is exact

   5 -> evaluation gave a faithful result and we know that it is
   inexact

   6 -> evaluation gave a correctly rounded result, we do not know if
   it is inexact

   7 -> evaluation gave a correctly rounded result and we know that it
   is inexact

*/
static inline int evaluateFaithfulWithCutOffFastInternalImplementation(mpfr_t result, node *func, node *deriv, mpfr_t x, mpfr_t cutoff, mp_prec_t startprecOrig, node *altX) {
  mp_prec_t p, prec, pTemp, pRes;
  sollya_mpfi_t yI, xI, cutoffI, dummyI;
  int okay;
  mpfr_t resUp, resDown, resUpFaith, resDownFaith, resDownTemp;
  mpfr_t cutoffLeft, cutoffRight;
  mpfr_t yILeft, yIRight, yILeftCheck, yIRightCheck;
  int testCutOff;
  int correctlyRounded;
  mp_prec_t startprec;
  int precisionIncreased;
  int res;
  int retry;
  mp_exp_t expoCutoff;
  int lazyHookUsed;
  int noLazyHooks, forceNoLazyHooks, numberLazyHookUsed;
  mp_prec_t lazyStartPrec;

  /* Determine some sensible starting precision */
  startprec = startprecOrig;
  pRes = mpfr_get_prec(result);
  if ((mpfr_sgn(cutoff) == 0) || (!mpfr_number_p(cutoff))) {
    /* The cutoff is zero or Nan */
    expoCutoff = mpfr_get_emin_min();

    /* If the cutoff is zero or NaN, we actually need at least as much
       working precision as the precision of the result we are to
       produce.
    */
    if (startprec < pRes + 10) startprec = pRes + 10;
  } else {
    /* The cutoff is a real and not zero */
    expoCutoff = mpfr_get_exp(cutoff);
    if (expoCutoff >= -1) {
      expoCutoff = -2;
    }
    expoCutoff--;
    if ((expoCutoff >= 0) || (expoCutoff < mpfr_get_emin_min())) expoCutoff = mpfr_get_emin_min();
  }

  /* Try a short-circuit for common case evaluations */
  if ((altX == NULL) && (((void *) x) != ((void *) result))) {
    if (firstTryEvaluateFaithfulWithCutOffFastInternalImplementation(&res, result, func, x, startprec, &expoCutoff)) {
      return res;
    }
  }

  /* Check if we have a constant expression to evaluate at and if so,
     check if it is constant
  */
  if ((altX != NULL) && (!isConstant(altX))) {
    /* Here the alternate abscissa expression is not constant */
    mpfr_set_nan(result);
    return 0;
  }

  /* We need a dummy interval for some evaluations */
  sollya_mpfi_init2(dummyI,12);
  sollya_mpfi_set_ui(dummyI,1);

  /* We test the cutoff only if it is not zero */
  testCutOff = 1;

  prec = mpfr_get_prec(cutoff);
  mpfr_init2(cutoffLeft, prec);
  mpfr_init2(cutoffRight, prec);
  sollya_mpfi_init2(cutoffI,prec);
  mpfr_abs(cutoffRight,cutoff,GMP_RNDU);

  /* We test the cutoff only if it is not zero */
  mpfr_set_si(cutoffLeft,0,GMP_RNDN);
  if (mpfr_cmp(cutoffRight,cutoffLeft) == 0) testCutOff = 0;

  mpfr_neg(cutoffLeft,cutoffRight,GMP_RNDD);
  sollya_mpfi_interv_fr(cutoffI,cutoffLeft,cutoffRight);
  mpfr_clear(cutoffLeft);
  mpfr_clear(cutoffRight);

  /* Initialize the variables that hold the possible roundings */
  prec = mpfr_get_prec(result);
  mpfr_init2(resUp,prec);
  mpfr_init2(resDown,prec);
  mpfr_init2(resUpFaith,prec);
  mpfr_init2(resDownFaith,prec);
  mpfr_init2(resDownTemp,prec);

  /* Determine a starting precision */
  if (startprec > prec) prec = startprec;

  /* Use up a little more memory with the first malloc
     The starting subsequent mpf*_set_prec will not malloc
  */
  sollya_mpfi_init2(yI,startprec*16);
  mpfr_init2(yILeft,startprec*16);
  mpfr_init2(yIRight,startprec*16);

  /* Initialize an interval for the abscissa point we evaluate at */
  if (altX == NULL) {
    /* If we have no alternate abscissa point, copy x into an interval
       with its own precision
    */
    p = mpfr_get_prec(x);
    sollya_mpfi_init2(xI,p);
    sollya_mpfi_interv_fr(xI,x,x);
  } else {
    /* Otherwise, initialize an interval with a little more than the
       starting precision and then set the precision to the starting
       precision. This allows for the first subsequent mpf*_set_prec
       not to malloc which saves a couple of mallocs and copies.
    */
    sollya_mpfi_init2(xI,startprec * 16);
    sollya_mpfi_set_prec(xI,startprec);
  }

  /* Start the rounding loop */
  p=startprec;
  okay = 0;
  retry = 0;
  numberLazyHookUsed = 0;
  forceNoLazyHooks = 0;
  lazyStartPrec = p;
  while (p < prec * 512) {
    correctlyRounded = 0;
    lazyHookUsed = 0;
    noLazyHooks = (forceNoLazyHooks || (retry < 4));

    sollya_mpfi_set_prec(yI,p);

    /* If we evaluate at a constant expression instead of a point,
       evaluate the constant expression to the current (loop)
       precision.
    */
    if (altX != NULL) {
      sollya_mpfi_set_prec(xI,p);
      evaluateIntervalInternalFast(xI, altX, NULL, dummyI, (retry < 2), noLazyHooks, &expoCutoff, &lazyHookUsed);
    }

    mpfr_set_prec(yILeft,p);
    mpfr_set_prec(yIRight,p);
    evaluateIntervalInternalFast(yI, func, deriv, xI, (retry < 2), noLazyHooks, &expoCutoff, &lazyHookUsed);
    sollya_mpfi_get_left(yILeft,yI);
    sollya_mpfi_get_right(yIRight,yI);
    mpfr_set(resDown,yILeft,GMP_RNDN);
    mpfr_set(resUp,yIRight,GMP_RNDN);
    /*
      No, two sollya_mpfi_get_left/right on resDown and resUp directly do not work:

      ----|-----------|-----------|----
      |--yI-|
      left        RN          right
    */
    if ((mpfr_number_p(resDown)) &&
	(mpfr_number_p(resUp)) &&
	(sollya_mpfi_bounded_p(yI))) {
      /* Check if correct rounding is already possible */
      if (mpfr_cmp(resDown,resUp) == 0) {
	/* Correct rounding is possible, resDown = resUp is that
	   correct rounding
	*/
	okay = 1;
	correctlyRounded = 1;
      } else {
	/* Check if faithful rounding is already possible because the
           proof interval is between two consecutive floating-point
           numbers in the target format.
	*/
	sollya_mpfi_get_left(resDownFaith,yI);
	sollya_mpfi_get_right(resUpFaith,yI);
	mpfr_set(resDownTemp, resDownFaith, GMP_RNDN); /* exact, same precision */
	mpfr_nextabove(resDownTemp);
	if (mpfr_number_p(resDownFaith) &&
	    mpfr_number_p(resUpFaith) &&
	    (mpfr_cmp(resDownTemp,resUpFaith) == 0)) {
	  /* A faithful rounding is already possible because the
	     proof interval is between two consecutive floating-point
	     numbers in the target format.

	     As the rounding process implemented here will start with
             some small precision (yielding large proof intervals at
             first) and increase that precision more and more until
             correct or faithful rounding becomes possible and since
             faithful rounding is easier than correct rounding, there
             are more chances to fall into this case of faithful
             rounding than into the case for correct rounding.

	     In order to allow for correct rounding in most cases
	     anyway, we are hence going to try to recompute the proof
	     interval at a slightly higher precision to see if we
	     can't get correct rounding "for just one dollar more."

          */
	  pTemp = p + 10;

	  /* Recompute at that slightly higher precision */
	  sollya_mpfi_set_prec(yI,pTemp);

	  /* If we evaluate at a constant expression instead of a point,
	     evaluate the constant expression to the current (loop)
	     precision.
	  */
	  if (altX != NULL) {
	    sollya_mpfi_set_prec(xI,pTemp);
	    evaluateIntervalInternalFast(xI, altX, NULL, dummyI, (retry < 2), (retry < 4), &expoCutoff, &lazyHookUsed);
	  }

	  mpfr_set_prec(yILeft,pTemp);
	  mpfr_set_prec(yIRight,pTemp);
	  evaluateIntervalInternalFast(yI, func, deriv, xI, (retry < 2), (retry < 4), &expoCutoff, &lazyHookUsed);
	  sollya_mpfi_get_left(yILeft,yI);
	  sollya_mpfi_get_right(yIRight,yI);
	  mpfr_set(resDown,yILeft,GMP_RNDN);
	  mpfr_set(resUp,yIRight,GMP_RNDN);

	  /* Check if we have a correct rounding now. */
	  if (mpfr_number_p(resDown) &&
	      mpfr_number_p(resUp) &&
	      sollya_mpfi_bounded_p(yI) &&
	      (mpfr_cmp(resDown,resUp) == 0)) {
	    /* We did find a correct rounding now, which we're gonna
	       use.
	    */
	    okay = 1;
	    correctlyRounded = 1;
	  } else {
	    /* We couldn't find a correct rounding with a
	       recomputation at a slightly higher precision. So we'll
	       go for our faithful rounding.
	    */
	    okay = 1;
	    mpfr_set(resDown, resDownFaith, GMP_RNDN); /* exact, same precision */
	    mpfr_set(resUp, resUpFaith, GMP_RNDN); /* exact, same precision */
	  }
	} else {
	  /* Check if faithful rounding is already possible because
	     the proof interval contains a single floating-point
	     number in the target format.
	  */
	  mpfr_set(resDownFaith,yILeft,GMP_RNDU);
	  mpfr_set(resUpFaith,yIRight,GMP_RNDD);
	  if (mpfr_number_p(resDownFaith) &&
	      mpfr_number_p(resUpFaith) &&
	      (mpfr_cmp(resDownFaith,resUpFaith) == 0)) {
	    /* A faithful rounding is already possible because the
	       proof interval contains a single floating-point number
	       in the target format.

	       As the rounding process implemented here will start with
	       some small precision (yielding large proof intervals at
	       first) and increase that precision more and more until
	       correct or faithful rounding becomes possible and since
	       faithful rounding is easier than correct rounding, there
	       are more chances to fall into this case of faithful
	       rounding than into the case for correct rounding.

	       In order to allow for correct rounding in most cases
	       anyway, we are hence going to try to recompute the proof
	       interval at a slightly higher precision to see if we
	       can't get correct rounding "for just one dollar more."

	    */
	    pTemp = p + 10;

	    /* Recompute at that slightly higher precision */
	    sollya_mpfi_set_prec(yI,pTemp);

	    /* If we evaluate at a constant expression instead of a point,
	       evaluate the constant expression to the current (loop)
	       precision.
	    */
	    if (altX != NULL) {
	      sollya_mpfi_set_prec(xI,pTemp);
	      evaluateIntervalInternalFast(xI, altX, NULL, dummyI, (retry < 2), (retry < 4), &expoCutoff, &lazyHookUsed);
	    }

	    mpfr_set_prec(yILeft,pTemp);
	    mpfr_set_prec(yIRight,pTemp);
	    evaluateIntervalInternalFast(yI, func, deriv, xI, (retry < 2), (retry < 4), &expoCutoff, &lazyHookUsed);
	    sollya_mpfi_get_left(yILeft,yI);
	    sollya_mpfi_get_right(yIRight,yI);
	    mpfr_set(resDown,yILeft,GMP_RNDN);
	    mpfr_set(resUp,yIRight,GMP_RNDN);

	    /* Check if we have a correct rounding now. */
	    if (mpfr_number_p(resDown) &&
		mpfr_number_p(resUp) &&
		sollya_mpfi_bounded_p(yI) &&
		(mpfr_cmp(resDown,resUp) == 0)) {
	      /* We did find a correct rounding now, which we're gonna
		 use.
	      */
	      okay = 1;
	      correctlyRounded = 1;
	    } else {
	      /* We couldn't find a correct rounding with a
		 recomputation at a slightly higher precision. So we'll
		 go for our faithful rounding.
	      */
	      okay = 1;
	      mpfr_set(resDown, resDownFaith, GMP_RNDN); /* exact, same precision */
	      mpfr_set(resUp, resUpFaith, GMP_RNDN); /* exact, same precision */
	    }
	  } else {
	    /* Check if the proof interval is already contained in the
	       cutoff interval
	    */
	    if (testCutOff && (okay == 0)) {
	      if (sollya_mpfi_is_inside(yI,cutoffI)) {
		sollya_mpfi_mid(resUp,yI);
		mpfr_set(resDown, resUp, GMP_RNDN); /* exact, same precision */
		okay = 2;
	      }
	    }
	  }
	}
      }
    } else {
      if (accessThruMemRef(func)->nodeType == DIV) {
	if ((accessThruMemRef(func)->child1->nodeType == MEMREF) &&
	    (accessThruMemRef(func)->child2->nodeType == MEMREF)) {
	  if ((accessThruMemRef(func)->child1->cache->evalCacheX != NULL) &&
	      (accessThruMemRef(func)->child1->cache->evalCacheY != NULL) &&
	      (accessThruMemRef(func)->child2->cache->evalCacheX != NULL) &&
	      (accessThruMemRef(func)->child2->cache->evalCacheY != NULL)) {
	    if (sollya_mpfi_is_point_and_real(xI) &&
		sollya_mpfi_equal_p(*(accessThruMemRef(func)->child1->cache->evalCacheX), xI) &&
		sollya_mpfi_equal_p(*(accessThruMemRef(func)->child2->cache->evalCacheX), xI)) {
	      if (sollya_mpfi_is_zero(*(accessThruMemRef(func)->child1->cache->evalCacheY)) &&
		  sollya_mpfi_is_zero(*(accessThruMemRef(func)->child2->cache->evalCacheY))) {
		break;
	      }
	    }
	  }
	}
      }
    }
    if (okay > 0) break;

    /* Now sensibly increase the working precision for the next
       round
    */
    retry++;
    precisionIncreased = 0;
    if (lazyHookUsed) {
      /* In the case when the result we got was computed with a lazy
	 hook evaluation, it makes sense to increase the precision
	 only in very small steps but not eternally
      */
      p = p + (7 + numberLazyHookUsed + ((numberLazyHookUsed * numberLazyHookUsed) >> 4)) + (p >> 6);
      precisionIncreased = 1;
      numberLazyHookUsed++;
      if (numberLazyHookUsed < 2) {
	lazyStartPrec = p;
      }
      if (numberLazyHookUsed > 15) {
	if ((p > (lazyStartPrec + ((lazyStartPrec >> 2) + 1))) || (numberLazyHookUsed > 1000)) {
	  forceNoLazyHooks = 1;
	}
      }
    } else {
      numberLazyHookUsed = 0;
      if (p < pRes + 10) {
	/* The working precision is less than the precision of the
	   result. This makes sense only if there is still hope that the
	   final result may fall below the cutoff.

	   There is still hope if the cutoff interval and the
	   current evaluation interval have a point in common.

	   When we loose faith that we still could do the cutoff with
	   low precision, we can directly jump to the current working
	   precision + precision of the result. The fact that we are
	   here means that we could not achieve enough precision for the
	   cutoff, hence that the expression cancelled on all bits (the
	   current precision).

	*/
	if (testCutOff) {
	  if (!sollya_mpfi_have_common_real_point(yI, cutoffI)) {
	    p = (p + 10 > pRes + p + 10? p + 10 : pRes + p + 10);
	    precisionIncreased = 1;
	  } else {
	    if ((p << 1) > (pRes + p + 10)) {
	      p = (p + 10 > pRes + p + 10? p + 10 : pRes + p + 10);
	      precisionIncreased = 1;
	    }
	  }
	} else {
	  p = (p + 10 > pRes + p + 10? p + 10 : pRes + p + 10);
	  precisionIncreased = 1;
	}
      }
    }

    if (!precisionIncreased) {
      /* Here, perform some default precision incrementation.

	 The right way to do so might be a matter of religion.

	 The formula below seems to be fine in practice.
      */
      p <<= 1;
    }
  }
  sollya_mpfi_clear(xI);
  mpfr_clear(yILeft);
  mpfr_clear(yIRight);

  if (okay > 0) {
    /* Now compute the correct rounding wherever possible and the
       correct rounding of the midpoint of both possible answers
       otherwise

       To do so, compute the midpoint of resUp and resDown and round
       it to the target precision.

       There are 3 possible cases:

       a) The loop above has been left because both endpoints of the
       proof interval rounded to the same value, resUp = resDown.
       In this case, the midpoint of resUp and resDown is
       resUp = resDown, which is also the correct rounding.

       b) The loop above has been left because resUp was the
       next floating-point number above resDown. resDown
       and resUp are the two possible answers for a faithful
       rounding. Their midpoint is the point where the rounding
       changes. The proof interval contained that midpoint.
       Returning the correct rounding of the midpoint to the
       nearest floating-point value hence chooses the even of
       both possible faithful roundings.

       c) The loop above has been left because the proof interval
       was contained in [-cutoff,cutoff]. In this case, resUp and
       resDown have both been set to the midpoint of the proof
       interval, which is certainly less than cutoff in magnitude.
       Computing the midpoint of resUp and resDown and rounding
       to target precision does not change anything to that value.
    */
    mpfr_div_2ui(resDown, resDown, 1, GMP_RNDN); /* exact, power of 2 */
    mpfr_div_2ui(resUp, resUp, 1, GMP_RNDN); /* exact, power of 2 */
    mpfr_add(result,resDown, resUp, GMP_RNDN); /* Computes the midpoint and rounds it to target precision */

    /* Now check if we are in a case when we can tell if the result is
       an exact or inexact one. We do this check only if we had a
       correct or faithful result; we do nothing for the case when
       the proof interval was below the threshold.
    */
    if (okay == 1) {
      /* We have a correctly or faithfully rounded result in result and we
	 want to check if the proof interval we had, yI allows anything
	 else to be said than just "it's a faithful rounding".

	 We won't do the additional check if there are quantities
	 involved that are no real numbers.
      */
      p = sollya_mpfi_get_prec(yI);
      mpfr_init2(yILeftCheck, p);
      mpfr_init2(yIRightCheck, p);
      sollya_mpfi_get_left(yILeftCheck,yI);
      sollya_mpfi_get_right(yIRightCheck,yI);
      if (mpfr_number_p(result) &&
	  mpfr_number_p(yILeftCheck) &&
	  mpfr_number_p(yIRightCheck)) {
	/* Here, we have real numbers to work on */
	if (mpfr_equal_p(yILeftCheck, yIRightCheck)) {
	  /* Here, the proof interval is reduced to a point */
	  if (mpfr_equal_p(result, yILeftCheck)) {
	    /* The result is exact as it is equal to both bounds of
	       the proof interval.
	    */
	    okay = 4;
	  } else {
	    /* The result is inexact as it is unequal to both bounds
	       of the proof interval which are themselves equal one to
	       another.
	    */
	    okay = 5;
	  }
	} else {
	  /* Given the correctly or faithfully rounded result, check
	     if the proof interval is sufficiently small to prove that
	     the result is not in the proof interval, which means that
	     the result must be inexact.
	  */
	  if (!((mpfr_cmp(yILeftCheck, result) <= 0) && (mpfr_cmp(result, yIRightCheck) <= 0))) {
	    /* Here, the result is not in the proof interval and
	       therefore it is necessarily inexact.
	    */
	    okay = 5;
	  }
	}
      }
      mpfr_clear(yILeftCheck);
      mpfr_clear(yIRightCheck);
    }

    /* If we got a correctly rounded result, we can express this fact
       with the return code
    */
    if (correctlyRounded) {
      if (okay == 1) {
	okay = 6;
      } else {
	if (okay == 5) {
	  okay = 7;
	}
      }
    }
  } else {
    mpfr_set_nan(result);
    if( (!mpfr_number_p(resUp)) || (!mpfr_number_p(resDown))) okay=3;
  }

  sollya_mpfi_clear(yI);
  mpfr_clear(resUp);
  mpfr_clear(resDown);
  mpfr_clear(resUpFaith);
  mpfr_clear(resDownFaith);
  mpfr_clear(resDownTemp);
  sollya_mpfi_clear(cutoffI);
  sollya_mpfi_clear(dummyI);
  return okay;
}

int evaluateFaithfulWithCutOffFast(mpfr_t result, node *func, node *deriv, mpfr_t x, mpfr_t cutoff, mp_prec_t startprec) {
  return evaluateFaithfulWithCutOffFastInternalImplementation(result, func, deriv, x, cutoff, startprec, NULL);
}

int evaluateFaithfulAtConstantExpression(mpfr_t result, node *func, node *deriv, node *x, mpfr_t cutoff, mp_prec_t startprec) {
  mpfr_t dummy;
  if (x == NULL) {
    mpfr_set_nan(result);
    return 0;
  }
  return evaluateFaithfulWithCutOffFastInternalImplementation(result, func, deriv, dummy, cutoff, startprec, x);
}

int evaluateFaithfulWithCutOff(mpfr_t result, node *func, mpfr_t x, mpfr_t cutoff, mp_prec_t startprec) {
  node *deriv;
  int res;

  if ((2*startprec) < (mpfr_get_prec(result) + 10)) {
    printMessage(12,SOLLYA_MSG_DIFFERENTIATING_FOR_EVAL_AS_START_PREC_LOW,"Information: Differentiating while evaluating because start precision (%d bits) too low.\n",
		 (int)startprec);
    deriv = differentiate(func);
  }else deriv = NULL;
  res = evaluateFaithfulWithCutOffFast(result, func, deriv, x, cutoff, startprec);
  if (res==3) res=0;
  if (deriv != NULL) free_memory(deriv);
  return res;
}


node *makeDoubleConstant(double d) {
  node *tempNode;

  tempNode = allocateNode();
  tempNode->nodeType = CONSTANT;
  tempNode->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(tempNode->value),64);
  mpfr_set_d(*(tempNode->value),d,GMP_RNDN);

  return tempNode;
}

int evaluateSign(int *s, node *func);

node *convertConstantToFunctionInPiInner(node *tree) {
  node *res;
  int a;

  if (tree->nodeType == MEMREF) return convertConstantToFunctionInPiInner(getMemRefChild(tree));

  if (tree->nodeType == PI_CONST) {
    return makeVariable();
  }

  a = arity(tree);
  switch (a) {
  case 0:
    res = copyTree(tree);
    break;
  case 1:
    res = allocateNode();
    res->nodeType = tree->nodeType;
    if (tree->nodeType == UNARY_BASE_FUNC) {
      res->baseFun = tree->baseFun;
    }
    if (tree->nodeType == LIBRARYFUNCTION) {
      res->libFun = tree->libFun;
      res->libFunDeriv = tree->libFunDeriv;
    }
    if (tree->nodeType == PROCEDUREFUNCTION) {
      res->libFunDeriv = tree->libFunDeriv;
      res->child2 = copyThing(tree->child1);
    }
    res->child1 = convertConstantToFunctionInPiInner(tree->child1);
    break;
  case 2:
    res = allocateNode();
    res->nodeType = tree->nodeType;
    res->child1 = convertConstantToFunctionInPiInner(tree->child1);
    res->child2 = convertConstantToFunctionInPiInner(tree->child2);
    break;
  default:
    sollyaFprintf(stderr,"Error: convertConstantToFunctionInPiInner: unknown arity (%d).\n",a);
    exit(1);
  }

  return res;
}

node *convertConstantToFunctionInPi(node *tree) {
  if (!isConstant(tree)) return NULL;
  return convertConstantToFunctionInPiInner(tree);
}


int containsPi(node *tree) {
  int a;

  if (tree->nodeType == MEMREF) return containsPi(getMemRefChild(tree));

  if (tree->nodeType == PI_CONST) return 1;

  a = arity(tree);
  switch (a) {
  case 0:
    return 0;
    break;
  case 1:
    return containsPi(tree->child1);
    break;
  case 2:
    return (containsPi(tree->child1) || containsPi(tree->child2));
  default:
    sollyaFprintf(stderr,"Error: containsPi: unknown arity (%d).\n",a);
    exit(1);
  }
  return 0;
}


int compareConstant(int *cmp, node *func1, node *func2, node *difference, int doNotEval) {
  node *diff, *rawDiff, *rawDiff2;
  int res, okay;
  mpfr_t value, dummyX;
  int okayA, okayB, signA, signB;
  node *tempNode;
  node **coefficients;
  int degree, i;
  int allZero, allOkay;
  int freeDiff = 1;

  okay = 0;
  res = -17; /* Makes compiler happy */

  if (difference == NULL) {
    rawDiff = addMemRef(makeSub(copyTree(func1),copyTree(func2)));
    rawDiff2 = simplifyRationalErrorfree(rawDiff);
    diff = simplifyTreeErrorfree(rawDiff2);
    freeDiff = 1;
    free_memory(rawDiff);
    free_memory(rawDiff2);
  } else {
    diff = difference;
    freeDiff = 0;
  }
  mpfr_init2(value,12);
  mpfr_init2(dummyX,12);
  mpfr_set_ui(dummyX,1,GMP_RNDN);
  if ((!doNotEval) &&
      evaluateFaithful(value, diff, dummyX, defaultprecision) &&
      mpfr_number_p(value)) {
    res = mpfr_sgn(value);
    okay = 1;
  } else {
    if ((accessThruMemRef(func1)->nodeType == DIV) &&
	evaluateSign(&signA,accessThruMemRef(func1)->child2) &&
	(signA != 0)) {
      tempNode = makeMul(copyTree(accessThruMemRef(func1)->child2),copyTree(func2));
      okayB = compareConstant(&signB, accessThruMemRef(func1)->child1, tempNode, NULL, 1);
      if (okayB) {
	okay = 1;
	res = signB;
      }
      free_memory(tempNode);
    }
    if (!okay) {
      if ((accessThruMemRef(func2)->nodeType == DIV) &&
	  evaluateSign(&signA,accessThruMemRef(func2)->child2) &&
	  (signA != 0)) {
	tempNode = makeMul(copyTree(accessThruMemRef(func2)->child2),copyTree(func1));
	okayB = compareConstant(&signB, tempNode, accessThruMemRef(func2)->child1, NULL, 1);
	if (okayB) {
	  okay = 1;
	  res = signB;
	}
	free_memory(tempNode);
      }
      if (!okay) {
	if ( (accessThruMemRef(func1)->nodeType == UNARY_BASE_FUNC) &&
             (accessThruMemRef(func2)->nodeType == UNARY_BASE_FUNC) &&
             (accessThruMemRef(func1)->baseFun == accessThruMemRef(func2)->baseFun) &&
             (accessThruMemRef(func1)->baseFun->isDefinedEverywhere) &&
             (accessThruMemRef(func1)->baseFun->monotonicity == INCREASING) ) {
	  okayA = compareConstant(&signA, accessThruMemRef(func1)->child1, accessThruMemRef(func2)->child1, NULL, 0);
	  if (okayA) {
	    okay = 1;
	    res = signA;
	  }
	}
	if (!okay) {
	  if ( ((accessThruMemRef(func1)->nodeType == NEG) && (accessThruMemRef(func2)->nodeType == NEG)) ||
               ( (accessThruMemRef(func1)->nodeType == UNARY_BASE_FUNC) &&
                 (accessThruMemRef(func2)->nodeType == UNARY_BASE_FUNC) &&
                 (accessThruMemRef(func1)->baseFun == accessThruMemRef(func2)->baseFun) &&
                 (accessThruMemRef(func1)->baseFun->isDefinedEverywhere) &&
                 (accessThruMemRef(func1)->baseFun->monotonicity == DECREASING) )
               ) {
	    okayA = compareConstant(&signA, accessThruMemRef(func1)->child1, accessThruMemRef(func2)->child1, NULL, 0);
	    if (okayA) {
	      okay = 1;
	      res = -signA;
	    }
	  }
	  if (!okay) {
	    if (containsPi(diff)) {
	      if ((tempNode = convertConstantToFunctionInPi(diff)) != NULL) {
		if (isPolynomial(tempNode)) {
		  /* Here we have tempNode(pi) = diff and tempNode a polynomial */
		  getCoefficients(&degree, &coefficients, tempNode);
		  if (degree >= 0) {
		    allZero = 1; allOkay = 1;
		    for (i=0;i<=degree;i++) {
		      if (coefficients[i] != NULL) {
			if (evaluateSign(&signA, coefficients[i])) {
			  if (signA != 0) {
			    allZero = 0;
			    break;
			  }
			} else {
			  allOkay = 0;
			  break;
			}
		      }
		    }
		    if (allOkay && allZero) {
		      okay = 1;
		      res = 0;
		    }
		    for (i=0;i<=degree;i++)
		      if (coefficients[i] != NULL) free_memory(coefficients[i]);
		    safeFree(coefficients);
		  }
		}
		free_memory(tempNode);
	      }
	    }
	    /* Put next case here */
	  }
	}
      }
    }
  }

  mpfr_clear(dummyX);
  mpfr_clear(value);
  if (freeDiff) free_memory(diff);

  if (okay) *cmp = res;
  return okay;
}

int evaluateSignFast(int *s, node *constFunc) {
  int okay;
  sollya_mpfi_t y;
  if (!isConstant(constFunc)) return 0;
  if (accessThruMemRef(constFunc)->nodeType == CONSTANT) {
    if (!mpfr_number_p(*(accessThruMemRef(constFunc)->value))) return 0;
    *s = mpfr_sgn(*(accessThruMemRef(constFunc)->value));
    return 1;
  }

  okay = 0;
  sollya_mpfi_init2(y, 12);
  evaluateConstantExpressionToInterval(y, constFunc);
  if (!sollya_mpfi_has_nan(y)) {
    if (sollya_mpfi_is_zero(y)) {
      okay = 1;
      *s = 0;
    } else {
      if (!sollya_mpfi_has_zero(y)) {
	okay = 1;
	if (sollya_mpfi_is_nonneg(y)) {
	  *s = 1;
	} else {
	  *s = -1;
	}
      }
    }
  }
  sollya_mpfi_clear(y);

  return okay;
}

int evaluateSign(int *s, node *rawFunc) {
  int sign, okay, okayA, okayB;
  mpfr_t value, dummyX;
  sollya_mpfi_t valueI;
  int signA, signB;
  node *tempNode;
  node *func, *rawFunc2;

  okay = 0;
  if (!isConstant(rawFunc)) return 0;

  if ((accessThruMemRef(rawFunc)->nodeType == CONSTANT) &&
      (!mpfr_number_p(*(accessThruMemRef(rawFunc)->value)))) return 0;

  mpfr_init2(value,12);
  mpfr_init2(dummyX,12);
  mpfr_set_ui(dummyX,1,GMP_RNDN);
  if (evaluateFaithful(value, rawFunc, dummyX, defaultprecision) &&
      mpfr_number_p(value)) {
    sign = mpfr_sgn(value);
    okay = 1;
  } else {
    rawFunc2 = simplifyRationalErrorfree(rawFunc);
    func = simplifyTreeErrorfree(rawFunc2);
    free_memory(rawFunc2);
    if (evaluateFaithful(value, func, dummyX, defaultprecision) &&
	mpfr_number_p(value)) {
      sign = mpfr_sgn(value);
      okay = 1;
    } else {
      switch (accessThruMemRef(func)->nodeType) {
      case CONSTANT:
	sign = mpfr_sgn(*(accessThruMemRef(func)->value));
	okay = 1;
	break;
      case ADD:
	tempNode = makeNeg(copyTree(accessThruMemRef(func)->child2));
	okay = compareConstant(&sign, accessThruMemRef(func)->child1, tempNode, NULL, 0);
	free_memory(tempNode);
	break;
      case SUB:
	okay = compareConstant(&sign, accessThruMemRef(func)->child1, accessThruMemRef(func)->child2, func, 0);
	break;
      case MUL:
	signA = 0;
	signB = 0;
	okay = (evaluateSign(&signA, accessThruMemRef(func)->child1) && evaluateSign(&signB, accessThruMemRef(func)->child2));
	sign = signA * signB;
	break;
      case DIV:
	okayA = (evaluateSign(&signA, accessThruMemRef(func)->child1) && evaluateSign(&signB, accessThruMemRef(func)->child2));
	if (okayA && (signB != 0)) {
	  okay = 1;
	  sign = signA * signB;
	}
	break;
      case NEG:
	okay = evaluateSign(&signA, accessThruMemRef(func)->child1);
	sign = -1 * signA;
	break;
      case UNARY_BASE_FUNC:
        okay = accessThruMemRef(func)->baseFun->evalsign(&sign, accessThruMemRef(func)->child1);
        break;
      case POW:
	okayA = evaluateSign(&signA, accessThruMemRef(func)->child1);
	okayB = evaluateSign(&signB, accessThruMemRef(func)->child2);
	if (okayA && okayB) {
	  if (signB == 0) {
	    okay = 1;
	    sign = 1;
	  } else {
	    if (signA == 0) {
	      if (signB > 0) {
		okay = 1;
		sign = 0;
	      }
	    } else {
	      if ((signA > 0) &&
		  (signB >= 0)) {
		okay = 1;
		sign = 1;
	      }
	    }
	  }
	}
	break;
      case LIBRARYFUNCTION:
	break;
      case PROCEDUREFUNCTION:
	break;
      case PI_CONST:
	okay = 1;
	sign = 1;
	break;
      case LIBRARYCONSTANT:
        /* By definition, a library constant is known with a relative error
           smaller that ~ 2^(-prec). So we can decide the sign, based on low
           approximation of the constant. */
        sollya_mpfi_init2(valueI, 12);
        libraryConstantToInterval(valueI, accessThruMemRef(func));
        if (sollya_mpfi_is_zero(valueI)) {
          okay = 1;
          sign = 0;
        }
        else {
          if (sollya_mpfi_has_zero(valueI)) {
            okay = 0;
            sign = 0;
          }
          else {
            okay = 1;
            sign = (sollya_mpfi_is_nonneg(valueI))?1:(-1);
          }
        }
        break;
      default:
	sollyaFprintf(stderr,"Error: evaluateSign: unknown identifier (%d) in the tree\n",func->nodeType);
	exit(1);
      }
    }
    free_memory(func);
  }
  mpfr_clear(value);
  mpfr_clear(dummyX);

  if (okay) *s = sign;
  return okay;
}
