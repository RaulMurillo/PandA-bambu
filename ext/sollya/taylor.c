/*

  Copyright 2006-2018 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France

  and by

  Sorbonne Université
  CNRS, Laboratoire d'Informatique de Paris 6, LIP6
  F - 75005 Paris
  France.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "taylor.h"
#include "expression.h"
#include "general.h"
#include "execute.h"



node *taylor(node* tree, int degree, node* point, mp_prec_t prec) {
  node *copy, *temp, *temp2, *fderiv, *fderivsubst, *denominator, *numerator, *expon, *variable, *term, *pointTemp;
  mpfr_t *value;
  mpz_t denominatorGMP;
  int i;

  if (!isConstant(point)) {
    printMessage(1,SOLLYA_MSG_DEVELOPMENT_POINT_NOT_CONSTANT,"Warning: the expression given for the development point is not constant.\n");
    printMessage(1,SOLLYA_MSG_CONTINUATION,"Will evaluate the expression in %s = 0 before using it as development point.\n",((variablename == NULL) ? "_x_" : variablename));
    temp = allocateNode();
    temp->nodeType = CONSTANT;
    value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,prec);
    mpfr_set_d(*value,0.0,GMP_RNDN);
    temp->value = value;
    temp2 = substitute(point,temp);
    pointTemp = simplifyTreeErrorfree(temp2);
    free_memory(temp);
    free_memory(temp2);
  } else {
    pointTemp = copyTree(point);
  }

  value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*value,prec);
  mpfr_set_d(*value,0.0,GMP_RNDN);
  copy = allocateNode();
  copy->nodeType = CONSTANT;
  copy->value = value;

  mpz_init(denominatorGMP);
  fderiv = copyTree(tree);
  for (i=0;i<=degree;i++) {
    temp = addMemRef(substitute(fderiv,pointTemp));
    fderivsubst = addMemRef(simplifyTreeErrorfree(temp));
    tryRepresentAsPolynomial(fderivsubst);
    free_memory(temp);
    mpz_fac_ui(denominatorGMP,(unsigned int) i);
    value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((prec > 74) ? (prec) : (74)));
    if(mpfr_set_z(*value,denominatorGMP,GMP_RNDN) != 0) {
      if (!noRoundingWarnings) {
	printMessage(1,SOLLYA_MSG_ROUNDING_ON_COMPUTATION_OF_TAYLOR_COEFFICIENT,"Warning: rounding occurred on computing a taylor constant factor.\n");
	printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the working precision.\n");
      }
    }
    denominator = allocateNode();
    denominator->nodeType = CONSTANT;
    denominator->value = value;
    value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,prec);
    mpfr_set_d(*value,1.0,GMP_RNDN);
    numerator = allocateNode();
    numerator->nodeType = CONSTANT;
    numerator->value = value;
    temp = allocateNode();
    temp->nodeType = DIV;
    temp->child1 = addMemRef(numerator);
    temp->child2 = addMemRef(denominator);
    temp2 = allocateNode();
    temp2->nodeType = MUL;
    temp2->child1 = addMemRef(temp);
    temp2->child2 = addMemRef(fderivsubst);
    variable = addMemRef(makeVariable());
    value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((prec > 74) ? (prec) : (74)));
    if(mpfr_set_si(*value,i,GMP_RNDN) != 0) {
      if (!noRoundingWarnings) {
	printMessage(1,SOLLYA_MSG_ROUNDING_ON_COMPUTATION_OF_TAYLOR_POWER,"Warning: rounding occurred on computing a taylor exponent.\n");
	printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the working precision.\n");
      }
    }
    expon = allocateNode();
    expon->nodeType = CONSTANT;
    expon->value = value;
    temp = allocateNode();
    temp->nodeType = POW;
    temp->child1 = variable;
    temp->child2 = addMemRef(expon);
    term = allocateNode();
    term->nodeType = MUL;
    term->child1 = addMemRef(temp2);
    term->child2 = addMemRef(temp);
    temp = allocateNode();
    temp->nodeType = ADD;
    temp->child1 = addMemRef(copy);
    temp->child2 = addMemRef(term);
    copy = addMemRef(temp);
    if (i < degree) {
      temp = differentiate(fderiv);
      free_memory(fderiv);
      fderiv = temp;
    }
  }
  mpz_clear(denominatorGMP);

  free_memory(fderiv);

  copy = addMemRef(copy);
  tryRepresentAsPolynomial(copy);
  temp = addMemRef(horner(copy));
  free_memory(copy);
  free_memory(pointTemp);
  return temp;
}
