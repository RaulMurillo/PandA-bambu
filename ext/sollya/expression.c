/*

  Copyright 2006-2018 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2),

  Centre de recherche INRIA Sophia Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France,

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

  CNRS, LIP6, UPMC
  Sorbonne Universités, UPMC Univ Paris 06,
  CNRS, LIP6 UMR 7606, 4 place Jussieu 75005 Paris.

  Contributors Ch. Lauter, S. Chevillard, M. Mezzarobba

  christoph.lauter@ens-lyon.org
  sylvain.chevillard@ens-lyon.org
  marc@mezzarobba.net

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

#include <limits.h>
#include <inttypes.h>
#include <mpfr.h>
#include "mpfi-compat.h"
#include <gmp.h>
#include "expression.h"
#include <stdio.h> /* fprinft, fopen, fclose, */
#include <stdlib.h> /* exit, free, mktemp */
#include <string.h>
#include <errno.h>
#include "general.h"
#include "double.h"
#include "chain.h"
#include "execute.h"
#include "infnorm.h"
#include "signalhandling.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define MAXDIFFSIMPLSIZE 100
#define MAXDIFFSIMPLDEGREE 25
#define MAXDIFFPOLYSPECIALDEGREE 300

#if (!(defined(HAVE_MP_BITCNT_T) && (HAVE_MP_BITCNT_T)))
typedef unsigned long int mp_bitcnt_t;
#endif

static inline void copyTreeAnnotations(node *new, node *old) {
  if (new == NULL) return;
  if (old == NULL) return;
  if (new->nodeType != MEMREF) return;
  if (old->nodeType != MEMREF) return;
  if (new == old) return;
  new->cache->isCorrectlyTyped = old->cache->isCorrectlyTyped;
  if ((old->cache->derivCache != NULL) && (new->cache->derivCache == NULL)) {
    new->cache->derivCache = copyThing(old->cache->derivCache);
  }
  if ((old->cache->derivUnsimplCache != NULL) && (new->cache->derivUnsimplCache == NULL)) {
    new->cache->derivUnsimplCache = copyThing(old->cache->derivUnsimplCache);
  }
  if ((old->cache->simplifyCache != NULL) && (new->cache->simplifyCache == NULL)) {
    new->cache->simplifyCache = copyThing(old->cache->simplifyCache);
  }
  if ((old->cache->simplifyCacheRationalMode >= 0) && (new->cache->simplifyCacheRationalMode < 0)) {
    new->cache->simplifyCacheRationalMode = old->cache->simplifyCacheRationalMode;
  }
  if ((old->cache->simplifyCacheDoesNotSimplify >= 0) && (new->cache->simplifyCacheDoesNotSimplify < 0)) {
    new->cache->simplifyCacheDoesNotSimplify = old->cache->simplifyCacheDoesNotSimplify;
  }
  addEvaluationHookFromCopy(&(new->cache->evaluationHook), old->cache->evaluationHook);
}

static inline void copyTreeAnnotationsNoSimplifications(node *new, node *old) {
  if (new == NULL) return;
  if (old == NULL) return;
  if (new->nodeType != MEMREF) return;
  if (old->nodeType != MEMREF) return;
  if (new == old) return;
  new->cache->isCorrectlyTyped = old->cache->isCorrectlyTyped;
  if ((old->cache->derivCache != NULL) && (new->cache->derivCache == NULL)) {
    new->cache->derivCache = copyThing(old->cache->derivCache);
  }
  if ((old->cache->derivUnsimplCache != NULL) && (new->cache->derivUnsimplCache == NULL)) {
    new->cache->derivUnsimplCache = copyThing(old->cache->derivUnsimplCache);
  }
  addEvaluationHookFromCopy(&(new->cache->evaluationHook), old->cache->evaluationHook);
}

void tryCopyTreeAnnotations(node *newTree, node *oldTree) {
  copyTreeAnnotations(newTree, oldTree);
}

int isSyntacticallyEqualCheap(node *tree1, node *tree2);

void simplifyMpfrPrec(mpfr_t rop, mpfr_t op) {
  mp_prec_t prec;

  prec = mpfr_min_prec(op);
  if (prec < 12) prec = 12;
  mpfr_set_prec(rop, prec);
  mpfr_set(rop, op, GMP_RNDN); /* exact */
}



void sollya_mpfr_from_mpfi(mpfr_t rop, mpfr_t op, int n, int (*mpfifun)(sollya_mpfi_t, sollya_mpfi_t, int)) {
  sollya_mpfi_t opI, ropI, ropItemp;

  sollya_mpfi_init2(opI,mpfr_get_prec(op));
  sollya_mpfi_init2(ropItemp,mpfr_get_prec(rop)+2);
  sollya_mpfi_set_fr(opI,op);

  enterExternalCode();
  mpfifun(ropItemp,opI,n);
  leaveExternalCode();

  sollya_init_and_convert_interval(ropI, ropItemp);

  sollya_mpfi_mid(rop,ropI);

  sollya_mpfi_clear(opI);
  sollya_mpfi_clear(ropI);
  sollya_mpfi_clear(ropItemp);
}

void sollya_mpfr_from_mpfi_data(mpfr_t rop, mpfr_t op, int n, int (*mpfifun)(sollya_mpfi_t, sollya_mpfi_t, int, void *), void *data) {
  sollya_mpfi_t opI, ropI, ropItemp;

  sollya_mpfi_init2(opI,mpfr_get_prec(op));
  sollya_mpfi_init2(ropItemp,mpfr_get_prec(rop)+2);
  sollya_mpfi_set_fr(opI,op);

  enterExternalCode();
  mpfifun(ropItemp,opI,n,data);
  leaveExternalCode();

  sollya_init_and_convert_interval(ropI, ropItemp);

  sollya_mpfi_mid(rop,ropI);

  sollya_mpfi_clear(opI);
  sollya_mpfi_clear(ropI);
  sollya_mpfi_clear(ropItemp);
}


void free_memory(node *tree) {
  if (tree == NULL) return;
  freeThing(tree);
}


void fprintHeadFunction(FILE *fd,node *tree, char *x, char *y) {
  int i;

  if (tree == NULL) return;
  switch (tree->nodeType) {
  case MEMREF:
    fprintHeadFunction(fd,getMemRefChild(tree), x, y);
    break;
  case VARIABLE:
    if (x != NULL) sollyaFprintf(fd,"%s",x); else sollyaFprintf(fd,"x");
    break;
  case CONSTANT:
    fprintValue(fd,*(tree->value));
    break;
  case ADD:
    sollyaFprintf(fd,"%s + %s",x,y);
    break;
  case SUB:
    sollyaFprintf(fd,"%s - %s",x,y);
    break;
  case MUL:
    sollyaFprintf(fd,"%s * %s",x,y);
    break;
  case DIV:
    sollyaFprintf(fd,"%s / %s",x,y);
    break;
  case NEG:
    sollyaFprintf(fd,"-%s",x);
    break;
  case UNARY_BASE_FUNC:
    sollyaFprintf(fd,"%s",tree->baseFun->functionName);
    sollyaFprintf(fd,"(%s)",x);
    break;
  case POW:
    sollyaFprintf(fd,"%s^%s",x,y);
    break;
  case LIBRARYFUNCTION:
    {
      sollyaFprintf(fd,"(");
      for (i=1;i<=tree->libFunDeriv;i++) {
	sollyaFprintf(fd,"diff(");
      }
      sollyaFprintf(fd,"%s",tree->libFun->functionName);
      for (i=1;i<=tree->libFunDeriv;i++) {
	sollyaFprintf(fd,")");
      }
      sollyaFprintf(fd,")(%s)",x);
    }
    break;
  case PROCEDUREFUNCTION:
    {
      sollyaFprintf(fd,"(");
      for (i=1;i<=tree->libFunDeriv;i++) {
	sollyaFprintf(fd,"diff(");
      }
      sollyaFprintf(fd,"function(");
      fPrintThing(fd,tree->child2);
      sollyaFprintf(fd,")");
      for (i=1;i<=tree->libFunDeriv;i++) {
	sollyaFprintf(fd,")");
      }
      sollyaFprintf(fd,")(%s)",x);
    }
    break;
  case PI_CONST:
    sollyaFprintf(fd,"pi");
    break;
  case LIBRARYCONSTANT:
    sollyaFprintf(fd,"%s",tree->libFun->functionName);
    break;
  default:
    sollyaFprintf(stderr,"fprintHeadFunction: unknown identifier (%d) in the tree\n",tree->nodeType);
    exit(1);
  }
  return;
}

int precedence(node *tree) {
  switch (tree->nodeType) {
  case MEMREF:
    return precedence(getMemRefChild(tree));
    break;
  case CONSTANT:
  case VARIABLE:
  case PI_CONST:
    return 1;
    break;
  case ADD:
  case SUB:
    return 2;
    break;
  case MUL:
    return 3;
    break;
  case DIV:
    return 4;
    break;
  case NEG:
    return 5;
    break;
  case POW:
    return 6;
    break;
  default:
    return 7;
  }
  return 0;
}


int isInfix(node *tree) {
  char *str;
  int res;
  switch(tree->nodeType) {
  case MEMREF:
    return isInfix(getMemRefChild(tree));
    break;
  case CONSTANT:
    if (mpfr_sgn(*(tree->value)) < 0) return 1;
    if ((dyadic == 2) || (dyadic == 3)) {
      str = sprintValue(tree->value);
      res = (sollya_strchr(str,'*') != NULL);
      safeFree(str);
      return res;
    }
    break;
  case PI_CONST:
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
  case NEG:
    return 1;
    break;
  default: return 0;
  }
  return 0;
}

void removeTrailingZeros(char *outbuf, char *inbuf) {
  char *temp, *temp2, *temp3;

  temp = inbuf; temp2 = outbuf; temp3 = outbuf;
  while ((temp != NULL) && (*temp != '\0')) {
    *temp2 = *temp;
    if (*temp2 != '0') {
      temp3 = temp2;
    }
    temp2++;
    temp++;
  }
  temp3++;
  *temp3 = '\0';
}

void printHexadecimalValue(mpfr_t x);

void printValue(mpfr_t *value) {
  char *str;

  str = sprintValue(value);
  sollyaPrintf("%s",str);
  safeFree(str);
}

char *sprintMidpointMode(mpfr_t a, mpfr_t b) {
  mp_exp_t e1, e2;
  char *str1, *str2, *str, *str3;
  mpfr_t aP, bP;
  int sign, len1, len2, len, i;
  mp_prec_t prec, p;
  char *s1, *s2;
  int possibleLength;

  if (mpfr_sgn(a) != mpfr_sgn(b)) return NULL;

  if (mpfr_zero_p(a)) {
    str = safeCalloc(7,sizeof(char));
    sprintf(str,"[0]");
    return str;
  }

  prec = mpfr_get_prec(a);
  p = mpfr_get_prec(b);

  if (p > prec) prec = p;

  mpfr_init2(aP,prec);
  mpfr_init2(bP,prec);

  sign = mpfr_sgn(a);
  if (sign > 0) {
    mpfr_set(aP,a,GMP_RNDN);
    mpfr_set(bP,b,GMP_RNDN);
  } else {
    mpfr_neg(aP,b,GMP_RNDN);
    mpfr_neg(bP,a,GMP_RNDN);
  }

  str1 = mpfr_get_str(NULL,&e1,10,0,aP,GMP_RNDD);
  str2 = mpfr_get_str(NULL,&e2,10,0,bP,GMP_RNDU);

  str3 = safeCalloc(strlen(str1) + 1, sizeof(char));
  removeTrailingZeros(str3,str1);
  safeFree(str1);
  str1 = str3;

  str3 = safeCalloc(strlen(str2) + 1, sizeof(char));
  removeTrailingZeros(str3,str2);
  safeFree(str2);
  str2 = str3;

  if (e1 == e2) {
    if (strcmp(str1,str2) == 0) {
      mpfr_set(aP,a,GMP_RNDN);
      str3 = sprintValue(&aP);
      str = (char *) safeCalloc(strlen(str3) + 3, sizeof(char));
      sprintf(str,"[%s]",str3);
      safeFree(str3);
    } else {

      if (str1[0] == str2[0]) {
	len1 = strlen(str1);
	len2 = strlen(str2);
	len = len1;
	if (len2 < len) len = len2;
	i = 0;
	while ((i < len) && (str1[i] == str2[i])) {
	  i++;
	}
	possibleLength = i;
	s1 = mpfr_get_str(NULL,&e1,10,possibleLength+1,aP,GMP_RNDD);
	s2 = mpfr_get_str(NULL,&e2,10,possibleLength+1,bP,GMP_RNDU);

	if (e1 == e2) {
	  if (s1[0] == s2[0]) {
	    len1 = strlen(s1);
	    len2 = strlen(s2);
	    len = len1;
	    if (len2 < len) len = len2;
	    str = (char *) safeCalloc(len+6,sizeof(char));
	    i = 0;
	    while ((i < len) && (s1[i] == s2[i])) {
	      str[i] = s1[i];
	      i++;
	    }
	    str[i] = '~';
	    if (sign > 0)
	      str[i+1] = s1[i];
	    else
	      str[i+1] = s2[i];
	    str[i+2] = '/';
	    if (sign > 0)
	      str[i+3] = s2[i];
	    else
	      str[i+3] = s1[i];
	    str[i+4] = '~';
	    str3 = (char *) safeCalloc(strlen(str)+1,sizeof(char));
	    removeTrailingZeros(str3,str);
	    safeFree(str);
	    str = str3;
	    str3 = (char *) safeCalloc(strlen(str)+69,sizeof(char));
	    if (sign < 0) {
	      if (e1 == 0) {
		sprintf(str3,"-0.%s",str);
	      } else {
		sprintf(str3,"-0.%se%d",str,(int)e1);
	      }
	    } else {
	      if (e1 == 0) {
		sprintf(str3,"0.%s",str);
	      } else {
		sprintf(str3,"0.%se%d",str,(int)e1);
	      }
	    }
	    safeFree(str);
	    str = str3;
	    str3 = (char *) safeCalloc(strlen(str)+1,sizeof(char));
	    sprintf(str3,"%s",str);
	    safeFree(str);
	    str = str3;
	  } else {
	    str = NULL;
	  }
	} else {
	  str = NULL;
	}
	safeFree(s1);
	safeFree(s2);
      } else {
	str = NULL;
      }
    }
  } else {
    str = NULL;
  }

  mpfr_free_str(str1);
  mpfr_free_str(str2);

  mpfr_clear(aP);
  mpfr_clear(bP);
  return str;
}

char *sPrintBinary(mpfr_t x) {
  mpfr_t xx;
  int negative;
  mp_prec_t prec;
  mp_exp_t expo;
  char *raw, *formatted, *temp1, *temp2, *str3;
  char *temp3=NULL;
  char *resultStr;
  int len;

  prec = mpfr_get_prec(x);
  mpfr_init2(xx,prec);
  mpfr_abs(xx,x,GMP_RNDN);
  negative = 0;
  if (mpfr_sgn(x) < 0) negative = 1;
  raw = mpfr_get_str(NULL,&expo,2,0,xx,GMP_RNDN);
  if (raw == NULL) {
    sollyaFprintf(stderr,"Error: unable to get a string for the given number.\n");
    exit(1);
  } else {
    formatted = safeCalloc(strlen(raw) + 3, sizeof(char));
    temp1 = raw; temp2 = formatted;
    if (negative) {
      *temp2 = '-';
      temp2++;
    }
    *temp2 = *temp1; temp2++; temp1++;
    if (*temp1 != '\0') {
      *temp2 = '.';
      temp2++;
    }
    while (*temp1 != '\0') {
      *temp2 = *temp1;
      temp2++; temp1++;
    }
    str3 = (char *) safeCalloc(strlen(formatted)+2,sizeof(char));
    removeTrailingZeros(str3,formatted);
    len = strlen(str3) - 1;
    if (str3[len] == '.') {
      str3[len] = '\0';
    }
    if (!mpfr_zero_p(x)) {
      if (mpfr_number_p(x)) {
	temp3 = (char *) safeCalloc(strlen(str3)+74,sizeof(char));
	if ((((int) expo)-1) != 0)
	  sprintf(temp3,"%s_2 * 2^(%d)",str3,((int)expo)-1);
	else
	  sprintf(temp3,"%s_2",str3);
      } else {
	temp3 = (char *) safeCalloc(strlen(raw) + 2,sizeof(char));
	if (negative)
	  sprintf(temp3,"-%s",raw);
	else
	  sprintf(temp3,"%s",raw);
      }
    }
    else {
      temp3 = (char *) safeCalloc(2,sizeof(char));
      sprintf(temp3,"0");
    }
    safeFree(formatted);
    safeFree(str3);
  }
  mpfr_free_str(raw);
  mpfr_clear(xx);
  resultStr = (char *) safeCalloc(strlen(temp3) + 1,sizeof(char));
  sprintf(resultStr,"%s",temp3);
  safeFree(temp3);
  return resultStr;
}

static inline int __sprintfValue_sprintf(char *str, const char *format, ...) {
  va_list varlist;
  int res;
  va_start(varlist, format);
  deferSignalHandling();
  res = vsprintf(str, format, varlist);
  resumeSignalHandling();
  va_end(varlist);
  return res;
}

char *sPrintHexadecimal(mpfr_t x) {
  char *res, *curr;
  mpfr_exp_t E;
  mpz_t m, n, z;
  int s;
  mp_bitcnt_t zerosRight;
  size_t q, w, a, i;
  long long int tmp;
  mpfr_exp_t Et;

  /* Display zero */
  if (mpfr_zero_p(x)) {
    res = safeCalloc(strlen("0") + 1, sizeof(char));
    strcpy(res, "0");
    return res;
  }

  /* Display NaN */
  if (mpfr_nan_p(x)) {
    res = safeCalloc(strlen("NaN") + 1, sizeof(char));
    strcpy(res, "NaN");
    return res;
  }

  /* Display +/- Infinity */
  if (mpfr_inf_p(x)) {
    if (mpfr_sgn(x) < 0) {
      /* Display -infty */
      res = safeCalloc(strlen("-infty") + 1, sizeof(char));
      strcpy(res, "-infty");
    } else {
      /* Display infty */
      res = safeCalloc(strlen("infty") + 1, sizeof(char));
      strcpy(res, "infty");
    }
    return res;
  }

  /* Here, the input is neither zero, nor NaN, nor an infinity

     We can therefore split it into a sign, an exponent and an odd
     integer mantissa.

  */
  mpz_init(m);
  E = mpfr_get_z_2exp(m, x);
  s = 0;
  if (mpz_sgn(m) < 0) {
    mpz_neg(m, m);
    s = 1;
  }
  zerosRight = mpz_scan1(m, 0);
  mpz_fdiv_q_2exp(m, m, zerosRight);
  E += zerosRight;

  /* Here, we have x = (-1)^s * 2^E * m

     where m = 2 * k + 1, k in N

     We now compute q = floor(log2(m)) and
     then

     n = m - 2^q.

  */
  q = mpz_sizeinbase(m, 2);
  q--;

  mpz_init(n);
  mpz_set_si(n, 1);
  if (q != ((size_t) 0)) {
    mpz_mul_2exp(n, n, (mp_bitcnt_t) q);
  }
  mpz_sub(n, m, n);

  /* Now, we have to display

     x = (-1)^s * 2^(E + q) * (1 + n * 2^(-q))

     We start by checking if n is zero. In this
     case m = 2^q + n is an integer power of 2.

  */
  if (mpz_sgn(n) == 0) {
    /* Here, we have n = 0, so we have to display

       x = (-1)^s * 2^(E + q) * 1

    */
    E += q;
    res = safeCalloc(1 + 2 + 1 + 8 * sizeof(mpfr_exp_t) + 1, sizeof(char));
    if (s) {
      res[0] = '-';
      curr = res + 1;
    } else {
      curr = res;
    }
    curr[0] = '0';
    curr[1] = 'x';
    curr[2] = '1';
    curr[3] = 'p';
    curr += 4;
    if (sizeof(mpfr_exp_t) <= sizeof(int)) {
      __sprintfValue_sprintf(curr, "%d", (int) E);
    } else {
      if (sizeof(mpfr_exp_t) <= sizeof(long int)) {
	__sprintfValue_sprintf(curr, "%ld", (long int) E);
      } else {
	if (sizeof(mpfr_exp_t) <= sizeof(long long int)) {
	  __sprintfValue_sprintf(curr, "%lld", (long long int) E);
	} else {
	  tmp = (long long int) E;
	  Et = (mpfr_exp_t) tmp;
	  if (E == Et) {
	    __sprintfValue_sprintf(curr, "%lld", tmp);
	  } else {
	    sollyaFprintf(stderr,"Error: sPrintHexadecimal: mpfr_exp_t variable too large to be displayed.\n");
	    exit(1);
	  }
	}
      }
    }
  } else {
    /* Here, we have n = 0, so we have to display

       x = (-1)^s * 2^(E + q) * (1 + z * 2^(-w))

       where

       w = 4 * ceil(q/4)

       which is an integer divisible by 4

       and

       z = n * 2^(w-q)

       Hence

       w - q = 4 * ceil(q / 4) - q
       = 4 * (q / 4 + delta) - q
       = 4 * delta

       where delta is bounded by 0 <= delta < 1

       Since w - q is integer, we have

       0 <= w - q <= 3.

       Further we have

       1 <= n <= 2^q - 1

       and hence

       1 <= z <= 2^w - 1.

    */
    w = q >> 2;
    if ((q & ((size_t) 3)) != ((size_t) 0)) w++;
    w <<= 2;
    mpz_init(z);
    if ((w - q) != ((size_t) 0)) {
      mpz_mul_2exp(z, n, (mp_bitcnt_t) (w - q));
    } else {
      mpz_set(z, n);
    }
    E += q;
    res = safeCalloc(1 + 2 + 1 + 1 + ((size_t) (w >> 2)) + 1 + 1 + 8 * sizeof(mpfr_exp_t) + 1, sizeof(char));
    if (s) {
      res[0] = '-';
      curr = res + 1;
    } else {
      curr = res;
    }
    curr[0] = '0';
    curr[1] = 'x';
    curr[2] = '1';
    curr[3] = '.';
    curr += 4;
    a = mpz_sizeinbase(z, 16);
    for (i=((size_t) 0);i<((size_t) ((w >> 2) - a));i++) {
      curr[i] = '0';
    }
    curr += ((size_t) ((w >> 2) - a));
    mpz_get_str(curr, 16, z);
    curr += a;
    curr[0] = 'p';
    curr++;
    if (sizeof(mpfr_exp_t) <= sizeof(int)) {
      __sprintfValue_sprintf(curr, "%d", (int) E);
    } else {
      if (sizeof(mpfr_exp_t) <= sizeof(long int)) {
	__sprintfValue_sprintf(curr, "%ld", (long int) E);
      } else {
	if (sizeof(mpfr_exp_t) <= sizeof(long long int)) {
	  __sprintfValue_sprintf(curr, "%lld", (long long int) E);
	} else {
	  tmp = (long long int) E;
	  Et = (mpfr_exp_t) tmp;
	  if (E == Et) {
	    __sprintfValue_sprintf(curr, "%lld", tmp);
	  } else {
	    sollyaFprintf(stderr,"Error: sPrintHexadecimal: mpfr_exp_t variable too large to be displayed.\n");
	    exit(1);
	  }
	}
      }
    }
    mpz_clear(z);
  }

  /* Clear temporaries */
  mpz_clear(n);
  mpz_clear(m);

  /* Return the result */
  return res;
}

void printBinary(mpfr_t x) {
  char *str;

  str = sPrintBinary(x);
  sollyaPrintf("%s",str);
  safeFree(str);
}

void printHexadecimalValue(mpfr_t x) {
  char *str;

  str = sPrintHexadecimal(x);
  sollyaPrintf("%s",str);
  safeFree(str);
}

static inline void __sprintfValue_init_value(mpfr_t x, mpfr_t *valPtr) {
  int ternary;

  /* Check if input is NaN or Inf */
  if (!mpfr_number_p(*valPtr)) {
    mpfr_init2(x, mpfr_get_prec(*valPtr));
    mpfr_set(x, *valPtr, GMP_RNDN);
    return;
  }

  /* Here the input is a real number */
  mpfr_init2(x, tools_precision);
  ternary = mpfr_set(x, *valPtr, GMP_RNDN);
  if ((ternary == 0) && (mpfr_number_p(x))) {
    /* The number holds on tools_precision bits */
    return;
  }

  /* The number does not hold on tools_precision bits, so use its own
     precision
  */
  mpfr_set_prec(x, mpfr_get_prec(*valPtr));
  mpfr_set(x, *valPtr, GMP_RNDN); /* exact */
}

static inline char *__sprintfValue_print_special(mpfr_t x) {
  char *res;

  /* Display zero */
  if (mpfr_zero_p(x)) {
    res = safeCalloc(strlen("0") + 1, sizeof(char));
    strcpy(res, "0");
    return res;
  }

  /* Display NaN */
  if (mpfr_nan_p(x)) {
    res = safeCalloc(strlen("NaN") + 1, sizeof(char));
    strcpy(res, "NaN");
    return res;
  }

  /* Display +/- Infinity */
  if (mpfr_inf_p(x)) {
    if (mpfr_sgn(x) < 0) {
      /* Display -infty */
      res = safeCalloc(strlen("-infty") + 1, sizeof(char));
      strcpy(res, "-infty");
    } else {
      /* Display infty */
      res = safeCalloc(strlen("infty") + 1, sizeof(char));
      strcpy(res, "infty");
    }
    return res;
  }

  /* Returning NULL on purpose for unhandled cases */
  return NULL;
}

static inline int __sprintfValue_is_small_integer(mpfr_t x) {
  if (!mpfr_number_p(x)) return 0;
  if (mpfr_zero_p(x)) return 0;
  if (!mpfr_integer_p(x)) return 0;
  return (mpfr_get_exp(x) <= 128);
}

static inline char *__sprintfValue_print_integer(mpfr_t x) {
  mpz_t z;
  char *res;

  /* Return NULL on purpose for unhandled cases */
  if (!mpfr_number_p(x)) return NULL;
  if (!mpfr_integer_p(x)) return NULL;

  /* x is an integer, convert it to MPZ */
  mpz_init(z);
  mpfr_get_z(z, x, GMP_RNDN); /* exact */

  /* Allocate output memory */
  res = safeCalloc(mpz_sizeinbase(z, 10) + 2 + 1, sizeof(char));
  mpz_get_str(res, 10, z);

  /* Clear temporary */
  mpz_clear(z);

  /* Return the result */
  return res;
}

static inline char *__sprintfValue_print_dyadic(mpfr_t x) {
  mpfr_exp_t E;
  mpz_t m;
  int s;
  mp_bitcnt_t zerosRight;
  char *res;
  char *curr;
  long long int tmp;
  mpfr_exp_t Et;

  /* Return NULL on purpose for unhandled cases */
  if (!mpfr_number_p(x)) return NULL;

  /* Handle the case of zero */
  if (mpfr_zero_p(x)) {
    res = safeCalloc(strlen("0") + 1, sizeof(char));
    strcpy(res, "0");
    return res;
  }

  /* Here the input is a non-zero real number

     Decompose x into (-1)^s * 2^E * m

  */
  mpz_init(m);
  E = mpfr_get_z_2exp(m, x);
  s = 0;
  if (mpz_sgn(m) < 0) {
    mpz_neg(m,m);
    s = 1;
  }
  zerosRight = mpz_scan1(m, 0);
  mpz_fdiv_q_2exp(m, m, zerosRight);
  E += zerosRight;

  /* Here, we have x = (-1)^s * 2^E * m

     where m = 2 * k + 1, k in Z

  */
  if (E == ((mpfr_exp_t) 0)) {
    /* We have no exponent to display */
    res = safeCalloc(mpz_sizeinbase(m, 10) + 2 + 1 + 1, sizeof(char));
    if (s) {
      res[0] = '-';
      curr = res + 1;
    } else {
      curr = res;
    }
    mpz_get_str(curr, 10, m);
  } else {
    /* We have the exponent to display */
    res = safeCalloc(mpz_sizeinbase(m, 10) + 2 + 1 + 1 + sizeof(mpfr_exp_t) * 8 + 1 + 1, sizeof(char));
    if (s) {
      res[0] = '-';
      curr = res + 1;
    } else {
      curr = res;
    }
    mpz_get_str(curr, 10, m);
    curr = res + strlen(res);
    curr[0] = 'b';
    curr++;
    if (sizeof(mpfr_exp_t) <= sizeof(int)) {
      __sprintfValue_sprintf(curr, "%d", (int) E);
    } else {
      if (sizeof(mpfr_exp_t) <= sizeof(long int)) {
	__sprintfValue_sprintf(curr, "%ld", (long int) E);
      } else {
	if (sizeof(mpfr_exp_t) <= sizeof(long long int)) {
	  __sprintfValue_sprintf(curr, "%lld", (long long int) E);
	} else {
	  tmp = (long long int) E;
	  Et = (mpfr_exp_t) tmp;
	  if (E == Et) {
	    __sprintfValue_sprintf(curr, "%lld", tmp);
	  } else {
	    sollyaFprintf(stderr,"Error: __sprintfValue_print_dyadic: mpfr_exp_t variable too large to be displayed.\n");
	    exit(1);
	  }
	}
      }
    }
  }

  /* Free temporary variable */
  mpz_clear(m);

  /* Return the result */
  return res;
}

static inline char *__sprintfValue_print_powers(mpfr_t x) {
  mpfr_exp_t E;
  mpz_t m;
  int s;
  mp_bitcnt_t zerosRight;
  char *res;
  char *curr;
  long long int tmp;
  mpfr_exp_t Et;

  /* Return NULL on purpose for unhandled cases */
  if (!mpfr_number_p(x)) return NULL;

  /* Handle the case of zero */
  if (mpfr_zero_p(x)) {
    res = safeCalloc(strlen("0") + 1, sizeof(char));
    strcpy(res, "0");
    return res;
  }

  /* Here the input is a non-zero real number

     Decompose x into (-1)^s * 2^E * m

  */
  mpz_init(m);
  E = mpfr_get_z_2exp(m, x);
  s = 0;
  if (mpz_sgn(m) < 0) {
    mpz_neg(m,m);
    s = 1;
  }
  zerosRight = mpz_scan1(m, 0);
  mpz_fdiv_q_2exp(m, m, zerosRight);
  E += zerosRight;

  /* Here, we have x = (-1)^s * 2^E * m

     where m = 2 * k + 1, k in Z

  */
  if (E == ((mpfr_exp_t) 0)) {
    /* We have no exponent to display */
    res = safeCalloc(mpz_sizeinbase(m, 10) + 2 + 1 + 1, sizeof(char));
    if (s) {
      res[0] = '-';
      curr = res + 1;
    } else {
      curr = res;
    }
    mpz_get_str(curr, 10, m);
  } else {
    /* We have the exponent to display */
    res = safeCalloc(mpz_sizeinbase(m, 10) + 2 + 1 + 1 + 7 + sizeof(mpfr_exp_t) * 8 + 1 + 1, sizeof(char));
    if (s) {
      res[0] = '-';
      curr = res + 1;
    } else {
      curr = res;
    }
    mpz_get_str(curr, 10, m);
    curr = res + strlen(res);
    if (sizeof(mpfr_exp_t) <= sizeof(int)) {
      if (E < ((mpfr_exp_t) 0)) {
	__sprintfValue_sprintf(curr, " * 2^(%d)", (int) E);
      } else {
	__sprintfValue_sprintf(curr, " * 2^%d", (int) E);
      }
    } else {
      if (sizeof(mpfr_exp_t) <= sizeof(long int)) {
	if (E < ((mpfr_exp_t) 0)) {
	  __sprintfValue_sprintf(curr, " * 2^(%ld)", (long int) E);
	} else {
	  __sprintfValue_sprintf(curr, " * 2^%ld", (long int) E);
	}
      } else {
	if (sizeof(mpfr_exp_t) <= sizeof(long long int)) {
	  if (E < ((mpfr_exp_t) 0)) {
	    __sprintfValue_sprintf(curr, " * 2^(%lld)", (long long int) E);
	  } else {
	    __sprintfValue_sprintf(curr, " * 2^%lld", (long long int) E);
	  }
	} else {
	  tmp = (long long int) E;
	  Et = (mpfr_exp_t) tmp;
	  if (E == Et) {
	    if (E < ((mpfr_exp_t) 0)) {
	      __sprintfValue_sprintf(curr, " * 2^(%lld)", tmp);
	    } else {
	      __sprintfValue_sprintf(curr, " * 2^%lld", tmp);
	    }
	  } else {
	    sollyaFprintf(stderr,"Error: __sprintfValue_print_dyadic: mpfr_exp_t variable too large to be displayed.\n");
	    exit(1);
	  }
	}
      }
    }
  }

  /* Free temporary variable */
  mpz_clear(m);

  /* Return the result */
  return res;
}

static inline mpfr_exp_t __sprintfValue_print_decimal_dec_exponent(mpfr_t x) {
  mp_prec_t prec;
  mpfr_t tmp;
  mpfr_exp_t F;
  intmax_t Fintmax;

  /* Return dummy value for unhandled cases */
  if ((!mpfr_number_p(x)) || (mpfr_zero_p(x))) return 0;

  /* Initialize temporary variable */
  prec = 2 * 8 * sizeof(mpfr_exp_t) + 12;
  if (((mp_prec_t) (2 * 8 * sizeof(intmax_t) + 12)) > prec) prec = 2 * 8 * sizeof(intmax_t) + 12;
  if (mpfr_get_prec(x) > prec) prec = mpfr_get_prec(x);
  mpfr_init2(tmp, prec);

  /* Compute absolute value of x */
  mpfr_abs(tmp, x, GMP_RNDN); /* exact as prec >= prec(x) */

  /* Compute log10(abs(x)) with rounding downwards, precision is
     enough to avoid double rounding
  */
  mpfr_log10(tmp, tmp, GMP_RNDD);

  /* Compute floor(log(abs(x))) with rounding downwards */
  mpfr_rint_floor(tmp, tmp, GMP_RNDD); /* no double rounding possible */

  /* Get machine integer (intmax_t) representation of that integer */
  Fintmax = mpfr_get_sj(tmp, GMP_RNDD); /* exact if not overflowing */

  /* Convert to mpfr_exp_t */
  F = (mpfr_exp_t) Fintmax;

  /* Clear temporary variable */
  mpfr_clear(tmp);

  /* Return the result */
  return F;
}

static inline void __sprintfValue_print_decimal_mul2(mpq_t q, mpfr_exp_t E) {
  mp_bitcnt_t tmp;

  if (E < ((mpfr_exp_t) 0)) {
    tmp = (mp_bitcnt_t) (-E);
    mpz_mul_2exp(mpq_denref(q), mpq_denref(q), tmp);
  } else {
    tmp = (mp_bitcnt_t) E;
    mpz_mul_2exp(mpq_numref(q), mpq_numref(q), tmp);
  }
  mpq_canonicalize(q);
}

static void __sprintfValue_print_decimal_pow5_fallback_pow(mpz_t z, mpz_t b, mp_bitcnt_t E) {
  unsigned int tmp;
  mp_bitcnt_t EE, H, L;
  mpz_t t1, t2;

  tmp = (unsigned int) E;
  EE = (mp_bitcnt_t) tmp;
  if (E == EE) {
    mpz_pow_ui(z,b,tmp);
    return;
  }

  /* This case should actually never happen */
  H = E >> 14;
  L = E - (H << 14);

  /* Here we have

     E = 2^14 * H + L

     and hence

     b^E = b^(2^14 * H + L) = (b^(2^14))^H * b^L

     where we are sure that 2^14 and L will hold on an unsigned int.

  */
  mpz_init(t1);
  mpz_init(t2);
  tmp = (unsigned int) (1 << 14);
  mpz_pow_ui(t1,b,tmp);
  __sprintfValue_print_decimal_pow5_fallback_pow(t1, t1, H);
  tmp = (unsigned int) L;
  mpz_pow_ui(t2,b,tmp);
  mpz_mul(z, t1, t2);
  mpz_clear(t2);
  mpz_clear(t1);
}

static inline void __sprintfValue_print_decimal_pow5_fallback(mpz_t z, mp_bitcnt_t E) {
  mpz_t tmp;

  mpz_init(tmp);
  mpz_set_si(tmp, 5);
  __sprintfValue_print_decimal_pow5_fallback_pow(z, tmp, E);
  mpz_clear(tmp);
}

static inline void __sprintfValue_print_decimal_pow5(mpz_t z, mp_bitcnt_t E) {
  unsigned int tmp;
  mp_bitcnt_t EE;

  tmp = (unsigned int) E;
  EE = (mp_bitcnt_t) tmp;
  if (E == EE) {
    mpz_ui_pow_ui(z,5,tmp);
    return;
  }

  /* Fall-back case that should never happen */
  __sprintfValue_print_decimal_pow5_fallback(z, E);
}

static inline void __sprintfValue_print_decimal_mul5(mpq_t q, mpfr_exp_t E) {
  mpz_t f;
  mp_bitcnt_t tmp;

  mpz_init(f);
  if (E < ((mpfr_exp_t) 0)) {
    tmp = (mp_bitcnt_t) (-E);
  } else {
    tmp = (mp_bitcnt_t) E;
  }
  __sprintfValue_print_decimal_pow5(f, tmp);
  if (E < ((mpfr_exp_t) 0)) {
    mpz_mul(mpq_denref(q), mpq_denref(q), f);
  } else {
    mpz_mul(mpq_numref(q), mpq_numref(q), f);
  }
  mpz_clear(f);
  mpq_canonicalize(q);
}

void sollya_mpq_nearestint(mpq_t rop, mpq_t op) {
  mpz_t q, r;

  /* Initialize two integer variables */
  mpz_init(q);
  mpz_init(r);

  /* Eucliadian division plus correction */
  mpz_fdiv_qr(q, r, mpq_numref(op), mpq_denref(op));
  mpz_mul_2exp(r, r, 1);

  /* Now we have

     op = q + 1/2 * r/b

     where b is the denominator of op.

     Hence

     nearestint(op) = q + nearestint(1/2 * r/b)

     We know that abs(1/2 * r/b) < 1 and
     that r and b have the same sign.

     Therefore 0 <= 1/2 * r/b < 1

     and

     / 1   if r > b
     nearestint(1/2 * r/b) = |
     \ 0   otherwise.


     So, finally,

     / 1    if r > b
     nearestint(op) = q + |
     \ 0    otherwise.

  */
  if (mpz_cmp(r,mpq_denref(op)) > 0) {
    mpz_add_ui(q, q, (unsigned long int) 1);
  }

  /* Here q = nearestint(op).

     Assign q to rop.

  */
  mpq_set_z(rop, q);
  mpq_canonicalize(rop);

  /* Clear the temporaries */
  mpz_clear(r);
  mpz_clear(q);
}

static inline char *__sprintfValue_print_decimal(mpfr_t x) {
  char *res;
  mpfr_exp_t E, F, G, H, Ht, K, D, Kt;
  mpz_t m, n;
  int s;
  mpq_t q;
  mp_prec_t k;
  int akFactor;
  mpq_t ak, bk, deltak, tmp, ten, nr;
  int mIsEven;
  char *curr;
  long long int tmpLongLongInt;
  size_t mantLen, d;

  /* Return NULL on purpose for unhandled cases */
  if (!mpfr_number_p(x)) return NULL;

  /* Handle the case of zero */
  if (mpfr_zero_p(x)) {
    res = safeCalloc(strlen("0") + 1, sizeof(char));
    strcpy(res, "0");
    return res;
  }

  /* Here the input is a non-zero real number

     Decompose x into (-1)^s * 2^E * m

  */
  mpz_init(m);
  E = mpfr_get_z_2exp(m, x);
  s = 0;
  if (mpz_sgn(m) < 0) {
    mpz_neg(m,m);
    s = 1;
  }

  /* Compute the decimal exponent

     F = floor(log10(2^E * m)) = floor(log10(abs(x)))

     We have

     1 <= 10^-F * 2^E * m < 10.

  */
  F = __sprintfValue_print_decimal_dec_exponent(x);

  /* Initialize a mpq_t variable to hold

     q = 2^E / 10^F = 2^G / 5^F

     with G = E - F.

  */
  G = E - F;
  mpq_init(q);
  mpq_set_si(q, 1, (unsigned long int) 1);
  mpq_canonicalize(q);
  __sprintfValue_print_decimal_mul2(q, G);
  __sprintfValue_print_decimal_mul5(q, -F);

  /* We know that

     10^(F - k + 1) * nearestint(2^E/10^F * m * 10^(k - 1))

     is the round-to-nearest decimal representation at k digits of
     our input x (modulo the sign).

     Now let alpha be defined as

     / -0.25    if m is an integer power of 2
     alpha = |
     \ -0.5     otherwise

     and let

     beta = 0.5

     Then 10^(F - k + 1) * nearestint(2^E/10^F * m * 10^(k - 1))

     is the round-to-nearest decimal representation of 2^E * m
     that will round back to 2^E * m iff

     2^E * (m + alpha) <(=) 10^(F - k + 1) * nearestint(2^E/10^F * m * 10^(k - 1)) <(=) 2^E * (m + beta)

     This is equivalent to

     a_k <(=) delta_k <(=) b_k          (1)

     with

     a_k = alpha * 10^(k - 1) * 2^E/10^F,

     b_k = beta * 10^(k - 1) * 2^E/10^F

     and

     delta_k = nearestint(2^E/10^F * m * 10^(k - 1)) - 2^E/10^F * m * 10^(k - 1).

     It is clear that

     a_(k+1) = 10 * a_k

     b_(k+1) = 10 * b_k

     and it can be shown that

     delta_(k+1) = 10 * delta_k - nearestint(10 * delta_k).

     So we can start with some k, say k = 1, and test if
     the condition (1) above is satisfied (in which case we
     know that smallest decimal precision k for which the
     decimal rounding rounds back to binary, recovering x), or
     increment k and update a_k, b_k and delta_k at low cost.

     We shall start with k = 4 and setup a_k, b_k and delta_k.

     For k = 1, we have

     / -25/100 * 2^E/10^F    if m is an integer power of 2
     a_k = |
     \ -50/100 * 2^E/10^F    otherwise

     and

     b_k = 50/100 * 2^E/10^F.

  */
  k = 1;
  mpq_init(ak);
  mpq_init(bk);
  mpq_init(deltak);
  mpq_init(tmp);
  mpq_init(ten);

  /* Test if m > 0 is an integer power of 2

     A positive integer is an integer power of 2 iff its binary
     representation has only one bit set to one.

  */
  if (mpz_popcount(m) == ((mp_bitcnt_t) 1)) {
    akFactor = -25;
  } else {
    akFactor = -50;
  }
  mpq_set_si(tmp, akFactor, (unsigned long int) 100);
  mpq_canonicalize(tmp);
  mpq_mul(ak, q, tmp);

  mpq_set_si(tmp, 50, (unsigned long int) 100);
  mpq_canonicalize(tmp);
  mpq_mul(bk, q, tmp);

  mpq_set_si(tmp, 100, (unsigned long int) 100);
  mpq_canonicalize(tmp);
  mpq_mul(deltak, q, tmp);
  mpq_set_z(tmp, m);
  mpq_canonicalize(tmp);
  mpq_mul(deltak, deltak, tmp);
  sollya_mpq_nearestint(tmp, deltak);
  mpq_sub(deltak, tmp, deltak);

  /* Now determine if m is even or odd. If m is even, we have to test
     condition (1) as

     a_k <= delta_k <= b_k

     otherwise we have to test

     a_k < delta_k < b_k.

  */
  mIsEven = (mpz_tstbit(m, 0) == 0);

  /* Initialize helper variable to 10 */
  mpq_set_si(ten, 10, (unsigned long int) 1);
  mpq_canonicalize(ten);

  /* Loop to determine smallest decimal precision k */
  while (!(mIsEven ? ((
		       mpq_cmp(ak, deltak) <= 0
		       ) &&
		      (
		       mpq_cmp(deltak, bk) <= 0
		       )
		      ) :
	   ((
	     mpq_cmp(ak, deltak) < 0
	     ) &&
	    (
	     mpq_cmp(deltak, bk) < 0
	     )
	    ))) {
    mpq_mul(ak, ak, ten);
    mpq_mul(bk, bk, ten);
    mpq_mul(deltak, deltak, ten);
    sollya_mpq_nearestint(tmp, deltak);
    mpq_sub(deltak, deltak, tmp);
    k++;
  }

  /* Here we know the smallest decimal precision k for which the
     decimal rounding to nearest will round back to binary (with
     round-to-nearest and the original binary precision), resulting
     in the original value of x.

     Now actually compute the decimal rounding:

     10^(F - k + 1) * n

     with

     n = nearestint(2^E/10^F * m * 10^(k - 1))

  */
  mpq_init(nr);
  mpz_init(n);
  mpq_set_z(nr, m);
  mpq_mul(nr, q, nr);
  __sprintfValue_print_decimal_mul2(nr, k - 1);
  __sprintfValue_print_decimal_mul5(nr, k - 1);
  sollya_mpq_nearestint(nr, nr);
  mpz_tdiv_q(n, mpq_numref(nr), mpq_denref(nr));

  /* Now construct the string representation */
  if (mpz_sgn(n) == 0) {
    /* In very rare cases 0 might be the shortest rounding. Produce
       0.0 as output in this case.
    */
    res = safeCalloc(strlen("0.0") + 1, sizeof(char));
    strcpy(res, "0.0");
  } else {
    /* Regular case

       We start by adjusting the exponent F: in some cases, the
       mantissa n exceeds 10^k - 1 because of rounding. In this case,
       n is divisible by 10

    */
    while (mpz_divisible_ui_p(n, (unsigned long int) 10)) {
      mpz_divexact_ui(n, n, (unsigned long int) 10);
      F++;
    }

    /* Now check if there is only one mantissa digit or if there are
       more than that
    */
    if (mpz_cmpabs_ui(n, (unsigned long int) 10) < 0) {
      /* There is only one mantissa digit.

	 Now distinguish the case when the corresponding exponent F -
	 k + 1 is zero or not.

      */
      if ((F - k + 1) == ((mpfr_exp_t) 0)) {
	/* There is only one mantissa digit and no exponent to be
	   displayed. In order to allow for visual differentiation of
	   this case (where there has been rounding) and an integer,
	   we display ".0" after the mantissa digit.
	*/
	res = safeCalloc(1 + mpz_sizeinbase(n, 10) + 2 + 1 + 2 + 1, sizeof(char));
	if (s) {
	  res[0] = '-';
	  curr = res + 1;
	} else {
	  curr = res;
	}
	mpz_get_str(curr, 10, n);
	curr = res + strlen(res);
	strcpy(curr, ".0");
      } else {
	if ((F - k + 1) == ((mpfr_exp_t) -1)) {
	  /* There is only one mantissa digit to be displayed and the exponent is -1.
	     Instead of displaying Xe-1, we display 0.X
	  */
	  res = safeCalloc(1 + mpz_sizeinbase(n, 10) + 2 + 1 + 2 + 1, sizeof(char));
	  if (s) {
	    res[0] = '-';
	    curr = res + 1;
	  } else {
	    curr = res;
	  }
	  strcpy(curr, "0.");
	  curr = res + strlen(res);
	  mpz_get_str(curr, 10, n);
	} else {
	  /* There is only one mantissa digit but an exponent to be
	     displayed
	  */
	  res = safeCalloc(1 + mpz_sizeinbase(n, 10) + 2 + 1 + 1 + sizeof(mpfr_exp_t) * 8 + 1 + 1, sizeof(char));
	  if (s) {
	    res[0] = '-';
	    curr = res + 1;
	  } else {
	    curr = res;
	  }
	  mpz_get_str(curr, 10, n);
	  curr = res + strlen(res);
	  H = F - k + 1;
	  if (sizeof(mpfr_exp_t) <= sizeof(int)) {
	    __sprintfValue_sprintf(curr, "e%d", (int) H);
	  } else {
	    if (sizeof(mpfr_exp_t) <= sizeof(long int)) {
	      __sprintfValue_sprintf(curr, "e%ld", (long int) H);
	    } else {
	      if (sizeof(mpfr_exp_t) <= sizeof(long long int)) {
		__sprintfValue_sprintf(curr, "e%lld", (long long int) H);
	      } else {
		tmpLongLongInt = (long long int) H;
		Ht = (mpfr_exp_t) tmpLongLongInt;
		if (H == Ht) {
		  __sprintfValue_sprintf(curr, "e%lld", tmpLongLongInt);
		} else {
		  sollyaFprintf(stderr,"Error: __sprintfValue_print_decimal: mpfr_exp_t variable too large to be displayed.\n");
		  exit(1);
		}
	      }
	    }
	  }
	}
      }
    } else {
      /* There is more than one mantissa digit to be displayed */
      res = safeCalloc(1 + mpz_sizeinbase(n, 10) + 2 + 1 + 1 + 1 + 1 + 1 + 1 + sizeof(mpfr_exp_t) * 8 + 1 + 1, sizeof(char));
      if (s) {
	res[0] = '-';
	curr = res + 1;
      } else {
	curr = res;
      }
      mpz_get_str(curr, 10, n);
      mantLen = strlen(curr);

      /* Place the decimal point *into* the string of digits unless
	 this makes the exponent field become -1, in which case we
	 place "0." ahead of the string.
      */
      H = F - k + 1;
      D = (mpfr_exp_t) (mantLen - ((size_t) 1));
      K = H + D;
      if ((((mpfr_exp_t) 0) < K) && (K <= ((mpfr_exp_t) 4))) {
	K = (mpfr_exp_t) 0;
      }
      D = K - H;
      if (D < ((mpfr_exp_t) 1)) {
	D = (mpfr_exp_t) 1;
      }
      if (D > ((mpfr_exp_t) (mantLen - ((size_t) 1)))) {
	D = (mpfr_exp_t) (mantLen - ((size_t) 1));
      }
      K = H + D;
      if (K == ((mpfr_exp_t) (-1))) {
	D++;
	K = H + D;
	d = (size_t) D;
	if (d == mantLen) {
	  sollya_memmove((void *) (((curr + mantLen) - d) + 2), (void *) ((curr + mantLen) - d), d * sizeof(char));
	  ((curr + mantLen) - d)[0] = '0';
	  ((curr + mantLen) - d)[1] = '.';
	  curr += mantLen;
	  curr += 2;
	} else {
	  sollya_memmove((void *) (((curr + mantLen) - d) + 1), (void *) ((curr + mantLen) - d), d * sizeof(char));
	  ((curr + mantLen) - d)[0] = '.';
	  curr += mantLen;
	  curr++;
	}
      } else {
	d = (size_t) D;
	sollya_memmove((void *) (((curr + mantLen) - d) + 1), (void *) ((curr + mantLen) - d), d * sizeof(char));
	((curr + mantLen) - d)[0] = '.';
	curr += mantLen;
	curr++;
      }

      /* Now display the exponent field if it needs to be displayed */
      if (K != ((mpfr_exp_t) 0)) {
	if (sizeof(mpfr_exp_t) <= sizeof(int)) {
	  __sprintfValue_sprintf(curr, "e%d", (int) K);
	} else {
	  if (sizeof(mpfr_exp_t) <= sizeof(long int)) {
	    __sprintfValue_sprintf(curr, "e%ld", (long int) K);
	  } else {
	    if (sizeof(mpfr_exp_t) <= sizeof(long long int)) {
	      __sprintfValue_sprintf(curr, "e%lld", (long long int) K);
	    } else {
	      tmpLongLongInt = (long long int) K;
	      Kt = (mpfr_exp_t) tmpLongLongInt;
	      if (K == Kt) {
		__sprintfValue_sprintf(curr, "e%lld", tmpLongLongInt);
	      } else {
		sollyaFprintf(stderr,"Error: __sprintfValue_print_decimal: mpfr_exp_t variable too large to be displayed.\n");
		exit(1);
	      }
	    }
	  }
	}
      }
    }
  }

  /* Free temporary variables */
  mpz_clear(m);
  mpq_clear(q);
  mpq_clear(ak);
  mpq_clear(bk);
  mpq_clear(deltak);
  mpq_clear(tmp);
  mpq_clear(ten);
  mpq_clear(nr);
  mpz_clear(n);

  /* Return the result */
  return res;
}


#define __SPRINTVALUE_DYADIC_DECIMAL      (0)
#define __SPRINTVALUE_DYADIC_DYADIC       (1)
#define __SPRINTVALUE_DYADIC_POWERS       (2)
#define __SPRINTVALUE_DYADIC_BINARY       (3)
#define __SPRINTVALUE_DYADIC_HEXADECIMAL  (4)

char *sprintValue(mpfr_t *valPtr) {
  char *res;
  mpfr_t x;

  /* Initialize our representation x of the value to be displayed */
  __sprintfValue_init_value(x, valPtr);

  /* Display NaNs, infinities and zero */
  if ((!mpfr_number_p(x)) || (mpfr_zero_p(x))) {
    /* Use an auxiliary function */
    res = __sprintfValue_print_special(x);

    /* Clear our representation x of the value to be displayed */
    mpfr_clear(x);

    /* Return the result */
    return res;
  }

  /* Here the input is a non-zero real number

     Now, if the display mode is 'binary numbers' or 'hexadecimal
     numbers' do not use a special code for small integers.

  */
  if (dyadic == __SPRINTVALUE_DYADIC_BINARY) {
    /* Use an auxiliary function */
    res = sPrintBinary(x);

    /* Clear our representation x of the value to be displayed */
    mpfr_clear(x);

    /* Return the result */
    return res;
  }
  if (dyadic == __SPRINTVALUE_DYADIC_HEXADECIMAL) {
    /* Use an auxiliary function */
    res = sPrintHexadecimal(x);

    /* Clear our representation x of the value to be displayed */
    mpfr_clear(x);

    /* Return the result */
    return res;
  }

  /* Here the input is a non-zero real number

     Continue by checking if the number is an integer that holds on
     128 bits.

  */
  if (__sprintfValue_is_small_integer(x)) {
    /* Use an auxiliary function */
    res = __sprintfValue_print_integer(x);

    /* Clear our representation x of the value to be displayed */
    mpfr_clear(x);

    /* Return the result */
    return res;
  }

  /* Here the input is a non-zero real number and the display mode
     should be taken into account
  */
  switch (dyadic) {
  case __SPRINTVALUE_DYADIC_HEXADECIMAL:
    res = sPrintHexadecimal(x);
    break;
  case __SPRINTVALUE_DYADIC_BINARY:
    res = sPrintBinary(x);
    break;
  case __SPRINTVALUE_DYADIC_DYADIC:
    res = __sprintfValue_print_dyadic(x);
    break;
  case __SPRINTVALUE_DYADIC_POWERS:
    res = __sprintfValue_print_powers(x);
    break;
  case __SPRINTVALUE_DYADIC_DECIMAL:
  default:
    res = __sprintfValue_print_decimal(x);
  }

  /* Clear our representation x of the value to be displayed */
  mpfr_clear(x);

  /* Return the result */
  return res;
}

void printMpfr(mpfr_t x) {
  mpfr_t tmp;
  mp_prec_t prec;

  prec = mpfr_get_prec(x);
  mpfr_init2(tmp,prec);
  mpfr_set(tmp,x,GMP_RNDN);

  printValue(&tmp);
  sollyaPrintf("\n");

  mpfr_clear(tmp);
}


void fprintValueWithPrintMode(FILE *fd, mpfr_t value) {
  char *str;
  mpfr_t temp;
  mp_prec_t p;

  p = mpfr_get_prec(value);
  mpfr_init2(temp,p);
  mpfr_set(temp,value,GMP_RNDN);
  str = sprintValue(&temp);
  mpfr_clear(temp);
  sollyaFprintf(fd,"%s",str);
  safeFree(str);

}

void fprintTreeWithPrintMode(FILE *fd, node *tree) {
  int pred, i;

  if (tree->nodeType == MEMREF) {
    fprintTreeWithPrintMode(fd, getMemRefChild(tree));
    return;
  }

  if (fullParentheses) pred = 100; else pred = precedence(tree);

  switch (tree->nodeType) {
  case VARIABLE:
    if (variablename != NULL) {
      sollyaFprintf(fd,"%s",variablename);
    } else {
      sollyaFprintf(fd,"_x_");
    }
    break;
  case CONSTANT:
    fprintValueWithPrintMode(fd,*(tree->value));
    break;
  case ADD:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," + ");
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaFprintf(fd,")");
    break;
  case SUB:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," - ");
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaFprintf(fd,")");
    break;
  case MUL:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," * ");
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaFprintf(fd,")");
    break;
  case DIV:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (tree->child1->nodeType != CONSTANT))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," / ");
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaFprintf(fd,")");
    break;
  case NEG:
    sollyaFprintf(fd,"-");
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaFprintf(fd,")");
    break;
  case UNARY_BASE_FUNC:
    sollyaFprintf(fd,"%s(",tree->baseFun->functionName);
    fprintTreeWithPrintMode(fd,tree->child1);
    sollyaFprintf(fd,")");
    break;
  case POW:
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaFprintf(fd,"(");
    fprintTreeWithPrintMode(fd,tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd,"^");
    if (isInfix(tree->child2) && ((precedence(tree->child2) <= pred)
				  || ((tree->child2->nodeType == CONSTANT)
				      && ((dyadic == 2) || (dyadic == 3))))) {
      sollyaFprintf(fd,"(");
    }
    fprintTreeWithPrintMode(fd,tree->child2);
    if (isInfix(tree->child2) && ((precedence(tree->child2) <= pred)
				  || ((tree->child2->nodeType == CONSTANT)
				      && ((dyadic == 2) || (dyadic == 3))))) {
      sollyaFprintf(fd,")");
    }
    break;
  case LIBRARYFUNCTION:
    {
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,"diff(");
	}
	sollyaFprintf(fd,"%s",tree->libFun->functionName);
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  sollyaFprintf(fd,"%s(",tree->libFun->functionName);
	  fprintTreeWithPrintMode(fd,tree->child1);
	  sollyaFprintf(fd,")");
	} else {
	  sollyaFprintf(fd,"(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,"diff(");
	  }
	  sollyaFprintf(fd,"%s",tree->libFun->functionName);
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,")");
	  }
	  sollyaFprintf(fd,")(");
	  fprintTreeWithPrintMode(fd,tree->child1);
	  sollyaFprintf(fd,")");
	}
      }
    }
    break;
  case PROCEDUREFUNCTION:
    {
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,"diff(");
	}
	sollyaFprintf(fd,"function(");
	fPrintThing(fd,tree->child2);
	sollyaFprintf(fd,")");
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  sollyaFprintf(fd,"(function(");
	  fPrintThing(fd,tree->child2);
	  sollyaFprintf(fd,"))(");
	  fprintTreeWithPrintMode(fd,tree->child1);
	  sollyaFprintf(fd,")");
	} else {
	  sollyaFprintf(fd,"(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,"diff(");
	  }
	  sollyaFprintf(fd,"function(");
	  fPrintThing(fd,tree->child2);
	  sollyaFprintf(fd,")");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,")");
	  }
	  sollyaFprintf(fd,")(");
	  fprintTreeWithPrintMode(fd,tree->child1);
	  sollyaFprintf(fd,")");
	}
      }
    }
    break;
  case PI_CONST:
    sollyaFprintf(fd,"pi");
    break;
  case LIBRARYCONSTANT:
    sollyaFprintf(fd,"%s",tree->libFun->functionName);
    break;
  default:
    sollyaFprintf(stderr,"Error: fprintTreeWithPrintMode: unknown identifier in the tree\n");
    exit(1);
  }
  return;
}






void fprintValue(FILE *fd, mpfr_t value) {
  mpfr_t y;
  char *str, *str2;
  mp_exp_t e, expo;
  mp_prec_t prec;

  if (mpfr_zero_p(value)) {
    sollyaFprintf(fd,"0");
  } else {
    prec = mpfr_get_prec(value);
    mpfr_init2(y,prec+10);
    mpfr_set(y,value,GMP_RNDN);
    if (mpfr_sgn(y) < 0) {
      sollyaFprintf(fd,"-"); mpfr_neg(y,y,GMP_RNDN);
    }
    if (!mpfr_number_p(value)) {
      str = mpfr_get_str(NULL,&e,10,0,y,GMP_RNDN);
      sollyaFprintf(fd,"%s",str);
    } else {
      expo = mpfr_get_exp(y);
      if (mpfr_set_exp(y,prec+10)) {
	printMessage(1,SOLLYA_MSG_OUT_OF_CURRENT_EXPONENT_RANGE,"\nWarning: upon printing to a file: %d is not in the current exponent range of a variable. Values printed may be wrong.\n",(int)(prec+10));
      }
      expo -= prec+10;
      while (mpfr_integer_p(y)) {
	mpfr_div_2ui(y,y,1,GMP_RNDN);
	expo += 1;
      }
      expo--;
      if (mpfr_mul_2ui(y,y,1,GMP_RNDN) != 0) {
	if (!noRoundingWarnings) {
	  printMessage(1,SOLLYA_MSG_INADVERTED_ROUNDING_WHILE_DISPLAYING,"\nWarning: upon printing to a file: rounding occurred. Values printed may be wrong.\n");
	}
      }
      str = mpfr_get_str(NULL,&e,10,0,y,GMP_RNDN);
      str2 = (char *) safeCalloc(strlen(str)+1,sizeof(char));
      strncpy(str2,str,e);
      sollyaFprintf(fd,"%sb%d",str2,(int)expo);
      safeFree(str2);
    }
    safeFree(str);
    mpfr_clear(y);
  }
}

void fprintValueForXml(FILE *fd, mpfr_t value) {
  mpfr_t y, h;
  char *str, *str2;
  mp_exp_t e, expo;
  mp_prec_t prec;
  int negate, val;

  if (mpfr_zero_p(value)) {
    sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> 0 </cn>\n");
  } else {
    prec = mpfr_get_prec(value);
    mpfr_init2(y,prec+10);
    mpfr_set(y,value,GMP_RNDN);
    val = mpfr_get_si(y,GMP_RNDN);
    mpfr_init2(h,prec);
    mpfr_set_si(h,val,GMP_RNDN);
    if (mpfr_number_p(y) && (mpfr_cmp(h,y) == 0)) {
      mpfr_clear(h);
      sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> %d </cn>\n",val);
    } else {
      mpfr_clear(h);
      negate = 0;
      if (mpfr_sgn(y) < 0) {
	negate = 1; mpfr_neg(y,y,GMP_RNDN);
      }
      if (!mpfr_number_p(value)) {
	str = mpfr_get_str(NULL,&e,10,0,y,GMP_RNDN);
	if (!negate)
	  sollyaFprintf(fd,"<cn type=\"real\"> %s </cn>\n",str);
	else
	  sollyaFprintf(fd,"<cn type=\"real\"> -%s </cn>\n",str);
      } else {
	expo = mpfr_get_exp(y);
	if (mpfr_set_exp(y,prec+10)) {
	  printMessage(1,SOLLYA_MSG_OUT_OF_CURRENT_EXPONENT_RANGE,"\nWarning: upon printing to a file: %d is not in the current exponent range of a variable. Values printed may be wrong.\n",(int)(prec+10));
	}
	expo -= prec+10;
	while (mpfr_integer_p(y)) {
	  mpfr_div_2ui(y,y,1,GMP_RNDN);
	  expo += 1;
	}
	expo--;
	if (mpfr_mul_2ui(y,y,1,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_INADVERTED_ROUNDING_WHILE_DISPLAYING,"\nWarning: upon printing to a file: rounding occurred. Values printed may be wrong.\n");
	  }
	}
	str = mpfr_get_str(NULL,&e,10,0,y,GMP_RNDN);
	str2 = (char *) safeCalloc(strlen(str)+1,sizeof(char));
	strncpy(str2,str,e);
	if (!negate) {
	  sollyaFprintf(fd,"<apply>\n");
	  sollyaFprintf(fd,"<times/>\n");
	  sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> %s </cn>\n",str2);
	  sollyaFprintf(fd,"<apply>\n");
	  sollyaFprintf(fd,"<power/>\n");
	  sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> 2 </cn>\n");
	  sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> %d </cn>\n",(int) expo);
	  sollyaFprintf(fd,"</apply>\n");
	  sollyaFprintf(fd,"</apply>\n");
	} else {
	  sollyaFprintf(fd,"<apply>\n");
	  sollyaFprintf(fd,"<times/>\n");
	  sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> -%s </cn>\n",str2);
	  sollyaFprintf(fd,"<apply>\n");
	  sollyaFprintf(fd,"<power/>\n");
	  sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> 2 </cn>\n");
	  sollyaFprintf(fd,"<cn type=\"integer\" base=\"10\"> %d </cn>\n",(int) expo);
	  sollyaFprintf(fd,"</apply>\n");
	  sollyaFprintf(fd,"</apply>\n");
	}
	safeFree(str2);
      }
    }
    mpfr_clear(y);
  }
}



void printTree(node *tree) {
  int pred, i;

  if (tree->nodeType == MEMREF) {
    printTree(getMemRefChild(tree));
    return;
  }

  if (fullParentheses) pred = 100; else pred = precedence(tree);

  switch (tree->nodeType) {
  case VARIABLE:
    if (variablename != NULL) {
      sollyaPrintf("%s",variablename);
    } else {
      sollyaPrintf("_x_");
    }
    break;
  case CONSTANT:
    printValue(tree->value);
    break;
  case ADD:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf("(");
    printTree(tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf(")");
    sollyaPrintf(" + ");
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaPrintf("(");
    printTree(tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaPrintf(")");
    break;
  case SUB:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf("(");
    printTree(tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf(")");
    sollyaPrintf(" - ");
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaPrintf("(");
    printTree(tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaPrintf(")");
    break;
  case MUL:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf("(");
    printTree(tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf(")");
    sollyaPrintf(" * ");
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaPrintf("(");
    printTree(tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      sollyaPrintf(")");
    break;
  case DIV:
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf("(");
    printTree(tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      sollyaPrintf(")");
    sollyaPrintf(" / ");
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaPrintf("(");
    printTree(tree->child2);
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      sollyaPrintf(")");
    break;
  case NEG:
    sollyaPrintf("-");
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaPrintf("(");
    printTree(tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaPrintf(")");
    break;
  case UNARY_BASE_FUNC:
    sollyaPrintf("%s(",tree->baseFun->functionName);
    printTree(tree->child1);
    sollyaPrintf(")");
    break;
  case POW:
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaPrintf("(");
    printTree(tree->child1);
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      sollyaPrintf(")");
    sollyaPrintf("^");
    if (isInfix(tree->child2) && ((precedence(tree->child2) <= pred)
				  || ((accessThruMemRef(tree->child2)->nodeType == CONSTANT)
				      && ((dyadic == 2) || (dyadic == 3))))) {
      sollyaPrintf("(");
    }
    printTree(tree->child2);
    if (isInfix(tree->child2) && ((precedence(tree->child2) <= pred)
				  || ((accessThruMemRef(tree->child2)->nodeType == CONSTANT)
				      && ((dyadic == 2) || (dyadic == 3))))) {
      sollyaPrintf(")");
    }
    break;
  case LIBRARYFUNCTION:
    {
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaPrintf("diff(");
	}
	sollyaPrintf("%s",tree->libFun->functionName);
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaPrintf(")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  sollyaPrintf("%s(",tree->libFun->functionName);
	  printTree(tree->child1);
	  sollyaPrintf(")");
	} else {
	  sollyaPrintf("(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaPrintf("diff(");
	  }
	  sollyaPrintf("%s",tree->libFun->functionName);
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaPrintf(")");
	  }
	  sollyaPrintf(")(");
	  printTree(tree->child1);
	  sollyaPrintf(")");
	}
      }
    }
    break;
  case PROCEDUREFUNCTION:
    {
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaPrintf("diff(");
	}
	sollyaPrintf("function(");
	printThing(tree->child2);
	sollyaPrintf(")");
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaPrintf(")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  sollyaPrintf("(function(");
	  printThing(tree->child2);
	  sollyaPrintf("))(");
	  printTree(tree->child1);
	  sollyaPrintf(")");
	} else {
	  sollyaPrintf("(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaPrintf("diff(");
	  }
	  sollyaPrintf("function(");
	  printThing(tree->child2);
	  sollyaPrintf(")");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaPrintf(")");
	  }
	  sollyaPrintf(")(");
	  printTree(tree->child1);
	  sollyaPrintf(")");
	}
      }
    }
    break;
  case PI_CONST:
    sollyaPrintf("pi");
    break;
  case LIBRARYCONSTANT:
    sollyaPrintf("%s",tree->libFun->functionName);
    break;
  default:
    sollyaFprintf(stderr,"Error: printTree: unknown identifier in the tree\n");
    exit(1);
  }
  return;
}

char *sprintTree(node *tree) {
  int pred, i;
  char *buffer, *buffer1, *buffer2, *finalBuffer, *tempBuf;

  if (tree->nodeType == MEMREF) {
    return sprintTree(getMemRefChild(tree));
  }

  buffer1 = NULL;
  buffer2 = NULL;
  if (fullParentheses) pred = 100; else pred = precedence(tree);
  switch (tree->nodeType) {
  case VARIABLE:
    if (variablename == NULL) {
      buffer = (char *) safeCalloc(4,sizeof(char));
      sprintf(buffer,"_x_");
    } else {
      buffer = (char *) safeCalloc(strlen(variablename)+1,sizeof(char));
      sprintf(buffer,"%s",variablename);
    }
    break;
  case CONSTANT:
    buffer = sprintValue(tree->value);
    break;
  case ADD:
    buffer1 = sprintTree(tree->child1);
    buffer2 = sprintTree(tree->child2);
    buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 9, sizeof(char));
    tempBuf = buffer;
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,")");
    tempBuf += sprintf(tempBuf," + ");
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer2);
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      tempBuf += sprintf(tempBuf,")");
    break;
  case SUB:
    buffer1 = sprintTree(tree->child1);
    buffer2 = sprintTree(tree->child2);
    buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 9, sizeof(char));
    tempBuf = buffer;
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,")");
    tempBuf += sprintf(tempBuf," - ");
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer2);
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      tempBuf += sprintf(tempBuf,")");
    break;
  case MUL:
    buffer1 = sprintTree(tree->child1);
    buffer2 = sprintTree(tree->child2);
    buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 9, sizeof(char));
    tempBuf = buffer;
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,")");
    tempBuf += sprintf(tempBuf," * ");
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer2);
    if (isInfix(tree->child2) && (precedence(tree->child2) < pred))
      tempBuf += sprintf(tempBuf,")");
    break;
  case DIV:
    buffer1 = sprintTree(tree->child1);
    buffer2 = sprintTree(tree->child2);
    buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 9, sizeof(char));
    tempBuf = buffer;
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer1);
    if (isInfix(tree->child1) && (precedence(tree->child1) < pred) && (accessThruMemRef(tree->child1)->nodeType != CONSTANT))
      tempBuf += sprintf(tempBuf,")");
    tempBuf += sprintf(tempBuf," / ");
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer2);
    if (isInfix(tree->child2) && (precedence(tree->child2) <= pred))
      tempBuf += sprintf(tempBuf,")");
    break;
  case NEG:
    buffer1 = sprintTree(tree->child1);
    buffer = (char *) safeCalloc(strlen(buffer1) + 4, sizeof(char));
    if (isInfix(tree->child1)) sprintf(buffer,"-(%s)",buffer1); else sprintf(buffer,"-%s",buffer1);
    break;
  case UNARY_BASE_FUNC:
    buffer1 = sprintTree(tree->child1);
    buffer = (char *) safeCalloc(strlen(buffer1) + strlen(tree->baseFun->functionName) + 3, sizeof(char));
    sprintf(buffer,"%s(%s)",tree->baseFun->functionName, buffer1);
    break;
  case POW:
    buffer1 = sprintTree(tree->child1);
    buffer2 = sprintTree(tree->child2);
    buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 9, sizeof(char));
    tempBuf = buffer;
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer1);
    if (isInfix(tree->child1) && (precedence(tree->child1) <= pred))
      tempBuf += sprintf(tempBuf,")");
    tempBuf += sprintf(tempBuf,"^");
    if (isInfix(tree->child2) && ((precedence(tree->child2) <= pred)
				  || ((accessThruMemRef(tree->child2)->nodeType == CONSTANT)
				      && ((dyadic == 2) || (dyadic == 3)))))
      tempBuf += sprintf(tempBuf,"(");
    tempBuf += sprintf(tempBuf,"%s",buffer2);
    if (isInfix(tree->child2) && ((precedence(tree->child2) <= pred)
				  || ((accessThruMemRef(tree->child2)->nodeType == CONSTANT)
				      && ((dyadic == 2) || (dyadic == 3)))))
      tempBuf += sprintf(tempBuf,")");
    break;
  case LIBRARYFUNCTION:
    {
      buffer1 = sprintTree(tree->child1);
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	buffer = (char *) safeCalloc(strlen(tree->libFun->functionName) + 6 * tree->libFunDeriv + 1, sizeof(char));
	tempBuf = buffer;
	for (i=1;i<=tree->libFunDeriv;i++) {
	  tempBuf += sprintf(tempBuf,"diff(");
	}
	tempBuf += sprintf(tempBuf,"%s",tree->libFun->functionName);
	for (i=1;i<=tree->libFunDeriv;i++) {
	  tempBuf += sprintf(tempBuf,")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  buffer = (char *) safeCalloc(strlen(tree->libFun->functionName) + strlen(buffer1) + 2 + 1, sizeof(char));
	  tempBuf = buffer;
	  tempBuf += sprintf(tempBuf,"%s(",tree->libFun->functionName);
	  tempBuf += sprintf(tempBuf,"%s",buffer1);
	  tempBuf += sprintf(tempBuf,")");
	} else {
	  buffer = (char *) safeCalloc(strlen(tree->libFun->functionName) + strlen(buffer1) + 6 * tree->libFunDeriv + 4 + 1, sizeof(char));
	  tempBuf = buffer;
	  tempBuf += sprintf(tempBuf,"(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    tempBuf += sprintf(tempBuf,"diff(");
	  }
	  tempBuf += sprintf(tempBuf,"%s",tree->libFun->functionName);
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    tempBuf += sprintf(tempBuf,")");
	  }
	  tempBuf += sprintf(tempBuf,")(");
	  tempBuf += sprintf(tempBuf,"%s",buffer1);
	  tempBuf += sprintf(tempBuf,")");
	}
      }
    }
    break;
  case PROCEDUREFUNCTION:
    {
      buffer1 = sprintTree(tree->child1);
      buffer2 = sPrintThing(tree->child2);
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	buffer = (char *) safeCalloc(strlen(buffer2) + 6 * tree->libFunDeriv + 10 + 1, sizeof(char));
	tempBuf = buffer;
	for (i=1;i<=tree->libFunDeriv;i++) {
	  tempBuf += sprintf(tempBuf,"diff(");
	}
	tempBuf += sprintf(tempBuf,"function(");
	tempBuf += sprintf(tempBuf,"%s",buffer2);
	tempBuf += sprintf(tempBuf,")");
	for (i=1;i<=tree->libFunDeriv;i++) {
	  tempBuf += sprintf(tempBuf,")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 14 + 1, sizeof(char));
	  tempBuf = buffer;
	  tempBuf += sprintf(tempBuf,"(function(");
	  tempBuf += sprintf(tempBuf,"%s",buffer2);
	  tempBuf += sprintf(tempBuf,"))(");
	  tempBuf += sprintf(tempBuf,"%s",buffer1);
	  tempBuf += sprintf(tempBuf,")");
	} else {
	  buffer = (char *) safeCalloc(strlen(buffer1) + strlen(buffer2) + 6 * tree->libFunDeriv + 14 + 1, sizeof(char));
	  tempBuf = buffer;
	  tempBuf += sprintf(tempBuf,"(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    tempBuf += sprintf(tempBuf,"diff(");
	  }
	  tempBuf += sprintf(tempBuf,"function(");
	  tempBuf += sprintf(tempBuf,"%s",buffer2);
	  tempBuf += sprintf(tempBuf,")");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    tempBuf += sprintf(tempBuf,")");
	  }
	  tempBuf += sprintf(tempBuf,")(");
	  tempBuf += sprintf(tempBuf,"%s",buffer1);
	  tempBuf += sprintf(tempBuf,")");
	}
      }
    }
    break;
  case PI_CONST:
    buffer = (char *) safeCalloc(3, sizeof(char));
    sprintf(buffer,"pi");
    break;
  case LIBRARYCONSTANT:
    buffer = (char *) safeCalloc(strlen(tree->libFun->functionName) + 1,sizeof(char));
    sprintf(buffer,"%s", tree->libFun->functionName);
    break;
  default:
    sollyaFprintf(stderr,"Error: sprintTree: unknown identifier in the tree\n");
    exit(1);
  }

  finalBuffer = (char *) safeCalloc(strlen(buffer)+1,sizeof(char));
  sprintf(finalBuffer,"%s",buffer);
  safeFree(buffer);
  if (buffer1 != NULL) safeFree(buffer1);
  if (buffer2 != NULL) safeFree(buffer2);
  return finalBuffer;
}


void fprintTree(FILE *fd, node *tree) {
  int i;

  if (tree == NULL) return;
  switch (tree->nodeType) {
  case MEMREF:
    fprintTree(fd, getMemRefChild(tree));
    break;
  case VARIABLE:
    if (variablename == NULL) {
      sollyaFprintf(fd,"_x_");
    } else {
      sollyaFprintf(fd,"%s",variablename);
    }
    break;
  case CONSTANT:
    fprintValue(fd,*(tree->value));
    break;
  case ADD:
    if (isInfix(tree->child1))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child1);
    if (isInfix(tree->child1))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," + ");
    if (isInfix(tree->child2))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child2);
    if (isInfix(tree->child2))
      sollyaFprintf(fd,")");
    break;
  case SUB:
    if (isInfix(tree->child1))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child1);
    if (isInfix(tree->child1))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," - ");
    if (isInfix(tree->child2))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child2);
    if (isInfix(tree->child2))
      sollyaFprintf(fd,")");
    break;
  case MUL:
    if (isInfix(tree->child1))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child1);
    if (isInfix(tree->child1))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," * ");
    if (isInfix(tree->child2))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child2);
    if (isInfix(tree->child2))
      sollyaFprintf(fd,")");
    break;
  case DIV:
    if (isInfix(tree->child1))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child1);
    if (isInfix(tree->child1))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd," / ");
    if (isInfix(tree->child2))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child2);
    if (isInfix(tree->child2))
      sollyaFprintf(fd,")");
    break;
  case NEG:
    sollyaFprintf(fd,"-");
    if (isInfix(tree->child1))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child1);
    if (isInfix(tree->child1))
      sollyaFprintf(fd,")");
    break;
  case UNARY_BASE_FUNC:
    sollyaFprintf(fd,"%s(",tree->baseFun->functionName);
    fprintTree(fd,tree->child1);
    sollyaFprintf(fd,")");
    break;
  case POW:
    if (isInfix(tree->child1))
      sollyaFprintf(fd,"(");
    fprintTree(fd,tree->child1);
    if (isInfix(tree->child1))
      sollyaFprintf(fd,")");
    sollyaFprintf(fd,"^(");
    fprintTree(fd,tree->child2);
    sollyaFprintf(fd,")");
    break;
  case LIBRARYFUNCTION:
    {
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,"diff(");
	}
	sollyaFprintf(fd,"%s",tree->libFun->functionName);
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  sollyaFprintf(fd,"%s(",tree->libFun->functionName);
	  fprintTree(fd,tree->child1);
	  sollyaFprintf(fd,")");
	} else {
	  sollyaFprintf(fd,"(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,"diff(");
	  }
	  sollyaFprintf(fd,"%s",tree->libFun->functionName);
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,")");
	  }
	  sollyaFprintf(fd,")(");
	  fprintTree(fd,tree->child1);
	  sollyaFprintf(fd,")");
	}
      }
    }
    break;
  case PROCEDUREFUNCTION:
    {
      if (accessThruMemRef(tree->child1)->nodeType == VARIABLE) {
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,"diff(");
	}
	sollyaFprintf(fd,"function(");
	fPrintThing(fd,tree->child2);
	sollyaFprintf(fd,")");
	for (i=1;i<=tree->libFunDeriv;i++) {
	  sollyaFprintf(fd,")");
	}
      } else {
	if (tree->libFunDeriv == 0) {
	  sollyaFprintf(fd,"(function(");
	  fPrintThing(fd,tree->child2);
	  sollyaFprintf(fd,"))(");
	  fprintTree(fd,tree->child1);
	  sollyaFprintf(fd,")");
	} else {
	  sollyaFprintf(fd,"(");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,"diff(");
	  }
	  sollyaFprintf(fd,"function(");
	  fPrintThing(fd,tree->child2);
	  sollyaFprintf(fd,")");
	  for (i=1;i<=tree->libFunDeriv;i++) {
	    sollyaFprintf(fd,")");
	  }
	  sollyaFprintf(fd,")(");
	  fprintTree(fd,tree->child1);
	  sollyaFprintf(fd,")");
	}
      }
    }
    break;
  case PI_CONST:
    sollyaFprintf(fd,"pi");
    break;
  case LIBRARYCONSTANT:
    sollyaFprintf(fd,"%s",tree->libFun->functionName);
    break;
  default:
    sollyaFprintf(stderr,"Error: fprintTree: unknown identifier in the tree\n");
    exit(1);
  }
  return;
}

node* copyTreeInner(node *tree);

node* copyTree(node *tree) {
  return addMemRef(copyTreeInner(tree));
}

node* copyTreeInner(node *tree) {
  node *copy;
  mpfr_t *value;
  mp_prec_t prec, p;
  mpfr_t temp;

  if (tree == NULL) return tree;

  if (tree->nodeType == MEMREF) {
    tree->libFunDeriv++;
    return tree;
  }

  switch (tree->nodeType) {
  case VARIABLE:
    copy = makeVariable();
    break;
  case CONSTANT:
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    prec = tools_precision;
    p = mpfr_get_prec(*(tree->value));
    if (p > prec) prec = p;
    mpfr_init2(temp,prec);
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    copy->value = value;
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = copyTreeInner(tree->child1);
    copy->child2 = copyTreeInner(tree->child2);
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = copyTreeInner(tree->child1);
    copy->child2 = copyTreeInner(tree->child2);
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = copyTreeInner(tree->child1);
    copy->child2 = copyTreeInner(tree->child2);
    break;
  case DIV:
    copy = allocateNode();
    copy->nodeType = DIV;
    copy->child1 = copyTreeInner(tree->child1);
    copy->child2 = copyTreeInner(tree->child2);
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = copyTreeInner(tree->child1);
    break;
  case UNARY_BASE_FUNC:
    copy = allocateNode();
    copy->nodeType = UNARY_BASE_FUNC;
    copy->baseFun = tree->baseFun;
    copy->child1 = copyTreeInner(tree->child1);
    break;
  case POW:
    copy = allocateNode();
    copy->nodeType = POW;
    copy->child1 = copyTreeInner(tree->child1);
    copy->child2 = copyTreeInner(tree->child2);
    break;
  case LIBRARYFUNCTION:
    copy = allocateNode();
    copy->nodeType = LIBRARYFUNCTION;
    copy->libFun = tree->libFun;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = copyTreeInner(tree->child1);
    break;
  case PROCEDUREFUNCTION:
    copy = allocateNode();
    copy->nodeType = PROCEDUREFUNCTION;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = copyTreeInner(tree->child1);
    copy->child2 = copyThing(tree->child2);
    break;
  case PI_CONST:
    copy = allocateNode();
    copy->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    copy = allocateNode();
    copy->nodeType = LIBRARYCONSTANT;
    copy->libFun = tree->libFun;
    break;
  default:
    sollyaFprintf(stderr,"Error: copyTreeInner: unknown identifier (%d) in the tree\n", tree->nodeType);
    exit(1);
  }
  return copy;
}

int sollya_mpfr_to_mpq( mpq_t y, mpfr_t x){
  mpz_t mant;
  mp_exp_t expo;
  mpq_t aux;
  if (mpfr_number_p(x)) {
    mpz_init(mant);
    expo = mpfr_get_z_exp(mant,x);
    mpq_init(aux);
    mpq_set_z(aux,mant);

    if (expo>=0)
      mpq_mul_2exp(y,aux,(unsigned int)expo);
    else
      mpq_div_2exp(y,aux,(unsigned int)(-expo));

    mpq_clear(aux);
    mpz_clear(mant);
    return 1;
  }
  else return 0;
}



mp_prec_t getMpzPrecision(mpz_t x) {
  mp_prec_t prec;
  int p, dyadicValue;

  if (mpz_cmp_si(x, 0) == 0) {
    prec = 12;
  } else {
    prec = mpz_sizeinbase(x, 2);
    dyadicValue = mpz_scan1(x, 0);
    p = prec - dyadicValue;
    if (p < 12) prec = 12; else prec = p;
  }

  return prec;
}

int sollya_mpfr_to_mpz( mpz_t y, mpfr_t x){
  mpz_t aux;
  mp_exp_t expo;
  int negative;
  mpfr_t t;

  mpfr_init2(t, mpfr_get_prec(x));
  mpfr_floor(t, x); /* No double rounding as precision the same */
  if (mpfr_number_p(t)) {
    negative = 0;
    if (mpfr_sgn(t) < 0) {
      negative = 1;
      mpfr_neg(t, t, GMP_RNDN); /* exact */
    }

    mpz_init(aux);
    expo = mpfr_get_z_exp(aux,t);

    if (expo>=0)
      mpz_mul_2exp(y,aux,(unsigned int)expo);
    else
      mpz_div_2exp(y,aux,(unsigned int)(-expo));

    if (negative) {
      mpz_neg(y, y);
    }

    mpz_clear(aux);
    mpfr_clear(t);
    return 1;
  } else {
    mpfr_clear(t);
    return 0;
  }
}

void mpz_to_mpfr(mpfr_t res, mpz_t a) {
  mpfr_set_prec(res, getMpzPrecision(a));
  mpfr_set_z(res, a, GMP_RNDN);
}

int mpqHoldsOnMpfr(mpfr_t res, mpq_t a) {
  mpfr_t num, den, tmp;
  mp_prec_t prec, p;
  int tern;

  mpfr_init2(num, 12);
  mpfr_init2(den, 12);
  mpz_to_mpfr(num, mpq_numref(a));
  mpz_to_mpfr(den, mpq_denref(a));
  prec = mpfr_get_prec(num);
  p = mpfr_get_prec(den);
  if (p > prec) prec = p;
  prec += 20;
  mpfr_init2(tmp, prec);
  tern = mpfr_div(tmp, num, den, GMP_RNDN);
  if (tern != 0) {
    mpfr_clear(tmp);
    mpfr_clear(den);
    mpfr_clear(num);
    return 0;
  }
  mpfr_set_prec(res, mpfr_get_prec(tmp));
  mpfr_set(res, tmp, GMP_RNDN);
  mpfr_clear(tmp);
  mpfr_clear(den);
  mpfr_clear(num);

  return 1;
}

static inline int __tryMpzPow(mpz_t res, mpz_t a, mpz_t b) {
  unsigned long int bui;
  uint64_t alog2;

  if (mpz_sgn(a) == 0) {
    if (mpz_sgn(b) == 0) {
      mpz_set_si(res, 1);
      return 1;
    }
    if (mpz_sgn(b) < 0) {
      return 0;
    }
    mpz_set_si(res, 0);
    return 1;
  }

  if (mpz_cmp_si(a, 1) == 0) {
    mpz_set_si(res, 1);
    return 1;
  }

  if (mpz_sgn(b) < 0) {
    return 0;
  }

  if (mpz_sgn(b) == 0) {
    mpz_set_si(res, 1);
    return 1;
  }

  if (!mpz_fits_ulong_p(b)) {
    return 0;
  }

  bui = mpz_get_ui(b);

  alog2 = (uint64_t) mpz_sizeinbase(a,2);
  if (alog2 >= (((uint64_t) 1) << 20)) {
    return 0;
  }

  if (((uint64_t) bui) >= (((uint64_t) 1) << 20)) {
    return 0;
  }

  if ((((uint64_t) bui) * alog2) >= (((uint64_t) 1) << 32)) {
    return 0;
  }

  mpz_pow_ui(res, a, bui);

  return 1;
}

/* res = a^(1/b) */
static inline int __tryMpzInvPow(mpz_t res, mpz_t a, mpz_t b) {
  mpz_t neg_a, tmp1, tmp2, tmp;
  int rr;
  unsigned long int bui;

  if (mpz_sgn(a) < 0) {
    mpz_init(neg_a);
    mpz_neg(neg_a, a);
    mpz_init(tmp1);
    rr = __tryMpzInvPow(tmp1, neg_a, b);
    if (!rr) {
      mpz_clear(neg_a);
      mpz_clear(tmp1);
      return 0;
    }
    mpz_clear(neg_a);
    mpz_init(tmp2);
    mpz_neg(tmp2, tmp1);
    rr = __tryMpzPow(tmp1, tmp2, b);
    if (!rr) {
      mpz_clear(tmp1);
      mpz_clear(tmp2);
      return 0;
    }
    if (mpz_cmp(tmp1, a) != 0) {
      mpz_clear(tmp1);
      mpz_clear(tmp2);
      return 0;
    }
    mpz_set(res, tmp2);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    return 1;
  }

  if (mpz_sgn(b) == 0) {
    return 0;
  }

  if (mpz_cmp_si(a, 1) == 0) {
    mpz_set_si(res, 1);
    return 1;
  }

  if (mpz_sgn(b) <= 0) {
    return 0;
  }

  if (mpz_sgn(a) == 0) {
    mpz_set_si(res, 0);
  }

  if (mpz_cmp_si(b, 1) == 0) {
    mpz_set(res, a);
    return 1;
  }

  /* a >= 2, b >= 2 */
  if (!mpz_fits_ulong_p(b)) {
    return 0;
  }

  bui = mpz_get_ui(b);

  mpz_init(tmp);
  rr = mpz_root(tmp, a, bui);

  if (!rr) {
    mpz_clear(tmp);
    return 0;
  }

  mpz_set(res, tmp);

  return 1;
}

static inline int __tryMpqPow(mpq_t res, mpq_t a, mpq_t b) {
  int rr;
  mpq_t rcp_a, neg_b;
  mpz_t p,q,r,s;

  /* 0^b */
  if (mpq_sgn(a) == 0) {
    /* 0^0 */
    if (mpq_sgn(b) == 0) {
      mpq_set_si(res, 1, 1u);
      return 1;
    }
    /* 0^-p/q */
    if (mpq_sgn(b) < 0) {
      return 0;
    }
    /* 0^p/q */
    mpq_set_si(res, 0, 1u);
    return 1;
  }

  /* a != 0, a^0 */
  if (mpq_sgn(b) == 0) {
    mpq_set_si(res, 1, 1u);
    return 1;
  }

  /* b < 0, a^b = (1/a)^-b */
  if (mpq_sgn(b) < 0) {
    mpq_init(rcp_a);
    mpq_init(neg_b);
    mpq_inv(rcp_a, a);
    mpq_neg(neg_b, b);
    rr = __tryMpqPow(res, rcp_a, neg_b);
    mpq_clear(rcp_a);
    mpq_clear(neg_b);
    return rr;
  }

  /* Here a != 0, b > 0 */
  mpz_init(p);
  mpq_canonicalize(a);
  mpq_canonicalize(b);
  if (!__tryMpzInvPow(p, mpq_numref(a), mpq_denref(b))) {
    mpz_clear(p);
    return 0;
  }
  mpz_init(q);
  if (!__tryMpzInvPow(q, mpq_denref(a), mpq_denref(b))) {
    mpz_clear(p);
    mpz_clear(q);
    return 0;
  }
  mpz_init(r);
  if (!__tryMpzPow(r, p, mpq_numref(b))) {
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(r);
    return 0;
  }
  mpz_init(s);
  if (!__tryMpzPow(s, q, mpq_numref(b))) {
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(s);
    return 0;
  }
  mpz_set(mpq_numref(res),r);
  mpz_set(mpq_denref(res),s);
  mpq_canonicalize(res);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(r);
  mpz_clear(s);

  return 1;
}

static inline int __tryEvaluateConstantTermToMpqSpecialBivariate(mpq_t res, mpq_t a, mpq_t b, int nodeType) {

  switch (nodeType) {
  case ADD:
    mpq_add(res, a, b);
    return 1;
    break;
  case SUB:
    mpq_sub(res, a, b);
    return 1;
    break;
  case MUL:
    mpq_mul(res, a, b);
    return 1;
    break;
  case DIV:
    if (mpq_sgn(b) == 0) return 0;
    mpq_div(res, a, b);
    return 1;
    break;
  case POW:
    return __tryMpqPow(res, a, b);
    break;
  default:
    return 0;
    break;
  }

  return 0;
}

static inline int __tryEvaluateConstantTermToMpqOtherUnivariate(mpq_t res, mpfr_t a, node *t, node *memreftree) {
  mpfr_t y;
  sollya_mpfi_t Y;
  node *tree;

  tree = t;
  if (memreftree != NULL) {
    if (memreftree->nodeType == MEMREF) {
      tree = memreftree;
    }
  }

  switch (t->nodeType) {
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    sollya_mpfi_init2(Y, mpfr_get_prec(a) + 10);
    evaluateConstantExpressionToSharpInterval(Y, tree);
    if (sollya_mpfi_is_point_and_real(Y)) {
      mpfr_init2(y, sollya_mpfi_get_prec(Y));
      sollya_mpfi_get_left(y, Y);
      sollya_mpfi_clear(Y);
      if (mpfr_number_p(y)) {
	sollya_mpfr_to_mpq(res, y);
        mpfr_clear(y);
        return 1;
      }
      mpfr_clear(y);
      return 0;
    }
    sollya_mpfi_clear(Y);
    return 0;
    break;
  default:
    return 0;
    break;
  }

  return 0;
}

static inline int __tryEvaluateConstantTermToMpqInner(mpq_t res, node *tree, node *memreftree) {
  mpq_t c1, c2;
  int r;
  mpfr_t t;
  sollya_mpfi_t tI;
  node *b;

  if (tree == NULL) return 0;

  switch (tree->nodeType) {
  case MEMREF:
    return __tryEvaluateConstantTermToMpqInner(res, getMemRefChild(tree), tree);
    break;
  case CONSTANT:
    if (mpfr_number_p(*(tree->value))) {
      sollya_mpfr_to_mpq(res,*(tree->value));
      return 1;
    }
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    mpq_init(c1);
    if (!tryEvaluateConstantTermToMpq(c1, tree->child1)) {
      mpq_clear(c1);
      return 0;
    }
    mpq_init(c2);
    if (!tryEvaluateConstantTermToMpq(c2, tree->child2)) {
      mpq_clear(c2);
      mpq_clear(c1);
      return 0;
    }
    r = __tryEvaluateConstantTermToMpqSpecialBivariate(res, c1, c2, tree->nodeType);
    mpq_clear(c2);
    mpq_clear(c1);
    return r;
    break;
  case NEG:
  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    mpq_init(c1);
    if (!tryEvaluateConstantTermToMpq(c1, tree->child1)) {
      mpq_clear(c1);
      return 0;
    }
    switch (tree->nodeType) {
    case NEG:
      mpq_neg(res, c1);
      mpq_clear(c1);
      return 1;
      break;
    case UNARY_BASE_FUNC:
      r = tree->baseFun->try_exact_rational_eval(res, c1);
      mpq_clear(c1);
      return r;
      break;
    case LIBRARYFUNCTION:
    case PROCEDUREFUNCTION:
      mpfr_init2(t, 12);
      if (!mpqHoldsOnMpfr(t, c1)) {
        mpfr_clear(t);
        mpq_clear(c1);
        return 0;
      }
      mpq_clear(c1);
      r = __tryEvaluateConstantTermToMpqOtherUnivariate(res, t, tree, memreftree);
      mpfr_clear(t);
      return r;
      break;
    }
    break;
  case PI_CONST:
    return 0;
    break;
  case LIBRARYCONSTANT:
    sollya_mpfi_init2(tI, tools_precision * 10 + 12);
    if (memreftree != NULL) {
      if (memreftree->nodeType == MEMREF) {
	b = memreftree;
      } else {
	b = tree;
      }
    } else {
      b = tree;
    }
    evaluateConstantExpressionToSharpInterval(tI, b);
    if (!sollya_mpfi_is_point_and_real(tI)) {
      sollya_mpfi_clear(tI);
      return 0;
    }
    mpfr_init2(t, sollya_mpfi_get_prec(tI));
    sollya_mpfi_get_left(t, tI);
    sollya_mpfi_clear(tI);
    if (!mpfr_number_p(t)) {
      mpfr_clear(t);
      return 0;
    }
    sollya_mpfr_to_mpq(res, t);
    mpfr_clear(t);
    return 1;
    break;
  default:
    return 0;
  }

  return 0;
}

int tryEvaluateConstantTermToMpq(mpq_t res, node *tree) {
  return __tryEvaluateConstantTermToMpqInner(res, tree, NULL);
}


node *dividePolynomialByPowerOfVariableUnsafe(node *tree, int alpha);

int containsNotANumbers(node * tree) {
  int numberChilds;
  int res;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->containsNotANumbersIsCached) {
      return tree->cache->containsNotANumbersCacheResult;
    }

    res = containsNotANumbers(getMemRefChild(tree));

    if (!tree->cache->containsNotANumbersIsCached) {
      tree->cache->containsNotANumbersCacheResult = res;
      tree->cache->containsNotANumbersIsCached = 1;
    }

    return res;
  }

  if (tree->nodeType == CONSTANT) {
    if (mpfr_nan_p(*(tree->value)))
      return 1;
    else
      return 0;
  }

  if (tree->nodeType == VARIABLE) {
    return 0;
  }

  numberChilds = arity(tree);
  switch (numberChilds) {
  case 0:
    return 0;
    break;
  case 1:
    return containsNotANumbers(tree->child1);
    break;
  case 2:
    return (containsNotANumbers(tree->child1) ||
	    containsNotANumbers(tree->child2));
    break;
  default:
    sollyaFprintf(stderr,"Error: containsNotANumbers: unknown arity of tree node symbol.\n");
    exit(1);
  }

  return 1;
}

int containsOnlyRealNumbers(node * tree) {
  int numberChilds;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) return 1;
    return containsOnlyRealNumbers(getMemRefChild(tree));
  }

  if (tree->nodeType == CONSTANT) {
    if (mpfr_number_p(*(tree->value)))
      return 1;
    else
      return 0;
  }

  if (tree->nodeType == VARIABLE) {
    return 1;
  }

  numberChilds = arity(tree);
  switch (numberChilds) {
  case 0:
    return 1;
    break;
  case 1:
    return containsOnlyRealNumbers(tree->child1);
    break;
  case 2:
    return (containsOnlyRealNumbers(tree->child1) &&
	    containsOnlyRealNumbers(tree->child2));
    break;
  default:
    sollyaFprintf(stderr,"Error: containsOnlyRealNumbers: unknown arity of tree node symbol.\n");
    exit(1);
  }

  return 1;
}

node* simplifyTreeErrorfreeInnerst(node *tree, int rec, int doRational);

node* simplifyTreeErrorfreeInnerNoHookHandling(node *tree, int rec, int doRational) {
  node *res;

  if ((tree != NULL) &&
      ((tree->nodeType == MEMREF) &&
       (tree->cache->simplifyCacheRationalMode >= doRational))) {
    if (tree->cache->simplifyCacheDoesNotSimplify == 1) {
      return copyTree(tree);
    } else {
      if (tree->cache->simplifyCache != NULL) {
	return copyTree(tree->cache->simplifyCache);
      }
    }
  }

  res = addMemRef(simplifyTreeErrorfreeInnerst(tree, rec, doRational));

  if ((tree != NULL) && (res != NULL) &&
      (tree->nodeType == MEMREF)) {
    if (isSyntacticallyEqualCheap(tree,res)) {
      free_memory(res);
      res = copyTree(tree);
      tree->cache->simplifyCacheDoesNotSimplify = 1;
      tree->cache->simplifyCacheRationalMode = doRational;
    } else {
      if (((tree->nodeType == MEMREF) &&
	   (tree->cache->evaluationHook != NULL)) &&
	  ((res->nodeType == MEMREF) &&
	   (res->cache->evaluationHook == NULL))) {
	res->cache->isCorrectlyTyped = tree->cache->isCorrectlyTyped;
	addEvaluationHookFromCopy(&(res->cache->evaluationHook), tree->cache->evaluationHook);
	if ((res->cache->derivCache == NULL) && (tree->cache->derivCache != NULL)) {
	  res->cache->derivCache = copyTree(tree->cache->derivCache);
	}
      }
      if (tree->cache->simplifyCache == NULL) {
	if (res->nodeType == MEMREF) {
	  tree->cache->simplifyCache = copyTree(res);
	  tree->cache->simplifyCacheRationalMode = doRational;
	  tree->cache->simplifyCacheDoesNotSimplify = 0;
	}
      } else {
	if ((tree->cache->simplifyCacheRationalMode >= 0) &&
	    (tree->cache->simplifyCacheRationalMode < doRational)) {
	  free_memory(tree->cache->simplifyCache);
	  tree->cache->simplifyCache = copyTree(res);
	  tree->cache->simplifyCacheRationalMode = doRational;
	  tree->cache->simplifyCacheDoesNotSimplify = 0;
	}
      }
      if (res->nodeType == MEMREF) {
	if ((tree->cache->derivCache != NULL) &&
	    (res->cache->derivCache == NULL)) {
	  res->cache->derivCache = copyTree(tree->cache->derivCache);
	}
	if ((tree->cache->derivUnsimplCache != NULL) &&
	    (res->cache->derivUnsimplCache == NULL)) {
	  res->cache->derivUnsimplCache = copyTree(tree->cache->derivUnsimplCache);
	}
      }
    }
  }

  return res;
}

node* simplifyTreeErrorfreeInner(node *tree, int rec, int doRational) {
  node *res;

  res = simplifyTreeErrorfreeInnerNoHookHandling(tree, rec, doRational);

  if ((res != tree) &&
      (tree->nodeType == MEMREF) &&
      (res->nodeType == MEMREF) &&
      treeContainsHooks(tree) &&
      (!treeContainsHooks(res))) {
    res = rewriteThingWithMemRefReuse(res, tree);
  }

  if ((res != tree) &&
      (tree->nodeType == MEMREF) &&
      (res->nodeType == MEMREF)) {
    copyTreeAnnotationsNoSimplifications(res, tree);
  }

  return res;
}

/*
  Let f the function represented by the expression tree.
  If isNotUniformlyInfinite returns true, f is defined over all real numbers except
  maybe on a discrete set of points (by defined on a point x, we mean that
  f(x) can be given a value in R U {-Inf, +Inf}) and the set of points x where
  |f(x)| = +Inf is discrete.
*/
int isNotUniformlyInfinite(node *tree) {
  sollya_mpfi_t y, x;
  int res;
  mpfr_t yPt;

  /* If f evaluates over the whole real to a closed subset of the
     reals, we know that it is bounded.
  */
  sollya_mpfi_init2(x, 64);
  sollya_mpfi_set_full_range(x);
  sollya_mpfi_init2(y, 64);
  evaluateInterval(y, tree, NULL, x);
  if (sollya_mpfi_has_infinity(y) ||
      sollya_mpfi_has_nan(y)) {
    res = 0;
  } else {
    res = 1;
  }
  sollya_mpfi_clear(y);
  sollya_mpfi_clear(x);
  if (res) return 1;

  /* Consider the set R of rational expressions formed only from x and finite
     constants.
     An easy structural induction on R shows that if r \in R evaluates to a
     finite non-null real number at some point x, it is so for all points x but
     a finite number of them.
  */
  if (isRationalFunction(tree) &&
      containsOnlyRealNumbers(tree)) {
    sollya_mpfi_init2(x, 64);
    sollya_mpfi_set_d(x, 0.7403649888628021);
    sollya_mpfi_init2(y, 64);
    evaluateInterval(y, tree, NULL, x);
    if (sollya_mpfi_has_infinity(y) ||
	sollya_mpfi_has_nan(y) ||
	sollya_mpfi_has_zero(y)) {
      res = 0;
    } else {
      res = 1;
    }
    sollya_mpfi_clear(y);
    sollya_mpfi_clear(x);
    if (res) return res;
  }

  /* If f is of the form f = g +/- h and both g and h are bounded by a
     real, f is bounded by a real.
  */
  switch (accessThruMemRef(tree)->nodeType) {
  case ADD:
  case SUB:
  case MUL:
    return (isNotUniformlyInfinite(accessThruMemRef(tree)->child1) &&
	    isNotUniformlyInfinite(accessThruMemRef(tree)->child2));
    break;
  default:
    break;
  }

  /* If f is of the form f = g / h and both g and h are not uniformly infinite
     and h is not uniformly zero, f is bounded by a real.
  */
  if (accessThruMemRef(tree)->nodeType == DIV) {
    return (isNotUniformlyInfinite(accessThruMemRef(tree)->child1) &&
	    isNotUniformlyZero(accessThruMemRef(tree)->child2));
  }

  /* f is of the form f = g(h), g is defined on the whole real line
     and h is not uniformly infinite, so is f.
  */
  if ( (accessThruMemRef(tree)->nodeType == NEG) ||
       ((accessThruMemRef(tree)->nodeType == UNARY_BASE_FUNC) && (accessThruMemRef(tree)->baseFun->isDefinedEverywhere)) )
    return isNotUniformlyInfinite(accessThruMemRef(tree)->child1);

  /* f is of the form f = g^h, g can be shown to be strictly positive
     and finite over the reals and h is not uniformly infinite, then f is
     not uniformly infinite.
  */
  if (accessThruMemRef(tree)->nodeType == POW) {
    if (!isNotUniformlyInfinite(accessThruMemRef(tree)->child2)) return 0;
    sollya_mpfi_init2(x, 64);
    sollya_mpfi_set_full_range(x);
    sollya_mpfi_init2(y, 64);
    evaluateInterval(y, accessThruMemRef(tree)->child1, NULL, x);
    if (sollya_mpfi_has_infinity(y) ||
	sollya_mpfi_has_nan(y) ||
	sollya_mpfi_has_zero(y) ||
	sollya_mpfi_has_negative_numbers(y)) {
      res = 0;
    } else {
      res = 1;
    }
    sollya_mpfi_clear(y);
    sollya_mpfi_clear(x);
    if (res) return res;
  }

  /* f is of the form f = g^k, g is not uniformly infinite and k is integer
     and non-zero then f is not uniformly infinite.
  */
  if (accessThruMemRef(tree)->nodeType == POW) {
    if (!isNotUniformlyInfinite(accessThruMemRef(tree)->child1)) return 0;
    if (accessThruMemRef(accessThruMemRef(tree)->child2)->nodeType == CONSTANT)
      return (mpfr_integer_p(*(accessThruMemRef(accessThruMemRef(tree)->child2)->value)) &&
	      (mpfr_sgn(*(accessThruMemRef(accessThruMemRef(tree)->child2)->value)) >= 0));
    if (!isConstant(accessThruMemRef(tree)->child2)) return 0;
    sollya_mpfi_init2(y, 64);
    evaluateConstantExpressionToInterval(y, accessThruMemRef(tree)->child2);
    if (sollya_mpfi_is_point_and_real(y)) {
      mpfr_init2(yPt, sollya_mpfi_get_prec(y));
      sollya_mpfi_get_left(yPt, y); /* exact */
      res = (mpfr_integer_p(yPt) && (mpfr_sgn(yPt) >= 0));
      mpfr_clear(yPt);
    } else {
      res = 0;
    }
    sollya_mpfi_clear(y);
    return res;
  }

  /* For all other expressions, we cannot exclude that might evaluate
     to infinity.
  */
  return 0;
}

/*
  Let f the function represented by the expression tree.
  If isNotUniformlyZero returns true, f is defined over all real numbers except
  maybe on a discrete set of points (in particular, f is not identically equal
  to NaN) and the set of its zeros is discrete.
*/
int isNotUniformlyZero(node *tree) {
  sollya_mpfi_t y, x;
  int res;
  mpfr_t yPt;

  /* If f evaluates over the whole real to a closed subset of the real
     not containing zero, we know that it is not uniformly zero.
  */
  sollya_mpfi_init2(x, 64);
  sollya_mpfi_set_full_range(x);
  sollya_mpfi_init2(y, 64);
  evaluateInterval(y, tree, NULL, x);
  if (sollya_mpfi_has_nan(y) ||
      sollya_mpfi_has_zero(y)) {
    res = 0;
  } else {
    res = 1;
  }
  sollya_mpfi_clear(y);
  sollya_mpfi_clear(x);
  if (res) return 1;

  /* Consider the set R of rational expressions formed only from x and finite
     constants.
     An easy structural induction on R shows that if r \in R evaluates to a
     finite non-null real number at some point x, it is so for all points x but
     a finite number of them.
  */
  if (isRationalFunction(tree) &&
      containsOnlyRealNumbers(tree)) {
    sollya_mpfi_init2(x, 64);
    sollya_mpfi_set_d(x, 0.7403649888628021);
    sollya_mpfi_init2(y, 64);
    evaluateInterval(y, tree, NULL, x);
    if (sollya_mpfi_has_infinity(y) ||
	sollya_mpfi_has_nan(y) ||
	sollya_mpfi_has_zero(y)) {
      res = 0;
    } else {
      res = 1;
    }
    sollya_mpfi_clear(y);
    sollya_mpfi_clear(x);
    return res;
  }

  /* If f is of the form f = g * h and both g and h are not uniformly
     zero, f is not uniformly zero.
  */
  if (accessThruMemRef(tree)->nodeType == MUL) {
    return (isNotUniformlyZero(accessThruMemRef(tree)->child1) &&
	    isNotUniformlyZero(accessThruMemRef(tree)->child2));
  }

  /* If f is of the form f = g / h, g is not uniformly zero and h
     stays bounded by a real number and is not uniformly zero, then f
     is not uniformly zero.
  */
  if (accessThruMemRef(tree)->nodeType == DIV) {
    return (isNotUniformlyZero(accessThruMemRef(tree)->child1) &&
	    isNotUniformlyZero(accessThruMemRef(tree)->child2) &&
	    isNotUniformlyInfinite(accessThruMemRef(tree)->child2));
  }

  /* f is of the form f = g(h), g is defined on the whole real line
     and has no zero over the reals or only has a zero at zero and h
     stays both bounded by a real and is not uniformly zero over any
     interval that is not reduced to a point, then f is not uniformly
     zero over any interval that is not reduced to a point.
  */
  if ( (accessThruMemRef(tree)->nodeType == UNARY_BASE_FUNC) &&
       (accessThruMemRef(tree)->baseFun->doesNotVanish) )
    return isNotUniformlyInfinite(accessThruMemRef(tree)->child1);

  if ((accessThruMemRef(tree)->nodeType == NEG) ||
      ((accessThruMemRef(tree)->nodeType == UNARY_BASE_FUNC) &&
       (accessThruMemRef(tree)->baseFun->onlyZeroIsZero) &&
       (accessThruMemRef(tree)->baseFun->isDefinedEverywhere)
       )
      )
    return isNotUniformlyZero(accessThruMemRef(tree)->child1);

  /* f is of the form f = g^h, g can be shown to be strictly positive
     over the reals and h stays both bounded by a real, then f is not
     uniformly zero over any interval that is not reduced to a point.
  */
  if (accessThruMemRef(tree)->nodeType == POW) {
    if (!isNotUniformlyInfinite(accessThruMemRef(tree)->child2)) return 0;
    sollya_mpfi_init2(x, 64);
    sollya_mpfi_set_full_range(x);
    sollya_mpfi_init2(y, 64);
    evaluateInterval(y, accessThruMemRef(tree)->child1, NULL, x);
    if (sollya_mpfi_has_infinity(y) ||
	sollya_mpfi_has_nan(y) ||
	sollya_mpfi_has_zero(y) ||
	sollya_mpfi_has_negative_numbers(y)) {
      res = 0;
    } else {
      res = 1;
    }
    sollya_mpfi_clear(y);
    sollya_mpfi_clear(x);
    if (res) return res;
  }

  /* f is of the form f = g^k, g is not uniformly zero over any
     interval and k is integer and non-zero then f is not uniformly
     zero over any interval that is not reduced to a point.
  */
  if (accessThruMemRef(tree)->nodeType == POW) {
    if (!isNotUniformlyZero(accessThruMemRef(tree)->child1)) return 0;
    if (accessThruMemRef(accessThruMemRef(tree)->child2)->nodeType == CONSTANT)
      return (mpfr_integer_p(*(accessThruMemRef(accessThruMemRef(tree)->child2)->value)) &&
	      (mpfr_sgn(*(accessThruMemRef(accessThruMemRef(tree)->child2)->value)) >= 0));
    if (!isConstant(accessThruMemRef(tree)->child2)) return 0;
    sollya_mpfi_init2(y, 64);
    evaluateConstantExpressionToInterval(y, accessThruMemRef(tree)->child2);
    if (sollya_mpfi_is_point_and_real(y)) {
      mpfr_init2(yPt, sollya_mpfi_get_prec(y));
      sollya_mpfi_get_left(yPt, y); /* exact */
      res = (mpfr_integer_p(yPt) && (mpfr_sgn(yPt) >= 0));
      mpfr_clear(yPt);
    } else {
      res = 0;
    }
    sollya_mpfi_clear(y);
    return res;
  }

  /* For all other expressions, we cannot exclude that they are not uniformly
     zero on any interval subset of the reals that is not reduced to a point.
  */
  return 0;
}

/*
  Let f the function represented by the expression tree.
  If canDoSimplificationSubtraction returns true, f-f is defined over all real
  numbers except maybe on a discrete set of points (in particular, f is not
  identically equal to NaN) and has value zero on all points where it is
  defined.
*/
int canDoSimplificationSubtraction(node *tree) {
  return isNotUniformlyInfinite(tree);
}

/*
  Let f the function represented by the expression tree.
  If canDoSimplificationDivision returns true, f/f is defined over all real
  numbers except maybe on a discrete set of points (in particular, f is not
  identically equal to NaN) and has value one on all points where it is defined.
*/
int canDoSimplificationDivision(node *tree) {
  return (isNotUniformlyInfinite(tree) && isNotUniformlyZero(tree));
}

int isIntegerConstant(node *tree) {
  mpq_t q;
  int s;

  if (tree == NULL) return 0;
  if (!isConstant(tree)) return 0;
  if (accessThruMemRef(tree)->nodeType == CONSTANT) {
    return (mpfr_number_p(*(accessThruMemRef(tree)->value)) && mpfr_integer_p(*(accessThruMemRef(tree)->value)));
  }
  mpq_init(q);
  if (tryEvaluateConstantTermToMpq(q, tree)) {
    if (mpz_divisible_p(mpq_numref(q), mpq_denref(q))) {
      mpq_clear(q);
      return 1;
    } else {
      mpq_clear(q);
      return 0;
    }
  }
  mpq_clear(q);

  switch (tree->nodeType) {
  case MEMREF:
    return isIntegerConstant(getMemRefChild(tree));
    break;
  case UNARY_BASE_FUNC:
    switch (tree->baseFun->baseFunctionCode) {
    case NEARESTINT:
    case FLOOR:
    case CEIL:
      if (isNotUniformlyInfinite(tree->child1)) {
        return 1;
      }
      break;
    default:
      break;
    }
    break;
  case ADD:
  case SUB:
  case MUL:
    if (isIntegerConstant(tree->child1) &&
	isIntegerConstant(tree->child2)) {
      return 1;
    }
    break;
  case POW:
    if (isIntegerConstant(tree->child1) &&
	isIntegerConstant(tree->child2)) {
      if (evaluateSign(&s, tree->child2)) {
	return (s > 0);
      }
    }
    break;
  default:
    break;
  }

  return 0;
}

node* simplifyTreeErrorfreeInnerst(node *tree, int rec, int doRational) {
  node *simplChild1, *simplChild2, *simplified, *recsimplified;
  mpfr_t *value;
  mpfr_t temp;
  mp_prec_t prec, p, pppp;
  int alpha, beta;
  node *temp1, *temp2, *temp3, *temp4;
  mpq_t resMpq;
  mpfr_t num, denom, resDiv, resA, resB;
  int numberChilds;
  int signOkay, sign;
  node *res;
  node *kind;

  if (tree == NULL) return NULL;
  if ((tree->nodeType == MEMREF) &&
      (tree->child1 != NULL) &&
      (tree->cache->polynomialRepresentation == NULL)) {
    kind = getMemRefChild(tree);
    if (accessThruMemRef(kind)->nodeType == CONSTANT) {
      return copyTree(tree);
    }
  }
  if (tree->nodeType == MEMREF) {
    if ((tree->arguments != NULL) &&
	(*((mp_prec_t *) tree->arguments->value) >= 12) &&
	sollya_mpfi_is_point_and_real(*((sollya_mpfi_t *) tree->arguments->next->value))) {
      mpfr_init2(temp, sollya_mpfi_get_prec(*((sollya_mpfi_t *) tree->arguments->next->value)));
      sollya_mpfi_get_left(temp, *((sollya_mpfi_t *) tree->arguments->next->value));
      res = makeConstant(temp);
      mpfr_clear(temp);
      return addMemRef(res);
    }
    if (tree->cache->isConstantIsCached &&
	tree->cache->isConstantCacheResult &&
	(tree->cache->evaluationHook != NULL)) {
      if ((tree->cache->pointEvalCacheY != NULL) &&
	  (tree->cache->pointEvalCacheResultType == POINT_EVAL_EXACT)) {
	res = makeConstant(*(tree->cache->pointEvalCacheY));
	return addMemRef(res);
      }
    }
    if (tree->cache->polynomialRepresentation != NULL) {
      if ((tree->child1 == NULL) || tree->cache->memRefChildFromPolynomial) return copyTree(tree);
      res = addMemRefEvenOnNull(NULL);
      if (res != NULL) {
	res->cache->polynomialRepresentation = polynomialFromCopy(tree->cache->polynomialRepresentation);
	copyTreeAnnotationsNoSimplifications(res, tree);
	return res;
      }
    }
    res = addMemRef(simplifyTreeErrorfreeInner(getMemRefChild(tree), rec, doRational));
    copyTreeAnnotationsNoSimplifications(res, tree);
    return res;
  }

  if ((tree->nodeType == CONSTANT) && (mpfr_nan_p(*(tree->value)))) return copyTree(tree);
  if (tree->nodeType != VARIABLE) {
    numberChilds = arity(tree);
    switch (numberChilds) {
    case 0:
      break;
    case 1:
      if ((accessThruMemRef(tree->child1)->nodeType == CONSTANT) && (mpfr_nan_p(*(accessThruMemRef(tree->child1)->value)))) return copyTree(tree->child1);
      break;
    case 2:
      if ((accessThruMemRef(tree->child1)->nodeType == CONSTANT) && (mpfr_nan_p(*(accessThruMemRef(tree->child1)->value)))) {
	if (isConstant(tree)) return copyTree(tree->child1);
	return copyTree(tree);
      }
      if ((accessThruMemRef(tree->child2)->nodeType == CONSTANT) && (mpfr_nan_p(*(accessThruMemRef(tree->child2)->value)))) {
	if (isConstant(tree)) return copyTree(tree->child2);
	return copyTree(tree);
      }
      break;
    default:
      sollyaFprintf(stderr,"Error: simplifyTreeErrorfreeInnerst: unknown arity of tree node symbol.\n");
      exit(1);
    }
  }

  if (doRational && isConstant(tree) && (tree->nodeType != CONSTANT)) {
    mpq_init(resMpq);
    if (tryEvaluateConstantTermToMpq(resMpq, tree)) {
      mpfr_init2(num,getMpzPrecision(mpq_numref(resMpq)));
      mpfr_init2(denom,getMpzPrecision(mpq_denref(resMpq)));
      mpfr_set_z(num,mpq_numref(resMpq),GMP_RNDN); /* exact */
      mpfr_set_z(denom,mpq_denref(resMpq),GMP_RNDN); /* exact */
      pppp = defaultprecision;
      if (mpfr_get_prec(num) > pppp) pppp = mpfr_get_prec(num);
      if (mpfr_get_prec(denom) > pppp) pppp = mpfr_get_prec(denom);
      mpfr_init2(resDiv,pppp);
      if ((mpfr_div(resDiv,num,denom,GMP_RNDN) == 0) &&
	  mpfr_number_p(resDiv)) {
	mpfr_init2(resA,mpfr_get_prec(resDiv)+10);
	mpfr_set(resA,resDiv,GMP_RNDN); /* exact */
	simplifyMpfrPrec(resA, resDiv);
	simplified = makeConstant(resA);
	mpfr_clear(resA);
      } else {
	mpfr_init2(resA,mpfr_get_prec(num)+10);
	mpfr_set(resA,num,GMP_RNDN); /* exact */
	simplifyMpfrPrec(resA, num);
	mpfr_init2(resB,mpfr_get_prec(denom)+10);
	mpfr_set(resB,denom,GMP_RNDN); /* exact */
	simplifyMpfrPrec(resB, denom);
	simplified = makeDiv(makeConstant(resA),makeConstant(resB));
	mpfr_clear(resA);
	mpfr_clear(resB);
      }

      mpfr_clear(num);
      mpfr_clear(denom);
      mpfr_clear(resDiv);
      mpq_clear(resMpq);
      return simplified;
    }
    mpq_clear(resMpq);
  }

  if ((tree->nodeType == DIV) &&
      (!containsNotANumbers(tree)) &&
      (isPolynomial(tree->child1)) &&
      (isPolynomial(tree->child2)) &&
      ((alpha = getMaxPowerDivider(tree->child1)) > 0) &&
      ((beta = getMaxPowerDivider(tree->child2)) > 0)) {
    if (alpha == beta) {
      temp1 = dividePolynomialByPowerOfVariableUnsafe(tree->child1, alpha);
      temp2 = dividePolynomialByPowerOfVariableUnsafe(tree->child2, alpha);
      temp3 = allocateNode();
      temp3->nodeType = DIV;
      temp3->child1 = temp1;
      temp3->child2 = temp2;
      temp4 = simplifyTreeErrorfreeInner(temp3,rec,doRational);
      free_memory(temp3);
      return temp4;
    } else {
      temp1 = allocateNode();
      temp1->nodeType = DIV;
      temp1->child1 = copyTree(tree->child1);
      temp2 = allocateNode();
      temp1->child2 = temp2;
      temp2->nodeType = POW;
      temp2->child1 = makeVariable();
      temp2->child2 = allocateNode();
      temp2->child2->nodeType = CONSTANT;
      temp2->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(temp2->child2->value),8 * sizeof(int) + 10);
      mpfr_set_si(*(temp2->child2->value),alpha,GMP_RNDN);
      temp3 = simplifyTreeErrorfreeInner(temp1, rec, doRational);
      free_memory(temp1);
      temp1 = allocateNode();
      temp1->nodeType = DIV;
      temp1->child1 = copyTree(tree->child2);
      temp2 = allocateNode();
      temp1->child2 = temp2;
      temp2->nodeType = POW;
      temp2->child1 = makeVariable();
      temp2->child2 = allocateNode();
      temp2->child2->nodeType = CONSTANT;
      temp2->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(temp2->child2->value),8 * sizeof(int) + 10);
      mpfr_set_si(*(temp2->child2->value),beta,GMP_RNDN);
      temp4 = simplifyTreeErrorfreeInner(temp1, rec, doRational);
      free_memory(temp1);
      temp1 = allocateNode();
      temp1->nodeType = DIV;
      temp1->child1 = temp3;
      temp1->child2 = temp4;
      temp2 = simplifyTreeErrorfreeInner(temp1, rec, doRational);
      free_memory(temp1);
      if (alpha > beta) {
	temp1 = allocateNode();
	temp1->nodeType = POW;
	temp1->child1 = makeVariable();
	temp1->child2 = allocateNode();
	temp1->child2->nodeType = CONSTANT;
	temp1->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*(temp1->child2->value),8 * sizeof(int) + 10);
	mpfr_set_si(*(temp1->child2->value),alpha - beta,GMP_RNDN);
      } else {
	temp3 = allocateNode();
	temp3->nodeType = POW;
	temp3->child1 = makeVariable();
	temp3->child2 = allocateNode();
	temp3->child2->nodeType = CONSTANT;
	temp3->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*(temp3->child2->value),8 * sizeof(int) + 10);
	mpfr_set_si(*(temp3->child2->value),beta - alpha,GMP_RNDN);
	temp1 = allocateNode();
	temp1->nodeType = DIV;
	temp1->child2 = temp3;
	temp1->child1 = allocateNode();
	temp1->child1->nodeType = CONSTANT;
	temp1->child1->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*(temp1->child1->value),12);
	mpfr_set_si(*(temp1->child1->value),1,GMP_RNDN);
      }
      temp4 = allocateNode();
      temp4->nodeType = MUL;
      temp4->child1 = temp2;
      temp4->child2 = simplifyTreeErrorfreeInner(temp1,rec, doRational);
      free_memory(temp1);
      return temp4;
    }
  }

  switch (tree->nodeType) {
  case VARIABLE:
    simplified = makeVariable();
    break;
  case CONSTANT:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,mpfr_get_prec(*(tree->value)) + 10);
    simplifyMpfrPrec(temp, *(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    simplified->value = value;
    break;
  case ADD:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplChild2 = simplifyTreeErrorfreeInner(tree->child2,rec, doRational);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      prec = 2 * tools_precision;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild1)->value));
      if (p > prec) prec = p;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild2)->value));
      if (p > prec) prec = p;
      prec += 10;
      if (prec > 256 * tools_precision) prec = 256 * tools_precision;
      mpfr_init2(*value,prec);
      simplified->value = value;
      if ((mpfr_add(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN) != 0) ||
	  (!mpfr_number_p(*value))) {
	simplified->nodeType = ADD;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
	mpfr_clear(*value);
	safeFree(value);
      } else {
	free_memory(simplChild1);
	free_memory(simplChild2);
      }
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	freeNode(simplified);
	simplified = simplChild2;
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  if (accessThruMemRef(simplChild1)->nodeType == NEG) {
	    simplified->nodeType = SUB;
	    simplified->child1 = simplChild2;
	    simplified->child2 = copyTree(accessThruMemRef(simplChild1)->child1);
	    free_memory(simplChild1);
	    if (rec > 0) {
	      recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
	      free_memory(simplified);
	      simplified = recsimplified;
	    }
	  } else {
	    if (accessThruMemRef(simplChild2)->nodeType == NEG) {
	      simplified->nodeType = SUB;
	      simplified->child1 = simplChild1;
	      simplified->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
	      free_memory(simplChild2);
	      if (rec > 0) {
		recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
		free_memory(simplified);
		simplified = recsimplified;
	      }
	    } else {
	      if (isSyntacticallyEqualCheap(simplChild1,simplChild2)) {
		simplified->nodeType = MUL;
		simplified->child1 = allocateNode();
		simplified->child1->nodeType = CONSTANT;
		simplified->child1->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
		mpfr_init2(*(simplified->child1->value),tools_precision);
		mpfr_set_d(*(simplified->child1->value),2.0,GMP_RNDN);
		simplified->child2 = simplChild1;
		free_memory(simplChild2);
		if (rec > 0) {
		  recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
		  free_memory(simplified);
		  simplified = recsimplified;
		}
	      } else {
		simplified->nodeType = ADD;
		simplified->child1 = simplChild1;
		simplified->child2 = simplChild2;
	      }
	    }
	  }
	}
      }
    }
    break;
  case SUB:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplChild2 = simplifyTreeErrorfreeInner(tree->child2,rec, doRational);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      prec = 2 * tools_precision;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild1)->value));
      if (p > prec) prec = p;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild2)->value));
      if (p > prec) prec = p;
      prec += 10;
      if (prec > 256 * tools_precision) prec = 256 * tools_precision;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,prec);
      simplified->value = value;
      if ((mpfr_sub(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN) != 0) ||
	  (!mpfr_number_p(*value))) {
	simplified->nodeType = SUB;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
	mpfr_clear(*value);
	safeFree(value);
      } else {
	free_memory(simplChild1);
	free_memory(simplChild2);
      }
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	simplified->nodeType = NEG;
	simplified->child1 = simplChild2;
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  if (isSyntacticallyEqualCheap(simplChild1,simplChild2) && canDoSimplificationSubtraction(simplChild1)) {
	    free_memory(simplChild1);
	    free_memory(simplChild2);
	    simplified->nodeType = CONSTANT;
	    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	    mpfr_init2(*value,tools_precision);
	    simplified->value = value;
	    mpfr_set_d(*value,0.0,GMP_RNDN);
	  } else {
	    if (accessThruMemRef(simplChild2)->nodeType == NEG) {
	      simplified->nodeType = ADD;
	      simplified->child1 = simplChild1;
	      simplified->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
	      free_memory(simplChild2);
	    } else {
	      if ((accessThruMemRef(simplChild1)->nodeType == UNARY_BASE_FUNC) &&
                  (accessThruMemRef(simplChild1)->baseFun->baseFunctionCode == EXP) &&
		  (accessThruMemRef(simplChild2)->nodeType == CONSTANT) &&
		  (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) &&
		  (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
		simplified->nodeType = UNARY_BASE_FUNC;
                simplified->baseFun = basefun_expm1;
		simplified->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
		free_memory(simplChild1);
		free_memory(simplChild2);
	      } else {
		simplified->nodeType = SUB;
		simplified->child1 = simplChild1;
		simplified->child2 = simplChild2;
	      }
	    }
	  }
	}
      }
    }
    break;
  case MUL:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplChild2 = simplifyTreeErrorfreeInner(tree->child2,rec, doRational);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      prec = 2 * tools_precision;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild1)->value));
      if (p > prec) prec = p;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild2)->value));
      if (p > prec) prec = p;
      prec += 10;
      if (prec > 256 * tools_precision) prec = 256 * tools_precision;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,prec);
      simplified->value = value;
      if ((mpfr_mul(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN) != 0) ||
	  (!mpfr_number_p(*value))) {
	simplified->nodeType = MUL;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
	mpfr_clear(*value);
	safeFree(value);
      } else {
	free_memory(simplChild1);
	free_memory(simplChild2);
      }
    } else {
      if ((((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) ||
	   ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value))))) &&
	  (isNotUniformlyInfinite(simplChild1) && isNotUniformlyInfinite(simplChild2))) {
	free_memory(simplChild1);
	free_memory(simplChild2);
	simplified->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	simplified->value = value;
	mpfr_set_d(*value,0.0,GMP_RNDN);
      } else {
	if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild1)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild1)->value)))) {
	  free_memory(simplChild1);
	  freeNode(simplified);
	  simplified = simplChild2;
	} else {
	  if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	    free_memory(simplChild2);
	    freeNode(simplified);
	    simplified = simplChild1;
	  } else {
	    if ((accessThruMemRef(simplChild1)->nodeType == DIV) &&
		(accessThruMemRef(accessThruMemRef(simplChild1)->child1)->nodeType == CONSTANT) &&
		(mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(simplChild1)->child1)->value),1.0) == 0) &&
		(!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(simplChild1)->child1)->value)))) {
	      simplified->nodeType = DIV;
	      simplified->child1 = simplChild2;
	      simplified->child2 = copyTree(accessThruMemRef(simplChild1)->child2);
	      free_memory(simplChild1);
	    } else {
	      if ((accessThruMemRef(simplChild2)->nodeType == DIV) &&
		  (accessThruMemRef(accessThruMemRef(simplChild2)->child1)->nodeType == CONSTANT) &&
		  (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(simplChild2)->child1)->value),1.0) == 0) &&
		  (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(simplChild2)->child1)->value)))) {
		simplified->nodeType = DIV;
		simplified->child1 = simplChild1;
		simplified->child2 = copyTree(accessThruMemRef(simplChild2)->child2);
		free_memory(simplChild2);
	      } else {
		if ((accessThruMemRef(simplChild1)->nodeType == NEG) &&
		    (accessThruMemRef(simplChild2)->nodeType == NEG)) {
		  simplified->nodeType = MUL;
		  simplified->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
		  simplified->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
		  free_memory(simplChild1);
		  free_memory(simplChild2);
		} else {
		  if (accessThruMemRef(simplChild1)->nodeType == NEG) {
		    simplified->nodeType = NEG;
		    simplified->child1 = allocateNode();
		    simplified->child1->nodeType = MUL;
		    simplified->child1->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
		    simplified->child1->child2 = simplChild2;
		    free_memory(simplChild1);
		  } else {
		    if (accessThruMemRef(simplChild2)->nodeType == NEG) {
		      simplified->nodeType = NEG;
		      simplified->child1 = allocateNode();
		      simplified->child1->nodeType = MUL;
		      simplified->child1->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
		      simplified->child1->child1 = simplChild1;
		      free_memory(simplChild2);
		    } else {
		      if ((accessThruMemRef(simplChild1)->nodeType == MUL) &&
			  (isConstant(simplChild2)) &&
			  (isConstant(accessThruMemRef(simplChild1)->child1))) {
			simplified->nodeType = MUL;
			simplified->child1 = allocateNode();
			simplified->child1->nodeType = MUL;
			simplified->child1->child2 = simplChild2;
			simplified->child1->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
			simplified->child2 = copyTree(accessThruMemRef(simplChild1)->child2);
			free_memory(simplChild1);
			if (rec > 0) {
			  recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
			  free_memory(simplified);
			  simplified = recsimplified;
			}
		      } else {
			if ((accessThruMemRef(simplChild1)->nodeType == MUL) &&
			    (isConstant(simplChild2)) &&
			    (isConstant(accessThruMemRef(simplChild1)->child2))) {
			  simplified->nodeType = MUL;
			  simplified->child1 = allocateNode();
			  simplified->child1->nodeType = MUL;
			  simplified->child1->child1 = simplChild2;
			  simplified->child1->child2 = copyTree(accessThruMemRef(simplChild1)->child2);
			  simplified->child2 = copyTree(accessThruMemRef(simplChild1)->child1);
			  free_memory(simplChild1);
			  if (rec > 0) {
			    recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
			    free_memory(simplified);
			    simplified = recsimplified;
			  }
			} else {
			  if ((accessThruMemRef(simplChild2)->nodeType == MUL) &&
			      (isConstant(simplChild1)) &&
			      (isConstant(accessThruMemRef(simplChild2)->child1))) {
			    simplified->nodeType = MUL;
			    simplified->child2 = allocateNode();
			    simplified->child2->nodeType = MUL;
			    simplified->child2->child2 = simplChild1;
			    simplified->child2->child1 = copyTree(accessThruMemRef(simplChild2)->child1);
			    simplified->child1 = copyTree(accessThruMemRef(simplChild2)->child2);
			    free_memory(simplChild2);
			    if (rec > 0) {
			      recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
			      free_memory(simplified);
			      simplified = recsimplified;
			    }
			  } else {
			    if ((accessThruMemRef(simplChild2)->nodeType == MUL) &&
				(isConstant(simplChild1)) &&
				(isConstant(accessThruMemRef(simplChild2)->child2))) {
			      simplified->nodeType = MUL;
			      simplified->child2 = allocateNode();
			      simplified->child2->nodeType = MUL;
			      simplified->child2->child1 = simplChild1;
			      simplified->child2->child2 = copyTree(accessThruMemRef(simplChild2)->child2);
			      simplified->child1 = copyTree(accessThruMemRef(simplChild2)->child1);
			      free_memory(simplChild2);
			      if (rec > 0) {
				recsimplified = simplifyTreeErrorfreeInner(simplified,rec-1, doRational);
				free_memory(simplified);
				simplified = recsimplified;
			      }
			    } else {
			      if ((accessThruMemRef(simplChild2)->nodeType == DIV) &&
				  (isSyntacticallyEqualCheap(simplChild1,accessThruMemRef(simplChild2)->child2) && canDoSimplificationDivision(simplChild1))) {
				freeNode(simplified);
				free_memory(simplChild1);
				simplified = copyTree(accessThruMemRef(simplChild2)->child1);
				free_memory(simplChild2);
			      } else {
				if ((accessThruMemRef(simplChild1)->nodeType == DIV) &&
				    (isSyntacticallyEqualCheap(simplChild2,accessThruMemRef(simplChild1)->child2) && canDoSimplificationDivision(simplChild2))) {
				  freeNode(simplified);
				  free_memory(simplChild2);
				  simplified = copyTree(accessThruMemRef(simplChild1)->child1);
				  free_memory(simplChild1);
				} else {
				  simplified->nodeType = MUL;
				  simplified->child1 = simplChild1;
				  simplified->child2 = simplChild2;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    break;
  case DIV:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplChild2 = simplifyTreeErrorfreeInner(tree->child2,rec, doRational);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      prec = 2 * tools_precision;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild1)->value));
      if (p > prec) prec = p;
      p = 2 * mpfr_get_prec(*(accessThruMemRef(simplChild2)->value));
      if (p > prec) prec = p;
      prec += 10;
      if (prec > 256 * tools_precision) prec = 256 * tools_precision;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,prec);
      simplified->value = value;
      if ((mpfr_div(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN) != 0) ||
	  (!mpfr_number_p(*value))) {
	simplified->nodeType = DIV;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
	mpfr_clear(*value);
	safeFree(value);
      } else {
	free_memory(simplChild1);
	free_memory(simplChild2);
      }
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	if ((!isConstant(simplChild2)) && (isNotUniformlyZero(simplChild2))) {
	  free_memory(simplChild1);
	  free_memory(simplChild2);
	  simplified->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,tools_precision);
	  simplified->value = value;
	  mpfr_set_d(*value,0.0,GMP_RNDN);
	} else {
	  signOkay = evaluateSign(&sign, simplChild2);
	  if (signOkay && (sign != 0)) {
	    free_memory(simplChild1);
	    free_memory(simplChild2);
	    simplified->nodeType = CONSTANT;
	    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	    mpfr_init2(*value,tools_precision);
	    simplified->value = value;
	    mpfr_set_d(*value,0.0,GMP_RNDN);
	  } else {
	    simplified->nodeType = DIV;
	    simplified->child1 = simplChild1;
	    simplified->child2 = simplChild2;
	  }
	}
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  if (isSyntacticallyEqualCheap(simplChild1,simplChild2) && canDoSimplificationDivision(simplChild1)) {
	    simplified->nodeType = CONSTANT;
	    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	    mpfr_init2(*value,tools_precision);
	    simplified->value = value;
	    mpfr_set_d(*value,1.0,GMP_RNDN);
	    free_memory(simplChild1);
	    free_memory(simplChild2);
	  } else {
	    if ((accessThruMemRef(simplChild1)->nodeType == NEG) &&
		(accessThruMemRef(simplChild2)->nodeType == NEG)) {
	      simplified->nodeType = DIV;
	      simplified->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
	      simplified->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
	      free_memory(simplChild1);
	      free_memory(simplChild2);
	    } else {
	      if (accessThruMemRef(simplChild1)->nodeType == NEG) {
		simplified->nodeType = NEG;
		simplified->child1 = allocateNode();
		simplified->child1->nodeType = DIV;
		simplified->child1->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
		simplified->child1->child2 = simplChild2;
		free_memory(simplChild1);
	      } else {
		if (accessThruMemRef(simplChild2)->nodeType == NEG) {
		  simplified->nodeType = NEG;
		  simplified->child1 = allocateNode();
		  simplified->child1->nodeType = DIV;
		  simplified->child1->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
		  simplified->child1->child1 = simplChild1;
		  free_memory(simplChild2);
		} else {
		  if (accessThruMemRef(simplChild2)->nodeType == DIV) {
		    simplified->nodeType = MUL;
		    simplified->child1 = simplChild1;
		    simplified->child2 = allocateNode();
		    simplified->child2->nodeType = DIV;
		    simplified->child2->child1 = copyTree(accessThruMemRef(simplChild2)->child2);
		    simplified->child2->child2 = copyTree(accessThruMemRef(simplChild2)->child1);
		    free_memory(simplChild2);
		    if (rec > 0) {
		      recsimplified = simplifyTreeErrorfreeInner(simplified,rec, doRational);
		      free_memory(simplified);
		      simplified = recsimplified;
		    }
		  } else {
		    if (accessThruMemRef(simplChild1)->nodeType == DIV) {
		      simplified->nodeType = DIV;
		      simplified->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
		      simplified->child2 = allocateNode();
		      simplified->child2->nodeType = MUL;
		      simplified->child2->child1 = simplChild2;
		      simplified->child2->child2 = copyTree(accessThruMemRef(simplChild1)->child2);
		      free_memory(simplChild1);
		      if (rec > 0) {
			recsimplified = simplifyTreeErrorfreeInner(simplified,rec, doRational);
			free_memory(simplified);
			simplified = recsimplified;
		      }
		    } else {
		      if ((simplChild1->nodeType == UNARY_BASE_FUNC) &&
                          (simplChild1->baseFun->baseFunctionCode == SIN) &&
			  (simplChild2->nodeType == UNARY_BASE_FUNC) &&
                          (simplChild2->baseFun->baseFunctionCode == COS) &&
                          (isSyntacticallyEqualCheap(accessThruMemRef(simplChild1)->child1,accessThruMemRef(simplChild2)->child1))) {
			simplified->nodeType = UNARY_BASE_FUNC;
			simplified->baseFun = basefun_tan;
			simplified->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
			free_memory(simplChild1);
			free_memory(simplChild2);
		      } else {
			simplified->nodeType = DIV;
			simplified->child1 = simplChild1;
			simplified->child2 = simplChild2;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    break;
  case NEG:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      prec = tools_precision;
      p = mpfr_get_prec(*(accessThruMemRef(simplChild1)->value));
      if (p > prec) prec = p;
      prec += 10;
      if (prec > 256 * tools_precision) prec = 256 * tools_precision;
      mpfr_init2(*value,prec);
      simplified->value = value;
      if ((mpfr_neg(*value, *(accessThruMemRef(simplChild1)->value), GMP_RNDN) != 0) ||
	  (mpfr_nan_p(*value))) {
	simplified->nodeType = NEG;
	simplified->child1 = simplChild1;
	mpfr_clear(*value);
	safeFree(value);
      } else {
	free_memory(simplChild1);
      }
    } else {
      if (accessThruMemRef(simplChild1)->nodeType == NEG) {
	freeNode(simplified);
	simplified = copyTree(accessThruMemRef(simplChild1)->child1);
	free_memory(simplChild1);
      } else {
	simplified->nodeType = NEG;
	simplified->child1 = simplChild1;
      }
    }
    break;
  case UNARY_BASE_FUNC:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplified = tree->baseFun->simplify(simplChild1);
    break;
  case POW:
    simplChild1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    simplChild2 = simplifyTreeErrorfreeInner(tree->child2,rec, doRational);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      if ((mpfr_pow(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN) != 0) ||
	  (!mpfr_number_p(*value)) || (!mpfr_number_p(*(accessThruMemRef(simplChild1)->value)))) {
	simplified->nodeType = POW;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
	mpfr_clear(*value);
	safeFree(value);
      } else {
	free_memory(simplChild1);
	free_memory(simplChild2);
      }
    } else {
      if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	freeNode(simplified);
	free_memory(simplChild2);
	simplified = simplChild1;
      } else {
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) &&
	    (accessThruMemRef(simplChild1)->nodeType == POW) &&
	    (accessThruMemRef(accessThruMemRef(simplChild1)->child2)->nodeType == CONSTANT) &&
	    (mpfr_mul(*value,*(accessThruMemRef(simplChild2)->value),*(accessThruMemRef(accessThruMemRef(simplChild1)->child2)->value),GMP_RNDN) == 0)) {
	  simplified->nodeType = POW;
	  simplified->child1 = copyTree(accessThruMemRef(simplChild1)->child1);
	  simplified->child2 = allocateNode();
	  simplified->child2->nodeType = CONSTANT;
	  simplified->child2->value = value;
	  free_memory(simplChild1);
	  free_memory(simplChild2);
	} else {
	  mpfr_clear(*value);
	  safeFree(value);
	  if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) &&
	      (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),2.0) == 0) &&
	      (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value))) &&
	      (accessThruMemRef(simplChild1)->nodeType == UNARY_BASE_FUNC) &&
              (accessThruMemRef(simplChild1)->baseFun->baseFunctionCode == SQRT)) {
            freeNode(simplified);
	    simplified = copyTree(accessThruMemRef(simplChild1)->child1);
	    free_memory(simplChild1);
	    free_memory(simplChild2);
	  } else {
	    simplified->nodeType = POW;
	    simplified->child1 = simplChild1;
	    simplified->child2 = simplChild2;
	  }
	}
      }
    }
    break;
  case LIBRARYFUNCTION:
    simplified = allocateNode();
    simplified->nodeType = LIBRARYFUNCTION;
    simplified->libFun = tree->libFun;
    simplified->libFunDeriv = tree->libFunDeriv;
    simplified->child1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    break;
  case PROCEDUREFUNCTION:
    simplified = allocateNode();
    simplified->nodeType = PROCEDUREFUNCTION;
    simplified->libFunDeriv = tree->libFunDeriv;
    simplified->child2 = copyThing(tree->child2);
    simplified->child1 = simplifyTreeErrorfreeInner(tree->child1,rec, doRational);
    break;
  case PI_CONST:
    simplified = allocateNode();
    simplified->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    simplified = allocateNode();
    simplified->nodeType = LIBRARYCONSTANT;
    simplified->libFun = tree->libFun;
    break;
  default:
    sollyaFprintf(stderr,"Error: simplifyTreeErrorfreeInnerst: unknown identifier in the tree\n");
    exit(1);
  }

  return simplified;
}

node *simplifyRationalErrorfree(node *tree) {
  return simplifyTreeErrorfreeInner(tree,1,1);
}

node *simplifyTreeErrorfree(node *tree) {
  node *temp;

  temp = simplifyTreeErrorfreeInner(tree,1,rationalMode);

  if (verbosity >= 7) {
    if (!isSyntacticallyEqualCheap(temp,tree)) {
      if (verbosity < 9) {
	printMessage(7,SOLLYA_MSG_EXPRESSION_HAS_BEEN_SIMPLIFIED,"Information: an expression has been simplified.\n");
      } else {
	printMessage(9,SOLLYA_MSG_EXPRESSION_HAS_BEEN_SIMPLIFIED_TO_ANOTHER_ONE,"Information: expression '%b' has been simplified to expression '%b'.\n",tree,temp);
      }
    }
  }
  return temp;
}

int isPolynomial(node *tree);
node *differentiatePolynomialUnsafe(node *tree);

node* differentiateUnsimplifiedInner(node *tree);

node* differentiateUnsimplified(node *tree) {
  node *res;

  if ((tree->nodeType == MEMREF) &&
      (tree->cache->derivCache != NULL)) {
    return copyTree(tree->cache->derivCache);
  }

  if ((tree->nodeType == MEMREF) &&
      (tree->cache->derivUnsimplCache != NULL)) {
    return copyTree(tree->cache->derivUnsimplCache);
  }

  res = addMemRef(differentiateUnsimplifiedInner(tree));

  if ((tree->nodeType == MEMREF) &&
      (tree->cache->derivUnsimplCache == NULL) &&
      (res->nodeType == MEMREF)) {
    tree->cache->derivUnsimplCache = copyTree(res);
  }

  return res;
}

int containsNonDifferentiableSubfunctions(node *tree) {

  /* Memory references */
  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      return 0; /* Polynomials are always differentiable */
    }
    return containsNonDifferentiableSubfunctions(getMemRefChild(tree));
  }

  /* Constant expressions are always differentiable */
  if (isConstant(tree)) {
    return 0;
  }

  /* Other cases */
  switch (tree->nodeType) {
  case VARIABLE:
    return 0;
    break;
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    return 0;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    return (containsNonDifferentiableSubfunctions(tree->child1) ||
	    containsNonDifferentiableSubfunctions(tree->child2));
    break;
  case UNARY_BASE_FUNC:
    if (tree->baseFun->isDifferentiableEverywhere)
      return containsNonDifferentiableSubfunctions(tree->child1);
    else return 1;
    break;
  case NEG:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return containsNonDifferentiableSubfunctions(tree->child1);
    break;
  case ABS:
  default:
    sollyaFprintf(stderr,"Error: containsNonDifferentiableSubfunctions: unknown identifier (%d) in the tree\n", tree->nodeType);
    exit(1);
  }
  return 0;
}

node* differentiateUnsimplifiedInner(node *tree) {
  node *derivative;
  mpfr_t *mpfr_temp;
  node *temp_node, *temp_node2, *temp_node3, *f_diff, *g_diff, *f_copy, *g_copy, *g_copy2;
  node *temp_node4, *f_copy2;
  node *temp;
  int deg;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      temp = addMemRefEvenOnNull(NULL);
      if (temp != NULL) {
	temp->cache->polynomialRepresentation = polynomialDeriv(tree->cache->polynomialRepresentation);
	return temp;
      }
    }
    return addMemRef(differentiateUnsimplifiedInner(getMemRefChild(tree)));
  }

  if (isConstant(tree)) {
    mpfr_temp = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*mpfr_temp,tools_precision);
    mpfr_set_d(*mpfr_temp,0.0,GMP_RNDN);
    temp_node = allocateNode();
    temp_node->nodeType = CONSTANT;
    temp_node->value = mpfr_temp;
    derivative = temp_node;
  } else {
    if (isPolynomial(tree) && ((deg = getDegreeSilent(tree)) <= MAXDIFFPOLYSPECIALDEGREE) && (deg >= 0)) {
      derivative = differentiatePolynomialUnsafe(tree);
    } else {

      switch (tree->nodeType) {
      case VARIABLE:
	mpfr_temp = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*mpfr_temp,tools_precision);
	mpfr_set_d(*mpfr_temp,1.0,GMP_RNDN);
	temp_node = allocateNode();
	temp_node->nodeType = CONSTANT;
	temp_node->value = mpfr_temp;
	derivative = temp_node;
	break;
      case CONSTANT:
	mpfr_temp = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*mpfr_temp,tools_precision);
	mpfr_set_d(*mpfr_temp,0.0,GMP_RNDN);
	temp_node = allocateNode();
	temp_node->nodeType = CONSTANT;
	temp_node->value = mpfr_temp;
	derivative = temp_node;
	break;
      case ADD:
	temp_node = allocateNode();
	temp_node->nodeType = ADD;
	temp_node->child1 = differentiateUnsimplified(tree->child1);
	temp_node->child2 = differentiateUnsimplified(tree->child2);
	derivative = temp_node;
	break;
      case SUB:
	temp_node = allocateNode();
	temp_node->nodeType = SUB;
	temp_node->child1 = differentiateUnsimplified(tree->child1);
	temp_node->child2 = differentiateUnsimplified(tree->child2);
	derivative = temp_node;
	break;
      case MUL:
	if (isConstant(tree->child1)) {
	  f_diff = differentiateUnsimplified(tree->child2);
	  g_copy = copyTree(tree->child1);
	  temp_node = allocateNode();
	  temp_node->nodeType = MUL;
	  temp_node->child1 = g_copy;
	  temp_node->child2 = f_diff;
	  derivative = temp_node;
	} else {
	  if (isConstant(tree->child2)) {
	    f_diff = differentiateUnsimplified(tree->child1);
	    g_copy = copyTree(tree->child2);
	    temp_node = allocateNode();
	    temp_node->nodeType = MUL;
	    temp_node->child2 = g_copy;
	    temp_node->child1 = f_diff;
	    derivative = temp_node;
	  } else {
	    f_copy = copyTree(tree->child1);
	    g_copy = copyTree(tree->child2);
	    f_diff = differentiateUnsimplified(tree->child1);
	    g_diff = differentiateUnsimplified(tree->child2);
	    temp_node = allocateNode();
	    temp_node->nodeType = ADD;
	    temp_node2 = allocateNode();
	    temp_node2->nodeType = MUL;
	    temp_node3 = allocateNode();
	    temp_node3->nodeType = MUL;
	    temp_node->child1 = temp_node2;
	    temp_node->child2 = temp_node3;
	    temp_node2->child1 = f_copy;
	    temp_node2->child2 = g_diff;
	    temp_node3->child1 = g_copy;
	    temp_node3->child2 = f_diff;
	    derivative = temp_node;
	  }
	}
	break;
      case DIV:
	f_copy = copyTree(tree->child1);
	g_copy = copyTree(tree->child2);
	f_diff = differentiateUnsimplified(tree->child1);
	g_diff = differentiateUnsimplified(tree->child2);
	temp_node = allocateNode();
	temp_node->nodeType = SUB;
	temp_node2 = allocateNode();
	temp_node2->nodeType = MUL;
	temp_node3 = allocateNode();
	temp_node3->nodeType = MUL;
	temp_node->child1 = temp_node2;
	temp_node->child2 = temp_node3;
	temp_node2->child1 = g_copy;
	temp_node2->child2 = f_diff;
	temp_node3->child1 = f_copy;
	temp_node3->child2 = g_diff;
	g_copy = copyTree(tree->child2);
	temp_node2 = allocateNode();
	temp_node2->nodeType = POW;
	temp_node2->child1 = g_copy;
	temp_node4 = allocateNode();
	temp_node4->nodeType = CONSTANT;
	mpfr_temp = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*mpfr_temp,tools_precision);
	mpfr_set_d(*mpfr_temp,2.0,GMP_RNDN);
	temp_node4->value = mpfr_temp;
	temp_node2->child2 = temp_node4;
	temp_node3 = allocateNode();
	temp_node3->nodeType = DIV;
	temp_node3->child1 = temp_node;
	temp_node3->child2 = temp_node2;
	derivative = temp_node3;
	break;
      case NEG:
	g_diff = differentiateUnsimplified(tree->child1);
	temp_node = allocateNode();
	temp_node->nodeType = NEG;
	temp_node->child1 = g_diff;
	derivative = temp_node;
	break;
      case UNARY_BASE_FUNC:
        derivative = tree->baseFun->diff_expr(tree->child1);
        break;
      case POW:
	if (isConstant(tree->child2)) {
	  g_copy = copyTree(tree->child2);
	  g_copy2 = copyTree(tree->child2);
	  f_diff = differentiateUnsimplified(tree->child1);
	  f_copy = copyTree(tree->child1);
	  temp_node2 = allocateNode();
	  temp_node2->nodeType = MUL;
	  temp_node2->child1 = g_copy;
	  temp_node2->child2 = f_diff;
	  temp_node3 = allocateNode();
	  temp_node3->nodeType = SUB;
	  temp_node4 = allocateNode();
	  temp_node4->nodeType = CONSTANT;
	  mpfr_temp = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*mpfr_temp,tools_precision);
	  mpfr_set_d(*mpfr_temp,1.0,GMP_RNDN);
	  temp_node4->value = mpfr_temp;
	  temp_node3->child2 = temp_node4;
	  temp_node3->child1 = g_copy2;
	  temp_node4 = allocateNode();
	  temp_node4->nodeType = POW;
	  temp_node4->child1 = f_copy;
	  temp_node4->child2 = temp_node3;
	  temp_node = allocateNode();
	  temp_node->nodeType = MUL;
	  temp_node->child1 = temp_node4;
	  temp_node->child2 = temp_node2;
	  derivative = temp_node;
	} else {
	  f_diff = differentiateUnsimplified(tree->child1);
	  f_copy = copyTree(tree->child1);
	  f_copy2 = copyTree(tree->child1);
	  g_copy = copyTree(tree->child2);
	  g_diff = differentiateUnsimplified(tree->child2);
          derivative = makeMul( makePow(copyTree(tree->child1), copyTree(tree->child2)),
                                makeAdd( makeDiv( makeMul(f_diff,  g_copy),
                                                  f_copy2
                                                  ),
                                         makeMul(g_diff, makeLog(f_copy))
                                         )
                                );
        }
	break;
      case LIBRARYFUNCTION:
	g_copy = copyTree(tree->child1);
	g_diff = differentiateUnsimplified(tree->child1);
	temp_node = allocateNode();
	temp_node->nodeType = MUL;
	temp_node->child2 = g_diff;
	temp_node2 = allocateNode();
	temp_node2->nodeType = LIBRARYFUNCTION;
	temp_node2->libFun = tree->libFun;
	temp_node2->libFunDeriv = tree->libFunDeriv + 1;
	temp_node->child1 = temp_node2;
	temp_node2->child1 = g_copy;
	derivative = temp_node;
	break;
      case PROCEDUREFUNCTION:
	g_copy = copyTree(tree->child1);
	g_diff = differentiateUnsimplified(tree->child1);
	temp_node = allocateNode();
	temp_node->nodeType = MUL;
	temp_node->child2 = g_diff;
	temp_node2 = allocateNode();
	temp_node2->nodeType = PROCEDUREFUNCTION;
        temp_node2->child2 = copyThing(tree->child2);
	temp_node2->libFunDeriv = tree->libFunDeriv + 1;
	temp_node->child1 = temp_node2;
	temp_node2->child1 = g_copy;
	derivative = temp_node;
	break;
      case PI_CONST:
      case LIBRARYCONSTANT:
	mpfr_temp = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*mpfr_temp,tools_precision);
	mpfr_set_d(*mpfr_temp,0.0,GMP_RNDN);
	temp_node = allocateNode();
	temp_node->nodeType = CONSTANT;
	temp_node->value = mpfr_temp;
	derivative = temp_node;
	break;
      default:
	sollyaFprintf(stderr,"Error: differentiateUnsimplified: unknown identifier in the tree\n");
	exit(1);
      }
    }
  }
  return addMemRef(derivative);
}

int isHorner(node *);
int isCanonical(node *);

node* differentiateInner(node *tree);

node* differentiate(node *tree) {
  node *res;

  if ((tree->nodeType == MEMREF) &&
      (tree->cache->derivCache != NULL)) {
    return copyTree(tree->cache->derivCache);
  }

  res = addMemRef(differentiateInner(tree));

  if ((tree->nodeType == MEMREF) &&
      (tree->cache->derivCache == NULL) &&
      (res->nodeType == MEMREF)) {
    tree->cache->derivCache = copyTree(res);
  }

  return res;
}

node* differentiateInner(node *tree) {
  node *temp, *temp2, *temp3;

  printMessage(10,SOLLYA_MSG_FORMALLY_DIFFERENTIATING_AN_EXPRESSION,"Information: formally differentiating a function.\n");

  printMessage(11,SOLLYA_MSG_FORMALLY_DIFFERENTIATING_A_PARTICULAR_EXPR,"Information: differentiating the expression '%b'\n",tree);

  if ((tree->nodeType == MEMREF) &&
      (tree->cache->polynomialRepresentation != NULL)) {
    temp = addMemRefEvenOnNull(NULL);
    if (temp != NULL) {
      temp->cache->polynomialRepresentation = polynomialDeriv(tree->cache->polynomialRepresentation);
      return temp;
    }
  }

  if (isPolynomial(tree) && (isHorner(tree) || isCanonical(tree))) {
    temp3 = differentiateUnsimplified(tree);
    temp = simplifyTreeErrorfree(temp3);
    free_memory(temp3);
  } else {
    if ((treeSize(tree) > MAXDIFFSIMPLSIZE) || (getDegree(tree) > MAXDIFFSIMPLDEGREE)) {
      printMessage(7,SOLLYA_MSG_EXPR_TOO_BIG_FOR_SIMPLIFICATION_BEFORE_DIFF,"Information: will not simplify the given expression before differentiating because it is too big.\n");
      temp = differentiateUnsimplified(tree);
    } else {
      temp3 = simplifyTreeErrorfree(tree);
      temp2 = differentiateUnsimplified(temp3);
      temp = simplifyTreeErrorfree(temp2);
      free_memory(temp3);
      free_memory(temp2);
    }
  }
  return temp;
}

node *gcd(node *a, node *b) {
  node *aSimpl, *bSimpl, *res;
  polynomial_t p, q, r;

  if ((a->nodeType == MEMREF) &&
      (b->nodeType == MEMREF) &&
      (a->cache->polynomialRepresentation != NULL) &&
      (b->cache->polynomialRepresentation != NULL)) {
    r = polynomialGcd(a->cache->polynomialRepresentation,
		      b->cache->polynomialRepresentation);
    res = polynomialGetExpression(r);
    polynomialFree(r);
    return res;
  }

  aSimpl = simplifyRationalErrorfree(a);
  bSimpl = simplifyRationalErrorfree(b);

  tryRepresentAsPolynomial(aSimpl);
  tryRepresentAsPolynomial(bSimpl);

  if (polynomialFromExpressionOnlyRealCoeffs(&p, aSimpl)) {
    if (polynomialFromExpressionOnlyRealCoeffs(&q, bSimpl)) {
      r = polynomialGcd(p, q);
      res = polynomialGetExpression(r);
      polynomialFree(r);
      polynomialFree(q);
    } else {
      res = addMemRef(makeConstantInt(1));
    }
    polynomialFree(p);
  } else {
    res = addMemRef(makeConstantInt(1));
  }

  free_memory(bSimpl);
  free_memory(aSimpl);

  return res;
}

node *polydiv(node *a, node *b) {
  node *aSimpl, *bSimpl, *res;
  polynomial_t p, q, rq, rr;

  if ((a->nodeType == MEMREF) &&
      (b->nodeType == MEMREF) &&
      (a->cache->polynomialRepresentation != NULL) &&
      (b->cache->polynomialRepresentation != NULL)) {
    polynomialDivExtended(&rq, &rr,
			  a->cache->polynomialRepresentation,
			  b->cache->polynomialRepresentation);
    res = polynomialGetExpression(rq);
    polynomialFree(rq);
    polynomialFree(rr);
    return res;
  }

  aSimpl = simplifyRationalErrorfree(a);
  bSimpl = simplifyRationalErrorfree(b);

  tryRepresentAsPolynomial(aSimpl);
  tryRepresentAsPolynomial(bSimpl);

  if (polynomialFromExpressionOnlyRealCoeffs(&p, aSimpl)) {
    if (polynomialFromExpressionOnlyRealCoeffs(&q, bSimpl)) {
      polynomialDivExtended(&rq, &rr, p, q);
      res = polynomialGetExpression(rq);
      polynomialFree(rq);
      polynomialFree(rr);
      polynomialFree(q);
    } else {
      res = addMemRef(makeConstantInt(0));
    }
    polynomialFree(p);
  } else {
    res = addMemRef(makeConstantInt(0));
  }

  free_memory(bSimpl);
  free_memory(aSimpl);

  return res;
}

node *polymod(node *a, node *b) {
  node *aSimpl, *bSimpl, *res;
  polynomial_t p, q, rq, rr;

  if ((a->nodeType == MEMREF) &&
      (b->nodeType == MEMREF) &&
      (a->cache->polynomialRepresentation != NULL) &&
      (b->cache->polynomialRepresentation != NULL)) {
    polynomialDivExtended(&rq, &rr,
			  a->cache->polynomialRepresentation,
			  b->cache->polynomialRepresentation);
    res = polynomialGetExpression(rr);
    polynomialFree(rq);
    polynomialFree(rr);
    return res;
  }

  aSimpl = simplifyRationalErrorfree(a);
  bSimpl = simplifyRationalErrorfree(b);

  tryRepresentAsPolynomial(aSimpl);
  tryRepresentAsPolynomial(bSimpl);

  if (polynomialFromExpressionOnlyRealCoeffs(&p, aSimpl)) {
    if (polynomialFromExpressionOnlyRealCoeffs(&q, bSimpl)) {
      polynomialDivExtended(&rq, &rr, p, q);
      res = polynomialGetExpression(rr);
      polynomialFree(rq);
      polynomialFree(rr);
      polynomialFree(q);
    } else {
      res = addMemRef(copyThing(a));
    }
    polynomialFree(p);
  } else {
    res = addMemRef(copyThing(a));
  }

  free_memory(bSimpl);
  free_memory(aSimpl);

  return res;
}

node *eucldiv(node *a, node *b) {
  node *aSimpl, *bSimpl, *res;
  polynomial_t p, q, rq, rr;

  if ((a->nodeType == MEMREF) &&
      (b->nodeType == MEMREF) &&
      (a->cache->polynomialRepresentation != NULL) &&
      (b->cache->polynomialRepresentation != NULL)) {
    polynomialDivExtended(&rq, &rr,
			  a->cache->polynomialRepresentation,
			  b->cache->polynomialRepresentation);
    res = polynomialGetExpression(rq);
    polynomialFree(rq);
    polynomialFree(rr);
    return res;
  }

  aSimpl = simplifyRationalErrorfree(a);
  bSimpl = simplifyRationalErrorfree(b);

  tryRepresentAsPolynomial(aSimpl);
  tryRepresentAsPolynomial(bSimpl);

  if (polynomialFromExpressionOnlyRealCoeffs(&p, aSimpl)) {
    if (polynomialFromExpressionOnlyRealCoeffs(&q, bSimpl)) {
      polynomialDivExtended(&rq, &rr, p, q);
      res = polynomialGetExpression(rq);
      polynomialFree(rq);
      polynomialFree(rr);
      polynomialFree(q);
    } else {
      res = addMemRef(makeConstantInt(0));
    }
    polynomialFree(p);
  } else {
    res = addMemRef(makeConstantInt(0));
  }

  free_memory(bSimpl);
  free_memory(aSimpl);

  return res;
}

node *euclmod(node *a, node *b) {
  node *aSimpl, *bSimpl, *res;
  polynomial_t p, q, rq, rr;

  if ((a->nodeType == MEMREF) &&
      (b->nodeType == MEMREF) &&
      (a->cache->polynomialRepresentation != NULL) &&
      (b->cache->polynomialRepresentation != NULL)) {
    polynomialDivExtended(&rq, &rr,
			  a->cache->polynomialRepresentation,
			  b->cache->polynomialRepresentation);
    res = polynomialGetExpression(rr);
    polynomialFree(rq);
    polynomialFree(rr);
    return res;
  }

  aSimpl = simplifyRationalErrorfree(a);
  bSimpl = simplifyRationalErrorfree(b);

  tryRepresentAsPolynomial(aSimpl);
  tryRepresentAsPolynomial(bSimpl);

  if (polynomialFromExpressionOnlyRealCoeffs(&p, aSimpl)) {
    if (polynomialFromExpressionOnlyRealCoeffs(&q, bSimpl)) {
      polynomialDivExtended(&rq, &rr, p, q);
      res = polynomialGetExpression(rr);
      polynomialFree(rq);
      polynomialFree(rr);
      polynomialFree(q);
    } else {
      res = addMemRef(copyThing(a));
    }
    polynomialFree(p);
  } else {
    res = addMemRef(copyThing(a));
  }

  free_memory(bSimpl);
  free_memory(aSimpl);

  return res;
}

int evaluateConstantExpression(mpfr_t result, node *tree, mp_prec_t prec) {
  mpfr_t cutoff;
  int res;

  if (!isConstant(tree)) return 0;
  mpfr_init2(cutoff, 12);
  mpfr_set_si(cutoff, 0, GMP_RNDN);
  res = evaluateFaithfulWithCutOffFast(result, tree, NULL, cutoff, cutoff, prec);
  if ((res == 0) || (res == 3)) {
    evaluate(result, tree, cutoff, prec);
  }
  mpfr_clear(cutoff);
  return 1;
}

node* simplifyTreeInnerst(node *tree);

node* simplifyTreeInner(node *tree) {
  node *res;

  res = addMemRef(simplifyTreeInnerst(tree));

  if ((tree != NULL) && (res != NULL) &&
      (tree->nodeType == MEMREF) &&
      (!((res->nodeType == MEMREF) && (((tree->child1 != NULL) && (!tree->cache->memRefChildFromPolynomial)) && (res->child1 == NULL)))) &&
      isSyntacticallyEqualCheap(tree,res)) {
    free_memory(res);
    res = copyTree(tree);
  }

  return res;
}

node* simplifyTreeInnerst(node *tree) {
  node *simplChild1, *simplChild2, *simplified;
  mpfr_t *value;
  mpfr_t temp, x, y;
  sollya_mpfi_t tempI;
  int numberChilds;
  node *res;
  node *kind;

  if ((tree->nodeType == MEMREF) &&
      (tree->child1 != NULL) &&
      (tree->cache->polynomialRepresentation == NULL)) {
    kind = getMemRefChild(tree);
    if (accessThruMemRef(kind)->nodeType == CONSTANT) {
      if (mpfr_get_prec(*(accessThruMemRef(kind)->value)) <= tools_precision) {
	return copyTree(tree);
      }
    }
  }
  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      if (((tree->child1 == NULL) || tree->cache->memRefChildFromPolynomial) && polynomialCoefficientsAreDyadic(tree->cache->polynomialRepresentation, 0)) {
	return copyTree(tree);
      }
      res = addMemRefEvenOnNull(NULL);
      if (res != NULL) {
	res->cache->polynomialRepresentation = polynomialRoundDyadic(tree->cache->polynomialRepresentation,
								     tools_precision);
	return res;
      }
    }
    return addMemRef(simplifyTreeInner(getMemRefChild(tree)));
  }

  if ((tree->nodeType == CONSTANT) && (mpfr_nan_p(*(tree->value)))) return copyTree(tree);
  if (tree->nodeType != VARIABLE) {
    numberChilds = arity(tree);
    switch (numberChilds) {
    case 0:
      break;
    case 1:
      if ((accessThruMemRef(tree->child1)->nodeType == CONSTANT) && (mpfr_nan_p(*(accessThruMemRef(tree->child1)->value)))) return copyTree(tree);
      break;
    case 2:
      if ((accessThruMemRef(tree->child1)->nodeType == CONSTANT) && (mpfr_nan_p(*(accessThruMemRef(tree->child1)->value)))) {
	if (isConstant(tree)) return copyTree(tree->child1);
	return copyTree(tree);
      }
      if ((accessThruMemRef(tree->child2)->nodeType == CONSTANT) && (mpfr_nan_p(*(accessThruMemRef(tree->child2)->value)))) {
	if (isConstant(tree)) return copyTree(tree->child2);
	return copyTree(tree);
      }
      break;
    default:
      sollyaFprintf(stderr,"Error: simplifyTreeInner: unknown arity of tree node symbol.\n");
      exit(1);
    }
  }

  if (isConstant(tree) && (tree->nodeType != CONSTANT)) {
    mpfr_init2(x, tools_precision);
    mpfr_init2(y, tools_precision);
    mpfr_set_si(x, 1, GMP_RNDN); /* exact */
    if (evaluateFaithful(y, tree, x, tools_precision)) {
      simplified = makeConstant(y);
      mpfr_clear(y);
      mpfr_clear(x);
      return simplified;
    }
    mpfr_clear(y);
    mpfr_clear(x);
  }

  switch (tree->nodeType) {
  case VARIABLE:
    simplified = makeVariable();
    break;
  case CONSTANT:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,tools_precision);
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    simplified->value = value;
    break;
  case ADD:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplChild2 = simplifyTreeInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_add(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	freeNode(simplified);
	simplified = simplChild2;
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  simplified->nodeType = ADD;
	  simplified->child1 = simplChild1;
	  simplified->child2 = simplChild2;
	}
      }
    }
    break;
  case SUB:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplChild2 = simplifyTreeInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_sub(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	simplified->nodeType = NEG;
	simplified->child1 = simplChild2;
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  simplified->nodeType = SUB;
	  simplified->child1 = simplChild1;
	  simplified->child2 = simplChild2;
	}
      }
    }
    break;
  case MUL:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplChild2 = simplifyTreeInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_mul(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if (((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) ||
	  ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value))))) {
	free_memory(simplChild1);
	free_memory(simplChild2);
	simplified->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	simplified->value = value;
	mpfr_set_d(*value,0.0,GMP_RNDN);
      } else {
	if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild1)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild1)->value)))) {
	  free_memory(simplChild1);
	  freeNode(simplified);
	  simplified = simplChild2;
	} else {
	  if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	    free_memory(simplChild2);
	    freeNode(simplified);
	    simplified = simplChild1;
	  } else {
	    simplified->nodeType = MUL;
	    simplified->child1 = simplChild1;
	    simplified->child2 = simplChild2;
	  }
	}
      }
    }
    break;
  case DIV:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplChild2 = simplifyTreeInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_div(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	free_memory(simplChild2);
	simplified->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	simplified->value = value;
	mpfr_set_d(*value,0.0,GMP_RNDN);
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0)) {
	  if (mpfr_nan_p(*(accessThruMemRef(simplChild2)->value))) {
	    simplified->nodeType = MUL;
	    simplified->child1 = simplChild1;
	    simplified->child2 = simplChild2;
	  } else {
	    free_memory(simplChild2);
	    freeNode(simplified);
	    simplified = simplChild1;
	  }
	} else {
	  simplified->nodeType = DIV;
	  simplified->child1 = simplChild1;
	  simplified->child2 = simplChild2;
	}
      }
    }
    break;
  case NEG:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_neg(*value, *(accessThruMemRef(simplChild1)->value), GMP_RNDN);
      free_memory(simplChild1);
    } else {
      simplified->nodeType = NEG;
      simplified->child1 = simplChild1;
    }
    break;
  case UNARY_BASE_FUNC:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      tree->baseFun->point_eval(*value, *(accessThruMemRef(simplChild1)->value));
      free_memory(simplChild1);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = tree->baseFun;
      simplified->child1 = simplChild1;
    }
    break;
  case POW:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplChild2 = simplifyTreeInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_pow(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	freeNode(simplified);
	free_memory(simplChild2);
	simplified = simplChild1;
      } else {
	simplified->nodeType = POW;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
      }
    }
    break;
  case LIBRARYFUNCTION:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      if (tree->libFun->hasData) {
	sollya_mpfr_from_mpfi_data(*value, *(accessThruMemRef(simplChild1)->value), tree->libFunDeriv, tree->libFun->code, tree->libFun->data);
      } else {
	sollya_mpfr_from_mpfi(*value, *(accessThruMemRef(simplChild1)->value), tree->libFunDeriv, tree->libFun->code);
      }
      free_memory(simplChild1);
    } else {
      simplified->nodeType = LIBRARYFUNCTION;
      simplified->child1 = simplChild1;
      simplified->libFun = tree->libFun;
      simplified->libFunDeriv = tree->libFunDeriv;
    }
    break;
  case PROCEDUREFUNCTION:
    simplChild1 = simplifyTreeInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      computeFunctionWithProcedureMpfr(*value, tree->child2, *(accessThruMemRef(simplChild1)->value), (unsigned int) tree->libFunDeriv);
      free_memory(simplChild1);
    } else {
      simplified->nodeType = PROCEDUREFUNCTION;
      simplified->child1 = simplChild1;
      simplified->child2 = copyThing(tree->child2);
      simplified->libFunDeriv = tree->libFunDeriv;
    }
    break;
  case PI_CONST:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    mpfr_const_pi(*value,GMP_RNDN);
    simplified->value = value;
    break;
  case LIBRARYCONSTANT:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    sollya_mpfi_init2(tempI, tools_precision);
    libraryConstantToInterval(tempI, tree);
    sollya_mpfi_get_fr(*value, tempI);
    sollya_mpfi_clear(tempI);
    simplified->value = value;
    break;
  default:
    sollyaFprintf(stderr,"Error: simplifyTreeInnerst: unknown identifier (%d) in the tree\n",tree->nodeType);
    exit(1);
  }

  return simplified;
}

node* simplifyAllButDivisionInnerst(node *tree);

node* simplifyAllButDivisionInner(node *tree) {
  node *res;

  res = addMemRef(simplifyAllButDivisionInnerst(tree));

  if ((tree != NULL) && (res != NULL) &&
      (tree->nodeType == MEMREF) &&
      isSyntacticallyEqualCheap(tree,res)) {
    free_memory(res);
    res = copyTree(tree);
  }

  return res;
}

node* simplifyAllButDivisionInnerst(node *tree) {
  node *simplChild1, *simplChild2, *simplified;
  mpfr_t *value;
  mpfr_t temp;
  sollya_mpfi_t tempI;
  node *res;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      if (((tree->child1 == NULL) || tree->cache->memRefChildFromPolynomial) && polynomialCoefficientsAreRational(tree->cache->polynomialRepresentation, 0))
	return copyTree(tree);
      res = addMemRefEvenOnNull(NULL);
      if (res != NULL) {
	res->cache->polynomialRepresentation = polynomialRoundRational(tree->cache->polynomialRepresentation,
								       tools_precision);
	return res;
      }
    }
    return addMemRef(simplifyAllButDivisionInner(getMemRefChild(tree)));
  }

  switch (tree->nodeType) {
  case VARIABLE:
    simplified = makeVariable();
    break;
  case CONSTANT:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,((tools_precision > mpfr_get_prec(*(tree->value))) ? (tools_precision) : (mpfr_get_prec(*(tree->value)))));
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    simplified->value = value;
    break;
  case ADD:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplChild2 = simplifyAllButDivisionInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_add(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	freeNode(simplified);
	simplified = simplChild2;
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  simplified->nodeType = ADD;
	  simplified->child1 = simplChild1;
	  simplified->child2 = simplChild2;
	}
      }
    }
    break;
  case SUB:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplChild2 = simplifyAllButDivisionInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_sub(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) {
	free_memory(simplChild1);
	simplified->nodeType = NEG;
	simplified->child1 = simplChild2;
      } else {
	if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value)))) {
	  free_memory(simplChild2);
	  freeNode(simplified);
	  simplified = simplChild1;
	} else {
	  simplified->nodeType = SUB;
	  simplified->child1 = simplChild1;
	  simplified->child2 = simplChild2;
	}
      }
    }
    break;
  case MUL:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplChild2 = simplifyAllButDivisionInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_mul(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if (((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild1)->value)))) ||
	  ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_zero_p(*(accessThruMemRef(simplChild2)->value))))) {
	free_memory(simplChild1);
	free_memory(simplChild2);
	simplified->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	simplified->value = value;
	mpfr_set_d(*value,0.0,GMP_RNDN);
      } else {
	if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild1)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild1)->value)))) {
	  free_memory(simplChild1);
	  freeNode(simplified);
	  simplified = simplChild2;
	} else {
	  if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	    free_memory(simplChild2);
	    freeNode(simplified);
	    simplified = simplChild1;
	  } else {
	    simplified->nodeType = MUL;
	    simplified->child1 = simplChild1;
	    simplified->child2 = simplChild2;
	  }
	}
      }
    }
    break;
  case DIV:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplChild2 = simplifyAllButDivisionInner(tree->child2);
    simplified = allocateNode();
    simplified->nodeType = DIV;
    simplified->child1 = simplChild1;
    simplified->child2 = simplChild2;
    break;
  case NEG:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_neg(*value, *(accessThruMemRef(simplChild1)->value), GMP_RNDN);
      free_memory(simplChild1);
    } else {
      simplified->nodeType = NEG;
      simplified->child1 = simplChild1;
    }
    break;
  case UNARY_BASE_FUNC:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      tree->baseFun->point_eval(*value, *(accessThruMemRef(simplChild1)->value));
      free_memory(simplChild1);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = tree->baseFun;
      simplified->child1 = simplChild1;
    }
    break;
  case POW:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplChild2 = simplifyAllButDivisionInner(tree->child2);
    simplified = allocateNode();
    if ((accessThruMemRef(simplChild1)->nodeType == CONSTANT) && (accessThruMemRef(simplChild2)->nodeType == CONSTANT)) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_pow(*value, *(accessThruMemRef(simplChild1)->value), *(accessThruMemRef(simplChild2)->value), GMP_RNDN);
      free_memory(simplChild1);
      free_memory(simplChild2);
    } else {
      if ((accessThruMemRef(simplChild2)->nodeType == CONSTANT) && (mpfr_cmp_d(*(accessThruMemRef(simplChild2)->value),1.0) == 0) && (!mpfr_nan_p(*(accessThruMemRef(simplChild2)->value)))) {
	freeNode(simplified);
	free_memory(simplChild2);
	simplified = simplChild1;
      } else {
	simplified->nodeType = POW;
	simplified->child1 = simplChild1;
	simplified->child2 = simplChild2;
      }
    }
    break;
  case LIBRARYFUNCTION:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      if (tree->libFun->hasData) {
	sollya_mpfr_from_mpfi_data(*value, *(accessThruMemRef(simplChild1)->value), tree->libFunDeriv, tree->libFun->code, tree->libFun->data);
      } else {
	sollya_mpfr_from_mpfi(*value, *(accessThruMemRef(simplChild1)->value), tree->libFunDeriv, tree->libFun->code);
      }
      free_memory(simplChild1);
    } else {
      simplified->nodeType = LIBRARYFUNCTION;
      simplified->child1 = simplChild1;
      simplified->libFun = tree->libFun;
      simplified->libFunDeriv = tree->libFunDeriv;
    }
    break;
  case PROCEDUREFUNCTION:
    simplChild1 = simplifyAllButDivisionInner(tree->child1);
    simplified = allocateNode();
    if (accessThruMemRef(simplChild1)->nodeType == CONSTANT) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      computeFunctionWithProcedureMpfr(*value, tree->child2, *(accessThruMemRef(simplChild1)->value), (unsigned int) tree->libFunDeriv);
      free_memory(simplChild1);
    } else {
      simplified->nodeType = PROCEDUREFUNCTION;
      simplified->child1 = simplChild1;
      simplified->child2 = copyThing(tree->child2);
      simplified->libFunDeriv = tree->libFunDeriv;
    }
    break;
  case PI_CONST:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    mpfr_const_pi(*value,GMP_RNDN);
    simplified->value = value;
    break;
  case LIBRARYCONSTANT:
    simplified = allocateNode();
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    sollya_mpfi_init2(tempI, tools_precision);
    libraryConstantToInterval(tempI, tree);
    sollya_mpfi_get_fr(*value, tempI);
    sollya_mpfi_clear(tempI);
    simplified->value = value;
    break;
  default:
    sollyaFprintf(stderr,"Error: simplifyAllButDivisionInnerst: unknown identifier (%d) in the tree\n",tree->nodeType);
    exit(1);
  }

  return simplified;
}

node *simplifyTree(node *tree) {
  node *temp, *temp2;

  temp = simplifyTreeErrorfree(tree);
  temp2 = simplifyTreeInner(temp);
  free_memory(temp);

  return temp2;
}

node *simplifyAllButDivision(node *tree) {
  node *temp, *temp2, *temp3, *temp4;

  temp = addMemRef(simplifyTreeErrorfree(tree));
  temp2 = addMemRef(simplifyRationalErrorfree(temp));
  temp3 = addMemRef(simplifyTreeErrorfree(temp2));
  temp4 = addMemRef(simplifyAllButDivisionInner(temp3));
  free_memory(temp);
  free_memory(temp2);
  free_memory(temp3);
  return addMemRef(temp4);
}


void evaluate(mpfr_t result, node *tree, mpfr_t x, mp_prec_t prec) {
  mpfr_t stack1, stack2;
  sollya_mpfi_t stackI, X, Y;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      polynomialEvalMpfr(result, tree->cache->polynomialRepresentation, x);
      return;
    }
    if (tree->cache->evaluationHook != NULL) {
      sollya_mpfi_init2(X,mpfr_get_prec(x));
      sollya_mpfi_set_fr(X, x);
      sollya_mpfi_init2(Y,mpfr_get_prec(result));
      if (evaluateWithEvaluationHook(Y, X, prec, 1, tree->cache->evaluationHook)) {
	if (sollya_mpfi_has_zero(Y)) {
	  mpfr_set_si(result, 0, GMP_RNDN);
	} else {
	  sollya_mpfi_mid(result, Y);
	}
	sollya_mpfi_clear(X);
	sollya_mpfi_clear(Y);
	return;
      }
      sollya_mpfi_clear(X);
      sollya_mpfi_clear(Y);
    }
    evaluate(result, getMemRefChild(tree), x, prec);
    return;
  }

  mpfr_init2(stack1, prec);
  mpfr_init2(stack2, prec);

  switch (tree->nodeType) {
  case VARIABLE:
    mpfr_set(result, x, GMP_RNDN);
    break;
  case CONSTANT:
    mpfr_set(result, *(tree->value), GMP_RNDN);
    break;
  case ADD:
    evaluate(stack1, tree->child1, x, prec);
    evaluate(stack2, tree->child2, x, prec);
    mpfr_add(result, stack1, stack2, GMP_RNDN);
    break;
  case SUB:
    evaluate(stack1, tree->child1, x, prec);
    evaluate(stack2, tree->child2, x, prec);
    mpfr_sub(result, stack1, stack2, GMP_RNDN);
    break;
  case MUL:
    evaluate(stack1, tree->child1, x, prec);
    evaluate(stack2, tree->child2, x, prec);
    mpfr_mul(result, stack1, stack2, GMP_RNDN);
    break;
  case DIV:
    evaluate(stack1, tree->child1, x, prec);
    evaluate(stack2, tree->child2, x, prec);
    mpfr_div(result, stack1, stack2, GMP_RNDN);
    break;
  case NEG:
    evaluate(stack1, tree->child1, x, prec);
    mpfr_neg(result, stack1, GMP_RNDN);
    break;
  case UNARY_BASE_FUNC:
    evaluate(stack1, tree->child1, x, prec);
    tree->baseFun->point_eval(result, stack1);
    break;
  case POW:
    evaluate(stack1, tree->child1, x, prec);
    evaluate(stack2, tree->child2, x, prec);
    mpfr_pow(result, stack1, stack2, GMP_RNDN);
    break;
  case LIBRARYFUNCTION:
    evaluate(stack1, tree->child1, x, prec);
    if (tree->libFun->hasData) {
      sollya_mpfr_from_mpfi_data(result, stack1, tree->libFunDeriv, tree->libFun->code, tree->libFun->data);
    } else {
      sollya_mpfr_from_mpfi(result, stack1, tree->libFunDeriv, tree->libFun->code);
    }
    break;
  case PROCEDUREFUNCTION:
    evaluate(stack1, tree->child1, x, prec);
    computeFunctionWithProcedureMpfr(result, tree->child2, stack1, (unsigned int) tree->libFunDeriv);
    break;
  case PI_CONST:
    mpfr_const_pi(result, GMP_RNDN);
    break;
  case LIBRARYCONSTANT:
    sollya_mpfi_init2(stackI, mpfr_get_prec(result));
    libraryConstantToInterval(stackI, tree);
    sollya_mpfi_get_fr(result, stackI);
    sollya_mpfi_clear(stackI);
    break;
  default:
    sollyaFprintf(stderr,"Error: evaluate: unknown identifier in the tree\n");
    exit(1);
  }

  mpfr_clear(stack1); mpfr_clear(stack2);
  return;
}

int arity(node *tree) {
  switch (tree->nodeType) {
  case MEMREF:
    return arity(getMemRefChild(tree));
    break;
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    return 0;
    break;

  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    return 2;
    break;

  case VARIABLE:
  case NEG:
  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return 1;
    break;

  default:
    sollyaFprintf(stderr,"Error: arity: unknown identifier in the tree\n");
    exit(1);
  }
}


int isSyntacticallyEqualCheap(node *tree1, node *tree2) {

  if (tree1 == tree2) return 1;

  if ((tree1->nodeType == MEMREF) &&
      (tree2->nodeType == MEMREF) &&
      (tree1->cache->polynomialRepresentation != NULL) &&
      (tree2->cache->polynomialRepresentation != NULL) &&
      (tree1->cache->polynomialRepresentation == tree2->cache->polynomialRepresentation)) return 1;

  if (tree1->nodeType == MEMREF) return isSyntacticallyEqualCheap(getMemRefChild(tree1), tree2);
  if (tree2->nodeType == MEMREF) return isSyntacticallyEqualCheap(tree1, getMemRefChild(tree2));

  if (tree1->nodeType != tree2->nodeType) return 0;
  if (tree1->nodeType == PI_CONST) return 1;
  if ((tree1->nodeType == UNARY_BASE_FUNC) &&
      (tree1->baseFun != tree2->baseFun)) return 0;
  if ((tree1->nodeType == LIBRARYFUNCTION) &&
      ((tree1->libFun != tree2->libFun) ||
       (tree1->libFunDeriv != tree2->libFunDeriv))) return 0;
  if (tree1->nodeType == LIBRARYCONSTANT)
    return (tree1->libFun == tree2->libFun);
  if ((tree1->nodeType == PROCEDUREFUNCTION) &&
      ((!isEqualThing(tree1->child2, tree2->child2)) ||
       (tree1->libFunDeriv != tree2->libFunDeriv))) return 0;
  if (tree1->nodeType == CONSTANT) {
    if (mpfr_equal_p(*(tree1->value),*(tree2->value)))
      return 1;
    else
      return 0;
  }
  if (tree1->nodeType == VARIABLE) return 1;

  if (arity(tree1) == 1) {
    if (!isSyntacticallyEqualCheap(tree1->child1,tree2->child1)) return 0;
  } else {
    if (!isSyntacticallyEqualCheap(tree1->child1,tree2->child1)) return 0;
    if (!isSyntacticallyEqualCheap(tree1->child2,tree2->child2)) return 0;
  }

  return 1;
}

int isSyntacticallyEqual(node *tree1, node *tree2) {
  if ((tree1->nodeType == MEMREF) &&
      (tree2->nodeType == MEMREF) &&
      (tree1->cache->polynomialRepresentation != NULL) &&
      (tree2->cache->polynomialRepresentation != NULL)) {
    return polynomialEqual(tree1->cache->polynomialRepresentation,
			   tree2->cache->polynomialRepresentation, 0);
  }

  return isSyntacticallyEqualCheap(tree1, tree2);
}

int isConstant(node *tree);

int isPolynomial(node *tree) {
  int res;
  node *temp;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) return 1;
    return isPolynomial(getMemRefChild(tree));
  }

  if (isConstant(tree)) return 1;

  switch (tree->nodeType) {
  case VARIABLE:
    res = 1;
    break;
  case CONSTANT:
    res = 1;
    break;
  case ADD:
    res = isPolynomial(tree->child1) && isPolynomial(tree->child2);
    break;
  case SUB:
    res = isPolynomial(tree->child1) && isPolynomial(tree->child2);
    break;
  case MUL:
    res = isPolynomial(tree->child1) && isPolynomial(tree->child2);
    break;
  case DIV:
    res = isPolynomial(tree->child1) && isConstant(tree->child2);
    break;
  case NEG:
    res = isPolynomial(tree->child1);
    break;

  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
    res = 0;
    break;

  case POW:
    {
      res = 0;
      if (isPolynomial(tree->child1)) {
	if (accessThruMemRef(tree->child2)->nodeType == CONSTANT)
	  temp = tree->child2;
	else
	  temp = simplifyTreeErrorfree(tree->child2);
	if (accessThruMemRef(temp)->nodeType == CONSTANT) {
	  if (mpfr_integer_p(*(accessThruMemRef(temp)->value))) {
	    if (mpfr_sgn(*(accessThruMemRef(temp)->value)) >= 0) {
	      res = 1;
	    }
	  }
	}
	if (accessThruMemRef(tree->child2)->nodeType != CONSTANT) free_memory(temp);
      }
    }
    break;

  case PROCEDUREFUNCTION:
    res = 0;
    break;
  case PI_CONST:
  case LIBRARYCONSTANT:
    res = 1;
    break;
  default:
    sollyaFprintf(stderr,"Error: isPolynomial: unknown identifier in the tree\n");
    exit(1);
  }
  return res;
}

int isRationalFunction(node *tree) {
  node *temp;
  int res;

  /* Handle all polynomials and constant expressions (constants, pi,
     library constants etc.)
  */
  if (isPolynomial(tree)) return 1;

  /* Handle memory reference counting */
  if (tree->nodeType == MEMREF) {
    return isRationalFunction(getMemRefChild(tree));
  }

  switch (tree->nodeType) {
  case ADD:
  case SUB:
  case MUL:
  case DIV:
    /* Sums, productions and quotients are rational functions iff
       their operands are rational functions.
    */
    return (isRationalFunction(tree->child1) &&
	    isRationalFunction(tree->child2));
    break;
  case NEG:
    /* -f is a rational function iff f is a rational function */
    return isRationalFunction(tree->child1);
    break;
  case POW:
    /* f^k is a rational function if f is a rational function and k is
       a constant integer
    */
    if (!isRationalFunction(tree->child1)) return 0;
    res = 0;
    if (isPolynomial(tree->child1)) {
      if (accessThruMemRef(tree->child2)->nodeType == CONSTANT)
	temp = tree->child2;
      else
	temp = simplifyTreeErrorfree(tree->child2);
      if (accessThruMemRef(temp)->nodeType == CONSTANT) {
	if (mpfr_integer_p(*(accessThruMemRef(temp)->value))) {
	  res = 1;
	}
      }
      if (accessThruMemRef(tree->child2)->nodeType != CONSTANT) free_memory(temp);
    }
    return res;
    break;
  default:
    break;
  }

  /* All other functions are deemed not to be rational functions */
  return 0;
}


int isAffine(node *tree) {
  int res;
  node *temp;

  if (tree->nodeType == MEMREF) {
    return isAffine(getMemRefChild(tree));
  }

  switch (tree->nodeType) {
  case VARIABLE:
    res = 1;
    break;
  case CONSTANT:
    res = 1;
    break;
  case ADD:
    res = isAffine(tree->child1) && isAffine(tree->child2);
    break;
  case SUB:
    res = isAffine(tree->child1) && isAffine(tree->child2);
    break;
  case MUL:
    res = isAffine(tree->child1) && isAffine(tree->child2);
    break;
  case NEG:
    res = isAffine(tree->child1);
    break;

  case DIV:
  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    res = 0;
    break;

  case POW:
    {
      res = 0;
      if (isAffine(tree->child1)) {
        if (accessThruMemRef(tree->child2)->nodeType == CONSTANT)
          temp = tree->child2;
        else
          temp = simplifyTreeErrorfree(tree->child2);
        if (accessThruMemRef(temp)->nodeType == CONSTANT) {
          if (mpfr_number_p(*(accessThruMemRef(temp)->value)) && mpfr_integer_p(*(accessThruMemRef(temp)->value))) {
            if (mpfr_sgn(*(accessThruMemRef(temp)->value)) > 0) {
              res = 1;
            }
          }
        }
        if (accessThruMemRef(tree->child2)->nodeType != CONSTANT) free_memory(temp);
      }
    }
    break;

  case PI_CONST:
  case LIBRARYCONSTANT:
    res = 1;
    break;
  default:
    sollyaFprintf(stderr,"Error: isAffine: unknown identifier in the tree\n");
    exit(1);
  }
  return res;
}


#define MAX_MACRO(a,b) (a) > (b) ? (a) : (b)
#define MIN_MACRO(a,b) (a) < (b) ? (a) : (b)

int getDegreeUnsafe(node *tree, int silent) {
  int l, r;
  mpfr_t temp;
  node *simplifiedExponent;
  unsigned int h;

  if (isConstant(tree)) return 0;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      return polynomialGetDegreeAsInt(tree->cache->polynomialRepresentation);
    }
    return getDegreeUnsafe(getMemRefChild(tree), silent);
  }

  switch (tree->nodeType) {
  case VARIABLE:
    return 1;
    break;
  case CONSTANT:
  case PI_CONST:
    return 0;
    break;
  case ADD:
    l = getDegreeUnsafe(tree->child1, silent);
    r = getDegreeUnsafe(tree->child2, silent);
    if ((l >= 0) && (r >= 0)) return MAX_MACRO(l,r); else return -1;
    break;
  case SUB:
    l = getDegreeUnsafe(tree->child1, silent);
    r = getDegreeUnsafe(tree->child2, silent);
    if ((l >= 0) && (r >= 0)) return MAX_MACRO(l,r); else return -1;
    break;
  case MUL:
    l = getDegreeUnsafe(tree->child1, silent);
    r = getDegreeUnsafe(tree->child2, silent);
    if ((l >= 0) && (r >= 0)) return l + r; else return -1;
    break;
  case DIV:
    return getDegreeUnsafe(tree->child1, silent);
    break;
  case POW:
    {
      l = getDegreeUnsafe(tree->child1, silent);
      if (accessThruMemRef(tree->child2)->nodeType != CONSTANT) {
        simplifiedExponent = simplifyRationalErrorfree(tree->child2);
        if ((accessThruMemRef(simplifiedExponent)->nodeType == CONSTANT) &&
            mpfr_integer_p(*(accessThruMemRef(simplifiedExponent)->value)) &&
            (mpfr_sgn(*(accessThruMemRef(simplifiedExponent)->value)) >= 0)) {
          h = mpfr_get_ui(*(accessThruMemRef(simplifiedExponent)->value),GMP_RNDN);
          mpfr_init2(temp,mpfr_get_prec(*(accessThruMemRef(simplifiedExponent)->value)) + 10);
          mpfr_set_ui(temp,h,GMP_RNDN);
          if (mpfr_cmp(*(accessThruMemRef(simplifiedExponent)->value),temp) != 0) {
            if (!silent) printMessage(1, SOLLYA_MSG_DEGREE_OF_POLYNOMIAL_DOESNT_HOLD_ON_MACHINE_INT, "Warning: tried to compute polynomial degree of an expression using a power operator with an exponent which cannot be represented on an integer variable.\n");
            mpfr_clear(temp);
            free_memory(simplifiedExponent);
            return -1;
          }
          mpfr_clear(temp);
          free_memory(simplifiedExponent);
	  r = (int) h;
          if ((l >= 0) && (r >= 0)) return l * r; else return -1;
        } else {
          sollyaFprintf(stderr,"Error: getDegreeUnsafe: an error occurred. The exponent in a power operator is not constant, not integer or not non-negative.\n");
          exit(1);
        }
        free_memory(simplifiedExponent);
      } else {
        if (!mpfr_integer_p(*(accessThruMemRef(tree->child2)->value))) {
          sollyaFprintf(stderr,"Error: getDegreeUnsafe: an error occurred. The exponent in a power operator is not integer.\n");
          exit(1);
        }
        if (mpfr_sgn(*(accessThruMemRef(tree->child2)->value)) < 0) {
          sollyaFprintf(stderr,"Error: getDegreeUnsafe: an error occurred. The exponent in a power operator is negative.\n");
          exit(1);
        }

        h = mpfr_get_ui(*(accessThruMemRef(tree->child2)->value),GMP_RNDN);
        mpfr_init2(temp,mpfr_get_prec(*(accessThruMemRef(tree->child2)->value)) + 10);
        mpfr_set_ui(temp,h,GMP_RNDN);
        if (mpfr_cmp(*(accessThruMemRef(tree->child2)->value),temp) != 0) {
          if (!silent) printMessage(1, SOLLYA_MSG_DEGREE_OF_POLYNOMIAL_DOESNT_HOLD_ON_MACHINE_INT, "Warning: tried to compute polynomial degree of an expression using a power operator with an exponent which cannot be represented on an integer variable.\n");
          mpfr_clear(temp);
          return -1;
        }
        mpfr_clear(temp);
	r = (int) h;
        if ((l >= 0) && (r >= 0)) return l * r; else return -1;
      }
    }
    break;
  case NEG:
    return getDegreeUnsafe(tree->child1, silent);
    break;
  default:
    sollyaFprintf(stderr,"Error: getDegreeUnsafe: an error occurred on handling the expression tree\n");
    exit(1);
  }
}

int getDegreeUnsafeMpz(mpz_t res, node *tree) {
  mpz_t left, right;
  int l, r;
  node *simplifiedExponent;

  if (isConstant(tree)) {
    mpz_set_si(res, 0);
    return 1;
  }

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      polynomialGetDegree(res, tree->cache->polynomialRepresentation);
      return 1;
    }
    return getDegreeUnsafeMpz(res, getMemRefChild(tree));
  }

  switch (tree->nodeType) {
  case VARIABLE:
    mpz_set_si(res, 1);
    return 1;
    break;
  case CONSTANT:
  case PI_CONST:
    mpz_set_si(res, 1);
    return 1;
    break;
  case ADD:
    mpz_init(left);
    mpz_init(right);
    l = getDegreeUnsafeMpz(left, tree->child1);
    r = getDegreeUnsafeMpz(right, tree->child2);
    if (l && r) {
      if (mpz_cmp(left, right) < 0) {
	mpz_set(res, right);
      } else {
	mpz_set(res, left);
      }
      mpz_clear(left);
      mpz_clear(right);
      return 1;
    }
    mpz_clear(left);
    mpz_clear(right);
    break;
  case SUB:
    mpz_init(left);
    mpz_init(right);
    l = getDegreeUnsafeMpz(left, tree->child1);
    r = getDegreeUnsafeMpz(right, tree->child2);
    if (l && r) {
      if (mpz_cmp(left, right) < 0) {
	mpz_set(res, right);
      } else {
	mpz_set(res, left);
      }
      mpz_clear(left);
      mpz_clear(right);
      return 1;
    }
    mpz_clear(left);
    mpz_clear(right);
    break;
  case MUL:
    mpz_init(left);
    mpz_init(right);
    l = getDegreeUnsafeMpz(left, tree->child1);
    r = getDegreeUnsafeMpz(right, tree->child2);
    if (l && r) {
      mpz_add(res, left, right);
      mpz_clear(left);
      mpz_clear(right);
      return 1;
    }
    mpz_clear(left);
    mpz_clear(right);
    break;
  case DIV:
    return getDegreeUnsafeMpz(res, tree->child1);
    break;
  case POW:
    {
      mpz_init(left);
      l = getDegreeUnsafeMpz(left, tree->child1);
      if (accessThruMemRef(tree->child2)->nodeType != CONSTANT) {
        simplifiedExponent = simplifyRationalErrorfree(tree->child2);
        if ((accessThruMemRef(simplifiedExponent)->nodeType == CONSTANT) &&
            mpfr_integer_p(*(accessThruMemRef(simplifiedExponent)->value)) &&
            (mpfr_sgn(*(accessThruMemRef(simplifiedExponent)->value)) >= 0)) {
	  mpz_init(right);
	  mpfr_get_z(right, *(accessThruMemRef(simplifiedExponent)->value), GMP_RNDN); /* exact */
	  if (l) {
	    mpz_mul(res, left, right);
	    mpz_clear(left);
	    mpz_clear(right);
	    free_memory(simplifiedExponent);
	    return 1;
	  }
	  mpz_init(right);
        } else {
          sollyaFprintf(stderr,"Error: getDegreeUnsafeMpz: an error occurred. The exponent in a power operator is not constant, not integer or not non-negative.\n");
          exit(1);
        }
        free_memory(simplifiedExponent);
      } else {
        if (!mpfr_integer_p(*(accessThruMemRef(tree->child2)->value))) {
          sollyaFprintf(stderr,"Error: getDegreeUnsafe: an error occurred. The exponent in a power operator is not integer.\n");
          exit(1);
        }
        if (mpfr_sgn(*(accessThruMemRef(tree->child2)->value)) < 0) {
          sollyaFprintf(stderr,"Error: getDegreeUnsafe: an error occurred. The exponent in a power operator is negative.\n");
          exit(1);
        }
	mpz_init(right);
	mpfr_get_z(right, *(accessThruMemRef(tree->child2)->value), GMP_RNDN); /* exact */
	if (l) {
	  mpz_mul(res, left, right);
	  mpz_clear(left);
	  mpz_clear(right);
	  return 1;
	}
	mpz_init(right);
      }
      mpz_init(left);
    }
    break;
  case NEG:
    return getDegreeUnsafeMpz(res, tree->child1);
    break;
  default:
    sollyaFprintf(stderr,"Error: getDegreeUnsafeMpz: an error occurred on handling the expression tree\n");
    exit(1);
  }

  return 0;
}

int getDegree(node *tree) {
  if (!isPolynomial(tree)) return -1;
  return getDegreeUnsafe(tree, 0);
}

int getDegreeSilent(node *tree) {
  if (!isPolynomial(tree)) return -1;
  return getDegreeUnsafe(tree, 1);
}

int getDegreeMpz(mpz_t res, node *tree) {
  if (!isPolynomial(tree)) {
    mpz_set_si(res, -1);
    return 0;
  }
  return getDegreeUnsafeMpz(res, tree);
}

int getDegreeMpzVerified(mpz_t res, node *tree) {
  int okay, k, gottaBreak;
  node *tempNode;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation == NULL) {
      tryRepresentAsPolynomial(tree);
    }
    if (tree->cache->polynomialRepresentation != NULL) {
      polynomialGetDegree(res, tree->cache->polynomialRepresentation);
      if (mpz_cmp_si(res, -1) == 0) {
	printMessage(1,SOLLYA_MSG_DEGREE_OF_POLYNOMIAL_LARGER_THAN_MULTIPRECISION_INT,
		     "Warning: the degree of the given polynomial is larger than the largest multiprecision integer that can be held in memory. The polynomial's degree will be returned as -1.\n");
      }
      return 1;
    }
  }


  okay = getDegreeMpz(res, tree);

  if (okay) {
    k = mpz_get_si(res);
    if ((mpz_cmp_si(res, k) == 0) && (k > 0)) {
      while (k > 0) {
	tempNode = getIthCoefficient(tree, k);
	gottaBreak = 1;
	if (accessThruMemRef(tempNode)->nodeType == CONSTANT) {
	  if (mpfr_zero_p(*(accessThruMemRef(tempNode)->value))) {
	    gottaBreak = 0;
	  }
	}
	free_memory(tempNode);
	if (gottaBreak) break;
	k--;
	mpz_set_si(res, k);
      }
    }
  }

  return okay;
}

int isPolynomialExtraSafe(node *tree) {
  return (isPolynomial(tree) && (getDegreeSilent(tree) >= 0));
}

int getMaxPowerDividerUnsafe(node *tree) {
  int l, r;
  mpfr_t temp;
  node *simplifiedNode;

  if (isConstant(tree)) return 0;

  if (tree->nodeType == MEMREF) return getMaxPowerDividerUnsafe(getMemRefChild(tree));

  switch (tree->nodeType) {
  case VARIABLE:
    return 1;
    break;
  case CONSTANT:
  case PI_CONST:
    return 0;
    break;
  case ADD:
    l = getMaxPowerDividerUnsafe(tree->child1);
    r = getMaxPowerDividerUnsafe(tree->child2);
    return MIN_MACRO(l,r);
    break;
  case SUB:
    l = getMaxPowerDividerUnsafe(tree->child1);
    r = getMaxPowerDividerUnsafe(tree->child2);
    return MIN_MACRO(l,r);
    break;
  case MUL:
    l = getMaxPowerDividerUnsafe(tree->child1);
    r = getMaxPowerDividerUnsafe(tree->child2);
    return l + r;
    break;
  case DIV:
    return getMaxPowerDividerUnsafe(tree->child1);
    break;
  case POW:
    {
      l = getMaxPowerDividerUnsafe(tree->child1);
      if (l == 0) return 0;
      simplifiedNode = simplifyRationalErrorfree(tree->child2);
      if (accessThruMemRef(simplifiedNode)->nodeType != CONSTANT) {
	printMessage(1,SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_IS_NOT_CONSTANT,
		     "Warning: an attempt was made to compute the degree of the maximal polynomial divider of a polynomial in an expression using a power operator with an exponent which is not a constant but a constant expression.\n");
	free_memory(simplifiedNode);
	return -1;
      }
      if (!mpfr_integer_p(*(accessThruMemRef(simplifiedNode)->value))) {
	printMessage(1,SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_IS_NOT_INTEGER,
		     "Warning: an attempt was made to compute the degree of the maximal polynomial divider of a polynomial in an expression using a power operator with an exponent which is not an integer.\n");
	free_memory(simplifiedNode);
	return -1;
      }
      if (mpfr_sgn(*(accessThruMemRef(simplifiedNode)->value)) < 0) {
	printMessage(1,SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_IS_NEGATIVE,
		     "Warning: an attempt was made to compute the degree of the maximal polynomial divider of a polynomial in an expression using a power operator with an exponent which is negative.\n");
	free_memory(simplifiedNode);
	return -1;
      }

      r = mpfr_get_si(*(accessThruMemRef(simplifiedNode)->value),GMP_RNDN);
      mpfr_init2(temp,mpfr_get_prec(*(accessThruMemRef(simplifiedNode)->value)) + 10);
      mpfr_set_si(temp,r,GMP_RNDN);
      if (mpfr_cmp(*(accessThruMemRef(simplifiedNode)->value),temp) != 0) {
	printMessage(1,SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_DOESNT_HOLD_ON_MACHINE_INT,
                     "Warning: tried to compute degree of maximal polynomial divider of a polynomial in an expression using a power operator with an exponent which cannot be represented on an integer variable.\n");
	mpfr_clear(temp);
	free_memory(simplifiedNode);
	return -1;
      }
      mpfr_clear(temp);
      free_memory(simplifiedNode);
      return l * r;
    }
    break;
  case NEG:
    return getMaxPowerDividerUnsafe(tree->child1);
    break;
  default:
    sollyaFprintf(stderr,"Error: getMaxPowerDividerUnsafe: an error occurred on handling the expression tree\n");
    exit(1);
  }
}

int getMaxPowerDivider(node *tree) {
  if (!isPolynomial(tree)) return -1;
  return getMaxPowerDividerUnsafe(tree);
}

node* makeBinomial(node *a, node *b, int n, int s) {
  node *tree, *coeff, *aPow, *bPow, *tempNode, *tempNode2;
  mpfr_t *coeffVal, *mpfr_temp;
  unsigned int i;
  mpz_t coeffGMP;
  mp_prec_t prec;

  tree = allocateNode();
  tree->nodeType = CONSTANT;
  mpfr_temp = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*mpfr_temp,tools_precision);
  mpfr_set_d(*mpfr_temp,0.0,GMP_RNDN);
  tree->value = mpfr_temp;
  mpz_init(coeffGMP);
  for (i=0;i<=((unsigned int) n);i++) {
    mpz_bin_uiui(coeffGMP,(unsigned int) n,i);
    prec = mpz_sizeinbase (coeffGMP, 2) + 10;
    if (prec < tools_precision) prec = tools_precision;
    coeffVal = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*coeffVal,prec);
    if(mpfr_set_z(*coeffVal,coeffGMP,GMP_RNDN) != 0) {
      if (!noRoundingWarnings) {
	printMessage(1,SOLLYA_MSG_ROUNDING_UPON_BINOMIAL_COEFFICIENT_COMPUTATION,"Warning: on expanding a power operator a rounding occurred when calculating a binomial coefficient.\n");
	printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the working precision.\n");
      }
    }
    if ((s < 0) && (((((unsigned int) n) - i) & 1) != 0)) { /* This is a modulo 2 to determine eveness */
      mpfr_neg(*coeffVal,*coeffVal,GMP_RNDN);
    }
    coeff = allocateNode();
    coeff->nodeType = CONSTANT;
    coeff->value = coeffVal;
    aPow = allocateNode();
    aPow->nodeType = POW;
    aPow->child1 = copyTree(a);
    tempNode = allocateNode();
    tempNode->nodeType = CONSTANT;
    mpfr_temp = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*mpfr_temp,tools_precision);
    if(mpfr_set_ui(*mpfr_temp,i,GMP_RNDN) != 0) {
      if (!noRoundingWarnings) {
	printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: on expanding a power operator a rounding occurred when calculating an exponent constant.\n");
	printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the working precision.\n");
      }
    }
    tempNode->value = mpfr_temp;
    aPow->child2 = tempNode;
    bPow = allocateNode();
    bPow->nodeType = POW;
    bPow->child1 = copyTree(b);
    tempNode = allocateNode();
    tempNode->nodeType = CONSTANT;
    mpfr_temp = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*mpfr_temp,tools_precision);
    if(mpfr_set_ui(*mpfr_temp,((unsigned int) n) - i,GMP_RNDN) != 0) {
      printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: on expanding a power operator a rounding occurred when calculating an exponent constant.\n");
      printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the working precision.\n");
    }
    tempNode->value = mpfr_temp;
    bPow->child2 = tempNode;
    tempNode = allocateNode();
    tempNode->nodeType = MUL;
    tempNode->child1 = coeff;
    tempNode->child2 = aPow;
    tempNode2 = allocateNode();
    tempNode2->nodeType = MUL;
    tempNode2->child1 = tempNode;
    tempNode2->child2 = bPow;
    tempNode = allocateNode();
    tempNode->nodeType = ADD;
    tempNode->child1 = tree;
    tempNode->child2 = tempNode2;
    tree = tempNode;
  }
  mpz_clear(coeffGMP);

  return tree;
}


node* expandPowerInPolynomialUnsafe(node *tree) {
  node *copy, *left, *tempTree, *tempTree2, *tempTree3;
  int r, i;
  mpfr_t temp;
  mpfr_t *value;

  if (isConstant(tree)) return copyTree(tree);

  if (tree->nodeType == MEMREF) return addMemRef(expandPowerInPolynomialUnsafe(getMemRefChild(tree)));

  switch (tree->nodeType) {
  case VARIABLE:
    return copyTree(tree);
    break;
  case CONSTANT:
  case PI_CONST:
    return copyTree(tree);
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = expandPowerInPolynomialUnsafe(tree->child1);
    copy->child2 = expandPowerInPolynomialUnsafe(tree->child2);
    return copy;
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = expandPowerInPolynomialUnsafe(tree->child1);
    copy->child2 = expandPowerInPolynomialUnsafe(tree->child2);
    return copy;
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = expandPowerInPolynomialUnsafe(tree->child1);
    copy->child2 = expandPowerInPolynomialUnsafe(tree->child2);
    return copy;
    break;
  case DIV:
    copy = allocateNode();
    copy->nodeType = DIV;
    copy->child1 = expandPowerInPolynomialUnsafe(tree->child1);
    copy->child2 = expandPowerInPolynomialUnsafe(tree->child2);
    return copy;
    break;
  case POW:
    {
      left = expandPowerInPolynomialUnsafe(tree->child1);
      if (accessThruMemRef(tree->child2)->nodeType != CONSTANT) {
	sollyaFprintf(stderr,"Error: expandPowerInPolynomialUnsafe: an error occurred. The exponent in a power operator is not constant.\n");
	exit(1);
      }
      if (!mpfr_integer_p(*(accessThruMemRef(tree->child2)->value))) {
	sollyaFprintf(stderr,"Error: expandPowerInPolynomialUnsafe: an error occurred. The exponent in a power operator is not integer.\n");
	exit(1);
      }
      if (mpfr_sgn(*(accessThruMemRef(tree->child2)->value)) < 0) {
	sollyaFprintf(stderr,"Error: expandPowerInPolynomialUnsafe: an error occurred. The exponent in a power operator is negative.\n");
	exit(1);
      }

      r = mpfr_get_si(*(accessThruMemRef(tree->child2)->value),GMP_RNDN);
      mpfr_init2(temp,mpfr_get_prec(*(accessThruMemRef(tree->child2)->value)) + 10);
      mpfr_set_si(temp,r,GMP_RNDN);
      if (mpfr_cmp(*(accessThruMemRef(tree->child2)->value),temp) != 0) {
	sollyaFprintf(stderr,"Error: expandPowerInPolynomialUnsafe: an error occurred. Tried to expand an expression using a power operator with an exponent ");
	sollyaFprintf(stderr,"which cannot be represented on an integer variable.\n");
	mpfr_clear(temp);
	exit(1);
      }
      mpfr_clear(temp);
      if (r > 1) {
	switch (accessThruMemRef(left)->nodeType) {
	case VARIABLE:
	case CONSTANT:
	case PI_CONST:
	  tempTree = copyTree(left);
	  for (i=1;i<r;i++) {
	    tempTree2 = allocateNode();
	    tempTree2->nodeType = MUL;
	    tempTree2->child1 = tempTree;
	    tempTree2->child2 = copyTree(left);
	    tempTree = tempTree2;
	  }
	  break;
	case ADD:
	  tempTree = makeBinomial(accessThruMemRef(left)->child1,accessThruMemRef(left)->child2,r,1);
	  break;
	case SUB:
	  tempTree = makeBinomial(accessThruMemRef(left)->child1,accessThruMemRef(left)->child2,r,-1);
	  break;
	case MUL:
	  tempTree = allocateNode();
	  tempTree->nodeType = MUL;
	  tempTree2 = allocateNode();
	  tempTree2->nodeType = POW;
	  tempTree2->child1 = copyTree(accessThruMemRef(left)->child1);
	  tempTree3 = allocateNode();
	  tempTree3->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,tools_precision);
	  mpfr_set_si(*value,r,GMP_RNDN);
	  tempTree3->value = value;
	  tempTree2->child2 = tempTree3;
	  tempTree->child1 = tempTree2;
	  tempTree2 = allocateNode();
	  tempTree2->nodeType = POW;
	  tempTree2->child1 = copyTree(accessThruMemRef(left)->child2);
	  tempTree3 = allocateNode();
	  tempTree3->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,tools_precision);
	  mpfr_set_si(*value,r,GMP_RNDN);
	  tempTree3->value = value;
	  tempTree2->child2 = tempTree3;
	  tempTree->child2 = tempTree2;
	  break;
	case DIV:
	  tempTree = allocateNode();
	  tempTree->nodeType = DIV;
	  tempTree2 = allocateNode();
	  tempTree2->nodeType = POW;
	  tempTree2->child1 = copyTree(accessThruMemRef(left)->child1);
	  tempTree3 = allocateNode();
	  tempTree3->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,tools_precision);
	  mpfr_set_si(*value,r,GMP_RNDN);
	  tempTree3->value = value;
	  tempTree2->child2 = tempTree3;
	  tempTree->child1 = tempTree2;
	  tempTree2 = allocateNode();
	  tempTree2->nodeType = POW;
	  tempTree2->child1 = copyTree(accessThruMemRef(left)->child2);
	  tempTree3 = allocateNode();
	  tempTree3->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,tools_precision);
	  mpfr_set_si(*value,r,GMP_RNDN);
	  tempTree3->value = value;
	  tempTree2->child2 = tempTree3;
	  tempTree->child2 = tempTree2;
	  break;
	case NEG:
	  if (r & 1) { /* This is a modulo 2 to determine eveness, not a logical test */
	    /* r is odd */
	    tempTree = allocateNode();
	    tempTree->nodeType = NEG;
	    tempTree2 = allocateNode();
	    tempTree2->nodeType = POW;
	    tempTree2->child1 = copyTree(accessThruMemRef(left)->child1);
	    tempTree3 = allocateNode();
	    tempTree3->nodeType = CONSTANT;
	    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	    mpfr_init2(*value,tools_precision);
	    mpfr_set_si(*value,r,GMP_RNDN);
	    tempTree3->value = value;
	    tempTree2->child2 = tempTree3;
	    tempTree->child1 = tempTree2;
	  } else {
	    tempTree = allocateNode();
	    tempTree->nodeType = POW;
	    tempTree->child1 = copyTree(accessThruMemRef(left)->child1);
	    tempTree3 = allocateNode();
	    tempTree3->nodeType = CONSTANT;
	    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	    mpfr_init2(*value,tools_precision);
	    mpfr_set_si(*value,r,GMP_RNDN);
	    tempTree3->value = value;
	    tempTree->child2 = tempTree3;
	  }
	  break;
	default:
	  if (isConstant(left)) return copyTree(tree);

	  sollyaFprintf(stderr,"Error: expandPowerInPolynomialUnsafe: an error occurred on handling the expanded expression subtree\n");
	  exit(1);
	}
	copy = expandPowerInPolynomialUnsafe(tempTree);
	free_memory(tempTree);
      } else {
	if (r == 1) {
	  copy = copyTree(left);
	} else {
	  copy = allocateNode();
	  copy->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,tools_precision);
	  mpfr_set_d(*value,1.0,GMP_RNDN);
	  copy->value = value;
	}
      }
      free_memory(left);
      return copy;
    }
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = expandPowerInPolynomialUnsafe(tree->child1);
    return copy;
    break;
  default:

    if (isConstant(tree)) return copyTree(tree);

    sollyaFprintf(stderr,"Error: expandPowerInPolynomialUnsafe: an error occurred on handling the expression tree\n");
    exit(1);
  }
}

node* expandPowerInPolynomial(node *tree) {
  if (getDegree(tree) < 0) return copyTree(tree);
  return expandPowerInPolynomialUnsafe(tree);
}



node* expandDivision(node *tree) {
  node *copy, *left, *right, *tempNode;
  mpfr_t *value;
  mpfr_t temp;

  if (tree->nodeType == MEMREF) {
    return addMemRef(expandDivision(getMemRefChild(tree)));
  }

  switch (tree->nodeType) {
  case VARIABLE:
    copy = makeVariable();
    break;
  case CONSTANT:
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,tools_precision);
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    copy->value = value;
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = expandDivision(tree->child1);
    copy->child2 = expandDivision(tree->child2);
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = expandDivision(tree->child1);
    copy->child2 = expandDivision(tree->child2);
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = expandDivision(tree->child1);
    copy->child2 = expandDivision(tree->child2);
    break;
  case DIV:
    left = expandDivision(tree->child1);
    right = expandDivision(tree->child2);
    if (accessThruMemRef(right)->nodeType == DIV) {
      tempNode = allocateNode();
      tempNode->nodeType = DIV;
      tempNode->child1 = copyTree(accessThruMemRef(right)->child2);
      tempNode->child2 = copyTree(accessThruMemRef(right)->child1);
      free_memory(right);
      copy = allocateNode();
      copy->nodeType = MUL;
      copy->child1 = left;
      copy->child2 = tempNode;
    } else {
      copy = allocateNode();
      copy->nodeType = DIV;
      copy->child1 = left;
      copy->child2 = right;
    }
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = expandDivision(tree->child1);
    break;
  case UNARY_BASE_FUNC:
    copy = allocateNode();
    copy->nodeType = UNARY_BASE_FUNC;
    copy->baseFun = tree->baseFun;
    copy->child1 = expandDivision(tree->child1);
    break;
  case POW:
    copy = allocateNode();
    copy->nodeType = POW;
    copy->child1 = expandDivision(tree->child1);
    copy->child2 = expandDivision(tree->child2);
    break;
  case LIBRARYFUNCTION:
    copy = allocateNode();
    copy->nodeType = LIBRARYFUNCTION;
    copy->libFun = tree->libFun;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = expandDivision(tree->child1);
    break;
  case PROCEDUREFUNCTION:
    copy = allocateNode();
    copy->nodeType = PROCEDUREFUNCTION;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = expandDivision(tree->child1);
    copy->child2 = copyThing(tree->child2);
    break;
  case PI_CONST:
    copy = allocateNode();
    copy->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    copy = allocateNode();
    copy->nodeType = LIBRARYCONSTANT;
    copy->libFun = tree->libFun;
    break;

  default:
    sollyaFprintf(stderr,"Error: expandDivision: unknown identifier in the tree\n");
    exit(1);
  }
  return copy;
}



node* expandPolynomialUnsafe(node *tree) {
  node *left, *right, *copy, *tempNode, *tempNode2, *tempNode3, *tempNode4;
  mpfr_t *value;

  if (isConstant(tree)) return copyTree(tree);

  if (tree->nodeType == MEMREF) return addMemRef(expandPolynomialUnsafe(getMemRefChild(tree)));

  switch (tree->nodeType) {
  case VARIABLE:
    return copyTree(tree);
    break;
  case CONSTANT:
  case PI_CONST:
    return copyTree(tree);
    break;
  case ADD:
    left = expandPolynomialUnsafe(tree->child1);
    right = expandPolynomialUnsafe(tree->child2);
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = left;
    copy->child2 = right;
    return copy;
    break;
  case SUB:
    left = expandPolynomialUnsafe(tree->child1);
    right = expandPolynomialUnsafe(tree->child2);
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = left;
    copy->child2 = right;
    return copy;
    break;
  case MUL:
    left = expandPolynomialUnsafe(tree->child1);
    right = expandPolynomialUnsafe(tree->child2);
    switch (accessThruMemRef(left)->nodeType) {
    case VARIABLE:
    case CONSTANT:
    case PI_CONST:
      if (isConstant(right)) {
	copy = allocateNode();
	copy->nodeType = MUL;
	copy->child1 = left;
	copy->child2 = right;
      } else {
	switch (accessThruMemRef(right)->nodeType) {
	case VARIABLE:
	case CONSTANT:
	case PI_CONST:
	  copy = allocateNode();
	  copy->nodeType = MUL;
	  copy->child1 = left;
	  copy->child2 = right;
	  break;
	default:
	  tempNode = allocateNode();
	  tempNode->nodeType = MUL;
	  tempNode->child1 = right;
	  tempNode->child2 = left;
	  copy = expandPolynomialUnsafe(tempNode);
	  free_memory(tempNode);
	}
      }
      break;
    case MUL:
      switch (accessThruMemRef(right)->nodeType) {
      case ADD:
      case SUB:
      case NEG:
      case DIV:
	tempNode = allocateNode();
	tempNode->nodeType = MUL;
	tempNode->child1 = right;
	tempNode->child2 = left;
	copy = expandPolynomialUnsafe(tempNode);
	free_memory(tempNode);
	break;
      default:
	copy = allocateNode();
	copy->nodeType = MUL;
	copy->child1 = left;
	copy->child2 = right;
	break;
      }
      break;
    case ADD:
      tempNode = allocateNode();
      tempNode->nodeType = ADD;
      tempNode2 = allocateNode();
      tempNode2->nodeType = MUL;
      tempNode3 = allocateNode();
      tempNode3->nodeType = MUL;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2->child2 = copyTree(right);
      tempNode3->child1 = copyTree(accessThruMemRef(left)->child2);
      tempNode3->child2 = right;
      tempNode->child1 = tempNode2;
      tempNode->child2 = tempNode3;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case SUB:
      tempNode = allocateNode();
      tempNode->nodeType = SUB;
      tempNode2 = allocateNode();
      tempNode2->nodeType = MUL;
      tempNode3 = allocateNode();
      tempNode3->nodeType = MUL;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2->child2 = copyTree(right);
      tempNode3->child1 = copyTree(accessThruMemRef(left)->child2);
      tempNode3->child2 = right;
      tempNode->child1 = tempNode2;
      tempNode->child2 = tempNode3;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case NEG:
      tempNode = allocateNode();
      tempNode->nodeType = NEG;
      tempNode2 = allocateNode();
      tempNode2->nodeType = MUL;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2->child2 = right;
      tempNode->child1 = tempNode2;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case DIV:
      if (isConstant(left)) {
	if (isConstant(right)) {
	  copy = allocateNode();
	  copy->nodeType = MUL;
	  copy->child1 = left;
	  copy->child2 = right;
	} else {
	  switch (accessThruMemRef(right)->nodeType) {
	  case ADD:
	  case SUB:
	  case NEG:
	  case DIV:
	    tempNode = allocateNode();
	    tempNode->nodeType = MUL;
	    tempNode->child1 = right;
	    tempNode->child2 = left;
	    copy = expandPolynomialUnsafe(tempNode);
	    free_memory(tempNode);
	    break;
	  default:
	    copy = allocateNode();
	    copy->nodeType = MUL;
	    copy->child1 = left;
	    copy->child2 = right;
	    break;
	  }
	}
      } else {
	tempNode = allocateNode();
	tempNode->nodeType = MUL;
	tempNode2 = allocateNode();
	tempNode2->nodeType = MUL;
	tempNode3 = allocateNode();
	tempNode3->nodeType = DIV;
	tempNode4 = allocateNode();
	tempNode4->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	mpfr_set_d(*value,1.0,GMP_RNDN);
	tempNode4->value = value;
	tempNode3->child1 = tempNode4;
	tempNode3->child2 = copyTree(accessThruMemRef(left)->child2);
	tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
	tempNode2->child2 = right;
	tempNode->child1 = tempNode3;
	tempNode->child2 = tempNode2;
	free_memory(left);
	copy = expandPolynomialUnsafe(tempNode);
	free_memory(tempNode);
      }
      break;
    default:
      if (isConstant(left)) {
	if (isConstant(right)) {
	  return copyTree(tree);
	} else {
	  tempNode = allocateNode();
	  tempNode->nodeType = MUL;
	  tempNode->child1 = right;
	  tempNode->child2 = left;
	  copy = expandPolynomialUnsafe(tempNode);
	  free_memory(tempNode);
	  return copy;
	}
      } else {
	sollyaFprintf(stderr,"Error: expandPolynomialUnsafe: an error occurred on handling the MUL left rewritten expression subtree\n");
	exit(1);
      }
    }
    return copy;
    break;
  case DIV:
    left = expandPolynomialUnsafe(tree->child1);
    right = expandPolynomialUnsafe(tree->child2);
    switch (accessThruMemRef(left)->nodeType) {
    case CONSTANT:
    case PI_CONST:
      copy = allocateNode();
      copy->nodeType = DIV;
      copy->child1 = left;
      copy->child2 = right;
      break;
    case VARIABLE:
      copy = allocateNode();
      copy->nodeType = MUL;
      tempNode = allocateNode();
      tempNode->nodeType = DIV;
      tempNode->child2 = right;
      tempNode2 = allocateNode();
      tempNode2->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      mpfr_set_d(*value,1.0,GMP_RNDN);
      tempNode2->value = value;
      tempNode->child1 = tempNode2;
      copy->child2 = left;
      copy->child1 = tempNode;
      break;
    case ADD:
      tempNode = allocateNode();
      tempNode->nodeType = ADD;
      tempNode2 = allocateNode();
      tempNode2->nodeType = DIV;
      tempNode3 = allocateNode();
      tempNode3->nodeType = DIV;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2->child2 = copyTree(right);
      tempNode3->child1 = copyTree(accessThruMemRef(left)->child2);
      tempNode3->child2 = right;
      tempNode->child1 = tempNode2;
      tempNode->child2 = tempNode3;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case SUB:
      tempNode = allocateNode();
      tempNode->nodeType = SUB;
      tempNode2 = allocateNode();
      tempNode2->nodeType = DIV;
      tempNode3 = allocateNode();
      tempNode3->nodeType = DIV;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2->child2 = copyTree(right);
      tempNode3->child1 = copyTree(accessThruMemRef(left)->child2);
      tempNode3->child2 = right;
      tempNode->child1 = tempNode2;
      tempNode->child2 = tempNode3;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case MUL:
      tempNode = allocateNode();
      tempNode->nodeType = MUL;
      tempNode2 = allocateNode();
      tempNode2->nodeType = DIV;
      tempNode3 = allocateNode();
      tempNode3->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      mpfr_set_d(*value,1.0,GMP_RNDN);
      tempNode3->value = value;
      tempNode4 = allocateNode();
      tempNode4->nodeType = MUL;
      tempNode2->child1 = tempNode3;
      tempNode2->child2 = right;
      tempNode->child1 = tempNode2;
      tempNode->child2 = tempNode4;
      tempNode4->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode4->child2 = copyTree(accessThruMemRef(left)->child2);
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case DIV:
      tempNode = allocateNode();
      tempNode->nodeType = DIV;
      tempNode->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2 = allocateNode();
      tempNode2->nodeType = MUL;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child2);
      tempNode2->child2 = right;
      tempNode->child2 = tempNode2;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    case NEG:
      tempNode = allocateNode();
      tempNode->nodeType = NEG;
      tempNode2 = allocateNode();
      tempNode2->nodeType = DIV;
      tempNode2->child1 = copyTree(accessThruMemRef(left)->child1);
      tempNode2->child2 = right;
      tempNode->child1 = tempNode2;
      free_memory(left);
      copy = expandPolynomialUnsafe(tempNode);
      free_memory(tempNode);
      break;
    default:
      if (isConstant(left)) {
	return copyTree(tree);
      } else {
	sollyaFprintf(stderr,"Error: expandPolynomialUnsafe: an error occurred on handling the DIV left rewritten expression subtree\n");
	exit(1);
      }
    }
    return copy;
    break;
  case NEG:
    left = expandPolynomialUnsafe(tree->child1);
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = left;
    return copy;
    break;
  default:
    if (isConstant(tree)) {
      return copyTree(tree);
    } else {
      sollyaFprintf(stderr,"Error: expandPolynomialUnsafe: an error occurred on handling the expression tree\n");
      exit(1);
    }
  }
}



node* expandPolynomial(node *tree) {
  node *temp, *temp2;
  if (getDegree(tree) < 0) return copyTree(tree);
  temp = expandPowerInPolynomialUnsafe(tree);
  temp2 = expandPolynomialUnsafe(temp);
  free_memory(temp);
  return temp2;
}

node* expandUnsimplified(node *tree) {
  node *copy;
  mpfr_t *value;
  mpfr_t temp;

  if (tree->nodeType == MEMREF) {
    return addMemRef(expandUnsimplified(getMemRefChild(tree)));
  }

  if (!isConstant(tree) && isPolynomial(tree)) return expandPolynomial(tree);

  switch (tree->nodeType) {
  case VARIABLE:
    copy = makeVariable();
    break;
  case CONSTANT:
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,tools_precision);
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    copy->value = value;
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = expand(tree->child1);
    copy->child2 = expand(tree->child2);
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = expand(tree->child1);
    copy->child2 = expand(tree->child2);
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = expand(tree->child1);
    copy->child2 = expand(tree->child2);
    break;
  case DIV:
    copy = allocateNode();
    copy->nodeType = DIV;
    copy->child1 = expand(tree->child1);
    copy->child2 = expand(tree->child2);
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = expand(tree->child1);
    break;
  case UNARY_BASE_FUNC:
    copy = allocateNode();
    copy->nodeType = UNARY_BASE_FUNC;
    copy->baseFun = tree->baseFun;
    copy->child1 = expand(tree->child1);
    break;
  case POW:
    copy = allocateNode();
    copy->nodeType = POW;
    copy->child1 = expand(tree->child1);
    copy->child2 = expand(tree->child2);
    break;
  case LIBRARYFUNCTION:
    copy = allocateNode();
    copy->nodeType = LIBRARYFUNCTION;
    copy->libFun = tree->libFun;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = expand(tree->child1);
    break;
  case PROCEDUREFUNCTION:
    copy = allocateNode();
    copy->nodeType = PROCEDUREFUNCTION;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = expand(tree->child1);
    copy->child2 = copyThing(tree->child2);
    break;
  case PI_CONST:
    copy = allocateNode();
    copy->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    copy = allocateNode();
    copy->nodeType = LIBRARYCONSTANT;
    copy->libFun = tree->libFun;
    break;
  default:
    sollyaFprintf(stderr,"Error: expand: unknown identifier in the tree\n");
    exit(1);
  }
  return copy;
}

node* expand(node *tree) {
  node *temp, *temp2, *temp3;
  temp3 = expandDivision(tree);
  temp = expandUnsimplified(temp3);
  temp2 = simplifyTreeErrorfree(temp);
  free_memory(temp);
  free_memory(temp3);
  return temp2;
}


int isConstant(node *tree) {
  int res;
  int r;

  switch (tree->nodeType) {
  case MEMREF:
    if (tree->cache->isConstantIsCached) {
      return tree->cache->isConstantCacheResult;
    }
    r = 0;
    if (tree->arguments != NULL) {
      r = 1;
    } else {
      if (tree->value != NULL) {
	r = 0;
      } else {
	if (tree->cache->polynomialRepresentation != NULL) {
	  res = (polynomialGetDegreeAsInt(tree->cache->polynomialRepresentation) == 0);
	} else {
	  res = isConstant(getMemRefChild(tree));
	}
	if (!res) {
	  tree->value = (mpfr_t *) (-1);
	}
	r = res;
      }
    }
    if (!tree->cache->isConstantIsCached) {
      tree->cache->isConstantCacheResult = r;
      tree->cache->isConstantIsCached = 1;
    }
    return r;
    break;
  case VARIABLE:
    return 0;
    break;
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    return 1;
    break;

  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    return (!(!((isConstant(tree->child1) && isConstant(tree->child2)))));
    break;

  case NEG:
  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return isConstant(tree->child1);
    break;
  default:
    sollyaFprintf(stderr,"Error: isConstant: unknown identifier in the tree\n");
    exit(1);
  }
}


int isMonomial(node *tree) {

  switch (tree->nodeType) {
  case MEMREF:
    return isMonomial(getMemRefChild(tree));
    break;
  case MUL:
    return (isMonomial(tree->child1) && isMonomial(tree->child2));
    break;
  case NEG:
    return isMonomial(tree->child1);
    break;
  case VARIABLE:
    return 1;
    break;
  case DIV:
    return (isConstant(tree->child2)) && isMonomial(tree->child1);
  default:
    return isConstant(tree);
  }
  return 0;
}



node* getCoefficientsInMonomialUnsafe(node *polynom) {
  node *leftSub, *rightSub, *coeffs;
  mpfr_t *value;

  if (isConstant(polynom)) return copyTree(polynom);

  if (polynom->nodeType == MEMREF) return getCoefficientsInMonomialUnsafe(getMemRefChild(polynom));

  if (polynom->nodeType == VARIABLE) return NULL;

  if (polynom->nodeType == MUL) {
    leftSub = getCoefficientsInMonomialUnsafe(polynom->child1);
    rightSub = getCoefficientsInMonomialUnsafe(polynom->child2);
    if ((leftSub == NULL) && (rightSub == NULL)) return NULL;
    if (leftSub == NULL) return rightSub;
    if (rightSub == NULL) return leftSub;
    coeffs = allocateNode();
    coeffs->nodeType = MUL;
    coeffs->child1 = leftSub;
    coeffs->child2 = rightSub;
    return coeffs;
  }

  if (polynom->nodeType == DIV) {
    leftSub = getCoefficientsInMonomialUnsafe(polynom->child1);
    if (leftSub == NULL) {
      coeffs = allocateNode();
      coeffs->nodeType = DIV;
      coeffs->child1 = allocateNode();
      coeffs->child1->nodeType = CONSTANT;
      coeffs->child1->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(coeffs->child1->value),10);
      mpfr_set_d(*(coeffs->child1->value),1.0,GMP_RNDN);
      coeffs->child2 = copyTree(polynom->child2);
    } else {
      coeffs = allocateNode();
      coeffs->nodeType = DIV;
      coeffs->child1 = leftSub;
      coeffs->child2 = copyTree(polynom->child2);
    }
    return coeffs;
  }


  if (polynom->nodeType == NEG) {
    leftSub = getCoefficientsInMonomialUnsafe(polynom->child1);
    rightSub = allocateNode();
    rightSub->nodeType = CONSTANT;
    value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    mpfr_set_d(*value,-1.0,GMP_RNDN);
    rightSub->value = value;
    if (leftSub == NULL) return rightSub;
    coeffs = allocateNode();
    coeffs->nodeType = MUL;
    coeffs->child1 = leftSub;
    coeffs->child2 = rightSub;
    return coeffs;
  }

  sollyaFprintf(stderr,"Error: getCoefficientsInMonomialUnsafe: an error occurred. The expression does not have the correct monomial form.\n");
  exit(1);
  return NULL;
}


void getCoefficientsUnsafe(node **monomials, node *polynom, int sign) {
  int degree;
  node *temp, *coeff, *temp2;
  mpfr_t *value;
  node *simplified, *simplifiedTemp;

  if (polynom->nodeType == MEMREF) {
    getCoefficientsUnsafe(monomials, getMemRefChild(polynom), sign);
    return;
  }

  if (isMonomial(polynom)) {
    degree = getDegree(polynom);
    coeff = getCoefficientsInMonomialUnsafe(polynom);
    if (coeff == NULL) {
      coeff = allocateNode();
      coeff->nodeType = CONSTANT;
      value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      mpfr_set_d(*value,1.0,GMP_RNDN);
      coeff->value = value;
    }
    temp = monomials[degree];
    if (temp == NULL) {
      if (sign < 0) {
	temp2 = allocateNode();
	temp2->nodeType = NEG;
	temp2->child1 = coeff;
	coeff = temp2;
      }
      temp = coeff;
    } else {
      temp2 = allocateNode();
      if (sign > 0) temp2->nodeType = ADD; else temp2->nodeType = SUB;
      temp2->child1 = temp;
      temp2->child2 = coeff;
      temp = temp2;
    }
    monomials[degree] = temp;
    return;
  }

  if (polynom->nodeType == ADD) {
    getCoefficientsUnsafe(monomials,polynom->child1,sign);
    getCoefficientsUnsafe(monomials,polynom->child2,sign);
    return;
  }

  if (polynom->nodeType == SUB) {
    getCoefficientsUnsafe(monomials,polynom->child1,sign);
    getCoefficientsUnsafe(monomials,polynom->child2,-sign);
    return;
  }

  if (polynom->nodeType == NEG) {
    getCoefficientsUnsafe(monomials,polynom->child1,-sign);
    return;
  }


  simplifiedTemp = expandPowerInPolynomialUnsafe(polynom);
  simplified = expandPolynomialUnsafe(simplifiedTemp);

  printMessage(7,SOLLYA_MSG_RECURSION_ON_POLY_COEFFICIENTS_EXTRACTION,"Warning: recursion on coefficients extraction: %b\ntransformed to\n%b\n",polynom,simplified);

  getCoefficientsUnsafe(monomials, simplified, sign);

  free_memory(simplifiedTemp);
  free_memory(simplified);

}

int isPowerOfVariable(node *);

void getCoefficientsHornerUnsafe(node **coefficients, node *poly, int offset, int sign) {
  int deg, newSign;
  node *newCoeff, *temp;

  if (poly->nodeType == MEMREF) {
    getCoefficientsHornerUnsafe(coefficients, getMemRefChild(poly), offset, sign);
    return;
  }

  if (isConstant(poly)) {
    newCoeff = copyTree(poly);
  } else {
    if (poly->nodeType == SUB) newSign = -1; else newSign = 1;
    newCoeff = copyTree(poly->child1);
    if ((accessThruMemRef(poly->child2)->nodeType == MUL) &&
	isConstant(accessThruMemRef(poly->child2)->child1) &&
	isPowerOfVariable(accessThruMemRef(poly->child2)->child2)) {
      deg = getDegree(accessThruMemRef(poly->child2)->child2);
      getCoefficientsHornerUnsafe(coefficients,accessThruMemRef(poly->child2)->child1,offset+deg,sign*newSign);
    } else {
      if ((accessThruMemRef(poly->child2)->nodeType == MUL) &&
	  (accessThruMemRef(accessThruMemRef(poly->child2)->child1)->nodeType == MUL) &&
	  isPowerOfVariable(accessThruMemRef(accessThruMemRef(poly->child2)->child1)->child1) &&
	  isConstant(accessThruMemRef(accessThruMemRef(poly->child2)->child1)->child2) &&
	  isConstant(accessThruMemRef(poly->child2)->child2)) {
	deg = getDegree(accessThruMemRef(accessThruMemRef(poly->child2)->child1)->child1);
	temp = allocateNode();
	temp->nodeType = MUL;
	temp->child1 = copyTree(accessThruMemRef(accessThruMemRef(poly->child2)->child1)->child2);
	temp->child2 = copyTree(accessThruMemRef(poly->child2)->child2);
	getCoefficientsHornerUnsafe(coefficients,temp,offset+deg,sign*newSign);
	free_memory(temp);
      } else {
	if (isPowerOfVariable(poly->child2)) {
	  deg = getDegree(poly->child2);
	  temp = allocateNode();
	  temp->nodeType = CONSTANT;
	  temp->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*(temp->value),17);
	  mpfr_set_d(*(temp->value),1.0,GMP_RNDN);
	  getCoefficientsHornerUnsafe(coefficients,temp,offset+deg,sign*newSign);
	  free_memory(temp);
	} else {
	  deg = getDegree(accessThruMemRef(poly->child2)->child1);
	  getCoefficientsHornerUnsafe(coefficients,accessThruMemRef(poly->child2)->child2,offset+deg,sign*newSign);
	}
      }
    }
  }

  if (coefficients[offset] == NULL) {
    if (sign == -1) {
      temp = allocateNode();
      temp->nodeType = NEG;
      temp->child1 = newCoeff;
      coefficients[offset] = temp;
    } else {
      coefficients[offset] = newCoeff;
    }
  } else {
    temp = allocateNode();
    if (sign == 1) temp->nodeType = ADD; else temp->nodeType = SUB;
    temp->child1 = coefficients[offset];
    temp->child2 = newCoeff;
    coefficients[offset] = temp;
  }

}

void getCoefficientsHorner(node **coefficients, node *poly) {
  int offset;

  if (poly->nodeType == MEMREF) {
    getCoefficientsHorner(coefficients, getMemRefChild(poly));
    return;
  }

  printMessage(7,SOLLYA_MSG_POLY_COEFF_EXTRACTION_SPECIAL_ALGO_FOR_HORNER,"Information: extraction of coefficient terms from a polynomial uses a special algorithm for Horner forms.\n");

  if (poly->nodeType == MUL) {
    offset = getDegree(poly->child1);
    getCoefficientsHornerUnsafe(coefficients,poly->child2,offset,1);
    return;
  }
  getCoefficientsHornerUnsafe(coefficients,poly,0,1);
}

int isPowerOfVariable(node *);
int isCanonicalMonomial(node *);

void getCoefficientsCanonicalUnsafe(node **coefficients, node *poly) {
  int deg, sign;
  node *newCoeff, *temp;

  if (poly->nodeType == MEMREF) {
    getCoefficientsCanonicalUnsafe(coefficients, getMemRefChild(poly));
    return;
  }

  if (isConstant(poly)) {
    sign = 1;
    deg = 0;
    newCoeff = copyTree(poly);
  } else {
    if (isCanonicalMonomial(poly)) {
      deg = getDegree(poly);
      sign = 1;
      if (isPowerOfVariable(poly)) {
        newCoeff = allocateNode();
        newCoeff->nodeType = CONSTANT;
        newCoeff->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(newCoeff->value),17);
        mpfr_set_d(*(newCoeff->value),1.0,GMP_RNDN);
      } else {
	newCoeff = copyTree(poly->child1);
      }
    } else {
      getCoefficientsCanonicalUnsafe(coefficients,poly->child1);
      if (poly->nodeType == ADD) sign = 1; else sign = -1;
      if (isConstant(poly->child2)) {
	deg = 0;
	newCoeff = copyTree(poly->child2);
      } else {
	deg = getDegree(poly->child2);
	if (isPowerOfVariable(poly->child2)) {
	  newCoeff = allocateNode();
	  newCoeff->nodeType = CONSTANT;
	  newCoeff->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*(newCoeff->value),17);
	  mpfr_set_d(*(newCoeff->value),1.0,GMP_RNDN);
	} else {
	  newCoeff = copyTree(accessThruMemRef(poly->child2)->child1);
	}
      }
    }
  }

  if (coefficients[deg] == NULL) {
    if (sign == -1) {
      temp = allocateNode();
      temp->nodeType = NEG;
      temp->child1 = newCoeff;
      coefficients[deg] = temp;
    } else {
      coefficients[deg] = newCoeff;
    }
  } else {
    temp = allocateNode();
    if (sign == 1) temp->nodeType = ADD; else temp->nodeType = SUB;
    temp->child1 = coefficients[deg];
    temp->child2 = newCoeff;
    coefficients[deg] = temp;
  }

}

void getCoefficientsCanonical(node **coefficients, node *poly) {

  printMessage(7,SOLLYA_MSG_POLY_COEFF_EXTRACT_SPECIAL_ALGO_FOR_CANONICAL,"Information: extraction of coefficient terms from a polynomial uses a special algorithm for canonical forms.\n");

  getCoefficientsCanonicalUnsafe(coefficients, poly);
}


int isHorner(node *);
int isCanonical(node *);
node *dividePolynomialByPowerOfVariableUnsafe(node *tree, int alpha);
node *makePowerOfConstant(node *constTree, int k);
node *makeBinomialCoefficient(unsigned int n, unsigned int k);

/* Computes the coefficients of the polynomial p^k where p is
// given by the coefficients in input
*/
void computePowerOfPolynomialCoefficients(int *degreeRes, node ***coeffRes,
                                          node **coeffs, int degree, int k) {
  int i, t;
  node **coeffsQ;
  int degreeQ;
  node *binom;
  node *constPow;
  node *factor;
  node *temp;

  *degreeRes = k * degree;
  *coeffRes = (node **) safeCalloc(*degreeRes+1,sizeof(node *));

  if (k == 0) {
    for (i=1;i<=*degreeRes;i++) {
      (*coeffRes)[i] = makeConstantDouble(0.0);
    }
    (*coeffRes)[0] = makeConstantDouble(1.0);
    return;
  }

  if (degree == 0) {
    for (i=1;i<=*degreeRes;i++) {
      (*coeffRes)[i] = makeConstantDouble(0.0);
    }
    if (coeffs[0] != NULL) {
      (*coeffRes)[0] = makePowerOfConstant(coeffs[0],k);
    } else {
      (*coeffRes)[0] = makeConstantDouble(0.0);
    }
    return;
  }

  if (k == 1) {
    for (i=0;i<=*degreeRes;i++) {
      if (coeffs[i] != NULL) {
        (*coeffRes)[i] = copyTree(coeffs[i]);
      } else {
        (*coeffRes)[i] = makeConstantDouble(0.0);
      }
    }
    return;
  }

  for (i=0;i<=*degreeRes;i++) {
    (*coeffRes)[i] = makeConstantDouble(0.0);
  }

  for (t=0;t<=k;t++) {
    if ((coeffs[0] != NULL) &&
        (!((accessThruMemRef(coeffs[0])->nodeType == CONSTANT) &&
           (mpfr_zero_p(*(accessThruMemRef(coeffs[0])->value)))))) {
      computePowerOfPolynomialCoefficients(&degreeQ, &coeffsQ,
                                           &(coeffs[1]), degree - 1, t);
      binom = makeBinomialCoefficient(k, t);
      constPow = makePowerOfConstant(coeffs[0],k-t);
      factor = makeMul(binom,constPow);
      for (i=t;i<=t+degreeQ;i++) {
        if (coeffsQ[i-t] != NULL) {
          (*coeffRes)[i] = makeAdd((*coeffRes)[i],
                                   makeMul(copyTree(factor),
                                           coeffsQ[i-t]));
        }
      }
      safeFree(coeffsQ);
      free_memory(factor);
    } else {
      if (k == t) {
        computePowerOfPolynomialCoefficients(&degreeQ, &coeffsQ,
                                             &(coeffs[1]), degree - 1, t);
        factor = makeBinomialCoefficient(k, t);
        for (i=t;i<=t+degreeQ;i++) {
          if (coeffsQ[i-t] != NULL) {
            (*coeffRes)[i] = makeAdd((*coeffRes)[i],
                                     makeMul(copyTree(factor),
                                             coeffsQ[i-t]));
          }
        }
        safeFree(coeffsQ);
        free_memory(factor);
      }
    }
  }

  for (i=0;i<=*degreeRes;i++) {
    if ((*coeffRes)[i] != NULL) {
      temp = simplifyTreeErrorfree((*coeffRes)[i]);
      free_memory((*coeffRes)[i]);
      (*coeffRes)[i] = temp;
    }
  }
}

int getCoefficientsInnerAlternate(int *degree, node ***coefficients, node *poly) {
  polynomial_t polynomial;
  int deg;
  unsigned int myDegree, t, k;
  node **myCoeffs;

  if (poly == NULL) return 0;
  if (!polynomialFromExpressionOnlyRealCoeffs(&polynomial, poly)) {
    return 0;
  }

  if (polynomialGetCoefficients(&myCoeffs, &myDegree, polynomial)) {
    deg = myDegree;
    t = deg;
    if ((deg >= 0) && (t == myDegree)) {
      for (k=0u;k<=myDegree;k++) {
	if ((accessThruMemRef(myCoeffs[k])->nodeType == CONSTANT) &&
	    mpfr_zero_p(*(accessThruMemRef(myCoeffs[k])->value))) {
	  free_memory(myCoeffs[k]);
	  myCoeffs[k] = NULL;
	}
      }
      *coefficients = myCoeffs;
      *degree = deg;
      polynomialFree(polynomial);
      return 1;
    } else {
      for (k=0u;k<=myDegree;k++) {
	free_memory(myCoeffs[k]);
      }
      safeFree(myCoeffs);
    }
  }

  polynomialFree(polynomial);

  return 0;
}

void getCoefficientsInner(int *degree, node ***coefficients, node *poly) {
  node *temp, *temp2, *temp3, *temp4;
  int i,k,j, mpd;
  node **coefficients1, **coefficients2;
  int degree1, degree2;
  mpfr_t y;

  if (getCoefficientsInnerAlternate(degree, coefficients, poly)) {
    return;
  }

  if (poly->nodeType == MEMREF) {
    getCoefficients(degree, coefficients, getMemRefChild(poly));
    return;
  }

  *degree = getDegree(poly);
  if (*degree < 0) {
    printMessage(1,SOLLYA_MSG_TRIED_TO_EXTRACT_COEFFS_OF_STH_NOT_POLYNOMIAL,"Warning: Tried to get coefficients of an expression that is not a polynomial.\n");
    return;
  }

  *coefficients = (node**) safeCalloc((*degree + 1),sizeof(node*));
  for (i=0;i<=*degree;i++) (*coefficients)[i] = NULL;


  if (isHorner(poly)) {
    getCoefficientsHorner(*coefficients,poly);
    return;
  }

  if (isCanonical(poly)) {
    getCoefficientsCanonical(*coefficients,poly);
    return;
  }


  if ((poly->nodeType == ADD) || (poly->nodeType == SUB)) {
    getCoefficients(&degree1, &coefficients1, poly->child1);
    getCoefficients(&degree2, &coefficients2, poly->child2);
    for (i=0;i<=degree1;i++) {
      (*coefficients)[i] = coefficients1[i];
    }
    safeFree(coefficients1);
    for (i=0;i<=degree2;i++) {
      if (coefficients2[i] != NULL) {
	if ((*coefficients)[i] == NULL) {
	  if (poly->nodeType == SUB) {
	    temp = allocateNode();
	    temp->nodeType = NEG;
	    temp->child1 = coefficients2[i];
	    (*coefficients)[i] = temp;
	  } else {
	    (*coefficients)[i] = coefficients2[i];
	  }

	} else {
	  temp = allocateNode();
	  temp->nodeType = poly->nodeType;
	  temp->child1 = (*coefficients)[i];
	  temp->child2 = coefficients2[i];
	  (*coefficients)[i] = temp;
	}
      }
    }
    safeFree(coefficients2);
    return;
  }

  if (poly->nodeType == MUL) {
    getCoefficients(&degree1, &coefficients1, poly->child1);
    getCoefficients(&degree2, &coefficients2, poly->child2);
    for (i=0;i<=degree1;i++) {
      for (k=0;k<=degree2;k++) {
	if ((coefficients1[i] != NULL) && (coefficients2[k] != NULL)) {
	  j = i + k;
	  temp = allocateNode();
	  temp->nodeType = MUL;
	  temp->child1 = copyTree(coefficients1[i]);
	  temp->child2 = copyTree(coefficients2[k]);
	  if ((*coefficients)[j] == NULL) {
	    (*coefficients)[j] = temp;
	  } else {
	    temp2 = allocateNode();
	    temp2->nodeType = ADD;
	    temp2->child1 = (*coefficients)[j];
	    temp2->child2 = temp;
	    (*coefficients)[j] = temp2;
	  }
	}
      }
    }
    for (i=0;i<=degree1;i++) free_memory(coefficients1[i]);
    for (i=0;i<=degree2;i++) free_memory(coefficients2[i]);
    safeFree(coefficients1);
    safeFree(coefficients2);
    return;
  }

  if ((poly->nodeType == POW) &&
      (accessThruMemRef(poly->child2)->nodeType == CONSTANT) &&
      (mpfr_integer_p(*(accessThruMemRef(poly->child2)->value)))) {
    k = mpfr_get_si(*(accessThruMemRef(poly->child2)->value),GMP_RNDN);
    mpfr_init2(y,8 * sizeof(int) + 10);
    mpfr_set_si(y,k,GMP_RNDN);
    if ((mpfr_cmp(y,*(accessThruMemRef(poly->child2)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(poly->child2)->value))) &&
	(k > 0)) {
      if ((mpd = getMaxPowerDivider(poly->child1)) > 0) {
        temp = dividePolynomialByPowerOfVariableUnsafe(poly->child1, mpd);
        temp2 = allocateNode();
        temp2->nodeType = POW;
        temp2->child1 = temp;
        temp2->child2 = copyTree(poly->child2);
        getCoefficients(&degree1, &coefficients1, temp2);
        free_memory(temp2);
        for (i=0;i<=degree1;i++)
          (*coefficients)[i + k * mpd] = coefficients1[i];
        safeFree(coefficients1);
        mpfr_clear(y);
        return;
      }

      if (k == 2) {
        getCoefficients(&degree1, &coefficients1, poly->child1);
        getCoefficients(&degree2, &coefficients2, poly->child1);
        for (i=0;i<=degree1;i++) {
          for (k=0;k<=degree2;k++) {
            if ((coefficients1[i] != NULL) && (coefficients2[k] != NULL)) {
              j = i + k;
              temp = allocateNode();
              temp->nodeType = MUL;
              temp->child1 = copyTree(coefficients1[i]);
              temp->child2 = copyTree(coefficients2[k]);
              if ((*coefficients)[j] == NULL) {
                (*coefficients)[j] = temp;
              } else {
                temp2 = allocateNode();
                temp2->nodeType = ADD;
                temp2->child1 = (*coefficients)[j];
                temp2->child2 = temp;
                (*coefficients)[j] = temp2;
              }
            }
          }
        }
        for (i=0;i<=degree1;i++) free_memory(coefficients1[i]);
        for (i=0;i<=degree2;i++) free_memory(coefficients2[i]);
        safeFree(coefficients1);
        safeFree(coefficients2);
        mpfr_clear(y);
        return;
      }

      getCoefficients(&degree1, &coefficients1, poly->child1);
      for (i=0;i<=degree1;i++) {
        if (coefficients1[i] == NULL)
          coefficients1[i] = makeConstantDouble(0.0);
      }

      computePowerOfPolynomialCoefficients(&degree2, &coefficients2,
                                           coefficients1, degree1,
                                           k);
      for (i=0;i<=degree2;i++) {
        if (coefficients2[i] != NULL) {
          temp = simplifyTreeErrorfree(coefficients2[i]);
          free_memory(coefficients2[i]);
          coefficients2[i] = temp;
        }
      }

      for (i=0;i<=degree2;i++) {
        if ((coefficients2[i] != NULL) &&
            (!((accessThruMemRef(coefficients2[i])->nodeType == CONSTANT) &&
               (mpfr_zero_p(*(accessThruMemRef(coefficients2[i])->value)))))) {
          (*coefficients)[i] = copyTree(coefficients2[i]);
        }
      }

      for (i=0;i<=degree1;i++) free_memory(coefficients1[i]);
      safeFree(coefficients1);
      for (i=0;i<=degree2;i++) free_memory(coefficients2[i]);
      safeFree(coefficients2);
      mpfr_clear(y);
      return;
    }
    mpfr_clear(y);
  }

  temp = simplifyTreeErrorfree(poly);
  temp2 = expandPowerInPolynomialUnsafe(temp);
  temp3 = expandPolynomialUnsafe(temp2);
  temp4 = simplifyTreeErrorfree(temp3);

  getCoefficientsUnsafe(*coefficients,temp4,1);

  free_memory(temp);
  free_memory(temp2);
  free_memory(temp3);
  free_memory(temp4);
}

void getCoefficients(int *degree, node ***coefficients, node *poly) {
  int i, deg;
  unsigned int myDegree, t, k;
  node **myCoeffs;

  if (poly->nodeType == MEMREF) {
    if (poly->cache->polynomialRepresentation == NULL) {
      tryRepresentAsPolynomial(poly);
    }
    if (poly->cache->polynomialRepresentation != NULL) {
      if (polynomialGetCoefficients(&myCoeffs, &myDegree, poly->cache->polynomialRepresentation)) {
	deg = myDegree;
	t = deg;
	if ((deg >= 0) && (t == myDegree)) {
	  for (k=0u;k<=myDegree;k++) {
	    if ((accessThruMemRef(myCoeffs[k])->nodeType == CONSTANT) &&
		mpfr_zero_p(*(accessThruMemRef(myCoeffs[k])->value))) {
	      free_memory(myCoeffs[k]);
	      myCoeffs[k] = NULL;
	    }
	  }
	  *coefficients = myCoeffs;
	  *degree = deg;
	  return;
	} else {
	  for (k=0u;k<=myDegree;k++) {
	    free_memory(myCoeffs[k]);
	  }
	  safeFree(myCoeffs);
	}
      }
    }
  }

  getCoefficientsInner(degree, coefficients, poly);

  if (*degree >= 0) {
    for (i=0;i<=*degree;i++) {
      if ((*coefficients)[i] != NULL) (*coefficients)[i] = addMemRef((*coefficients)[i]);
    }
  }
}


node* hornerPolynomialUnsafe(node *tree) {
  node *copy, *temp, *temp2, *temp3, *temp4, *simplified;
  node **monomials;
  int degree, i, k, e;
  mpfr_t *value;

  simplified = simplifyTreeErrorfree(tree);

  if (isHorner(simplified)) {
    degree = getDegree(simplified);
    monomials = (node**) safeCalloc((degree + 1),sizeof(node*));
    for (i=0;i<=degree;i++) monomials[i] = NULL;
    getCoefficientsHorner(monomials,simplified);
  } else {
    if (isCanonical(simplified)) {
      degree = getDegree(simplified);
      monomials = (node**) safeCalloc((degree + 1),sizeof(node*));
      for (i=0;i<=degree;i++) monomials[i] = NULL;
      getCoefficientsCanonical(monomials,simplified);
    } else {
      getCoefficients(&degree,&monomials,simplified);
    }
  }

  while ((degree >= 0) && (monomials[degree] == NULL)) degree--;
  if ((degree < 0) || (monomials[degree] == NULL)) {
    for (i=0;i<=degree;i++) {
      if (monomials[i] != NULL) free_memory(monomials[i]);
    }
    safeFree(monomials);
    return makeConstantInt(0);
  }

  copy = copyTree(monomials[degree]);

  for (i=degree-1;i>=0;i--) {
    if (monomials[i] == NULL) {
      if (i == 0) {
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp->child1 = temp2;
	temp->child2 = copy;
	copy = temp;
      } else {
	for (k=i-1;((monomials[k]==NULL) && (k > 0));k--);
	e = (i - k) + 1;
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp3 = allocateNode();
	temp3->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	if (mpfr_set_si(*value,e,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: rounding occurred on representing a monomial power exponent with %d bits.\n",
			 (int) tools_precision);
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	}
	temp3->value = value;
	temp4 = allocateNode();
	temp4->nodeType = POW;
	temp4->child1 = temp2;
	temp4->child2 = temp3;
	temp->child1 = temp4;
	temp->child2 = copy;
	copy = temp;
	if (monomials[k] != NULL) {
	  temp = allocateNode();
	  temp->nodeType = ADD;
	  temp->child1 = copyTree(monomials[k]);
	  temp->child2 = copy;
	  copy = temp;
	}
	i = k;
      }
    } else {
      temp = allocateNode();
      temp->nodeType = MUL;
      temp2 = makeVariable();
      temp->child1 = temp2;
      temp->child2 = copy;
      copy = temp;
      temp = allocateNode();
      temp->nodeType = ADD;
      temp->child1 = copyTree(monomials[i]);
      temp->child2 = copy;
      copy = temp;
    }
  }


  for (i=0;i<=degree;i++) {
    if (monomials[i] != NULL) free_memory(monomials[i]);
  }
  safeFree(monomials);

  free_memory(simplified);
  return copy;
}

node* dividePolynomialByPowerOfVariableUnsafe(node *tree, int alpha) {
  node *copy, *temp, *temp2, *temp3, *temp4, *simplified;
  node **monomials;
  int degree, i, k, e;
  mpfr_t *value;

  simplified = simplifyTreeErrorfree(tree);

  getCoefficients(&degree,&monomials,simplified);

  if (alpha > 0) {
    for (i=0;i<alpha;i++)
      if (monomials[i] != NULL) free_memory(monomials[i]);
    for (i=alpha;i<=degree;i++) {
      monomials[i-alpha] = monomials[i];
    }
    degree = degree - alpha;
  }

  while ((degree >= 0) && (monomials[degree] == NULL)) degree--;
  if ((degree < 0) || (monomials[degree] == NULL)) {
    for (i=0;i<=degree;i++) {
      if (monomials[i] != NULL) free_memory(monomials[i]);
    }
    safeFree(monomials);
    free_memory(simplified);
    return makeConstantInt(0);
  }

  copy = copyTree(monomials[degree]);

  for (i=degree-1;i>=0;i--) {
    if (monomials[i] == NULL) {
      if (i == 0) {
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp->child1 = temp2;
	temp->child2 = copy;
	copy = temp;
      } else {
	for (k=i-1;((monomials[k]==NULL) && (k > 0));k--);
	e = (i - k) + 1;
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp3 = allocateNode();
	temp3->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	if (mpfr_set_si(*value,e,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: rounding occurred on representing a monomial power exponent with %d bits.\n",
			 (int) tools_precision);
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	}
	temp3->value = value;
	temp4 = allocateNode();
	temp4->nodeType = POW;
	temp4->child1 = temp2;
	temp4->child2 = temp3;
	temp->child1 = temp4;
	temp->child2 = copy;
	copy = temp;
	if (monomials[k] != NULL) {
	  temp = allocateNode();
	  temp->nodeType = ADD;
	  temp->child1 = copyTree(monomials[k]);
	  temp->child2 = copy;
	  copy = temp;
	}
	i = k;
      }
    } else {
      temp = allocateNode();
      temp->nodeType = MUL;
      temp2 = makeVariable();
      temp->child1 = temp2;
      temp->child2 = copy;
      copy = temp;
      temp = allocateNode();
      temp->nodeType = ADD;
      temp->child1 = copyTree(monomials[i]);
      temp->child2 = copy;
      copy = temp;
    }
  }


  for (i=0;i<=degree;i++) {
    if (monomials[i] != NULL) free_memory(monomials[i]);
  }
  safeFree(monomials);

  free_memory(simplified);
  return copy;
}



node* hornerPolynomial(node *tree) {
  node *temp;

  if (isConstant(tree)) return copyTree(tree);
  if (getDegree(tree) < 0) return copyTree(tree);
  if (isHorner(tree)) {
    printMessage(7,SOLLYA_MSG_EXPR_NOT_HORNERIZED_AS_ALREADY_HORNERIZED,"Information: no Horner simplification will be performed because the given tree is already in Horner form.\n");
    return copyTree(tree);
  }

  temp = hornerPolynomialUnsafe(tree);

  return temp;
}


node* hornerUnsimplified(node *tree) {
  node *copy, *res;
  mpfr_t *value;
  mpfr_t temp;
  polynomial_t p;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      if (polynomialIsHornerized(tree->cache->polynomialRepresentation)) return copyTree(tree);
      if (tree->child1 == NULL) {
	p = polynomialHornerize(tree->cache->polynomialRepresentation);
	polynomialFree(tree->cache->polynomialRepresentation);
	tree->cache->polynomialRepresentation = p;
	return copyTree(tree);
      }
      res = addMemRefEvenOnNull(NULL);
      if (res != NULL) {
	res->cache->polynomialRepresentation = polynomialHornerize(tree->cache->polynomialRepresentation);
	copyTreeAnnotationsNoSimplifications(res, tree);
	return res;
      }
    }
    res = addMemRef(hornerUnsimplified(getMemRefChild(tree)));
    copyTreeAnnotationsNoSimplifications(res, tree);
    return res;
  }

  if (isPolynomial(tree)) return hornerPolynomial(tree);

  switch (tree->nodeType) {
  case VARIABLE:
    copy = makeVariable();
    break;
  case CONSTANT:
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,tools_precision);
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    copy->value = value;
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = horner(tree->child1);
    copy->child2 = horner(tree->child2);
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = horner(tree->child1);
    copy->child2 = horner(tree->child2);
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = horner(tree->child1);
    copy->child2 = horner(tree->child2);
    break;
  case DIV:
    copy = allocateNode();
    copy->nodeType = DIV;
    copy->child1 = horner(tree->child1);
    copy->child2 = horner(tree->child2);
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = horner(tree->child1);
    break;
  case UNARY_BASE_FUNC:
    copy = allocateNode();
    copy->nodeType = UNARY_BASE_FUNC;
    copy->baseFun = tree->baseFun;
    copy->child1 = horner(tree->child1);
    break;
  case POW:
    copy = allocateNode();
    copy->nodeType = POW;
    copy->child1 = horner(tree->child1);
    copy->child2 = horner(tree->child2);
    break;
  case LIBRARYFUNCTION:
    copy = allocateNode();
    copy->nodeType = LIBRARYFUNCTION;
    copy->libFun = tree->libFun;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = horner(tree->child1);
    break;
  case PROCEDUREFUNCTION:
    copy = allocateNode();
    copy->nodeType = PROCEDUREFUNCTION;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = horner(tree->child1);
    copy->child2 = copyThing(tree->child2);
    break;
  case PI_CONST:
    copy = allocateNode();
    copy->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    copy = allocateNode();
    copy->nodeType = LIBRARYCONSTANT;
    copy->libFun = tree->libFun;
    break;
  default:
    sollyaFprintf(stderr,"Error: horner: unknown identifier in the tree\n");
    exit(1);
  }
  return copy;
}


int isPowerOfVariable(node *tree) {
  if (tree->nodeType == MEMREF) return isPowerOfVariable(getMemRefChild(tree));
  if (tree->nodeType == VARIABLE) return 1;
  if ((tree->nodeType == POW) &&
      (accessThruMemRef(tree->child1)->nodeType == VARIABLE) &&
      (accessThruMemRef(tree->child2)->nodeType == CONSTANT) &&
      mpfr_integer_p(*(accessThruMemRef(tree->child2)->value))) return 1;
  return 0;
}

int isHornerUnsafe(node *tree) {
  if (isConstant(tree)) return 1;
  if (tree->nodeType == MEMREF) return isHornerUnsafe(getMemRefChild(tree));
  if (((tree->nodeType == ADD) || (tree->nodeType == SUB)) &&
      isConstant(tree->child1) &&
      (accessThruMemRef(tree->child2)->nodeType == MUL) &&
      isPowerOfVariable(accessThruMemRef(tree->child2)->child1) &&
      isHornerUnsafe(accessThruMemRef(tree->child2)->child2)) return 1;
  if (((tree->nodeType == ADD) || (tree->nodeType == SUB)) &&
      isConstant(tree->child1) &&
      isPowerOfVariable(tree->child2)) return 1;
  if (((tree->nodeType == ADD) || (tree->nodeType == SUB)) &&
      isConstant(tree->child1) &&
      (accessThruMemRef(tree->child2)->nodeType == MUL) &&
      (accessThruMemRef(accessThruMemRef(tree->child2)->child1)->nodeType == MUL) &&
      isPowerOfVariable(accessThruMemRef(accessThruMemRef(tree->child2)->child1)->child1) &&
      isConstant(accessThruMemRef(accessThruMemRef(tree->child2)->child1)->child2) &&
      isConstant(accessThruMemRef(tree->child2)->child2)) return 1;
  if (((tree->nodeType == ADD) || (tree->nodeType == SUB)) &&
      isConstant(tree->child1) &&
      (accessThruMemRef(tree->child2)->nodeType == MUL) &&
      isConstant(accessThruMemRef(tree->child2)->child1) &&
      isPowerOfVariable(accessThruMemRef(tree->child2)->child2)) return 1;
  return 0;
}

int isHorner(node *tree) {
  if (tree->nodeType == MEMREF) {
    if (((tree->child1 == NULL) || (tree->cache->memRefChildFromPolynomial)) && (tree->cache->polynomialRepresentation != NULL)) {
      return polynomialIsHornerized(tree->cache->polynomialRepresentation);
    }
    return isHorner(getMemRefChild(tree));
  }
  if ((tree->nodeType == ADD) || (tree->nodeType == SUB))
    return isHornerUnsafe(tree);
  if (tree->nodeType == MUL) {
    return isPowerOfVariable(tree->child1) && isHornerUnsafe(tree->child2);
  }
  return 0;
}

node* hornerInner(node *);

node* hornerWork(node *tree) {
  node *res;
  polynomial_t p;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      if (polynomialIsHornerized(tree->cache->polynomialRepresentation)) return copyTree(tree);
      if (tree->child1 == NULL) {
	p = polynomialHornerize(tree->cache->polynomialRepresentation);
	polynomialFree(tree->cache->polynomialRepresentation);
	tree->cache->polynomialRepresentation = p;
	return copyTree(tree);
      }
      res = addMemRefEvenOnNull(NULL);
      if (res != NULL) {
	res->cache->polynomialRepresentation = polynomialHornerize(tree->cache->polynomialRepresentation);
	copyTreeAnnotationsNoSimplifications(res, tree);
	return res;
      }
    }
  }

  res = addMemRef(hornerInner(tree));

  if ((tree != NULL) && (res != NULL) && (tree != res) &&
      (tree->nodeType == MEMREF) &&
      isSyntacticallyEqualCheap(tree,res)) {
    free_memory(res);
    res = copyTree(tree);
  }

  if (((tree->nodeType == MEMREF) &&
       (tree->cache->evaluationHook != NULL)) &&
      ((res->nodeType == MEMREF) &&
       (res->cache->evaluationHook == NULL))) {
    res->cache->isCorrectlyTyped = tree->cache->isCorrectlyTyped;
    addEvaluationHookFromCopy(&(res->cache->evaluationHook), tree->cache->evaluationHook);
    if ((res->cache->derivCache == NULL) && (tree->cache->derivCache != NULL)) {
      res->cache->derivCache = copyTree(tree->cache->derivCache);
    }
  }

  return res;
}

node* horner(node *tree) {
  node *res;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->hornerCache != NULL) {
      res = copyTree(tree->cache->hornerCache);
    } else {
      res = hornerWork(tree);
      if ((tree->cache->hornerCache != NULL) &&
	  (res->nodeType == MEMREF)) {
	tree->cache->hornerCache = copyTree(res);
      }
    }
  } else {
    res = hornerWork(tree);
  }

  if (((tree->nodeType == MEMREF) &&
       (tree->cache->evaluationHook != NULL)) &&
      ((res->nodeType == MEMREF) &&
       (res->cache->evaluationHook == NULL))) {
    res->cache->isCorrectlyTyped = tree->cache->isCorrectlyTyped;
    addEvaluationHookFromCopy(&(res->cache->evaluationHook), tree->cache->evaluationHook);
    if ((res->cache->derivCache == NULL) && (tree->cache->derivCache != NULL)) {
      res->cache->derivCache = copyTree(tree->cache->derivCache);
    }
  }

  return res;
}

node* hornerInner(node *tree) {
  node *temp, *temp2, *temp3;
  int i;

  if (isHorner(tree) || isPowerOfVariable(tree)) {
    printMessage(7,SOLLYA_MSG_EXPR_NOT_HORNERIZED_AS_ALREADY_HORNERIZED,"Information: no Horner simplification will be performed because the given tree is already in Horner form.\n");
    return copyTree(tree);
  }

  temp3 = simplifyTreeErrorfree(tree);

  i = 0;
  do {
    temp = hornerUnsimplified(temp3);
    temp2 = simplifyTreeErrorfree(temp);
    free_memory(temp);
    free_memory(temp3);
    temp3 = temp2;
    i++;
  } while ((!(isHorner(temp3) || (!isPolynomial(temp3)))) && (i < 3));

  return temp3;
}


node *differentiatePolynomialHornerUnsafe(node *tree) {
  int degree, i, k, e;
  node **monomials;
  node *temp, *temp2, *temp3, *temp4, *copy;
  mpfr_t *value;
  mp_prec_t prec;

  copy = NULL;

  getCoefficients(&degree,&monomials,tree);

  if (degree == 0) {
    for (i=0;i<=degree;i++) {
      free_memory(monomials[i]);
    }
    safeFree(monomials);
    return addMemRef(makeConstantInt(0));
  }

  if (monomials[0] != NULL) free_memory(monomials[0]);

  for (i=1;i<=degree;i++) {
    if (monomials[i] != NULL) {
      if (accessThruMemRef(monomials[i])->nodeType == CONSTANT) {
	value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,mpfr_get_prec(*(accessThruMemRef(monomials[i])->value))+(sizeof(int)*8));
	if (mpfr_mul_si(*value,*(accessThruMemRef(monomials[i])->value),i,GMP_RNDN) != 0)
	  printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_HORNER_POLY,"Warning: rounding occurred while differentiating a polynomial in Horner form.\n");
	free_memory(monomials[i]);
	monomials[i] = allocateNode();
	monomials[i]->nodeType = CONSTANT;
	monomials[i]->value = value;
	temp = monomials[i];
      } else {
	temp = allocateNode();
	temp->nodeType = MUL;
	temp->child1 = allocateNode();
	temp->child1->nodeType = CONSTANT;
	temp->child1->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	prec = 8 * sizeof(int) + 10; if (tools_precision > prec) prec = tools_precision;
	mpfr_init2(*(temp->child1->value),prec);
	if (mpfr_set_si(*(temp->child1->value),i,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_HORNER_POLY,"Warning: on differentiating a polynomial in Horner form rounding occurred while representing the degree of a monomial on a constant of the given precision\n");
	  }
	}
	temp->child2 = monomials[i];
      }
      monomials[i-1] = temp;
    } else {
      monomials[i-1] = NULL;
    }
  }

  degree--;

  copy = copyTree(monomials[degree]);

  for (i=degree-1;i>=0;i--) {
    if (monomials[i] == NULL) {
      if (i == 0) {
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp->child1 = temp2;
	temp->child2 = copy;
	copy = temp;
      } else {
	for (k=i-1;((monomials[k]==NULL) && (k > 0));k--);
	e = (i - k) + 1;
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp3 = allocateNode();
	temp3->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,tools_precision);
	if (mpfr_set_si(*value,e,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: rounding occurred on representing a monomial power exponent with %d bits.\n",
			 (int) tools_precision);
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	}
	temp3->value = value;
	temp4 = allocateNode();
	temp4->nodeType = POW;
	temp4->child1 = temp2;
	temp4->child2 = temp3;
	temp->child1 = temp4;
	temp->child2 = copy;
	copy = temp;
	if (monomials[k] != NULL) {
	  temp = allocateNode();
	  temp->nodeType = ADD;
	  temp->child1 = copyTree(monomials[k]);
	  temp->child2 = copy;
	  copy = temp;
	}
	i = k;
      }
    } else {
      temp = allocateNode();
      temp->nodeType = MUL;
      temp2 = makeVariable();
      temp->child1 = temp2;
      temp->child2 = copy;
      copy = temp;
      temp = allocateNode();
      temp->nodeType = ADD;
      temp->child1 = copyTree(monomials[i]);
      temp->child2 = copy;
      copy = temp;
    }
  }


  for (i=0;i<=degree;i++) {
    if (monomials[i] != NULL) free_memory(monomials[i]);
  }
  safeFree(monomials);

  if (copy == NULL) {
    copy = addMemRef(makeConstantInt(0));
  }

  return copy;
}

node *differentiatePolynomialUnsafe(node *tree) {
  node *copy, *temp, *temp2, *temp3, *temp4, *temp5;
  int degree, i;
  node **monomials;
  mpfr_t *value;

  copy = NULL;

  if (isHorner(tree)) {
    printMessage(25,SOLLYA_MSG_DIFFERENTIATION_USES_SPECIAL_ALGO_FOR_HORNER,"Information: differentiating a polynomial in Horner form uses a special algorithm.\n");
    return differentiatePolynomialHornerUnsafe(tree);
  }

  degree = getDegree(tree);

  if (degree == 0) {
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
    mpfr_set_d(*value,0.0,GMP_RNDN);
    copy->value = value;
  } else {
    getCoefficients(&degree,&monomials,tree);

    if (monomials[degree] == NULL) {
      copy = allocateNode();
      copy->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
      mpfr_set_d(*value,0.0,GMP_RNDN);
      copy->value = value;
      monomials[degree] = copy;
    }

    if (degree >= 2) {
      if (degree > 2) {
	temp = copyTree(monomials[degree]);
	temp2 = allocateNode();
	temp2->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
	if (mpfr_set_si(*value,degree,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_POLYNOMIAL,"Warning: rounding occurred on differentiating a polynomial. A constant could not be written on %d bits.\n",
			 (int) ((tools_precision >= 128) ? (tools_precision) : (128)));
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	}
	temp2->value = value;
	temp3 = allocateNode();
	temp3->nodeType = MUL;
	temp3->child1 = temp2;
	temp3->child2 = temp;
	temp2 = allocateNode();
	temp2->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
	if (mpfr_set_si(*value,degree-1,GMP_RNDN) != 0) {
	  printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_POLYNOMIAL,
                       "Warning: rounding occurred on differentiating a polynomial. An exponent constant could not be written on %d bits.\n",
                       (int) ((tools_precision >= 128) ? (tools_precision) : (128)));
	  printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	}
	temp2->value = value;
	temp = makeVariable();
	temp4 = allocateNode();
	temp4->nodeType = POW;
	temp4->child1 = temp;
	temp4->child2 = temp2;
	temp5 = allocateNode();
	temp5->nodeType = MUL;
	temp5->child1 = temp3;
	temp5->child2 = temp4;
	copy = temp5;
      } else {
	temp = copyTree(monomials[degree]);
	temp2 = allocateNode();
	temp2->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
	if (mpfr_set_si(*value,degree,GMP_RNDN) != 0) {
	  printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_POLYNOMIAL,"Warning: rounding occurred on differentiating a polynomial. A constant could not be written on %d bits.\n",
                       (int) ((tools_precision >= 128) ? (tools_precision) : (128)));
	  printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	}
	temp2->value = value;
	temp3 = allocateNode();
	temp3->nodeType = MUL;
	temp3->child1 = temp2;
	temp3->child2 = temp;
	temp = makeVariable();
	temp4 = allocateNode();
	temp4->nodeType = MUL;
	temp4->child1 = temp3;
	temp4->child2 = temp;
	copy = temp4;
      }

      for (i=degree-1;i>1;i--) {
	if (monomials[i] != NULL) {
	  temp = copyTree(monomials[i]);
	  temp2 = allocateNode();
	  temp2->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
	  if (mpfr_set_si(*value,i,GMP_RNDN) != 0) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_POLYNOMIAL,"Warning: rounding occurred on differentiating a polynomial. A constant could not be written on %d bits.\n",
                         (int) ((tools_precision >= 128) ? (tools_precision) : (128)));
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	  temp2->value = value;
	  temp3 = allocateNode();
	  temp3->nodeType = MUL;
	  temp3->child1 = temp2;
	  temp3->child2 = temp;
	  temp2 = allocateNode();
	  temp2->nodeType = CONSTANT;
	  value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*value,((tools_precision >= 128) ? (tools_precision) : (128)));
	  if (mpfr_set_si(*value,i-1,GMP_RNDN) != 0) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_POLYNOMIAL,
                         "Warning: rounding occurred on differentiating a polynomial. An exponent constant could not be written on %d bits.\n",
                         (int) ((tools_precision >= 128) ? (tools_precision) : (128)));
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	  temp2->value = value;
	  temp = makeVariable();
	  temp4 = allocateNode();
	  temp4->nodeType = POW;
	  temp4->child1 = temp;
	  temp4->child2 = temp2;
	  temp5 = allocateNode();
	  temp5->nodeType = MUL;
	  temp5->child1 = temp3;
	  temp5->child2 = temp4;
	  temp = copy;
	  temp2 = allocateNode();
	  temp2->nodeType = ADD;
	  temp2->child1 = temp5;
	  temp2->child2 = temp;
	  copy = temp2;
	}
      }

      if (monomials[1] != NULL) {
	temp = copyTree(monomials[1]);
	temp2 = allocateNode();
	temp2->nodeType = ADD;
	temp2->child1 = temp;
	temp2->child2 = copy;
	copy = temp2;
      }
    } else {
      if (degree >= 1) {
	copy = copyTree(monomials[1]);
      } else {
	copy = addMemRef(makeConstantInt(0));
      }
    }

    for (i=0;i<=degree;i++) {
      if (monomials[i] != NULL) free_memory(monomials[i]);
    }
    safeFree(monomials);
  }

  if (copy == NULL) {
    copy = addMemRef(makeConstantInt(0));
  }

  return copy;
}


int getNumeratorDenominator(node **numerator, node **denominator, node *tree) {
  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      if (polynomialGetDegreeAsInt(tree->cache->polynomialRepresentation) != 0) {
	*numerator = copyTree(tree);
	*denominator = NULL;
	return 0;
      }
    }
    return getNumeratorDenominator(numerator, denominator, getMemRefChild(tree));
  }
  if (tree->nodeType == DIV) {
    *numerator = copyTree(tree->child1);
    *denominator = copyTree(tree->child2);
    return 1;
  }
  else {
    *numerator = copyTree(tree);
    *denominator = NULL;
    return 0;
  }
}

node *makeBinomialCoefficient(unsigned int n, unsigned int k) {
  mpz_t coeffGMP;
  mp_prec_t prec;
  mpfr_t *coeffVal;
  node *res;

  mpz_init(coeffGMP);
  mpz_bin_uiui(coeffGMP,n,k);
  prec = mpz_sizeinbase(coeffGMP, 2) + 10;
  if (prec < tools_precision) prec = tools_precision;
  coeffVal = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*coeffVal,prec);
  if(mpfr_set_z(*coeffVal,coeffGMP,GMP_RNDN) != 0) {
    if (!noRoundingWarnings) {
      printMessage(1,SOLLYA_MSG_ROUNDING_UPON_BINOMIAL_COEFFICIENT_COMPUTATION,"Warning: rounding occurred when calculating a binomial coefficient.\n");
      printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the working precision.\n");
    }
  }
  mpz_clear(coeffGMP);
  res = allocateNode();
  res->nodeType = CONSTANT;
  res->value = coeffVal;
  return res;
}

node *makePowerOfConstant(node *constTree, int k) {
  node *temp, *res;

  if (k == 1) {
    return copyTree(constTree);
  }
  temp = allocateNode();
  temp->nodeType = CONSTANT;
  temp->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(temp->value),8 * sizeof(k) + 10);
  mpfr_set_si(*(temp->value),k,GMP_RNDN);
  res = allocateNode();
  res->nodeType = POW;
  res->child1 = copyTree(constTree);
  res->child2 = temp;

  return res;
}

/* polynomialShiftAndScaleAbscissaUnsafe(poly, a, b) returns
//
// p(a + b * x)
//
// for a polynomial p and a, b constant nodes
//
// If poly is not a polynomial the tool is terminated.
// If a or b is not a constant expression, the shifting and
// scaling is performed as if they were.
*/
node *polynomialShiftAndScaleAbscissaUnsafe(node *poly, node *a, node *b) {
  node *res;
  node **coeffs;
  node **coeffsRes;
  node *temp;
  int degree;
  int i,k;

  getCoefficients(&degree, &coeffs, poly);
  if (degree < 0) {
    sollyaFprintf(stderr,"Error: polynomialShiftAndScaleAbscissaUnsafe: the given expression is not a polynomial\n");
    exit(1);
  }
  for (i=0;i<=degree;i++) {
    if (coeffs[i] == NULL) {
      coeffs[i] = makeConstantDouble(0.0);
    }
  }

  coeffsRes = (node **) safeCalloc(degree+1,sizeof(node *));
  for (i=0;i<=degree;i++) {
    coeffsRes[i] = makeConstantDouble(0.0);
  }

  for (i=0;i<=degree;i++) {
    for (k=0;k<=i;k++) {
      temp = makeMul(copyTree(coeffs[i]),
                     makeMul(makeBinomialCoefficient(i,k),
                             makeMul(makePowerOfConstant(a, i-k),
                                     makePowerOfConstant(b, k))));
      coeffsRes[k] = makeAdd(coeffsRes[k],temp);
    }
  }

  for (i=0;i<=degree;i++) {
    if (coeffsRes[i] != NULL) {
      temp = simplifyTreeErrorfree(coeffsRes[i]);
      free_memory(coeffsRes[i]);
      coeffsRes[i] = temp;
    }
  }

  res = makePolynomialConstantExpressions(coeffsRes, degree);

  for (i=0;i<=degree;i++) {
    if (coeffs[i] != NULL) free_memory(coeffs[i]);
    if (coeffsRes[i] != NULL) free_memory(coeffsRes[i]);
  }
  safeFree(coeffs);
  safeFree(coeffsRes);

  return res;
}

/* Returns p(q)
 */
node *substitutePolynomialUnsafe(node *p, node *q) {
  node *res;
  node **coeffsP, **coeffsQ, **coeffs, **coeffsQPi;
  int degP, degQ, deg, i, k, degQPi;
  node *temp;

  getCoefficients(&degP, &coeffsP, p);
  if (degP < 0) {
    sollyaFprintf(stderr,"Error: substitutePolynomialUnsafe: the given expression is not a polynomial\n");
    exit(1);
  }

  getCoefficients(&degQ, &coeffsQ, q);
  if (degQ < 0) {
    sollyaFprintf(stderr,"Error: substitutePolynomialUnsafe: the given expression is not a polynomial\n");
    exit(1);
  }
  for (i=0;i<=degQ;i++) {
    if (coeffsQ[i] == NULL) {
      coeffsQ[i] = makeConstantDouble(0.0);
    }
  }

  deg = degP * degQ;
  coeffs = (node **) safeCalloc(deg+1,sizeof(node *));
  for (i=0;i<=deg;i++) {
    coeffs[i] = makeConstantDouble(0.0);
  }

  for (i=0;i<=degP;i++) {
    if (coeffsP[i] != NULL) {
      computePowerOfPolynomialCoefficients(&degQPi, &coeffsQPi,
                                           coeffsQ, degQ, i);
      for (k=0;k<=degQPi;k++) {
        if (coeffsQPi[k] != NULL) {
          coeffs[k] = makeAdd(coeffs[k],
                              makeMul(copyTree(coeffsP[i]),
                                      coeffsQPi[k]));
        }
      }
      safeFree(coeffsQPi);
    }
  }

  for (i=0;i<=deg;i++) {
    if (coeffs[i] != NULL) {
      temp = simplifyTreeErrorfree(coeffs[i]);
      free_memory(coeffs[i]);
      coeffs[i] = temp;
    }
  }

  res = makePolynomialConstantExpressions(coeffs, deg);

  for (i=0;i<=degP;i++) {
    if (coeffsP[i] != NULL) free_memory(coeffsP[i]);
  }
  safeFree(coeffsP);

  for (i=0;i<=degQ;i++) {
    if (coeffsQ[i] != NULL) free_memory(coeffsQ[i]);
  }
  safeFree(coeffsQ);

  for (i=0;i<=deg;i++) {
    if (coeffs[i] != NULL) free_memory(coeffs[i]);
  }
  safeFree(coeffs);

  return res;
}

node *substituteInner(node* tree, node *t, int doNotEvaluate, int maySimplify);

static inline node *substituteEnhancedInner(node* tree, node *t, int doNotEvaluate, int maySimplify) {
  node *res;

  if (maySimplify) {
    if ((tree->nodeType == MEMREF) &&
	(t->nodeType == MEMREF)) {
      if (tree->cache->polynomialRepresentation == NULL) {
	tryRepresentAsPolynomial(tree);
      }
      if ((tree->cache->polynomialRepresentation != NULL) &&
	  (t->cache->polynomialRepresentation == NULL)) {
	tryRepresentAsPolynomial(t);
      }
      if ((tree->cache->polynomialRepresentation != NULL) &&
	  (t->cache->polynomialRepresentation != NULL)) {
	res = addMemRefEvenOnNull(NULL);
	if (res != NULL) {
	  res->cache->polynomialRepresentation = polynomialCompose(tree->cache->polynomialRepresentation,
								   t->cache->polynomialRepresentation);
	  return res;
	}
      }
    }
  }

  return substituteInner(tree, t, doNotEvaluate, maySimplify);
}

node *substituteEnhanced(node* tree, node *t, int doNotEvaluate, int maySimplify) {
  node *res;

  res = addMemRef(substituteEnhancedInner(tree, t, doNotEvaluate, maySimplify));

  return res;
}

node *substitute(node* tree, node *t) {
  node *tt, *res;
  int freeTT;

  freeTT = 0;
  if ((t != NULL) && (t->nodeType == MEMREF)) {
    tt = t;
    freeTT = 0;
  } else {
    if (t == NULL) {
      tt = t;
      freeTT = 0;
    } else {
      tt = addMemRef(copyTree(t));
      freeTT = 1;
    }
  }

  res = substituteEnhanced(tree, tt, 0, 1);

  if (freeTT) {
    free_memory(tt);
  }

  return res;
}

static inline int treeContainsHooksInner(node *tree, unsigned int call) {
  int res;

  if (tree == NULL) return 0;
  if (tree->nodeType == MEMREF) {
    if (tree->cache->containsHooksIsCached &&
	(tree->cache->containsHooksCall == call)) {
      return tree->cache->containsHooksCacheResult;
    }
    res = 0;
    if (tree->cache->evaluationHook != NULL) {
      res = 1;
    } else {
      if ((tree->child1 == NULL) &&
	  (tree->cache->polynomialRepresentation != NULL)) {
	res = 0;
      } else {
	res = treeContainsHooksInner(getMemRefChild(tree),call);
      }
    }
    if (!(tree->cache->containsHooksIsCached &&
	  (tree->cache->containsHooksCall == call))) {
      tree->cache->containsHooksCacheResult = res;
      tree->cache->containsHooksCall = call;
      tree->cache->containsHooksIsCached = 1;
    }
    return res;
  }

  switch (tree->nodeType) {
  case VARIABLE:
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    return 0;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    if (treeContainsHooksInner(tree->child1,call)) return 1;
    if (treeContainsHooksInner(tree->child2,call)) return 1;
    return 0;
    break;
  case UNARY_BASE_FUNC:
  case NEG:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return treeContainsHooksInner(tree->child1,call);
    break;
  default:
    sollyaFprintf(stderr,"Error: treeContainsHooksInner: unknown identifier in the tree\n");
    exit(1);
  }

  return 0;
}

unsigned int treeContainsHooksGlobalCallCounter = 0u;

int treeContainsHooks(node *tree) {
  int res;

  treeContainsHooksGlobalCallCounter++;

  res = treeContainsHooksInner(tree, treeContainsHooksGlobalCallCounter);

  return res;
}

node *substituteInner(node* tree, node *t, int doNotEvaluate, int maySimplify) {
  node *copy;
  mpfr_t *value;
  mpfr_t temp;
  node **coeffs;
  int degree;
  int i;
  sollya_mpfi_t tEval, treeEval;
  mp_prec_t treeEvalPrec;
  mpfr_t tEl, tEr;
  sollya_mpfi_t aPrioriBoundForConstantExpr;
  int haveAPrioriBoundForConstantExpr;
  mp_prec_t *precPtr;
  sollya_mpfi_t *intervalPtr;
  node *tempDerivCache;


  haveAPrioriBoundForConstantExpr = 0;

  if (isPolynomial(tree) &&
      isPolynomial(t)) {
    if ((getDegree(t) == 1) &&
        (getDegree(tree) >= 2)) {
      getCoefficients(&degree, &coeffs, t);
      if (degree == 1) {
        for (i=0;i<=degree;i++) {
          if (coeffs[i] == NULL) {
            coeffs[i] = makeConstantDouble(0.0);
          }
        }

        copy = polynomialShiftAndScaleAbscissaUnsafe(tree, coeffs[0], coeffs[1]);

        for (i=0;i<=degree;i++) {
          if (coeffs[i] != NULL) free_memory(coeffs[i]);
        }
        safeFree(coeffs);
        return copy;
      }
      for (i=0;i<=degree;i++) {
        if (coeffs[i] != NULL) free_memory(coeffs[i]);
      }
      safeFree(coeffs);
    }

    if ((getDegree(t) >= 2) &&
        (getDegree(tree) >= 2)) {
      copy = substitutePolynomialUnsafe(tree,t);
      return copy;
    }
  }

  if (isConstant(t) && (!isConstant(tree)) && (!doNotEvaluate)) {
    copy = NULL;
    sollya_mpfi_init2(tEval, tools_precision + 10);
    sollya_mpfi_init2(treeEval, tools_precision + 10);

    evaluateConstantExpressionToInterval(tEval, t);
    evaluateInterval(treeEval, tree, NULL, tEval);

    treeEvalPrec = sollya_mpfi_get_prec(treeEval);
    mpfr_init2(tEl, treeEvalPrec);
    mpfr_init2(tEr, treeEvalPrec);
    sollya_mpfi_get_left(tEl, treeEval);
    sollya_mpfi_get_right(tEr, treeEval);

    if (mpfr_number_p(tEr) &&
	mpfr_number_p(tEl)) {
      sollya_mpfi_init2(aPrioriBoundForConstantExpr,sollya_mpfi_get_prec(treeEval));
      sollya_mpfi_set(aPrioriBoundForConstantExpr,treeEval);
      haveAPrioriBoundForConstantExpr = 1;
      if (mpfr_equal_p(tEr, tEl)) {
	copy = addMemRef(makeConstant(tEr));
      }
    }

    mpfr_clear(tEl);
    mpfr_clear(tEr);
    sollya_mpfi_clear(tEval);
    sollya_mpfi_clear(treeEval);
    if (copy != NULL) {
      if (haveAPrioriBoundForConstantExpr) {
	sollya_mpfi_clear(aPrioriBoundForConstantExpr);
	haveAPrioriBoundForConstantExpr = 0;
      }
      return copy;
    }
  }

  if (tree->nodeType == MEMREF) {
    if ((tree->cache->substituteCacheY != NULL) &&
	(tree->cache->substituteCacheX != NULL) &&
	(tree->cache->substituteCacheX == t) &&
	((!!tree->cache->substituteCacheMaySimplify) == (!!maySimplify))) {
      copy = addMemRef(copyTree(tree->cache->substituteCacheY));
    } else {
      copy = addMemRef(substituteInner(getMemRefChild(tree), t, 1, maySimplify));

      if ((copy != NULL) &&
	  (t != NULL) &&
	  (t != copy) &&
	  (tree != copy) &&
	  (tree != t)) {
	if ((copy->nodeType == MEMREF) &&
	    (t->nodeType == MEMREF)) {
	  if (tree->cache->substituteCacheY != NULL) {
	    freeThing(tree->cache->substituteCacheY);
	    tree->cache->substituteCacheY = NULL;
	  }
	  if (tree->cache->substituteCacheX != NULL) {
	    freeThing(tree->cache->substituteCacheX);
	    tree->cache->substituteCacheX = NULL;
	  }
	  tree->cache->substituteCacheX = addMemRef(copyTree(t));
	  tree->cache->substituteCacheY = addMemRef(copyTree(copy));
	  tree->cache->substituteCacheMaySimplify = !!maySimplify;
	}
      }

      if ((copy->nodeType == MEMREF) && (tree->cache->evaluationHook != NULL)) {
	addEvaluationHookFromComposition(&(copy->cache->evaluationHook), tree->cache->evaluationHook, t);

	if (((copy->cache->derivCache == NULL) &&
	     ((tree->cache->derivCache != NULL) &&
	      ((tree->cache->derivCache->nodeType == MEMREF) &&
	       (treeContainsHooks(tree->cache->derivCache))))) &&
	    (!isConstant(copy))) {
	  tempDerivCache = addMemRef(makeMul(addMemRef(substituteInner(tree->cache->derivCache, t, 1, maySimplify)),
					     addMemRef(differentiate(t))));
	  copy->cache->derivCache = simplifyTreeErrorfree(tempDerivCache);
	  free_memory(tempDerivCache);
	}
      }

      if (haveAPrioriBoundForConstantExpr) {
	if (copy->nodeType == MEMREF) {
	  if (copy->arguments == NULL) {
	    precPtr = (mp_prec_t *) safeMalloc(sizeof(mp_prec_t));
	    *precPtr = sollya_mpfi_get_prec(aPrioriBoundForConstantExpr);
	    intervalPtr = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	    sollya_mpfi_init2(*intervalPtr, sollya_mpfi_get_prec(aPrioriBoundForConstantExpr));
	    sollya_mpfi_set(*intervalPtr, aPrioriBoundForConstantExpr);
	    copy->arguments = addElement(addElement(NULL, intervalPtr), precPtr);
	  }
	}
	sollya_mpfi_clear(aPrioriBoundForConstantExpr);
	haveAPrioriBoundForConstantExpr = 0;
      }
    }
    if (haveAPrioriBoundForConstantExpr) {
      sollya_mpfi_clear(aPrioriBoundForConstantExpr);
      haveAPrioriBoundForConstantExpr = 0;
    }
    return copy;
  }

  switch (tree->nodeType) {
  case VARIABLE:
    copy = copyTree(t);
    break;
  case CONSTANT:
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,mpfr_get_prec(*(tree->value)));
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    copy->value = value;
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    copy->child2 = substituteEnhanced(tree->child2,t,doNotEvaluate,maySimplify);
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    copy->child2 = substituteEnhanced(tree->child2,t,doNotEvaluate,maySimplify);
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    copy->child2 = substituteEnhanced(tree->child2,t,doNotEvaluate,maySimplify);
    break;
  case DIV:
    copy = allocateNode();
    copy->nodeType = DIV;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    copy->child2 = substituteEnhanced(tree->child2,t,doNotEvaluate,maySimplify);
    break;
  case UNARY_BASE_FUNC:
    copy = allocateNode();
    copy->nodeType = UNARY_BASE_FUNC;
    copy->baseFun = tree->baseFun;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    break;
  case POW:
    copy = allocateNode();
    copy->nodeType = POW;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    copy->child2 = substituteEnhanced(tree->child2,t,doNotEvaluate,maySimplify);
    break;
  case LIBRARYFUNCTION:
    copy = allocateNode();
    copy->nodeType = LIBRARYFUNCTION;
    copy->libFun = tree->libFun;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    break;
  case PROCEDUREFUNCTION:
    copy = allocateNode();
    copy->nodeType = PROCEDUREFUNCTION;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = substituteEnhanced(tree->child1,t,doNotEvaluate,maySimplify);
    copy->child2 = copyThing(tree->child2);
    break;
  case PI_CONST:
    copy = allocateNode();
    copy->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    copy = allocateNode();
    copy->nodeType = LIBRARYCONSTANT;
    copy->libFun = tree->libFun;
    break;

  default:
    sollyaFprintf(stderr,"Error: substituteInner: unknown identifier in the tree\n");
    exit(1);
  }
  copy = addMemRef(copy);

  if (haveAPrioriBoundForConstantExpr) {
    if (copy->nodeType == MEMREF) {
      if (copy->arguments == NULL) {
	precPtr = (mp_prec_t *) safeMalloc(sizeof(mp_prec_t));
	*precPtr = sollya_mpfi_get_prec(aPrioriBoundForConstantExpr);
	intervalPtr = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
	sollya_mpfi_init2(*intervalPtr, sollya_mpfi_get_prec(aPrioriBoundForConstantExpr));
	sollya_mpfi_set(*intervalPtr, aPrioriBoundForConstantExpr);
	copy->arguments = addElement(addElement(NULL, intervalPtr), precPtr);
      }
    }
    sollya_mpfi_clear(aPrioriBoundForConstantExpr);
    haveAPrioriBoundForConstantExpr = 0;
  }

  return copy;
}

void composePolynomialsInner(sollya_mpfi_t *res, int degR, sollya_mpfi_t *p, int degP, sollya_mpfi_t *q, int degQ, mp_prec_t prec) {
  sollya_mpfi_t *r, *s, *t;
  int i, j, k, l;
  sollya_mpfi_t temp;

  /* Initialize a temporary
   */
  sollya_mpfi_init2(temp, prec);

  /* Allocate two scratch arrays of the length of
  // the result polynomial.
  */
  r = safeCalloc(degR+1,sizeof(sollya_mpfi_t));
  s = safeCalloc(degR+1,sizeof(sollya_mpfi_t));
  for (i=0;i<=degR;i++) {
    sollya_mpfi_init2(r[i],prec);
    sollya_mpfi_init2(s[i],prec);
    sollya_mpfi_set_si(r[i],0);
    sollya_mpfi_set_si(s[i],0);
  }

  /* Perform a stupid Horner:
  //
  // p(q(x)) = a_0 + q(x) * (.... a_i * q(x) * pTilde(q(x)) ...)
  */
  sollya_mpfi_set(r[0],p[degP]);
  for (i=degP-1;i>=0;i--) {
    /* Multiply pTilde(q) with q, result in s
       while clearing up r
    */
    for (j=0;j<=degQ*(degP-1-i);j++) {
      for (k=0;k<=degQ;k++) {
	l = j + k;
	sollya_mpfi_mul(temp,r[j],q[k]);
	sollya_mpfi_add(s[l],s[l],temp);
      }
      sollya_mpfi_set_si(r[j],0);
    }
    /* Add in a_i */
    sollya_mpfi_add(s[0],s[0],p[i]);

    /* Swap r and s for next step */
    t = s;
    s = r;
    r = t;
  }

  /* Here, the result is in r
  //
  // Copy it into the result array res
  */
  for (i=0;i<=degR;i++) {
    sollya_mpfi_set(res[i],r[i]);
  }

  /* Clear and free the scratch arrays
   */
  for (i=0;i<=degR;i++) {
    sollya_mpfi_clear(r[i]);
    sollya_mpfi_clear(s[i]);
  }
  safeFree(r);
  safeFree(s);

  /* Clear the temporary
   */
  sollya_mpfi_clear(temp);

}

void composePolynomials(node **poly, chain **radii, node *p, node *q, mp_prec_t prec) {
  int degP, degQ, i, degR;
  node **coeffsP, **coeffsQ;
  sollya_mpfi_t *polyP, *polyQ, *polyR;
  mpfr_t *coeffsR;
  sollya_mpfi_t **radiiArray;

  /* Check if both p and q are polynomials
  // and implement a fallback solution
  // for when this is not the case
  */
  if ((!isPolynomial(p)) || (!isPolynomial(q))) {
    *radii = NULL;
    *poly = substitute(p,q);
    return;
  }

  /* Here, we know that p and q are polynomials
  //
  // Continue by getting all coefficients of p and q
  */
  getCoefficients(&degP, &coeffsP, p);
  getCoefficients(&degQ, &coeffsQ, q);

  /* Evaluate all coefficients of p and q to
  // small intervals in two arrays polyP and polyQ.
  */
  polyP = (sollya_mpfi_t *) safeCalloc(degP+1,sizeof(sollya_mpfi_t));
  for (i=0;i<=degP;i++) {
    sollya_mpfi_init2(polyP[i],prec);
    if (coeffsP[i] != NULL) {
      evaluateConstantExpressionToSharpInterval(polyP[i], coeffsP[i]);
    } else {
      sollya_mpfi_set_si(polyP[i],0);
    }
  }
  polyQ = (sollya_mpfi_t *) safeCalloc(degQ+1,sizeof(sollya_mpfi_t));
  for (i=0;i<=degQ;i++) {
    sollya_mpfi_init2(polyQ[i],prec);
    if (coeffsQ[i] != NULL) {
      evaluateConstantExpressionToSharpInterval(polyQ[i], coeffsQ[i]);
    } else {
      sollya_mpfi_set_si(polyQ[i],0);
    }
  }

  /* Free the arrays with the unevaluated coefficients
   */
  for (i=0;i<=degP;i++) {
    if (coeffsP[i] != NULL) {
      free_memory(coeffsP[i]);
    }
  }
  safeFree(coeffsP);
  for (i=0;i<=degQ;i++) {
    if (coeffsQ[i] != NULL) {
      free_memory(coeffsQ[i]);
    }
  }
  safeFree(coeffsQ);

  /* Allocate and initialize an array of interval coefficients for the result
  // polynomial
  */
  degR = degP * degQ;
  polyR = (sollya_mpfi_t *) safeCalloc(degR+1,sizeof(sollya_mpfi_t));
  for (i=0;i<=degR;i++) {
    sollya_mpfi_init2(polyR[i],prec);
  }

  /* Have an auxiliary function do the real work of composing p and q
   */
  composePolynomialsInner(polyR,degR,polyP,degP,polyQ,degQ,prec);

  /* Free the arrays with the evaluated coefficients
   */
  for (i=0;i<=degP;i++) {
    sollya_mpfi_clear(polyP[i]);
  }
  safeFree(polyP);
  for (i=0;i<=degQ;i++) {
    sollya_mpfi_clear(polyQ[i]);
  }
  safeFree(polyQ);

  /* Allocate and compute an array of centerpoints
  //
  // Allocate also pointers to intervals and
  // initialize intervals that will hold the radii.
  //
  // Clear the intervals for the output coefficients
  // as well.
  */
  coeffsR = (mpfr_t *) safeCalloc(degR+1,sizeof(mpfr_t));
  radiiArray = (sollya_mpfi_t **) safeCalloc(degR+1,sizeof(sollya_mpfi_t *));
  for (i=0;i<=degR;i++) {
    mpfr_init2(coeffsR[i],prec+1);
    sollya_mpfi_mid(coeffsR[i],polyR[i]);
    radiiArray[i] = (sollya_mpfi_t *) safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*(radiiArray[i]),prec);
    sollya_mpfi_sub_fr(*(radiiArray[i]),polyR[i],coeffsR[i]);
    sollya_mpfi_clear(polyR[i]);
  }

  /* Free the array of intervals for the output coefficients
   */
  safeFree(polyR);

  /* Convert the array of centerpoints to a tree.
   */
  *poly = makePolynomial(coeffsR, degR);

  /* Free the arrays with the centerpoints
   */
  for (i=0;i<=degR;i++) {
    mpfr_clear(coeffsR[i]);
  }
  safeFree(coeffsR);

  /* Convert the array of radii to a list of radii
   */
  *radii = NULL;
  for (i=0;i<=degR;i++) {
    *radii = addElement(*radii,(void *) (radiiArray[i]));
  }

  /* Free the array holding the pointers to the
  // intervals of radii.
  */
  safeFree(radiiArray);

}

int readHexadecimal(mpfr_t rop, char *c) {
  mpfr_t vrd, vru;
  mp_prec_t p;
  int res, resA, resB;
  char *c2;

  c2 = (char *) safeCalloc(strlen(c) + 2, sizeof(char));
  strcpy(c2, c);

  if ((c2[strlen(c2) - 1] == 'p') ||
      (c2[strlen(c2) - 1] == 'P')) c2[strlen(c2)] = '0';


  p = mpfr_get_prec(rop);

  mpfr_init2(vrd, p);
  mpfr_init2(vru, p);

  resA = mpfr_set_str(vrd, c2, 16, GMP_RNDD);
  resB = mpfr_set_str(vru, c2, 16, GMP_RNDU);

  if (!resA && !resB) {
    if (mpfr_cmp(vrd,vru) == 0) {
      mpfr_set(rop,vrd,GMP_RNDN);
      res = 1;
    } else {
      resA = mpfr_set_str(vrd, c2, 16, GMP_RNDN);
      if (!resA) {
	mpfr_set(rop, vrd, GMP_RNDN);
	res = 0;
      } else {
	mpfr_set_nan(rop);
	res = 0;
      }
    }
  } else {
    mpfr_set_nan(rop);
    res = 0;
  }

  mpfr_clear(vrd);
  mpfr_clear(vru);
  safeFree(c2);

  return res;
}


int readDyadic(mpfr_t res, char *c) {
  char *mantissa, *exponent, *curr, *curr2;
  mpfr_t mant, expo, temp1, temp2;
  mp_prec_t prec;
  int rounding;

  mantissa = (char *) safeCalloc(strlen(c)+1,sizeof(char));
  exponent = (char *) safeCalloc(strlen(c)+1,sizeof(char));
  curr = c; curr2 = mantissa;
  while ((*curr != '\0') && (*curr != 'b') && (*curr != 'B')) {
    *curr2 = *curr;
    curr2++;
    curr++;
  }
  if (*curr != '\0') curr++;
  curr2 = exponent;
  while (*curr != '\0') {
    *curr2 = *curr;
    curr2++;
    curr++;
  }

  rounding = 1;

  prec = mpfr_get_prec(res);
  mpfr_init2(mant,prec);
  mpfr_init2(expo,prec);
  mpfr_init2(temp1,prec);
  mpfr_init2(temp2,prec);

  mpfr_set_str(temp1,mantissa,10,GMP_RNDU);
  mpfr_set_str(temp2,mantissa,10,GMP_RNDD);
  if (mpfr_cmp(temp1,temp2) != 0) {
    rounding = 0;
    mpfr_set_str(temp1,mantissa,10,GMP_RNDN);
  }
  if (mpfr_set(mant,temp1,GMP_RNDN) != 0) {
    rounding = 0;
  }
  mpfr_set_str(temp1,exponent,10,GMP_RNDU);
  mpfr_set_str(temp2,exponent,10,GMP_RNDD);
  if (mpfr_cmp(temp1,temp2) != 0) {
    rounding = 0;
    mpfr_set_str(temp1,exponent,10,GMP_RNDN);
  }
  if (mpfr_exp2(expo,temp1,GMP_RNDN) != 0) {
    rounding = 0;
  }
  if (mpfr_mul(res,mant,expo,GMP_RNDN) != 0) {
    rounding = 0;
  }

  if (!mpfr_number_p(res)) rounding = 1;

  mpfr_clear(mant);
  mpfr_clear(expo);
  mpfr_clear(temp1);
  mpfr_clear(temp2);
  safeFree(mantissa);
  safeFree(exponent);
  return rounding;
}

node *makePolynomialConstantExpressions(node **coeffs, int deg) {
  node *copy;
  int i, degree, e, k;
  node *temp;
  node *temp2;
  node *temp3;
  mpfr_t *value;
  node *temp4;
  node *res;
  unsigned int degU;
  polynomial_t p;

  if (deg < 0) {
    sollyaFprintf(stderr,"Error: makePolynomialConstantExpressions: degree of polynomial to be built is negative\n");
    exit(1);
  }

  res = addMemRefEvenOnNull(NULL);
  if (res != NULL) {
    degU = deg;
    if (polynomialFromConstantExpressionCoefficients(&p, coeffs, degU)) {
      res->cache->polynomialRepresentation = p;
      return res;
    } else {
      res->cache->polynomialRepresentation = polynomialFromIntConstant(1);
      free_memory(res);
    }
  }

  degree = deg;
  while ((degree > 0) &&
         ((coeffs[degree] == NULL) ||
          ((accessThruMemRef(coeffs[degree])->nodeType == CONSTANT) &&
           (mpfr_zero_p(*(accessThruMemRef(coeffs[degree])->value)))))) degree--;

  if (degree == 0) {
    if (coeffs[0] == NULL) return makeConstantDouble(0.0);
    return copyTree(coeffs[0]);
  }

  copy = copyTree(coeffs[degree]);
  for (i=degree-1;i>=0;i--) {
    if ((coeffs[i] == NULL) ||
        ((accessThruMemRef(coeffs[i])->nodeType == CONSTANT) &&
         (mpfr_zero_p(*(accessThruMemRef(coeffs[i])->value))))) {
      if (i == 0) {
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp->child1 = temp2;
	temp->child2 = copy;
	copy = temp;
      } else {
	for (k=i-1;(((coeffs[k] == NULL) ||
                     ((accessThruMemRef(coeffs[k])->nodeType == CONSTANT) &&
                      (mpfr_zero_p(*(accessThruMemRef(coeffs[k])->value))))) && (k > 0));k--);
	e = (i - k) + 1;
	temp = allocateNode();
	temp->nodeType = MUL;
	temp2 = makeVariable();
	temp3 = allocateNode();
	temp3->nodeType = CONSTANT;
	value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*value,(mp_prec_t) ((mp_prec_t) tools_precision > (mp_prec_t) (8 * sizeof(e) + 10) ? (mp_prec_t) tools_precision : (mp_prec_t) (8 * sizeof(e) + 10)));
	if (mpfr_set_si(*value,e,GMP_RNDN) != 0) {
	  if (!noRoundingWarnings) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: rounding occurred on representing a monomial power exponent with %d bits.\n",
			 (int) mpfr_get_prec(*value));
	    printMessage(1,SOLLYA_MSG_CONTINUATION,"Try to increase the precision.\n");
	  }
	}
	temp3->value = value;
	temp4 = allocateNode();
	temp4->nodeType = POW;
	temp4->child1 = temp2;
	temp4->child2 = temp3;
	temp->child1 = temp4;
	temp->child2 = copy;
	copy = temp;
	if (!((coeffs[k] == NULL) ||
              ((accessThruMemRef(coeffs[k])->nodeType == CONSTANT) &&
               (mpfr_zero_p(*(accessThruMemRef(coeffs[k])->value)))))) {
	  temp = allocateNode();
	  temp->nodeType = ADD;
	  temp->child1 = copyTree(coeffs[k]);
	  temp->child2 = copy;
	  copy = temp;
	}
	i = k;
      }
    } else {
      temp = allocateNode();
      temp->nodeType = MUL;
      temp2 = makeVariable();
      temp->child1 = temp2;
      temp->child2 = copy;
      copy = temp;
      temp = allocateNode();
      temp->nodeType = ADD;
      temp->child1 = copyTree(coeffs[i]);
      temp->child2 = copy;
      copy = temp;
    }
  }
  return addMemRef(copy);
}

/* Builds a polynomial expression, written in Horner form of the polynomial
   sum_{i=0}^degree  coefficients[i]*x^i
   Moreover, while writting it in Horner form, it accounts for the sparsity of
   the polynomial.
*/
node *makePolynomial(mpfr_t *coefficients, int degree) {
  node **coeffs;
  int i;
  node *poly;

  /* Allocate and build an array of constant expressions
     representing the coefficients.
     Do not represent zeros.
  */
  coeffs = (node **) safeCalloc(degree+1,sizeof(node *));
  for (i=0;i<=degree;i++) {
    if (!mpfr_zero_p(coefficients[i])) {
      coeffs[i] = makeConstant(coefficients[i]);
    } else {
      /* The coefficient is zero. Do not represent it. */
      coeffs[i] = NULL;
    }
  }

  /* Build the polynomial */
  poly = makePolynomialConstantExpressions(coeffs, degree);

  /* Free the array of constant expressions representing the
     coefficients
  */
  for (i=0;i<=degree;i++) {
    if (coeffs[i] != NULL) free_memory(coeffs[i]);
  }
  safeFree(coeffs);

  /* Return the polynomial */
  return poly;
}


int treeSize(node *tree) {
  int size;

  if (tree == NULL) return 0;
  switch (tree->nodeType) {
  case MEMREF:
    if (tree->cache->treeSizeCacheFilled) {
      return tree->cache->treeSizeCache;
    }

    size = treeSize(getMemRefChild(tree));

    tree->cache->treeSizeCache = size;
    tree->cache->treeSizeCacheFilled = 1;

    return size;
    break;

  case VARIABLE:
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    return 1;
    break;

  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    return treeSize(tree->child1) + treeSize(tree->child2) + 1;
    break;

  case NEG:
  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    return treeSize(tree->child1) + 1;
    break;

  default:
    sollyaFprintf(stderr,"Error: treeSize: unknown identifier (%d) in the tree\n",tree->nodeType);
    exit(1);
  }
  return -1;
}


int highestDegreeOfPolynomialSubexpression(node *tree) {
  int l, r;

  if (tree->nodeType == MEMREF) return highestDegreeOfPolynomialSubexpression(getMemRefChild(tree));

  if (isPolynomial(tree)) return getDegree(tree);

  switch (arity(tree)) {
  case 2:
    l = highestDegreeOfPolynomialSubexpression(tree->child1);
    r = highestDegreeOfPolynomialSubexpression(tree->child2);
    return l > r ? l : r;
    break;
  case 1:
    return highestDegreeOfPolynomialSubexpression(tree->child1);
    break;
  case 0:
    return getDegree(tree);
    break;
  default:
    sollyaFprintf(stderr,"Error: unknown arity of an operator.\n");
    exit(1);
  }

  return -1;
}

node *makeAddSimplified(node *a, node *b) {

  if ((accessThruMemRef(a)->nodeType == CONSTANT) &&
      mpfr_zero_p(*(accessThruMemRef(a)->value))) {
    free_memory(a);
    return b;
  }

  if ((accessThruMemRef(b)->nodeType == CONSTANT) &&
      mpfr_zero_p(*(accessThruMemRef(b)->value))) {
    free_memory(b);
    return a;
  }

  return makeAdd(a, b);
}

node *makeSubSimplified(node *a, node *b) {

  if ((accessThruMemRef(a)->nodeType == CONSTANT) &&
      mpfr_zero_p(*(accessThruMemRef(a)->value))) {
    free_memory(a);
    return b;
  }

  if ((accessThruMemRef(b)->nodeType == CONSTANT) &&
      mpfr_zero_p(*(accessThruMemRef(b)->value))) {
    free_memory(b);
    return a;
  }

  return makeSub(a, b);
}

node *makeMulSimplified(node *a, node *b) {

  if ((accessThruMemRef(a)->nodeType == CONSTANT) &&
      mpfr_zero_p(*(accessThruMemRef(a)->value))) {
    free_memory(b);
    return a;
  }

  if ((accessThruMemRef(b)->nodeType == CONSTANT) &&
      mpfr_zero_p(*(accessThruMemRef(b)->value))) {
    free_memory(a);
    return b;
  }

  if ((accessThruMemRef(a)->nodeType == CONSTANT) &&
      mpfr_number_p(*(accessThruMemRef(a)->value)) &&
      (mpfr_cmp_si(*(accessThruMemRef(a)->value),1) == 0)) {
    free_memory(a);
    return b;
  }

  if ((accessThruMemRef(b)->nodeType == CONSTANT) &&
      mpfr_number_p(*(accessThruMemRef(b)->value)) &&
      (mpfr_cmp_si(*(accessThruMemRef(b)->value),1) == 0)) {
    free_memory(b);
    return a;
  }

  return makeMul(a, b);
}

int tryGetIthCoefficientSparseUnsafe(node **res, node *poly, int i) {
  node *resLeft, *resRight;
  mpfr_t iAsMpfr, tAsMpfr;
  mpz_t tMpz, kMpz;
  int t, k;

  resLeft = NULL;
  resRight = NULL;

  if (i < 0) return 0;

  if (isConstant(poly)) {
    if (i == 0) {
      *res = copyTree(poly);
    } else {
      *res = makeConstantDouble(0.0);
    }
    return 1;
  }

  switch (poly->nodeType) {
  case MEMREF:
    return tryGetIthCoefficientSparseUnsafe(res, getMemRefChild(poly), i);
    break;
  case CONSTANT:
  case PI_CONST:
  case LIBRARYCONSTANT:
    /* p = c */
    if (i == 0) {
      *res = copyTree(poly);
    } else {
      *res = makeConstantDouble(0.0);
    }
    return 1;
    break;
  case VARIABLE:
    /* p = x */
    if (i == 1) {
      *res = makeConstantDouble(1.0);
    } else {
      *res = makeConstantDouble(0.0);
    }
    return 1;
    break;
  case ADD:
    /* p = q + r */
    if (tryGetIthCoefficientSparseUnsafe(&resLeft, poly->child1, i) &&
	tryGetIthCoefficientSparseUnsafe(&resRight, poly->child2, i)) {
      *res = makeAddSimplified(resLeft, resRight);
      resLeft = NULL;
      resRight = NULL;
      return 1;
    }
    break;
  case SUB:
    /* p = q - r */
    if (tryGetIthCoefficientSparseUnsafe(&resLeft, poly->child1, i) &&
	tryGetIthCoefficientSparseUnsafe(&resRight, poly->child2, i)) {
      *res = makeSubSimplified(resLeft, resRight);
      resLeft = NULL;
      resRight = NULL;
      return 1;
    }
    break;
  case MUL:
    /* p = c * q */
    if (isConstant(poly->child1)) {
      if (tryGetIthCoefficientSparseUnsafe(&resRight, poly->child2, i)) {
	*res = makeMulSimplified(copyTree(poly->child1), resRight);
	resRight = NULL;
	return 1;
      }
    }
    /* p = q * c */
    if (isConstant(poly->child2)) {
      if (tryGetIthCoefficientSparseUnsafe(&resLeft, poly->child1, i)) {
	*res = makeMulSimplified(resLeft, copyTree(poly->child2));
	resLeft = NULL;
	return 1;
      }
    }
    /* p = x * q */
    if (accessThruMemRef(poly->child1)->nodeType == VARIABLE) {
      if (tryGetIthCoefficientSparseUnsafe(res, poly->child2, i - 1)) return 1;
    }
    /* p = q * x */
    if (accessThruMemRef(poly->child2)->nodeType == VARIABLE) {
      if (tryGetIthCoefficientSparseUnsafe(res, poly->child1, i - 1)) return 1;
    }
    /* p = x^t * q */
    if ((accessThruMemRef(poly->child1)->nodeType == POW) &&
	(accessThruMemRef(accessThruMemRef(poly->child1)->child1)->nodeType == VARIABLE) &&
	((accessThruMemRef(accessThruMemRef(poly->child1)->child2)->nodeType == CONSTANT) &&
	 (mpfr_number_p(*(accessThruMemRef(accessThruMemRef(poly->child1)->child2)->value))))) {
      t = mpfr_get_si(*(accessThruMemRef(accessThruMemRef(poly->child1)->child2)->value),GMP_RNDN);
      mpfr_init2(tAsMpfr, 8 * sizeof(int) + 10);
      mpfr_set_si(tAsMpfr, t, GMP_RNDN); /* exact */
      if (mpfr_cmp(*(accessThruMemRef(accessThruMemRef(poly->child1)->child2)->value), tAsMpfr) == 0) {
	if (t <= i) {
	  if (tryGetIthCoefficientSparseUnsafe(res, poly->child2, i - t)) {
	    mpfr_clear(tAsMpfr);
	    return 1;
	  }
	} else {
	  mpfr_clear(tAsMpfr);
	  *res = makeConstantDouble(0.0);
	  return 1;
	}
      }
      mpfr_clear(tAsMpfr);
    }
    /* p = q * x^t */
    if ((accessThruMemRef(poly->child2)->nodeType == POW) &&
	(accessThruMemRef(accessThruMemRef(poly->child2)->child1)->nodeType == VARIABLE) &&
	((accessThruMemRef(accessThruMemRef(poly->child2)->child2)->nodeType == CONSTANT) &&
	 (mpfr_number_p(*(accessThruMemRef(accessThruMemRef(poly->child2)->child2)->value))))) {
      t = mpfr_get_si(*(accessThruMemRef(accessThruMemRef(poly->child2)->child2)->value),GMP_RNDN);
      mpfr_init2(tAsMpfr, 8 * sizeof(int) + 10);
      mpfr_set_si(tAsMpfr, t, GMP_RNDN); /* exact */
      if (mpfr_cmp(*(accessThruMemRef(accessThruMemRef(poly->child2)->child2)->value), tAsMpfr) == 0) {
	if (t <= i) {
	  if (tryGetIthCoefficientSparseUnsafe(res, poly->child1, i - t)) {
	    mpfr_clear(tAsMpfr);
	    return 1;
	  }
	} else {
	  mpfr_clear(tAsMpfr);
	  *res = makeConstantDouble(0.0);
	  return 1;
	}
      }
      mpfr_clear(tAsMpfr);
    }
    /* p = q * r, degree(q) = t, degree(r) = k, t * k = i */
    mpz_init(tMpz);
    if (getDegreeMpz(tMpz, poly->child1)) {
      t = mpz_get_si(tMpz);
      if (mpz_cmp_si(tMpz, t) == 0) {
	mpz_init(kMpz);
	if (getDegreeMpz(kMpz, poly->child2)) {
	  k = mpz_get_si(kMpz);
	  if (mpz_cmp_si(kMpz, k) == 0) {
	    if ((t * k == i) &&
                (t >= 0) && (t <= i) &&
		(k >= 0) && (k <= i) &&
		((((k == 0) || (t == 0)) && (i == 0)) || (t == i / k))
		) {
	      if (tryGetIthCoefficientSparseUnsafe(&resLeft, poly->child1, t) &&
		  tryGetIthCoefficientSparseUnsafe(&resRight, poly->child2, k)) {
		*res = makeMulSimplified(resLeft, resRight);
		resLeft = NULL;
		resRight = NULL;
		return 1;
	      }
	    }
	  }
	}
	mpz_clear(kMpz);
      }
    }
    mpz_clear(tMpz);
    /* Continue with other optimized ways to get the i-th coefficient
       of a sparse polynomial that is a product here.
    */
    break;
  case DIV:
    /* p = q / c */
    if (isConstant(poly->child2)) {
      if (tryGetIthCoefficientSparseUnsafe(&resLeft, poly->child1, i)) {
	*res = makeDiv(resLeft, copyTree(poly->child2));
	resLeft = NULL;
	return 1;
      }
    }
    break;
  case POW:
    /* p = x^k */
    if ((accessThruMemRef(poly->child1)->nodeType == VARIABLE) &&
	(accessThruMemRef(poly->child2)->nodeType == CONSTANT) &&
	mpfr_number_p(*(accessThruMemRef(poly->child2)->value))) {
      mpfr_init2(iAsMpfr, 8 * sizeof(int) + 10);
      mpfr_set_si(iAsMpfr, i, GMP_RNDN); /* exact */
      if (mpfr_cmp(*(accessThruMemRef(poly->child2)->value), iAsMpfr) == 0) {
	*res = makeConstantDouble(1.0);
	mpfr_clear(iAsMpfr);
	return 1;
      } else {
	if (mpfr_integer_p(*(accessThruMemRef(poly->child2)->value))) {
	  *res = makeConstantDouble(0.0);
	  mpfr_clear(iAsMpfr);
	  return 1;
	}
      }
      mpfr_clear(iAsMpfr);
    }
    /* p = q^t with degree(q) = k and k * t = i */
    if ((accessThruMemRef(poly->child2)->nodeType == CONSTANT) &&
	mpfr_number_p(*(accessThruMemRef(poly->child2)->value))) {
      t = mpfr_get_si(*(accessThruMemRef(poly->child2)->value),GMP_RNDN);
      mpfr_init2(tAsMpfr, 8 * sizeof(int) + 10);
      mpfr_set_si(tAsMpfr, t, GMP_RNDN); /* exact */
      if (mpfr_cmp(*(accessThruMemRef(poly->child2)->value), tAsMpfr) == 0) {
	if ((t != 0) && (i % t == 0)) {
	  k = getDegreeSilent(poly->child1);
	  if ((k > 0) && (k <= t) && (k * t == i)) {
	    if (tryGetIthCoefficientSparseUnsafe(res, poly->child1, k)) {
	      mpfr_clear(tAsMpfr);
	      return 1;
	    }
	  }
	}
      }
      mpfr_clear(tAsMpfr);
    }
    /* Continue with other optimized ways to get the i-th coefficient
       of a sparse polynomial that is a power here.
    */
    break;
  default:
    return 0;
  }

  if (resLeft != NULL) free_memory(resLeft);
  if (resRight != NULL) free_memory(resRight);

  return 0;
}

int tryGetIthCoefficientSparse(node **res, node *poly, int i) {
  node *myres;

  if (!isPolynomial(poly)) return 0;

  myres = NULL;
  if (tryGetIthCoefficientSparseUnsafe(&myres, poly, i)) {
    if (myres == NULL) return 0;
    *res = simplifyTreeErrorfree(myres);
    free_memory(myres);
    return 1;
  }

  return 0;
}

node *getIthCoefficient(node *poly, int i) {
  node *tempNode;
  node **coefficients;
  int degree, k;

  if (poly->nodeType == MEMREF) {
    if (poly->cache->polynomialRepresentation == NULL) {
      tryRepresentAsPolynomial(poly);
    }
    if (poly->cache->polynomialRepresentation != NULL) {
      return polynomialGetIthCoefficientIntIndex(poly->cache->polynomialRepresentation, i);
    }
  }

  if ((!isPolynomial(poly)) || (i < 0)) {
    tempNode = allocateNode();
    tempNode->nodeType = CONSTANT;
    tempNode->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*(tempNode->value),10);
    mpfr_set_d(*(tempNode->value),0.0,GMP_RNDN);
    return tempNode;
  }

  if (tryGetIthCoefficientSparse(&tempNode, poly, i)) {
    printMessage(8, SOLLYA_MSG_SPECIAL_ALGORITHM_USED_FOR_COEFF, "Information: a special algorithm is used to extract the i-th coefficient of an expression that is a polynomial.\n");
    return tempNode;
  }

  getCoefficients(&degree, &coefficients, poly);

  if ((i > degree) || (coefficients[i] == NULL)) {
    tempNode = allocateNode();
    tempNode->nodeType = CONSTANT;
    tempNode->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*(tempNode->value),10);
    mpfr_set_d(*(tempNode->value),0.0,GMP_RNDN);
  } else {
    tempNode = copyTree(coefficients[i]);
  }

  for (k=0;k<=degree;k++) {
    if (coefficients[k] != NULL) free_memory(coefficients[k]);
  }

  safeFree(coefficients);

  return tempNode;
}


node *getSubpolynomial(node *poly, chain *monomials, int fillDegrees, mp_prec_t prec) {
  node *tempNode, *tempNode2, *tempNode3;
  node **coefficients;
  int degree, k, currDeg, maxDegree;
  chain *curr;

  tempNode = allocateNode();
  tempNode->nodeType = CONSTANT;
  tempNode->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(tempNode->value),prec);
  mpfr_set_d(*(tempNode->value),0.0,GMP_RNDN);

  if (!isPolynomial(poly)) {
    return tempNode;
  }

  getCoefficients(&degree, &coefficients, poly);

  curr = monomials;

  maxDegree = -1;

  while (curr != NULL) {
    currDeg = *((int *) (curr->value));
    if (currDeg > maxDegree) maxDegree = currDeg;
    if ((currDeg >= 0) && (currDeg <= degree) && (coefficients[currDeg] != NULL)) {
      tempNode2 = allocateNode();
      tempNode2->nodeType = POW;
      tempNode3 = allocateNode();
      tempNode3->nodeType = CONSTANT;
      tempNode3->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(tempNode3->value),prec);
      if (mpfr_set_si(*(tempNode3->value),currDeg,GMP_RNDN) != 0) {
	printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: during subpolynomial extraction, the exponent of a power could not be represented exactly on with the given precision.\n");
      }
      tempNode2->child2 = tempNode3;
      tempNode3 = makeVariable();
      tempNode2->child1 = tempNode3;
      tempNode3 = allocateNode();
      tempNode3->nodeType = MUL;
      tempNode3->child2 = tempNode2;
      tempNode3->child1 = copyTree(coefficients[currDeg]);
      tempNode2 = allocateNode();
      tempNode2->nodeType = ADD;
      tempNode2->child2 = tempNode3;
      tempNode2->child1 = tempNode;
      tempNode = tempNode2;
    }
    curr = curr->next;
  }

  if (fillDegrees) {
    for (k=maxDegree+1;k<=degree;k++) {
      if (coefficients[k] != NULL) {
	tempNode2 = allocateNode();
	tempNode2->nodeType = POW;
	tempNode3 = allocateNode();
	tempNode3->nodeType = CONSTANT;
	tempNode3->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	mpfr_init2(*(tempNode3->value),prec);
	if (mpfr_set_si(*(tempNode3->value),k,GMP_RNDN) != 0) {
	  printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: during subpolynomial extraction, the exponent of a power could not be represented exactly on with the given precision.\n");
	}
	tempNode2->child2 = tempNode3;
	tempNode3 = makeVariable();
	tempNode2->child1 = tempNode3;
	tempNode3 = allocateNode();
	tempNode3->nodeType = MUL;
	tempNode3->child2 = tempNode2;
	tempNode3->child1 = copyTree(coefficients[k]);
	tempNode2 = allocateNode();
	tempNode2->nodeType = ADD;
	tempNode2->child2 = tempNode3;
	tempNode2->child1 = tempNode;
	tempNode = tempNode2;
      }
    }
  }


  for (k=0;k<=degree;k++) {
    if (coefficients[k] != NULL) free_memory(coefficients[k]);
  }

  safeFree(coefficients);

  tempNode2 = horner(tempNode);

  free_memory(tempNode);

  return tempNode2;
}

node *makeCanonicalPolyUnsafe(node *poly, mp_prec_t prec) {
  node **coefficients;
  int degree, k;
  node *tempNode, *tempNode2, *tempNode3;

  getCoefficients(&degree, &coefficients, poly);

  tempNode = allocateNode();
  tempNode->nodeType = CONSTANT;
  tempNode->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(tempNode->value),prec);
  mpfr_set_d(*(tempNode->value),0.0,GMP_RNDN);
  for (k=0;k<=degree;k++) {
    if (coefficients[k] != NULL) {
      if (k == 0) {
	tempNode2 = allocateNode();
        tempNode2->nodeType = ADD;
        tempNode2->child2 = coefficients[k];
        tempNode2->child1 = tempNode;
        tempNode = tempNode2;
      } else {
	if (k == 1) {
	  tempNode3 = makeVariable();
	  tempNode2 = allocateNode();
	  tempNode2->nodeType = MUL;
	  tempNode2->child2 = tempNode3;
	  tempNode2->child1 = coefficients[k];
	  tempNode3 = allocateNode();
	  tempNode3->nodeType = ADD;
	  tempNode3->child2 = tempNode2;
	  tempNode3->child1 = tempNode;
	  tempNode = tempNode3;
	} else {
	  tempNode2 = allocateNode();
	  tempNode2->nodeType = POW;
	  tempNode3 = allocateNode();
	  tempNode3->nodeType = CONSTANT;
	  tempNode3->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
	  mpfr_init2(*(tempNode3->value),prec);
	  if (mpfr_set_si(*(tempNode3->value),k,GMP_RNDN) != 0) {
	    printMessage(1,SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION,"Warning: during transformation to canonical form, the exponent of a power could not be represented exactly on with the given precision.\n");
	  }
	  tempNode2->child2 = tempNode3;
	  tempNode3 = makeVariable();
	  tempNode2->child1 = tempNode3;
	  tempNode3 = allocateNode();
	  tempNode3->nodeType = MUL;
	  tempNode3->child2 = tempNode2;
	  tempNode3->child1 = coefficients[k];
	  tempNode2 = allocateNode();
	  tempNode2->nodeType = ADD;
	  tempNode2->child2 = tempNode3;
	  tempNode2->child1 = tempNode;
	  tempNode = tempNode2;
	}
      }
    }
  }

  safeFree(coefficients);

  tempNode2 = simplifyTreeErrorfree(tempNode);

  free_memory(tempNode);

  return tempNode2;
}


int isCanonicalMonomial(node *tree) {

  if (tree->nodeType == MEMREF) return isCanonicalMonomial(getMemRefChild(tree));

  if (isConstant(tree)) return 1;

  if (isPowerOfVariable(tree)) return 1;

  if ((tree->nodeType == MUL) &&
      isConstant(tree->child1) &&
      isPowerOfVariable(tree->child2)) return 1;

  return 0;
}


int isCanonicalUnsafe(node *tree) {
  int deg1, deg2;

  if (tree->nodeType == MEMREF) {
    if (((tree->child1 == NULL) || (tree->cache->memRefChildFromPolynomial)) && (tree->cache->polynomialRepresentation != NULL)) {
      return polynomialIsCanonicalized(tree->cache->polynomialRepresentation);
    }
    return isCanonicalUnsafe(getMemRefChild(tree));
  }

  if (isConstant(tree) || isCanonicalMonomial(tree)) return 1;

  if ((tree->nodeType == ADD) || (tree->nodeType == SUB)) {
    if (!isCanonicalUnsafe(tree->child1)) return 0;
    if (!isCanonicalMonomial(tree->child2)) return 0;
    deg1 = getDegree(tree->child1);
    deg2 = getDegree(tree->child2);
    if (deg1 >= deg2) return 0;
    return 1;
  }

  return 0;
}


int isCanonical(node *tree) {
  return isCanonicalUnsafe(tree);
}



node *makeCanonical(node *tree, mp_prec_t prec) {
  node *copy, *res;
  mpfr_t *value;
  mpfr_t temp;
  polynomial_t p;

  if (tree->nodeType == MEMREF) {
    if (tree->cache->polynomialRepresentation != NULL) {
      if (polynomialIsCanonicalized(tree->cache->polynomialRepresentation)) return copyTree(tree);
      if (tree->child1 == NULL) {
	p = polynomialCanonicalize(tree->cache->polynomialRepresentation);
	polynomialFree(tree->cache->polynomialRepresentation);
	tree->cache->polynomialRepresentation = p;
	return copyTree(tree);
      }
      res = addMemRefEvenOnNull(NULL);
      if (res != NULL) {
	res->cache->polynomialRepresentation = polynomialCanonicalize(tree->cache->polynomialRepresentation);
	copyTreeAnnotationsNoSimplifications(res, tree);
	return res;
      }
    }
    res = addMemRef(makeCanonical(getMemRefChild(tree), prec));
    copyTreeAnnotationsNoSimplifications(res, tree);
    return res;
  }

  if (isCanonical(tree)) {
    printMessage(7,SOLLYA_MSG_EXPR_NOT_CANONICALIZED_AS_ALREADY_CANONICAL,"Information: no canonical form simplification will be performed because the given tree is already canonical.\n");
    return copyTree(tree);
  }

  if (isPolynomial(tree)) return makeCanonicalPolyUnsafe(tree,prec);

  switch (tree->nodeType) {
  case VARIABLE:
    copy = makeVariable();
    break;
  case CONSTANT:
    copy = allocateNode();
    copy->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(temp,tools_precision);
    simplifyMpfrPrec(temp,*(tree->value));
    mpfr_init2(*value,mpfr_get_prec(temp));
    mpfr_set(*value,temp,GMP_RNDN);
    mpfr_clear(temp);
    copy->value = value;
    break;
  case ADD:
    copy = allocateNode();
    copy->nodeType = ADD;
    copy->child1 = makeCanonical(tree->child1,prec);
    copy->child2 = makeCanonical(tree->child2,prec);
    break;
  case SUB:
    copy = allocateNode();
    copy->nodeType = SUB;
    copy->child1 = makeCanonical(tree->child1,prec);
    copy->child2 = makeCanonical(tree->child2,prec);
    break;
  case MUL:
    copy = allocateNode();
    copy->nodeType = MUL;
    copy->child1 = makeCanonical(tree->child1,prec);
    copy->child2 = makeCanonical(tree->child2,prec);
    break;
  case DIV:
    copy = allocateNode();
    copy->nodeType = DIV;
    copy->child1 = makeCanonical(tree->child1,prec);
    copy->child2 = makeCanonical(tree->child2,prec);
    break;
  case NEG:
    copy = allocateNode();
    copy->nodeType = NEG;
    copy->child1 = makeCanonical(tree->child1,prec);
    break;
  case UNARY_BASE_FUNC:
    copy = allocateNode();
    copy->nodeType = UNARY_BASE_FUNC;
    copy->baseFun = tree->baseFun;
    copy->child1 = makeCanonical(tree->child1,prec);
    break;
  case POW:
    copy = allocateNode();
    copy->nodeType = POW;
    copy->child1 = makeCanonical(tree->child1,prec);
    copy->child2 = makeCanonical(tree->child2,prec);
    break;
  case LIBRARYFUNCTION:
    copy = allocateNode();
    copy->nodeType = LIBRARYFUNCTION;
    copy->libFun = tree->libFun;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = makeCanonical(tree->child1,prec);
    break;
  case PROCEDUREFUNCTION:
    copy = allocateNode();
    copy->nodeType = PROCEDUREFUNCTION;
    copy->libFunDeriv = tree->libFunDeriv;
    copy->child1 = makeCanonical(tree->child1,prec);
    copy->child2 = copyThing(tree->child2);
    break;
  case PI_CONST:
    copy = allocateNode();
    copy->nodeType = PI_CONST;
    break;
  case LIBRARYCONSTANT:
    copy = allocateNode();
    copy->nodeType = LIBRARYCONSTANT;
    copy->libFun = tree->libFun;
    break;
  default:
    sollyaFprintf(stderr,"Error: makeCanonical: unknown identifier in the tree\n");
    exit(1);
  }
  return copy;
}


static inline node *__makeVariable() {
  node *res;

  res = allocateNode();
  res->nodeType = VARIABLE;

  return addMemRef(res);
}

node *__makeVariableCache = NULL;
int __makingAVariable = 0;

void freeVariableCache() {
  if (__makeVariableCache != NULL) {
    free_memory(__makeVariableCache);
    __makeVariableCache = NULL;
  }
  __makingAVariable = 0;
}

void resetVariableCacheHandling() {
  __makingAVariable = 0;
}

node *makeVariable() {
  node *res;

  if (__makingAVariable) {
    return __makeVariable();
  }
  __makingAVariable = 1;

  if (__makeVariableCache != NULL) {
    return addMemRef(copyThing(__makeVariableCache));
  }

  res = __makeVariable();
  if ((__makeVariableCache == NULL) &&
      ((res != NULL) &&
       (res->nodeType == MEMREF))) {
    __makeVariableCache = res;
    res = addMemRef(copyThing(__makeVariableCache));
  }

  __makingAVariable = 0;

  return res;
}

node *makeConstantDouble(double d) {
  node *res;

  res = allocateNode();
  res->nodeType = CONSTANT;
  res->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res->value),53);
  mpfr_set_d(*(res->value),d,GMP_RNDN);

  return res;
}

node *makeConstantInt(int a) {
  node *res;

  res = allocateNode();
  res->nodeType = CONSTANT;
  res->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res->value),sizeof(int)*8);
  mpfr_set_si(*(res->value), a, GMP_RNDN);

  return res;
}

node *makeConstantUnsignedInt(unsigned int a) {
  node *res;

  res = allocateNode();
  res->nodeType = CONSTANT;
  res->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res->value),sizeof(unsigned int)*8);
  mpfr_set_ui(*(res->value), a, GMP_RNDN);

  return res;
}

node *makeConstantMpz(mpz_t x) {
  node *res;
  mp_prec_t prec;

  res = allocateNode();
  res->nodeType = CONSTANT;
  if (mpz_sgn(x) == 0) {
    prec = 12;
  } else {
    prec = (mp_prec_t) mpz_sizeinbase(x, 2);
    if (prec < 12) prec = 12;
  }
  res->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res->value),prec);
  mpfr_set_z(*(res->value),x,GMP_RNDN); /* exact */

  return res;
}

node *makeConstant(mpfr_t x) {
  node *res;

  res = allocateNode();
  res->nodeType = CONSTANT;
  res->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*(res->value),mpfr_get_prec(x));
  mpfr_set(*(res->value),x,GMP_RNDN);

  return res;
}

node *makeBinary(node *op1, node *op2, int opType) {
  node *res;

  res = allocateNode();
  res->nodeType = opType;
  res->child1 = op1;
  res->child2 = op2;

  return res;
}

node *makeUnary(node *op1, baseFunction const *f) {
  node *res;

  res = allocateNode();
  res->nodeType = UNARY_BASE_FUNC;
  res->baseFun = f;
  res->child1 = op1;

  return res;
}

node *makeAdd(node *op1, node *op2) {
  return makeBinary(op1,op2,ADD);
}

node *makeSub(node *op1, node *op2) {
  return makeBinary(op1,op2,SUB);
}

node *makeMul(node *op1, node *op2) {
  return makeBinary(op1,op2,MUL);
}

node *makeDiv(node *op1, node *op2) {
  return makeBinary(op1,op2,DIV);
}

node *makeNeg(node *op1) {
  node *res;

  res = allocateNode();
  res->nodeType = NEG;
  res->child1 = op1;

  return res;
}

node *makePow(node *op1, node *op2) {
  return makeBinary(op1,op2,POW);
}

node *makePi() {
  node *res;

  res = allocateNode();
  res->nodeType = PI_CONST;

  return res;
}

int readDecimalConstant(mpfr_t result, char *str) {
  mpfr_t a,b;
  int ternary;

  mpfr_init2(a,tools_precision);
  mpfr_init2(b,tools_precision);

  mpfr_set_str(a,str,10,GMP_RNDD);
  mpfr_set_str(b,str,10,GMP_RNDU);
  if (mpfr_cmp(a,b) != 0) {
    if (!noRoundingWarnings) {
      printMessage(1,SOLLYA_MSG_ROUNDING_OCCURRED_WHILE_READING_A_CONSTANT,
		   "Warning: Rounding occurred when converting the constant \"%s\" to floating-point with %d bits.\n",
		   str,(int) tools_precision);
      printMessage(1,SOLLYA_MSG_CONTINUATION,"If safe computation is needed, try to increase the precision.\n");
    }
    ternary = mpfr_set_str(a,str,10,GMP_RNDN);
  } else {
    ternary = 0;
  }

  mpfr_set_prec(result, tools_precision);
  mpfr_set(result,a,GMP_RNDN);

  mpfr_clear(a);
  mpfr_clear(b);

  return ternary;
}

