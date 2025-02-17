/*

  Copyright 2007-2017 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2),

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France

  and by

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  Boîte Courrier 169
  4, place Jussieu
  F-75252 Paris Cedex 05
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if (defined(HAVE_DLADDR) && (HAVE_DLADDR) && defined(HAVE_DL_INFO_TYPE) && (HAVE_DL_INFO_TYPE))
#define _GNU_SOURCE 1
#define _DARWIN_C_SOURCE 1
#endif

#include <dlfcn.h>


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <mpfr.h>
#include "mpfi-compat.h"
#include "expression.h"
#include "general.h"
#include "library.h"
#include "chain.h"
#include "execute.h"

#ifndef RTLD_LOCAL
#define RTLD_LOCAL 0
#endif

#define LIBRARYFUNCTION_CASE 0
#define LIBRARYCONSTANT_CASE 1
#define LIBRARYPROC_CASE 2

/* List of Sollya keywords we must not bind dynamic objects to.

   The last entry must be NULL.

*/
static const char* sollyaKeywords[] = { "D",
					"DD",
					"DE",
					"Pi",
					"RD",
					"RN",
					"RU",
					"RZ",
					"SG",
					"TD",
					"abs",
					"absolute",
					"accurateinfnorm",
					"acos",
					"acosh",
					"annotatefunction",
					"asciiplot",
					"asin",
					"asinh",
					"atan",
					"atanh",
					"autodiff",
					"autosimplify",
					"bashexecute",
					"begin",
					"binary",
					"bind",
					"boolean",
					"by",
					"canonical",
					"ceil",
					"chebyshevform",
					"checkinfnorm",
					"coeff",
					"composepolynomials",
					"constant",
					"cos",
					"cosh",
					"decimal",
					"default",
					"degree",
					"denominator",
					"diam",
					"dieonerrormode",
					"diff",
					"dirtyfindzeros",
					"dirtyinfnorm",
					"dirtyintegral",
					"dirtysimplify",
					"display",
					"div",
					"do",
					"double",
					"doubledouble",
					"doubleextended",
					"dyadic",
					"else",
					"end",
					"erf",
					"erfc",
					"error",
					"evaluate",
					"execute",
					"exp",
					"expand",
					"expm1",
					"exponent",
					"externalplot",
					"externalproc",
					"false",
					"file",
					"findzeros",
					"fixed",
					"floating",
					"floor",
					"for",
					"fpfindzeros",
					"fpminimax",
					"from",
					"fullparentheses",
					"function",
					"gcd",
					"getsuppressedmessages",
					"getbacktrace",
					"guessdegree",
					"head",
					"hexadecimal",
					"honorcoeffprec",
					"hopitalrecursions",
					"horner",
					"if",
					"implementpoly",
					"implementconstant",
					"in",
					"inf",
					"infnorm",
					"integer",
					"integral",
					"isbound",
					"isevaluable",
					"length",
					"library",
					"libraryconstant",
					"list",
					"log",
					"log10",
					"log1p",
					"log2",
					"mantissa",
					"max",
					"mid",
					"midpointmode",
					"min",
					"mod",
					"nearestint",
					"numberroots",
					"nop",
					"numerator",
					"object",
					"objectname",
					"of",
					"off",
					"on",
					"parse",
					"perturb",
					"pi",
					"plot",
					"points",
					"postscript",
					"postscriptfile",
					"powers",
					"prec",
					"precision",
					"print",
					"printbinary",
					"printdouble",
					"printexpansion",
					"printfloat",
					"printhexa",
					"printsingle",
					"printxml",
					"proc",
					"procedure",
					"quit",
					"range",
					"rationalapprox",
					"rationalmode",
					"readfile",
					"readxml",
					"relative",
					"remez",
					"rename",
					"restart",
					"return",
					"revert",
					"round",
					"roundcoefficients",
					"roundcorrectly",
					"roundingwarnings",
					"searchgal",
					"showmessagenumbers",
					"simplify",
					"simplifysafe",
					"sin",
					"single",
					"sinh",
					"sort",
					"sqrt",
					"string",
					"subpoly",
					"substitute",
					"sup",
					"supnorm",
					"suppressmessage",
					"tail",
					"tan",
					"tanh",
					"taylor",
					"taylorform",
					"taylorrecursions",
					"then",
					"time",
					"timing",
					"to",
					"tripledouble",
					"true",
					"unsuppressmessage",
					"var",
					"verbosity",
					"version",
					"void",
					"while",
					"worstcase",
					"write",
					"zerodenominators",
					"_x_",
					"quad",
					"halfprecision",
					"QD",
					"HP",
					"bashevaluate",
					"help",
					"match",
					"with",
					NULL }; /* The last entry must be NULL */

chain *openedFunctionLibraries = NULL;
chain *openedConstantLibraries = NULL;
chain *openedProcLibraries = NULL;

chain *globalLibraryFunctions = NULL;
chain *globalLibraryConstants = NULL;
chain *globalLibraryProcedures = NULL;

void copyIdentifierSymbols(char *ptr, const char *src) {
  char *currPtr;
  const char *currSrc;

  currPtr = ptr;
  currSrc = src;
  while ((*currPtr = *currSrc) != '\0') {
    if (currPtr == ptr) {
      if (((*currSrc >= 'a') && (*currSrc <= 'z')) ||
	  ((*currSrc >= 'A') && (*currSrc <= 'Z')) ||
	  (*currSrc == '_') ||
	  (*currSrc == '$')) {
	currPtr++;
      }
    } else {
      if (((*currSrc >= 'a') && (*currSrc <= 'z')) ||
	  ((*currSrc >= 'A') && (*currSrc <= 'Z')) ||
	  ((*currSrc >= '0') && (*currSrc <= '9')) ||
	  (*currSrc == '_') ||
	  (*currSrc == '$')) {
	currPtr++;
      }
    }
    currSrc++;
  }
}

char *getBaseFunctionName(void *func, const char *base, int hasData, void *data) {
  char *myFuncName, *funcName;
#if (defined(HAVE_DLADDR) && (HAVE_DLADDR) && defined(HAVE_DL_INFO_TYPE) && (HAVE_DL_INFO_TYPE))
  int errorOccurred;
  Dl_info myInfo;
  char *temp;
  errorOccurred = 0;
  if (dladdr(func, &myInfo) != 0) {
    if ((myInfo.dli_sname != NULL) && (myInfo.dli_saddr != NULL)) {
      if (myInfo.dli_saddr == func) {
	if (strlen(myInfo.dli_sname) > 0) {
	  temp = (char *) safeCalloc(strlen(myInfo.dli_sname) + 8 * sizeof(void *) + 1 + 1, sizeof(char));
	  if (hasData) {
	    sprintf(temp, "%s_%p", myInfo.dli_sname, data);
	  } else {
	    sprintf(temp, "%s", myInfo.dli_sname);
	  }
	  myFuncName = (char *) safeCalloc(strlen(myInfo.dli_sname) + 8 * sizeof(void *) + 1 + 1, sizeof(char));
	  copyIdentifierSymbols(myFuncName, temp);
	  safeFree(temp);
	} else {
	  errorOccurred = 1;
	}
      } else {
	errorOccurred = 1;
      }
    } else {
      errorOccurred = 1;
    }
  } else {
    errorOccurred = 1;
  }
  if (errorOccurred) {
    myFuncName = (char *) safeCalloc(strlen(base) + 3 + 8 * sizeof(void *) + 1 + 8 * sizeof(void *) + 1 + 8 * sizeof(void *) + 1,sizeof(char));
    if (hasData) {
      sprintf(myFuncName, "%s_%p_%p", base, func, data);
    } else {
      sprintf(myFuncName, "%s_%p", base, func);
    }
  }
#else
  myFuncName = (char *) safeCalloc(strlen(base) + 3 + 8 * sizeof(void *) + 1 + 8 * sizeof(void *) + 1 + 8 * sizeof(void *) + 1,sizeof(char));
  if (hasData) {
    sprintf(myFuncName, "%s_%p_%p", base, func, data);
  } else {
    sprintf(myFuncName, "%s_%p", base, func);
  }
#endif
  funcName = (char *) safeCalloc(strlen(myFuncName)+1,sizeof(char));
  strcpy(funcName, myFuncName);
  safeFree(myFuncName);
  return funcName;
}

char *filterSymbolName(char *basename) {
  char *res;

  res = (char *) safeCalloc(strlen(basename)+1,sizeof(char));

  copyIdentifierSymbols(res, basename);

  return res;
}

int symbolNameIsKeyword(char *name) {
  int i;

  for (i=0; sollyaKeywords[i] != NULL; i++) {
    if (!strcmp(name, sollyaKeywords[i])) return 1;
  }

  return 0;
}

int isValidIdentifier(char *name) {
  char *tmp;
  int equal;
  
  if (symbolNameIsKeyword(name)) return 0;
  tmp = (char *) safeCalloc(strlen(name) + 1, sizeof(char));
  copyIdentifierSymbols(tmp, name);
  equal = (strcmp(tmp, name) == 0);
  safeFree(tmp);
  return equal;
}

char *unifySymbolName(char *basename) {
  char *res, *str;
  uint64_t number;
  int found;
  size_t lenbasename;

  if (!(symbolNameAlreadyUsed(basename) || symbolNameIsKeyword(basename))) {
    res = (char *) safeCalloc(strlen(basename)+1,sizeof(char));
    strcpy(res, basename);
    return res;
  }

  lenbasename = strlen(basename);
  number = 0;
  found = 0;
  str = NULL;
  do {
    if (str != NULL) {
      safeFree(str);
      str = NULL;
    }
    str = (char *) safeCalloc(lenbasename + 1 + 8 * sizeof(uint64_t) + 1, sizeof(char));
    sprintf(str,"%s_%" PRIu64, basename, number);
    if (!(symbolNameAlreadyUsed(str) || symbolNameIsKeyword(str))) {
      found = 1;
    } else {
      number++;
    }
  } while (!found);

  res = safeCalloc(strlen(str)+1,sizeof(char));
  strcpy(res, str);
  safeFree(str);

  return res;
}

libraryHandle *getLibraryHandle(char *libraryName, int type) {
  chain *curr;
  chain *openedLibraries;
  libraryHandle *currHandle;
  void *dlfcnHandle, *otherDlfcnHandle;
  int shallowCopy, foundElsewhere, i;
  chain *otherLibraries[2];
  
  openedLibraries = NULL;

  switch(type) {
  case LIBRARYFUNCTION_CASE: openedLibraries = openedFunctionLibraries; break;
  case LIBRARYCONSTANT_CASE: openedLibraries = openedConstantLibraries; break;
  case LIBRARYPROC_CASE: openedLibraries = openedProcLibraries; break;
  }

  curr = openedLibraries;
  while (curr != NULL) {
    currHandle = (libraryHandle *) curr->value;
    if (strcmp(currHandle->libraryName,libraryName) == 0)
      return currHandle;
    curr = curr->next;
  }

  switch (type) {
  case LIBRARYFUNCTION_CASE:
    otherLibraries[0] = openedConstantLibraries;
    otherLibraries[1] = openedProcLibraries;
    break;
  case LIBRARYCONSTANT_CASE:
    otherLibraries[0] = openedFunctionLibraries;
    otherLibraries[1] = openedProcLibraries;
    break;
  case LIBRARYPROC_CASE:
    otherLibraries[0] = openedConstantLibraries;
    otherLibraries[1] = openedFunctionLibraries;
    break;
  default:
    sollyaFprintf(stderr,"Error: getLibraryHandle: unknown type (%d)\n", type);
    exit(1);
  }
  otherDlfcnHandle = NULL;
  for (i=0,foundElsewhere=0;(i<2)&&(!foundElsewhere);i++) {
    curr = otherLibraries[i];
    while (curr != NULL) {
      currHandle = (libraryHandle *) curr->value;
      if (strcmp(currHandle->libraryName,libraryName) == 0) {
	otherDlfcnHandle = currHandle->libraryDescriptor;
	foundElsewhere = 1;
	break;
      }
      curr = curr->next;
    }
  }

  if (!foundElsewhere) {
    dlerror();
    dlfcnHandle = dlopen(libraryName, RTLD_LOCAL | RTLD_NOW);
    if (dlfcnHandle == NULL)
      return NULL;
    shallowCopy = 0;
  } else {
    dlfcnHandle = otherDlfcnHandle;
    shallowCopy = 1;
  }

  currHandle = (libraryHandle *) safeMalloc(sizeof(libraryHandle));
  currHandle->libraryName = (char *) safeCalloc(strlen(libraryName)+1,sizeof(char));
  strcpy(currHandle->libraryName,libraryName);
  currHandle->libraryDescriptor = dlfcnHandle;
  currHandle->functionList = NULL;
  currHandle->shallowCopy = shallowCopy;

  openedLibraries = addElement(openedLibraries,currHandle);

  switch(type) {
  case LIBRARYFUNCTION_CASE: openedFunctionLibraries = openedLibraries; break;
  case LIBRARYCONSTANT_CASE: openedConstantLibraries = openedLibraries; break;
  case LIBRARYPROC_CASE: openedProcLibraries = openedLibraries; break;
  }

  return currHandle;
}

libraryHandle *getLibraryFunctionHandle(char *libraryName) {
  return getLibraryHandle(libraryName, LIBRARYFUNCTION_CASE);
}

libraryHandle *getConstantLibraryHandle(char *libraryName) {
  return getLibraryHandle(libraryName, LIBRARYCONSTANT_CASE);
}

libraryHandle *getProcLibraryHandle(char *libraryName) {
  return getLibraryHandle(libraryName, LIBRARYPROC_CASE);
}


/* Functions related to library functions */
libraryFunction *bindFunction(char* libraryName, char *functionName) {
  libraryHandle *libHandle;
  libraryFunction *currFunct;
  int (*myFunction)(sollya_mpfi_t, sollya_mpfi_t, int);
  char *error;
  sollya_mpfi_t op, rop;

  currFunct = getFunction(functionName);
  if (currFunct != NULL) {
    printMessage(1,SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_ALREADY_BOUND,"Warning: a function named \"%s\" has already been bound.\n",functionName);
    return currFunct;
  }

  libHandle = getLibraryFunctionHandle(libraryName);
  if (libHandle == NULL) {
    printMessage(1,SOLLYA_MSG_COULD_NOT_OPEN_LIBRARY_WITH_EXTERN_FUNC_OR_PROC,"Error: could not open library \"%s\" for binding \"%s\": %s\n",libraryName,functionName,dlerror());
    return NULL;
  }

  dlerror();
  myFunction = (int (*)(sollya_mpfi_t, sollya_mpfi_t, int)) dlsym(libHandle->libraryDescriptor, functionName);
  if ((error = dlerror()) != NULL) {
    printMessage(1, SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_NOT_FOUND_IN_LIBRARY, "Error: could not find function \"%s\" in library \"%s\" for binding: %s\n",functionName,libraryName,error);
    return NULL;
  }

  sollya_mpfi_init2(rop,12);
  sollya_mpfi_init2(op,12);
  sollya_mpfi_set_d(op,1.0);

  enterExternalCode();
  myFunction(rop,op,0);
  leaveExternalCode();

  sollya_mpfi_clear(rop);
  sollya_mpfi_clear(op);

  currFunct = (libraryFunction *) safeMalloc(sizeof(libraryFunction));
  currFunct->functionName = (char *) safeCalloc(strlen(functionName)+1,sizeof(char));
  strcpy(currFunct->functionName,functionName);
  currFunct->code = (void *) myFunction;
  currFunct->hasData = 0;
  currFunct->data = NULL;
  currFunct->dealloc = NULL;

  libHandle->functionList = addElement(libHandle->functionList,currFunct);

  return currFunct;
}

libraryFunction *getFunction(char *functionName) {
  chain *currLibList, *currFunList;
  libraryFunction *currFunct;
  libraryHandle *currLibHandle;

  for (currFunList=globalLibraryFunctions;currFunList!=NULL;currFunList=currFunList->next) {
    currFunct = (libraryFunction *) currFunList->value;
    if (strcmp(currFunct->functionName,functionName) == 0)
      return currFunct;
  }

  currLibList = openedFunctionLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    currFunList = currLibHandle->functionList;
    while (currFunList != NULL) {
      currFunct = (libraryFunction *) currFunList->value;
      if (strcmp(currFunct->functionName,functionName) == 0)
        return currFunct;
      currFunList = currFunList->next;
    }
    currLibList = currLibList->next;
  }

  return NULL;
}

libraryFunction *getFunctionByPtr(void *func, int hasData, void *data) {
  chain *currLibList, *currFunList;
  libraryFunction *currFunct;
  libraryHandle *currLibHandle;

  for (currFunList=globalLibraryFunctions;currFunList!=NULL;currFunList=currFunList->next) {
    currFunct = (libraryFunction *) currFunList->value;
    if ((currFunct->code == func) && ((!(currFunct->hasData && hasData)) || (currFunct->data == data)))
      return currFunct;
  }

  currLibList = openedFunctionLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    currFunList = currLibHandle->functionList;
    while (currFunList != NULL) {
      currFunct = (libraryFunction *) currFunList->value;
      if ((currFunct->code == func) && ((!(currFunct->hasData && hasData)) || (currFunct->data == data)))
        return currFunct;
      currFunList = currFunList->next;
    }
    currLibList = currLibList->next;
  }

  return NULL;
}

static libraryFunction *__bindFunctionByPtrImpl(char *suggestedName, void *func, int hasData, void *data, void *dealloc) {
  libraryFunction *res;
  char *unifiedName, *basename, *filteredSuggestedName, *filteredBaseName;

  res = getFunctionByPtr(func, hasData, data);
  if (res != NULL) {
    if (hasData) {
      if (res->hasData) {
	if (res->dealloc != dealloc) {
	  res->dealloc = dealloc;
	}
      }
    }
    return res;
  }
    
  if (suggestedName != NULL) {
    filteredSuggestedName = filterSymbolName(suggestedName);
    if (filteredSuggestedName[0] == '\0') {
      basename = getBaseFunctionName(func, "func", hasData, data);
      filteredBaseName = filterSymbolName(basename);
      safeFree(basename);
      if (filteredBaseName[0] == '\0') {
	unifiedName = unifySymbolName("func");
      } else {
	unifiedName = unifySymbolName(filteredBaseName);
      }
      safeFree(filteredBaseName);
    } else {
      unifiedName = unifySymbolName(filteredSuggestedName);
    }
    safeFree(filteredSuggestedName);
  } else {
    basename = getBaseFunctionName(func, "func", hasData, data);
    filteredBaseName = filterSymbolName(basename);
    safeFree(basename);
    if (filteredBaseName[0] == '\0') {
      unifiedName = unifySymbolName("func");
    } else {
      unifiedName = unifySymbolName(filteredBaseName);
    }
    safeFree(filteredBaseName);
  }

  res = (libraryFunction *) safeMalloc(sizeof(libraryFunction));
  res->functionName = unifiedName;
  res->code = func;
  res->hasData = hasData;
  res->data = data;
  res->dealloc = dealloc;

  globalLibraryFunctions = addElement(globalLibraryFunctions, res);

  return res;
}

libraryFunction *bindFunctionByPtr(char *suggestedName, int (*func)(mpfi_t, mpfi_t, int)) {
  return __bindFunctionByPtrImpl(suggestedName, ((void *) func), 0, NULL, NULL);
}

libraryFunction *bindFunctionByPtrWithData(char *suggestedName, int (*func)(mpfi_t, mpfi_t, int, void *), void *data, void (*dealloc)(void *)) {
  return __bindFunctionByPtrImpl(suggestedName, ((void *) func), 1, data, ((void *) dealloc));
}

void freeFunctionLibraries() {
  chain *currLibList, *currFunList, *prevFunList, *prevLibList;
  libraryFunction *currFunct;
  libraryHandle *currLibHandle;
  int (*myFunction)(void);
  int res;

  currFunList = globalLibraryFunctions;
  while (currFunList != NULL) {
    currFunct = (libraryFunction *) currFunList->value;
    if (currFunct->dealloc != NULL) {
      if (currFunct->hasData) {
	enterExternalCode();
	((void (*)(void *)) currFunct->dealloc)(currFunct->data);
	leaveExternalCode();
      } else {
	enterExternalCode();
	((void (*)()) currFunct->dealloc)();
	leaveExternalCode();
      }
    }
    safeFree(currFunct->functionName);
    safeFree(currFunList->value);
    prevFunList = currFunList;
    currFunList = currFunList->next;
    safeFree(prevFunList);
  }

  currLibList = openedFunctionLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    if (!(currLibHandle->shallowCopy)) {
      dlerror();
      myFunction = (int (*)(void)) dlsym(currLibHandle->libraryDescriptor, "sollya_external_lib_close");
      if (dlerror() == NULL) {
	enterExternalCode();
	res = myFunction();
	leaveExternalCode();
	if (res) {
	  printMessage(1,SOLLYA_MSG_LIBRARY_CLOSER_ERROR,"Warning: while closing library \"%s\", the function \"sollya_external_lib_close\" was found and called but it signaled the error %d\n",currLibHandle->libraryName,res);
	}
      }
    }
    currFunList = currLibHandle->functionList;
    while (currFunList != NULL) {
      currFunct = (libraryFunction *) currFunList->value;
      if (currFunct->dealloc != NULL) {
	if (currFunct->hasData) {
	  enterExternalCode();
	  ((void (*)(void *)) currFunct->dealloc)(currFunct->data);
	  leaveExternalCode();
	} else {
	  enterExternalCode();
	  ((void (*)()) currFunct->dealloc)();
	  leaveExternalCode();
	}
      }
      safeFree(currFunct->functionName);
      safeFree(currFunList->value);
      prevFunList = currFunList;
      currFunList = currFunList->next;
      safeFree(prevFunList);
    }
    if (!(currLibHandle->shallowCopy)) {
      dlerror();
      if (dlclose(currLibHandle->libraryDescriptor) != 0)
	printMessage(1,SOLLYA_MSG_COULD_NOT_CLOSE_LIBRARY,"Warning: could not close library \"%s\": %s\n",currLibHandle->libraryName,dlerror());
    }
    safeFree(currLibHandle->libraryName);
    safeFree(currLibHandle);
    prevLibList = currLibList;
    currLibList = currLibList->next;
    safeFree(prevLibList);
  }
  openedFunctionLibraries = NULL;
}



/* Functions related to library constants */
libraryFunction *bindConstantFunction(char* libraryName, char *functionName) {
  libraryHandle *libHandle;
  libraryFunction *currFunct;
  void (*myFunction)(mpfr_t, mp_prec_t);
  char *error;
  mpfr_t op;

  currFunct = getConstantFunction(functionName);
  if (currFunct != NULL) {
    printMessage(1,SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_ALREADY_BOUND,"Warning: a function named \"%s\" has already been bound.\n",functionName);
    return currFunct;
  }

  libHandle = getConstantLibraryHandle(libraryName);
  if (libHandle == NULL) {
    printMessage(1,SOLLYA_MSG_COULD_NOT_OPEN_LIBRARY_WITH_EXTERN_FUNC_OR_PROC,"Error: could not open library \"%s\" for binding \"%s\": %s\n",libraryName,functionName,dlerror());
    return NULL;
  }

  dlerror();
  myFunction = (void (*)(mpfr_t, mp_prec_t)) dlsym(libHandle->libraryDescriptor, functionName);
  if ((error = dlerror()) != NULL) {
    printMessage(1, SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_NOT_FOUND_IN_LIBRARY, "Error: could not find function \"%s\" in library \"%s\" for binding: %s\n",functionName,libraryName,error);
    return NULL;
  }

  mpfr_init2(op,20);

  enterExternalCode();
  myFunction(op,5);
  leaveExternalCode();

  mpfr_clear(op);

  currFunct = (libraryFunction *) safeMalloc(sizeof(libraryFunction));
  currFunct->functionName = (char *) safeCalloc(strlen(functionName)+1,sizeof(char));
  strcpy(currFunct->functionName,functionName);
  currFunct->code = (void *) myFunction;
  currFunct->hasData = 0;
  currFunct->data = NULL;
  currFunct->dealloc = NULL;

  libHandle->functionList = addElement(libHandle->functionList,currFunct);

  return currFunct;
}

libraryFunction *getConstantFunction(char *functionName) {
  chain *currLibList, *currFunList;
  libraryFunction *currFunct;
  libraryHandle *currLibHandle;

  for (currFunList=globalLibraryConstants;currFunList!=NULL;currFunList=currFunList->next) {
    currFunct = (libraryFunction *) currFunList->value;
    if (strcmp(currFunct->functionName,functionName) == 0)
      return currFunct;
  }

  currLibList = openedConstantLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    currFunList = currLibHandle->functionList;
    while (currFunList != NULL) {
      currFunct = (libraryFunction *) currFunList->value;
      if (strcmp(currFunct->functionName,functionName) == 0)
        return currFunct;
      currFunList = currFunList->next;
    }
    currLibList = currLibList->next;
  }

  return NULL;
}

libraryFunction *getConstantFunctionByPtr(void *func, int hasData, void *data) {
  chain *currLibList, *currFunList;
  libraryFunction *currFunct;
  libraryHandle *currLibHandle;

  for (currFunList=globalLibraryConstants;currFunList!=NULL;currFunList=currFunList->next) {
    currFunct = (libraryFunction *) currFunList->value;
    if ((currFunct->code == func) && ((!(currFunct->hasData && hasData)) || (currFunct->data == data)))
      return currFunct;
  }

  currLibList = openedConstantLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    currFunList = currLibHandle->functionList;
    while (currFunList != NULL) {
      currFunct = (libraryFunction *) currFunList->value;
      if ((currFunct->code == func) && ((!(currFunct->hasData && hasData)) || (currFunct->data == data)))
        return currFunct;
      currFunList = currFunList->next;
    }
    currLibList = currLibList->next;
  }

  return NULL;
}

static libraryFunction *__bindConstantFunctionByPtrImpl(char *suggestedName, void *func, int hasData, void *data, void *dealloc) {
  libraryFunction *res;
  char *unifiedName, *basename, *filteredBaseName, *filteredSuggestedName;

  res = getConstantFunctionByPtr(func, hasData, data);
  if (res != NULL) {
    if (hasData) {
      if (res->hasData) {
	if (res->dealloc != dealloc) {
	  res->dealloc = dealloc;
	}
      }
    }
    return res;
  }

  if (suggestedName != NULL) {
    filteredSuggestedName = filterSymbolName(suggestedName);
    if (filteredSuggestedName[0] == '\0') {
      basename = getBaseFunctionName(func, "const", hasData, data);
      filteredBaseName = filterSymbolName(basename);
      safeFree(basename);
      if (filteredBaseName[0] == '\0') {
	unifiedName = unifySymbolName("const");
      } else {
	unifiedName = unifySymbolName(filteredBaseName);
      }
      safeFree(filteredBaseName);
    } else {
      unifiedName = unifySymbolName(filteredSuggestedName);
    }
    safeFree(filteredSuggestedName);
  } else {
    basename = getBaseFunctionName(func, "const", hasData, data);
    filteredBaseName = filterSymbolName(basename);
    safeFree(basename);
    if (filteredBaseName[0] == '\0') {
      unifiedName = unifySymbolName("const");
    } else {
      unifiedName = unifySymbolName(filteredBaseName);
    }
    safeFree(filteredBaseName);
  }

  res = (libraryFunction *) safeMalloc(sizeof(libraryFunction));
  res->functionName = unifiedName;
  res->code = func;
  res->hasData = hasData;
  res->data = data;
  res->dealloc = dealloc;
  

  globalLibraryConstants = addElement(globalLibraryConstants, res);

  return res;
}

libraryFunction *bindConstantFunctionByPtr(char *suggestedName, void (*func)(mpfr_t, mp_prec_t)) {
  return __bindConstantFunctionByPtrImpl(suggestedName, ((void *) func), 0, NULL, NULL);
}

libraryFunction *bindConstantFunctionByPtrWithData(char *suggestedName, void (*func)(mpfr_t, mp_prec_t, void *), void *data, void (*dealloc)(void *)) {
  return __bindConstantFunctionByPtrImpl(suggestedName, ((void *) func), 1, data, ((void *) dealloc));
}

void freeConstantLibraries() {
  chain *currLibList, *currFunList, *prevFunList, *prevLibList;
  libraryFunction *currFunct;
  libraryHandle *currLibHandle;
  int res;
  int (*myFunction)(void);

  currFunList = globalLibraryConstants;
  while (currFunList != NULL) {
    currFunct = (libraryFunction *) currFunList->value;
    if (currFunct->dealloc != NULL) {
      if (currFunct->hasData) {
	enterExternalCode();
	((void (*)(void *)) currFunct->dealloc)(currFunct->data);
	leaveExternalCode();
      } else {
	enterExternalCode();
	((void (*)()) currFunct->dealloc)();
	leaveExternalCode();
      }
    }
    safeFree(currFunct->functionName);
    safeFree(currFunList->value);
    prevFunList = currFunList;
    currFunList = currFunList->next;
    safeFree(prevFunList);
  }

  currLibList = openedConstantLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    if (!(currLibHandle->shallowCopy)) {
      dlerror();
      myFunction = (int (*)(void)) dlsym(currLibHandle->libraryDescriptor, "sollya_external_lib_close");
      if (dlerror() == NULL) {
	enterExternalCode();
	res = myFunction();
	leaveExternalCode();
	if (res) {
	  printMessage(1,SOLLYA_MSG_LIBRARY_CLOSER_ERROR,"Warning: while closing library \"%s\", the function \"sollya_external_lib_close\" was found and called but it signaled the error %d\n",currLibHandle->libraryName,res);
	}
      }
    }
    currFunList = currLibHandle->functionList;
    while (currFunList != NULL) {
      currFunct = (libraryFunction *) currFunList->value;
      if (currFunct->dealloc != NULL) {
	if (currFunct->hasData) {
	  enterExternalCode();
	  ((void (*)(void *)) currFunct->dealloc)(currFunct->data);
	  leaveExternalCode();
	} else {
	  enterExternalCode();
	  ((void (*)()) currFunct->dealloc)();
	  leaveExternalCode();
	}
      }
      safeFree(currFunct->functionName);
      safeFree(currFunList->value);
      prevFunList = currFunList;
      currFunList = currFunList->next;
      safeFree(prevFunList);
    }
    if (!(currLibHandle->shallowCopy)) {
      dlerror();
      if (dlclose(currLibHandle->libraryDescriptor) != 0)
	printMessage(1,SOLLYA_MSG_COULD_NOT_CLOSE_LIBRARY,"Warning: could not close library \"%s\": %s\n",currLibHandle->libraryName,dlerror());
    }
    safeFree(currLibHandle->libraryName);
    safeFree(currLibHandle);
    prevLibList = currLibList;
    currLibList = currLibList->next;
    safeFree(prevLibList);
  }
  openedConstantLibraries = NULL;
}

/* Functions related to external procedures */
libraryProcedure *bindProcedure(char* libraryName, char *procedureName, chain *signature) {
  libraryHandle *libHandle;
  libraryProcedure *currProc;
  char *error;
  void *myFunction;

  currProc = getProcedure(procedureName);
  if (currProc != NULL) {
    printMessage(1,SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_ALREADY_BOUND,"Warning: a function named \"%s\" has already been bound.\n",procedureName);
    return currProc;
  }

  libHandle = getProcLibraryHandle(libraryName);
  if (libHandle == NULL) {
    printMessage(1,SOLLYA_MSG_COULD_NOT_OPEN_LIBRARY_WITH_EXTERN_FUNC_OR_PROC,"Error: could not open library \"%s\" for binding \"%s\": %s\n",libraryName,procedureName,dlerror());
    return NULL;
  }

  dlerror();
  myFunction = dlsym(libHandle->libraryDescriptor, procedureName);
  if ((error = dlerror()) != NULL) {
    printMessage(1, SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_NOT_FOUND_IN_LIBRARY, "Error: could not find function \"%s\" in library \"%s\" for binding: %s\n",procedureName,libraryName,error);
    return NULL;
  }

  currProc = (libraryProcedure *) safeMalloc(sizeof(libraryProcedure));
  currProc->procedureName = (char *) safeCalloc(strlen(procedureName)+1,sizeof(char));
  strcpy(currProc->procedureName,procedureName);
  currProc->code = myFunction;
  currProc->signature = copyChainWithoutReversal(signature, copyIntPtrOnVoid);
  currProc->hasData = 0;
  currProc->data = NULL;
  currProc->dealloc = NULL;


  libHandle->functionList = addElement(libHandle->functionList,currProc);

  return currProc;
}

libraryProcedure *getProcedure(char *procedureName) {
  chain *currLibList, *currProcList;
  libraryProcedure *currProc;
  libraryHandle *currLibHandle;

  for (currProcList=globalLibraryProcedures;currProcList!=NULL;currProcList=currProcList->next) {
    currProc = (libraryProcedure *) currProcList->value;
    if (strcmp(currProc->procedureName,procedureName) == 0)
      return currProc;
  }

  currLibList = openedProcLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    currProcList = currLibHandle->functionList;
    while (currProcList != NULL) {
      currProc = (libraryProcedure *) currProcList->value;
      if (strcmp(currProc->procedureName,procedureName) == 0)
	return currProc;
      currProcList = currProcList->next;
    }
    currLibList = currLibList->next;
  }

  return NULL;
}

libraryProcedure *getProcedureByPtr(void *ptr, int hasData, void *data) {
  chain *currLibList, *currProcList;
  libraryProcedure *currProc;
  libraryHandle *currLibHandle;

  for (currProcList=globalLibraryProcedures;currProcList!=NULL;currProcList=currProcList->next) {
    currProc = (libraryProcedure *) currProcList->value;
    if ((currProc->code == ptr) && ((!(currProc->hasData && hasData)) || (currProc->data == data)))
      return currProc;
  }

  currLibList = openedProcLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    currProcList = currLibHandle->functionList;
    while (currProcList != NULL) {
      currProc = (libraryProcedure *) currProcList->value;
      if ((currProc->code == ptr) && ((!(currProc->hasData && hasData)) || (currProc->data == data)))
	return currProc;
      currProcList = currProcList->next;
    }
    currLibList = currLibList->next;
  }

  return NULL;
}

static libraryProcedure *__bindProcedureByPtrImpl(int resType, int *argTypes, int arity, char *suggestedName, void *func, int hasData, void *data, void *dealloc) {
  libraryProcedure *res;
  char *unifiedName, *basename, *filteredBaseName, *filteredSuggestedName;
  chain *signature, *temp;
  int i;
  int *t;

  if (arity < 0) return NULL;
  if (arity > 1) {
    for (i=1;i<arity;i++) {
      if (argTypes[i] == VOID_TYPE)
	return NULL;
    }
  }
  
  res = getProcedureByPtr(func, hasData, data);
  if (res != NULL) {
    if (hasData) {
      if (res->hasData) {
	if (res->dealloc != dealloc) {
	  res->dealloc = dealloc;
	}
      }
    }
    return res;
  }

  signature = NULL;
  t = (int *) safeMalloc(sizeof(int));
  *t = resType;
  signature = addElement(signature, t);
  if (arity == 0) {
    t = (int *) safeMalloc(sizeof(int));
    *t = VOID_TYPE;
    signature = addElement(signature, t);
  } else {
    for (i=0;i<arity;i++) {
      t = (int *) safeMalloc(sizeof(int));
      *t = argTypes[i];
      signature = addElement(signature, t);
    }
  }
  temp = copyChain(signature, copyIntPtrOnVoid);
  freeChain(signature,freeIntPtr);
  signature = temp;
  
  if (suggestedName != NULL) {
    filteredSuggestedName = filterSymbolName(suggestedName);
    if (filteredSuggestedName[0] == '\0') {
      basename = getBaseFunctionName(func, "proc", hasData, data);
      filteredBaseName = filterSymbolName(basename);
      safeFree(basename);
      if (filteredBaseName[0] == '\0') {
	unifiedName = unifySymbolName("proc");
      } else {
	unifiedName = unifySymbolName(filteredBaseName);
      }
      safeFree(filteredBaseName);
    } else {
      unifiedName = unifySymbolName(filteredSuggestedName);
    }
    safeFree(filteredSuggestedName);
  } else {
    basename = getBaseFunctionName(func, "proc", hasData, data);
    filteredBaseName = filterSymbolName(basename);
    safeFree(basename);
    if (filteredBaseName[0] == '\0') {
      unifiedName = unifySymbolName("proc");
    } else {
      unifiedName = unifySymbolName(filteredBaseName);
    }
    safeFree(filteredBaseName);
  }

  res = (libraryProcedure *) safeMalloc(sizeof(libraryProcedure));
  res->procedureName = unifiedName;
  res->code = func;
  res->signature = signature;
  res->hasData = hasData;
  res->data = data;
  res->dealloc = dealloc;

  globalLibraryProcedures = addElement(globalLibraryProcedures, res);

  return res;
}

libraryProcedure *bindProcedureByPtr(int resType, int *argTypes, int arity, char *suggestedName, void *func) {
  return __bindProcedureByPtrImpl(resType, argTypes, arity, suggestedName, func, 0, NULL, NULL);
}

libraryProcedure *bindProcedureByPtrWithData(int resType, int *argTypes, int arity, char *suggestedName, void *func, void *data, void (*dealloc)(void *)) {
  return __bindProcedureByPtrImpl(resType, argTypes, arity, suggestedName, func, 1, data, ((void *) dealloc));
}

void freeProcLibraries() {
  chain *currLibList, *currProcList, *prevProcList, *prevLibList;
  libraryProcedure *currProc;
  libraryHandle *currLibHandle;
  int res;
  int (*myFunction)(void);

  currProcList = globalLibraryProcedures;
  while (currProcList != NULL) {
    currProc = (libraryProcedure *) currProcList->value;
    if (currProc->dealloc != NULL) {
      if (currProc->hasData) {
	enterExternalCode();
	((void (*)(void *)) currProc->dealloc)(currProc->data);
	leaveExternalCode();
      } else {
	enterExternalCode();
	((void (*)()) currProc->dealloc)();
	leaveExternalCode();
      }
    }
    safeFree(currProc->procedureName);
    freeChain(currProc->signature,freeIntPtr);
    safeFree(currProcList->value);
    prevProcList = currProcList;
    currProcList = currProcList->next;
    safeFree(prevProcList);
  }
  
  currLibList = openedProcLibraries;
  while (currLibList != NULL) {
    currLibHandle = (libraryHandle *) currLibList->value;
    if (!(currLibHandle->shallowCopy)) {
      dlerror();
      myFunction = (int (*)(void)) dlsym(currLibHandle->libraryDescriptor, "sollya_external_lib_close");
      if (dlerror() == NULL) {
	enterExternalCode();
	res = myFunction();
	leaveExternalCode();
	if (res) {
	  printMessage(1,SOLLYA_MSG_LIBRARY_CLOSER_ERROR,"Warning: while closing library \"%s\", the function \"sollya_external_lib_close\" was found and called but it signaled the error %d\n",currLibHandle->libraryName,res);
	}
      }
    }
    currProcList = currLibHandle->functionList;
    while (currProcList != NULL) {
      currProc = (libraryProcedure *) currProcList->value;
      if (currProc->dealloc != NULL) {
	if (currProc->hasData) {
	  ((void (*)(void *)) currProc->dealloc)(currProc->data);
	} else {
	  ((void (*)()) currProc->dealloc)();
	}
      }
      safeFree(currProc->procedureName);
      freeChain(currProc->signature,freeIntPtr);
      safeFree(currProcList->value);
      prevProcList = currProcList;
      currProcList = currProcList->next;
      safeFree(prevProcList);
    }
    if (!(currLibHandle->shallowCopy)) {
      dlerror();
      if (dlclose(currLibHandle->libraryDescriptor) != 0)
	printMessage(1,SOLLYA_MSG_COULD_NOT_CLOSE_LIBRARY,"Warning: could not close library \"%s\": %s\n",currLibHandle->libraryName,dlerror());
    }
    safeFree(currLibHandle->libraryName);
    safeFree(currLibHandle);
    prevLibList = currLibList;
    currLibList = currLibList->next;
    safeFree(prevLibList);
  }
  openedProcLibraries = NULL;
}

