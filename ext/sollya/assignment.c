/*

  Copyright 2006-2015 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

  and by

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

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
#include "chain.h"
#include "assignment.h"
#include "expression.h"
#include "execute.h"
#include "general.h"


chain *addEntry(chain *symTbl, char *name, void *value, void * (*copyValue) (void *)) {
  entry *newEntry;

  if (containsEntry(symTbl,name)) return symTbl;

  newEntry = (entry *) safeMalloc(sizeof(entry));
  newEntry->name = (char *) safeCalloc(strlen(name)+1,sizeof(char));
  strcpy(newEntry->name,name);
  newEntry->value = copyValue(value);
  symTbl = addElement(symTbl,newEntry);
  return symTbl;
}

int containsEntry(chain *symTbl, char *name) {
  chain *curr;

  curr = symTbl;
  while (curr != NULL) {
    if (strcmp(((entry *) (curr->value))->name,name) == 0) return 1;
    curr = curr->next;
  }

  return 0;
}

void *getEntry(chain *symTbl, char *name, void * (*copyValue) (void *)) {
  chain *curr;
  void *result;

  result = NULL;
  curr = symTbl;
  while (curr != NULL) {
    if (strcmp(((entry *) (curr->value))->name,name) == 0) {
      result = copyValue(((entry *) curr->value)->value);
      break;
    }
    curr = curr->next;
  }

  return result;
}

void freeEntry(void *e, void (*f) (void *)) {
  f(((entry *) e)->value);
  safeFree(((entry *) e)->name);
  safeFree(e);
}


chain *removeEntry(chain *symTbl, char *name, void (*f) (void *)) {
  chain *curr, *prev, *newSymTbl;

  curr = symTbl; prev = NULL;
  while (curr != NULL) {
    if (strcmp(((entry *) (curr->value))->name,name) == 0) {
      if ((prev == NULL) && (curr->next == NULL)) {
	newSymTbl = NULL;
      } else {
	if (prev == NULL) {
	  newSymTbl = curr->next;
	} else {
	  prev->next = curr->next;
	  newSymTbl = symTbl;
	}
      }
      freeEntry(((entry *) curr->value),f);
      safeFree(curr);
      return newSymTbl;
    }
    prev = curr;
    curr = curr->next;
  }

  return symTbl;
}


void freeSymbolTable(chain *symTbl, void (*f) (void *)) {
  if (symTbl != NULL) {
    if (symTbl->next != NULL) freeSymbolTable(symTbl->next,f);
    freeEntry(symTbl->value,f);
    safeFree(symTbl);
  }
}


void freeNothing(void *thing) {
  UNUSED_PARAM(thing); return;
}

void freeDeclaredSymbolTable(chain *declSymTbl, void (*f) (void *)) {
  chain *curr;

  curr = declSymTbl;
  while (curr != NULL) {
    freeSymbolTable((chain *) (curr->value), f);
    curr->value = NULL;
    curr = curr->next;
  }

  freeChain(declSymTbl, freeNothing);
}

char *getEntryName(chain *symTbl, chain *declSymTbl, void *value, int (*f)(void *, void *)) {
  chain *curr, *currFrame, *currDecl, *currFrameAbove, *currDeclAbove;
  int shadowed;
  char *res;

  /* Go over global symbol table first, for each hit check if the
     identifier isn't shadowed by the declared symbol table. 
  */
  for (curr=symTbl;curr!=NULL;curr=curr->next) {
    if (f(value, ((entry *) (curr->value))->value)) {
      /* We have a hit. Check if its name isn't shadowed */
      shadowed = 0;
      for (currFrame=declSymTbl;((currFrame!=NULL)&&(!shadowed));currFrame=currFrame->next) {
	for (currDecl=((chain *) (currFrame->value));((currDecl!=NULL)&&(!shadowed));currDecl=currDecl->next) {
	  if (strcmp(((entry *) (curr->value))->name, ((entry *) (currDecl->value))->name) == 0) {
	    shadowed = 1;
	  }
	}
      }
      if (!shadowed) {
	res = (char *) safeCalloc(strlen(((entry *) (curr->value))->name) + 1, sizeof(char));
	strcpy(res, ((entry *) (curr->value))->name);
	return res;
      }
    }
  }

  /* Now go over the declared symbol table, frame by frame. For each
     hit check if the identifier isn't shadowed by an identifier in an
     above frame of the declared symbol table.
  */
  for (currFrame=declSymTbl;currFrame!=NULL;currFrame=currFrame->next) {
    for (currDecl=((chain *)(currFrame->value));currDecl!=NULL;currDecl=currDecl->next) {
      if (f(value, ((entry *) (currDecl->value))->value)) {
	/* We have a hit. Check if its name isn't shadowed in a frame
	   above. 
	*/
	shadowed = 0;
	for (currFrameAbove=declSymTbl;
	     ((currFrameAbove!=NULL) &&
	      (currFrameAbove!=currFrame) &&
	      (!shadowed));
	     currFrameAbove=currFrameAbove->next) {
	  for (currDeclAbove=((chain *) (currFrameAbove->value));((currDeclAbove!=NULL)&&(!shadowed));currDeclAbove=currDeclAbove->next) {
	    if (strcmp(((entry *) (currDecl->value))->name, ((entry *) (currDeclAbove->value))->name) == 0) {
	      shadowed = 1;
	    }
	  }
	}
	if (!shadowed) {
	  res = (char *) safeCalloc(strlen(((entry *) (currDecl->value))->name) + 1, sizeof(char));
	  strcpy(res, ((entry *) (currDecl->value))->name);
	  return res;
	}
      }
    }
  }

  /* We didn't find anything or everything we found was shadowed. */
  return NULL;
}

chain *pushFrame(chain *declSymTbl) {
  return addElement(declSymTbl, NULL);
}

chain *popFrame(chain *declSymTbl, void (*f) (void *)) {
  chain *newDeclSymTbl;

  if (declSymTbl == NULL) return NULL;

  newDeclSymTbl = declSymTbl->next;

  freeSymbolTable((chain *) (declSymTbl->value), f);

  safeFree(declSymTbl);

  return newDeclSymTbl;
}

chain *declareNewEntry(chain *declSymTbl, char *name, void *value, void * (*copyValue) (void *)) {
  chain *newValue;

  if (declSymTbl == NULL) return NULL;

  if (containsEntry((chain *) (declSymTbl->value), name)) return declSymTbl;

  newValue = addEntry((chain *) (declSymTbl->value), name, value, copyValue);

  declSymTbl->value = newValue;

  return declSymTbl;
}

chain *replaceDeclaredEntry(chain *declSymTbl, char *name, void *value, void * (*copyValue) (void *), void (*freeValue) (void *)) {
  chain *curr;
  chain *newValue;

  if (declSymTbl == NULL) return NULL;

  curr = declSymTbl;
  while (curr != NULL) {
    if (containsEntry((chain *) (curr->value), name)) {
      newValue = removeEntry((chain *) (curr->value), name, freeValue);
      curr->value = newValue;
      newValue = addEntry((chain *) (curr->value), name, value, copyValue);
      curr->value = newValue;
      break;
    }
    curr = curr->next;
  }

  return declSymTbl;
}

int containsDeclaredEntry(chain *declSymTbl, char *name) {
  chain *curr;

  curr = declSymTbl;
  while (curr != NULL) {
    if (containsEntry((chain *) (curr->value), name)) return 1;
    curr = curr->next;
  }

  return 0;
}

void *getDeclaredEntry(chain *declSymTbl, char *name, void * (*copyValue) (void *)) {
  chain *curr;

  curr = declSymTbl;
  while (curr != NULL) {
    if (containsEntry((chain *) (curr->value), name)) return getEntry((chain *) (curr->value), name, copyValue);
    curr = curr->next;
  }

  return NULL;
}



chain *assignDeclaredEntry(chain *declSymTbl, char *name, void *value, void * (*copyValue) (void *), void (*freeValue) (void *)) {
  chain *newDeclSymTbl;

  if (containsDeclaredEntry(declSymTbl, name))
    newDeclSymTbl = replaceDeclaredEntry(declSymTbl, name, value, copyValue, freeValue);
  else
    newDeclSymTbl = declareNewEntry(declSymTbl, name, value, copyValue);

  return newDeclSymTbl;
}

int performListPrependOnEntry(chain *symTbl, char *ident, node *tree) {
  chain *curr, *newArgs;
  node *oldNode, *newNode;
  int okay;
  int oldChecked, addedChecked, newChecked;
  size_t neededSize;

  oldChecked = 0;
  newChecked = 0;
  addedChecked = 0;
  if ((tree->nodeType == MEMREF) && (tree->cache->isCorrectlyTyped)) addedChecked = 1;

  okay = 0;
  curr = symTbl;
  while (curr != NULL) {
    if (strcmp(((entry *) (curr->value))->name,ident) == 0) {
      oldNode = (node *) (((entry *) curr->value)->value);
      if (oldNode->nodeType == MEMREF) {
	oldNode->cache->hashComputed = 0;
      }
      if ((oldNode->nodeType == MEMREF) && (oldNode->cache->isCorrectlyTyped)) oldChecked = 1;
      newChecked = oldChecked && addedChecked;
      while (1) {
	if (oldNode->nodeType != MEMREF) break; else { if (!newChecked) oldNode->cache->isCorrectlyTyped = 0; }
	if (oldNode->libFunDeriv > 1) break;
	oldNode = getMemRefChild(oldNode);
      }
      switch (oldNode->nodeType) {
      case MEMREF:
	if (!newChecked) oldNode->cache->isCorrectlyTyped = 0;
	if (oldNode->libFunDeriv > 1) {
	  if ((getMemRefChild(oldNode)->nodeType == LIST) ||
	      (getMemRefChild(oldNode)->nodeType == FINALELLIPTICLIST)) {
	    newArgs = addElement(copyChainWithoutReversal(getMemRefChild(oldNode)->arguments, copyThingOnVoid), tree);
	    newNode = allocateNode();
	    newNode->nodeType = getMemRefChild(oldNode)->nodeType;
	    newNode->argArray = NULL;
	    newNode->argArraySize = 0;
	    newNode->argArrayAllocSize = 0;
	    newNode->arguments = newArgs;
	    newNode = addMemRef(newNode);
	    if (newChecked && (newNode->nodeType == MEMREF)) {
	      newNode->cache->isCorrectlyTyped = 1;
	    }
	    ((entry *) curr->value)->value = newNode;
	    freeThing(oldNode);
	    okay = 1;
	  } else {
	    newNode = deepCopyThing(oldNode);
	    if ((newNode->nodeType == LIST) ||
		(newNode->nodeType == FINALELLIPTICLIST)) {
	      freeThing(oldNode);
	      newNode->arguments = addElement(newNode->arguments, tree);
	      newNode->argArray = NULL;
	      newNode->argArraySize = 0;
	      newNode->argArrayAllocSize = 0;
	      newNode = addMemRef(newNode);
	      /* We might want to re-use the re-usable memref'ed parts
		 of the original expression here 
	      */
	      ((entry *) curr->value)->value = newNode;
	      if (newChecked && (newNode->nodeType == MEMREF)) {
		newNode->cache->isCorrectlyTyped = 1;
	      }
	      okay = 1;
	    } else {
	      freeThing(newNode);
	    }
	  }
	}
	break;
      case LIST:
      case FINALELLIPTICLIST:
	oldNode->arguments = addElement(oldNode->arguments, tree);
	if (oldNode->argArray != NULL) {
	  neededSize = ((size_t) (oldNode->argArraySize + 1)) * sizeof(node *);
	  if (neededSize <= (oldNode->argArrayAllocSize)) {
	    oldNode->argArraySize++;
	    (oldNode->argArray)[(oldNode->argArraySize - 1) - 0] = tree;
	  } else {
	    if ((neededSize <= ((size_t) (2 * oldNode->argArrayAllocSize))) &&
		(((size_t) (2 * oldNode->argArrayAllocSize)) <= ((size_t) SOLLYA_MAX_ARG_ARRAY_ALLOC_SIZE)) &&
		(((size_t) (2 * oldNode->argArrayAllocSize)) > ((size_t) 0))) {
	      oldNode->argArrayAllocSize = (size_t) (2 * oldNode->argArrayAllocSize);
	      oldNode->argArray = safeRealloc(oldNode->argArray, oldNode->argArrayAllocSize);
	      oldNode->argArraySize++;
	      (oldNode->argArray)[(oldNode->argArraySize - 1) - 0] = tree;
	    } else {
	      safeFree(oldNode->argArray);
	      oldNode->argArray = NULL;
	      oldNode->argArraySize = 0;
	      oldNode->argArrayAllocSize = 0;
	    }
	  }
	}
	okay = 1;
	break;
      default:
	okay = 0;
	break;
      }
      break;
    }
    curr = curr->next;
  }

  return okay;
}

int performListPrependOnDeclaredEntry(chain *declSymTbl, char *name, node *tree) {
  chain *curr;

  curr = declSymTbl;
  while (curr != NULL) {
    if (containsEntry((chain *) (curr->value), name)) return performListPrependOnEntry((chain *) (curr->value), name, tree);
    curr = curr->next;
  }

  return 0;
}

int performListTailOnEntry(chain *symTbl, char *ident) {
  chain *curr, *newArgs;
  node *oldNode, *newNode;
  int okay;
  int oldChecked;

  oldChecked = 0;

  okay = 0;
  curr = symTbl;
  while (curr != NULL) {
    if (strcmp(((entry *) (curr->value))->name,ident) == 0) {
      oldNode = (node *) (((entry *) curr->value)->value);
      if (oldNode->nodeType == MEMREF) {
	oldNode->cache->hashComputed = 0;
      }
      if ((oldNode->nodeType == MEMREF) && (oldNode->cache->isCorrectlyTyped)) oldChecked = 1;
      while (1) {
	if (oldNode->nodeType != MEMREF) break;
	if (oldNode->libFunDeriv > 1) break;
	oldNode = getMemRefChild(oldNode);
      }
      switch (oldNode->nodeType) {
      case MEMREF:
	if (oldNode->libFunDeriv > 1) {
	  if (((getMemRefChild(oldNode)->nodeType == LIST) ||
	       (getMemRefChild(oldNode)->nodeType == FINALELLIPTICLIST)) &&
	      ((getMemRefChild(oldNode)->arguments != NULL) &&
	       (getMemRefChild(oldNode)->arguments->next != NULL) &&
	       (getMemRefChild(oldNode)->arguments->next != NULL))) {
	    newArgs = copyChainWithoutReversal(getMemRefChild(oldNode)->arguments->next, copyThingOnVoid);
	    newNode = allocateNode();
	    newNode->nodeType = getMemRefChild(oldNode)->nodeType;
	    newNode->arguments = newArgs;
	    newNode->argArray = NULL;
	    newNode->argArraySize = 0;
	    newNode->argArrayAllocSize = 0;
	    newNode = addMemRef(newNode);
	    if (oldChecked && (newNode->nodeType == MEMREF)) {
	      newNode->cache->isCorrectlyTyped = 1;
	    }
	    ((entry *) curr->value)->value = newNode;
	    freeThing(oldNode);
	    okay = 1;
	  } else {
	    newNode = deepCopyThing(oldNode);
	    if (((newNode->nodeType == LIST) ||
		 (newNode->nodeType == FINALELLIPTICLIST)) &&
		((getMemRefChild(oldNode)->arguments != NULL) &&
		 (getMemRefChild(oldNode)->arguments->next != NULL) &&
		 (getMemRefChild(oldNode)->arguments->next->next != NULL))) {
	      freeThing(oldNode);
	      freeThing((node *) (newNode->arguments->value));
	      newArgs = newNode->arguments->next;
	      safeFree(newNode->arguments);
	      newNode->arguments = newArgs;
	      newNode->argArray = NULL;
	      newNode->argArraySize = 0;
	      newNode->argArrayAllocSize = 0;
	      newNode = addMemRef(newNode);
	      /* We might want to re-use the re-usable memref'ed parts
		 of the original expression here 
	      */
	      if (oldChecked && (newNode->nodeType == MEMREF)) {
		newNode->cache->isCorrectlyTyped = 1;
	      }
	      ((entry *) curr->value)->value = newNode;
	      okay = 1;
	    } else {
	      freeThing(newNode);
	    }
	  }
	}
	break;
      case LIST:
      case FINALELLIPTICLIST:
	if ((oldNode->arguments != NULL) &&
	    (oldNode->arguments->next != NULL) &&
	    (oldNode->arguments->next->next != NULL)) {
	  freeThing((node *) (oldNode->arguments->value));
	  newArgs = oldNode->arguments->next;
	  safeFree(oldNode->arguments);
	  oldNode->arguments = newArgs;
	  if (oldNode->argArray != NULL) {
	    if (oldNode->argArraySize >= 2) {
	      oldNode->argArraySize--;
	    } else {
	      safeFree(oldNode->argArray);
	      oldNode->argArray = NULL;
	      oldNode->argArraySize = 0;
	      oldNode->argArrayAllocSize = 0;
	    }
	  }
	  okay = 1;
	} else {
	  okay = 0;
	}
	break;
      default:
	okay = 0;
	break;
      }
      break;
    }
    curr = curr->next;
  }

  return okay;
}

int performListTailOnDeclaredEntry(chain *declSymTbl, char *name) {
  chain *curr;

  curr = declSymTbl;
  while (curr != NULL) {
    if (containsEntry((chain *) (curr->value), name)) return performListTailOnEntry((chain *) (curr->value), name);
    curr = curr->next;
  }

  return 0;
}
