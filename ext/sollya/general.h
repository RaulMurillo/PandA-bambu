/*

  Copyright 2007-2018 by

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

#ifndef GENERAL_H
#define GENERAL_H


#include <mpfr.h>
#include "expression.h"
#include "chain.h"
#include "library.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "sollya-messaging.h"
#include "bitfields.h"

#define DEFAULTPOINTS 501
#define DEFAULTPRECISION 165
#define DEFAULTDIAM 1e-4
#define DEFAULTDIAM2 2e-3
#define DEFAULTTAYLORRECURSIONS 0
#define DEFAULTHOPITALRECURSIONS 2
#define MAXHORNERTREESIZE 5000
#define MAXHORNERDEGREE 5000
#define MAXAUTOSIMPLSIZE 5500
#define CHEAPSIMPLIFYSIZE 4000
#define PARIMEMSIZE 300000000

#define UNUSED_PARAM(_unused_param_x) ((void)(_unused_param_x))

#define PRINT_MODE_LEGACY            0
#define PRINT_MODE_WARNING_TO_STDERR 1
#define PRINT_MODE_WARNING_TO_FILE   2

extern int oldAutoPrint;
extern int oldVoidPrint;
extern int oldExternalProcedurePrint;
extern int printMode;
extern FILE *warnFile;
extern char *variablename;
extern bitfield suppressedMessages;
extern mp_prec_t defaultprecision;
extern mp_prec_t tools_precision;
extern int defaultpoints;
extern int taylorrecursions;
extern int dyadic;
extern int verbosity;
extern int activateMessageNumbers;
extern int canonical;
extern int fileNumber;
extern int autosimplify;
extern int timecounting;
extern chain *timeStack;
extern int fullParentheses;
extern int midpointMode;
extern int dieOnErrorMode;
extern int rationalMode;
extern int noRoundingWarnings;
extern int hopitalrecursions;
extern int eliminatePromptBackup;
extern int libraryMode;
extern void *scanner;
extern node *temp_node;
extern mpfr_t *mpfr_temp;
extern mpfr_t *mpfr_temp2;
extern node *minitree;
extern int miniparserCharactersRead;
extern int miniparserSemicolonAtEnd;
extern int miniparserEofReached;
extern char *constBuffer;
extern char *constBuffer2;
extern node* parsedThing;
extern int promptToBePrinted;
extern int lastWasSyntaxError;
extern int lastCorrectlyExecuted;
extern int helpNotFinished;
extern char *tempString;
extern char *tempString2;
extern node *tempNode;
extern int tempInteger;
extern chain *symbolTable;
extern chain *declaredSymbolTable;
extern chain *backtraceStack;
extern mpfr_t statediam;
extern node *parsedThingIntern;
extern int *tempIntPtr;
extern FILE *inputFile;
extern int inputFileOpened;
extern int flushOutput;
extern int executingExternalCode;

extern node *memRefChainStart;

extern void *rpl_malloc(size_t n);
extern void *rpl_realloc(void *, size_t n);

extern int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used;
extern int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized;
extern int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized;
extern int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized;
extern sollya_mpfi_t __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x;
extern sollya_mpfi_t __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y;
extern mpfr_t __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp;
extern int    __sparsePolynomialEvalMpfr_var_used;
extern int    __sparsePolynomialEvalMpfr_scratch_initialized;
extern mpfr_t __sparsePolynomialEvalMpfr_scratch;
extern int    __sparsePolynomialEvalMpfi_var_used;
extern int    __sparsePolynomialEvalMpfi_scratch_initialized;
extern sollya_mpfi_t __sparsePolynomialEvalMpfi_scratch;

struct __sollya_time_struct_t {
  int64_t seconds;
  int64_t nano_seconds;
};
typedef struct __sollya_time_struct_t sollya_time_t;

void safeFree(void *);

void printPrompt(void);
void recoverFromError(void);
void demaskString(char*, char*);
char *maskString(char *);
void *safeCalloc (size_t nmemb, size_t size);
void *safeMalloc (size_t size);
void *safeRealloc (void *ptr, size_t size);
int printMessage(int verb, int msgNum, const char *format, ...);
int messageHasEnoughVerbosityAndIsNotSuppressed(int verb, int msgNum);
int sollyaVfprintf(FILE *fd, const char *format, va_list varlist);
int sollyaPrintf(const char *format, ...);
int sollyaFprintf(FILE *fd, const char *format, ...);
char *mpfr_to_binary_str(mpfr_t x);
int removeSpaces(char *outbuf, char *inbuf);
int removeMidpointMode(char *outbuf, char *inbuf);
void cutMidpointStringIntoTwo(char *str1, char *str2, char *str);
void freeCounter(void);
void pushTimeCounter(void);
void popTimeCounter(char *s);
void newReadFileStarted();
void carriageReturnLexed();
void newTokenLexed();
void considerDyingOnError();
void restartTool();
void normalMode();
void redMode();
void blueMode();
void blinkMode();
void unblinkMode();
void parseMode();
void outputMode();
void warningMode();
void saveMode();
void restoreMode();
void changeToWarningMode();
int  getDisplayColor();
void setDisplayColor(int);
int installMessageCallback(int (*msgHandler) (sollya_msg_t, void *), void *);
int uninstallMessageCallback();
void getMessageCallback(int (**)(sollya_msg_t, void *), void **);
int initializeLibraryMode(void *(*)(size_t),
			  void *(*)(size_t, size_t),
			  void *(*)(void *, size_t),
			  void (*)(void*),
			  void *(*)(void *, size_t, size_t),
			  void (*)(void *, size_t),
			  int,
			  char **,
			  void (*)(void *(*)(size_t),
				   void *(*)(void *, size_t, size_t),
				   void (*)(void *, size_t)),
			  void (*)(void *(**)(size_t),
				   void *(**)(void *, size_t, size_t),
				   void (**)(void *, size_t)));
int finalizeLibraryMode();
mp_prec_t getToolPrecision();
void setToolPrecision(mp_prec_t prec);
int getToolPoints();
void setToolPoints(int p);
int getToolTaylorRecursions();
void setToolTaylorRecursions(int i);
int getToolHopitalRecursions();
void setToolHopitalRecursions(int i);
int getToolDiameter(mpfr_t rop);
void setToolDiameter(mpfr_t op);
char *getTempDir();
char *getUniqueId();
int getDisplayMode();
int setDisplayMode(int);
int getVerbosity();
int setVerbosity(int);
int getCanonical();
void setCanonical(int);
int getAutosimplify();
void setAutosimplify(int);
int getFullParentheses();
void setFullParentheses(int);
int getMidpointMode();
void setMidpointMode(int);
int getDieOnErrorMode();
void setDieOnErrorMode(int);
int getTimecounting();
void setTimecounting(int);
int getRoundingWarnings();
void setRoundingWarnings(int);
int getRationalMode();
void setRationalMode(int);
double sollya_mpfr_get_d(mpfr_srcptr, mpfr_rnd_t);
sollya_mpfi_t *getReusedGlobalMPFIVars(unsigned int, mp_prec_t);
void returnReusedGlobalMPIVars(unsigned int);
mpfr_t *getReusedGlobalMPFRVars(unsigned int, mp_prec_t);
void returnReusedGlobalMPFRVars(unsigned int);
int sollyaLibPrintmessage(int, int, const char *, ...);
int sollya_getc(FILE *);
int sollya_feof(FILE *);
int sollya_ferror(FILE *);
size_t sollya_fread(void *, size_t, size_t, FILE *);
size_t sollya_fwrite(const void *, size_t, size_t, FILE *);
int sollya_gettime(sollya_time_t *);
void readManipulate();
void parserFlushInput();
void enterExternalCode();
void leaveExternalCode();
char *getGnuplotName();

/* Compatibility macros */

#if defined(HAVE_MEMMOVE) && HAVE_MEMMOVE
#define sollya_memmove(dest, src, n) memmove((dest),(src),(n))
#else
void *sollya_memmove_impl(void *, const void *, size_t);
#define sollya_memmove(dest, src, n) sollya_memmove_impl((dest),(src),(n))
#endif

#if defined(HAVE_MEMSET) && HAVE_MEMSET
#define sollya_memset(s, c, n) memset((s),(c),(n))
#else
void *sollya_memset_impl(void *, int, size_t);
#define sollya_memset(s, c, n) sollya_memset_impl((s),(c),(n))
#endif

#if defined(HAVE_STRCHR) && HAVE_STRCHR
#define sollya_strchr(s, c) strchr((s),(c))
#else
char *sollya_strchr_impl(const char *, int);
#define sollya_strchr(s, c) sollya_strchr_impl((s),(c))
#endif

#if defined(HAVE_STRRCHR) && HAVE_STRRCHR
#define sollya_strrchr(s, c) strrchr((s),(c))
#else
char *sollya_strrchr_impl(const char *, int);
#define sollya_strrchr(s, c) sollya_strrchr_impl((s),(c))
#endif

#if defined(HAVE_STRTOL) && HAVE_STRTOL
#define sollya_strtol(nptr, endptr, base) strtol((nptr),(endptr),(base))
#else
long int sollya_strtol_impl(const char *, char **, int);
#define sollya_strtol(nptr, endptr, base) sollya_strtol_impl((nptr),(endptr),(base))
#endif

#if defined(HAVE_DUP2) && HAVE_DUP2
#define sollya_dup2(oldfd, newfd) dup2((oldfd), (newfd))
#else
int sollya_dup2_impl(int, int);
#define sollya_dup2(oldfd, newfd) sollya_dup2_impl((oldfd), (newfd))
#endif

#if defined(HAVE_STRSTR) && HAVE_STRSTR
#define sollya_strstr(haystack, needle) strstr((haystack),(needle))
#else
char *sollya_strstr_impl(const char *, const char *);
#define sollya_strstr(haystack, needle) sollya_strstr_impl((haystack),(needle))
#endif

#endif /* ifdef GENERAL_H*/
