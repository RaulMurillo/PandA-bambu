/*

Copyright 2007-2017 by

Laboratoire de l'Informatique du Parallelisme,
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2)

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

%{

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "expression.h"
#include "internparser.h"
#include "general.h"
#include "chain.h"

#define YY_NO_UNPUT 1
#define YY_NO_INPUT 1

/* Mess with the mallocs used by the lexer */
extern void *parserCalloc(size_t, size_t);
extern void *parserMalloc(size_t);
extern void *parserRealloc(void *, size_t);
extern void parserFree(void *);
#undef malloc
#undef realloc
#undef calloc
#undef free
#define malloc parserMalloc
#define realloc parserRealloc
#define calloc parserCalloc
#define free UNUSED_PARAM(yyscanner); parserFree
/* End of the malloc mess */

/* Mess with fprintf to use sollyaFprintf instead. */
#undef  fprintf
#define fprintf UNUSED_PARAM(yyscanner); sollyaFprintf
/* End of fprintf mess */

/* Mess with getc... and fread to use sollya_getc, sollya_feof, sollya_ferror, sollya_fwrite and sollya_fread instead. */
#undef  getc
#undef  fread
#undef  feof
#undef  ferror
#undef  fwrite
#define getc sollya_getc
#define fread sollya_fread
#define feof sollya_feof
#define ferror sollya_ferror
#define fwrite sollya_fwrite
/* End of getc, feof, ferror, fwrite and fread mess */

/* const mess */
#if defined(const)
#undef const
#endif
/* End of the const mess */

%}

%option noyyalloc noyyrealloc noyyfree

%option noyywrap
%option always-interactive
%option prefix="internyy"
%option reentrant
%option bison-bridge
%option nounput

%x commentstate

CHAR		[a-zA-Z]
NUMBER		[0-9]
HEXNUMBER       (([0-9])|([ABCDEFabcdef]))
SPECIAL         [_$§]

ACONSTANT        ((({NUMBER}+|({NUMBER}*"."{NUMBER}+))(((((" ")*)(([eE]))([+-])?{NUMBER}+)?)))|("infty"|"nan"|"NAN"|"NaN"|"@Inf@"|"@NaN@"))
CONSTANT        ({ACONSTANT})|((("%")({NUMBER})({NUMBER}+)("%"))({ACONSTANT}))
DYADICCONSTANT  ({NUMBER}+)([bB])([+-]?)({NUMBER}+)
HEXCONSTANT     ("0x")({HEXNUMBER}{16}|{HEXNUMBER}{8})

HEXADECIMALCONSTANT ("0x"|"0X")(({HEXNUMBER}*"."{HEXNUMBER}+)|({HEXNUMBER}+((".")?)))(("p"|"P")(([+-])?)({NUMBER}*))

BINARYCONSTANT  (([0-1])+|(([0-1])*"."([0-1])+))"_2"

AMIDPOINTCONSTANT        (({NUMBER}+|({NUMBER}*"."{NUMBER}+))(("~"({NUMBER})+("/")({NUMBER})+("~")))(((((" ")*)(([eE]))([+-])?{NUMBER}+)?)))
MIDPOINTCONSTANT         (({AMIDPOINTCONSTANT})|(("-")((" ")*)({AMIDPOINTCONSTANT})))

IDENTIFIER        ({CHAR}|{SPECIAL})({CHAR}|{NUMBER}|{SPECIAL})*

DOT             "."

LPAR            "("
RPAR            ")"

LBRACKET        "["
RBRACKET        "]"

PI              ("pi")|("Pi")

PLUS            "+"
MINUS           "-"
APPROX          "~"
MUL             "*"
DIV             "/"
POW             "^"
SQRT            "sqrt"
EXP             "exp"
FREEVAR         "_x_"
LOG             "log"
LOG2            "log2"
LOG10           "log10"
SIN             "sin"
COS             "cos"
TAN             "tan"
ASIN            "asin"
ACOS            "acos"
ATAN            "atan"
SINH            "sinh"
COSH            "cosh"
TANH            "tanh"
ASINH           "asinh"
ACOSH           "acosh"
ATANH           "atanh"
ABS             "abs"
ERF             "erf"
ERFC            "erfc"
LOG1P           "log1p"
EXPM1           "expm1"

EQUAL           "="
COMMA           ","
PREC            "prec"
POINTS          "points"
EXCLAMATION     "!"
DOUBLECOLON     "::"
COLON           ":"
DOTCOLON        ".:"
COLONDOT        ":."
ASSIGNEQUAL     ":="

SEMICOLON       ";"
QUIT            "quit"
PRINT           "print"
SUPPRESSMESSAGE "suppressmessage"
UNSUPPRESSMESSAGE "unsuppressmessage"
DIFF            "diff"
DIRTYSIMPLIFY   "dirtysimplify"
CANONICAL       "canonical"
PLOT            "plot"
REMEZ           "remez"
ANNOTATEFUNCTION "annotatefunction"
FPMINIMAX       "fpminimax"
INFNORM         "infnorm"
SUPNORM         "supnorm"
DIAM            "diam"

DOUBLELONG       "double"
SINGLELONG       "single"
QUADLONG         "quad"
HALFPRECISIONLONG "halfprecision"
DOUBLEDOUBLELONG "doubledouble"
TRIPLEDOUBLELONG "tripledouble"
DOUBLEEXTENDEDLONG "doubleextended"


DOUBLESHORT      "D"
SINGLESHORT      "SG"
QUADSHORT        "QD"
HALFPRECISIONSHORT "HP"
DOUBLEDOUBLESHORT "DD"
TRIPLEDOUBLESHORT "TD"
DOUBLEEXTENDEDSHORT "DE"


DOUBLE          ({DOUBLELONG}|{DOUBLESHORT})
SINGLE          ({SINGLELONG}|{SINGLESHORT})
QUAD            ({QUADLONG}|{QUADSHORT})
HALFPRECISION   ({HALFPRECISIONLONG}|{HALFPRECISIONSHORT})
DOUBLEDOUBLE    ({DOUBLEDOUBLELONG}|{DOUBLEDOUBLESHORT})
TRIPLEDOUBLE    ({TRIPLEDOUBLELONG}|{TRIPLEDOUBLESHORT})
DOUBLEEXTENDED  ({DOUBLEEXTENDEDLONG}|{DOUBLEEXTENDEDSHORT})

HORNER          "horner"
DEGREE          "degree"
EXPAND          "expand"

SIMPLIFYSAFE1   "simplify"
SIMPLIFYSAFE2   "simplifysafe"

SIMPLIFYSAFE    ({SIMPLIFYSAFE1}|{SIMPLIFYSAFE2})

TAYLOR          "taylor"
TAYLORFORM      "taylorform"
CHEBYSHEVFORM   "chebyshevform"
AUTODIFF        "autodiff"

FINDZEROS       "findzeros"
FPFINDZEROS     "fpfindzeros"
DIRTYINFNORM    "dirtyinfnorm"
GCD             "gcd"
EUCLDIV         "div"
EUCLMOD         "mod"
NUMBERROOTS     "numberroots"

EVALUATE        "evaluate"


NUMERATOR       "numerator"
DENOMINATOR     "denominator"

DYADIC          "dyadic"
DISPLAY         "display"
ON              "on"
OFF             "off"
POWERS          "powers"
BINARY          "binary"
HEXADECIMAL     "hexadecimal"

INTEGRAL        "integral"
DIRTYINTEGRAL   "dirtyintegral"

STRINGDELIMITER [\"]
OCTALCHAR       [01234567]
OCTAL           ({OCTALCHAR})(({OCTALCHAR})?)(({OCTALCHAR})?)
HEXACHAR        [0123456789ABCDEFabcdef]
HEXA            ({HEXACHAR})(({HEXACHAR})?)
STRING          ({STRINGDELIMITER})((("\\\\")|(("\\")[\"\'\?ntabfrv])|(("\\")({OCTAL}))|(("\\x")({HEXA}))|([^\"\\]))*)({STRINGDELIMITER})


VERBOSITY       "verbosity"
SHOWMESSAGENUMBERS "showmessagenumbers"

WORSTCASE       "worstcase"

LEFTANGLE       "<"
RIGHTANGLEUNDERSCORE ">_"
RIGHTANGLEDOT ">."
RIGHTANGLE      ">"
STARLEFTANGLE   "*<"
RIGHTANGLESTAR  ">*"

COMPAREEQUAL    "=="
EXCLAMATIONEQUAL "!="

AND             "&&"
OR              "||"

SUBSTITUTE      "substitute"
COMPOSEPOLYNOMIALS "composepolynomials"

DOTS            "..."

IMPLEMENTPOLY   "implementpoly"
IMPLEMENTCONST  "implementconstant"

CHECKINFNORM    "checkinfnorm"

TAYLORRECURSIONS "taylorrecursions"

PRINTHEXA       ("printhexa"|"printdouble")
PRINTFLOAT      ("printfloat"|"printsingle")
PRINTBINARY     "printbinary"

ROUNDCOEFFICIENTS "roundcoefficients"


RESTART         "restart"


ZERODENOMINATORS "zerodenominators"
ISEVALUABLE     "isevaluable"
HONORCOEFFPREC  "honorcoeffprec"

ACCURATEINFNORM  "accurateinfnorm"


FILE            "file"
POSTSCRIPT      "postscript"
POSTSCRIPTFILE  "postscriptfile"

PRINTEXPANSION  "printexpansion"

BASHEXECUTE     "bashexecute"
BASHEVALUATE    "bashevaluate"
GETSUPPRESSEDMESSAGES "getsuppressedmessages"
GETBACKTRACE    "getbacktrace"
EXTERNALPLOT    "externalplot"
PERTURB         "perturb"

COEFF           "coeff"
SUBPOLY         "subpoly"

QUESTIONMARK    "?"

SEARCHGAL       "searchgal"


RATIONALAPPROX  "rationalapprox"

READ            "read"

COMMENTSTART    "/*"
COMMENTEND      "*/"

ONELINECOMMENT  ("//")([^\n])*"\n"

WRITE           "write"

ASCIIPLOT       "asciiplot"


ROUNDTOFORMAT   "round"
MINUSWORD       "RD"
PLUSWORD        "RU"
ZEROWORD        "RZ"
NEAREST         "RN"


GUESSDEGREE     "guessdegree"

PARSE           "parse"

AUTOSIMPLIFY    "autosimplify"

TIMING          "timing"
TIME            "time"
FULLPARENTHESES "fullparentheses"
MIDPOINTMODE    "midpointmode"
DIEONERRORMODE  "dieonerrormode"
SUPPRESSWARNINGS "roundingwarnings"
RATIONALMODE    "rationalmode"

LIBRARY         "library"
LIBRARYCONSTANT "libraryconstant"

HOPITALRECURSIONS "hopitalrecursions"

HELP            "help"

DIRTYFINDZEROS  "dirtyfindzeros"

CEIL            "ceil"
FLOOR           "floor"
NEARESTINT      "nearestint"

HEAD            "head"
TAIL            "tail"

VERTBAR         "|"
AT              "@"

IF              "if"
THEN            "then"
ELSE            "else"
FOR             "for"
IN              "in"
FROM            "from"
TO              "to"
BY              "by"
DO              "do"
BEGIN           "begin"
END             "end"
LBRACE          "{"
RBRACE          "}"

WHILEDEF        "while"

TRUE            "true"
FALSE           "false"
DEFAULT         "default"

MATCH           "match"
WITH            "with"

RENAME          "rename"
OBJECTNAME      "objectname"
BIND            "bind"

LENGTH          "length"

ABSOLUTE        "absolute"
RELATIVE        "relative"
DECIMAL         "decimal"
FIXED           "fixed"
FLOATING        "floating"

ERROR           "error"

PRINTXML        "printxml"

INF             "inf"
MID             "mid"
SUP             "sup"

MIN             "min"
MAX             "max"

READXML         "readxml"

MANTISSA        "mantissa"
EXPONENT        "exponent"
PRECISION       "precision"

REVERT          "revert"
SORT            "sort"

READFILE        "readfile"

ROUNDCORRECTLY  "roundcorrectly"

EXECUTE         "execute"

ISBOUND         "isbound"

EXTERNALPROC    "externalproc"
VOID            "void"
CONSTANTTYPE    "constant"
FUNCTION        "function"
OBJECT          "object"
RANGE           "range"
INTEGER         "integer"
STRINGTYPE      "string"
BOOLEAN         "boolean"
LIST            "list"
OF              "of"

VAR             "var"

PROC            "proc"
PROCEDURE       "procedure"
RETURN          "return"
NOP             "nop"


%%

%{

%}




{COMMENTSTART}  {     BEGIN(commentstate); }

<commentstate>{COMMENTEND} { BEGIN(INITIAL); }

<commentstate>. { /* Eat up comments */

                 }

<commentstate>[\n] { /* Eat up newlines in comments */

		}



{ONELINECOMMENT} {  /* Eat up comments */

                 }

<<EOF>>         {
                      yyterminate();
                }
{MIDPOINTCONSTANT}                          {
					      constBuffer2 = (char *) safeCalloc(yyleng+1,sizeof(char));
					      if (removeSpaces(constBuffer2,yytext)) {
						printMessage(2,SOLLYA_MSG_SPACES_REMOVED_FROM_CONSTANT_IN_SCIENTIF_NOTAT,"Information: removed spaces in scientific notation constant \"%s\", it will be considered as \"%s\"\n",yytext,constBuffer2);
					      }
					      yylval->value = constBuffer2;
                                              return MIDPOINTCONSTANTTOKEN; }
{CONSTANT}                                  {
					      constBuffer2 = (char *) safeCalloc(yyleng+1,sizeof(char));
					      if (removeSpaces(constBuffer2,yytext)) {
						printMessage(2,SOLLYA_MSG_SPACES_REMOVED_FROM_CONSTANT_IN_SCIENTIF_NOTAT,"Information: removed spaces in scientific notation constant \"%s\", it will be considered as \"%s\"\n",yytext,constBuffer2);
					      }
					      yylval->value = constBuffer2;
                                               return CONSTANTTOKEN; }
{DYADICCONSTANT} 			    {
                                              constBuffer = (char *) safeCalloc(yyleng+1,sizeof(char));
					      strncpy(constBuffer,yytext,yyleng);
					      yylval->value = constBuffer;
                                               return DYADICCONSTANTTOKEN; }
{HEXADECIMALCONSTANT}     		    { constBuffer = (char *) safeCalloc(yyleng+1,sizeof(char));
					      strncpy(constBuffer,yytext,yyleng);
					      yylval->value = constBuffer;
                                              return HEXADECIMALCONSTANTTOKEN; }
{HEXCONSTANT}     			    { constBuffer = (char *) safeCalloc(yyleng+1,sizeof(char));
					      strncpy(constBuffer,yytext,yyleng);
					      yylval->value = constBuffer;
                                               return HEXCONSTANTTOKEN; }
{BINARYCONSTANT}  			    { constBuffer = (char *) safeCalloc(yyleng-1,sizeof(char));
					      strncpy(constBuffer,yytext,yyleng-2);
					      yylval->value = constBuffer;
                                               return BINARYCONSTANTTOKEN; }

{PI}              			    {  return PITOKEN; }

{STRING}          			    {
					      constBuffer = (char *) safeCalloc(yyleng-1,sizeof(char));
					      constBuffer2 = (char *) safeCalloc(yyleng-1,sizeof(char));
					      strncpy(constBuffer2,yytext+1,yyleng-2);
					      demaskString(constBuffer,constBuffer2);
					      safeFree(constBuffer2);
					      yylval->value = constBuffer;
                                               return STRINGTOKEN; }

{LPAR}            			    {  return LPARTOKEN; }
{RPAR}            			    {  return RPARTOKEN; }
{LBRACKET}        			    {  return LBRACKETTOKEN; }
{RBRACKET}        			    {  return RBRACKETTOKEN; }
{EXCLAMATIONEQUAL}			    {  return EXCLAMATIONEQUALTOKEN; }
{COMPAREEQUAL}                              {  return COMPAREEQUALTOKEN; }
{EQUAL}           			    {  return EQUALTOKEN; }
{ASSIGNEQUAL}          			    {  return ASSIGNEQUALTOKEN; }
{COMMA}           			    {  return COMMATOKEN; }
{EXCLAMATION}				    {  return EXCLAMATIONTOKEN; }
{SEMICOLON}       			    {  return SEMICOLONTOKEN; }
{LEFTANGLE}       			    {  return LEFTANGLETOKEN; }
{STARLEFTANGLE}   			    {  return STARLEFTANGLETOKEN; }
{RIGHTANGLEUNDERSCORE}			    {  return RIGHTANGLEUNDERSCORETOKEN; }
{RIGHTANGLEDOT} 			    {  return RIGHTANGLEDOTTOKEN; }
{RIGHTANGLESTAR}			    {  return RIGHTANGLESTARTOKEN; }
{RIGHTANGLE}      			    {  return RIGHTANGLETOKEN; }
{DOTS}            			    {  return DOTSTOKEN; }
{DOT}            			    {  return DOTTOKEN; }
{QUESTIONMARK}				    {  return QUESTIONMARKTOKEN; }
{AND}         				    {  return ANDTOKEN; }
{OR}         				    {  return ORTOKEN; }
{VERTBAR}				    {  return VERTBARTOKEN; }
{VAR}					    {  return VARTOKEN; }
{PROC}					    {  return PROCTOKEN; }
{PROCEDURE}				    {  return PROCEDURETOKEN; }
{RETURN}				    {  return RETURNTOKEN; }
{NOP}	         			    {  return NOPTOKEN; }
{AT}					    {  return ATTOKEN; }
{DOUBLECOLON}				    {  return DOUBLECOLONTOKEN; }
{COLON}  				    {  return COLONTOKEN; }
{DOTCOLON}				    {  return DOTCOLONTOKEN; }
{COLONDOT}				    {  return COLONDOTTOKEN; }


{PLUS}            			    {  return PLUSTOKEN; }
{MINUS}           			    {  return MINUSTOKEN; }
{MUL}          				    {  return MULTOKEN; }
{DIV}            			    {  return DIVTOKEN; }
{POW}            			    {  return POWTOKEN; }

{SQRT}            			    {  return SQRTTOKEN; }
{EXP}             			    {  return EXPTOKEN; }
{FREEVAR}             			    {  return FREEVARTOKEN; }
{LOG}             			    {  return LOGTOKEN; }
{LOG2}            			    {  return LOG2TOKEN; }
{LOG10}           			    {  return LOG10TOKEN; }
{SIN}             			    {  return SINTOKEN; }
{COS}             			    {  return COSTOKEN; }
{TAN}             			    {  return TANTOKEN; }
{ASIN}            			    {  return ASINTOKEN; }
{ACOS}            			    {  return ACOSTOKEN; }
{ATAN}            			    {  return ATANTOKEN; }
{SINH}            			    {  return SINHTOKEN; }
{COSH}            			    {  return COSHTOKEN; }
{TANH}            			    {  return TANHTOKEN; }
{ASINH}           			    {  return ASINHTOKEN; }
{ACOSH}           			    {  return ACOSHTOKEN; }
{ATANH}           			    {  return ATANHTOKEN; }
{ABS}             			    {  return ABSTOKEN; }
{ERF}             			    {  return ERFTOKEN; }
{ERFC}            			    {  return ERFCTOKEN; }
{LOG1P}           			    {  return LOG1PTOKEN; }
{EXPM1}           			    {  return EXPM1TOKEN; }
{DOUBLE}          			    {  return DOUBLETOKEN; }
{SINGLE}          			    {  return SINGLETOKEN; }
{QUAD}          			    {  return QUADTOKEN; }
{HALFPRECISION}          		    {  return HALFPRECISIONTOKEN; }
{DOUBLEDOUBLE}				    {  return DOUBLEDOUBLETOKEN; }
{TRIPLEDOUBLE}    			    {  return TRIPLEDOUBLETOKEN; }
{DOUBLEEXTENDED}  			    {  return DOUBLEEXTENDEDTOKEN; }
{CEIL}            			    {  return CEILTOKEN; }
{FLOOR}           			    {  return FLOORTOKEN; }
{NEARESTINT}           			    {  return NEARESTINTTOKEN; }

{PREC}            			    {  return PRECTOKEN; }
{POINTS}          			    {  return POINTSTOKEN; }
{DIAM}            			    {  return DIAMTOKEN; }
{DISPLAY}          			    {  return DISPLAYTOKEN; }
{VERBOSITY}       			    {  return VERBOSITYTOKEN; }
{SHOWMESSAGENUMBERS}       		    {  return SHOWMESSAGENUMBERSTOKEN; }
{CANONICAL}       			    {  return CANONICALTOKEN; }
{AUTOSIMPLIFY}    			    {  return AUTOSIMPLIFYTOKEN; }
{TAYLORRECURSIONS}			    {  return TAYLORRECURSIONSTOKEN; }
{TIMING}          			    {  return TIMINGTOKEN; }
{TIME}          			    {  return TIMETOKEN; }
{FULLPARENTHESES} 			    {  return FULLPARENTHESESTOKEN; }
{MIDPOINTMODE}    			    {  return MIDPOINTMODETOKEN; }
{DIEONERRORMODE}    			    {  return DIEONERRORMODETOKEN; }
{SUPPRESSWARNINGS}			    {  return SUPPRESSWARNINGSTOKEN; }
{RATIONALMODE}			            {  return RATIONALMODETOKEN; }
{HOPITALRECURSIONS}			    {  return HOPITALRECURSIONSTOKEN; }

{ON}              			    {  return ONTOKEN; }
{OFF}             			    {  return OFFTOKEN; }
{DYADIC}				    {  return DYADICTOKEN; }
{POWERS}          			    {  return POWERSTOKEN; }
{BINARY}          			    {  return BINARYTOKEN; }
{HEXADECIMAL}          			    {  return HEXADECIMALTOKEN; }
{FILE}            			    {  return FILETOKEN; }
{POSTSCRIPT}      			    {  return POSTSCRIPTTOKEN; }
{POSTSCRIPTFILE}  			    {  return POSTSCRIPTFILETOKEN; }
{PERTURB}         			    {  return PERTURBTOKEN; }
{MINUSWORD}       			    {  return MINUSWORDTOKEN; }
{PLUSWORD}        			    {  return PLUSWORDTOKEN; }
{ZEROWORD}        			    {  return ZEROWORDTOKEN; }
{NEAREST}         			    {  return NEARESTTOKEN; }
{HONORCOEFFPREC} 			    {  return HONORCOEFFPRECTOKEN; }
{TRUE}					    {  return TRUETOKEN; }
{FALSE}					    {  return FALSETOKEN; }
{DEFAULT}				    {  return DEFAULTTOKEN; }
{MATCH} 				    {  return MATCHTOKEN; }
{WITH}  				    {  return WITHTOKEN; }
{HEAD}   				    {  return HEADTOKEN; }
{ROUNDCORRECTLY}   		       	    {  return ROUNDCORRECTLYTOKEN; }
{REVERT}   				    {  return REVERTTOKEN; }
{SORT}   				    {  return SORTTOKEN; }
{MANTISSA}   				    {  return MANTISSATOKEN; }
{EXPONENT}   				    {  return EXPONENTTOKEN; }
{PRECISION}   				    {  return PRECISIONTOKEN; }
{TAIL}   				    {  return TAILTOKEN; }
{INF}   				    {  return INFTOKEN; }
{MID}   				    {  return MIDTOKEN; }
{SUP}   				    {  return SUPTOKEN; }
{MIN}   				    {  return MINTOKEN; }
{MAX}   				    {  return MAXTOKEN; }
{READXML}   				    {  return READXMLTOKEN; }
{LENGTH}   				    {  return LENGTHTOKEN; }
{ABSOLUTE}   				    {  return ABSOLUTETOKEN; }
{RELATIVE}   				    {  return RELATIVETOKEN; }
{DECIMAL}   				    {  return DECIMALTOKEN; }
{FIXED}   				    {  return FIXEDTOKEN; }
{FLOATING}   				    {  return FLOATINGTOKEN; }

{ERROR}   				    {  return ERRORTOKEN; }
{READFILE}   				    {  return READFILETOKEN; }
{ISBOUND}   				    {  return ISBOUNDTOKEN; }

{QUIT}             			    {
			                        return FALSEQUITTOKEN;
                                            }
{RESTART}         			    { return FALSERESTARTTOKEN; }
{LIBRARY}         			    {  return LIBRARYTOKEN; }
{LIBRARYCONSTANT}         		    {  return LIBRARYCONSTANTTOKEN; }
{DIFF}            			    {  return DIFFTOKEN; }
{DIRTYSIMPLIFY}				    {  return DIRTYSIMPLIFYTOKEN; }
{REMEZ}           			    {  return REMEZTOKEN; }
{ANNOTATEFUNCTION}           		    {  return ANNOTATEFUNCTIONTOKEN; }
{FPMINIMAX}           			    {  return FPMINIMAXTOKEN; }
{HORNER}          			    {  return HORNERTOKEN; }
{EXPAND}          			    {  return EXPANDTOKEN; }
{SIMPLIFYSAFE}				    {  return SIMPLIFYSAFETOKEN; }
{TAYLOR}         			    {  return TAYLORTOKEN; }
{TAYLORFORM}         			    {  return TAYLORFORMTOKEN; }
{CHEBYSHEVFORM}         		    {  return CHEBYSHEVFORMTOKEN; }
{AUTODIFF}         			    {  return AUTODIFFTOKEN; }
{DEGREE}          			    {  return DEGREETOKEN; }
{NUMERATOR}       			    {  return NUMERATORTOKEN; }
{DENOMINATOR}     			    {  return DENOMINATORTOKEN; }
{SUBSTITUTE}      			    {  return SUBSTITUTETOKEN; }
{COMPOSEPOLYNOMIALS}      		    {  return COMPOSEPOLYNOMIALSTOKEN; }
{COEFF}           			    {  return COEFFTOKEN; }
{SUBPOLY}         			    {  return SUBPOLYTOKEN; }
{ROUNDCOEFFICIENTS} 			    {  return ROUNDCOEFFICIENTSTOKEN; }
{RATIONALAPPROX}  			    {  return RATIONALAPPROXTOKEN; }
{APPROX}           			    {  return APPROXTOKEN; }
{ACCURATEINFNORM}  			    {  return ACCURATEINFNORMTOKEN; }
{ROUNDTOFORMAT}   			    {  return ROUNDTOFORMATTOKEN; }
{EVALUATE}        			    {  return EVALUATETOKEN; }
{PARSE}           			    {  return PARSETOKEN; }
{PRINT}           			    {  return PRINTTOKEN; }
{SUPPRESSMESSAGE}           		    {  return SUPPRESSMESSAGETOKEN; }
{UNSUPPRESSMESSAGE}           		    {  return UNSUPPRESSMESSAGETOKEN; }
{PRINTXML}           			    {  return PRINTXMLTOKEN; }
{EXECUTE}           			    {  return EXECUTETOKEN; }
{PLOT}            			    {  return PLOTTOKEN; }
{PRINTHEXA}       			    {  return PRINTHEXATOKEN; }
{PRINTFLOAT}       			    {  return PRINTFLOATTOKEN; }
{PRINTBINARY}     			    {  return PRINTBINARYTOKEN; }
{PRINTEXPANSION}  			    {  return PRINTEXPANSIONTOKEN; }
{BASHEXECUTE}     			    {  return BASHEXECUTETOKEN; }
{BASHEVALUATE}     			    {  return BASHEVALUATETOKEN; }
{GETSUPPRESSEDMESSAGES}     		    {  return GETSUPPRESSEDMESSAGESTOKEN; }
{GETBACKTRACE}     		            {  return GETBACKTRACETOKEN; }
{EXTERNALPLOT}    			    {  return EXTERNALPLOTTOKEN; }
{WRITE}           			    {  return WRITETOKEN; }
{ASCIIPLOT}       			    {  return ASCIIPLOTTOKEN; }
{RENAME}         			    {  return RENAMETOKEN; }
{OBJECTNAME}         			    {  return OBJECTNAMETOKEN; }
{BIND}         			            {  return BINDTOKEN; }


{INFNORM}         			    {  return INFNORMTOKEN; }
{SUPNORM}         			    {  return SUPNORMTOKEN; }
{FINDZEROS}       			    {  return FINDZEROSTOKEN; }
{FPFINDZEROS}     			    {  return FPFINDZEROSTOKEN; }
{DIRTYINFNORM}    			    {  return DIRTYINFNORMTOKEN; }
{GCD}                  			    {  return GCDTOKEN; }
{EUCLDIV}                  		    {  return EUCLDIVTOKEN; }
{EUCLMOD}                  		    {  return EUCLMODTOKEN; }
{NUMBERROOTS}    			    {  return NUMBERROOTSTOKEN; }
{INTEGRAL}        			    {  return INTEGRALTOKEN; }
{DIRTYINTEGRAL}				    {  return DIRTYINTEGRALTOKEN; }
{WORSTCASE}       			    {  return WORSTCASETOKEN; }
{IMPLEMENTPOLY}				    {  return IMPLEMENTPOLYTOKEN; }
{IMPLEMENTCONST}                            {  return IMPLEMENTCONSTTOKEN; }
{CHECKINFNORM}    			    {  return CHECKINFNORMTOKEN; }
{ZERODENOMINATORS}			    {  return ZERODENOMINATORSTOKEN; }
{ISEVALUABLE}     			    {  return ISEVALUABLETOKEN; }
{SEARCHGAL}       			    {  return SEARCHGALTOKEN; }
{GUESSDEGREE}     			    {  return GUESSDEGREETOKEN; }
{DIRTYFINDZEROS}  			    {  return DIRTYFINDZEROSTOKEN; }

{IF}					    {  return IFTOKEN; }
{THEN}					    {  return THENTOKEN; }
{ELSE}					    {  return ELSETOKEN; }
{FOR}					    {  return FORTOKEN; }
{IN}					    {  return INTOKEN; }
{FROM}					    {  return FROMTOKEN; }
{TO}					    {  return TOTOKEN; }
{BY}					    {  return BYTOKEN; }
{DO}					    {  return DOTOKEN; }
{BEGIN}					    {  return BEGINTOKEN; }
{END}					    {  return ENDTOKEN; }
{LBRACE}				    {  return LEFTCURLYBRACETOKEN; }
{RBRACE}				    {  return RIGHTCURLYBRACETOKEN; }
{WHILEDEF}				    {  return WHILETOKEN; }

{EXTERNALPROC}                              {  return EXTERNALPROCTOKEN; }
{VOID}                                      {  return VOIDTOKEN; }
{CONSTANTTYPE}                              {  return CONSTANTTYPETOKEN; }
{FUNCTION}                                  {  return FUNCTIONTOKEN; }
{OBJECT}                                    {  return OBJECTTOKEN; }
{RANGE}                                     {  return RANGETOKEN; }
{INTEGER}                                   {  return INTEGERTOKEN; }
{STRINGTYPE}                                {  return STRINGTYPETOKEN; }
{BOOLEAN}                                   {  return BOOLEANTOKEN; }
{LIST}                                      {  return LISTTOKEN; }
{OF}                                        {  return OFTOKEN; }



[ \t]		{ /* Eat up spaces and tabulators */
		}

[\r]            { /* Eat up carriage return for compatibility with DOS-like systems */
                }

[\n]		{

		}

{IDENTIFIER}        			    { constBuffer = (char *) safeCalloc(yyleng+1,sizeof(char));
					      strncpy(constBuffer,yytext,yyleng);
					      yylval->value = constBuffer;
					       return IDENTIFIERTOKEN; }


.		{ /* otherwise */
			printMessage(1,SOLLYA_MSG_A_CHARACTER_COULD_NOT_BE_RECOGNIZED,"The character \"%s\" cannot be recognized. Will ignore it.\n",
				yytext);
}

%%


void *yyalloc (yy_size_t  size , yyscan_t yyscanner)
{
        UNUSED_PARAM(yyscanner);
        return (void *) malloc( size );
}

void *yyrealloc  (void * ptr, yy_size_t  size , yyscan_t yyscanner)
{
        /* The cast to (char *) in the following accommodates both
         * implementations that use char* generic pointers, and those
         * that use void* generic pointers.  It works with the latter
         * because both ANSI C and C++ allow castless assignment from
         * any pointer type to void*, and deal with argument conversions
         * as though doing an assignment.
         */
        UNUSED_PARAM(yyscanner);
        return (void *) realloc( (char *) ptr, size );
}

void yyfree (void * ptr , yyscan_t yyscanner)
{
        UNUSED_PARAM(yyscanner);
        free( (char *) ptr );   /* see yyrealloc() for (char *) cast */
}
