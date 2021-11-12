/**
 * Porting of the libm library to the PandA framework
 * starting from the original FDLIBM 5.3 (Freely Distributable LIBM) developed by SUN
 * plus the newlib version 1.19 from RedHat and plus uClibc version 0.9.32.1 developed by Erik Andersen.
 * The author of this port is Fabrizio Ferrandi from Politecnico di Milano.
 * The porting fall under the LGPL v2.1, see the files COPYING.LIB and COPYING.LIBM_PANDA in this directory.
 * Date: September, 11 2013.
 */
/* sf_tanh.c -- float version of s_tanh.c.
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */

/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

#include "math_privatef.h"

static const float one = 1.0, two = 2.0, tiny = 1.0e-30;

float tanhf(float x)
{
   float t, z;
   int jx, ix;

   GET_FLOAT_WORD(jx, x);
   ix = jx & 0x7fffffff;

   /* x is INF or NaN */
   if(!FLT_UWORD_IS_FINITE(ix))
   {
      /* tanh(+-inf)=+-1 */
      /* tanh(NaN) = NaN */
      return FLT_UWORD_IS_INFINITE(ix) ? copysignf(one, x) : __builtin_nanf("");
   }

   /* |x| < 22 */
   if(ix < 0x41b00000)
   {                           /* |x|<22 */
      if(ix < 0x24000000)      /* |x|<2**-55 */
         return x * (one + x); /* tanh(small) = small */
      if(ix >= 0x3f800000)
      { /* |x|>=1  */
         t = expm1f(two * fabsf(x));
         z = one - two / (t + two);
      }
      else
      {
         t = expm1f(-two * fabsf(x));
         z = -t / (t + two);
      }
      /* |x| > 22, return +-1 */
   }
   else
   {
      z = one - tiny; /* raised inexact flag */
   }
   return (jx >= 0) ? z : -z;
}

#ifdef _DOUBLE_IS_32BITS

double tanh(double x)
{
   return (double)tanhf((float)x);
}

#endif /* defined(_DOUBLE_IS_32BITS) */
