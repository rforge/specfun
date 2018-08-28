/*
 *  Copyright (C) 2004--2018 Martin Maechler

 *  Providing an R interface to the new (in March 2004)
 *  non-API  qchisq_appr()  in R's src/nmath/qgamma.c
 *
 * Goal: Investigate use of becoming it part of RMathlib -- for other packages to use
 */

#include <math.h>
#include <float.h> /* DBL_MIN etc */

#include <Rmath.h>

#include "DPQpkg.h"


/* For now copy-pasted here, version of qgamma.c unchanged since Aug.2015 [as of 2018-08-23]
 * Source at https://svn.r-project.org/R/trunk/src/nmath/qgamma.c
 *           ====================================================
 * MM: at ~/R/D/r-devel/R/src/nmath/qgamma.c
 *        ----------------------------------
 *
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000--2015 The R Core Team
 *  Copyright (C) 2004--2015 The R Foundation
 *  based on AS 91 (C) 1979 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *	Compute the quantile function of the gamma distribution.
 *
 *  NOTES
 *
 *	This function is based on the Applied Statistics
 *	Algorithm AS 91 ("ppchi2") and via pgamma(.) AS 239.
 *
 *	R core improvements:
 *	o  lower_tail, log_p
 *      o  non-trivial result for p outside [0.000002, 0.999998]
 *	o  p ~ 1 no longer gives +Inf; final Newton step(s)
 *
 *  REFERENCES
 *
 *	Best, D. J. and D. E. Roberts (1975).
 *	Percentage Points of the Chi-Squared Distribution.
 *	Applied Statistics 24, page 385.  */

#ifdef USING_R
#include "nmath.h"
#include "dpq.h"
#else // MM_R our substitute for the above R/src/nmath-internal headers :

// -------- excerpt from  nmath.h ---------------------------------

#include <R_ext/Error.h>
# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) warning(fmt,x,x2,x3,x4,x5)

#include <R_ext/Arith.h>
#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) gettext (String)
#else
#define _(String) (String)
#endif


#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

/* FIXME: 'ERR' and 'ERROR' below are misnomers, both in R and stand-alone Rmathlib:
   -----  they only *warn* after all; they are private headers, so can change anytime
*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }

/* For a long time prior to R 2.3.0 ML_ERROR did nothing.
   We don't report ME_DOMAIN errors as the callers collect ML_NANs into
   a single warning.
 */
#define ML_ERROR(x, s) { \
   if(x > ME_DOMAIN) { \
       char *msg = ""; \
       switch(x) { \
       case ME_DOMAIN: \
	   msg = _("argument out of domain in '%s'\n");	\
	   break; \
       case ME_RANGE: \
	   msg = _("value out of range in '%s'\n");	\
	   break; \
       case ME_NOCONV: \
	   msg = _("convergence failed in '%s'\n");	\
	   break; \
       case ME_PRECISION: \
	   msg = _("full precision may not have been achieved in '%s'\n"); \
	   break; \
       case ME_UNDERFLOW: \
	   msg = _("underflow occurred in '%s'\n");	\
	   break; \
       } \
       MATHLIB_WARNING(msg, s); \
   } \
}
#endif
// MM_R  end_of { #include <nmath.h> substitute }


#ifdef DEBUG_qgamma
# define DEBUG_q
# include <R_ext/Print.h>
#endif

// MM_R attribute_hidden
double qchisq_appr(double p, double nu, double g /* = log Gamma(nu/2) */,
		   logical lower_tail, logical log_p, double tol /* EPS1 */)
{
#define C7	4.67
#define C8	6.66
#define C9	6.73
#define C10	13.32

    double alpha, a, c, ch, p1;
    double p2, q, t, x;

    /* test arguments and initialise */

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(nu))
	return p + nu;
#endif
    R_Q_P01_check(p);
    if (nu <= 0) ML_ERR_return_NAN;

    alpha = 0.5 * nu;/* = [pq]gamma() shape */
    c = alpha-1;

    if(nu < (-1.24)*(p1 = R_DT_log(p))) {	/* for small chi-squared */
	/* log(alpha) + g = log(alpha) + log(gamma(alpha)) =
	 *        = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
	 *  catastrophic cancellation when alpha << 1
	 */
	double lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g);
	ch = exp((lgam1pa + p1)/alpha + M_LN2);
#ifdef DEBUG_qgamma
	REprintf(" small chi-sq., ch0 = %g\n", ch);
#endif

    } else if(nu > 0.32) {	/*  using Wilson and Hilferty estimate */

	x = qnorm(p, 0, 1, lower_tail, log_p);
	p1 = 2./(9*nu);
	ch = nu*pow(x*sqrt(p1) + 1-p1, 3);

#ifdef DEBUG_qgamma
	REprintf(" nu > .32: Wilson-Hilferty; x = %7g\n", x);
#endif
	/* approximation for p tending to 1: */
	if( ch > 2.2*nu + 6 )
	    ch = -2*(R_DT_Clog(p) - c*log(0.5*ch) + g);

    } else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */

	ch = 0.4;
	a = R_DT_Clog(p) + g + c*M_LN2;
#ifdef DEBUG_qgamma
	REprintf(" nu <= .32: a = %7g\n", a);
#endif
	do {
	    q = ch;
	    p1 = 1. / (1+ch*(C7+ch));
	    p2 = ch*(C9+ch*(C8+ch));
	    t = -0.5 +(C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2;
	    ch -= (1- exp(a+0.5*ch)*p2*p1)/t;
	} while(fabs(q - ch) > tol * fabs(ch));
    }

    return ch;
}
// end of {copy-paste from R/src/nmath/qgamma.c}

//--------- a simple vectorized .C() interface --------------------

// FIXME: Use a .Call() and then vectorize in both main args (p, nu)
void qchisq_appr_v(double *P, int *n, double *nu, double *tol,
		   logical *lower_tail, logical *log_p,
		   /* result: */ double *q)
{
    double g = lgammafn(0.5* *nu);
    for(int i = 0; i < *n; i++)
	q[i] = qchisq_appr(P[i], *nu, g, *lower_tail, *log_p, *tol);
    return;
}
