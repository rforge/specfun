#include "DPQpkg.h" // before stdio (for MINGW_...)

#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _typ)/sizeof(name ## _typ[0]), name ##_typ}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


/* .C calls */

// qchisq_appr.c :
static R_NativePrimitiveArgType qchisq_appr_v_typ[7] = {
    REALSXP, INTSXP, REALSXP, REALSXP,
    LGLSXP,  LGLSXP, REALSXP };
// pnchisq-it.c :
static R_NativePrimitiveArgType Pnchisq_it_typ[] = {
    /* x, f, theta : */ REALSXP, REALSXP, REALSXP,
    /* errmax, reltol, itrmax, verbose: */ REALSXP, REALSXP, INTSXP, INTSXP,
    /* i_0, n_terms: */ INTSXP,  INTSXP,
    /* terms, prob : */ REALSXP, REALSXP };
static R_NativePrimitiveArgType ncbeta_typ[] = {
    REALSXP, REALSXP, REALSXP, REALSXP,
    /* n: */ INTSXP, LGLSXP,
    /* errmax: */ REALSXP, INTSXP, INTSXP,
    /* res: */ REALSXP };


// wienergerm_nchisq.c :
static R_NativePrimitiveArgType pchisqV_typ[] = {
    REALSXP, INTSXP,  REALSXP, REALSXP,
    LGLSXP,  LGLSXP,  INTSXP };



static const R_CMethodDef CEntries[] = {
    CDEF(qchisq_appr_v),
    CDEF(Pnchisq_it),
    CDEF(ncbeta),
    CDEF(pchisqV),

    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(R_algdiv,   2), // <-- ./algdiv.c
    CALLDEF(chk_LDouble,3), // <-- ./ppois-direct.c
    CALLDEF(ppoisD,     4), // <-- ./ppois-direct.c
    CALLDEF(Pnchisq_R, 14), // <-- ./pnchisq-it.c

    CALLDEF(R_log1pmx,  1), // <--> DPQ-misc.c
    CALLDEF(R_log1pexp, 1),
    CALLDEF(R_log1mexp, 1),
    CALLDEF(R_lgamma1p, 1),
    CALLDEF(R_frexp, 1),
    CALLDEF(R_ldexp, 2),

    {NULL, NULL, 0}
};

static R_FortranMethodDef FortEntries[] = {
    {"noncechi", (DL_FUNC) &F77_SUB(noncechi), 6}, // <- ./wienergerm_nchisq_F.f
    {NULL, NULL, 0}
};


void R_init_DPQ(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
