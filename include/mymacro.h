/* ----------------- */
/* --- mymacro.h --- */
/* ----------------- */

#ifndef __MY_MACRO_H__
#define __MY_MACRO_H__

// activer ou desactiver le define ci-dessous pour passer du mode de mise au point
// au mode benchmark

#define ZERO zero_f32matrix(Z, 0, n-1, 0, n-1)

#define ENABLE_BENCHMARK

// -------------------------------------------
// -- ne rien ecrire en dessous de cette ligne
// -------------------------------------------

#ifdef ENABLE_BENCHMARK

#pragma message("ENABLE_BENCHMARK is ON")

// check predef
// https://sourceforge.net/p/predef/wiki/Architectures/

#ifdef __INTEL_COMPILER
#pragma message("icc compiler")
#include <ia32intrin.h>
#endif

#ifdef __GNUC__
#pragma message("gcc compiler ...")
#include <x86intrin.h>
#endif

//#define CHRONO(X,t)  tmin = 1e38; for(run=0; run<nrun; run++) { t0 = (double) _rdtsc(); for(iter=0; iter<niter; iter++) { X; } t1 = (double) _rdtsc(); dt=t1-t0; if(dt<tmin) tmin = dt; } t = tmin / (double) niter

#define CHRONO(X,t)                       \
    tmin = 1e38;                          \
    for(run=0; run<nrun; run++) {         \
        t0 = (double) _rdtsc();           \
        for(iter=0; iter<niter; iter++) { \
            X;                            \
        }                                 \
        t1 = (double) _rdtsc();           \
        dt=t1-t0; if(dt<tmin) tmin = dt;  \
    }                                     \
    t = tmin / (double) niter

#define BENCH(X) X
#define DEBUG(X)

#else

#pragma message("ENABLE_BENCHMARK is OFF")

#define CHRONO(X,t)  X
#define BENCH(X)
#define DEBUG(X) X
#endif

#endif // __MY_MACRO_H__
