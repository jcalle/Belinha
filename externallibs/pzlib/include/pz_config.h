// This file is a template that will be used during INSTALL process
// It holds all definitions (-D**, libraries, etc) passed to compiler when building this library


#ifndef PZSOURCEDIR
#define PZSOURCEDIR "/Users/labmec/Documents/GitHub/neopz"
#endif

#ifndef REFPATTERNDIR
#define REFPATTERNDIR "/Users/labmec/Documents/GitHub/neopz/Refine/RefPatterns"
#endif

#ifndef PZ_BRANCH
#define PZ_BRANCH "master"
#endif

#ifndef PZ_REVISION
#define PZ_REVISION "afd0b2314"
#endif

#ifndef PZ_REVISION_DATE
#define PZ_REVISION_DATE "Tue Jan 22 17:39:59 2019"
#endif

#ifndef LOG4CXX
/* #undef LOG4CXX */
#endif

#ifndef _AUTODIFF
/* #undef _AUTODIFF */
#endif

#ifndef USING_BOOST
/* #undef USING_BOOST */
#endif

#ifndef USING_IPO
/* #undef USING_IPO */
#endif

#ifndef USING_METIS
/* #undef USING_METIS */
#endif

#ifndef USING_OPENSSL
/* #undef USING_OPENSSL */
#endif

#ifndef USING_TBB
/* #undef USING_TBB */
#endif

#ifndef USING_OPENMP
/* #undef USING_OPENMP */
#endif

#ifndef USING_LIKWID
/* #undef USING_LIKWID */
#endif

#ifndef USING_LIBNUMA
/* #undef USING_LIBNUMA */
#endif

#ifndef USING_MATLAB_ENGINE
/* #undef USING_MATLAB_ENGINE */
#endif

#ifndef USING_BLAS
/* #undef USING_BLAS */
#endif

#ifndef USING_LAPACK
/* #undef USING_LAPACK */
#endif

#ifndef USING_NEW_SKYLMAT
/* #undef USING_NEW_SKYLMAT */
#endif

#ifndef USING_PAPI
/* #undef USING_PAPI */
#endif

#ifndef USING_HWLOC
/* #undef USING_HWLOC */
#endif

#ifndef USING_MKL
/* #undef USING_MKL */
#endif

#ifdef USING_MKL
#define USING_LAPACK
#define USING_BLAS
#endif

//WIN32 defs
#ifndef VC
/* #undef VC */
#endif

//Type defs
#ifndef REALdouble
#define REALdouble
#endif

#ifndef STATEdouble
#define STATEdouble
#endif

#ifndef HDIVPIOLA
/* #undef HDIVPIOLA */
#endif

#ifndef STATE_COMPLEX
/* #undef STATE_COMPLEX */
#endif
