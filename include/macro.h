#ifndef MACRO_H
#define MACRO_H

#include"../config.h"

//#define DEBUG  // activate some debug tests

#define CSTAR_BC   //C* bounday conditions

//#define TEMPORAL_GAUGE  // fix theta=0 on temporal links

//#define LORENZ_GAUGE // add (alpha/2)(\sum_{mu} \partial_{mu} theta_{mu})^2 to the energy

//#define LINKS_FIXED_TO_ONE  //to remove gauge fields

//#define COMPACT_MASS // (m^2/2)theta^2 -> -m^2 cos(theta)

#ifdef TEMPORAL_GAUGE
  #ifndef CSTAR_BC
    #error "TEMPORAL_GAUGE can be defined only with C* boundary conditions"
  #endif

  #ifdef LORENZ_GAUGE
    #error "TEMPORAL_GAUGE and LORENZ_GAUGE can not be used togheter"
  #endif
#endif

#ifdef LORENZ_GAUGE
  #ifdef TEMPORAL_GAUGE
    #error "TEMPORAL_GAUGE and LORENZ_GAUGE can not be used togheter"
  #endif
#endif

// function to access matrix elements
#define m(X,Y) ((X)*NFLAVOUR + (Y))

#define MIN_VALUE 1.0e-13

#define INT_ALIGN 16
#define DOUBLE_ALIGN 32

static const double PI=3.141592653589793238462643383279502884197169399375105820974944;
static const double PI2=6.283185307179586476925286766559005768394338798750211641949889;
static const double HALF_PI=1.570796326794896619231321691639751442098584699687552910487472;

#define STD_STRING_LENGTH 50 // standarg lenght of unknown strings

// way to print a macro: if
// #define val1 val2
// then QUOTEME(val1) give the string "val2"
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

// to activate posix_memalign in stdlib.h
#define _POSIX_C_SOURCE 200809L

#endif
