AC_PREREQ([2.59])
AC_INIT([cpn_nc], [0.4.7], [claudio.bonati82@gmail.com])
AM_INIT_AUTOMAKE(subdir-objects)
AC_CONFIG_SRCDIR([./include/macro.h])
AC_CONFIG_HEADERS([config.h])


# overwrite default values
if test x${ac_cv_env_CFLAGS_set} != xset; then
  CFLAGS=''
fi 
if test x${ac_cv_env_CPPFLAGS_set} != xset; then
  CPPFLAGS=''
fi 



### GENERAL CHECKS


# Checks for programs.
AC_PROG_AWK
AC_PROG_CC


# Checks for header files.
AC_CHECK_HEADERS([limits.h], [], [AC_MSG_ERROR([Required header limits.h not found.])])

AC_CHECK_HEADERS([math.h], [], [AC_MSG_ERROR([Required header math.h not found.])])

AC_CHECK_HEADERS([openssl/md5.h], [], [AC_MSG_ERROR([Required openssl/md5.h not found.])])

AC_CHECK_HEADERS([stdio.h], [], [AC_MSG_ERROR([Required header stdio.h not found.])])

AC_CHECK_HEADERS([stdlib.h], [], [AC_MSG_ERROR([Required header stdlib.h not found.])])

AC_CHECK_HEADERS([time.h], [], [AC_MSG_ERROR([Required header time.h not found.])])


# Checks for libraries.
AC_SEARCH_LIBS([MD5_Init], [crypto], [], [AC_MSG_ERROR([Required crypto library not found.])])

AC_SEARCH_LIBS([sqrt], [m], [], [AC_MSG_ERROR([Required c mathematical library not found.])]) 

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T



### ENABLING & CO


### FURTHER CHECKS



### VARIABLES

AC_ARG_VAR(N_f, the number of flavours (default 2))

AC_ARG_VAR(ST_dim, spacetime dimensionality (default 3))


### DEFINE & CO


### TEST

if test x${N_f} = x; then
  AC_MSG_WARN([ 
   *****
   No chosen value for N_f: using the default N_f=2
   ***** 
   ])
  AC_DEFINE_UNQUOTED(NFLAVOUR, 2, [Number of flavours])
else
  AC_DEFINE_UNQUOTED(NFLAVOUR, ${N_f}, [Number of flavours])
fi


if test x${ST_dim} = x; then
  AC_MSG_WARN([ 
   *****
   No chosen value for ST_dim: using the default ST_dim=3
   ***** 
   ])
  AC_DEFINE_UNQUOTED(STDIM, 3, [Spacetime dimensionality])
else
  AC_DEFINE_UNQUOTED(STDIM, ${ST_dim}, [Spacetime dimensionality])
fi


### CONDITIONALS


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
