#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(SPH2GRID, 0.1, arm2arm@gmail.com)
AM_INIT_AUTOMAKE(SPH2GRID, 0.1)
AC_CONFIG_SRCDIR([src/sphvolume.h])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt strrchr strstr])

AC_OUTPUT(Makefile src/Makefile doc/Makefile idl/Makefile)
