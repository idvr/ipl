AC_PREREQ([2.50])
AC_INIT([libipl], [1.01], [andrea.lagala@gmail.com])
AC_CONFIG_AUX_DIR([config])
AM_CONFIG_HEADER([config_auto.h:config/config.h.in])
AC_CONFIG_SRCDIR([src/corner.c])

AC_CONFIG_MACRO_DIR([m4])
LT_INIT

AM_INIT_AUTOMAKE
AC_LANG(C)

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for arguments.
AC_ARG_WITH([libfftw3], AS_HELP_STRING([--without-libfftw3], [do not include libfftw3 support]))
AC_ARG_WITH([libfuzzylite], AS_HELP_STRING([--without-libfuzzylite], [do not include libfuzzylite support]))

# Checks for headers.

# Checks for libraries.
AC_CHECK_LIB([lept], [pixCreate],
    AC_SUBST([LIBLEPT_LIBS], [-llept]), [
	echo "Error! You need to have liblept installed."
	exit -1
	])

AS_IF([test "x$with_libfftw3" != xno],
  AC_CHECK_LIB([fftw3], [fftw_malloc],
    AC_DEFINE([HAVE_LIBFFTW3], 1, [Define to 1 if you have libfftw3.]) AC_SUBST([LIBFTTW3_LIBS], [-lfftw3]),
    AS_IF([test "x$with_libfftw3" = xyes], AC_MSG_ERROR([libfftw3 support requested but library not found]))
  )
)

AS_IF([test "x$with_libfuzzylite" != xno],
  AC_CHECK_LIB([fuzzylite], [deflate],
    AC_DEFINE([HAVE_LIBFUZZYLITE], 1, [Define to 1 if you have libfuzzylite.]) AC_SUBST([LIBFUZZYLITE_LIBS], [-lfuzzylite]),
    AS_IF([test "x$with_libfuzzylite" = xyes], AC_MSG_ERROR([libfuzzylite support requested but library not found]))
  )
)

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
