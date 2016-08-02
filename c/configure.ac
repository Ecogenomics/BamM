dnl Autoconf for BamM
AC_PREREQ(2.61)
AC_INIT(BamM, 0.1, mike@mikeimelfort.com)
AC_CONFIG_MACRO_DIR([m4])

dnl Keep helpers in build-aux/ subdirectory to not leave too much junk.
#AC_CONFIG_AUX_DIR([build-aux])

# test for a c compiler
AC_PROG_CC
AC_PROG_CXX

# test for libm
m_flags=
AC_CHECK_LIB([m],[lgamma,log,exp],
             [m_flags="-lm"
              AC_DEFINE([HAVE_M],[1],[Defines to 1 if the m library (-lm) is found])],
             [AC_MSG_ERROR([libm not found])])

# test for libcfu
AX_LIBCFU
if test $HAVE_LIBCFU = no; then
AC_MSG_ERROR([Cannot find libcfu])
else
AM_CONDITIONAL([LIBCFU_FOUND], [test $HAVE_LIBCFU = yes])
fi

# test for libhts
AX_LIBHTS
if test $HAVE_LIBHTS = no; then
AC_MSG_ERROR([Cannot find libhts])
else
AM_CONDITIONAL([LIBHTS_FOUND], [test $HAVE_LIBHTS = yes])
fi

dnl Add the default location for libcfu and htslib
AC_CONFIG_SUBDIRS([libcfu-0.03 htslib-1.3.1])

dnl Write #defines to config.h file.
dnl You need to generate config.h.in, for example using autoheader.
AC_CONFIG_HEADERS([config.h])
dnl Generate Makefile from Makefile.in.
AC_CONFIG_FILES([Makefile])
dnl This generates fore-mentioned files.
AC_OUTPUT
