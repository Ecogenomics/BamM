AC_CONFIG_MACRO_DIR([m4])    
AC_INIT([libcfu], [0.03], mike@mikeimelfort.com, [libcfu])
AM_INIT_AUTOMAKE

AC_PROG_LIBTOOL
AC_LANG_CPLUSPLUS
AC_PROG_CC
AC_PROG_CXX
SAVECFLAGS="$CFLAGS"
SAVELIBS="$LIBS"
CFLAGS="-pthread $CFLAGS"
LIBS="-pthread $LIBS"
AC_CHECK_FUNC(pthread_mutex_init,,
    CFLAGS="$SAVECFLAGS"
    LIBS="$SAVELIBS"
    AC_CHECK_LIB(pthread, pthread_mutex_init,,
        AC_CHECK_LIB(pthreads, pthread_mutex_init,,
            AC_MSG_ERROR([No pthreads support found.]))))

AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])

AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_CHECK_FUNCS([gettimeofday memset strcasecmp strncasecmp])
AC_OUTPUT(Makefile src/Makefile)

