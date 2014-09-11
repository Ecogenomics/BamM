dnl SYNOPSIS
dnl
dnl AX_LIBHTS
dnl
dnl DESCRIPTION
dnl
dnl This macro provides tests of availability for libhts. This macros
dnl checks for libhts Parser headers and libraries and defines compilation flags
dnl
dnl Macro supports following options and their values:
dnl
dnl --with-libhts-inc - path to base directory with libhts headers
dnl --with-libhts-lib - linker flags for libhts
dnl
dnl This macro calls:
dnl
dnl AC_SUBST(LIBHTS_CPPFLAGS)
dnl AC_SUBST(LIBHTS_LDFLAGS)
dnl AC_SUBST([LIBHTS_LIBS])
dnl AC_SUBST(LIBHTS_LIB_DIR)
dnl
dnl And sets:
dnl
dnl HAVE_LIBHTS
dnl
dnl LICENSE
dnl
dnl Copyright (c) 2012,2013 Connor Skennerton, Michael Imelfort
dnl
dnl Copying and distribution of this file, with or without modification, are
dnl permitted in any medium without royalty provided the copyright notice
dnl and this notice are preserved. This file is offered as-is, without any
dnl warranty.

AC_DEFUN([AX_LIBHTS],
[
    AC_ARG_WITH([libhts-inc],
        AS_HELP_STRING([--with-libhts-inc=@<:@ARG@:>@],
        [libhts headers are at this location (ARG=path)]
        ),
        [libhts_include_dir="$withval"],
        [
        dnl Default behavior is implicit yes
        if test -d /usr/local/include/htslib ; then
            libhts_include_dir=/usr/local/include
        elif test -d /usr/include/htslib; then
            libhts_include_dir=/usr/include
        else
            libhts_include_dir=""
        fi
        ]
    )
    
    AC_ARG_WITH([libhts-lib],
        AS_HELP_STRING([--with-libhts-lib=@<:@ARG@:>@],
        [libhts library is at this location (ARG=path)]
        ),
        [libhts_lib_dir="$withval"],
        [
        dnl Default behavior is implicit yes
        if test -f /usr/local/lib/libhts.so ; then
            libhts_lib_dir=/usr/local/lib
        elif test -f /usr/lib/libhts.so; then
            libhts_lib_dir=/usr/lib
        else
            libhts_lib_dir=""
        fi
        ]
    )

    LIBHTS_CPPFLAGS=""
    LIBHTS_LDFLAGS=""

    dnl
    dnl Collect include/lib paths and flags
    dnl

    if test -n "$libhts_lib_dir"; then
        libhts_ldflags="-L$libhts_lib_dir"
    else
        libhts_lib_found="no"
        AC_MSG_RESULT([not found])
    fi
    
    libhts_libs="-lhts"

    dnl
    dnl Check libhts files
    dnl

    saved_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS -I$libhts_include_dir"

    saved_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $libhts_ldflags $PTHREAD_LDFLAGS"

    saved_LIBS="$LIBS"
    LIBS="$libhts_libs $PTHREAD_LIBS $LIBS"
    
    dnl echo "libs: ${LIBS}"
    dnl echo "cppflags: ${CPPFLAGS}"
    dnl echo "ldflags: ${LDFLAGS}"
    dnl
    dnl Check libhts headers
    dnl
    AC_MSG_CHECKING([for libhts headers in $libhts_include_dir])

    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM(
            [[
@%:@include <$libhts_include_dir/htslib/sam.h>
            ]],
            [[]]
        )],
        [
        LIBHTS_CPPFLAGS="-I$libhts_include_dir"
        libhts_header_found="yes"
        AC_MSG_RESULT([found])
        ],
        [
        libhts_header_found="no"
        AC_MSG_RESULT([not found])
        ]
    )

    dnl
    dnl Check libhts libraries
    dnl
    if test "$libhts_header_found" = "yes"; then

        AC_MSG_CHECKING([for libhts libraries in $libhts_lib_dir])

        AC_LINK_IFELSE([
            AC_LANG_PROGRAM(
                [[
@%:@include <$libhts_include_dir/htslib/sam.h>
                ]],
                [[
bam_hdr_t *h = 0;
                ]]
            )],
            [
            LIBHTS_LDFLAGS="$libhts_ldflags"
            LIBHTS_LIBS="$libhts_libs"
            libhts_lib_found="yes"
            AC_MSG_RESULT([found])
            ],
            [
            libhts_lib_found="no"
            AC_MSG_RESULT([not found])
            ]
        )
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"
 
    if test "$libhts_header_found" = "yes" -a "$libhts_lib_found" = "yes"; then

        LIBHTS_LIB_DIR="$libhts_lib_dir"
    
        AC_SUBST([LIBHTS_CPPFLAGS])
        AC_SUBST([LIBHTS_LDFLAGS])
        AC_SUBST([LIBHTS_LIBS])
        AC_SUBST([LIBHTS_LIB_DIR])

        HAVE_LIBHTS="yes"
    else
        HAVE_LIBHTS="no"
    fi
])