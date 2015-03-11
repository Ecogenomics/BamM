# generated automatically by aclocal 1.14.1 -*- Autoconf -*-

# Copyright (C) 1996-2013 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# AM_CONDITIONAL                                            -*- Autoconf -*-

# Copyright (C) 1997-2013 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[AC_PREREQ([2.52])dnl
 m4_if([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
       [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])dnl
AC_SUBST([$1_FALSE])dnl
_AM_SUBST_NOTMAKE([$1_TRUE])dnl
_AM_SUBST_NOTMAKE([$1_FALSE])dnl
m4_define([_AM_COND_VALUE_$1], [$2])dnl
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([[conditional "$1" was never defined.
Usually this means the macro was only invoked conditionally.]])
fi])])

# Copyright (C) 2006-2013 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_SUBST_NOTMAKE(VARIABLE)
# ---------------------------
# Prevent Automake from outputting VARIABLE = @VARIABLE@ in Makefile.in.
# This macro is traced by Automake.
AC_DEFUN([_AM_SUBST_NOTMAKE])

# AM_SUBST_NOTMAKE(VARIABLE)
# --------------------------
# Public sister of _AM_SUBST_NOTMAKE.
AC_DEFUN([AM_SUBST_NOTMAKE], [_AM_SUBST_NOTMAKE($@)])

dnl SYNOPSIS
dnl
dnl AX_LIBCFU
dnl
dnl DESCRIPTION
dnl
dnl This macro provides tests of availability for libcfu. This macros
dnl checks for libcfu Parser headers and libraries and defines compilation flags
dnl
dnl Macro supports following options and their values:
dnl
dnl --with-libcfu-inc - path to base directory with libcfu headers
dnl --with-libcfu-lib - linker flags for libcfu
dnl
dnl This macro calls:
dnl
dnl AC_SUBST(LIBCFU_CPPFLAGS)
dnl AC_SUBST(LIBCFU_LDFLAGS)
dnl AC_SUBST([LIBCFU_LIBS])
dnl AC_SUBST(LIBCFU_LIB_DIR)
dnl
dnl And sets:
dnl
dnl HAVE_LIBCFU
dnl
dnl LICENSE
dnl
dnl Copyright (c) 2012,2013 Connor Skennerton, Michael Imelfort
dnl
dnl Copying and distribution of this file, with or without modification, are
dnl permitted in any medium without royalty provided the copyright notice
dnl and this notice are preserved. This file is offered as-is, without any
dnl warranty.

AC_DEFUN([AX_LIBCFU],
[
    AC_ARG_WITH([libcfu-inc],
        AS_HELP_STRING([--with-libcfu-inc=@<:@ARG@:>@],
        [libcfu headers are at this location (ARG=path)]
        ),
        [libcfu_include_dir="$withval"],
        [
        dnl Default behavior is implicit yes
if test -f /usr/local/include/cfuhash.h ; then
            libcfu_include_dir=/usr/local/include
        elif test -f /usr/include/cfuhash.h; then
            libcfu_include_dir=/usr/include
        else
            libcfu_include_dir=""
        fi
        ]
    )

    AC_ARG_WITH([libcfu-lib],
        AS_HELP_STRING([--with-libcfu-lib=@<:@ARG@:>@],
        [libcfu library is at this location (ARG=path)]
        ),
        [libcfu_lib_dir="$withval"],
        [
        dnl Default behavior is implicit yes
if test -f /usr/local/lib/libcfu.a ; then
            libcfu_lib_dir=/usr/local/lib
        elif test -f /usr/lib/libcfu.a; then
            libcfu_lib_dir=/usr/lib
        else
            libcfu_lib_dir=""
        fi
        ]
    )

    LIBCFU_CPPFLAGS=""
    LIBCFU_LDFLAGS=""

    dnl
    dnl Collect include/lib paths and flags
    dnl

    if test -n "$libcfu_lib_dir"; then
        libcfu_ldflags="-L$libcfu_lib_dir"
    else
        saved_dir=$(pwd)
        CFU_SRC_DIR=$saved_dir/"libcfu-0.03"
        AC_MSG_RESULT([libcfu location not specified - building from local version])
        cd $CFU_SRC_DIR
        ./configure
        make clean && make && make install
        cd $saved_dir
        libcfu_ldflags="-L$CFU_SRC_DIR/lib"
        libcfu_lib_dir="$CFU_SRC_DIR/lib"
        libcfu_include_dir="$CFU_SRC_DIR/include"
    fi

    libcfu_libs="-lcfu"

    dnl
    dnl Check libcfu files
    dnl

    saved_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS -I$libcfu_include_dir"

    saved_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $libcfu_ldflags $PTHREAD_LDFLAGS"

    saved_LIBS="$LIBS"
    LIBS="$libcfu_libs $PTHREAD_LIBS $LIBS"

    dnl echo "libs: ${LIBS}"
    dnl echo "cppflags: ${CPPFLAGS}"
    dnl echo "ldflags: ${LDFLAGS}"
    dnl
    dnl Check libcfu headers
    dnl
    AC_MSG_CHECKING([for libcfu headers in $libcfu_include_dir])

    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM(
            [[
@%:@include <$libcfu_include_dir/cfuhash.h>
            ]],
            [[]]
        )],
        [
        LIBCFU_CPPFLAGS="-I$libcfu_include_dir"
        libcfu_header_found="yes"
        AC_MSG_RESULT([found])
        ],
        [
        libcfu_header_found="no"
        AC_MSG_RESULT([not found])
        ]
    )

    dnl
    dnl Check libcfu libraries
    dnl
    if test "$libcfu_header_found" = "yes"; then

        AC_MSG_CHECKING([for libcfu libraries in $libcfu_lib_dir])

        AC_LINK_IFELSE([
            AC_LANG_PROGRAM(
                [[
@%:@include <$libcfu_include_dir/cfuhash.h>
                ]],
                [[
cfuhash_table_t *links = cfuhash_new_with_initial_size(30); cfuhash_destroy(links);
                ]]
            )],
            [
            LIBCFU_LDFLAGS="$libcfu_ldflags"
            LIBCFU_LIBS="$libcfu_libs"
            libcfu_lib_found="yes"
            AC_MSG_RESULT([found])
            ],
            [
            libcfu_lib_found="no"
            AC_MSG_RESULT([not found])
            ]
        )
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

    if test "$libcfu_header_found" = "yes" -a "$libcfu_lib_found" = "yes"; then

        AC_SUBST([LIBCFU_CPPFLAGS])
        AC_SUBST([LIBCFU_LDFLAGS])
        AC_SUBST([LIBCFU_LIBS])
        AC_SUBST([LIBCFU_LIB_DIR])

        HAVE_LIBCFU="yes"
    else
        HAVE_LIBCFU="no"
    fi
])
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
        elif test -d htslib; then
            libhts_include_dir=htslib
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
        saved_dir=$(pwd)
        HTS_SRC_DIR=$saved_dir/"htslib-1.2.1"
        AC_MSG_RESULT([htslib location not specified - building from local version])
        cd $HTS_SRC_DIR
        ./configure
        make clean && make
        cd $saved_dir
        libhts_ldflags="-L$HTS_SRC_DIR"
        libhts_lib_dir="$HTS_SRC_DIR"
        libhts_include_dir="$HTS_SRC_DIR/htslib"
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
@%:@include <$libhts_include_dir/sam.h>
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
@%:@include <$libhts_include_dir/sam.h>
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
