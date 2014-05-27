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
        libcfu_lib_found="no"
        AC_MSG_RESULT([not found])
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