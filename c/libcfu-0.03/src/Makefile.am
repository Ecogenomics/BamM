__FINAL_LIB_dir = $(top_srcdir)/lib
__FINAL_LIB__LIBRARIES = libcfu.a
libcfu_a_SOURCES = cfuhash.c cfutime.c cfustring.c cfulist.c cfuconf.c cfuthread_queue.c cfu.c cfuopt.c
AM_CFLAGS = -I../include -std=gnu99 -W -pthread -g -fPIC -pthread -O2 -Wall -Wno-unused-but-set-variable -Wno-unused-function
LIB_FLAGS = -static-libgcc -shared -Wl,-rpath

AM_CPPFLAGS = -I$(top_srcdir)/include

default: libcfu.a

clean:
	rm -f *.o
	rm -f ../lib/libcfu.a

