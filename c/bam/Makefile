CC = gcc

LIBHTS_CPPFLAGS = -I/home/uqmimelf/local/opt/htslib
LIBHTS_LDFLAGS = -L/home/uqmimelf/local/opt/htslib
LIBHTS_LIBS = -lhts
LIBHTS_LIB_DIR = /home/uqmimelf/local/opt/htslib

LIBCFU_CPPFLAGS = -I/home/uqmimelf/local/opt/libcfu-0.03/include
LIBCFU_LDFLAGS = -L/home/uqmimelf/local/opt/libcfu-0.03/lib
LIBCFU_LIBS = -lcfu
LIBCFU_LIB_DIR = 

CFLAGS = -g -fPIC -pthread -O2 -Wall $(LIBHTS_CPPFLAGS) $(LIBCFU_CPPFLAGS)
LIB_FLAGS = -static-libgcc -shared -Wl,-rpath,$(LIBHTS_LIB_DIR),-soname,libPMBam.so.0
LIBS =  -lm $(LIBCFU_LDFLAGS) $(LIBCFU_LIBS) $(LIBHTS_LDFLAGS) $(LIBHTS_LIBS)
EXECUTABLE = bamParser
BMM_BAM_LIB = libPMBam.a

TEST_SOURCES = example.c bamParser.c pairedLink.c
LIB_SOURCES = bamParser.c pairedLink.c

LIBPMBAM_OBJS = \
        bamParser.o \
        pairedLink.o

all: test library
        
.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(BMM_BAM_LIB) : $(LIBPMBAM_OBJS)
	@-rm -f $@
	$(CC) $(CFLAGS) $(LIB_FLAGS) -o $@ $(LIB_SOURCES) $(LIBS)

$(EXECUTABLE): $(TEST_SOURCES)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

test: $(EXECUTABLE)

library: $(BMM_BAM_LIB)

clean:
	$(RM) $(EXECUTABLE)
	$(RM) *.o
	$(RM) $(BMM_BAM_LIB)
