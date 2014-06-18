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
LIB_FLAGS = -static-libgcc -shared -Wl,-rpath,$(LIBHTS_LIB_DIR),-soname,libBamM.so.0
LIBS =  -lm $(LIBCFU_LDFLAGS) $(LIBCFU_LIBS) $(LIBHTS_LDFLAGS) $(LIBHTS_LIBS)

BMM_PARSER = bamParser
BMM_EXTRACTOR = bamExtractor
BMM_LIB = libBamM.a

PARSER_SOURCES = parser.c bamParser.c pairedLink.c
EXTRACTOR_SOURCES = extractor.c bamParser.c pairedLink.c
LIBBAMM_SOURCES = bamParser.c pairedLink.c

LIBBAMM_OBJS = bamParser.o pairedLink.o

#all: parser extractor library
all: parser library
       
.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(BMM_LIB) : $(LIBBAMM_OBJS)
	@-rm -f $@
	$(CC) $(CFLAGS) $(LIB_FLAGS) -o $@ $(LIBBAMM_SOURCES) $(LIBS)

$(BMM_PARSER): $(PARSER_SOURCES)
	$(CC) $(CFLAGS) -o $@ $(PARSER_SOURCES) $(LIBS)

$(BMM_EXTRACTOR): $(EXTRACTOR_SOURCES)
	$(CC) $(CFLAGS) -o $@ $(EXTRACTOR_SOURCES) $(LIBS)

extractor: $(BMM_EXTRACTOR)
	@echo "Building extractor"

parser: $(BMM_PARSER)
	@echo "Building parser"

library: $(BMM_LIB)
	@echo "Building library"

clean:
	$(RM) $(BMM_PARSER)
	$(RM) $(BMM_EXTRACTOR)
	$(RM) *.o
	$(RM) $(BMM_LIB)
