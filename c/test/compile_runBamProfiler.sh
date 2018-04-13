gcc -std=gnu99 -g -fPIC -pthread -O2 -Wall -I/home/tlamberton/git/BamM/c/htslib-1.3.1/htslib -I/home/tlamberton/git/BamM/c/libcfu-0.03/include -static-libgcc -Wl,-rpath,/home/tlamberton/git/BamM/c/htslib-1.3.1 -o runBamProfiler runBamProfiler.c ./bamFilter.c ./bamProfiler.c -lm -L/home/tlamberton/git/BamM/c/libcfu-0.03/lib -lcfu -L/home/tlamberton/git/BamM/c/htslib-1.3.1 -lhts

