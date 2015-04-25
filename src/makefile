#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

PROGS = PKA
# sequence_feature

# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++11 -O2

# define any directories containing header files other than /usr/include
#
INCLUDES = -I../include


.PHONY: depend clean

all: $(PROGS) install

PKA: sequence.o markov.o positional_kmer.o container.o utility.o text.o stat.o ushuffle.o PKA.o
	$(CC) $(CFLAGS) $^ -o $@

#sequence_feature: sequence_feature.o utility.o sequence.o ushuffle.o  positional_kmer.o stat.o container.o
#	$(CC) $(CFLAGS) $^ -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

# ushuffle.c
.c.o:
	gcc -O3 $(INCLUDES) -c $< -o $@	

install:
	cp $(PROGS) ../bin
	mv $(PROGS) ../web/$(PROGS)

clean:
	$(RM) *.o *~ ../bin/*

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
