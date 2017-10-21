#
# Makefile for mxlu1fac
#
# Copyright (C) 2005 Yin Zhang yzhang@cs.utexas.edu
#

CC = gcc
MX = mex
RM = rm -f

OBJS 	= commonlib.o hbio.o lusol.o lusolio.o lusolmain.o mmio.o myblas.o

MXSRCS  = mxlu1fac.c commonlib.c lusol.c myblas.c

CFLAGS 	= -DYZHANG -DNDEBUG -Wall -fexceptions -pthread -ansi -fPIC -O3 # -g

MXFLAGS = -DMATLAB -DYZHANG -O # -g

#BLAS 	= /p/lib/prescott_gotoblas.so
BLAS	= /usr/local/matlab/bin/glnx86/libmwrefblas.so

LIBS 	= $(BLAS) -lm -lpthread -ldl

#--------------------------------------------------------------------------

all: mxlu1fac lusol

mxlu1fac: $(MXSRCS)
	$(MX) $(MXFLAGS) $(MXSRCS) $(LIBS)

lusol: $(OBJS)
	$(CC) -o lusol $(CFLAGS) $(OBJS) $(LIBS)

clean: 
	$(RM) -f $(OBJS) *.dll 

distclean: clean
	$(RM) -f *.mex* lusol



