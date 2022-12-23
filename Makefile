# SWIPE
# 
# Smith-Waterman database searches with Inter-sequence Parallel Execution
# 
# Copyright (C) 2008-2013 Torbjorn Rognes, University of Oslo, 
# Oslo University Hospital and Sencel Bioinformatics AS
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contact: Torbjorn Rognes <torognes@ifi.uio.no>, 
# Department of Informatics, University of Oslo, 
# PO Box 1080 Blindern, NO-0316 Oslo, Norway

# Makefile for SWIPE

# BLAST libs: 
# old (2.2.29): /scratch/cube/tuechler/swipe_LIBS
# new (2.13.0): /scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++
# with compo_thresholds: /scratch/cube/tuechler/swipe_LIBS/custom_lib/3

# set LD_LIBRARY_PATH:
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/cube/tuechler/swipe_LIBS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++/ReleaseMT/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/ReleaseMT/lib

# MPI_COMPILE=`mpicxx --showme:compile`
# MPI_LINK=`mpicxx --showme:link`

COMMON=-g

# COMPILE OPTIONS:
# -DTIME_PAIRCOUNT: print runtimes of alignment phases + print number of computed pairs in each alignment phase
# -DCOMPO_THRESHOLDS: allows to set custom thresholds for compositional matrix adjustment (needs to be compiled with modified BLAST library: /scratch/cube/tuechler/swipe_LIBS/custom_lib/3)
COMPILEOPT=-Wall -O3 -DCOMPO_ADJUSTMENT -DTIME_PAIRCOUNT #-DSWLIB_8BIT #-DCOMPO_THRESHOLDS

LIBS=-lpthread -lcomposition_adjustment -lxblast -lxncbi

# GNU options
CXX=g++  # works with newest version gcc/12.2.0

# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I../swipe_LIBS/ncbi-tools++/ -I../swipe_LIBS/src -I../swipe_LIBS/src/algo/blast/core
CXXFLAGS=$(COMPILEOPT) $(COMMON) -I/scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++/include -I/scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++/src -I/scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++/src/algo/blast/core -I/scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++/ReleaseMT/inc
# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/include -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/src -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/src/algo/blast/core -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/ReleaseMT/inc

# LINKFLAGS=$(COMMON) -L../swipe_LIBS/lib
LINKFLAGS=$(COMMON) -L/scratch/cube/tuechler/ncbi_blast2/ncbi-blast-2.13.0+-src/c++/ReleaseMT/lib
# LINKFLAGS=$(COMMON) -L/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/ReleaseMT/lib

PROG=swipe_sswlib # mpiswipe

all : $(PROG)

clean :
	rm -f *.o *.ii *.s *.i *~ $(PROG) gmon.out

OBJS = database.o asnparse.o align.o matrices.o \
	stats.o hits.o query.o \
	search63.o search16.o search16s.o search7.o search7_ssse3.o \
        fasta.o adjusted.o SSW/src/ssw.o

DEPS = swipe.h Makefile

.SUFFIXES:.o .cc

swipe_sswlib : swipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) -g

# mpiswipe : mpiswipe.o $(OBJS)
# 	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(MPI_LINK)

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o : %.c $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# mpiswipe.o : swipe.cc $(DEPS)
# 	$(CXX) $(CXXFLAGS) -DMPISWIPE $(MPI_COMPILE) -c -o $@ swipe.cc

search7_ssse3.o : search7.cc $(DEPS)
	$(CXX) -mssse3 $(CXXFLAGS) -DSWIPE_SSSE3 -c -o $@ search7.cc
