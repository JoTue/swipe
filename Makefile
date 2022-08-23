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

export LD_LIBRARY_PATH=/apps/compat/8.6/lib64:/apps/compat/8.6/lib:/apps/python3/3.10.4/lib:/apps/lua/5.4.2/lib:/apps/gcc/12.1.0/lib64:/apps/gcc/12.1.0/lib:/apps/java/14.0.2/lib/server:/apps/R/4.2.0/lib64:/apps/slurm/22.05.2/lib/slurm:/apps/slurm/22.05.2/lib:/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/ReleaseMT/lib
# export LD_LIBRARY_PATH=/apps/compat/8.6/lib64:/apps/compat/8.6/lib:/apps/python3/3.10.4/lib:/apps/lua/5.4.2/lib:/apps/gcc/12.1.0/lib64:/apps/gcc/12.1.0/lib:/apps/java/14.0.2/lib/server:/apps/R/4.2.0/lib64:/apps/slurm/22.05.2/lib/slurm:/apps/slurm/22.05.2/lib:/scratch/cube/tuechler/swipe_LIBS/lib

# BLAST lib with compo_thresholds: /scratch/cube/tuechler/swipe_LIBS/custom_lib/3

MPI_COMPILE=`mpicxx --showme:compile`
MPI_LINK=`mpicxx --showme:link`

COMMON=-g
#COMMON=-pg -g

COMPILEOPT=-Wall -O3 -DCOMPO_ADJUSTMENT -DCOMPO_THRESHOLDS #-DSWLIB_8BIT #-DCOMPO_THRESHOLDS

LIBS=-lpthread -lcomposition_adjustment -lxblast -lxncbi

# Intel options
#CXX=icpc
#CXXFLAGS=$(COMPILEOPT) $(COMMON) -Wno-missing-declarations -fast
#LINKFLAGS=$(COMMON)

# GNU options
CXX=g++  # works with gcc/7.4.0
#CXXFLAGS=$(COMPILEOPT) $(COMMON) -I../ncbi-blast-2.2.29+-src/c++/src -I../ncbi-blast-2.2.29+-src/c++/include -I../ncbi-blast-2.2.29+-src/c++/DebugMT/inc -I../ncbi-blast-2.2.29+-src/c++/src/algo/blast/core
# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I../ncbi_blast/ncbi-blast-2.13.0+-src/c++/src -I../ncbi_blast/ncbi-blast-2.13.0+-src/c++/include -I../ncbi_blast/ncbi-blast-2.13.0+-src/c++/src/algo/blast/core
# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I../ncbi_blast/ncbi-blast-2.2.29+-src/c++/src -I../ncbi_blast/ncbi-blast-2.2.29+-src/c++/include -I../ncbi_blast/ncbi-blast-2.2.29+-src/c++/src/algo/blast/core
# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I/apps/ncbiblastplus/2.11.0/include/ncbi-tools++/ -I../ncbi_blast/ncbi-blast-2.2.29+-src/c++/src -I../ncbi_blast/ncbi-blast-2.2.29+-src/c++/src/algo/blast/core

# works with newest available libraries
# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I../swipe_LIBS/ncbi-tools++/ -I../swipe_LIBS/src -I../swipe_LIBS/src/algo/blast/core
# CXXFLAGS=$(COMPILEOPT) $(COMMON) -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/ncbi-blast-2.13.0+-src/c++/include -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/ncbi-blast-2.13.0+-src/c++/src -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/ncbi-blast-2.13.0+-src/c++/src/algo/blast/core -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/ncbi-blast-2.13.0+-src/c++/ReleaseMT/inc
CXXFLAGS=$(COMPILEOPT) $(COMMON) -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/include -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/src -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/src/algo/blast/core -I/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/ReleaseMT/inc

# LINKFLAGS=$(COMMON) -L../ncbi-blast-2.2.29+-src/c++/DebugMT/lib
# LINKFLAGS=$(COMMON) -L../swipe_LIBS/lib
# LINKFLAGS=$(COMMON) -L/scratch/cube/tuechler/simap2/string2020/clip/lib

LINKFLAGS=$(COMMON) -L/scratch/cube/tuechler/swipe_LIBS/custom_lib/3/ncbi-blast-2.13.0+-src/c++/ReleaseMT/lib
# LINKFLAGS=$(COMMON) -L../swipe_LIBS/lib

# PROG=swipe # mpiswipe
PROG=swipe_test # mpiswipe

all : $(PROG)

clean :
	rm -f *.o *.ii *.s *.i *~ $(PROG) gmon.out

OBJS = database.o asnparse.o align.o matrices.o \
	stats.o hits.o query.o \
	search63.o search16.o search16s.o search7.o search7_ssse3.o \
        fasta.o adjusted.o ssw.c

DEPS = swipe.h Makefile

.SUFFIXES:.o .cc

swipe_test : swipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) -g

mpiswipe : mpiswipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(MPI_LINK)

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mpiswipe.o : swipe.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -DMPISWIPE $(MPI_COMPILE) -c -o $@ swipe.cc

search7_ssse3.o : search7.cc $(DEPS)
	$(CXX) -mssse3 $(CXXFLAGS) -DSWIPE_SSSE3 -c -o $@ search7.cc
