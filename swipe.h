/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2013 Torbjorn Rognes, University of Oslo,
    Oslo University Hospital and Sencel Bioinformatics AS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include <iostream>
#include <chrono>  // for high_resolution_clock

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <inttypes.h>
#include <ctype.h>
#include <sys/stat.h>
#ifdef _WIN32
  #include <time.h>
#else
  #include <sys/times.h>
#endif
#include <fcntl.h>
#include <unistd.h>
#ifdef _WIN32
  #include "mman.h"
  #include <winsock.h>
#else
  #include <sys/mman.h>
  #include <arpa/inet.h>
#endif
#include <pthread.h>
#include <getopt.h>
#include <math.h>
#include <x86intrin.h>
#include <vector>

#ifdef MPISWIPE
#include <mpi.h>
#endif

#ifdef __APPLE__
#include <libkern/OSByteOrder.h>
#define bswap_32 OSSwapInt32
#define bswap_64 OSSwapInt64
#elif defined(_WIN32)
  #include "byteswap.h"
  #include "getpagesize.h"
#else
  #include <byteswap.h>
#endif

#if defined(__MINGW32__)
  #define free  _aligned_free
  #if !defined(__MINGW64_VERSION_MAJOR)
    #define _aligned_malloc __mingw_aligned_malloc
    #define _aligned_free  __mingw_aligned_free
  #endif // __MINGW64_VERSION_MAJOR
#endif // __MINGW32__

#ifdef COMPO_ADJUSTMENT
#include <algo/blast/core/blast_encoding.h> // for BLASTAA_SIZE
#include <algo/blast/core/blast_stat.h>
#include <algo/blast/composition_adjustment/composition_adjustment.h>
#endif // COMPO_ADJUSTMENT

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define SWIPE_VERSION "2.0.11_s1_sswlib"

// Should be 32bits integer
typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

#define WIDTH 32
#define WIDTH_SHIFT 5
#define BLOCKWIDTH 32

#define ext1 ".ssq"
#define ext2 ".ssi"
#define ext3 ".shd"
#define ext4 ".shi"

extern char BIAS;


//#define BIASED

#ifdef BIASED
#define ZERO 0x00
#else
#define ZERO 0x80
#endif

void vector_print(BYTE * vector);
void vector_print_word(WORD * vector);

void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);
void xfree(void *ptr);


extern long cpu_feature_ssse3;
extern long cpu_feature_sse41;

extern const char * queryname;
extern const char * matrixname;
extern long gapopen;
extern long gapextend;
extern long gapopenextend;
extern long symtype;
extern long matchscore;
extern long mismatchscore;
extern long totalhits;
extern const char * gencode_names[];
extern long querystrands;
extern double minexpect;
extern double expect;
extern long maxmatches;
extern unsigned int threads;
extern const char * databasename;
extern long alignments;
extern long queryno;
extern long compute7;
extern long show_taxid;
extern long effdbsize;
extern long mask;

extern char map_ncbi_nt4[];
extern char map_ncbi_nt16[];
extern char map_ncbi_aa[];
extern char map_sound[];

extern const char * sym_ncbi_nt4;
extern const char * sym_ncbi_nt16;
extern const char * sym_ncbi_nt16u;
extern const char * sym_ncbi_aa;
extern const char * sym_sound;

extern char ntcompl[];
extern char d_translate[];

extern FILE * out;

extern const char mat_blosum45[];
extern const char mat_blosum50[];
extern const char mat_blosum62[];
extern const char mat_blosum80[];
extern const char mat_blosum90[];
extern const char mat_pam30[];
extern const char mat_pam70[];
extern const char mat_pam250[];

extern long SCORELIMIT_7;
extern long SCORELIMIT_8;
extern long SCORELIMIT_16;
extern long SCORELIMIT_32;
extern long SCORELIMIT_63;
extern char BIAS;

extern char * score_matrix_7;
extern unsigned char * score_matrix_8;
extern int16_t * score_matrix_16;
extern int32_t * score_matrix_32;
extern int64_t * score_matrix_63;

struct sequence
{
  char * seq;
  char * seq_unmasked;
  long len;
};

struct query_s
{
  struct sequence nt[2]; /* 2 strands */
  struct sequence aa[6]; /* 6 frames */
  char * description;
  long dlen;
  long symtype;
  long strands;
  char * map;
  const char * sym;
#ifdef COMPO_ADJUSTMENT
  Blast_AminoAcidComposition composition;
  Blast_AminoAcidComposition composition_unmasked;
#endif // COMPO_ADJUSTMENT
};

extern struct query_s query;
//extern char * qseq;
//extern long qlen;

struct db_thread_s;

struct time_info
{
  time_t t1, t2;
#ifndef _WIN32
  struct tms times1, times2;
#endif // _WIN32
  clock_t wc1, wc2;
  long clk_tck;

  char * starttime;
  char * endtime;
  double elapsed;
  double speed;
};

extern struct time_info ti;

void fatal(const char * message);
void fatal(const char * format, const char * message);

void search7(BYTE * * q_start,
	     BYTE gap_open_penalty,
	     BYTE gap_extend_penalty,
	     BYTE * score_matrix,
	     BYTE * dprofile,
	     BYTE * hearray,
	     struct db_thread_s * dbt,
	     long sequences,
	     long * seqnos,
	     long * scores,
	     long qlen);

void search7_ssse3(BYTE * * q_start,
		   BYTE gap_open_penalty,
		   BYTE gap_extend_penalty,
		   BYTE * score_matrix,
		   BYTE * dprofile,
		   BYTE * hearray,
		   struct db_thread_s * dbt,
		   long sequences,
		   long * seqnos,
		   long * scores,
		   long qlen);

void search16(WORD * * q_start,
	      WORD gap_open_penalty,
	      WORD gap_extend_penalty,
	      WORD * score_matrix,
	      WORD * dprofile,
	      WORD * hearray,
	      struct db_thread_s * dbt,
	      long sequences,
	      long * seqnos,
	      long * scores,
	      long * bestpos,
	      int qlen);

void search16s(WORD * * q_start,
	       WORD gap_open_penalty,
	       WORD gap_extend_penalty,
	       WORD * score_matrix,
	       WORD * dprofile,
	       WORD * hearray,
	       struct db_thread_s * * dbta,
	       long sequences,
	       long * seqnos,
	       long * scores,
	       long * bestpos,
	       long * bestq,
	       int qlen);

long fullsw(char * dseq,
	    char * dend,
	    char * qseq,
	    char * qend,
	    long * hearray,
	    int64_t * score_matrix,
	    WORD gap_open_penalty,
	    WORD gap_extend_penalty);

void align(char * a_seq,
	   char * b_seq,
	   long M,
	   long N,
	   int64_t * scorematrix,
	   long q,
	   long r,
	   long * a_begin,
	   long * b_begin,
	   long * a_end,
	   long * b_end,
	   char ** alignment,
	   long * s);

void query_init(const char * queryname, long symtype, long strands);
void query_exit();
int query_read();
void query_show();

void score_matrix_init();
void score_matrix_free();
void score_matrix_dump();

void translate_init(long qtableno, long dtableno);
char * revcompl(char * seq, long len);
void translate(char * dna, long dlen,
               long strand, long frame, long table,
               char ** protp, long * plenp);

struct asnparse_info;
typedef struct asnparse_info * apt;

apt parser_create();
void parser_destruct(apt p);

long parse_header(apt p, unsigned char * buf, long len, long memb, long (*f)(long),
		  long show_gis, long indent, long maxlen,
		  long linelen, long maxdeflines, long show_descr);

void parse_getdeflines(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, long * deflines, char *** deflinetable);
void parse_gettitle(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, char ** title);

long parse_getdeflinecount(apt p, unsigned char * buf, long len,
                           long memb, long(*f_checktaxid)(long));

void db_open(long symtype, const char * basename, char * taxidfilename);
void db_close();
long db_getseqcount();
long db_getseqcount_masked();
long db_getsymcount();
long db_getsymcount_masked();
long db_getlongest();
char* db_gettitle();
char* db_gettime();
long db_getvolumecount();
long db_getseqcount_volume(long v);
long db_getseqcount_volume_masked(long v);
long db_ismasked();
long db_getversion();

long db_getvolume(long seqno);

struct db_thread_s * db_thread_create();
void db_thread_destruct(struct db_thread_s * t);

long db_check_taxid(long taxid);

void db_parse_header(struct db_thread_s * t, char * address, long length,
		     long show_gis,
		     long * deflines, char *** deflinetable);

void db_showheader(struct db_thread_s * t, char * address, long length,
		   long show_gis, long indent,
		   long maxlen, long linelen, long maxdeflines, long show_descr);
void db_getshowheader(struct db_thread_s * t, long seqno,
		      long show_gis, long indent,
		      long maxlen, long linelen, long maxdeflines);

void db_show_fasta(struct db_thread_s * t, long seqno,
		   long strand, long frame, long split);

long db_check_inclusion(struct db_thread_s * t, long seqno);

void db_mapsequences(struct db_thread_s * t, long firstseqno, long lastseqno);
void db_mapheaders(struct db_thread_s * t, long firstseqno, long lastseqno);

void db_getsequence(struct db_thread_s * t, long seqno, long strand, long frame,
		    char ** addressp, long * lengthp, long * ntlenp, int c);
void db_getheader(struct db_thread_s * t, long seqno, char ** address,
		  long * length);
void db_print_seq_map(char * address, long length, const char * map);


void hits_init(long descriptions, long alignments, long minscore,
	       long maxscore, double minexpect, double expect, int show_nostats);
void hits_enter(long seqno, long score, long qstrand, long qframe,
		long dstrand, long dframe, long align_hint, long bestq);
long * hits_sort();
long hits_getcount();
void hits_align(struct db_thread_s * t, long i);
void hits_show_begin(long view);
void hits_show_end(long view);
void hits_show(long view, long show_gis);
void hits_show_score_only();
void hits_empty();
void hits_exit();
void hits_gethit(long i, long * seqno, long * score,
		 long * qstrand, long * qframe,
		 long * dstrand, long * dframe);
void hits_getfull(long i,
		  long * seqno,
		  long * score,
		  long * align_q_start,
		  long * align_q_end,
		  long * align_d_start,
		  long * align_d_end,
		  char ** header, long * header_len,
		  char ** seq, long * seq_len,
		  char ** align, long * align_len);
void hits_enter_align_hint(long i, long q_end, long d_end);
void hits_enter_header(long i, char * header, long header_len);
#ifdef SWLIB_8BIT
void hits_enter_mat(long i, int8_t* mat, long aa_size);
int8_t* hits_getmat(long i);
#else
void hits_enter_mat(long i, int16_t* mat, long aa_size);
int16_t* hits_getmat(long i);
#endif // SWLIB_8BIT
char* hits_getseq(long i);
long hits_getdlen(long i);
void hits_enter_seq(long hitno, char* buffer, long len);
void hits_enter_align_coord(long i,
			    long align_q_start,
			    long align_q_end,
			    long align_d_start,
			    long align_d_end,
			    long dlennt);
void hits_enter_align_string(long hitno, char * align, long align_len);
void hits_defline_split(char * defline,
			long * gi,
			char ** link, int * linklen,
			char ** rest);

long stats_getparams_nt(long matchscore,
			long mismatchscore,
			long gopen,
			long gextend,
			double * lambda,
			double * K,
			double * H,
			double * alpha,
			double * beta);

long stats_getparams(const char * matrix,
		     long gopen,
		     long gextend,
		     double * lambda,
		     double * K,
		     double * H,
		     double * alpha,
		     double * beta);

long stats_getprefs(const char * matrix,
		    long * gopen,
		    long * gextend);


typedef int Int4;
typedef long Int8;
typedef double Nlm_FloatHi;

#ifdef COMPO_ADJUSTMENT

#define COMPOSITIONAL_MASK_COMP 1
#define COMPOSITIONAL_MASK_BOTH 2
#define COMPOSITIONAL_MASK_NONE 3
#define COMPOSITIONAL_MASK_SYMM 4
#define COMPOSITIONAL_MASK_BOTH_MATRIXONLY 5
#define COMPOSITIONAL_MASK_MATRIXONLY_SYMM 6

#define HIT_SUBJECT_QUERY_BEST_BL50 0b0001
#define HIT_SUBJECT_QUERY_BEST_BL62 0b0010
#define HIT_LARGE_SCORE             0b1000

#ifdef SWLIB_8BIT
static const int scaling_factor = 3;
#else
static const int scaling_factor = 32;
#endif // SWLIB_8BIT
static const int scaling_factor_BL62 = 32;
void compo_init(const char *matrixName, BlastScoreBlk **sbp, Blast_MatrixInfo **scaledMatrixInfo, int scaling_factor);
void compo_done(BlastScoreBlk **sbp, Blast_MatrixInfo **scaledMatrixInfo);
int compo_align(long *score_out, Blast_CompositionWorkspace * NRrecord, BlastScoreBlk *sbp, Blast_MatrixInfo *scaledMatrixInfo, int unmask, const Uint1 *data, int subject_length, long gapopen, long gapextend, int *matchStart, int *queryStart, int *matchEnd, int *queryEnd);
int compo_adjusted_matrix(Blast_CompositionWorkspace * NRrecord, BlastScoreBlk *sbp, Blast_MatrixInfo *scaledMatrixInfo, const Blast_AminoAcidComposition* query_composition, int query_length, const Blast_AminoAcidComposition* subject_composition, int subject_length, ECompoAdjustModes compo_adjust_mode, int stage, EMatrixAdjustRule *matrix_adjust_rule, double thresh_length, double thresh_distance, double thresh_angle);
using namespace std;
bool readFastaSequences(const char* dbFilePath, vector< vector<unsigned char> >* seqs);
void hits_set_align_string(long hitno, char * align, long score_align);
void hits_enter_score(long i, long score);
void hits_enter_adjusted_score(long i, long score, long score_blast, long score_blast_rev, long flags);
void count_align_matrix(long i, int64_t * score_matrix, const char *q_seq, long q_len, const char *d_seq, long d_len);
void show_align(long i);
//static bool preliminaryTestNearIdentical(int queryLength, int queryStart, int queryEnd, int matchStart, int matchEnd, int score, double cutoff);
// count matrix_adjust_rule
extern long mar2_0;
extern long mar2_4;
extern long mar2_other;
extern long mar3_0;
extern long mar3_4;
extern long mar3_other;
// db_composition functions
void enter_subject(long seqno, Blast_AminoAcidComposition subject_composition, char *subject_sequence_masked, long subject_length);
void get_subject_composition(long seqno, Blast_AminoAcidComposition* subject_composition);
char* get_subject_sequence_masked(long seqno);
void get_subject_length(long seqno, long* subject_length);
void compute_subject_compositions();
#endif // COMPO_ADJUSTMENT

#include "blastkar_partial.h"
