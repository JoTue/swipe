/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution
    Copyright (C) 2008-2014 Torbjorn Rognes, University of Oslo, 
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
// for align_adjusted3():
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#include "swipe.h"
#include "SSW/src/ssw.h"
#include <string>
#include <fstream>
#include <sstream>

#ifdef TIME_PAIRCOUNT
// time measuring
double time_stage1 = 0.0;
double time_stage2 = 0.0;
double time_stage2b = 0.0;
double time_stage2c = 0.0;
double time_writing = 0.0;
double time_writing_matrix = 0.0;
double time_writing_subject = 0.0;
double time_writing_query = 0.0;
double time_stage2d = 0.0;
double time_stage2_init = 0.0;
double time_stage2_sw = 0.0;
double time_stage2e = 0.0;
double time_stage2f = 0.0;
double time_stage3 = 0.0;
double time_work = 0.0;
// count number of pairs computed in each stage
long pairs1 = 0;
long pairs2 = 0;
long pairs2b = 0;
long pairs2c = 0;
long pairs3 = 0;
long pairs3b = 0;
long pairs3c = 0;
long pairs3_rev = 0;
long skipped_stage2 = 0;
#endif //TIME_PAIRCOUNT

/* ARGUMENTS AND THEIR DEFAULTS */

#define DEFAULT_MAXMATCHES 250
#define DEFAULT_ALIGNMENTS 100
#define DEFAULT_MINSCORE 1
#define DEFAULT_MAXSCORE (LONG_MAX)
#define DEFAULT_QUERYNAME "-"
#define DEFAULT_DATABASENAME ""
#define DEFAULT_GAPOPEN 0
#define DEFAULT_GAPEXTEND 0
#define DEFAULT_MATRIXNAME "BLOSUM62"
#define DEFAULT_MATCHSCORE 1
#define DEFAULT_MISMATCHSCORE (-3)
#define DEFAULT_THREADS 1
#define DEFAULT_VIEW 0
#define DEFAULT_SYMTYPE 1
#define DEFAULT_SHOW_GIS 0
#define DEFAULT_SHOW_TAXID 0
#define DEFAULT_EXPECT 10.0
#define DEFAULT_MINEXPECT 0.0
#define DEFAULT_QUERYSTRANDS 3
#define DEFAULT_QUERY_GENCODE 1
#define DEFAULT_DB_GENCODE 1
#define DEFAULT_SUBALIGNMENTS 1
#define DEFAULT_DUMP 0
#define DEFAULT_OUT stdout
#define DEFAULT_EFFDBSIZE 0
#define DEFAULT_MASK 0
#define DEFAULT_COMPO_ADJUST_MODE eCompositionMatrixAdjust
#define DEFAULT_NO_COMPO_STATS 0
#define DEFAULT_SKIP_STAGE2_SCORE 128
#define DEFAULT_THRESH_LENGTH 3.0
#define DEFAULT_THRESH_DISTANCE 0.16
#define DEFAULT_THRESH_ANGLE 70.0
#define DEFAULT_TEMPORARY_DIR ""
#define DEFAULT_SSW_BINARY_DIR "/scratch/cube/tuechler/swipe/SSW/bin"
#define DEFAULT_ALIGN_ADJUSTED_VERSION 1

char * progname;
const char * matrixname;
const char * databasename;
const char * queryname;
char * taxidfilename;
char * outfile = NULL;

double expect;
double minexpect;
long minscore, minscore2;
long maxscore;
long alignments;
long maxmatches;
long gapopen;
long gapextend;
unsigned int threads;
long view;
long symtype;
long show_gis;
long show_taxid;
long matchscore;
long mismatchscore;
long gapopenextend;
long querystrands;
long query_gencode;
long db_gencode;
long subalignments;
long dump;
long effdbsize;
long mask;
ECompoAdjustModes compo_adjust_mode;
long no_compo_stats;
long skip_stage2_score;
double thresh_length;
double thresh_distance;
double thresh_angle;
const char * temporary_dir;
const char * ssw_binary_dir;
long align_adjusted_version;

/* Other variables */

long queryno;

long cpu_feature_sse2;
long cpu_feature_ssse3;
long cpu_feature_sse41;

#define MAX_THREADS 256

pthread_mutex_t countmutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;

pthread_t pthread_id[MAX_THREADS];

#ifdef MPISWIPE
int mpirank, mpisize;
#endif

long maxchunksize;
long volnext;
long seqnext;
long * volchunks;
long * volseqs;

long compute7;
long compute16;
long compute32;
long compute63;

long rounds7;
long rounds16;
long rounds32;
long rounds63;

long totalhits;

FILE * out = stdout;

struct time_info ti;

struct search_data
{
  struct db_thread_s * dbt;
  struct db_thread_s * dbta[8];

  BYTE * dprofile;
  BYTE * hearray;
  BYTE ** qtable[6];

  long * scores;
  long * bestpos;
  long * bestq;
  long * start_list;
  long * in_list;
  long * out_list;
  long * tmp_list;

  long qlen[6];

  long start_count;
  long in_count;
  long out_count;

  long * start_hits;

  long seqfirst, seqlast;

  long qstrand1, qstrand2, qframe1, qframe2;
  long dstrand1, dstrand2, dframe1, dframe2;
};

void fatal(const char * message)
{
  if (message)
    fprintf(stderr, "%s\n", message);
  exit(1);
}

void fatal(const char * format, const char * message)
{
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;

  void * t = NULL;
#ifdef _WIN32
  t = _aligned_malloc(size, alignment);
#else

  posix_memalign(& t, alignment, size);
#endif // _WIN32

  if (t==NULL)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

void xfree(void *ptr)
{
#ifdef _WIN32
  _aligned_free(ptr);
#else
  free(ptr);
#endif // _WIN32
}

long alignedhits;
long alignedall;
long * hits_sorted;

long align_volnext;
long align_seqnext;

long * align_volseqs;
long * align_volchunks;

void align_init(struct search_data * sdp)
{
  sdp->dbt = db_thread_create();

  for(int i=0; i<8; i++)
    sdp->dbta[i] = db_thread_create();

  sdp->dprofile = (BYTE*) xmalloc(4*16*32);
  long qlen = 0;
  long hearraylen = 0;

  for(int i = 0; i < 6; i++)
  {
    sdp->qtable[i] = 0;
  }

  if (symtype == 0)
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
      {
	qlen = query.nt[s].len;
	sdp->qlen[3*s] = qlen;
	sdp->qtable[3*s] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	for(int i=0; i<qlen; i++)
	{
	  sdp->qtable[3*s][i] = sdp->dprofile + 16*query.nt[s].seq[i];
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
  }
  else if ((symtype == 1) || (symtype == 3) || (symtype == 5))
  {
    qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
    for(int i=0; i<qlen; i++)
    {
      sdp->qtable[0][i] = sdp->dprofile + 16*query.aa[0].seq[i];
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((symtype == 2) || (symtype == 4))
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
	for(int f=0; f<3; f++)
	{
	  qlen = query.aa[3*s+f].len;
	  sdp->qlen[3*s+f] = qlen;
	  sdp->qtable[3*s+f] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	  for(int i=0; i<qlen; i++)
	  {
	    sdp->qtable[3*s+f][i] = sdp->dprofile + 16*query.aa[3*s+f].seq[i];
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
  }

  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray = (BYTE*) xmalloc(hearraylen*32);

  long listsize = maxchunksize * sizeof(long);
  //  if ((symtype == 3) || (symtype == 4))
  //    listsize *= 6;

  sdp->start_list = (long*) xmalloc(listsize);
  sdp->start_hits = (long*) xmalloc(listsize);
  sdp->in_list = (long*) xmalloc(listsize);
  sdp->out_list = (long*) xmalloc(listsize);
  sdp->scores = (long*) xmalloc(listsize);
  sdp->bestpos = (long*) xmalloc(listsize);
  sdp->bestq = (long*) xmalloc(listsize);

  if (symtype == 0)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 2)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 3)
  {
    sdp->qstrand1 = 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = 0;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 1;
    sdp->dframe2 = 2;
  }
  else if (symtype == 4)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 1;
    sdp->dframe2 = 2;
  }
  else
  {
    sdp->qstrand1 = 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = 0;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
}

void align_chunk(struct search_data * sdp, long hitfirst, long hitlast)
{
  if (hitlast < alignments)
  {

    for (long qstrand = sdp->qstrand1; qstrand <= sdp->qstrand2; qstrand++)
      for(long qframe = sdp->qframe1; qframe <= sdp->qframe2; qframe++)
      {
	sdp->start_count = 0;

	for(long hitno = hitfirst; hitno <= hitlast; hitno++)
	{
	  long hs = hits_sorted[hitno];
	  long seqno, score, hqstrand, hqframe, hdstrand, hdframe;

	  hits_gethit(hs, & seqno, & score, & hqstrand, & hqframe,
		      & hdstrand, & hdframe);

	  if ((qstrand == hqstrand) && (qframe == hqframe))
	  {
	    sdp->start_hits[sdp->start_count] = hs;
	    sdp->start_list[sdp->start_count] =
	      (seqno << 3) | (hdstrand << 2) | hdframe;
	    sdp->start_count++;
	  }
	}

	if (sdp->start_count)
	{
	  //	  printf("Aligning %ld sequences.\n", sdp->start_count);


	  BYTE ** qtable = sdp->qtable[3*qstrand+qframe];
	  long qlen = sdp->qlen[3*qstrand+qframe];

	  /* 16-bit search, 8x1 db symbols, with alignment end */

	  pthread_mutex_lock(&countmutex);
	  compute32 += sdp->in_count;
	  rounds32++;
	  pthread_mutex_unlock(&countmutex);

	  search16s((WORD**)qtable,
		    gapopenextend,
		    gapextend,
		    (WORD*)(score_matrix_16),
		    (WORD*)(sdp->dprofile),
		    (WORD*)(sdp->hearray),
		    sdp->dbta,
		    sdp->start_count,
		    sdp->start_list,
		    sdp->scores,
		    sdp->bestpos,
		    sdp->bestq,
		    qlen);

	  for (int i=0; i<sdp->start_count; i++)
	  {
	    long pos = sdp->bestpos[i];
	    long bestq = sdp->bestq[i];

	    //	  fprintf(out, "seqno=%ld score=%ld bestpos=%ld\n", seqno, score, pos);

	    long hitno = sdp->start_hits[i];

	    if (sdp->scores[i] < SCORELIMIT_16)
	      hits_enter_align_hint(hitno, bestq, pos);
	  }
	}
      }

  }

  for(long hitno = hitfirst; hitno <= hitlast; hitno++)
    hits_align(sdp->dbt, hits_sorted[hitno]);
}

void align_done(struct search_data * sdp)
{
  for(int i = 0; i < 6; i++)
  {
    if (sdp->qtable[i])
      free(sdp->qtable[i]);
  }

  free(sdp->dprofile);
  free(sdp->hearray);
  free(sdp->scores);
  free(sdp->bestpos);
  free(sdp->bestq);
  free(sdp->start_list);
  free(sdp->start_hits);
  free(sdp->in_list);
  free(sdp->out_list);

  for(int i=0; i<8; i++)
    db_thread_destruct(sdp->dbta[i]);

  db_thread_destruct(sdp->dbt);
}


void calc_chunks(long volcount,
		 long par,
		 long channels,
		 long * volseqs,
		 long * volchunks,
		 long * totalchunks,
		 long * biggestchunk)
{
#ifdef DEBUG
  fprintf(out, "Calculating chunk distribution among volumes\n");
#endif

  long volsused = 0;
  long chunksizes[volcount];
  long totalseqs = 0;
  long maxchunksize = 0;
  long vv = 0;
  for(long v = 0; v < volcount; v++)
  {
    if (volseqs[v])
    {
      totalseqs += volseqs[v];
      volsused++;
      volchunks[v] = 1;
      chunksizes[v] = volseqs[v];
      if (chunksizes[v] > maxchunksize)
      {
	maxchunksize = chunksizes[v];
	vv = v;
      }
    }
    else
    {
      chunksizes[v] = 0;
      volchunks[v] = 0;
    }
  }

  long upper = channels;
  if (totalseqs >= 4 * channels * par)
    upper *= (long) (floor(sqrt((1.0 * totalseqs) / (channels * par))));

  long chunks = volsused;
  long minchunks = totalseqs < par ? totalseqs : par;

  while((maxchunksize > upper) || (chunks < minchunks))
  {
    volchunks[vv]++;
    chunks++;
    chunksizes[vv] = (volseqs[vv] + volchunks[vv] - 1) / volchunks[vv];

    maxchunksize = 0;
    vv = 0;
    for(long v=0; v < volcount; v++)
    {
      if (chunksizes[v] > maxchunksize)
      {
	vv = v;
	maxchunksize = chunksizes[v];
      }
    }
  }

#ifdef DEBUG
#ifdef MPISWIPE
  if (!mpirank)
#endif
  {
    fprintf(out, "\n");
    fprintf(out, "Chunk distribution:\n");
    fprintf(out, "Total seqs: %ld\n", totalseqs);
    fprintf(out, "Total chunks: %ld  Upper size limit: %ld\n", chunks, upper);
    fprintf(out, "Vol Sequences Chunks Maxchunksize\n");
    for(long v = 0; v < volcount; v++)
    {
      fprintf(out, "%3ld %9ld %6ld %12ld\n",
	      v, volseqs[v], volchunks[v], chunksizes[v]);
    }
    fprintf(out, "\n");
  }
#endif

  *biggestchunk = maxchunksize;
  *totalchunks = chunks;
}

void align_threads_init()
{
  long hits = hits_getcount();

  hits_sorted = hits_sort();

  long bins = 7;

  align_volseqs = (long*) xmalloc(bins*sizeof(long));
  align_volchunks = (long*) xmalloc(bins*sizeof(long));

  for(long i = 0; i<bins; i++)
    align_volseqs[i] = 0;

  for(long i = 0; i<hits; i++)
  {
    long seqno;
    long score;
    long qstrand;
    long qframe;
    long dstrand;
    long dframe;

    if (i>=alignments)
      align_volseqs[6]++;
    else
    {
      hits_gethit(i, & seqno, & score,
		  & qstrand, & qframe,
		  & dstrand, & dframe);

      align_volseqs[3*qstrand+qframe]++;
    }
  }

  long totalchunks;

  calc_chunks(bins,
	      threads,
	      8,
	      align_volseqs,
	      align_volchunks,
	      & totalchunks,
	      & maxchunksize);

  alignedhits = 0;
  align_volnext = 0;

  while ((align_volnext < bins) && (align_volchunks[align_volnext] == 0))
    align_volnext++;
}

void align_threads_done()
{
  free(hits_sorted);
  free(align_volchunks);
  free(align_volseqs);
}

int align_getwork(long * first, long * last)
{
  int status = 0;
  long bins = 7;
  long volcount = bins;

  pthread_mutex_lock(&workmutex);
  if (align_volnext < volcount)
  {
    long seqcount = align_volseqs[align_volnext];
    long chunks = align_volchunks[align_volnext];
    long chunksize = ((seqcount+chunks-1) / chunks);

    * first = alignedhits;
    * last = alignedhits + chunksize - 1;

    alignedhits += chunksize;
    status = 1;

    align_volseqs[align_volnext] -= chunksize;
    align_volchunks[align_volnext]--;

    while ((align_volnext < bins) && (align_volchunks[align_volnext] == 0))
      align_volnext++;
  }
  pthread_mutex_unlock(&workmutex);
  return status;
}

void * align_worker(void *)
{
  search_data sd;
  align_init(&sd);

  long i, j;
  while(align_getwork(&i, &j))
    align_chunk(&sd, i, j);

  align_done(&sd);
  return 0;
}

void align_threads()
{
  unsigned int t;
  void * status;

  align_threads_init();

  for(t=0; t<threads; t++)
    {
      if (pthread_create(pthread_id + t, 0, align_worker, &t)) // no thread safety issues because the worker does not access arg
	fatal("Cannot create thread.");
    }

  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      fatal("Cannot join thread.");
  }

  align_threads_done();
}

//#ifdef COMPO_ADJUSTMENT
BlastScoreBlk *sbp;
Blast_MatrixInfo *scaledMatrixInfo;
Blast_CompositionWorkspace *NRrecord;

BlastScoreBlk *sbpBL62;
Blast_MatrixInfo *scaledMatrixInfoBL62;
Blast_CompositionWorkspace *NRrecordBL62;

using namespace std;
vector< vector<unsigned char> >* dbSequences = new vector< vector<unsigned char> >();
int64_t *adjusted_score_matrix_63, *adjusted_score_matrix_63_rev, *adjusted_score_matrix_63_fullsw;
#ifdef SWLIB_8BIT
uint8_t minScoreSWlib;
#else
uint16_t minScoreSWlib;
#endif // SWLIB_8BIT

int matrix_min = 0;
int matrix_max = 0;

int16_t* mata;
int8_t* mat8;
s_align* result = NULL;
alignment_end ae = {0};
alignment_end* bests;
BYTE* hearray = NULL;
long hesize = query.aa[0].len*32;

void printSequence(const unsigned char* sequence, int len) {
  for (int i = 0; i < len; i++)
    putchar(NCBISTDAA_TO_AMINOACID[sequence[i]]);
  putchar('\n');
  putchar('\n');
}

void print_matrix(int64_t* matrix) {
    for (int a = 0; a < 32; a++) { // ms should be sbps->matrix->nrows
        for (int b = 0; b < 32; b++)
            fprintf(stderr, "%4ld ", matrix[(a<<5) + b]);
        fprintf(stderr, "\n");
    }
}

void align_adjusted_init() {
  if (mask) {
    readFastaSequences(databasename, dbSequences);
  //if(dump)
/*      for (unsigned int i = 0; i < dbSequences->size(); i++) {
        unsigned char* seq = (*dbSequences)[i].data();
        printSequence(seq, (*dbSequences)[i].size());
      }
*/}
  NRrecord = Blast_CompositionWorkspaceNew();
  int status_code = Blast_CompositionWorkspaceInit(NRrecord, matrixname);
  if (status_code != 0) {
    fprintf(stderr, "Blast_CompositionWorkspaceInit error: %d\n", status_code);
    exit(1);
  }
  compo_init(matrixname, &sbp, &scaledMatrixInfo, scaling_factor);

  // for BLOSUM62
  NRrecordBL62 = Blast_CompositionWorkspaceNew();
  status_code = Blast_CompositionWorkspaceInit(NRrecordBL62, "BLOSUM62");
  if (status_code != 0) {
    fprintf(stderr, "Blast_CompositionWorkspaceInit error: %d\n", status_code);
    exit(1);
  }
  compo_init("BLOSUM62", &sbpBL62, &scaledMatrixInfoBL62, scaling_factor_BL62);

  adjusted_score_matrix_63 = (int64_t *) xmalloc(32*32*sizeof(int64_t));
  memset(adjusted_score_matrix_63, -1, 32*32*8);
  adjusted_score_matrix_63_rev = (int64_t *) xmalloc(32*32*sizeof(int64_t));
  adjusted_score_matrix_63_fullsw = (int64_t *) xmalloc(32*32*sizeof(int64_t));
  //print_matrix(score_matrix_63);
  
  //fprintf(out, "SCORELIMIT_7: %ld\n", SCORELIMIT_7); // 113 for BL50 and 117 for BL62
  mata = (int16_t*)calloc(BLASTAA_SIZE*BLASTAA_SIZE, sizeof(int16_t));
  mat8 = (int8_t*)calloc(BLASTAA_SIZE*BLASTAA_SIZE, sizeof(int8_t));

//  int8_t* read_reverse = 0;
//  long reverse_len = 0;
  
//  int err = 0;
//  s_align* result = NULL;
  result = (s_align*)calloc(1, sizeof(s_align));

#ifdef SWLIB_8BIT
  minScoreSWlib = 255;
#else
  minScoreSWlib = 32767;
#endif // SWLIB_8BIT
  long tmpMinScore = minscore2 * scaling_factor;
  if (tmpMinScore > scaling_factor && tmpMinScore < minScoreSWlib)
    minScoreSWlib = tmpMinScore;
 
  if (mask == COMPOSITIONAL_MASK_SYMM)
    hearray = (BYTE*)xmalloc(hesize);
}

// -i ../swimd/test_data/db/uniprot_sprot15.fasta -d ../swimd/test_data/db/uniprot_sprot15.fasta -m 88 -b 1000 -v 1000 -e 10000 -M BLOSUM50 -G 13 -E 2 -s
// -i test/test_c.fa -d test/test_mcw_masked.fa -e 10000000 -m 88 -b 10 -v 1 -M BLOSUM50 -G 13 -E 2 -B 38 -s
// -i test/RFCL_PYRNV.fa -d test/YOCA_CAEEL.fa -e 10000000 -m 88 -b 10 -v 1 -M BLOSUM50 -G 13 -E 2 -c 50 -B 80 -s4

char map_ncbi_int2aa[28] = 
  {
    '-', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
    'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'U', '*', 'O', 'J'
  };

string int2aa(char* num_seq, long seqlen) 
{
  string aa_seq = "";
  for (int i = 0; i < seqlen; i++) {
    aa_seq = aa_seq + map_ncbi_int2aa[num_seq[i]];
  }
  return aa_seq;
}

// default and fastest align_adjusted version:
void align_adjusted() {
  long hits = hits_getcount();

  struct db_thread_s * dbt = db_thread_create();

  char *subject_sequence_masked, *subject_sequence_unmasked;
  long subject_length;
  char *query_sequence_masked = query.aa[0].seq;
  char *query_sequence_unmasked = (mask ? query.aa[0].seq_unmasked : query.aa[0].seq);
  long query_length = query.aa[0].len;
  char *qseq, *dseq;
  long qlen, dlen;
  long ntlen;
  long seqno;
  long score;
  long qstrand;
  long qframe;
  long dstrand;
  long dframe;
  long adjusted_score;
  long matchStart, queryStart, matchEnd, queryEnd, score_align;
  //long matchStart_r, queryStart_r, matchEnd_r, queryEnd_r, score_align_r;
  char * alignment;
  long gap_open = gapopen * scaling_factor;
  long gap_extend = gapextend * scaling_factor;
  long gap_open_extend = gap_open + gap_extend;
  long gapopen_BlastDef = 11 * scaling_factor_BL62;
  long gapextend_BlastDef = 1 * scaling_factor_BL62;
  long gapopenextend_BlastDef = gapopen_BlastDef + gapextend_BlastDef;
  EMatrixAdjustRule matrix_adjust_rule;
  ECompoAdjustModes compo_adjust_mode_stage3 = eCompositionMatrixAdjust;
  int16_t* mat;
  char dseq2;
  s_profile* p = NULL;
  int32_t maskLen = query_length / 2;
  
  for(long i = 0; i<hits; i++)
  {
    long flags = 0;
    bool symm_check = false;
    long adjusted_score_blast = 0, adjusted_score_blast_rev = 0;
    hits_gethit(i, & seqno, & score,
		  & qstrand, & qframe,
		  & dstrand, & dframe);

    if (score < minscore)
      continue;

    // stage 2
#ifdef TIME_PAIRCOUNT
#ifdef COMPO_THRESHOLDS
    printf("stage2\n");
#endif //COMPO_THRESHOLDS
    pairs2++;
    auto begin2 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    matchStart = queryStart = matchEnd = queryEnd = 0;
    Blast_AminoAcidComposition subject_composition, subject_composition_unmasked;

    // // OLD:
    // db_mapsequences(dbt, seqno, seqno);
    // db_getsequence(dbt, seqno, dstrand, dframe, &subject_sequence_masked, &subject_length, &ntlen, 0);
    // subject_length--;
    // if (mask)
    //   subject_sequence_unmasked = (char*)(*dbSequences)[seqno].data();
    // else
    //   subject_sequence_unmasked = subject_sequence_masked;

    // Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, (const Uint1*)subject_sequence_masked, subject_length);
    // if (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
    //   Blast_ReadAaComposition(&subject_composition_unmasked, BLASTAA_SIZE, (const Uint1*)subject_sequence_unmasked, subject_length);
    //   symm_check = (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) && (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
    // }

    // NEW:
    get_subject_composition(seqno, &subject_composition);
    get_subject_composition_unmasked(seqno, &subject_composition_unmasked);
    subject_sequence_masked = get_subject_sequence_masked(seqno);
    get_subject_length(seqno, &subject_length);
    if (mask)
      subject_sequence_unmasked = (char*)(*dbSequences)[seqno].data();
    else
      subject_sequence_unmasked = subject_sequence_masked;

#ifdef TIME_PAIRCOUNT
    // time init, aa composition
    auto end2b = std::chrono::high_resolution_clock::now();
    time_stage2b += (end2b - begin2).count();
#endif //TIME_PAIRCOUNT

    if (compo_adjust_mode == eNoCompositionBasedStats || score >= skip_stage2_score) {

#ifdef TIME_PAIRCOUNT
       // time stage 2
      auto end2 = std::chrono::high_resolution_clock::now();
      time_stage2 += (end2 - begin2).count();
      skipped_stage2++;
#endif //TIME_PAIRCOUNT

      adjusted_score = score * scaling_factor;
    } else {

#ifdef TIME_PAIRCOUNT
    auto begin2c = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT
    
    matrix_adjust_rule = eDontAdjustMatrix;  // TODO REMOVE?
    if (compo_adjust_mode != eNoCompositionBasedStats){
      compo_adjusted_matrix(NRrecord, sbp, scaledMatrixInfo,
        mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
        mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode, 2, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
    }

    if (matrix_adjust_rule == 0 && compo_adjust_mode == eCompositionMatrixAdjust && no_compo_stats == 1) {
       // skip stage 2 alignment

#ifdef TIME_PAIRCOUNT
      auto end2c = std::chrono::high_resolution_clock::now();
      time_stage2c += (end2c - begin2c).count();
      auto end2 = std::chrono::high_resolution_clock::now();
      time_stage2 += (end2 - begin2).count();
#endif //TIME_PAIRCOUNT

      adjusted_score = score * scaling_factor;
    } else {

    // convert adjusted score matrix to sswlib matrix
    for (int a = 0; a < BLASTAA_SIZE; a++)
      for (int b = 0; b < BLASTAA_SIZE; b++) {
#ifdef SWLIB_8BIT
        mat8[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
        // std::cout << sbp->matrix->data[b][a] << " ";
#else
        mata[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
        // std::cout << sbp->matrix->data[b][a] << " ";
#endif // SWLIB_8BIT
    }

#ifdef TIME_PAIRCOUNT
    // time 1st comp. adjustment
    auto end2c = std::chrono::high_resolution_clock::now();
    time_stage2c += (end2c - begin2c).count();
    auto begin2d = std::chrono::high_resolution_clock::now();
    auto begin2_init = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

#ifdef SWLIB_8BIT
    p = ssw_init((const int8_t*)(mask != COMPOSITIONAL_MASK_BOTH ? query_sequence_unmasked : query_sequence_masked), query_length, mat8, BLASTAA_SIZE, 0);

#ifdef TIME_PAIRCOUNT
    // time ssw_init_word
    auto end2_init = std::chrono::high_resolution_clock::now();
    time_stage2_init += (end2_init - begin2_init).count();
    auto begin2_sw = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT
    
    sw_sse2_byte2((const int8_t*)(mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? subject_sequence_unmasked : subject_sequence_masked), 0, subject_length, p->readLen, gap_open_extend, gap_extend, p->profile_byte, minScoreSWlib, p->bias, maskLen, &ae);

#ifdef TIME_PAIRCOUNT
    // time sw
    auto end2_sw = std::chrono::high_resolution_clock::now();
    time_stage2_sw += (end2_sw - begin2_sw).count();
#endif //TIME_PAIRCOUNT
    
    result->score1 = ae.score;
    result->ref_end1 = ae.ref;
    result->read_end1 = ae.read;
#else
    p = ssw_init_word((const int8_t*)(mask != COMPOSITIONAL_MASK_BOTH ? query_sequence_unmasked : query_sequence_masked), query_length, mata, BLASTAA_SIZE);

#ifdef TIME_PAIRCOUNT
    // time ssw_init_word
    auto end2_init = std::chrono::high_resolution_clock::now();
    time_stage2_init += (end2_init - begin2_init).count();
    auto begin2_sw = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    sw_sse2_word2((const int8_t*)(mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? subject_sequence_unmasked : subject_sequence_masked), 0, subject_length, p->readLen, gap_open_extend, gap_extend, p->profile_word, minScoreSWlib, maskLen, &ae);

#ifdef TIME_PAIRCOUNT    
    // time sw
    auto end2_sw = std::chrono::high_resolution_clock::now();
    time_stage2_sw += (end2_sw - begin2_sw).count();
#endif //TIME_PAIRCOUNT

    result->score1 = ae.score;
    result->ref_end1 = ae.ref;
    result->read_end1 = ae.read;
#endif // SWLIB_8BIT
    if (p)
       init_destroy(p);

#ifdef TIME_PAIRCOUNT
    // time 1st alignment
    auto end2d = std::chrono::high_resolution_clock::now();
    time_stage2d += (end2d - begin2d).count();
    auto begin2e = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT
    
    if (symm_check && result->score1 < minscore2 * scaling_factor) {

#ifdef TIME_PAIRCOUNT
      pairs2b++;
#endif //TIME_PAIRCOUNT

      matrix_adjust_rule = eDontAdjustMatrix;
      compo_adjusted_matrix(NRrecord, sbp, scaledMatrixInfo, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
      // convert adjusted score matrix to swipe matix
      for (int a = 0; a < BLASTAA_SIZE; a++)
        for (int b = 0; b < BLASTAA_SIZE; b++) {
#ifdef SWLIB_8BIT
          mat8[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
#else
          mata[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
#endif // SWLIB_8BIT
        }

#ifdef TIME_PAIRCOUNT
      // time reverse comp. adjustment
      auto end2e = std::chrono::high_resolution_clock::now();
      time_stage2e += (end2e - begin2e).count();
      auto begin2f = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

#ifdef SWLIB_8BIT
      p = ssw_init((const int8_t*)subject_sequence_unmasked, subject_length, mat8, BLASTAA_SIZE, 0);
      sw_sse2_byte2((const int8_t*)(mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? query_sequence_unmasked : query_sequence_masked), 0, query_length, p->readLen, gap_open_extend, gap_extend, p->profile_byte, minScoreSWlib, p->bias, maskLen, &ae);
#else
      p = ssw_init_word((const int8_t*)subject_sequence_unmasked, subject_length, mata, BLASTAA_SIZE);
      sw_sse2_word2((const int8_t*)(mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? query_sequence_unmasked : query_sequence_masked), 0, query_length, p->readLen, gap_open_extend, gap_extend, p->profile_word, minScoreSWlib, maskLen, &ae);
#endif // SWLIB_8BIT
      if (p)
        init_destroy(p);
      if (ae.score > result->score1) {

#ifdef TIME_PAIRCOUNT
        pairs2c++;
#endif //TIME_PAIRCOUNT

        result->score1 = ae.score;
        flags |= HIT_SUBJECT_QUERY_BEST_BL50;
      }

#ifdef TIME_PAIRCOUNT
      // time reverse alignment
      auto end2f = std::chrono::high_resolution_clock::now();
      time_stage2f += (end2e - begin2e).count();
#endif //TIME_PAIRCOUNT

    }
    adjusted_score = result->score1;

#ifdef TIME_PAIRCOUNT
    // time stage 2
    auto end2 = std::chrono::high_resolution_clock::now();
    time_stage2 += (end2 - begin2).count();
#endif //TIME_PAIRCOUNT

    }
    }
    
    if (adjusted_score >= minscore2 * scaling_factor) {
      // stage 3
#ifdef TIME_PAIRCOUNT
#ifdef COMPO_THRESHOLDS
      printf("stage3\n");
#endif //COMPO_THRESHOLDS
      pairs3++;
      auto begin3 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT
      
      matrix_adjust_rule = eDontAdjustMatrix;
      if (compo_adjust_mode_stage3 != eNoCompositionBasedStats){
        compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, 
          mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
          mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode_stage3, 3, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
      }

      // convert adjusted score matrix to swipe matix
      for (int a = 0; a < BLASTAA_SIZE; a++)
        for (int b = 0; b < BLASTAA_SIZE; b++) {
          mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
          adjusted_score_matrix_63[(a<<5) + b] = sbpBL62->matrix->data[a][b];
          adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
        }
      p = ssw_init_word((const int8_t*)query_sequence_unmasked, query_length, mata, BLASTAA_SIZE);
      sw_sse2_word2((const int8_t*)subject_sequence_unmasked, 0, subject_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);
    
      if (p)
        init_destroy(p);
      result->score1 = ae.score;
      result->ref_end1 = ae.ref;
      result->read_end1 = ae.read;
      if (result->score1 == 32767) {
        flags |= HIT_LARGE_SCORE;
      }
      
      matchStart = queryStart = 0;
      qseq = query_sequence_unmasked;
      qlen = query_length;
      dseq = subject_sequence_unmasked;
      dlen = subject_length;
      hits_enter_seq(i, dseq, dlen);
      
      if (symm_check) {

#ifdef TIME_PAIRCOUNT
        pairs3b++;
#endif //TIME_PAIRCOUNT

        if (result->score1 == 32767) {
          if (hesize < dlen * 32) {
            hesize = dlen*32;
            hearray = (BYTE*)realloc(hearray, hesize);
          }
          //TODO: replace fullsw by an int version of swlib
          result->score1 = fullsw(dseq, dseq + dlen, qseq, qseq + qlen, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
        }
                  
        matrix_adjust_rule = eDontAdjustMatrix;
        compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode_stage3, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
        // convert adjusted score matrix to swipe matix
        for (int a = 0; a < BLASTAA_SIZE; a++)
          for (int b = 0; b < BLASTAA_SIZE; b++) {
            mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
            adjusted_score_matrix_63_rev[(a<<5) + b] = sbpBL62->matrix->data[a][b];
            adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
          }
        p = ssw_init_word((const int8_t*)subject_sequence_unmasked, subject_length, mata, BLASTAA_SIZE);
        sw_sse2_word2((const int8_t*)query_sequence_unmasked, 0, query_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);
        if (p)
          init_destroy(p);
        result->score2 = ae.score;
        result->ref_end2 = ae.ref;
        result->read_end2 = ae.read;
        if (result->score2 == 32767) {
          flags |= HIT_LARGE_SCORE;
          char *qseq2 = subject_sequence_unmasked;
          long qlen2 = subject_length;
          char *dseq2 = query_sequence_unmasked;
          long dlen2 = query_length;
          if (hesize < dlen2 * 32) {
            hesize = dlen2*32;
            hearray = (BYTE*)realloc(hearray, hesize);
          }
          //TODO: replace fullsw by an int version of swlib
          result->score2 = fullsw(dseq2, dseq2 + dlen2, qseq2, qseq2 + qlen2, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
        }
//        if (round(result->score2/scaling_factor_BL62) > round(result->score1/scaling_factor_BL62)) {
//        if (round((float)result->score2/scaling_factor_BL62) > round((float)result->score1/scaling_factor_BL62)) {
        if (result->score2 > result->score1) {

#ifdef TIME_PAIRCOUNT
          pairs3c++;
#endif //TIME_PAIRCOUNT

          flags |= HIT_SUBJECT_QUERY_BEST_BL62;
          //fprintf(stderr, "symmetrical score is larger than adjusted blastp score: %d vs %d\n", (int)round(result->score2/scaling_factor_BL62), (int)round(result->score1/scaling_factor_BL62));
        }
        adjusted_score_blast_rev = result->score2;
      }
      
      int64_t *used_score_matrix_63 = adjusted_score_matrix_63;
      
      if (flags & HIT_LARGE_SCORE)
        matchEnd = queryEnd = score_align = 0;
      else {
        if (flags & HIT_SUBJECT_QUERY_BEST_BL62) {
          adjusted_score_blast = result->score2;
          matchEnd = result->ref_end2;
          queryEnd = result->read_end2;
          adjusted_score_blast_rev = result->score1;
          used_score_matrix_63 = adjusted_score_matrix_63_rev;
          qseq = subject_sequence_unmasked;
          qlen = subject_length;
          dseq = query_sequence_unmasked;
          dlen = query_length;
        } else {
          adjusted_score_blast = result->score1;
          matchEnd = result->ref_end1;
          queryEnd = result->read_end1;
        }
        score_align = adjusted_score_blast;
      }
      //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
/*      fprintf(out, "query sequence: \n");
      db_print_seq_map(qseq, qlen, sym_ncbi_aa);
      fprintf(out, "subject sequence: \n");
      db_print_seq_map(dseq, dlen, sym_ncbi_aa);
*/      align(qseq,
        dseq,
        qlen,
        dlen,
        used_score_matrix_63,
        gapopen_BlastDef,
        gapextend_BlastDef,
        & queryStart,
        & matchStart,
        & queryEnd,
        & matchEnd,
        & alignment,
        & score_align);
        adjusted_score_blast = score_align;
//          fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
/*        } else {
        score_align = adjusted_score_blast;
        // reverse query sequence; queryEnd is 0-based alignment ending position
        int32_t start = 0, end = queryEnd;
        if (end > reverse_len) {
          if (reverse_len == 0)
            read_reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));
          else
            read_reverse = (int8_t*)realloc(read_reverse, (end + 1) * sizeof(int8_t));
          reverse_len = end;
        }
        while (start <= end) {
      read_reverse[start] = qseq[end];
          read_reverse[end] = qseq[start];
          ++ start;
          -- end;
        }
        __m128i* vP = qP_word16(read_reverse, mata, queryEnd + 1, p->n);
        alignment_end* bests_reverse = sw_sse2_word2((const int8_t*)subject_sequence, 1, matchEnd + 1, queryEnd + 1, gapopen_BlastDef + gapextend_BlastDef, gapextend_BlastDef, vP, score_align, maskLen);
        free(vP);
        matchStart = bests_reverse[0].ref;
        queryStart = queryEnd - bests_reverse[0].read;
        result->score1 = bests_reverse[0].score;
        free(bests_reverse);
        matchStart_r = matchEnd_r = queryStart_r = queryEnd_r = score_align_r = 0;
      align(&qseq[queryStart],
        &subject_sequence[matchStart],
        queryEnd - queryStart + 1,
        matchEnd - matchStart + 1,
        adjusted_score_matrix_63,
        gapopen_BlastDef,
        gapextend_BlastDef,
        & queryStart_r,
        & matchStart_r,
        & queryEnd_r,
        & matchEnd_r,
        & alignment,
        & score_align_r);
        //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", queryStart, queryEnd, matchStart, matchEnd, score_align_r, result->score1);
      }*/
      hits_set_align_string(i, alignment, score_align);
      hits_enter_align_coord(i, queryStart, queryEnd, matchStart, matchEnd, 0);
//        score_matrix_63 = used_score_matrix_63;
//        score_matrix_dump();
      count_align_matrix(i, used_score_matrix_63, qseq, qlen, dseq, dlen);
//        show_align(i);

#ifdef TIME_PAIRCOUNT
      // time stage 3
      auto end3 = std::chrono::high_resolution_clock::now();
      time_stage3 += (end3 - begin3).count();
#endif //TIME_PAIRCOUNT

    } else
      hits_enter_align_coord(i, queryStart, queryEnd, matchStart, matchEnd, 0);
    hits_enter_adjusted_score(i, adjusted_score, adjusted_score_blast, adjusted_score_blast_rev, flags);
  }
  db_thread_destruct(dbt);
}

// align all hits for one query at once, calling SSW via command line
void align_adjusted4() {
  long hits = hits_getcount();

  struct db_thread_s * dbt = db_thread_create();

  char *subject_sequence_masked, *subject_sequence_unmasked;
  long subject_length;
  char *query_sequence_masked = query.aa[0].seq;
  char *query_sequence_unmasked = (mask ? query.aa[0].seq_unmasked : query.aa[0].seq);
  long query_length = query.aa[0].len;
  char *qseq, *dseq;
  long qlen, dlen;
  long ntlen;
  long seqno;
  long score;
  long qstrand;
  long qframe;
  long dstrand;
  long dframe;
  long adjusted_score;
  long matchStart, queryStart, matchEnd, queryEnd, score_align;
  //long matchStart_r, queryStart_r, matchEnd_r, queryEnd_r, score_align_r;
  char * alignment;
  long gap_open = gapopen * scaling_factor;
  long gap_extend = gapextend * scaling_factor;
  long gap_open_extend = gap_open + gap_extend;
  long gapopen_BlastDef = 11 * scaling_factor_BL62;
  long gapextend_BlastDef = 1 * scaling_factor_BL62;
  long gapopenextend_BlastDef = gapopen_BlastDef + gapextend_BlastDef;
  EMatrixAdjustRule matrix_adjust_rule;
  ECompoAdjustModes compo_adjust_mode_stage3 = eCompositionMatrixAdjust;
#ifdef SWLIB_8BIT
  int8_t* mat;
#else
  int16_t* mat;
#endif // SWLIB_8BIT
  char *dseq2;
  string aa_seq;
  std::vector< int > stage3_indices;
  Blast_AminoAcidComposition subject_composition, subject_composition_unmasked;
  
  s_profile* p = NULL;
  int32_t maskLen = query_length / 2;

  const int8_t* query_sequence = (const int8_t*)(mask != COMPOSITIONAL_MASK_BOTH ? query_sequence_unmasked : query_sequence_masked);
  string file_db_str(temporary_dir);
  file_db_str +=  "/db.fa";
  std::ofstream file_db(file_db_str);
  string file_matrix_str(temporary_dir);
  file_matrix_str +=  "/matrix.txt";
  std::ofstream file_matrix(file_matrix_str);
  string stage2_res_str(temporary_dir);
  stage2_res_str +=  "/stage2.res";

#ifdef TIME_PAIRCOUNT
#ifdef COMPO_THRESHOLDS
  printf("stage2\n");
#endif //COMPO_THRESHOLDS
  auto begin2 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

  for(long i = 0; i<hits; i++)
  {
    bool symm_check = false;
    hits_gethit(i, & seqno, & score,
		  & qstrand, & qframe,
		  & dstrand, & dframe);

    if (score < minscore)
      continue;

    // stage 2

#ifdef TIME_PAIRCOUNT
    pairs2++;
    auto begin2 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    // if stage1-score >= skip_stage2_score, skip (bypass) stage2 for this hit and pass it on to stage3
    if (compo_adjust_mode == eNoCompositionBasedStats || score >= skip_stage2_score) {
#ifdef TIME_PAIRCOUNT
       // time stage 2
      auto end2 = std::chrono::high_resolution_clock::now();
      time_stage2 += (end2 - begin2).count();
      skipped_stage2++;
#endif //TIME_PAIRCOUNT
      adjusted_score = score * scaling_factor;
      if (adjusted_score >= minscore2 * scaling_factor) {
        stage3_indices.push_back(i);
      }
    } else {

    matchStart = queryStart = matchEnd = queryEnd = 0;

    // // OLD:
    // db_mapsequences(dbt, seqno, seqno);
    // db_getsequence(dbt, seqno, dstrand, dframe, &subject_sequence_masked, &subject_length, &ntlen, 0);
    // subject_length--;
    // if (mask)
    //   subject_sequence_unmasked = (char*)(*dbSequences)[seqno].data();
    // else
    //   subject_sequence_unmasked = subject_sequence_masked;

    // Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, (const Uint1*)subject_sequence_masked, subject_length);
    // if (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
    //   Blast_ReadAaComposition(&subject_composition_unmasked, BLASTAA_SIZE, (const Uint1*)subject_sequence_unmasked, subject_length);
    //   symm_check = (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) && (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
    // }

    // NEW:
    get_subject_composition(seqno, &subject_composition);
    get_subject_composition_unmasked(seqno, &subject_composition_unmasked);
    subject_sequence_masked = get_subject_sequence_masked(seqno);
    subject_sequence_unmasked = get_subject_sequence_unmasked(seqno);
    get_subject_length(seqno, &subject_length);

    if (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
        symm_check = (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) && (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
      }

#ifdef TIME_PAIRCOUNT
    // time init, aa composition
    auto end2b = std::chrono::high_resolution_clock::now();
    time_stage2b += (end2b - begin2).count();
    auto begin2c = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT
    
    matrix_adjust_rule = eDontAdjustMatrix;
    if (compo_adjust_mode != eNoCompositionBasedStats){
      compo_adjusted_matrix(NRrecord, sbp, scaledMatrixInfo,
        mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
        mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode, 2, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
    }

    // in the conditional matrix adj. mode, use stage1 score instead of the compo based statistics score if no_compo_stats=true
    if (matrix_adjust_rule == 0 && compo_adjust_mode == eCompositionMatrixAdjust && no_compo_stats == 1) {
       // skip stage 2 alignment
#ifdef TIME_PAIRCOUNT
      auto end2c = std::chrono::high_resolution_clock::now();
      time_stage2c += (end2c - begin2c).count();
      auto end2 = std::chrono::high_resolution_clock::now();
      time_stage2 += (end2 - begin2).count();
#endif //TIME_PAIRCOUNT
      adjusted_score = score * scaling_factor;
      if (adjusted_score >= minscore2 * scaling_factor) {
        stage3_indices.push_back(i);
      }
    } else {
    // convert adjusted score matrix to sswlib matrix
    for (int a = 0; a < BLASTAA_SIZE; a++)
      for (int b = 0; b < BLASTAA_SIZE; b++) {
#ifdef SWLIB_8BIT
        mat8[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
#else
        mata[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
#endif // SWLIB_8BIT
    }

#ifdef TIME_PAIRCOUNT
    // time 1st comp. adjustment
    auto end2c = std::chrono::high_resolution_clock::now();
    time_stage2c += (end2c - begin2c).count();
    auto begin_writing = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

//     // write matrix to file
//     file_matrix << "#" << i << "\n ";
//     for (int a = 0; a < (BLASTAA_SIZE); a++){
//       file_matrix << map_ncbi_int2aa[a] << " ";
//     }
//     file_matrix << "\n";
//     for (int a = 0; a < (BLASTAA_SIZE); a++){
//       file_matrix << map_ncbi_int2aa[a] << " ";
//       for (int b = 0; b < (BLASTAA_SIZE); b++){
// #ifdef SWLIB_8BIT
//         file_matrix << (int32_t)mat8[a*BLASTAA_SIZE + b] << " ";
// #else
//         file_matrix << (int32_t)mata[a*BLASTAA_SIZE + b] << " ";
// #endif // SWLIB_8BIT
//       }
//       file_matrix << "\n";
//     }

    // write matrix to file
    std::stringstream ss_matrix;
    ss_matrix << "#" << i << "\n ";
    for (int a = 0; a < (BLASTAA_SIZE); a++){
      ss_matrix << map_ncbi_int2aa[a] << " ";
    }
    ss_matrix << "\n";
    for (int a = 0; a < (BLASTAA_SIZE); a++){
      ss_matrix << map_ncbi_int2aa[a] << " ";
      for (int b = 0; b < (BLASTAA_SIZE); b++){
#ifdef SWLIB_8BIT
        ss_matrix << (int32_t)mat8[a*BLASTAA_SIZE + b] << " ";
#else
        ss_matrix << (int32_t)mata[a*BLASTAA_SIZE + b] << " ";
#endif // SWLIB_8BIT
      }
      ss_matrix << "\n";
    }
    file_matrix << ss_matrix.rdbuf();

#ifdef TIME_PAIRCOUNT
    // time file writing
    auto end_writing_matrix = std::chrono::high_resolution_clock::now();
    time_writing_matrix += (end_writing_matrix - begin_writing).count();
    auto begin_writing_subject = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    // write subject sequence to file
    aa_seq = int2aa(mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? subject_sequence_unmasked : subject_sequence_masked, subject_length);
    // use i as header 
    file_db << ">" << i << "\n" << aa_seq << "\n";

#ifdef TIME_PAIRCOUNT
    // time file writing
    auto end_writing = std::chrono::high_resolution_clock::now();
    time_writing += (end_writing - begin_writing).count();
    time_writing_subject += (end_writing - begin_writing_subject).count();
#endif //TIME_PAIRCOUNT

    }
    }
  } // Stage 2 comp. adjustement finished
  
  file_db.close();
  file_matrix.close();

#ifdef TIME_PAIRCOUNT
  auto begin_writing = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

  // write query_sequence to file
  aa_seq = int2aa((char*)query_sequence, query_length);
  string file_query_str(temporary_dir);
  file_query_str +=  "/query.fa";
  std::ofstream file_query(file_query_str);
  file_query << ">" << query.description << "\n" << aa_seq << "\n";
  file_query.close();

#ifdef TIME_PAIRCOUNT
  // time file writing
  auto end_writing = std::chrono::high_resolution_clock::now();
  time_writing += (end_writing - begin_writing).count();
  time_writing_query += (end_writing - begin_writing).count();
  auto begin2d = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

  string system_str = ssw_binary_dir;
#ifdef SWLIB_8BIT
  system_str +=  "/ssw_main_sw_sse2_byte_multimatrix";
#else
  system_str +=  "/ssw_main_sw_sse2_word_multimatrix";
#endif // SWLIB_8BIT
  system_str += " -o ";
  system_str += std::to_string(gap_open_extend);
  system_str += " -e ";
  system_str += std::to_string(gap_extend);
  system_str += " -f ";
  system_str += std::to_string(minScoreSWlib);
  system_str += " -p ";
  system_str += file_db_str;
  system_str += " ";
  system_str += file_query_str;
  system_str += " ";
  system_str += file_matrix_str;
  system_str += " > ";
  system_str +=  stage2_res_str;
  system(system_str.c_str());

#ifdef TIME_PAIRCOUNT
  // time 1st alignment
  auto end2d = std::chrono::high_resolution_clock::now();
  time_stage2d += (end2d - begin2d).count();
  // time stage 2
  auto end2 = std::chrono::high_resolution_clock::now();
  time_stage2 += (end2 - begin2).count();
  auto begin3 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

  // stage 3
  // parse stage2 results to find i´s for stage3
  char line[256];
  FILE *stage2_fp = fopen(stage2_res_str.c_str(), "r");
  while (fgets(line, 256, stage2_fp)) {
    string str(line);
    if (str.find("tar") != std::string::npos) {
      // parse i
      stringstream ss;  
      ss << str.substr(13, strlen(line));
      long i;
      ss >> i;
      // store i in array
      stage3_indices.push_back(i);
    }
  }

  // loop over stage3_indices
  for(const auto& i: stage3_indices) {

#ifdef TIME_PAIRCOUNT
#ifdef COMPO_THRESHOLDS
    printf("stage3\n");
#endif //COMPO_THRESHOLDS
    pairs3++;
#endif //TIME_PAIRCOUNT

    long flags = 0;
    bool symm_check = false;
    long adjusted_score_blast = 0, adjusted_score_blast_rev = 0;
    Blast_AminoAcidComposition subject_composition, subject_composition_unmasked;
    
    // get seqno
    hits_gethit(i, & seqno, & score,
    & qstrand, & qframe,
    & dstrand, & dframe);
    // get subject sequence/composition/length
    get_subject_composition(seqno, &subject_composition);
    get_subject_composition_unmasked(seqno, &subject_composition_unmasked);
    subject_sequence_masked = get_subject_sequence_masked(seqno);
    subject_sequence_unmasked = get_subject_sequence_unmasked(seqno);
    get_subject_length(seqno, &subject_length);

    if (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
      symm_check = (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) && (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
    }
    matrix_adjust_rule = eDontAdjustMatrix;
    if (compo_adjust_mode_stage3 != eNoCompositionBasedStats){
      compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, 
        mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
        mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode_stage3, 3, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
    }

    // convert adjusted score matrix to swipe matix
    for (int a = 0; a < BLASTAA_SIZE; a++)
      for (int b = 0; b < BLASTAA_SIZE; b++) {
        mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
        adjusted_score_matrix_63[(a<<5) + b] = sbpBL62->matrix->data[a][b];
        adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
      }
    p = ssw_init_word((const int8_t*)query_sequence_unmasked, query_length, mata, BLASTAA_SIZE);
    sw_sse2_word2((const int8_t*)subject_sequence_unmasked, 0, subject_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);

    if (p)
      init_destroy(p);
    result->score1 = ae.score;
    result->ref_end1 = ae.ref;
    result->read_end1 = ae.read;
    if (result->score1 == 32767) {
      flags |= HIT_LARGE_SCORE;
    }
    
    matchStart = queryStart = 0;
    qseq = query_sequence_unmasked;
    qlen = query_length;
    dseq = subject_sequence_unmasked;
    dlen = subject_length;
    hits_enter_seq(i, dseq, dlen);
    
    if (symm_check) {

#ifdef TIME_PAIRCOUNT
      pairs3b++;
#endif //TIME_PAIRCOUNT

      if (result->score1 == 32767) {
        if (hesize < dlen * 32) {
          hesize = dlen*32;
          hearray = (BYTE*)realloc(hearray, hesize);
        }
        //TODO: replace fullsw by an int version of swlib
        result->score1 = fullsw(dseq, dseq + dlen, qseq, qseq + qlen, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
      }
                
      matrix_adjust_rule = eDontAdjustMatrix;
      compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode_stage3, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
      // convert adjusted score matrix to swipe matix
      for (int a = 0; a < BLASTAA_SIZE; a++)
        for (int b = 0; b < BLASTAA_SIZE; b++) {
          mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
          adjusted_score_matrix_63_rev[(a<<5) + b] = sbpBL62->matrix->data[a][b];
          adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
        }
      p = ssw_init_word((const int8_t*)subject_sequence_unmasked, subject_length, mata, BLASTAA_SIZE);
      sw_sse2_word2((const int8_t*)query_sequence_unmasked, 0, query_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);
      if (p)
        init_destroy(p);
      result->score2 = ae.score;
      result->ref_end2 = ae.ref;
      result->read_end2 = ae.read;
      if (result->score2 == 32767) {
        flags |= HIT_LARGE_SCORE;
        char *qseq2 = subject_sequence_unmasked;
        long qlen2 = subject_length;
        char *dseq2 = query_sequence_unmasked;
        long dlen2 = query_length;
        if (hesize < dlen2 * 32) {
          hesize = dlen2*32;
          hearray = (BYTE*)realloc(hearray, hesize);
        }
        //TODO: replace fullsw by an int version of swlib
        result->score2 = fullsw(dseq2, dseq2 + dlen2, qseq2, qseq2 + qlen2, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
      }
//        if (round(result->score2/scaling_factor_BL62) > round(result->score1/scaling_factor_BL62)) {
//        if (round((float)result->score2/scaling_factor_BL62) > round((float)result->score1/scaling_factor_BL62)) {
      if (result->score2 > result->score1) {

#ifdef TIME_PAIRCOUNT
        pairs3c++;
#endif //TIME_PAIRCOUNT

        flags |= HIT_SUBJECT_QUERY_BEST_BL62;
        //fprintf(stderr, "symmetrical score is larger than adjusted blastp score: %d vs %d\n", (int)round(result->score2/scaling_factor_BL62), (int)round(result->score1/scaling_factor_BL62));
      }
      adjusted_score_blast_rev = result->score2;
    }
    
    int64_t *used_score_matrix_63 = adjusted_score_matrix_63;
    
    if (flags & HIT_LARGE_SCORE)
      matchEnd = queryEnd = score_align = 0;
    else {
      if (flags & HIT_SUBJECT_QUERY_BEST_BL62) {
        adjusted_score_blast = result->score2;
        matchEnd = result->ref_end2;
        queryEnd = result->read_end2;
        adjusted_score_blast_rev = result->score1;
        used_score_matrix_63 = adjusted_score_matrix_63_rev;
        qseq = subject_sequence_unmasked;
        qlen = subject_length;
        dseq = query_sequence_unmasked;
        dlen = query_length;
      } else {
        adjusted_score_blast = result->score1;
        matchEnd = result->ref_end1;
        queryEnd = result->read_end1;
      }
      score_align = adjusted_score_blast;
    }
    //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
/*      fprintf(out, "query sequence: \n");
    db_print_seq_map(qseq, qlen, sym_ncbi_aa);
    fprintf(out, "subject sequence: \n");
    db_print_seq_map(dseq, dlen, sym_ncbi_aa);
*/      align(qseq,
      dseq,
      qlen,
      dlen,
      used_score_matrix_63,
      gapopen_BlastDef,
      gapextend_BlastDef,
      & queryStart,
      & matchStart,
      & queryEnd,
      & matchEnd,
      & alignment,
      & score_align);
      adjusted_score_blast = score_align;
//          fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
/*        } else {
      score_align = adjusted_score_blast;
      // reverse query sequence; queryEnd is 0-based alignment ending position
      int32_t start = 0, end = queryEnd;
      if (end > reverse_len) {
        if (reverse_len == 0)
          read_reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));
        else
          read_reverse = (int8_t*)realloc(read_reverse, (end + 1) * sizeof(int8_t));
        reverse_len = end;
      }
      while (start <= end) {
    read_reverse[start] = qseq[end];
        read_reverse[end] = qseq[start];
        ++ start;
        -- end;
      }
      __m128i* vP = qP_word16(read_reverse, mata, queryEnd + 1, p->n);
      alignment_end* bests_reverse = sw_sse2_word2((const int8_t*)subject_sequence, 1, matchEnd + 1, queryEnd + 1, gapopen_BlastDef + gapextend_BlastDef, gapextend_BlastDef, vP, score_align, maskLen);
      free(vP);
      matchStart = bests_reverse[0].ref;
      queryStart = queryEnd - bests_reverse[0].read;
      result->score1 = bests_reverse[0].score;
      free(bests_reverse);
      matchStart_r = matchEnd_r = queryStart_r = queryEnd_r = score_align_r = 0;
    align(&qseq[queryStart],
      &subject_sequence[matchStart],
      queryEnd - queryStart + 1,
      matchEnd - matchStart + 1,
      adjusted_score_matrix_63,
      gapopen_BlastDef,
      gapextend_BlastDef,
      & queryStart_r,
      & matchStart_r,
      & queryEnd_r,
      & matchEnd_r,
      & alignment,
      & score_align_r);
      //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", queryStart, queryEnd, matchStart, matchEnd, score_align_r, result->score1);
    }*/
    hits_set_align_string(i, alignment, score_align);
    hits_enter_align_coord(i, queryStart, queryEnd, matchStart, matchEnd, 0);
//        score_matrix_63 = used_score_matrix_63;
//        score_matrix_dump();
    count_align_matrix(i, used_score_matrix_63, qseq, qlen, dseq, dlen);
//        show_align(i);

    hits_enter_adjusted_score(i, adjusted_score, adjusted_score_blast, adjusted_score_blast_rev, flags);
  } 

#ifdef TIME_PAIRCOUNT
  // time stage 3
  auto end3 = std::chrono::high_resolution_clock::now();
  time_stage3 += (end3 - begin3).count();
#endif //TIME_PAIRCOUNT

  if (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
    std::ofstream file_db(file_db_str);
    std::ofstream file_matrix(file_matrix_str);
    for(long i = 0; i<hits; i++) {
      // skip if i was already processed (result->score1 >= minscore2 * scaling_factor)
      if(std::find(stage3_indices.begin(), stage3_indices.end(), i) != stage3_indices.end()) {
        continue;
      }

#ifdef TIME_PAIRCOUNT
      auto begin2 = std::chrono::high_resolution_clock::now();
      auto begin2e = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

      hits_gethit(i, & seqno, & score,
        & qstrand, & qframe,
        & dstrand, & dframe);
      
      get_subject_composition(seqno, &subject_composition);
      get_subject_composition_unmasked(seqno, &subject_composition_unmasked);
      subject_sequence_masked = get_subject_sequence_masked(seqno);
      subject_sequence_unmasked = get_subject_sequence_unmasked(seqno);
      get_subject_length(seqno, &subject_length);
      
      // TODO: loop only over i with valid symm_check + not already processed
      bool symm_check = (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
      if (!symm_check) {
        continue;
      }

#ifdef TIME_PAIRCOUNT
      pairs2b++;
#endif //TIME_PAIRCOUNT

      matrix_adjust_rule = eDontAdjustMatrix;
      compo_adjusted_matrix(NRrecord, sbp, scaledMatrixInfo, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
      // convert adjusted score matrix to swipe matix
      for (int a = 0; a < BLASTAA_SIZE; a++)
        for (int b = 0; b < BLASTAA_SIZE; b++) {
#ifdef SWLIB_8BIT
          mat8[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
#else
          mata[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
#endif // SWLIB_8BIT
        }

      // write matrix to file
      file_matrix << "#" << i << "\n ";
      for (int a = 0; a < (BLASTAA_SIZE); a++){
        file_matrix << map_ncbi_int2aa[a] << " ";
      }
      file_matrix << "\n";
      for (int a = 0; a < (BLASTAA_SIZE); a++){
        file_matrix << map_ncbi_int2aa[a] << " ";
        for (int b = 0; b < (BLASTAA_SIZE); b++){
  #ifdef SWLIB_8BIT
          file_matrix << (int32_t)mat8[a*BLASTAA_SIZE + b] << " ";
  #else
          file_matrix << (int32_t)mata[a*BLASTAA_SIZE + b] << " ";
  #endif // SWLIB_8BIT
        }
        file_matrix << "\n";
      }

      // write subject sequence to file
      aa_seq = int2aa(subject_sequence_unmasked, subject_length);
      // use i as header 
      file_db << ">" << i << "\n" << aa_seq << "\n";

#ifdef TIME_PAIRCOUNT
      // time reverse comp. adjustment
      auto end2e = std::chrono::high_resolution_clock::now();
      time_stage2e += (end2e - begin2e).count();
#endif //TIME_PAIRCOUNT

    }
    file_db.close();
    file_matrix.close();

    // write query_sequence to file
    const int8_t* query_sequence = (const int8_t*)(mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? query_sequence_unmasked : query_sequence_masked);
    aa_seq = int2aa((char*)query_sequence, query_length);
    string file_query_str(temporary_dir);
    file_query_str +=  "/query.fa";
    std::ofstream file_query(file_query_str);
    file_query << ">" << query.description << "\n" << aa_seq << "\n";
    file_query.close();

#ifdef TIME_PAIRCOUNT
    auto begin2f = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    string system_str = ssw_binary_dir;
  #ifdef SWLIB_8BIT
    system_str +=  "/ssw_main_sw_sse2_byte_multimatrix";
  #else
    system_str +=  "/ssw_main_sw_sse2_word_multimatrix";
  #endif // SWLIB_8BIT
    system_str += " -o ";
    system_str += std::to_string(gap_open_extend);
    system_str += " -e ";
    system_str += std::to_string(gap_extend);
    system_str += " -f ";
    system_str += std::to_string(minScoreSWlib);
    system_str += " -p ";
    system_str += file_query_str;
    system_str += " ";
    system_str += file_db_str;
    system_str += " ";
    system_str += file_matrix_str;
    system_str += " > ";
    system_str +=  stage2_res_str;
    system(system_str.c_str());

#ifdef TIME_PAIRCOUNT
    // time reverse alignment
    auto end2f = std::chrono::high_resolution_clock::now();
    time_stage2f += (end2f - begin2f).count();
    // time stage 2
    auto end2 = std::chrono::high_resolution_clock::now();
    time_stage2 += (end2 - begin2).count();
    auto begin3 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    // stage 3
    char line[256];
    FILE *stage2_fp = fopen(stage2_res_str.c_str(), "r");
    while (fgets(line, 256, stage2_fp)) {
      string str(line);
      if (str.find("que") != std::string::npos) {

#ifdef TIME_PAIRCOUNT
#ifdef COMPO_THRESHOLDS
        printf("stage3\n");
#endif //COMPO_THRESHOLDS
        pairs3_rev++;
        pairs3++;
#endif //TIME_PAIRCOUNT

        long flags = 0;
        flags |= HIT_SUBJECT_QUERY_BEST_BL50;
        bool symm_check = false;
        long adjusted_score_blast = 0, adjusted_score_blast_rev = 0;
        Blast_AminoAcidComposition subject_composition, subject_composition_unmasked;

        // parse i
        stringstream ss;  
        ss << str.substr(12, strlen(line));
        long i;
        ss >> i;
        // get seqno
        hits_gethit(i, & seqno, & score,
        & qstrand, & qframe,
        & dstrand, & dframe);
        // get subject sequence/composition/length
        get_subject_composition(seqno, &subject_composition);
        get_subject_composition_unmasked(seqno, &subject_composition_unmasked);
        subject_sequence_masked = get_subject_sequence_masked(seqno);
        subject_sequence_unmasked = get_subject_sequence_unmasked(seqno);
        get_subject_length(seqno, &subject_length);

        if (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
          symm_check = (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) && (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
        }

        matrix_adjust_rule = eDontAdjustMatrix;
        if (compo_adjust_mode_stage3 != eNoCompositionBasedStats){
          compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, 
            mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
            mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode_stage3, 3, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
        }

        // convert adjusted score matrix to swipe matix
        for (int a = 0; a < BLASTAA_SIZE; a++)
          for (int b = 0; b < BLASTAA_SIZE; b++) {
            mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
            adjusted_score_matrix_63[(a<<5) + b] = sbpBL62->matrix->data[a][b];
            adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
          }

        p = ssw_init_word((const int8_t*)query_sequence_unmasked, query_length, mata, BLASTAA_SIZE);
        sw_sse2_word2((const int8_t*)subject_sequence_unmasked, 0, subject_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);

        if (p)
          init_destroy(p);
        result->score1 = ae.score;
        result->ref_end1 = ae.ref;
        result->read_end1 = ae.read;
        if (result->score1 == 32767) {
          flags |= HIT_LARGE_SCORE;
        }
        
        matchStart = queryStart = 0;
        qseq = query_sequence_unmasked;
        qlen = query_length;
        dseq = subject_sequence_unmasked;
        dlen = subject_length;
        hits_enter_seq(i, dseq, dlen);
        
        if (symm_check) {

#ifdef TIME_PAIRCOUNT
          pairs3b++;
#endif //TIME_PAIRCOUNT

          if (result->score1 == 32767) {
            if (hesize < dlen * 32) {
              hesize = dlen*32;
              hearray = (BYTE*)realloc(hearray, hesize);
            }
            //TODO: replace fullsw by an int version of swlib
            result->score1 = fullsw(dseq, dseq + dlen, qseq, qseq + qlen, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
          }
                    
          matrix_adjust_rule = eDontAdjustMatrix;
          compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode_stage3, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
          // convert adjusted score matrix to swipe matix
          for (int a = 0; a < BLASTAA_SIZE; a++)
            for (int b = 0; b < BLASTAA_SIZE; b++) {
              mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
              adjusted_score_matrix_63_rev[(a<<5) + b] = sbpBL62->matrix->data[a][b];
              adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
            }
          p = ssw_init_word((const int8_t*)subject_sequence_unmasked, subject_length, mata, BLASTAA_SIZE);
          sw_sse2_word2((const int8_t*)query_sequence_unmasked, 0, query_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);
          if (p)
            init_destroy(p);
          result->score2 = ae.score;
          result->ref_end2 = ae.ref;
          result->read_end2 = ae.read;
          if (result->score2 == 32767) {
            flags |= HIT_LARGE_SCORE;
            char *qseq2 = subject_sequence_unmasked;
            long qlen2 = subject_length;
            char *dseq2 = query_sequence_unmasked;
            long dlen2 = query_length;
            if (hesize < dlen2 * 32) {
              hesize = dlen2*32;
              hearray = (BYTE*)realloc(hearray, hesize);
            }
            //TODO: replace fullsw by an int version of swlib
            result->score2 = fullsw(dseq2, dseq2 + dlen2, qseq2, qseq2 + qlen2, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
          }
  //        if (round(result->score2/scaling_factor_BL62) > round(result->score1/scaling_factor_BL62)) {
  //        if (round((float)result->score2/scaling_factor_BL62) > round((float)result->score1/scaling_factor_BL62)) {
          if (result->score2 > result->score1) {

#ifdef TIME_PAIRCOUNT
            pairs3c++;
#endif //TIME_PAIRCOUNT

            flags |= HIT_SUBJECT_QUERY_BEST_BL62;
            //fprintf(stderr, "symmetrical score is larger than adjusted blastp score: %d vs %d\n", (int)round(result->score2/scaling_factor_BL62), (int)round(result->score1/scaling_factor_BL62));
          }
          adjusted_score_blast_rev = result->score2;
        }
        
        int64_t *used_score_matrix_63 = adjusted_score_matrix_63;
        
        if (flags & HIT_LARGE_SCORE)
          matchEnd = queryEnd = score_align = 0;
        else {
          if (flags & HIT_SUBJECT_QUERY_BEST_BL62) {
            adjusted_score_blast = result->score2;
            matchEnd = result->ref_end2;
            queryEnd = result->read_end2;
            adjusted_score_blast_rev = result->score1;
            used_score_matrix_63 = adjusted_score_matrix_63_rev;
            qseq = subject_sequence_unmasked;
            qlen = subject_length;
            dseq = query_sequence_unmasked;
            dlen = query_length;
          } else {
            adjusted_score_blast = result->score1;
            matchEnd = result->ref_end1;
            queryEnd = result->read_end1;
          }
          score_align = adjusted_score_blast;
        }
        //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
  /*      fprintf(out, "query sequence: \n");
        db_print_seq_map(qseq, qlen, sym_ncbi_aa);
        fprintf(out, "subject sequence: \n");
        db_print_seq_map(dseq, dlen, sym_ncbi_aa);
  */      align(qseq,
          dseq,
          qlen,
          dlen,
          used_score_matrix_63,
          gapopen_BlastDef,
          gapextend_BlastDef,
          & queryStart,
          & matchStart,
          & queryEnd,
          & matchEnd,
          & alignment,
          & score_align);
          adjusted_score_blast = score_align;
  //          fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
  /*        } else {
          score_align = adjusted_score_blast;
          // reverse query sequence; queryEnd is 0-based alignment ending position
          int32_t start = 0, end = queryEnd;
          if (end > reverse_len) {
            if (reverse_len == 0)
              read_reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));
            else
              read_reverse = (int8_t*)realloc(read_reverse, (end + 1) * sizeof(int8_t));
            reverse_len = end;
          }
          while (start <= end) {
        read_reverse[start] = qseq[end];
            read_reverse[end] = qseq[start];
            ++ start;
            -- end;
          }
          __m128i* vP = qP_word16(read_reverse, mata, queryEnd + 1, p->n);
          alignment_end* bests_reverse = sw_sse2_word2((const int8_t*)subject_sequence, 1, matchEnd + 1, queryEnd + 1, gapopen_BlastDef + gapextend_BlastDef, gapextend_BlastDef, vP, score_align, maskLen);
          free(vP);
          matchStart = bests_reverse[0].ref;
          queryStart = queryEnd - bests_reverse[0].read;
          result->score1 = bests_reverse[0].score;
          free(bests_reverse);
          matchStart_r = matchEnd_r = queryStart_r = queryEnd_r = score_align_r = 0;
        align(&qseq[queryStart],
          &subject_sequence[matchStart],
          queryEnd - queryStart + 1,
          matchEnd - matchStart + 1,
          adjusted_score_matrix_63,
          gapopen_BlastDef,
          gapextend_BlastDef,
          & queryStart_r,
          & matchStart_r,
          & queryEnd_r,
          & matchEnd_r,
          & alignment,
          & score_align_r);
          //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", queryStart, queryEnd, matchStart, matchEnd, score_align_r, result->score1);
        }*/
        hits_set_align_string(i, alignment, score_align);
        hits_enter_align_coord(i, queryStart, queryEnd, matchStart, matchEnd, 0);
  //        score_matrix_63 = used_score_matrix_63;
  //        score_matrix_dump();
        count_align_matrix(i, used_score_matrix_63, qseq, qlen, dseq, dlen);
  //        show_align(i);

        hits_enter_adjusted_score(i, adjusted_score, adjusted_score_blast, adjusted_score_blast_rev, flags);
      } 
    }

#ifdef TIME_PAIRCOUNT
    // time stage 3
    auto end3 = std::chrono::high_resolution_clock::now();
    time_stage3 += (end3 - begin3).count();
#endif //TIME_PAIRCOUNT

    }
  db_thread_destruct(dbt);
}

// // trying to speed up align_adjusted (without using SSW executable):
// void align_adjusted5() {
//   long hits = hits_getcount();

//   struct db_thread_s * dbt = db_thread_create();

//   char *subject_sequence_masked, *subject_sequence_unmasked;
//   long subject_length;
//   char *query_sequence_masked = query.aa[0].seq;
//   char *query_sequence_unmasked = (mask ? query.aa[0].seq_unmasked : query.aa[0].seq);
//   long query_length = query.aa[0].len;
//   char *qseq, *dseq;
//   long qlen, dlen;
//   long ntlen;
//   long seqno;
//   long score;
//   long qstrand;
//   long qframe;
//   long dstrand;
//   long dframe;
//   long adjusted_score;
//   long matchStart, queryStart, matchEnd, queryEnd, score_align;
//   //long matchStart_r, queryStart_r, matchEnd_r, queryEnd_r, score_align_r;
//   char * alignment;
//   long gap_open = gapopen * scaling_factor;
//   long gap_extend = gapextend * scaling_factor;
//   long gap_open_extend = gap_open + gap_extend;
//   long gapopen_BlastDef = 11 * scaling_factor_BL62;
//   long gapextend_BlastDef = 1 * scaling_factor_BL62;
//   long gapopenextend_BlastDef = gapopen_BlastDef + gapextend_BlastDef;
//   EMatrixAdjustRule matrix_adjust_rule;
//   ECompoAdjustModes compo_adjust_mode_stage3 = eCompositionMatrixAdjust;
//   int16_t* mat;
//   char dseq2;
//   s_profile* p = NULL;
//   int32_t maskLen = query_length / 2;
  
//   for(long i = 0; i<hits; i++)
//   {
//     long flags = 0;
//     bool symm_check = false;
//     long adjusted_score_blast = 0, adjusted_score_blast_rev = 0;
//     hits_gethit(i, & seqno, & score,
// 		  & qstrand, & qframe,
// 		  & dstrand, & dframe);

//     if (score < minscore)
//       continue;

//     // stage 2
// #ifdef TIME_PAIRCOUNT
// #ifdef COMPO_THRESHOLDS
//     printf("stage2\n");
// #endif //COMPO_THRESHOLDS
//     pairs2++;
//     auto begin2 = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT

//     matchStart = queryStart = matchEnd = queryEnd = 0;
//     Blast_AminoAcidComposition subject_composition, subject_composition_unmasked;

//     // OLD:
//     db_mapsequences(dbt, seqno, seqno);
//     db_getsequence(dbt, seqno, dstrand, dframe, &subject_sequence_masked, &subject_length, &ntlen, 0);
//     subject_length--;
//     if (mask)
//       subject_sequence_unmasked = (char*)(*dbSequences)[seqno].data();
//     else
//       subject_sequence_unmasked = subject_sequence_masked;

//     Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, (const Uint1*)subject_sequence_masked, subject_length);
//     if (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) {
//       Blast_ReadAaComposition(&subject_composition_unmasked, BLASTAA_SIZE, (const Uint1*)subject_sequence_unmasked, subject_length);
//       symm_check = (mask == COMPOSITIONAL_MASK_SYMM || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM) && (subject_composition.numTrueAminoAcids < subject_length || query.composition.numTrueAminoAcids < query_length);
//     }

//     // // NEW:
//     // get_subject_composition(seqno, &subject_composition);
//     // get_subject_composition_unmasked(seqno, &subject_composition_unmasked);
//     // subject_sequence_masked = get_subject_sequence_masked(seqno);
//     // get_subject_length(seqno, &subject_length);
//     // if (mask)
//     //   subject_sequence_unmasked = (char*)(*dbSequences)[seqno].data();
//     // else
//     //   subject_sequence_unmasked = subject_sequence_masked;

// #ifdef TIME_PAIRCOUNT
//     // time init, aa composition
//     auto end2b = std::chrono::high_resolution_clock::now();
//     time_stage2b += (end2b - begin2).count();
// #endif //TIME_PAIRCOUNT

//     if (compo_adjust_mode == eNoCompositionBasedStats || score >= skip_stage2_score) {

// #ifdef TIME_PAIRCOUNT
//        // time stage 2
//       auto end2 = std::chrono::high_resolution_clock::now();
//       time_stage2 += (end2 - begin2).count();
//       skipped_stage2++;
// #endif //TIME_PAIRCOUNT

//       adjusted_score = score * scaling_factor;
//     } else {

// #ifdef TIME_PAIRCOUNT
//     auto begin2c = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT
    
//     matrix_adjust_rule = eDontAdjustMatrix;  // TODO REMOVE?
//     if (compo_adjust_mode != eNoCompositionBasedStats){
//       compo_adjusted_matrix(NRrecord, sbp, scaledMatrixInfo,
//         mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
//         mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode, 2, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
//     }

//     if (matrix_adjust_rule == 0 && compo_adjust_mode == eCompositionMatrixAdjust && no_compo_stats == 1) {
//        // skip stage 2 alignment

// #ifdef TIME_PAIRCOUNT
//       auto end2c = std::chrono::high_resolution_clock::now();
//       time_stage2c += (end2c - begin2c).count();
//       auto end2 = std::chrono::high_resolution_clock::now();
//       time_stage2 += (end2 - begin2).count();
// #endif //TIME_PAIRCOUNT

//       adjusted_score = score * scaling_factor;
//     } else {

//     // convert adjusted score matrix to sswlib matrix
//     for (int a = 0; a < BLASTAA_SIZE; a++)
//       for (int b = 0; b < BLASTAA_SIZE; b++) {
// #ifdef SWLIB_8BIT
//         mat8[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
//         // std::cout << sbp->matrix->data[b][a] << " ";
// #else
//         mata[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
//         // std::cout << sbp->matrix->data[b][a] << " ";
// #endif // SWLIB_8BIT
//     }

// #ifdef TIME_PAIRCOUNT
//     // time 1st comp. adjustment
//     auto end2c = std::chrono::high_resolution_clock::now();
//     time_stage2c += (end2c - begin2c).count();
//     auto begin2d = std::chrono::high_resolution_clock::now();
//     auto begin2_init = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT

// #ifdef SWLIB_8BIT
//     p = ssw_init((const int8_t*)(mask != COMPOSITIONAL_MASK_BOTH ? query_sequence_unmasked : query_sequence_masked), query_length, mat8, BLASTAA_SIZE, 0);

// #ifdef TIME_PAIRCOUNT
//     // time ssw_init_word
//     auto end2_init = std::chrono::high_resolution_clock::now();
//     time_stage2_init += (end2_init - begin2_init).count();
//     auto begin2_sw = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT
    
//     sw_sse2_byte2((const int8_t*)(mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? subject_sequence_unmasked : subject_sequence_masked), 0, subject_length, p->readLen, gap_open_extend, gap_extend, p->profile_byte, minScoreSWlib, p->bias, maskLen, &ae);

// #ifdef TIME_PAIRCOUNT
//     // time sw
//     auto end2_sw = std::chrono::high_resolution_clock::now();
//     time_stage2_sw += (end2_sw - begin2_sw).count();
// #endif //TIME_PAIRCOUNT
    
//     result->score1 = ae.score;
//     result->ref_end1 = ae.ref;
//     result->read_end1 = ae.read;
// #else
//     p = ssw_init_word((const int8_t*)(mask != COMPOSITIONAL_MASK_BOTH ? query_sequence_unmasked : query_sequence_masked), query_length, mata, BLASTAA_SIZE);

// #ifdef TIME_PAIRCOUNT
//     // time ssw_init_word
//     auto end2_init = std::chrono::high_resolution_clock::now();
//     time_stage2_init += (end2_init - begin2_init).count();
//     auto begin2_sw = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT

//     int32_t m, s1 = 67108864;
//     int8_t* ref_num = (int8_t*)malloc(s1);
//     while (subject_length > s1) {
//       ++s1;
//       kroundup32(s1);
//       ref_num = (int8_t*)realloc(ref_num, s1);
//     }
//     for (m = 0; m < subject_length; ++m) ref_num[m] = (mask == COMPOSITIONAL_MASK_NONE || mask == COMPOSITIONAL_MASK_BOTH_MATRIXONLY || mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? subject_sequence_unmasked : subject_sequence_masked)[m];

//     sw_sse2_word2(ref_num, 0, subject_length, p->readLen, gap_open_extend, gap_extend, p->profile_word, minScoreSWlib, maskLen, &ae);

// #ifdef TIME_PAIRCOUNT    
//     // time sw
//     auto end2_sw = std::chrono::high_resolution_clock::now();
//     time_stage2_sw += (end2_sw - begin2_sw).count();
// #endif //TIME_PAIRCOUNT

//     result->score1 = ae.score;
//     result->ref_end1 = ae.ref;
//     result->read_end1 = ae.read;
// #endif // SWLIB_8BIT
//     if (p)
//        init_destroy(p);

// #ifdef TIME_PAIRCOUNT
//     // time 1st alignment
//     auto end2d = std::chrono::high_resolution_clock::now();
//     time_stage2d += (end2d - begin2d).count();
//     auto begin2e = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT
    
//     if (symm_check && result->score1 < minscore2 * scaling_factor) {

// #ifdef TIME_PAIRCOUNT
//       pairs2b++;
// #endif //TIME_PAIRCOUNT

//       matrix_adjust_rule = eDontAdjustMatrix;
//       compo_adjusted_matrix(NRrecord, sbp, scaledMatrixInfo, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
//       // convert adjusted score matrix to swipe matix
//       for (int a = 0; a < BLASTAA_SIZE; a++)
//         for (int b = 0; b < BLASTAA_SIZE; b++) {
// #ifdef SWLIB_8BIT
//           mat8[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
// #else
//           mata[a * BLASTAA_SIZE + b] = sbp->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
// #endif // SWLIB_8BIT
//         }

// #ifdef TIME_PAIRCOUNT
//       // time reverse comp. adjustment
//       auto end2e = std::chrono::high_resolution_clock::now();
//       time_stage2e += (end2e - begin2e).count();
//       auto begin2f = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT

// #ifdef SWLIB_8BIT
//       p = ssw_init((const int8_t*)subject_sequence_unmasked, subject_length, mat8, BLASTAA_SIZE, 0);
//       sw_sse2_byte2((const int8_t*)(mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? query_sequence_unmasked : query_sequence_masked), 0, query_length, p->readLen, gap_open_extend, gap_extend, p->profile_byte, minScoreSWlib, p->bias, maskLen, &ae);
// #else
//       p = ssw_init_word((const int8_t*)subject_sequence_unmasked, subject_length, mata, BLASTAA_SIZE);
//       sw_sse2_word2((const int8_t*)(mask == COMPOSITIONAL_MASK_MATRIXONLY_SYMM ? query_sequence_unmasked : query_sequence_masked), 0, query_length, p->readLen, gap_open_extend, gap_extend, p->profile_word, minScoreSWlib, maskLen, &ae);
// #endif // SWLIB_8BIT
//       if (p)
//         init_destroy(p);
//       if (ae.score > result->score1) {

// #ifdef TIME_PAIRCOUNT
//         pairs2c++;
// #endif //TIME_PAIRCOUNT

//         result->score1 = ae.score;
//         flags |= HIT_SUBJECT_QUERY_BEST_BL50;
//       }

// #ifdef TIME_PAIRCOUNT
//       // time reverse alignment
//       auto end2f = std::chrono::high_resolution_clock::now();
//       time_stage2f += (end2e - begin2e).count();
// #endif //TIME_PAIRCOUNT

//     }
//     adjusted_score = result->score1;

// #ifdef TIME_PAIRCOUNT
//     // time stage 2
//     auto end2 = std::chrono::high_resolution_clock::now();
//     time_stage2 += (end2 - begin2).count();
// #endif //TIME_PAIRCOUNT

//     }
//     }
    
//     if (adjusted_score >= minscore2 * scaling_factor) {
//       // stage 3
// #ifdef TIME_PAIRCOUNT
// #ifdef COMPO_THRESHOLDS
//       printf("stage3\n");
// #endif //COMPO_THRESHOLDS
//       pairs3++;
//       auto begin3 = std::chrono::high_resolution_clock::now();
// #endif //TIME_PAIRCOUNT
      
//       matrix_adjust_rule = eDontAdjustMatrix;
//       if (compo_adjust_mode_stage3 != eNoCompositionBasedStats){
//         compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, 
//           mask != COMPOSITIONAL_MASK_BOTH && mask != COMPOSITIONAL_MASK_BOTH_MATRIXONLY ? &query.composition_unmasked : &query.composition, query_length,
//           mask == COMPOSITIONAL_MASK_NONE ? &subject_composition_unmasked : &subject_composition, subject_length, compo_adjust_mode_stage3, 3, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
//       }

//       // convert adjusted score matrix to swipe matix
//       for (int a = 0; a < BLASTAA_SIZE; a++)
//         for (int b = 0; b < BLASTAA_SIZE; b++) {
//           mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
//           adjusted_score_matrix_63[(a<<5) + b] = sbpBL62->matrix->data[a][b];
//           adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
//         }
//       p = ssw_init_word((const int8_t*)query_sequence_unmasked, query_length, mata, BLASTAA_SIZE);
//       sw_sse2_word2((const int8_t*)subject_sequence_unmasked, 0, subject_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);
    
//       if (p)
//         init_destroy(p);
//       result->score1 = ae.score;
//       result->ref_end1 = ae.ref;
//       result->read_end1 = ae.read;
//       if (result->score1 == 32767) {
//         flags |= HIT_LARGE_SCORE;
//       }
      
//       matchStart = queryStart = 0;
//       qseq = query_sequence_unmasked;
//       qlen = query_length;
//       dseq = subject_sequence_unmasked;
//       dlen = subject_length;
//       hits_enter_seq(i, dseq, dlen);
      
//       if (symm_check) {

// #ifdef TIME_PAIRCOUNT
//         pairs3b++;
// #endif //TIME_PAIRCOUNT

//         if (result->score1 == 32767) {
//           if (hesize < dlen * 32) {
//             hesize = dlen*32;
//             hearray = (BYTE*)realloc(hearray, hesize);
//           }
//           //TODO: replace fullsw by an int version of swlib
//           result->score1 = fullsw(dseq, dseq + dlen, qseq, qseq + qlen, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
//         }
                  
//         matrix_adjust_rule = eDontAdjustMatrix;
//         compo_adjusted_matrix(NRrecordBL62, sbpBL62, scaledMatrixInfoBL62, &subject_composition_unmasked, subject_length, &query.composition, query_length, compo_adjust_mode_stage3, -1, &matrix_adjust_rule, thresh_length, thresh_distance, thresh_angle);
//         // convert adjusted score matrix to swipe matix
//         for (int a = 0; a < BLASTAA_SIZE; a++)
//           for (int b = 0; b < BLASTAA_SIZE; b++) {
//             mata[a * BLASTAA_SIZE + b] = sbpBL62->matrix->data[b][a]; // because of asymmetric matrix, [b][a] is correct for sswlib and fullsw
//             adjusted_score_matrix_63_rev[(a<<5) + b] = sbpBL62->matrix->data[a][b];
//             adjusted_score_matrix_63_fullsw[(a<<5) + b] = sbpBL62->matrix->data[b][a];
//           }
//         p = ssw_init_word((const int8_t*)subject_sequence_unmasked, subject_length, mata, BLASTAA_SIZE);
//         sw_sse2_word2((const int8_t*)query_sequence_unmasked, 0, query_length, p->readLen, gapopenextend_BlastDef, gapextend_BlastDef, p->profile_word, 32767, maskLen, &ae);
//         if (p)
//           init_destroy(p);
//         result->score2 = ae.score;
//         result->ref_end2 = ae.ref;
//         result->read_end2 = ae.read;
//         if (result->score2 == 32767) {
//           flags |= HIT_LARGE_SCORE;
//           char *qseq2 = subject_sequence_unmasked;
//           long qlen2 = subject_length;
//           char *dseq2 = query_sequence_unmasked;
//           long dlen2 = query_length;
//           if (hesize < dlen2 * 32) {
//             hesize = dlen2*32;
//             hearray = (BYTE*)realloc(hearray, hesize);
//           }
//           //TODO: replace fullsw by an int version of swlib
//           result->score2 = fullsw(dseq2, dseq2 + dlen2, qseq2, qseq2 + qlen2, (long*) hearray, adjusted_score_matrix_63_fullsw, gapopenextend_BlastDef, gapextend_BlastDef);
//         }
// //        if (round(result->score2/scaling_factor_BL62) > round(result->score1/scaling_factor_BL62)) {
// //        if (round((float)result->score2/scaling_factor_BL62) > round((float)result->score1/scaling_factor_BL62)) {
//         if (result->score2 > result->score1) {

// #ifdef TIME_PAIRCOUNT
//           pairs3c++;
// #endif //TIME_PAIRCOUNT

//           flags |= HIT_SUBJECT_QUERY_BEST_BL62;
//           //fprintf(stderr, "symmetrical score is larger than adjusted blastp score: %d vs %d\n", (int)round(result->score2/scaling_factor_BL62), (int)round(result->score1/scaling_factor_BL62));
//         }
//         adjusted_score_blast_rev = result->score2;
//       }
      
//       int64_t *used_score_matrix_63 = adjusted_score_matrix_63;
      
//       if (flags & HIT_LARGE_SCORE)
//         matchEnd = queryEnd = score_align = 0;
//       else {
//         if (flags & HIT_SUBJECT_QUERY_BEST_BL62) {
//           adjusted_score_blast = result->score2;
//           matchEnd = result->ref_end2;
//           queryEnd = result->read_end2;
//           adjusted_score_blast_rev = result->score1;
//           used_score_matrix_63 = adjusted_score_matrix_63_rev;
//           qseq = subject_sequence_unmasked;
//           qlen = subject_length;
//           dseq = query_sequence_unmasked;
//           dlen = query_length;
//         } else {
//           adjusted_score_blast = result->score1;
//           matchEnd = result->ref_end1;
//           queryEnd = result->read_end1;
//         }
//         score_align = adjusted_score_blast;
//       }
//       //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
// /*      fprintf(out, "query sequence: \n");
//       db_print_seq_map(qseq, qlen, sym_ncbi_aa);
//       fprintf(out, "subject sequence: \n");
//       db_print_seq_map(dseq, dlen, sym_ncbi_aa);
// */      align(qseq,
//         dseq,
//         qlen,
//         dlen,
//         used_score_matrix_63,
//         gapopen_BlastDef,
//         gapextend_BlastDef,
//         & queryStart,
//         & matchStart,
//         & queryEnd,
//         & matchEnd,
//         & alignment,
//         & score_align);
//         adjusted_score_blast = score_align;
// //          fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", queryStart, queryEnd, matchStart, matchEnd, adjusted_score_blast, adjusted_score_blast/scaling_factor_BL62);
// /*        } else {
//         score_align = adjusted_score_blast;
//         // reverse query sequence; queryEnd is 0-based alignment ending position
//         int32_t start = 0, end = queryEnd;
//         if (end > reverse_len) {
//           if (reverse_len == 0)
//             read_reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));
//           else
//             read_reverse = (int8_t*)realloc(read_reverse, (end + 1) * sizeof(int8_t));
//           reverse_len = end;
//         }
//         while (start <= end) {
//       read_reverse[start] = qseq[end];
//           read_reverse[end] = qseq[start];
//           ++ start;
//           -- end;
//         }
//         __m128i* vP = qP_word16(read_reverse, mata, queryEnd + 1, p->n);
//         alignment_end* bests_reverse = sw_sse2_word2((const int8_t*)subject_sequence, 1, matchEnd + 1, queryEnd + 1, gapopen_BlastDef + gapextend_BlastDef, gapextend_BlastDef, vP, score_align, maskLen);
//         free(vP);
//         matchStart = bests_reverse[0].ref;
//         queryStart = queryEnd - bests_reverse[0].read;
//         result->score1 = bests_reverse[0].score;
//         free(bests_reverse);
//         matchStart_r = matchEnd_r = queryStart_r = queryEnd_r = score_align_r = 0;
//       align(&qseq[queryStart],
//         &subject_sequence[matchStart],
//         queryEnd - queryStart + 1,
//         matchEnd - matchStart + 1,
//         adjusted_score_matrix_63,
//         gapopen_BlastDef,
//         gapextend_BlastDef,
//         & queryStart_r,
//         & matchStart_r,
//         & queryEnd_r,
//         & matchEnd_r,
//         & alignment,
//         & score_align_r);
//         //fprintf(out, "%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", queryStart, queryEnd, matchStart, matchEnd, score_align_r, result->score1);
//       }*/
//       hits_set_align_string(i, alignment, score_align);
//       hits_enter_align_coord(i, queryStart, queryEnd, matchStart, matchEnd, 0);
// //        score_matrix_63 = used_score_matrix_63;
// //        score_matrix_dump();
//       count_align_matrix(i, used_score_matrix_63, qseq, qlen, dseq, dlen);
// //        show_align(i);

// #ifdef TIME_PAIRCOUNT
//       // time stage 3
//       auto end3 = std::chrono::high_resolution_clock::now();
//       time_stage3 += (end3 - begin3).count();
// #endif //TIME_PAIRCOUNT

//     } else
//       hits_enter_align_coord(i, queryStart, queryEnd, matchStart, matchEnd, 0);
//     hits_enter_adjusted_score(i, adjusted_score, adjusted_score_blast, adjusted_score_blast_rev, flags);
//   }
//   db_thread_destruct(dbt);
// }

void align_adjusted_done() {
  if (result)
    align_destroy(result);

  align_end_destroy(&ae);
  free(mata);
  free(mat8);
  if (hearray)
    free(hearray);
//  free(read_reverse);

  Blast_CompositionWorkspaceFree(&NRrecord);
  compo_done(&sbp, &scaledMatrixInfo);

  Blast_CompositionWorkspaceFree(&NRrecordBL62);
  compo_done(&sbpBL62, &scaledMatrixInfoBL62);
  if (mask)
    delete(dbSequences);
  free(adjusted_score_matrix_63);
  free(adjusted_score_matrix_63_rev);
  free(adjusted_score_matrix_63_fullsw);
}
//#endif

void args_getstring(int i, int argc, char **argv, char ** result, char * error)
{
  if (i+1 < argc)
    *result = argv[i+1];
  else
    fatal(error);
}

void args_getnum(int i, int argc, char **argv, long * result, char * error)
{
  if (i+1 < argc)
    *result = atol(argv[i+1]);
  else
    fatal(error);
}

// convert arg (as word or number) into ECompoAdjustModes
ECompoAdjustModes convert(const std::string& arg)
{
    if(arg == "eNoCompositionBasedStats" || arg == "0") return eNoCompositionBasedStats;
    else if(arg == "eCompositionBasedStats" || arg == "1") return eCompositionBasedStats;
    else if(arg == "eCompositionMatrixAdjust" || arg == "2") return eCompositionMatrixAdjust;
    else if(arg == "eCompoForceFullMatrixAdjust" || arg == "3") return eCompoForceFullMatrixAdjust;
    else fprintf(stderr, "Unknown ECompoAdjustModes value");
}

void args_show()
{
  if (view == 0)
  {

#ifndef MPISWIPE
    if (! cpu_feature_ssse3)
    {
      fprintf(out, "The performance is reduced because this CPU lacks SSSE3.\n\n");
    }
#endif

    const char * symtypestring[] = { "Nucleotide", "Amino acid", "Translated query", "Translated database", "Both translated", "Sound" };

    //      char * viewtypestring[] = { "plain", 0, 0, 0, 0, 0, 0, "xml",
    //			  "tab-separated", "tab-separated with comments" };

    fprintf(out, "Database file:     %s\n", databasename);
    fprintf(out, "Database title:    %s\n", db_gettitle());
    fprintf(out, "Database time:     %s\n", db_gettime());

    if (db_ismasked())
      {
	fprintf(out, "Database size:     %ld residues", db_getsymcount_masked());
	fprintf(out, " in %ld sequences\n", db_getseqcount_masked());
      }
      else
      {
	fprintf(out, "Database size:     %ld residues", db_getsymcount());
	fprintf(out, " in %ld sequences\n", db_getseqcount());
      }

      fprintf(out, "Longest db seq:    %ld residues\n", db_getlongest());

      if (effdbsize > 0)
	fprintf(out, "Effecive db size:  %ld\n", effdbsize);

      fprintf(out, "Query file name:   %s\n", queryname);

      long qlen;
      if ((symtype == 0) || (symtype == 2) || (symtype == 4))
	qlen = query.nt[0].len;
      else
	qlen = query.aa[0].len;

      fprintf(out, "Query length:      %ld residues\n", qlen);

      query_show();

      if (symtype == 0)
      {
	fprintf(out, "Query strands:     ");
	switch (querystrands)
	{
	case 1:
	  fprintf(out, "Plus");
	  break;
	case 2:
	  fprintf(out, "Minus");
	  break;
	case 3:
	  fprintf(out, "Plus and minus");
	  break;
	}
	fprintf(out, "\n");
	fprintf(out, "Score matrix:      %ld/%ld\n", matchscore, mismatchscore);
      }
      else
      {
	fprintf(out, "Score matrix:      %s\n", matrixname);
      }

      fprintf(out, "Gap penalty:       %ld+%ldk\n", gapopen, gapextend);
      fprintf(out, "Max expect shown:  %-g\n", expect);
      fprintf(out, "Min score shown:   %ld\n", minscore);
      fprintf(out, "Max matches shown: %ld\n", maxmatches);
      fprintf(out, "Alignments shown:  %ld\n", alignments);
      fprintf(out, "Show gi's:         %ld\n", show_gis);
      fprintf(out, "Show taxid's:      %ld\n", show_taxid);
      fprintf(out, "Threads:           %d\n", threads);
      fprintf(out, "Symbol type:       %s\n", symtypestring[symtype]);
      if ((symtype == 2) || (symtype == 4))
	fprintf(out, "Query genetic code:%s (%ld)\n", gencode_names[query_gencode-1], query_gencode);
      if ((symtype == 3) || (symtype == 4))
	fprintf(out, "DB genetic code:   %s (%ld)\n", gencode_names[db_gencode-1], db_gencode);

#if 0
      if ((symtype == 2) || (symtype == 4))
      {
	for(long s=0; s<2; s++)
	  for(long f=0; f<3; f++)
	  {
	    fprintf(out, "Translation of query, frame %c%ld:\n", s?'-':'+', f+1);
	    //	    translate(query_, qlen, s, f, 0, & prot, & plen);

	    long plen = query.aa[3*s+f].len;
	    for (int j=0; j<plen; j+=60)
	    {
	      for(int k=0; k<60; k++)
	      {
		if (j+k < plen)
		  putc(sym_ncbi_aa[query.aa[3*s+f].seq[j+k]], out);
		else
		  break;
	      }
	      fprintf(out, "\n");
	    }
	    fprintf(out, "\n");
	  }
      }
#endif

      // fprintf(out, "View:              %s\n", viewtypestring[view]);
      if (taxidfilename)
	fprintf(out, "Taxid filename:    %s\n", taxidfilename);
      fprintf(out, "\n");
    }
}

void args_usage()
{
  /* options unused by BLAST: chkuxHN */
  /* options used by SWIPE:   chkuxHN  */

  fprintf(out, "Usage: %s [OPTIONS]\n", progname);
  fprintf(out, "  -h, --help                 show help\n");
  fprintf(out, "  -d, --db=FILE              sequence database base name (required)\n");
  fprintf(out, "  -i, --query=FILE           query sequence filename (stdin)\n");
  fprintf(out, "  -M, --matrix=NAME/FILE     score matrix name or filename (BLOSUM62)\n");
  fprintf(out, "  -q, --penalty=NUM          penalty for nucleotide mismatch (-3)\n");
  fprintf(out, "  -r, --reward=NUM           reward for nucleotide match (1)\n");
  fprintf(out, "  -G, --gapopen=NUM          gap open penalty (11)\n");
  fprintf(out, "  -E, --gapextend=NUM        gap extension penalty (1)\n");
  fprintf(out, "  -v, --num_descriptions=NUM sequence descriptions to show (250)\n");
  fprintf(out, "  -b, --num_alignments=NUM   sequence alignments to show (100)\n");
  fprintf(out, "  -e, --evalue=REAL          maximum expect value of sequences to show (10.0)\n");
  fprintf(out, "  -k, --minevalue=REAL       minimum expect value of sequences to show (0.0)\n");
  fprintf(out, "  -c, --min_score=NUM        minimum score of sequences to show (1)\n");
  fprintf(out, "  -B, --min_score2=NUM       minimum score of sequences to compute unmasked Blast score\n");
  fprintf(out, "  -u, --max_score=NUM        maximum score of sequences to show (inf.)\n");
  fprintf(out, "  -a, --num_threads=NUM      number of threads to use [1-%d] (1)\n", MAX_THREADS);
  fprintf(out, "  -m, --outfmt=NUM           output format [0,7-9=plain,xml,tsv,tsv+,88=score only] (0)\n");
  fprintf(out, "  -I, --show_gis             show gi numbers in results (no)\n");
  fprintf(out, "  -p, --symtype=NAME/NUM     symbol type/translation [0-4] (1)\n");
  fprintf(out, "  -S, --strand=NAME/NUM      query strands to search [1-3] (3)\n");
  fprintf(out, "  -Q, --query_gencode=NUM    query genetic code [1-23] (1)\n");
  fprintf(out, "  -D, --db_gencode=NUM       database genetic code [1-23] (1)\n");
  fprintf(out, "  -x, --taxidlist=FILE       taxid list filename (none)\n");
  fprintf(out, "  -N, --dump=NUM             dump database [0-2=no,yes,split headers] (0)\n");
  fprintf(out, "  -H, --show_taxid           show taxid etc in results (no)\n");
  fprintf(out, "  -o, --out=FILE             output file (stdout)\n");
  fprintf(out, "  -z, --dbsize=NUM           set effective database size (0)\n");
  fprintf(out, "  -s, --soft_masking         lower case regions of sequences are masked before alignment\n");
  fprintf(out, "  -A, --compo_adjust_mode=NAME/NUM  mode of composition adjustment\n");
  fprintf(out, "  -n, --no_compo_stats       replace composition based statistics by no adjustment (only applies to compo_adjust_mode=2 (conditional))\n");
  fprintf(out, "  -X, --skip_stage2_score=NUM       stage1 scores (max. 127) greater than this threshold go directly to stage3 (bypassing stage2)\n");
  fprintf(out, "  -U, --thresh_length=REAL    length_ratio_threshold used in compo_mode_condition.c\n");
  fprintf(out, "  -V, --thresh_distance=REAL  query_match_distance_threshold used in compo_mode_condition.c\n");
  fprintf(out, "  -W, --thresh_angle=REAL     angle_degree_threshold used in compo_mode_condition.c\n");
  fprintf(out, "  -t, --temporary_dir=DIR     temporary directory path\n");
  fprintf(out, "  -w, --ssw_binary_dir=DIR    directory path of SSW binaries (ssw_main_sw_sse2_byte_multimatrix and ssw_main_sw_sse2_word_multimatrix)\n");
  fprintf(out, "  -R, --align_adjusted_version=NUM  align_adjusted_version (4)\n");
}

void args_help()
{
  char title[] = "SWIPE " SWIPE_VERSION;
  char ref[] = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
  fprintf(out, "%s [%s %s]\n\n%s\n\n", title, __DATE__, __TIME__, ref);

  args_usage();
}

void args_init(int argc, char **argv)
{
  /* Set defaults */
  gapopen = DEFAULT_GAPOPEN;
  gapextend = DEFAULT_GAPEXTEND;
  matrixname = "";
  queryname = DEFAULT_QUERYNAME;
  databasename = DEFAULT_DATABASENAME;
  minscore = DEFAULT_MINSCORE;
  minscore2 = DEFAULT_MINSCORE;
  maxscore = DEFAULT_MAXSCORE;
  maxmatches = DEFAULT_MAXMATCHES;
  alignments = DEFAULT_ALIGNMENTS;
  threads = DEFAULT_THREADS;
  view = DEFAULT_VIEW;
  symtype = DEFAULT_SYMTYPE;
  show_gis = DEFAULT_SHOW_GIS;
  show_taxid = DEFAULT_SHOW_TAXID;
  expect = DEFAULT_EXPECT;
  minexpect = DEFAULT_MINEXPECT;
  taxidfilename = NULL;
  matchscore = DEFAULT_MATCHSCORE;
  mismatchscore = DEFAULT_MISMATCHSCORE;
  querystrands = DEFAULT_QUERYSTRANDS;
  query_gencode = DEFAULT_QUERY_GENCODE;
  db_gencode = DEFAULT_DB_GENCODE;
  subalignments = DEFAULT_SUBALIGNMENTS;
  dump = DEFAULT_DUMP;
  effdbsize = DEFAULT_EFFDBSIZE;
  mask = DEFAULT_MASK;
  compo_adjust_mode = DEFAULT_COMPO_ADJUST_MODE;
  no_compo_stats = DEFAULT_NO_COMPO_STATS;
  skip_stage2_score = DEFAULT_SKIP_STAGE2_SCORE;
  thresh_length = DEFAULT_THRESH_LENGTH;
  thresh_distance = DEFAULT_THRESH_DISTANCE;
  thresh_angle = DEFAULT_THRESH_ANGLE;
  temporary_dir = DEFAULT_TEMPORARY_DIR;
  ssw_binary_dir = DEFAULT_SSW_BINARY_DIR;
  align_adjusted_version = DEFAULT_ALIGN_ADJUSTED_VERSION;

  progname = argv[0];

  opterr = 1;
  char short_options[] = "d:i:M:q:r:G:E:S:v:b:c:B:u:e:k:a:m:p:x:C:Q:D:F:K:N:o:z:A:X:U:V:W:t:w:R:IHhns::";

  static struct option long_options[] =
  {
    {"db",               required_argument, NULL, 'd' },
    {"query",            required_argument, NULL, 'i' },
    {"matrix",           required_argument, NULL, 'M' },
    {"penalty",          required_argument, NULL, 'q' },
    {"reward",           required_argument, NULL, 'r' },
    {"gapopen",          required_argument, NULL, 'G' },
    {"gapextend",        required_argument, NULL, 'E' },
    {"strand",           required_argument, NULL, 'S' },
    {"num_descriptions", required_argument, NULL, 'v' },
    {"num_alignments",   required_argument, NULL, 'b' },
    {"min_score",        required_argument, NULL, 'c' },
    {"min_score2",       required_argument, NULL, 'B' },
    {"max_score",        required_argument, NULL, 'u' },
    {"evalue",           required_argument, NULL, 'e' },
    {"minevalue",        required_argument, NULL, 'k' },
    {"num_threads",      required_argument, NULL, 'a' },
    {"outfmt",           required_argument, NULL, 'm' },
    {"symtype",          required_argument, NULL, 'p' },
    {"taxid",            required_argument, NULL, 'x' },
    {"comp_based_stats", required_argument, NULL, 'C' },
    {"query_gencode",    required_argument, NULL, 'Q' },
    {"db_gencode",       required_argument, NULL, 'D' },
    {"filter",           required_argument, NULL, 'F' },
    {"subalignments",    required_argument, NULL, 'K' },
    {"dump",             required_argument, NULL, 'N' },
    {"out",              required_argument, NULL, 'o' },
    {"dbsize",           required_argument, NULL, 'z' },
    {"show_gis",         no_argument,       NULL, 'I' },
    {"show_taxid",       no_argument,       NULL, 'H' },
    {"soft_masking",     optional_argument, NULL, 's' },
    {"compo_adjust_mode",required_argument, NULL, 'A' },
    {"no_compo_stats",   no_argument,       NULL, 'n' },
    {"skip_stage2_score",required_argument, NULL, 'X' },
    {"thresh_length",    required_argument, NULL, 'U' },
    {"thresh_distance",  required_argument, NULL, 'V' },
    {"thresh_angle",     required_argument, NULL, 'W' },
    {"temporary_dir",    required_argument, NULL, 't' },
    {"ssw_binary_dir",   required_argument, NULL, 'w' },
    {"align_adjusted_version",    required_argument, NULL, 'R' },
    {"help",             no_argument,       NULL, 'h' },
    { 0, 0, 0, 0 }
  };

  int option_index = 0;
  int c;

  while (1)
    {
      c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
	break;

      switch(c)
	{
	case 'a':
	  /* threads */
	  threads = atol(optarg);
	  break;

	case 'b':
	  /* alignments */
	  alignments = atol(optarg);
	  break;

	case 'c':
	  /* min score threshold */
	  minscore = atol(optarg);
	  break;

	case 'B':
	  /* min score threshold */
	  minscore2 = atol(optarg);
	  break;

	case 'C':
	  /* composition-based adjustments */
	  if ( (strcasecmp(optarg, "F") != 0) && (strcmp(optarg, "0") != 0) )
	    fatal("Composition-based score adjustments not supported.");
	  break;

	case 'd':
	  /* database */
	  databasename = optarg;
	  break;

	case 'D':
	  /* database genetic code */
	  db_gencode = atol(optarg);
	  break;

	case 'e':
	  /* evalue */
	  expect = atof(optarg);
	  break;

	case 'E':
	  /* gap extend */
	  gapextend = atol(optarg);
	  break;

	case 'F':
	  /* filter */
	  if ( (strlen(optarg) != 0) && (strcasecmp(optarg, "F") != 0) )
	    fatal("Query sequence filtering not supported.");
	  break;

	case 'G':
	  /* gap open */
	  gapopen = atol(optarg);
	  break;

	case 'h':
	  args_help();
	  exit(1);
	  break;

	case 'H':
	  /* show_taxid */
	  show_taxid = 1;
	  break;

	case 'i':
	  /* query */
	  queryname = optarg;
	  break;

	case 'I':
	  /* show_gis */
	  show_gis = 1;
	  break;

	case 'k':
	  /* min evalue threshold */
	  minexpect = atof(optarg);
	  break;

	case 'K':
	  /* subalignments */
	  subalignments = atol(optarg);
	  break;

	case 'm':
	  /* view */
	  view = atol(optarg);
	  break;

	case 'M':
	  /* matrix */
	  matrixname = optarg;
	  break;

	case 'N':
	  /* dump */
	  dump = atol(optarg);
	  break;

	case 'o':
	  /* output file */
	  outfile = optarg;
	  break;

	case 'p':
	  /* symtype */
	  if (strcmp(optarg, "blastn") == 0)
	    symtype = 0;
	  else if (strcmp(optarg, "blastp") == 0)
	    symtype = 1;
	  else if (strcmp(optarg, "blastx") == 0)
	    symtype = 2;
	  else if (strcmp(optarg, "tblastn") == 0)
	    symtype = 3;
	  else if (strcmp(optarg, "tblastx") == 0)
	    symtype = 4;
	  else if (strcmp(optarg, "sound") == 0)
	    symtype = 5;
	  else
	    symtype = atol(optarg);
	  break;

	case 'q':
	  /* penalty */
	  mismatchscore = atol(optarg);
	  break;

	case 'Q':
	  /* query genetic code */
	  query_gencode = atol(optarg);
	  break;

	case 'r':
	  /* reward */
	  matchscore = atol(optarg);
	  break;

	case 'S':
	  if (strcmp(optarg, "plus") == 0)
	    querystrands = 1;
	  else if (strcmp(optarg, "minus") == 0)
	    querystrands = 2;
	  else if (strcmp(optarg, "both") == 0)
	    querystrands = 3;
	  else
	    querystrands = atol(optarg);
	  break;

	case 'u':
	  /* maxscore */
	  maxscore = atol(optarg);
	  break;

	case 'v':
	  /* max matches shown */
	  maxmatches = atol(optarg);
	  break;

	case 'x':
	  /* taxid filename */
	  taxidfilename = optarg;
	  break;

	case 'z':
	  /* effective db size */
	  effdbsize = atol(optarg);
	  break;

	case 's':
	  /* soft_masked */
	  if (optarg)
        mask = atol(optarg);
      else
	    mask = COMPOSITIONAL_MASK_COMP;
	  break;

  case 'A':
	  /* compo_adjust_mode 
    eNoCompositionBasedStats       = 0, 
    eCompositionBasedStats         = 1, 
    eCompositionMatrixAdjust       = 2, // default
    eCompoForceFullMatrixAdjust    = 3*/
    compo_adjust_mode = convert(optarg);
	  break;
  
  case 'n':
	  /* no_compo_stats */
	  no_compo_stats = 1;
	  break;
  
  case 'X':
	  /* skip_stage2_score*/
    skip_stage2_score = atol(optarg);
	  break;

  case 'U':
	  /* thresh_length*/
    thresh_length = atof(optarg);
	  break;

  case 'V':
	  /* thresh_distance*/
    thresh_distance = atof(optarg);
	  break;

  case 'W':
	  /* thresh_angle*/
    thresh_angle = atof(optarg);
	  break;

  case 't':
	  /* temporary_dir */
	  temporary_dir = optarg;
	  break;

  case 'w':
	  /* ssw_binary_dir */
	  ssw_binary_dir = optarg;
	  break;

  case 'R':
	  /* align_adjusted_version */
	  align_adjusted_version = atol(optarg);
	  break;

	case '?':
	default:
	  args_usage();
	  exit(1);
	  break;
	}
    }

  if (outfile)
  {
    FILE * f = fopen(outfile, "w");
    if (! f)
      fatal("Unable to open output file for writing.");
    out = f;
  }

  long gopen_default;
  long gextend_default;

  if (symtype == 0)
  {
    if (gapopen == 0)
      gapopen = 5;
    if (gapextend == 0)
      gapextend = 2;
  }
  else if (symtype < 5)
  {
    if (strlen(matrixname) == 0)
      matrixname = DEFAULT_MATRIXNAME;

    if (stats_getprefs(matrixname, & gopen_default, & gextend_default))
    {
      if (gapopen == 0)
	gapopen = gopen_default;
      if (gapextend == 0)
	gapextend = gextend_default;
    }
    else
    {
      if ((gapopen == 0) && (gapextend == 0))
	fatal("Unknown score matrix. Gap penalties must be specified (-G and -E).");
    }
  }
  else if (symtype == 5)
  {
    if (strlen(matrixname) == 0)
      matrixname = "IDENTITY_5_1";
    if (gapopen == 0)
      gapopen = 15;
    if (gapextend == 0)
      gapextend = 5;
  }

  gapopenextend = gapopen + gapextend;

  if (effdbsize < 0)
    fatal("Illegal effective db size specified");

  if ((threads < 1) || (threads > MAX_THREADS))
    fatal("Illegal number of threads specified");

  if (strlen(databasename) == 0)
    fatal("No database specified.");

  if (!((view==0)||(view==7)||(view==8)||(view==9)||(view==88)||(view==99)))
    fatal("Illegal view type.");

  if ((gapopen < 0) || (gapextend < 0) || ((gapopen + gapextend) < 1))
    fatal("Illegal gap penalties.");

  if ((symtype < 0) || (symtype > 5))
    fatal("Illegal symbol type.");

  if ((querystrands < 1) || (querystrands > 3))
    fatal("Illegal query strands specified.");

  if ((querystrands == 2) && ((symtype == 1) || (symtype == 3) || (symtype == 4)))
    fatal("Illegal strand specified for protein query.");

  if ((query_gencode < 1)  || (query_gencode > 23) || (! gencode_names[query_gencode-1]))
    fatal("Illegal query genetic code specified.");

  if ((db_gencode < 1) || (db_gencode > 23) || (! gencode_names[db_gencode-1]))
    fatal("Illegal database genetic code specified.");

  if ((dump<0) || (dump>2))
    fatal("Illegal dump mode.");

  translate_init(query_gencode, db_gencode);
}


void vector_fill(BYTE * vector, BYTE value)
{
  memset(vector, value, 16);
}

void vector_print(BYTE * vector)
{
  for(int i=0; i<16; i++)
    fprintf(out, " %02x", vector[i]);
}

void vector_print_word(WORD * vector)
{
  for(int i=0; i<8; i++)
    fprintf(out, " %04x", vector[i]);
}


void search_init(struct search_data * sdp)
{
  sdp->dbt = db_thread_create();
  sdp->dprofile = (BYTE*) xmalloc(4*16*32);
  long qlen = 0;
  long hearraylen = 0;

  for(int i = 0; i < 6; i++)
  {
    sdp->qtable[i] = 0;
  }

  if (symtype == 0)
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
      {
	qlen = query.nt[s].len;
	sdp->qlen[3*s] = qlen;
	sdp->qtable[3*s] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	for(int i=0; i<qlen; i++)
	{
	  sdp->qtable[3*s][i] = sdp->dprofile + 64*query.nt[s].seq[i];
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
  }
  else if ((symtype == 1) || (symtype == 3) || (symtype == 5))
  {
    qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
    for(int i=0; i<qlen; i++)
    {
      sdp->qtable[0][i] = sdp->dprofile + 64*query.aa[0].seq[i];
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((symtype == 2) || (symtype == 4))
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
	for(int f=0; f<3; f++)
	{
	  qlen = query.aa[3*s+f].len;
	  sdp->qlen[3*s+f] = qlen;
	  sdp->qtable[3*s+f] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	  for(int i=0; i<qlen; i++)
	  {
	    sdp->qtable[3*s+f][i] = sdp->dprofile + 64*query.aa[3*s+f].seq[i];
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
  }

  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray = (BYTE*) xmalloc(hearraylen*32);

  long listsize = maxchunksize * sizeof(long);
  if ((symtype == 3) || (symtype == 4))
    listsize *= 6;

  sdp->start_list = (long*) xmalloc(listsize);
  sdp->in_list = (long*) xmalloc(listsize);
  sdp->out_list = (long*) xmalloc(listsize);
  sdp->scores = (long*) xmalloc(listsize);
  sdp->bestpos = (long*) xmalloc(listsize);
  sdp->bestq = (long*) xmalloc(listsize);

  if (symtype == 0)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 2)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 3)
  {
    sdp->qstrand1 = 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = 0;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 1;
    sdp->dframe2 = 2;
  }
  else if (symtype == 4)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 1;
    sdp->dframe2 = 2;
  }
  else
  {
    sdp->qstrand1 = 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = 0;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }

}

void search_done(struct search_data * sdp)
{
  for(int i = 0; i < 6; i++)
  {
    if (sdp->qtable[i])
      free(sdp->qtable[i]);
  }

  free(sdp->dprofile);
  free(sdp->hearray);
  free(sdp->scores);
  free(sdp->bestpos);
  free(sdp->bestq);
  free(sdp->start_list);
  free(sdp->in_list);
  free(sdp->out_list);
  db_thread_destruct(sdp->dbt);
}

int search_getwork(long * first, long * last)
{
  int status = 0;
  long volcount = db_getvolumecount();

  pthread_mutex_lock(&workmutex);
  if (volnext < volcount)
  {
    long seqcount = volseqs[volnext];
    long chunks = volchunks[volnext];
    long chunksize = ((seqcount+chunks-1) / chunks);

    * first = seqnext;
    * last = seqnext + chunksize - 1;
    seqnext += chunksize;
    status = 1;

    //    fprintf(out, "Processing sequences %d to %d (%d sequences) in volume %ld.\n", *first, *last, *last - * first + 1, volnext);

    volseqs[volnext] -= chunksize;
    volchunks[volnext]--;

    while ((volnext < volcount) && (volchunks[volnext] == 0))
      volnext++;
  }
  pthread_mutex_unlock(&workmutex);
  return status;
}


void search_chunk(struct search_data * sdp)
{
  //  fprintf(out, "Searching seqnos %ld to %ld\n", sdp->seqfirst, sdp->seqlast);

  if(taxidfilename)
    db_mapheaders(sdp->dbt, sdp->seqfirst, sdp->seqlast);

  sdp->start_count = 0;
  for(long seqno = sdp->seqfirst; seqno <= sdp->seqlast; seqno++)
  {
    if (db_check_inclusion(sdp->dbt, seqno))
    {
      if ((symtype == 3) || (symtype == 4))
      {
	for(long dstrand = sdp->dstrand1; dstrand <= sdp->dstrand2; dstrand++)
	  for(long dframe = sdp->dframe1; dframe <= sdp->dframe2; dframe++)
	  {
	    sdp->start_list[sdp->start_count++] =
	      (seqno << 3) | (dstrand << 2) | dframe;
	  }
      }
      else
      {
	sdp->start_list[sdp->start_count++] = seqno << 3;
      }
    }
  }

  if (sdp->start_count == 0)
    return;

  long s1 = sdp->start_list[0] >> 3;
  long s2 = sdp->start_list[sdp->start_count-1] >> 3;

  // fprintf(out, "Mapping seqnos %ld to %ld\n", s1, s2);

  db_mapsequences(sdp->dbt, s1, s2);

  for (long qstrand = sdp->qstrand1; qstrand <= sdp->qstrand2; qstrand++)
    for(long qframe = sdp->qframe1; qframe <= sdp->qframe2; qframe++)
    {
      long dstrand, dframe;

      sdp->out_count = sdp->start_count;
      memcpy(sdp->out_list, sdp->start_list, sdp->start_count * sizeof(long));

      BYTE ** qtable = sdp->qtable[3*qstrand+qframe];
      long qlen = sdp->qlen[3*qstrand+qframe];

#if 1

      /* 7-bit search */

      sdp->tmp_list = sdp->in_list;
      sdp->in_list = sdp->out_list;
      sdp->out_list = sdp->tmp_list;
      sdp->in_count = sdp->out_count;

      if (sdp->in_count > 0)
      {
	pthread_mutex_lock(&countmutex);
	compute7 += sdp->in_count;
	rounds7++;
	pthread_mutex_unlock(&countmutex);

	// fprintf(out, "Searching seqnos %ld to %ld\n", sdp->in_list[0], sdp->in_list[sdp->in_count-1]);

	if (cpu_feature_ssse3)
	  search7_ssse3(qtable,
			gapopenextend,
			gapextend,
			(BYTE*) score_matrix_7,
			sdp->dprofile,
			sdp->hearray,
			sdp->dbt,
			sdp->in_count,
			sdp->in_list,
			sdp->scores,
			qlen);
	else
	  search7(qtable,
		  gapopenextend,
		  gapextend,
		  (BYTE*) score_matrix_7,
		  sdp->dprofile,
		  sdp->hearray,
		  sdp->dbt,
		  sdp->in_count,
		  sdp->in_list,
		  sdp->scores,
		  qlen);

	sdp->out_count = 0;

	for (int i=0; i<sdp->in_count; i++)
	{

#ifdef TIME_PAIRCOUNT
    pairs1++;
#endif //TIME_PAIRCOUNT

	  long seqnosf = sdp->in_list[i];
	  long score = sdp->scores[i];

	  if (score >= minscore)  // if (score >= minscore && score < SCORELIMIT_7) 
	  {
	    long seqno = seqnosf >> 3;
	    dstrand = (seqnosf >> 2) & 1;
	    dframe = seqnosf & 3;

	    if ((symtype == 0) && qstrand)
	      hits_enter(seqno, score, 0, 0, 1, 0, -1, -1);
	    else
	      hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, -1, -1);
	  }
	  // else
	  // {
	  //   sdp->out_list[sdp->out_count++] = seqnosf;
	  // }
	}
      }
  return;
#endif

#if 1

      /* 16-bit search */

      sdp->tmp_list = sdp->in_list;
      sdp->in_list = sdp->out_list;
      sdp->out_list = sdp->tmp_list;
      sdp->in_count = sdp->out_count;

      if (sdp->in_count > 0)
      {
	pthread_mutex_lock(&countmutex);
	compute16 += sdp->in_count;
	rounds16++;
	pthread_mutex_unlock(&countmutex);

	search16((WORD**)qtable,
		 gapopenextend,
		 gapextend,
		 (WORD*)(score_matrix_16),
		 (WORD*)(sdp->dprofile),
		 (WORD*)(sdp->hearray),
		 sdp->dbt,
		 sdp->in_count,
		 sdp->in_list,
		 sdp->scores,
		 sdp->bestpos,
		 qlen);

	sdp->out_count = 0;

	for (int i=0; i<sdp->in_count; i++)
	{
	  long seqnosf = sdp->in_list[i];
	  long score = sdp->scores[i];
	  if (score < SCORELIMIT_16)
	  {
	    long seqno = seqnosf >> 3;
	    dstrand = (seqnosf >> 2) & 1;
	    dframe = seqnosf & 3;

	    long pos = sdp->bestpos[i];

	    //	    fprintf(out, "seqno=%ld score=%ld bestpos=%ld\n", seqno, score, pos);

	    if ((symtype == 0) && qstrand)
	      hits_enter(seqno, score, 0, 0, 1, 0, pos, -1);
	    else
	      hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, pos, -1);
	  }
	  else
	  {
	    sdp->out_list[sdp->out_count++] = seqnosf;
	  }
	}
      }
#endif

#if 1

      /* 63-bit search */

      sdp->tmp_list = sdp->in_list;
      sdp->in_list = sdp->out_list;
      sdp->out_list = sdp->tmp_list;
      sdp->in_count = sdp->out_count;

      if (sdp->in_count > 0)
      {
	pthread_mutex_lock(&countmutex);
	compute63 += sdp->in_count;
	rounds63++;
	pthread_mutex_unlock(&countmutex);

	for (int i=0; i<sdp->in_count; i++)
	{
	  long seqnosf = sdp->in_list[i];
	  long seqno = seqnosf >> 3;
	  dstrand = (seqnosf >> 2) & 1;
	  dframe = seqnosf & 3;

	  char * address;
	  long length;
	  long ntlen;
	  db_getsequence(sdp->dbt, seqno, dstrand, dframe,
			 & address, & length, & ntlen, 0);
	  char * dbegin = address;
	  char * dend = address + length - 1;

	  char * q;
	  if (symtype == 0)
	    q = query.nt[qstrand].seq;
	  else
	    q = query.aa[3*qstrand+qframe].seq;

	  long score = fullsw(dbegin,
			      dend,
			      (char*) q,
			      (char*) q + qlen,
			      (long*) sdp->hearray,
			      score_matrix_63,
			      gapopenextend,
			      gapextend);

	  if ((symtype == 0) && qstrand)
	    hits_enter(seqno, score, 0, 0, 1, 0, -1, -1);
	  else
	    hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, -1, -1);
	}
      }

#endif
    }
}


void * worker(void *)
{
  struct search_data sd;
  search_init(&sd);

  while(search_getwork(&sd.seqfirst, &sd.seqlast))
    search_chunk(&sd);

  search_done(&sd);
  return 0;
}


void prepare_search(long par)
{
  volnext = 0;
  seqnext = 0;

#if 1

  long volcount = db_getvolumecount();
  for(long v = 0; v < volcount; v++)
    volseqs[v] = db_getseqcount_volume(v);

  long totalchunks;

  calc_chunks(volcount,
	      par,
	      16,
	      volseqs,
	      volchunks,
	      & totalchunks,
	      & maxchunksize);

  while ((volnext < volcount) && (volchunks[volnext] == 0))
    volnext++;

#else
  long seqcount = db_getseqcount();
  long chunkcount;
  if (seqcount >= 16*par)
  {
    chunkcount = par * (long) (floor(sqrt((1.0 * seqcount) / (16 * par))));
  }
  else if (seqcount >= par)
  {
    chunkcount = par;
  }
  else
  {
    chunkcount = seqcount;
  }
//  fprintf(out, "Chunkcount: %ld\n", chunkcount);
  long volcount = db_getvolumecount();
  long rest_chunks = chunkcount;
  long rest_seqcount = seqcount;
  //  fprintf(out, "Seqcount, chunkcount, ratio (total): %ld, %ld, %f\n", seqcount, chunkcount, 1.0 * seqcount / chunkcount);
  maxchunksize = 0;
  for(long i = 0; i < volcount; i++)
  {
    long volsize = db_getseqcount_volume(i);
    long volchunk = ((rest_chunks * volsize) + rest_seqcount - 1) / rest_seqcount;
    volseqs[i] = volsize;
    volchunks[i] = volchunk;
    long chunksize = (volsize + volchunk -1) / volchunk;
    if (chunksize > maxchunksize)
      maxchunksize = chunksize;
    // fprintf(out, "Seqcount, chunkcount, ratio (%ld): %ld, %ld, %f\n", i, volsize, volchunk, 1.0 * volsize / volchunk);
    rest_chunks -= volchunk;
    rest_seqcount -= volsize;
  }
#endif
}

void run_threads()
{
  unsigned int t;
  void * status;

  for(t=0; t<threads; t++)
    {
      if (pthread_create(pthread_id + t, 0, worker, &t))
	fatal("Cannot create thread.");
    }

  //TODO: create the realigning thread heare
  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      fatal("Cannot join thread.");
  }
}

#define cpuid(l1,l2,a,b,c,d)						\
  __asm__ __volatile__							\
    ("cpuid": "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (l1), "c" (l2));

void cpu_features()
{
  unsigned int a __attribute__ ((unused));
  unsigned int b __attribute__ ((unused));
  unsigned int c,d;
  cpuid(1,0,a,b,c,d);
  cpu_feature_sse2  = (d >> 26) & 1;
  cpu_feature_ssse3 = (c >>  9) & 1;
  cpu_feature_sse41 = (c >> 19) & 1;
}

void clock_start(struct time_info * tip)
{
  time(& tip->t1);                 /* time(2)   */
#ifdef _WIN32
  tip->clk_tck = 1000;
  tip->wc1 = GetTickCount();
#else
  tip->clk_tck = sysconf(_SC_CLK_TCK);
  tip->wc1 = times(& tip->times1); /* times (2) */
#endif // _WIN32
}

void clock_stop(struct time_info * tip)
{
  struct tm tms;
  char buf[30];
  char timeformat[] = "%a, %e %b %Y %T UTC";

#ifdef _WIN32
  tip->wc2 = GetTickCount();
#else
  tip->wc2 = times(& tip->times2);
#endif // _WIN32
  time(& tip->t2);

  gmtime_r(&tip->t1, & tms);
  strftime(buf, 30, timeformat, & tms);
  tip->starttime = (char*) xmalloc(30);
  strcpy(tip->starttime, buf);

  gmtime_r(&tip->t2, & tms);
  strftime(buf, 30, timeformat, & tms);
  tip->endtime = (char*) xmalloc(30);
  strcpy(tip->endtime, buf);

  tip->elapsed = ((double)(tip->wc2 - tip->wc1)) / tip->clk_tck;

  double speed = ((double)db_getsymcount_masked());

  if (symtype == 0)
  {
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  else if (symtype == 1)
  {
    speed *= query.aa[0].len;
  }
  else if (symtype == 2)
  {
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  else if (symtype == 3)
  {
    speed *= 2;
    speed *= query.aa[0].len;
  }
  else if (symtype == 4)
  {
    speed *= 2;
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  speed /= tip->elapsed;
  tip->speed = speed;

  if (view == 0)
  {
    fprintf(out, "Search started:    %s\n", tip->starttime);
    fprintf(out, "Search completed:  %s\n", tip->endtime);
    fprintf(out, "Elapsed:           %.2fs\n", tip->elapsed);
    fprintf(out, "Speed:             %.3f GCUPS\n", tip->speed / 1e9);
    fprintf(out, "\n");
  }

  free(tip->starttime);
  tip->starttime = 0;
  free(tip->endtime);
  tip->endtime = 0;
}

struct subject_list
{
  char *subject_sequence_masked;
  char *subject_sequence_unmasked;
  Blast_AminoAcidComposition subject_composition;
  Blast_AminoAcidComposition subject_composition_unmasked;
  long subject_length;
} * subject_list;

pthread_mutex_t hitsmutex2 = PTHREAD_MUTEX_INITIALIZER;

void enter_subject(long seqno, char *subject_sequence_masked, char *subject_sequence_unmasked, Blast_AminoAcidComposition subject_composition, Blast_AminoAcidComposition subject_composition_unmasked, long subject_length)
{
  pthread_mutex_lock(&hitsmutex2);

  subject_list[seqno].subject_composition = subject_composition;
  subject_list[seqno].subject_composition_unmasked = subject_composition_unmasked;
  subject_list[seqno].subject_sequence_masked = (char*)calloc(subject_length, sizeof(char));
  memcpy(subject_list[seqno].subject_sequence_masked, subject_sequence_masked, subject_length*sizeof(char));
  subject_list[seqno].subject_sequence_unmasked = (char*)calloc(subject_length, sizeof(char));
  memcpy(subject_list[seqno].subject_sequence_unmasked, subject_sequence_unmasked, subject_length*sizeof(char));
  subject_list[seqno].subject_length = subject_length;

  pthread_mutex_unlock(&hitsmutex2);
}

void get_subject_composition(long seqno, Blast_AminoAcidComposition* subject_composition)
{
  struct subject_list * h = subject_list + seqno;
  *subject_composition = h->subject_composition;
}

void get_subject_composition_unmasked(long seqno, Blast_AminoAcidComposition* subject_composition_unmasked)
{
  struct subject_list * h = subject_list + seqno;
  *subject_composition_unmasked = h->subject_composition_unmasked;
}

void get_subject_length(long seqno, long* subject_length)
{
  struct subject_list * h = subject_list + seqno;
  *subject_length = h->subject_length;
}

char* get_subject_sequence_masked(long seqno)
{
  return subject_list[seqno].subject_sequence_masked;
}

char* get_subject_sequence_unmasked(long seqno)
{
  return subject_list[seqno].subject_sequence_unmasked;
}

void compute_subject_compositions()
{
  struct db_thread_s * dbt = db_thread_create();
  long seqno;
  long dstrand;
  long dframe;
  char *subject_sequence_masked, *subject_sequence_unmasked;
  long subject_length;
  long ntlen;
  Blast_AminoAcidComposition subject_composition, subject_composition_unmasked;
  long db_seqs;

  db_seqs = db_getseqcount();
  subject_list = (struct subject_list *) xmalloc(db_seqs * sizeof(struct subject_list));
  for(long seqno = 0; seqno < db_seqs; seqno++)
  {
    db_mapsequences(dbt, seqno, seqno);
    db_getsequence(dbt, seqno, dstrand, dframe, &subject_sequence_masked, &subject_length, &ntlen, 0);
    subject_length--;
    if (mask)
      subject_sequence_unmasked = (char*)(*dbSequences)[seqno].data();
    else
      subject_sequence_unmasked = subject_sequence_masked;
    Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, (const Uint1*)subject_sequence_masked, subject_length);
    Blast_ReadAaComposition(&subject_composition_unmasked, BLASTAA_SIZE, (const Uint1*)subject_sequence_unmasked, subject_length);
    enter_subject(seqno, subject_sequence_masked, subject_sequence_unmasked, subject_composition, subject_composition_unmasked, subject_length);
  }
}

void work()
{
  args_show();
  hits_init(maxmatches, alignments, minscore, maxscore, minexpect, expect, view==0);

  compute7 = 0;
  compute16 = 0;
  compute32 = 0;
  compute63 = 0;
  rounds7 = 0;
  rounds16 = 0;
  rounds32 = 0;
  rounds63 = 0;

  //  totalhits = 0;

  prepare_search(threads);

  if (view==0)
  {
    fprintf(out, "Searching...");
    fflush(out);
  }

  clock_start(&ti);

#ifdef TIME_PAIRCOUNT
  auto begin1 = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

  run_threads();

#ifdef TIME_PAIRCOUNT
  // time stage 1
  auto end1 = std::chrono::high_resolution_clock::now();
  time_stage1 += (end1 - begin1).count();
#endif //TIME_PAIRCOUNT

#if 1
  if (view == 0)
    fprintf(out, "...............................................done\n\n");
#endif

  clock_stop(&ti);

#if 0
  if (view == 0)
  {
    fprintf(out, "Computed (7bit):   %ld sequences in %ld rounds\n", compute7, rounds7);
    fprintf(out, "Computed (16bit):  %ld sequences in %ld rounds\n", compute16, rounds16);
    fprintf(out, "Computed (63bit):  %ld sequences in %ld rounds\n", compute63, rounds63);
    //    fprintf(out, "Total hits:        %ld\n", totalhits);
    fprintf(out, "\n");
  }
#endif

#if 0
  if (view==0)
  {
    fprintf(out, "Aligning...");
    fflush(out);
  }
#endif

  //  if (view == 0)
  //    clock_start(&ti);

  if (view!=88){
    align_threads();
  }
#ifdef COMPO_ADJUSTMENT
  else {
    if (align_adjusted_version == 1)
      align_adjusted();
    else if (align_adjusted_version == 4)
      align_adjusted4();
    else if (align_adjusted_version == 5)
      align_adjusted5();
    else {
      fprintf(stderr, "Invalid align_adjusted version (%ld)\n", align_adjusted_version);
      exit(1);
    }

//	fprintf(out, "# Fields: Query id, Subject id, masked score, mskAdjScrBL50, umskAdjScrBL62, q. start, q. end, s. start, s. end, %% identity, %% positives, alignment length, mismatches, gap openings\n");
  }
#endif // COMPO_ADJUSTMENT

#if 0
  if (view == 0)
    fprintf(out, "...............................................done\n\n");
#endif

  //  if (view == 0)
  //    clock_stop(&ti);

  hits_show(view, show_gis);
  hits_exit();
}

int main(int argc, char**argv)
{
  // Record start time
  auto start_main = std::chrono::high_resolution_clock::now();

#ifdef MPISWIPE
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS)
    fatal("Unable to initialize MPI.");

  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  MPI_Comm_size(MPI_COMM_WORLD, & size);

  mpirank = rank;
  mpisize = size;

#endif

  cpu_features();

#ifndef MPISWIPE
  if (! cpu_feature_sse2)
    fatal("This program requires a processor with SSE2.\n");
#endif

  args_init(argc,argv);

  db_open(symtype, databasename, taxidfilename);

  volchunks = (long*) xmalloc(db_getvolumecount() * sizeof(long));
  volseqs   = (long*) xmalloc(db_getvolumecount() * sizeof(long));

  if(dump)
  {
    struct db_thread_s * t = db_thread_create();
    long seqcount = db_getseqcount();
    for (long i=0; i < seqcount; i++)
      db_show_fasta(t, i, 0, 0, dump-1);
    db_thread_destruct(t);
  }
  else
  {
    score_matrix_init();

    queryno = 0;

    query_init(queryname, symtype, querystrands);

#ifdef MPISWIPE
    if (!mpirank)
#endif
    {
      hits_show_begin(view);
    }

#ifdef COMPO_ADJUSTMENT
    align_adjusted_init();

#ifdef TIME_PAIRCOUNT
    auto begin2b = std::chrono::high_resolution_clock::now();
#endif //TIME_PAIRCOUNT

    compute_subject_compositions();

#ifdef TIME_PAIRCOUNT
    auto end2b = std::chrono::high_resolution_clock::now();
    time_stage2b += (end2b - begin2b).count();
#endif //TIME_PAIRCOUNT

#endif //COMPO_ADJUSTMENT

    while (query_read())
    {
    //TODO: implement masking
    auto start_work = std::chrono::high_resolution_clock::now();

#ifdef MPISWIPE
      if (rank == 0)
	master(size);
      else
	slave(rank,size);
#else
      work();
#endif

      queryno++;
#ifdef TIME_PAIRCOUNT
      auto end_work = std::chrono::high_resolution_clock::now();
      time_work += (end_work - start_work).count();
#endif //TIME_PAIRCOUNT
    }

#ifdef MPISWIPE
    if (!mpirank)
#endif
    {
      hits_show_end(view);
    }

#ifdef COMPO_ADJUSTMENT
    align_adjusted_done();
#endif

    query_exit();

#ifdef DEBUG
    fprintf(out, "Freeing score matrix\n");
#endif
    score_matrix_free();
  }

  free(volchunks);
  free(volseqs);

#ifdef DEBUG
  fprintf(out, "Closing db\n");
#endif
  db_close();

#ifdef MPISWIPE
#ifdef DEBUG
  fprintf(out, "Finalizing (rank %d).\n", rank);
#endif
  MPI_Finalize();
#ifdef DEBUG
  fprintf(out, "Finalized (rank %d).\n", rank);
#endif
#endif

 // fprintf(stderr, "matrix_min = %d, matrix_max = %d\n", matrix_min, matrix_max);

  if (outfile)
    fclose(out);

#ifdef TIME_PAIRCOUNT
  // print time, counts
  auto end_main = std::chrono::high_resolution_clock::now();
  std::cout << "time_main: " << (end_main - start_main).count() * 1e-9 << " s\n";
  std::cout << "time_work: " << time_work * 1e-9 << " s\n";
  std::cout << "time_stage1: " << time_stage1 * 1e-9 << " s\n";
  std::cout << "time_stage2: " << time_stage2 * 1e-9 << " s\n";
  std::cout << "time_stage2b: " << time_stage2b * 1e-9 << " s\n";
  std::cout << "time_stage2c: " << time_stage2c * 1e-9 << " s\n";
  std::cout << "time_writing: " << time_writing * 1e-9 << " s\n";
  std::cout << "time_writing_matrix: " << time_writing_matrix * 1e-9 << " s\n";
  std::cout << "time_writing_subject: " << time_writing_subject * 1e-9 << " s\n";
  std::cout << "time_writing_query: " << time_writing_query * 1e-9 << " s\n";
  std::cout << "time_stage2d: " << time_stage2d * 1e-9 << " s\n";
  std::cout << "time_stage2_init: " << time_stage2_init * 1e-9 << " s\n";
  std::cout << "time_stage2_sw: " << time_stage2_sw * 1e-9 << " s\n";
  std::cout << "time_stage2e: " << time_stage2e * 1e-9 << " s\n";
  std::cout << "time_stage2f: " << time_stage2f * 1e-9 << " s\n";
  std::cout << "time_stage3: " << time_stage3 * 1e-9 << " s\n";
  std::cout << "pairs1: " << pairs1 << "\n";
  std::cout << "pairs2: " << pairs2 << "\n";
  std::cout << "pairs2b: " << pairs2b << "\n";
  std::cout << "pairs2c: " << pairs2c << "\n";
  std::cout << "pairs3: " << pairs3 << "\n";
  std::cout << "pairs3b: " << pairs3b << "\n";
  std::cout << "pairs3c: " << pairs3c << "\n";
  std::cout << "pairs3_rev: " << pairs3_rev << "\n";
  std::cout << "mar2_0: " << mar2_0 << "\n";
  std::cout << "mar2_4: " << mar2_4 << "\n";
  std::cout << "mar2_other: " << mar2_other << "\n";
  std::cout << "mar3_0: " << mar3_0 << "\n";
  std::cout << "mar3_4: " << mar3_4 << "\n";
  std::cout << "mar3_other: " << mar3_other << "\n";
  std::cout << "skipped_stage2: " << skipped_stage2 << "\n";
#endif //TIME_PAIRCOUNT
}