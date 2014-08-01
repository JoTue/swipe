#include <cstdio>
#include <vector>

using namespace std;

#define NA 123
#define EL 125
#define ES 126
#define AAMASK 127

/*       0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15	*/
/* 32	    !  "  #  $  %  &  '  (  )  *  +  ,  -  .  / 	*/
/* 48	 0  1  2  3  4  5  6  7  8  9  :  ;  <  =  >  ? 	*/
/* 64	 @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O		*/
/* 80	 P  Q  R  S  T  U  V  W  X  Y  Z  [  \  ]  ^  _		*/
/* 96	 `  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o 	*/
/*112	 p  q  r  s  t  u  v  w  x  y  z  {  |  }  ~  ^?	*/

int ncbistdaa[128]={
        EL,NA,NA,NA,NA,NA,NA,NA,NA,NA,EL,NA,NA,EL,NA,NA,	/* 15 */
        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,	/* 31 */
        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,25,NA,NA,NA,NA,NA,	/* 47 */
        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,	/* 63 */
        NA, 1, 2, 3, 4, 5, 6, 7, 8, 9,27,10,11,12,13,26,	/* 79 */
        14,15,16,17,18,24,19,20,21,22,23,NA,NA,NA,NA,NA,	/* 95 */
        NA, 1, 2, 3, 4, 5, 6, 7, 8, 9,27,10,11,12,13,26,	/*111 */
        14,15,16,17,18,24,19,20,21,22,23,NA,NA,NA,NA,NA};	/*127 */

/** Reads sequences from fasta file.
 * @param [in] file File pointer to database.
 * @param [out] seqs Sequences will be stored here, each sequence as vector of indexes from alphabet.
 * @return true if reached end of file, otherwise false.
 */
bool readFastaSequences(const char* dbFilePath, vector< vector<unsigned char> >* seqs) {

    FILE* file = fopen(dbFilePath, "r");
    if (file == 0) {
        printf("Error: There is no file with name %s\n", dbFilePath);
        return 1;
    }
    
    seqs->clear();

    long numResiduesRead = 0;
    bool inHeader = false;
    bool inSequence = false;
    int buffSize = 4096;
    unsigned char buffer[buffSize];
    while (!feof(file)) {
        int read = fread(buffer, sizeof(char), buffSize, file);
        for (int i = 0; i < read; ++i) {
            unsigned char c = buffer[i];
            if (inHeader) { // I do nothing if in header
                if (c == '\n')
                    inHeader = false;
            } else {
                if (c == '>') {
                    inHeader = true;
                    inSequence = false;
                } else {
                    if (c == '\r' || c == '\n')
                        continue;
                    // If starting new sequence, initialize it.
                    if (inSequence == false) {
                        if (seqs->size() > 0) {
                            numResiduesRead += seqs->back().size();
                        }
                        inSequence = true;
                        seqs->push_back(vector<unsigned char>());
                    }

                    seqs->back().push_back(ncbistdaa[c & AAMASK]);
                }
            }
        }
    }
    fclose(file);
    return true;
}
