/*
    compositional adjusted scoring matrix aligment
*/

#include "swipe.h"
#include <algo/blast/composition_adjustment/smith_waterman.h>
#include <algo/blast/composition_adjustment/composition_constants.h>
#include <algo/blast/composition_adjustment/redo_alignment.h>
#include <algo/blast/api/blast_setup.hpp>
#include <algo/blast/core/matrix_freq_ratios.h>
//#include "blast_kappa.h"

/** pseudocounts for relative-entropy-based score matrix adjustment */
int kReMatrixAdjustmentPseudocounts = 20;

/**
 * A callback routine: compute lambda for the given score
 * probabilities.
 * (@sa calc_lambda_type).
 */
static double s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
{

    int i;                 /* loop index */
    int score_range;       /* range of possible scores */
    double avg;            /* expected score of aligning two characters */
    Blast_ScoreFreq freq;  /* score frequency data */

    score_range = max_score - min_score + 1;
    avg = 0.0;
    for (i = 0;  i < score_range;  i++) {
        avg += (min_score + i) * probs[i];
    }
    freq.score_min = min_score;
    freq.score_max = max_score;
    freq.obs_min = min_score;
    freq.obs_max = max_score;
    freq.sprob0 = probs;
    freq.sprob = &probs[-min_score];
    freq.score_avg = avg;

    return Blast_KarlinLambdaNR(&freq, lambda0);
}

void compo_init(const char *matrixName, BlastScoreBlk **bsb, Blast_MatrixInfo **matrixInfo) {
//    const char *matrixName = "BLOSUM50";

    BlastScoreBlk *sbp = BlastScoreBlkNew(BLASTAA_SEQ_CODE, 1);
    sbp->name = strdup(matrixName);
    Blast_ScoreBlkMatrixFill(sbp, &ncbi::blast::BlastFindMatrixPath);
    Blast_ScoreBlkKbpIdealCalc(sbp);
    Blast_MatrixInfo *scaledMatrixInfo = Blast_MatrixInfoNew(BLASTAA_SIZE, BLASTAA_SIZE, 0);
    scaledMatrixInfo->matrixName = strdup(matrixName);
    scaledMatrixInfo->ungappedLambda = sbp->kbp_ideal->Lambda / scaling_factor;

    /* Frequency ratios for the matrix */
    SFreqRatios * stdFreqRatios = NULL;

    stdFreqRatios = _PSIMatrixFrequencyRatiosNew(scaledMatrixInfo->matrixName);
    for (int i = 0;  i < BLASTAA_SIZE;  i++) {
      for (int j = 0;  j < BLASTAA_SIZE;  j++) {
          scaledMatrixInfo->startFreqRatios[i][j] = stdFreqRatios->data[i][j];
      }
    }
    _PSIMatrixFrequencyRatiosFree(stdFreqRatios);

    Blast_Int4MatrixFromFreq(scaledMatrixInfo->startMatrix, scaledMatrixInfo->cols,
                           scaledMatrixInfo->startFreqRatios, scaledMatrixInfo->ungappedLambda);

    *bsb = sbp;
    *matrixInfo = scaledMatrixInfo;
}

void compo_done(BlastScoreBlk **sbp, Blast_MatrixInfo **scaledMatrixInfo) {
    Blast_MatrixInfoFree(scaledMatrixInfo);
    BlastScoreBlkFree(*sbp);
}

int compo_align(long *score_out, Blast_CompositionWorkspace * NRrecord, BlastScoreBlk *sbp, Blast_MatrixInfo *scaledMatrixInfo, const Uint1 *data, int nData, long gapopen, long gapextend,
                int *matchStart, int *queryStart, int *matchEnd, int *queryEnd) {

    int score = -1;
    int err = 0;
    *score_out = -1;

    /* adjust_search_failed is true only if Blast_AdjustScores
     * is called and returns a nonzero value */
    int adjust_search_failed = FALSE;


	Blast_ForbiddenRanges forbidden = {0,};
	forbidden.isEmpty = TRUE;

    double pvalueForThisPair = (-1); /* p-value for this match for composition; -1 == no adjustment*/
    double LambdaRatio; /*lambda ratio*/
    /* which test function do we use to see if a composition-adjusted
       p-value is desired; value needs to be passed in eventually*/
    int compositionTestIndex = 0;
    /* which mode of composition adjustment is actually used? */
    EMatrixAdjustRule matrix_adjust_rule = eDontAdjustMatrix;

    Blast_AminoAcidComposition subject_composition;
    Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, data, nData);
//    Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, &data[res.matchSeqStart], res.matchSeqEnd-res.matchSeqStart);

    adjust_search_failed =
      Blast_AdjustScores(sbp->matrix->data,
                         (mask ? &query.composition_unmasked : &query.composition), query.aa[0].len,
                         &subject_composition, nData,
                         scaledMatrixInfo, eCompositionMatrixAdjust,
                         kReMatrixAdjustmentPseudocounts, NRrecord,
                         &matrix_adjust_rule, &s_CalcLambda,
                         &pvalueForThisPair,
                         compositionTestIndex,
                         &LambdaRatio);

    if (adjust_search_failed < 0)
        return adjust_search_failed;  // Score adjustment error

    const char *seq = (mask ? query.aa[0].seq_unmasked : query.aa[0].seq);
    err = Blast_SmithWatermanScoreOnly( &score, matchEnd, queryEnd,
                                        data, nData, (const Uint1*)seq, query.aa[0].len, sbp->matrix->data, gapopen, gapextend, false, &forbidden );
    if (err)
        return err;
    
    *score_out = score;

/*    if (matchStart) {
        int updatedScore;            // score found by the SW algorithm run in reverse
        err = Blast_SmithWatermanFindStart(&updatedScore, matchStart, queryStart, data, nData, (const Uint1*)seq, sbp->matrix->data, gapopen, gapextend, *matchEnd, *queryEnd, score, false, &forbidden);
        // the redone alignment
        BlastCompo_Alignment * newAlign;
        BlastCompo_GappingParams gapping_params = {0};

        gapping_params.gap_open = gapopen;
        gapping_params.gap_extend = gapextend;

        err = s_NewAlignmentUsingXdrop(&newAlign, queryEnd, matchEnd,
                                        *queryStart, *matchStart, score,
                                        &query, &window->query_range,
                                        0,
                                        &subject, &window->subject_range,
                                        0,
                                        gapping_params, matrix_adjust_rule);
    }
*/

    return err;
}
