#include <cstdio>
#include "edlib.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector edlibNW(StringVector query, NumericVector queryLength,
                         StringVector target, NumericVector targetLength) {
    int n = query.size();
    NumericVector scores(n);
    
    for(int i = 0; i < n; ++i) {
        EdlibAlignResult result = edlibAlign(query[i], queryLength[i],
                                             target[i], targetLength[i],
                                             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0)
                                             );
        scores[i] = result.editDistance;
        edlibFreeAlignResult(result);
    }
    
    return scores;
    
}
