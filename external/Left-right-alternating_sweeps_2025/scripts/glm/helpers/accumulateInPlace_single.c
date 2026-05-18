#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <stdint.h>

/* Simple helper MEX function for performing C-style "in-place" 
 * accumulation of values from a matrix, with row-wise weights.
 */

void mexFunction(int nOUT, mxArray *pOUT[], int nINP, const mxArray *pINP[]) {
    
    //mexPrintf("nINP = %u\n", nINP);
    if (nINP != 7) {
        mexErrMsgTxt("Usage: accumulateInPlace(YIn, YOut colIndsIn, colIndsOut, rowIndsIn, weights, rezero)");
    }
    if (nOUT != 0) {
        mexErrMsgTxt("This mex function returns no outputs.");
    }

    // usage:
    // accumulateInPlace(YIn, YOut colIndsIn, colIndsOut, rowIndsIn, weights)
    float*      YIn         = mxGetPr(pINP[0]);
    float*      YOut        = mxGetPr(pINP[1]);
    uint32_t*   colIndsIn   = (uint32_t*) mxGetData(pINP[2]);
    uint32_t*   colIndsOut  = (uint32_t*) mxGetData(pINP[3]);
    uint32_t*   rowIndsIn   = (uint32_t*) mxGetData(pINP[4]);
    float*      weights     = mxGetPr(pINP[5]);
    int         rezero      = (int) mxGetScalar(pINP[6]);
    
    int nRowsIn = mxGetM(pINP[0]);
    int nRowsOut = mxGetM(pINP[1]);
    int nCols = mxGetN(pINP[2]);
    
    int colIdxOut, colIdxIn;
    uint32_t idxOut0, idxIn0, idxOut, idxIn;
    
    // Loop through columns
    for (uint32_t c=0; c<nCols; c++) {
        idxOut0 = colIndsOut[c]*nRowsOut;
        idxIn0 = colIndsIn[c]*nRowsIn;
        for (uint32_t r=0; r<nRowsOut; r++) {
            idxOut = idxOut0 + r;
            idxIn = idxIn0 + rowIndsIn[r];
            if (rezero) {
                YOut[idxOut] = YIn[idxIn]*weights[r];
            } else {
                YOut[idxOut] += YIn[idxIn]*weights[r];
            }
        }
    }
}