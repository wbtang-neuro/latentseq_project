#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <stdint.h>

/* 
 * In-place matrix-vector multiplication using subset of matrix columns
 */

void mexFunction(int nOUT, mxArray *pOUT[], int nINP, const mxArray *pINP[]) {
    
    /* usage:
     * mvmInPlace(X, y, v, ciX)
     * X: input matrix (dims [m, n])
     * y: output vector (numel = m)
     * v: vector to multiply with X
     * ciX: indices of subset of X columns to multiply with v
     */
    
    
    float*     X           = mxGetPr(pINP[0]);
    float*     y           = mxGetPr(pINP[1]);
    double*    v           = mxGetPr(pINP[2]);
    uint16_t*   colIndsIn   = (uint16_t*) mxGetData(pINP[3]);
    
    size_t nCols = mxGetNumberOfElements(pINP[2]);
    size_t nRows = mxGetM(pINP[0]);
    size_t idxIn0, idxIn;
    
    // Loop through columns
    for (size_t c=0; c<nCols; c++) {
        idxIn0 = colIndsIn[c]*nRows;
        // Loop through rows
        for (size_t r=0; r<nRows; r++) {
            idxIn = idxIn0 + r;
            y[r] += X[idxIn]*v[c];
        }
    }
}