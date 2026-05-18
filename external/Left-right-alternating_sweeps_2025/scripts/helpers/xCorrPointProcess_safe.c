/* xCorrPointProcess.c
 *
 * Calculates crosscorrelation histogram for two point process time series.
 * Usage: [C, B] = xCorrPointProcess(t1, t2, binSize, nBins, normalize)
 *
 * INPUTS
 * t1, t2    - double column vectors containing sorted spike times
 * binSize   - scalar double specifying the bin size to use
 * bBins     - scalar double specifying the number of bins (should be odd)
 * normalize - use "1" to normalize output by number f t1 events; use "0"
 *             to output raw counts
 *
 * OUTPUTS
 * C         - crosscorrelation histogram
 * B         - crosscorrelation bin centre values (N.B. MClust uses bin edges!)
 *
 * The algorithm differs from that of the MClust CrossCorr function in that
 * it does not iterate over bins centred around events in t1.  Rather, it 
 * iterates over events in t2 that fall within the current t1-centered time window, 
 * assigning them to the appropriate bin based on the time difference from 
 * the current t1 event. For cases where concidence of t1 and t2 is sparse
 * (i.e. for a single t1 event, most of the perievent bins are empty), 
 * this is much more efficient than iterating over bins. Computation time
 * increases with the width of the time window (not the number of bins).
 *
 * The MClust algorithm may be marginally faster in cases where bin 
 * occupancy is very dense.
 *
 * Rich 2017-01-29
 *
 * Updated 2024-08-09 to fix inconsistent behavior between compilers
 *
 */

#include "mex.h"
#include <matrix.h>

#define min(a, b) a<b ? a : b
#define TOLERANCE 0.00001 // to avoid rounding errors

int myRound(double num) {
    if (num >= 0) {
        return (int)(num + 0.5);
    } else {
        return (int)(num - 0.5);
    }
}

void mexFunction(int nOUT, mxArray *pOUT[], int nINP, const mxArray *pINP[]) {
    
    double *t1;
    double *t2;
    double binsize;
    double *C, *B;
    double windowSize;
    
    int norm, ignoreExact;
    mwSize nbins, nt1, nt2, iMidBin;
    mwIndex i1, i2, ii2;
    int t2endReached = 0;
    
    /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
    if (nINP < 4)
        mexErrMsgTxt("Call with at least t1, t2, binsize and nbins as inputs");
    if (nOUT < 1 || nOUT > 2)
        mexErrMsgTxt("Requires one or two outputs.");
    
    /* check validity of inputs */
    if (mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
        mexErrMsgTxt("t1 must be a row or column vector");
    if (mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
        mexErrMsgTxt("t2 must be a row or column vector");
    if (mxGetM(pINP[2]) * mxGetN(pINP[2]) != 1)
        mexErrMsgTxt("binsize must be scalar");
    if (mxGetM(pINP[3]) * mxGetN(pINP[3]) != 1)
        mexErrMsgTxt("nbins must be scalar");
    if (nINP==5 && mxGetM(pINP[4]) * mxGetN(pINP[4]) != 1)
        mexErrMsgTxt("norm must be scalar");
    if (nINP>5)
        mexErrMsgTxt("Usage: xCorrPointProcess(t1, t2, nBins, binSize, normalize)");
    
    /* unpack inputs */
    nt1 = mxGetM(pINP[0]) * mxGetN(pINP[0]);
    t1 = mxGetPr(pINP[0]);
    nt2 = mxGetM(pINP[1]) * mxGetN(pINP[1]);
    t2 = mxGetPr(pINP[1]);
    binsize = mxGetScalar(pINP[2]);
    nbins = (int)mxGetScalar(pINP[3]);
    
    // Determine normalization
    if (nINP==5) {
        norm = mxGetScalar(pINP[4]);
        if (norm < 0 || norm > 1) {
            mexErrMsgTxt("norm must be either 0 or 1");
        }
    }
    else {
        norm = 1;
    }
    
    // Determine value of "ignoreExact"
    if (nINP==6) {
        ignoreExact = (int)mxGetScalar(pINP[5]);
        if (ignoreExact < 0 || ignoreExact > 1) {
            mexErrMsgTxt("ignoreExact must be either 0 or 1");
        }
    }
    else {
        ignoreExact = 0;
    }
    
    /* we want nbins to be odd */
    if ((nbins / 2) * 2 == nbins)
        nbins++;
    
    iMidBin = (nbins-1)/2;
    
    pOUT[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
    C = mxGetPr(pOUT[0]);
    
    if(nOUT == 2) {
        double m;
        
        pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
        B =  mxGetPr(pOUT[1]);
        m = - binsize * (nbins / 2);
        
        for(mwIndex j = 0; j < nbins; j++)
            B[j] = m + j * binsize;
    }
    
    windowSize = ((double)(nbins)/2) * binsize;
    
    if(nt2 > 0) { // Handle possibility of t2 being empty
        
        i2 = 0;
        
        // For each time in t1
        for (i1 = 0; i1 < nt1; i1++) {
            
            // Scan along t2 until reaching a time inside the window
            while (t2[i2] <= t1[i1]-windowSize) {
                // If end of t2 reached
                if (i2 == nt2) {
                    t2endReached = 1;
                    break;
                }
                i2++;
            }
            
            // If end of t2 reached, we're done
            if (t2endReached) break;
            
            // Initialize the t2 moving cursor position
            ii2 = i2;
            
            // Scan along from the current time in t2 until reaching a time
            // no longer in the window.  Assign each time to the relevant
            // bin, according to the time difference between t1 and t2
            mwIndex idx;
            
            // Iterate through events in t2 until reaching the end of t2 or 
            // the right-hand side of the time window
            while ( (ii2<nt2) && (t2[ii2]) <= t1[i1]+windowSize) {
                
                // For the current t2-t1 lag, calculate the appropriate lag bin index
                idx = myRound((t2[ii2]-t1[i1]) / binsize) + iMidBin;
                idx = min(idx, nbins-1);
                C[idx]++;  // Increment bin count

                // Move to next t2 event
                ii2++;
                
            }
        }
    }
    
    if (norm) {
        for(mwIndex j = 0; j < nbins; j++) {
            C[j] /= nt1 * binsize;
        }
    }
    
}