/* extractGroupedDataMex.cpp
 *
 * C = EXTRACTGROUPEDDATA(GROUPS, X1, X2, ... , XN)
 *
 * Efficient extraction of group-indexed data from a numeric matrix into a
 * cell array of numeric vectors, with each cell element corresponding to
 * one group.
 *
 * Must be compiled with the R2018a MEX API by using the "-R2018a" flag:
 * mex -R2018a extractGroupedDataMex.cpp
 *
 * R.J. Gardner 2019/05/17
 *
 */

#include "mex.h"
#include <math.h>
#include "matrix.h"

void *getData(const mxArray *arr);

// template <class Tg> void extractSub1(mxArray **pLHS[], int nINP, const mxArray *pRHS[]);

template <class Tg> void extractSub(mxArray **ppLHS[], int nINP, const mxArray *pRHS[]);

template <class Tx, class Tg> void extractOneArray(mxArray **ppLHS[], int xIndex, const mxArray *pRHS[], mwSize nGroups, mwSize *groupHist[], Tg *groups[]);

void mexFunction(int nOUT, mxArray *pLHS[], int nINP, const mxArray *pRHS[]) {
    
    // Parse inputs
    if (nINP < 2 || nOUT != 1)
        mexErrMsgTxt("Syntax: C = extractGroupedData(groups, X1, ... , XN).");

    mwSize    ng  = mxGetNumberOfElements(pRHS[0]);
    for (int i=1; i<nINP; i++) {
//         mwSize    nx = mxGetNumberOfElements(pRHS[i]);
        mwSize nx = mxGetM(pRHS[i]);
        if (nx != ng) {
            mexErrMsgTxt("Every input X array must have same number of rows as the number of elements in GROUPS.");
        }
    }
    
    mxClassID gType = mxGetClassID(pRHS[0]);
    switch (gType) {
        case mxDOUBLE_CLASS:
            extractSub<mxDouble>(&pLHS, nINP, pRHS);
            break;
        case mxUINT8_CLASS:
            extractSub<mxUint8>(&pLHS, nINP, pRHS);
            break;
        case mxUINT16_CLASS:
            extractSub<mxUint16>(&pLHS, nINP, pRHS);
            break;
        case mxUINT32_CLASS:
            extractSub<mxUint32>(&pLHS, nINP, pRHS);
            break;
        case mxUINT64_CLASS:
            extractSub<mxUint64>(&pLHS, nINP, pRHS);
            break;
        default:
            mexErrMsgTxt("Input 'groups' is of an invalid type. Valid types are double, or any unsigned integer type.");
    }
    
}

void *getData(const mxArray *arr) {
    mxClassID cls = mxGetClassID(arr);
    switch (cls) {
        case mxDOUBLE_CLASS:
            return mxGetDoubles(arr);
        case mxSINGLE_CLASS:
            return mxGetSingles(arr);
        case mxUINT8_CLASS:
            return mxGetUint8s(arr);
        case mxUINT16_CLASS:
            return mxGetUint16s(arr);
        case mxUINT32_CLASS:
            return mxGetUint32s(arr);
        case mxUINT64_CLASS:
            return mxGetUint64s(arr);
            
        case mxINT8_CLASS:
            return mxGetInt8s(arr);
        case mxINT16_CLASS:
            return mxGetInt16s(arr);
        case mxINT32_CLASS:
            return mxGetInt32s(arr);
        case mxINT64_CLASS:
            return mxGetInt64s(arr);
            
        case mxLOGICAL_CLASS:
            return mxGetLogicals(arr);
            
        default:
            mexErrMsgTxt("Invalid type.");
    }
}

template <class Tg> void extractSub(mxArray **ppLHS[], int nINP, const mxArray *pRHS[]) {
    // Calculate group index histogram and pass to x-type-specific
    // subroutine for each input array
    
    int nXArrays = nINP-1;
//     mwSize numelArray = mxGetNumberOfElements(pRHS[0]);
    mwSize nRows = mxGetNumberOfElements(pRHS[0]);
    Tg *groups = (Tg *) getData(pRHS[0]);
//     mwSize nGroups = getNGroups<Tg>(groups, numelArray);
    mwSize nGroups = getNGroups<Tg>(groups, nRows);
    
    // Calculate group histogram
    mwSize *groupHist = (mwSize *) mxCalloc(nGroups+1, sizeof(mwSize));
    mwIndex ig, iig;
    for (mwIndex i=0; i<nRows; i++) {
        ig = groups[i];
        groupHist[ig] += 1;
    }
    //mexPrintf("Calculated group histogram\n");
    
    // Allocate the output cell array
    mwSize dims[2] = {nGroups, nXArrays};
    //mexPrintf("Creating cell array with length %u...\n", nGroups);
    mxArray *pout = mxCreateCellArray(2, dims);
    (*ppLHS)[0] = pout;
    
    for (int ix=0; ix<nXArrays; ix++) {
        mxClassID xtype = mxGetClassID(pRHS[ix+1]);
        switch (xtype) {
            case mxDOUBLE_CLASS:
                extractOneArray<mxDouble, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxSINGLE_CLASS:
                extractOneArray<mxSingle, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxUINT8_CLASS:
                extractOneArray<mxUint8, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxUINT16_CLASS:
                extractOneArray<mxUint16, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxUINT32_CLASS:
                extractOneArray<mxUint32, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxUINT64_CLASS:
                extractOneArray<mxUint64, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
                
            case mxINT8_CLASS:
                extractOneArray<mxInt8, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxINT16_CLASS:
                extractOneArray<mxInt16, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxINT32_CLASS:
                extractOneArray<mxInt32, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
            case mxINT64_CLASS:
                extractOneArray<mxInt64, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
                
            case mxLOGICAL_CLASS:
                extractOneArray<mxLogical, Tg>(&pout, ix, pRHS, nGroups, &groupHist, &groups);
                break;
        }
    }
    
}

template <class Tx, class Tg> void extractOneArray(mxArray **ppout, int xIndex, const mxArray *pRHS[], mwSize nGroups, mwSize *groupHist[], Tg *groups[]) {
    
    //mexPrintf("Extracting input index %u...\n", xIndex);
    
    const mxArray *px = pRHS[xIndex + 1];
    Tx **cellElements;
    mxArray *arr;
    mxClassID xType = mxGetClassID(px);
    
    /*
     * Allocate a pointer for the numeric data in each cell element. We allocate for
     * nGroups+1, so that we can use the MATLAB 1-based indices unaltered.
     * We also create an index counter for each cell element, to keep track
     * of the current index when distributing the values of x.
     */
    
    //mexPrintf("Allocating cell element arrays...\n", nGroups);
    
    mwSize nrows = mxGetNumberOfElements(pRHS[0]);
    mwSize ncols = mxGetNumberOfElements(px) / nrows;
    
    cellElements = (Tx **) (mxCalloc(nGroups+1, sizeof(Tx **)));
    
    mwIndex cellIdx = xIndex*nGroups;
    for (mwIndex i=0; i<nGroups; i++) {
        // For the current cell, allocate the appropriate mxArray
//         arr = mxCreateNumericMatrix((*groupHist)[i+1], 1, xType, mxREAL);
        arr = mxCreateNumericMatrix((*groupHist)[i+1], ncols, xType, mxREAL);
        mxSetCell(*ppout, cellIdx++, arr);
        cellElements[i+1] = (Tx *) getData(arr);
        //mexPrintf("Allocated %u elements to group %u\n", (*groupHist)[i+1], i);
    }
    
    // Distribute x values into cells according to indices
    mwIndex *groupCounters = (mwIndex *) mxCalloc(nGroups+1, sizeof(mwIndex));
    Tx *x = reinterpret_cast<Tx *>(getData(px));
    
    mwIndex ig, iig, icolx;
    
    for (mwIndex i=0; i<nrows; i++) {
        ig = (*groups)[i];
        if (ig != 0) {
            iig = groupCounters[ig]++;
            for (mwIndex c=0; c<ncols; c++) {
                icolx = c * nrows;
                cellElements[ig][iig] = x[i+icolx];
                iig += (*groupHist)[ig];
            }
        }
    }
}

template <class T> mwSize getNGroups(T *groups, mwSize n) {
    mwSize nGroups = 0;
    double ig;
    for (mwIndex i=0; i<n; i++) {
        ig = groups[i];
        if (ig > nGroups)
            nGroups = ig;
        // Check validity of all group identifiers. Groups are specified by
        // 1-based indexing, and must be integer-valued. Zero is permitted
        // as an indicator of excluded data.
        if (ig < 0 || fmod(ig, 1.0) != 0.0 )
            mexErrMsgTxt("Input 'groups' must contain only non-negative integers.");
    }
    return nGroups;
}
