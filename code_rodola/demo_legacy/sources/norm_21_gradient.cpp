/*
 *  
 * mex norm_21_gradient.cpp COMPFLAGS="/openmp $COMPFLAGS"
 */
 
#include <limits>
#include <iostream>
#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{ 
    if (nrhs != 3)
		mexErrMsgTxt("Usage: [S] = norm_21_gradient(C,A,B) where ||CA-B||_{2,1}.");
    if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 

	// input
	const double *C = (double*)mxGetPr(prhs[0]);	
    const double *A = (double*)mxGetPr(prhs[1]);	
    const double *B = (double*)mxGetPr(prhs[2]);	
    
    //matrix MxN
    const int m = int( mxGetM(prhs[1]));
    const int n = int( mxGetN(prhs[1]));
    
    
    // output
	plhs[0] = mxCreateDoubleMatrix( (mwSize)m, (mwSize)m, mxREAL);
	double* gC = mxGetPr(plhs[0]);
  
    std::vector<double> sums_i(n);
    
    #pragma omp parallel for
    for(int j=0; j<n; ++j)
    {
        double sum_i = 0;
        for(int i=0; i<m; ++i)
        {
            double sum_r = 0;
            for(int r=0; r<m; ++r)
                sum_r += C[i + r*m]*A[r + j*m];
            sum_r -= B[i + j*m];
            sum_i += sum_r*sum_r;
        }
        sums_i[j] = 1/std::sqrt(sum_i);
    }
    
	// run	
    #pragma omp parallel for
    for(int p=0; p<m; ++p)
    for(int q=0; q<m; ++q)
    {             
        double sum_j=0;
        for(int j=0; j<n; ++j)
        {
            double sum_r = 0;
            for(int r=0; r<m; ++r)
               sum_r += C[p + r*m]*A[r + j*m];
            sum_r -= B[p + j*m];
            sum_j += A[q + j*m]*sums_i[j]*sum_r;
        }
                            
         gC[p + q*m] = sum_j;
    }
    
    //delete[] sums1;
    return;
}
