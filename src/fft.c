/****************************************************************************
 * Copyright (C) 2016 by Ludovic Delchambre                                 *
 *                                                                          *
 * This file is part of wcorrQRL.                                           *
 *                                                                          *
 *   wcorrQRL is free software: you can redistribute it and/or modify       *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation, version 3 of the License.                *
 *                                                                          *
 *   wcorrQRL is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with wcorrQRL.  If not, see <http://www.gnu.org/licenses/>.        *
 ****************************************************************************/
 
/**
* @file fft.c
* @brief Implementation of FFTW3 interface functions.
* @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
* @version 1.0
* @date 11 March 2016
*/

#include <fftw3.h>
#include <fft.h>


unsigned long fftr_size(const unsigned long nvar, const unsigned long nobs) {
  return nobs*((nvar & ~1UL) + 2); // nobs*(2*floor(nvar/2)+2)
}

void *fft_alloc(const unsigned long size) {
  return fftw_malloc(size);
}

void fft_free(void *ptr) {
  fftw_free(ptr);
}

double *fftr(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs) {
  const unsigned long nvarout = (nvar>>1)+1; // Since FFT of real data are symmetric, we only store nvar/2+1 fftw_complex elements as output (symmetric part is not duplicated)
  const int n = nvar; // Store nvar within an int since we will pass it as a pointer
  
  // Create the FFTW plan 
  fftw_plan plan = fftw_plan_many_dft_r2c(1, // [int rank]     Rank 1 DFT
                                          &n, // [const int *n] Number of variables within input array
                                          nobs, // [int howmany] Number of observations (number of DFT to perform)
                                          (double *) X, // [double *in] Input array is X (it is cast to non-const but will not be modified since FFTW_DESTROY_INPUT is not set)
                                          NULL, // [const int *inembed] Distance between each rank in input array (Not used since rank=1)
                                          1, // [int istride] Distance between successive variables in input array (in unit of double)
                                          n, // [int idist]   Distance between 2 observations in input array (in unit of double)
                                          (fftw_complex *) Y, // [fftw_complex *out] Output array is Y
                                          NULL, // [const int *onembed] Distance between each rank in output array (Not used since rank=1)
                                          1, // [int ostride] Distance between successive variables in output array (in unit of fftw_complex)
                                          nvarout, // [int odist] Distance between 2 observations in output array (in unit of fftw_complex)
                                          FFTW_ESTIMATE // [unsigned flags] Quickly choose a plan without performing full benchmarks (maybe sub-optimal but take less time)
                                          );
                                          
  // If plan building fails, quit
  if(!plan)
    return NULL;
  
  // Execute FFTW plan
  fftw_execute(plan);
      
  
  // Destroy FFTW plan
  fftw_destroy_plan(plan);
  
  return Y;
}

double *ifftr(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs) {
  const unsigned long nvarout = (nvar>>1)+1; // Since FFT of real data are symmetric, we only store nvar/2+1 fftw_complex elements as output (symmetric part is not duplicated)
  const int n = nvar;
  unsigned long i,nelem;
  
  // Create the FFTW plan 
  fftw_plan plan = fftw_plan_many_dft_c2r(1, // [int rank] Rank 1 DFT
                                          &n, // [const int *n] Number of variables within input array
                                          nobs, // [int howmany] Number of observations (number of DFT to perform)
                                          (fftw_complex *) X, // [fftw_complex *in] Input array is X (it is cast to non-const but will not be modified since FFTW_DESTROY_INPUT is not set)
                                          NULL, // [const int *inembed] Distance between each rank in input array (Not used since rank=1)
                                          1, // [int istride] Distance between successive variables in input array (in unit of fftw_complex)
                                          nvarout, // [int idist] Distance between 2 observations in input array (in unit of fftw_complex)
                                          Y, // [double *out] Output array is Y
                                          NULL, // [const int *onembed] Distance between each rank in output array (Not used since rank=1)
                                          1, // [int ostride] Distance between successive variables in output array (in unit of double)
                                          n, // [int odist] Distance between 2 observations in output array (in unit of double)
                                          FFTW_ESTIMATE // [unsigned flags] Quickly choose a plan without performing full benchmarks (maybe sub-optimal but take less time)
                                          );
                                          
  // If plan building fails, quit
  if(!plan)
    return NULL;
  
  // Execute FFTW plan
  fftw_execute(plan);

  // Normalize result by nvar (FFTW compute unormalized transform)
  nelem = nvar * nobs;
  for(i = 0; i < nelem; i++)
    Y[i] /= nvar;
      
  
  // Destroy FFTW plan
  fftw_destroy_plan(plan);
  
  return Y;
}

/**
 * @brief Complex DFT and inverse DFT calculation wrapper
 *
 * Compute the DFT or inverse DFT of the complex matrix X of size (nvar x nobs) 
 * into the complex matrix Y of size (nvar x nobs) depending on the sign of sign
 * parameter.
 * 
 * @param[out] Y Output matrix of complex numbers (double[2]) of size (nvar x nobs)
 * @param[in] X Input matrix of complex numbers (double[2]) of size (nvar x nobs) [can be the same as Y]
 * @param[in] nvar Number of variables (rows) within X and Y
 * @param[in] nobs Number of observations (columns) within X and Y
 * @param[in] sign If -1 computes the DFT, if 1 computes the inverse DFT of X
 *
 * @return Pointer to Y or NULL if DFT fails
 */
double *_fft(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs, int sign) {
  const int n = nvar;
  unsigned long i,nelem;
  
  fftw_plan plan = fftw_plan_many_dft(1, // [int rank] Rank 1 DFT
  																	&n, // [const int *n] Number of variables within input array
  																	nobs, // [int howmany] Number of observations (number of DFT to perform)
                                    (fftw_complex *) X, // [fftw_complex *in] Input array is X (it is cast to non-const but will not be modified since FFTW_DESTROY_INPUT is not set)
                                    NULL, // [const int *inembed] Distance between each rank in input array (Not used since rank=1)
                                    1, // [int istride] Distance between successive variables in input array (in unit of fftw_complex)
                                    nvar, // [int idist] Distance between 2 observations in input array (in unit of fftw_complex)
                                    (fftw_complex *) Y, // [fftw_complex *out] Output array is Y
                                    NULL, // [const int *onembed] Distance between each rank in output array (Not used since rank=1)
                                    1, // [int ostride] Distance between successive variables in output array (in unit of fftw_complex)
                                    nvar, // [int odist] Distance between 2 observations in output array (in unit of fftw_complex)
                                    sign, // sign of the exponent in the formula that defines the Fourier transform (-1 or +1)
                                    FFTW_ESTIMATE); // [unsigned flags] Quickly choose a plan without performing full benchmarks (maybe sub-optimal but take less time)
  // If plan building fails, quit
  if(!plan)
    return NULL;
  
  // Execute FFTW plan
  fftw_execute(plan);
  
  // If we compute the inverse transform (sign == 1), normalize result by nvar (FFTW compute unnormalized transform)
  if(sign == 1) {
    nelem = 2 * nvar * nobs;
    for(i = 0; i < nelem; i++)
      Y[i] /= nvar;
  }
      
  
  // Destroy FFTW plan
  fftw_destroy_plan(plan);
  
  return Y;
}

double *fft(double *F, const double *X, const unsigned long nvar, const unsigned long nobs) {
  return _fft(F,X,nvar,nobs,-1);
}

double *ifft(double *F, const double *X, const unsigned long nvar, const unsigned long nobs) {
  return _fft(F,X,nvar,nobs,1);
}
