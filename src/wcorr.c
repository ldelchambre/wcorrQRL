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
* @file wcorr.c
* @brief Implementation of the wcorrQRL core functions.
* @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
* @version 1.0
* @date 11 March 2016
*/
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <fft.h>
#include <wcorr.h>

double *wcorr_lookup_alloc(const unsigned long nvar, const unsigned long ncomp) {
  return fft_alloc(sizeof(double)*fftr_size(nvar,ncomp*(ncomp+3)/2)); // ncomp*(ncomp+3)/2 is the number of columns of the lookup tables
}

void wcorr_lookup_free(double *L) {
  fft_free(L);
}

double *wcorr_lookup(double *L, const double *T, const unsigned long nvar, const unsigned long ncomp) {
  const unsigned long ncol = ncomp*(ncomp+3)/2; // Number of column in the lookup table
  const double (*Tp)[nvar] = (const double (*)[nvar]) T; // Reinterpret T as being an array of observations
  double (*Lp)[nvar] = (double (*)[nvar]) L; // Reinterpret L as being an array of observations
  unsigned long i,j,k,l;// General usage counters

  // Fill the lookup table with the element-wise product of the templates plus
  // hte individuals templates 
  i = 0;
  for(j = 0; j < ncomp; j++,i++) {
    for(k = j; k < ncomp; k++,i++)
      for(l = 0; l < nvar; l++)
        Lp[i][l] = Tp[j][l] * Tp[k][l]; // L(:,i) = T(:,j) .* T(:,k)
    for(l = 0; l < nvar; l++)
      Lp[i][l] = Tp[j][l]; // L(:,i) = T(:,j)
  }
  
  return fftr(L,L,nvar,ncol); // Return the DFT of the extended array of templates
}

/**
 * Compute the weighted cross correlation function of an observation s having
 * weights w against a set of templates T and for which the QR inner product 
 * table is given in A.
 *
 * @param[out] y Pointer to the array where to store the resulting cross correlation function.
 * @param[in,out] A Pointer to the inner product lookup table.
 * @param[in] w Pointer to the weight vector associated with the observations s.
 * @param[in] s Pointer to the observation.
 * @param[in] T Pointer to the templates observations to use. Size of (nvar x ncomp) in column-major order.
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 *
 * @see Delchambre L., 2016, MNRAS, "Redshift determination through weighted phase correlation: a linearithmic implementation"
 */
void _wcorr(double *y, double *A, const double *w, const double *s, const double *T, const unsigned long nvar, const unsigned long ncomp) {
	// Matrices used in the QR orthogonalization
	double X[ncomp*(ncomp+1)]; // X = [(W.*T)(1:ncomp,:) (W.*S)(1:ncomp)]
	
	// General usage variable
	double alpha, gamma, temp;
  unsigned long i,j,k,l,il,z;
	
	// Reinterpret some matrices in order to ease understanding
	double (*Xp)[ncomp] = (double (*)[ncomp]) X;
	const double (*Tp)[nvar] = (const double (*)[nvar]) T;
	double (*Ap)[nvar] = (double (*)[nvar]) A;
	
	// Browse all shift
	for(z = 0; z < nvar; z++) {
		// Build X and Y matrices
		for(i = 0; i < ncomp; i++) {
		  for(j = 0; j < ncomp; j++)
		    Xp[i][j] = w[j] * Tp[i][(j+z)%nvar];
		  Xp[ncomp][i] = w[i] * s[i];
		}
		
		// Reduce all rows of X
		il = 0; // The column index within the lookup table
		for(i = 0; i < ncomp; i++) {
			alpha = Ap[il++][z];
			// Check that the ith columns has a non zero norm
			if(0. < alpha && DBL_EPSILON * Ap[i][z] < alpha) {
		    // Annihilate the ith column of X thanks to factorized QR algorithm
			  alpha = (Xp[i][i] < 0) ? -sqrt(alpha) : sqrt(alpha);
			  temp = (Xp[i][i] += alpha);
			  for(j = i + 1; j <= ncomp; j++) {
				  gamma = (Xp[j][i] + Ap[il++][z]/alpha) / temp;
				  for(k = i; k < ncomp; k++)
				    Xp[j][k] -= gamma * Xp[i][k];
		  	}
			
			  // Update the lookup table
			  l = il;
			  for(j = i + 1; j < ncomp; j++)
			    for(k = j; k <= ncomp; k++, l++)
			      Ap[l][z] -= Xp[j][i]*Xp[k][i];
			} else {
				// If the ith column of X has a zero norm then X is rank-deficient
				Xp[ncomp][i] = 0;
				il += ncomp - i;
			}
		}
		
		// Compute the resulting weighted cross correlation function
		temp = 0;
		for(i = 0; i < ncomp; i++)
		  temp += Xp[ncomp][i] * Xp[ncomp][i];
		y[z] = temp;
	}
}

/**
 * Compute the weighted cross correlation function of an observation having 
 * ncomp starting weights equals to zero and for which the inner product table 
 * is given by A. Note that the resulting cross correlation function will be
 * shifted according to the value of zpos.
 *
 * @param[out] y Pointer to the array where to store the resulting cross correlation function.
 * @param[in,out] A Pointer to the inner product lookup table.
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 * @param[in] zpos Shift that was applied to the observation in order to get leading zeros weights.
 *
 * @see Delchambre L., 2016, MNRAS, "Redshift determination through weighted phase correlation: a linearithmic implementation"
 */
void _wcorr_zopt(double *y, double *A, const unsigned long nvar, const unsigned long ncomp, const unsigned long zpos) {
	// Reinterpret A matrices in order to facilitate understanding
	double (*Ap)[nvar] = (double (*)[nvar]) A;
	
	// General usage variable
  double alpha, norm2, res;
	unsigned long i,j,k,l,next,il,z,iy;
	
	// Compute the cross correlation function insertion index by taking into 
	// account the shift that was applied to the observation.
	iy = nvar - zpos;
	
	// Browse all shift
	for(z = 0; z < nvar; z++) {
		res = 0;
		for(i = il = 0; i < ncomp; i++) {
			// Annihilate ith column
			norm2 = Ap[il][z];
			next = il + ncomp - i;
			// If ith column has a non zero length then this problem is rank deficient.
			if(0. < norm2 && DBL_EPSILON * Ap[i][z] < norm2) {
			  // Compute the resulting weighted cross correlation function
			  res += Ap[next][z]*Ap[next][z]/norm2;
			  // Update the lookup table
			  for(j = 1, l = next + 1; j < ncomp - i; j++) {
				  alpha = Ap[il+j][z] / norm2;
				  for(k = j; k <= ncomp - i; k++, l++)
			  	    Ap[l][z] -= alpha * Ap[il+k][z];
		      }
		  }
		  il = next + 1;
		}
		y[iy] = res;
		iy = (iy + 1) % nvar;
	}
}

double *wcorr(double *corr, const double *L, const double *T, const unsigned long nvar, const unsigned long ncomp, const double *W, const double *S, const unsigned long n, const unsigned long nobs, const int search_zopt) {
  // Reinterpret some matrices in order to ease understanding
  const unsigned long nfft = fftr_size(nvar,1);
  const double (*Lp)[nfft] = (const double (*)[nfft]) L;
  double (*corrp)[nvar] = (double (*)[nvar]) corr;
  
  const double *Sp; // Pointer to the spectrum of the current observation
  const double *Wp; // Pointer to the weights associated with Sp
	double (*FA)[nfft]; // Pointer to the inner product table in the Fourier domain
	double (*A)[nvar];// Pointer to the inner product table in the time domain
	double (*FY)[nfft]; // Pointer to the image matrix in the Fourier domain
  double (*Y)[nvar]; // Pointer to the image matrix in the time domain
  unsigned long obs; // The current observation index
  
  // Zero weights optimization variables
  unsigned long zpos;
  int zopt;
  
  // General usage variables
  unsigned long i,j,k,l;
  double temp;
	
	// Check that the number of variable within the observations is less or equal than the number of variables within the templates
	if(nvar < n)
	  return NULL;
	
	// Allocate matrices
	FA = (double (*)[nfft]) wcorr_lookup_alloc(nvar, ncomp);
	FY = (double (*)[nfft]) fft_alloc(sizeof(double) * 2 * nfft);
	A = (double (*)[nvar]) FA;
	Y = (double (*)[nvar]) FY;
	if(!FA || !FY) {
	  wcorr_lookup_free((double *) FA);
	  fft_free(FY);
	  return NULL;
	}
	
	// For each observation
	for(obs = 0; obs < nobs; obs++) {
	  // Initialize Sp and Wp for this observation
	  Sp = &S[obs*n];
	  Wp = &W[obs*n];

		// Check if we can shift Sp such that we have ncomp consecutive zeros within Wp
    zopt = 0;
    zpos = nvar;
    if(ncomp <= nvar - n) { // By considering the padding we have to set
      zopt = 1;
      zpos = n; // We have ncomp consecutive zeros beyond Sp
    } else if(search_zopt) { // Otherwise, try to find ncomp consecutive zeros within Wp (if search_zopt != 0)
      i = 0;
      while(i < n && !zopt) {
        // Browse all zeros from i
        j = i;
        while(j < n && Wp[j] == 0 && j - i < ncomp)
          j++;
        // We find ncomp consecutive zeros
        if(j - i == ncomp) {
          zopt = 1;
          zpos = i;
        } else
          i = j + 1;
      }
    }
    // Here:
    //   if zopt == 1 => Wp[zpos:zpos+ncomp-1] == 0
    
    // Copy the square weights and weighted observation into Y by taking into
    // account the potential shift leading to the 'zero-optimized' mode
    if(zopt) {
      for(i = nvar - zpos, j = 0; j < zpos; i++, j++) {
				temp = Wp[j]*Wp[j];
				Y[0][i] = temp;
				Y[1][i] = temp*Sp[j];
			}
			for(i = 0; i < ncomp; i++)
			  Y[0][i] = Y[1][i] = 0;
		  j += ncomp;
    } else
      i = j = 0;
		for(; j < n; i++, j++) {
			temp = Wp[j]*Wp[j];
			Y[0][i] = temp;
			Y[1][i] = temp*Sp[j];
		}
		
		// Pad Y with 0
		for(; j < nvar; i++, j++)
		  Y[0][i] = Y[1][i] = 0;
      
    // Convert Y to frequency domain
    fftr((double *) FY, (const double *) Y, nvar, 2);
    
    // Build the inner product table A in the frequency domain
    i = 0;
    for(j = 0; j < ncomp; j++) {
      for(k = j; k < ncomp; k++) {
        for(l = 0; l < nfft; l += 2) {
          // (FA[i][l] + FA[i][l+1] i) = conj(FY[0][l]+FY[0][l+1] i) * (Lp[i][l] * Lp[i][l+1] i)
          FA[i][l]   = FY[0][l]*Lp[i][l]   + FY[0][l+1]*Lp[i][l+1];
          FA[i][l+1] = FY[0][l]*Lp[i][l+1] - FY[0][l+1]*Lp[i][l];
        }
        i++;
      }
      for(l = 0; l < nfft; l += 2) {
        // (FA[i][l] + FA[i][l+1] i) = conj(FY[1][l]+FY[1][l+1] i) * (Lp[i][l] * Lp[i][l+1] i)
        FA[i][l]   = FY[1][l]*Lp[i][l]   + FY[1][l+1]*Lp[i][l+1];
        FA[i][l+1] = FY[1][l]*Lp[i][l+1] - FY[1][l+1]*Lp[i][l];
      }
      i++;
    }

    // Convert the inner product table into the time domain
    ifftr((double *) A, (const double *) FA, nvar, i);
    
    // Execute the factorized QR algorithm with lookup table algorithm
    if(zopt)
      _wcorr_zopt(corrp[obs], (double *) A, nvar, ncomp, zpos); // In the 'zero-optimized' mode
    else
      _wcorr(corrp[obs], (double *) A, Wp, Sp, T, nvar, ncomp); // In the classical mode
	}
	
	// Free temporary arrays
	wcorr_lookup_free((double *) FA);
	fft_free(FY);
	
	return corr;
}
