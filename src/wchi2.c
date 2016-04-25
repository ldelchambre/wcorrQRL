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
 * @file wchi2.c
 * @brief Implementation of the normal equation approach to the weighted phase correlation problem.
 * @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
 * @version 1.0
 * @date 11 March 2016
 */
#include <ldl.h>
#include <wchi2.h>

int wchi2_normeq(double *y, const double *T, const unsigned long nvar, const unsigned long ncomp,  const double *W, const double *S, const unsigned long n, const unsigned long nobs) {
  // Pointer to array
  const double (*Tp)[nvar] = (const double (*)[nvar]) T; // Pointer to T
  const double (*Wp)[n] = (const double (*)[n]) W; // Pointer to W
  const double (*Sp)[n] = (const double (*)[n]) S; // Pointer to S
  double (*yp)[nvar] = (double (*)[nvar]) y; // Pointer to result matrix y
  const double *aw, *as; // Temporary array to rows of W and S
  
  // Normal equation matrices
  double A[ncomp*ncomp];
  double (*Ap)[ncomp] = (double (*)[ncomp]) A;
  double b[ncomp];
  
  // Non zero weight related variables
  unsigned long Wnz[n];
  unsigned long nz;
  
  // Result for one observation
  double res[nvar];
  
  // General usage variables
  unsigned long obs,i,j,k,l,z;
  double temp, sum, w2;
  
  // For each observation
  for(obs = 0; obs < nobs; obs++) {
    // Initialize aw and as for this observation
    aw = Wp[obs];
    as = Sp[obs];
    
    // Find non-zero weights variables
    for(i = nz = 0; i < n; i++) {
      if(aw[i] != 0.)
        Wnz[nz++] = i;
    }
    
    // Browse each shift
    for(z = 0; z < nvar; z++) {
    
      // Build the normal equation design matrices A and b in Ax=b
      for(i = 0; i < ncomp; i++){
        for(j = 0; j <= i; j++) {
          sum = 0;
          for(k = 0; k < nz; k++) {
            l = Wnz[k];
            w2 = aw[l] * aw[l];
            l = (l + z) % nvar;
            sum += Tp[i][l] * w2 * Tp[j][l];
          }
          Ap[j][i] = sum;
        }
        
        // Build the normal equation design matrices b in Ax=b
        sum = 0;
        for(k = 0; k < nz; k++) {
          l = Wnz[k];
          w2 = aw[l] * aw[l] * as[l];
          l = (l + z) % nvar;
          sum += Tp[i][l] * w2;
        }
        b[i] = sum;
      }
      
      // LDL decompose A
      if(ldl_decomposition(A, ncomp))
        return -1;
      
      // Find the solution vector
      if(!ldl_solve(b, A, b,ncomp))
        return -1;
      
      // Compute chi square
      sum = 0;
      for(i = 0; i < nz; i++) {
        k = Wnz[i];
        l = (k + z) % nvar;
        temp = 0;
        for(j = 0; j < ncomp; j++)
          temp += b[j] * Tp[j][l];
        temp = aw[k] * (as[k] - temp);
        sum += temp * temp;
      }
      
      // Assign result
      res[z] = sum;
    }
    
    // Copy result within y (such that y can be S or W [but not T!])
    for(i = 0; i < nvar; i++)
      yp[obs][i] = res[i];
    
  }
  
  return 0;
}
