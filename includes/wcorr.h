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
 * @file wcorr.h
 * @brief Prototypes of the wcorrQRL core functions.
 *
 * Given a set of templates stored in T of size (nvar x ncomp) in column-major 
 * order and a set of observations S having weights W each of size (n x nobs) 
 * and also stored in column major order, we have that the cross correlation
 * function resulting from the weighted phase correlation of T with S whose
 * weights are W will be stored in corr [size of (nvar x nobs) in column-major
 * order] by the following code :
 * @code{.c}
 *   double *L = wcorr_lookup_alloc(nvar,ncomp); // Allocated the lookup table
 *   wcorr_lookup(L,T,nvar,ncomp); // Build the lookup table
 *   wcorr(corr,L,T,nvar,ncomp,W,S,n,nobs,0); // Compute the weighted cross correlation function
 *   wcorr_lookup_free(L); // Delete the lookup table
 * @endcode
 *
 * @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
 * @version 1.0
 * @date 11 March 2016
 * @see Delchambre L., 2016, MNRAS, "Redshift determination through weighted phase correlation: a linearithmic implementation"
 */
#ifndef _WCORR_H_
#define _WCORR_H_

/**
 * @brief Lookup table allocation
 *
 * Allocate memory designed to contains a lookup table associated with ncomp 
 * templates each consisting in nvar variables.
 *
 * @param[in] nvar The number of variable within each template.
 * @param[in] ncomp The number of template that will be used.
 *
 * @return Pointer to the allocated memory or NULL if allocation fails.
 * @warning Memory allocated with #wcorr_lookup_alloc must be freed with #wcorr_lookup_free.
 */
double *wcorr_lookup_alloc(const unsigned long nvar, const unsigned long ncomp);

/**
 * @brief Lookup table freeing
 * 
 * Free the memory previously allocated with #wcorr_lookup_alloc.
 *
 * @param[in] L Pointer to the memory previously allocated with #wcorr_lookup_alloc.
 * @note If L is NULL, then nothing happens.
 */
void wcorr_lookup_free(double *L);

/**
 * @brief Lookup table building
 *
 * Build the initial lookup table associated with the ncomp templates of size 
 * nvar stored in T. The resulting initial lookup table will be stored in L.
 *
 *   The lookup table is defined as the DFT of each template, plus the DFT of
 * the element-wise multiplication of each individual template with all other
 * templates (including itself).
 *
 * @param[out] L The pointer to the memory where to store the resulting lookup table (should be allocated with #wcorr_lookup_alloc).
 * @param[in] T Pointer to the templates observations to use. Size of (nvar x ncomp) in column-major order.
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 *
 * @return Pointer to L if DFT succeed, NULL otherwise.
 */
double *wcorr_lookup(double *L, const double *T, const unsigned long nvar, const unsigned long ncomp);

/**
 * @brief Computation of the weighted cross correlation function through factorized QR algorithm with lookup table
 *
 * Compute the weighted cross correlation function of a set of observations S of size 
 * (n x nobs), whose weights are given by W [also of size (n x nobs)], 
 * against a set of template T of size (nvar x ncomp) by using the lookup 
 * table L through a factorized QR algorithm with lookup table as described in
 * Delchambre L., 2016, MNRAS, "Redshift determination through weighted phase 
 * correlation: a linearithmic implementation".
 *
 * @param[out] corr Pointer to the array where to store the resulting cross correlation function. Size of (nvar x nobs) in column-major order.
 * @param[in] L Pointer to the lookup table to use as returned by #wcorr_lookup (L, T, nvar, ncomp).
 * @param[in] T Pointer to the templates observations to use. Size of (nvar x ncomp) in column-major order.
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 * @param[in] W Pointer to the weight matrix associated with the observations S. Size of (n x nobs) in column-major order.
 * @param[in] S Pointer to the matrix of observations. Size of (n x nobs) in column-major order.
 * @param[in] n The number of variables (rows) within S and W.
 * @param[in] nobs The number of observations (columns) within S and W.
 * @param[in] search_zopt If search_zopt is different than 0, the method will seek within each observation for a set of ncomp consecutive zero weights in order to optimize processing. Otherwise, optimization will only occurs if ncomp <= nvar - n. Recommended value is 0.
 *
 * @return Pointer to corr upon success. NULL if nvar < n or temporary arrays allocation fails.
 * @note Setting the search_zopt parameter to a value different than 0 will slow down the algorithm if most of the observations within S, don't have any ncomp consecutive zero weights within W.
 */
double *wcorr(double *corr, const double *L, const double *T, const unsigned long nvar, const unsigned long ncomp, const double *W, const double *S, const unsigned long n, const unsigned long nobs, const int search_zopt);

#endif
