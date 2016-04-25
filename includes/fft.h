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
* @file fft.h
* @brief Prototypes of FFTW3 interface functions.
*
*  The functions designed within this file aims at facilitate and shorten
* the calls to the FFTW3 library functions.
*
* @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
* @version 1.0
* @date 11 March 2016
*/

#ifndef _FFT_H_
#define _FFT_H_
/**
 * @brief Real DFT size calculation
 *
 * Return the size in unit of double of the array designed to contains the
 * DFT of a real data matrix of size (nvar x nobs). Due to the fact that DFTs 
 * of real data are symmetric, the transform of a matrix of size (nvar x nobs)
 * will be a complex matrix of size (nvar/2+1 x nobs).
 *
 * Example:
 * @code{.c}
 *   double *X = fft_alloc(nvar*sizeof(double));
 *   double *Y = fft_alloc(fftr_size(nvar)*sizeof(double));
 *   fftr(Y,X,nvar);
 *   ifftr(X,Y,nvar);
 *   fft_free(X);fft_free(Y);
 * @endcode
 *
 * @param[in] nvar The number of variable (rows) within matrix
 * @param[in] nobs The number of observations (columns) within matrix
 *
 * @return The size (in unit of double) of the array designed to contains the DFT of the real data matrix.
 */
unsigned long fftr_size(const unsigned long nvar, const unsigned long nobs);

/**
 * @brief DFT memory allocation
 *
 * Allocate array in a way that is optimal in order to perform DFT though
 * #fftr, #ifftr, #fft and #ifft functions (FFTW3 use "Single Instruction 
 * Multiple Data" [SIMD] instructions such as to fasten data processing).
 *
 * @param[in] size The size in byte of the array to allocate
 *
 * @return Pointer to the allocated memory or NULL if allocation fails.
 * @note Any memory allocated with #fft_alloc must be freed using #fft_free.
 */
void *fft_alloc(const unsigned long size);

/**
 * @brief DFT memory freeing
 *
 * Free memory allocated by #fft_alloc
 *
 * @param[in] ptr Pointer returned by #fft_alloc
 */
void fft_free(void *ptr);

/**
 * @brief Real DFT calculation
 *
 * Compute the DFT of each observation of the real matrix X of size (nvar x nobs), and
 * store result within the complex matrix Y of size (nvar/2+1 x nobs).
 * Matrix Y is of type double[2], first number being real part, second being imaginary part.
 * The fact that only nvar/2+1 rows-doublets are returned within Y is due to the fact that DFTs of
 * real data are symmetric and are consequently not duplicated.
 *
 * @param[out] Y Output matrix of complex number (double[2]) of size (nvar/2+1 x nobs) [must be at least of size #fftr_size (nvar,nobs)]
 * @param[in] X Input matrix of real number (double) of size (nvar x nobs) [can be the same as Y]
 * @param[in] nvar Number of variables (rows) within X
 * @param[in] nobs Number of observations (columns) within X
 *
 * @return Pointer to Y or NULL if DFT fails
 */
double *fftr(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs);

/**
 * @brief Real inverse DFT calculation
 *
 * Compute the inverse DFT of each observation of the complex matrix X of size 
 * (nvar/2+1 x nobs), and store (real) result within Y of size (nvar x nobs). 
 * Matrix X is of type double[2], first number being real part, second being imaginary part.
 * Matrix X must be thought as a symmetric matrix as defined if #fftr function.
 *
 * @param[out] Y Output matrix of real number of size (nvar x nobs)
 * @param[in] X Input matrix of complex number (double[2]) of size (nvar/2+1 x nobs) [can be the same as Y]
 * @param[in] nvar Number of variables (rows) within Y
 * @param[in] nobs Number of observations (columns) within Y
 *
 * @return Pointer to Y or NULL if DFT fails
 */
double *ifftr(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs);

/**
 * @brief Complex DFT calculation
 *
 * Compute the DFT of the complex matrix X of size (nvar x nobs) into the complex matrix
 * Y of size (nvar x nobs).
 * 
 * @param[out] Y Output matrix of complex numbers (double[2]) of size (nvar x nobs)
 * @param[in] X Input matrix of complex numbers (double[2]) of size (nvar x nobs) [can be the same as Y]
 * @param[in] nvar Number of variables (rows) within X and Y
 * @param[in] nobs Number of observations (columns) within X and Y
 *
 * @return Pointer to Y or NULL if DFT fails
 */
double *fft(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs);

/**
 * @brief Complex inverse DFT calculation
 *
 * Compute the inverse DFT of the complex matrix X of size (nvar x nobs) into the complex matrix
 * Y of size (nvar x nobs).
 * 
 * @param[out] Y Output matrix of complex numbers (double[2]) of size (nvar x nobs)
 * @param[in] X Input matrix of complex numbers (double[2]) of size (nvar x nobs) [can be the same as Y]
 * @param[in] nvar Number of variables (rows) within X and Y
 * @param[in] nobs Number of observations (columns) within X and Y
 *
 * @return Pointer to Y or NULL if DFT fails
 */
double *ifft(double *Y, const double *X, const unsigned long nvar, const unsigned long nobs);
#endif
