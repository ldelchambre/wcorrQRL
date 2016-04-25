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
 * @file ldl.h
 * @brief Prototypes of the LDL matrix decomposition and solving functions.
 *
 *  LDL decomposition is a Cholesky-related decomposition of the form A = L*D*L'
 * where A is the symmetric matrix to be decomposed; L is a lower unit triangular
 * matrix and D is a diagonal matrix.
 *  Its main advantage over the Cholesky decomposition is that it does not 
 * require to compute any square roots and allow negative diagonal elements.
 *
 * @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
 * @version 1.0
 * @date 11 March 2016
 * @see Golub G.H., Van Loan C.F., 1996, Matrix Computations, 3rd edn. The Johns Hopkins Univ. Press, London
 */

#ifndef LDL_H_
#define LDL_H_

/**
 * @brief Perform a LDL decomposition of the matrix A
 *
 *  Compute the decomposition A = L*D*L' where L is a lower unit triangular
 * matrix and D is a diagonal matrix.
 *
 *	This method browse the lower left part of A and perform an in-place
 * LDL decomposition. After call, A will contains the L matrix within its lower
 * left part with the unit diagonal elements being replaced with D.
 *
 * @param[in,out] A Pointer to the symmetric matrix to be decomposed in column-major order.
 * @param[in] size The number of columns and rows of A
 * 
 * @return 0 if the LDL decomposition succeed, -1 if A is rank deficient.
 * @note This decomposition has an algorithmic complexity of O(N^3) flops.
 */
int ldl_decomposition(double *A, const unsigned long size);

/**
 * @brief Solve the linear system of equations based on its LDL decomposition
 *
 *  Solve the linear system of equation A * x = y where A comes from a LDL 
 * decomposition of a symmetric matrix (see #ldl_decomposition).
 *
 * @param[out] x Pointer to the solution vector
 * @param[in] A Pointer to a LDL decomposed matrix
 * @param[in] y Pointer to the image vector [can be the same as x]
 * @param[in] size The size of y 
 *
 * @return Pointer to x
 * @note This method has an algorithmic complexity of O(N^2) flops
 */
double *ldl_solve(double *x, const double *A, const double *y,	const unsigned long size);

#endif
