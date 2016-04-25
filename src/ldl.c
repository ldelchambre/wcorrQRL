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
 * @file ldl.c
 * @brief Implementation of the LDL matrix decomposition and solving functions.
 * @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
 * @version 1.0
 * @date 11 March 2016
 */
#include <ldl.h>

int ldl_decomposition(double *A, const unsigned long size) {
  double (*M)[size] = (double (*)[size]) A;
	unsigned long i, j, k;
	double sum;

	// Compute decomposition of each column
	for (j = 0; j < size; j++) {

		// Compute diagonal element
		sum = M[j][j];
		for (k = 0; k < j; k++)
			sum -= M[k][k] * M[k][j] * M[k][j];

		// If sum == 0, then A is rank deficient
		if (sum == 0.)
			return -1;
		M[j][j] = sum;

		// Compute jth column
		for (i = j + 1; i < size; i++) {
			sum = M[j][i];
			for (k = 0; k < j; k++)
				sum -= M[k][k] * M[k][i] * M[k][j];
			M[j][i] = sum / M[j][j];
		}
	}
	
	return 0;
}

double *ldl_solve(double *x, const double *A, const double *y,
		const unsigned long size) {
	const double (*L)[size] = (const double (*)[size]) A;
	unsigned long i, j;
	double sum;

	// Perform forward substitution
	for (i = 0; i < size; i++) {
		sum = y[i];
		for (j = 0; j < i; j++)
			sum -= L[j][i] * x[j];
		x[i] = sum;
	}

	// Divide answer by diagonal element
	for (i = 0; i < size; i++) {
	  if(L[i][i] != 0.)
		  x[i] /= L[i][i]; 
		else
		  x[i] = 0;
  }

	// Backward substitution
	for (i = size; i--;) {
		sum = x[i];
		for (j = size; i < --j;)
			sum -= L[i][j] * x[j];
		x[i] = sum;
	}

	return x;
}
