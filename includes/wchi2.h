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
 * @file wchi2.h
 * @brief Prototypes of the normal equation approach to the weighted phase correlation problem.
 * @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
 * @version 1.0
 * @date 11 March 2016
 */
#ifndef _WCHI2_H
#define _WCHI2_H

/**
 * @param W Weights (in inverse sigma unit) 
 */
 
 /**
 * @brief Computation of the weighted phase correlation through the solution of the normal equations problems.
 *
 * Compute the minimal chi-squares of a set of observations S of size 
 * (n x nobs), whose weights are given by W [also of size (n x nobs)], 
 * against a set of template T of size (nvar x ncomp) by using the solution 
 * the the normal equations.
 *
 * @param[out] y Pointer to the array containing the minimal chi-square associated with each phase. Size of (nvar x nobs) in column-major order.
 * @param[in] T Pointer to the templates observations to use. Size of (nvar x nobs) in column-major order.
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 * @param[in] W Pointer to the weight matrix associated with the observations S. Size of (n x nobs) in column-major order. [If n == nvar, then y can be either S or W]
 * @param[in] S Pointer to the matrix of observations. Size of (n x nobs) in column-major order.
 * @param[in] n The number of variables (rows) within S and W.
 * @param[in] nobs The number of observations (columns) within S and W.
 *
 * @return 0 upon success, -1 if we failed to retrieve the solution to the normal equations for any phase.
 */
int wchi2_normeq(double *y, const double *T, const unsigned long nvar, const unsigned long ncomp,  const double *W, const double *S, const unsigned long n, const unsigned long nobs);

#endif
