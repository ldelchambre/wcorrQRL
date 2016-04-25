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
 * @file demo.c
 * @author Ludovic Delchambre <ldelchambre@ulg.ac.be>
 * @version 1.0
 * @date 11 March 2016
 * @brief wcorrQRL demonstration program.
 *
 *  This demonstration program compute mock templates upon which it build 
 * shifted mock observation coming along with their weight. It further compute
 * the phase at which the observation match at best the templates in a weighted
 * least-squares sense. This last point is performed by using two methods,
 * namely the solutions of the normal equations associated with each of the 
 * weighted least-squares problem through the use of a Cholesky-related matrix 
 * decomposition algorithm; and the factorized QR decomposition with lookup 
 * tables.
 *
 *   This program print to the standard output the observation index; the 
 * effective shift that was applied to the observation; the execution time of
 * the factorized QR decomposition with lookup tables; its estimated shift;
 * the execution time of the solution to the normal equations and its estimated
 * shift.
 */
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <wchi2.h>
#include <wcorr.h>

#define DEFAULT_NVAR 10000 /**< Default number of templates variables to use */
#define DEFAULT_NCOMP 10   /**< Default number of templates to use */
#define DEFAULT_N 10000    /**< Default number of observation variables to use */
#define DEFAULT_NOBS 32    /**< Default number of observation to use */
#define DEFAULT_NOISE 1e-2 /**< Default maximal noise to add to the mock data [sigma unit] */

/**
 * @brief Generate normally distributed random variables of mean mu and of standard deviation sigma.
 * 
 * @param[in] mu The mean of the normal distribution to use
 * @param[in] sigma the standard deviation of the normal distribution to use
 *
 * @return Normally distributed random variables of mean mu and of standard deviation sigma
 * @note This code was adapted from Press W.H., Tuekolsky S.A., Vetterling W.T., Flannery B.P., 2002, Numerical recipes in C++: The Art of Scientific Computing, 2nd edn. Cambridge Univ. Press, New York
 */
double randn(const double mu, const double sigma) {
  static int isset = 0;
  static double z0 = 0, z1 = 0; // Dummy initialization to avoid compilation warning
  double u0,u1;
  
  // If z was already computed
  if(isset) {
    isset = 0; // Say that we have to re-compute z0, z1
    return z1*sigma+mu;
  }
  
  // Pick up u1,u2 in [DBL_EPSILON, 1]
  do {
    u0 = 1./RAND_MAX*rand();
    u1 = 1./RAND_MAX*rand();
  } while(u0 <= DBL_EPSILON);
  
  // Compute the Box-Muller transform of u1, u2
  z0 = sqrt(-2.0 * log(u0)) * cos(2.*M_PI*u1);
	z1 = sqrt(-2.0 * log(u0)) * sin(2.*M_PI*u1);
	
	// Say that z1 is still usable
	isset = 1;
	
	return z0 * sigma + mu;
}

/**
 * @brief Build random templates
 *
 *   Build ncomp random template having nvar variables uniformly distributed in
 * the interval [-1,1].
 *
 * @param[out] T Pointer to the templates observations in column-major order.
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 */
void build_mock_templates(double *T, const unsigned long nvar, const unsigned long ncomp) {
  double (*Tp)[nvar] = (double (*)[nvar]) T;
  unsigned long i,j;
  
  for(i = 0; i < ncomp; i++)
    for(j = 0; j < nvar; j++)
      Tp[i][j] = 2.*rand()/RAND_MAX-1.;
}

/**
 * @brief Build noisy and shifted observations based on a set of templates
 *
 * 	Build an observation, s, of size n along with its associated weights, w, 
 * and taken as the shifted linear combination of the templates T 
 * (of size nvar x ncomp) to which we add a normal noise.
 *
 * @param[out] s Pointer to the resulting observation
 * @param[out] w Pointer to the resulting weights
 * @param[in] n The size of s and w
 * @param[in] T Pointer to the templates
 * @param[in] nvar The number of variables (rows) within T.
 * @param[in] ncomp The number of templates (columns) within T.
 * @param[in] shift The shift to apply to the generated observation (shift <= nvar).
 * @param[in] noise The maximal noise to add to the observation [in sigma unit]. 
 */
void build_mock_observation(double *s, double *w, const unsigned long n, const double *T, const unsigned long nvar, const unsigned long ncomp, const unsigned long shift, const double noise) {
  const double (*Tp)[nvar] = (const double (*)[nvar]) T;
  double a, sigma;
  unsigned long i,j;
  
  // Initialize observation
  for(i = 0; i < n; i++)
    s[i] = 0;
  
  // Compute a random linear combination of the templates at the observation
  for(i = 0; i < ncomp; i++) {
    // Pick-up the random linear coefficient to apply to the template
    a = (2.* rand()/RAND_MAX-1.);
    // Add the shifted template 
    for(j = 0; j < n; j++)
      s[j] += Tp[i][(j + shift) % nvar] * a;
  }
  
  // Add noise to the observation
  for(i = 0; i < n; i++) {
    sigma = noise * (1. + rand()) / (1. + RAND_MAX); // sigma = ]0, noise]
    s[i] += randn(0.,sigma);      
    w[i] = 1./sigma;
  }
}

/**
 * @brief Compute the minimal value between a and b
 * @param[in] a A number
 * @param[in] b Another number
 * @return 1 if a < b; 0 otherwise
 */
int min(const double a, const double b) {
  return (a < b);
}

/**
 * @brief Compute the maximal value between a and b
 * @param[in] a A number
 * @param[in] b Another number
 * @return 1 if a > b; 0 otherwise
 */
int max(const double a, const double b) {
  return (a > b);
}

/**
 * @brief Find the index of the value contained within x that is minimal according to a given ordering sequence.
 * @param[in] x Pointer to the array where to search for the minimal ordered value
 * @param[in] nvar The number of elements within x
 * @param[in] compar Pointer to a function designed to define the ordering sequence to use
 * @return The index of the minimal ordered value defined according to the compar function
 */
unsigned long find_shift(const double *x, unsigned long nvar, int (*compar)(const double a, const double b)) {
  unsigned long imin = 0;
  double cmin = (nvar) ? x[0] : 0;
  unsigned long i;
  
  for(i = 1; i < nvar; i++) {
    if(compar(x[i], cmin)) {
      imin = i;
      cmin = x[i];
    }
  }
  
  return imin;
}

/**
 * @brief Compute the time difference between two CPU clock count.
 * @param[in] t0 First CPU clock count (farthest)
 * @param[in] t1 Second CPU clock count (nearest)
 * @return The time difference (in seconds) between t0 and t1 or -1 either if t0 or t1 is equal to (clock_t) -1 or that a wrap around was detected
 * @warning The CPU clock count can wrap around, on 32 bits system this happens each 72 minutes, on 64 bits systems this happens each 584540 years.
 */
double clock_diff(const clock_t t0, const clock_t t1) {
  if(t0 != (clock_t) -1 && t1 != (clock_t) -1 && t0 < t1)
    return ((double) (t1 - t0)) / CLOCKS_PER_SEC; // Return the time difference (in seconds) between the two clock time
  else
    return -1.; // If either t0 = -1 or t1 = -1 or the clock time wrap around, return -1
}

/**
 * @brief Main function of the demonstration program
 * 
 * Call as ./demo [nvar n ncomp nobs noise]
 *
 * @param[in] argc The Number of input arguments
 * @param[in] argv The array of string argument(s)
 * @return 0 upon success and a number different than 0 otherwise
 */
int main(const int argc, const char **argv) {
  // Get program input parameter
  const unsigned long nvar = (argc > 1) ? strtoul(argv[1],NULL,10) : DEFAULT_NVAR;
  const unsigned long n = (argc > 2) ? strtoul(argv[2],NULL,10) : DEFAULT_N;
  const unsigned long ncomp = (argc > 3) ? strtoul(argv[3],NULL,10) : DEFAULT_NCOMP;
  const unsigned long nobs = (argc > 4) ? strtoul(argv[4],NULL,10) : DEFAULT_NOBS;
  const double noise = (argc > 5) ? strtod(argv[5],NULL) : DEFAULT_NOISE;
  double *T = NULL; // Pointer to templates
  double *s = NULL; // Pointer to the current observation
  double *w = NULL; // Pointer to the weights
  double *corr = NULL; // Pointer to the weighted cross correlation function
  double *L = NULL; // Pointer to the lookup table
  clock_t t0, t1; // Clock times
  double tqr,tneq; // Execution times respectively for the QRL algorithm and for the normal equation
  unsigned long shift, shiftqr, shiftneq; // Real shift and shift found by the two algorithms
  unsigned long i;
  
  // Allocate arrays
  T = (double *) calloc(nvar*ncomp,sizeof(double));
  s = (double *) calloc(n,sizeof(double));
  w = (double *) calloc(n,sizeof(double));
  corr = (double *) calloc(nvar,sizeof(double));
  L = wcorr_lookup_alloc(nvar, ncomp);
  
  // Check that all arrays were correctly allocated
  if(T == NULL || s == NULL || w == NULL || corr == NULL || L == NULL) {
    free(T);free(s);free(w);free(corr);free(L);
    fprintf(stderr, "An error occurs during array allocation\n");
    return -1;
  }
  
  // Init the random number generator
  srand(1);
  
  // Build mock templates
  build_mock_templates(T, nvar, ncomp);
  
  // Build the initial lookup table associated with T
  wcorr_lookup(L, T, nvar, ncomp);
  
  // Print the usage we made
  fprintf(stderr, "Usage: %s nvar=%lu n=%lu ncomp=%lu nobs=%lu noise=%g\n", argv[0], nvar, n, ncomp, nobs, noise);
  
  // Check if we will run in zero-optimized mode?
  if(nvar - n < ncomp)
    fprintf(stderr, "Running demonstration in classical mode\n");
  else
    fprintf(stderr, "Running demonstration in zero-optimized mode\n");

  // Print header line
  fprintf(stdout, "Observation | Real shift | QRL CPU time | QRL shift | Neq. CPU time | Neq. shift \n");
  
  // Compute weighted phase correlation for each observation
  for(i = 0; i < nobs; i++) {
    // Select a random shift
    shift = rand() % nvar;
    
    // Build a mock observation
    build_mock_observation(s,w,n,T,nvar,ncomp,shift,noise); 
    
    // Compute the weighted phase correlation function through factorized QR algorithm with lookup tables
    t0 = clock();
    wcorr(corr, L, T, nvar, ncomp, w, s, n, 1, 0.);
    t1 = clock();
    
    // Get the found shift
    shiftqr = find_shift(corr, nvar, max);
    
    // Get the CPU execution time
    tqr = clock_diff(t0, t1);
      
    // Compute the weighted phase correlation function through normal equations
    t0 = clock();
    wchi2_normeq(corr, T, nvar, ncomp,  w, s, n, 1);
    t1 = clock();
    
    // Get the found shift
    shiftneq = find_shift(corr, nvar, min);
    
    // Get the CPU execution time
    tneq = clock_diff(t0, t1);
    
    // Print result
    fprintf(stdout, "%11lu | %10lu |  %10.5fs |  %8lu |   %10.5fs | %10lu \n", i, shift, tqr, shiftqr, tneq, shiftneq);
    
    // 
    fflush(stdout);
  }
  
  // Free arrays
  free(T);
  free(s);
  free(w);
  free(corr);
  wcorr_lookup_free(L);
  
  return 0;
}
