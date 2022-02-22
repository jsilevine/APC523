#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <chrono>
#include <math.h>
#include <omp.h>

#define NDEF    64
#define NGHOST  1
#define MAXITER 100000
#define OUTFREQ 1000
#define FOUT    1

int main(int argc, char *argv[]) {
  // Initialize variables
  int n, nwg, np = 1; //np is number of cells in parallel block
  int t,i,j    = 0;
  int nthreads = 1;
  double h   = 1.;
  double eps = 1.e-6;

  // Parse input
  if (argc == 1) {
    n = NDEF;
  } else if (argc == 2) {
    n = atoi(argv[1]);
  } else if (argc == 3) {
    n = atoi(argv[1]);
    eps = atof(argv[2]);
  } else if (argc == 4) {
    n = atoi(argv[1]);
    eps = atof(argv[2]);
    nthreads = atoi(argv[3]); //take N threads as input to program
  } else {
    std::cout << "Usage: ./laplace [n eps nthreads]" << std::endl;
    exit(1);
  }
  nwg = n + 2*NGHOST;
  h   = 1./n; //gridcell size
  np = (n/nthreads); // size of parallel block

  // Allocate and fill the working array and previous step array
  double* arr_curr = (double*) calloc(nwg*nwg, sizeof(double));
  double* arr_prev = (double*) calloc(nwg*nwg, sizeof(double));

  // Fill in the source function
  double  x, y = 0.;
  double* fval = (double*) malloc(n*n*sizeof(double));
  for (i = 0; i < n; ++i) {
    x = (i + .5) * h;
    for (j = 0; j < n; ++j) {
      y = (j + .5) * h;
      fval[i*n+j] =
        pow(M_PI,2.) * (16.*pow(y,2.)+81.*pow(x,4.))
          * cos(2*M_PI*pow(y,2.)) * cos(3*M_PI*pow(x,3.))
        + 18.*M_PI*x * sin(3*M_PI*pow(x,3.)) * cos(2*M_PI*pow(y,2.))
        + 4. *M_PI   * cos(3*M_PI*pow(x,3.)) * sin(2*M_PI*pow(y,2.));
    }
  }

  // Main loop
  int     iter = 0;
  double  err  = std::numeric_limits<double>::max();
  double  terr = 0.;
  double  ffac = .25 * pow(h, 2.);
  double* tarr = NULL;
  auto    strt = std::chrono::steady_clock::now();

  // initialize threads
  omp_set_num_threads(nthreads);

  while (err > eps && iter < MAXITER) {
    // Jacobi method
    err = 0.;

#pragma omp parallel firstprivate(arr_prev, terr, np) shared(arr_curr) reduction(max:err) {
    #pragma omp for
    for (t = 0; t < nthreads; ++t) { // loop over parallel blocks
      for (i = NGHOST+t*np; i < NGHOST+(t+1)*np; ++i) { // only calculate values for interior cells
        for (j = NGHOST+t*np; j < NGHOST+(t+1)*np; ++j) {
          arr_curr[(i*nwg)+j] =
            0.25 * (arr_prev[(i-1)*nwg+j]     + arr_prev[(i+1)*nwg+j] +
                    arr_prev[i    *nwg+(j-1)] + arr_prev[i    *nwg+(j+1)]) -
            ffac * fval[(i-NGHOST)*n+(j-NGHOST)];
          terr = std::abs(arr_curr[i*nwg+j] - arr_prev[i*nwg+j]);
          err  = (terr > err) ? terr: err;
        }
      }
    }
    #pragma omp barrier // make sure threads have finished running before
                          // dealing with boundary conditions

  #pragma omp for
    for (t = 0; t < nthreads; ++t) { // loop over parallel blocks
      for (i = NGHOST+t*np; i < NGHOST+(t+1)*np; ++i) {
        for (j = 0; j < NGHOST; ++j) {
          arr_curr[i*nwg+j]              = arr_curr[i*nwg+NGHOST];
          arr_curr[i*nwg+nwg-NGHOST+j]   = arr_curr[i*nwg+nwg-NGHOST-1];
          arr_curr[j*nwg+i]              = arr_curr[NGHOST*nwg+i];
          arr_curr[(nwg-NGHOST+j)*nwg+i] = arr_curr[(nwg-NGHOST-1)*nwg+i];
        }
      }
    }

    // Swap pointers
    tarr     = arr_prev;
    arr_prev = arr_curr;
    arr_curr = tarr;

    iter++;
    if (iter%OUTFREQ == 0 || err < eps) {
      std::cout << "Iter. " << std::setw(8) << iter <<
        ", err = " << std::scientific <<  err << std::endl;
    }
  }

  // Print timing this doesn't work on my mac so commenting out until moving to adroit
  auto stop = std::chrono::steady_clock::now();
  std::chrono::duration<double>  durr = stop - strt;
  std::cout << "Finished in " << durr.count() << " s" << std::endl;

  // Final time step output
  if (FOUT) {
    std::ofstream arr_file("final_phi.csv"); // modified to export .csv bc easier4me
    if (arr_file.is_open()) {
     // arr_file << n << ", " << NGHOST << ", " << eps << ", " << iter
      //  << ", " << err << std::endl;
      for (i = NGHOST; i < nwg-NGHOST; ++i) {
        for (j = NGHOST; j < nwg-NGHOST; ++j) {
          arr_file << arr_prev[i*nwg+j] << ", ";
        }
        arr_file << std::endl;
      }
      arr_file.close();
    } else {
      std::cout << "Can't open file for final phi output";
    }
    // write runtime output
    std::ofstream timing_file("time.csv", std::ios::app);
    if (timing_file.is_open()) {

      timing_file << nthreads << ", " << durr.count() << ", " << std::endl;
      timing_file.close();

    } else {
      std::cout << "Can't open file for timing output";
    }
  }

  // Free arrays
  free(arr_curr);
  free(arr_prev);

  return 0;
}
