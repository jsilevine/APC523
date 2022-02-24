#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <chrono>
#include <math.h>

#define NDEF    64
#define NGHOST  2
#define MAXITER 100000
#define OUTFREQ 1000
#define FOUT    1

int main(int argc, char *argv[]) {
  // Initialize variables
  int n = 1;
  int i,j    = 0;
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
  } else {
    std::cout << "format is program [n]"] << std::endl;
}
nwg = n + 2*NGHOST;
h   = 1./n; //gridcell size

// initialize arrays
double* grad[2][n*n]; //gradient

  // Fill in the source function
  double  x, y = 0.;
  for (i = 0; i < n; ++i) {
    x = (i + .5) * h;
    for (j = 0; j < n; ++j) {
      y = (j + .5) * h;
      grad[1][] =
        pow(M_PI,2.) * (16.*pow(y,2.)+81.*pow(x,4.))
          * cos(2*M_PI*pow(y,2.)) * cos(3*M_PI*pow(x,3.))
        + 18.*M_PI*x * sin(3*M_PI*pow(x,3.)) * cos(2*M_PI*pow(y,2.))
        + 4. *M_PI   * cos(3*M_PI*pow(x,3.)) * sin(2*M_PI*pow(y,2.));
    }
  }

  // Print timing this doesn't work on my mac so commenting out until moving to adroit
  //auto stop = std::chrono::steady_clock::now();
  //std::chrono::duration<double>  durr = stop - strt;
  //std::cout << "Finished in " << durr.count() << " s" << std::endl;

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
  }

  // Free arrays
  free(arr_curr);
  free(arr_prev);

  return 0;
}
