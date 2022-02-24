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
  double x,y = 0.;
  double h   = 1.;
  double eps = 1.e-6;

  // Parse input
  if (argc == 1) {
    h = NDEF;
  } else if (argc == 2) {
    h = atof(argv[1]);
  } else {
    std::cout << "format is program [n]" << std::endl;
  }

  //calculate n
  n = 1./h;

  // initialize arrays
  double grad[2][(n-4)*(n-4)]; //gradient
  double true_grad[2][(n-4)*(n-4)];
  double phi[n*n];
  double L2[2];
  L2[0] = 0.;
  L2[1] = 0.;

// generate phi
for (i = 0; i<n; ++i) {
  x = (i + .5) * h;
  for (j=0; j<n; ++j) {
    y = (j + .5) * h;
    phi[i*n+j] = sin(4. * M_PI * pow(x,2.)) * cos(2. * M_PI * pow(y, 3.));
  }
 }
 
  // loop over interior
  for (int i = 2; i < n-2; ++i) {
    for (int j = 2; j < n-2; ++j) { 
      grad[0][(i-2)*(n-4) + (j-2)] =
	(-phi[(i+2)*n + j] + (8. * phi[(i+1)*n + j]) - (8. * phi[(i-1)*n + j]) + phi[(i-2)*n + j]) / (12. * h);
      grad[1][(i-2)*(n-4) + (j-2)] =
	(-phi[(i*n) + j-2] + (8. * phi[(i*n) + j-1]) - (8. * phi[(i*n) + j+1]) + phi[(i*n)+ j+2]) / (12. * h);
      }
  } 
 
  // Calculate true gradient
  for (i = 0; i < (n-4); ++i) {
    x = ((i+2) + .5) * h;
    for (j = 0; j < (n-4); ++j) {
      y = ((j+2) + .5) * h;
      true_grad[0][i*(n-4)+j] =
	2. * M_PI * (4. * x * cos(4. * M_PI * pow(x, 2.)) * cos(-2. * M_PI * pow(y, 3.)));
      true_grad[1][i*(n-4)+j] =
	2. * M_PI * (3. * pow(y,2.) * sin(4. * M_PI * pow(x,2.)) * sin(-2. * M_PI * pow(y,3.)));
    }
  }
  
// Calculate L2 norm
for (i = 0; i<n-4; ++i) {
  for (j=0; j<n-4; ++j) {

    L2[0] += pow(std::abs(true_grad[0][i*(n-4)+j] - grad[0][i*(n-4)+j]), 2.);
    L2[1] += pow(std::abs(true_grad[1][i*(n-4)+j] - grad[1][i*(n-4)+j]), 2.);

  }
 }

 L2[0] = sqrt(L2[0] * pow(h, 2.));
 L2[1] = sqrt(L2[1] * pow(h, 2.));
 
 // Final time step output
 if (FOUT) {
   std::ofstream arr_file("L2_h.csv", std::ios::app); // modified to export .csv bc easier4me
   if (arr_file.is_open()) {
     arr_file << h << ", " << L2[0] << ", " << L2[1] << std::endl;
     arr_file.close();
   } else {
     std::cout << "Can't open file for final phi output";
   }
 }
 
 return 0;
}
