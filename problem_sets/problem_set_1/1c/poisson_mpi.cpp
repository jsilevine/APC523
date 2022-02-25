#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <chrono>
#include <math.h>
#include <mpi.h>
#include <algorithm>
#include <iterator>

#define NDEF    64
#define NGHOST  1
#define MAXITER 1000000
#define OUTFREQ 100
#define FOUT    1

int main(int argc, char *argv[]) {
  // Initialize variables
  int n, nwg, np, npwg = 1; //np is number of cells in parallel block
  int i, j, t = 0;
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

  auto    strt = std::chrono::steady_clock::now();
  
  // initialize MPI process
  MPI_Init(NULL, NULL); 
  
  nwg = n + 2*NGHOST;
  h   = 1./n; //gridcell size
  np = (n/nthreads); // size of parallel block
  npwg = np + 2*NGHOST;


  // Main loop
  double  err  = std::numeric_limits<double>::max();
  double  ffac = .25 * pow(h, 2.);
  
  //  double* arr_final = (double*) calloc(n*n, sizeof(double));
  double* arr_final = (double*) calloc(n*n, sizeof(double));
  double* arr_curr = (double*) calloc(npwg*nwg, sizeof(double));
  double* arr_prev = (double*) calloc(npwg*nwg, sizeof(double));      
  double* agg_arr = (double*) calloc(np*n, sizeof(double));
  double* agg_send = (double*) calloc(np*n, sizeof(double));
  
  double* tosendup = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* tosenddown = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* torecvup = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* torecvdown = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* tarr = (double*) calloc(npwg*nwg, sizeof(double));
  double terr = 0.;
  double err_m = 0.;
  int tn = 0;
  int world_size = 0;
  int iter = 0;
  
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &tn);

  std::cout << world_size << std::endl;
  
  // Fill in the source function
  double  x, y = 0.;
  double* fval = (double*) malloc(n*np*sizeof(double));
  for (i = (tn*np); i < (tn*np)+np; ++i) {
	x = (i + .5) * h;
	for (j = 0; j < n; ++j) {
	  y = (j + .5) * h;
	  fval[(i-(tn*np))*n+j] =
	    pow(M_PI,2.) * (16.*pow(y,2.)+81.*pow(x,4.))
	    * cos(2*M_PI*pow(y,2.)) * cos(3*M_PI*pow(x,3.))
	    + 18.*M_PI*x * sin(3*M_PI*pow(x,3.)) * cos(2*M_PI*pow(y,2.))
	    + 4. *M_PI   * cos(3*M_PI*pow(x,3.)) * sin(2*M_PI*pow(y,2.));
	}
  }
    while (err > eps && iter < MAXITER) {
      err_m = 0.;
      
	for (i = NGHOST; i < npwg-NGHOST; ++i) { // only calculate values for interior cells
	  for (j = NGHOST; j < nwg-NGHOST; ++j) {
	    arr_curr[(i*nwg)+j] =
	      0.25 * (arr_prev[(i-1)*nwg+j]     + arr_prev[(i+1)*nwg+j] +
		      arr_prev[i    *nwg+(j-1)] + arr_prev[i    *nwg+(j+1)]) -
	      ffac * fval[(i-NGHOST)*n+(j-NGHOST)];
	    terr = std::abs(arr_curr[i*nwg+j] - arr_prev[i*nwg+j]);
	    err_m = (terr > err_m) ? terr: err_m;
	    if (err_m == 0) {
	      std::cout << "err_m = 0 and arr_curr =  " << arr_curr[i*nwg+j] << std::endl;
	    }
	  }
	}
	
      for (i = NGHOST; i < npwg-NGHOST; ++i) { 
	  for (j = 0; j < NGHOST; ++j) {
	    arr_curr[i*nwg+j]              = arr_curr[i*nwg+NGHOST];
	    arr_curr[i*nwg+nwg-NGHOST+j]   = arr_curr[i*nwg+nwg-NGHOST-1];
	    
	    // only fill in top and bottom boundary conditions if top or bottom cell
	    if (tn == 0) {
	      arr_curr[j*nwg+i]              = arr_curr[NGHOST*nwg+i];
	    } else if (tn == world_size-1) {
	      arr_curr[(npwg-NGHOST+j)*nwg+i] = arr_curr[(npwg-NGHOST-1)*nwg+i];
	    }
	    }
	}

      if (world_size > 1) { 
	for (i = 0; i<nwg; ++i) {
	  tosenddown[i] = arr_curr[(nwg*NGHOST) + i];
	  tosendup[i] = arr_curr[((npwg-NGHOST-1)*nwg)+i];
	}      

	MPI_Barrier(MPI_COMM_WORLD);
      
      // communications
      if (tn == 0) {
	MPI_Send(tosendup, nwg, MPI_LONG_DOUBLE, tn+1, 0, MPI_COMM_WORLD);
	MPI_Recv(torecvup, nwg, MPI_LONG_DOUBLE, tn+1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
      } else if (tn == world_size-1) {
	MPI_Recv(torecvdown, nwg, MPI_LONG_DOUBLE, tn-1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
	MPI_Send(tosenddown, nwg, MPI_LONG_DOUBLE, tn-1, 0, MPI_COMM_WORLD);
      } else {
	//	receive and send to both top and bottom if in the middle
	MPI_Recv(torecvdown, nwg, MPI_LONG_DOUBLE, tn-1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
	MPI_Send(tosendup, nwg, MPI_LONG_DOUBLE, tn+1, 0, MPI_COMM_WORLD);

	MPI_Recv(torecvup, nwg, MPI_LONG_DOUBLE, tn+1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
	MPI_Send(tosenddown, nwg, MPI_LONG_DOUBLE, tn-1, 0, MPI_COMM_WORLD);
	
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(&err_m, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //collect errors into max
      MPI_Barrier(MPI_COMM_WORLD);
      
      for (i = 0; i<nwg*NGHOST; ++i) {

	if (tn == world_size-1) {
	  arr_curr[i] = torecvdown[i];
	} else if (tn == 0) {
	  arr_curr[(npwg-NGHOST)*nwg+i] = torecvup[i];
	} else {
	  arr_curr[i] = torecvdown[i];
	  arr_curr[(npwg-NGHOST)*nwg+i] = torecvup[i];
	}	
      }
    }

      if (world_size == 1) {
	err = err_m;
      }
      
      tarr     = arr_prev;
      arr_prev = arr_curr;
      arr_curr = tarr;
      
      MPI_Barrier(MPI_COMM_WORLD);
      iter++;

      if(tn == 0) {
	     
	if (iter%OUTFREQ == 0 || err < eps) {
	std::cout << "Iter. " << std::setw(8) << iter <<
	  ", err = " << std::scientific <<  err << std::endl;
      }
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
  }
    
    for (i = NGHOST; i<npwg-NGHOST; ++i) {
      for (j = NGHOST; j<nwg-NGHOST; ++j) {
	  agg_send[(i-NGHOST)*n+(j-NGHOST)] = arr_curr[i*nwg+j];
	}
    }

    if (tn != 0) {
     MPI_Send(agg_send, np*n, MPI_LONG_DOUBLE, 0, tn, MPI_COMM_WORLD);
    }
    
     MPI_Barrier(MPI_COMM_WORLD);
    
    if(tn == 0){
      for (j=0; j<np*n; ++j) {
	  arr_final[j] = agg_send[j];	  
	}
      for (i=1; i<world_size; ++i) {
	MPI_Recv(agg_arr, np*n, MPI_LONG_DOUBLE, i, i, MPI_COMM_WORLD, NULL);
	    for (j = i*np*n; j<(i*np*n)+(np*n); ++j) {
	      arr_final[j] = agg_arr[j-(i*np*n)];	  
	    }	    
      }
      
  auto stop = std::chrono::steady_clock::now();
  std::chrono::duration<double>  durr = stop - strt;
  std::cout << "Finished in " << durr.count() << " s" << std::endl;

  // Final time step output
  if (FOUT) {
    // write runtime output
    std::ofstream timing_file("time.csv", std::ios::app);
    if (timing_file.is_open()) {

      timing_file << nthreads << ", " << durr.count() << std::endl;
      timing_file.close();

    } else {
      std::cout << "Can't open file for timing output";
    }
  }
    }
    
  MPI_Finalize();
 
  return 0;
}
