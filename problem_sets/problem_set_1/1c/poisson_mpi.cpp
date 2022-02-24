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
#define MAXITER 100000
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
  nwg = n + 2*NGHOST;
  h   = 1./n; //gridcell size
  np = (n/nthreads); // size of parallel block
  npwg = np + 2*NGHOST;

  std::cout << "npwg = "<< npwg << " np = " << np << " nwg = " << nwg << std::endl;
  
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
  double  err  = std::numeric_limits<double>::max();
  double  ffac = .25 * pow(h, 2.);
  auto    strt = std::chrono::steady_clock::now();

  // initialize MPI process
  MPI_Init(NULL, NULL); 
  
  // Allocate and fill the working array and previous step array
  double* arr_final = (double*) calloc(n*n, sizeof(double));
  double* arr_curr = (double*) calloc(npwg*nwg, sizeof(double));
  double* arr_prev = (double*) calloc(npwg*nwg, sizeof(double));      
  double* agg_arr = (double*) calloc(np*n, sizeof(double));
  double* agg_send = (double*) calloc(np*n, sizeof(double));
  double* tosendup = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* tosenddown = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* torecvup = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* torecvdown = (double*) calloc(nwg*NGHOST, sizeof(double));
  double* tarr = NULL;
  double terr = 0.;
  double err_m = 0.;
  int tn = 0;
  int world_size = 0;
  int iter = 0;
  
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &tn);
  
  //  std::cout << "I am thread: " << tn << std::endl;
  
    while (err > eps && iter < MAXITER) {
 
      err_m = 0.;
     
	for (i = NGHOST; i < npwg-NGHOST; ++i) { // only calculate values for interior cells
	  for (j = NGHOST; j < nwg-NGHOST; ++j) {
	    arr_curr[(i*nwg)+j] =
	      0.25 * (arr_prev[(i-1)*nwg+j]     + arr_prev[(i+1)*nwg+j] +
		      arr_prev[i    *nwg+(j-1)] + arr_prev[i    *nwg+(j+1)]) -
	      ffac * fval[tn*np*n+(i-NGHOST)*n+(j-NGHOST)];
	    terr = std::abs(arr_curr[i*nwg+j] - arr_prev[i*nwg+j]);
	    err_m = (terr > err_m) ? terr: err_m;

	  }
	}
      
      for (i = NGHOST; i < npwg-NGHOST; ++i) { // only calculate values for interior cells
	  for (j = 0; j < NGHOST; ++j) {
	    arr_curr[i*nwg+j]              = arr_curr[i*nwg+NGHOST];
	    arr_curr[i*nwg+nwg-NGHOST+j]   = arr_curr[i*nwg+nwg-NGHOST-1];
	    
	    // only fill in top and bottom boundary conditions if top or bottom cell
	    if (tn == 0) {
	      arr_curr[j*nwg+i]              = arr_curr[NGHOST*nwg+i];
	    } else if (tn != 0 & tn == world_size-1) {
	      arr_curr[(nwg-NGHOST+j)*nwg+i] = arr_curr[(nwg-NGHOST-1)*nwg+i];
	    }
	    }
	}

      if (nthreads > 1) {
      
      for (i = 0; i<nwg*NGHOST; ++i) {

	tosenddown[i] = arr_curr[(nwg*NGHOST) + i];
	tosendup[i] = arr_curr[((npwg-NGHOST-1)*nwg)+i];
		
      }

      MPI_Barrier(MPI_COMM_WORLD);
      
      // communications
      if (tn == 0) {
	//receive from and send to top only if bottom cell
	MPI_Send(tosendup, nwg, MPI_LONG_DOUBLE, tn+1, tn*10+(tn+1), MPI_COMM_WORLD);
	std::cout << "I, thread: " << tn << "am sending rows: " << nwg*NGHOST << " through " << nwg*NGHOST + nwg*NGHOST << " to thread " << tn+1 << std::endl;
	MPI_Recv(torecvup, nwg, MPI_LONG_DOUBLE, tn+1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
	std::cout << "I, thread: " << tn << "received from thread: " << tn+1 << std::endl;
      } else if (tn == tn !=0 & world_size-1) {
	MPI_Send(tosenddown, nwg, MPI_LONG_DOUBLE, tn-1, tn*10+(tn-1), MPI_COMM_WORLD);
	std::cout << "I, thread: " << tn << " am sending rows: " << (npwg-NGHOST-1)*nwg << " through " << (npwg-NGHOST-1)*nwg + nwg*NGHOST << " to thread " << tn-1 << std::endl; 
	MPI_Recv(torecvdown, nwg, MPI_LONG_DOUBLE, tn-1, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
	std::cout << "I, thread: " << tn << " received from thread: " << tn-1 << std::endl;
      } else {
	//receive and send to both top and bottom if in the middle
	//	MPI_Recv(torecvup, 1, MPI_LONG_DOUBLE, tn+1, (tn+1)*10+tn, MPI_COMM_WORLD, NULL);
	//	MPI_Send(tosendup, 1, MPI_LONG_DOUBLE, tn+1, tn*10+(tn+1), MPI_COMM_WORLD);
	//	MPI_Recv(torecvdown, 1, MPI_LONG_DOUBLE, tn-1, (tn-1)*10+tn, MPI_COMM_WORLD, NULL);
	//	MPI_Send(tosenddown, 1, MPI_LONG_DOUBLE, tn-1, tn*10+(tn-1), MPI_COMM_WORLD);
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
      std::cout << "crashing after barrier " << std::endl;
      
      MPI_Allreduce(&err_m, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //collect errors into max
      std::cout << "crashing after this " << std::endl;

      if (tn == 0) {
      } else if (tn == world_size-1) {
      }
      //reassign edge values
      for (i = 0; i<nwg*NGHOST; ++i) {

	if (tn == world_size-1) {
	  arr_curr[i] = torecvdown[i];
	} else if (tn == 0) {
	  arr_curr[(npwg-NGHOST)*nwg+i] = torecvup[i];
	} else {
	  arr_curr[i] = torecvdown[i];
	  arr_curr[(npwg-NGHOST)*nwg+i] = torecvdown[i];
	}
	
      }
      }
      std::cout << "crashing after this " << std::endl;
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
    
    // MPI_Gather(agg_send, np*n, MPI_LONG_DOUBLE, arr_final, np*n, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(tn == 0){

      for (j=0; j<np*n; ++j) {

	  arr_final[j] = agg_send[j];
	  
	}

      for (i=1; i<world_size; ++i) {
	    MPI_Recv(agg_arr, np*n, MPI_LONG_DOUBLE, i, i, MPI_COMM_WORLD, NULL);
	    
	    for (j = i*n; j<(i*n)+(np*n); ++j) {
	      //	      std::cout << "the issue is accessing arr_final" << agg_arr[90] << std::endl;
	      arr_final[j] = agg_arr[j-(i*n)];
	  
	    }
	    
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
  MPI_Finalize();

  std::cout << "Made it past finalize, arr_final contains: " << arr_final[90] << std::endl;
  
  auto stop = std::chrono::steady_clock::now();
  std::chrono::duration<double>  durr = stop - strt;
  std::cout << "Finished in " << durr.count() << " s" << std::endl;

  // Final time step output
  if (FOUT) {
    std::ofstream arr_file("final_phi.csv"); // modified to export .csv bc easier4me
    if (arr_file.is_open()) {
     // arr_file << n << ", " << NGHOST << ", " << eps << ", " << iter
      //  << ", " << err << std::endl;
      for (i = 0; i < np; ++i) {
        for (j = 0; j < np; ++j) {
          arr_file << arr_final[i*np+j] << ", ";
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

      timing_file << nthreads << ", " << durr.count() << std::endl;
      timing_file.close();

    } else {
      std::cout << "Can't open file for timing output";
    }
  }
  
  return 0;
}
