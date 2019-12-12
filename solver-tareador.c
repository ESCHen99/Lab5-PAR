#include "heat.h"
#include <stdio.h>
#include <tareador.h>

/*
 * Function to copy one matrix into another
 */

void copy_mat (double *u, double *v, unsigned sizex, unsigned sizey)
{
    for (int i=1; i<= sizex-2; i++)
        for (int j=1; j<= sizey-2; j++) 
            v[ i*sizey+j ] = u[ i*sizey+j ];
}

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
    double diff, sum=0.0;
  
    int howmany=1;
 
    for (int blockid = 0; blockid < howmany; ++blockid) {
  //    tareador_start_task("relax_jacobi");
	  int i_start = lowerb(blockid, howmany, sizex);
      int i_end = upperb(blockid, howmany, sizex);
      for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
//        tareador_start_task("relax_jacobi_i");
		tareador_disable_object(&sum);
		for (int j=1; j<= sizey-2; j++) {
	     tareador_start_task("relax_jacobi_j");
		 utmp[i*sizey+j]= 0.25 * ( u[ i*sizey     + (j-1) ]+  // left
	                               u[ i*sizey     + (j+1) ]+  // right
				       u[ (i-1)*sizey + j     ]+  // top
				       u[ (i+1)*sizey + j     ]); // bottom
	     diff = utmp[i*sizey+j] - u[i*sizey + j];
	     sum += diff * diff; 
		 tareador_end_task("relax_jacobi_j");
	 	}
		tareador_enable_object(&sum);
	//	tareador_end_task("relax_jacobi_i");
      }
	 // tareador_end_task("relax_jacobi");
    }

    return sum;
}

/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
	tareador_disable_object(&sum);
    int howmany=1;
    for (int blockid = 0; blockid < howmany; ++blockid) {
      //tareador_start_task("relax_gauss");	
      int i_start = lowerb(blockid, howmany, sizex);
      int i_end = upperb(blockid, howmany, sizex);
      for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
        tareador_start_task("relax_gauss_i");
		for (int j=1; j<= sizey-2; j++) {
	    //tareador_start_task("relax_gauss_j");
		unew= 0.25 * ( u[ i*sizey	+ (j-1) ]+  // left
			   u[ i*sizey	+ (j+1) ]+  // right
			   u[ (i-1)*sizey	+ j     ]+  // top
			   u[ (i+1)*sizey	+ j     ]); // bottom
	    diff = unew - u[i*sizey+ j];
	    sum += diff * diff; 
	    u[i*sizey+j]=unew; 
		//tareador_end_task("relax_gauss_j");
        }
       tareador_end_task("relax_gauss_i");
	  }
	  //tareador_end_task("relax_gauss");
    }
	tareador_enable_object(&sum);
    return sum;
}
