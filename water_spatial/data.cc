#ifndef __TSP_DATA_IMPL__
#define __TSP_DATA_IMPL__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

#include "data.h"


int NPROCS 			= 2;			// current Number of Processors
int Tmk_page_size  		= 16384;		// default TMK_page_size is 16K

// pthread pool
	// pthread context pool
	pthread_t pthreads[MAX_PROCS];
	
	// pthread lock pool
	pthread_mutex_t locks[NLOCKS] = {PTHREAD_MUTEX_INITIALIZER};
	
	// pthread condition variable pool
	pthread_cond_t conds[NCONDS] = {PTHREAD_COND_INITIALIZER};

// pthread data used for barrier
pthread_mutex_t mutex_arrival[NBARRIERS] = {PTHREAD_MUTEX_INITIALIZER};
pthread_cond_t cond_arrival[NBARRIERS]	 = {PTHREAD_COND_INITIALIZER };
int arrived[NBARRIERS]  	  = {0};	// this is pthread_barrier()'s internal data

//pthread_mutex_t print_lock    		= PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_t mol_lock[MAXMOLS]   	= {PTHREAD_MUTEX_INITIALIZER};
//pthread_mutex_t IntrafVirLock		= PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_t InterfVirLock		= PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_t KinetiSumLock		= PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_t PotengSumLock		= PTHREAD_MUTEX_INITIALIZER;

//pthread_mutex_t gl_lock				= PTHREAD_MUTEX_INITIALIZER;	// lock to guard all write access to gl (GlobalMemory)
//pthread_mutex_t var_lock			= PTHREAD_MUTEX_INITIALIZER;	// lock to guard all write access to global memory VAR


struct timeval start, end;


///////////////////////////////////////////////////////////////
//
//	a default barrier, without taking a barrier id
//
//////////////////////////////////////////////////////////////
void pthread_barrier(int barrier_index){

	// check the barrier_index is valid: within range
	if(barrier_index > NBARRIERS || barrier_index < 0){
		fprintf(stderr," barrier index %d out of range [0..%d]\b",
				barrier_index, NBARRIERS);
		fflush(stderr);
		exit(-1);					
	}

	pthread_mutex_lock(&mutex_arrival[barrier_index]);

	arrived[barrier_index]++;
	if(arrived[barrier_index] < NPROCS){
		pthread_cond_wait(&cond_arrival[barrier_index], &mutex_arrival[barrier_index]);
	}
	else{
		pthread_cond_broadcast(&cond_arrival[barrier_index]);
		arrived[barrier_index] = 0; // ready for next barrier
	}

	pthread_mutex_unlock(&mutex_arrival[barrier_index]);

}

#endif //__TSP_DATA_IMPL__
