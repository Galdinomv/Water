#ifndef __TSP_DATA_HEADER__
#define __TSP_DATA_HEADER__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

#include "global.h"

#define MAX_PROCS 		16

//#define NPROCS 		2
#define NBARRIERS 	32			// max 32 Barriers				
#define NCONDS 		4			// max 4  Conditions			
#define NLOCKS      4 * 1024 	// max 4K Locks					

extern int  NPROCS 	  ;
extern int	Tmk_page_size;		// default TMK_page_size is 16K 			

// pthread pool
	// pthread context pool
	extern pthread_t pthreads[MAX_PROCS];
	
	// pthread lock pool
	extern pthread_mutex_t locks[NLOCKS];
	
	// pthread condition variable pool
	extern pthread_cond_t conds[NCONDS];

// pthread data used for barrier
extern pthread_mutex_t mutex_arrival[NBARRIERS];
extern pthread_cond_t cond_arrival[NBARRIERS]  ;

//extern pthread_mutex_t print_lock   ;
//extern pthread_mutex_t mol_lock[MAXMOLS];
//extern pthread_mutex_t IntrafVirLock;
//extern pthread_mutex_t InterfVirLock;
//extern pthread_mutex_t KinetiSumLock;
//extern pthread_mutex_t PotengSumLock;

// global (shared) memory access lock
//extern pthread_mutex_t gl_lock		;	// lock to guard all write access to gl (GlobalMemory)
//extern pthread_mutex_t var_lock		;	// lock to guard all write access to global memory VAR


typedef struct pthread_param_struct_t{
	int threadID;

	// place more thread data here
	// ...
	
}ParamStruct;

extern struct timeval start, end;		// start time, stop time
extern int debug ;							// debugging flag



// declarations of function prototypes
extern void pthread_barrier(int);









#endif //__TSP_DATA_HEADER__
