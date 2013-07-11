#ifndef __WATER_KINETI_IMPL__
#define __WATER_KINETI_IMPL__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "data.h"

#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

/* 	this routine computes kinetic energy in each of the three
  	spatial dimensions, and puts the computed values in the
	SUM array */ 
//	Notes on PThread parallelization
//	- all non-threadID params (NMOL, SUM[], HMAS, OMAS) are global shard variables
//	- 
//	- all shared data access have been guarded with locks (good)
//	- 
void KINETI(int NMOL, tm_double SUM[], double HMAS, double OMAS, int threadID)
{
    int dir, mol;
    double S;
    tm_double *tempptr; 
    
    // loop over the three directions
    for (dir = XDIR; dir <= ZDIR; dir++) 
    {
        S = 0.0;
        // loop over the molecules
        for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
        {
            tempptr = VAR[mol].F[VEL][dir]; 
            S += ( tempptr[H1] * tempptr[H1] + tempptr[H2] * tempptr[H2] ) * HMAS
                + (tempptr[O] * tempptr[O]) * OMAS;
        }// for mol

        int coun = 0;
        INIT_TRANSACTIONS();
        BEGIN_TRANSACTION();
        if(coun!=0)
        {	
            printf("aborted 2\n");fflush(stdout);
        }
        coun++;
        //pthread_mutex_lock(&KinetiSumLock);
        SUM[dir] += S;
        //pthread_mutex_unlock(&KinetiSumLock);
        COMMIT_TRANSACTION();
        coun =0;
    } // for "dir"

} // end of subroutine KINETI

#endif //__WATER_KINETI_IMPL__
