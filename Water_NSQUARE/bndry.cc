#ifndef __WATER_BNDRY_IMPL__
#define __WATER_BNDRY_IMPL__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include "data.h"

#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"

// 	this routine puts the molecules back inside the box if they are out
//
//	Note: on pthread parallelization
//	BNDRY is called within redsync only, and the callsite is protected with pthread_mutex_lock
//	as a result,the uses of VAR (global variable) within BNDRY is safe
//
void BNDRY(int threadID)
{
    int mol, dir; 
    
    // for each molecule
    for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
    {
        // for each direction
        for ( dir = XDIR; dir <= ZDIR; dir++ ) 
        {
            // if the oxygen atom is out of the box
            if (VAR[mol].F[DISP][dir][O] > BOXL) 
            {
                // move all three atoms back in the box
                int coun = 0;
                INIT_TRANSACTIONS();
                BEGIN_TRANSACTION();
                if(coun!=0)
                {	
                    printf("aborted 10\n");fflush(stdout);
                }
                coun++;

                //pthread_mutex_lock(&var_lock);
                VAR[mol].F[DISP][dir][H1] -= BOXL;
                VAR[mol].F[DISP][dir][O]  -= BOXL;
                VAR[mol].F[DISP][dir][H2] -= BOXL;                
                //pthread_mutex_unlock(&var_lock);

                COMMIT_TRANSACTION();
                coun =0;
            }
            else 
                if (VAR[mol].F[DISP][dir][O] < 0.00) 
                {
                    int coun = 0;
                    INIT_TRANSACTIONS();
                    BEGIN_TRANSACTION();
                    if(coun!=0)
                    {	
                        printf("aborted 11\n");fflush(stdout);
                    }
                    coun++;	

                    //pthread_mutex_lock(&var_lock);
                    VAR[mol].F[DISP][dir][H1] += BOXL;
                    VAR[mol].F[DISP][dir][O]  += BOXL;
                    VAR[mol].F[DISP][dir][H2] += BOXL;
                    //pthread_mutex_unlock(&var_lock);

                    COMMIT_TRANSACTION();
                    coun =0;
                }// end of if
        } // for dir
    } // for mol

} // end of subroutine BNDRY


#endif //__WATER_BNDRY_IMPL__
