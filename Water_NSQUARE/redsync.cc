#ifndef __WATER_REDSYNC_IMPL__
#define __WATER_REDSYNC_IMPL__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
//using namespace std;
#include "data.h"

#define PR "%30.10"

#include "comments.h"	// moves leading comments to another file for now

#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "split.h"
#include "global.h"
#include "include/src/tm_basics.h"

//
// routine that implements the time-steps. Called by main routine and calls others
//
double  MDMAIN(int NFSV, int NFRST, int NSTEP, int NRST, int NPRINT, int NSAVE, int LKT, int NORD1, int threadID)
{
    double RES =0.0;
    int i =0;
    double POTA = 0.0, POTR =0.0, POTRF =0.0;
    double XVIR = 0.0, AVGT =0.0, TEN = 0.0;

    // wait till everyone gets to beginning; not necessary
    pthread_barrier(0);

    // MOLECULAR DYNAMICS LOOP OVER ALL TIME-STEPS	
    for (i=1;i <= NSTEP; i++) 
    {
        TTMV[threadID] = TTMV[threadID] + 1.00;		// ? on this
        
        // initialize various shared sums
        if (threadID == 0) 
        {
            int dir;
            INIT_TRANSACTIONS();
            BEGIN_TRANSACTION();		
            //pthread_mutex_lock(&gl_lock);	// gl data initialization
            gl->VIR 	= 0.0;
            gl->POTA 	= 0.0;
            gl->POTR 	= 0.0;
            gl->POTRF 	= 0.0;
            for (dir = XDIR; dir <= ZDIR; dir++){
                gl->SUM[dir] = 0.0;
            }
            //pthread_mutex_unlock(&gl_lock);
            COMMIT_TRANSACTION();
        }

        pthread_barrier(gl->start);
        
        PREDIC(TLC, NORD1, threadID);
        INTRAF(&gl->VIR, threadID);
        
        pthread_barrier(gl->start);
        
        INTERF(FORCES,&gl->VIR,  threadID);
        
        pthread_barrier(gl->start);

        //----- some hardcoding from INTERF() -----
        int mol,dir;
        for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
        {
            for ( dir = XDIR; dir  <= ZDIR; dir++) 
            {
                int coun = 0;
                INIT_TRANSACTIONS();
                BEGIN_TRANSACTION();
                if(coun!=0)
                {	
                    printf("aborted 14\n");fflush(stdout);
                }
                coun++;
                
                //pthread_mutex_lock(&var_lock);                
                VAR[mol].F[FORCES][dir][H1] = VAR[mol].F[FORCES][dir][H1] * FHM;
                VAR[mol].F[FORCES][dir][H2] = VAR[mol].F[FORCES][dir][H2] * FHM;
                VAR[mol].F[FORCES][dir][O]  = VAR[mol].F[FORCES][dir][O]  * FOM;
                //pthread_mutex_unlock(&var_lock);
                
                COMMIT_TRANSACTION();
                coun=0;
            } // for dir
        } // for mol
        //----- Some hardcoding from INTERF() -----


		
        CORREC(PCC, NORD1, threadID);
        BNDRY(threadID);
        KINETI(NMOL, gl->SUM, HMAS, OMAS, threadID); // this (access to gl->SUM) is good
        
        pthread_barrier(gl->start);

        // --> looking at this now
        TKIN[threadID] = TKIN[threadID] + gl->SUM[0] + gl->SUM[1] + gl->SUM[2];
        TVIR[threadID] = TVIR[threadID] - gl->VIR;
        
        pthread_barrier(0);
        
        // check if potential energy is to be computed, and if
        // printing and/or saving is to be done, this time step.
        // Note that potential energy is computed once every NPRINT
        // time-steps
        if (((i % NPRINT) == 0) || ( (NSAVE > 0) && ((i % NSAVE) == 0)))
        {
            //  call potential energy computing routine
            POTENG(&gl->POTA, &gl->POTR, &gl->POTRF, threadID);	// this used of gl->... data is good
            
            pthread_barrier(gl->start);	// this barrier can cause problems on non-entry threads

            // modify computed sums
            POTA  	= gl->POTA * FPOT;
            POTR  	= gl->POTR * FPOT;
            POTRF 	= gl->POTRF* FPOT;
            
            // compute some values to print
            XVIR 	= TVIR[threadID] * FPOT * 0.50 / TTMV[threadID];
            AVGT 	= TKIN[threadID] * FKIN * TEMP * 2.00/(3.00 * TTMV[threadID]);
            TEN 	= (gl->SUM[0] + gl->SUM[1] + gl->SUM[2]) * FKIN;
            RES 	= POTA + POTR + POTRF + TEN;
            
            pthread_barrier(gl->start);	
            
            if ((i % NPRINT == 0) && (threadID == 0)) 
            {
                pthread_mutex_lock(&print_lock);
                fprintf(six,"\t%5d "PR"lf "PR"lf "PR"lf "PR"lf\n %16.3lf "PR"lf  "PR"lf\n",
                    		i,TEN,POTA,POTR,POTRF,RES,AVGT,XVIR);
                fflush(six);
                pthread_mutex_unlock(&print_lock);
            }// if (i % NPRINT)
	}// end big if


        // wait for everyone to finish time-step
        pthread_barrier(gl->start);

    }// for loop on "i"

    tempXTT[threadID] = RES;
    
    return RES;
    
}// end of subroutine MDMAIN

#endif //__WATER_REDSYNC_IMPL__
