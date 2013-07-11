#ifndef __WATER_PREDCOR_IMPL__
#define __WATER_PREDCOR_IMPL__

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

//
// predicts new values for displacement and its five derivatives
// this routine calculates predicted F(X), F'(X), F''(X), ...
// NOR1 = NORDER + 1 = 7 (for a sixth-order method)
//
//	Note: on PThread parallellelization
//	there is only one callsite to PREDIC, inside redsync.c
//	and, the callsite is within a critical section
//	
//	So, it is safe
//
void PREDIC(double C[], int NOR1, int threadID)
{
	int JIZ;
	int  JI;
	int  L;
	double S;
        int func, mol, dir, atom;

        JIZ = 2;
        // .....loop over F(X), F'(X), F''(X), .....
	for (func = 0; func < NORDER; func++) 
	{
		for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++)
	    	for ( dir = 0; dir < NDIR; dir++)
				for ( atom = 0; atom < NATOM; atom++ ) 
				{
		    		        JI = JIZ;
					// sum over Taylor Series
					S = 0.0;
					for ( L = func; L < NORDER; L++) 
					{
						S += C[JI] * VAR[mol].F[L+1][dir][atom];
						JI++;
		    		} // for on L
				
					// updating VAR shared memory, need a lock
				int coun = 0;
				INIT_TRANSACTIONS();
				BEGIN_TRANSACTION();
				if(coun!=0)
				{	
					printf("aborted 6\n");fflush(stdout);
				}
				coun++;
					//pthread_mutex_lock(&var_lock);	
			    		VAR[mol].F[func][dir][atom] += S;
					//pthread_mutex_unlock(&var_lock);
				COMMIT_TRANSACTION();
				coun =0;
				} // for on atom
				JIZ += NOR1;
    } // for func

}// end of subroutine PREDIC

//
// corrects the predicted values, based on forces etc. computed in the interim
//
// double PCC [] the predictor-corrector constants */
// int NOR NORDER + 1 = 7 for a sixth-order method) */
//
/*
.....this routine calculates corrected F(X), F'(X), F"(X), ....
     from corrected F(X) = predicted F(X) + PCC(1)*(FR-SD)
     where SD is predicted accl. F"(X) and FR is computed 
     accl. (force/mass) at predicted position
*/
//	Note: on PThread parallellelization
//	
//
//	
//
void CORREC(double PCC[], int NOR1, int threadID)
{
	double Y;
	int mol, dir, atom, func;

	for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
	{
        //for (mol = 0; mol < NMOL; mol++){
		for (dir = 0; dir < NDIR; dir++) 
		{
		    for (atom = 0; atom < NATOM; atom++) 
			{
				Y = VAR[mol].F[FORCES][dir][atom] - VAR[mol].F[ACC][dir][atom];
				for ( func = 0; func < NOR1; func++)
				{
	
					int coun = 0;
					INIT_TRANSACTIONS();
					BEGIN_TRANSACTION();
						if(coun!=0)
					{	
						printf("aborted 13\n");fflush(stdout);
					}
					coun++;
					//pthread_mutex_lock(&var_lock);	
		    			VAR[mol].F[func][dir][atom] += PCC[func] * Y;   
					//pthread_mutex_unlock(&var_lock);
					COMMIT_TRANSACTION();
					coun =0;
		    		
				}
	    	} // for atom
		} // for dir
	} // for mol

}// end of subroutine CORREC

#endif //__WATER_PREDCOR_IMPL__
