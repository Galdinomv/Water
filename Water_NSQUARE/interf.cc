#ifndef __WATER_INTERF_IMPL__
#define __WATER_INTERF_IMPL__

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include "data.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"
#include "fileio.h"

#define	LOCAL_COMP

void UPDATE_FORCES(int , int , int , double [], double [], double [], double [], non_tm_molecule_type []);
void INTERF(int, tm_double *, int);

/* 	This routine gets called both from main() and from mdmain().
	When called from main(), it is used to estimate the initial
	accelerations by computing intermolecular forces.  When called
	from mdmain(), it is used to compute intermolecular forces.
	The parameter DEST specifies whether results go into the 
	accelerations or the forces. Uses routine UPDATE_FORCES in this
	file, and routine CSHIFT in file cshift.U 
*/

/*
	this routine calculates inter-molecular interaction forces
	the distances are arranged in the order  M-M, M-H1, M-H3, H1-M,
	H3-M, H1-H3, H1-H1, H3-H1, H3-H3, O-O, O-H1, O-H3, H1-O, H3-O, 
	where the M are "centers" of the molecules.
*/
//
//	Note: on pthread parallelization
//		there are 2 calls for INTERF, from within redsync and water.c
//	in redsync.c: INTERF is within critical section,
//	in water.c:  INTERF is within the data initialization phase
//	so, in both cases, VAR accesses inside INTERF are safe
//
void INTERF(int DEST, tm_double * VIR, int threadID)
{
        // privatize localVAR[]
        #ifdef	LOCAL_COMP
            non_tm_molecule_type	localVAR[MAXMOLS] = {0};
        #endif	// LOCAL_COMP
        
        int mol, comp, dir, icomp;
        int comp_last, half_mol;
        int KC, K;
        double XL[15], YL[15], ZL[15], RS[15], FF[15], RL[15]; 
        // per- interaction arrays that hold some computed distances
        double FTEMP;
        double LVIR = 0.0;
        
        #ifdef	LOCAL_COMP
            memset(localVAR, 0, NMOL*sizeof(non_tm_molecule_type));
        #endif	// LOCAL_COMP
        
        half_mol = NMOL / 2;
        for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
        {
            comp_last = mol + half_mol;
            
            if (NMOL%2==0 && ((!(mol%2) && (mol < half_mol)) || ((mol%2) && mol > half_mol))) 
            {
                comp_last -= 1;
            }
            
            for (icomp = mol+1; icomp <= comp_last; icomp++) 
            {
                comp = icomp;
                if (comp > NMOL1) 
                {
                    comp = comp % NMOL;
                }
			
                // compute some intermolecular distances
                // only write into XL[], which is local
                CSHIFT( VAR[mol].F[DISP][XDIR], VAR[comp].F[DISP][XDIR],
                        VAR[mol].VM[XDIR], VAR[comp].VM[XDIR], XL, BOXH, BOXL);
                
                // only write into YL[], which is local
                CSHIFT( VAR[mol].F[DISP][YDIR],VAR[comp].F[DISP][YDIR],
                        VAR[mol].VM[YDIR],VAR[comp].VM[YDIR], YL, BOXH, BOXL);
                
                // only write into ZL[], which is local
                CSHIFT( VAR[mol].F[DISP][ZDIR],VAR[comp].F[DISP][ZDIR],
                        VAR[mol].VM[ZDIR],VAR[comp].VM[ZDIR],ZL,BOXH,BOXL);
                
                KC=0;
                for (K = 0; K < 9; K++)
                {
                    RS[K] = XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K];
                    if (RS[K] > CUT2){
                        KC++;
                    }
                } // for loop K

                if (KC != 9)
                {
                    for (K = 0; K < 14; K++)
                    {
                        FF[K]=0.0;
                    }
                    if (RS[0] < CUT2) 
                    {
                        FF[0]=QQ4/(RS[0]*sqrt(RS[0]))+REF4;
                        LVIR = LVIR + FF[0]*RS[0];
                    } // if
      		
                    for (K = 1; K < 5; K++) 
                    {
                        if (RS[K] < CUT2) 
                        { 
                            FF[K] = -QQ2/(RS[K]*sqrt(RS[K]))-REF2;
                            LVIR  = LVIR + FF[K]*RS[K];
                        } // if
                        
                        if (RS[K+4] <= CUT2) 
                        { 
                            RL[K+4] = sqrt(RS[K+4]);
                            FF[K+4] = QQ/(RS[K+4]*RL[K+4])+REF1;
                            LVIR    = LVIR + FF[K+4]*RS[K+4];
                        } // if
                    } // for loop on "K"
    
                    if (KC == 0) 
                    {
                        RS[9] = XL[9]*XL[9]+YL[9]*YL[9]+ZL[9]*ZL[9];
                        RL[9] = sqrt(RS[9]);
                        FF[9] = AB1*exp(-B1*RL[9])/RL[9];
                        LVIR  = LVIR + FF[9]*RS[9];
		    	
                        for (K = 10; K < 14; K++) 
                        { 
                            FTEMP 	= AB2*exp(-B2*RL[K-5])/RL[K-5];
                            FF[K-5] = FF[K-5]+FTEMP;
                            LVIR	= LVIR+FTEMP*RS[K-5];
                            RS[K]	= XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K];
                            RL[K]	= sqrt(RS[K]);
                            FF[K]	= (AB3*exp(-B3*RL[K])-AB4*exp(-B4*RL[K]))/RL[K];
                            LVIR 	= LVIR + FF[K]*RS[K];
                        } // for K
      		
                    } // if KC == 0

                    UPDATE_FORCES(DEST, mol, comp, XL, YL, ZL, FF, localVAR);
                
                }  // if KC != 9
            } // for icomp loop
        } // for mol loop

#ifdef	LOCAL_COMP
    {
        int i, start = StartMol[threadID];	// <- this is a place needs to change
        
        for (i = 0; i < NMOL; i++)
        {
            static int	times = 0;            
            mol = (start + i) % NMOL;
	    
            if (localVAR[mol].F[DEST][XDIR][O]  ||
                localVAR[mol].F[DEST][XDIR][H1] ||
                localVAR[mol].F[DEST][XDIR][H2] ||
                localVAR[mol].F[DEST][YDIR][O]  ||
                localVAR[mol].F[DEST][YDIR][H1] ||
                localVAR[mol].F[DEST][YDIR][H2] ||
                localVAR[mol].F[DEST][ZDIR][O]  ||
                localVAR[mol].F[DEST][ZDIR][H1] ||
                localVAR[mol].F[DEST][ZDIR][H2])
            {
                int coun = 0;
                INIT_TRANSACTIONS();
                BEGIN_TRANSACTION();
                if(coun!=0)
                {	
                    printf("aborted 5\n");fflush(stdout);
                }
                coun++;
    
                //pthread_mutex_lock(&mol_lock[mol]);
                VAR[mol].F[DEST][XDIR][O]  += localVAR[mol].F[DEST][XDIR][O];
                VAR[mol].F[DEST][XDIR][H1] += localVAR[mol].F[DEST][XDIR][H1];
                VAR[mol].F[DEST][XDIR][H2] += localVAR[mol].F[DEST][XDIR][H2];
                VAR[mol].F[DEST][YDIR][O]  += localVAR[mol].F[DEST][YDIR][O];
                VAR[mol].F[DEST][YDIR][H1] += localVAR[mol].F[DEST][YDIR][H1];
                VAR[mol].F[DEST][YDIR][H2] += localVAR[mol].F[DEST][YDIR][H2];
                VAR[mol].F[DEST][ZDIR][O]  += localVAR[mol].F[DEST][ZDIR][O];
                VAR[mol].F[DEST][ZDIR][H1] += localVAR[mol].F[DEST][ZDIR][H1];
                VAR[mol].F[DEST][ZDIR][H2] += localVAR[mol].F[DEST][ZDIR][H2];
                //pthread_mutex_unlock(&mol_lock[mol]);
    
                COMMIT_TRANSACTION();
                coun =0;
            }// end of long conditional if
        }// end of for "i" loop
    }
#endif	//LOCAL_COMP

    //  accumulate the running sum from private per-interaction partial sums
    int coun = 0;
    INIT_TRANSACTIONS();
    BEGIN_TRANSACTION();
	
    if(coun!=0)
    {	
        printf("aborted 1\n");fflush(stdout);
    }
    
    coun++;	
    //pthread_mutex_lock(&InterfVirLock);
    *VIR = *VIR + LVIR;
    //pthread_mutex_unlock(&InterfVirLock);
    COMMIT_TRANSACTION();
    coun =0;
    // wait till all forces are updated
}// end of subroutine INTERF


/* 	
 * from the computed distances etc., compute the intermolecular forces and 
 * update the force (or acceleration) locations 
 * 
 */
void UPDATE_FORCES(int DEST, int mol, int comp, double XL[], double YL[], double ZL[], double FF[], non_tm_molecule_type	localVAR[MAXMOLS]){
	int K;
	double G110[3], G23[3], G45[3], TT[3], TT1[3], TT2[3];
	double GG[15][3];

	// CALCULATE X-COMPONENT FORCES
	for (K = 0; K < 14; K++)  
	{
		GG[K+1][XDIR] = FF[K]*XL[K];
		GG[K+1][YDIR] = FF[K]*YL[K];
		GG[K+1][ZDIR] = FF[K]*ZL[K];
	}

	G110[XDIR] 	= 	GG[10][XDIR]+GG[1][XDIR]*C1;
	G110[YDIR] 	= 	GG[10][YDIR]+GG[1][YDIR]*C1;
	G110[ZDIR] 	= 	GG[10][ZDIR]+GG[1][ZDIR]*C1;
	G23[XDIR] 	= 	GG[2][XDIR]+GG[3][XDIR];
	G23[YDIR] 	= 	GG[2][YDIR]+GG[3][YDIR];
	G23[ZDIR] 	= 	GG[2][ZDIR]+GG[3][ZDIR];
	G45[XDIR]	=	GG[4][XDIR]+GG[5][XDIR];
	G45[YDIR]	=	GG[4][YDIR]+GG[5][YDIR];
	G45[ZDIR]	=	GG[4][ZDIR]+GG[5][ZDIR];
	TT1[XDIR] 	=	GG[1][XDIR]*C2;
	TT1[YDIR] 	=	GG[1][YDIR]*C2;
	TT1[ZDIR] 	=	GG[1][ZDIR]*C2;
	TT[XDIR] 	=	G23[XDIR]*C2+TT1[XDIR];
	TT[YDIR] 	=	G23[YDIR]*C2+TT1[YDIR];
	TT[ZDIR] 	=	G23[ZDIR]*C2+TT1[ZDIR];
	TT2[XDIR]	=	G45[XDIR]*C2+TT1[XDIR];
	TT2[YDIR]	=	G45[YDIR]*C2+TT1[YDIR];
	TT2[ZDIR]	=	G45[ZDIR]*C2+TT1[ZDIR];
    
	// lock locations for the molecule to be updated
	#ifdef	LOCAL_COMP
		// if LOCAL_COMP is defined
		localVAR[mol].F[DEST][XDIR][O]  += G110[XDIR] + GG[11][XDIR] +GG[12][XDIR]+C1*G23[XDIR];

		localVAR[mol].F[DEST][XDIR][H1] += GG[6][XDIR]+GG[7][XDIR]+GG[13][XDIR]+TT[XDIR]+GG[4][XDIR];

	    localVAR[mol].F[DEST][XDIR][H2] += GG[8][XDIR]+GG[9][XDIR]+GG[14][XDIR]+TT[XDIR]+GG[5][XDIR];

	    localVAR[mol].F[DEST][YDIR][O]  += G110[YDIR]+GG[11][YDIR]+GG[12][YDIR]+C1*G23[YDIR];

	    localVAR[mol].F[DEST][YDIR][H1] += GG[6][YDIR]+GG[7][YDIR]+GG[13][YDIR]+TT[YDIR]+GG[4][YDIR];

	    localVAR[mol].F[DEST][YDIR][H2] += GG[8][YDIR]+GG[9][YDIR]+GG[14][YDIR]+TT[YDIR]+GG[5][YDIR];

	    localVAR[mol].F[DEST][ZDIR][O]  += G110[ZDIR]+GG[11][ZDIR]+GG[12][ZDIR]+C1*G23[ZDIR];

	    localVAR[mol].F[DEST][ZDIR][H1] += GG[6][ZDIR]+GG[7][ZDIR]+GG[13][ZDIR]+TT[ZDIR]+GG[4][ZDIR];

	    localVAR[mol].F[DEST][ZDIR][H2] += GG[8][ZDIR]+GG[9][ZDIR]+GG[14][ZDIR]+TT[ZDIR]+GG[5][ZDIR];
	    
	    localVAR[comp].F[DEST][XDIR][O] += -G110[XDIR]-GG[13][XDIR]-GG[14][XDIR]-C1*G45[XDIR];

	    localVAR[comp].F[DEST][XDIR][H1] += -GG[6][XDIR]-GG[8][XDIR]-GG[11][XDIR]-TT2[XDIR]-GG[2][XDIR];

	    localVAR[comp].F[DEST][XDIR][H2] +=
	    	-GG[7][XDIR]-GG[9][XDIR]-GG[12][XDIR]-TT2[XDIR]-GG[3][XDIR];

	    localVAR[comp].F[DEST][YDIR][O] += 
	    	-G110[YDIR]-GG[13][YDIR]-GG[14][YDIR]-C1*G45[YDIR];

	    localVAR[comp].F[DEST][YDIR][H1] += 
	    	-GG[6][YDIR]-GG[8][YDIR]-GG[11][YDIR]-TT2[YDIR]-GG[2][YDIR];

	    localVAR[comp].F[DEST][YDIR][H2] += 
	    	-GG[7][YDIR]-GG[9][YDIR]-GG[12][YDIR]-TT2[YDIR]-GG[3][YDIR];

	    localVAR[comp].F[DEST][ZDIR][O] +=
	    	-G110[ZDIR]-GG[13][ZDIR]-GG[14][ZDIR]-C1*G45[ZDIR];

	    localVAR[comp].F[DEST][ZDIR][H1] +=
	    	-GG[6][ZDIR]-GG[8][ZDIR]-GG[11][ZDIR]-TT2[ZDIR]-GG[2][ZDIR];

	    localVAR[comp].F[DEST][ZDIR][H2] +=
	    	-GG[7][ZDIR]-GG[9][ZDIR]-GG[12][ZDIR]-TT2[ZDIR]-GG[3][ZDIR];	
	
	#else
		// LOCAL_COMP is not defined
		
		BEGIN_TRANSACTION();
		if(coun!=0)
		{	
			printf("aborted 7\n");fflush(stdout);
		}
		coun++;	
		//pthread_mutex_lock(&mol_lock[mol]);
			VAR[mol].F[DEST][XDIR][O] +=
		    	G110[XDIR] + GG[11][XDIR] +GG[12][XDIR]+C1*G23[XDIR];
		    VAR[mol].F[DEST][XDIR][H1] += 
		    	GG[6][XDIR]+GG[7][XDIR]+GG[13][XDIR]+TT[XDIR]+GG[4][XDIR];
		    VAR[mol].F[DEST][XDIR][H2] +=
		    	GG[8][XDIR]+GG[9][XDIR]+GG[14][XDIR]+TT[XDIR]+GG[5][XDIR];
		    VAR[mol].F[DEST][YDIR][O]  += 
		    	G110[YDIR]+GG[11][YDIR]+GG[12][YDIR]+C1*G23[YDIR];
		    VAR[mol].F[DEST][YDIR][H1] += 
		    	GG[6][YDIR]+GG[7][YDIR]+GG[13][YDIR]+TT[YDIR]+GG[4][YDIR];
		    VAR[mol].F[DEST][YDIR][H2] += 
		    	GG[8][YDIR]+GG[9][YDIR]+GG[14][YDIR]+TT[YDIR]+GG[5][YDIR];
		    VAR[mol].F[DEST][ZDIR][O]  +=
		    	G110[ZDIR]+GG[11][ZDIR]+GG[12][ZDIR]+C1*G23[ZDIR];
		    VAR[mol].F[DEST][ZDIR][H1] +=
		    	GG[6][ZDIR]+GG[7][ZDIR]+GG[13][ZDIR]+TT[ZDIR]+GG[4][ZDIR];
		    VAR[mol].F[DEST][ZDIR][H2] +=
		    	GG[8][ZDIR]+GG[9][ZDIR]+GG[14][ZDIR]+TT[ZDIR]+GG[5][ZDIR];
		//pthread_mutex_unlock(&mol_lock[mol]);
		
		//pthread_mutex_lock(&mol_lock[comp]);
			VAR[comp].F[DEST][XDIR][O] +=
		    	-G110[XDIR]-GG[13][XDIR]-GG[14][XDIR]-C1*G45[XDIR];
		    VAR[comp].F[DEST][XDIR][H1] +=
		    	-GG[6][XDIR]-GG[8][XDIR]-GG[11][XDIR]-TT2[XDIR]-GG[2][XDIR];
		    VAR[comp].F[DEST][XDIR][H2] +=
		    	-GG[7][XDIR]-GG[9][XDIR]-GG[12][XDIR]-TT2[XDIR]-GG[3][XDIR];
		    VAR[comp].F[DEST][YDIR][O] += 
		    	-G110[YDIR]-GG[13][YDIR]-GG[14][YDIR]-C1*G45[YDIR];
		    VAR[comp].F[DEST][YDIR][H1] += 
		    	-GG[6][YDIR]-GG[8][YDIR]-GG[11][YDIR]-TT2[YDIR]-GG[2][YDIR];
		    VAR[comp].F[DEST][YDIR][H2] += 
		    	-GG[7][YDIR]-GG[9][YDIR]-GG[12][YDIR]-TT2[YDIR]-GG[3][YDIR];
		    VAR[comp].F[DEST][ZDIR][O] +=
		    	-G110[ZDIR]-GG[13][ZDIR]-GG[14][ZDIR]-C1*G45[ZDIR];
		    VAR[comp].F[DEST][ZDIR][H1] +=
		    	-GG[6][ZDIR]-GG[8][ZDIR]-GG[11][ZDIR]-TT2[ZDIR]-GG[2][ZDIR];
		    VAR[comp].F[DEST][ZDIR][H2] +=
		    	-GG[7][ZDIR]-GG[9][ZDIR]-GG[12][ZDIR]-TT2[ZDIR]-GG[3][ZDIR];
		//pthread_mutex_unlock(&mol_lock[comp]);
		COMMIT_TRANSACTION();
		coun =0;
	#endif

}// end of subroutine UPDATE_FORCES

#endif //__WATER_INTERF_IMPL__
