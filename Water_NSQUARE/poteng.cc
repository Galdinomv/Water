#ifndef __WATER_POTENG_IMPL__
#define __WATER_POTENG_IMPL__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "data.h"

#include "mdvar.h"
#include "frcnst.h"
#include "water.h"
#include "wwpot.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

//
// some shared sums computed by POTENG
//
/*
    this routine calculates the potential energy of the system.
     FC11 ,FC12, FC13, and FC33 are the quardratic force constants
*/
//	Note: on PThread parallelization, for update global shared memory VAR
//	there is a single callsite to POTENG inside redsync.c
//	since there is no protection on this call, var_lock Critical Section is added
//
//	it is now safe.
//
void POTENG(tm_double * POTA, tm_double * POTR, tm_double * PTRF, int threadID)
{
    int    mol=0, comp=0;
    int    half_mol=0;
    int    KC=0, K=0;
    double R1=0.0, R2=0.0, RX=0.0, COS=0.0, DT=0.0, DR1=0.0, DR2=0.0, DR1S=0.0, DR2S=0.0, DRP=0.0;
    double XL[15]={0.0}, YL[15]={0.0}, ZL[15]={0.0}, RS[15]={0.0}, RL[15]={0.0};
    double DTS=0.0;
    double LPOTA=0.0, LPOTR=0.0, LPTRF=0.0;
    double tempa=0.0, tempb=0.0, tempc=0.0;

    //  compute intra-molecular potential energy
    LPOTA=0.0;
    
    for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
    {
        int coun = 0;
        INIT_TRANSACTIONS();
        BEGIN_TRANSACTION();
        if(coun!=0)
        {	
        	printf("aborted 4\n");fflush(stdout);
        }
        coun++;
        //pthread_mutex_lock(&var_lock);
        VAR[mol].VM[XDIR] =   C1 * VAR[mol].F[DISP][XDIR][ O] 
        					+ C2 * (VAR[mol].F[DISP][XDIR][H1] 
        					+ VAR[mol].F[DISP][XDIR][H2] );
        
        VAR[mol].VM[YDIR] = C1*VAR[mol].F[DISP][YDIR][O] 
        					+ C2*(VAR[mol].F[DISP][YDIR][H1] 
        					+ VAR[mol].F[DISP][YDIR][H2] );
        					
        VAR[mol].VM[ZDIR] = C1*VAR[mol].F[DISP][ZDIR][O] +
        		    C2*(VAR[mol].F[DISP][ZDIR][H1] + 
        			VAR[mol].F[DISP][ZDIR][H2] );
        //pthread_mutex_unlock(&var_lock);
        COMMIT_TRANSACTION();
        coun =0;
        tempa = VAR[mol].F[DISP][XDIR][O]-VAR[mol].F[DISP][XDIR][H1];
        tempb = VAR[mol].F[DISP][YDIR][O]-VAR[mol].F[DISP][YDIR][H1];
        tempc = VAR[mol].F[DISP][ZDIR][O]-VAR[mol].F[DISP][ZDIR][H1]; 
            R1 = tempa * tempa + tempb * tempb + tempc * tempc;

        tempa = VAR[mol].F[DISP][XDIR][O]-VAR[mol].F[DISP][XDIR][H2];
        tempb = VAR[mol].F[DISP][YDIR][O]-VAR[mol].F[DISP][YDIR][H2];
        tempc = VAR[mol].F[DISP][ZDIR][O]-VAR[mol].F[DISP][ZDIR][H2]; 
            R2 = tempa * tempa + tempb * tempb + tempc * tempc;

        RX = ( (VAR[mol].F[DISP][XDIR][O]-VAR[mol].F[DISP][XDIR][H1]) *
               (VAR[mol].F[DISP][XDIR][O]-VAR[mol].F[DISP][XDIR][H2])  ) +
             ( (VAR[mol].F[DISP][YDIR][O]-VAR[mol].F[DISP][YDIR][H1]) *
               (VAR[mol].F[DISP][YDIR][O]-VAR[mol].F[DISP][YDIR][H2])  ) +
             ( (VAR[mol].F[DISP][ZDIR][O]-VAR[mol].F[DISP][ZDIR][H1]) *
               (VAR[mol].F[DISP][ZDIR][O]-VAR[mol].F[DISP][ZDIR][H2])  );

        R1 	= 	sqrt(R1);
        R2 	= 	sqrt(R2);
        COS	= 	RX/(R1*R2);
        DT 	= 	(acos(COS)-ANGLE)*ROH;
        DR1	= 	R1-ROH;
        DR2	= 	R2-ROH;
        DR1S=	DR1*DR1;
        DR2S=	DR2*DR2;
        DRP	=	DR1+DR2;
        DTS	=	DT*DT;
        
        LPOTA += (FC11*(DR1S+DR2S)+FC33*DTS)*0.5
           		+ FC12*DR1*DR2+FC13*DRP*DT
               	+ (FC111*(DR1S*DR1+DR2S*DR2)+FC333*DTS*DT+FC112*DRP*DR1*DR2
               	+ FC113*(DR1S+DR2S)*DT+FC123*DR1*DR2*DT+FC133*DRP*DTS)*ROHI;
        
        LPOTA += (FC1111*(DR1S*DR1S+DR2S*DR2S)+FC3333*DTS*DTS
        		+ FC1112*(DR1S+DR2S)*DR1*DR2+FC1122*DR1S*DR2S
        		+ FC1113*(DR1S*DR1+DR2S*DR2)*DT+FC1123*DRP*DR1*DR2*DT
        		+ FC1133*(DR1S+DR2S)*DTS+FC1233*DR1*DR2*DTS
        		+ FC1333*DRP*DTS*DT)*ROHI2;
    }// for mol



    pthread_barrier(gl->PotengBar);	// this barrier has to be removed, due to the calling feature from within MSMAIN()
	
    //compute inter-molecular potential energy
    LPOTR=0.0;
    LPTRF=0.0;
    half_mol = NMOL/2;

    for (mol = StartMol[threadID]; mol < StartMol[threadID+1]; mol++) 
    {
        int comp_last = mol + half_mol;
        int icomp;
        if (NMOL%2==0 && ((!(mol%2) && (mol < half_mol)) || ((mol%2) && mol > half_mol))) 
        {
        	comp_last -= 1;
        }
        
        for (icomp = mol+1; icomp <= comp_last; icomp++) 
        {
            comp = icomp;
            if (comp > NMOL1) 
            {
            comp = comp%NMOL;
            }
			
            CSHIFT(VAR[mol].F[DISP][XDIR],VAR[comp].F[DISP][XDIR],
               VAR[mol].VM[XDIR], VAR[comp].VM[XDIR],XL,BOXH,BOXL);
            
            CSHIFT(VAR[mol].F[DISP][YDIR],VAR[comp].F[DISP][YDIR],
               VAR[mol].VM[YDIR], VAR[comp].VM[YDIR],YL,BOXH,BOXL);
            
            CSHIFT(VAR[mol].F[DISP][ZDIR],VAR[comp].F[DISP][ZDIR],
               VAR[mol].VM[ZDIR], VAR[comp].VM[ZDIR],ZL,BOXH,BOXL);
            
            KC=0;
            for (K = 0; K < 9; K++) 
            {
                RS[K] = XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K];
                if (RS[K] > CUT2)
                {
                    KC = KC+1;
                }
            } // for k

            if (KC != 9) 
            {
                for (K = 0; K < 9; K++) 
                {
                    if (RS[K] <= CUT2)
                    {
                        RL[K]=sqrt(RS[K]);
                    }
                    else 
                    {
                        RL[K]=CUTOFF;
                        RS[K]=CUT2;
                    } // else
                } // for K

                LPOTR= LPOTR-QQ2/RL[1]-QQ2/RL[2]-QQ2/RL[3]-QQ2/RL[4]
                        +QQ /RL[5]+QQ /RL[6]+QQ /RL[7]+QQ /RL[8]
                        +QQ4/RL[0];
                
                LPTRF= LPTRF-REF2*RS[0]-REF1*((RS[5]+RS[6]+RS[7]+RS[8])*0.5
                                -RS[1]-RS[2]-RS[3]-RS[4]);

                if (KC <= 0) 
                {
                    for (K = 9; K <  14; K++)
                    {
                        RL[K]=sqrt(XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K]);
                    }
					
                    LPOTR = LPOTR+A1* exp(-B1*RL[9])
                                + A2*(exp(-B2*RL[ 5])+exp(-B2*RL[6])
                                + exp(-B2*RL[ 7])+exp(-B2*RL[ 8]))
                                + A3*(exp(-B3*RL[10])+exp(-B3*RL[11])
                                + exp(-B3*RL[12])+exp(-B3*RL[13]))
                                - A4*(exp(-B4*RL[10])+exp(-B4*RL[11])
                                + exp(-B4*RL[12])+exp(-B4*RL[13]));
                } // if KC <= 0
            } // if KC != 9
	} // for icomp
    } // for mol

    // update shared sums from computed  private sums
    int coun = 0;	
    INIT_TRANSACTIONS();
    BEGIN_TRANSACTION();
    if(coun!=0)
    {	
        printf("aborted\n");fflush(stdout);
    }
    coun++;
    //pthread_mutex_lock(&PotengSumLock);
        *POTA = *POTA + LPOTA;
        *POTR = *POTR + LPOTR;
        *PTRF = *PTRF + LPTRF;
    //pthread_mutex_unlock(&PotengSumLock);
    COMMIT_TRANSACTION();
    coun = 0;
}// end of subroutine POTENG

#endif //__WATER_POTENG_IMPL__
