#include <stdio.h>

#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "split.h"
#include "global.h"
#include "pthread.h"
# include "data.h"

/************************************************************************/

/* double  MDMAIN(NFSV,NFRST,NSTEP,NRST,NPRINT,NSAVE,LKT,NORD1,threadID) */
/*   int NFSV,NFRST,NSTEP,NRST,NPRINT,NSAVE,LKT,NORD1; */
/*   unsigned threadID; */
double MDMAIN(int NFSV, int NFRST, int NSTEP, int NRST, int NPRINT, int NSAVE, int LKT, int NORD1, unsigned threadID) 
{


	
    
    double TVIR = 0.0;
    double TTMV = 0.0;
    double TKIN = 0.0;
    double XTT;
    int i,j,k;
    double POTA,POTR,POTRF;
    double XVIR,AVGT,TEN;
    struct link *curr_ptr;
    struct list_of_boxes *new_box, *curr_box;

    //cout << threadID << " was here " << endl;
    //cout << threadID << " thinks start_end is " << start_end << endl;
    //cout << threadID << " was here " << endl;

    /* calculate initial value for acceleration */
    INTRAF(threadID);

    //cout << threadID << " finished intraf " << endl;

    //cout << threadID << " has arrived at barrier 10" << endl;
	pthread_barrier(10);
    //cout << threadID << " went past barrier 10" << endl;

    INTERF(ACC,threadID);

    //cout << threadID << " has finished interf " << endl;
    pthread_barrier(11);
 
    /* MOLECULAR DYNAMICS LOOP */
    
    for (i=1;i <= NSTEP; i++) {
      
        TTMV=TTMV+1.00;
        
        /* POSSIBLE ENHANCEMENT:  Here's where one start measurements to avoid 
           cold-start effects.  Recommended to do this at the beginning of the
           second timestep; i.e. if (i == 2).
           */
        
        /* initialize various shared sums */
        if (threadID == 0) {
            int dir;
            if (i >= 2) {
	      {/*long time();*/ (gl->trackstart) = time(0);};
            }                
            gl->VIR = 0.0;
            gl->POTA = 0.0;
            gl->POTR = 0.0;
            gl->POTRF = 0.0;
            for (dir = XDIR; dir <= ZDIR; dir++)
                gl->SUM[dir] = 0.0;
        }
        
        if ((threadID == 0) && (i >= 2)) {
	  {/*long time();*/ (gl->intrastart) = time(0);};
        }
        
        pthread_barrier(10);

        PREDIC(TLC,NORD1,threadID);
        INTRAF(threadID);
        
	pthread_barrier(10);
        
        if ((threadID == 0) && (i >= 2)) {
	  {/*long time();*/ (gl->intraend) = time(0);};
            gl->intratime = gl->intratime + gl->intraend - gl->intrastart;
        }
        
        if ((threadID == 0) && (i >= 2)) {
	  {/*long time();*/ (gl->interstart) = time(0);};
        }
        
        INTERF(FORCES,threadID); 
        
        if ((threadID == 0) && (i >= 2)) {
	  {/*long time();*/ (gl->interend) = time(0);};
            gl->intertime = gl->intertime + gl->interend - gl->interstart;
        }
        
        CORREC(NORD1,threadID);
        
        BNDRY(threadID);
//printf("gl->SUM[0] = %f\n", gl->SUM[0]);
//printf("gl->SUM[1] = %f\n", gl->SUM[1]);
//printf("gl->SUM[2] = %f\n", gl->SUM[2]);        

        KINETI(gl->SUM,HMAS,OMAS,threadID);
//printf("gl->SUM[0] = %f\n", gl->SUM[0]);
//printf("gl->SUM[1] = %f\n", gl->SUM[1]);
//printf("gl->SUM[2] = %f\n", gl->SUM[2]);        
      
        pthread_barrier(10);
        
        if ((threadID == 0) && (i >= 2)) {
	  {/* long time();*/ (gl->intraend) = time(0);};
            gl->intratime = gl->intratime + gl->intraend - gl->interend;
        }
        
        TKIN=TKIN+gl->SUM[0]+gl->SUM[1]+gl->SUM[2];
//printf("TKIN = %f\n", TKIN);

        TVIR=TVIR-gl->VIR;
        
        /* CHECK if  PRINTING AND/OR SAVING IS TO BE DONE */
        
        if ( ((i % NPRINT) == 0) || ((NSAVE > 0) && ((i % NSAVE) == 0))) {
            
            /* if so, call poteng to compute potential energy.  Note
               that we are attributing all the time in poteng to intermolecular
               computation although some of it is intramolecular (see poteng.C) */
            
            if ((threadID == 0) && (i >= 2)) {
	      {/*long time();*/ (gl->interstart) = time(0);};
            }
            
            POTENG(&gl->POTA,&gl->POTR,&gl->POTRF,threadID);
            
            pthread_barrier(10);
            
            if ((threadID == 0) && (i >= 2)) {
	      {/*long time();*/ (gl->interend) = time(0);};
                gl->intertime = gl->intertime + gl->interend - gl->interstart;
            }
            
            POTA=gl->POTA*FPOT;
            POTR=gl->POTR*FPOT;
            POTRF=gl->POTRF*FPOT;
            XVIR=TVIR*FPOT*0.50/TTMV;
            AVGT=TKIN*FKIN*TEMP*2.00/(3.00*TTMV);
            TEN=(gl->SUM[0]+gl->SUM[1]+gl->SUM[2])*FKIN;
            XTT=POTA+POTR+POTRF+TEN;
            
            /* if it is time to print output as well ... */      
            if ((i % NPRINT) == 0 && threadID == 0) {

                printf("     %5d %14.5lf %12.5lf %12.5lf %12.5lf \n"
                        ,i,TEN,POTA,POTR,POTRF);
                printf(" %16.3lf %16.5lf %16.5lf\n",XTT,AVGT,XVIR);
                fflush(stdout);

            }
            
        }
        
        pthread_barrier(10);
        
        if ((threadID == 0) && (i >= 2)) {
	  {/*long time();*/ (gl->trackend) = time(0);};
            gl->tracktime = gl->tracktime + gl->trackend - gl->trackstart;
        }
        
    } /* for i */
if(threadID==0){
    tempXTT = XTT;
}
    return(XTT);
    
} /* mdmain.c */


