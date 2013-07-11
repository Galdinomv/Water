
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "mdvar.h"
#include "global.h"
#include "mddata.h"
#include "water.h"
#include "data.h"
#include "cnst.h"
#include "fileio.h"
#include "frcnst.h"
#include "parameters.h"
#include "randno.h"
#include "split.h"
#include "wwpot.h"

struct GlobalMemory *gl;        /* pointer to the Global Memory
                                   structure, which contains the lock,
                                   barrier, and some scalar variables */

double BOX_LENGTH = 0.0;
int BPS_SQRD = 0;
double  TLC[100] = {0.0}, ELPST = 0.0, TKIN[MAX_PROCS]={0.0}, TVIR[MAX_PROCS]={0.0}, TTMV[MAX_PROCS]={0.0}, FPOT=0.0, FKIN=0.0;
int IX[3*MXOD2+1]={0}, IRST =0, NVAR=0, NXYZ=0, NXV=0, IXF=0, IYF=0, IZF=0, IMY=0, IMZ=0;
int NumProcs = 0;
int NumBoxes = 0;
struct GlobalSharedMemory *gsm;
int NSTEP =0;
int NSAVE =0;
int NRST  =0;
int NPRINT=0;
int NFMC  =0;
int NORD1 =0;
double UNITT = 0.0;
double UNITL = 0.0;
double UNITM = 0.0;
double BOLTZ = 0.0;
double AVGNO = 0.0;
double PCC[11] = {0.0};
FILE *one 	= NULL;	
FILE *five	= NULL;
FILE *six	= NULL;
FILE *nfmc	= NULL;
double FC11 = 0.0;
double FC12 = 0.0;
double FC13 = 0.0;
double FC33 = 0.0;
double FC111 = 0.0;
double FC333 = 0.0;
double FC112 = 0.0;
double FC113 = 0.0;
double FC123 = 0.0;
double FC133 = 0.0;
double FC1111 = 0.0;
double FC3333 = 0.0;
double FC1112 = 0.0;
double FC1122 = 0.0;
double FC1113 = 0.0;
double FC1123 = 0.0;
double FC1133 = 0.0;
double FC1233 = 0.0;
double FC1333 = 0.0;
double TEMP=0.0, RHO=0.0, TSTEP=0.0, BOXL=0.0, BOXH=0.0, CUTOFF=0.0, CUT2=0.0;
int    NMOL=0, NORDER=0, NATMO=0, NATMO3=0, NMOL1=0;
double R3[128]={0.0}, R1 =0.0;
int I2=0;
double OMAS=0.0, HMAS=0.0, WTMOL=0.0, ROH=0.0, ANGLE=0.0, FHM=0.0, FOM=0.0, ROHI=0.0,ROHI2=0.0;
int NATOMS=0;
double QQ=0.0, A1=0.0, B1=0.0, A2=0.0, B2=0.0, A3=0.0, B3=0.0;
double A4=0.0, B4=0.0, AB1=0.0, AB2=0.0, AB3=0.0, AB4=0.0;
double C1=0.0, C2=0.0, QQ2=0.0, QQ4=0.0;
double REF1=0.0, REF2=0.0, REF4=0.0;

extern char	*optarg;
double tempXTT = 0.0;

int II = 0;                         /*  variables explained in common.h */
int i = 0;
int NDATA =  0;
int NFRST=11;
int NFSV=10;
int LKT=0;


//first_last_array **start_end; /* ptr to array of start/end box #s *\/ */
//int NumProcs;                 /* number of processors being used;
//                                 run-time input           */

void * thread_work(void * param_in);

//***************************************************************************
int main(int argc, char *argv[])
{
    struct link *curr_ptr;
    unsigned threadID;
    FILE *fp;
    char *input_file = "sample.in";
    int c,tsteps = 0;
    ParamStruct param[MAX_PROCS] = {0};
    int rc, status;
    int i;
	    
    // process arguments
    while ((c = getopt(argc, argv, "i:t:hH?:")) != -1)
        switch (c) 
        {
            case 'i':	// "i" for input file
                input_file = optarg;
                break;
            case 't':	// "t" for t number of steps
                tsteps = atoi(optarg);
                break;
            case 'h':
            case 'H':
            case '?':
                exit(0);
                break;
            default:	// no argument is fine
                break;			  
        }

    // default values for the control parameters of the driver
    // are in parameters.h
    six = stderr;

    // input file for particle displacements
    if (!(fp = fopen(input_file,"r"))) {
        fprintf(stderr, "Unable to open '%s'\n", input_file); 
        fflush(stderr);
        exit(-1);
    }
			
    fscanf(fp, "%lf%d%d%d%d%d%d%d",
                &TSTEP, &NMOL, &NSTEP, &NORDER, 
                &NSAVE, &NRST, &NPRINT, &NFMC
            );

    if (tsteps) {
        NSTEP = tsteps;
    }
    if (NMOL > MAXMOLS) {
        fprintf(stderr, "Lock array in global.H has size %d < %d (NMOL)\n", MAXMOLS, NMOL);
        fflush(stderr);
        exit(-1);
    }

    NumProcs = NPROCS;
    printf("*******NumProcs=%d************\n", NumProcs);
    six = stdout;
    double temp, rho;
    temp = 298.0;
    rho = 0.9980;

    TEMP = temp;
    RHO  = rho;
            	
    // read input
    printf("*** HACK FOR NOW ***\n");
    printf("Initial NMOL %d\n",NMOL);
    TSTEP = ((double) 1.5e-16);
    NMOL = (int) 512;
    NSTEP = (int) 100;		
    NORDER = (int) 6;
    NSAVE = (int) -1;
    NRST = (int) 3000;
    NPRINT = (int) 1;
    NFMC = (int) 0;
    //NumProcs = (int) Tmk_nprocs;
    CUTOFF = ((double) 6.212752);
       
    //printf("Using %d\n",NumProcs); 
    //printf(" procs on %d steps of %d mols\n",NSTEP,NMOL); 
    //printf("Other parameters: \n"); 
    //printf("TSTEP = %19.18lf\n",TSTEP);
    //printf("NORDER = %d\n",NORDER);
    //printf("NSAVE = %d\n",NSAVE);
    //printf("NRST = %d\n",NRST);
    //printf("NPRINT = %d\n",NPRINT);
    //printf("NFMC = %d\n",NFMC);
    //printf("CUTOFF = %f\n",CUTOFF);

    // set up scaling factors and constants
    NORD1 = NORDER + 1;
    CNSTNT(NORD1,TLC);  // sub. call to set up constants
    SYSCNS();    // sub. call to initialize system constants

    printf("%d boxes with %d threads",BOX_PER_SIDE * BOX_PER_SIDE * BOX_PER_SIDE,NumProcs);
    
    if (NumProcs > (BOX_PER_SIDE * BOX_PER_SIDE * BOX_PER_SIDE)) {
        printf("ERROR: less boxes %d than processors\n ",BOX_PER_SIDE * BOX_PER_SIDE * BOX_PER_SIDE,NumProcs);
        exit(-1);
    }
    printf("\nTEMPERATURE                  = %f\n",TEMP);
    printf("DENSITY                      = %f\n",RHO);
    printf("NUMBER OF MOLECULES          = %d\n",NMOL);
    printf("NUMBER OF PROCESSORS         = %d\n",NumProcs);
    printf("TIME STEP                    = %19.18lf\n",TSTEP);
    printf("ORDER USED TO SOLVE F=MA     = %d\n",NORDER);
    printf("NO. OF TIME STEPS            = %d\n",NSTEP);
    printf("FREQUENCY OF DATA SAVING     = %d\n",NSAVE);
    printf("FREQUENCY TO WRITE RST FILE  = %d\n",NRST);  
    
    //----- Starting threads -----------------------------	
    for(i=0; i < NumProcs; i++){
        //printf("Starting thread %d\n\n\n\n",i);fflush(stdout);
        param[i].threadID = i;
        pthread_create(&pthreads[i],NULL, thread_work, (void *) (&param[i]));
    }
   
    //----- Processing -----------------------------------	
    gettimeofday(&start, NULL);
    // Waiting for all threads to join before continue
    for(i=0; i < NumProcs; i++){
        rc = pthread_join(pthreads[i], (void **)&status);
        if (rc){
            printf("ERROR; return code from pthread_join() of thread %d is %d\n", i, rc);
            exit(-1);
        }
    }
    gettimeofday(&end, NULL);
	
    //----- Joined threads ------------------------------
    fprintf(stderr, "\nExited Happily with XTT %g\n", tempXTT);
    fprintf(stderr, "Elapsed time: %.2f seconds\n",
            (((end.tv_sec * 1000000.0) +  end.tv_usec) -
            ((start.tv_sec * 1000000.0) + start.tv_usec)) / 1000000.0);
    fflush(stderr);    
    
    return 0;  
} // main()


// 0 = use the first context on all processors before using the second context on the first processor
// 1 = use both contexts on first processor before using the second processor
#define SMT_FIRST 1
#define MAX_THREAD_CONTEXTS 4

//***************************************************************************
void* thread_work(void* param_in)
{

    ParamStruct *theParam = (ParamStruct *) (param_in);
    int threadID = theParam->threadID;
    double XTT;	
    //cout << threadID << " of " << NumProcs << " just before goto" << endl;

    // bind threads to processors
    printf("&&&&&&&&&&&&&&&NUMBER OF PROCESSORS=%d\n", NumProcs);
    
    int mask;
    if (!SMT_FIRST){
        /* change this code so that it avoids binding to SMT contexts
        and instead uses MP contexts first */
        mask = 1 << ((threadID % NumProcs)*2 + threadID / NumProcs); 
    } else {
        mask = 1 << (threadID % MAX_THREAD_CONTEXTS); 
    }
    sched_setaffinity(0, sizeof(mask), (cpu_set_t*)&mask);
    
    //----- Start threads at the same time --------------
    pthread_barrier(0);

    //----- Main thread does init for shared mem --------
    if( threadID == 0) 
    {    
        //gl = (struct GlobalMemory *) malloc(sizeof(struct GlobalMemory));
        gl = new GlobalMemory;
        gsm = (struct GlobalSharedMemory *) malloc(sizeof(struct GlobalSharedMemory));
        gsm->global_ns_count = 0;
        gsm->global_boxes_count = 0;
    	
        gl->start = 0;		// used as a barrier id
        gl->InterfBar = 1;	// used as a barrier id
        gl->PotengBar = 2;	// used as a barrier id
        //gl->IntrafVirLock = 1;	// used as a pthread_mutex_lock index
        //gl->InterfVirLock = 2;	// used as a pthread_mutex_lock index
       	 	    
        int kk = 0;
        for (kk = 0; kk < NMOL; kk++){
            gl->MolLock[kk] = kk + 8;
        }
        
        //gl->KinetiSumLock = 6;
        //gl->PotengSumLock = 7;  
    }

    if(threadID == 0)
    { // do memory initializations
    
        int pid, procnum, i, j, k, l;
        struct list_of_boxes *new_box;
        struct list_of_boxes *temp_box;
        int xprocs, yprocs, zprocs;
        int x_inc, y_inc, z_inc;
        int x_ct, y_ct, z_ct;
        int x_left, y_left, z_left;
        int x_first, y_first, z_first;
        int x_last, y_last, z_last;
        double proccbrt;
        unsigned gmem_size = sizeof(struct GlobalMemory);
        
        {;};  // macro call to initialize shared memory etc.
        
        
        // Calculate start and finish box numbers for processors
        xprocs = 0;
        yprocs = 0;
        proccbrt = (double) pow((double) NumProcs, 1.0/3.0) + 0.00000000000001;
        j = (int) proccbrt;
        if (j<1) j = 1;
        while ((xprocs == 0) && (j>0)) 
        {
            k = (int) sqrt((double) (NumProcs / j));
            if (k<1) k=1;
            while ((yprocs == 0) && (k>0)) 
            {
                l = NumProcs/(j*k);
                if ((j*k*l) == NumProcs) 
                {
                    xprocs = j;
                    yprocs = k;
                    zprocs = l;
                } // if
                k--;
            } // while yprocs && k
            j--;
        } // while xprocs && j

        //printf("xprocs = %d yprocs = %d  zprocs = %d\n",xprocs,yprocs,zprocs);
        //printf("BOX_PER_SIDE %d\n",BOX_PER_SIDE);
        // Fill in start_end array values        

        procnum = 0;
        x_inc = BOX_PER_SIDE/xprocs;
        y_inc = BOX_PER_SIDE/yprocs;
        z_inc = BOX_PER_SIDE/zprocs;
        
        x_left = BOX_PER_SIDE - (xprocs*x_inc);
        y_left = BOX_PER_SIDE - (yprocs*y_inc);
        z_left = BOX_PER_SIDE - (zprocs*z_inc);
        //printf("x_inc = %d y_inc = %d z_inc = %d \n ",x_inc,y_inc,z_inc);
        //printf("x_left = %d y_left = %d z_left = %d \n ",x_left,y_left,z_left);        
        fflush(stdout);

        x_first = 0;
        x_ct = x_left;
        x_last = -1;
        x_inc++;
        for (i=0; i<xprocs; i++) 
        {
            y_ct = y_left;
            if (x_ct == 0) x_inc--;
            x_last += x_inc;
            y_first = 0;
            y_last = -1;
            y_inc++;
            for (j=0; j<yprocs; j++) 
            {
                z_ct = z_left;
                if (y_ct == 0) y_inc--;
                y_last += y_inc;
                z_first = 0;
                z_last = -1;
                z_inc++;
                for (k=0; k<zprocs; k++) 
                {
                    if (z_ct == 0) z_inc--;
                    z_last += z_inc;
                    start_end[procnum]->box[XDIR][FIRST] = x_first;
                    start_end[procnum]->box[XDIR][LAST] = 
                        min((int) x_last, (int) BOX_PER_SIDE - 1);
                    start_end[procnum]->box[YDIR][FIRST] = y_first;
                    start_end[procnum]->box[YDIR][LAST] = 
                        min((int) y_last, (int) BOX_PER_SIDE - 1);
                    start_end[procnum]->box[ZDIR][FIRST] = z_first;
                    start_end[procnum]->box[ZDIR][LAST] = 
                        min((int) z_last, (int) BOX_PER_SIDE - 1);
                    z_first = z_last + 1;
                    z_ct--;
                    procnum++;
                }
                y_first = y_last + 1;
                y_ct--;
            }
            x_first = x_last + 1;
            x_ct--;
        }
        
        // Set all box ptrs to null
        for (i=0; i<NumProcs; i++) my_boxes[i] = NULL;
        
        // Set up links for all boxes for initial interf and intraf
        temp_box = my_boxes[0];
        while (temp_box) {
            temp_box = temp_box->next_box;
        }
        
        // Allocate space for BOX array
        int bps = BOX_PER_SIDE;
        
        for (i=0; i < BOX_PER_SIDE; i++) 
        {  
            for (j=0; j < BOX_PER_SIDE; j++) 
            {        
                for (k=0; k < BOX_PER_SIDE; k++) 
                {
                    BOX[i][j][k].list = NULL;
                    {;};
                }
            }
        } // for i
	
        // macro calls to initialize synch variables
        {;};
        {;};
        {;};
        {;};
        {;};
        {;};
        {;};
        {;};
        {;};
    }


    //pthread_barrier(0);
    //printf("threadID = %d\n\n\n\n\n\n", threadID);fflush(stdout);
    if(threadID != 0) goto skipinit;

    // cout << "SPHERICAL CUTOFF RADIUS    = " << CUTOFF << "  ANGSTROM" << endl;
    
    IRST=0;
    
    // call initialization routine
    INITIA();
     
    struct list_of_boxes *new_box, *curr_box;
    threadID = 0;


    for(threadID=0; threadID < NumProcs; threadID++) 
    {
        int i=0, j=0, k=0;
        
        for (i=start_end[threadID]->box[XDIR][FIRST]; i<=start_end[threadID]->box[XDIR][LAST]; i++) 
        {
            for (j=start_end[threadID]->box[YDIR][FIRST]; j<=start_end[threadID]->box[YDIR][LAST]; j++) 
            {
                for (k=start_end[threadID]->box[ZDIR][FIRST]; k<=start_end[threadID]->box[ZDIR][LAST]; k++) 
                {    
                    new_box = &gsm->GLOBAL_BOXES[gsm->global_boxes_count];
                    gsm->global_boxes_count++;
                    new_box->coord[XDIR] = i;
                    new_box->coord[YDIR] = j;
                    new_box->coord[ZDIR] = k;
                    new_box->next_box = NULL;
                    curr_box = my_boxes[threadID];
                    if (curr_box == NULL)
                        my_boxes[threadID] = new_box;
                    else {
                        while (curr_box->next_box != NULL)
                            curr_box = curr_box->next_box;
                        curr_box->next_box = new_box;
                    } // else
                }
            }
        }
    }
    threadID = 0;
    gl->tracktime = 0;
    gl->intratime = 0;
    gl->intertime = 0;

    // initialize Index to 1 so that the first created child gets id 1, not 0    
    gl->Index = 1;
    
    // spawn helper processes   
    {;};
    printf("CREATETIME = %u\n", gl->createend-gl->createstart);
    {;};
    
    if (NSAVE > 0) {  // not true for input decks provided
        {;};
        fprintf(six,"COLLECTING X AND V DATA AT EVERY %4d TIME STEPS \n",NSAVE);
        {;};
    }
    
    // CALL ROUTINE TO DO THE TIMESTEPS, with own id passed as 0    
    {/*long time();*/ (gl->computestart) = time(0);};


    //----- Shared region executed by all threads -----------------     
skipinit:
    //----- Barrier 1 ------------
    pthread_barrier(0);
    //cout << threadID << " thinks gl = " << gl << " and gsm = " << gsm << endl;
    //cerr << threadID << " thinks gsm->start_end is " << start_end << " and gsm->BOX is " << BOX << endl;
    //cout << threadID << " went past the first barrier" << endl;    
    if(threadID == 0) 
    {    
        //cout<< threadID << " has called MDMAIN" << endl;
        //printf("bjshgfskgfklsdgf\n\n\n\n");
        XTT = MDMAIN(NFSV,NFRST,NSTEP,NRST,NPRINT,NSAVE,LKT,NORD1,threadID); 
    } else {
        //printf("here kmn,m.b,\n\n\n\n");fflush(stdout);
        //cout << threadID << " has called MDMAIN" << endl;
        XTT = MDMAIN(NFSV,NFRST,NSTEP,NRST,NPRINT,NSAVE,LKT,NORD1,threadID); 
    }
    //cout << threadID << " is just before the second barrier" << endl;

    //----- Barrier 2 ------------
    pthread_barrier(1);
    if(threadID == 0) 
    {
        // macro to make main process wait for all others to finish
        {;};
        {/*long time();*/ (gl->computeend) = time(0);};
    
        printf("COMPUTESTART (after initialization) = %u\n",gl->computestart);
        printf("COMPUTEEND = %u\n",gl->computeend);
        printf("COMPUTETIME (after initialization) = %u\n",gl->computeend-gl->computestart);
        printf("Measured Time (2nd timestep onward) = %u\n",gl->tracktime);
        printf("Intramolecular time only (2nd timestep onward) = %u\n",gl->intratime);
        printf("Intermolecular time only (2nd timestep onward) = %u\n",gl->intertime);
        printf("Other time (2nd timestep onward) = %u\n",gl->tracktime - gl->intratime - gl->intertime);
        //printf("Exited Happily with XTT = %f\n",XTT);
    }
     
    //----- Barrier 3 ------------
    pthread_barrier(0);
    //cout << threadID << " is exiting" << endl;
}
