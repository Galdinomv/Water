#ifndef __WATER_WATER_IMPL__
#define __WATER_WATER_IMPL__

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "data.h"
#include "comments.h"	// moves leading comments to another file
#include "split.h"
#include "global.h"

//include files for declarations
//#define extern
#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "frcnst.h"
#include "randno.h"
#include "global.h"
#include "include/src/tm_basics.h" 
//#undef extern

const int NFRST	=	11;
const int NFSV	=	10;
const int LKT	=	0;

// pointer to the Global Memory structure, which contains the lock, barrier, 
// and some scalar variables
struct GlobalMemory *gl = 	NULL;
molecule_type *VAR= NULL;
// this section explicitly collects all the externally declared data                                
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

// scalar -> 1D Array privatization also happens here
double TLC[100] = {0.0}, ELPST = 0.0, TKIN[MAX_PROCS]={0.0}, TVIR[MAX_PROCS]={0.0}, TTMV[MAX_PROCS]={0.0}, FPOT=0.0, FKIN=0.0;
int IX[3*MXOD2+1]={0}, IRST =0, NVAR=0, NXYZ=0, NXV=0, IXF=0, IYF=0, IZF=0, IMY=0, IMZ=0;
// number of the first molecule to be handled by this process; used for static scheduling
int StartMol[MAX_PROCS+1]={0};
int StartPos[MAX_PROCS+1]={0};
int MolsPerProc=0;	// number of mols per processor


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

// XTT is promoted from main() local to global level, with a reduction on it
double XTT = 0.0;
double tempXTT[MAX_PROCS] = {0.0};	// a reduction on XTT is needed

extern char	*optarg;

void usage(void);	// display usage(help) options
int main(int, char **);	// main function


// display usage(help) options
void usage(void)
{
	fprintf(stderr, "water:\n\t%s \n\t%s \n\t%s \n","-i input_file","-t steps","-h: usage");	
}

void InitData(void)
{	
}

//////////////////////////////////////////////////////////
//
//	The parallel entry of work that each thread
//	will take after being created successfully
//
///////////////////////////////////////////////////////////

// 0 = use the first context on all processors before using the second context on the first processor
// 1 = use both contexts on first processor before using the second processor
#define SMT_FIRST 1
#define MAX_THREAD_CONTEXTS 4

void * thread_work(void * param_in)
{
    ParamStruct *theParam = (ParamStruct *) (param_in);
    int threadID = theParam->threadID;
    tm_double VIR = 0.0;

    // bind threads to processors
    printf("&&&&&&&&&&&&&&&NUMBER OF PROCESSORS=%d\n", NPROCS);

    int mask;
    if (!SMT_FIRST){
        /* change this code so that it avoids binding to SMT contexts
        and instead uses MP contexts first */
        mask = 1 << ((threadID % NPROCS)*2 + threadID / NPROCS); 
    } else {
        mask = 1 << (threadID % MAX_THREAD_CONTEXTS); 
    }
    sched_setaffinity(0, sizeof(mask), (cpu_set_t*)&mask);
    
    // make sure all threads have arrived here, before continue	
    pthread_barrier(0);
        
    if (threadID == 0) // Do memory initializations
    { 
        int kk;
        unsigned mol_size 	= sizeof(molecule_type) * NMOL;
        unsigned gmem_size 	= sizeof(struct GlobalMemory);
        
        // allocate space for main (VAR) data structure as well as
        // synchronization variables
        //VAR = (molecule_type *) malloc(mol_size);
        VAR = new molecule_type[NMOL];	
        if(VAR == NULL)
        {
            printf("Unable to allocate meory for VAR\n");
            exit(-1);	
        }
        //memset(VAR, '\0', mol_size);
        int g,e,r,t;
        for(g=0; g<NMOL; g++)
        {
            VAR[g].VM[1] = 0.0;
            VAR[g].VM[2] = 0.0;
            VAR[g].VM[3] = 0.0;
            for(e=0;e<MXOD2;e++)
            {
                for(r=0;r<NDIR;r++)
                {
                    for(t=0;t<NATOM;t++)
                    {
                        VAR[g].F[e][r][t] = 0.0;
                    }
                }
            }
        }		

        gl = new GlobalMemory;
        //gl = (struct GlobalMemory *) malloc(gmem_size);
        if(gl == NULL){
            printf("Unable to allocate memory for gl\n");
            exit(-1);	
        }
        //memset(gl, '\0', gmem_size);
		
        // macro calls to initialize synch varibles, and other fields inside GlobalMemory
        gl->start = 0;		// used as a barrier id
        gl->InterfBar = 1;	// used as a barrier id
        gl->PotengBar = 2;	// used as a barrier id
        
        // TRANSACTIONAL => ELIM LOCKS
        //gl->IntrafVirLock = 1;	// used as a pthread_mutex_lock index
        //gl->InterfVirLock = 2;	// used as a pthread_mutex_lock index        
        //for (kk = 0; kk < NMOL; kk++){
        //	gl->MolLock[kk] = kk + 8;
        //}
        //gl->KinetiSumLock = 6;
        //gl->PotengSumLock = 7;

        gl->VIR = 0.0;
        gl->SUM[1] = 0.0 ;
        gl->POTA = 0.0;
        gl->POTR = 0.0;
        gl->POTRF = 0.0;
    }// end if (threadID==0)

    // set up control for static scheduling
    // some more initialization for StartMol array, done by thread 0
    if(threadID == 0)
    {
        int pid;
        MolsPerProc = NMOL / NPROCS;
        StartMol[0] = 0;
	    
        for (pid = 1; pid < NPROCS; pid += 1) 
        {
            // StartMol[pid] = StartMol[pid-1] + MolsPerProc;
            StartMol[pid] = (pid * NMOL) / NPROCS;
        }
        StartMol[NPROCS] = NMOL;
    }
	
    if(threadID == 0)
    {
        // sub. call to initialize system constants
        SYSCNS();
        
        pthread_mutex_lock(&print_lock);
            fprintf(six,	"SPHERICAL CUTOFF RADIUS    = %8.4f ANGSTROM\n",CUTOFF);
            fflush(six);
        pthread_mutex_unlock(&print_lock);
		
        IRST=0;
	
        // if there is no input displacement file, and we are to
        // initialize to a regular lattice
        if (NFMC == 0) 
        {
            fclose(nfmc);
            nfmc = NULL;
        }

        // initialization routine; also reads displacements and
        //	sets up random velocities
        INITIA(nfmc);
        
        // .......ESTIMATE ACCELERATION FROM F/M
        // note that these initial calls to the force-computing 
        // routines  (INTRAF and INTERF) use only 1 process since 
        // others haven't been created yet 
        {
            int tmp_smol 	= StartMol[1];
            StartMol[0] 	= 0;
            StartMol[1] 	= NMOL;
            INTRAF(&VIR, threadID);
            INTERF(ACC, &VIR, threadID);

            //----- some hardcoding moved from INTERF() -----

            int mal,diro,coun = 0;;
            for (mal = StartMol[threadID]; mal < StartMol[threadID+1]; mal++) 
            {
                for ( diro = XDIR; diro  <= ZDIR; diro++) 
                {
                
                    INIT_TRANSACTIONS();
                    BEGIN_TRANSACTION();
                    if(coun!=0)
                    {	
                    printf("aborted 9\n");fflush(stdout);
                    }
                    coun++;	
                    //pthread_mutex_lock(&var_lock);
                    
                    VAR[mal].F[ACC][diro][H1] = VAR[mal].F[ACC][diro][H1] * FHM;
                    VAR[mal].F[ACC][diro][H2] = VAR[mal].F[ACC][diro][H2] * FHM;
                    VAR[mal].F[ACC][diro][O]  = VAR[mal].F[ACC][diro][O]  * FOM;
                    
                    //pthread_mutex_unlock(&var_lock);
                    COMMIT_TRANSACTION();
                    coun =0;
                } // for diro
            } // for mal

            //----- some hardcoding moved from INTERF() -----



            StartMol[1] 	= tmp_smol;
        }// estimate of acceleration
            
        NFMC = -1;
        //.....START MOLECULAR DYNAMIC LOOP
        if (NFMC < 0) 
        {
            ELPST	=	0.00;
            //TKIN	=	0.00;
            memset(&TKIN, 0, sizeof(TKIN));
            //TVIR	=	0.00;
            memset(&TVIR, 0, sizeof(TVIR));
            //TTMV    =   0.00;
            memset(&TTMV, 0, sizeof(TTMV));
        }
		
        if(NSAVE > 0){ // not true for input decks provided
            fprintf(six,"COLLECTING X AND V DATA AT EVERY %4d TIME STEPS \n",NSAVE);
        }

    } // end if (threadID==0)
    // end of data initialization with threadID == 0
		
    pthread_barrier(0);
    
    // this is the last step toward a successful PThread parallelized Water code
    MDMAIN(NFSV, NFRST, NSTEP, NRST, NPRINT, NSAVE, LKT, NORD1, threadID);
    // Note:
    // value(s) have been stored into tempXXT[] for a later on reduction
    //
	
    pthread_barrier(0);	
    XTT = tempXTT[0];
    pthread_barrier(0);	
	
    //return NULL;
}

int main(int argc, char **argv)
{
    FILE *fp;
    char *input_file = "sample.in";
    int mol, pid, func, c, dir, atom, tsteps = 0;
    //double XTT = 0.0;
    tm_double VIR = 0.0;
    ParamStruct param[MAX_PROCS] = {0};
    int rc, status;
    //struct timeval start, finish;
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
            usage();
            exit(0);
            break;
        default:	// no argument is fine
            break;			  
    }

    // default values for the control parameters of the driver
    // are in parameters.h
    six = stderr;
    
    // input file for particle displacements
    nfmc = fopen("LWI12","r"); 
    if(nfmc == NULL)
    {
        printf("Unable to open LWI12 file\n");
        exit(-1);	
    }
	
    TEMP  	=	298.0;
    RHO   	=	0.9980;
    CUTOFF	=	0.0;
    
    if (!(fp = fopen(input_file,"r"))) 
    {
        fprintf(stderr, "Unable to open '%s'\n", input_file);
        fflush(stderr);
        exit(-1);
    }

    // READ INPUT
    /*
    *   TSTEP = time interval between steps
    *   NMOL  = number of molecules to be simulated
    *   NSTEP = number of time-steps to be simulated
    *   NORDER = order of the predictor-corrector method. 6 by default
    *   NSAVE = frequency of data saving.  -1 by default
    *   NRST  = no longer used
    *   NPRINT = frequency (in time-steps) of computing pot. energy and
    *            writing output file.  setting it larger than NSTEP
    *            means not doing any I/O until the very end.
    *   NFMC = file number to read initial displacements from.  set to 
    *          0 if program should generate a regular lattice initially.
    *
    *   NPROCS = number of processors to be used.
    */

    // this allows the control on the number of CPUS (NPROSC)
    fscanf(fp, "%lf%d%d%d%d%d%d%d",&TSTEP, &NMOL, &NSTEP, &NORDER, &NSAVE, &NRST, &NPRINT, &NFMC);
	
    if (tsteps) 
    {
        NSTEP = tsteps;
    }

    if (NMOL > MAXMOLS) 
    {
        fprintf(stderr, "Lock array in global.H has size %d < %d (NMOL)\n", MAXMOLS, NMOL);
        fflush(stderr);
        exit(-1);
    }

    if (NPROCS > MAX_PROCS) 
    {
        fprintf(stderr, "**ERROR**\nCannot process %d threads, maximum of %d threads allowed\n", NPROCS, MAX_PROCS);
        fflush(stderr);
        exit(-1);
    }
    printf("Using %d procs on %d steps of %d mols\n", NPROCS, NSTEP, NMOL);

    // Data Initialization Process
    // SET UP SCALING FACTORS AND CONSTANTS
    NORD1 = NORDER + 1;
    
    // Sub. call to set up constants
    CNSTNT(NORD1, TLC);
    
    fprintf(six, "\nTEMPERATURE              = %8.2f K\n",		TEMP);
    fprintf(six, "DENSITY                    = %8.5f G/C.C.\n",		RHO);
    fprintf(six, "NUMBER OF MOLECULES        = %8d\n",			NMOL);
    fprintf(six, "NUMBER OF PROCESSORS       = %8d\n",			NPROCS);
    fprintf(six, "TIME STEP                  = %8.2e SEC\n",		TSTEP);
    fprintf(six, "ORDER USED TO SOLVE F=MA   = %8d \n",			NORDER);
    fprintf(six, "NO. OF TIME STEPS          = %8d \n",			NSTEP);
    fprintf(six, "FREQUENCY OF DATA SAVING   = %8d \n",			NSAVE);
    fprintf(six, "FREQUENCY TO WRITE RST FILE= %8d \n",			NRST);

    // Create up to NPROCS number of threads
    for(i=0; i < NPROCS; i++)
    {
        param[i].threadID = i;
        pthread_create(&pthreads[i],NULL, thread_work, (void *) (&param[i]));
    }

    // Collect start time
    gettimeofday(&start, NULL);

        // Waiting for all threads to join before continue
        for(i=0; i < NPROCS; i++){
            rc = pthread_join(pthreads[i], (void **)&status);
            if (rc){
                printf("ERROR; return code from pthread_join() of thread %d is %d\n", i, rc);
                exit(-1);
            }
        }

    // Collect end time when all threads have returned
    gettimeofday(&end, NULL);
    
    fprintf(stderr, "\nExited Happily with XTT %g\n", XTT);
    fprintf(stderr, "Elapsed time: %.2f seconds\n",
    		(((end.tv_sec * 1000000.0) +  end.tv_usec) -
    		((start.tv_sec * 1000000.0) + start.tv_usec)) / 1000000.0);
    fflush(stderr);
	
    return 0;	
}// main.c

#endif //__WATER_WATER_IMPL__
