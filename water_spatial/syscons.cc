
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

 
#include <stdio.h>
#include <math.h>


#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"

int    BOX_PER_SIDE;

void SYSCNS()                    /* sets up some system constants */
{
    TSTEP=TSTEP/UNITT;        /* time between steps */
    NATMO=NATOMS*NMOL;        /* total number of atoms in system */
    NATMO3=NATMO*3; /* number of atoms * number of spatial dimensions */
    double fpot, unitm, unitl, unitt, boltz, temp, natmo;
    unitm = UNITM; 
    unitl = UNITL;
    unitt = UNITT;
    boltz = BOLTZ;
    temp = TEMP;
    natmo = NATMO;
    FPOT= unitm * pow((unitl/unitt),2.0) / (boltz*temp*natmo);
    //printf("ANMOL ANMOL natmo = %28.27lf\n\n\n\n", natmo);
    //printf("ANMOL ANMOL unitm = %28.27lf\n\n\n\n", unitm);	
    //printf("ANMOL ANMOL FPOT = %20.19lf\n\n\n\n",FPOT);
    FKIN=FPOT*0.50/(TSTEP*TSTEP);
    
    /* computed length of the cubical "box".  Note that box size is
     * computed as being large enough to handle the input number of 
     * water molecules
     */
    
    //printf("NMOL = %d\n",NMOL);
    //printf("WTMOL = %f\n",WTMOL);
    //printf("UNITM = %26.25lf\n",UNITM);
    //printf("RHO = %f\n",RHO);

    double nmol, wtmol, rho;
    nmol =NMOL;
    wtmol = WTMOL;
    rho = RHO;
    unitm = UNITM;
    BOXL= pow( (nmol*wtmol*unitm/rho),(1.00/3.00));  
    
    /* normalized length of computational box (in Angstroms) */
    
    BOXL=BOXL/UNITL;    
    
    /* # of boxes per side */
    
    BOXH = BOXL*0.50;
   
    /* set cutoff radius if it was not already read in nonzero from
       the input file in water.C.  If it was defined in the input file,
       then that definition is retained. */
    
    if (CUTOFF == 0.0) {
        CUTOFF=max(BOXH,CUTOFF);    /* cutoff radius is max of BOXH 
                                       and default (= 0); i.e. CUTOFF
                                       radius is set to half the normalized
                                       box length */
// printf("CUTOFFFFF = %26.25lf\n",CUTOFF);
    }
    
    if (CUTOFF > 11.0) CUTOFF = 11.0; /* cutoff never greater than 11 
                                         Angstrom */
   
    BOX_PER_SIDE = (int)(BOXL/CUTOFF);      
 
    /* BOX_PER_SIDE is always >=1 */
    if (BOX_PER_SIDE == 0) BOX_PER_SIDE = 1;

    printf("BOX_PER_SIDE = %d\n",BOX_PER_SIDE);
    printf("BOXL = %f\n",BOXL);
    printf("CUTOFF = %f\n",CUTOFF);
    
    /* Length of a box in Angstroms */
    
    BOX_LENGTH = BOXL / ( (double) BOX_PER_SIDE);
    
    BPS_SQRD = BOX_PER_SIDE * BOX_PER_SIDE;
    
    NumBoxes = BOX_PER_SIDE * BOX_PER_SIDE * BOX_PER_SIDE;
    
    REF1= -QQ/(CUTOFF*CUTOFF*CUTOFF);
    REF2=2.00*REF1;
    REF4=2.00*REF2;
    CUT2=CUTOFF*CUTOFF;        /* square of cutoff radius,  used 
                                  to actually decide whether an 
                                  interaction should be computed in 
                                  INTERF and POTENG */
    
    FHM=(TSTEP*TSTEP*0.50)/HMAS;
    FOM=(TSTEP*TSTEP*0.50)/OMAS;
    NMOL1=NMOL-1;
    
}       /* end of subroutine SYSCNS */


