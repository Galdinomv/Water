
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
#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"
#include "data.h"

/* KINETI(NMOL,SUM,HMAS,OMAS,ProcID) */
/*   int NMOL; */
/*   double HMAS,OMAS; */
/*   double SUM[]; */
/*   unsigned ProcID; */
  
/*   /\* this routine computes kinetic energy in each of the three */
/*      spatial dimensions, and puts the computed values in the */
/*      SUM array *\/  */
void KINETI( tm_double SUM[3], double HMAS, double OMAS, unsigned ProcID)
{
    int dir, i, j, k;
    double S;
    p_tm_link curr_ptr;
    struct list_of_boxes *curr_box;
    double *tempptr;
    
    /* Loop over three directions */
    
    for (dir = XDIR; dir <= ZDIR; dir++) {
        S=0.0;
        
        /* loop over processor's boxes */
        curr_box = my_boxes[ProcID];
        
        while (curr_box) {
            
            i = curr_box->coord[XDIR];  /* X coordinate of box */
            j = curr_box->coord[YDIR];  /* Y coordinate of box */
            k = curr_box->coord[ZDIR];  /* Z coordinate of box */
            
            /* loop over the molecules */
            
            curr_ptr = BOX[i][j][k].list;
            while (curr_ptr) {
	      tempptr = curr_ptr->mol.F[VEL][dir];
                S += (tempptr[H1] * tempptr[H1] +
                      tempptr[H2] * tempptr[H2] ) * HMAS +
                          (tempptr[O] * tempptr[O]) * OMAS;
                curr_ptr = curr_ptr->next_mol;
            } /* while curr_ptr */
            
            curr_box = curr_box->next_box;
            
        } /* while curr_box */ 
        
	INIT_TRANSACTIONS();
	BEGIN_TRANSACTION();
	//pthread_mutex_lock(&KinetiSumLock);
        SUM[dir]  = SUM[dir] + S;
	//pthread_mutex_unlock(&KinetiSumLock);
	COMMIT_TRANSACTION();
        
    } /* for dir */
    
} /* end of subroutine KINETI */
