/*
 * Power Crust software, by Nina Amenta, Sunghee Choi and Ravi Krishna Kolluri.
 * Copyright (c) 2000 by the University of Texas
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee under the GNU Public License is hereby granted, 
 * provided that this entire notice  is included in all copies of any software 
 * which is or includes a copy or modification of this software and in all copies 
 * of the supporting documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include "hull.h"


extern struct heap_array *heap_A;
extern int heap_length,heap_size; 
extern struct polelabel *adjlist;
extern struct plist **opplist;
extern double theta, deep;
extern int defer;

int propagate()
{
    int pid;

    pid = extract_max();
    if (adjlist[pid].in > adjlist[pid].out) 
        adjlist[pid].label = IN;
    else adjlist[pid].label = OUT;
    /*  fprintf(DFILE,"pole %d in %f out %f label %d\n",pid, adjlist[pid].in, adjlist[pid].out, adjlist[pid].label);  */
    
    if (pid != -1) {
        /*    fprintf(DFILE,"propagating pole %d..\n",pid);  */
        opp_update(pid);
        sym_update(pid);
    }
    return pid;
}

void opp_update(int pi)
{
    struct plist *pindex;
    int npi, nhi;
    double temp;

    pindex = opplist[pi];
    while (pindex!=NULL) { 
        npi = pindex->pid;
        if (defer) {
            if (adjlist[npi].bad == BAD_POLE) {
                /*  fprintf(DFILE,"found bad pole.. defer its labeling\n"); */
                pindex = pindex->next;
                continue; 
            } 
        }
        if (adjlist[npi].label == INIT) { /* not yet labeled */
            if (adjlist[npi].hid == 0) { /* not in the heap */
                if (adjlist[pi].in > adjlist[pi].out) {
                    /* propagate in*cos to out */
                    adjlist[npi].out = (-1.0) * adjlist[pi].in * pindex->angle;
                    /*  fprintf(DFILE,"pole %d.out = %f\n",npi,adjlist[npi].out); */
                    insert_heap(npi,adjlist[npi].out);
                }
                else if (adjlist[pi].in < adjlist[pi].out) {
                    /* propagate out*cos to in */
                    adjlist[npi].in = (-1.0) * adjlist[pi].out * pindex->angle;
                    /* fprintf(DFILE,"pole %d.in = %f\n",npi,adjlist[npi].in); */
                    insert_heap(npi,adjlist[npi].in);
                }
            }
            else { /* in the heap */
                nhi = adjlist[npi].hid;
                if (adjlist[pi].in > adjlist[pi].out) {
                    /* propagate in*cos to out */
                    temp = (-1.0) * adjlist[pi].in * pindex->angle;
                    if (temp > adjlist[npi].out) { 
                        adjlist[npi].out = temp;
                        update_pri(nhi,npi);
                    }
                }
                else if (adjlist[pi].in < adjlist[pi].out) {
                    /* propagate out*cos to in */
                    temp = (-1.0) * adjlist[pi].out * pindex->angle;
                    if (temp > adjlist[npi].in) {
                        adjlist[npi].in = temp;
                        update_pri(nhi,npi);
                    }
                }
            }
        }
        pindex = pindex->next;
    }
}

void sym_update(int pi)
{
    struct edgesimp *eindex;
    int npi, nhi;
    double temp;

    eindex = adjlist[pi].eptr;
    while (eindex!=NULL) {
        npi = eindex->pid;
        if (defer) {
            if (adjlist[npi].bad == BAD_POLE) {
                eindex = eindex->next;
                /*fprintf(DFILE,"found bad pole.. defer its labeling\n");*/
                continue;
            }
        }

        /* try to label deeply intersecting unlabeled neighbors */
        if  ((adjlist[npi].label==INIT) && (eindex->angle > theta)) { 
            /* not yet labeled */
            if (adjlist[npi].hid == 0) { /* not in the heap */
                if (adjlist[pi].in > adjlist[pi].out) {
                    /* propagate in*cos to in */
                    adjlist[npi].in = adjlist[pi].in * eindex->angle;
                    insert_heap(npi,adjlist[npi].in);
                }
                else if (adjlist[pi].in < adjlist[pi].out) {
                    /* propagate out*cos to out */
                    adjlist[npi].out = adjlist[pi].out * eindex->angle;
                    insert_heap(npi,adjlist[npi].out);
                }
            }
            else { /* in the heap */
                if (heap_A[adjlist[npi].hid].pid != npi) fprintf(DFILE,"ERROR\n");
                nhi = adjlist[npi].hid;
                if (adjlist[pi].in > adjlist[pi].out) {
                    /* propagate in*cos to in */
                    temp = adjlist[pi].in * eindex->angle;
                    if (temp > adjlist[npi].in) { 
                        adjlist[npi].in = temp;
                        update_pri(nhi,npi);
                    }
                }
                else if (adjlist[pi].in < adjlist[pi].out) {
                    /* propagate out*cos to out */
                    temp = adjlist[pi].out * eindex->angle;
                    if (temp > adjlist[npi].out) {
                        adjlist[npi].out = temp;
                        update_pri(nhi,npi);
                    }
                }
            }
        }
        eindex = eindex->next;
    }
}

void update_pri(int hi, int pi) 
    /* update heap_A[hi].pri using adjlist[pi].in/out */
{
    double pr;

    if ((heap_A[hi].pid != pi)||(adjlist[pi].hid != hi)) {
        fprintf(DFILE,"Error update_pri!\n");
        return;
    }
    if (adjlist[pi].in==0.0) {
        pr = adjlist[pi].out;
    }
    else if (adjlist[pi].out == 0.0) {
        pr = adjlist[pi].in;
    }
    else { /* both in/out nonzero */
        if (adjlist[pi].in > adjlist[pi].out) {
            pr =  adjlist[pi].in - adjlist[pi].out - 1;
        }
        else {
            pr = adjlist[pi].out - adjlist[pi].in - 1;
        }
    }
    update(hi,pr);
}

void label_unlabeled(int num)
{
    struct plist *pindex;
    struct edgesimp *eindex;
    int npi,i, opplabel;
    int tlabel;
    double tangle, tangle1;


    for (i=0;i<num;i++) { 
        if (adjlist[i].label == INIT) { /* pole i is unlabeled.. try to label now */
            tlabel = INIT;
            opplabel = INIT;
            pindex = opplist[i];
            if ((pindex == NULL)&&(adjlist[i].eptr==NULL)) {
                fprintf(DFILE,"no opp pole, no adjacent pole!\n");
                continue;
            }
            /* check whether there is opp pole */
            while (pindex!=NULL) { /* opp pole */
                npi = pindex->pid;
                if (adjlist[npi].label != INIT) {
                    fprintf(DFILE,"opp is labeled\n");
                    if (opplabel == INIT) opplabel = adjlist[npi].label;
                    else if (opplabel != adjlist[npi].label) {
                        /* opp poles have different labels ... inconsistency! */
                        fprintf(DFILE,"opp poles have inconsistent labels\n");
                        opplabel = INIT; /* ignore the label of opposite poles */
                    } 
                }
                pindex = pindex->next;
            }

            tangle = -3.0;
            tangle1 = -3.0;
            eindex = adjlist[i].eptr;
            while (eindex != NULL) {
                npi = eindex->pid; 
                if (adjlist[npi].label == IN) {
                    if (tangle < eindex->angle) {
                        tangle = eindex->angle;
                    }
                }
                else if (adjlist[npi].label == OUT) {
                    if (tangle1 < eindex->angle) {
                        tangle1 = eindex->angle;
                    }
                }
                eindex = eindex->next;
            }
            /* now tangle, tangle 1 are angles of most deeply interesecting in, out poles */
            if (tangle == -3.0) { /* there was no in poles */
                if (tangle1 == -3.0) { /* there was no out poles */
                    if (opplabel == INIT) { /* cannot trust opp pole or no labeled opp pole */
                        fprintf(DFILE, "1: cannot label pole %d\n", i);
                    }
                    else if (opplabel == IN) {
                        adjlist[i].label = OUT;
                    }
                    else  adjlist[i].label = IN;
                }
                else if (tangle1 > deep) { /* interesecting deeply only out poles */
                    adjlist[i].label = OUT;
                }
                else { /* no deeply intersecting poles . use opp pole */
                    if (opplabel == INIT) { /* cannot trust opp pole or no labeled opp pole */
                        fprintf(DFILE, "2: cannot label pole %d\n", i);
                    }
                    else if (opplabel == IN) {
                        adjlist[i].label = OUT;
                    }
                    else  adjlist[i].label = IN;
                }
            }
            else if (tangle1 == -3.0) { /* there are in pole but no out pole */
                if (tangle > deep) { /* interesecting deeply only in poles */
                    adjlist[i].label = IN;
                }
                else { /* no deeply intersecting poles . use opp pole */
                    if (opplabel == INIT) { /* cannot trust opp pole or no labeled opp pole */
                        fprintf(DFILE, "2: cannot label pole %d\n", i);
                    }
                    else if (opplabel == IN) {
                        adjlist[i].label = OUT;
                    }
                    else  adjlist[i].label = IN;
                }
            }
            else { /* there are both in/out poles */
                if (tangle > deep) {
                    if (tangle1 > deep) { /* intersecting both deeply */
                        /* use opp */
                        if (opplabel == INIT) { /* cannot trust opp pole or no labeled opp pole */
                            /* then give label with bigger angle */
                            fprintf(DFILE,"intersect both deeply but no opp,in %f out %f.try to label more deeply intersected label.\n", tangle, tangle1);
                            if (tangle > tangle1) {
                                adjlist[i].label = IN;
                            }
                            else adjlist[i].label = OUT;
                        }
                        else if (opplabel == IN) {
                            adjlist[i].label = OUT;
                        }
                        else  adjlist[i].label = IN;
                    }
                    else { /* intersecting only in deeply */
                        adjlist[i].label = IN;
                    }
                }
                else if (tangle1 > deep) { /* intersecting only out deeply */
                    adjlist[i].label = OUT;
                }
                else { /* no deeply intersecting poles . use opp pole */
                    if (opplabel == INIT) { /* cannot trust opp pole or no labeled opp pole */
                        fprintf(DFILE, "3: cannot label pole %d\n", i); 
                    }
                    else if (opplabel == IN) {
                        adjlist[i].label = OUT;
                    }
                    else  adjlist[i].label = IN;
                }
            } 



 
            /* no labeled opp pole - label pole same as the most deeply intersecting labeled pole ... no longer needed because opp values are already propagated..
               tangle = -3.0;
               eindex = adjlist[i].eptr;
               while (eindex != NULL) {
               npi = eindex->pid; 
               if (adjlist[npi].label == IN) {
               if (tangle < eindex->angle) {
               tangle = eindex->angle;
               tlabel = IN;
               }
               }
               else if (adjlist[npi].label == OUT) {
               if (tangle < eindex->angle) {
               tangle = eindex->angle;
               tlabel = OUT;
               }
               }
               eindex = eindex->next;
               }
               fprintf(DFILE,"pole %d  max angle %f label %d\n", i,tangle,tlabel);
               adjlist[i].label = tlabel;    
            */


        }
    }
}
