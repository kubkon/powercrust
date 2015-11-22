/* io.c : input-output */


/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
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


double mult_up = 1.0;
Coord mins[MAXDIM]
= {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX},
    maxs[MAXDIM]
    = {-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX};

void panic(char *fmt, ...) {
    va_list args;

    va_start(args, fmt);
    vfprintf(DFILE, fmt, args);
    fflush(DFILE);
    va_end(args);

    exit(1);
}

/*
FILE * popen(char*, char*);
void pclose(FILE*);
*/

char tmpfilenam[L_tmpnam];

FILE* efopen(char *file, char *mode) {
    FILE* fp;
    if ((fp = fopen(file, mode))!=NULL) return fp;
    fprintf(DFILE, "couldn't open file %s mode %s\n",file,mode);
    exit(1);
    return NULL;
}

void efclose(FILE* file) {
    fclose(file);
    file = NULL;
}

FILE* epopen(char *com, char *mode) {
    FILE* fp;
    if ((fp = popen(com, mode))!=NULL) return fp;
    fprintf(stderr, "couldn't open stream %s mode %s\n",com,mode);
    exit(1);
    return 0;
}



void print_neighbor_snum(FILE* F, neighbor *n){

    assert(site_num!=NULL);
    if (n->vert)
        fprintf(F, "%ld ", (*site_num)(n->vert));
    else
        fprintf(F, "NULL vert ");
    fflush(stdout);
}

void print_basis(FILE *F, basis_s* b) {
    if (!b) {fprintf(F, "NULL basis ");fflush(stdout);return;}
    if (b->lscale<0) {fprintf(F, "\nbasis computed");return;}
    fprintf(F, "\n%p  %d \n b=",(void*)b,b->lscale);
    print_point(F, rdim, b->vecs);
    fprintf(F, "\n a= ");
    print_point_int(F, rdim, b->vecs+rdim); fprintf(F, "   ");
    fflush(F);
}

void print_simplex_num(FILE *F, simplex *s) {
    fprintf(F, "simplex ");
    if(!s) fprintf(F, "NULL ");
    else fprintf(F, "%p  ", (void*)s);
}

void print_neighbor_full(FILE *F, neighbor *n){

    if (!n) {fprintf(F, "null neighbor\n");return;}

    print_simplex_num(F, n->simp);
    print_neighbor_snum(F, n);fprintf(F, ":  ");fflush(F);
    if (n->vert) {
        /*      if (n->basis && n->basis->lscale <0) fprintf(F, "trans ");*/
        /* else */ print_point(F, pdim,n->vert);
        fflush(F);
    }
    print_basis(F, n->basis);
    fflush(F);
    fprintf(F, "\n");
}

void *print_facet(FILE *F, simplex *s, print_neighbor_f *pnfin) {
    int i;
    neighbor *sn = s->neigh;

    /*  fprintf(F, "%d ", s->mark);*/
    for (i=0; i<cdim;i++,sn++) (*pnfin)(F, sn);
    fprintf(F, "\n");
    fflush(F);
    return NULL;
}




void *print_simplex_f(simplex *s, FILE *F, print_neighbor_f *pnfin){

    static print_neighbor_f *pnf;

    if (pnfin) {pnf=pnfin; if (!s) return NULL;}

    print_simplex_num(F, s);
    fprintf(F, "\n");
    if (!s) return NULL;
    fprintf(F, "normal ="); print_basis(F, s->normal); fprintf(F, "\n");
    fprintf(F, "peak ="); (*pnf)(F, &(s->peak));
    fprintf (F, "facet =\n");fflush(F);
    return print_facet(F, s, pnf);
}

void *print_simplex(simplex *s, void *Fin) {
    static FILE *F;

    if (Fin) {F=(FILE*)Fin; if (!s) return NULL;}

    return print_simplex_f(s, F, 0);

}


void print_triang(simplex *root, FILE *F, print_neighbor_f *pnf) {
    print_simplex(0,F);
    print_simplex_f(0,0,pnf);
    visit_triang(root, print_simplex);
}

void *p_peak_test(simplex *s) {return (s->peak.vert==p) ? (void*)s : (void*)NULL;}



void *check_simplex(simplex *s, void *dum){

    int i,j,k,l;
    neighbor *sn, *snn, *sn2;
    simplex *sns;
    site vn;

    for (i=-1,sn=s->neigh-1;i<cdim;i++,sn++) {
        sns = sn->simp;
        if (!sns) {
            fprintf(DFILE, "check_triang; bad simplex\n");
            print_simplex_f(s, DFILE, &print_neighbor_full);
            fprintf(DFILE, "site_num(p)=%ld\n",  site_num(p));
            return s;
        }
        if (!s->peak.vert && sns->peak.vert && i!=-1) {
            fprintf(DFILE, "huh?\n");
            print_simplex_f(s, DFILE, &print_neighbor_full);
            print_simplex_f(sns, DFILE, &print_neighbor_full);
            exit(1);
        }
        for (j=-1,snn=sns->neigh-1; j<cdim && snn->simp!=s; j++,snn++);
        if (j==cdim) {
            fprintf(DFILE, "adjacency failure:\n");
            DEBEXP(-1,site_num(p))
                print_simplex_f(sns, DFILE, &print_neighbor_full);
            print_simplex_f(s, DFILE, &print_neighbor_full);
            exit(1);
        }
        for (k=-1,snn=sns->neigh-1; k<cdim; k++,snn++){
            vn = snn->vert;
            if (k!=j) {
                for (l=-1,sn2 = s->neigh-1;
                     l<cdim && sn2->vert != vn;
                     l++,sn2++);
                if (l==cdim) {
                    fprintf(DFILE, "cdim=%d\n",cdim);
                    fprintf(DFILE, "error: neighboring simplices with incompatible vertices:\n");
                    print_simplex_f(sns, DFILE, &print_neighbor_full);
                    print_simplex_f(s, DFILE, &print_neighbor_full);
                    exit(1);
                }
            }
        }
    }
    return NULL;
}

int p_neight(simplex *s, int i, void *dum) {return s->neigh[i].vert !=p;}

void check_triang(simplex *root){visit_triang(root, &check_simplex);}

void check_new_triangs(simplex *s){visit_triang_gen(s, check_simplex, p_neight);}


/* outfuncs: given a list of points, output in a given format */

void vv_out(point *v, int vdim, FILE *Fin, int amble) {
    /* sunghee */
    static FILE *F;
    int i,j;

    if (Fin) {F=Fin; if (!v) return;}
    for (j=0;j<vdim;j++) {
        for (i=0;i<3;i++) {
            fprintf(F, "%G ", v[j][i]/mult_up);
        }
        fprintf(F, " | ");
    }
    fprintf(F, "\n");
    return;
}

void vlist_out(point *v, int vdim, FILE *Fin, int amble) {

    static FILE *F;
    int j;

    if (Fin) {F=Fin; if (!v) return;}

    for (j=0;j<vdim;j++) fprintf(F, "%ld ", (site_num)(v[j]));
    fprintf(F,"\n");

    return;
}

void off_out(point *v, int vdim, FILE *Fin, int amble) {

    static FILE *F, *G;
    static FILE *OFFFILE;
    static char offfilenam[L_tmpnam];
    char comst[100], buf[200];
    int j,i;

    if (Fin) {F=Fin;}

    if (pdim!=3) { warning(-10, off apparently for 3d points only); return;}

    if (amble==0) {
        for (i=0;i<vdim;i++) if (v[i]==coordsAtInfinity) return;
        fprintf(OFFFILE, "%d ", vdim);
        for (j=0;j<vdim;j++) fprintf(OFFFILE, "%ld ", (site_num)(v[j]));
        fprintf(OFFFILE,"\n");
    } else if (amble==-1) {
        OFFFILE = efopen(tmpnam(offfilenam), "w");
    } else {
        efclose(OFFFILE);

        fprintf(F, "    OFF\n");

        sprintf(comst, "wc %s", tmpfilenam);
        G = epopen(comst, "r");
        fscanf(G, "%d", &i);
        fprintf(F, " %d", i);
        pclose(G);

        sprintf(comst, "wc %s", offfilenam);
        G = epopen(comst, "r");
        fscanf(G, "%d", &i);
        fprintf(F, " %d", i);
        pclose(G);

        fprintf (F, " 0\n");

        G = efopen(tmpfilenam, "r");
        while (fgets(buf, sizeof(buf), G)) fprintf(F, "%s", buf);
        efclose(G);

        G = efopen(offfilenam, "r");


        while (fgets(buf, sizeof(buf), G)) fprintf(F, "%s", buf);
        efclose(G);
    }

}



void mp_out(point *v, int vdim, FILE *Fin, int amble) {


    /* should fix scaling */

    static int figno=1;
    static FILE *F;

    if (Fin) {F=Fin;}

    if (pdim!=2) { warning(-10, mp for planar points only); return;}
    if (amble==0) {
        int i;
        if (!v) return;
        for (i=0;i<vdim;i++) if (v[i]==coordsAtInfinity) {
            point t=v[i];
            v[i]=v[vdim-1];
            v[vdim-1] = t;
            vdim--;
            break;
        }
        fprintf(F, "draw ");
        for (i=0;i<vdim;i++)
            fprintf(F,
                    (i+1<vdim) ? "(%Gu,%Gu)--" : "(%Gu,%Gu);\n",
                    v[i][0]/mult_up,v[i][1]/mult_up
                );
    } else if (amble==-1) {
        if (figno==1) fprintf(F, "u=1pt;\n");
        fprintf(F , "beginfig(%d);\n",figno++);
    } else if (amble==1) {
        fprintf(F , "endfig;\n");
    }
}


void ps_out(point *v, int vdim, FILE *Fin, int amble) {

    static FILE *F;
    static double scaler;

    if (Fin) {F=Fin;}

    if (pdim!=2) { warning(-10, ps for planar points only); return;}

    if (amble==0) {
        int i;
        if (!v) return;
        for (i=0;i<vdim;i++) if (v[i]==coordsAtInfinity) {
            point t=v[i];
            v[i]=v[vdim-1];
            v[vdim-1] = t;
            vdim--;
            break;
        }
        fprintf(F,
                "newpath %G %G moveto\n",
                v[0][0]*scaler,v[0][1]*scaler);
        for (i=1;i<vdim;i++)
            fprintf(F,
                    "%G %G lineto\n",
                    v[i][0]*scaler,v[i][1]*scaler
                );
        fprintf(F, "stroke\n");
    } else if (amble==-1) {
        float len[2], maxlen;
        fprintf(F, "%%!PS\n");
        len[0] = maxs[0]-mins[0]; len[1] = maxs[1]-mins[1];
        maxlen = (len[0]>len[1]) ? len[0] : len[1];
        scaler = 216/maxlen;

        fprintf(F, "%%%%BoundingBox: %G %G %G %G \n",
                mins[0]*scaler,
                mins[1]*scaler,
                maxs[0]*scaler,
                maxs[1]*scaler);
        fprintf(F, "%%%%Creator: hull program\n");
        fprintf(F, "%%%%Pages: 1\n");
        fprintf(F, "%%%%EndProlog\n");
        fprintf(F, "%%%%Page: 1 1\n");
        fprintf(F, " 0.5 setlinewidth [] 0 setdash\n");
        fprintf(F, " 1 setlinecap 1 setlinejoin 10 setmiterlimit\n");
    } else if (amble==1) {
        fprintf(F , "showpage\n %%%%EOF\n");
    }
}

void cpr_out(point *v, int vdim, FILE *Fin, int amble) {

    static FILE *F;
    int i;

    if (Fin) {F=Fin; if (!v) return;}

    if (pdim!=3) { warning(-10, cpr for 3d points only); return;}

    for (i=0;i<vdim;i++) if (v[i]==coordsAtInfinity) return;

    fprintf(F, "t %G %G %G %G %G %G %G %G %G 3 128\n",
            v[0][0]/mult_up,v[0][1]/mult_up,v[0][2]/mult_up,
            v[1][0]/mult_up,v[1][1]/mult_up,v[1][2]/mult_up,
            v[2][0]/mult_up,v[2][1]/mult_up,v[2][2]/mult_up
        );
}


/* vist_funcs for different kinds of output: facets, alpha shapes, etc. */



void *facets_print(simplex *s, void *p) {

    static out_func *out_func_here;
    point v[MAXDIM];
    int j;

    if (p) {out_func_here = (out_func*)p; if (!s) return NULL;}

    for (j=0;j<cdim;j++) v[j] = s->neigh[j].vert;

    out_func_here(v,cdim,0,0);

    return NULL;
}

void *ridges_print(simplex *s, void *p) {

    static out_func *out_func_here;
    point v[MAXDIM];
    int j,k,vnum;

    if (p) {out_func_here = (out_func*)p; if (!s) return NULL;}

    for (j=0;j<cdim;j++) {
        vnum=0;
        for (k=0;k<cdim;k++) {
            if (k==j) continue;
            v[vnum++] = (s->neigh[k].vert);
        }
        out_func_here(v,cdim-1,0,0);
    }
    return NULL;
}


void *afacets_print(simplex *s, void *p) {

    static out_func *out_func_here;
    point v[MAXDIM];
    int j,k,vnum;

    if (p) {out_func_here = (out_func*)p; if (!s) return NULL;}

    for (j=0;j<cdim;j++) { /* check for ashape consistency */
        for (k=0;k<cdim;k++) if (s->neigh[j].simp->neigh[k].simp==s) break;
        if (alph_test(s,j,0)!=alph_test(s->neigh[j].simp,k,0)) {
            DEB(-10,alpha-shape not consistent)
                DEBTR(-10)
                print_simplex_f(s,DFILE,&print_neighbor_full);
            print_simplex_f(s->neigh[j].simp,DFILE,&print_neighbor_full);
            fflush(DFILE);
            exit(1);
        }
    }
    for (j=0;j<cdim;j++) {
        vnum=0;
        if (alph_test(s,j,0)) continue;
        for (k=0;k<cdim;k++) {
            if (k==j) continue;
            v[vnum++] = s->neigh[k].vert;
        }
        out_func_here(v,cdim-1,0,0);
    }
    return NULL;
}
