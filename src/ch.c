/* ch.c : numerical functions for hull computation */

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



#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string.h>


#include "hull.h"

extern short power_diagram;

short check_overshoot_f=0;


simplex *ch_root;

#define NEARZERO(d) ((d) < FLT_EPSILON && (d) > -FLT_EPSILON)
#define SMALL (100*FLT_EPSILON)*(100*FLT_EPSILON)

#define SWAP(X,a,b) {X t; t = a; a = b; b = t;}

#define DMAX

double Huge;

#define check_overshoot(x)                          \
    {if (CHECK_OVERSHOOT && check_overshoot_f && ((x)>9e15))        \
        warning(-20, overshot exact arithmetic)}            \


#define DELIFT 0
int basis_vec_size;


#define lookupshort(a,b,whatb,c,whatc)                  \
{                                   \
    int i;                              \
    neighbor *x;                            \
    c = NULL;                           \
    for (i=-1, x = a->neigh-1; (x->whatb != b) && (i<cdim) ; i++, x++);\
    if (i<cdim) c = x->whatc;                   \
}                                   \


Coord Vec_dot(point x, point y) {
    int i;
    Coord sum = 0;
    for (i=0;i<rdim;i++) sum += x[i] * y[i];
    return sum;
}

Coord Vec_dot_pdim(point x, point y) {
    int i;
    Coord sum = 0;
    for (i=0;i<pdim;i++) sum += x[i] * y[i];
    /*  check_overshoot(sum); */
    return sum;
}

Coord Norm2(point x) {
    int i;
    Coord sum = 0;
    for (i=0;i<rdim;i++) sum += x[i] * x[i];
    return sum;
}

void Ax_plus_y(Coord a, point x, point y) {
    int i;
    for (i=0;i<rdim;i++) {
        *y++ += a * *x++;
    }
}

void Ax_plus_y_test(Coord a, point x, point y) {
    int i;
    for (i=0;i<rdim;i++) {
        check_overshoot(*y + a * *x);
        *y++ += a * *x++;
    }
}

void Vec_scale(int n, Coord a, Coord *x)
{
    register Coord *xx = x,
        *xend = xx + n;
    while (xx!=xend) *xx++ *= a;
}

void Vec_scale_test(int n, Coord a, Coord *x)
{
    register Coord *xx = x,
        *xend = xx + n  ;

    check_overshoot(a);

    while (xx!=xend) {
        *xx *= a;
        check_overshoot(*xx);
        xx++;
    }
}



int exact_bits;
float b_err_min, b_err_min_sq;

double logb(double); /* on SGI machines: returns floor of log base 2 */

static short vd;
static basis_s  tt_basis = {0,1,-1,0,0,{0}};
static basis_s *tt_basisp = &tt_basis, *infinity_basis;


STORAGE(basis_s)

    typedef Coord site_struct;

    Coord   coordsAtInfinity[10]={57.2,0,0,0,0}; /* point at infinity for vd; value not used */

    void print_site(site p, FILE* F)
{print_point(F, pdim,p);fprintf(F, "\n");}

#define VA(x) ((x)->vecs+rdim)
#define VB(x) ((x)->vecs)




/* tables for runtime stats */
int A[100]={0}, B[100] ={0}, C[100] = {0}, D[100] ={0};
int tot =0, totinf=0, bigt=0;

#define two_to(x) ( ((x)<20) ? 1<<(x) : ldexp(1,(x)) )



double sc(basis_s *v,simplex *s, int k, int j) {
    /* amount by which to scale up vector, for reduce_inner */

    double      labound;
    static int  lscale;
    static double   max_scale,
        ldetbound,
        Sb;

    if (j<10) {
        labound = logb(v->sqa)/2;
        max_scale = exact_bits - labound - 0.66*(k-2)-1  -DELIFT;
        if (max_scale<1) {
            warning(-10, overshot exact arithmetic);
            max_scale = 1;

        }

        if (j==0) {
            int i;
            neighbor *sni;
            basis_s *snib;

            ldetbound = DELIFT;

            Sb = 0;
            for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
                snib = sni->basis;
                Sb += snib->sqb;
                ldetbound += logb(snib->sqb)/2 + 1;
                ldetbound -= snib->lscale;
            }
        }
    }
    if (ldetbound - v->lscale + logb(v->sqb)/2 + 1 < 0) {
        DEBS(-2)
            DEBTR(-2) DEBEXP(-2, ldetbound)
            print_simplex_f(s, DFILE, &print_neighbor_full);
        print_basis(DFILE,v);
        EDEBS
            return 0;
    } else {
        lscale = logb(2*Sb/(v->sqb + v->sqa*b_err_min))/2;
        if (lscale > max_scale) {
            lscale = max_scale;
        } else if (lscale<0) lscale = 0;
        v->lscale += lscale;
        return two_to(lscale);
    }
}


double lower_terms(basis_s* v) {

    point vp = v->vecs;
    int i,j,h,hh=0;
    int facs[6] = {2,3,5,7,11,13};
    double out = 1;

    DEBTR(-10) print_basis(DFILE, v); printf("\n");
    DEBTR(0)

        for (j=0;j<6;j++) do {
            for (i=0; i<2*rdim && facs[j]*floor(vp[i]/facs[j])==vp[i];i++);
            if ((h = (i==2*rdim))!=0) {
                hh=1;
                out *= facs[j];
                for (i=0;i<2*rdim; i++) vp[i]/=facs[j];
            }
        } while (h);
    /*  if (hh) {DEBTR(-10)  print_basis(DFILE, v);} */
    return out;
}

double lower_terms_point(point vp) {

    int i,j,h,hh=0;
    int facs[6] = {2,3,5,7,11,13};
    double out = 1;

    for (j=0;j<6;j++) do {
        for (i=0; i<2*rdim && facs[j]*floor(vp[i]/facs[j])==vp[i];i++);
        if ((h = (i==2*rdim))!=0) {
            hh=1;
            out *= facs[j];
            for (i=0;i<2*rdim; i++) vp[i]/=facs[j];
        }
    } while (h);
    return out;
}


int reduce_inner(basis_s *v, simplex *s, int k) {

    point   va = VA(v),
        vb = VB(v);
    int i,j;
    double  dd;
    double  scale;
    basis_s *snibv;
    neighbor *sni;
    static int failcount;

    /*  lower_terms(v); */
    v->sqa = v->sqb = Norm2(vb);
    if (k<=1) {
        memcpy(vb,va,basis_vec_size);
        return 1;
    }
    /*  if (vd) {
        snibv = s->neigh[1].basis;
        scale = floor(sqrt(snibv->sqa/v->sqa));
        if (scale > 1) Vec_scale(rdim,scale,va);
        dd = Vec_dot(VA(snibv),va)/snibv->sqa;
        dd = -floor(0.5+dd);
        Ax_plus_y( dd, VA(snibv), va);
        }
    */
    for (j=0;j<250;j++) {

        memcpy(vb,va,basis_vec_size);
        for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
            snibv = sni->basis;
            dd = -Vec_dot(VB(snibv),vb)/ snibv->sqb;
            Ax_plus_y( dd, VA(snibv), vb);
        }
        v->sqb = Norm2(vb);
        v->sqa = Norm2(va);

        if (2*v->sqb >= v->sqa) {B[j]++; return 1;}

        Vec_scale_test(rdim,scale = sc(v,s,k,j),va);

        for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
            snibv = sni->basis;
            dd = -Vec_dot(VB(snibv),va)/snibv->sqb;
            dd = floor(dd+0.5);
            Ax_plus_y_test( dd, VA(snibv), va);
        }
    }
    if (failcount++<10) {
        DEB(-8, reduce_inner failed on:)
            DEBTR(-8) print_basis(DFILE, v);
        print_simplex_f(s, DFILE, &print_neighbor_full);
    }
    return 0;
}

#define trans(z,p,q) {int i; for (i=0;i<pdim;i++) z[i+rdim] = z[i] = p[i] - q[i];}
#define lift(z,s) {if (vd) z[2*rdim-1] =z[rdim-1]= ldexp(Vec_dot_pdim(z,z), -DELIFT);}
/*not scaling lift to 2^-DELIFT */



int reduce(basis_s **v, point p, simplex *s, int k) {

    point   z;
    point   tt = s->neigh[0].vert;

    if (!*v) NEWLRC(basis_s,(*v))
                 else (*v)->lscale = 0;
    z = VB(*v);
    if (vd || power_diagram ) {
        if (p==coordsAtInfinity) memcpy(*v,infinity_basis,basis_s_size);
        else {trans(z,p,tt); lift(z,s);}
    } else trans(z,p,tt);
    return reduce_inner(*v,s,k);
}


void get_basis_sede(simplex *s) {

    int k=1;
    neighbor *sn = s->neigh+1,
        *sn0 = s->neigh;

    if (( vd || power_diagram) && sn0->vert == coordsAtInfinity && cdim >1) {
        SWAP(neighbor, *sn0, *sn );
        NULLIFY(basis_s,sn0->basis);
        sn0->basis = tt_basisp;
        tt_basisp->ref_count++;
    } else {
        if (!sn0->basis) {
            sn0->basis = tt_basisp;
            tt_basisp->ref_count++;
        } else while (k < cdim && sn->basis) {k++;sn++;}
    }
    while (k<cdim) {
        NULLIFY(basis_s,sn->basis);
        reduce(&sn->basis,sn->vert,s,k);
        k++; sn++;
    }
}


int out_of_flat(simplex *root, point p) {

    static neighbor p_neigh={0,0,0};

    if (!p_neigh.basis) p_neigh.basis = (basis_s*) malloc(basis_s_size);

    p_neigh.vert = p;
    cdim++;
    root->neigh[cdim-1].vert = root->peak.vert;
    NULLIFY(basis_s,root->neigh[cdim-1].basis);
    get_basis_sede(root);
    if ((vd || power_diagram)&& root->neigh[0].vert == coordsAtInfinity) return 1;
    reduce(&p_neigh.basis,p,root,cdim);
    if (p_neigh.basis->sqa != 0) return 1;
    cdim--;
    return 0;
}


double cosangle_sq(basis_s* v,basis_s* w)  {
    double dd;
    point   vv=v->vecs,
        wv=w->vecs;
    dd = Vec_dot(vv,wv);
    return dd*dd/Norm2(vv)/Norm2(wv);
}


int check_perps(simplex *s) {

    static basis_s *b = NULL;
    point   z,y;
    point   tt;
    double  dd;
    int i,j;

    for (i=1; i<cdim; i++) if (NEARZERO(s->neigh[i].basis->sqb)) return 0;
    if (!b) assert(b = (basis_s*)malloc(basis_s_size));
    else b->lscale = 0;
    z = VB(b);
    tt = s->neigh[0].vert;
    for (i=1;i<cdim;i++) {
        y = s->neigh[i].vert;
        if ( (vd || power_diagram)&& y==coordsAtInfinity) memcpy(b, infinity_basis, basis_s_size);
        else {trans(z,y,tt); lift(z,s);}
        if (s->normal && cosangle_sq(b,s->normal)>b_err_min_sq) {DEBS(0)
                                                                     DEB(0,bad normal) DEBEXP(0,i) DEBEXP(0,dd)
                                                                     print_simplex_f(s, DFILE, &print_neighbor_full);
        EDEBS
            return 0;
        }
        for (j=i+1;j<cdim;j++) {
            if (cosangle_sq(b,s->neigh[j].basis)>b_err_min_sq) {
                DEBS(0)
                    DEB(0,bad basis)DEBEXP(0,i) DEBEXP(0,j) DEBEXP(0,dd)
                    DEBTR(-8) print_basis(DFILE, b);
                print_simplex_f(s, DFILE, &print_neighbor_full);
                EDEBS
                    return 0;
            }
        }
    }
    return 1;
}


void get_normal_sede(simplex *s) {

    neighbor *rn;
    int i,j;

    get_basis_sede(s);
    if (rdim==3 && cdim==3) {
        point   c,
            a = VB(s->neigh[1].basis),
            b = VB(s->neigh[2].basis);
        NEWLRC(basis_s,s->normal);
        c = VB(s->normal);
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
        s->normal->sqb = Norm2(c);
        for (i=cdim+1,rn = ch_root->neigh+cdim-1; i; i--, rn--) {
            for (j = 0; j<cdim && rn->vert != s->neigh[j].vert;j++);
            if (j<cdim) continue;
            if (rn->vert==coordsAtInfinity) {
                if (c[2] > -b_err_min) continue;
            } else  if (!sees(rn->vert,s)) continue;
            c[0] = -c[0]; c[1] = -c[1]; c[2] = -c[2];
            break;
        }
        DEBS(-1) if (!check_perps(s)) exit(1); EDEBS
                                                   return;
    }

    for (i=cdim+1,rn = ch_root->neigh+cdim-1; i; i--, rn--) {
        for (j = 0; j<cdim && rn->vert != s->neigh[j].vert;j++);
        if (j<cdim) continue;
        reduce(&s->normal,rn->vert,s,cdim);
        if (s->normal->sqb != 0) break;
    }

    DEBS(-1) if (!check_perps(s)) {DEBTR(-1) exit(1);} EDEBS

                                                           }


void get_normal(simplex *s) {get_normal_sede(s); return;}

int sees(site p, simplex *s) {

    static basis_s *b = NULL;
    point   tt,zz;
    double  dd,dds;
    int i;


    if (!b) assert(b = (basis_s*)malloc(basis_s_size));
    else b->lscale = 0;
    zz = VB(b);
    if (cdim==0) return 0;
    if (!s->normal) {
        get_normal_sede(s);
        for (i=0;i<cdim;i++) NULLIFY(basis_s,s->neigh[i].basis);
    }
    tt = s->neigh[0].vert;
    if (vd || power_diagram) {
        if (p==coordsAtInfinity) memcpy(b,infinity_basis,basis_s_size);
        else {trans(zz,p,tt); lift(zz,s);}
    } else trans(zz,p,tt);
    for (i=0;i<3;i++) {
        dd = Vec_dot(zz,s->normal->vecs);
        if (dd == 0.0) {
            DEBS(-7) DEB(-6,degeneracy:); DEBEXP(-6,site_num(p));
            print_site(p, DFILE); print_simplex_f(s, DFILE, &print_neighbor_full); EDEBS
                                                                                       return 0;
        }
        dds = dd*dd/s->normal->sqb/Norm2(zz);
        if (dds > b_err_min_sq) return (dd<0);
        get_basis_sede(s);
        reduce_inner(b,s,cdim);
    }
    DEBS(-7) if (i==3) {
        DEB(-6, looped too much in sees);
        DEBEXP(-6,dd) DEBEXP(-6,dds) DEBEXP(-6,site_num(p));
        print_simplex_f(s, DFILE, &print_neighbor_full); exit(1);}
    EDEBS
        return 0;
}





double radsq(simplex *s) {

    point n;
    neighbor *sn;
    int i;

    /* square of ratio of circumcircle radius to
       max edge length for Delaunay tetrahedra */


    for (i=0,sn=s->neigh;i<cdim;i++,sn++)
        if (sn->vert == coordsAtInfinity) return Huge;

    if (!s->normal) get_normal_sede(s);

    /* compute circumradius */
    n = s->normal->vecs;

    if (NEARZERO(n[rdim-1])) {return Huge;}

    return Vec_dot_pdim(n,n)/4/n[rdim-1]/n[rdim-1];
}


void *zero_marks(simplex * s, void *dum) { s->mark = 0; return NULL; }

void *one_marks(simplex * s, void *dum) {s->mark = 1; return NULL;}

void *show_marks(simplex * s, void *dum) {printf("%d",s->mark); return NULL;}


#define swap_points(a,b) {point t; t=a; a=b; b=t;}

int alph_test(simplex *s, int i, void *alphap) {
    /*returns 1 if not an alpha-facet */

    simplex *si;
    double  rs,rsi,rsfi;
    neighbor *scn, *sin;
    int k, nsees, ssees;
    static double alpha;

    if (alphap) {alpha=*(double*)alphap; if (!s) return 1;}
    if (i==-1) return 0;

    si = s->neigh[i].simp;
    scn = s->neigh+cdim-1;
    sin = s->neigh+i;
    nsees = 0;

    for (k=0;k<cdim;k++) if (s->neigh[k].vert==coordsAtInfinity && k!=i) return 1;
    rs = radsq(s);
    rsi = radsq(si);

    if (rs < alpha &&  rsi < alpha) return 1;

    swap_points(scn->vert,sin->vert);
    NULLIFY(basis_s, s->neigh[i].basis);
    cdim--;
    get_basis_sede(s);
    reduce(&s->normal,coordsAtInfinity,s,cdim);
    rsfi = radsq(s);

    for (k=0;k<cdim;k++) if (si->neigh[k].simp==s) break;

    ssees = sees(scn->vert,s);
    if (!ssees) nsees = sees(si->neigh[k].vert,s);
    swap_points(scn->vert,sin->vert);
    cdim++;
    NULLIFY(basis_s, s->normal);
    NULLIFY(basis_s, s->neigh[i].basis);

    if (ssees) return alpha<rs;
    if (nsees) return alpha<rsi;

    assert(rsfi<=rs+FLT_EPSILON && rsfi<=rsi+FLT_EPSILON);

    return alpha<=rsfi;
}


void *conv_facetv(simplex *s, void *dum) {
    int i;
    for (i=0;i<cdim;i++) if (s->neigh[i].vert==coordsAtInfinity) {return s;}
    return NULL;
}

short mi[MAXPOINTS], mo[MAXPOINTS];

void *mark_points(simplex *s, void *dum) {
    int i, snum;
    neighbor *sn;

    for  (i=0,sn=s->neigh;i<cdim;i++,sn++) {
        if (sn->vert==coordsAtInfinity) continue;
        snum = site_num(sn->vert);
        if (s->mark) mo[snum] = 1;
        else mi[snum] = 1;
    }
    return NULL;
}

void* visit_outside_ashape(simplex *root, visit_func visit) {
    return visit_triang_gen(visit_hull(root, conv_facetv), visit, alph_test);
}

int check_ashape(simplex *root, double alpha) {

    int i;

    for (i=0;i<MAXPOINTS;i++) {mi[i] = mo[i] = 0;}

    visit_hull(root, zero_marks);

    alph_test(0,0,&alpha);
    visit_outside_ashape(root, one_marks);

    visit_hull(root, mark_points);

    for (i=0;i<MAXPOINTS;i++) if (mo[i] && !mi[i]) {return 0;}

    return 1;
}

double find_alpha(simplex *root) {

    int i;
    float al=0,ah,am;

    for (ah=i=0;i<pdim;i++) ah += (maxs[i]-mins[i])*(maxs[i]-mins[i]);
    assert(check_ashape(root,ah));
    for (i=0;i<17;i++) {
        if (check_ashape(root, am = (al+ah)/2)) ah = am;
        else al = am;
        if ((ah-al)/ah<.5) break;
    }
    return 1.1*ah;
}






void vols(fg *f, Tree *t, basis_s* n, int depth) {

    static simplex *s;
    static neighbor *sn;
    int tdim = cdim;
    basis_s *nn = 0;
    int signum;
    point nnv;
    double sqq;


    if (!t) {return;}

    if (!s) {NEWL(simplex,s); sn = s->neigh;}
    cdim = depth;
    s->normal = n;
    if (depth>1 && sees(t->key,s)) signum = -1; else signum = 1;
    cdim = tdim;

    if (t->fgs->dist == 0) {
        sn[depth-1].vert = t->key;
        NULLIFY(basis_s,sn[depth-1].basis);
        cdim = depth; get_basis_sede(s); cdim = tdim;
        reduce(&nn, coordsAtInfinity, s, depth);
        nnv = nn->vecs;
        if (t->key==coordsAtInfinity || f->dist==Huge || NEARZERO(nnv[rdim-1]))
            t->fgs->dist = Huge;
        else
            t->fgs->dist = Vec_dot_pdim(nnv,nnv)
                /4/nnv[rdim-1]/nnv[rdim-1];
        if (!t->fgs->facets) t->fgs->vol = 1;
        else vols(t->fgs, t->fgs->facets, nn, depth+1);
    }

    assert(f->dist!=Huge || t->fgs->dist==Huge);
    if (t->fgs->dist==Huge || t->fgs->vol==Huge) f->vol = Huge;
    else {
        sqq = t->fgs->dist - f->dist;
        if (NEARZERO(sqq)) f->vol = 0;
        else f->vol += signum
                 *sqrt(sqq)
                 *t->fgs->vol
                 /(cdim-depth+1);
    }
    vols(f, t->left, n, depth);
    vols(f, t->right, n, depth);

    return;
}


void find_volumes(fg *faces_gr, FILE *F) {
    if (!faces_gr) return;
    vols(faces_gr, faces_gr->facets, 0, 1);
    print_fg(faces_gr, F);
}


gsitef *get_site;
site_n *site_num;


void set_ch_root(simplex *s) {ch_root = s; return;}
/* set root to s, for purposes of getting normals etc. */


simplex *build_convex_hull(gsitef *get_s, site_n *site_numm, short dim, short vdd) {

    /*
      get_s     returns next site each call;
      hull construction stops when NULL returned;
      site_numm returns number of site when given site;
      dim       dimension of point set;
      vdd       if (vdd) then return Delaunay triangulation


    */

    simplex *s, *root;

    if (!Huge) Huge = DBL_MAX*DBL_MAX;

    cdim = 0;
    get_site = get_s;
    site_num = site_numm;
    pdim = dim;
    vd = vdd;

    exact_bits = DBL_MANT_DIG*log(FLT_RADIX)/log(2);
    b_err_min = DBL_EPSILON*MAXDIM*(1<<MAXDIM)*MAXDIM*3.01;
    b_err_min_sq = b_err_min * b_err_min;

    assert(get_site!=NULL); assert(site_num!=NULL);

    rdim = vd ? pdim+1 : pdim;
    if (rdim > MAXDIM)
        panic("dimension bound MAXDIM exceeded; rdim=%d; pdim=%d\n",
              rdim, pdim);
    /*  fprintf(DFILE, "rdim=%d; pdim=%d\n", rdim, pdim); fflush(DFILE);*/

    point_size = site_size = sizeof(Coord)*pdim;
    basis_vec_size = sizeof(Coord)*rdim;
    basis_s_size = sizeof(basis_s)+ (2*rdim-1)*sizeof(Coord);
    simplex_size = sizeof(simplex) + (rdim-1)*sizeof(neighbor);
    Tree_size = sizeof(Tree);
    fg_size = sizeof(fg);

    root = NULL;
    if (vd || power_diagram ) {
        p = coordsAtInfinity;
        NEWLRC(basis_s, infinity_basis);
        infinity_basis->vecs[2*rdim-1]
            = infinity_basis->vecs[rdim-1]
            = 1;
        infinity_basis->sqa
            = infinity_basis->sqb
            = 1;
    } else if (!(p = (*get_site)())) return 0;

    NEWL(simplex,root);

    ch_root = root;

    copy_simp(s,root);
    root->peak.vert = p;
    root->peak.simp = s;
    s->peak.simp = root;

    buildhull(root);
    return root;
}


void free_hull_storage(void) {
    free_basis_s_storage();
    free_simplex_storage();
    free_Tree_storage();
    free_fg_storage();
}
