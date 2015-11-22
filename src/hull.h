/* hull.h */
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

/*
 * This file is a significant modification of Ken Clarkson's file hull.h
 * We include his copyright notice in accordance with its terms.
 *                                                                     - Nina, Sunghee and Ravi
 */

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

#ifndef HDR
#define HDR 1

#include "points.h"
#include "stormacs.h"

#define MAXDIM 8
#define BLOCKSIZE 100000
#define MAXBLOCKS 1000
#define DEBUG -7
#define CHECK_OVERSHOOT 1
#define EXACT 1 /* sunghee */
#define NRAND 5  /* number of random points chosen to check orientation */
#define SMALL_ENOUGH .0001
#define MAXNF 100 /* max number of conv hull triangles adjacent to a vertex */
#define MAXTA 100000
#define MAXTEMPA 100000
#define CNV 0 /* sunghee : status of simplex, if it's on convex hull */
#define VV 1 /* sunghee :    if it's regular simplex  */
#define SLV -1 /* sunghee : if orient3d=0, sliver simplex */
#define AV 2 /* if av contains the averaged pole vector */
#define PSLV -2 /* probably sliver */
#define POLE_OUTPUT 3 /* VV is pole and it's ouput */
#define SQ(a) ((a)*(a)) /* sunghee */

#define BAD_POLE -1

#define IN 2
#define OUT 1
#define INIT 0
#define NONE -1

#define FIRST 0
#define NO 1
#define DEG -1
#define NORM 2
#define VOR 3
#define VOR_NORM 4
#define SLVT 7
#define PSLVT 8
#define SURF 5
#define OPP 6

#define FIRST_EDGE 0
#define POW 1
#define NOT_POW 2
#define VISITED 3

/*RAVI */

#define VALIDEDGE 24
#define INVALIDEDGE 23
#define INEDGE 25
#define OUTEDGE 26
#define ADDAXIS 13
#define PRESENT 19
#define FIXED 20
#define REMOVED 21  /* for the thinning  stuff */

/* for priority queue */
#define LEFT(i) ((i)*2)
#define RIGHT(i) ((i)*2+1)
#define PARENT(i) ((i)/2)


extern char tmpfilenam[L_tmpnam];

extern short check_overshoot_f;

FILE* efopen(char *, char *);
void  efclose(FILE* file);



extern FILE *DFILE;

#define DEBS(qq)  {if (DEBUG>qq) {
#define EDEBS }}
#define DEBOUT DFILE
#define DEB(ll,mes)  DEBS(ll) fprintf(DEBOUT,#mes "\n");fflush(DEBOUT); EDEBS
#define DEBEXP(ll,exp) DEBS(ll) fprintf(DEBOUT,#exp "=%G\n", (double) exp); fflush(DEBOUT); EDEBS
#define DEBTR(ll) DEBS(ll) fprintf(DEBOUT, __FILE__ " line %d \n" ,__LINE__);fflush(DEBOUT); EDEBS
#define warning(lev, x)                     \
    {static int messcount;                  \
        if (++messcount<=10) {DEB(lev,x) DEBTR(lev)}    \
        if (messcount==10) DEB(lev, consider yourself warned) \
    }                           \


#define SBCHECK(s) /*                               \
                                                    {double Sb_check=0;                             \
                                                    int i;                                      \
                                                    for (i=1;i<cdim;i++) if (s->neigh[i].basis)             \
                                                    Sb_check+=s->neigh[i].basis->sqb;       \
                                                    if ((float)(Sb_check - s->Sb) !=0.0)                            \
                                                    {DEBTR DEB(bad Sb); DEBEXP(s->Sb) DEBEXP(Sb_check);print_simplex(s); exit(1);}}*/\




typedef point site;

extern site p;          /* the current site */

extern Coord coordsAtInfinity[10];  /* point at infinity for Delaunay triang */

extern int
rdim,   /* region dimension: (max) number of sites specifying region */
    cdim,   /* number of sites currently specifying region */
    site_size, /* size of malloc needed for a site */
    point_size;  /* size of malloc needed for a point */



typedef struct basis_s {
    struct basis_s *next; /* free list */
    int ref_count;  /* storage management */
    int lscale;    /* the log base 2 of total scaling of vector */
    Coord sqa, sqb; /* sums of squared norms of a part and b part */
    Coord vecs[1]; /* the actual vectors, extended by malloc'ing bigger */
} basis_s;
STORAGE_GLOBALS(basis_s)


    typedef struct neighbor {
        site vert; /* vertex of simplex */
        /*        short edgestatus[3];  FIRST_EDGE if not visited
                  NOT_POW if not dual to powercrust faces
                  POW if dual to powercrust faces */
        struct simplex *simp; /* neighbor sharing all vertices but vert */
        basis_s *basis; /* derived vectors */
    } neighbor;

typedef struct simplex {
    struct simplex *next;   /* used in free list */
    short mark;
    site vv; /* Voronoi vertex of simplex ; sunghee */
    double sqradius; /* squared radius of Voronoi ball */
    /*        site av; */ /* averaged pole */
    /*        double cond; */
    /*    float Sb; */
    short status;/* sunghee : 0(CNV) if on conv hull so vv contains normal vector;
                    1(VV) if vv points to circumcenter of simplex;
                    -1(SLV) if cond=0 so vv points to hull
                    2(AV) if av contains averaged pole */
    long poleindex; /* for 1st DT, if status==POLE_OUTPUT, contains poleindex; for 2nd, contains vertex index for powercrust output for OFF file format */
    short edgestatus[6]; /* edge status :(01)(02)(03)(12)(13)(23)
                            FIRST_EDGE if not visited
                            VISITED
                            NOT_POW if not dual to powercrust faces
                            POW if dual to powercrust faces */
    /*  short tristatus[4];   triangle status :
        FIRST if not visited
        NO   if not a triangle
        DEG  if degenerate triangle
        SURF if surface triangle
        NORM if fails normal test
        VOR  if falis voronoi edge test
        VOR_NORM if fails both test */
    /* NOTE!!! neighbors has to be the LAST field in the simplex stucture,
       since it's length gets altered by some tricky Clarkson-move.
       Also peak has to be the one before it.
       Don't try to move these babies!! */
    long visit;     /* number of last site visiting this simplex */
    basis_s* normal;    /* normal vector pointing inward */
    neighbor peak;      /* if null, remaining vertices give facet */
    neighbor neigh[1];  /* neighbors of simplex */
} simplex;
STORAGE_GLOBALS(simplex)


    /* Ravi:  for the thinning stuff */

/* represent a node in the graph */

    typedef  struct  spole { /* simple version to rep neighbors */
        long index;
        struct spole *next;
    } snode;

typedef   struct vpole{
    long index; /* index of the node */
    long pindex; /* index in the actual list of poles */
    double px;
    double py;
    double pz;
    double pr;  /* the radius of the ball centered here */
    double perpma; /* perpendicular distance from sample to medial axis */
    double pw;
    snode  *adj;
    int status;  /* for thinning */
    int label;  /* might be useful for computing the crust again */
    long substitute; /* if removed points to the substitute node */
    double estlfs; /* the estimated lfs of each ball */
} vnode ;

/* edges in the powerface */

typedef struct enode {
    long sindex;
    long dindex;
}   edge;

typedef struct fnode {
    long index1;
    long index2;
    long index3;
} face;



/* end defn for medial axis thinning */


/* structure for list of opposite poles, opplist. */
typedef struct plist {
    long pid;
    double angle;
    struct plist *next;
} plist;

/* regular triangulation edge, between pole pid to center of simp? */
typedef struct edgesimp {
    short kth;
    double angle;   /* angle between balls */
    struct simplex *simp;
    long pid;
    struct edgesimp *next;
} edgesimp;

/* additional info about poles: label for pole, pointer to list of regular
   triangulation edges, squared radius of  polar ball. adjlist is an
   array of polelabels. */
typedef struct polelabel {
    struct edgesimp *eptr;
    short bad;
    short label;
    double in; /* 12/7/99 Sunghee for priority queue */
    double out; /* 12/7/99 Sunghee for priority queue */
    int hid; /* 0 if not in the heap, otherwise heap index 1..heap_size*/
    double sqradius;
    double oppradius; /* minimum squared radius of this or any opposite ball */
    double samp_distance;
    int grafindex; /* index in thinning graph data structure */
} polelabel;

typedef struct queuenode {
    long pid;
    struct queuenode *next;
} queuenode;

typedef struct temp {
    struct simplex *simp;
    int vertptr[3];
    int novert;
    /* 0,1,2,3 : 3 vertices but ts->neigh[ti].vert are vertices of triangle */
} temp;

typedef struct tarr {
    int tempptr[50];
    int num_tempptr;
    long vert;
} tarr;

/*
typedef struct tlist {
  int tempptr;
  struct tlist *next;
} tlist;
*/

typedef struct fg_node fg;
typedef struct tree_node Tree;
struct tree_node {
    Tree *left, *right;
    site key;
    int size;   /* maintained to be the number of nodes rooted here */
    fg *fgs;
    Tree *next; /* freelist */
};

STORAGE_GLOBALS(Tree)


    typedef struct fg_node {
        Tree *facets;
        double dist, vol;   /* of Voronoi face dual to this */
        fg *next;       /* freelist */
        short mark;
        int ref_count;
    } fg_node;

STORAGE_GLOBALS(fg)


    typedef void* visit_func(simplex *, void *);
    typedef int test_func(simplex *, int, void *);
typedef void out_func(point *, int, FILE*, int);


/* Ravi thin axis */

void thinaxis();
void printaxis();
void initialize();

/* from driver, e.g., hullmain.c */

typedef site gsitef(void);

extern gsitef *get_site;

typedef long site_n(site);
extern site_n *site_num;

site get_site_offline(long); /* sunghee */

extern double bound[8][3];

void read_bounding_box(long);

extern double mult_up;

extern Coord mins[MAXDIM], maxs[MAXDIM];

typedef short zerovolf(simplex *);

extern double Huge;

extern double bound[8][3];

void read_bounding_box(long);
void construct_face(simplex *, short);

/* from segt.c or ch.c */

simplex *build_convex_hull(gsitef*, site_n*, short, short);

void free_hull_storage(void);

int sees(site, simplex *);

void get_normal(simplex *s);

int out_of_flat(simplex*, site);

void set_ch_root(simplex*);

void print_site(site, FILE*);

void print_normal(simplex*);

visit_func check_marks;

double find_alpha(simplex*);
test_func alph_test;
void* visit_outside_ashape(simplex*, visit_func);

void get_basis_sede(simplex *);

void compute_distance(simplex**,int,double*);

    /* for debugging */
int check_perps(simplex*);

void find_volumes(fg*, FILE*);

#define MAXPOINTS 10000
extern short mi[MAXPOINTS], mo[MAXPOINTS];


/* from hull.c */


void *visit_triang_gen(simplex *, visit_func, test_func);
void *visit_triang(simplex *, visit_func);
void* visit_hull(simplex *, visit_func);

neighbor *op_simp(simplex *a, simplex *b);

neighbor *op_vert(simplex *a, site b);

simplex *new_simp(void);

void buildhull(simplex *);


/* from io.c */

void panic(char *fmt, ...);

typedef void print_neighbor_f(FILE*, neighbor*);
extern print_neighbor_f
print_neighbor_full,
    print_neighbor_snum;

void check_triang(simplex*);

short is_bound(simplex *);

void check_new_triangs(simplex *);

void print_extra_facets(void);

void *print_facet(FILE*, simplex*, print_neighbor_f*);

void print_basis(FILE*, basis_s*);

void *print_simplex_f(simplex*, FILE*, print_neighbor_f*);

void *print_simplex(simplex*, void*);

void print_triang(simplex*, FILE*, print_neighbor_f*);


out_func vlist_out, ps_out, cpr_out, mp_out, off_out, vv_out;
/* sunghee : added vlist_out */
/* functions for different formats */

/* added compute axis RAVI */
visit_func facets_print, afacets_print, ridges_print,
    compute_vv, compute_pole1, compute_pole2, test_surface,
    compute_2d_power_vv, compute_3d_power_vv, compute_3d_power_edges,compute_axis;
/* to print facets, alpha facets, ridges */
/* Sunghee added compute_cc, compute_pole1, compute_pole2, test_surface */

void test_temp();

/* Nina's functions in crust.c */
short is_bound(simplex *);
int close_pole(double*,double*,double);
int antiLabel(int);
int cantLabelAnything(int);
void labelPole(int,int);
void newOpposite(int, int, double);
double computePoleAngle(simplex*, simplex*, double*);
void outputPole(FILE*, FILE*, simplex*, int, double*, int*,double);

void print_edge_dat(fg *, FILE *);


/* from pointops.c */

void print_point(FILE*, int, point);
void print_point_int(FILE*, int, point);
Coord maxdist(int,point p1, point p2);



/* from rand.c */

double double_rand(void);
void init_rand(long seed);


/* from fg.c, for face graphs */

fg *build_fg(simplex*);

void print_fg(fg*, FILE *);

void print_fg_alt(fg*, FILE *, int);

void print_hist_fg(simplex *, fg*, FILE *);

/*  void arena_check(void); */  /* from hobby's debugging malloc  */

/* from predicates.c, math.c */
void normalize(double*);
double sqdist(double*, double*);
void dir_and_dist(double*, double*, double*, double*);
double dotabac(double*, double*, double*);
double maxsqdist(double*, double*, double*, double*);
double dotabc(double*, double*, double*);
void crossabc(double*, double*, double*, double*);
void tetcircumcenter(double*, double*, double*, double*, double*,double*);
void tricircumcenter3d(double*, double*, double*, double*,double*);
void exactinit();
double orient3d(double*, double*, double*, double*);
double orient2d(double*, double*, double*);

/* heap.c */
typedef struct heap_array {
    int pid;
    double pri;
} heap_array;

void init_heap(int);
void heapify(int);
int extract_max();
int insert_heap(int , double);
void update(int , double);

/* label.c */
void opp_update(int);
void sym_update(int);
void update_pri(int,int);
int propagate();
void label_unlabeled(int);


/*power.c */
int correct_orientation(double*,double*,double*,double*,double*);

#endif
