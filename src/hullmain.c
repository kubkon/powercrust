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
 * This file is a significant modification of Ken Clarkson's file hullmain.c.
 * We include his copyright notice in accordance with its terms.
 *                                                   - Nina, Sunghee and Ravi
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

#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>

#define POINTSITES 1

#include "hull.h"

double bound[8][3], omaxs[3], omins[3];  /* 8 vertices for bounding box */
point site_blocks[MAXBLOCKS];
int num_blocks;
extern int numfaces;
struct queuenode *queue;
struct queuenode *qend;
int num_vtxs = 0;
int num_faces = 0;

/* Data structures for poles */
struct simplex **pole1, **pole2;  /* arrays of poles - per sample*/
struct polelabel *adjlist;  /* array of additional info - per pole */
struct plist **opplist; /* opposite pid and angle between poles - per pole */
double* lfs_lb;  /*  array of lower bounds for lfs of each sample */
double est_r = 0.6;   /* estimated value of r - user input */


int num_poles = 0, num_axedgs = 0, num_axfaces = 0;
double *pole1_distance, *pole2_distance;

/* for priority queue */
extern int heap_size;

/* int  getopt(int, char**, char*); */
extern char *optarg;
extern int optind;
extern int opterr;
extern int scount;
extern int v1[6], v2[6], v3[6], v4[6];
long num_sites;
static short vd = 1;

short power_diagram = 0; /* 1 if power diagram */

static int dim;
static long s_num = 0; /* site number */

double theta = 0.0; /* input argument - angle defining deep intersection */
double deep = 0.0; /* input argument.. same as theta for labeling unlabled pole */
int defer = 0; /* input argument -D 1 if you don't want to propagate bad poles */

int poleInput=0; /* are the poles given as input */

FILE *INFILE, *OUTFILE, *DFILE, *TFILE, *SPFILE, *POLE, *PC, *PNF, *INPOLE,
     *INPBALL, *INVBALL, *AXIS, *AXISFACE;

int *rverts;

int* select_random_points(int Nv) /* for orientation testing */
{ /* Nv : Number of vertices (sites) */
  int i, j;
  int *rverts;

  rverts = (int*) malloc (NRAND * sizeof(int));

  srandom(Nv);  /* seed the random number generator */
  for (i = 0; i < NRAND; i++)
  {
    j = random() % Nv;
    rverts[i] = j;
  }
  return rverts;
}

long site_numm(site p)
{
  long i,j;

  if (( vd || power_diagram) && p==coordsAtInfinity)
  {
    return -1;
  }

  if (!p)
  {
    return -2;
  }

  for (i = 0; i < num_blocks; i++)
  {
    if ((j = p - site_blocks[i]) >= 0 && j < BLOCKSIZE * dim)
    {
      return j / dim + BLOCKSIZE * i;
    }
  }

  return -3;
}


site new_site(site p, long j)
{
  assert(num_blocks + 1 < MAXBLOCKS);
  if (0 == (j % BLOCKSIZE))
  {
    assert(num_blocks < MAXBLOCKS);
    return site_blocks[num_blocks++]=(site)malloc(BLOCKSIZE*site_size);
  }
  else
  {
    return p + dim;
  }
}

void read_bounding_box(long j)
{
  int i, k;
  double center[3], width;

  omaxs[0] = maxs[0];
  omins[0] = mins[0];
  omaxs[1] = maxs[1];
  omins[1] = mins[1];
  omaxs[2] = maxs[2];
  omins[2] = mins[2];

  center[0] = (maxs[0] - mins[0])/2;
  center[1] = (maxs[1] - mins[1])/2;
  center[2] = (maxs[2] - mins[2])/2;

  if ((maxs[0] - mins[0]) > (maxs[1] - mins[1]))
  {
    if ((maxs[2] - mins[2]) > (maxs[0] - mins[0]))
    {
      width = maxs[2] - mins[2];
    }
    else
    {
      width = maxs[0] - mins[0];
    }
  }
  else
  {
    if ((maxs[1] - mins[1]) > (maxs[2] - mins[2]))
    {
      width = maxs[1] - mins[1];
    }
    else
    {
      width = maxs[2] - mins[2];
    }
  }

  width = width * 4;

  bound[0][0] = center[0] + width;
  bound[1][0] = bound[0][0];
  bound[2][0] = bound[0][0];
  bound[3][0] = bound[0][0];
  bound[0][1] = center[1] + width;
  bound[1][1] = bound[0][1];
  bound[4][1] = bound[0][1];
  bound[5][1] = bound[0][1];
  bound[0][2] = center[2] + width;
  bound[2][2] = bound[0][2];
  bound[4][2] = bound[0][2];
  bound[6][2] = bound[0][2];
  bound[4][0] = center[0] - width;
  bound[5][0] = bound[4][0];
  bound[6][0] = bound[4][0];
  bound[7][0] = bound[4][0];
  bound[2][1] = center[1] - width;
  bound[3][1] = bound[2][1];
  bound[6][1] = bound[2][1];
  bound[7][1] = bound[2][1];
  bound[1][2] = center[2] - width;
  bound[3][2] = bound[1][2];
  bound[5][2] = bound[1][2];
  bound[7][2] = bound[1][2];

  for (i = 0; i < 8; i++)
  {
    fprintf(DFILE,
            "%f %f %f\n",
            bound[i][0]/mult_up,
            bound[i][1]/mult_up,
            bound[i][2]/mult_up);
  }

  for (k = 0; k < 3; k++)
  {
    p[k] = bound[0][k];
  }

  for (i=1;i<8;i++)
  {
    p = new_site(p, j+i);
    for (k = 0; k < 3; k++)
    {
      p[k] = bound[i][k];
    }
  }
  maxs[0] = bound[0][0];
  mins[0] = bound[4][0];
  maxs[1] = bound[0][1];
  mins[1] = bound[2][1];
  maxs[2] = bound[0][2];
  mins[2] = bound[1][2];
}

site read_next_site(long j)
{
  int i=0, k=0;
  static char buf[100], *s;

  if (j!=-1)
  {
    p = new_site(p,j);
  }

  if (j!=0)
  {
    while ((s = fgets(buf, sizeof(buf), INFILE)))
    {
      if (buf[0] == '%')
      {
        continue;
      }
      for (k=0; buf[k] && isspace(buf[k]); k++) {};
      if (buf[k])
      {
        break;
      }
    }
  }

  if (!s)
  {
    return 0;
  }

  if (j!=0)
  {
    assert(TFILE != NULL);
    fprintf(TFILE, "%s", &(buf[k]));
    fflush(TFILE);
  }
  while (buf[k])
  {
    while (buf[k] && isspace(buf[k]))
    {
      k++;
    }
    if (buf[k] && j!=-1)
    {
      if (sscanf(buf+k,"%lf",p+i)==EOF)
      {
        fprintf(DFILE, "bad input line: %s\n", buf);
        exit(1);
      }
      p[i] = floor(mult_up * p[i] + 0.5);
      mins[i] = (mins[i] < p[i]) ? mins[i] : p[i];
      maxs[i] = (maxs[i] > p[i]) ? maxs[i] : p[i];
    }
    if (buf[k])
    {
      i++;
    }
    while (buf[k] && !isspace(buf[k]))
    {
      k++;
    }
  }

  if (!dim)
  {
    dim = i;
  }
  if (i!=dim)
  {
    DEB(-10,inconsistent input);
    DEBTR(-10);
    exit(1);
  }

  return p;
}

/* reads a site from storage we're managing outselves */
site get_site_offline(long i)
{
  if (i>=num_sites)
  {
    return NULL;
  }
  else {
    return site_blocks[i/BLOCKSIZE]+(i%BLOCKSIZE)*dim;
  }
}

long *shufmat;

void make_shuffle(void)
{
  long i,t,j;
  static long mat_size = 0;

  if (mat_size<=num_sites)
  {
    mat_size = num_sites+1;
    shufmat = (long*)malloc(mat_size*sizeof(long));
  }
  for (i = 0; i <= num_sites; i++)
  {
    shufmat[i] = i;
  }
  for (i = 0; i < num_sites; i++)
  {
    t = shufmat[i];
    j = i + (num_sites-i)*double_rand();
    shufmat[i] = shufmat[j];
    shufmat[j] = t;
  }
}

static long (*shuf)(long);

long noshuffle(long i)
{
  return i;
}

long shufflef(long i)
{
  return shufmat[i];
}

static site (*get_site_n)(long);

/* returns shuffled, offline sites or reads an unshuffled site, depending on
   how get_site_n and shuf are set up. */
site get_next_site(void)
{
  /*  static long s_num = 0; */
  return (*get_site_n)((*shuf)(s_num++));
}


void errline(char *s)
{
  fprintf(stderr, s);
  fprintf(stderr,"\n");
  return;
}

void tell_options(void)
{
  errline("options:");
  errline( "-m mult  multiply by mult before rounding;");
  errline( "-s seed  shuffle with srand(seed);");
  errline( "-i<name> read input from <name>;");
  errline( "-X<name> chatter to <name>;");
  errline( "-oF<name>  prefix of output files is <name>;");
  errline( "-t min cosine of allowed dihedral angle btwn polar balls");
  errline( "-w same as -t, but for trying to label unlabled poles, the second time around.");
  errline( "-D no propagation for 1st pole of non-manifold cells");
  errline( "-B throw away both poles for non-manifold cells");
  errline( "-R guess for value of r, used to eliminate bad second poles");
}

void echo_command_line(FILE *F, int argc, char **argv)
{
  fprintf(F,"%%");
  while (--argc >= 0)
  {
    fprintf(F, "%s%s", *argv++, (argc>0) ? " " : "");
  }
  fprintf(F,"\n");
}

char *output_forms[] = {"vn", "ps", "mp", "cpr", "off"};

out_func *out_funcs[] = {&vlist_out, &ps_out, &mp_out, &cpr_out, &off_out};

int set_out_func(char *s)
{
  int i;
  for (i=0;i< sizeof(out_funcs)/(sizeof (out_func*)); i++)
  {
    if (strcmp(s,output_forms[i])==0)
    {
      return i;
    }
  }
  tell_options();
  return 0;
}

void make_output(simplex *root,
                 void *(*visit_gen)(simplex*, visit_func* visit),
                 visit_func* visit,
                 out_func* out_funcp,
                 FILE *F)
{
    out_funcp(0,0,F,-1);
    visit(0, out_funcp);
    visit_gen(root, visit);
    out_funcp(0,0,F,1);
    /*  efclose(F); */
}

int main(int argc, char **argv) {
  long seed = 0, poleid = 0;
  short shuffle = 1,
        output = 1,
        hist = 0,
        vol = 0,
        ofn = 0,
        ifn = 0,
        bad = 0; /* for -B */
  int option, num_poles = 0;
  double pole_angle;
  char ofile[50] = "",
       ifile[50] = "",
       ofilepre[50] = "";
  FILE *INPOLE, *OUTPOLE, *HEAD,*POLEINFO;
  int main_out_form=0, i,k;

  simplex *root;

  struct edgesimp *eindex;
  double samp[3];
  double tmp_pt[3];
  int numbadpoles=0;
  double x,y,z,r,d;
  int l;
  out_func *mof;
  visit_func *pr;

  /* some default values */
  mult_up = 100000;
  est_r = 1;
  DFILE = stderr;

  while ((option = getopt(argc, argv, "i:m:rs:DBo:X::f:t:w:R:p")) != EOF) {
    switch (option)
    {
      case 'm' :
        sscanf(optarg,"%lf",&mult_up);
        DEBEXP(-4,mult_up);
        break;
      case 's':
        seed = atol(optarg);
        shuffle = 1;
        break;
      case 'D':
        defer = 1;
        break;
      case 'B':
        bad = 1;
        break;
      case 'i' :
        strcpy(ifile, optarg);
        break;
      case 'X' :
        DFILE = efopen(optarg, "w");
        break;
      case 'f' :
        main_out_form = set_out_func(optarg);
        break;
      case 'o':
        switch (optarg[0])
        {
          case 'o':
            output=1;
            break;
          case 'N':
            output=1;
            break; /* output is never set to zero */
          case 'v':
            vd = vol = 1;
            break;
          case 'h':
            hist = 1;
            break;
          case 'F':
            strcpy(ofile, optarg+1);
            break;
          default:
            errline("illegal output option");
            exit(1);
        }
        break;
      case 't':
        sscanf(optarg,"%lf",&theta);
        break;
      case 'w':
        sscanf(optarg,"%lf",&deep);
        break;
      case 'R':
        sscanf(optarg,"%lf",&est_r);
        break;
      case 'p':
        poleInput=1;
        break;
      default :
        tell_options();
        exit(1);
    }
  }
  AXIS=fopen("axis","w");
  AXISFACE=fopen("axisface","w");
  POLE=fopen("pole","w");
  TFILE = efopen(tmpnam(tmpfilenam), "w");

  if (!poleInput)
  {
    ifn = (strlen(ifile)!=0);
    INFILE = ifn ? efopen(ifile, "r") : stdin;
    fprintf(DFILE, "reading from %s\n", ifn ? ifile : "stdin");

    ofn = (strlen(ofile)!=0);

    strcpy(ofilepre, ofn ? ofile : (ifn ? ifile : "hout") );

    if (output)
    {
      if (ofn && main_out_form > 0)
      {
        strcat(ofile, ".");
        strcat(ofile, output_forms[main_out_form]);
      }
      OUTFILE = ofn ? efopen(ofile, "w") : stdout;
      fprintf(DFILE, "main output to %s\n", ofn ? ofile : "stdout");
    }
    else
    {
      fprintf(DFILE, "no main output\n");
    }

    read_next_site(-1);
    fprintf(DFILE,"dim=%d\n",dim);
    fflush(DFILE);
    if (dim > MAXDIM)
    {
      panic("dimension bound MAXDIM exceeded");
    }

    point_size = site_size = sizeof(Coord)*dim;

    fprintf(DFILE, "reading sites...");
    for (num_sites=0; read_next_site(num_sites); num_sites++);
    fprintf(DFILE,"done; num_sites=%ld\n", num_sites);
    fflush(DFILE);
    read_bounding_box(num_sites);
    num_sites += 8;
    fprintf(DFILE,"shuffling...");
    init_rand(seed);
    make_shuffle();
    shuf = &shufflef;
    get_site_n = get_site_offline;

    /* Step 1: compute DT of input point set */
    root = build_convex_hull(get_next_site, site_numm, dim, vd);

    /* Step 2: Find poles */
    pole1 = (struct simplex **) calloc(num_sites, sizeof(struct simplex *));
    pole2 = (struct simplex **) calloc(num_sites, sizeof(struct simplex *));
    lfs_lb = (double*) calloc(num_sites, sizeof(double));

    fprintf(DFILE, "done\n");
    fflush(DFILE);

    mof = out_funcs[main_out_form];
    pr = facets_print;

    if (main_out_form==0)
    {
      echo_command_line(OUTFILE,argc,argv);
    }

    exactinit(); /* Shewchuk's exact arithmetic initialization */

    pr = compute_vv;
    fprintf(DFILE, "Computing Voronoi vertices and 1st poles....\n");
    make_output(root, visit_hull, pr, mof, OUTFILE);

    pr = compute_pole2;
    fprintf(DFILE, "\n\n\ncomputing 2nd poles....\n");
    make_output(root, visit_hull, pr, mof, OUTFILE);

    /* poles with weights. Input to regular triangulation */
    SPFILE = fopen("sp","w");

    /* initialize the sample distance info for the poles */
    pole1_distance=(double *) malloc(num_sites*sizeof(double));
    pole2_distance=(double *) malloc(num_sites*sizeof(double));

    compute_distance(pole1,num_sites-8,pole1_distance);
    compute_distance(pole2,num_sites-8,pole2_distance);

    /* intialize list of lists of pointers to opposite poles */
    opplist = (struct plist**) calloc(num_sites*2, sizeof(struct plist *));

    /* data about poles; adjacencies, labels, radii */
    adjlist = (struct polelabel *) calloc(num_sites*2, sizeof(struct polelabel));

    /* loop through sites, writing out poles */
    for (i=0;i<num_sites-8;i++)
    {
      /* rescale the sample to real input coordinates */
      for (k=0; k<3; k++)
      {
        samp[k] = get_site_offline(i)[k]/mult_up;
      }

      /* output poles, both to debugging file and for weighted DT */
      /* remembers sqaured radius */
      if ((pole1[i]!=NULL)&&(pole1[i]->status != POLE_OUTPUT))
      {
        /* if second pole is closer than we think it should be... */
        if ((pole2[i]!=NULL) &&
             bad &&
             close_pole(samp,pole2[i]->vv,lfs_lb[i]))
        {
          numbadpoles++;
        }
        else
        {
          outputPole(POLE,SPFILE,pole1[i],poleid++,samp,&num_poles,pole1_distance[i]);
        }
      }

      if ( (pole2[i]!=NULL) && (pole2[i]->status != POLE_OUTPUT))
      {
        /* if pole is closer than we think it should be... */
        if (close_pole(samp,pole2[i]->vv,lfs_lb[i]))
        {
          /* remember opposite bad for late labeling */
          if (!bad)
          {
            adjlist[pole1[i]->poleindex].bad = BAD_POLE;
          }
          numbadpoles++;
          continue;
        }
        /* otherwise... */
        outputPole(POLE,SPFILE,pole2[i],poleid++,samp,&num_poles,pole2_distance[i]);
      }

      /* keep list of opposite poles for later coloring */
      if ((pole1[i]!=NULL)&&(pole2[i]!=NULL)&&
          (pole1[i]->status == POLE_OUTPUT) &&
          (pole2[i]->status == POLE_OUTPUT))
      {
        pole_angle = computePoleAngle(pole1[i],pole2[i],samp);
        newOpposite(pole1[i]->poleindex, pole2[i]->poleindex, pole_angle);
        newOpposite(pole2[i]->poleindex, pole1[i]->poleindex, pole_angle);
      }
    }
    efclose(POLE);
    efclose(SPFILE);
    fprintf(DFILE,"bad poles=%d\n",numbadpoles);

    free_hull_storage();
  } /* do this if the input was a set of samples not poles */
  else
  {
    /* read data from the input file and put it in sp */
    /* initialize adjlist to the right size and initialize labels */
    INFILE=fopen(ifile,"r");
    SPFILE=fopen("sp","w");

    fprintf(DFILE,"%s",ifile);
    mof = out_funcs[main_out_form];
    pr = facets_print;
    OUTFILE=stdout;

    fscanf(INFILE,"%d",&num_poles);/* get the number of poles ..*/
    adjlist = (struct polelabel *) calloc(num_poles, sizeof(struct polelabel));
    for (i = 0; i < num_poles; i++)
    {
      fscanf(INFILE,"%le ",&x);
      fscanf(INFILE,"%le ",&y);
      fscanf(INFILE,"%le ",&z);
      fscanf(INFILE,"%le ",&r);
      fscanf(INFILE,"%d ",&l);
      fscanf(INFILE,"%le ",&d);

      /* we  have the square of the radius */
      fprintf(SPFILE,"%f %f %f %f\n",x,y,z,SQ(x)+SQ(y)+SQ(z)-r);
      fprintf(POLE,"%f %f %f \n",x,y,z);
      /*  fprintf(DFILE,"%f %f %f %d\n",x,y,z,num_poles);*/
      adjlist[i].sqradius=r;
      adjlist[i].label=l;
    }

    efclose(INFILE);
    efclose(SPFILE);
    efclose(POLE);
  }

  power_diagram = 1;
  vd = 0;
  dim = 4;

  INFILE = fopen("sp","r");
  fprintf(DFILE,"num_blocks = %d\n",num_blocks);

  num_blocks = 0;
  s_num = 0;
  scount = 0;
  read_next_site(-1);
  fprintf(DFILE,"dim=%d\n",dim);
  fflush(DFILE);

  point_size = site_size = sizeof(Coord)*dim;
  /* save points in order read */
  for (num_sites=0; read_next_site(num_sites); num_sites++);

  fprintf(DFILE,"done; num_sites=%ld\n", num_sites);fflush(DFILE);
  efclose(INFILE);

  /* set up the shuffle */
  fprintf(DFILE,"shuffling...");
  init_rand(seed);
  make_shuffle();
  shuf = &shufflef;
  get_site_n = get_site_offline;  /* returns stored points, unshuffled */

  /* Compute weighted DT  */
  root = build_convex_hull(get_next_site, site_numm, dim, vd);

  fprintf(DFILE,"scount=%d, s_num=%ld\n",scount,s_num);
  fprintf(DFILE, "done\n"); fflush(DFILE);

  /* file of faces */
  PNF = fopen("pnf","w");
  /* file of points */
  PC = fopen("pc","w");

  /* compute adjacencies and find angles of ball intersections */
  queue = NULL;
  pr = compute_3d_power_vv;
  make_output(root, visit_hull, pr, mof, OUTFILE);

  /* Begin by labeling everything outside a big bounding box as outside */
  /* labeling */
  if(!poleInput)
  { /* if we dont have the labels */
    fprintf(DFILE,"num_poles=%d\n", num_poles);
    init_heap(num_poles);
    for (i = 0; i < num_poles; i++)
    {
      if ((get_site_offline(i)[0]>(2*omaxs[0]-omins[0]))||
          (get_site_offline(i)[0]<(2*omins[0]-omaxs[0]))||
          (get_site_offline(i)[1]>(2*omaxs[1]-omins[1]))||
          (get_site_offline(i)[1]<(2*omins[1]-omaxs[1]))||
          (get_site_offline(i)[2]>(2*omaxs[2]-omins[2]))||
          (get_site_offline(i)[2]<(2*omins[2]-omaxs[2])))
      {
        adjlist[i].hid = insert_heap(i,1.0);
        adjlist[i].out = 1.0;
        adjlist[i].label = OUT;
      }
    }

    while (heap_size != 0)
    {
      propagate();
    }
    label_unlabeled(num_poles);
  }

  /* Enough labeling; let's look at the poles and output a crust!  */
  INPOLE = fopen("inpole","w");
  OUTPOLE = fopen("outpole","w");

  /* for visualization of polar balls: */
  INPBALL = fopen("inpball","w");  /* inner poles with radii */
  POLEINFO = fopen("tpoleinfo","w");

  for (i=0;i<num_poles;i++)
  {
    for (k=0; k<3; k++)
    {
      tmp_pt[k] = get_site_offline(i)[k]/mult_up;
    }

    fprintf(POLEINFO,"%f %f %f %f %d %f \n ",tmp_pt[0],tmp_pt[1],tmp_pt[2],
            adjlist[i].sqradius,adjlist[i].label,adjlist[i].samp_distance);

    if ((adjlist[i].label != IN) && (adjlist[i].label != OUT))
    {
        fprintf(DFILE,"pole %d label %d\n",i,adjlist[i].label);
    }
    else
    {
      if (adjlist[i].label == IN)
      {
          fprintf(INPOLE,"%f %f %f\n",tmp_pt[0],tmp_pt[1],tmp_pt[2]);
          fprintf(INPBALL,"%f %f %f %f \n",
                  tmp_pt[0],tmp_pt[1],tmp_pt[2],sqrt(adjlist[i].sqradius));
      }
      else if (adjlist[i].label == OUT)
      {
        fprintf(OUTPOLE,"%f %f %f\n",tmp_pt[0],tmp_pt[1],tmp_pt[2]);
      }

      eindex = adjlist[i].eptr;
      while (eindex!=NULL)
      {
        if ((i < eindex->pid) &&
          (antiLabel(adjlist[i].label) == adjlist[eindex->pid].label))
        {
          construct_face(eindex->simp,eindex->kth);
        }
        eindex = eindex->next;
      }
    }
  }

  efclose(PC);
  efclose(PNF);
  efclose(POLEINFO);

  /* powercrust output done... */
  HEAD = fopen("head","w");
  fprintf(HEAD,"OFF\n");
  fprintf(HEAD,"%d %d %d\n",num_vtxs,num_faces,0);
  efclose(HEAD);
  system("cat head pc pnf > pc.off");
  system("rm head pc pnf");

  /* compute the medial axis */
  pr=compute_axis;
  fprintf(DFILE,"\n\n computing the medial axis ....\n");
  make_output(root,visit_hull,pr,mof,OUTFILE);

  HEAD = fopen("head","w");
  fprintf(HEAD,"OFF\n");
  fprintf(HEAD,"%d %d %d\n",num_poles,num_axedgs,0);
  efclose(HEAD);
  efclose(AXIS);
  system("cat head pole axis > axis.off");

  HEAD=fopen("head","w");
  fprintf(HEAD,"%d %d \n", num_poles,num_axedgs);
  efclose(HEAD);

  system("cat head tpoleinfo axis > poleinfo");

  HEAD = fopen("head","w");
  fprintf(HEAD,"OFF\n");
  fprintf(HEAD,"%d %d %d\n",num_poles,num_axfaces,0);
  efclose(HEAD);
  efclose(AXISFACE);
  system("cat head pole axisface > axisface.off");
  system("rm -f head pole axis axisface tpoleinfo sp");
  /* power shape output done */

  efclose(INPOLE);
  efclose(OUTPOLE);
  efclose(INPBALL);
  efclose(TFILE);
  free(adjlist);

  free_hull_storage();
  efclose(DFILE);
  exit(0);
}

/* for each pole array, compute the maximum of the distances on the sample */
void compute_distance(simplex** poles, int size, double* distance)
{
  int i,j,k,l;
  double indices[4][3]; /* the coords of the four vertices of the simplex*/
  point v[MAXDIM];
  simplex* currSimplex;

  double maxdistance=0;
  double currdistance;

  for(l=0;l<size;l++)
  {  /* for each pole do*/
    if(poles[l]!=NULL)
    {
      currSimplex=poles[l];

      /* get the coordinates of the  four endpoints */
      for(j=0;j<4;j++)
      {
        v[j]=currSimplex->neigh[j].vert;
        for(k=0;k<3;k++)
        {
          indices[j][k]=v[j][k]/mult_up;
        }
      }

      /* now compute the actual distance  */
      maxdistance=0;

      for(i=0;i<4;i++)
      {
        for(j=i+1;j<4;j++)
        {
          currdistance= SQ(indices[i][0]-indices[j][0]) +
              SQ(indices[i][1]-indices[j][1])+ SQ(indices[i][2]-indices[j][2]);
          currdistance=sqrt(currdistance);
          if(maxdistance<currdistance)
          {
            maxdistance=currdistance;
          }
        }
      }
      distance[l]=maxdistance;
    }
  }
}
