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
#include <stdlib.h>
#include <stdio.h>

#include "points.h"

int	pdim;	/* point dimension */

#define NEARZERO(d)	((d) < FLT_EPSILON && (d) > -FLT_EPSILON)
Coord maxdist(int dim, point p1, point p2) {
	Coord	x,y,
		d = 0;
	int i = dim;


	while (i--) {
		x = *p1++;
		y = *p2++;
		d += (x<y) ? y-x : x-y ;
	}

	return d;
}

void print_point(FILE *F, int dim, point p) {
	int j;
	if (!p) {
		fprintf(F, "NULL");
		return;
	}
	for (j=0;j<dim;j++) fprintf(F, "%g  ", *p++);
}

void print_point_int(FILE *F, int dim, point p) {
	int j;
	if (!p) {
		fprintf(F, "NULL");
		return;
	}
	for (j=0;j<dim;j++) fprintf(F, "%.20g  ", *p++);
}


int scale(int dim, point p) {
	Coord max = 0;
	int i;
	Coord abs,val;
	for (i=0;i<dim;i++) {
		val = p[i];
		abs = (val > 0) ? val: -val;
		max = (abs > max) ? abs : max;
	}


	if (max< 100*DBL_EPSILON) {
		fprintf(stderr, "fails to scale: ");
		print_point(stderr, dim,p);fflush(stderr);
		fprintf(stderr, "\n");
		return 1;
	}

	for (i=0;i<dim;i++) p[i] /= max;

	return 0;
}


/*
  int normalize(int dim, point p) {
  Coord norm;
  int i;

  norm = norm2(dim,p);

  if (norm < 3*FLT_EPSILON) {
  fprintf(stderr, "fails to normalize: ");
  print_point(dim,p);fflush(stdout);
  fprintf(stderr, "\n");
  return 1;
  }


  for (i=0;i<dim;i++) p[i] /= norm;

  return 0;
  }

*/
