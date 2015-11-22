/* fg.c : face graph of hull, and splay trees */

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

Tree* insert(site, Tree*);

Tree* delete(site, Tree*);

void printtree(Tree*,int);

void printtree_flat(Tree*);


/* splay tree code */

/*
Fri Oct 21 21:15:01 EDT 1994
    style changes, removed Sedgewickized...
    Ken Clarkson

*/
/*
           An implementation of top-down splaying with sizes
             D. Sleator <sleator@cs.cmu.edu>, January 1994.

  This extends top-down-splay.c to maintain a size field in each node.
  This is the number of nodes in the subtree rooted there.  This makes
  it possible to efficiently compute the rank of a key.  (The rank is
  the number of nodes to the left of the given key.)  It it also
  possible to quickly find the node of a given rank.  Both of these
  operations are illustrated in the code below.  The remainder of this
  introduction is taken from top-down-splay.c.

  "Splay trees", or "self-adjusting search trees" are a simple and
  efficient data structure for storing an ordered set.  The data
  structure consists of a binary tree, with no additional fields.  It
  allows searching, insertion, deletion, deletemin, deletemax,
  splitting, joining, and many other operations, all with amortized
  logarithmic performance.  Since the trees adapt to the sequence of
  requests, their performance on real access patterns is typically even
  better.  Splay trees are described in a number of texts and papers
  [1,2,3,4].

  The code here is adapted from simple top-down splay, at the bottom of
  page 669 of [2].  It can be obtained via anonymous ftp from
  spade.pc.cs.cmu.edu in directory /usr/sleator/public.

  The chief modification here is that the splay operation works even if the
  item being splayed is not in the tree, and even if the tree root of the
  tree is NULL.  So the line:

                              t = splay(i, t);

  causes it to search for item with key i in the tree rooted at t.  If it's
  there, it is splayed to the root.  If it isn't there, then the node put
  at the root is the last one before NULL that would have been reached in a
  normal binary search for i.  (It's a neighbor of i in the tree.)  This
  allows many other operations to be easily implemented, as shown below.

  [1] "Data Structures and Their Algorithms", Lewis and Denenberg,
       Harper Collins, 1991, pp 243-251.
  [2] "Self-adjusting Binary Search Trees" Sleator and Tarjan,
       JACM Volume 32, No 3, July 1985, pp 652-686.
  [3] "Data Structure and Algorithm Analysis", Mark Weiss,
       Benjamin Cummins, 1992, pp 119-130.
  [4] "Data Structures, Algorithms, and Performance", Derick Wood,
       Addison-Wesley, 1993, pp 367-375
*/


STORAGE(Tree)


#define compare(i,j) (site_num(i)-site_num(j))
    /* This is the comparison.                                       */
    /* Returns <0 if i<j, =0 if i=j, and >0 if i>j                   */

#define node_size(x) ((x) ? ((x)->size) : 0 )
    /* This macro returns the size of a node.  Unlike "x->size",     */
    /* it works even if x=NULL.  The test could be avoided by using  */
    /* a special version of NULL which was a real node with size 0.  */

    Tree * splay (site i, Tree *t)
    /* Splay using the key i (which may or may not be in the tree.) */
    /* The starting root is t, and the tree used is defined by rat  */
    /* size fields are maintained */
{
    Tree N, *l, *r, *y;
    int comp, root_size, l_size, r_size;

    if (!t) return t;
    N.left = N.right = NULL;
    l = r = &N;
    root_size = node_size(t);
    l_size = r_size = 0;

    for (;;) {
        comp = compare(i, t->key);
        if (comp < 0) {
            if (!t->left) break;
            if (compare(i, t->left->key) < 0) {
                y = t->left;                           /* rotate right */
                t->left = y->right;
                y->right = t;
                t->size = node_size(t->left) + node_size(t->right) + 1;
                t = y;
                if (!t->left) break;
            }
            r->left = t;                               /* link right */
            r = t;
            t = t->left;
            r_size += 1+node_size(r->right);
        } else if (comp > 0) {
            if (!t->right) break;
            if (compare(i, t->right->key) > 0) {
                y = t->right;                          /* rotate left */
                t->right = y->left;
                y->left = t;
                t->size = node_size(t->left) + node_size(t->right) + 1;
                t = y;
                if (!t->right) break;
            }
            l->right = t;                              /* link left */
            l = t;
            t = t->right;
            l_size += 1+node_size(l->left);
        } else break;
    }
    l_size += node_size(t->left);  /* Now l_size and r_size are the sizes of */
    r_size += node_size(t->right); /* the left and right trees we just built.*/
    t->size = l_size + r_size + 1;

    l->right = r->left = NULL;

    /* The following two loops correct the size fields of the right path  */
    /* from the left child of the root and the right path from the left   */
    /* child of the root.                                                 */
    for (y = N.right; y != NULL; y = y->right) {
        y->size = l_size;
        l_size -= 1+node_size(y->left);
    }
    for (y = N.left; y != NULL; y = y->left) {
        y->size = r_size;
        r_size -= 1+node_size(y->right);
    }

    l->right = t->left;                                /* assemble */
    r->left = t->right;
    t->left = N.right;
    t->right = N.left;

    return t;
}

Tree * insert(site i, Tree * t) {
    /* Insert key i into the tree t, if it is not already there. */
    /* Return a pointer to the resulting tree.                   */
    Tree * new;

    if (t != NULL) {
        t = splay(i,t);
        if (compare(i, t->key)==0) {
            return t;  /* it's already there */
        }
    }
    NEWL(Tree,new)
        if (!t) {
            new->left = new->right = NULL;
        } else if (compare(i, t->key) < 0) {
            new->left = t->left;
            new->right = t;
            t->left = NULL;
            t->size = 1+node_size(t->right);
        } else {
            new->right = t->right;
            new->left = t;
            t->right = NULL;
            t->size = 1+node_size(t->left);
        }
    new->key = i;
    new->size = 1 + node_size(new->left) + node_size(new->right);
    return new;
}

Tree * delete(site i, Tree *t) {
    /* Deletes i from the tree if it's there.               */
    /* Return a pointer to the resulting tree.              */
    Tree * x;
    int tsize;

    if (!t) return NULL;
    tsize = t->size;
    t = splay(i,t);
    if (compare(i, t->key) == 0) {               /* found it */
        if (!t->left) {
            x = t->right;
        } else {
            x = splay(i, t->left);
            x->right = t->right;
        }
        FREEL(Tree,t);
        if (x) x->size = tsize-1;
        return x;
    } else {
        return t;                         /* It wasn't there */
    }
}

Tree *find_rank(int r, Tree *t) {
    /* Returns a pointer to the node in the tree with the given rank.  */
    /* Returns NULL if there is no such node.                          */
    /* Does not change the tree.  To guarantee logarithmic behavior,   */
    /* the node found here should be splayed to the root.              */
    int lsize;
    if ((r < 0) || (r >= node_size(t))) return NULL;
    for (;;) {
        lsize = node_size(t->left);
        if (r < lsize) {
            t = t->left;
        } else if (r > lsize) {
            r = r - lsize -1;
            t = t->right;
        } else {
            return t;
        }
    }
}

void printtree_flat_inner(Tree * t) {
    if (!t) return;

    printtree_flat_inner(t->right);
    printf("%f ", *(t->key));
    fflush(stdout);
    printtree_flat_inner(t->left);
}

void printtree_flat(Tree * t) {
    if (!t) {
        printf("<empty tree>");
        return;
    }
    printtree_flat_inner(t);
}


void printtree(Tree * t, int d) {
    int i;
    if (!t) return;

    printtree(t->right, d+1);
    for (i=0; i<d; i++) printf("  ");
    printf("%f(%d)\n", *(t->key), t->size);
    fflush(stdout);
    printtree(t->left, d+1);
}






fg *faces_gr_t;

STORAGE(fg)

#define snkey(x) site_num((x)->vert)

    fg *find_fg(simplex *s,int q) {

    fg *f;
    neighbor *si, *sn = s->neigh;
    Tree *t;

    if (q==0) return faces_gr_t;
    if (!faces_gr_t) NEWLRC(fg, faces_gr_t);
    f = faces_gr_t;
    for (si=sn; si<sn+cdim; si++) if (q & (1<<(si-sn))) {
        t = f->facets = insert(si->vert,f->facets);
        if (!t->fgs) NEWLRC(fg, (t->fgs))
                         f = t->fgs;
    }
    return f;
}

void *add_to_fg(simplex *s, void *dum) {

    neighbor t, *si, *sj, *sn = s->neigh;
    fg *fq;
    int q,m,Q=1<<cdim;
    /* sort neigh by site number */
    for (si=sn+2;si<sn+cdim;si++)
        for (sj=si; sj>sn+1 && snkey(sj-1) > snkey(sj); sj--)
        {t=*(sj-1); *(sj-1) = *sj; *sj = t;}

    NULLIFY(basis_s,s->normal);
    NULLIFY(basis_s,s->neigh[0].basis);

    /* insert subsets */
    for (q=1; q<Q; q++) find_fg(s,q);

    /* include all superset relations */
    for (q=1; q<Q; q++) {
        fq = find_fg(s,q);
        assert(fq);
        for (m=1,si=sn;si<sn+cdim;si++,m<<=1) if (!(q&m)) {
            fq->facets = insert(si->vert,fq->facets);
            fq->facets->fgs = find_fg(s, q|m);
        }
    }
    return NULL;
}

fg *build_fg(simplex *root) {
    faces_gr_t= 0;
    visit_hull(root, add_to_fg);
    return faces_gr_t;
}

void visit_fg_i(   void (*v_fg)(Tree *, int, int),
                   Tree *t, int depth, int vn, int boundary) {
    int boundaryc = boundary;

    if (!t) return;

    assert(t->fgs);
    if (t->fgs->mark!=vn) {
        t->fgs->mark = vn;
        if (t->key!=coordsAtInfinity && !mo[site_num(t->key)]) boundaryc = 0;
        v_fg(t,depth, boundaryc);
        visit_fg_i(v_fg, t->fgs->facets,depth+1, vn, boundaryc);
    }
    visit_fg_i(v_fg, t->left,depth,vn, boundary);
    visit_fg_i(v_fg, t->right,depth,vn,boundary);
}

void visit_fg(fg *faces_gr, void (*v_fg)(Tree *, int, int)) {
    static int fg_vn;
    visit_fg_i(v_fg, faces_gr->facets, 0, ++fg_vn, 1);
}

int visit_fg_i_far(void (*v_fg)(Tree *, int),
                   Tree *t, int depth, int vn) {
    int nb = 0;

    if (!t) return 0;

    assert(t->fgs);
    if (t->fgs->mark!=vn) {
        t->fgs->mark = vn;
        nb = (t->key==coordsAtInfinity) || mo[site_num(t->key)];
        if (!nb && !visit_fg_i_far(v_fg, t->fgs->facets,depth+1,vn))
            v_fg(t,depth);
    }
    nb = visit_fg_i_far(v_fg, t->left,depth,vn) || nb;
    nb = visit_fg_i_far(v_fg, t->right,depth,vn) || nb;
    return nb;
}

void visit_fg_far(fg *faces_gr, void (*v_fg)(Tree *, int)) {
    static int fg_vn;
    visit_fg_i_far(v_fg,faces_gr->facets, 0, --fg_vn);
}



FILE *FG_OUT;

void p_fg(Tree* t, int depth, int bad) {
    static int fa[MAXDIM];
    int i;
    static double mults[MAXDIM];

    if (mults[0]==0) {
        mults[pdim] = 1;
        for (i=pdim-1; i>=0; i--) mults[i] = mult_up*mults[i+1];
    }

    fa[depth] = site_num(t->key);
    for (i=0;i<=depth;i++)
        fprintf(FG_OUT, "%d ", fa[i]);
    fprintf(FG_OUT, "   %G\n", t->fgs->vol/mults[depth]);
}

int p_fg_x_depth;

void p_fg_x(Tree*t, int depth, int bad) {

    static int fa[MAXDIM];
    static point fp[MAXDIM];
    int i;

    fa[depth] = site_num(t->key);
    fp[depth] = t->key;

    if (depth==p_fg_x_depth) for (i=0;i<=depth;i++)
        fprintf(FG_OUT, "%d%s", fa[i], (i==depth) ? "\n" : " ");
}

void print_fg_alt(fg *faces_gr, FILE *F, int fd) {
    FG_OUT=F;
    if (!faces_gr) return;
    p_fg_x_depth = fd;
    visit_fg(faces_gr, p_fg_x);
    fclose(FG_OUT);
}


void print_fg(fg *faces_gr, FILE *F) {FG_OUT=F; visit_fg(faces_gr, p_fg);}


double fg_hist[100][100], fg_hist_bad[100][100],fg_hist_far[100][100];

void h_fg(Tree *t, int depth, int bad) {
    if (!t->fgs->facets) return;
    if (bad) {
        fg_hist_bad[depth][t->fgs->facets->size]++;
        return;
    }
    fg_hist[depth][t->fgs->facets->size]++;
}

void h_fg_far(Tree* t, int depth) {
    if (t->fgs->facets) fg_hist_far[depth][t->fgs->facets->size]++;
}


void print_hist_fg(simplex *root, fg *faces_gr, FILE *F) {
    int i,j,k;
    double tot_good[100], tot_bad[100], tot_far[100];
    for (i=0;i<20;i++) {
        tot_good[i] = tot_bad[i] = tot_far[i] = 0;
        for (j=0;j<100;j++) {
            fg_hist[i][j]= fg_hist_bad[i][j]= fg_hist_far[i][j] = 0;
        }
    }
    if (!root) return;

    find_alpha(root);

    if (!faces_gr) faces_gr = build_fg(root);

    visit_fg(faces_gr, h_fg);
    visit_fg_far(faces_gr, h_fg_far);

    for (j=0;j<100;j++) for (i=0;i<20;i++) {
        tot_good[i] += fg_hist[i][j];
        tot_bad[i] += fg_hist_bad[i][j];
        tot_far[i]  += fg_hist_far[i][j];
    }

    for (i=19;i>=0 && !tot_good[i] && !tot_bad[i]; i--);
    fprintf(F,"totals   ");
    for (k=0;k<=i;k++) {
        if (k==0) fprintf(F, "  ");
        else fprintf(F,"            ");
        fprintf(F, "%d/%d/%d",
                (int)tot_far[k], (int)tot_good[k], (int)tot_good[k] + (int)tot_bad[k]);
    }


    for (j=0;j<100;j++) {
        for (i=19; i>=0 && !fg_hist[i][j] && !fg_hist_bad[i][j]; i--);
        if (i==-1) continue;
        fprintf(F, "\n%d    ",j);fflush(F);

        for (k=0;k<=i;k++) {
            if (k==0) fprintf(F, "  ");
            else fprintf(F,"            ");
            if (fg_hist[k][j] || fg_hist_bad[k][j])
                fprintf(F,
                        "%2.1f/%2.1f/%2.1f",
                        tot_far[k] ? 100*fg_hist_far[k][j]/tot_far[k]+.05 : 0,
                        tot_good[k] ? 100*fg_hist[k][j]/tot_good[k]+.05 : 0,
                        100*(fg_hist[k][j]+fg_hist_bad[k][j])/(tot_good[k]+tot_bad[k])+.05
                    );
        }
    }
    fprintf(F, "\n");
}
