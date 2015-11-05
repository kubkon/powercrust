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

#define max_blocks 10000
#define Nobj 10000

#define STORAGE_GLOBALS(X)		\
					\
extern size_t X##_size;			\
extern X *X##_list;			\
extern X *new_block_##X(int);		\
extern void flush_##X##_blocks(void);	\
void free_##X##_storage(void);		\


#define INCP(X,p,k) ((X*) ( (char*)p + (k) * X##_size)) /* portability? */


#define STORAGE(X)						\
								\
size_t	X##_size;						\
X	*X##_list = 0;						\
								\
X *new_block_##X(int make_blocks)				\
{	int i;							\
	static	X *X##_block_table[max_blocks];			\
	X *xlm, *xbt;					\
	static int num_##X##_blocks;				\
	if (make_blocks) {					\
		assert(num_##X##_blocks<max_blocks);		\
        DEB(0, before) DEBEXP(0, Nobj * X##_size)			\
								\
		xbt = X##_block_table[num_##X##_blocks++] =	(X*)malloc(Nobj * X##_size); \
 		memset(xbt,0,Nobj * X##_size);	\
		if (!xbt) {					\
			DEBEXP(-10,num_##X##_blocks)		\
		}						\
		assert(xbt);					\
								\
		xlm = INCP(X,xbt,Nobj);				\
		for (i=0;i<Nobj; i++) {				\
			xlm = INCP(X,xlm,(-1));			\
			xlm->next = X##_list;			\
			X##_list = xlm;				\
		}						\
								\
		return X##_list;				\
	}							\
								\
	for (i=0; i<num_##X##_blocks; i++)			\
		free(X##_block_table[i]);			\
	num_##X##_blocks = 0;					\
	X##_list = 0;						\
	return 0;						\
}								\
								\
void free_##X##_storage(void) {new_block_##X(0);}		\
/*end of STORAGE*/

#define NEWL(X,p)						\
{								\
 	p = X##_list ? X##_list : new_block_##X(1);		\
	assert(p);						\
 	X##_list = p->next;					\
}								\



#define NEWLRC(X,p)						\
{								\
	p = X##_list ? X##_list : new_block_##X(1);		\
	assert(p);						\
	X##_list = p->next;					\
	p->ref_count = 1;					\
}								\


#define FREEL(X,p)						\
{								\
	memset((p),0,X##_size);					\
	(p)->next = X##_list;					\
	X##_list = p;						\
}								\


#define dec_ref(X,v)	{if ((v) && --(v)->ref_count == 0) FREEL(X,(v));}
#define inc_ref(X,v)	{if (v) v->ref_count++;}
#define NULLIFY(X,v)	{dec_ref(X,v); v = NULL;}



#define mod_refs(op,s)					\
{							\
	int i;						\
	neighbor *mrsn;					\
							\
	for (i=-1,mrsn=s->neigh-1;i<cdim;i++,mrsn++)	\
		op##_ref(basis_s, mrsn->basis);		\
}

#define free_simp(s)				\
{	mod_refs(dec,s);			\
	FREEL(basis_s,s->normal);		\
	FREEL(simplex, s);			\
}						\


#define copy_simp(new,s)			\
{	NEWL(simplex,new);			\
	memcpy(new,s,simplex_size);		\
	mod_refs(inc,s);			\
}						\





#if 0
STORAGE_GLOBALS(type)
    STORAGE(type)
    NEWL(type,xxx)
    FREEL(type,xxx)
    dec_ref(type,xxxx)
    inc_ref(type,xxxx)
    NULLIFY(type,xxxx)
#endif
