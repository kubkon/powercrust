#	type "make".

CC	= gcc
CPLUS  = g++
AR	= ar
CFLAGS	= -g -Wall
OBJS	= hull.o ch.o io.o crust.o power.o rand.o pointops.o fg.o math.o predicates.o heap.o label.o
HDRS	= hull.h points.h pointsites.h stormacs.h
SRC	= hull.c ch.c io.c crust.c power.c rand.c pointops.c fg.c math.c predicates.c heap.c label.c
PROG	= powercrust
LIB	= lib$(PROG).a


all	: $(PROG) simplify orient

$(OBJS) : $(HDRS)

hullmain.o	: $(HDRS)

$(PROG)	: $(OBJS) hullmain.o
	mkdir target
	$(CC) $(CFLAGS) $(OBJS) hullmain.o -o target/$(PROG) -lm
	$(AR) rcv target/$(LIB) $(OBJS)

simplify: powershape.C sdefs.h
	$(CPLUS) -o target/simplify powershape.C -lm

orient: setNormals.C ndefs.h
	$(CPLUS) -o target/orient setNormals.C -lm

clean	:
	-rm -f $(OBJS) hullmain.o
	-rm -rf target
