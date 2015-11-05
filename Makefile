CC     = gcc
CPLUS  = g++
AR     = ar
CFLAGS = -g -Wall
_OBJS  = hull.o ch.o io.o crust.o power.o rand.o pointops.o fg.o math.o \
			   predicates.o heap.o label.o
OBJS   = $(patsubst %,src/%, $(_OBJS))
_HDRS  = hull.h points.h pointsites.h stormacs.h
HDRS   = $(patsubst %,src/%, $(_HDRS))
_SRC   = hull.c ch.c io.c crust.c power.c rand.c pointops.c fg.c math.c \
			   predicates.c heap.c label.c
SRC    = $(patsubst %,src/%, $(_SRC))
PROG   = powercrust
LIB    = lib$(PROG).a


all	: $(PROG) simplify orient

$(OBJS) : $(HDRS)

hullmain.o	: $(HDRS)

$(PROG)	: $(OBJS) src/hullmain.o
	mkdir target
	$(CC) $(CFLAGS) $(OBJS) src/hullmain.o -o target/$(PROG) -lm
	$(AR) rcv target/$(LIB) $(OBJS)

simplify: src/powershape.C src/sdefs.h
	$(CPLUS) -o target/simplify src/powershape.C -lm

orient: src/setNormals.C src/ndefs.h
	$(CPLUS) -o target/orient src/setNormals.C -lm

clean	:
	-rm -f $(OBJS) src/hullmain.o
	-rm -rf target
