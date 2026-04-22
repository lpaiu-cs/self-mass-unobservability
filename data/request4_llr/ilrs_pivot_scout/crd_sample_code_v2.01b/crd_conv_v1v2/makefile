# CRD Sample code
#   common_c

#---------------------------------------------------------------------------

CONSRC= crd_conv.c
CONOBJ= crd_conv.o

crd_conv: $(CONOBJ)
	gcc -g $(CONOBJ) ../common_c/read_crd.o ../common_c/write_crd.o \
		-lm -o crd_conv

#---------------------------------------------------------------------------
CRTSRC= crd_create.c
CRTOBJ= crd_create.o

crd_create: $(CRTOBJ)
	gcc -g $(CRTOBJ) ../common_c/read_crd.o ../common_c/write_crd.o \
		../common_c/grtodoy.o \
		-lm -o crd_create

#---------------------------------------------------------------------------
all:
	make crd_conv
	make crd_create


#crd_conv.o: crd_conv.c
#	gcc -g -c -o crd_conv.o crd_conv.c

clean:
	rm -f *.o crd_conv crd_create

