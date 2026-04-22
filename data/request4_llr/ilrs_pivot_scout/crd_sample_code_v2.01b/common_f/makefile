# CRD Sample code
#   common_f
#
# NOTE: It may be necessary to either replace "f77" with the name of the
# FORTRAN compiler being used, or to link that compiler to "f77".


#---------------------------------------------------------------------------

WRTSRC= wr_testf.c
WRTOBJ= wr_testf.o

wr_testf: $(WRTOBJ)
	f77 -g -fno-align-commons $(WRTOBJ) \
        read_crdf.o write_crdf.o trimlen.o -lm \
        -o wr_testf

#---------------------------------------------------------------------------
all: read_crdf.o write_crdf.o trimlen.o wr_testf

read_crdf.o: read_crdf.f
	f77 -g -c -fno-align-commons -o read_crdf.o read_crdf.f

write_crdf.o: write_crdf.f
	f77 -g -c -fno-align-commons -o write_crdf.o write_crdf.f

trimlen.o: trimlen.f
	f77 -g -c -fno-align-commons -o trimlen.o trimlen.f

clean:
	rm -f *.o
