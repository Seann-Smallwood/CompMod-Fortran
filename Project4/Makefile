
OBJS1 = numtype.o prob1.o 

PROG1 = run

F90 = gfortran

F90FLAGS = -O3 -funroll-loops  -fexternal-blas

LIBS = -framework Accelerate

LDFLAGS = $(LIBS)

all: $(PROG1) 

$(PROG1): $(OBJS1)
	$(F90) $(LDFLAGS) -o $@ $(OBJS1) 

clean:
	rm -f $(PROG1) *.{o,mod} fort.*

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<



