
FC = gfortran
FFLAGS = -O3 -Wall -Wextra -fcheck=all -ffree-line-length-none

TARGET = solver
OBJECTS = vardec.o sub.o exact.o main.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

vardec.o: vardec.F90
	$(FC) $(FFLAGS) -c $<

sub.o: sub.F90 vardec.o
	$(FC) $(FFLAGS) -c $<

exact.o: exact.F90 vardec.o
	$(FC) $(FFLAGS) -c $<

main.o: main.F90 vardec.o sub.o exact.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.mod $(TARGET)