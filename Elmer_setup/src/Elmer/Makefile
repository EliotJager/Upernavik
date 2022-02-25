LIBS=PASSIVE RAMP SlipCoef_LR Scalar_OUTPUT DISTANCE ComputeNormal SSASolver

all: $(LIBS)

% : %.F90
	elmerf90 $^ -o $@

clean:
	rm -f $(LIBS)
