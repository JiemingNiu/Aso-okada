FC		= mpiifort

FORFILES	= nlokada_main.f nlokada_beta_phi_theta.f nlokada_phi_theta.f DC3Dfortran.f

OBJECTS		= $(FORFILES:%.f=%.o)

FFLAGS		= -O3 -mkl -static-intel

PROGRAM 	= mpi_nlokada

$(PROGRAM): 	$(OBJECTS)
		$(FC) $(FFLAGS) $(OBJECTS) -o $@

clean:
	rm *.o 
