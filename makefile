# ----------------------------------------------------------------------------------------------------------------------------------
# --- Quarkwood Magneto 8.1 --------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# --- The Makefile -----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------



# COMPILER OPTIONS
# ----------------

# GDB with traps for floating point exceptions
# OPTS = -g -ffpe-trap=invalid,overflow -fbounds-check

# GDB
# OPTS = -g

# GDB with profiling
# OPTS = -g -pg

# Optimised
  OPTS = -O3



# ----------------------------------------------------------------------------------------------------------------------------------

OBJS 	= ./objects/
FFLAGS 	= -I$(OBJS) -J$(OBJS) $(OPTS)



MAGNETO:	$(OBJS)A.o  $(OBJS)B.o  $(OBJS)C.o  $(OBJS)D.o	\
                $(OBJS)E.o  $(OBJS)F.o  $(OBJS)G.o  $(OBJS)H.o	\
		magneto.f90

		gfortran $(FFLAGS) -o MAGNETO			\
                $(OBJS)A.o  $(OBJS)B.o  $(OBJS)C.o  $(OBJS)D.o 	\
		$(OBJS)E.o  $(OBJS)F.o  $(OBJS)G.o  $(OBJS)H.o	\
		magneto.f90



$(OBJS)H.o:	H-error-mod.f90 $(OBJS)A.o
		gfortran $(FFLAGS) -c -o $(OBJS)H.o		\
		H-error-mod.f90 

$(OBJS)G.o:	G-evolve-mod.f90 $(OBJS)A.o $(OBJS)B.o 		\
		$(OBJS)C.o $(OBJS)D.o $(OBJS)E.o $(OBJS)F.o
		gfortran $(FFLAGS) -c -o $(OBJS)G.o		\
		G-evolve-mod.f90	

$(OBJS)F.o:	F-solve-mod.f90 $(OBJS)A.o $(OBJS)E.o
		gfortran $(FFLAGS) -c -o $(OBJS)F.o		\
		F-solve-mod.f90	

$(OBJS)E.o:	E-recon-mod.f90 $(OBJS)A.o $(OBJS)D.o
		gfortran $(FFLAGS) -c -o $(OBJS)E.o		\
		E-recon-mod.f90	

$(OBJS)D.o:	D-grids-mod.f90 $(OBJS)A.o $(OBJS)B.o 			
		gfortran $(FFLAGS) -c -o $(OBJS)D.o		\
		D-grids-mod.f90	

$(OBJS)C.o:	C-tests-mod.f90 $(OBJS)A.o $(OBJS)B.o			
		gfortran $(FFLAGS) -c -o $(OBJS)C.o		\
		C-tests-mod.f90

$(OBJS)B.o:	B-fluid-mod.f90 $(OBJS)A.o			
		gfortran $(FFLAGS) -c -o $(OBJS)B.o		\
		B-fluid-mod.f90

$(OBJS)A.o:	A-setup-mod.f90
		gfortran $(FFLAGS) -c -o $(OBJS)A.o		\
		A-setup-mod.f90



clean:
		rm -f *~

fullclean:
		rm -f *~ $(OBJS)* ./output/* MAGNETO
