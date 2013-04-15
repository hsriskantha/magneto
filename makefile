# ----------------------------------------------------------------------------------------------------------------------------------
# --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# --- The Makefile -----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
#
#     Copyright 2012, 2013 Hari Sriskantha.
#     This file is part of Magneto.
#
#     Magneto is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#     Magneto is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#   
#     A copy of the GNU General Public License can be found in the folder 'readmes', or at <http://www.gnu.org/licenses/>.
#
# ----------------------------------------------------------------------------------------------------------------------------------


# COMPILER
# --------

  COMPILER = gfortran


# DEBUGGING OPTIONS
# -----------------
#   Uncomment just one of these options...

# GDB with traps for floating point exceptions
# OPTS = -g -ffpe-trap=invalid,overflow -fbounds-check

# GDB
# OPTS = -g

# GDB with profiling
# OPTS = -g -pg

# Optimised
  OPTS = -O3



# ----------------------------------------------------------------------------------------------------------------------------------

OBJS 	 = ./objects/
FFLAGS 	 = -I$(OBJS) -J$(OBJS) $(OPTS)



MAGNETO:	$(OBJS)A.o $(OBJS)B.o  $(OBJS)C.o  $(OBJS)D.o	\
                $(OBJS)E.o $(OBJS)F.o  $(OBJS)G.o  $(OBJS)H.o	\
		magneto.f90

		$(COMPILER) $(FFLAGS) -o MAGNETO			\
                $(OBJS)A.o $(OBJS)B.o  $(OBJS)C.o  $(OBJS)D.o	\
		$(OBJS)E.o $(OBJS)F.o  $(OBJS)G.o  $(OBJS)H.o	\
		magneto.f90



$(OBJS)H.o:	H-error-mod.f90 $(OBJS)A.o
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)H.o		\
		H-error-mod.f90 

$(OBJS)G.o:	G-evolve-mod.f90 $(OBJS)A.o $(OBJS)B.o 		\
		$(OBJS)C.o $(OBJS)D.o $(OBJS)E.o $(OBJS)F.o
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)G.o		\
		G-evolve-mod.f90	

$(OBJS)F.o:	F-solve-mod.f90 $(OBJS)A.o $(OBJS)E.o
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)F.o		\
		F-solve-mod.f90	

$(OBJS)E.o:	E-recon-mod.f90 $(OBJS)A.o $(OBJS)D.o
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)E.o		\
		E-recon-mod.f90	

$(OBJS)D.o:	D-grids-mod.f90 $(OBJS)A.o $(OBJS)B.o 			
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)D.o		\
		D-grids-mod.f90	

$(OBJS)C.o:	C-tests-mod.f90 $(OBJS)A.o $(OBJS)B.o			
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)C.o		\
		C-tests-mod.f90

$(OBJS)B.o:	B-fluid-mod.f90 $(OBJS)A.o			
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)B.o		\
		B-fluid-mod.f90

$(OBJS)A.o:	A-setup-mod.f90
		$(COMPILER) $(FFLAGS) -c -o $(OBJS)A.o		\
		A-setup-mod.f90



clean:
		rm -f *~ $(OBJS)*

fullclean:
		rm -f *~ $(OBJS)* ./output/* MAGNETO
