VPATH = .
sources = data.f90 input.f90 init.f90 props.f90 calcu.f90 calcv.f90 calcp.f90 calct.f90 promod.f90 main.f90   \
	  writeq.f90 lisolv.f90 force.f90 solvex.f90 geom.f90 bforce.f90 setbc.f90

OBJ := $(addsuffix .o, $(basename ${sources}))

#  IBM option #
#fc = xlf
#opt = -O -qautodbl=dblpad
#opt = -O

# gfortran option #
fc = gfortran
#opt = -fdefault-real-8 
opt  = 

teach: 	$(OBJ)
	$(fc) ${opt} -o teach ${OBJ}

${OBJ}: $(sources)
	$(fc) $(opt) -c ${sources}

clean:
ifeq ($(OS), Windows_NT)
	del /Q $(OBJ) teach
else
	rm -f $(OBJ) teach
endif
