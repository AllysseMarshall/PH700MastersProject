# This is an commentary line in a makefile
# Start of the makefile
Run: dc_decsol.o decsol.o n2_ch4.o params.o radau5.o reaction_scheme.o
	f95 -o Run dc_decsol.o decsol.o n2_ch4.o params.o radau5.o reaction_scheme.o

konc.dat: parsereact.py
	python parsereact.py reakce_new_19.txt
dc_decsol.o: dc_decsol.f
	f95 -c dc_decsol.f
decsol.o: decsol.f
	f95 -c decsol.f
n2_ch4.o: params.o reaction_scheme.o n2_ch4.f90
	f95 -c n2_ch4.f90
params.o: params.f90
	f95 -c params.f90
params.mod: params.o params.f90
	f95 -c params.f90
radau5.o: radau5.f
	f95 -c radau5.f
reaction_scheme.o: reaction_scheme.f90
	f95 -c reaction_scheme.f90
reaction_scheme.mod: reaction_scheme.o reaction_scheme.f90
	f95 -c reaction_scheme.f90
clean:
	rm dc_decsol.o decsol.o n2_ch4.o params.o radau5.o reaction_scheme.o params.mod reaction_scheme.mod Run konc.dat n2_ch4.png n2_ch4_N1.png n2_ch4_N2.png n2_ch4_O1.png n2_ch4_O2.png react_tmp
# End of the makefile
