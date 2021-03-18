           Instruction manual

The list of species and reactions is in file reakce_new_19.txt.
There are also the rate coefficients for each reaction.
This file is a text file, it can be edited by a text editor (e.g. vi).
First, run the Python program parsereact.py using the command

python parsereact.py reakce_new_19.txt

Then compile the Fortran codes using command

make

and run the program using command 

./a.out

The program runs several minutes and produces the file konc.dat,
which is again the text file.
The files with extension .plt plot the time dependence of different
concentrations. The files with extension .plt have be set as executable.

I use Ubuntu 16.04 Linux distribution. Programs python (version 2.7),
gfortran and gnuplot have to be instaled.



