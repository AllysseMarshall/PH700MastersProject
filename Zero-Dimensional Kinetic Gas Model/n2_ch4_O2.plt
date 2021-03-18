#!/usr/bin/gnuplot
set logscale x
set logscale y
set style data line
set yrange [100:]
set xrange [1e-12:]
#unset key
#set format x "10^{%T}"
#set format y "10^{%T}"
set xlabel "t (s)"
set ylabel "n (cm^-3)"
#set terminal postscript eps enhanced
#set output "konc1.eps"
plot "konc.dat" using 1:197 t "H2CO", "" u 1:202 t "CH3O", "" u 1:203 t "C2H5O", "" u 1:200 t "CH3CO",\
      "" u 1:207 t "CH3OH", "" u 1:208 t "C2H5OH", "" u 1:221 t "NO","" u 1:218 t "N2O"
quit
