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
plot "konc.dat" using 1:173 t "O", "" u 1:175 t "O2", "" u 1:178 t "O3", "" u 1:185 t "CO",\
      "" u 1:186 t "CO2", "" u 1:189 t "H2O", "" u 1:191 t "H2O2"
quit
