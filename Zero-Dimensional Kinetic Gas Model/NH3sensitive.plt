#!/usr/bin/gnuplot
set logscale x
set style data line
set yrange [1e13:]
set xrange [1e-12:]
#unset key
#set format x "10^{%T}"
#set format y "10^{%T}"
set xlabel "t (s)"
set ylabel "n (cm^-3)"
#set terminal postscript eps enhanced
#set output "konc1.eps"
plot "konc.dat" using 1:96 t "NH2", "" u 1:97 t "NH3", "" u 1:238 t "NH2(s)"
quit
