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
plot "konc.dat" using 1:104 t "CN", "" u 1:114 t "CH3CN", "" u 1:123 t "C2N2", "" u 1:135 t "C2H5CN"
quit
