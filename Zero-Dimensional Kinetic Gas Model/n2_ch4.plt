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
plot "konc.dat" using 1:86 t "N(4S)", "" u 1:106 t "HCN", "" u 1:123 t "C2N2", "" u 1:18 t "C2H2",\
      "" u 1:22 t "C2H6", "" u 1:27 t "C3H8", "" u 1:4 t "H2", "" u 1:16 t "CH4"
quit
