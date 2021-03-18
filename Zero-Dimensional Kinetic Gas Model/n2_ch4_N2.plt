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
plot "konc.dat" using 1:86 t "N(4S)", "" u 1:88 t "N2(A)", "" u 1:90 t "N2(B)", "" u 1:91 t "N2(C)",\
      "" u 1:97 t "NH3"
quit
