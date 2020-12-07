#!/usr/local/bin/gnuplot

set terminal postscript
set out "image.eps"
set multiplot
plot "data/lines.dat" with lines
plot "data/points.dat"
