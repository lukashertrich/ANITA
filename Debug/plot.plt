reset
set size square
set title "Density Traversed - Hemispherical Sweep"
set xlabel "x [pi/2]" rotate parallel
set ylabel "y [pi/2]" rotate parallel
set zlabel "Density traversed [kgm/m^3]" rotate parallel offset 10,10,0
set pm3d depthorder hidden3d
set pm3d implicit
set pm3d noborder
splot "diagnostic.dat" with lines