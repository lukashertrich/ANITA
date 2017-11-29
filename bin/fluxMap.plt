reset
set pm3d depthorder hidden3d
set pm3d implicit
set pm3d interpolate 1,1
set xlabel 'phi'
set ylabel 'theta'
splot "fluxMap.dat" w pm3d