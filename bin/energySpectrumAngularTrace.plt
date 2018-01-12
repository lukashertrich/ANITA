reset
set size square
set logscale x
set xlabel 'Energy [eV]'
set ylabel 'Angle below horizon [rad]'
set zlabel 'Transmittance [Fraction of initial flux]'
set pm3d interpolate 2,2	
set pm3d map

splot  'energySpectrumAngularTrace.dat' using 1:2:3 with pm3d
