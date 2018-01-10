reset
set xlabel 'Energy [eV]'
set ylabel 'Angle below horizon [rad]'
set zlabel 'Transmittance [Fraction of initial flux]'
set pm3d implicit
splot  'energySpectrumAngularTrace.dat' using 1:2:3 with pm3d
