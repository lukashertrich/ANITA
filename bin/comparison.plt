reset
set key top left
set title "Transmittance of 10^{19} eV neutrino flux at south pole"
set xlabel "Nadir [deg]"
set ylabel "Transmittance [Fraction of initial flux]"
#plot 'angularTrace.dat' using 1 : 2 with lines title 'CC' lw 2, \
#'angularTrace.dat' using 1 : 3 with lines title 'NC' lw 2, \
#'legacyTransmittance.dat' using 3 : 5 with lines title 'Legacy Earth Model [NC & CC average]' lw 2,\
#'angularTrace.dat' using 1 : 4 with lines title 'Total' lw 2,\
#'angularTrace.dat' using 1 : 5 with lines title 'Antineutrino CC' lw 2,\
#'angularTrace.dat' using 1 : 6 with lines title 'Antineutrino NC' lw 2,\
#'angularTrace.dat' using 1 : 7 with lines title 'Antineutrino Total' lw 2,\

plot 'legacyTransmittance.dat' using 3 : ($5 - $6) with lines title 'Difference of outputs' lw 2
