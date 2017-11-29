reset
set title "Transmittance of 10^{19} eV neutrino flux at south pole"
set xlabel "Angle below horizon [rad]"
set ylabel "Transmittance [Fraction of initial flux]"
plot 'angularTrace.dat' using 1 : 2 with lines title 'CC', \
'angularTrace.dat' using 1 : 3 with lines title 'NC', \
'angularTrace.dat' using 1 : 4 with lines title 'Total',\
'angularTrace.dat' using 1 : 5 with lines title 'Antineutrino CC',\
'angularTrace.dat' using 1 : 6 with lines title 'Antineutrino NC',\
'angularTrace.dat' using 1 : 7 with lines title 'Antineutrino Total',\
