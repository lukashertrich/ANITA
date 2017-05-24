
/*
 * gcc anita_earthmodel2.c -o ~/bin/anita_earthmodel2 -lm -O2
 * 
anita_earthmodel: monte carlo neutrino interactions coming into
antarctica, anita-style.
This code determines the distribution of neutrinos of a given
energy that survive to produce a detectable event near the
surface.
P. Gorham March 2003, from earlier code dating back to 1999.
Based on code for GLUE.

Modified September 2004 to limit tracks that are too downgoing or
upcoming, until anitamcE4.c can handle the upcoming...

fixed a stupid bug on earth density, added firn density 2/17/2005 -PG

more fixes to earth density model  6/5/2005  --PG

10/11/2005, added better random seed generator --PG

11/26/2005:  checked some of the geometric formulas, fixed an
error in integrate_Lint() which made no difference 
( formula was correct by symmetry, though formally incorrect)  --PG

7/31/2007, minor changes to prelim. earth model, interaction length
   integration  --PG

11/15/2010:  updated the cross sections --PG

10/15/13 bug in firn density call fixed

*/
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h> 
#include <unistd.h>
#include <signal.h>

#define PI    M_PI

#define DRef 600.e3	/*  reference distance to "horizon"*/
#define Remean   6371.0e3     /* earth mean radius m */
#define Re   6356.755e3   /* polar earth radius m */
#define RHOearth  5.515  // gm cm^-3, mean density
double Recm = Re*100.0;

#define NNY  1000
#define DNY (double)NNY

#define Nref 1.78  //n for ice

double Unif[NNY],Yref[NNY];

double Enu=1.e20, E00=1.e20;

#define C0       2.9979e8    /* speed of light */
#define C        3.336 	    /* inverse vacuum speed of light, ns/m */

double  MAXDEPTH=3000.0;       /* maximum shower depth, m*/

// Globals for chord calculations and integrations;
double Ltot=0.0, K1=0.0, CHORDdepth=0.0, Rr=0.0;

#define MAXDIST1  2000.0e5   /* cm, distance from entry to allow vertices for downgoing*/
#define MAXDIST2  2000.0e5   /* cm, distance from exit to allow vertices from emerging */


//---------------------------------------------------------------------------------------

/* these are coefficients for a polynomial fit to ESS over 1e17-1e21 eV */
static double c[6] = { 1.661142492611783e+04,
		     -4.616049469872646e+03,
                      5.104845489878782e+02,
                     -2.812808939782252e+01,
                      7.727397573928863e-01,
                     -8.473959043935221e-03 };

double get_ESS_flux(lEnu)  // input: log10 of  Enu in eV, returns flux in cm^-2 s^-1 sr^-1
double lEnu;
{
   double lflux,flux;

  // first calculate log10(E * F(E)) which is what the fit gives,then tak antilog, dividing out lE
    lflux = c[0] + c[1]*lEnu + c[2]*pow(lEnu,2) + c[3]*pow(lEnu,3) + c[4]*pow(lEnu,4) + c[5]*pow(lEnu,5);
    flux = pow(10.,(lflux-lEnu));
  return(flux);
}

#define lESS_Emin 16.45
#define lFlux0   -1.688140237260814e+01   // log10(E*F(E)) at lESS_Emin

double get_Enu()  // return neutrino energy according to ESS spectrum		
{
      double x,y,wx,E,flux0;

      flux0 =  get_ESS_flux(lESS_Emin);
    //fprintf(stderr,"%e\n", flux0); sleep(1);

	wx = 0; y=1;
	while(y>=wx) {
		x=1.+ 31622.*drand48();  // gives up to 20.95 eV
		y=drand48();
		wx = get_ESS_flux(log10(x)+lESS_Emin)/flux0;
	 // fprintf(stderr,"%e %e %e\n", x, y, wx); sleep(1);
	     }
      E = pow(10.,lESS_Emin)*x;
      return(E);
}



/*>>>>>>>>>>>>>>>>>>  start main <<<<<<<<<<<<<<<<<<<<<<<<<*/
 	
main(argc,argv)
int argc; char **argv;
{
 double a,r,R,d,gam,gmin,gmax,dr,rmin,rmax,dgam,pirfac,Fratio,rho,rho2,theta;
 double beta(),acos(),cos(),sin(),exp(),sqrt(),dtime=0, thist[10000],dF,btot;
 double Gratio, lam1, lam2, lam3, lamrat, Aeff, b410,b480,b550,Ef;
 int idtime, i,j, ilam, iangle, Nray, Nevt,nabs, nearint(),ntir,iaz,iel;
 double  mcos, angle, drand48(), fabs(), Labs, Lscat,phi,psi, thetaC,cross_section_factor;
 double getLabs(),getLscat(),gettheta(),getphi(),getpsi(),get_shower_theta(),D;
 void rotate(), propagate(), get_fscat(), get_shower_dist();
 float gasdev();
 double SNRm, nSNR;
 double thet_n0, thet_out,depth, Stheta, Sphi, Spsi,pvel,depthm;
 double nbar,tandel,kw0,lambda,elev,Ve,Vn,SNR,Vn1,Vobs;
 double get_dens(),Ecr,Esh,upangle,Omega,get_shower_y(),yshower,Ntrials;
 float *Itheta, *Isum, *Ntheta, Ev_detected;
 double get_Esh(),get_Lint_cc(),get_Lint_nc(),get_depth(),get_nadir_angle();
 double get_updown(), get_thet(),Lint,Jy,Nev,get_dens(), Lintcc,Lintnc, get_Enu();
  int fp,fpsum, nue,evno, iLint, ESS;
  double nu0,numin,numax,integrate_Lint(),maxdepth,upang,updown,nadA,ang_corr,Dist,thet;
   struct timeval seed_time;

/*___________________________end of declarations__________________________*/

//USAGE:      

	if(argc < 3){
	fprintf(stderr,
"usage: anita_earthmodel2 [Nevt][ spectrum='ESS' or 'Enu'][Enu(eV)][maxdepth(m)][cross section factor] \n");
	fprintf(stderr, "If single energy 'Enu' then argv[3] is the neutrino energy, else ignored\n");
		exit(0);
		}
	ESS=0;
	Nevt = atoi(argv[1]);
	if(strcmp(argv[2],"ESS")== 0) ESS=1;
      fprintf(stderr,"ESS= %d\n", ESS); fflush(stdout);
	Enu = atof(argv[3]);
	maxdepth = atof(argv[4]);
	cross_section_factor = atof(argv[5]);
	MAXDEPTH = maxdepth*100.0;

 	nabs =0;
	Sphi = 0.0; Spsi=0.0;
	Nev = (double)Nevt;
	
	thetaC = acos(1./Nref);

	// a large prime seed; modified October 2005 --PG
	gettimeofday(&seed_time, (struct timezone*)0);
	//srand48(1299811);
	srand48(seed_time.tv_usec);
	fprintf(stderr,"#Seed value= %d\n", seed_time.tv_usec);


 fprintf(stdout,"#Enu= %e eV\t maxdepth= %e m\n", Enu, maxdepth);

/* --------------start the main loop on neutrino events HERE ----------- */

printf("#Ntrials\tevno\tnadA(d)\tupang(d)\tLint(kmwe)\tLtot(km)\tChord(km)\tdepth(m)\tcc/nc\tEnu(EeV)\n");
  Ntrials = 0.0; evno = 1;


 for(i=0;i<Nevt;i++){

 /* we loop on input angles, interaction length and depth
     until we get a final depth of less than MAXDEPTH */

	depth = 1.e8;

	while(depth > MAXDEPTH){   

	  // get the input nadir angle, from horizon to nadir
	  // the entry angle of the track takes care of downgoing events
	  // exit angle for upcoming events
	  nadA = get_nadir_angle();  // just throw the dice for cos(nadA)

	  // get the length of the chord:
  	  K1 = 2.*Recm*cos(nadA);

	  if(ESS==1) Enu=get_Enu();

	  //get neutrino interaction length (gm cm^2):
	  Lintcc = get_Lint_cc(Enu,cross_section_factor);
	  Lintnc = get_Lint_nc(Enu,cross_section_factor);
	  Lint = Lintcc>Lintnc? Lintnc: Lintcc;
	  iLint = Lintcc>Lintnc? 0 : 1;   // 0 for nc 1 for cc

	  // this is the depth of the middle of the chord:
	  //CHORDdepth = Recm*(1.0-sqrt(1.0-K1*K1/(4.0*Recm*Recm)));  NOT USED
	  // this is correct, but is not used  11/26/2005  --PG

	  depth = integrate_Lint(nadA,Lint); // integrate along track to get gm/cm^2 at interaction
	  				// this function returns the depth below nearest surface

	  // get the complement of upangle
	   thet = get_thet();  // also uses globals Ltot, K1
	  // upangle is the emergence angle wrt local surface ==> negative for downgoing
	   upang = thet-PI/2.0;
	   
	 // calculate depth at shower maximum ==> assume 3m further along track
	   depthm = depth - 300.0*upang;  // use cm here

	// exclude escaping tracks, those too close to surface, and those too large: 
	// >22 deg upangle is a problem 
	if( depthm<5.0 || upang<-.1 || upang >= PI/2.-thetaC-0.2)depth = MAXDEPTH+1.0; 
	// last part dumps this track
	
	Ntrials += 2.0; // downgoing on entrance and upcoming on exit, 2 trials
	}
	/* if we got out from above, it means we have a valid interaction, within allowed depth range*/

	Dist = (K1-Ltot);

	printf("%8d\t %6d\t %5.2f\t    %6.3f\t    %8.4f\t    %8.4f\t    %8.4f\t    %8.4f\t   %d\t   %8.4f\n",
		 (int)Ntrials,evno, nadA*180.0/PI,(upang)*180.0/PI,Lint/1.e5,Ltot/1.e5,K1/1.e5,depthm/1.e2,iLint,Enu/1e18);
	fflush(stdout);
	evno++;


    }

}  /* >>>>>>>>>>>>>>>>> end main <<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/* functions follow: */


 
#define HALF (0.500000000000)

int nearint(x)
double x; 
{ 
        double a, rem;
        int sign, i, newi;
        sign = (x<0) ? -1 : 1;
        a = fabs(x);
        i = floor(a);
        rem = a-(double)i;
        newi = (i + ((rem > HALF) ? 1 : 0)) * sign;
        return(newi);
}  




//---------------------------------------------------------------------------------------

double get_nadir_angle()  // returns angle that is flat in cosine distribution (solid angle)
		

{
  double drand48(), acos();

  return (acos(drand48()) );
  
  //return(asin(sqrt(drand48()))); test of Brian Mercurio's version 11/26/05
}

//---------------------------------------------------------------------------------------

// cross sections originally from Gandhi et al. 1998 hep-ph/9807264
// updated to Sarkar & Sarkar 2008, gives only cc cross section, scaled the others
/* produces an exponentially distributed interaction length, charged current  */
double get_Lint_cc(E,cross_section_factor)
double E,cross_section_factor;
{
   double sigcc,signc,sigbcc,sigbnc,avgsigcc,avgsignc,Na,L,sigbar,sigtot,expdev(),cc2nc_ratio;

	//sigcc = 5.53e-36*pow(E/1.e9,0.363); // Old Gandhi et al 1996 cross sections
	//signc = 2.31e-36*pow(E/1.e9,0.363);
	//sigbcc = 5.52e-36*pow(E/1.e9,0.363);
	//sigbnc = 2.29e-36*pow(E/1.e9,0.363);

	cc2nc_ratio = 2.39;
	//cc_plus_nc_totalfactor = 1.414;  // combined total cross section compared to cc

	/* neutrino to antineutrino ratio is very close to 1 */
	sigcc = 1.e-36 * exp ( 82.893 - 98.8*( pow( log(E/1.e9),-0.0964 ) ) );
	//sigtot = 7.82e-36*pow(E/1.e9,0.363);
	//sigbar = avgsigcc;

	sigcc *= cross_section_factor;

	Na = 6.022e23;


	L = 1./(Na*sigcc);  /* interaction length now in cm water equiv */


	return(L*expdev()); // returns total gm cm^2 interaction length
}

/* produces an exponentially distributed interaction length, neutral current  */
double get_Lint_nc(E,cross_section_factor)
double E,cross_section_factor;
{
   double sigcc,signc,sigbcc,sigbnc,avgsigcc,avgsignc,Na,L,sigbar,sigtot, expdev(),cc2nc_ratio;

	//sigcc = 5.53e-36*pow(E/1.e9,0.363);
	//signc = 2.31e-36*pow(E/1.e9,0.363);
	//sigbcc = 5.52e-36*pow(E/1.e9,0.363);
	//sigbnc = 2.29e-36*pow(E/1.e9,0.363);
	
	cc2nc_ratio = 2.39;
	signc = 1.e-36 * exp ( 82.893 - 98.8*( pow( log(E/1.e9),-0.0964 ) ) ) / cc2nc_ratio;
	
	//sigtot = 7.82e-36*pow(E/1.e9,0.363);
	//sigbar = avgsignc;

	signc *= cross_section_factor;

	Na = 6.022e23;


	L = 1./(Na*signc);  /* interaction length now in cm water equiv */


	return(L*expdev()); // returns total gm cm^2 interaction length
}

//----------------------------------------------------------------------------------------

// This function is based on the Preliminary Earth Model
// uses  polar radius now -- PG 9/25/04
// polar radius Re   6356.755e3 m  
// stupid bug fixed 2/27/2005: part of this  was using km instead of cm!!--PG
// another density profile bug fixed 6/5/2005
// modify last layers slightly 7/31/07


double get_dens(double rr) // return density as a function of radius, gm/cm^3
{
double x, get_firn_dens();
	x = rr/Recm;
	if(rr<1221.5e5) return(13.0885-8.8381*x*x);
	if(rr>=1221.5e5 && rr<3480.0e5)
	    return(12.5815-1.2638*x-3.6426*x*x-5.5281*x*x*x);
	if(rr>=3480.e5 && rr<5701.e5)
	    return(7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x);
	if(rr>=5701.e5 && rr<5761.e5)return(5.3197-1.4836*x);
	if(rr>=5761.e5 && rr<5960.e5)return(11.2494-8.0298*x);
	if(rr>=5960.e5 && rr<6140.e5)return(7.1089-3.8045*x);
	if(rr>=6140.e5 && rr<6335.e5)return(2.691+0.6924*x);
	if(rr>=6335.e5 && rr<6341.8e5)return(2.9);
	if(rr>=6341.8e5 && rr<6353.e5)return(2.6);
	if(rr>=6353.e5 && rr<6356.655e5) return(0.9); // ice above 3.755km deep
	if(rr>=6356.655e5 &&  rr<6356.755e5) return get_firn_dens((Recm-rr)); // this part for the firn
	if(rr>=6356.755e5) return(0.0);
}


//-----------------------------------------------------------------------------------------
//

double integrate_Lint(uA, LL)   // uA is the nadir angle
     double uA, LL;             // LL in interaction depth in gm cm^2
{
  double gmtot, L1,dL,d1,AA,get_dens(),DeltaL;

	if(uA>=PI/2.0){
		AA = PI-uA;
		}else{
		AA=uA;
		}

//fprintf(stderr,"chord=%e \n", K1/1.e5);
  gmtot = 0.0;
  L1 = 0.0;
  dL = 200.0;  // these are cm ==> 2m increments
  while(gmtot<LL){
    DeltaL = K1-L1; // added 11/26/05
    /*  the original line--did it work by symmetry, even though wrong?? */
    /* Rr = sqrt(Recm*Recm+L1*L1-2.0*Recm*L1*cos(AA)); original version?? wrong? */
    // above version does not follow cosine law directly , rather works by accident
    // correct version here 11/26/05, but the old is algebraically equal.
      Rr = sqrt(Recm*Recm + DeltaL*DeltaL - 2.0*Recm*DeltaL*cos(AA));
    if(Rr>Recm){
      Ltot = L1;
      d1 = -999999.0;
      break;
    }
    dL= (Rr>0.999*Recm)? 200.0 : 100000.0;
    gmtot += get_dens(Rr)*dL;  // total grammage along path
    L1 += dL;  // L1 is total path length to interaction
    d1 = Recm-Rr; // this is the distance below earth's surface for interaction
  }
  Ltot = L1;
//fprintf(stderr,"r1=%e  L1=%e  Ltot=%e\n", Rr/1e5,L1/1e5,Ltot/1e5);
  return d1;
}

//------------------------------------------------------------------------------------------------

// the complement of the angle between radius vector to interaction Rr and the incoming ray:
// checked again 11/26/2005  --PG
double get_thet()
{
	double thet = 0.0;
	return thet = acos(  (Rr*Rr + (K1-Ltot)*(K1-Ltot) - Recm*Recm )/(2.0*Rr*(K1-Ltot)) );
}

//----------------------------------------------
double get_firn_dens(double depth)  // depth in cm
{
// new routine added 2/27/2005 for firn density
// fit from gnuplot
//f(x) = 0.68 + a1*(1.0-exp(-b1*x));
//a1 = .22; b1 = .0131
// 
	double depth1;
	
	depth1 = depth/100.0;
	return (0.68 + 0.22*(1.0-exp(-0.0131*depth1)));
}

//------------------------------------------------------------------------------------------------
//exponential deviate

double expdev()
{
	double drand48();
	return -log(drand48());
}

//------------------------------------------------------------------------------------------------

//----------------------------END-----------------------------------------------------------------
