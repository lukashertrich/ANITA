/*
 * gcc anita3mcE10.c -o ~/bin/anita3mcE10 -lm -O3
 * 
anitamcE10: monte carlo radio emission from GZK neutrino ice showers.
This version used field strengths.
	
new version August 2003 based on GLUE MC

modified Feb. 2004 to generate 4 different frequency bands  --PG

modified September 2004 to make some updates & corrections.  --PG


March, 2005, added firn dependence in angles and field strength,
added different parameterization for high energies to avoid
energy saturaton, latter is small correction  --PG

Late March 2005, some corrections to how surface refraction is
handled. Added term for firn index matching, changed surface
refraction to use analytic forms rather than lookup.  --PG

March 26 2005, fixed incorrect handling of some neutral current 
events, y factor was reversed--PG

5/20/05:  error fixed on call get_field (was using nue instead of nut)

6/29/2005:  modified testing mode for sety, yin, etc.  Also, testing
indicates that 64x128 images are probably not fine-grained enough
for nue LPM showers at 3e18 eV. Tried 128x256 for tests, appears
OK at 3e18. --PG

6/29/2005:  added temporary fixes to failure of ZHS angle width param
at high energies above 1e19 eV  --PG

7/1/2005:  changed y distribution from lookup to von Neuman method
using parameterization from fits of Gandhi et al 1995.

10/13/2005: better random seed generator, corrected input BW;
improved angular width generation, added/debugged more firn corrections.  

7/31/07:  is there a mistake in rotation matrix, y-convention?
   seems that Goldstein p. 607 may be wrong, now using Wolfram
   mathworld version...does it make any difference?

8/1/07  added some handling for ray info output, ran some checks on
       ray tracing.

9/15/09:  This version now modified to account for ANITA-2 parameters

11/17/2010: this version for 10 frequency bands, used to synthesize full band

10/16/13, fixed firn density bugs....
--PG

*/
	
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>

//#include <values.h>

#define PI    M_PI

/* structure for raybundle & associated parameters */
typedef struct {
		 double x,y,z,t0,	/* position vector, start time */
			cx,cy,cz,	/* velocity  "     */
			L,		/* total distance  traveled */
			el,	/* elevation angle out wrt surface */
			az,	/* az angle out wrt surface */
			T,      /* total transmission coeff incl: absorp. & surface transm.,  */
			E,     /* magnitude of field*/
			Epar, 
			Eperp;
				} raybundle; 

raybundle ray;


typedef struct {
                 double had, // hadronic part of energy
		         em;  // electromagnetic (leptonic) energy
		 } showerE;
		 

#define DRef 600.e3	/* 600 km reference distance */
#define Re   6356.755e3     /* polar earth radius m */

#define RHOearth 2.75        /* most chords less than 100 km, not used currently 9/04*/

//#define NA   7200		/* number of angular bins to use */
//#define NA   14400		/* number of angular bins to use */
#define NA   28800		/* number of angular bins to use */
#define DNA  (double)NA
#define NNY   1000		/* number of bjorken y values */
#define DNY  (double)NNY
#define NEL  128 		/* size of angular intensity array */
#define NAZ  256 		/* size of angular intensity array */
#define NI   NEL*NAZ
#define NFREQ 10

showerE  Dtheta[NFREQ];   // global for keeping track of angular width of Cherenkov cone, 9/04--PG
double Sham[NFREQ][NA+1];  /* array of cumulative shower angular dist fn */

float *Thetout;
float *Transout; 
int NFIRN = 346409; // firn ray trace lookup file
#define NTOUT = 350000;

/* FOR ANITA-2: elog note 381  Lowband=200-350MHz, Mid1=330-600, Mid2=630-1100, Full=200-1200. */



static double NUBAR[NFREQ]= { 250.e6,350e6,450e6,550e6,650e6,750e6,850e6,950e6,1050e6,1150e6 }; // center frequencies, unweighted
static double BW[NFREQ] = {100e6,100e6,100e6,100e6,100e6,100e6,100e6,100e6,100e6,100e6};  // 3 dB bandwidths

  /* when the frequency bands are weighted by the antenna effective height, we have to take the
     weighed average of f*heff, ==> we get <f>_w = (fmax-fmin)/ln(fmax/fmin) */
/* with the change to narrower bands the difference is negligible */
static double NUBARW[NFREQ]= {250.e6,350e6,450e6,550e6,650e6,750e6,850e6,950e6,1050e6,1150e6}; 

static double alpha0[NFREQ] = {5.793480e-04 , 6.644813e-04 , 7.393419e-04 , 7.960876e-04 , 
				8.406714e-04 , 8.782108e-04 , 9.123559e-04 , 9.460117e-04 , 9.796675e-04 , 1.013323e-03 };

static double dnufac[NFREQ] = {5,5,5,5,5,5,5,5,5,5}; // number of frequency integration intervals per band

//globals for y, brems, pairs, and photonuclears
double Unif[NNY],Yref[NNY];  // for Bjorken y -distribution calcs
double mubremUnif[120],mupairsUnif[120], mupnUnif[120];  // these all for CDFs for taus & mus
double mubremYref[120],mupairsYref[120],mupnYref[120];
double taubremUnif[120],taupairsUnif[120], taupnUnif[120];
double taubremYref[120],taupairsYref[120],taupnYref[120];
double Nmubrem,Nmupairs,Nmupn;
int mubremwc, mupairswc, mupnwc;
double Ntaubrem,Ntaupairs,Ntaupn;
int taubremwc, taupairswc, taupnwc;

#define K          5.9375   /* inverse speed of light in ice in ns/m, Nref=1.78 */
#define C0         2.9979e8    /* speed of light */
#define C          3.336 	    /* inverse vacuum speed of light, ns/m */
#define H          8.0       /* ~longitudinal size (m) of cascade shower max */
#define Nref       1.78      /* index of refraction of ice */
#define Nfirn      1.325     /* index of refraction at the surface */
#define Nair       1.000186   /* radio index of refraction of dry air, 3km altitude, -20C*/
#define TIR        0.604594    /* total internal reflectance angle wrt normal */
#define NU0        1.15e9     /* reference CR output frequ */
#define ROUGH      0.0280   /* surface local slope rms roughness, 1 deg=0.0175 rad std, 1.7 degrees from Taylor data*/
#define MEANSLOPE  0.0175   /* mean surface slope ~km scales, 1 deg=0.0175 rad.std*/


double  MAXDEPTH=3000.0;       /* maximum shower depth */

double NUMIN[NFREQ],NUMAX[NFREQ];

double E00 = 1.e20;
double Enu = 1.e20;
#define MAXDIST  700e3   /* meters */
#define ZHSFAC   0.707   // the dreaded ZHS field strength correction

/*>>>>>>>>>>>>>>>>>>  start main <<<<<<<<<<<<<<<<<<<<<<<<<*/
 	
main(argc,argv)
int argc; char **argv;
{

/* ______________________local variables ________________________________*/
 double  a,r,R,d,gam,gmin,gmax,dr,rmin,rmax,dgam,pirfac,Fratio,rho,rho2,theta;
 double  beta(),acos(),cos(),sin(),exp(),sqrt(),dtime=0, dF,btot;
 double  Gratio, lam1, lam2, lam3, lamrat, Aeff, b410,b480,b550,Ef;
 int     idtime, i,j, ilam, iangle, Nray, Nevt,nabs, nearint(),ntir,iaz,iel,ESS;
 double  mcos, angle, drand48(), fabs(), Labs, Lscat,phi,psi;
 double  getLabs(),getLscat(),gettheta(),getphi(),getpsi(),get_shower_theta(),D;
 void    rotatey(), propagate(), get_fscat(), get_shower_dist();
 void    set_ydist(),set_mu_ydist(),set_tau_ydist(), set_shower_dist();
 double  gasdev(), trialn=0.0,evno,Ltot,Chord,upang;
 double  SNRm, nSNR,rp,tp, Esum1,radperpix,Elangle, Ev_detected, Ptot[NFREQ];
 double  thet_n0, thet_out,depth, Stheta, Sphi, Spsi,pvel;
 double  nbar,tandel,kw0[NFREQ],lambda[NFREQ],elev,Ve,Vn,SNR,Vn1,Vobs;
 double  get_dens(),Ecr,upangle,Omega,get_shower_y(),yshower,sinp,tfirn;
 float   *Etheta[NFREQ],*EthetaV[NFREQ], *EthetaH[NFREQ],*EthetaP[NFREQ],*Ntheta[NFREQ];
 double  get_Esh(),get_Lint(),get_depth(),get_shower_upangle(),Lint,Jy,Nev, get_Ptot();
 showerE Esh, get_shower_E(); 
 //extern  double Sham[][];
 char    stmp[80],imgname[80], *get_E(); 
 int     fp[NFREQ],fP[NFREQ],fpsum, nue, nut,sety,nutin,current, roughON, slopeON, RayON,firstray=1;
 double  nu0,numin[NFREQ],numax[NFREQ],get_field(),heff,gain,polfac,lam,maxdepth, Efmax[NFREQ];
 double  nadA,ang_corr,Dist,thet,get_thet(), Ehad, Eem,yin, cr,ci, thet_firn,tpar,tperp; 
 double  get_thetout(),get_transout(),get_avg_alpha(),thetamax,Enu1;
 char    *getE(), *Energy, *Energy1;
 FILE    *fray, *fopen();  // these are for ray outputs
  struct timeval seed_time;
  
/*___________________________end of declarations__________________________*/

USAGE:      

	if(argc < 9){
	fprintf(stderr,
"usage: anita3mcE10 [N rays][Nevt][Enu][maxdepth][Image basename][roughness? 1/0][slopes? 1/0][rays out?(1/0)] <eventfile >logfile\n");
	fprintf(stderr,"\t\tfor single energy Enu= energy in eV, for ESS spectrum, energy=0.0\n");
		exit(0);
		}



	// a large prime seed; modified October 2005 --PG
	gettimeofday(&seed_time, (struct timezone*)0);
	//srand48(1299811);
	srand48(seed_time.tv_usec);
	fprintf(stderr,"#Seed value= %d\n", seed_time.tv_usec);
	
	ESS=0;
	Nray = atoi(argv[1]);
	Nevt = atoi(argv[2]);
	Enu = atof(argv[3]);
	  if(Enu==0.0) {
	      ESS=1; 
	    fprintf(stderr,"ESS= %d\n", ESS); fflush(stdout);
	  }
	maxdepth = atof(argv[4]);
	roughON = atoi(argv[6]);
	slopeON = atoi(argv[7]);
	RayON = atoi(argv[8]);
	
	
	
	if(ESS != 1) {
	    Energy = get_E(Enu);  // get a string for the energy-dependent CDF file handles
	    fprintf(stderr,"Energy = %s\n", Energy);
	    fflush(stderr);
	}
	
	// old method for y distribution--was inaccurate! --PG 7/1/2005
	/* set_ydist();   read in Bjorken y cumulative file */
	//----------------------------------
	if(ESS != 1){
	  set_mu_ydist(Energy);
	  set_tau_ydist(Energy);
	fprintf(stderr,"all lepton files read...\n");
	fflush(stderr);
	}

	
	//set_thetout();  /* read in and set the angle lookup tables for firn*/

	sety = 0;  // 1 for diagnostics, 0 for normal
	yin = 1.0;
	//nutin = 0;  // electron neutrino for diagnostics
	nutin = 1;  // muon neutrino, for diagnostics
	//nutin = 2;  // muon neutrino, for diagnostics
	
	 MAXDEPTH = maxdepth;

	// some useful arrays of frequencies
	for(j=0;j<NFREQ;j++){
	  NUMIN[j]=(NUBAR[j]-BW[j]/2.);
	  NUMAX[j]=(NUBAR[j]+BW[j]/2.);
	  numin[j] = NUMIN[j]/1.e6;
	  numax[j] = NUMAX[j]/1.e6;
	  }
	
	nu0 = NU0/1.e6;  /* Cherenkov reference freq in MHz */
	  
 	nabs =0;
	Sphi = 0.0; Spsi=0.0;
	Nev = (double)Nevt;


	/* create the NFREQ ouput image files for each frequency band, and plane-of-pol*/
	for(j=0;j<NFREQ;j++){
          sprintf(imgname,"%s%d.img", argv[5],(int)(NUBAR[j]/1.e6));
          if((fp[j] = creat(imgname,0644))<=0){
          fprintf(stdout,"open failure on %s\n", imgname, j);
           }
	  sprintf(imgname,"%s%dP.img", argv[5],(int)(NUBAR[j]/1.e6));
	  if((fP[j] = creat(imgname,0644))<=0){
          fprintf(stdout,"open failure on %s\n", imgname, j);
           }
	 }
	 
	 /* open the ray output files */
	if(RayON) fray = fopen("rays.out","w");
		

	Ev_detected = 0.0;
    //  read over the two event file header lines to position for the read
   	fscanf(stdin,"%*s%*s%*s%*s%*s%*s");
  	fscanf(stdin,"%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s");

/*-------------------------------------------------------------------------*/
/* --------------start the main loop on neutrino events HERE ----------- */
/*-------------------------------------------------------------------------*/
while(Nevt > 0){

  /* allocate the intensity array and hit counter*/
      for(j=0;j<NFREQ;j++){
	Etheta[j] = (float*)calloc(NI, sizeof(float));
	EthetaV[j] = (float*)calloc(NI, sizeof(float));
	EthetaH[j] = (float*)calloc(NI, sizeof(float));
	EthetaP[j] = (float*)calloc(NI, sizeof(float));
	Ntheta[j] = (float*)calloc(NI, sizeof(float));
	}

	
   // read in the event geometry parameters from the event file, written as:
   // printf("#nadA(d)\tupang(d)\tLint(kmwe)\tLtot(km)\tChord(km)\tdepth(m)\n");


	fscanf(stdin,"%lf%lf%lf%lf%lf%lf%lf%lf%d%lf",&trialn,&evno,&nadA,&upang,&Lint,&Ltot,&Chord,&depth,&current,&Enu1);
//   fprintf(stderr,"Esh.em, Esh.had,nut,nadA,upang,Lint,Ltot,Chord,depth,current,Enu\nn");
	
  // get neutrino type: 0 for e, 1 for mu, 2 for tau
        nut = (sety == 1) ? nutin : (int)(3.0*drand48());  // truncate a random number
	
	if(ESS==1){
   
	  Enu = Enu1*1e18;  // Enu1 is in EeV
// fprintf(stderr, "Enu = %e\n", Enu);
	  if(Enu>=3e17){
	    Energy1 = get_E(Enu);
	    }else{
	    Energy1 = get_E(3e17); // since we only have CDFs down to 3e17
	    }
//fprintf(stderr, "Enu = %s\n", Energy1);
	  set_mu_ydist(Energy1);
	  set_tau_ydist(Energy1);
//	fprintf(stderr,"all lepton files read...\n");
//	fflush(stderr);
        }
  // get the Bjorken y-factor, returns a structure with hadronic & EM energy
        Esh = get_shower_E(Enu, nut, sety, yin, current);

   /*set up the cherenkov angular distribution functions-- Esh is a structure here*/
     
	set_shower_dist(Esh,nut,depth);

//	for(i=0;i<NA;i++)printf("%d  %e\n", i,Sham[0][i]);
//	exit(0);

	/*  ADD A RANDOM SLOPE COMPONENT TO UPANG */
	if(slopeON) upang += MEANSLOPE*180./PI*((double)gasdev());  // convert mean slope to degrees, add random to upang
	fprintf(stdout,"%9.0f %e %e %d %e  %e  %e  %e  %e %e %d %e\n",
	trialn,Esh.em,Esh.had,nut,nadA,upang,Lint,Ltot,Chord,depth,current,Enu1);
	fflush(stdout);
	
	nadA *= PI/180.0; // now go back to radians for the calcs
	upang *= PI/180.0;


   // this angle used in calculations below
	Stheta = PI/2. - upang;


/* determine the RF absorption */
	nbar = Nref;
	for(i=0;i<NFREQ;i++){
		 lambda[i] = C0/NUBAR[i];
		 kw0[i] = get_avg_alpha(alpha0[i],depth);
		 kw0[i] = alpha0[i];
//		 fprintf(stderr,"depth= %f, avg alpha=%e\n",depth, kw0[i]);
		 }

/* this is the field attenuation coefficient vs. frequency
   for ice at average -50 C (eg Vostok upper 2km) */
/* the attenuation will be higher on the ice shelves or nearer coast 
	kw0[0] = 1./1500.0;
	kw0[1] = 1./1200.0;
	kw0[2] = 1./1050.0;
	kw0[3] = 1./900.0; */

/*	fprintf(stderr,"tandel = %7.5f,  Latten= %7.5f\n", tandel, 1./kw0); */

/* =========== loop over different frequency bands============================ : */
   for(j=0;j<NFREQ;j++){

	Efmax[j] = -99;
	thetamax=-99;
	for(i=0;i<NI;i++){
	   Etheta[j][i] = 0.0;
	   EthetaV[j][i] = 0.0;
	   EthetaH[j][i] = 0.0;
	   EthetaP[j][i] = 0.0;	   
	   Ntheta[j][i] = 0.0;
	   }
  /* -----------loop over Nray rays in the shower : */  
    for(i=0;i<Nray;i++){

   /* initialize ray along track given by shower on rotated axis*/
	ray.t0 = 0.0; 
        ray.x = 0; ray.cx = 0;
	ray.y = 0; ray.cy = 0;
	ray.z = -depth; 
	ray.cx = 0;
	ray.cz = (1./K); /*ray along theta direction, m/ns*/
	ray.L = 0; 


//	ray.cz = 1./K*0.746;	// puts the shower at the cherenkov angle 
//	ray.cy = 1./K*0.666;  


   /* now re-aim the ray according to expected shower distribution: */
	theta = get_shower_theta(j);
//	theta = acos(1./Nref);

	phi = 0.0 ; //psi = -PI/2; 
//	phi = getphi();  //phi is not used in these coordinates...
	psi = getpsi(); 
	rotatey(&ray,theta,phi,psi);/* rotate through Euler angles, y-conv */

   /* now rotate by shower angles */
	rotatey(&ray,-Stheta, Sphi,Spsi); 
	
/*	get field strength at 1 m for this angle, insert ZHS factor here */ 
	Ef = ZHSFAC*get_field(Esh,1.0,nu0,numin[j],numax[j],theta,nut,depth,dnufac[j]);

	if(Ef>Efmax[j]){
		Efmax[j] = Ef;
		thetamax=theta;
		} 


  /* now calculate distance to surface along ray trajectory*/
	D = (ray.cz < 0.0) ? 1e4 : (depth/ray.cz)/K;
	if(D==1e4) continue;


  /* ray velocity as a check */
  /*	pvel = sqrt(ray.cx*ray.cx+ray.cy*ray.cy+ray.cz*ray.cz);

       printf("Rph = %6.4f,%6.4f,%6.4f,%e\n", ray.x, ray.y, ray.z, D); 
       printf("%6.4f  %6.4f  %6.4f  %e  %e\n", ray.cx, ray.cy, ray.cz, D,pvel*K); 
	continue; */

	// ----------this loop propagates the individual rays
    while(ray.z <= 0.0){ 

	 propagate(&ray,D); /* to get exit posn in x-y plane */


  //    printf("Rph= %6.4f %6.4f %6.4f\n", ray.x, ray.y, ray.z); 
/*   printf("Cph = %6.4f,%6.4f,%6.4f\n", ray.cx, ray.cy, ray.cz); */

  	/* the  angle wrt normal just below FIRN underside surface */
	/* and for now we ignore the ray displacement in the firn */
	thet_n0 = atan2(hypot(ray.cx,ray.cy), ray.cz);
	 
	/* now add a random component to represent surface roughness */
	/* this approximates a facet model, but including only slopes and
	   not piston offsets of the surface. The effect is that all facets which
	   yield output rays in the same direction are assumed to add coherently in
	   the far field. The scale length of the facets is effectively the mean
	   spacing of the number of rays traced. Probably the most important effect
	   is that it allows rays to escape from slightly beyond TIR */
	if(roughON) thet_n0 += ROUGH*(double)gasdev();
		
	
	/* compute the ray polarization components:

  	/* the refracted angle */
	thet_out = (Nref/Nair*sin(thet_n0) <= 1.0) ? 
			asin(Nref/Nair*sin(thet_n0)) : -99.; 
			
		if(thet_out == -99.){ /* this ray suffers TIR */
/*			printf("%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.6f\n",
			  ray.az,ray.el,ray.x,ray.y,ray.L,ray.T); */
			ray.T = 0.0;
			ntir++;
			break;
			}		
				
	ray.az = asin(ray.y/hypot(ray.y,ray.x));
	ray.el = PI/2.0 - thet_out;
	
	// Polarization angles seen from payload:
	sinp = cos(-upang*PI/180.)*sin(ray.az)/sin(theta); 
	ray.Eperp = Ef*sinp;
	ray.Epar = Ef*sqrt(1.-sinp*sinp);


	if(thet_n0<0.0)break;  // hack to fix rare unknown problem
	//thet_out = get_thetout(thet_n0)

	// need angle just before firn exit for the modified Fresnel coeff.: 
	thet_firn = asin(Nair/Nfirn*sin(thet_out));

	/* now weight down by absorp*transmission coefficient */

	// These terms are from Lehtinen et al. 2004 (FORTE paper):
	// Modified Fresnel coefficient for ice-->firn-->air
	// But there is an additional term for the index matching
	// due to continuous gradient of firn (below).
	cr = cos(thet_out);
	ci = cos(thet_firn);
	tpar = 2.*cr/(Nfirn*cr + ci);
	tperp = 2.*cr/(Nfirn*ci + cr);
	// This adds a correction for index matching of the firn
	//  (from N. Lehtinen note, ANITA elog 2005)
	tfirn = sqrt(Nref*cos(thet_n0)/(Nfirn*cos(thet_firn)));
	
//	fprintf(stdout, "Nrays=%d cr=%e ci=%e tp=%e\n",Nrays, cr,ci,tp);
//	fflush(stdout);
//	fprintf(stderr,"tp=%e  exp=%e ray.T = %e\n", tp, exp(-kw0[j]*D), ray.T);
	
	ray.T = exp(-kw0[j]*D); // just include attenuation first
	
	if(ray.T <= 0.000001 || thet_out == -99.) break;  /* too much absorption */ 
	
	// Epar==parallel to plane of incidence == Vertical pol.
	//Eperp == perp to plane-of-inc.== Hpol
	ray.Epar  *= tfirn*tpar*ray.T;
	ray.Eperp *= tfirn*tperp*ray.T;
	
	//ray.Epar  *= ray.T;
	//ray.Eperp *= ray.T;
	
	ray.E = hypot(ray.Eperp,ray.Epar);//use the magnitude of E for sky map
	
	
			
		iaz = nearint((ray.az+PI/2.)/PI * NAZ);
		iel = (int)((ray.el)/(PI/2.) * NEL);

/* fill in with peak field strength at this distance, weighted by
transmission & absorption factor, units are V/m/MHz  */

		Etheta[j][iaz+NAZ*iel] += (float)(ray.E/DRef);
		EthetaV[j][iaz+NAZ*iel] += (float)(ray.Epar/DRef);
		EthetaH[j][iaz+NAZ*iel] += (float)(ray.Eperp/DRef);
		Ntheta[j][iaz+NAZ*iel] += 1.0;
		
//	if(j==0 && fabs(ray.az)<1e-2 && fabs(theta-acos(1./Nref))<1e-3)
//	fprintf(stdout,"Ef=%e Etheta=%e Ntheta=%e ray.Epar=%e ray.az=%e ray.el=%e tpar=%4.3f tperp=%4.3f tfirn=%4.3f theta=%e\n", 
//		Ef*BW[0]/1e6,Etheta[j][iaz+NAZ*iel],Ntheta[j][iaz+NAZ*iel],ray.Epar*BW[0]/1e6,ray.az,ray.el,tpar,tperp,tfirn,theta);


OUTPUT1:		// used to get ray patterns out
		 if(RayON){
		 	if(firstray){
					 fprintf(fray,"#AZ,d     EL,d      X,m      Y,m       D,m      T\n");
					 firstray=0;
				  }
			if(j==2){  // set j to the frequency band you want
			fprintf(fray,"%+7.4f  %+7.4f  %+7.4f  %+7.4f  %+7.4f  %+7.6f\n",
			  ray.az*180/PI,ray.el*180/PI,ray.x,ray.y,ray.L,ray.T);
			fflush(stdout);
			}  
		   }

		   break;    /* get out and start another ray */


      }    /* end of ray while */

     }    /* end of Nray loop */
     
     	//printf("Efmax = %e V/m/MHz at 1m, thetamax= %f\n", Efmax[j],thetamax);
	
   }  /* end of frequency bin loop */
   
 
/*%%%%%%%%%%%%%%%%%%%%% print out final results %%%%%%%%%%%%%%%%%%%%%%%%*/

   //  normalize the field strength by the number of ray bundles & the angular area of the pixel--
   //  the latter is necessary since we are using a very simple projection on the sky,
   // and the pixels are getting effectively smaller and smaller as we approach the zenith

	radperpix = (double)(PI/2.0/NEL);
	for(j=0;j<NFREQ;j++){
//		fprintf(stderr,"Efmax%d = %e V/m\n", j+1, Efmax[j]*BW[j]/1e6);
		}

        for(j=0;j<NFREQ;j++){
	   for(i=0;i<NI;i++){
		Elangle = (double)(0.5 + i/NAZ)*radperpix;
		if(Etheta[j][i] == 0.0) continue;
		Etheta[j][i] = Etheta[j][i]/(Ntheta[j][i]);  
		Etheta[j][i] = (Etheta[j][i]/cos(Elangle));// correct for Mercator projection
		// Vertical pol
		EthetaV[j][i] = EthetaV[j][i]/Ntheta[j][i];  
		EthetaV[j][i] = (EthetaV[j][i]/cos(Elangle));// correct for Mercator projection
		// H pol
		EthetaH[j][i] = EthetaH[j][i]/Ntheta[j][i];  
		EthetaH[j][i] = (EthetaH[j][i]/cos(Elangle));// correct for Mercator projection
		
		EthetaP[j][i] = atan2(EthetaV[j][i], EthetaH[j][i]);
	    }
	 }	
   	Nevt--;

	for(j=0;j<NFREQ;j++){
        write(fp[j],(float*)Etheta[j], NI*sizeof(float));
        write(fP[j],(float*)EthetaP[j], NI*sizeof(float));	  
	  free(Etheta[j]);
	  free(EthetaV[j]);
	  free(EthetaH[j]);
	  free(EthetaP[j]);
	  free(Ntheta[j]);
	  
         }


 }   
/*-------------------------------------------------------------------------*/
/* --------------end the main loop on neutrino events HERE ----------- */
/*-------------------------------------------------------------------------*/

  	if(RayON)fclose(fray);

  	for(j=0;j<NFREQ;j++)close(fp[j]); 


	exit(0);

}  /* >>>>>>>>>>>>>>>>> end main <<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/* functions follow: */
//-------------------------------------------------------------------	
void propagate(fot, L)
raybundle *fot;
double L;  /* distance to propagate along cx,cy,cz */
{
	double t;
	t=L*K;   /* L*K is now in ns */

	fot->L += sqrt(fot->cx*fot->cx +
		  fot->cy*fot->cy +
		  fot->cz*fot->cz)*t;
	fot->x +=  fot->cx*t;  
	fot->y +=  fot->cy*t;
	fot->z +=  fot->cz*t;

}

//-------------------------------------------------------------------

double getLabs(a)   /* a is the abs coeff in m^-1, this just for testing*/
double a;
{
	double drand48();
 /* equal prob per unit length=> exp*/
	return (1.0/a)*(-log(drand48()));
}

//-------------------------------------------------------------------
void rotateSB(fot,theta,phi,psi)  /* uses Euler angles, space to body axes */
raybundle *fot; double theta,phi,psi;
{
    double tmpx,tmpy,tmpz, cos(),sin();

	tmpx = fot->cx*(cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi))
		+ fot->cy*(cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi))
		+ fot->cz*(sin(psi)*sin(theta));
	tmpy = fot->cx*(-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi))
		+fot->cy*(-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi))
		+fot->cz*(cos(psi)*sin(theta));
	tmpz = fot->cx*sin(theta)*sin(phi)
		-fot->cy*sin(theta)*cos(phi)
			+fot->cz*cos(theta);
	fot->cx= tmpx;   fot->cy = tmpy;    fot->cz = tmpz;
}
//-------------------------------------------------------------------
void rotateBS(fot,theta,phi,psi)  /* uses Euler angles, body to space axes */
raybundle *fot; double theta,phi,psi;
{
    double tmpx,tmpy,tmpz, cos(),sin();

	tmpx = fot->cx*(cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi))
		+ fot->cy*(-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi))
		+ fot->cz*(sin(theta)*sin(phi));
	tmpy = fot->cx*(cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi))
		+fot->cy*(-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi))
		+fot->cz*(-sin(theta)*cos(phi));
	tmpz = fot->cx*sin(theta)*sin(psi)
		+fot->cy*sin(theta)*cos(psi)
			+fot->cz*cos(theta);
	fot->cx= tmpx;   fot->cy = tmpy;    fot->cz = tmpz;
}
//-------------------------------------------------------------------
void rotatey(fot,theta,phi,psi)  /* Euler angles, space to body axes, y conv */
raybundle *fot; double theta,phi,psi;
{
    double tmpx,tmpy,tmpz,ctheta,stheta,cphi,sphi,cpsi,spsi;
	ctheta=cos(theta); stheta=sin(theta); 
	cphi=cos(phi); sphi=sin(phi);
	cpsi=cos(psi);  spsi=sin(psi);

	tmpx = fot->cx*(-spsi*sphi+ctheta*cphi*cpsi)
		+ fot->cy*(spsi*cphi+ctheta*sphi*cpsi)
		+ fot->cz*(-cpsi*stheta);
	tmpy = fot->cx*(-cpsi*sphi-ctheta*cphi*spsi)
		+fot->cy*(cpsi*cphi-ctheta*sphi*spsi)
		+fot->cz*(spsi*stheta);
	tmpz = fot->cx*(stheta*cpsi)  // according to Goldstein, p. 607--wrong??
		+fot->cy*(stheta*spsi)
			+fot->cz*ctheta;
//	tmpz = fot->cx*(stheta*cphi)  // from Wolfram mathworld
//		+fot->cy*(stheta*sphi)
//			+fot->cz*ctheta;
	fot->cx= tmpx;   fot->cy = tmpy;    fot->cz = tmpz;
}

//-------------------------------------------------------------------
double get_shower_theta(int ifreq)  /* from cumulative shower function */
{
	/* this function takes a ray, should be aligned with
           shower axis, and returns angle in a direction drawn according
	   to the cascade cherenkov light distribution */
   double dr,angle, drand48(), theta,phi,psi;
   int iangle, nearint();
   double  getphi(), getpsi();

	dr = (1.-drand48())*DNA;
	iangle = (int)nearint(dr);
	theta = PI*(Sham[ifreq][iangle]);
	return(theta);

	
}
//-------------------------------------------------------------------

double getphi()   /* uniform distribution */
{
   double angle, drand48();  

	angle = drand48()*2.0*PI;
	return(angle);
}
//-------------------------------------------------------------------
double getpsi()   /* uniform distribution */
{
   double angle, drand48();  

	angle = drand48()*2.0*PI;
	return(angle);
}

//-------------------------------------------------------------------
 
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
//-------------------------------------------------------------------
double gasdev()
{
/* returns gaussian with zero mean unit variance */
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;
        double drand48();

        if  (iset == 0) {
                do {
                        v1=2.0*drand48()-1.0;
                        v2=2.0*drand48()-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
} 
//-------------------------------------------------------------------
showerE get_shower_E(double Enu, int nutype, int sety, double yin, int current)  
/* uses tabulated cumulative Bjorken y distribution */
{

   double Bjy, By, dny, drand48(),flavor,i,expdev(),pow();
   double Bjy1,Ltau,Etau,eps=1.e-10, dr48;
   double kpoi(),Bjyem1,Bjyem2,Bjypn,Emu;
   int iy,iymax,ibrem,ipairs,ipn, nearint(),inue=0,inumu, inutau=0, decay;
   double tgamma, Ptau;
   extern double Yref[];
   double get_yvalue();
   showerE Esh;

	/* first determine if its an electron, muon or tau neutrino */
	flavor = nutype;
	inue = flavor == 0 ? 1 : 0;
	inumu = flavor == 1 ? 1 : 0 ;
	inutau = flavor >= 2 ? 1 : 0;
	
	if(sety == 1){   // for diagnostics
	
	    if(current == 1) {   // charged current test 
	    
		if(inue){
		  Bjy = yin;
		  Esh.em =   (1.-Bjy+eps)*Enu;  
		  Esh.had = (Bjy+eps)*Enu;
		  }
		if(inumu || inutau){
		  Bjy = yin;
		  Esh.em =   1.e8;  
		  Esh.had = (Bjy+eps)*Enu;
		  }
		  
		}else{   //  neutral current test
		
		if(inue){
		  Bjy = yin;
		  Esh.em =  1.e8;  
		  Esh.had = (Bjy+eps)*Enu;
		  }
		if(inumu || inutau){
		  Bjy = yin;
		  Esh.em =   1.e8;  
		  Esh.had = (Bjy+eps)*Enu;
		  }
	     }
	     return(Esh);
	 }
		  

	// if we got here, we are not in testing mode
	// OLD METHOD FOR GETTING Y DISTRIBUTION
	/* dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 )*DNY;
	iy = (int)floor(dny);
	Bjy = Yref[iy];  */
	
	//--------------------------------
	//  new method  7/1/2005
	Bjy = get_yvalue();
	//-----------------------------

	 /* for neutral current events (current=0), just make hadronic shower only,
	     for any lepton flavor */
	   /* for charged current we use cumulative distribution functions from Ped Miocinovic,
	       derived from running the AMANDA monte carlo mmc, these give shower probabilities
	       for all of the different processes, brems, pairs, and photonuclears*/
	if(inue==1){
		Esh.em = (current==1) ? Enu*(1.0-Bjy+eps): 1.e8;
		Esh.had = Enu*(Bjy+eps);
		}
	if(inumu == 1){
	    	if(current == 1) { // charged current
			// now step through brems, pairs, & pn to see if a shower occurs with the muon
			Emu = Enu*(1.-Bjy);
			iymax = -99;
			//brems:
			ibrem = (int)lrint(kpoi(Nmubrem)); // Poisson on the number of showers in this process
			for(i=0;i<ibrem;i++){
				dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 );
				if(dny<mubremUnif[0]) continue;  //since some of the CDFs have y cutoffs
				dny *= mubremwc-1.; // gives the bin number in distribution
		 		iy = (int)floor(dny);
				if(iy>iymax){
					iymax = iy;
				}
			}
			Bjyem1 = iymax>0? mubremYref[iymax]: eps;
			//end of brems
			iymax = -99;			
			//pairs:
			ipairs = (int)lrint(kpoi(Nmupairs)); // Poisson on the number of showers in this process
			for(i=0;i<ipairs;i++){
				dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 );
				if(dny<mupairsUnif[0]) continue;  //since some of the CDFs have y cutoffs
				dny *= mupairswc-1.; // gives the bin number in distribution
		 		iy = (int)floor(dny);
				if(iy>iymax){
					iymax = iy;
				}
			}
		  	Bjyem2 = iymax>0? mupairsYref[iymax]: eps;
			//end of pairs
			iymax = -99;			
			//photonuclear:
			// this is a fit to the mean number of showers vs. energy
			Nmupn = pow(10.,1.80*(1.-exp(-1.17*(log10(Emu)-17.87)))-1.);
			//printf("Emu=%e Nmupn= %e\n", Emu, Nmupn);
			ipn = (int)lrint(kpoi(Nmupn)); // Poisson on the number of showers in this process
			for(i=0;i<ipn;i++){
				dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 );
				if(dny<mupnUnif[0]) continue;  //since some of the CDFs have y cutoffs
				dny *= mupnwc-1.; // gives the bin number in distribution
		 		iy = (int)floor(dny);
				if(iy>iymax){
					iymax = iy;
				}
			}
			Bjypn = iymax>0? mupnYref[iymax]: eps;	
			// end of photonuclears
			//printf("Bjyem1=%e  Bjyem2=%e  Bjypn=%e  Emu=%e\n", Bjyem1, Bjyem2,Bjypn,Emu);
			Esh.em = Bjyem1>Bjyem2 ? Emu*Bjyem1: Emu*Bjyem2;  
			Esh.had =  Enu*(Bjy+eps) > Bjypn*Emu ? Enu*(Bjy+eps) : Bjypn*Emu;
			if(Esh.had >= Esh.em){
				Esh.em = eps*Enu;
				}else{
				Esh.had = eps*Enu;
				}
		   }else{     // neutral current
			Esh.em =   1.e8;  
			Esh.had = (Bjy)*Enu;
		       }
		 }
		
        if(inutau == 1){
	
	    	if(current == 1) { // charged current
			// now step through brems, pairs, & pn to see if a shower occurs with the tauon
			Etau = Enu*(1.-Bjy);
			iymax = -99;
			//brems:
			ibrem = (int)lrint(kpoi(Ntaubrem)); // Poisson on the number of showers in this process
			for(i=0;i<ibrem;i++){
				dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 );
				if(dny<taubremUnif[0]) continue;  //since some of the CDFs have y cutoffs
				dny *= taubremwc-1.; // gives the bin number in distribution
		 		iy = (int)floor(dny);
				if(iy>iymax){
					iymax = iy;
				}
			}
			Bjyem1 = iymax>0? taubremYref[iymax]: eps;
			//end of brems
			iymax = -99;			
			//pairs:
			ipairs = (int)lrint(kpoi(Ntaupairs)); // Poisson on the number of showers in this process
			for(i=0;i<ipairs;i++){
				dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 );
				if(dny<taupairsUnif[0]) continue;  //since some of the CDFs have y cutoffs
				dny *= taupairswc-1.; // gives the bin number in distribution
		 		iy = (int)floor(dny);
				if(iy>iymax){
					iymax = iy;
				}
			}
		  	Bjyem2 = iymax>0? taupairsYref[iymax]: eps;
			//end of pairs
			iymax = -99;			
			//photonuclear:
			// this is a fit to the mean number of showers vs. energy
			Ntaupn = pow(10.,1.82*(1.-exp(-2.13*(log10(Etau)-17.95)))-1.);
			//printf("Etau=%e Ntaupn= %e\n", Etau, Ntaupn);
			ipn = (int)lrint(kpoi(Ntaupn)); // Poisson on the number of showers in this process
			for(i=0;i<ipn;i++){
				dny = ( (dr48=drand48())==1.0 ? dr48-eps : dr48 );
				if(dny<taupnUnif[0]) continue;  //since some of the CDFs have y cutoffs
				dny *= taupnwc-1.; // gives the bin number in distribution
		 		iy = (int)floor(dny);
				if(iy>iymax){
					iymax = iy;
				}
			}
			Bjypn = iymax>0? taupnYref[iymax]: eps;	
			// end of photonuclears
			//printf("Bjyem1=%e  Bjyem2=%e  Bjypn=%e  Emu=%e\n", Bjyem1, Bjyem2,Bjypn,Etau);			
			Esh.em = Bjyem1>Bjyem2 ? Etau*Bjyem1: Etau*Bjyem2;  
			Esh.had =  Enu*(Bjy+eps) > Bjypn*Etau ? Enu*(Bjy+eps) : Bjypn*Etau;
			if(Esh.had >= Esh.em){
				Esh.em = eps*Enu;
				}else{
				Esh.had = eps*Enu;
				}			
			
		    }else{     // neutral current
			Esh.em =   1.e8;  
			Esh.had = (Bjy)*Enu;
		     }
		}

/*printf(">>tau event: Bjy=%f Bjy1=%f Yref=%f Etau=%f Ltau=%f\n", Bjy, Bjy1, Yref[iy], Etau, Ltau); */


	return(Esh);

}
//-------------------------------------------------------------------
void set_shower_dist(showerE Es, int nut, double depth) /*sets up shower angle function cumulative array*/
{
  double gam,sham[NFREQ][NA+1],tmp[NFREQ][NA+1],nlocal,nu,dnu,bw;    
  double thetaC,nu0,nubar[NFREQ],Elpm,deltathet1[NFREQ],deltathet2[NFREQ],Efac2,gamrad;
  static char *cp = {"e"};
  double s0,s1,slope,E1e,E1h,fEh, fEe,trlh,trle,localdensity,get_firn_dens();
  int i,j,k,i0,i1;
  showerE E;
  double epss = 1.e-10;

	E = Es;
	
	localdensity = get_firn_dens(100.*depth);  // depth in meters originally, convert to cm for this call
	//nlocal = get_n_ice(depth);
	
	Elpm = 2.e15*0.90/localdensity;  /* for ice*/
	
	// thetaC = acos(1./nlocal);  this for later when ray tracing is ready
	thetaC = acos(1./Nref);
	nu0 = 500.e6;  /* this is different than the coherence NU0 */
	for(j=0;j<NFREQ;j++){
	   nubar[j] = NUBAR[j];
	}
	
	/* from Alvarez-Muniz & Zas, Phys Lett B 411 (1997) we have
	 Efield(theta) = Efield(thetaC)*exp(-ln2*((theta-thetaC)/sigma)**2)
	 where sigma = 2.7deg * nu0/nu * (Elpm/(0.14*Esh+Elpm))**0.3 for Esh>1PeV
	 for electron showers  IN ICE */
	 
	/* from Alvarez-Muniz & Zas, Phys Lett B 434 (1998) we have
	 sigma = 1deg * nu0/nu * (1.71 - 0.0121*log10(E/1TeV) ) for 0.1PeV <E<100PeV
	 sigma = 1deg * nu0/nu * (4.23 - 0.785*log10(E/1TeV)+0.055* (log10(E/1TeV))**2) for E>100PeV
	 for hadronic showers in ICE*/

	E1h = log10(E.had/1.e12); /* now in log units TeV */
	E1h = (E1h > 7.0)? 7.0: E1h;  // kluge to avoid divergence of parameterization
				      // above 1e19 eV.
	E1e = log10(E.em/1.e12);

	fEh = -1.27e-2-4.76e-2*(E1h+3)-2.07e-3*(E1h+3)*(E1h+3)+0.52*sqrt(E1h+3);
	trlh = 6.25*fEh*E.had/1.e9 *0.90/localdensity;  /* proj. tracklength in m */

	fEe = -1.27e-2-4.76e-2*(E1e+3)-2.07e-3*(E1e+3)*(E1e+3)+0.52*sqrt(E1e+3);
	trle = 6.25*fEe*E.em/1.e9 *0.90/localdensity;  /* proj. tracklength in m */
	

	/* these are the std. deviations for gaussian Cherenkov angle width */

       for(j=0;j<NFREQ;j++){         // for each frequency band
       	    deltathet1[j] = 0.;
	    deltathet2[j] = 0.;
	    dnu = BW[j]/200.0;  // interval for new integration 10/13/2005
	    bw = 0.0;
	    
	if(nut==0){              // electron neutrinos

	  // first width is for EM shower, 2nd for hadronic
	  // added integration over nu for all flavors--10/13/2005 PG
	  
	   for(nu=NUBAR[j]-BW[j]/2.;nu<=NUBAR[j]+BW[j]/2.;nu+= dnu){
	     deltathet1[j] += (2.7*PI/180.0)*nu0/nu * 
	   	pow(Elpm/(0.14*E.em+Elpm),0.3) * localdensity/0.9 * dnu;
	     deltathet2[j] += (1.0*PI/180.0)*nu0/nu * 
	   	(4.23-0.785*E1h+5.5e-2*E1h*E1h) * localdensity/0.9 * dnu;
		bw += dnu;
	     }
	    deltathet1[j] /= bw;
	    deltathet2[j] /= bw;
	    
	   }else if(nut == 1){  // muon neutrinos

	   //deltathet1[j] = (1.0*PI/180.0)*nu0/nubar[j] * (4.23-0.785*E1e+5.5e-2*E1e*E1e);
	   for(nu=(NUBAR[j]-BW[j]/2.);nu<=(NUBAR[j]+BW[j]/2.);nu+= dnu){
	     deltathet1[j] += (2.7*PI/180.0)*nu0/nu * 
	   	pow(Elpm/(0.14*E.em+Elpm),0.3) * localdensity/0.9 * dnu;
	     deltathet2[j] += (1.0*PI/180.0)*nu0/nu * 
	   	(4.23-0.785*E1h+5.5e-2*E1h*E1h) * localdensity/0.9 * dnu;
		bw += dnu;
	     }
	    deltathet1[j] /= bw;
	    deltathet2[j] /= bw; 
	    
	   }else if(nut == 2){  // tau neutrinos

	   //deltathet1[j] = (1.0*PI/180.0)*nu0/nubar[j] * (4.23-0.785*E1e+5.5e-2*E1e*E1e);
	   for(nu=NUBAR[j]-BW[j]/2.;nu<=NUBAR[j]+BW[j]/2.;nu+= dnu){	   
	   deltathet1[j] += (2.7*PI/180.0)*nu0/nu * 
	   	pow(Elpm/(0.14*E.em+Elpm),0.3) * localdensity/0.9 * dnu;
	   deltathet2[j] += (1.0*PI/180.0)*nu0/nu * 
	   	(4.23-0.785*E1h+5.5e-2*E1h*E1h) * localdensity/0.9 * dnu;
		bw += dnu;
	     }
	    deltathet1[j] /= bw;
	    deltathet2[j] /= bw;
	    
	  }
        }
	
	//update the global variables for the angular width
	for(j=0;j<NFREQ;j++){
	   Dtheta[j].em = deltathet1[j];
	   Dtheta[j].had = deltathet2[j];
//	  fprintf(stderr,"f=%e nut=%d  dthet1= %e  dthet2= %e\n",nubar[j],nut,  deltathet1[j], deltathet2[j]);
	   }

      for(j=0;j<NFREQ;j++){
       for(i=0;i<=NA;i++){
	sham[j][i] = 0.0;
	Sham[j][i] = 0.0;
	}	
       }

     for(j=0;j<NFREQ;j++){
       for(i=0;i<=NA;i++){
	gam = (((double)i)/DNA) *PI;    /* polar angle in rad, 0-PI */
    /* this uses Jaime's function, scaled for field */
	sham[j][i] = sin(gam)*E.em*exp(-0.693*pow((gam-thetaC)/deltathet1[j],2.0))
	           + sin(gam)*E.had*exp(-0.693*pow((gam-thetaC)/deltathet2[j],2.0));

	if(i>0) for(k=i-1;k>=0;k--)Sham[j][i] += sham[j][k]; // cumulate it
	  }
        }
	  // this normalizes the cumulative distribution:
      for(j=0;j<NFREQ;j++){
	for(i=0;i<=NA;i++){
	Sham[j][i] /= Sham[j][NA];
	tmp[j][i] = 0.0;
//	fprintf(stdout,"%d  %f\n", i, Sham[i]);
	}
       }
/* now reflect around y=x to create inverse CDF*/
      for(j=0;j<NFREQ;j++){
	for(i=0;i<=NA;i++){
	   tmp[j][(int)nearint(Sham[j][i]*DNA)] = (double)i/DNA;
	}
       }
	/* this vector will have some holes because the map is not 1-1 */
	/* so we need to linearly interpolate over the holes : */
        /* linear interpolation: */
      for(j=0;j<NFREQ;j++){
	s0 = 0.0;  i0 = 0; 
	for(i=1;i<NA;i++){
		if(tmp[j][i] == 0.0){
			continue;
		}else{
			s1= tmp[j][i];
			i1 = i;
			slope = (s1-s0)/(double)(i1-i0);
			for(k=i0+1;k<i1;k++){
				tmp[j][k] = tmp[j][i0]+slope*(double)(k-i0);
			}
		}
		i0=i1;
		s0 = s1;
	} 
	tmp[j][NA] = Sham[j][NA];
      }
/* end of linear interpolation */

      for(j=0;j<NFREQ;j++){
	for(i=1;i<NA;i++){
		Sham[j][i] = tmp[j][i]; 
	}
       }

//	for(i=0;i<=NA;i++)printf("%d %e\n",i,Sham[0][i]);

}
//-------------------------------------------------------------------


double get_upang()  // returns angle that is flat is cosine distribution (solid angle)
{
  double drand48();
  return acos(drand48());
}

//-------------------------------------------------------------------
double get_nadir_angle()  // returns angle that is flat in cosine distribution (solid angle)
			  // distributed only over first quadrant: 0=nadir, pi/2=horizon
{
  double drand48();
  return acos(drand48());
}

//-------------------------------------------------------------------

/* this for a more realistic earth parametrization */
double get_shower_upangle(L,d)
double L,d;  /* should be in m  */
{
  double ang,upang,drand48();
	if(L < 0.53e8){
		ang = PI/2.0 - (0.832522e8-L)/0.53e8;
		} else {
		ang = PI/2.0 - 0.56*acos(L/1.1e8);
		}
	upang = ang- d/L;  /* corrects for depth of interaction */
	return(upang);
}  
//-------------------------------------------------------------------

double get_dens(h) // not used currently
double h;
{
double p0,k,z,p;

	p0 = 1.27;
	k = 0.121;	
	z = 100.*h;
	return p = p0 + k*log(z+1.0);
}

//-------------------------------------------------------------------
double expdev()
{
	double drand48();
	return (-log(drand48()+1.e-12));
}
//-------------------------------------------------------------------
double kpoi(double m)          /* m=mean; returns kran according to poisson dist */
{
double drand48();
double z = 0.0,rnd, log();
double kran = 0.;
int k;

        for(k = 1; k <= 100 ; k++) {
                rnd  = drand48();
                z += (-log( rnd ))/m ;
                if (z > 1.0 )
                        break;
                kran += 1.0 ;
        }

        return(kran);
}

//-------------------------------------------------------------------

void set_ydist()  // Bjorken y distribution, not used currently
{
	int i=0,maxi;
	FILE *fy, *fopen();
	extern double Unif[], Yref[];

	
	fy = fopen("ymc.dat", "r");

	for(i=0;i<NNY;i++){
	   fscanf(fy, "%lf%lf", Unif+i,Yref+i);
			}
	fclose(fy);

/*	for(i=0;i<150;i++)printf("%f\n", Yref[i]); */
}

//-----------------------------------------------------------------
// Uses von Neuman method to get y values distributed according
// to a fitted parameterization  7/1/2005
// fixed a minor bug 8/1/2007 =--PG
// verfied mean of output distribution --PG 8/1/2007
double get_yvalue()
{
	double y,W,fy,fymax;
	
	fymax = (1.e10/0.301e-31*1e-38 * pow(10.,-1.92));  //fit at 1e19 eV

	W = fymax;
	fy = 0.0;
	while(W==W){
		W = fymax*(1.0-drand48());
		y = drand48();
		fy =  1.e10/0.301e-31*1e-38 * pow(10.,(-2.92 + exp(-y/0.0124) - 2.56*y + 1.46*y*y));
		if(W<fy) break;
	}
	return(y);
}
//-------------------------------------------------------------------
void set_mu_ydist(char *energy)  // read in CDFs for mu lepton brem, pairs, photonuclears
{
	int i=0,maxi;
	FILE *fmubrem, *fmupair,*fmupn, *fopen();
	extern double mubremUnif[],mupairsUnif[], mupnUnif[];
	extern double mubremYref[],mupairsYref[],mupnYref[];
	char muPN[80];
	
	
	fmubrem = fopen("mubremsCDF.dat", "r");
	fmupair = fopen("mupairsCDF.dat", "r");
	sprintf(muPN,"mupn%sCDF.dat",energy);
	//printf("%s\n", muPN); fflush(stdout);
	fmupn = fopen(muPN, "r");
	
	fscanf(fmubrem, "%lf%d", &Nmubrem,&mubremwc);
	fscanf(fmupair, "%lf%d",&Nmupairs,&mupairswc);
	fscanf(fmupn, "%lf%d", &Nmupn,&mupnwc);	
	

	
	for(i=0;i<mubremwc-1;i++)
	   fscanf(fmubrem, "%lf%lf", mubremUnif+i,mubremYref+i);
	for(i=0;i<mupairswc-1;i++)
	   fscanf(fmupair, "%lf%lf", mupairsUnif+i,mupairsYref+i);
	for(i=0;i<mupnwc-1;i++)
	   fscanf(fmupn, "%lf%lf", mupnUnif+i,mupnYref+i);	
	   		
	/*printf("&Nmubrem,&mubremwc = %e  %d\n",Nmubrem,mubremwc);
	printf("&Nmupairs,&mupairswc = %e  %d\n",Nmupairs,mupairswc);
	printf("&Nmupn,&mupnwc = %e  %d\n",Nmupn,mupnwc);*/	
			
	fclose(fmubrem);
	fclose(fmupair);
	fclose(fmupn);
/*	for(i=0;i<150;i++)printf("%f\n", Yref[i]); */
}

//-------------------------------------------------------------------
void set_tau_ydist(char *energy)  // read in CDFs for tau lepton brem, pairs, photonuclears
{
	int i=0,maxi;
	FILE *ftaubrem,*ftaupair, *ftaupn,*fopen();
	extern double taubremUnif[],taupairsUnif[], taupnUnif[];
	extern double taubremYref[],taupairsYref[],taupnYref[];
	char tauPN[80];
	
	
	ftaubrem = fopen("taubremsCDF.dat", "r");
	ftaupair = fopen("taupairsCDF.dat", "r");
	sprintf(tauPN,"taupn%sCDF.dat",energy);
	ftaupn = fopen(tauPN, "r");
	

	fscanf(ftaubrem, "%lf%d", &Ntaubrem,&taubremwc);
	fscanf(ftaupair, "%lf%d",&Ntaupairs,&taupairswc );
	fscanf(ftaupn, "%lf%d", &Ntaupn,&taupnwc);	
	/*printf("&Ntaubrem,&taubremwc = %e  %d\n",Ntaubrem,taubremwc);
	printf("&Ntaupairs,&taupairswc = %e  %d\n",Ntaupairs,taupairswc);
	printf("&Ntaupn,&taupnwc = %e  %d\n",Ntaupn,taupnwc);*/

	for(i=0;i<taubremwc-1;i++)
	   fscanf(ftaubrem, "%lf%lf", taubremUnif+i,taubremYref+i);
	for(i=0;i<taupairswc-1;i++)
	   fscanf(ftaupair, "%lf%lf", taupairsUnif+i,taupairsYref+i);
	for(i=0;i<taupnwc-1;i++)
	   fscanf(ftaupn, "%lf%lf", taupnUnif+i,taupnYref+i);			
		
		
	fclose(ftaubrem);
	fclose(ftaupair);
	fclose(ftaupn);
/*	for(i=0;i<150;i++)printf("%f\n", Yref[i]); */
}

//-------------------------------------------------------------------
double get_dist()   /* uniform distance distribution up to MAXDIST */
/* just throw the dice for x and y, calculate radius */
{
   double d, drand48(),dx,dy;  

	d = 2*MAXDIST;
    while(d>MAXDIST){
	dx = drand48()*MAXDIST;
	dy = drand48()*MAXDIST;
	d = hypot(dx,dy);
	}
	return(d);
}



//-------------------------------------------------------------------
/* field strength from Alvarez-muniz et al PRD 62,63001 (2000) */
/* Ec = cascade energy, D is distance in m, nu0 is rolloff freq,
numin, numax limits of bandpass in MHz, result  */

double get_field(Ec,D,nu0,numin,numax,thet,nut,depth,dnufac)
showerE Ec; double D,nu0,numin,numax,thet,depth,dnufac;
int nut;
{
	double EE,dfac,nu,Fac,field,dnu,nn,ange,angmu,angtau,Elfac,Emfac,Etfac;
	double fE,trl,Elpm,E1,thetaC,E,deltathet,bw,dthet,nn0, Es[2];
	double nice,localdensity,Firnfactor,n_ice,get_n_ice(),get_firn_dens();
	int i;


	bw = numax-numin;
	dnu = (numax-numin)/dnufac;
	
	Es[0] = Ec.had;  // hadronic part
	Es[1] = Ec.em;   // electromagnetic part
	field = 0.0;
	
	// these four lines added 2/27/2005
	// accounts for extra path length, lower refractive index in firn
	n_ice = get_n_ice(depth);
	localdensity = get_firn_dens(100.*depth); // convert meters to cm for firn call
	
	Firnfactor = 0.90/localdensity* sqrt(cos(1./get_n_ice(depth)))/sqrt(cos(1./1.79));
	//printf("n_ice= %f localdensity= %f Firnfactor = %e\n",n_ice,localdensity, Firnfactor);		
	
  for(nu=numin;nu<numax;nu+=dnu){	
   for(i=0;i<2;i++){


	EE = (Es[i]/1e12);
	dfac = 1./D;


  /* Fac contains various terms including the differential */
	// Fac = 2.53e-7*dfac*EE*Firnfactor*dnu; origal ZHS; below is a
	// modified empirical formula which rolls the energy off
	// slowly above 1e23 eV to account for energy saturation
	
	Fac = dfac*Firnfactor*2.53e-7*(EE*sqrt((1.-exp(-1.e23/Es[i]))) + 
		sqrt(EE*(1.-exp(-Es[i]/1.e23)))) ;
		
	//Fac = 2.53e-7*dfac*EE*Firnfactor;

	E = Es[i];

	Elpm = 2.e15 * 0.9/localdensity;  /* electromagnetic, for ice*/

	thetaC = acos(1./Nref);

	E1 = log10(E/1.e12); /* now in log units TeV, USED ONLY FOR ANGLES */
	if(nut > 0 && E1 > 7.0) E1 = 7.0;  // kluge to fix divergence of ZHS param.!!

	fE = -1.27e-2-4.76e-2*(E1+3)-2.07e-3*(E1+3)*(E1+3)+0.52*sqrt(E1+3);
	trl = 6.25*fE*Es[i]/1.e9;  /* proj. tracklength in m */

	dthet = thet-thetaC;
	
	/* from Alvarez-Muniz & Zas, Phys Lett B 411 (1997) we have
	 Efield(theta) = Efield(thetaC)*exp(-ln2*((theta-thetaC)/sigma)**2)
	 where sigma = 2.7deg * nu0/nu * (Elpm/(0.14*Esh+Elpm))**0.3 for Esh>1PeV
	 for electron showers  IN ICE */
	 
	/* from Alvarez-Muniz & Zas, Phys Lett B 434 (1998) we have
	 sigma = 1deg * nu0/nu * (1.71 - 0.0121*log10(E/1TeV) ) for 0.1PeV <E<100PeV
	 sigma = 1deg * nu0/nu * (4.23 - 0.785*log10(E/1TeV)+0.055* (log10(E/1TeV))**2) for E>100PeV
	 for hadronic showers in ICE*/
	

 	nn = (nu/nu0);
	nn0 = (nu/500.0); /* for the angles not the coherence */

        if(nut == 0){  /* electron neutrino; 2.6 factor for decrease in shower length, etc.*/

	   ange = (2.7*PI/180.0);
	   Elfac =  pow(Elpm/(0.14*Es[i]+Elpm),0.3);
	   deltathet = ange/nn0 * Elfac*localdensity/0.90;
	   if(i==0){                          // hadronic part of e- shower
	  	ange =  (1.0*PI/180.0);
	   	Emfac = (4.23-0.785*E1+5.5e-2*E1*E1);
 	   	deltathet = ange/nn0 * Emfac *localdensity/0.90;
		}

	   }else if(nut == 1){  /* muon neutrino */

	   angmu =  (1.0*PI/180.0);
	   Emfac = (4.23-0.785*E1+5.5e-2*E1*E1);
 	   deltathet = angmu/nn0 * Emfac *localdensity/0.90;
		
	   }else if(nut >= 2){ // tau neutrino, similar to muon

	   angtau =  (1.0*PI/180.0);
	   Etfac = (4.23-0.785*E1+5.5e-2*E1*E1);
 	   deltathet = angtau/nn0 * Etfac *localdensity/0.90;
	}

	field += Fac*nn*(1./(1.+pow(nn,1.44)))
		 *exp(-0.693*pow(dthet/deltathet,2.0)) * dnu;

       }
     }
	return(field/bw);  /* V/m/MHz */
}
//-------------------------------------------------------------------

double get_n_ice(double depthm) // depth in meters
{
	return (1.325 + 0.463*(1.0-exp(-.0140*depthm)));
}
//-------------------------------------------------------------------
void set_thetout()
{
FILE *fthet, *fopen();
float tmp1, att;
extern float *Thetout, *Transout;
extern int NFIRN;
int i;

	fthet= fopen("firntrace1.dat", "r");
	for(i=0; i<NFIRN;i++){
		fscanf(fthet, "%f  %f  %f  %f", 
		&tmp1 , Thetout+i, Transout+i, &att);
//		Transout[i] *= att;
//fprintf(stderr,"%f  %f\n", Thetout[i], Transout[i]);
		}
	fclose(fthet);
}
//-------------------------------------------------------------------

double get_thetout(tfirn)
double tfirn;
{
	int tindex;
	tindex = (int)floor((double)NFIRN*(tfirn/TIR));
	
	if(tindex>=NFIRN){
		 return(-99.);
	}else{
	return( Thetout[tindex] );
	}
}
//-------------------------------------------------------------------
double get_transout(tfirn)
double tfirn;
{
	int tindex;
	tindex = (int)floor((double)NFIRN*(tfirn/TIR));
	if(tindex>=NFIRN){
		 return(0.0);
	}else{
	return( Transout[tindex] );
	}
}

//----------------------------------------------
double get_firn_dens(double depth)  // depth in cm
{
// new routine added 2/27/2005 for firn density
// fit from gnuplot using rice data
//f(x) = 0.68 + a1*(1.0-exp(-b1*x));
//a1 = .22; b1 = .0131
// 
	double depth1;
	
	depth1 = depth/100.0;
	return (0.68 + 0.22*(1.0-exp(-0.0131*depth1)));
}


//-------------------------------------------------------------------
/* this is a polynomial fit to the change in RF loss coefficient at 300 MHz
vs. depth for data from Bogorodsky's book */
/* the average from depth to 0 is computed and returned*/
double get_avg_alpha(double alpha, double depth)  // depth in m
{
	double alpha0_300=6.6e-4, coeff=6.24e-17;
	double depth1,losstot,avgloss,dD;
	dD = depth/100.;
	losstot = 1.0;
	for(depth1=depth;depth1>=0.0; depth1 -= dD){
		losstot *= exp(-(alpha + coeff*pow(depth1,4.0))*dD);
		}
	avgloss = -log(losstot)/depth;
	return (avgloss);  //loss factor per meter
}

//------------------------------------------------------------------

char *get_E(double Enu)
{
double lEnu;
int iEnu;
static char *E[] = {"3e16" "1e17","3e17","1e18","3e18","1e19","3e19","1e20","3e20","1e21","3e21","1e22"};
static char E0[] = {"ESS"};
	if(Enu==0.0) return(E0);
	lEnu = log10(Enu);
	iEnu = (int)lrint(2.*lEnu - 33.5);
	return(E[iEnu]);
}



