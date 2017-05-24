/*
anitaVeffE10R: monte carlo radio emission air showers detected by 
antarctic balloon, using field strength array 

gcc anita3VeffE10R.c -o ~/bin/anita3VeffE10R -lm -O2

This version for 10 frequency bands
2/15/2004  -PG

modified 9/2004 for updates to anitamcE.c --PG

major updates in March/April 2005
Changed integration to pure MC for final section.

11/21/05:  added hooks for bottom dropring array --PG

12/12/2007:  updates for ZHS corrections, etc.

September 2009 major changes to adapt for ANITA-2

11/5/2010:  take this version for ANITA-3, lower the threshold

11/17/2010:  this version uses 10 frequency bands, prints out
 	
 10/15/2013, modified hit printout --PG
*/
	
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <values.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>


#define PI    M_PI


#define horb   36.0     /* orbit altitude, km, above 3km high ice sheet*/
#define Re     6359.755  /* polar radius, km, including 3 km ice sheet */
#define RH     (sqrt((Re+horb)*(Re+horb)-Re*Re))  // distance to horizon from payload */
#define ThetH  (atan2(RH,Re))

#define C0     2.9979e8

double  MAXANG = 0.1;

double MAXDIST= 690.00;   /* km */


double NBTHR;
int NTRIG0 = 2;    /* basic ring coincidence level */
double angperpix;

#define NEL  128 		/* size of angular intensity array */
#define NAZ  256 		/* size of angular intensity array */
int NI;
#define IRANGE   600.0000	/* km initial range for input image scaling */

#define Navo  6.022e23
#define MINANG  0.00    /* minimum nadir angle to use */

#define NUPANGLE  100   /* used for input elevation angle acceptance array */

#define NFREQ 10

static double NUBAR[NFREQ]= {250.e6,350e6,450e6,550e6,650e6,750e6,850e6,950e6,1050e6,1150e6  }; // center frequencies, unweighted
static double BW[NFREQ] = {100e6,100e6,100e6,100e6,100e6,100e6,100e6,100e6,100e6,100e6};  // 3 dB bandwidths
static double dnufac[NFREQ] = {5,5,5,5,5,5,5,5,5,5}; // number of frequency integration intervals per band
//static double GAIN0[NFREQ] = { 8.85 ,  9.15 ,  9.46 ,  9.81 , 10.21 , 10.61 , 10.92 , 11.02 , 10.90 , 10.65 };  // gains at NUBAR[i] OLD
static double GAIN0V[NFREQ] = { 7.0 ,  8.0 ,  8.5 ,  11.0 , 9.2 , 9.5 , 10.5 , 11.0 , 13.6 , 14.5 };  // gains at NUBAR[i]
static double GAIN0H[NFREQ] = { 6.7 ,  7.0 ,  7.7 ,  11.0 , 9.0 , 9.0 , 10.0 , 9.0 , 12.0 , 13.0 };  // gains at NUBAR[i]
// NOTE: although we have two different gains for V and Hpol, an average beamwidth is used for both -- this should be fixed later!!
//static double GAIN22[NFREQ] = {5.5 ,6.5, 6.,6.};  // gains at 22.5 degrees off boresight


#define NANT 27   // out of 40 total

/* antenna elevation and azimuth (deg.) of boresight, for 27 antennas */
static double ANTEL[NANT] = {-10.,-10.,-10.,-10.,-10.,\
                                -10.,-10.,-10.,-10.,\
			   -10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,\
			   -10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.}; // last line for 2nd lower tier

static double ANTAZ[NANT] = {-90.,-67.5,-45.,-22.5,0.,22.5,45.,67.5,90.,\
			  -90.,-67.5,-45.,-22.5,0.,22.5,45.,67.5,90.,\
			  -90.,-67.5,-45.,-22.5,0.,22.5,45.,67.5,90.}; // last line for 2nd lower tier

/*        X X X X X  5 of top 8
           X X X X   4 of next 8   
	  XXXXXXXXX  9 of 1st lower tier 16, adjacent 
	  XXXXXXXXX  9 of 2nd lower tier 16, adjacent 
	 
This represents a half array */
static int PHISECTOR[NANT] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9};

/* these are ANITA-2 sytem temperatures assuming 90K system temperature, 150-240K ice*/
/*  actual measured average values correspond to about 270K  in quiet conditions */
static double TSYSANT[NANT]={270.,270.,270.,270.,270.,\
                                270.,270.,270.,270.,\
			   270.,270.,270.,270.,270.,270.,270.,270.,270.,\
			   270.,270.,270.,270.,270.,270.,270.,270.,270.}; // last line for 2nd lower tier

/* PIDgoals for all but Nadirs were 14MHz, 8MHZ, 8MHZ, 1MHz */
/* for the nadirs they were 200KHz,200KHz,200KHz,40KHz */ 
//static double BTHR[NFREQ] = { 1.5, 1.8, 1.8, 2.45};  // original threshold in sigmas for 4 bands
//static double BTHR[NFREQ] = { 1.5, 1.8, 1.8, 2.92};  // threshold in sigmas for 4 bands, 2.92 based on 100Hz raw rate

static double BTHR = 2.4;  // threshold in voltage sigmas for full band
static double BWtot= 1000.;  // MHz

/*>>>>>>>>>>>>>>>>>>  start main <<<<<<<<<<<<<<<<<<<<<<<<<*/
 	
main(int argc,char **argv)
{
   int Nevt, NHIT1,NHIT2,NHIT3,N2trig0,N2trig1,N2trig2,N2trig3, Nimg,Coinc,i,j,k,m,n,q;
   int kp,km,mp,mm,np,nm, dropringON, Goodpixel;
   double kpslope,kmslope,kangle,EjfV[NFREQ],EjfH[NFREQ];
   int iant,jf,Nrange, iangle,ii,fim,npix,p, NADIRHIT, TOTHIT,*Nut, Ievt, NEVT, NuImg;
   double Esh,angle=0,Nev,R,angfac, dthet, theta,thetamax, Thetevt, Sevt, Revt,Phinu;
   double angle0=.01,dist,E,VPMM,Volts_per_m_per_MHz;
   double get_shower_angle(),get_shower_E(),get_dist(),exp();
   double time,rate,Nhit,beamfac,polfac=1.0;
   double hitchan, chanv,noise, gasdev();
   double thetae, Rx,Ry,beta,thetaf, MAXarea, Km3sr,rbeta;
   double dbeta,dphi,Domega,phi,range,depthmax,attenlength;
   double Dr,x,ang,angx,angy,rho, V0, OM0,meanangle0,sigangle0;
   double Vn[NANT],AeffV[NFREQ],heffV[NFREQ],AeffH[NFREQ],heffH[NFREQ],lam[NFREQ],thresh[NANT];
   double nhitV[NANT],nhitH[NANT],V_rcvrV[NFREQ][NANT],V_rcvrH[NFREQ][NANT],Vobs_V[NANT],Vobs_H[NANT];
   int    L1_V[NANT], L1_H[NANT], L2[NANT], L3, NeutrinosUsed, *IMGHIT;
   float  *Img[NFREQ], *PImg[NFREQ];
   double beta0,beta1,rang0,rang1,theta0,theta1,get_theta();
   double Km2sr,km2percm2,sigma_nu(),Enu,fluxlim,maxdepth;
   double Tsys,Z=50.0, Z0=376.99, K=1.38e-23, km3sr, beamfacV,beamfacH;
   double X1, Thet, Dtheta,nadangle,xel,X2,X0,gasdev(),TrueSNR_V[NANT],TrueSNR_H[NANT];
   double Vobs_lcp, Vobs_rcp, AtmpV, AtmpH;
   double *Ftrials,*EE1,*EEem,*EEhad,*LL1,*dd1,*Up1,*Enu0;
   double gain,Gain,f0,qrange, rfthet, rfaz, dnu[NFREQ],nu, get_rician();
   char   Evsum[80], Evsumout[80], Imbasename[NFREQ][80],PImbasename[NFREQ][80];
   double MHIT, Ntrials, anglesum, angle2sum, meanangle, sigangle, nangle,dropringfac;
   double Upangle0[NUPANGLE+1],Upangle[NUPANGLE+1], NUPperbin,uprat, omega;
   double minangle,maxangle,anglerange,dangle, beampower(), INtrials, Ntrials1, cross_section_factor;

	FILE *Flog, *Fhit,  *fopen();
	int Fimage[NFREQ], Pimage[NFREQ], ESS;
  struct timeval seed_time;

/*___________________________end of declarations__________________________*/

USAGE:      

	if(argc < 8){
	fprintf(stderr,
"usage: anita3VeffE10R [Nimg][NEVENT][maxdepth(m)][image file basename][log file][Enu][output hitfile][xsection factor]\n");
		exit(0);
		}
		
	// a large prime seed; modified October 2005 --PG
	gettimeofday(&seed_time, (struct timezone*)0);
	//srand48(1299811);
	srand48(seed_time.tv_usec);
	fprintf(stdout,"#Seed value= %d\n", seed_time.tv_usec);
		

	Nimg = atoi(argv[1]);  /* images in archive file*/
	NEVT = atoi(argv[2]);  /* number of trial interactions*/
	maxdepth = atof(argv[3])/1000.; /* maximum depth in km used for images */
	Flog = fopen(argv[5],"r");
	Enu = atof(argv[6]);
	Fhit = fopen(argv[7], "w");
	cross_section_factor = atof(argv[8]);
	dropringON = 1;  // always ON for this one
	
	ESS=0;
	if(Enu==0.0) {
	      ESS=1;
	      Enu= 1e18; // approx. average ESS energy
	    fprintf(stderr,"ESS= %d\n", ESS); fflush(stdout);
	  }
	
	for(i=0;i<NFREQ;i++){
		sprintf(Imbasename[i],"%s%d.img", argv[4],(int)(NUBAR[i]/1.e6));
		sprintf(PImbasename[i],"%s%dP.img", argv[4],(int)(NUBAR[i]/1.e6));
		if( (Fimage[i] = open(Imbasename[i],O_RDONLY)) != -1 ){
			fprintf(stderr,"Image file %s opened...\n", Imbasename[i]);
		    }else{ 
			fprintf(stderr,"Error opening %s\n", Imbasename[i]);
			exit(0);
		   }
		if( (Pimage[i] = open(PImbasename[i],O_RDONLY)) != -1 ){
			fprintf(stderr,"PImage file %s opened...\n", PImbasename[i]);
		    }else{ 
			fprintf(stderr,"Error opening %s\n", PImbasename[i]);
			exit(0);
		   }
		}

	angperpix = PI/(double)NAZ;

	MAXDIST = sqrt((Re+horb)*(Re+horb) - Re*Re);

	/* this counter array histograms the number of times an image contributes to
		the hit total, makes sure we have a decent number of contributors */
	IMGHIT = (int*)calloc(Nimg, sizeof(int));

	N2trig0 = 1;         // single antenna trigger
	N2trig1 = NTRIG0;    // coinc level within upper 2 rings
	
	/* polarization factor for linear-->Circular pol */
	polfac = 0.90;  /* 10/2/2008: includes factor of 0.9 for some Hpol content */
	
	/*  this estimate uses the mean frequency 
	 for(jf=0;jf<NFREQ;jf++){
	  lam[jf] = C0/NUBAR[jf];
	  Aeff[jf] = lam[jf]*lam[jf]/(4.*PI)*GAIN0[jf];
	  heff[jf] = 2.*sqrt(Z*Aeff[jf]/378.0);
	   //fprintf(stderr,"Aeff[%4.0f] =%e  heff[%4.0f] =%e\n", NUBAR[jf]/1e6,Aeff[jf],NUBAR[jf]/1e6,heff[jf]);
	 } */
	 
	 
	 /* this estimate uses the weighted average Aeff and heff */
	for(jf=0;jf<NFREQ;jf++){
	   dnu[jf]= BW[jf]/dnufac[jf];
	   AeffV[jf] = 0.0;
	   heffV[jf] = 0.0;
	   AeffH[jf] = 0.0;
	   heffH[jf] = 0.0;
	   for(nu=NUBAR[jf]-BW[jf]/2.;nu<= NUBAR[jf]+BW[jf]/2.; nu+=dnu[jf]){
	      lam[jf] = C0/NUBAR[jf];
	      AtmpV = lam[jf]*lam[jf]/(4.*PI)*GAIN0V[jf];
	      AeffV[jf] += AtmpV*dnu[jf];
	      heffV[jf] += 2.*sqrt(Z*AtmpV/Z0)*dnu[jf]; /* Z0=377ohms*/
	      AtmpH = lam[jf]*lam[jf]/(4.*PI)*GAIN0H[jf];
	      AeffH[jf] += AtmpH*dnu[jf];
	      heffH[jf] += 2.*sqrt(Z*AtmpH/Z0)*dnu[jf]; /* Z0=377ohms*/	      
	    }
	     AeffV[jf]/= BW[jf];
	     heffV[jf]/= BW[jf];
	     AeffH[jf]/= BW[jf];
	     heffH[jf]/= BW[jf];
	     //	 fprintf(stderr,"Aeff[%4.0f]=%e  heff[%4.0f]=%e\n", NUBAR[jf]/1e6,Aeff[jf],NUBAR[jf]/1e6,heff[jf]);
	 }
	 
	 

	Nev = (double)Nevt;
	rho = 0.90;   /* density of ice */

 	Nhit=0.0;
	NI = NEL*NAZ;
	
	// allocate the image & polarization array
	for(jf=0;jf<NFREQ;jf++){
		 Img[jf] = (float*)calloc(Nimg*NI,sizeof(float));
		 PImg[jf] = (float*)calloc(Nimg*NI,sizeof(float));
		}
	
	
   /* first read in image & polarization files */
	for(jf=0;jf<NFREQ;jf++){
		 read(Fimage[jf], Img[jf], Nimg*NI*sizeof(float));
		 read(Pimage[jf], PImg[jf], Nimg*NI*sizeof(float));
		}
  /*  the loaded image is 256 x 128 by Nimg in memory */

	Ftrials=(double*)calloc(Nimg,sizeof(double));
	EEem  = (double*)calloc(Nimg,sizeof(double));
	EEhad = (double*)calloc(Nimg,sizeof(double));
	Nut   = (int*)calloc(Nimg,sizeof(int));
	Up1   = (double*)calloc(Nimg,sizeof(double));
	LL1   = (double*)calloc(Nimg,sizeof(double));
	dd1   = (double*)calloc(Nimg,sizeof(double));
	Enu0   = (double*)calloc(Nimg,sizeof(double));
  
   /* and the event info from the summary file */
	for(i=0;i<Nimg;i++){
		
       		fscanf(Flog,"%lf %lf %lf %d %*f %lf %lf %*f %*f %lf %*d %lf\n", 
         	Ftrials+i,EEem+i,EEhad+i,Nut+i,Up1+i,LL1+i,dd1+i, Enu0+i);      		
		//printf("%e %e %e %d %e %e %e\n",Ftrials[i],EEem[i],EEhad[i],Nut[i],Up1[i],LL1[i],dd1[i]); 
    		//sleep(1);
	  }

	
	/* number of total trials: events tried here* (trials/image)*/
	Ntrials1 = ((double)(NEVT))*(Ftrials[Nimg-1]/(double)Nimg);
	
//fprintf(Fhit,"#TrialN  Ecas-em,eV Ecas-had,eV Nut(012/emt) depth,m  range,km  upang,d  nadangle,d  az,d  E1 .... E10(V/m/MHz),Pavg(deg), beta(deg), Enu0(EeV)\n");


/* --------------start the main MC sequence on cascade events HERE ----------- */

     Nhit = 0;
     Km3sr = 0.0;
     MHIT = 0.0; NADIRHIT=0; TOTHIT=0;
     NeutrinosUsed = 0;
     
     for(Ievt=0;Ievt<NEVT;Ievt++){
     
       /* first generate random sample variables */
       
       	while((Thetevt=acos(drand48())) > ThetH){  
	     // these trials just to place the interaction unformly within horizon
	      }	
	Revt = sqrt(2.*Re*Re*(1.-cos(Thetevt)));  // chord distance from nadir to event
	// range from payload to point above event at surface
	range = sqrt((Re+horb)*(Re+horb) + Re*Re - 2.*(Re+horb)*Re*cos(Thetevt));
	qrange = range;
	
	Sevt = Re*Thetevt;   // arc distance on surface

	Phinu = PI*drand48();  // azimuthal angle relative to boresight to interaction
	
	while( (NuImg = (int)(Nimg*drand48())) > Nimg){};  // which neutrino to use for this event
	
	/* these are some necessary geometric terms: */
	X1 = Re*sin(Thetevt);  // cylindrical radius to event from axis
	// nadir angle:
	nadangle = asin(X1/range);
	beta = PI/2.- nadangle - Thetevt ;  // elevation angle of payload from event
	
	xel = nadangle*180.0/PI;  // nadir angle in degrees

	rfthet = nadangle - PI/2.;// angle below horizon of RF arrival at payload

	rfaz = PI/8.*(2.*(drand48()-0.50));  // random +/- 22.5 degree offset relative
					      // to center antenna boresight

	// indexing variables:
	k = (int)(floor(beta/(PI/2.)*NEL));
	kp = k+1; km = k-1;
	j = (int)(floor(Phinu/PI*NAZ));
        m = k*NAZ + j;
	mp = kp*NAZ + j; mm = km*NAZ + j;
	n = NuImg*NI + m;
	np = NuImg*NI + mp; nm = NuImg*NI + mm;
	
	
	// this short section interpolates the image pixel values
	kpslope=kmslope=0;
	Goodpixel = 1;

	  for(jf=0;jf<NFREQ;jf++){
	   	kangle = ((double)k+0.5)*angperpix;
		if(Img[jf][n] == 0.0) Goodpixel = 0;
		kpslope = (double)(Img[jf][np]-Img[jf][n])/angperpix;
/* next line includes sine of polarization angle to get Vpol  --PG 9/17/09*/
//		EjfV[jf] = kpslope*(beta-kangle) + Img[jf][n]*sin(PImg[jf][n]); // for Vpol
//		EjfH[jf] = kpslope*(beta-kangle) + Img[jf][n]*fabs(cos(PImg[jf][n])); // for Hpol
		EjfV[jf] =  Img[jf][n]*sin(PImg[jf][n]); // for Vpol
		EjfH[jf] =  Img[jf][n]*fabs(cos(PImg[jf][n])); // for Hpol		
		}

	if(Goodpixel == 0) continue;  // go to the next event

//fprintf(stdout,"r= %5.1f Phinu= %5.2f NuImg= %d X1= %5.1f nad= %5.3f beta= %5.3f k= %d j=%d m=%d n=%d\n",
//                  range,Phinu*180./PI,NuImg,X1,nadangle*180./PI,beta*180./PI,k,j,m,n);
//		  fflush(stdout);

	   

	     NHIT1=NHIT2=NHIT3= 0;  // counters for hits within single trigger
	     for(iant=0;iant<NANT;iant++){
	     	nhitV[iant]=0; // antenna hit counters
		nhitH[iant] = 0;
		}

            /*----------------------------------------------------------------- */
	        /*>>>>>>>>>>>>>> antenna loop <<<<<<<<<<<<*/

//	fprintf(stdout,"el= %f  az= %f  Img[n]= %e qrange=%e  E=%e\n", 
//		beta,phi,Ejf[jf],qrange,V_rcvr[jf]); 
  	 /* generate the voltages for each of 18 antennas, 2 pols, 4 freq each, 
	     rings have 5 antennas * 4 freq * 2 pol = 40 channels*/
                for(iant=0;iant<NANT;iant++){
		     Vn[iant] = sqrt(K*TSYSANT[iant]*Z*BWtot*1.e6); /* noise voltage rms in the load */
	             thresh[iant] = BTHR*Vn[iant];
		     TrueSNR_V[iant]=TrueSNR_H[iant]=0.0;
		     for(jf=0;jf<NFREQ;jf++){
			/* this is the signal amplitude sum:  heff/2*Efield/R*df */
		    	beamfacV=sqrt(beampower(jf,PI*ANTAZ[iant]/180.+rfaz,rfthet-PI*ANTEL[iant]/180.));
		    	beamfacH=sqrt(beampower(jf,rfthet-PI*ANTEL[iant]/180.,PI*ANTAZ[iant]/180.+rfaz));
			V_rcvrV[jf][iant] = sqrt(Z/Z0*AeffV[jf]) * EjfV[jf]*IRANGE/qrange*BW[jf]/1.e6 * beamfacV;
 			TrueSNR_V[iant] += V_rcvrV[jf][iant];
			V_rcvrH[jf][iant] = sqrt(Z/Z0*AeffH[jf]) * EjfH[jf]*IRANGE/qrange*BW[jf]/1.e6 * beamfacH;
 			TrueSNR_H[iant] += V_rcvrH[jf][iant];
			}
			TrueSNR_V[iant] /= Vn[iant];
			TrueSNR_H[iant] /= Vn[iant];
			Vobs_V[iant] = Vn[iant]*get_rician(TrueSNR_V[iant]); // proper rician
			Vobs_H[iant] = Vn[iant]*get_rician(TrueSNR_H[iant]); // proper rician
 
			if(fabs(Vobs_V[iant]) >=thresh[iant] && TrueSNR_V[iant]>0.5 ) nhitV[iant]++;
			if(fabs(Vobs_H[iant]) >=thresh[iant] && TrueSNR_H[iant]>0.5 ) nhitH[iant]++;
//	fprintf(stdout,"el=%f az=%f jf= %d ia= %d\tElcp= %g Ercp= %g thr=%g\tHIT= %d\n",
//	beta,phi,jf,iant,Vobs_lcp,Vobs_rcp,thresh[iant], NHIT);		       
		  }   /* end of  antenna loop */

		    



  /* now count up the L1 triggers */
    	       for(iant=0;iant<NANT;iant++){
			L1_V[iant] = 0;
			L1_H[iant] = 0;
			if(nhitV[iant] > 0 ) L1_V[iant]=1; 
			if(nhitH[iant] > 0 ) L1_H[iant]=1; 
				}

	/* now get the L2s */
		for(iant=0;iant<NANT;iant++) L2[iant] = 0;
		for(iant=0;iant<8;iant++){  // L2s for top ring
			if( (L1_V[iant]>0 && L1_V[iant+1]>0) || (L1_H[iant]>0 && L1_H[iant+1]>0) ){ 
				L2[iant]=1;
				L2[iant+1]=1;
			}
		  }

		for(iant=9;iant<17;iant++){  // L2s for middle ring
			if( (L1_V[iant]>0 && L1_V[iant+1]>0) || (L1_H[iant]>0 && L1_H[iant+1]>0) ){ 
				L2[iant]=1;
				L2[iant+1]=1;
			}
		  }
		for(iant=18;iant<26;iant++){  // L2s for bottom ring
			if((L1_V[iant]>0 && L1_V[iant+1]>0) || (L1_H[iant]>0 && L1_H[iant+1]>0) ){ 
				L2[iant]=1;
				L2[iant+1]=1;
			}
		  }
 

 
/* now call a routine with the L1 pattern to determine whether we got an L2 and L3 */
	       if( (L3 = getL3(L2))==1 ) { 
			MHIT += 1.0;  /* this is the main L3 trigger counter, */
			IMGHIT[NuImg] = 1;
			Nhit++;
// 		double PolAvg=0;
// 	        for(int ip=0;ip<NFREQ;ip++){  
// 		   PolAvg += 180/PI*PImg[ip][n];
// 		}
// 		PolAvg /= NFREQ;
  fprintf(Fhit,"%.0f %e %e %d %6.2f %5.1f %5.3f %5.3f %5.3f %5.3f %5.3f\n", 
  	Ftrials[NuImg],EEem[NuImg],EEhad[NuImg],Nut[NuImg],dd1[NuImg],qrange,Up1[NuImg],nadangle*180./PI,(Phinu-PI/2.)*180./PI,beta*180./PI, Enu0[NuImg]); 
	  for(iant=0;iant<NANT/3;iant++){
	    for(jf=0;jf<NFREQ;jf++){
	      fprintf(Fhit,"   %.2e %.2e",V_rcvrV[jf][iant],V_rcvrH[jf][iant]); 
	    }
	    fprintf(Fhit,"  %.2e  %.2e  %4.3f  %4.3f\n",Vobs_V[iant], Vobs_H[iant],TrueSNR_V[iant],TrueSNR_H[iant] );
	  } 
	    fprintf(Fhit,"\n");	

	  }   // end of  L3 if loop

	} // end of Event loop
 
           // trial volume of ice below payload (formula for spherical cap):
	V0 = 4.*PI*Re*Re*sin(ThetH/2.)*sin(ThetH/2.)*maxdepth*rho;  
	Km3sr = (MHIT/Ntrials1)*V0*2.*PI; // 2pi is the solid angle of trial neutrino directions
					  // because we only try events pointing in the hemisphere
					  // toward the payload, not away from it...
	

	for(i=0;i<Nimg;i++) NeutrinosUsed += IMGHIT[i];

   /*%%%%%%%%%%%%%%%%%%%%% print out final results %%%%%%%%%%%%%%%%%%%%%%%%*/

	km2percm2 = 1e-10;
	double km_per_cm = 1.e-5;
	double Lintnu_km = 1./(Navo*sigma_nu(Enu,cross_section_factor))  * km_per_cm;
	Km2sr = (Km3sr)/ Lintnu_km ;


	fprintf(stdout,"%d different Neutrinos used out of %d hits\n", NeutrinosUsed, (int)MHIT);

	fprintf(stdout,"Nhit= %6.1f 	Veff= %e  km^3sr  Lint=  %7.1f  Aeff= %e km2sr\n", 
			Nhit, Km3sr, Lintnu_km, Km2sr);

	fluxlim = 1./Km2sr*km2percm2/(30.*86400.);
	fprintf(stdout,"90pctCL Flux limit in one month live: < %e  per cm2 per s per sr\n",
			2.3*fluxlim);
	fprintf(stdout,"90pctCL Flux limit in one month live: < %e  per km2 per month per sr\n",
			2.3*fluxlim*1e10*(30.*86400.));

	exit(0);

}  /* >>>>>>>>>>>>>>>>> end main <<<<<<<<<<<<<<<<<<<<<<<<<<<*/

 
// for range r to the payload, returns angle relative to center of earth
double get_theta(r)  
double r;
{

	double theta;

	theta = acos((Re*Re + (Re+horb)*(Re+horb)-r*r)/(2.*Re*(Re+horb)));

	return(theta);
}	



double sigma_nu(E,cross_section_factor)  // 2nd arg is multiplier for SM cross section
double E,cross_section_factor;
{
// cross section from Sarkar & Sarkar 2008
	double sigcc,signc,cc2nc_ratio,sigtotal;

	sigcc = 1.e-36 * exp ( 82.893 - 98.8*( pow( log(E/1.e9),-0.0964 ) ) );
	cc2nc_ratio = 2.39;
	signc = 1.e-36 * exp ( 82.893 - 98.8*( pow( log(E/1.e9),-0.0964 ) ) ) / cc2nc_ratio;
	sigtotal = (sigcc+signc)*cross_section_factor;

	return(sigtotal);
}	

/* this function is a gaussian fit to the antenna beam profile, based on seavey data */
/* --PG Feb. 2004 */
// corrections to angle conversions, March 2005  --PG


double beampower( int ifreq, double thet, double phi)  // ifreq is the index in NUBAR
{
/* static double sigfE[4] = {24.4, 14.9, 23.4, 18.9};  these are for the seavey data values
static double sigfH[4] = {26.3, 17.1, 14.5, 15.8};	
static double sigfreq[4] = {300.e6,600.e6,900.e6,1200.e6}; */
//   sigma = FWHM/(2sqrt(2ln2))

// spline interpolated values for 250-1150 MHz in 100 MHz steps  OLD values
//static double sigE[NFREQ] = {26.93 , 21.90 , 17.52 , 15.04 , 15.63 , 18.77 , 22.25 , 23.84 , 22.90 , 20.41};
//static double sigH[NFREQ] = { 28.08 , 24.52 , 21.14 , 18.25 , 16.19 , 15.01 , 14.54 , 14.56 , 14.93 , 15.49};
// Guesstimates of new values based on Seavey data:
static double sigE[NFREQ] = { 33 ,  34 ,  21 , 15 ,  25 , 24 , 22 , 22 , 20 , 20};
static double sigH[NFREQ] = { 26 , 25,  25,  25,  16,  14,  18, 20, 15,  10};
	double degx,degy,powerfac;
	
	degx = 180./PI*thet;
	degy = 180./PI*phi;
	 // set to have a -18dB floor  3/13/2005 --PG
	
	powerfac = 0.01 + exp(-(degx*degx/(2.*sigH[ifreq]*sigH[ifreq]) + degy*degy/(2.*sigE[ifreq]*sigE[ifreq])) );
			
	return(powerfac);
}



double get_rician(double Sre)
{
/* returns Rician with signal level Sre measured in rms sigmas */
 	double re,im,A, gasdev();
    re = gasdev();
    im = gasdev();
    
    A = sqrt((Sre+re)*(Sre+re) + im*im);
    return(A);
} 
    
    

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

int getL3( int L2pat[] )
{
	int i, L3, iant, jant;

	for(iant=0;iant<NANT;iant++){
		for(jant=0;jant<NANT;jant++){
		if(        L2pat[iant]==1 
			&& L2pat[jant]==1 
			&& PHISECTOR[iant]==PHISECTOR[jant]){
		return(L3=1);
		}
	     }
	}
	return(L3=0);
}
