/*
 This is modification of SIGPROC barycentric correction code. The original code does not do
 frequency correction which is implimented in this code for SETI application

 Vishal Gajjar 
 Nov 4th, 2020
*/   

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/* 
   barycentre.c - refer a filterbank/timeseries file to the rest
   frame of the solar system barycentre by computing an appropriate
   set of polynomial coefficients for the data sampling interval,
   then keeping track of the difference between elapsed time (topo)
   versus the time computed from the coefficients (barycentric) and
   adding or removing time samples so that the two timescales stay
   fixed. Added samples are set to zero.

   Modification history:
   March 18, 2007, drl added -mypar option to override position in header
                   this is useful for observations where there may be a small
                   position offset between the telescope position and the true
                   pulsar position and the resulting difference can cause
                   the barycentric time series to drift. Also included refdm
                   in the calculation.
*/
#include "dedisperse.h"
#include <string.h>
char polyco_filename[80];
double fcent;

struct restruct {
	double mjdbary;
	double velrel;
};

/* subroutine to call TEMPO to calculate a polyco.bar file for barycentering */
char *make_polycofile(char ra[],char dec[],char topo[], char site,
			    double mjdtopo, double tsamp)
{
  FILE *resid2,*parfile, *tzfile;
  float junk;
  double mjdbary;
  char *polycofilename;
  polycofilename=(char *) malloc(80);

  parfile=fopen("tssb.par","w");
  fprintf(parfile,"PSR 0000+00\n");
  fprintf(parfile,"RAJ %s\n",ra);
  fprintf(parfile,"DECJ %s\n",dec);
  fprintf(parfile,"F0 1.0\n");
  fprintf(parfile,"DM %f\n",refdm);
  fprintf(parfile,"PEPOCH %s\n",topo);
  fclose(parfile);
  tzfile=fopen("tz.in","w");
  fprintf(tzfile,"%c    2  30  9 %lf\n",site,fcent);
  fprintf(tzfile,"\n \n");
  fprintf(tzfile,"0000+00 60 9 12 %lf\n",fcent);
  fclose(tzfile);
  tzfile=fopen("runtempo.csh","w");
  fprintf(tzfile,"#!/bin/csh\n",site);
  fprintf(tzfile,"tempo -z -f tssb.par << EOD\n");
  fprintf(tzfile,"%f %f\n",mjdtopo-1.0,mjdtopo+1.0);
  fprintf(tzfile,"EOD");
  fclose(tzfile);
  system("csh runtempo.csh > /dev/null");
  system("mv polyco.dat polyco.bar");
  system("rm -f tssb.par tz.in tz.tmp");
  system("rm -f fort.22 tempo.lis runtempo.csh ");
  strcpy(polycofilename,"polyco.bar");
  return(polycofilename);
}

/* subroutine to call TEMPO to calculate the barycentric MJD */
struct restruct barycentric_time(char ra[],char dec[],char topo[], char site,
			    double mjdtopo)
{
  FILE *resid2,*parfile, *timfile;
  float junk;
  struct restruct baryval;
  double velrel;
  double mjdbary;
  parfile=fopen("tssb.par","w");
  fprintf(parfile,"PSR 0000+00\n");
  fprintf(parfile,"RAJ %s\n",ra);
  fprintf(parfile,"DECJ %s\n",dec);
  fprintf(parfile,"F0 1.0\n");
  fprintf(parfile,"DM 0.0\n");
  fprintf(parfile,"PEPOCH %s\n",topo);
  fclose(parfile);
  timfile=fopen("tssb.tim","w");
  fprintf(timfile,"%c    0  0000+00 %.3lf %.13f     1.00\n",site,fcent,mjdtopo);
  fprintf(timfile,"%c    0  0000+00 %.3lf %.13f     1.00\n",site,fcent,mjdtopo);
  fprintf(timfile,"%c    0  0000+00 %.3lf %.13f     1.00\n",site,fcent,mjdtopo);
  fclose(timfile);
  system("tempo tssb.tim > /dev/null");
  resid2=fopen("resid2.tmp","r");
  // Orig
  //fread(&junk,4,1,resid2);
  //fread(&mjdbary,8,1,resid2);

  //Adding read of observed frequency to get the relative velocity 	
  int ii;
  double dd,femit;
  long long ll;
  static double d[9];
  fread(&ii, sizeof(int), 1,resid2);
  fread(&d, sizeof(double), 9, resid2);
  //fprintf(stderr,"%.12f\t%.12f\n",d[0],d[4]);
  femit = d[4];
  velrel = femit / fcent - 1.0;
 
  //fprintf(stderr,"%.12f\t%.12f\n",fobs,velrel);
  rewind(resid2);
  fread(&junk,4,1,resid2);
  fread(&mjdbary,8,1,resid2);
  baryval.mjdbary = mjdbary;
  baryval.velrel  = velrel;
  //fprintf(stderr,"%.12f %.12f\n",mjdbary,velrel);
  fclose(resid2);
  system("rm -f tssb.tim tssb.par 0000+00.par");
  system("rm -f matrix.tmp tempo.lis resid2.tmp ");
  return  baryval;
  //return(mjdbary);
}

void leftRotatebyOne(float* arr, int sizearr);
void rightRotatebyOne(float* arr, int sizearr);

/*Function to left rotate arr[] of size n by d*/
void leftRotate(float* arr, int shiftl, int sizearr)
{
    int shiftj=0;
    for (shiftj = 0; shiftj < shiftl; shiftj++){
        leftRotatebyOne(arr, sizearr);
    }	
}

void leftRotatebyOne(float* arr, int sizearr)
{
    float temp = arr[0];
    int shifti=0;
    
    //fprintf(stderr,"arr %d %d\n",shifti,sizearr);
    for (shifti = 0; shifti < sizearr - 1; shifti++){
        arr[shifti] = arr[shifti + 1];
	//fprintf(stderr,"arr %d %f\n",shifti,arr[shifti]);
        //printf("arr %c\n",arr[i]);
    }

    arr[shifti] = temp;
}

float* squeeze(float *arr, int nchans, double fch1, double foff, double velrel,int sqchan)
{	
   int seqi,seqj,seqk;
   double seqfreq1,seqfreq2;
   float *seqchanblk;
   seqchanblk = (float *)malloc(nchans*sizeof(float));
   memset(seqchanblk,0,nchans);
   seqi=seqj=seqk=0;
   seqfreq1=(fch1*(1 + velrel))*1000000;
   for(seqj=0;seqj<nchans;seqj++){
        seqfreq2 = ((fch1+(seqj+1)*foff)*(1 + velrel)*1000000); // subsequent channel frequency in barycentric frame
        if(seqfreq2-seqfreq1<(foff*1000000/2.0 && seqj+1<nchans)) {
        	fch1 = fch1+foff;
		//If two channels come closer than half of chan, add them
		seqchanblk[seqk] = (arr[seqj]+arr[seqj+1])/2.0;
		seqi+=1;
        }
	else {
		seqchanblk[seqk] = arr[seqj];	
		seqk+=1;
	}		
   	seqfreq1 = seqfreq1+foff*1000000; // Real frquency how they will get organize in the file 
   }
   //These number should match 
   //fprintf(stderr,"%d %d\n",seqi,seqj-seqk);
   //return seqi;
   sqchan=seqi;
   return seqchanblk;
}

/*Function to right rotate arr[] of size n by d*/
void rightRotate(float* arr, int shiftr, int sizearr)
{
    int shiftj=0;
    for (shiftj = 0; shiftj < shiftr; shiftj++)
        rightRotatebyOne(arr, sizearr);
}

void rightRotatebyOne(float* arr, int sizearr)
{
    int shifti=0;
    float temp = arr[sizearr-1];
    //fprintf(stderr,"arr %d %d\n",shifti,sizearr);
    for (shifti = sizearr-1; shifti > 0; shifti--){
        arr[shifti] = arr[shifti - 1];
    }
    arr[shifti] = temp;
}


char inpfile[80], outfile[80];
main (int argc, char *argv[]) 
{
  int drop,add,chandrop,chanadd,i,j,n,ntim,headersize,rah,ram,ded,dem;
  int ndropped=0,nadded=0,nbytes_per_sample,verbose=0,mypolyco=0;
  unsigned char *rawdata, *dummy;
  float *chanblk;
  double mjd, elapsed_time, barycentre_time,mjdtopostart,newtopo;
  double ras,des;
  char ra[80], dec[80], topo[80], sra[6], sde[6], site;
  char message[80], myparfile[80], line[80], key[80];
  struct POLYCO polyco;
  FILE *polycofile,*parfile;
  double velrel,mjdbary,nfreq1,nfreq2;
  struct restruct baryval;

  FILE *testfile;
  //testfile=fopen("FreqShift","w");

  if (argc<2 || help_required(argv[1])) {
    barycentre_help();
    exit(0);
  }
  print_version(argv[0],argv[1]);
  if (!file_exists(argv[1]))
    error_message("input file does not exist!");

  strcpy(inpfile,argv[1]);
  input=open_file(inpfile,"r");
  strcpy(outfile,"stdout");
  output=stdout;
  strcpy(myparfile,"");

  i=2;
  while (i<argc) {
    if (strings_equal(argv[i],"-verbose")) 
      verbose=1;
    if (strings_equal(argv[i],"-mypolyco")) 
      mypolyco=1;
    if (strings_equal(argv[i],"-myparfile")) 
      strcpy(myparfile,argv[++i]);
    i++;
  }

  if ((headersize=read_header(input))) {

    /* calculate centre frequency for use in TEMPO files */
    fcent=fch1+(double)(nchans/2)*foff;
    float origfch1=fch1;
    /* parse the header parameters for RA */
    angle_split(src_raj,&rah,&ram,&ras);
    if (ras<10.0) {
      sprintf(sra,"0%.1f",ras);
    } else {
      sprintf(sra,"%.1f",ras);
    }
    sprintf(ra,"%02d:%02d:%s",rah,ram,sra);
    /* parse the header parameters for DEC */
    angle_split(src_dej,&ded,&dem,&des);
    if (des<10.0) {
      sprintf(sde,"0%.1f",des);
    } else {
      sprintf(sde,"%.1f",des);
    }
    sprintf(dec,"%02d:%02d:%s",ded,dem,sde);
    /* now call TEMPO to calculate the barycentric MJD */
    sprintf(topo,"%.12f",tstart);
    site=tempo_site(telescope_id);
    if (verbose) 
      fprintf(stderr,"Telescope: %s TEMPO site code: %c\n",
	      telescope_name(telescope_id),site);
    if (!strings_equal(myparfile,"")) {
      parfile=open_file(myparfile,"r");
      while (fgets(line,80,parfile) != NULL) {
	strcpy(key,strtok(line," "));
	if (strings_equal(key,"RAJ"))
	  strcpy(ra,strtok(NULL," "));
	if (strings_equal(key,"DECJ"))
	  strcpy(dec,strtok(NULL," "));
      }
    }

    if (mypolyco) {
      strcpy(polyco_filename,"polyco.bar");
      if (verbose)
      fprintf(stderr,"Using barycentric polyco file: %s\n",polyco_filename);
    } else {
      strcpy(polyco_filename,make_polycofile(ra,dec,topo,site,tstart,tsamp));
      if (verbose)
      fprintf(stderr,"Created barycentric polyco file: %s\n",polyco_filename);
    }
      
    barycentric=1;
    
    double origtstart=tstart;
    baryval = barycentric_time(ra,dec,topo,site,tstart);
    mjdbary=baryval.mjdbary;
    float origmjdbary=mjdbary;
    velrel=baryval.velrel;
    mjdtopostart=tstart;
    double origbaryfch1 = origfch1*(1 + velrel);

    if (verbose) {
      fprintf(stderr,"Topocentric MJD %.12f\n",tstart);
      fprintf(stderr,"Barycentric MJD %.12f\n",mjdbary);
      fprintf(stderr,"Relative velocity %.12f\n",velrel);
    }
    /* write out header with barycentric MJD if required */
    send_string("HEADER_START");
    send_int("telescope_id",telescope_id); 
    send_int("machine_id",machine_id);
    send_coords(src_raj,src_dej,az_start,za_start);
    send_int("data_type",data_type);
    send_int("barycentric",1);
    send_int("pulsarcentric",0);
    if (nchans==1) send_double("refdm",refdm);
    if (fch1 == 0.0) 
      send_double("fch1",frequency_table[0]);
    else
      // Send barycentric frequency to header	    
      send_double("fch1",fch1*(1 + baryval.velrel));	      
      //send_double("fch1",fch1);
    send_int("nchans",nchans);
    if (nchans>1) send_double("foff",foff);
    send_int("nbits",nbits);  
    send_int("nifs",nifs);  
    send_double ("tstart",mjdbary); 
    send_double("tsamp",tsamp);
    send_int("nbeams",nbeams);
    send_int("ibeam",ibeam);
    send_string("HEADER_END");
    open_log("barycentre.monitor");

    ntim=nsamples(inpfile,headersize,nbits,nifs,nchans);
    //fprintf(stderr,"%d\n",ntim);
    mjd=tstart;
    polycofile=open_file(polyco_filename,"r");
    if (!read_polycoset(polycofile,&polyco)) {
      error_message("depolyco: error reading polyco file");
    } else {
      get_nearest_polyco(polyco_filename,mjd,&polyco);
    }

    i=n=drop=add=0;
    nbytes_per_sample=nchans*nbits*nifs/8;
    dummy=(char *) malloc(nbytes_per_sample);
    for (j=0;j<nbytes_per_sample;j++) dummy[j]=0;
    elapsed_time=barycentre_time=0.0;
    nfreq1=nfreq2=0;

    while (i<ntim) {
      n++;
      elapsed_time+=tsamp;
      //newtopo = mjdtopostart+elapsed_time/(86400); 
      baryval = barycentric_time(ra,dec,topo,site,mjd);
      barycentre_time+=tsamp*polyco_period(mjd,polyco);
      nfreq1 = (fch1*(1 + baryval.velrel))*1000000; //Emitted first channel frequency 
      
      if(i==1){
	for(j=0;j<nchans;j++){      
	nfreq2 = ((fch1+(j+1)*foff)*(1 + baryval.velrel)*1000000); // subsequent channel frequency in barycentric frame
	//fprintf(stderr,"%.12f\n",((fch1+j*foff) - (fch1+j*foff)/(1 + baryval.velrel))*1000000);	
	//nfreq1 = ((fch1/(1 + baryval.velrel)+j*foff)*1000000); // Real frquency how they will get organize in the file 
	//nfreq1 = ((fch1/(1 + baryval.velrel)+j*foff)*1000000); // Real frquency how they will get organize in the file 
	if(foff>0.0)
		if(baryval.velrel<0.0){
			chandrop=1; //Channels need to be squeezed 
			if(nfreq2-nfreq1<(foff*1000000/2.0)) {
				//nfreq1=nfreq1+foff*1000000;
				fch1 = fch1+foff;
				//fprintf(testfile,"%.12f\t %.12f\n",nfreq1,nfreq2-nfreq1);
			}	
		}	
	//fprintf(testfile,"nfreq1 %.12f nfreq2 %.12f foff %.12f nf2-nf1 %.12f\n",nfreq1,nfreq2,foff*1000000,nfreq2-nfreq1);
	//if(foff<0.0)
	//if(nfreq1-nfreq2<foff*1000000) nfreq2=nfreq2+foff*1000000;
	//fprintf(testfile,"%.12f\t %.12f\n",nfreq1,nfreq2-nfreq1);
	nfreq1 = nfreq1+foff*1000000; // Real frquency how they will get organize in the file 
	}
      	// How first frequency and center frequency change as function of time
      	//fprintf(stderr,"%.12f %.12f %.12f %.12f\n",fch1,fch1*(1 + baryval.velrel),(fch1+nchans*foff),(fch1+nchans*foff)*(1 + baryval.velrel));
      }		

      //fprintf(stderr,"%.10f %.10f\n",polyco_period(mjd,polyco),polyco);
      //fprintf(stderr,"%d %.12f %.12f %.12f %.12f\n",i,mjd,baryval.mjdbary,barycentre_time,baryval.velrel);
      //How first frequency and center frequency change as function of time
      //fprintf(stderr,"%.12f %.12f %.12f %.12f\n",fch1,fch1*(1 + baryval.velrel),(fch1+nchans*foff),(fch1+nchans*foff)*(1 + baryval.velrel));
      if (elapsed_time-barycentre_time>tsamp) {
	add=1;
	elapsed_time-=tsamp;
      } else if (barycentre_time-elapsed_time>tsamp) {
	drop=1;
	elapsed_time+=tsamp;
      }
      if (drop || add) {
	rawdata=(char *) malloc(n*nbytes_per_sample);
	fread(rawdata,1,n*nbytes_per_sample,input);
	if (drop) {
	  ndropped++;
	  fwrite(rawdata,1,(n-1)*nbytes_per_sample,output);
	} else {
	  nadded++;
	  fwrite(rawdata,1,n*nbytes_per_sample,output);
	  fwrite(dummy,1,nbytes_per_sample,output);
	}
	n=0;
	free(rawdata);
	drop=add=0;
        sprintf(message,"time:%.1fs",elapsed_time);
        update_log(message);
      }
      mjd+=tsamp/86400.0;
      get_nearest_polyco(polyco_filename,mjd,&polyco);
      i++;
    }
    chanblk = (float *)malloc(nchans*sizeof(float));
    //n=2;	
    baryval = barycentric_time(ra,dec,topo,site,origtstart);
    nfreq1 = (origfch1*(1 + baryval.velrel));
    fprintf(stderr,"nfreq1 %.12f origbaryfch1 %.12f\n",nfreq1,origbaryfch1);

    if (n) {
      //rawdata=(char *) malloc(n*nbytes_per_sample);
      //rawdata=(unsigned char*) malloc(nchans*nbits);
      //chanblk = (float *)malloc(nchans*sizeof(float));
      for(j=0;j<nchans;j++) chanblk[j] = 0.0;
      //fprintf(stderr,"%d\n",sizeof(chanblk[10]));
      mjd=origtstart;
      int lshift=0;
      int sqchan=0;
      for(i=0;i<n;i++){
	        baryval = barycentric_time(ra,dec,topo,site,mjd);
	        //barycentre_time+=tsamp*polyco_period(mjd,polyco);
		nfreq1 = (origfch1*(1 + baryval.velrel));
	        for(j=0;j<nchans;j++) fread(&chanblk[j],1,sizeof(float),input);
		//rightRotate(chanblk,i,nchans);
		if(baryval.velrel<0 && foff>0){
				if(nfreq1-origbaryfch1>foff) {
					lshift = ceil((nfreq1-origbaryfch1)/foff);
					leftRotate(chanblk,lshift,nchans);
				}
				//sqchan=squeeze(chanblk, nchans, origfch1, foff, baryval.velrel);
				chanblk=squeeze(chanblk, nchans, origfch1, foff, baryval.velrel,sqchan);
				fprintf(stderr,"At %lf freq diff  %lf Hz shifting by %d channels and squeeze by %d channels\n",baryval.mjdbary,(nfreq1-origbaryfch1)*1000000,lshift,sqchan);
		}		
		//leftRotate(rawdata,3*i,sizeof(rawdata));
		for(j=0;j<nchans;j++) fwrite(&chanblk[j],1,sizeof(float),output);
		mjd+=tsamp/86400.0;
      }	      
      //for(j=0;j<nchans;j++) free(chanblk[i]);
    }
    free(chanblk);
    free(dummy);
    update_log("finished");
    close_log("barycentre.monitor");
    if (verbose && nadded)
      fprintf(stderr,"added %d samples\n",nadded);
    if (verbose && ndropped)
      fprintf(stderr,"dropped %d samples\n",ndropped);
  }
  //fclose(testfile);
}


