#include "su.h"
#include "cwp.h"
#include "segy.h"
#include "stdio.h"

char *sdoc[] = {
"  Forward modling                             \n",
NULL};

float *broad_riker( float f1,float f2,float dt,int ntt );
float *broad_wavelet(char *path,float dt,float maxf,int ntt);
void convmtx(float *w,int lw,int nt,float **a0);
float maxfsinc (float x,float maxf);
void intt(float *sqa, float ct, float dt, float *va, int nt1);
void halfdir(float *data, float dt, int nt);

segy tr, tr1;

int main(int argc, char **argv)
{   
	FILE *fp, *fp1, *fp2;           
    int i, j, ix, it, ik, ntt;
	char path[200],filename[200]; 
	
	float *riker, **data;
    float f=30;
    int npx, ns;
    float v, dt, half_triker;
	
	float  tmid, tmin, tmax, t, va;
	int tmin1, tmax1;
	
	float sx, gx;
	float xx, tt, dsx, dgx, ts, tg, qtmp, qtmp1;
	float *s,*g;

	float **cmt, **cd;
	float off, sum;

	initargs (argc, argv);
   
   /**** model(signal dot) information ************************************/
	npx=182;  /*input sx gx number*/
	
	xx=500.0; /*single point position*/  
	tt=300.0;
	
	ns=500;   /* sampling */  
    v=2000;  /* velocity */
	off=100; /* offset   */
 
   
    /******************* riker *******************************/  
	dt=0.002; 
    riker=alloc1float(101);
    memset(riker, (int)'\0', sizeof(float)*101);
    for(i=0;i<101;i++)
    {
		t=(i-50)*dt; 
		riker[i]=(1-2*PI*PI*t*t*f*f)*exp(-PI*PI*t*t*f*f);
    }
	
    sprintf(filename,"../result/riker.bin");
    fp=fopen(filename,"wb");
	
    for(i=0;i<101;i++)
	{
		fwrite(&riker[i],4,1,fp);
	}
    fclose(fp);
	
	/************** mod_riker ******************************/
    half_triker=100*dt/2;
	halfdir(riker,dt,101); 
	
	sprintf(filename,"../result/mod_riker.bin");
	fp2=fopen(filename,"wb");
	
	for(j=0;j<101;j++)
	{ 
		//warn("riker[%d]=%f",j,riker[j]);
		fwrite(&riker[j],4,1,fp2);
	}
	fclose(fp2);

	/***** get geosystem **********************************/
    fp=fopen("sx.bin","rb");
    fp1=fopen("gx.bin","rb");
	warn("Get System ok!");
	
	s=alloc1float(npx);
	memset(s,(int)'\0',sizeof(float)*npx);
	g=alloc1float(npx);
	memset(g,(int)'\0',sizeof(float)*npx);
	
	data=alloc2float(ns,npx);
	memset(data[0], (int)'\0', sizeof(float)*ns*npx);
	
	for(i=0;i<npx;i++)	    
	{          
		//get every CDP's parameters  
		fread(&sx,4,1,fp);
		fread(&gx,4,1,fp1);
		s[i]=sx;
		g[i]=gx;
		
		dsx=fabs(sx-xx);
		dgx=fabs(gx-xx);
		ts=sqrt(dsx*dsx+tt*tt)/v;
		tg=sqrt(dgx*dgx+tt*tt)/v;
		tmid=(ts+tg);
		
		qtmp =tt/v/v/sqrt(2*PI*ts)/ts;
		qtmp1=tt/v/v/sqrt(2*PI*tg)/tg;
		//warn("qtmp=%f,qtmp1=%f\t",qtmp,qtmp1);
		
		tmin=(tmid-half_triker)/dt;
		tmax=(tmid+half_triker)/dt;
		warn("tmin=%f, tmax=%5f, tmid=%f,half_triker=%f\t",tmin,tmax,tmid,half_triker);
		tmin1=(int)(tmin);
       	tmax1=(int)(tmax);

       	if(tmin1<0) 
		{
			tmin1=0;     //timn=1
		}
       	if(tmax1>ns-100)
		{
			tmax1=ns-1;
		}			                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
			
		for(it=tmin1;it<=tmax1;it++)
		{
			t=it*dt-tmid; 
			t+=half_triker;        
	
			intt(riker, t, dt, &va, 101);   //interpolation
			
			data[i][it]=va*qtmp*qtmp1;
       	}
		
	}
	
	fclose(fp);
	fclose(fp1);
   
	sprintf(filename,"../result/data.su");
	fp1=fopen(filename,"wb");

	for(i=0;i<npx;i++)
	{
		tr1.cdp=i+1;
		tr1.ns=ns;
		tr1.dt=dt*1.0E6;
		tr1.sx=s[i];
		tr1.gx=g[i];
		tr1.offset=off;
		
		memcpy((void *)tr1.data, (const void *)data[i], sizeof(float)*ns); 
		fputtr(fp1, &tr1);	
	}  
	fclose(fp1);
 
 
 /************************* test ***********************************************************/
	
	/*
	//sinc_data
	sprintf(path,"../data/");
	ntt=101;
	warn("ntt=%d,max_f=%f",ntt,2*f);
	//riker=broad_wavelet(path,dt,2*f,ntt);
	riker=broad_riker(0,100,dt,ntt );
	//conv
	cmt=alloc2float(ns,ns);
	memset(cmt[0], (int)'\0', sizeof(float)*ns*ns);

	cd=alloc2float(ns,npx);
	memset(cd[0], (int)'\0', sizeof(float)*ns*npx);
	
	convmtx(riker,ntt,ns,cmt);

	//data wcm*r
	for(ix=0;ix<npx;ix++)
	{
		for(it=0;it<ns;it++)
		{
			sum=0.0;
			for(ik=0;ik<ns;ik++)
			{
				sum+=data[ix][ik]*cmt[ik][it];
			}
			cd[ix][it]=sum;
		}
		
	}

	sprintf(filename,"../result/data_sinc.su");
	fp1=fopen(filename,"wb");

	for(i=0;i<npx;i++)
	{
		tr1.cdp=i+1;
		tr1.ns=ns;
		tr1.dt=dt*1.0E6;
		tr1.sx=s[i];
		tr1.gx=g[i];
		tr1.offset=off;
		memcpy((void *)tr1.data, (const void *)cd[i], sizeof(float)*ns); 
		fputtr(fp1, &tr1);	
	}  
	fclose(fp1);

	free1float(riker);
	free2float(data);
	
	return 0;
	*/
}



//wavelet convolution matrix
void convmtx(float *w,int lw,int nt,float **a0) //a~(ld+lw-1)*ld
{
	int ix,iy;
	float **a;
	int lx,ly,lw2;

	lx=nt;
	lw2=lw/2;
	ly=nt+lw-1;
	
	a=alloc2float(ly,lx);
	
	memset(a[0], (int)'\0', sizeof(float)*ly*lx);	
	
	for(iy=0;iy<lw;iy++)
	{
		a[0][iy]=w[iy];
	}
	
	for(ix=1;ix<lx;ix++)
	{
		for(iy=ix;iy<ix+lw;iy++)
		{
			a[ix][iy]=w[iy-ix];	
		}
	}
	
	for(ix=0;ix<lx;ix++)
	{
		for(iy=lw2;iy<ly-lw2;iy++)
		{
			a0[ix][iy-lw2]=a[ix][iy];
		}
	}
	free2float(a);
}


float *broad_wavelet(char *path,float dt,float maxf,int ntt)
{
	int it;
	char filename[100];
	FILE *fp;
	float t;
	float *wtmp;
	
	wtmp=alloc1float(ntt);
	memset(wtmp, (int)'\0', sizeof(float)*ntt);
	
	for(it=0;it<ntt;it++)
	{
		t=(it-ntt/2)*dt; 
		wtmp[it]=maxfsinc(t,maxf);
	}
	
	sprintf(filename, "%s/sinc.bin", path);
  	fp=fopen(filename, "wb");
	
	for(it=0;it<ntt;it++)
	{
		fwrite(&wtmp[it],sizeof(float),1,fp);
	}		
	fclose(fp);
	
	return wtmp;
	free1float(wtmp);
	
}

float maxfsinc (float x,float maxf) //max_f
{
	float pix;

	if (x==0.0) 
	{
		return 1.0;
	}
	else
	{
		pix=2*PI*maxf*x;
		return sin(pix)/pix;
	}
}

void halfdir(float *data, float dt, int nt)
{
  	int nw, ntfft, iw, it;
  	float dw, tmpw, *datafft, w, t;
  	complex *dataw, ctmp;

  	ntfft=npfar(nt);
  	nw=ntfft/2+1;
  	dw=2.0*PI/(ntfft*dt); 

  	dataw=alloc1complex(nw);
	datafft=alloc1float(ntfft);
	
  	for(it=0; it<ntfft; it++)
	{
		datafft[it]=0.0;
	}
  	for(it=0, t=0; it<nt; it++, t+=dt)
	{
		datafft[it]=data[it];
    }
	
  	pfarc (-1, ntfft, datafft, dataw);
	
  	for(iw=0, w=0; iw<nw; iw++, w+=dw)
	{ 
      	tmpw=w;//根号
		//ctmp=cmplx(cos(PI*1/4.0)*tmpw,-sin(PI*1/4.0)*tmpw);
		ctmp=cmplx(cos(PI*1/2.0)*tmpw,sin(PI*1/2.0)*tmpw);  
		dataw[iw]=cmul(dataw[iw],ctmp);
		//dataw[iw]=crmul(dataw[iw],tmpw);
    }
	
	for(it=0; it<ntfft; it++)
	{
		datafft[it]=0.0;
	}

  	pfacr(1, ntfft, dataw, datafft);

  	for (it=0, t=0; it<ntfft; it++, t+=dt)
	{
		datafft[it]/=ntfft;
    }
  	for (it=0; it<nt; it++)
	{
		data[it]=datafft[it];
	}

  	free1float(datafft);
  	free1complex(dataw);
}

void intt(float *sqa, float ct, float dt, float *va, int nt1)
{
	float at;
	int i2;
	float  es, as, bs;
	float f1, f2, f3, f4;

	at=ct/dt;
	i2=(int)at;
	f1 = sqa[i2-1] ;
	f2 = sqa[i2] ;
	f3 = sqa[i2+1] ;
	f4 = sqa[i2+2];
	es = at - i2;
	as = es*(1 - es);
	bs = 3*as+6;
	*va=0.1666666666666666667*(bs*(f2+es*(f3-f2))-as*(f1+f1+es*(f4-f1)+f4));
}
    
/*** broad_riker ***/
float *broad_riker( float f1,float f2,float dt,int ntt )
{
	int ntt2,ix;
	float *w;
	float t;

	w=alloc1float(ntt);
	memset(w, 0, sizeof(float)*ntt);

	ntt2=ntt/2;
	
	for(ix=0;ix<ntt;ix++)
    {
		t=(ix-ntt2)*dt; 

		w[ix]=f2*exp(-PI*PI*t*t*f2*f2)-f1*exp(-PI*PI*t*t*f1*f1);
		w[ix]=w[ix]/(f2-f1);
    }

	return w;
	free1float(w);	
}
