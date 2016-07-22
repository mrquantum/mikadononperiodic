#include <iostream>
#include <math.h>
#include "nrutil.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>


#define ITMAX 1e6
#define EPS 1.0e-10
#define TOL 2.0e-8
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#include "structs.h"
#include "cgmethod.h"

using namespace std;
using namespace Eigen;

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double,networkinfo),networkinfo info);


double *dvector( int size )
{
    double *p = new double[size];
    if (!p) {
        nrerror("Allocation of dvector failed!");
    }else{
        return p;
    }
}



double dbrent(double ax, double bx, double cx, double (*f)(double,networkinfo),
              double (*df)(double,networkinfo),double tol, double *xmin,networkinfo info)
{
    int iter, ok1, ok2;
    double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
    double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
    
    a=(ax<cx?ax:cx);
    b=(ax>cx?ax:cx);
    
    x=w=v=bx;
    fw=fv=fx=(*f)(x,info);
    dw=dv=dx=(*df)(x,info);
    
    for(iter=1;iter<=ITMAX;iter++){
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+ZEPS;
        tol2=2.0*tol1;
        if(fabs(x-xm)<=tol2-0.5*(b-a)){
            *xmin=x;
            return fx;
        }
        if(fabs(e)>tol1){
            d1=2.0*(b-a);
            d2=d1;
            if(dw!=dx) d1=(w-x)*dx/(dx-dw); //secant method
            if(dv!=dx) d2=(v-x)*dx/(dx-dv);
            
            u1=x+d1;
            u2=x+d2;
            
            ok1=(a-u1)*(u1-b)>0.0 && dx*d1<=0.0;
            ok2=(a-u2)*(u2-b)>0.0 && dx*d2<=0.0;
            olde=e;
            e=d;
            if(ok1||ok2){
                if(ok1&&ok2){
                    d=(fabs(d1)<fabs(d2)?d1:d2);
                }else if(ok1){
                    d=d1;
                }else{
                    d=d2;
                }
                if(fabs(d)<=fabs(0.5*olde)){
                    u=x+d;
                    if(u-a<tol2||b-u<tol2) d=SIGN(tol1,xm-x);
                } else{
                    d=0.5*(e=(dx>=0.0?a-x:b-x));
                }
            }else{
                d=0.5*(e=(dx>=0.0?a-x:b-x));
            }
        }else{
            d=0.5*(e=(dx>0.0?a-x:b-x));
        }
        
        if(fabs(d)>=tol1){
            u=x+d;
            fu=(*f)(u,info);
        }else{
            u=x+SIGN(tol1,d);
            fu=(*f)(u,info);
            if(fu>fx){
                *xmin=x;
                return fx;
            }
        }
        du=(*df)(u,info);
        if(fu<=fx){
            if(u>=x) a=x; else b=x;
            MOV3(v,fv,dv,w,fw,dw)
            MOV3(w,fw,dw,x,fx,dx)
            MOV3(x,fx,dx,u,fu,du)
        }else{
            if(u<x) a=u; else b=u;
            if(fu<=fw || w==x){
                MOV3(v,fv,dv,w,fw,dw)
                MOV3(w,fw,dw,u,fu,du)
            }else if(fu<fv||v==x||v==w){
                MOV3(v,fv,dv,u,fu,du)
            }
        }
    }
    nrerror("Too many iterations in dbrent");
    return 0.0;
}

double brent(double ax,double bx,double cx,double (*f)(double,networkinfo info),double tol,
            double *xmin,networkinfo info)
{
    int iter;
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0;
    
    a=(ax<cx    ?       ax      :cx);
    b=(ax>cx    ?       ax      :cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x,info);
    
    for(iter=0;iter<ITMAX;iter++){
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if(fabs(x-xm)<=(tol2-0.5*(b-a))){
           *xmin=x;
           return fx;
        }
        if(fabs(e)>tol1){
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if(q>0.0) p=-p;
            q=fabs(q);
            etemp=e;
            e=d;
            if(fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x)){
                d=CGOLD*(e=(x>=xm?a-x:b-x));
            }else{
                d=p/q;
                u=x+d;
                if(u-a<tol2||b-u<tol2)
                    d=SIGN(tol1,xm-x);
            }
            
        } else {
            d=CGOLD*(e=(x>=xm?a-x:b-x));
        }
        u=(fabs(d)>=tol1?x+d:x+SIGN(tol1,d));
        fu=(*f)(u,info);
        if(fu<=fx){
            if(u>=x) a=x; else b=x;
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
        } else{
            if(u<x) a=u; else b=u;
            if(fu<=fw || w==x){
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if(fu<=fv || v==x||v==w){
                v=u;
                fv=fu;
            }
        }
        
    }
    nrerror("Too many iterations in Brent");
    *xmin=x;
    return fx;
}




int ncom;
double *pcom, *xicom, (*nrfunc)(double [],networkinfo info);

void linmin(double p[],double xi[], int n,double *fret,double (*func)(double [],networkinfo info),networkinfo info)
{
    
    double brent(double ax, double bx, double cx,
                double (*f)(double,networkinfo info), double tol, double *xmin,networkinfo info);
    double f1dim(double x,networkinfo info);
//     void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb, 
//                 double *fc,double (*func)(double));
//     
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    
    ncom=n;
    pcom=dvector(n);
    xicom=dvector(n);
    nrfunc=func;
    
    for(j=0;j<n;j++){
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0;
    xx=1.0;
    
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim,info);
    *fret=brent(ax,xx,bx,f1dim,TOL,&xmin,info);
    
    double *xcopy = dvector(n);
    
    
    
    for(j=0;j<n;j++){ //x=x+alpha*d
        xcopy[j] = xi[j];
        xi[j]*=xmin;
        p[j]+=xi[j];
    }
    /*
    for(j=0; j < n; ++j ){
        xi[j] = xcopy[j];
    }
    */
    delete pcom;
    delete xicom;
    delete xcopy;
    
}

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double [],networkinfo info);

double f1dim(double alpha,networkinfo info)
{
    int j;
    double f, *xt;
    
    xt=dvector(ncom);
    for(j=0;j<ncom;j++) xt[j]=pcom[j]+alpha*xicom[j];
    f=(*nrfunc)(xt,info);
    delete xt;
    return f;
}

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double,networkinfo info),networkinfo info)
//Given a function funct and given intitial points ax and bx this routine
//searches in the downhill direction and returns new poits ax bx cx that bracket 
//a minimum of this function. The values at those points are also returned
{
    double ulim,u,r,q,fu,dum;
    
    *fa=(*func)(*ax,info);
    *fb=(*func)(*bx,info);
    
    if(*fb>*fa){
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx,info);
    while(*fb>*fc){
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
            (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        
        if((*bx-u)*(u-*cx)>0.0){
            fu=(*func)(u,info);
            if(fu<*fc){
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            } else if(fu>*fb){
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u,info);
        }else if((*cx-u)*(u-ulim)>0.0){
            fu=(*func)(u,info);
            if(fu<*fc){
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                SHFT(*fb,*fc,fu,(*func)(u,info))
            }
        }else if((u-ulim)*(ulim-*cx)>0.0){
            u=ulim;
            fu=(*func)(u,info);
        }else{
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u,info);
        }
        SHFT(*ax,*bx,*cx,u);
        SHFT(*fa,*fb,*fc,fu)
    }
}

void frprmn(double p[],int n, double ftol, int *iter, double *fret,
            double (*func)(double [],networkinfo), void (*dfunc)(double [],double [],networkinfo),networkinfo &info,ofstream &EN)
{
    
    
    
    //void dlinmin(double p[],double xi[],int n, double *fret,
    //            double (*func)(double [],networkinfo),void (*dfunc)(double[], double [],networkinfo info));
    int j,its;
    double gg,gam,fp,dgg;
    double *g,*h,*xi;
    //initialisation
    g=dvector(n);
    h=dvector(n);
    xi=dvector(n);

    fp=(*func)(p,info);
    (*dfunc)(p,xi,info); //calculates the gradient vector xi at position p

    for(j=0;j<n;j++){
        g[j]=-xi[j];
        xi[j]=h[j]=g[j];
    }
    
    //the loop itself
    for(its=0;its<=ITMAX;its++){
        *iter=its; //how many cg steps
        dlinmin(p,xi,n,fret,func,dfunc,info);
        if(2.0*fabs(*fret-fp) <=ftol*(fabs(*fret)+fabs(fp)+EPS) && its>100){
            delete g;
            delete h;
            delete xi;
            return;
        }
        fp=*fret;
        double gsq=0.0;
        
        for(int ii=0;ii<n;ii++){
            gsq=gsq+xi[ii]*xi[ii];
        }
//         EN<<info.sheardeformation<<"\t"<<fp<<"\t"<<gsq<<endl;
        (*dfunc)(p,xi,info);
        dgg=gg=0.0;
        for(j=0;j<n;j++){
            gg+=g[j]*g[j];
            //dgg+=xi[j]*xi[j];
            dgg+=(xi[j]+g[j])*xi[j];
        }
        if(gg==0.0){
            delete g;
            delete h;
            delete xi;
            return;
        }
        gam=dgg/gg;
        for(j=0;j<n;j++){
            g[j]=-xi[j];
            xi[j]=h[j]=g[j]+gam*h[j];
        }
    }

    nrerror("Too many iterations in frprmn");
}

void (*nrdfun)(double [], double  [],networkinfo info);

void dlinmin(double p[], double xi[],int n, double *fret, double (*func)(double [],networkinfo),
             void (*dfunc)(double [], double [],networkinfo),networkinfo info)
{
    double dbrent(double ax, double bx, double cx,
                  double (*f)(double,networkinfo), double (*df)(double,networkinfo), double tol, double *xmin,networkinfo info);
    double f1dim(double x,networkinfo info);
    double df1dim(double x,networkinfo info);
    void mnbrak(double *ax, double *bx, double *cx, double *fa, double *a,
                double *fc, double (*func)(double,networkinfo),networkinfo info);
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    
    
    ncom=n;
    pcom=dvector(n);
    xicom=dvector(n);
    nrfunc=func;
    nrdfun=dfunc;
    
    for(j=0;j<n;j++){
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0;
    xx=1.0;
    
    mnbrak(&ax,&xx,&bx, &fa,&fx,&fb,f1dim,info);
    *fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin,info);
    for(j=0;j<n;j++){
        xi[j]*=xmin; 
        p[j]+=xi[j];
    }
    delete xicom;
    delete pcom;
    
}

extern void (*nrdfun)(double [], double [],networkinfo);
double df1dim(double x,networkinfo info)
{
    int j;
    double df1=0.0;
    double *xt, *df;
    
    xt=dvector(ncom);
    df=dvector(ncom);
    for(j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    (*nrdfun)(xt,df,info);
    for(j=0;j<ncom;j++) df1+=df[j]*xicom[j];
    delete df;
    delete xt;
    return df1;
}

