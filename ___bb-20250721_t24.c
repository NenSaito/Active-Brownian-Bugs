// ___bb-20250721_t24.c
// % gcc ___bb-20250721_t24.c -o out -lm -O3
// % time ./out seed n_Dr n_v0 e_Dth factor > _bb-20250721_t24-data-seed-n_Dr-n_v0-e_Dth-factor.txt &
/*
- factor=0,1,...,7
- th0=atan2(cos(ce1[i].th),sin(ce1[i].th))+M_PI*factor/4.0

# if n==tau0[9]-1
- calculate global-direction
- output
# if n==tau0[9]
- update tau0 and T0 and M0
- set ce0
- p_[i][k]=p(ce0[i].N)
- nq[k]=1
# if n<tau0[k]
- update ce0
- p_[i][k]+=p(ce0[i].N)
- nq[k]+=1
# if n==tau0[k]
- p_[i][k]+=0
- nq[k]+=0
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "MT.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

//define
#define dt (0.01)
#define L (1.0)
#define R (0.1)
#define Ns (50)
#define p0 (0.85)
#define q0 (0.15)

//structure
typedef struct{ 
    double x;
    double y;
    double th;
    int N;      
}Cell;              
                    
//function              
double p(int);      
double bc_r(double);
double bc_l(double);
int bc_b(int);
int bN(int,int,int);
double Nrand(int);
double ConvertToParameter(int);

int main(int argc,char* argv[]){
    time_t t_start=time(NULL);

    //command line argument
    int seed=strtol(argv[1],NULL,10);
    double Dr=2.0*exp((-6.0+strtol(argv[2],NULL,10)/4.0)*log(10)),v0=2.0*exp((-4.0+strtol(argv[3],NULL,10)/4.0)*log(10)),Dth=ConvertToParameter(strtol(argv[4],NULL,10));
    int factor=strtol(argv[5],NULL,10);

    //define
    int tau=(int)(3.0/Dth/dt),tau_th=(int)(1.0/Dth/dt),Mj=3000,nB=(int)(L/R);
    //declaration
    int tau0[10],nq[10];
    double p_[Mj][10],th0;
    Cell ce1[Mj],ce0[Mj];
    //initialisation
    double T0=1.0/Dth;
    int M=1500,M0=0,J0=0;
    for(int k=0;k<10;k++){
        tau0[k]=tau_th;
        nq[k]=1;
    }
    init_genrand(seed);

    //action
    for(int n=0;n<=tau;n++){
        //initialisation
        if(n==0){
            for(int i=0;i<M;i++){
                ce1[i].x=L*genrand_real2();
                ce1[i].y=L*genrand_real2();
                ce1[i].th=2.0*M_PI*genrand_real2();
                ce1[i].N=0;
            }
        }
        //update
        else if(n>0){
            //number of neighborhoods
            int cB[nB*nB],mB[nB*nB][4*Ns];
            for(int b=0;b<nB*nB;b++)
                cB[b]=0;
            for(int i=0;i<M;i++){
                ce1[i].N=0;
                int b=(int)(ce1[i].y/R)*nB+(int)(ce1[i].x/R);
                mB[b][cB[b]]=i;
                cB[b]+=1;
            }
            for(int i=0;i<M;i++){
                for(int k=0;k<9;k++){
                    int b=bN(k,(int)(ce1[i].y/R),(int)(ce1[i].x/R));
                    for(int c=0;c<cB[b];c++){
                        int j=mB[b][c];
                        if(i<j && sqrt(bc_l(ce1[i].x-ce1[j].x)*bc_l(ce1[i].x-ce1[j].x)+bc_l(ce1[i].y-ce1[j].y)*bc_l(ce1[i].y-ce1[j].y))<R){
                            ce1[i].N+=1;
                            ce1[j].N+=1;
                        }
                        else;
                    }
                }
            }
            //ce0.N
            for(int i=0;i<M0 && n>tau_th;i++){
                ce0[i].N=0;
                for(int k=0;k<9;k++){
                    int b=bN(k,(int)(ce0[i].y/R),(int)(ce0[i].x/R));
                    for(int c=0;c<cB[b];c++){
                        int j=mB[b][c];
                        if(sqrt(bc_l(ce0[i].x-ce1[j].x)*bc_l(ce0[i].x-ce1[j].x)+bc_l(ce0[i].y-ce1[j].y)*bc_l(ce0[i].y-ce1[j].y))<R)
                            ce0[i].N+=1;
                        else;
                    }
                }
            }

            //reproduction and death
            int j=0;
            Cell ce2[Mj];
            for(int i=0;i<M;i++){
                if(genrand_real2()<(p(ce1[i].N)+q0)*dt){
                    //p
                    if(genrand_real2()<p(ce1[i].N)/(p(ce1[i].N)+q0)){
                        ce2[j]=ce1[i];
                        ce2[j+1]=ce1[i];
                        j+=2;
                        //n==tau0[9]
                        if(J0==1){
                            ce0[M0]=ce1[i];
                            ce0[M0].th=th0;
                            M0+=1;
                        }
                        else;
                    }
                    //q
                    else;
                }
                //r
                else{
                    ce2[j]=ce1[i];
                    j+=1;
                }
            }
            M=j;
            for(int i=0;i<M;i++)
                ce1[i]=ce2[i];

            //time-development
            for(int i=0;i<M;i++){
                ce1[i].x=bc_r(ce1[i].x+sqrt(2.0*Dr*dt)*Nrand(0)+v0*cos(ce1[i].th)*dt);
                ce1[i].y=bc_r(ce1[i].y+sqrt(2.0*Dr*dt)*Nrand(1)+v0*sin(ce1[i].th)*dt);
                ce1[i].th+=sqrt(2.0*Dth*dt)*Nrand(0);
            }
        }
        else;

        //n==tau0[9]
        for(int k=0;k<10 && J0==1;k++){
            for(int i=0;i<M0;i++)
                p_[i][k]=0.0;
            nq[k]=0;
        }
        //n<tau0[9]
        for(int i=0;i<M0 && n>=tau_th && J0==0;i++){
            ce0[i].x=bc_r(ce0[i].x+v0*cos(ce0[i].th)*dt);
            ce0[i].y=bc_r(ce0[i].y+v0*sin(ce0[i].th)*dt);
            for(int k=0;k<10;k++){
                if(n<tau0[k]){
                    p_[i][k]+=p(ce0[i].N);
                    if(i==0)
                        nq[k]+=1;
                    else;
                }
                else;
            }
        }

        if(n==tau0[9]-1){
            //output, 12l
            for(int i=0;i<M0 && n>=tau_th;i++)
                printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",factor,seed,p_[i][0]/nq[0],p_[i][1]/nq[1],p_[i][2]/nq[2],p_[i][3]/nq[3],p_[i][4]/nq[4],p_[i][5]/nq[5],p_[i][6]/nq[6],p_[i][7]/nq[7],p_[i][8]/nq[8],p_[i][9]/nq[9]);

            //global direction
            double th_x=0.0,th_y=0.0;
            for(int i=0;i<M;i++){
                th_x+=cos(ce1[i].th);
                th_y+=sin(ce1[i].th);
            }
            th0=atan2(th_y,th_x)+M_PI*factor/4.0;

            //update tau0 and T0
            for(int k=0;k<10;k++)
                tau0[k]=(int)((T0+(k+1)/q0)/dt);
            T0+=10.0/q0;

            //initialisation
            M0=0;
            J0=1;
        }
        else
            J0=0;
    }

    time_t t_end=time(NULL);
    printf("%lf\n",difftime(t_end,t_start));

    return 0;
}

//function
double p(int N){
    double p1=p0*(1.0-(double)N/Ns);
    if(p1>0.0)
        return p1;
    else
        return 0.0;
}
double bc_r(double u){
    if(u<0.0)
        return u+L;
    else if(u>=L)
        return u-L;
    else
        return u;
}
double bc_l(double l){
    if(l>=0.5*L)
        return l-L;
    else if(l<-0.5*L)
        return l+L;
    else
        return l;
}
int bc_b(int b){
    int nB=(int)(L/R);
    if(b==-1)
        return nB-1;
    else if(b==nB)
        return 0;
    else
        return b;
}
int bN(int k,int j,int i){
    int nB=(int)(L/R);
         if(k==0) return bc_b(j-1)*nB+bc_b(i-1);
    else if(k==1) return bc_b(j-1)*nB+bc_b(i*1);
    else if(k==2) return bc_b(j-1)*nB+bc_b(i+1);
    else if(k==3) return bc_b(j*1)*nB+bc_b(i-1);
    else if(k==4) return bc_b(j*1)*nB+bc_b(i*1);
    else if(k==5) return bc_b(j*1)*nB+bc_b(i+1);
    else if(k==6) return bc_b(j+1)*nB+bc_b(i-1);
    else if(k==7) return bc_b(j+1)*nB+bc_b(i*1);
    else if(k==8) return bc_b(j+1)*nB+bc_b(i+1);
    else          return -1;
}
double Nrand(int n){
    if(n%2==0)
        return sqrt(-2.0*log(genrand_real3()))*cos(2.0*M_PI*genrand_real3());
    else
        return sqrt(-2.0*log(genrand_real3()))*sin(2.0*M_PI*genrand_real3());
}
double ConvertToParameter(int f){
    double e=(int)(f/100);
    for(int i=0;i<f%100;i++)
        e*=0.1;
    return e;
}

