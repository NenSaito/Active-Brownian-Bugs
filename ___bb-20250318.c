//___bb-20250318.c
//% gcc ___bb-20250318.c -o out -lm -O3
//% time ./out seed e_L e_Dr e_v0 e_Dth > _bb-20250318-data-'seed'-'e_L'-'e_Dr'-'e_v0'-'e_Dth'.txt
//$ gcc -std=c99 ___bb-20250318.c -o out -lm -O3
//$ for a in {0..4};do for b in 10;do for c in 206 506 105 205 505 104 204;do for d in 204 504 103 203 503 102 202;do for e in 203;do bsub -q "jobLimit52" "./out $a $b $c $d $e > _bb-20250419-data-'$a'-'$b'-'$c'-'$d'-'$e'.txt";done;done;done;done;done
//e_L in {10,15,20,30,40}
//e_Dr in {206,506,105,205,505,104,204}
//e_v0 in {204,504,103,203,503,102,202}
//e_Dth in {104,204,404,103,203,403,102}

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
#define Tau (1000000)
#define Tau_relax (500000)
#define interval (1000)
#define R_inter (0.1)
#define Ns (50)
#define p0 (0.85)
#define q0 (0.15)
#define rh0 (1500.0)
#define dr (0.001)

typedef struct{
    double x;
    double y;
    double theta;
    int N;
    int loc;
    double th0;
    int fl;
}Cell;

//function
double p(int);
double BC_u(double,double);
double BC_l(double,double);
int U_loc_bool(int,int,double);
double N_rand(int);
double ConvertToParameter(int);

int main(int argc,char* argv[]){
    //command line argument
    int seed=strtol(argv[1],NULL,10);
    init_genrand(seed);
    double L_square=0.1*strtol(argv[2],NULL,10);
    double Dr=ConvertToParameter(strtol(argv[3],NULL,10)),v0=ConvertToParameter(strtol(argv[4],NULL,10)),Dth=ConvertToParameter(strtol(argv[5],NULL,10));

    //define
    int M0=(int)(rh0*L_square*L_square),Mj=2*M0,Nr=(int)(0.5*L_square/dr);
    int M_now=M0;
    double M_=0.0,g_[Nr],S_=0.0,s_=0.0,n_relax=0.0;
    for(int k=0;k<Nr;k++)
        g_[k]=0.0;

    int flag_malloc=0;
    Cell cell_1[Mj],cell_2[Mj];
//    //malloc
//    int flag_malloc=0;
//    Cell* cell_1[Mj];
//    Cell* cell_2[Mj];
//    for(int i=0;i<Mj && flag_malloc==0;i++){
//        cell_1[i]=(Cell*)malloc(sizeof(Cell));
//        cell_2[i]=(Cell*)malloc(sizeof(Cell));
//        if(cell_1[i]==NULL || cell_2[i]==NULL)
//            flag_malloc=1;
//        else;
//    }
    if(flag_malloc==0){

        //output
        printf("%d %lf %lf %lf %lf %d %d %d %lf %lf %d %lf %lf %d %lf %d\n"
                ,seed,Dr,v0,Dth,dt,Tau,Tau_relax,interval,L_square,R_inter,Ns,p0,q0,M0,dr,Nr);

        //time-evolution
        for(int tau=0;tau<=Tau;tau++){
            //initialisation
            if(tau==0){
                for(int i=0;i<M_now;i++){
                    cell_1[i].x=L_square*genrand_real2();
                    cell_1[i].y=L_square*genrand_real2();
                    cell_1[i].theta=2.0*M_PI*genrand_real2();
                    cell_1[i].N=0;
                    cell_1[i].loc=(int)(cell_1[i].x/L_square)+(int)(L_square/R_inter)*(int)(cell_1[i].y/L_square);
                    cell_1[i].th0=cell_1[i].theta;
                    cell_1[i].fl=i;
                }
            }
            else if(tau>0){
                //cell number
                int M_new=0;
                for(int i=0;i<M_now;i++)
                    cell_1[i].N=0;
                for(int i=0;i<M_now;i++){
                    for(int j=i+1;j<M_now;j++){
                        if(U_loc_bool(cell_1[j].loc,cell_1[i].loc,L_square)){
                            if(BC_l(cell_1[j].x-cell_1[i].x,L_square)*BC_l(cell_1[j].x-cell_1[i].x,L_square)+BC_l(cell_1[j].y-cell_1[i].y,L_square)*BC_l(cell_1[j].y-cell_1[i].y,L_square)<R_inter*R_inter){
                                cell_1[i].N+=1;
                                cell_1[j].N+=1;
                            }
                            else;
                        }
                        else;
                    }
                }
                for(int i=0;i<M_now;i++){
                    if(genrand_real2()<(p(cell_1[i].N)+q0)*dt){
                        if(genrand_real2()<p(cell_1[i].N)/(p(cell_1[i].N)+q0)){
                            cell_2[M_new].x=cell_2[M_new+1].x=cell_1[i].x;
                            cell_2[M_new].y=cell_2[M_new+1].y=cell_1[i].y;
                            cell_2[M_new].theta=cell_2[M_new+1].theta=cell_1[i].theta;
                            cell_2[M_new].N=cell_2[M_new+1].N=cell_1[i].N;
                            cell_2[M_new].loc=cell_2[M_new+1].loc=cell_1[i].loc;
                            cell_2[M_new].th0=cell_2[M_new+1].th0=cell_1[i].th0;
                            cell_2[M_new].fl=cell_2[M_new+1].fl=cell_1[i].fl;
                            M_new+=2;
                        }
                        else;
                    }
                    else{
                        cell_2[M_new].x=cell_1[i].x;
                        cell_2[M_new].y=cell_1[i].y;
                        cell_2[M_new].theta=cell_1[i].theta;
                        cell_2[M_new].N=cell_1[i].N;
                        cell_2[M_new].loc=cell_1[i].loc;
                        cell_2[M_new].th0=cell_1[i].th0;
                        cell_2[M_new].fl=cell_1[i].fl;
                        M_new+=1;
                    }
                }
                M_now=M_new;

                //variable
                for(int i=0;i<M_now;i++){
                    cell_1[i].x=BC_u(cell_2[i].x+sqrt(2.0*Dr*dt)*N_rand(1)+v0*cos(cell_2[i].theta)*dt,L_square);
                    cell_1[i].y=BC_u(cell_2[i].y+sqrt(2.0*Dr*dt)*N_rand(2)+v0*sin(cell_2[i].theta)*dt,L_square);
                    cell_1[i].theta=cell_2[i].theta+sqrt(2.0*Dth*dt)*N_rand(0);
                    cell_1[i].N=cell_2[i].N;
                    cell_1[i].loc=(int)(cell_2[i].x/L_square)+(int)(L_square/R_inter)*(int)(cell_2[i].y/L_square);
                    cell_1[i].th0=cell_2[i].th0;
                    cell_1[i].fl=cell_2[i].fl;
                }

            }
            else;

            if(tau>=Tau_relax)
                n_relax+=1.0;
            else;

            //cell number
            if(tau>=Tau_relax)
                M_+=M_now;
            else;

            //pair correlation
            double g[Nr],C_pair=0.0;
            for(int k=0;k<Nr;k++)
                g[k]=0.0;
            for(int i=0;i<M_now;i++){
                for(int j=i+1;j<M_now;j++){
                    double r=sqrt(BC_l(cell_1[j].x-cell_1[i].x,L_square)*BC_l(cell_1[j].x-cell_1[i].x,L_square)+BC_l(cell_1[j].y-cell_1[i].y,L_square)*BC_l(cell_1[j].y-cell_1[i].y,L_square));
                    if(r<0.5*L_square){
                        g[(int)(r/dr)]+=1.0;
                        C_pair+=1.0;
                    }
                    else;
                }
            }
            double g_max=1.0;
            int J_g=0;
            for(int k=0;k<Nr;k++){
                g[k]*=0.25*L_square*L_square/(C_pair*(2.0*k+1.0)*dr*dr);
                if(J_g==0 && g[k]<1.0)
                    J_g=1;
                else if(J_g==1 && g_max<=g[k])
                    g_max=g[k];
                else;
                if(tau>=Tau_relax)
                    g_[k]+=g[k];
                else;
            }

            //orientational order
            double x_S,y_S,x_s,y_s;
            x_S=y_S=x_s=y_s=0.0;
            for(int i=0;i<M_now;i++){
                x_S+=cos(cell_1[i].theta);
                y_S+=sin(cell_1[i].theta);
                x_s+=cos(cell_1[i].th0);
                y_s+=sin(cell_1[i].th0);
            }
            double S1=sqrt(x_S*x_S+y_S*y_S)/M_now;
            double s1=sqrt(x_s*x_s+y_s*y_s)/M_now;
            if(tau>=Tau_relax){
                S_+=S1;
                s_+=s1;
            }
            else;

            //family line
            int Nf=(int)(M_now>0);
            for(int i=1;i<M_now && Nf>0;i++){
                int J_f=0;
                for(int j=i-1;j>=0 && J_f==0;j--){
                    if(cell_1[j].fl==cell_1[i].fl)
                        J_f=1;
                    else;
                }
                Nf+=1-J_f;
            }

            //output
            if(tau%interval==0){
                printf("%d %d %lf %lf %lf %d %d\n",tau,M_now,g_max,S1,s1,Nf,-1);
                for(int i=0;i<M_now;i++){
                    printf("%d %lf %lf %lf %d %lf %d\n",tau,cell_1[i].x,cell_1[i].y,cell_1[i].theta,cell_1[i].N,cell_1[i].th0,cell_1[i].fl);
                }
            }
            else;
        }

        //output
        printf("%lf %lf %lf",M_/n_relax,S_/n_relax,s_/n_relax);
        double g_max=1.0;
        int J_g=0;
        for(int k=0;k<Nr;k++){
            printf(" %lf",g_[k]/n_relax);
            if(J_g==0 && g_[k]/n_relax<1.0)
                J_g=1;
            else if(J_g==1 && g_max<=g_[k]/n_relax)
                g_max=g_[k]/n_relax;
            else;
        }
        printf(" %lf\n",g_max);

    }
    else
        printf("fail in \"malloc\"\n");

//    //free
//    for(int i=0;i<Mj;i++){
//        free(cell_1[i]);
//        free(cell_2[i]);
//    }

    return 0;
}

//function
double p(int N){
    double p1=p0*(1.0-(double)(N)/Ns);
    if(p1>0.0)
        return p1;
    else
        return 0.0;
}
double BC_u(double u,double L){
    if(u<0.0)
        return u+L;
    else if(u>=L)
        return u-L;
    else
        return u;
}
double BC_l(double l,double L){
    if(l>0.5*L)
        return L-l;
    else if(l<-0.5*L)
        return L+l;
    else
        return l;
}
int U_loc_bool(int loc_j,int loc_i,double L){
    int N=(int)(L/R_inter);
    if(loc_j==+(loc_i> N-1 && loc_i%N!=0)*(loc_i-N-1)
              +(loc_i<=N-1 && loc_i%N!=0)*(loc_i+N*N-N-1)
              +(loc_i> N-1 && loc_i%N==0)*(loc_i-1)
              +(loc_i<=N-1 && loc_i%N==0)*(N*N-1))
        return 1;
    else if(loc_j==+(loc_i> N-1)*(loc_i-N)
                   +(loc_i<=N-1)*(loc_i+N*N-N))
        return 1;
    else if(loc_j==+(loc_i> N-1 && loc_i%N!=N-1)*(loc_i-N+1)
                   +(loc_i<=N-1 && loc_i%N!=N-1)*(loc_i+N*N-N+1)
                   +(loc_i> N-1 && loc_i%N==N-1)*(loc_i-2*N+1)
                   +(loc_i<=N-1 && loc_i%N==N-1)*(N*N-N))
        return 1;
    else if(loc_j==+(loc_i%N!=0)*(loc_i-1)
                   +(loc_i%N==0)*(loc_i+N-1))
        return 1;
    else if(loc_j==loc_i)
        return 1;
    else if(loc_j==+(loc_i%N!=N-1)*(loc_i+1)
                   +(loc_i%N==N-1)*(loc_i-N+1))
        return 1;
    else if(loc_j==+(loc_i< N*N-N && loc_i%N!=0)*(loc_i+N-1)
                   +(loc_i>=N*N-N && loc_i%N!=0)*(loc_i-N*N+N-1)
                   +(loc_i< N*N-N && loc_i%N==0)*(loc_i+2*N-1)
                   +(loc_i>=N*N-N && loc_i%N==0)*(N-1))
        return 1;
    else if(loc_j==+(loc_i< N*N-N)*(loc_i+N)
                   +(loc_i>=N*N-N)*(loc_i-N*N+N))
        return 1;
    else if(loc_j==+(loc_i< N*N-N && loc_i%N!=N-1)*(loc_i+N+1)
                   +(loc_i>=N*N-N && loc_i%N!=N-1)*(loc_i-N*N+N+1)
                   +(loc_i< N*N-N && loc_i%N==N-1)*(loc_i+1)
                   +(loc_i>=N*N-N && loc_i%N==N-1)*(0))
        return 1;
    else
        return 0;
}
double N_rand(int n){
    if(n%2==1)
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