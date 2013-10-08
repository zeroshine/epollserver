#include <iostream>
#include<fstream>
#include <complex>
#include <stdio.h>
#include <math.h>
using namespace std;

double getM(complex<double> Vsa,complex<double> Vsb,complex<double> Vsc,complex<double> Isa,complex<double> Isb,complex<double> Isc,complex<double> Vra,complex<double> Vrb,complex<double> Vrc,complex<double> Ira,complex<double> Irb,complex<double> Irc,double R0,double R1,double L0,double L1,double C0,double C1,double Lm) {
    
    double T[3][3],iT[3][3],R[3][3],C[3][3],L[3][3];
    int i, j, k;
    double w=120*3.1415926535;
           
    double Cm=(C0-C1)/3,Cs=C0-2*Cm;
    double Rm=(R0-R1)/3,Rs=R0-2*Rm;
    double Lmu=(L0-L1)/3,Ls=L0-2*Lmu;
    complex<double> Z[3][3],Y[3][3],Z012[3][3],Y012[3][3],r[3][3];    
    complex<double> r0,r1,r2,Zc0,Zc1,Zc2;
    complex<double> invY012[3][3],Zc[3][3];
    complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2,temp;
    
    // 0,1,2 represents for 0,alpha,beta
    complex<double> Vs0=(1/sqrt(6))*(Vsa+Vsb+Vsc);
    complex<double> Vs1=(sqrt(2)/sqrt(3))*Vsa+(-1/sqrt(6))*Vsb+(-1/sqrt(6))*Vsc;
    complex<double> Vs2=(1/sqrt(2))*Vsb+(-1/sqrt(2))*Vsc;
    complex<double> Is0=(1/sqrt(6))*(Isa+Isb+Isc);
    complex<double> Is1=(sqrt(2)/sqrt(3))*Isa+(-1/sqrt(6))*Isb+(-1/sqrt(6))*Isc;
    complex<double> Is2=(1/sqrt(2))*Isb+(-1/sqrt(2))*Isc;
    complex<double> Vr0=(1/sqrt(6))*(Vra+Vrb+Vrc);
    complex<double> Vr1=(sqrt(2)/sqrt(3))*Vra+(-1/sqrt(6))*Vrb+(-1/sqrt(6))*Vrc;
    complex<double> Vr2=(1/sqrt(2))*Vrb+(-1/sqrt(2))*Vrc;
    complex<double> Ir0=(1/sqrt(6))*(Ira+Irb+Irc);
    complex<double> Ir1=(sqrt(2)/sqrt(3))*Ira+(-1/sqrt(6))*Irb+(-1/sqrt(6))*Irc;
    complex<double> Ir2=(1/sqrt(2))*Irb+(-1/sqrt(2))*Irc;
    ///////////////////////////////////////////////////////////////////////////////////
    
    for(i=0;i<=2;i++){
          for(j=0;j<=2;j++){
               if(i==j){
                        R[i][j]=Rs;
                        L[i][j]=Ls;
                        C[i][j]=Cs;}
               else {
                    R[i][j]=Rm;
                    L[i][j]=Lmu;
                    C[i][j]=Cm;}            
          }                                
    }
    
    //Z=R+j*w*L ; Y=j*w*C
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
              (Z[i][j]).real()= R[i][j];
              (Z[i][j]).imag()= w*L[i][j];  
              (Y[i][j]).real()= 0;
              (Y[i][j]).imag()= w*C[i][j];    
        }       
    }    
    
    //Clarke Matrix
    T[0][0]=sqrt(2)/sqrt(3);T[0][1]=sqrt(2)/sqrt(3);T[0][2]=0;
    T[1][0]=sqrt(2)/sqrt(3);T[1][1]=-1/sqrt(6);T[1][2]=1/sqrt(2);
    T[2][0]=sqrt(2)/sqrt(3);T[2][1]=-1/sqrt(6);T[2][2]=-1/sqrt(2);

    //inverse Clarke Matrix
    iT[0][0]=1/sqrt(6);iT[0][1]=1/sqrt(6);iT[0][2]=1/sqrt(6);
    iT[1][0]=sqrt(2)/sqrt(3);iT[1][1]=-1/sqrt(6);iT[1][2]=-1/sqrt(6);
    iT[2][0]=0;iT[2][1]=1/sqrt(2);iT[2][2]=-1/sqrt(2);
    
    //Z012=iT*Z*T
    complex<double> Z_X_T[3][3];
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z_X_T[i][j]+=Z[i][k]*T[k][j];
                           
               }            
        }      
    }   
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z012[i][j]+=iT[i][k]*Z_X_T[k][j];
                           
               }            
        }       
    } 
    
    //Y012=iT*Y*T
    complex<double> Y_X_T[3][3];
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Y_X_T[i][j]+=Y[i][k]*T[k][j];
                           
               }            
        }      
    }
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Y012[i][j]+=iT[i][k]*Y_X_T[k][j];
                           
               }            
        }       
    }
    
    //r=sqrt(Z012*Y012);r0=r(1,1),r1=r(2,2),r2=r(3,3) % r1=r2
    complex<double> Z012_X_Y012[3][3];
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z012_X_Y012[i][j]+=Z012[i][k]*Y012[k][j];
                           
               }            
        }      
    }
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {         
                           r[i][j]=sqrt(Z012_X_Y012[i][j]);
                                                
        }      
    }
    r0=r[0][0];
    r1=r[1][1];
    r2=r[2][2];  
    
    //Zc=sqrt(Z0dq/Y0dq);Zc0=Zc(1,1),Zcd=Zc(2,2),Zcq=Zc(3,3) %Zcd=Zcq    
    complex<double> delta;
    delta=(Y012[0][0])*(Y012[1][1])*(Y012[2][2]);
    delta+=(Y012[0][1])*(Y012[1][2])*(Y012[2][0]);
    delta+=(Y012[0][2])*(Y012[1][0])*(Y012[2][1]);
    delta-=(Y012[0][2])*(Y012[1][1])*(Y012[2][0]);
    delta-=(Y012[0][0])*(Y012[1][2])*(Y012[2][1]);
    delta-=(Y012[0][1])*(Y012[2][2])*(Y012[1][0]);
    
    invY012[0][0]=((Y012[1][1])*(Y012[2][2]))/delta;
    invY012[0][0]-=((Y012[1][2])*(Y012[2][1]))/delta;
    
    invY012[0][1]=((Y012[1][2])*(Y012[2][0]))/delta;
    invY012[0][1]-=((Y012[1][0])*(Y012[2][2]))/delta;
    
    invY012[0][2]=((Y012[1][0])*(Y012[2][1]))/delta;
    invY012[0][2]-=((Y012[1][1])*(Y012[2][0]))/delta;
    
    invY012[1][0]=((Y012[2][1])*(Y012[0][2]))/delta;
    invY012[1][0]-=((Y012[2][2])*(Y012[0][1]))/delta;
    
    invY012[1][1]=((Y012[0][0])*(Y012[2][2]))/delta;
    invY012[1][1]-=((Y012[0][2])*(Y012[2][0]))/delta;
    
    invY012[1][2]=((Y012[0][1])*(Y012[2][0]))/delta;
    invY012[1][2]-=((Y012[0][0])*(Y012[2][1]))/delta;
    
    invY012[2][0]=((Y012[0][1])*(Y012[1][2]))/delta;
    invY012[2][0]-=((Y012[0][2])*(Y012[1][1]))/delta;
    
    invY012[2][1]=((Y012[1][0])*(Y012[0][2]))/delta;
    invY012[2][1]-=((Y012[0][0])*(Y012[1][2]))/delta;
    
    invY012[2][2]=((Y012[0][0])*(Y012[1][1]))/delta;
    invY012[2][2]-=((Y012[0][1])*(Y012[1][0]))/delta;
    
    complex<double> Z012_div_Y012[3][3];
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z012_div_Y012[i][j]+=Z012[i][k]*invY012[k][j];                           
               }            
        }       
    }
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {                        
                           Zc[i][j]=sqrt(Z012_div_Y012[i][j]);                                                   
        }       
    }
    
    Zc0=Zc[0][0]; 
    Zc1=Zc[1][1]; 
    Zc2=Zc[2][2];
             
    //M0=(Vs0+Zc0*Is0)*(exp(gamma0*Lm))/2-(Vr0-Zc0*Ir0)/2;
    //Lm=0.001;
    temp=exp((r0*Lm));
    
    
    temp*=(Vs0+Zc0*Is0); 
    temp-=(Vr0-Zc0*Ir0);
    temp/=2;
    M0=temp;
    
    //N0=(Vr0+Zc0*Ir0)/2-(Vs0-Zc0*Is0)*exp(-gamma0*Lm)/2;
    temp=exp(-r0*Lm);
    temp*=-(Vs0-Zc0*Is0);
    temp+=Vr0+Zc0*Ir0;
    temp/=2;
    N0=temp;
    
    //D0=(log(M0/N0)/(2*gamma0*Lm));
    D0=log(M0/N0);
    D0/=(r0*Lm);
    D0/=2;

    
    //M1=(Vs1+Zcd*Is1)*(exp(gammad*Lm))/2-(Vr1-Zcd*Ir1)/2;
    temp=exp(r1*Lm);
    temp*=(Vs1+Zc1*Is1); 
    temp-=(Vr1-Zc1*Ir1);
    temp/=2;
    M1=temp;
    
    //N1=(Vr1+Zcd*Ir1)/2-(Vs1-Zcd*Is1)*exp(-gammad*Lm)/2;
    temp=exp(-r1*Lm);
    temp*=-(Vs1-Zc1*Is1);
    temp+=Vr1+Zc1*Ir1;
    temp/=2;
    N1=temp;
    
    //D1=(log(Md/Nd)/(2*gammad*Lm));
    D1=log(M1/N1);
    D1/=(r1*Lm);
    D1/=2;

    //M2=(Vs1+Zcd*Is1)*(exp(gammad*Lm))/2-(Vr1-Zcd*Ir1)/2;
    temp=exp(r2*Lm);
    temp*=(Vs2+Zc2*Is2); 
    temp-=(Vr2-Zc2*Ir2);
    temp/=2;
    M2=temp;
    
    //N2=(Vr1+Zcd*Ir1)/2-(Vs1-Zcd*Is1)*exp(-gammad*Lm)/2;
    temp=exp(-r2*Lm);
    temp*=-(Vs2-Zc2*Is2);
    temp+=Vr2+Zc2*Ir2;
    temp/=2;
    N2=temp;
    
    //D2=(log(Md/Nd)/(2*gammad*Lm));
    D2=log(M2/N2);
    D2/=(r2*Lm);
    D2/=2;  
        
    
    if(abs(M0)>=abs(M1) && abs(M0)>=abs(M2)){
        return abs(M0);
    }
    else if(abs(M1)>=abs(M2)){
        return abs(M1);
    }
    else return abs(M2);
}


















double getD(complex<double> Vsa,complex<double> Vsb,complex<double> Vsc,complex<double> Isa,complex<double> Isb,complex<double> Isc,complex<double> Vra,complex<double> Vrb,complex<double> Vrc,complex<double> Ira,complex<double> Irb,complex<double> Irc,double R0,double R1,double L0,double L1,double C0,double C1,double Lm) {
    
    
    
    
    double T[3][3],iT[3][3],R[3][3],C[3][3],L[3][3];
    int i, j, k;
    double pi=3.14159265;
    double w=120*pi;
    complex<double> ang;
    ang.real()=0;
    ang.imag()=2*pi/3.0;
    complex<double> a=exp(ang);
    complex<double> Z0,Z1,Z2,Y0,Y1,Y2;
    complex<double> r0,r1,r2,Zc0,Zc1,Zc2;
    complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2,temp;
    
        
    complex<double> Vs0=(Vsa+Vsb+Vsc)/3.0;
    complex<double> Vs1=(Vsa+a*Vsb+a*a*Vsc)/3.0;
    complex<double> Vs2=(Vsa+a*a*Vsb+a*Vsc)/3.0;
    complex<double> Is0=(Isa+Isb+Isc)/3.0;
    complex<double> Is1=(Isa+a*Isb+a*a*Isc)/3.0;
    complex<double> Is2=(Isa+a*a*Isb+a*Isc)/3.0;
    complex<double> Vr0=(Vra+Vrb+Vrc)/3.0;
    complex<double> Vr1=(Vra+a*Vrb+a*a*Vrc)/3.0;
    complex<double> Vr2=(Vra+a*a*Vrb+a*Vrc)/3.0;
    complex<double> Ir0=(Ira+Irb+Irc)/3.0;
    complex<double> Ir1=(Ira+a*Irb+a*a*Irc)/3.0;
    complex<double> Ir2=(Ira+a*a*Irb+a*Irc)/3.0;
    
    ///////////////////////////////////////////////////////////////////////////////////
    Z0.real()=R0;
    Z0.imag()=w*L0;
    Z1.real()=R1;
    Z1.imag()=w*L1;
    Z2.real()=R1;
    Z2.imag()=w*L1;    
    Y0.real()=0;
    Y0.imag()=w*C0;
    Y1.real()=0;
    Y1.imag()=w*C1;
    Y2.real()=0;
    Y2.imag()=w*C1;
    
    r0=sqrt(Z0*Y0);
    r1=sqrt(Z1*Y1);
    r2=sqrt(Z1*Y1);
    
    Zc0=sqrt(Z0/Y0);
    printf("--------------   %f +%f i  <<<<< Zc0\n",real(Zc0),imag(Zc0));
    
    Zc0=sqrt(((Vs0*Vs0)-(Vr0*Vr0))/((Is0*Is0)-(Ir0*Ir0)));
    printf("--------------   %f +%f i  <<<<< Zc0\n",real(Zc0),imag(Zc0));
    Zc1=sqrt(Z1/Y1); 
    printf("--------------   %f +%f i  <<<<< Zc1\n",real(Zc1),imag(Zc1));
    Zc1=sqrt(((Vs1*Vs1)-(Vr1*Vr1))/((Is1*Is1)-(Ir1*Ir1))); 
    printf("--------------   %f +%f i  <<<<< Zc1\n",real(Zc1),imag(Zc1));
    Zc2=sqrt(Z1/Y1);  
    
    M0=-((Vr0+Zc0*Ir0)-(Vs0-Zc0*Is0)*exp(r0*Lm))/2.0;
    N0=-((Vs0+Zc0*Is0)/exp(r0*Lm)-(Vr0-Zc0*Ir0))/2.0;
    D0=log(M0/N0)/(2.0*r0*Lm);
    
    N1=-((Vs1+Zc1*Is1)/exp(r1*Lm)-(Vr1-Zc1*Ir1))/2.0;    
    M1=-((Vr1+Zc1*Ir1)-(Vs1-Zc1*Is1)*exp(r1*Lm))/2.0;
    D1=log(M1/N1)/(2.0*r1*Lm);

    M2=-((Vr2+Zc2*Ir2)-(Vs2-Zc2*Is2)*exp(r2*Lm))/2.0;
    N2=-((Vs2+Zc2*Is2)/exp(r2*Lm)-(Vr2-Zc2*Ir2))/2.0;
    D2=log(M2/N2)/(2.0*r1*Lm);
    
    printf("%f +%f i  <<<<< Vs1\n",real(Vs1),imag(Vs1));
    printf("%f +%f i  <<<<< Vr1\n",real(Vr1),imag(Vr1));   
    printf("%f +%f i  <<<<< Is1\n",real(Is1),imag(Is1)); 
    printf("%f +%f i  <<<<< Ir1\n",real(Ir1),imag(Ir1));
    printf("%f +%f i  <<<<< Zc1\n",real(Zc1),imag(Zc1));
    printf("%f +%f i  <<<<< r1\n",real(r1),imag(r1));
    printf("%f +%f i  <<<<< M1\n",real(M1),imag(M1));
    printf("%f +%f i  <<<<< N1\n",real(N1),imag(N1));
    printf("%f +%f i  <<<<< D0\n",real(D0),imag(D0)); 
    printf("%f +%f i  <<<<< D1\n",real(D1),imag(D1));   
    printf("%f +%f i  <<<<< D2\n",real(D2),imag(D2));
    
    if(abs(M0)<=abs(M1) && abs(M0)<=abs(M2)){
        return abs(D0);
    }
    else if(abs(M1)<=abs(M2)){
        return abs(D1);
    }
    else return abs(D2);
}


















string getType(complex<double> Vsa,complex<double> Vsb,complex<double> Vsc,complex<double> Isa,complex<double> Isb,complex<double> Isc,complex<double> Vra,complex<double> Vrb,complex<double> Vrc,complex<double> Ira,complex<double> Irb,complex<double> Irc,double R0,double R1,double L0,double L1,double C0,double C1,double Lm) {
    
    double T[3][3],iT[3][3],R[3][3],C[3][3],L[3][3];
    int i, j, k;
    double w=120*3.1415926535;
           
    double Cm=(C0-C1)/3,Cs=C0-2*Cm;
    double Rm=(R0-R1)/3,Rs=R0-2*Rm;
    double Lmu=(L0-L1)/3,Ls=L0-2*Lmu;
    complex<double> Z[3][3],Y[3][3],Z012[3][3],Y012[3][3],r[3][3];    
    complex<double> r0,r1,r2,Zc0,Zc1,Zc2;
    complex<double> invY012[3][3],Zc[3][3];
    complex<double> Ma0,Ma1,Ma2,Mb0,Mb1,Mb2,Mc0,Mc1,Mc2,temp;
    double Threshold=300;
    
    //A-basis   0,1,2 represents for 0,alpha,beta
    complex<double> Vs_a0=(1/sqrt(6))*(Vsa+Vsb+Vsc);
    complex<double> Vs_a1=(sqrt(2)/sqrt(3))*Vsa+(-1/sqrt(6))*Vsb+(-1/sqrt(6))*Vsc;
    complex<double> Vs_a2=(1/sqrt(2))*Vsb+(-1/sqrt(2))*Vsc;
    complex<double> Is_a0=(1/sqrt(6))*(Isa+Isb+Isc);
    complex<double> Is_a1=(sqrt(2)/sqrt(3))*Isa+(-1/sqrt(6))*Isb+(-1/sqrt(6))*Isc;
    complex<double> Is_a2=(1/sqrt(2))*Isb+(-1/sqrt(2))*Isc;
    complex<double> Vr_a0=(1/sqrt(6))*(Vra+Vrb+Vrc);
    complex<double> Vr_a1=(sqrt(2)/sqrt(3))*Vra+(-1/sqrt(6))*Vrb+(-1/sqrt(6))*Vrc;
    complex<double> Vr_a2=(1/sqrt(2))*Vrb+(-1/sqrt(2))*Vrc;
    complex<double> Ir_a0=(1/sqrt(6))*(Ira+Irb+Irc);
    complex<double> Ir_a1=(sqrt(2)/sqrt(3))*Ira+(-1/sqrt(6))*Irb+(-1/sqrt(6))*Irc;
    complex<double> Ir_a2=(1/sqrt(2))*Irb+(-1/sqrt(2))*Irc;
    ///////////////////////////////////////////////////////////////////////////////////
    
    for(i=0;i<=2;i++){
          for(j=0;j<=2;j++){
               if(i==j){
                        R[i][j]=Rs;
                        L[i][j]=Ls;
                        C[i][j]=Cs;}
               else {
                    R[i][j]=Rm;
                    L[i][j]=Lmu;
                    C[i][j]=Cm;}            
          }                                
    }
    
    //Z=R+j*w*L ; Y=j*w*C
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
              (Z[i][j]).real()= R[i][j];
              (Z[i][j]).imag()= w*L[i][j];  
              (Y[i][j]).real()= 0;
              (Y[i][j]).imag()= w*C[i][j];    
        }       
    }    
    
    //Clarke Matrix
    T[0][0]=sqrt(2)/sqrt(3);T[0][1]=sqrt(2)/sqrt(3);T[0][2]=0;
    T[1][0]=sqrt(2)/sqrt(3);T[1][1]=-1/sqrt(6);T[1][2]=1/sqrt(2);
    T[2][0]=sqrt(2)/sqrt(3);T[2][1]=-1/sqrt(6);T[2][2]=-1/sqrt(2);

    //inverse Clarke Matrix
    iT[0][0]=1/sqrt(6);iT[0][1]=1/sqrt(6);iT[0][2]=1/sqrt(6);
    iT[1][0]=sqrt(2)/sqrt(3);iT[1][1]=-1/sqrt(6);iT[1][2]=-1/sqrt(6);
    iT[2][0]=0;iT[2][1]=1/sqrt(2);iT[2][2]=-1/sqrt(2);
    
    //Z012=iT*Z*T
    complex<double> Z_X_T[3][3];
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z_X_T[i][j]+=Z[i][k]*T[k][j];
                           
               }            
        }      
    }   
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z012[i][j]+=iT[i][k]*Z_X_T[k][j];
                           
               }            
        }       
    } 
    
    //Y012=iT*Y*T
    complex<double> Y_X_T[3][3];
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Y_X_T[i][j]+=Y[i][k]*T[k][j];
                           
               }            
        }      
    }
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Y012[i][j]+=iT[i][k]*Y_X_T[k][j];
                           
               }            
        }       
    }
    
    //r=sqrt(Z012*Y012);r0=r(1,1),r1=r(2,2),r2=r(3,3) % r1=r2
    complex<double> Z012_X_Y012[3][3];
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z012_X_Y012[i][j]+=Z012[i][k]*Y012[k][j];
                           
               }            
        }      
    }
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {         
                           r[i][j]=sqrt(Z012_X_Y012[i][j]);
                                                
        }      
    }
    r0=r[0][0];
    r1=r[1][1];
    r2=r[2][2];  
    
    //Zc=sqrt(Z0dq/Y0dq);Zc0=Zc(1,1),Zcd=Zc(2,2),Zcq=Zc(3,3) %Zcd=Zcq    
    complex<double> delta;
    delta=(Y012[0][0])*(Y012[1][1])*(Y012[2][2]);
    delta+=(Y012[0][1])*(Y012[1][2])*(Y012[2][0]);
    delta+=(Y012[0][2])*(Y012[1][0])*(Y012[2][1]);
    delta-=(Y012[0][2])*(Y012[1][1])*(Y012[2][0]);
    delta-=(Y012[0][0])*(Y012[1][2])*(Y012[2][1]);
    delta-=(Y012[0][1])*(Y012[2][2])*(Y012[1][0]);
    
    invY012[0][0]=((Y012[1][1])*(Y012[2][2]))/delta;
    invY012[0][0]-=((Y012[1][2])*(Y012[2][1]))/delta;
    
    invY012[0][1]=((Y012[1][2])*(Y012[2][0]))/delta;
    invY012[0][1]-=((Y012[1][0])*(Y012[2][2]))/delta;
    
    invY012[0][2]=((Y012[1][0])*(Y012[2][1]))/delta;
    invY012[0][2]-=((Y012[1][1])*(Y012[2][0]))/delta;
    
    invY012[1][0]=((Y012[2][1])*(Y012[0][2]))/delta;
    invY012[1][0]-=((Y012[2][2])*(Y012[0][1]))/delta;
    
    invY012[1][1]=((Y012[0][0])*(Y012[2][2]))/delta;
    invY012[1][1]-=((Y012[0][2])*(Y012[2][0]))/delta;
    
    invY012[1][2]=((Y012[0][1])*(Y012[2][0]))/delta;
    invY012[1][2]-=((Y012[0][0])*(Y012[2][1]))/delta;
    
    invY012[2][0]=((Y012[0][1])*(Y012[1][2]))/delta;
    invY012[2][0]-=((Y012[0][2])*(Y012[1][1]))/delta;
    
    invY012[2][1]=((Y012[1][0])*(Y012[0][2]))/delta;
    invY012[2][1]-=((Y012[0][0])*(Y012[1][2]))/delta;
    
    invY012[2][2]=((Y012[0][0])*(Y012[1][1]))/delta;
    invY012[2][2]-=((Y012[0][1])*(Y012[1][0]))/delta;
    
    complex<double> Z012_div_Y012[3][3];
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {
               for(k=0;k<=2;k++) {           
                           Z012_div_Y012[i][j]+=Z012[i][k]*invY012[k][j];                           
               }            
        }       
    }
    
    for (i =0 ; i<=2; i++) {
        for(j=0;j<=2;j++) {                        
                           Zc[i][j]=sqrt(Z012_div_Y012[i][j]);                                                   
        }       
    }
    
    Zc0=Zc[0][0]; 
    Zc1=Zc[1][1]; 
    Zc2=Zc[2][2];
             
    //Ma0=(Vs_a0+Zc0*Is_a0)*(exp(gamma0*Lm))/2-(Vr_a0-Zc0*Ir_a0)/2;
    temp=exp((r0*Lm));    
    temp*=(Vs_a0+Zc0*Is_a0); 
    temp-=(Vr_a0-Zc0*Ir_a0);
    temp/=2;
    Ma0=temp;
        
    //Ma1=(Vs_a1+Zcd*Is_a1)*(exp(gammad*Lm))/2-(Vr_a1-Zcd*Ir_a1)/2;
    temp=exp(r1*Lm);
    temp*=(Vs_a1+Zc1*Is_a1); 
    temp-=(Vr_a1-Zc1*Ir_a1);
    temp/=2;
    Ma1=temp;

    //Ma2=(Vs_a1+Zcd*Is_a1)*(exp(gammad*Lm))/2-(Vr_a1-Zcd*Ir_a1)/2;
    temp=exp(r2*Lm);
    temp*=(Vs_a2+Zc2*Is_a2); 
    temp-=(Vr_a2-Zc2*Ir_a2);
    temp/=2;
    Ma2=temp;
    
    //B-basis   0,1,2 represents for 0,alpha,beta
    complex<double> Vs_b0=(1/sqrt(6))*(Vsa+Vsb+Vsc);
    complex<double> Vs_b1=(sqrt(2)/sqrt(3))*Vsb+(-1/sqrt(6))*Vsc+(-1/sqrt(6))*Vsa;
    complex<double> Vs_b2=(1/sqrt(2))*Vsc+(-1/sqrt(2))*Vsa;
    complex<double> Is_b0=(1/sqrt(6))*(Isa+Isb+Isc);
    complex<double> Is_b1=(sqrt(2)/sqrt(3))*Isb+(-1/sqrt(6))*Isc+(-1/sqrt(6))*Isa;
    complex<double> Is_b2=(1/sqrt(2))*Isc+(-1/sqrt(2))*Isa;
    complex<double> Vr_b0=(1/sqrt(6))*(Vra+Vrb+Vrc);
    complex<double> Vr_b1=(sqrt(2)/sqrt(3))*Vrb+(-1/sqrt(6))*Vrc+(-1/sqrt(6))*Vra;
    complex<double> Vr_b2=(1/sqrt(2))*Vrc+(-1/sqrt(2))*Vra;
    complex<double> Ir_b0=(1/sqrt(6))*(Ira+Irb+Irc);
    complex<double> Ir_b1=(sqrt(2)/sqrt(3))*Irb+(-1/sqrt(6))*Irc+(-1/sqrt(6))*Ira;
    complex<double> Ir_b2=(1/sqrt(2))*Irc+(-1/sqrt(2))*Ira;
    
    //Mb0=(Vs_b0+Zc0*Is_b0)*(exp(gamma0*Lm))/2-(Vr_b0-Zc0*Ir_b0)/2;   
    temp=exp((r0*Lm));       
    temp*=(Vs_b0+Zc0*Is_b0); 
    temp-=(Vr_b0-Zc0*Ir_b0);
    temp/=2;
    Mb0=temp;
        
    //Mb1=(Vs_b1+Zcd*Is_b1)*(exp(gammad*Lm))/2-(Vr_b1-Zcd*Ir_b1)/2;
    temp=exp(r1*Lm);
    temp*=(Vs_b1+Zc1*Is_b1); 
    temp-=(Vr_b1-Zc1*Ir_b1);
    temp/=2;
    Mb1=temp;

    //Mb2=(Vs_b2+Zcd*Is_b2)*(exp(gammad*Lm))/2-(Vr_b2-Zcd*Ir_b2)/2;
    temp=exp(r2*Lm);
    temp*=(Vs_b2+Zc2*Is_b2); 
    temp-=(Vr_b2-Zc2*Ir_b2);
    temp/=2;
    Mb2=temp;
    
    //C-basis   0,1,2 represents for 0,alpha,beta
    complex<double> Vs_c0=(1/sqrt(6))*(Vsa+Vsb+Vsc);
    complex<double> Vs_c1=(sqrt(2)/sqrt(3))*Vsc+(-1/sqrt(6))*Vsa+(-1/sqrt(6))*Vsb;
    complex<double> Vs_c2=(1/sqrt(2))*Vsa+(-1/sqrt(2))*Vsb;
    complex<double> Is_c0=(1/sqrt(6))*(Isa+Isb+Isc);
    complex<double> Is_c1=(sqrt(2)/sqrt(3))*Isc+(-1/sqrt(6))*Isa+(-1/sqrt(6))*Isb;
    complex<double> Is_c2=(1/sqrt(2))*Isa+(-1/sqrt(2))*Isb;
    complex<double> Vr_c0=(1/sqrt(6))*(Vra+Vrb+Vrc);
    complex<double> Vr_c1=(sqrt(2)/sqrt(3))*Vrc+(-1/sqrt(6))*Vra+(-1/sqrt(6))*Vrb;
    complex<double> Vr_c2=(1/sqrt(2))*Vra+(-1/sqrt(2))*Vrb;
    complex<double> Ir_c0=(1/sqrt(6))*(Ira+Irb+Irc);
    complex<double> Ir_c1=(sqrt(2)/sqrt(3))*Irc+(-1/sqrt(6))*Ira+(-1/sqrt(6))*Irb;
    complex<double> Ir_c2=(1/sqrt(2))*Ira+(-1/sqrt(2))*Irb;
    
    //Mc0=(Vs_c0+Zc0*Is_c0)*(exp(gamma0*Lm))/2-(Vr_c0-Zc0*Ir_c0)/2;   
    temp=exp((r0*Lm));       
    temp*=(Vs_c0+Zc0*Is_c0); 
    temp-=(Vr_c0-Zc0*Ir_c0);
    temp/=2;
    Mc0=temp;
        
    //Mc1=(Vs_c1+Zcd*Is_c1)*(exp(gammad*Lm))/2-(Vr_c1-Zcd*Ir_c1)/2;
    temp=exp(r1*Lm);
    temp*=(Vs_c1+Zc1*Is_c1); 
    temp-=(Vr_c1-Zc1*Ir_c1);
    temp/=2;
    Mc1=temp;

    //Mc2=(Vs_c2+Zcd*Is_c2)*(exp(gammad*Lm))/2-(Vr_c2-Zcd*Ir_c2)/2;
    temp=exp(r2*Lm);
    temp*=(Vs_c2+Zc2*Is_c2); 
    temp-=(Vr_c2-Zc2*Ir_c2);
    temp/=2;
    Mc2=temp;    
 
    if(abs(Ma0)>Threshold && abs(Mb0)>Threshold && abs(Mc0)>Threshold){
         if(abs(Ma1)<Threshold){return("AG");}            
         else if(abs(Mb1)<Threshold){return("BG");} 
         else if(abs(Mc1)<Threshold){return("CG");} 
         else {return("2PG");}           
    }
    else if(abs(Ma2)<Threshold){return("BCS");} 
    else if(abs(Mb2)<Threshold){return("CAS");} 
    else if(abs(Mc2)<Threshold){return("ABS");} 
    else{return("3PG");}
}














int main() {
      complex<double> Vsa[1000],Vsb[1000],Vsc[1000],Isa[1000],Isb[1000],Isc[1000],Vra[1000],Vrb[1000],Vrc[1000],Ira[1000],Irb[1000],Irc[1000];
      complex<double> Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_;
      double R0,R1,L0,L1,C0,C1;
      int i=0,t=0;
      double sum_D=0,ave_D;
      double M,D;
      //  complex<double> Vs[3],complex<double> Is[3],complex<double> Vr[3],complex<double> Ir[3],double Rin[2],double Lin[2],double Cin[2]
      double Lm=200;
      (Vsa[0]).real()=-20980.6538363469;//
      (Vsa[0]).imag()=-162259.169314497;//
      (Vsb[0]).real()=-236901.204229647;//
      (Vsb[0]).imag()=174131.421283562;//
      (Vsc[0]).real()=254346.542770219;//
      (Vsc[0]).imag()=168579.687420837;//
      (Isa[0]).real()=-15657.2683014278;//
      (Isa[0]).imag()=-6029.31600886684;//
      (Isb[0]).real()=-1421.86433569417;//
      (Isb[0]).imag()=136.515977262915;//
      (Isc[0]).real()=384.350730430943;//
      (Isc[0]).imag()=968.386891206540;//
      (Vra[0]).real()=-48010.9889656387;//
      (Vra[0]).imag()=-251670.145121785;//
      (Vrb[0]).real()=-220506.417270621;//
      (Vrb[0]).imag()=172446.973311729;//
      (Vrc[0]).real()=263863.633315267;//
      (Vrc[0]).imag()=103133.848753520;//
      (Ira[0]).real()=-4323.58186579215;//
      (Ira[0]).imag()=-726.924754326065;//
      (Irb[0]).real()=805.347409646477;//
      (Irb[0]).imag()=-738.309340390671;//
      (Irc[0]).real()=-888.640107359199;//
      (Irc[0]).imag()=-160.470574633819;//
      R0=0.3479;
      R1=0.0321;
      L0=0.00137;
      L1=0.000473;
      C0=3.80000000000000*0.00000001;
      C1=3.80000000000000*0.00000001;
      
      /*
      string line;
      ifstream S ("S.txt");
      ifstream vsa_r ("Vsa_real.txt");
      ifstream vsa_i ("Vsa_imag.txt");
      ifstream vsb_r ("Vsb_real.txt");
      ifstream vsb_i ("Vsb_imag.txt");
      ifstream vsc_r ("Vsc_real.txt");
      ifstream vsc_i ("Vsc_imag.txt");
   
      ifstream vra_r ("Vra_real.txt");
      ifstream vra_i ("Vra_imag.txt");
      ifstream vrb_r ("Vrb_real.txt");
      ifstream vrb_i ("Vrb_imag.txt");
      ifstream vrc_r ("Vrc_real.txt");
      ifstream vrc_i ("Vrc_imag.txt");
      
      ifstream isa_r ("Isa_real.txt");
      ifstream isa_i ("Isa_imag.txt");
      ifstream isb_r ("Isb_real.txt");
      ifstream isb_i ("Isb_imag.txt");
      ifstream isc_r ("Isc_real.txt");
      ifstream isc_i ("Isc_imag.txt");
   
      ifstream ira_r ("Ira_real.txt");
      ifstream ira_i ("Ira_imag.txt");
      ifstream irb_r ("Irb_real.txt");
      ifstream irb_i ("Irb_imag.txt");
      ifstream irc_r ("Irc_real.txt");
      ifstream irc_i ("Irc_imag.txt");
      if (S.is_open())
      {
         
         while ( S.good() )
         {
              getline (S,line);
              (Vsa[i]).real()=atof(line.c_str());
              (Vsa[i]).imag()=atof(line.c_str());
              //printf("%f <<<-------------- \n",(Vsa[i]).real());
              i++;
         }
   
         S.close();
      }
      */
      
      //for(i=0;i<=385;i++){printf("%f <<<-------------- \n",(Vsa[i]).imag());}
      
      Vsa_=Vsa[i]; 
      Vsb_=Vsb[i];  
      Vsc_=Vsc[i];
      Isa_=Isa[i]; 
      Isb_=Isb[i]; 
      Isc_=Isc[i]; 
      Vra_=Vra[i];  
      Vrb_=Vrb[i];
      Vrc_=Vrc[i];
      Ira_=Ira[i];    
      Irb_=Irb[i]; 
      Irc_=Irc[i];                                                                                                      
      
      
      //complex<double> a(1,2+50000*3.1415926535);
      //complex<double> b;
      //b=exp(a);
      //printf("%f <<< b.real\n",b.real());
      //printf("%f <<< b.imag\n",b.imag());
      //printf("%f <<< M\n",getM(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm));
      printf("%f <<< D\n",getD(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm));
      //printf("%f <<< D*Lm\n",Lm*getD(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm));
      //printf("%s <<< Type\n",getType(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm));
      
      /*
      for(;;){
            if(getM(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm)>300){
                        sum_D+=getD(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm); 
                        t++;  
                        ave_D=sum_D/(double)t;                                                                                                                               
                        i=0;
                        for(i=0;i<=2;i++){
                                          Vsa_=Vsa[t];   
                                          Vsb_=Vsb[t];
                                          Vsc_=Vsc[t];
                                          Isa_=Isa[t]; 
                                          Isb_=Isb[t]; 
                                          Isc_=Isc[t]; 
                                          Vra_=Vra[t]; 
                                          Vrb_=Vrb[t]; 
                                          Vrc_=Vrc[t];
                                          Ira_=Ira[t];  
                                          Irb_=Irb[t]; 
                                          Irc_=Irc[t];                                                                                                          
                        }
                        double d=( getD(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm)-ave_D )*Lm;
                        if(  d>=100 || d<=-100  ){
                             break;          
                        } 
                        else if(t>=10){
                              break;    
                        }
                        else {
                             continue;
                        }                                             
            }
            else break;
      }*/
      
  //printf("%f <<< M\n",getM(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm));
  //printf("%f <<< D\n",getD(Vsa_,Vsb_,Vsc_,Isa_,Isb_,Isc_,Vra_,Vrb_,Vrc_,Ira_,Irb_,Irc_,R0,R1,L0,L1,C0,C1,Lm));



  
  
  
   system("pause");
   
}




