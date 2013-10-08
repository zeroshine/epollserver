#include <iostream>
#include<fstream>
#include <complex>
#include <stdio.h>
#include <math.h>
#include <string.h>
using namespace std;


double getZc1_mag(double Vsa_mag,double Vsa_ang,double Vsb_mag,double Vsb_ang,double Vsc_mag,double Vsc_ang,double Isa_mag,double Isa_ang,double Isb_mag,double Isb_ang,double Isc_mag,double Isc_ang,double Vra_mag,double Vra_ang,double Vrb_mag,double Vrb_ang,double Vrc_mag,double Vrc_ang,double Ira_mag,double Ira_ang,double Irb_mag,double Irb_ang,double Irc_mag,double Irc_ang){

    complex<double> Vsa,Vsb,Vsc,Isa,Isb,Isc,Vra,Vrb,Vrc,Ira,Irb,Irc,Zc0,Zc1,r0,r1;
    double pi=3.14159265;
    double w=120*pi;
    complex<double> ang;
    ang.real()=0;
    ang.imag()=2*pi/3.0;
    complex<double> a=exp(ang);
    complex<double> Z0,Z1,Z2,Y0,Y1,Y2;
    
    complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2;
    
    Vsa.real()=Vsa_mag*cos(Vsa_ang);  Vsa.imag()=Vsa_mag*sin(Vsa_ang);
    Vsb.real()=Vsb_mag*cos(Vsb_ang);  Vsb.imag()=Vsb_mag*sin(Vsb_ang);
    Vsc.real()=Vsc_mag*cos(Vsc_ang);  Vsc.imag()=Vsc_mag*sin(Vsc_ang);
    Vra.real()=Vra_mag*cos(Vra_ang);  Vra.imag()=Vra_mag*sin(Vra_ang);
    Vrb.real()=Vrb_mag*cos(Vrb_ang);  Vrb.imag()=Vrb_mag*sin(Vrb_ang);
    Vrc.real()=Vrc_mag*cos(Vrc_ang);  Vrc.imag()=Vrc_mag*sin(Vrc_ang);
    Isa.real()=Isa_mag*cos(Isa_ang);  Isa.imag()=Isa_mag*sin(Isa_ang);
    Isb.real()=Isb_mag*cos(Isb_ang);  Isb.imag()=Isb_mag*sin(Isb_ang);
    Isc.real()=Isc_mag*cos(Isc_ang);  Isc.imag()=Isc_mag*sin(Isc_ang);
    Ira.real()=Ira_mag*cos(Ira_ang);  Ira.imag()=Ira_mag*sin(Ira_ang);
    Irb.real()=Irb_mag*cos(Irb_ang);  Irb.imag()=Irb_mag*sin(Irb_ang);
    Irc.real()=Irc_mag*cos(Irc_ang);  Irc.imag()=Irc_mag*sin(Irc_ang); 
    
    complex<double> Vs0=0.408248290463863*(Vsa+Vsb+Vsc);
    complex<double> Vs1=(0.816496580927726*Vsa-0.408248290463863*Vsb-0.408248290463863*Vsc);
    complex<double> Vs2=(0.707106781186548*Vsb-0.707106781186548*Vsc);
    complex<double> Is0=0.408248290463863*(Isa+Isb+Isc);
    complex<double> Is1=(0.816496580927726*Isa-0.408248290463863*Isb-0.408248290463863*Isc);
    complex<double> Is2=(0.707106781186548*Isb-0.707106781186548*Isc);
    complex<double> Vr0=0.408248290463863*(Vra+Vrb+Vrc);
    complex<double> Vr1=(0.816496580927726*Vra-0.408248290463863*Vrb-0.408248290463863*Vrc);
    complex<double> Vr2=(0.707106781186548*Vrb-0.707106781186548*Vrc);
    complex<double> Ir0=0.408248290463863*(Ira+Irb+Irc);
    complex<double> Ir1=(0.816496580927726*Ira-0.408248290463863*Irb-0.408248290463863*Irc);
    complex<double> Ir2=(0.707106781186548*Irb-0.707106781186548*Irc);
    
    Zc1=sqrt(((Vs1*Vs1)-(Vr1*Vr1))/((Is1*Is1)-(Ir1*Ir1))); 
    return abs(Zc1);
}







double getZc1_imag(complex<double> Vsa,complex<double> Vsb,complex<double> Vsc,complex<double> Isa,complex<double> Isb,complex<double> Isc,complex<double> Vra,complex<double> Vrb,complex<double> Vrc,complex<double> Ira,complex<double> Irb,complex<double> Irc){

    double pi=3.14159265;
    double w=120*pi;
    complex<double> ang;
    ang.real()=0;
    ang.imag()=2*pi/3.0;
    complex<double> a=exp(ang);
    complex<double> Z0,Z1,Z2,Y0,Y1,Y2;
    complex<double> r0,r1,r2,Zc0,Zc1,Zc2;
    complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2;

    complex<double> Vs0=0.408248290463863*(Vsa+Vsb+Vsc);
    complex<double> Vs1=(0.816496580927726*Vsa-0.408248290463863*Vsb-0.408248290463863*Vsc);
    complex<double> Vs2=(0.707106781186548*Vsb-0.707106781186548*Vsc);
    complex<double> Is0=0.408248290463863*(Isa+Isb+Isc);
    complex<double> Is1=(0.816496580927726*Isa-0.408248290463863*Isb-0.408248290463863*Isc);
    complex<double> Is2=(0.707106781186548*Isb-0.707106781186548*Isc);
    complex<double> Vr0=0.408248290463863*(Vra+Vrb+Vrc);
    complex<double> Vr1=(0.816496580927726*Vra-0.408248290463863*Vrb-0.408248290463863*Vrc);
    complex<double> Vr2=(0.707106781186548*Vrb-0.707106781186548*Vrc);
    complex<double> Ir0=0.408248290463863*(Ira+Irb+Irc);
    complex<double> Ir1=(0.816496580927726*Ira-0.408248290463863*Irb-0.408248290463863*Irc);
    complex<double> Ir2=(0.707106781186548*Irb-0.707106781186548*Irc);
    
    Zc1=sqrt(((Vs1*Vs1)-(Vr1*Vr1))/((Is1*Is1)-(Ir1*Ir1))); 
    return imag(Zc1);
}








double getD(double Vsa_mag,double Vsa_ang,double Vsb_mag,double Vsb_ang,double Vsc_mag,double Vsc_ang,double Isa_mag,double Isa_ang,double Isb_mag,double Isb_ang,double Isc_mag,double Isc_ang,double Vra_mag,double Vra_ang,double Vrb_mag,double Vrb_ang,double Vrc_mag,double Vrc_ang,double Ira_mag,double Ira_ang,double Irb_mag,double Irb_ang,double Irc_mag,double Irc_ang,double Zc0_mag,double Zc0_ang,double Zc1_mag,double Zc1_ang,double r0_mag,double r0_ang,double r1_mag,double r1_ang,double Lm) {
    
    complex<double> Vsa,Vsb,Vsc,Isa,Isb,Isc,Vra,Vrb,Vrc,Ira,Irb,Irc,Zc0,Zc1,r0,r1;
    double pi=3.141592653589793;
    double w=120*pi;
    complex<double> ang;
    ang.real()=0;
    ang.imag()=2*pi/3.0;
    complex<double> a=exp(ang);
    complex<double> Zc2=Zc1;
    complex<double> r2=r1;
    complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2;
    
    Vsa.real()=Vsa_mag*cos(Vsa_ang);  Vsa.imag()=Vsa_mag*sin(Vsa_ang);
    Vsb.real()=Vsb_mag*cos(Vsb_ang);  Vsb.imag()=Vsb_mag*sin(Vsb_ang);
    Vsc.real()=Vsc_mag*cos(Vsc_ang);  Vsc.imag()=Vsc_mag*sin(Vsc_ang);
    Vra.real()=Vra_mag*cos(Vra_ang);  Vra.imag()=Vra_mag*sin(Vra_ang);
    Vrb.real()=Vrb_mag*cos(Vrb_ang);  Vrb.imag()=Vrb_mag*sin(Vrb_ang);
    Vrc.real()=Vrc_mag*cos(Vrc_ang);  Vrc.imag()=Vrc_mag*sin(Vrc_ang);
    Isa.real()=Isa_mag*cos(Isa_ang);  Isa.imag()=Isa_mag*sin(Isa_ang);
    Isb.real()=Isb_mag*cos(Isb_ang);  Isb.imag()=Isb_mag*sin(Isb_ang);
    Isc.real()=Isc_mag*cos(Isc_ang);  Isc.imag()=Isc_mag*sin(Isc_ang);
    Ira.real()=Ira_mag*cos(Ira_ang);  Ira.imag()=Ira_mag*sin(Ira_ang);
    Irb.real()=Irb_mag*cos(Irb_ang);  Irb.imag()=Irb_mag*sin(Irb_ang);
    Irc.real()=Irc_mag*cos(Irc_ang);  Irc.imag()=Irc_mag*sin(Irc_ang); 
       
    Zc0.real()=Zc0_mag*cos(Zc0_ang);  Zc0.imag()=Zc0_mag*sin(Zc0_ang);
    Zc1.real()=Zc1_mag*cos(Zc1_ang);  Zc1.imag()=Zc1_mag*sin(Zc1_ang);
    
    r0.real()=r0_mag*cos(r0_ang);  r0.imag()=r0_mag*sin(r0_ang);
    r1.real()=r1_mag*cos(r1_ang);  r1.imag()=r1_mag*sin(r1_ang);
    
      
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
    
    /*
    complex<double> Vs0=0.408248290463863*(Vsa+Vsb+Vsc);
    complex<double> Vs1=(0.816496580927726*Vsa-0.408248290463863*Vsb-0.408248290463863*Vsc);
    complex<double> Vs2=(0.707106781186548*Vsb-0.707106781186548*Vsc);
    complex<double> Is0=0.408248290463863*(Isa+Isb+Isc);
    complex<double> Is1=(0.816496580927726*Isa-0.408248290463863*Isb-0.408248290463863*Isc);
    complex<double> Is2=(0.707106781186548*Isb-0.707106781186548*Isc);
    complex<double> Vr0=0.408248290463863*(Vra+Vrb+Vrc);
    complex<double> Vr1=(0.816496580927726*Vra-0.408248290463863*Vrb-0.408248290463863*Vrc);
    complex<double> Vr2=(0.707106781186548*Vrb-0.707106781186548*Vrc);
    complex<double> Ir0=0.408248290463863*(Ira+Irb+Irc);
    complex<double> Ir1=(0.816496580927726*Ira-0.408248290463863*Irb-0.408248290463863*Irc);
    complex<double> Ir2=(0.707106781186548*Irb-0.707106781186548*Irc);
    */
    ///////////////////////////////////////////////////////////////////////////////////    
    
    M0=((Vr0+Zc0*Ir0)-(Vs0-Zc0*Is0)*exp(r0*Lm))/2.0;
    N0=((Vs0+Zc0*Is0)/exp(r0*Lm)-(Vr0-Zc0*Ir0))/2.0;
    D0=log(M0/N0)/(2.0*r0*Lm);
    
    N1=((Vs1+Zc1*Is1)/exp(r1*Lm)-(Vr1-Zc1*Ir1))/2.0;    
    M1=((Vr1+Zc1*Ir1)-(Vs1-Zc1*Is1)*exp(r1*Lm))/2.0;
    D1=log(M1/N1)/(2.0*r1*Lm);

    M2=-((Vr2+Zc1*Ir2)-(Vs2-Zc1*Is2)*exp(r1*Lm))/2.0;
    N2=-((Vs2+Zc1*Is2)/exp(r1*Lm)-(Vr2-Zc1*Ir2))/2.0;
    D2=log(M2/N2)/(2.0*r1*Lm);
    
    printf("%f +%f i  <<<<< Vs1\n",real(Vs1),imag(Vs1));
    printf("%f +%f i  <<<<< Vr1\n",real(Vr1),imag(Vr1));   
    printf("%f +%f i  <<<<< Is1\n",real(Is1),imag(Is1)); 
    printf("%f +%f i  <<<<< Ir1\n",real(Ir1),imag(Ir1));
    printf("%f +%f i  <<<<< Zc1\n",real(Zc1),imag(Zc1));
    printf("%f +%f i  <<<<< Zc0\n",real(Zc0),imag(Zc0));
    printf("%f +%f i  <<<<< r1\n",real(r1),imag(r1));
    printf("%f +%f i  <<<<< r0\n",real(r0),imag(r0));
    printf("%f +%f i  <<<<< M1\n",real(M1),imag(M1));
    printf("%f +%f i  <<<<< N1\n",real(N1),imag(N1));
    printf("%f +%f i  <<<<< D0\n",real(D0),imag(D0)); 
    printf("%f +%f i  <<<<< D1\n",real(D1),imag(D1));   
    printf("%f +%f i  <<<<< D2\n",real(D2),imag(D2));
    
    if(real(D0)<1 && real(D0)>0 && abs(imag(D0))<abs(imag(D1)) && abs(imag(D0))<abs(imag(D2))){
        return abs(D0);
    }
    else if(real(D1)<1 && real(D1)>0 && abs(imag(D1))<abs(imag(D0)) && abs(imag(D1))<abs(imag(D2))){
        return abs(D1);
    }
    else if(real(D2)<1 && real(D2)>0 && abs(imag(D2))<abs(imag(D0)) && abs(imag(D2))<abs(imag(D1))){
        return abs(D2);
    }
    else return -1;
}











int getType(double Vsa_mag,double Vsa_ang,double Vsb_mag,double Vsb_ang,double Vsc_mag,double Vsc_ang,double Isa_mag,double Isa_ang,double Isb_mag,double Isb_ang,double Isc_mag,double Isc_ang,double Vra_mag,double Vra_ang,double Vrb_mag,double Vrb_ang,double Vrc_mag,double Vrc_ang,double Ira_mag,double Ira_ang,double Irb_mag,double Irb_ang,double Irc_mag,double Irc_ang,double Zc0_mag,double Zc0_ang,double Zc1_mag,double Zc1_ang,double r0_mag,double r0_ang,double r1_mag,double r1_ang,double Lm) {
    
    complex<double> Vsa,Vsb,Vsc,Isa,Isb,Isc,Vra,Vrb,Vrc,Ira,Irb,Irc,Zc0,Zc1,r0,r1;
    double pi=3.14159265;
    double w=120*pi;
    complex<double> ang;
    ang.real()=0;
    ang.imag()=2*pi/3.0;
    complex<double> a=exp(ang);
    complex<double> Zc2=Zc1;
    complex<double> r2=r1;
    complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2;
    complex<double> Mb0,Mb1,Mb2,Mc0,Mc1,Mc2;
    double Threshold=100;
    
    Vsa.real()=Vsa_mag*cos(Vsa_ang);  Vsa.imag()=Vsa_mag*sin(Vsa_ang);
    Vsb.real()=Vsb_mag*cos(Vsb_ang);  Vsb.imag()=Vsb_mag*sin(Vsb_ang);
    Vsc.real()=Vsc_mag*cos(Vsc_ang);  Vsc.imag()=Vsc_mag*sin(Vsc_ang);
    Vra.real()=Vra_mag*cos(Vra_ang);  Vra.imag()=Vra_mag*sin(Vra_ang);
    Vrb.real()=Vrb_mag*cos(Vrb_ang);  Vrb.imag()=Vrb_mag*sin(Vrb_ang);
    Vrc.real()=Vrc_mag*cos(Vrc_ang);  Vrc.imag()=Vrc_mag*sin(Vrc_ang);
    Isa.real()=Isa_mag*cos(Isa_ang);  Isa.imag()=Isa_mag*sin(Isa_ang);
    Isb.real()=Isb_mag*cos(Isb_ang);  Isb.imag()=Isb_mag*sin(Isb_ang);
    Isc.real()=Isc_mag*cos(Isc_ang);  Isc.imag()=Isc_mag*sin(Isc_ang);
    Ira.real()=Ira_mag*cos(Ira_ang);  Ira.imag()=Ira_mag*sin(Ira_ang);
    Irb.real()=Irb_mag*cos(Irb_ang);  Irb.imag()=Irb_mag*sin(Irb_ang);
    Irc.real()=Irc_mag*cos(Irc_ang);  Irc.imag()=Irc_mag*sin(Irc_ang); 
       
    Zc0.real()=Zc0_mag*cos(Zc0_ang);  Zc0.imag()=Zc0_mag*sin(Zc0_ang);
    Zc1.real()=Zc1_mag*cos(Zc1_ang);  Zc1.imag()=Zc1_mag*sin(Zc1_ang);
    
    r0.real()=r0_mag*cos(r0_ang);  r0.imag()=r0_mag*sin(r0_ang);
    r1.real()=r1_mag*cos(r1_ang);  r1.imag()=r1_mag*sin(r1_ang);
    
    complex<double> Vs0=0.408248290463863*(Vsa+Vsb+Vsc);
    complex<double> Vs1=(0.816496580927726*Vsa-0.408248290463863*Vsb-0.408248290463863*Vsc);
    complex<double> Vs2=(0.707106781186548*Vsb-0.707106781186548*Vsc);
    complex<double> Is0=0.408248290463863*(Isa+Isb+Isc);
    complex<double> Is1=(0.816496580927726*Isa-0.408248290463863*Isb-0.408248290463863*Isc);
    complex<double> Is2=(0.707106781186548*Isb-0.707106781186548*Isc);
    complex<double> Vr0=0.408248290463863*(Vra+Vrb+Vrc);
    complex<double> Vr1=(0.816496580927726*Vra-0.408248290463863*Vrb-0.408248290463863*Vrc);
    complex<double> Vr2=(0.707106781186548*Vrb-0.707106781186548*Vrc);
    complex<double> Ir0=0.408248290463863*(Ira+Irb+Irc);
    complex<double> Ir1=(0.816496580927726*Ira-0.408248290463863*Irb-0.408248290463863*Irc);
    complex<double> Ir2=(0.707106781186548*Irb-0.707106781186548*Irc);
    
    //B-basis
    complex<double> Vsb0=0.408248290463863*(Vsa+Vsb+Vsc);
    complex<double> Vsb1=(0.816496580927726*Vsb-0.408248290463863*Vsc-0.408248290463863*Vsa);
    complex<double> Vsb2=(0.707106781186548*Vsc-0.707106781186548*Vsa);
    complex<double> Isb0=0.408248290463863*(Isa+Isb+Isc);
    complex<double> Isb1=(0.816496580927726*Isb-0.408248290463863*Isc-0.408248290463863*Isa);
    complex<double> Isb2=(0.707106781186548*Isc-0.707106781186548*Isa);
    complex<double> Vrb0=0.408248290463863*(Vra+Vrb+Vrc);
    complex<double> Vrb1=(0.816496580927726*Vrb-0.408248290463863*Vrc-0.408248290463863*Vra);
    complex<double> Vrb2=(0.707106781186548*Vrc-0.707106781186548*Vra);
    complex<double> Irb0=0.408248290463863*(Ira+Irb+Irc);
    complex<double> Irb1=(0.816496580927726*Irb-0.408248290463863*Irc-0.408248290463863*Ira);
    complex<double> Irb2=(0.707106781186548*Irc-0.707106781186548*Ira);
    
    
    //C-basis
    complex<double> Vsc0=0.408248290463863*(Vsa+Vsb+Vsc);
    complex<double> Vsc1=(0.816496580927726*Vsc-0.408248290463863*Vsa-0.408248290463863*Vsb);
    complex<double> Vsc2=(0.707106781186548*Vsa-0.707106781186548*Vsb);
    complex<double> Isc0=0.408248290463863*(Isa+Isb+Isc);
    complex<double> Isc1=(0.816496580927726*Isc-0.408248290463863*Isa-0.408248290463863*Isb);
    complex<double> Isc2=(0.707106781186548*Isa-0.707106781186548*Isb);
    complex<double> Vrc0=0.408248290463863*(Vra+Vrb+Vrc);
    complex<double> Vrc1=(0.816496580927726*Vrc-0.408248290463863*Vra-0.408248290463863*Vrb);
    complex<double> Vrc2=(0.707106781186548*Vra-0.707106781186548*Vrb);
    complex<double> Irc0=0.408248290463863*(Ira+Irb+Irc);
    complex<double> Irc1=(0.816496580927726*Irc-0.408248290463863*Ira-0.408248290463863*Irb);
    complex<double> Irc2=(0.707106781186548*Ira-0.707106781186548*Irb);
    ///////////////////////////////////////////////////////////////////////////////////
           
    M0=-((Vr0+Zc0*Ir0)-(Vs0-Zc0*Is0)*exp(r0*Lm))/2.0;      
    M1=-((Vr1+Zc1*Ir1)-(Vs1-Zc1*Is1)*exp(r1*Lm))/2.0;
    M2=-((Vr2+Zc1*Ir2)-(Vs2-Zc1*Is2)*exp(r1*Lm))/2.0;
    Mb0=-((Vrb0+Zc0*Irb0)-(Vsb0-Zc0*Isb0)*exp(r0*Lm))/2.0;
    Mb1=-((Vrb1+Zc1*Irb1)-(Vsb1-Zc1*Isb1)*exp(r0*Lm))/2.0;
    Mb2=-((Vrb2+Zc1*Irb2)-(Vsb2-Zc1*Isb2)*exp(r0*Lm))/2.0;
    Mc0=-((Vrc0+Zc0*Irc0)-(Vsc0-Zc0*Isc0)*exp(r0*Lm))/2.0;
    Mc1=-((Vrc1+Zc1*Irc1)-(Vsc1-Zc1*Isc1)*exp(r0*Lm))/2.0;
    Mc2=-((Vrc2+Zc1*Irc2)-(Vsc2-Zc1*Isc2)*exp(r0*Lm))/2.0;
    
    printf("%f  <<<<< |M0|\n",abs(M0));
    printf("%f  <<<<< |M1|\n",abs(M1));
    printf("%f  <<<<< |M2|\n",abs(M2));
    printf("%f  <<<<< |Mb0|\n",abs(Mb0));
    printf("%f  <<<<< |Mb1|\n",abs(Mb1));
    printf("%f  <<<<< |Mb2|\n",abs(Mb2));
    printf("%f  <<<<< |Mc0|\n",abs(Mc0));
    printf("%f  <<<<< |Mc1|\n",abs(Mc1));
    printf("%f  <<<<< |Mc2|\n",abs(Mc2));
    
    if(abs(M0)<Threshold && abs(M1)<Threshold && abs(M2)<Threshold && abs(Mb0)<Threshold && abs(Mb1)<Threshold && abs(Mb2)<Threshold && abs(Mc0)<Threshold && abs(Mc1)<Threshold && abs(Mc2)<Threshold)
    {return 0;}
    if(abs(M0)>Threshold && abs(Mb0)>Threshold && abs(Mc0)>Threshold){
         if(abs(M2)<Threshold){return(1);} //ag           
         else if(abs(Mb2)<Threshold){return(2);} //bg
         else if(abs(Mc2)<Threshold){return(3);} //cg
         else {return(20);}  //2pg         
    }
    else if(abs(M1)<Threshold){return(23);} //bcs
    else if(abs(Mb1)<Threshold){return(31);} //cas
    else if(abs(Mc1)<Threshold){return(12);} //abs
    else {return(30);}//3p
}







int main() {
    
      double Vsa_mag,Vsa_ang,Vsb_mag,Vsb_ang,Vsc_mag,Vsc_ang,Isa_mag,Isa_ang,Isb_mag,Isb_ang,Isc_mag,Isc_ang,Vra_mag,Vra_ang,Vrb_mag,Vrb_ang,Vrc_mag,Vrc_ang,Ira_mag,Ira_ang,Irb_mag,Irb_ang,Irc_mag,Irc_ang,Zc0_mag,Zc0_ang,Zc1_mag,Zc1_ang,r0_mag,r0_ang,r1_mag,r1_ang;       
      double Lm=26.744;
      Zc0_mag=582.981302;
      Zc0_ang=-0.106477;  
      Zc1_mag=240.427786;
      Zc1_ang=-0.031633;
      r0_mag=0.001685;
      r0_ang=1.464319; 
      r1_mag=0.001284;
      r1_ang=1.539164; 
      
      //abg
      
      Vsa_mag=130.306267;  Vsa_ang=0.461748;
      Vsb_mag=133.614316;  Vsb_ang=-1.561165;
      Vsc_mag=324.175052;  Vsc_ang=2.442678;
      Vra_mag=101.101228;  Vra_ang=0.538367;
      Vrb_mag=110.167216;  Vrb_ang=-1.500249;
      Vrc_mag=329.230219;  Vrc_ang=2.502403;      
      Isa_mag=18.234762;  Isa_ang=-0.694024;
      Isb_mag=18.366729;  Isb_ang=2.9398;
      Isc_mag=1.893461;  Isc_ang=-0.431052;
      Ira_mag=45.282593;  Ira_ang=-0.588117;
      Irb_mag=44.336124;  Irb_ang=3.051381;
      Irc_mag=1.904676;  Irc_ang=2.732502;      
      
      
      //abcg
      /*
      Vsa_mag=115.701479;  Vsa_ang=1.739284;
      Vsb_mag=115.601396;  Vsb_ang=-0.353623;
      Vsc_mag=115.500884;  Vsc_ang=-2.449506;
      Vra_mag=144.123805;  Vra_ang=1.836625;
      Vrb_mag=143.954618;  Vrb_ang=-0.256512;
      Vrc_mag=143.875349;  Vrc_ang=-2.352561;      
      Isa_mag=27.762963;  Isa_ang=0.267669;
      Isb_mag=27.792904;  Isb_ang=-1.828136;
      Isc_mag=27.809657;  Isc_ang=2.362245;
      Ira_mag=34.768055;  Ira_ang=0.354269;
      Irb_mag=34.813105;  Irb_ang=-1.741347;
      Irc_mag=34.823824;  Irc_ang=2.449132;
      */
      
      
      printf("%f <<<<<<D\n",getD(Vsa_mag,Vsa_ang,Vsb_mag,Vsb_ang,Vsc_mag,Vsc_ang,Isa_mag,Isa_ang,Isb_mag,Isb_ang,Isc_mag,Isc_ang,Vra_mag,Vra_ang,Vrb_mag,Vrb_ang,Vrc_mag,Vrc_ang,Ira_mag,Ira_ang,Irb_mag,Irb_ang,Irc_mag,Irc_ang,Zc0_mag,Zc0_ang,Zc1_mag,Zc1_ang,r0_mag,r0_ang,r1_mag,r1_ang,Lm));
      printf("%d <<<<<<D\n",getType(Vsa_mag,Vsa_ang,Vsb_mag,Vsb_ang,Vsc_mag,Vsc_ang,Isa_mag,Isa_ang,Isb_mag,Isb_ang,Isc_mag,Isc_ang,Vra_mag,Vra_ang,Vrb_mag,Vrb_ang,Vrc_mag,Vrc_ang,Ira_mag,Ira_ang,Irb_mag,Irb_ang,Irc_mag,Irc_ang,Zc0_mag,Zc0_ang,Zc1_mag,Zc1_ang,r0_mag,r0_ang,r1_mag,r1_ang,Lm));

  
     system("pause");
     
     
      /*
      Vsa_mag=;  Vsa_ang=1;
      Vsb_mag=;  Vsb_ang=;
      Vsc_mag=;  Vsc_ang=;
      Vra_mag=;  Vra_ang=;
      Vrb_mag=;  Vrb_ang=;
      Vrc_mag=;  Vrc_ang=;      
      Isa_mag=;  Isa_ang=;
      Isb_mag=;  Isb_ang=;
      Isc_mag=;  Isc_ang=;
      Ira_mag=;  Ira_ang=;
      Irb_mag=;  Irb_ang=;
      Irc_mag=;  Irc_ang=;
      */
   
}




