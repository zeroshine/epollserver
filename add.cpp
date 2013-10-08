#include "add.h"
#include <math.h>
using namespace std;
const double Alg:: pi=3.1415926535;
const double Alg:: Lm=26.744;
const double Alg:: T=100;
double Alg:: getD() {
    if(abs(M0)<=abs(M1) && abs(M0)<=abs(M2)){
	return abs(D0);
    }
    else if(abs(M1)<=abs(M2)){
	return abs(D1);
    }
    else return abs(D2);
}

string Alg:: getType(){
    if(abs(Ma0)<T&&abs(Ma1)<T&&abs(Ma2)<T&&abs(Mb0)<T&&abs(Mb1)<T&&abs(Mb2)<T&&abs(Mc0)<T&&abs(Mc1)<T&& abs(Mc2)<T){
	return "safe";
    }
    if(abs(M0)>T&&abs(Mb0)>T&&abs(Mc0)>T){
	if(abs(M2)<T){
	    return("ag");
	} //ag           
	else if(abs(Mb2)<T){
	    return("bg");
	} //bg
	else if(abs(Mc2)<T){
	    return("cg");
	} //cg
	else {
	    return("2pg");
	}  //2pg         
    }else if(abs(M1)<T){
	return("bcs");
    } //bcs
    else if(abs(Mb1)<T){
	return("cas");
    } //cas
    else if(abs(Mc1)<T){
	return("abs");
    } //abs
    else {
	return("3p");
    }//3p
}
void Alg::parse(string s,bool isR){
    istringstream ss(s);
    string word;
    int i=1;
    if(!isR){
	while(ss>>word){
	    cout<<word<<endl;
	    if(i==1){
		Vsam=atof(word.c_str());        
	    }else if(i==2){
		Vsap=atof(word.c_str()); 
	    }else if(i==3){
		Vsbm=atof(word.c_str()); 
	    }else if(i==4){
		Vsbp=atof(word.c_str());
	    }else if(i==5){
		Vscm=atof(word.c_str()); 
	    }else if(i==6){
		Vscp=atof(word.c_str());
	    }else if(i==7){
		Isam=atof(word.c_str());
	    }else if(i==8){
		Isap=atof(word.c_str());
	    }else if(i==9){
		Isbm=atof(word.c_str());
	    }else if(i==10){
		Isbp=atof(word.c_str());
	    }else if(i==11){
		Iscm=atof(word.c_str());
	    }else if(i==12){
		Iscp=atof(word.c_str());
	    }else if(i==13){
		location_index=atof(word.c_str());
	    }else if(i==14){
		timestamp=atol(word.c_str());
	    }else if(i==15){
		_time=word;
	    }
	    i++;
	}
    }else{
	while(ss>>word){
	    cout<<word<<endl;
	    if(i==1){
		Vram=atof(word.c_str());
	    }else if(i==2){
		Vrap=atof(word.c_str()); 
	    }else if(i==3){
		Vrbm=atof(word.c_str()); 
	    }else if(i==4){
		Vrbp=atof(word.c_str());
	    }else if(i==5){
		Vrcm=atof(word.c_str()); 
	    }else if(i==6){
		Vrcp=atof(word.c_str());
	    }else if(i==7){
		Iram=atof(word.c_str());
	    }else if(i==8){
		Irap=atof(word.c_str());
	    }else if(i==9){
		Irbm=atof(word.c_str());
	    }else if(i==10){
		Irbp=atof(word.c_str());
	    }else if(i==11){
		Ircm=atof(word.c_str());
	    }else if(i==12){
		Ircp=atof(word.c_str());
	    }else if(i==13){
		location_index=atof(word.c_str());
	    }else if(i==14){
		timestamp=atol(word.c_str());
	    }else if(i==15){
		_time=word;
	    }
	    i++;
	}
    }
}

void Alg:: compute(){
    complex<double> one(1,0);
    complex<double> a(cos(2*pi/3),sin(2*pi/3));
    complex<double> Zc1=(582.982151*cos(-0.106478),582.982151*sin(-0.106478));
    complex<double> Zc0=(240.430117*cos(-0.031632),240.430117*sin(-0.031632));
    complex<double> r1=(0.001685*cos(1.464318),0.001685*sin(1.464318));
    complex<double> r0=(0.001284*cos(1.539164),0.001284*sin(1.539164 ));

    Vsa=(Vsam*cos(Vsap),Vsam*sin(Vsap));
    Vsb=(Vsbm*cos(Vsbp),Vsbm*sin(Vsbp));
    Vsc=(Vscm*cos(Vscp),Vscm*sin(Vscp));
    Vra=(Vram*cos(Vrap),Vram*sin(Vrap));
    Vrb=(Vrbm*cos(Vrbp),Vrbm*sin(Vrbp));
    Vrc=(Vrcm*cos(Vrcp),Vrcm*sin(Vrcp));
    Isa=(Isam*cos(Isap),Isam*sin(Isap));
    Isb=(Isbm*cos(Isbp),Isbm*sin(Isbp));
    Isc=(Iscm*cos(Iscp),Iscm*sin(Iscp));
    Ira=(Iram*cos(Irap),Iram*sin(Irap));
    Irb=(Irbm*cos(Irbp),Irbm*sin(Irbp));
    Irc=(Ircm*cos(Ircp),Ircm*sin(Ircp));

    Vs0=(Vsa+Vsb+Vsc)/3.0;
    Vs1=(Vsa+Vsb*a+Vsc*a*a)/3.0;
    Vs2=(Vsa+Vsb*a*a+Vsc*a)/3.0;
    Vr0=(Vra+Vrb+Vrc)/3.0;
    Vr1=(Vra+Vrb*a+Vrc*a*a)/3.0;
    Vr2=(Vra+Vrb*a*a+Vrc*a)/3.0;
    Is0=(Isa+Isb+Isc)/3.0;
    Is1=(Isa+Isb*a+Isc*a*a)/3.0;
    Is2=(Isa+Isb*a*a+Isc*a)/3.0;
    Ir0=(Ira+Irb+Irc)/3.0;
    Ir1=(Ira+Irb*a+Irc*a*a)/3.0;
    Ir2=(Ira+Irb*a*a+Irc*a)/3.0;
    //x1=(Vs1*Is1+Vr1*Ir1)/(Vs1*Ir1+Is1*Vr1);
    //x0=(Vs0*Is0+Vr0*Ir0)/(Vs0*Ir0+Is0*Vr0);
    //r1=log(x1+sqrt(x1*x1-one));
    //r0=log(x0+sqrt(x0*x0-one));
    //Zc1=sqrt(((Vs1*Vs1-Vr1*Vr1)/(Is1*Is1-Ir1*Ir1)));
    //Zc0=sqrt(((Vs0*Vs0-Vr0*Vr0)/(Is0*Is0-Ir0*Ir0)));
    //z1=r1*Zc1;
    //y1=r1/Zc1;
    //z0=r0*Zc0;
    //y0=r0/Zc0;
    M0=-((Vr0+Zc0*Ir0)-(Vs0-Zc0*Is0)*exp(r0))/2.0;
    N0=-((Vs0+Zc0*Is0)/exp(r0*Lm)-(Vr0-Zc0*Ir0))/2.0;
    if(abs(N0)!=0){    
	D0=log(M0/N0)/(2.0*r0*Lm);
    }
    // cout<<"D0"<<z1.imag()<<endl;

    N1=-((Vs1+Zc1*Is1)/exp(r1*Lm)-(Vr1-Zc1*Ir1))/2.0;    
    M1=-((Vr1+Zc1*Ir1)-(Vs1-Zc1*Is1)*exp(r1*Lm))/2.0;
    if(abs(N1)!=0){    
	D1=log(M1/N1)/(2.0*r1*Lm);
    }

    M2=-((Vr2+Zc1*Ir2)-(Vs2-Zc1*Is2)*exp(r1*Lm))/2.0;
    N2=-((Vs2+Zc1*Is2)/exp(r1*Lm)-(Vr2-Zc1*Ir2))/2.0;
    if(abs(N2)!=0){    
	D2=log(M2/N2)/(2.0*r1*Lm);
    }


    Vsa0=0.408248290463863*(Vsa+Vsb+Vsc);
    Vsa1=(0.816496580927726*Vsa-0.408248290463863*Vsb-0.408248290463863*Vsc);
    Vsa2=(0.707106781186548*Vsb-0.707106781186548*Vsc);
    Isa0=0.408248290463863*(Isa+Isb+Isc);
    Isa1=(0.816496580927726*Isa-0.408248290463863*Isb-0.408248290463863*Isc);
    Isa2=(0.707106781186548*Isb-0.707106781186548*Isc);
    Vra0=0.408248290463863*(Vra+Vrb+Vrc);
    Vra1=(0.816496580927726*Vra-0.408248290463863*Vrb-0.408248290463863*Vrc);
    Vra2=(0.707106781186548*Vrb-0.707106781186548*Vrc);
    Ira0=0.408248290463863*(Ira+Irb+Irc);
    Ira1=(0.816496580927726*Ira-0.408248290463863*Irb-0.408248290463863*Irc);
    Ira2=(0.707106781186548*Irb-0.707106781186548*Irc);

    //cout<<"A-basis"<<endl;
    Vsb0=0.408248290463863*(Vsa+Vsb+Vsc);
    Vsb1=(0.816496580927726*Vsb-0.408248290463863*Vsc-0.408248290463863*Vsa);
    Vsb2=(0.707106781186548*Vsc-0.707106781186548*Vsa);
    Isb0=0.408248290463863*(Isa+Isb+Isc);
    Isb1=(0.816496580927726*Isb-0.408248290463863*Isc-0.408248290463863*Isa);
    Isb2=(0.707106781186548*Isc-0.707106781186548*Isa);
    Vrb0=0.408248290463863*(Vra+Vrb+Vrc);
    Vrb1=(0.816496580927726*Vrb-0.408248290463863*Vrc-0.408248290463863*Vra);
    Vrb2=(0.707106781186548*Vrc-0.707106781186548*Vra);
    Irb0=0.408248290463863*(Ira+Irb+Irc);
    Irb1=(0.816496580927726*Irb-0.408248290463863*Irc-0.408248290463863*Ira);
    Irb2=(0.707106781186548*Irc-0.707106781186548*Ira);


    //cout<<"B-basis"<<endl;
    Vsc0=0.408248290463863*(Vsa+Vsb+Vsc);
    Vsc1=(0.816496580927726*Vsc-0.408248290463863*Vsa-0.408248290463863*Vsb);
    Vsc2=(0.707106781186548*Vsa-0.707106781186548*Vsb);
    Isc0=0.408248290463863*(Isa+Isb+Isc);
    Isc1=(0.816496580927726*Isc-0.408248290463863*Isa-0.408248290463863*Isb);
    Isc2=(0.707106781186548*Isa-0.707106781186548*Isb);
    Vrc0=0.408248290463863*(Vra+Vrb+Vrc);
    Vrc1=(0.816496580927726*Vrc-0.408248290463863*Vra-0.408248290463863*Vrb);
    Vrc2=(0.707106781186548*Vra-0.707106781186548*Vrb);
    Irc0=0.408248290463863*(Ira+Irb+Irc);
    Irc1=(0.816496580927726*Irc-0.408248290463863*Ira-0.408248290463863*Irb);
    Irc2=(0.707106781186548*Ira-0.707106781186548*Irb);
  //  cout<<"C-basis"<<endl;
    ///////////////////////////////////////////////////////////////////////////////////

    Ma0=-((Vr0+Zc0*Ir0)-(Vs0-Zc0*Is0)*exp(r0*Lm))/2.0;      
    Ma1=-((Vr1+Zc1*Ir1)-(Vs1-Zc1*Is1)*exp(r1*Lm))/2.0;
    Ma2=-((Vr2+Zc1*Ir2)-(Vs2-Zc1*Is2)*exp(r1*Lm))/2.0;
    Mb0=-((Vrb0+Zc0*Irb0)-(Vsb0-Zc0*Isb0)*exp(r0*Lm))/2.0;
    Mb1=-((Vrb1+Zc1*Irb1)-(Vsb1-Zc1*Isb1)*exp(r0*Lm))/2.0;
    Mb2=-((Vrb2+Zc1*Irb2)-(Vsb2-Zc1*Isb2)*exp(r0*Lm))/2.0;
    Mc0=-((Vrc0+Zc0*Irc0)-(Vsc0-Zc0*Isc0)*exp(r0*Lm))/2.0;
    Mc1=-((Vrc1+Zc1*Irc1)-(Vsc1-Zc1*Isc1)*exp(r0*Lm))/2.0;
    Mc2=-((Vrc2+Zc1*Irc2)-(Vsc2-Zc1*Isc2)*exp(r0*Lm))/2.0;
//    cout<<"MN done"<<endl;

}

long long int Alg::gettimestamp(){
    return timestamp;
}

int Alg::getlocation_index(){
    return location_index;
}
/*
   void Alg::comMatrix(){
   double T[3][3],iT[3][3],R[3][3],C[3][3],L[3][3];
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

for(int i=0;i<=2;++i){
for(int j=0;j<=2;++j){
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
for (int i =0;i<=2;++i) {
for(int j=0;j<=2;++j) {
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
for (int i =0 ; i<=2; ++i){
    for(int j=0;j<=2;++j){
	for(int k=0;k<=2;++k){           
	    Y_X_T[i][j]+=Y[i][k]*T[k][j];

	}            
    }      
}

for (int i=0;i<=2;++i){
    for(int j=0;j<=2;++j){
	for(int k=0;k<=2;++k){           
	    Y012[i][j]+=iT[i][k]*Y_X_T[k][j];
	}            
    }       
}

//r=sqrt(Z012*Y012);r0=r(1,1),r1=r(2,2),r2=r(3,3) % r1=r2
complex<double> Z012_X_Y012[3][3];
for (int i=0;i<=2;++i){
    for(int j=0;j<=2;++j){
	for(int k=0;k<=2;++k){           
	    Z012_X_Y012[i][j]+=Z012[i][k]*Y012[k][j];

	}            
    }      
}

for (int i=0;i<=2;++i){
    for(int j=0;j<=2;++j){         
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

for (int i=0;i<=2;++i){
    for(int j=0;j<=2;++j){
	for(int k=0;k<=2;k++){           
	    Z012_div_Y012[i][j]+=Z012[i][k]*invY012[k][j];                           
	}            
    }       
}

for (int i=0 ;i<=2;++i){
    for(int j=0;j<=2;++j){                        
	Zc[i][j]=sqrt(Z012_div_Y012[i][j]);                                                   
    }       
}

Zc0=Zc[0][0]; 
Zc1=Zc[1][1]; 
Zc2=Zc[2][2];

}
*/

//string Alg:: compute(){
//    converseCpx();
//    char s[1024];
//    sprintf(s,"%-012.4f %-012.4f %-012.4f %-012.4f %-012.4f %-012.4f %-012.4f %-012.4f",abs(z1)*cos(arg(z1)),abs(z1)*sin(arg(z1)),abs(y1)*cos(arg(y1)),abs(y1)*sin(arg(y1)),abs(z0)*cos(arg(z0)),abs(z0)*sin(arg(z0)),abs(y0)*cos(arg(y0)),abs(y0)*sin(arg(y0)));
//    return s;         

//}


