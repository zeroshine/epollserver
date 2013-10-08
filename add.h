#ifndef ADD_H
#define ADD_H
#include <string>
#include <complex>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
class Alg{
    public:
	void compute();
	void parse(string s,bool isR); 
        long long int gettimestamp();
	int getlocation_index();
        string getdata();
        string getType();
	double getD();
    private:
	//void converseCpx();
	//void comMatrix();
	double Vsam,Vsap,Vsbm,Vsbp,Vscm,Vscp;
	double Isam,Isap,Isbm,Isbp,Iscm,Iscp;
	double Vram,Vrap,Vrbm,Vrbp,Vrcm,Vrcp;
	double Iram,Irap,Irbm,Irbp,Ircm,Ircp;
	complex<double> Vsa,Vsb,Vsc,Isa,Isb,Isc,Vra,Vrb,Vrc,Ira,Irb,Irc;
        complex<double> Vs0,Vs1,Vs2,Is0,Is1,Is2,Vr0,Vr1,Vr2,Ir0,Ir1,Ir2;
        complex<double> Vsa0,Vsa1,Vsa2,Vsb0,Vsb1,Vsb2,Vsc0,Vsc1,Vsc2;
	complex<double> Isa0,Isa1,Isa2,Isb0,Isb1,Isb2,Isc0,Isc1,Isc2;
        complex<double> Vra0,Vra1,Vra2,Vrb0,Vrb1,Vrb2,Vrc0,Vrc1,Vrc2;
	complex<double> Ira0,Ira1,Ira2,Irb0,Irb1,Irb2,Irc0,Irc1,Irc2;
        complex<double> M0,N0,D0,M1,N1,D1,M2,N2,D2,temp;
	complex<double> Ma0,Mb0,Mc0,Ma1,Mb1,Mc1,Ma2,Mb2,Mc2;
	string _time;
	static const double pi;
	static const double Lm;
	static const double T;
        long long int timestamp;
	int location_index;

};

#endif
