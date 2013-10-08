#include <iostream>
#include <fstream>
#include "add.h"
using namespace std;
int main(){
    ifstream R_file("R.txt");
    ifstream S_file("S.txt");
    string line_R,line_S;
    while(getline(R_file,line_R)&&getline(S_file,line_S)){
        Alg* com = new Alg;
        com->parse(line_R,true);
        com->parse(line_S,false);
        com->compute();
        cout<<com->getType()<<endl;	
	delete com;
    }
    return 0;
}
