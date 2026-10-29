#include <cstdlib>
#include <iostream>
#include <math.h>

// LEO interception simulation

using namespace std;

double r, s, a, ba, mu;

int i, j, k;

int main(int argc, char *argv[])
{
    double r0[3];
    double br0[3];
    double v0[3];
    double bv0[3];
    double v[9];
    double vc[9];
    double f[9];
    double Jf[9][9];
    double b[9];
    double id[9][9];
    double c[9];
    double prod[9];

    // Initial conditions
    
    r0[0]=1;
    r0[1]=2;
    r0[2]=3;
    
    br0[0]=4;
    br0[1]=5;
    br0[2]=6;
    
    v0[0]=1;
    v0[1]=1;
    v0[2]=1;
    
    bv0[0]=2;
    bv0[1]=2;
    bv0[2]=2;
    
    a=1;
    ba=2;
    mu=1;

// Solution method
    
    v[0]=1;
    v[1]=1;
    v[2]=1;
    v[3]=1;
    v[4]=1;
    v[5]=1;
    v[6]=1;
    v[7]=1;
    v[8]=1;
    
    vc[0]=0;
    vc[1]=0;
    vc[2]=0;
    vc[3]=0;
    vc[4]=0;
    vc[5]=0;
    vc[6]=0;
    vc[7]=0;
    vc[8]=0;
    
    while (sqrt(((vc[0]-v[0])*(vc[0]-v[0]))+((vc[1]-v[1])*(vc[1]-v[1]))+((vc[2]-v[2])*(vc[2]-v[2]))+((vc[3]-v[3])*(vc[3]-v[3]))+((vc[4]-v[4])*(vc[4]-v[4]))+((vc[5]-v[5])*(vc[5]-v[5]))+((vc[6]-v[6])*(vc[6]-v[6]))+((vc[7]-v[7])*(vc[7]-v[7]))+((vc[8]-v[8])*(vc[8]-v[8])))>0.00001){
    
    for(i=0; i<9; i++){
       vc[i]=v[i];
    }
    
    f[0]=(ba*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))) + (ba/sqrt(mu)*((br0[0]*bv0[0])+(br0[1]*bv0[1])+(br0[2]*bv0[2]))*(1-(cos(v[0]/sqrt(ba))))) + (sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sqrt(ba)*sin(v[0]/sqrt(ba))) - (sqrt(mu)*v[2]);
    f[1]=(a*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))) + (a/sqrt(mu)*((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))*(1-(cos(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sqrt(a)*sin(v[1]/sqrt(a))) - (sqrt(mu)*v[2]);
    f[2]=((1-(ba/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*(1-cos(v[0]/sqrt(ba)))))*br0[0]) + ((v[2]-(ba/sqrt(mu)*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))))*bv0[0]) - ((1-(a/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a)))))*r0[0]) - ((v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*(v0[0]+v[3]));
    f[3]=((1-(ba/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*(1-cos(v[0]/sqrt(ba)))))*br0[1]) + ((v[2]-(ba/sqrt(mu)*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))))*bv0[1]) - ((1-(a/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a)))))*r0[1]) - ((v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*(v0[1]+v[4]));
    f[4]=((1-(ba/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*(1-cos(v[0]/sqrt(ba)))))*br0[2]) + ((v[2]-(ba/sqrt(mu)*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))))*bv0[2]) - ((1-(a/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a)))))*r0[2]) - ((v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*(v0[1]+v[5]));
    f[5]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*((v[6]*br0[0])+(v[7]*br0[1])+(v[8]*br0[2]))) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*((v[6]*bv0[0])+(v[7]*bv0[1])+(v[8]*bv0[2]))) + (sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))) + ((a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v0[0]-v[3]))+(v[7]*(bv0[1]-v0[1]-v[4]))+(v[8]*(bv0[2]-v0[2]-v[5]))));
    f[6]=v[3] + ((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))))*((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[0]))) - ( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) ) - ((v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*v[6]);
    f[7]=v[4] + ((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))))*((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[1]))) - ( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) ) - ((v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*v[7]);
    f[8]=v[5] + ((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))))*((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[2]))) - ( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) ) - ((v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*v[8]);
    
    cout<<f[0]<<endl;
    
    Jf[0][0]=(ba+(sqrt(ba/mu)*((br0[0]*bv0[0])+(br0[1]*bv0[1])+(br0[2]*bv0[2]))*sin(v[0]/sqrt(ba)))+((sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))-ba)*cos(v[0]/sqrt(ba))));
    Jf[0][1]=0;
    Jf[0][2]=-sqrt(mu);
    Jf[0][3]=0;
    Jf[0][4]=0;
    Jf[0][5]=0;
    Jf[0][6]=0;
    Jf[0][7]=0;
    Jf[0][8]=0;
    
    Jf[1][0]=0;
    Jf[1][1]=(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a))));
    Jf[1][2]=-sqrt(mu);
    Jf[1][3]=( (2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[0]) );
    Jf[1][4]=( (2*a*a/mu*(v0[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[1]) );
    Jf[1][5]=( (2*a*a/mu*(v0[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[2]) );
    Jf[1][6]=0;
    Jf[1][7]=0;
    Jf[1][8]=0;
    
    Jf[2][0]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[0]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[0]);
    Jf[2][1]=-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[0]+v[3]));
    Jf[2][2]=bv0[0]-v0[0]-v[3];
    Jf[2][3]=-( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) ) - (v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))));
    Jf[2][4]=-( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[2][5]=-( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[2][6]=0;
    Jf[2][7]=0;
    Jf[2][8]=0;
    
    Jf[3][0]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[1]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[1]);
    Jf[3][1]=-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[1]+v[4]));
    Jf[3][2]=bv0[1]-v0[1]-v[4];
    Jf[3][3]=-( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[3][4]=-( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) ) - (v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))));
    Jf[3][5]=-( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[3][6]=0;
    Jf[3][7]=0;
    Jf[3][8]=0;
    
    Jf[4][0]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[2]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[2]);
    Jf[4][1]=-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[2]+v[5]));
    Jf[4][2]=bv0[0]-v0[0]-v[3];
    Jf[4][3]=-( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[4][4]=-( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[4][5]=-( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) ) - (v[2]-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))));
    Jf[4][6]=0;
    Jf[4][7]=0;
    Jf[4][8]=0;
    
    Jf[5][0]=-(1/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*cos(v[0]/sqrt(ba))*((v[6]*br0[0])+(v[7]*br0[1])+(v[8]*br0[2]))) - (sqrt(ba/mu)*sin(v[0]/sqrt(ba))*((v[6]*bv0[0])+(v[7]*bv0[1])+(v[8]*bv0[2])));
    Jf[5][1]=(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*cos(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + (sqrt(a/mu)*sin(v[1]/sqrt(a))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))) + (((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v0[0]-v[3]))+(v[7]*(bv0[1]-v0[1]-v[4]))+(v[8]*(bv0[2]-v0[2]-v[5]))));
    Jf[5][2]=0;
    Jf[5][3]=( (((((1/(2*sqrt(a)*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*sin(v[1]/sqrt(a)))-(v[1]/(2*a*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + ((1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))) + ((1+(((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v0[0]-v[3]))+(v[7]*(bv0[1]-v0[1]-v[4]))+(v[8]*(bv0[2]-v0[2]-v[5])))))*(2*a*a*(v0[0]+v[3])/mu)) + (((a-(a*cos(v[1]/sqrt(a))))-(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*v[6]/sqrt(mu)) );
    Jf[5][4]=( (((((1/(2*sqrt(a)*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*sin(v[1]/sqrt(a)))-(v[1]/(2*a*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + ((1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))) + ((1+(((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v0[0]-v[3]))+(v[7]*(bv0[1]-v0[1]-v[4]))+(v[8]*(bv0[2]-v0[2]-v[5])))))*(2*a*a*(v0[1]+v[4])/mu)) + (((a-(a*cos(v[1]/sqrt(a))))-(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*v[7]/sqrt(mu)) );
    Jf[5][5]=( (((((1/(2*sqrt(a)*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*sin(v[1]/sqrt(a)))-(v[1]/(2*a*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + ((1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))) + ((1+(((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v0[0]-v[3]))+(v[7]*(bv0[1]-v0[1]-v[4]))+(v[8]*(bv0[2]-v0[2]-v[5])))))*(2*a*a*(v0[2]+v[5])/mu)) + (((a-(a*cos(v[1]/sqrt(a))))-(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*v[8]/sqrt(mu)) );
    Jf[5][6]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[0]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[0]) + (sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[0]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[0]+v[3])) + ((a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a))))/sqrt(mu)*(bv0[0]-v0[0]-v[3]));
    Jf[5][7]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[1]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[1]) + (sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[1]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[1]+v[4])) + ((a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a))))/sqrt(mu)*(bv0[1]-v0[1]-v[4]));
    Jf[5][8]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[2]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[2]) + (sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[2]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[2]+v[5])) + ((a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a))))/sqrt(mu)*(bv0[2]-v0[2]-v[5]));
    
    Jf[6][0]=0;
    Jf[6][1]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(((((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a)))) - (1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + ((((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) - (sqrt(a/mu)*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) *((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[0]))) + ((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))))*((2*a*a/mu*(v0[0]+v[3])*( 1 - cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*(v[1]*cos(v[1]/sqrt(a))/(sqrt(a)*sqrt(a)*sqrt(a)))) )) + (sqrt(a/mu)*(sin(v[1]/sqrt(a)))*r0[0]))) - ( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (1/sqrt(mu)*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) ) - (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v[6]);
    Jf[6][2]=-v[6];
    Jf[6][3]=0;
    Jf[6][4]=0;
    Jf[6][5]=0;
    Jf[6][6]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[0]+v[3])))*((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[0]))) - ( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[6][7]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[1]+v[4])))*((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[0]))) - ( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[1]+v[4]))) );
    Jf[6][8]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[2]+v[5])))*((2*a*a/mu*(v0[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[0]))) - ( 2*a*a/mu*(v0[0]+v[3])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[2]+v[5]))) );
    
    Jf[7][0]=0;
    Jf[7][1]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(((((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a)))) - (1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + ((((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) - (sqrt(a/mu)*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) *((2*a*a/mu*(v0[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[1]))) + ((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))))*((2*a*a/mu*(v0[1]+v[4])*( 1 - cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*(v[1]*cos(v[1]/sqrt(a))/(sqrt(a)*sqrt(a)*sqrt(a)))) )) + (sqrt(a/mu)*(sin(v[1]/sqrt(a)))*r0[1]))) - ( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (1/sqrt(mu)*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) ) - (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v[7]);
    Jf[7][2]=-v[7];
    Jf[7][3]=0;
    Jf[7][4]=0;
    Jf[7][5]=0;
    Jf[7][6]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[0]+v[3])))*((2*a*a/mu*(v0[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[1]))) - ( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[7][7]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[1]+v[4])))*((2*a*a/mu*(v0[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[1]))) - ( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[1]+v[4]))) );
    Jf[7][8]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[2]+v[5])))*((2*a*a/mu*(v0[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[1]))) - ( 2*a*a/mu*(v0[1]+v[4])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[2]+v[5]))) );

    Jf[8][0]=0;
    Jf[8][1]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(((((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a)))) - (1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) + ((((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*((((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[1]/sqrt(a)))-((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) - (sqrt(a/mu)*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) *((2*a*a/mu*(v0[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[2]))) + ((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5])))))*((2*a*a/mu*(v0[2]+v[5])*( 1 - cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*(v[1]*cos(v[1]/sqrt(a))/(sqrt(a)*sqrt(a)*sqrt(a)))) )) + (sqrt(a/mu)*(sin(v[1]/sqrt(a)))*r0[2]))) - ( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))*((v[6]*r0[0])+(v[7]*r0[1])+(v[8]*r0[2]))) - (1/sqrt(mu)*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*(v0[0]+v[3]))+(v[7]*(v0[1]+v[4]))+(v[8]*(v0[2]+v[5]))))) ) - (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v[8]);
    Jf[8][2]=-v[8];
    Jf[8][3]=0;
    Jf[8][4]=0;
    Jf[8][5]=0;
    Jf[8][6]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[0]+v[3])))*((2*a*a/mu*(v0[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[2]))) - ( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[0]+v[3]))) );
    Jf[8][7]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[1]+v[4])))*((2*a*a/mu*(v0[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[2]))) - ( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[1]+v[4]))) );
    Jf[8][8]=((1/(a+(sqrt(a/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[1]/sqrt(a)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sin(v[1]/sqrt(a))*r0[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v0[2]+v[5])))*((2*a*a/mu*(v0[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r0[2]))) - ( 2*a*a/mu*(v0[2]+v[5])*(-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r0[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v0[2]+v[5]))) );   
    
   for(i=0; i<9; i++){
      for(j=0; j<9; j++){
         cout<<Jf[i][j]<<"  ";
      }
      cout<<endl;
   }
    
    for(i=0; i<9; i++){
      for(j=0; j<9; j++){
          if(i==j){
             id[i][j]=1;
          }else{
             id[i][j]=0;
          }
      }
    }
    
    for(i=0; i<9; i++){
             
          if(Jf[i][i]==0){
             for(k=0; k<9; k++){
                if(k!=i && Jf[k][i]>0){
                   for(j=0; j<9; j++){
                      b[j]=Jf[i][j];
                         c[j]=id[i][j];
                      Jf[i][j]=Jf[k][j];
                         id[i][j]=id[k][j];
                      Jf[k][j]=b[j];
                         id[k][j]=c[j]; 
                   }
                }
             }
          }       
            
             
       s=Jf[i][i];
       for(j=0; j<9; j++){
          Jf[i][j]=Jf[i][j]/s;
             id[i][j]=id[i][j]/s;
          b[j]=Jf[i][j];
             c[j]=id[i][j];  
       }
       
       for(k=0; k<9; k++){
          if(k!=i){   
             r=Jf[k][i];
             for(j=0; j<9; j++){
                Jf[k][j]=Jf[k][j]-(r*b[j]);
                   id[k][j]=id[k][j]-(r*c[j]);
             }
          }
       }           
    }
    
    for(i=0; i<9; i++){
       prod[i]=0;
    }
    
    for(i=0; i<9; i++){
       for(j=0; j<9; j++){
          prod[i]=prod[i]+id[i][j]*f[j];
       }
    }
    
    for(i=0; i<9; i++){
       v[i]=v[i]-prod[i];
    }
    
}

// Results

for(i=0; i<9; i++){
   cout<<v[i]<<endl;
}
       
    system("PAUSE");
    return EXIT_SUCCESS;
}
