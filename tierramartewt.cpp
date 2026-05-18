#include <cstdlib>
#include <iostream>
#include <math.h>

// Simulación viaje Tierra-Marte

using namespace std;

using namespace std;

double r, s, a, ba, mu, el[3];

int i, j, k, l;

double a0, f0, g0, df0, dg0, r_1[3], ra1, ra2, v_1[3];

double dr1_dx1, dra_dx, deta3_da1, df1_da1, dg1_da1, ddf0_dx1, ddg0_dx1, da1_dx1, dvr1_dx1[3],dr_dx1[3], f1, g1, dv1_dx1[3];

double da1_dv1[3], dr_dt1[3], deta3_dx1, deta3_dt1, dr_dx[3], dra_da1, dra_dx1, dra_dv1[3];

double ddeta3_dx1x, ddeta3_dt1x, ddr_dxx[3], ddr_dx1x[3], ddr_dt1x[3], ddr_dx1t[3], ddeta3_da1a1, ddeta3_dx1a1;

double ddeta3_dt1a1, ddr_dxa1[3], ddf1_da1a1, ddg1_da1a1, ddr_dx1a1[3], ddr_dt1a1[3], ddeta3_dx1v1[3], ddeta3_dt1v1[3], ddr_dxv1[3][3], ddr_dx1v1[3][3];

double ddr_dt1v1[3][3], ddr1_dx1x1, ddvr1_dx1x1[3], dddf0_dx1x1, dddg0_dx1x1, ddv1_dx1x1[3], ddeta3_dx1x1, ddeta3_dt1x1, ddr_dx1x1[3], dda1_dx1x1;

double ddr_dt1x1[3], ddr_dxx1[3], ddeta3_dv1x1[3], ddr_dv1a1[3][3], ddr_dv1x1[3][3], df9_da1;
   
    double r0[3];
    double br0[3];
    double v0[3];
    double bv0[3];
    double v[11];
    double vc[11];
    double f[11];
    double Jf[11][11];
    double b[11];
    double id[11][11];
    double c[11];
    double prod[11];
    
    double h[3];
    double e[3];
    double p, bo;

int main(int argc, char *argv[])
{

// Condiciones iniciales
    
    cout.precision(15);
    
    r0[0]=-31176357.946332/149597870.7;
    r0[1]=143712337.858123/149597870.7;
    r0[2]=0/149597870.7;
    
    br0[0]=-237096018.783728/149597870.7;
    br0[1]=75813928.3880478/149597870.7;
    br0[2]=7373977.04244497/149597870.7;
    
    v0[0]=-29.6088726252273*(5.0227*pow(10,6))/149597870.7;
    v0[1]=-6.41400258755387*(5.0227*pow(10,6))/149597870.7;
    v0[2]=0*(5.0227*pow(10,6))/149597870.7;
    
    bv0[0]=-6.40409080768176*(5.0227*pow(10,6))/149597870.7;
    bv0[1]=-21.0389054159367*(5.0227*pow(10,6))/149597870.7;
    bv0[2]=-0.292202535710839*(5.0227*pow(10,6))/149597870.7;
    
    
    ba=(1.524);
    mu=1.32715*pow(10,11)*((5.0227*pow(10,6))*(5.0227*pow(10,6)))/(149597870.7*149597870.7*149597870.7);
    
// Elementos orbitas iniciales
    
    el[0]=((((((v0[0])*(v0[0]))+((v0[1])*(v0[1]))+((v0[2])*(v0[2])))/mu)-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*r0[0]) - (((r0[0]*(v0[0]))+(r0[1]*(v0[1]))+(r0[2]*(v0[2])))*v0[0]/mu);
    el[1]=((((((v0[0])*(v0[0]))+((v0[1])*(v0[1]))+((v0[2])*(v0[2])))/mu)-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*r0[1]) - (((r0[0]*(v0[0]))+(r0[1]*(v0[1]))+(r0[2]*(v0[2])))*v0[1]/mu);
    el[2]=((((((v0[0])*(v0[0]))+((v0[1])*(v0[1]))+((v0[2])*(v0[2])))/mu)-(1/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*r0[2]) - (((r0[0]*(v0[0]))+(r0[1]*(v0[1]))+(r0[2]*(v0[2])))*v0[2]/mu);

    for(i=0; i<3; i++){
       cout<<"1e"<<i<<"= "<<el[i]<<endl;
    }

    for(i=0; i<3; i++){
       cout<<"1d"<<i<<"= "<<-el[i]*1<<endl;
    }

    el[0]=((((((bv0[0])*(bv0[0]))+((bv0[1])*(bv0[1]))+((bv0[2])*(bv0[2])))/mu)-(1/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))))*br0[0]) - (((br0[0]*(bv0[0]))+(br0[1]*(bv0[1]))+(br0[2]*(bv0[2])))*bv0[0]/mu);
    el[1]=((((((bv0[0])*(bv0[0]))+((bv0[1])*(bv0[1]))+((bv0[2])*(bv0[2])))/mu)-(1/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))))*br0[1]) - (((br0[0]*(bv0[0]))+(br0[1]*(bv0[1]))+(br0[2]*(bv0[2])))*bv0[1]/mu);
    el[2]=((((((bv0[0])*(bv0[0]))+((bv0[1])*(bv0[1]))+((bv0[2])*(bv0[2])))/mu)-(1/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))))*br0[2]) - (((br0[0]*(bv0[0]))+(br0[1]*(bv0[1]))+(br0[2]*(bv0[2])))*bv0[2]/mu);

    for(i=0; i<3; i++){
       cout<<"2e"<<i<<"= "<<el[i]<<endl;
    }

    for(i=0; i<3; i++){
       cout<<"2d"<<i<<"= "<<-el[i]*1.524<<endl;
    }

// Semilla
    
    v[0]=3.96;
    v[1]=4.42;
    v[2]=5.97;
    v[3]=-0.06;
    v[4]=-0.06;
    v[5]=0.04;
    
    v[9]=0.4;
    v[10]=0.39;
    
    a0=(-1)/(((((v0[0])*(v0[0]))+((v0[1])*(v0[1]))+((v0[2])*(v0[2])))/mu)-(2/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))));
    f0=1-((a0/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*(1-cos(v[9]/sqrt(a0))));
    g0=v[10] - ((a0/sqrt(mu))*(v[9]-(sqrt(a0)*sin(v[9]/sqrt(a0)))));
    r_1[0]=(f0*r0[0])+(g0*v0[0]);
    r_1[1]=(f0*r0[1])+(g0*v0[1]);
    r_1[2]=(f0*r0[2])+(g0*v0[2]);
    ra1=sqrt(pow(r_1[0],2)+pow(r_1[1],2)+pow(r_1[2],2));
    ra1=(a0+(sqrt(a0/mu)*((r0[0]*(v0[0]))+(r0[1]*(v0[1]))+(r0[2]*(v0[2])))*sin(v[9]/sqrt(a0)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a0)*cos(v[9]/sqrt(a0))));
    df0=-((sqrt(mu*a0)/(ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*sin(v[9]/sqrt(a0)));
    dg0=1-((a0/ra1)*(1-cos(v[9]/sqrt(a0))));
    v_1[0]=(df0*r0[0])+(dg0*v0[0]);
    v_1[1]=(df0*r0[1])+(dg0*v0[1]);
    v_1[2]=(df0*r0[2])+(dg0*v0[2]);
    a=(-1)/(((((v_1[0]+v[3])*(v_1[0]+v[3]))+((v_1[1]+v[4])*(v_1[1]+v[4]))+((v_1[2]+v[5])*(v_1[2]+v[5])))/mu)-(2/sqrt((r_1[0]*r_1[0])+(r_1[1]*r_1[1])+(r_1[2]*r_1[2]))));
    
//calculo Hohmann

    cout<<"Av="<<sqrt(2*((mu/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))-(mu/(sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))+6968.1363)))) - sqrt(mu/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))<<"  "<<"t="<< 3.14159265358979323846*sqrt((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))+6968.1363)*(sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))+6968.1363)*(sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))+6968.1363)/(8*mu)) <<endl;

// Calculo primeros multiplicadores 

    double phi[3][3];
    double b1[3];
    double id1[3][3];
    double c1[3]; 
    double r1, s1, prod1[3];  
    
    
    phi[0][0]=( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[0]+v[3]))) ) + (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a)))))))  -  ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])))*((2*a*a/mu*(v_1[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0])));
    phi[0][1]=( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[0]+v[3]))) ) - ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])))*((2*a*a/mu*(v_1[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1])));
    phi[0][2]=( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[0]+v[3]))) ) - ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])))*((2*a*a/mu*(v_1[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2])));
    
    phi[1][0]=( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[1]+v[4]))) ) - ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])))*((2*a*a/mu*(v_1[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0])));
    phi[1][1]=( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[1]+v[4]))) ) + (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a)))))))  -  ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])))*((2*a*a/mu*(v_1[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1])));
    phi[1][2]=( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[1]+v[4]))) ) - ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])))*((2*a*a/mu*(v_1[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2])));
    
    phi[2][0]=( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[2]+v[5]))) ) - ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])))*((2*a*a/mu*(v_1[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0])));
    phi[2][1]=( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[2]+v[5]))) ) - ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])))*((2*a*a/mu*(v_1[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1])));
    phi[2][2]=( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*cos(v[1]/sqrt(a)))+(v[1]/2*sin(v[1]/sqrt(a))))*(v_1[2]+v[5]))) ) + (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a)))))))  -  ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])))*((2*a*a/mu*(v_1[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) ))+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2])));

// inversa
    
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
          if(i==j){
             id1[i][j]=1;
          }else{
             id1[i][j]=0;
          }
      }
    }
    
    for(i=0; i<3; i++){
             
          if(fabs(phi[i][i])<pow(10,-8)){
             for(k=0; k<3; k++ ){
                if(k!=i && fabs(phi[k][i])>pow(10,-8)){
                   for(j=0; j<3; j++){
                      b1[j]=phi[i][j];
                         c1[j]=id1[i][j];
                      phi[i][j]=phi[k][j];
                         id1[i][j]=id1[k][j];
                      phi[k][j]=b1[j];
                         id1[k][j]=c1[j]; 
                   }
                }
             }
          }       
            
             
       s1=phi[i][i];
       for(j=0; j<3; j++){
          phi[i][j]=phi[i][j]/s1;
             id1[i][j]=id1[i][j]/s1;
          b1[j]=phi[i][j];
             c1[j]=id1[i][j];  
       }
       
       for(k=0; k<3; k++){
          if(k!=i){   
             r1=phi[k][i];
             for(j=0; j<3; j++){
                phi[k][j]=phi[k][j]-(r1*b1[j]);
                   id1[k][j]=id1[k][j]-(r1*c1[j]);
             }
          }
       }           
    }
    
     for(i=0; i<3; i++){
       prod1[i]=0;
      for(j=0; j<3; j++){
         prod1[i]=prod1[i]+(v[j+3]*id1[j][i]);   
      }
    }
        
    v[6]=prod1[0];
    v[7]=prod1[1];
    v[8]=prod1[2];
    
//    for(i=0; i<3; i++){
//       cout<<v[i+6]<<endl;
//    }
    
    vc[0]=0;
    vc[1]=0;
    vc[2]=0;
    vc[3]=0;
    vc[4]=0;
    vc[5]=0;
    vc[6]=0;
    vc[7]=0;
    vc[8]=0;
    vc[9]=0;
    vc[10]=0;  
    
    l=0;

// Metodo de Newton
    
  while (sqrt(((vc[0]-v[0])*(vc[0]-v[0]))+((vc[1]-v[1])*(vc[1]-v[1]))+((vc[2]-v[2])*(vc[2]-v[2]))+((vc[3]-v[3])*(vc[3]-v[3]))+((vc[4]-v[4])*(vc[4]-v[4]))+((vc[5]-v[5])*(vc[5]-v[5]))+((vc[6]-v[6])*(vc[6]-v[6]))+((vc[7]-v[7])*(vc[7]-v[7]))+((vc[8]-v[8])*(vc[8]-v[8]))+((vc[9]-v[9])*(vc[9]-v[9]))+((vc[10]-v[10])*(vc[10]-v[10])))>pow(10,-12)){
    
    for(i=0; i<11; i++){
       vc[i]=v[i];
    }

// Variables Especiales    
    
    f0=1-((a0/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*(1-cos(v[9]/sqrt(a0))));
    g0=v[10] - ((a0/sqrt(mu))*(v[9]-(sqrt(a0)*sin(v[9]/sqrt(a0)))));
    r_1[0]=(f0*r0[0])+(g0*v0[0]);
    r_1[1]=(f0*r0[1])+(g0*v0[1]);
    r_1[2]=(f0*r0[2])+(g0*v0[2]);
    //ra1=sqrt(pow(r_1[0],2)+pow(r_1[1],2)+pow(r_1[2],2));
    ra1=(a0+(sqrt(a0/mu)*((r0[0]*(v0[0]))+(r0[1]*(v0[1]))+(r0[2]*(v0[2])))*sin(v[9]/sqrt(a0)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a0)*cos(v[9]/sqrt(a0))));
    df0=-((sqrt(mu*a0)/(ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*sin(v[9]/sqrt(a0)));
    dg0=1-((a0/ra1)*(1-cos(v[9]/sqrt(a0))));
    v_1[0]=(df0*r0[0])+(dg0*v0[0]);
    v_1[1]=(df0*r0[1])+(dg0*v0[1]);
    v_1[2]=(df0*r0[2])+(dg0*v0[2]);
    a=(-1)/(((((v_1[0]+v[3])*(v_1[0]+v[3]))+((v_1[1]+v[4])*(v_1[1]+v[4]))+((v_1[2]+v[5])*(v_1[2]+v[5])))/mu)-(2/sqrt((r_1[0]*r_1[0])+(r_1[1]*r_1[1])+(r_1[2]*r_1[2]))));
    
    f1=1-(a/ra1*(1-cos(v[1]/sqrt(a))));
    g1=(v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a)))));
    
    ra2=a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)));
    
    dr1_dx1=(((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(mu)*cos(v[9]/sqrt(a0))) - (sqrt(a0)*((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/a0)-1)*sin(v[9]/sqrt(a0)));
	
    dra_dx=(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[3]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a))) - (sqrt(a)*((ra1/a)-1)*sin(v[1]/sqrt(a)));
    
    deta3_da1=( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) );
    
    df1_da1=-(1/ra1)*(1-cos(v[1]/sqrt(a))-((v[1]/(2*sqrt(a)))*sin(v[1]/sqrt(a))));
    dg1_da1=-(1/sqrt(mu))*(v[1]-((3*sqrt(a)/2)*sin(v[1]/sqrt(a)))+((v[1]/2)*cos(v[1]/sqrt(a))));
    
    ddf0_dx1=-(df0*dr1_dx1/ra1)-(sqrt(mu)*cos(v[9]/sqrt(a0))/(ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))));
    ddg0_dx1=((1-dg0)*dr1_dx1/ra1)-(sqrt(a0)*sin(v[9]/sqrt(a0))/ra1);
    
    dvr1_dx1[0]=-((sqrt(a0)*sin(v[9]/sqrt(a0))/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*r0[0])-(a0/sqrt(mu)*(1-cos(v[9]/sqrt(a0)))*v0[0]);
    dvr1_dx1[1]=-((sqrt(a0)*sin(v[9]/sqrt(a0))/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*r0[1])-(a0/sqrt(mu)*(1-cos(v[9]/sqrt(a0)))*v0[1]);
    dvr1_dx1[2]=-((sqrt(a0)*sin(v[9]/sqrt(a0))/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*r0[2])-(a0/sqrt(mu)*(1-cos(v[9]/sqrt(a0)))*v0[2]);
    
    dv1_dx1[0]=(ddf0_dx1*r0[0])+(ddg0_dx1*v0[0]);
    dv1_dx1[1]=(ddf0_dx1*r0[1])+(ddg0_dx1*v0[1]);
    dv1_dx1[2]=(ddf0_dx1*r0[2])+(ddg0_dx1*v0[2]);
	
    
    da1_dx1=(2*a*a/mu)*((((v_1[0]+v[3])*dv1_dx1[0])+((v_1[1]+v[4])*dv1_dx1[1])+((v_1[2]+v[5])*dv1_dx1[2]))+(mu*dr1_dx1/(ra1*ra1)));
    
    dr_dx1[0]=(((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*(da1_dx1))+(f1*dvr1_dx1[0])+(g1*dv1_dx1[0])+(a*(1-cos(v[1]/sqrt(a)))*dr1_dx1*r_1[0]/(ra1*ra1));
    dr_dx1[1]=(((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*(da1_dx1))+(f1*dvr1_dx1[1])+(g1*dv1_dx1[1])+(a*(1-cos(v[1]/sqrt(a)))*dr1_dx1*r_1[1]/(ra1*ra1));
    dr_dx1[2]=(((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*(da1_dx1))+(f1*dvr1_dx1[2])+(g1*dv1_dx1[2])+(a*(1-cos(v[1]/sqrt(a)))*dr1_dx1*r_1[2]/(ra1*ra1));
    
    da1_dv1[0]=(2*a*a/mu)*(v_1[0]+v[3]);
    da1_dv1[1]=(2*a*a/mu)*(v_1[1]+v[4]);
    da1_dv1[2]=(2*a*a/mu)*(v_1[2]+v[5]);
    
    deta3_dx1=(deta3_da1*da1_dx1)+(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*(((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2]))+(((v_1[0]+v[3])*dvr1_dx1[0])+((v_1[1]+v[4])*dvr1_dx1[1])+((v_1[2]+v[5])*dvr1_dx1[2]))))+(sqrt(a)*sin(v[1]/sqrt(a))*dr1_dx1);
    
    deta3_dt1=(a/mu*(1-cos(v[1]/sqrt(a)))*(((v_1[0]+v[3])*v0[0])+((v_1[1]+v[4])*v0[1])+((v_1[2]+v[5])*v0[2])))+sqrt(mu);
    
    dr_dt1[0]=f1*v0[0]-(v_1[0]+v[3]);
    dr_dt1[1]=f1*v0[1]-(v_1[1]+v[4]);
    dr_dt1[2]=f1*v0[2]-(v_1[2]+v[5]);
    
    dr_dx[0]=-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0])-(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*(v_1[0]+v[3]));
    dr_dx[1]=-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1])-(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*(v_1[1]+v[4]));
    dr_dx[2]=-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2])-(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*(v_1[2]+v[5]));
    
    dra_da1=1 + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*((sin(v[1]/sqrt(a))/(2*sqrt(a*mu)))-(v[1]*cos(v[1]/sqrt(a))/(2*a*sqrt(mu))))) - cos(v[1]/sqrt(a)) + ((ra1-a)*sin(v[1]/sqrt(a))/(2*a*sqrt(a)));
    
    dra_dx1=(dra_da1*da1_dx1) + (a/sqrt(a*mu)*sin(v[1]/sqrt(a))*(((dvr1_dx1[0]*(v_1[0]+v[3]))+(dvr1_dx1[1]*(v_1[1]+v[4]))+(dvr1_dx1[2]*(v_1[2]+v[5])))+((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2])))) + (dr1_dx1*cos(v[1]/sqrt(a)));
    
    dra_dv1[0]=(dra_da1*da1_dv1[0]) + (sin(v[1]/sqrt(a))*sqrt(a/mu)*r_1[0]);
    dra_dv1[1]=(dra_da1*da1_dv1[1]) + (sin(v[1]/sqrt(a))*sqrt(a/mu)*r_1[1]);
    dra_dv1[2]=(dra_da1*da1_dv1[2]) + (sin(v[1]/sqrt(a))*sqrt(a/mu)*r_1[2]);
    
    
    ddeta3_dx1x=((1 -cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*((sin(v[1]/sqrt(a))/(2*sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (ra1/2*(v[1]*sin(v[1]/sqrt(a))/(sqrt(a)*a))) )*da1_dx1) + (a/sqrt(mu)*(sin(v[1]/sqrt(a))/sqrt(a))*(((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2]))+(((v_1[0]+v[3])*dvr1_dx1[0])+((v_1[1]+v[4])*dvr1_dx1[1])+((v_1[2]+v[5])*dvr1_dx1[2])))) + (cos(v[1]/sqrt(a))*dr1_dx1);
    
    ddeta3_dt1x=(1/sqrt(mu))*sin(v[1]/sqrt(a))*(((v_1[0]+v[3])*v0[0])+((v_1[1]+v[4])*v0[1])+((v_1[2]+v[5])*v0[2]));
    
    ddr_dxx[0]=-((1/ra1)*cos(v[1]/sqrt(a))*r_1[0])-(sqrt(a/mu)*sin(v[1]/sqrt(a))*(v_1[0]+v[3]));
    ddr_dxx[1]=-((1/ra1)*cos(v[1]/sqrt(a))*r_1[1])-(sqrt(a/mu)*sin(v[1]/sqrt(a))*(v_1[1]+v[4]));
    ddr_dxx[2]=-((1/ra1)*cos(v[1]/sqrt(a))*r_1[2])-(sqrt(a/mu)*sin(v[1]/sqrt(a))*(v_1[2]+v[5]));
    
    ddr_dx1x[0]=((( (-(1/ra1)*((sin(v[1]/sqrt(a))/(2*sqrt(a)))-((v[1]/(2*a))*cos(v[1]/sqrt(a))))) *r_1[0])+( (-(1/sqrt(mu))*(1-cos(v[1]/sqrt(a))-((v[1]/(2*sqrt(a)))*sin(v[1]/sqrt(a))))) *(v_1[0]+v[3])))*(da1_dx1)) - ( (sqrt(a)/ra1*sin(v[1]/sqrt(a))) *dvr1_dx1[0])+( (a/sqrt(mu)*cos(v[1]/sqrt(a))) *dv1_dx1[0])+(sqrt(a)*sin(v[1]/sqrt(a))*dr1_dx1*r_1[0]/(ra1*ra1));
    ddr_dx1x[1]=((( (-(1/ra1)*((sin(v[1]/sqrt(a))/(2*sqrt(a)))-((v[1]/(2*a))*cos(v[1]/sqrt(a))))) *r_1[1])+( (-(1/sqrt(mu))*(1-cos(v[1]/sqrt(a))-((v[1]/(2*sqrt(a)))*sin(v[1]/sqrt(a))))) *(v_1[1]+v[4])))*(da1_dx1)) - ( (sqrt(a)/ra1*sin(v[1]/sqrt(a))) *dvr1_dx1[1])+( (a/sqrt(mu)*cos(v[1]/sqrt(a))) *dv1_dx1[1])+(sqrt(a)*sin(v[1]/sqrt(a))*dr1_dx1*r_1[1]/(ra1*ra1));
    ddr_dx1x[2]=((( (-(1/ra1)*((sin(v[1]/sqrt(a))/(2*sqrt(a)))-((v[1]/(2*a))*cos(v[1]/sqrt(a))))) *r_1[2])+( (-(1/sqrt(mu))*(1-cos(v[1]/sqrt(a))-((v[1]/(2*sqrt(a)))*sin(v[1]/sqrt(a))))) *(v_1[2]+v[5])))*(da1_dx1)) - ( (sqrt(a)/ra1*sin(v[1]/sqrt(a))) *dvr1_dx1[2])+( (a/sqrt(mu)*cos(v[1]/sqrt(a))) *dv1_dx1[2])+(sqrt(a)*sin(v[1]/sqrt(a))*dr1_dx1*r_1[2]/(ra1*ra1));
    
    ddr_dt1x[0]=-sqrt(a)/ra1*sin(v[1]/sqrt(a))*v0[0];
    ddr_dt1x[1]=-sqrt(a)/ra1*sin(v[1]/sqrt(a))*v0[1];
    ddr_dt1x[2]=-sqrt(a)/ra1*sin(v[1]/sqrt(a))*v0[2];
    
    ddr_dx1t[0]=dv1_dx1[0];
    ddr_dx1t[1]=dv1_dx1[1];
    ddr_dx1t[2]=dv1_dx1[2];
    
    ddeta3_da1a1= -(3/(4*sqrt(a))*sin(v[1]/sqrt(a))) + (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((v[1]*v[1])/(4*a*sqrt(a))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*a*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*(-(sin(v[1]/sqrt(a))/(2*a*sqrt(a)))-(v[1]*v[1]*sin(v[1]/sqrt(a))/(2*a*a*sqrt(a)))+(v[1]/(2*a*a)*cos(v[1]/sqrt(a))))) ;
    
    ddeta3_dx1a1=(ddeta3_da1a1*da1_dx1) + (deta3_da1*da1_dx1*(2/a)) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))) - (v[1]*sin(v[1]/sqrt(a))/(2*sqrt(a*mu))) )*(((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2]))+(((v_1[0]+v[3])*dvr1_dx1[0])+((v_1[1]+v[4])*dvr1_dx1[1])+((v_1[2]+v[5])*dvr1_dx1[2])))) + (((sin(v[1]/sqrt(a))/(2*sqrt(a))) - (v[1]*cos(v[1]/sqrt(a))/(2*a)))*dr1_dx1);
    
    ddeta3_dt1a1=(((v_1[0]+v[3])*v0[0])+((v_1[1]+v[4])*v0[1])+((v_1[2]+v[5])*v0[2]))*((1-cos(v[1]/sqrt(a)))-(a*v[1]*sin(v[1]/sqrt(a))/(2*a*sqrt(a))))/sqrt(mu);
    
    ddr_dxa1[0]=((-(1/(2*ra1*sqrt(a))*sin(v[1]/sqrt(a))) + (v[1]*cos(v[1]/sqrt(a))/(2*ra1*a)))*r_1[0])-(((1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))) - (v[1]*sin(v[1]/sqrt(a))/(2*sqrt(a*mu))))*(v_1[0]+v[3]));
    ddr_dxa1[1]=((-(1/(2*ra1*sqrt(a))*sin(v[1]/sqrt(a))) + (v[1]*cos(v[1]/sqrt(a))/(2*ra1*a)))*r_1[1])-(((1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))) - (v[1]*sin(v[1]/sqrt(a))/(2*sqrt(a*mu))))*(v_1[1]+v[4]));
    ddr_dxa1[2]=((-(1/(2*ra1*sqrt(a))*sin(v[1]/sqrt(a))) + (v[1]*cos(v[1]/sqrt(a))/(2*ra1*a)))*r_1[2])-(((1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))) - (v[1]*sin(v[1]/sqrt(a))/(2*sqrt(a*mu))))*(v_1[2]+v[5]));
    
    ddf1_da1a1=-(1/ra1)*((-v[1]*sin(v[1]/sqrt(a))/(4*a*sqrt(a))) + (v[1]*v[1]*cos(v[1]/sqrt(a))/(4*a*a)));
    ddg1_da1a1=-(1/sqrt(mu))*((v[1]*v[1]*sin(v[1]/sqrt(a))/(4*a*sqrt(a))) - (3*sin(v[1]/sqrt(a))/(4*sqrt(a))) + (3*v[1]*cos(v[1]/sqrt(a))/(4*a)));
    
    ddr_dx1a1[0]=(((ddf1_da1a1*r_1[0]) + (ddg1_da1a1*(v_1[1]+v[3])))*da1_dx1) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*(da1_dx1*2/a)) + (df1_da1*dvr1_dx1[0]) + (dg1_da1*dv1_dx1[0]) + ((1-cos(v[1]/sqrt(a)))*dr1_dx1*r_1[0]/(ra1*ra1)) - (v[1]*sin(v[1]/sqrt(a))*dr1_dx1*r_1[0]/(2*ra1*ra1*sqrt(a)));
    ddr_dx1a1[1]=(((ddf1_da1a1*r_1[1]) + (ddg1_da1a1*(v_1[1]+v[4])))*da1_dx1) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*(da1_dx1*2/a)) + (df1_da1*dvr1_dx1[1]) + (dg1_da1*dv1_dx1[1]) + ((1-cos(v[1]/sqrt(a)))*dr1_dx1*r_1[1]/(ra1*ra1)) - (v[1]*sin(v[1]/sqrt(a))*dr1_dx1*r_1[1]/(2*ra1*ra1*sqrt(a)));
    ddr_dx1a1[2]=(((ddf1_da1a1*r_1[2]) + (ddg1_da1a1*(v_1[2]+v[5])))*da1_dx1) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*(da1_dx1*2/a)) + (df1_da1*dvr1_dx1[2]) + (dg1_da1*dv1_dx1[2]) + ((1-cos(v[1]/sqrt(a)))*dr1_dx1*r_1[2]/(ra1*ra1)) - (v[1]*sin(v[1]/sqrt(a))*dr1_dx1*r_1[2]/(2*ra1*ra1*sqrt(a)));
    
    ddr_dt1a1[0]=(df1_da1*v0[0])-(v_1[0]-v[3]); 
    ddr_dt1a1[1]=(df1_da1*v0[1])-(v_1[1]-v[4]);
    ddr_dt1a1[2]=(df1_da1*v0[2])-(v_1[2]-v[5]);
    
    ddeta3_dx1v1[0]=(ddeta3_dx1a1*da1_dv1[0]) + ((1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*da1_dx1*r_1[0]/sqrt(mu)) + (deta3_da1*(2*a*a*dv1_dx1[0])/mu) +(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dvr1_dx1[0]);
    ddeta3_dx1v1[1]=(ddeta3_dx1a1*da1_dv1[1]) + ((1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*da1_dx1*r_1[1]/sqrt(mu)) + (deta3_da1*(2*a*a*dv1_dx1[1])/mu) +(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dvr1_dx1[1]);
    ddeta3_dx1v1[2]=(ddeta3_dx1a1*da1_dv1[2]) + ((1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*da1_dx1*r_1[2]/sqrt(mu)) + (deta3_da1*(2*a*a*dv1_dx1[2])/mu) +(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dvr1_dx1[2]);
    
    ddeta3_dt1v1[0]=(ddeta3_dt1a1*da1_dv1[0]) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v0[0]);
    ddeta3_dt1v1[1]=(ddeta3_dt1a1*da1_dv1[1]) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v0[1]);
    ddeta3_dt1v1[2]=(ddeta3_dt1a1*da1_dv1[2]) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v0[2]);
    
    ddr_dxv1[0][0]=(ddr_dxa1[0]*da1_dv1[0])-(a*(1-cos(v[1]/sqrt(a)))/sqrt(mu));
    ddr_dxv1[0][1]=(ddr_dxa1[0]*da1_dv1[1]);
    ddr_dxv1[0][2]=(ddr_dxa1[0]*da1_dv1[2]);
    ddr_dxv1[1][0]=(ddr_dxa1[1]*da1_dv1[0]);
    ddr_dxv1[1][1]=(ddr_dxa1[1]*da1_dv1[1])-(a*(1-cos(v[1]/sqrt(a)))/sqrt(mu));
    ddr_dxv1[1][2]=(ddr_dxa1[1]*da1_dv1[2]);
    ddr_dxv1[2][0]=(ddr_dxa1[2]*da1_dv1[0]);
    ddr_dxv1[2][1]=(ddr_dxa1[2]*da1_dv1[1]);
    ddr_dxv1[2][2]=(ddr_dxa1[2]*da1_dv1[2])-(a*(1-cos(v[1]/sqrt(a)))/sqrt(mu));
    
    ddr_dx1v1[0][0]=(ddr_dx1a1[0]*da1_dv1[0]) + (2*a*a/mu*((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*dv1_dx1[0]) + (dg1_da1*da1_dx1);
    ddr_dx1v1[0][1]=(ddr_dx1a1[0]*da1_dv1[1]) + (2*a*a/mu*((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*dv1_dx1[1]);
    ddr_dx1v1[0][2]=(ddr_dx1a1[0]*da1_dv1[2]) + (2*a*a/mu*((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*dv1_dx1[2]);
    ddr_dx1v1[1][0]=(ddr_dx1a1[1]*da1_dv1[0]) + (2*a*a/mu*((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*dv1_dx1[0]);
    ddr_dx1v1[1][1]=(ddr_dx1a1[1]*da1_dv1[1]) + (2*a*a/mu*((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*dv1_dx1[1]) + (dg1_da1*da1_dx1);
    ddr_dx1v1[1][2]=(ddr_dx1a1[1]*da1_dv1[2]) + (2*a*a/mu*((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*dv1_dx1[2]);
    ddr_dx1v1[2][0]=(ddr_dx1a1[2]*da1_dv1[0]) + (2*a*a/mu*((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*dv1_dx1[0]);
    ddr_dx1v1[2][1]=(ddr_dx1a1[2]*da1_dv1[1]) + (2*a*a/mu*((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*dv1_dx1[1]);
    ddr_dx1v1[2][2]=(ddr_dx1a1[2]*da1_dv1[2]) + (2*a*a/mu*((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*dv1_dx1[2]) + (dg1_da1*da1_dx1);
    
    ddr_dt1v1[0][0]=(ddr_dt1a1[0]*da1_dv1[0])-1;
    ddr_dt1v1[0][1]=(ddr_dt1a1[0]*da1_dv1[1]);
    ddr_dt1v1[0][2]=(ddr_dt1a1[0]*da1_dv1[2]);
    ddr_dt1v1[1][0]=(ddr_dt1a1[1]*da1_dv1[0]);
    ddr_dt1v1[1][1]=(ddr_dt1a1[1]*da1_dv1[1])-1;
    ddr_dt1v1[1][2]=(ddr_dt1a1[1]*da1_dv1[2]);
    ddr_dt1v1[2][0]=(ddr_dt1a1[2]*da1_dv1[0]);
    ddr_dt1v1[2][1]=(ddr_dt1a1[2]*da1_dv1[1]);
    ddr_dt1v1[2][2]=(ddr_dt1a1[2]*da1_dv1[2])-1;
    
    ddr1_dx1x1=-(((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))/sqrt(a*mu)*sin(v[9]/sqrt(a0))) - (((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))/a0)-1)*cos(v[9]/sqrt(a0)));
    
    ddvr1_dx1x1[0]=-((cos(v[9]/sqrt(a0))/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*r0[0])-(sqrt(a/mu)*sin(v[9]/sqrt(a0))*v0[0]);
    ddvr1_dx1x1[1]=-((cos(v[9]/sqrt(a0))/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*r0[1])-(sqrt(a/mu)*sin(v[9]/sqrt(a0))*v0[1]);
    ddvr1_dx1x1[2]=-((cos(v[9]/sqrt(a0))/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))*r0[2])-(sqrt(a/mu)*sin(v[9]/sqrt(a0))*v0[2]);
    
    dddf0_dx1x1=-(2*sqrt(a0*mu)*sin(v[9]/sqrt(a0))*(dr1_dx1*dr1_dx1)/(ra1*ra1*ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))) + (sqrt(mu)*cos(v[9]/sqrt(a0))*dr1_dx1/(ra1*ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))) + (sqrt(a0*mu)*sin(v[9]/sqrt(a0))*ddr1_dx1x1/(ra1*ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))) + (sqrt(mu)*cos(v[9]/sqrt(a0))*dr1_dx1/(ra1*ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2])))) + (sqrt(mu)*sin(v[9]/sqrt(a0))/(sqrt(a0)*ra1*sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))));
    
    dddg0_dx1x1=-(2*a0*(1-cos(v[9]/sqrt(a0)))*(dr1_dx1*dr1_dx1)/(ra1*ra1*ra1)) + (sqrt(a0)*sin(v[9]/sqrt(a0))*dr1_dx1/(ra1*ra1)) + (a0*(1-cos(v[9]/sqrt(a0)))*ddr1_dx1x1/(ra1*ra1)) + (sqrt(a0)*sin(v[9]/sqrt(a0))*dr1_dx1/(ra1*ra1)) - (cos(v[9]/sqrt(a0))/(ra1));
    
    ddv1_dx1x1[0]=(dddf0_dx1x1*r0[0])+(dddg0_dx1x1*v0[0]);
    ddv1_dx1x1[1]=(dddf0_dx1x1*r0[1])+(dddg0_dx1x1*v0[1]);
    ddv1_dx1x1[2]=(dddf0_dx1x1*r0[2])+(dddg0_dx1x1*v0[2]);
    
    dda1_dx1x1=((2*a*a/mu)*((((v_1[0]+v[3])*ddv1_dx1x1[0])+((v_1[1]+v[4])*ddv1_dx1x1[1])+((v_1[2]+v[5])*ddv1_dx1x1[2]))+((dv1_dx1[0]*dv1_dx1[0])+(dv1_dx1[1]*dv1_dx1[1])+(dv1_dx1[2]*dv1_dx1[2]))-(2*mu*dr1_dx1*dr1_dx1/(ra1*ra1*ra1))+(mu*ddr1_dx1x1/(ra1*ra1)))) + (da1_dx1*(2/a)*da1_dx1);
    
    ddeta3_dx1x1=(((ddeta3_da1a1*da1_dx1) + (deta3_da1*da1_dx1*(2/a)) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))) - (v[1]*sin(v[1]/sqrt(a))/(2*sqrt(a*mu))))*(((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2]))+(((v_1[0]+v[3])*dvr1_dx1[0])+((v_1[1]+v[4])*dvr1_dx1[1])+((v_1[2]+v[5])*dvr1_dx1[2])))) + (1/(2*sqrt(a))*sin(v[1]/sqrt(a))*dr1_dx1) - (v[1]*cos(v[1]/sqrt(a))*dr1_dx1/(2*a)))*da1_dx1) + (( ((((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2])) + ((dvr1_dx1[0]*(v_1[0]+v[3]))+(dvr1_dx1[1]*(v_1[1]+v[4]))+(dvr1_dx1[2]*(v_1[2]+v[5]))))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (dr1_dx1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )*da1_dx1) + (deta3_da1*dda1_dx1x1) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*(((r_1[0]*ddv1_dx1x1[0])+(r_1[1]*ddv1_dx1x1[1])+(r_1[2]*ddv1_dx1x1[2]))+(2*((dvr1_dx1[0]*dv1_dx1[0])+(dvr1_dx1[1]*dv1_dx1[1])+(dvr1_dx1[2]*dv1_dx1[2])))+(((v_1[0]+v[3])*ddvr1_dx1x1[0])+((v_1[1]+v[4])*ddvr1_dx1x1[1])+((v_1[2]+v[5])*ddvr1_dx1x1[2]))))+(sqrt(a)*sin(v[1]/sqrt(a))*ddr1_dx1x1);
    
    ddeta3_dt1x1=(ddeta3_dt1a1*da1_dx1) + (a/sqrt(mu)*(1-cos(1-cos(v[1]/sqrt(a))))*((v0[0]*dv1_dx1[0])+(v0[1]*dv1_dx1[1])+(v0[2]*dv1_dx1[2])));
    
    ddr_dx1x1[0]=(ddr_dx1a1[0]*da1_dx1)+((-(dr1_dx1/ra1*df1_da1*r_1[0])+(df1_da1*dvr1_dx1[0])+(dg1_da1*dv1_dx1[0]))*da1_dx1)+(((df1_da1*r_1[0])+(dg1_da1*(v_1[0]*v[3])))*dda1_dx1x1) + (f1*ddvr1_dx1x1[0]) + (g1*ddv1_dx1x1[0]) + (2*a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*dvr1_dx1[0]) - (2*a/(ra1*ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*dr1_dx1*r_1[0]) + (a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*ddr1_dx1x1*r_1[0]);
    ddr_dx1x1[1]=(ddr_dx1a1[1]*da1_dx1)+((-(dr1_dx1/ra1*df1_da1*r_1[1])+(df1_da1*dvr1_dx1[1])+(dg1_da1*dv1_dx1[1]))*da1_dx1)+(((df1_da1*r_1[1])+(dg1_da1*(v_1[1]*v[4])))*dda1_dx1x1) + (f1*ddvr1_dx1x1[1]) + (g1*ddv1_dx1x1[1]) + (2*a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*dvr1_dx1[1]) - (2*a/(ra1*ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*dr1_dx1*r_1[1]) + (a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*ddr1_dx1x1*r_1[1]);
    ddr_dx1x1[2]=(ddr_dx1a1[2]*da1_dx1)+((-(dr1_dx1/ra1*df1_da1*r_1[2])+(df1_da1*dvr1_dx1[2])+(dg1_da1*dv1_dx1[2]))*da1_dx1)+(((df1_da1*r_1[2])+(dg1_da1*(v_1[2]*v[5])))*dda1_dx1x1) + (f1*ddvr1_dx1x1[2]) + (g1*ddv1_dx1x1[2]) + (2*a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*dvr1_dx1[2]) - (2*a/(ra1*ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*dr1_dx1*r_1[2]) + (a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*ddr1_dx1x1*r_1[2]);
    
    ddr_dt1x1[0]=(ddr_dt1a1[0]*da1_dx1) + (a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*v0[0])-dv1_dx1[0];
    ddr_dt1x1[1]=(ddr_dt1a1[1]*da1_dx1) + (a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*v0[1])-dv1_dx1[1];
    ddr_dt1x1[2]=(ddr_dt1a1[2]*da1_dx1) + (a/(ra1*ra1)*(1-cos(v[1]/sqrt(a)))*dr1_dx1*v0[2])-dv1_dx1[2];
    
    ddr_dxx1[0]=(ddr_dxa1[0]*da1_dx1) + (sqrt(a)/(ra1*ra1)*sin(v[1]/sqrt(a))*dr1_dx1*r_1[0]) - (sqrt(a)/ra1*sin(v[1])*dvr1_dx1[0]) - (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dv1_dx1[0]);
    ddr_dxx1[1]=(ddr_dxa1[1]*da1_dx1) + (sqrt(a)/(ra1*ra1)*sin(v[1]/sqrt(a))*dr1_dx1*r_1[1]) - (sqrt(a)/ra1*sin(v[1])*dvr1_dx1[1]) - (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dv1_dx1[1]);
    ddr_dxx1[2]=(ddr_dxa1[2]*da1_dx1) + (sqrt(a)/(ra1*ra1)*sin(v[1]/sqrt(a))*dr1_dx1*r_1[2]) - (sqrt(a)/ra1*sin(v[1])*dvr1_dx1[2]) - (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dv1_dx1[2]);
    
    ddeta3_dv1x1[0]=(((ddeta3_da1a1*da1_dv1[0])+(deta3_da1*(4*a/mu*(v_1[0]+v[3])))+(1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*r_1[0])-(v[1]/sqrt(a*mu)*sin(v[1]/sqrt(a))*r_1[0]))*da1_dx1) - (( ((((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2])) + ((dvr1_dx1[0]*(v_1[0]+v[3]))+(dvr1_dx1[1]*(v_1[1]+v[4]))+(dvr1_dx1[2]*(v_1[2]+v[5]))))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (dr1_dx1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )*da1_dv1[0]) + (deta3_da1*(2*a*a/mu*dv1_dx1[0])) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dvr1_dx1[0]);
    ddeta3_dv1x1[1]=(((ddeta3_da1a1*da1_dv1[1])+(deta3_da1*(4*a/mu*(v_1[1]+v[4])))+(1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*r_1[1])-(v[1]/sqrt(a*mu)*sin(v[1]/sqrt(a))*r_1[1]))*da1_dx1) - (( ((((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2])) + ((dvr1_dx1[0]*(v_1[0]+v[3]))+(dvr1_dx1[1]*(v_1[1]+v[4]))+(dvr1_dx1[2]*(v_1[2]+v[5]))))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (dr1_dx1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )*da1_dv1[1]) + (deta3_da1*(2*a*a/mu*dv1_dx1[1])) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dvr1_dx1[1]);
    ddeta3_dv1x1[2]=(((ddeta3_da1a1*da1_dv1[2])+(deta3_da1*(4*a/mu*(v_1[2]+v[5])))+(1/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*r_1[2])-(v[1]/sqrt(a*mu)*sin(v[1]/sqrt(a))*r_1[2]))*da1_dx1) - (( ((((r_1[0]*dv1_dx1[0])+(r_1[1]*dv1_dx1[1])+(r_1[2]*dv1_dx1[2])) + ((dvr1_dx1[0]*(v_1[0]+v[3]))+(dvr1_dx1[1]*(v_1[1]+v[4]))+(dvr1_dx1[2]*(v_1[2]+v[5]))))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (dr1_dx1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )*da1_dv1[2]) + (deta3_da1*(2*a*a/mu*dv1_dx1[2])) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*dvr1_dx1[2]);
    
    ddr_dv1a1[0][0]=(((ddf1_da1a1*r_1[0])+(ddg1_da1a1*(v_1[0]+v[3])))*da1_dv1[0]) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*(4*a/mu*(v_1[0]+v[3]))) + dg1_da1;
    ddr_dv1a1[0][1]=(((ddf1_da1a1*r_1[0])+(ddg1_da1a1*(v_1[0]+v[3])))*da1_dv1[1]) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*(4*a/mu*(v_1[1]+v[4])));
    ddr_dv1a1[0][2]=(((ddf1_da1a1*r_1[0])+(ddg1_da1a1*(v_1[0]+v[3])))*da1_dv1[2]) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*(4*a/mu*(v_1[2]+v[5])));
    ddr_dv1a1[1][0]=(((ddf1_da1a1*r_1[1])+(ddg1_da1a1*(v_1[1]+v[4])))*da1_dv1[0]) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*(4*a/mu*(v_1[0]+v[3])));
    ddr_dv1a1[1][1]=(((ddf1_da1a1*r_1[1])+(ddg1_da1a1*(v_1[1]+v[4])))*da1_dv1[1]) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*(4*a/mu*(v_1[1]+v[4]))) + dg1_da1;
    ddr_dv1a1[1][2]=(((ddf1_da1a1*r_1[1])+(ddg1_da1a1*(v_1[1]+v[4])))*da1_dv1[2]) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*(4*a/mu*(v_1[2]+v[5])));
    ddr_dv1a1[2][0]=(((ddf1_da1a1*r_1[2])+(ddg1_da1a1*(v_1[2]+v[5])))*da1_dv1[0]) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*(4*a/mu*(v_1[0]+v[3])));
    ddr_dv1a1[2][1]=(((ddf1_da1a1*r_1[2])+(ddg1_da1a1*(v_1[2]+v[5])))*da1_dv1[1]) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*(4*a/mu*(v_1[1]+v[4])));
    ddr_dv1a1[2][2]=(((ddf1_da1a1*r_1[2])+(ddg1_da1a1*(v_1[2]+v[5])))*da1_dv1[2]) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*(4*a/mu*(v_1[2]+v[5]))) + dg1_da1;
    
    ddr_dv1x1[0][0]=(ddr_dv1a1[0][0]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[0])+(df1_da1*dvr1_dx1[0])+(dg1_da1*dv1_dx1[0]))*da1_dv1[0]) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]*v[3])))*(2*a*a/mu*dv1_dx1[0]));
    ddr_dv1x1[0][1]=(ddr_dv1a1[0][1]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[0])+(df1_da1*dvr1_dx1[0])+(dg1_da1*dv1_dx1[0]))*da1_dv1[1]) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]*v[3])))*(2*a*a/mu*dv1_dx1[1]));
    ddr_dv1x1[0][2]=(ddr_dv1a1[0][2]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[0])+(df1_da1*dvr1_dx1[0])+(dg1_da1*dv1_dx1[0]))*da1_dv1[2]) + (((df1_da1*r_1[0])+(dg1_da1*(v_1[0]*v[3])))*(2*a*a/mu*dv1_dx1[2]));
    
    ddr_dv1x1[1][0]=(ddr_dv1a1[1][0]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[1])+(df1_da1*dvr1_dx1[1])+(dg1_da1*dv1_dx1[1]))*da1_dv1[0]) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]*v[4])))*(2*a*a/mu*dv1_dx1[0]));
    ddr_dv1x1[1][1]=(ddr_dv1a1[1][1]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[1])+(df1_da1*dvr1_dx1[1])+(dg1_da1*dv1_dx1[1]))*da1_dv1[1]) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]*v[4])))*(2*a*a/mu*dv1_dx1[1]));
    ddr_dv1x1[1][2]=(ddr_dv1a1[1][2]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[1])+(df1_da1*dvr1_dx1[1])+(dg1_da1*dv1_dx1[1]))*da1_dv1[2]) + (((df1_da1*r_1[1])+(dg1_da1*(v_1[1]*v[4])))*(2*a*a/mu*dv1_dx1[2]));
    
    ddr_dv1x1[2][0]=(ddr_dv1a1[2][0]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[2])+(df1_da1*dvr1_dx1[2])+(dg1_da1*dv1_dx1[2]))*da1_dv1[0]) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]*v[5])))*(2*a*a/mu*dv1_dx1[0]));
    ddr_dv1x1[2][1]=(ddr_dv1a1[2][1]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[2])+(df1_da1*dvr1_dx1[2])+(dg1_da1*dv1_dx1[2]))*da1_dv1[1]) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]*v[5])))*(2*a*a/mu*dv1_dx1[1]));
    ddr_dv1x1[2][2]=(ddr_dv1a1[2][2]*da1_dx1) + ((-(dr1_dx1/ra1*df1_da1*r_1[2])+(df1_da1*dvr1_dx1[2])+(dg1_da1*dv1_dx1[2]))*da1_dv1[2]) + (((df1_da1*r_1[2])+(dg1_da1*(v_1[2]*v[5])))*(2*a*a/mu*dv1_dx1[2]));
    
    df9_da1=((-(1/(ra2*ra2)*dra_da1*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)) + (1/ra2*((sqrt(mu)/ra1*ddeta3_dx1a1)+ddeta3_dt1a1)))*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))) + (1/ra2*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((ddr_dxa1[0]*v[6])+(ddr_dxa1[1]*v[7])+(ddr_dxa1[2]*v[8]))) - ((sqrt(mu)/ra1*((ddr_dx1a1[0]*v[6])+(ddr_dx1a1[1]*v[7])+(ddr_dx1a1[2]*v[8]))) + ((ddr_dt1a1[0]*v[6])+(ddr_dt1a1[1]*v[7])+(ddr_dt1a1[2]*v[8])));
    
// 9 ecuaciones (9 raices) 
    
    f[0]=(ba*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))) + (ba/sqrt(mu)*((br0[0]*bv0[0])+(br0[1]*bv0[1])+(br0[2]*bv0[2]))*(1-(cos(v[0]/sqrt(ba))))) + (sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sqrt(ba)*sin(v[0]/sqrt(ba))) - (sqrt(mu)*v[2]);
    f[1]=(a*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))) + (a/sqrt(mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*(1-(cos(v[1]/sqrt(a))))) + (ra1*sqrt(a)*sin(v[1]/sqrt(a))) - (sqrt(mu)*(v[2]-v[10]));
    f[2]=((1-(ba/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*(1-cos(v[0]/sqrt(ba)))))*br0[0]) + ((v[2]-(ba/sqrt(mu)*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))))*bv0[0]) - ((1-(a/ra1*(1-cos(v[1]/sqrt(a)))))*r_1[0]) - (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*(v_1[0]+v[3]));
    f[3]=((1-(ba/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*(1-cos(v[0]/sqrt(ba)))))*br0[1]) + ((v[2]-(ba/sqrt(mu)*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))))*bv0[1]) - ((1-(a/ra1*(1-cos(v[1]/sqrt(a)))))*r_1[1]) - (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*(v_1[1]+v[4]));
    f[4]=((1-(ba/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*(1-cos(v[0]/sqrt(ba)))))*br0[2]) + ((v[2]-(ba/sqrt(mu)*(v[0]-(sqrt(ba)*sin(v[0]/sqrt(ba))))))*bv0[2]) - ((1-(a/ra1*(1-cos(v[1]/sqrt(a)))))*r_1[2]) - (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*(v_1[2]+v[5]));
    f[5]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*((v[6]*br0[0])+(v[7]*br0[1])+(v[8]*br0[2]))) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*((v[6]*bv0[0])+(v[7]*bv0[1])+(v[8]*bv0[2]))) + (sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))) + ((a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5]))));
    
    f[6]=v[3] + ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[0]+v[3])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) ) - (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*v[6]);
    f[7]=v[4] + ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[1]+v[4])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) ) - (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*v[7]);
    f[8]=v[5] + ((1/(a+(sqrt(a/mu)*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))*sin(v[1]/sqrt(a)))+((ra1-a)*cos(v[1]/sqrt(a)))))*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[2]+v[5])*( v[1]-(sqrt(a)*sin(v[1]/sqrt(a))) + (a/2*((v[1]/a*cos(v[1]/sqrt(a)))-(sin(v[1]/sqrt(a))/sqrt(a)))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))) + (ra1/2*((sin(v[1]/sqrt(a))/sqrt(a))-(v[1]*cos(v[1]/sqrt(a))/a))) )) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) ) - (((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))))*v[8]);
    
    f[9]=(a0*(v[9]-(sqrt(a0)*sin(v[9]/sqrt(a0))))) + (a0/sqrt(mu)*((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))*(1-(cos(v[9]/sqrt(a0))))) + (sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))*sqrt(a0)*sin(v[9]/sqrt(a0))) - (sqrt(mu)*v[10]);
    f[10]=((1/ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((v[6]*dr_dx[0])+(v[7]*dr_dx[1])+(v[8]*dr_dx[2])))-((sqrt(mu)/ra1*((v[6]*dr_dx1[0])+(v[7]*dr_dx1[1])+(v[8]*dr_dx1[2])))+((v[6]*dr_dt1[0])+(v[7]*dr_dt1[1])+(v[8]*dr_dt1[2])));
    
// Jacobiano de las ecuaciones 
    
    Jf[0][0]=(ba+(sqrt(ba/mu)*((br0[0]*bv0[0])+(br0[1]*bv0[1])+(br0[2]*bv0[2]))*sin(v[0]/sqrt(ba)))+((sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))-ba)*cos(v[0]/sqrt(ba))));
    Jf[0][1]=0;
    Jf[0][2]=-sqrt(mu);
    Jf[0][3]=0;
    Jf[0][4]=0;
    Jf[0][5]=0;
    Jf[0][6]=0;
    Jf[0][7]=0;
    Jf[0][8]=0;
    Jf[0][9]=0;
    Jf[0][10]=0;
    
    
    Jf[1][0]=0;
    Jf[1][1]=ra2;
    Jf[1][2]=-sqrt(mu);
    Jf[1][3]=(deta3_da1*da1_dv1[0])+(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*r_1[0]);
    Jf[1][4]=(deta3_da1*da1_dv1[1])+(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*r_1[1]);
    Jf[1][5]=(deta3_da1*da1_dv1[2])+(a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*r_1[2]);
    Jf[1][6]=0;
    Jf[1][7]=0;
    Jf[1][8]=0;
    Jf[1][9]=deta3_dx1;
    Jf[1][10]=deta3_dt1;
    
    
    Jf[2][0]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[0]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[0]);
    Jf[2][1]=(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3]));
    Jf[2][2]=bv0[0]-v_1[0]-v[3];
    Jf[2][3]=-(((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*da1_dv1[0])-g1;
    Jf[2][4]=-(((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*da1_dv1[1]);
    Jf[2][5]=-(((df1_da1*r_1[0])+(dg1_da1*(v_1[0]+v[3])))*da1_dv1[2]);
    Jf[2][6]=0;
    Jf[2][7]=0;
    Jf[2][8]=0;
    Jf[2][9]=-dr_dx1[0];
    Jf[2][10]=-dr_dt1[0];
    
    
    Jf[3][0]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[1]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[1]);
    Jf[3][1]=(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4]));
    Jf[3][2]=bv0[1]-v_1[1]-v[4];
    Jf[3][3]=-(((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*da1_dv1[0]);
    Jf[3][4]=-(((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*da1_dv1[1])-g1;
    Jf[3][5]=-(((df1_da1*r_1[1])+(dg1_da1*(v_1[1]+v[4])))*da1_dv1[2]);
    Jf[3][6]=0;
    Jf[3][7]=0;
    Jf[3][8]=0;
    Jf[3][9]=-dr_dx1[1];
    Jf[3][10]=-dr_dt1[1];
    
    
    Jf[4][0]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[2]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[2]);
    Jf[4][1]=(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5]));
    Jf[4][2]=bv0[2]-v_1[2]-v[5];
    Jf[4][3]=-(((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*da1_dv1[0]);
    Jf[4][4]=-(((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*da1_dv1[1]);
    Jf[4][5]=-(((df1_da1*r_1[2])+(dg1_da1*(v_1[2]+v[5])))*da1_dv1[2])-g1;
    Jf[4][6]=0;
    Jf[4][7]=0;
    Jf[4][8]=0;
    Jf[4][9]=-dr_dx1[2];
    Jf[4][10]=-dr_dt1[2];
    
    
    Jf[5][0]=-(1/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*cos(v[0]/sqrt(ba))*((v[6]*br0[0])+(v[7]*br0[1])+(v[8]*br0[2]))) - (sqrt(ba/mu)*sin(v[0]/sqrt(ba))*((v[6]*bv0[0])+(v[7]*bv0[1])+(v[8]*bv0[2])));
    Jf[5][1]=(1/ra1*cos(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + (sqrt(a/mu)*sin(v[1]/sqrt(a))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))) + (((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5]))));
    Jf[5][2]=0;
    Jf[5][3]=( (((((1/(2*sqrt(a)*ra1)*sin(v[1]/sqrt(a)))-(v[1]/(2*a*ra1)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))) + ((1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5])))))*(2*a*a*(v_1[0]+v[3])/mu)) + (((a-(a*cos(v[1]/sqrt(a))))-ra2)*v[6]/sqrt(mu)) ) + (1/sqrt(mu)*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[0])*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5]))));
    Jf[5][4]=( (((((1/(2*sqrt(a)*ra1)*sin(v[1]/sqrt(a)))-(v[1]/(2*a*ra1)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))) + ((1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5])))))*(2*a*a*(v_1[1]+v[4])/mu)) + (((a-(a*cos(v[1]/sqrt(a))))-ra2)*v[7]/sqrt(mu)) ) + (1/sqrt(mu)*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[1])*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5]))));
    Jf[5][5]=( (((((1/(2*sqrt(a)*ra1)*sin(v[1]/sqrt(a)))-(v[1]/(2*a*ra1)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((1-(cos(v[1]/sqrt(a)))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))) + ((1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))/sqrt(mu)*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5])))))*(2*a*a*(v_1[2]+v[5])/mu)) + (((a-(a*cos(v[1]/sqrt(a))))-ra2)*v[8]/sqrt(mu)) ) + (1/sqrt(mu)*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[2])*((v[6]*(bv0[0]-v_1[0]-v[3]))+(v[7]*(bv0[1]-v_1[1]-v[4]))+(v[8]*(bv0[2]-v_1[2]-v[5]))));
    
    Jf[5][6]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[0]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[0]) + (sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])) + (ra2/sqrt(mu)*(bv0[0]-v_1[0]-v[3]));
    Jf[5][7]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[1]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[1]) + (sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])) + (ra2/sqrt(mu)*(bv0[1]-v_1[1]-v[4]));
    Jf[5][8]=-(sqrt(ba)/sqrt((br0[0]*br0[0])+(br0[1]*br0[1])+(br0[2]*br0[2]))*sin(v[0]/sqrt(ba))*br0[2]) - (ba/sqrt(mu)*(1-(cos(v[0]/sqrt(ba))))*bv0[2]) + (sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])) + (ra2/sqrt(mu)*(bv0[2]-v_1[2]-v[5]));
    
    Jf[5][9]=-((ddr_dxx1[0]*v[6])+(ddr_dxx1[1]*v[7])+(ddr_dxx1[2]*v[8])) + (dra_dx1/sqrt(mu)*((v[6]*(br0[0]-v_1[0]-v[3]))+(v[7]*(br0[1]-v_1[1]-v[4]))+(v[8]*(br0[2]-v_1[2]-v[5])))) - (ra2/sqrt(mu)*((v[6]*dv1_dx1[0])+(v[7]*dv1_dx1[1])+(v[8]*dv1_dx1[2])));
    Jf[5][10]=sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*v0[0])+(v[7]*v0[1])+(v[8]*v0[2]));
    
    
    Jf[6][0]=0;
    Jf[6][1]=((1/ra2)*(((((1/ra2)*((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) - (1/ra1*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) - (sqrt(a/mu)*sin(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) +  ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[0]+v[3])*( 1 - cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (ra1/2*(v[1]*sin(v[1]/sqrt(a))/(sqrt(a)*sqrt(a)*sqrt(a)))) )) + (sqrt(a/mu)*(sin(v[1]/sqrt(a)))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) ) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v[6]);
    Jf[6][2]=-v[6];
    Jf[6][3]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[0]+v[3])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[0]))) - ( (v_1[0]+v[3])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[6])  )*(2*a*a/mu*(v_1[0]+v[3]))) + (1 - ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[0])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[6]))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*deta3_da1))) - ( 2*a*a/mu*(v_1[0]+v[3])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[6])) ) - ( 2*a*a/mu*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) )) );
    Jf[6][4]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[0]+v[3])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[0]))) - ( (v_1[0]+v[3])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[6])  )*(2*a*a/mu*(v_1[1]+v[4]))) + (- ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[1])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[7]))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[7])) ) ) );
    Jf[6][5]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[0]+v[3])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[0]))) - ( (v_1[0]+v[3])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[6])  )*(2*a*a/mu*(v_1[2]+v[5]))) + (- ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[2])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[8]))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[8])) ) ) );
    
    Jf[6][6]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[0]+v[3]))) ) - ((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))));
    Jf[6][7]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[1]+v[4]))) );
    Jf[6][8]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])))*((2*a*a/mu*(v_1[0]+v[3])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[0]))) - ( 2*a*a/mu*(v_1[0]+v[3])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[2]+v[5]))) );
    
    Jf[6][9]= - (1/(ra2*ra2)*dra_dx1*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))*Jf[1][3]) + (1/ra2*((ddr_dxx1[0]*v[6])+(ddr_dxx1[1]*v[7])+(ddr_dxx1[2]*v[8]))*Jf[1][3]) + (1/ra2*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))*ddeta3_dv1x1[0]) - ((ddr_dv1x1[0][0]*v[6])+(ddr_dv1x1[1][0]*v[7])+(ddr_dv1x1[2][0]*v[8]));
    
    Jf[6][10]=v[6];
    
    Jf[7][0]=0;
    Jf[7][1]=((1/ra2)*(((((1/ra2)*((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) - (1/ra1*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) - (sqrt(a/mu)*sin(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[1]+v[4])*( 1 - cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (ra1/2*(v[1]*sin(v[1]/sqrt(a))/(sqrt(a)*sqrt(a)*sqrt(a)))) )) + (sqrt(a/mu)*(sin(v[1]/sqrt(a)))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) ) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v[7]);
    Jf[7][2]=-v[7];
    Jf[7][3]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[1]+v[4])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[1]))) - ( (v_1[1]+v[4])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[7])  )*(2*a*a/mu*(v_1[0]+v[3]))) + (- ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[0])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[6]))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[6])) ) ) );
    Jf[7][4]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[1]+v[4])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[1]))) - ( (v_1[1]+v[4])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[7])  )*(2*a*a/mu*(v_1[1]+v[4]))) + (1 - ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[1])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[7]))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*deta3_da1))) - ( 2*a*a/mu*(v_1[1]+v[4])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[7])) ) - ( 2*a*a/mu*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) )) );
    Jf[7][5]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[1]+v[4])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[1]))) - ( (v_1[1]+v[4])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[7])  )*(2*a*a/mu*(v_1[2]+v[5]))) + (- ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[2])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[8]))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[8])) ) ) );
    
    Jf[7][6]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[0]+v[3]))) );
    Jf[7][7]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[1]+v[4]))) ) - ((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))));
    Jf[7][8]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])))*((2*a*a/mu*(v_1[1]+v[4])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[1]))) - ( 2*a*a/mu*(v_1[1]+v[4])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[2]+v[5]))) );
    
    Jf[7][9]= - (1/(ra2*ra2)*dra_dx1*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))*Jf[1][4]) + (1/ra2*((ddr_dxx1[0]*v[6])+(ddr_dxx1[1]*v[7])+(ddr_dxx1[2]*v[8]))*Jf[1][4]) + (1/ra2*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))*ddeta3_dv1x1[1]) - ((ddr_dv1x1[0][1]*v[6])+(ddr_dv1x1[1][1]*v[7])+(ddr_dv1x1[2][1]*v[8]));
    
    Jf[7][10]=v[7];
    
    Jf[8][0]=0;
    Jf[8][1]=((1/ra2)*(((((1/ra2)*((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) - (1/ra1*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*((((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*cos(v[1]/sqrt(a)))-((ra1-a)/sqrt(a)*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) - (sqrt(a/mu)*sin(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[2]+v[5])*( 1 - cos(v[1]/sqrt(a)) - (v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))) + (ra1/2*(v[1]*sin(v[1]/sqrt(a))/(sqrt(a)*sqrt(a)*sqrt(a)))) )) + (sqrt(a/mu)*(sin(v[1]/sqrt(a)))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*((1/(2*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]/(2*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) ) + (a/sqrt(mu)*(1-cos(v[1]/sqrt(a)))*v[8]);
    Jf[8][2]=-v[8];
    Jf[8][3]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[2]+v[5])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[2]))) - ( (v_1[2]+v[5])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[8])  )*(2*a*a/mu*(v_1[0]+v[3]))) + (- ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[0])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[6]))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[6])) ) ) );
    Jf[8][4]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[2]+v[5])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[2]))) - ( (v_1[2]+v[5])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[8])  )*(2*a*a/mu*(v_1[1]+v[4]))) + (- ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[1])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[7]))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[7])) ) ) );
    Jf[8][5]=( (  ( ((1/ra2)*(((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(sqrt(a)/ra1*sin(v[1]/sqrt(a)))) + (1/(2*ra1)*(-(sin(v[1]/sqrt(a))/sqrt(a))+(v[1]*cos(v[1]/sqrt(a))/a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) + ((((1/ra2)*(1+(((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))-(v[1]*((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/(2*a*sqrt(mu))*cos(v[1]/sqrt(a)))-cos(v[1]/sqrt(a))+((ra1-a)*v[1]/(2*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a))))*(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a)))))) + (-(1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))+(v[1]/(2*sqrt(a*mu))*sin(v[1]/sqrt(a)))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) *((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(((v_1[2]+v[5])*((( (3*v[1]/(4*a)*cos(v[1]/sqrt(a))) + ((-(3/(4*sqrt(a)))+(v[1]*v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a))) + (((r_1[0]*(v_1[0]+v[3]))+(r_1[1]*(v_1[1]+v[4]))+(r_1[2]*(v_1[2]+v[5])))/sqrt(mu)*(-(v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))+(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))) + (ra1/2*((v[1]/(2*a*a)*cos(v[1]/sqrt(a)))-((1+(v[1]*v[1]/a))*sin(v[1]/sqrt(a))/(2*sqrt(a)*sqrt(a)*sqrt(a))))) ) *(2*a*a/mu)) + (deta3_da1*(4*a/mu))) ) + (((1/sqrt(mu)*(1-cos(v[1]/sqrt(a))))-(v[1]/(2*sqrt(a*mu))*(sin(v[1]/sqrt(a)))))*r_1[2]))) - ( (v_1[2]+v[5])*(((-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(4*a/mu)) + (( (1/ra1*((v[1]/(4*sqrt(a)*sqrt(a)*sqrt(a))*sin(v[1]/sqrt(a)))-(v[1]*v[1]/(4*a*a)*cos(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/(4*sqrt(mu))*((((3/sqrt(a))-(v[1]*v[1]/(sqrt(a)*sqrt(a)*sqrt(a))))*sin(v[1]/sqrt(a)))-(3*v[1]/a*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*(2*a*a/mu))) ) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[8])  )*(2*a*a/mu*(v_1[2]+v[5]))) + (1 - ((1/pow(ra2,2))*(sqrt(a/mu)*sin(v[1]/sqrt(a))*r_1[2])*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1 ) + (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*v[8]))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) + ((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5])))))*((2*a*a/mu*deta3_da1))) - ( 2*a*a/mu*(v_1[2]+v[5])*(- (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*v[8])) ) - ( 2*a*a/mu*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*((v[6]*r_1[0])+(v[7]*r_1[1])+(v[8]*r_1[2]))) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*((v[6]*(v_1[0]+v[3]))+(v[7]*(v_1[1]+v[4]))+(v[8]*(v_1[2]+v[5]))))) )) );
    
    Jf[8][6]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[0]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[0]+v[3])))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[0]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[0]+v[3]))) );
    Jf[8][7]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[1]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[1]+v[4])))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[1]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[1]+v[4]))) );
    Jf[8][8]=((1/ra2)*(-(sqrt(a)/ra1*sin(v[1]/sqrt(a))*r_1[2]) - (a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*(v_1[2]+v[5])))*((2*a*a/mu*(v_1[2]+v[5])*deta3_da1)+(a/sqrt(mu)*(1-(cos(v[1]/sqrt(a))))*r_1[2]))) - ( 2*a*a/mu*(v_1[2]+v[5])*(-(1/ra1*(1-cos(v[1]/sqrt(a))-(v[1]/(2*sqrt(a))*sin(v[1]/sqrt(a))))*r_1[2]) - (1/sqrt(mu)*(v[1]-(3*sqrt(a)/2*sin(v[1]/sqrt(a)))+(v[1]/2*cos(v[1]/sqrt(a))))*(v_1[2]+v[5]))) ) - ((v[2]-v[10])-(a/sqrt(mu)*(v[1]-(sqrt(a)*sin(v[1]/sqrt(a))))));   
      
    
    Jf[8][9]= - (1/(ra2*ra2)*dra_dx1*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))*Jf[1][5]) + (1/ra2*((ddr_dxx1[0]*v[6])+(ddr_dxx1[1]*v[7])+(ddr_dxx1[2]*v[8]))*Jf[1][5]) + (1/ra2*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))*ddeta3_dv1x1[2]) - ((ddr_dv1x1[0][2]*v[6])+(ddr_dv1x1[1][2]*v[7])+(ddr_dv1x1[2][2]*v[8]));
    
    Jf[8][10]=v[8];
    
    
    Jf[9][0]=0;
    Jf[9][1]=((-(1/(ra2*ra2)*dra_dx*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1))+(1/ra2*((sqrt(mu)/ra1*ddeta3_dx1x)+ddeta3_dt1x)))*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))) + (1/ra2*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((ddr_dxx[0]*v[6])+(ddr_dxx[1]*v[7])+(ddr_dxx[2]*v[8]))) - (sqrt(mu)/ra1*((ddr_dx1x[0]*v[6])+(ddr_dx1x[1]*v[7])+(ddr_dx1x[2]*v[8]))) - ((ddr_dt1x[0]*v[6])+(ddr_dt1x[1]*v[7])+(ddr_dt1x[2]*v[8]));
    Jf[9][2]=-(sqrt(mu)/ra1*((ddr_dx1t[0]*v[6])+(ddr_dx1t[1]*v[7])+(ddr_dx1t[2]*v[8])));
    
    Jf[9][3]= ((-(1/(ra2*ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*dra_dv1[0])+(1/ra2*((sqrt(mu)/ra1*ddeta3_dx1v1[0])+ddeta3_dt1v1[0])))*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))) + (1/ra2*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((ddr_dxv1[0][0]*v[6])+(ddr_dxv1[1][0]*v[7])+(ddr_dxv1[2][0]*v[8]))) - (sqrt(mu)/ra1*((ddr_dx1v1[0][0]*v[6])+(ddr_dx1v1[1][0]*v[7])+(ddr_dx1v1[2][0]*v[8]))) - ((ddr_dt1v1[0][0]*v[6])+(ddr_dt1v1[1][0]*v[7])+(ddr_dt1v1[2][0]*v[8]));
    Jf[9][4]= ((-(1/(ra2*ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*dra_dv1[1])+(1/ra2*((sqrt(mu)/ra1*ddeta3_dx1v1[1])+ddeta3_dt1v1[1])))*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))) + (1/ra2*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((ddr_dxv1[0][1]*v[6])+(ddr_dxv1[1][1]*v[7])+(ddr_dxv1[2][1]*v[8]))) - (sqrt(mu)/ra1*((ddr_dx1v1[0][1]*v[6])+(ddr_dx1v1[1][1]*v[7])+(ddr_dx1v1[2][1]*v[8]))) - ((ddr_dt1v1[0][1]*v[6])+(ddr_dt1v1[1][1]*v[7])+(ddr_dt1v1[2][1]*v[8]));;
    Jf[9][5]= ((-(1/(ra2*ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*dra_dv1[2])+(1/ra2*((sqrt(mu)/ra1*ddeta3_dx1v1[2])+ddeta3_dt1v1[2])))*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))) + (1/ra2*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((ddr_dxv1[0][2]*v[6])+(ddr_dxv1[1][2]*v[7])+(ddr_dxv1[2][2]*v[8]))) - (sqrt(mu)/ra1*((ddr_dx1v1[0][2]*v[6])+(ddr_dx1v1[1][2]*v[7])+(ddr_dx1v1[2][2]*v[8]))) - ((ddr_dt1v1[0][2]*v[6])+(ddr_dt1v1[1][2]*v[7])+(ddr_dt1v1[2][2]*v[8]));;
    
    Jf[9][6]=((1/ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*dr_dx[0])-((sqrt(mu)/ra1*dr_dx1[0])+dr_dt1[0]);
    Jf[9][7]=((1/ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*dr_dx[1])-((sqrt(mu)/ra1*dr_dx1[1])+dr_dt1[1]);
    Jf[9][8]=((1/ra2)*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*dr_dx[2])-((sqrt(mu)/ra1*dr_dx1[2])+dr_dt1[2]);
    
    Jf[9][9]= ((-(1/(ra2*ra2)*dra_dx1*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1))+(1/ra2*(-(sqrt(mu)/(ra1*ra1)*dr1_dx1*deta3_dx1)+(sqrt(mu)/ra1*ddeta3_dx1x1)+ddeta3_dt1x1)))*((dr_dx[0]*v[6])+(dr_dx[1]*v[7])+(dr_dx[2]*v[8]))) + (1/ra2*((sqrt(mu)/ra1*deta3_dx1)+deta3_dt1)*((ddr_dxx1[0]*v[6])+(ddr_dxx1[1]*v[7])+(ddr_dxx1[2]*v[8]))) + (sqrt(mu)/(ra1*ra1)*dr1_dx1*((dr_dx1[0]*v[6])+(dr_dx1[1]*v[7])+(dr_dx1[2]*v[8]))) - (sqrt(mu)/ra1*((ddr_dx1x1[0]*v[6])+(ddr_dx1x1[1]*v[7])+(ddr_dx1x1[2]*v[8]))) - ((ddr_dt1x1[0]*v[6])+(ddr_dt1x1[1]*v[7])+(ddr_dt1x1[2]*v[8]));
    
    Jf[9][10]=-(sqrt(mu)/ra1)*((dv1_dx1[0]*v[6])+(dv1_dx1[1]*v[7])+(dv1_dx1[2]*v[8]));
    
    
    Jf[10][0]=0;
    Jf[10][1]=0;
    Jf[10][2]=0;
    Jf[10][3]=0;
    Jf[10][4]=0;
    Jf[10][5]=0;
    Jf[10][6]=0;
    Jf[10][7]=0;
    Jf[10][8]=0;
    Jf[10][9]=(a0+(sqrt(a0/mu)*((r0[0]*v0[0])+(r0[1]*v0[1])+(r0[2]*v0[2]))*sin(v[9]/sqrt(a0)))+((sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))-a0)*cos(v[9]/sqrt(ba))));
    Jf[10][10]=-sqrt(mu);
    
/*    for(i=0; i<11; i++){
      for(j=0; j<11; j++){
         cout<<Jf[2][j]<<"  ";
      }
      cout<<endl;
  }
*/   

// Cálculo inversa

    for(i=0; i<11; i++){
      for(j=0; j<11; j++){
          if(i==j){
             id[i][j]=1;
          }else{
             id[i][j]=0;
          }
      }
    }
    
    for(i=0; i<11; i++){
             
          if(fabs(Jf[i][i])<pow(10,-14)){
             for(k=0; k<11; k++){
                if(k!=i && fabs(Jf[k][i])>pow(10,-14)){
                   for(j=0; j<11; j++){
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
       for(j=0; j<11; j++){
          Jf[i][j]=Jf[i][j]/s;
             id[i][j]=id[i][j]/s;
          b[j]=Jf[i][j];
             c[j]=id[i][j];  
       }
       
       for(k=0; k<11; k++){
          if(k!=i){   
             r=Jf[k][i];
             for(j=0; j<11; j++){
                Jf[k][j]=Jf[k][j]-(r*b[j]);
                   id[k][j]=id[k][j]-(r*c[j]);
             }
          }
       }           
    }
    
/*   for(i=0; i<11; i++){
      for(j=0; j<11; j++){
         cout<<id[10][j]<<"  ";
      }
      cout<<endl;
   }
*/  

// Operaciones Newton

    for(i=0; i<11; i++){
       prod[i]=0;
    }
    
    for(i=0; i<11; i++){
       for(j=0; j<11; j++){
          prod[i]=prod[i]+id[i][j]*f[j];
       }
    }
    
    for(i=0; i<11; i++){
       v[i]=v[i]-prod[i];
    }
    
    l=l+1;
    cout<<l<<endl;

    for(i=0; i<11; i++){
        cout<<"v"<<i<<"= "<<v[i]<<"   ";
    }cout<<endl;
 }

// Presentacion del resultado

    cout<<"---------------------"<<endl;
    cout<<"La solucion optima es"<<endl;
    cout<<"---------------------"<<endl;

    for(i=0; i<11; i++){
       cout<<"v"<<i<<"= "<<v[i]<<"   "<<endl;
    }
    
    cout<<"---------------------------"<<endl;   
    cout<<"t = "<<v[2]*(5.0227*pow(10,6))/86400 <<" dias"<<endl;
    cout<<"---------------------------"<<endl;   
    cout<<"t1 = "<<v[10]*(5.0227*pow(10,6))/86400 <<" dias"<<endl;     
    cout<<"--------------------------"<<endl;
    cout<<"Av = "<< sqrt((v[3]*v[3])+(v[4]*v[4])+(v[5]*v[5]))*149597870.7/(5.0227*pow(10,6)) <<" km/s"<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"vx = "<<v[3]*149597870.7/(5.0227*pow(10,6))<<" km/s"<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"vy = "<<v[4]*149597870.7/(5.0227*pow(10,6))<<" km/s"<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"vz = "<<v[5]*149597870.7/(5.0227*pow(10,6))<<" km/s"<<endl; 
    cout<<"---------------------------"<<endl;

// Elementos orbita de interseccion
/*
    double h[3];
    double e[3];
    double p, bo;

    h[0]=(r0[1]*(v0[2]+v[5]))-(r0[2]*(v0[1]+v[4]));
    h[1]=-((r0[0]*(v0[2]+v[5]))-(r0[2]*(v0[0]+v[3])));
    h[2]=(r0[0]*(v0[1]+v[4]))-(r0[1]*(v0[0]+v[3]));

    p=((h[0]*h[0])+(h[1]*h[1])+(h[2]*h[2]))/mu;
    a=(-1)/(((((v0[0]+v[3])*(v0[0]+v[3]))+((v0[1]+v[4])*(v0[1]+v[4]))+((v0[2]+v[5])*(v0[2]+v[5])))/mu)-(2/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))));
    bo=sqrt(p)*sqrt(a);

    cout<<"a= "<<a<<endl;
    cout<<"b= "<<bo<<endl;

    e[0]=(1/mu)*( ((((((v0[0]+v[3])*(v0[0]+v[3]))+((v0[1]+v[4])*(v0[1]+v[4]))+((v0[2]+v[5])*(v0[2]+v[5]))))-(mu/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*r0[0]) - (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))*(v0[0]+v[3])));
    e[1]=(1/mu)*( ((((((v0[0]+v[3])*(v0[0]+v[3]))+((v0[1]+v[4])*(v0[1]+v[4]))+((v0[2]+v[5])*(v0[2]+v[5]))))-(mu/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*r0[1]) - (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))*(v0[1]+v[4])));
    e[2]=(1/mu)*( ((((((v0[0]+v[3])*(v0[0]+v[3]))+((v0[1]+v[4])*(v0[1]+v[4]))+((v0[2]+v[5])*(v0[2]+v[5]))))-(mu/sqrt((r0[0]*r0[0])+(r0[1]*r0[1])+(r0[2]*r0[2]))))*r0[2]) - (((r0[0]*(v0[0]+v[3]))+(r0[1]*(v0[1]+v[4]))+(r0[2]*(v0[2]+v[5])))*(v0[2]+v[5])));

    for(i=0; i<3; i++){
       cout<<"h"<<i<<"= "<<h[i]<<endl;
    }

    for(i=0; i<3; i++){
       cout<<"e"<<i<<"= "<<e[i]<<endl;
    }

    for(i=0; i<3; i++){
       cout<<"d"<<i<<"= "<<-e[i]*a<<endl;
    }
*/
    system("PAUSE");
    return EXIT_SUCCESS;
}
