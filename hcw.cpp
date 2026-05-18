#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

int main(int argc, char *argv[])
{
    double x0, y0, z0, t, t1, vx0, vy0, vz0, w, mu, orb; 
    double rx, ry, rz, vx, vy, vz;
    double dvx, dvy, dvz;
    double ay, az, X[9], X1[9], rot2[3][3], rot3[3][3];
    int i,j;

// Calculos en Hill de r, v y Av del interceptor
    
    t=2400;
    t1=300;
    mu=398600.5;
    orb=590;    
    w=sqrt(mu/((6371+orb)*(6371+orb)*(6371+orb)));
    
    x0=0;
    y0=0;   
    z0=0;
    vx0=-0.1;
    vy0=-0.04;
    vz0=-0.02;
    
// r
    
    rx=(vx0/w*sin(w*t))-(((3*x0)+(2*vy0/w))*cos(w*t))+((4*x0)+(2*vy0/w));
       cout<<"x("<<t<<") = "<<rx<<endl;
    ry=(((6*x0)+(4*vy0/w))*sin(w*t))+((2*vx0)/w*cos(w*t))-(((6*w*x0)+(3*vy0))*t)+((y0)-(2*vx0/w));
       cout<<"y("<<t<<") = "<<ry<<endl;
    rz=(z0*cos(w*t))+(vz0/w*sin(w*t));
       cout<<"z("<<t<<") = "<<rz<<endl;   
       
// v
       
    vx=(vx0*cos(w*t))+(((3*x0*w)+(2*vy0))*sin(w*t));
       cout<<"vx("<<t<<") = "<<vx<<endl;
    vy=(((6*x0*w)+(4*vy0))*cos(w*t))-(2*vx0*sin(w*t))-((6*w*x0)+(3*vy0));
       cout<<"vy("<<t<<") = "<<vy<<endl;
    vz=-(z0*w*sin(w*t))+(vz0*cos(w*t));
       cout<<"vz("<<t<<") = "<<vz<<endl;    
       
// Av   
         
    dvy=((((6*rx*((w*t1)-sin(w*t1)))-ry)*w*sin(w*t1))-(2*w*rx*(4-3*cos(w*t1))*(1-cos(w*t1))))/(((4*sin(w*t1)-(3*w*t1))*sin(w*t1))+(4*(1-cos(w*t1))*(1-cos(w*t1))));
       cout<<"dvy0("<<t1<<") = "<<dvy<<endl;
    dvx=-((w*rx*(4-(3*cos(w*t1))))+((2*dvy)*(1-cos(w*t1))))/sin(w*t1);
       cout<<"dvx0("<<t1<<") = "<<dvx<<endl;
    dvz=-(rz*w*cos(w*t1)/sin(w*t1));
       cout<<"dvz0("<<t1<<") = "<<dvz<<endl;  

// Conversion de hill (RSW) a geocentricas (IJK) de r, v y Av del interceptor ([R|S|W]=Id) (mirar Av)
    
    X[0]=6961+rx;
    X[1]=0;
    X[2]=0;
    X[3]=0+vx;
    X[4]=7.57+vy;
    X[5]=0+vz;
    X[6]=0+dvx;
    X[7]=7.57+dvy;
    X[8]=0+dvz;
       
    ay=atan2(ry,6961);
    az=atan2(rz,6961);
    
    rot2[0][0]=cos(-ay);
    rot2[0][1]=0;
    rot2[0][2]=-sin(-ay);
    rot2[1][0]=0;
    rot2[1][1]=1;
    rot2[1][2]=0;
    rot2[2][0]=sin(-ay);
    rot2[2][1]=0;
    rot2[2][2]=cos(-ay);
    
    rot3[0][0]=cos(az);
    rot3[0][1]=sin(az);
    rot3[0][2]=0;
    rot3[1][0]=-sin(az);
    rot3[1][1]=cos(az);
    rot3[1][2]=0;
    rot3[2][0]=0;
    rot3[2][1]=0;
    rot3[2][2]=1;
    
    if(fabs(ry)>0){

// r
                   
       for(i=0; i<3; i++){
          X1[i]=0;      
          for(j=0; j<3; j++){
             X1[i]=X1[i]+rot2[i][j]*X[j];
          }          
       }
       for(i=0; i<3; i++){
          X[i]=X1[i];
       }
   
// v   
       
       for(i=0; i<3; i++){
          X1[i+3]=0;      
          for(j=0; j<3; j++){
             X1[i+3]=X1[i+3]+rot2[i][j]*X[j+3];
          }          
       }
       for(i=0; i<3; i++){
          X[i+3]=X1[i+3];
       }
   
// Av   
       
       for(i=0; i<3; i++){
          X1[i+6]=0;      
          for(j=0; j<3; j++){
             X1[i+6]=X1[i+6]+rot2[i][j]*X[j+6];
          }          
       }
       for(i=0; i<3; i++){
          X[i+6]=X1[i+6];
       }
    }
    
    if(fabs(rz)>0){
                   
// r
                   
       for(i=0; i<3; i++){
          X1[i]=0;      
          for(j=0; j<3; j++){
             X1[i]=X1[i]+rot3[i][j]*X[j];
          }          
       }
       for(i=0; i<3; i++){
          X[i]=X1[i];
       }
       
// v
       
       for(i=0; i<3; i++){
          X1[i+3]=0;      
          for(j=0; j<3; j++){
             X1[i+3]=X1[i+3]+rot3[i][j]*X[j+3];
          }          
       }
       for(i=0; i<3; i++){
          X[i+3]=X1[i+3];
       }
       
// Av
       
       for(i=0; i<3; i++){
          X1[i+6]=0;      
          for(j=0; j<3; j++){
             X1[i+6]=X1[i+6]+rot3[i][j]*X[j+6];
          }          
       }
       for(i=0; i<3; i++){
          X[i+6]=X1[i+6];
       }
    }
    
    for(i=0; i<3; i++){
       cout<<"r"<<i<<"= "<<X[i]<<endl;
    }
    
    for(i=0; i<3; i++){
       cout<<"v"<<i+3<<"= "<<X[i+3]<<endl;
    }
    
    for(i=0; i<3; i++){
       cout<<"dv"<<i+6<<"= "<<X[i+6]<<endl;
    }
     
    system("PAUSE");
    return EXIT_SUCCESS;
}
