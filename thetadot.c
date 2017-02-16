/*
=========================================================================
thetadot.c
=========================================================================
Here are implemented the routines to calculated the derivative of
the expansion.


*************************************************************
* Copyright (C) 2017					    *
*							    *
* Edison Montoya				       	    *
* eamonto@gmail.com		                            *
*							    *
* Up to date: Feb 16 2017				    *
* 							    *
*************************************************************
*/

#include <header.h>


// DERIVATIVE OF THE CLASSICAL EXPANSION
double class_thetadot(void)
{
  double val=0.0;

  //Raychaudhuri Equation
  val = -0.5*expansion*expansion -shear -16.0*M_PI*G*density; 

  return val;
}


// DERIVATIVE OF THE EXPANSION (EDWARD)
double Edward_thetadot(void)
{ 
  double val;
  double dev_mu1c1,dev_mu2c2,dev_mu3c3;
  
  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  dev_mu1c1 = -mu1*H1*c1 + mu1*sc1;
  dev_mu2c2 = -mu2*H2*c2 + mu2*sc2;
  dev_mu3c3 = -mu3*H3*c3 + mu3*sc3;
  
  val = 1.0/(2.0*gamma*lambda)
    *( (dev_mu1c1+dev_mu2c2)*cos(mu1*c1+mu2*c2)
     + (dev_mu1c1+dev_mu3c3)*cos(mu1*c1+mu3*c3)
     + (dev_mu3c3+dev_mu2c2)*cos(mu3*c3+mu2*c2)
     + Delta*L0/(mu1*p1)*( dev_mu1c1*sin(mu1*c1) + (1.5*sp1/p1-0.5*sp2/p2-0.5*sp3/p3)*cos(mu1*c1) )
     + Delta*L0/(mu2*p2)*( dev_mu2c2*sin(mu2*c2) + (1.5*sp2/p2-0.5*sp1/p1-0.5*sp3/p3)*cos(mu2*c2) )
     + Delta*L0/(mu3*p3)*( dev_mu3c3*sin(mu3*c3) + (1.5*sp3/p3-0.5*sp1/p1-0.5*sp2/p2)*cos(mu3*c3) )
       );
  
  return val;
}


// DERIVATIVE OF THE EXPANSION (ASIEH, NOT FRW)
double Asieh_notFRW_thetadot(void)
{
  double s1,s2,s3;
  double A_V,h_V,dev_Av,dev_hv;

  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  mu1c1 = mu1*c1;
  mu2c2 = mu2*c2;
  mu3c3 = mu3*c3;
  
  A_V = func_A(vol);
  h_V = func_h(vol);

  dev_Av = dev_A(vol,0);
  dev_hv = dev_h(vol,0);

  s1 = vol*vol*vol*expansion*pow(h_V/V_crit,5)/(V_crit*gamma*lambda)
    *cos(mu1c1)*( vol*h_V*dev_Av + 6.0*vol*A_V*dev_hv + 3.0*A_V*h_V)
    *(sin(mu2c2) + sin(mu3c3))
    +vol*vol*vol*A_V*pow(h_V/V_crit,6)/(gamma*lambda)
    *( cos(mu1c1)*cos(mu2c2)*(-mu2*H2*c2 + mu2*sc2) 
      +cos(mu1c1)*cos(mu3c3)*(-mu3*H3*c3 + mu3*sc3) 
      -sin(mu1c1)*(-mu1*H1*c1 + mu1*sc1)*( sin(mu2c2) + sin(mu3c3) )
       )
    -2.0*L0*pow(h_V/V_crit,3)/(V_crit*gamma)*p2*p3*sqrt(p2*p3/p1)*cos(mu1c1)
    *( vol*expansion*dev_Av*h_V + 4.0*vol*expansion*A_V*dev_hv 
       +0.5*A_V*h_V*( 3.0*expansion-4.0*sp1/p1) 
       )
    +2.0*L0*A_V*pow(h_V/V_crit,4)/(gamma)*p2*p3*sqrt(p2*p3/p1)*sin(mu1c1)*
    (-mu1*H1*c1 + mu1*sc1);

  s2 = vol*vol*vol*expansion*pow(h_V/V_crit,5)/(V_crit*gamma*lambda)
    *cos(mu2c2)*( vol*h_V*dev_Av + 6.0*vol*A_V*dev_hv + 3.0*A_V*h_V)
    *(sin(mu1c1) + sin(mu3c3))
    +vol*vol*vol*A_V*pow(h_V/V_crit,6)/(gamma*lambda)
    *( cos(mu2c2)*cos(mu1c1)*(-mu1*H1*c1 + mu1*sc1) 
      +cos(mu2c2)*cos(mu3c3)*(-mu3*H3*c3 + mu3*sc3) 
      -sin(mu2c2)*(-mu2*H2*c2 + mu2*sc2)*( sin(mu1c1) + sin(mu3c3) )
       )
    -2.0*L0*pow(h_V/V_crit,3)/(V_crit*gamma)*p1*p3*sqrt(p1*p3/p2)*cos(mu2c2)
    *( vol*expansion*dev_Av*h_V + 4.0*vol*expansion*A_V*dev_hv 
       +0.5*A_V*h_V*( 3.0*expansion-4.0*sp2/p2) 
       )
    +2.0*L0*A_V*pow(h_V/V_crit,4)/(gamma)*p1*p3*sqrt(p1*p3/p2)*sin(mu2c2)*
    (-mu2*H2*c2 + mu2*sc2);

  s3 = vol*vol*vol*expansion*pow(h_V/V_crit,5)/(V_crit*gamma*lambda)
    *cos(mu3c3)*( vol*h_V*dev_Av + 6.0*vol*A_V*dev_hv + 3.0*A_V*h_V)
    *(sin(mu2c2) + sin(mu1c1))
    +vol*vol*vol*A_V*pow(h_V/V_crit,6)/(gamma*lambda)
    *( cos(mu3c3)*cos(mu2c2)*(-mu2*H2*c2 + mu2*sc2) 
      +cos(mu3c3)*cos(mu1c1)*(-mu1*H1*c1 + mu1*sc1) 
      -sin(mu3c3)*(-mu3*H3*c3 + mu3*sc3)*( sin(mu2c2) + sin(mu1c1) )
       )
    -2.0*L0*pow(h_V/V_crit,3)/(V_crit*gamma)*p2*p1*sqrt(p2*p1/p3)*cos(mu3c3)
    *( vol*expansion*dev_Av*h_V + 4.0*vol*expansion*A_V*dev_hv 
       +0.5*A_V*h_V*( 3.0*expansion-4.0*sp3/p3) 
       )
    +2.0*L0*A_V*pow(h_V/V_crit,4)/(gamma)*p2*p1*sqrt(p2*p1/p3)*sin(mu3c3)*
    (-mu3*H3*c3 + mu3*sc3);

  return 0.5*(s1+s2+s3);
}


// DERIVATIVE OF THE EXPANSION (ASIEH,FRW)
double Asieh_FRW_thetadot(void)
{
  double r1,r2,r3;
  double A_V,dev_Av;
  
  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  mu1c1 = mu1*c1;
  mu2c2 = mu2*c2;
  mu3c3 = mu3*c3;

  A_V = func_A(vol);
  dev_Av = dev_A(vol,0);

  r1 = 1.0/(gamma*lambda)*( vol*expansion*dev_Av*cos(mu1c1) 
			   -A_V*(-mu1*H1*c1 + mu1*sc1)*sin(mu1c1) )
    *( sin(mu2c2) + sin(mu3c3) -2.0*L0*lambda*sqrt(p2*p3/p1)/p1 )
    + 1.0/(gamma*lambda)*A_V*cos(mu1c1)
    *( (-mu2*H2*c2 + mu2*sc2)*cos(mu2c2) + (-mu3*H3*c3 + mu3*sc3)*cos(mu3c3) 
       + 2.0*L0*lambda*sqrt(p2*p3/p1)*(sp1/p1 - expansion) );

  r2 = 1.0/(gamma*lambda)*( vol*expansion*dev_Av*cos(mu2c2) 
			   -A_V*(-mu2*H2*c2 + mu2*sc2)*sin(mu2c2) )
    *( sin(mu1c1) + sin(mu3c3) -2.0*L0*lambda*sqrt(p1*p3/p2)/p2 )
    + 1.0/(gamma*lambda)*A_V*cos(mu2c2)
    *( (-mu1*H1*c1 + mu1*sc1)*cos(mu1c1) + (-mu3*H3*c3 + mu3*sc3)*cos(mu3c3) 
       + 2.0*L0*lambda*sqrt(p1*p3/p2)*(sp2/p2 - expansion) );

  r3 = 1.0/(gamma*lambda)*( vol*expansion*dev_Av*cos(mu3c3) 
			   -A_V*(-mu3*H3*c3 + mu3*sc3)*sin(mu3c3) )
    *( sin(mu2c2) + sin(mu1c1) -2.0*L0*lambda*sqrt(p2*p1/p3)/p3 )
    + 1.0/(gamma*lambda)*A_V*cos(mu3c3)
    *( (-mu2*H2*c2 + mu2*sc2)*cos(mu2c2) + (-mu1*H1*c1 + mu1*sc1)*cos(mu1c1) 
       + 2.0*L0*lambda*sqrt(p2*p1/p3)*(sp3/p3 - expansion) );

  return 0.5*(r1+r2+r3);
}

