/*
=========================================================================
initialize.c
=========================================================================

Initialization of all quantities, including the parameters for
the adaptive methods (Fehlberg and Cash-Karp).

*********************************************************
*Copyright (C) 2017					*
*							*
*Edison Montoya						*
*eamonto@gmail.com					*
*							*
*Up to date: Feb 16 2017				*
*							*
*********************************************************
*/

#include <header.h>


////// INITIALIZATION OF ALL QUANTITIES
int initialize_all(void)
{
  //Observables
  vol = 1.0;
  H1 = 0.0;
  H2 = 0.0;
  H3 = 0.0;
  a1 = 0.0;
  a2 = 0.0;
  a3 = 0.0;
  a_prom = 0.0;
  omega = 0.0;
  sigma2 = 0.0;
  density = 0.0;
  constraint = 0.0;
  expansion = 0.0;
  dev_expansion = 0.0;
  shear = 0.0;
  Ricci = 0.0;
  k1 = 0.0;
  k2 = 0.0;
  k3 = 0.0;
  curvature_param = 0.0;
  x1 = 0.0;
  x2 = 0.0;
  x3 = 0.0;
  
  //Time
  run_time = initial_time;
  
  //Auxiliar Constants
  Lp = sqrt(hbar*G/(C*C*C)); 
  V0 = 2.0*M_PI*M_PI*r0*r0*r0;                
  Delta = 4.0*M_PI*sqrt(3.0)*gamma*Lp*Lp; //5.170045537718 
  lambda = sqrt(Delta);                   //2.273773413891       
  density_crit = sqrt(3.0)/(32.0*M_PI*M_PI*gamma*gamma*gamma*G*G*hbar);
  V_crit = 2.0*M_PI*gamma*lambda*Lp*Lp;
  inv_V_crit = 1.0/V_crit;  
  onegamma2 = 1.0+gamma*gamma;
  invgamma = 1.0/gamma;

  L0 = epsilon*pow(V0,1.0/3.0); //L0 for Edward, B-S, Classical

  //Parameter for the dynamical equations
  if(bianchi_switch==1) //Bianchi I
    {
      L0 = 0.0;
    }
  else if(bianchi_switch==2) //Bianchi IX 
    {
      if(Ed_A_BS_switch == 1) //Asieh's Equations
	L0 = -L0/2.0; //sigma for Asieh
    }

  //Step
  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  //Dynamical variables
  /* c1 = mu1c1; */
  /* c2 = mu2c2; */
  /* c3 = mu3c3; */

  /* mu1c1 = mu1*c1; */
  /* mu2c2 = mu2*c2; */
  /* mu3c3 = mu3*c3; */

  c1 = mu1c1/mu1;
  c2 = mu2c2/mu2;
  c3 = mu3c3/mu3;

  //Initial Field Momentum
  P_phi = field_momentum();  

  if(std_out==1) printf("\n Field Momentum -> %lf",P_phi);    

  //Auxiliar for integration
  c1_p = 0.0;
  c2_p = 0.0;  
  c3_p = 0.0;
  p1_p = 0.0;
  p2_p = 0.0;
  p3_p = 0.0;
  phi_p = 0.0;
  P_phi_p = 0.0;
 
  c1_a = 0.0;
  c2_a = 0.0; 
  c3_a = 0.0;
  p1_a = 0.0; 
  p2_a = 0.0; 
  p3_a = 0.0;
  phi_a = 0.0;
  P_phi_a = 0.0;
  
  sc1 = 0.0;
  sc2 = 0.0;
  sc3 = 0.0;
  sp1 = 0.0;
  sp2 = 0.0;
  sp3 = 0.0;
  sphi = 0.0;
  sP_phi = 0.0;

  //Auxiliar for potential
  pot_switch=0;
  pot_select=0;
  
  return 0;
}
