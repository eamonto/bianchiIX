/*
=========================================================================
observables.c
=========================================================================
Here are compute all the quantities that give relevant 
information about the system dynamics.


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


////// COMPUTE OBSERVABLES
int compute_obs(void)
{
  //  double x1,x2,x3;
  double P_phi_aux;
  double h_V;
  
  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  mu1c1 = mu1*c1;
  mu2c2 = mu2*c2;
  mu3c3 = mu3*c3;
  
  if(L0==0.0) {
    a1 = sqrt(p2*p3/p1)/L1;
    a2 = sqrt(p1*p3/p2)/L2;
    a3 = sqrt(p1*p2/p3)/L3;
  }else{
    a1 = sqrt(p2*p3/p1)/fabs(L0);
    a2 = sqrt(p1*p3/p2)/fabs(L0);
    a3 = sqrt(p1*p2/p3)/fabs(L0);
    
    if(Ed_A_BS_switch == 1) //Asieh's Equations
      {
	a1 = a1/2.0;
	a2 = a2/2.0;
	a3 = a3/2.0;
      }
  }

  a_prom = pow(a1*a2*a3,1.0/3.0);

  vol = sqrt(p1*p2*p3);
  
  rhs();  //Sources in the cosmic time
  
  H1 = (sp2/p2 + sp3/p3 - sp1/p1)/2.0;
  H2 = (sp1/p1 + sp3/p3 - sp2/p2)/2.0;
  H3 = (sp1/p1 + sp2/p2 - sp3/p3)/2.0;
  
  expansion = H1+H2+H3;
  
  shear = ((H1-H2)*(H1-H2)+(H2-H3)*(H2-H3)+(H3-H1)*(H3-H1))/3.0;
  
  P_phi_aux = field_momentum();
  
  constraint = fabs(P_phi_aux-P_phi)/P_phi_aux;
  
  if(Ed_A_BS_switch == 0) //Edward's Equations
    density = P_phi*P_phi/(2.0*p1*p2*p3) + potential();

  if(Ed_A_BS_switch == 1) //Asieh's Equations
    {
      h_V = func_h(vol);
      density = P_phi*P_phi*vol*pow(h_V/V_crit,6)/2.0 + potential();
    }

  // "Gravitational potential" for Edward's equations
  x1 = sqrt(p2*p3/(p1*p1*p1));
  x2 = sqrt(p1*p3/(p2*p2*p2));
  x3 = sqrt(p1*p2/(p3*p3*p3));
  
  if(fabs(expansion) > 1.0e-15) //Avoids division by zero
    {
      omega = 24.0*M_PI*G*density/(expansion*expansion);

      // omega = 24.0*M_PI*G*density/(expansion*expansion+0.25*pow(x1*x2*x3,2.0/3.0));

      sigma2 = 3.0*shear/(2.0*expansion*expansion);
            
      // Kasner Exponents
      k1 = H1/fabs(expansion);
      k2 = H2/fabs(expansion);
      k3 = H3/fabs(expansion);

      //Curvature Parameter for Edward's equations
      curvature_param = 3.0*L0*L0/(4.0*expansion*expansion)*
	(x1*x1 + x2*x2 + x3*x3 - 2.0*(x1*x2 + x1*x3+ x2*x3));

      if(Ed_A_BS_switch == 1) //Asieh's Equations
	curvature_param = 4.0*curvature_param; //This will change

    }
  

  //Derivative of the Expansion
  if(class_effec_switch == 0) //Classical
    {
      dev_expansion = class_thetadot();
    }
  else
    { //Effective
      if(Ed_A_BS_switch == 0) //Edward's Equations
	{
	  dev_expansion = Edward_thetadot();
	}
      else if(Ed_A_BS_switch == 1) //Asieh's Equations
	{
	  if(effec_effec_switch == 0) //Does not reduce to FRW
	    {
	      dev_expansion = Asieh_notFRW_thetadot();
	    }
	  else //Reduces to FRW
	    {
	      dev_expansion = Asieh_FRW_thetadot();
	    }
	}
      else //Exit
	{
	  printf("\n \t Not implemented\n");
	  exit(0);
	}
    }
  
  Ricci = 2.0*dev_expansion + (sp1/p1)*(sp1/p1)+(sp2/p2)*(sp2/p2)+(sp3/p3)*(sp3/p3);
  
  return 0;
}


