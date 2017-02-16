/*
=========================================================================
sources.c
=========================================================================
Here are implemented the right hand side of the equations 
that we want to integrate.


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


///// RIGHT HAND SIDE /////
int rhs(void)
{  

  if(class_effec_switch == 0) //Classical
    {
      class_bianchi_rhs();
    }
  else
    { //Effective
      
      if(Ed_A_BS_switch == 0) //Edward's Equations
	{
	  Edward_rhs();
	}
      else if(Ed_A_BS_switch == 1) //Asieh's Equations
	{
	  if(effec_effec_switch == 0) //Does not reduce to FRW
	    {
	      Asieh_notFRW_rhs();
	    }
	  else //Reduces to FRW
	    {	 
	      Asieh_FRW_rhs();
	    }
	}
      else  //Exit
	{
	  printf("\n \t Not implemented\n");
	  exit(0);
	}
    }
  
  return 0;
}


// CLASSICAL BIANCHI IX
int class_bianchi_rhs(void)
{
  double inv_vol;

  vol = sqrt(p1*p2*p3);

  inv_vol = 1.0/vol;

  sc1 = -invgamma*inv_vol*( p2*c1*c2 + p3*c1*c3 + L0*(p2*c3+p3*c2) 
		    + 0.5*L0*L0*onegamma2/p1*
		    ( 2.0*p1*p1
		      + (p2*p3/p1)*(p2*p3/p1)
		      - (p1*p2/p3)*(p1*p2/p3) 
		      - (p1*p3/p2)*(p1*p3/p2)
		      ) 
		    );
  
  sc2 = -invgamma*inv_vol*( p1*c2*c1 + p3*c2*c3 + L0*(p1*c3+p3*c1) 
		    + 0.5*L0*L0*onegamma2/p2*
		    ( 2.0*p2*p2 
		      + (p1*p3/p2)*(p1*p3/p2) 
		      - (p1*p2/p3)*(p1*p2/p3) 
		      - (p2*p3/p1)*(p2*p3/p1) 
		      ) 
		    );
  
  sc3 = -invgamma*inv_vol*( p1*c3*c1 + p2*c3*c2 + L0*(p2*c1+p1*c2) 
		    + 0.5*L0*L0*onegamma2/p3*
		    ( 2.0*p3*p3 
		      + (p2*p1/p3)*(p2*p1/p3) 
		      - (p3*p2/p1)*(p3*p2/p1) 
		      - (p1*p3/p2)*(p1*p3/p2) 
		      ) 
		    );
  
  
  sp1 = invgamma*inv_vol*(p1*p2*c2 + p1*p3*c3 + L0*p2*p3);
  
  sp2 = invgamma*inv_vol*(p2*p1*c1 + p2*p3*c3 + L0*p1*p3);
  
  sp3 = invgamma*inv_vol*(p3*p1*c1 + p3*p2*c2 + L0*p1*p2);
  
  if(pot_switch==1)
    {
      sphi   = P_phi*inv_vol;
      sP_phi =-vol*dev_pot();
    }
      
  return 0;
}


// EDWARD EQUATIONS
int Edward_rhs(void)
{
  double inv_vol;

  vol = sqrt(p1*p2*p3);

  inv_vol = 1.0/vol;

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;


  sc1 = -invgamma*inv_vol*( p2*p3/(Delta)*
		    ( sin(mu1*c1)*sin(mu2*c2)  
		      + sin(mu1*c1)*sin(mu3*c3) 
		      + sin(mu2*c2)*sin(mu3*c3) 
		      + 0.5*mu1*c1*cos(mu1*c1)*( sin(mu2*c2)+sin(mu3*c3) ) 
		      - 0.5*mu2*c2*cos(mu2*c2)*( sin(mu1*c1)+sin(mu3*c3) ) 
		      - 0.5*mu3*c3*cos(mu3*c3)*( sin(mu1*c1)+sin(mu2*c2) ) 
		      ) 
		    + L0*( 1.5/mu1*( p1*p2*sin(mu3*c3)/p3 
				     + p1*p3*sin(mu2*c2)/p2 
				     - p2*p3*sin(mu1*c1)/(3.0*p1) 
				     ) 
			   + 0.5*p2*p3*c1*cos(mu1*c1)/p1 
			   - 0.5*p2*c3*cos(mu3*c3)
			     - 0.5*p3*c2*cos(mu2*c2) 
			   ) 
		    + 0.5*L0*L0*onegamma2*p1*
		    ( 2.0 -(p2/p3)*(p2/p3) -(p3/p2)*(p3/p2) + (p2*p3/p1)*(p2*p3/p1)/(p1*p1) ) 
		    ) + 8.0*M_PI*G*gamma*p2*p3*potential();
  
  sc2 = -invgamma*inv_vol*( p1*p3/(Delta)*
		    ( sin(mu1*c1)*sin(mu2*c2)  
		      + sin(mu1*c1)*sin(mu3*c3) 
		      + sin(mu2*c2)*sin(mu3*c3) 
		      - 0.5*mu1*c1*cos(mu1*c1)*( sin(mu2*c2)+sin(mu3*c3) ) 
		      + 0.5*mu2*c2*cos(mu2*c2)*( sin(mu1*c1)+sin(mu3*c3) ) 
		      - 0.5*mu3*c3*cos(mu3*c3)*( sin(mu1*c1)+sin(mu2*c2) ) 
		      ) 
		    + L0*( 1.5/mu2*( p2*p3*sin(mu1*c1)/p1 
				     + p1*p2*sin(mu3*c3)/p3 
				     - p1*p3*sin(mu2*c2)/(3.0*p2) 
				     ) 
			   +0.5*p1*p3*c2*cos(mu2*c2)/p2 
			   -0.5*p3*c1*cos(mu1*c1) 
			   -0.5*p1*c3*cos(mu3*c3) 
			   )
		    + 0.5*L0*L0*onegamma2*p2*
		    ( 2.0 -(p1/p3)*(p1/p3) -(p3/p1)*(p3/p1) + (p1*p3/p2)*(p1*p3/p2)/(p2*p2) ) 
		    )  + 8.0*M_PI*G*gamma*p1*p3*potential();
  
  sc3 = -invgamma*inv_vol*( p1*p2/(Delta)*
		    ( sin(mu1*c1)*sin(mu2*c2)  
		      + sin(mu1*c1)*sin(mu3*c3) 
		      + sin(mu2*c2)*sin(mu3*c3) 
		      - 0.5*mu1*c1*cos(mu1*c1)*( sin(mu2*c2)+sin(mu3*c3) ) 
		      - 0.5*mu2*c2*cos(mu2*c2)*( sin(mu1*c1)+sin(mu3*c3) ) 
		      + 0.5*mu3*c3*cos(mu3*c3)*( sin(mu1*c1)+sin(mu2*c2) ) 
		      ) 
		    + L0*( 1.5/mu3*( p2*p3*sin(mu1*c1)/p1  
				     + p1*p3*sin(mu2*c2)/p2 
				     - p1*p2*sin(mu3*c3)/(3.0*p3) 
				     ) 
			   + 0.5*p1*p2*c3*cos(mu3*c3)/p3 
			   - 0.5*p2*c1*cos(mu1*c1) 
			   - 0.5*p1*c2*cos(mu2*c2) 
			   )
		    + 0.5*L0*L0*onegamma2*p3*
		    (2.0 -(p1/p2)*(p1/p2) -(p2/p1)*(p2/p1) + (p1*p2/p3)*(p1*p2/p3)/(p3*p3) ) 
		    ) + 8.0*M_PI*G*gamma*p1*p2*potential();
  

  sp1 = invgamma*inv_vol*( p1*p1/mu1*(sin(mu2*c2) + sin(mu3*c3)) + L0*p2*p3 )*cos(mu1*c1);
      
  sp2 = invgamma*inv_vol*( p2*p2/mu2*(sin(mu1*c1) + sin(mu3*c3)) + L0*p1*p3 )*cos(mu2*c2);
  
  sp3 = invgamma*inv_vol*( p3*p3/mu3*(sin(mu1*c1) + sin(mu2*c2)) + L0*p1*p2 )*cos(mu3*c3);
  
  if(pot_switch==1)
    {
      sphi   = P_phi*inv_vol;
      sP_phi =-vol*dev_pot();
    }
    
  return 0;
}


// ASIEH EQUATIONS, NOT REDUCE TO FRW
int Asieh_notFRW_rhs(void)
{
  double A_V,h_V,dev_A1,dev_A2,dev_A3,dev_h1,dev_h2,dev_h3;

  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;
  
  A_V = func_A(vol);
  h_V = func_h(vol);

  dev_A1 = dev_A(vol,1);
  dev_A2 = dev_A(vol,2);
  dev_A3 = dev_A(vol,3);

  dev_h1 = dev_h(vol,1);
  dev_h2 = dev_h(vol,2);
  dev_h3 = dev_h(vol,3);

  sp1 = A_V*vol*pow(h_V/V_crit,4)/gamma*cos(mu1*c1)
    *( pow(vol*h_V/V_crit,2)*p1/lambda*( sin(mu2*c2)+sin(mu3*c3) ) - 2.0*L0*p2*p3 );

  sp2 = A_V*vol*pow(h_V/V_crit,4)/gamma*cos(mu2*c2)
    *( pow(vol*h_V/V_crit,2)*p2/lambda*( sin(mu1*c1)+sin(mu3*c3) ) -2.0*L0*p1*p3 );

  sp3 = A_V*vol*pow(h_V/V_crit,4)/gamma*cos(mu3*c3)
    *( pow(vol*h_V/V_crit,2)*p3/lambda*( sin(mu1*c1)+sin(mu2*c2) ) -2.0*L0*p1*p2 );

  ///////////////////////////////////C1
  sc1 = -pow(h_V/V_crit,5)/(V_crit*gamma*Delta)
    *( 2.0*p2*p2*p3*p3*p1*A_V*h_V + pow(vol,4)*(dev_A1*h_V+6.0*A_V*dev_h1) )
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu1*c1)*sin(mu3*c3) + sin(mu2*c2)*sin(mu3*c3) )
    +2.0*L0*pow(h_V/V_crit,3)/(V_crit*gamma*lambda)
    *( sin(mu3*c3)*(2.0*p1*p2*p2*A_V*h_V + p1*p1*p2*p2*(dev_A1*h_V+4.0*A_V*dev_h1) )
      +sin(mu2*c2)*(2.0*p1*p3*p3*A_V*h_V + p1*p1*p3*p3*(dev_A1*h_V+4.0*A_V*dev_h1) )
      +sin(mu1*c1)*(                       p2*p2*p3*p3*(dev_A1*h_V+4.0*A_V*dev_h1) )
       )
    -A_V*pow(h_V/V_crit,4)*c1*cos(mu1*c1)/(2.0*gamma)
    *( vol*pow(vol*h_V/V_crit,2)/lambda*( sin(mu2*c2) +sin(mu3*c3) ) 
       -2.0*L0*p2*p3*sqrt(p2*p3/p1)
       )
    +A_V*pow(h_V/V_crit,4)*c2*cos(mu2*c2)/(2.0*gamma)
    *( p2*p2*p3*vol*h_V*h_V/(V_crit*V_crit*lambda)*( sin(mu1*c1) +sin(mu3*c3) )
       -2.0*L0*p3*vol
       )
    +A_V*pow(h_V/V_crit,4)*c3*cos(mu3*c3)/(2.0*gamma)
    *( p3*p3*p2*vol*h_V*h_V/(V_crit*V_crit*lambda)*( sin(mu1*c1) +sin(mu2*c2) )
       -2.0*L0*p2*vol
       )
    -L0*L0*onegamma2*pow(h_V/V_crit,3)/(V_crit*gamma)
    *( 4.0*p1*A_V*vol*h_V + (p1*p1 +p2*p2 + p3*p3)
       *( sqrt(p2*p3/p1)*A_V*h_V + 8.0*A_V*vol*dev_h1 + 2.0*dev_A1*vol*h_V )
       -4.0*V_crit*pow(h_V/V_crit,7.0)*A_V*pow(p1,3.0)
       *( pow(p2,4.0) + pow(p3,4.0) )
       -pow(h_V/V_crit,6.0)*( 10.0*dev_h1*A_V + h_V*dev_A1 )
       *( pow(p1*p2,4.0) + pow(p1*p3,4.0) + pow(p2*p3,4.0) )
       )
    +4.0*M_PI*G*gamma*sqrt(p2*p3/p1)*potential()
    +4.0*M_PI*G*gamma*P_phi*P_phi*pow(h_V/V_crit,5)/V_crit
    *( p2*p3*h_V + 6.0*vol*vol*dev_h1 );

  ///////////////////////////////////C2
  sc2 = -pow(h_V/V_crit,5)/(V_crit*gamma*Delta)
    *( 2.0*p1*p1*p3*p3*p2*A_V*h_V + pow(vol,4)*(dev_A2*h_V+6.0*A_V*dev_h2) )
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu1*c1)*sin(mu3*c3) + sin(mu2*c2)*sin(mu3*c3) )
    +2.0*L0*pow(h_V/V_crit,3)/(V_crit*gamma*lambda)
    *( sin(mu3*c3)*(2.0*p2*p1*p1*A_V*h_V + p1*p1*p2*p2*(dev_A2*h_V+4.0*A_V*dev_h2) )
      +sin(mu1*c1)*(2.0*p2*p3*p3*A_V*h_V + p2*p2*p3*p3*(dev_A2*h_V+4.0*A_V*dev_h2) )
      +sin(mu2*c2)*(                       p1*p1*p3*p3*(dev_A2*h_V+4.0*A_V*dev_h2) )
       )
    -A_V*pow(h_V/V_crit,4)*c2*cos(mu2*c2)/(2.0*gamma)
    *( vol*pow(vol*h_V/V_crit,2)/lambda*( sin(mu1*c1) +sin(mu3*c3) ) 
       -2.0*L0*p1*p3*sqrt(p1*p3/p2)
       )
    +A_V*pow(h_V/V_crit,4)*c1*cos(mu1*c1)/(2.0*gamma)
    *( p1*p1*p3*vol*h_V*h_V/(V_crit*V_crit*lambda)*( sin(mu2*c2) +sin(mu3*c3) )
       -2.0*L0*p3*vol
       )
    +A_V*pow(h_V/V_crit,4)*c3*cos(mu3*c3)/(2.0*gamma)
    *( p3*p3*p1*vol*h_V*h_V/(V_crit*V_crit*lambda)*( sin(mu1*c1) +sin(mu2*c2) )
       -2.0*L0*p1*vol
       )
    -L0*L0*onegamma2*pow(h_V/V_crit,3)/(V_crit*gamma)
    *( 4.0*p2*A_V*vol*h_V + (p1*p1 +p2*p2 + p3*p3)
       *( sqrt(p1*p3/p2)*A_V*h_V + 8.0*A_V*vol*dev_h2 + 2.0*dev_A2*vol*h_V )
       -4.0*V_crit*pow(h_V/V_crit,7.0)*A_V*pow(p2,3.0)
       *( pow(p1,4.0) + pow(p3,4.0) )
       -pow(h_V/V_crit,6.0)*( 10.0*dev_h2*A_V + h_V*dev_A2 )
       *( pow(p1*p2,4.0) + pow(p1*p3,4.0) + pow(p2*p3,4.0) )
       )
    +4.0*M_PI*G*gamma*sqrt(p1*p3/p2)*potential()
    +4.0*M_PI*G*gamma*P_phi*P_phi*pow(h_V/V_crit,5)/V_crit
    *( p1*p3*h_V + 6.0*vol*vol*dev_h2 );

  ///////////////////////////////////C3
  sc3 = -pow(h_V/V_crit,5)/(V_crit*gamma*Delta)
    *( 2.0*p2*p2*p1*p1*p3*A_V*h_V + pow(vol,4)*(dev_A3*h_V+6.0*A_V*dev_h3) )
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu1*c1)*sin(mu3*c3) + sin(mu2*c2)*sin(mu3*c3) )
    +2.0*L0*pow(h_V/V_crit,3)/(V_crit*gamma*lambda)
    *( sin(mu1*c1)*(2.0*p3*p2*p2*A_V*h_V + p2*p2*p3*p3*(dev_A3*h_V+4.0*A_V*dev_h3) )
      +sin(mu2*c2)*(2.0*p3*p1*p1*A_V*h_V + p1*p1*p3*p3*(dev_A3*h_V+4.0*A_V*dev_h3) )
      +sin(mu3*c3)*(                       p1*p1*p2*p2*(dev_A3*h_V+4.0*A_V*dev_h3) )
       )
    -A_V*pow(h_V/V_crit,4)*c3*cos(mu3*c3)/(2.0*gamma)
    *( vol*pow(vol*h_V/V_crit,2)/lambda*( sin(mu1*c1) +sin(mu2*c2) ) 
       -2.0*L0*p1*p2*sqrt(p1*p2/p3)
       )
    +A_V*pow(h_V/V_crit,4)*c2*cos(mu2*c2)/(2.0*gamma)
    *( p2*p2*p1*vol*h_V*h_V/(V_crit*V_crit*lambda)*( sin(mu1*c1) +sin(mu3*c3) )
       -2.0*L0*p1*vol
       )
    +A_V*pow(h_V/V_crit,4)*c1*cos(mu1*c1)/(2.0*gamma)
    *( p1*p1*p2*vol*h_V*h_V/(V_crit*V_crit*lambda)*( sin(mu2*c2) +sin(mu3*c3) )
       -2.0*L0*p2*vol
       )
    -L0*L0*onegamma2*pow(h_V/V_crit,3)/(V_crit*gamma)
    *( 4.0*p3*A_V*vol*h_V + (p1*p1 +p2*p2 + p3*p3)
       *( sqrt(p1*p2/p3)*A_V*h_V + 8.0*A_V*vol*dev_h3 + 2.0*dev_A3*vol*h_V )
       -4.0*V_crit*pow(h_V/V_crit,7.0)*A_V*pow(p3,3.0)
       *( pow(p1,4.0) + pow(p2,4.0) )
       -pow(h_V/V_crit,6.0)*( 10.0*dev_h3*A_V + h_V*dev_A3 )
       *( pow(p1*p2,4.0) + pow(p1*p3,4.0) + pow(p2*p3,4.0) )
       )
    +4.0*M_PI*G*gamma*sqrt(p1*p2/p3)*potential()
    +4.0*M_PI*G*gamma*P_phi*P_phi*pow(h_V/V_crit,5)/V_crit
    *( p1*p2*h_V + 6.0*vol*vol*dev_h3 );

  if(pot_switch==1)
    {
      sphi   = P_phi*vol*vol*pow(h_V/V_crit,6);
      sP_phi =-vol*dev_pot();
    }

  return 0;
}


// ASIEH EQUATIONS, REDUCE TO FRW
int Asieh_FRW_rhs(void)
{
  double A_V,h_V,dev_A1,dev_A2,dev_A3,dev_h1,dev_h2,dev_h3;

  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;
  
  A_V = func_A(vol);
  h_V = func_h(vol);

  dev_A1 = dev_A(vol,1);
  dev_A2 = dev_A(vol,2);
  dev_A3 = dev_A(vol,3);

  dev_h1 = dev_h(vol,1);
  dev_h2 = dev_h(vol,2);
  dev_h3 = dev_h(vol,3);

  sp1 = A_V*cos(mu1*c1)/(gamma*lambda)
    *( p1*( sin(mu2*c2) + sin(mu3*c3) ) -2.0*L0*lambda*sqrt(p2*p3/p1) );

  sp2 = A_V*cos(mu2*c2)/(gamma*lambda)
    *( p2*( sin(mu1*c1) + sin(mu3*c3) ) -2.0*L0*lambda*sqrt(p1*p3/p2) );

  sp3 = A_V*cos(mu3*c3)/(gamma*lambda)
    *( p3*( sin(mu1*c1) + sin(mu2*c2) ) -2.0*L0*lambda*sqrt(p1*p2/p3) );

  ///////////////////////////////////C1
  sc1 = -1.0/(gamma*Delta)
    *( 0.5*sqrt(p2*p3/p1)*A_V + dev_A1*vol )
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu1*c1)*sin(mu3*c3) + sin(mu2*c2)*sin(mu3*c3) )
    -2.0*L0/(gamma*lambda)
    *( ( A_V*p2*p3/(p1*p1) - dev_A1*p2*p3/p1 )*sin(mu1*c1)
      -( A_V*p3/p2 + dev_A1*p1*p3/p2 )*sin(mu2*c2)
      -( A_V*p2/p3 + dev_A1*p1*p2/p3 )*sin(mu3*c3)
       )
    -A_V*c1*cos(mu1*c1)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p2*p3/p1)/p1 + sin(mu2*c2) +sin(mu3*c3) )
    +A_V*c2*cos(mu2*c2)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p3/(p1*p2)) + p2/p1*( sin(mu1*c1)+sin(mu3*c3) ) )
    +A_V*c3*cos(mu3*c3)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p2/(p1*p3)) + p3/p1*( sin(mu1*c1)+sin(mu2*c2) ) )
    -A_V*L0*L0*onegamma2/gamma
    *( 5.0/2.0*pow(p2*p3/p1,1.5)/(p1*p1) - (p2*p2+p3*p3)/(p1*vol) 
       -1.5*sqrt(p1)*( pow(p2/p3,1.5)/p3+pow(p3/p2,1.5)/p2 ) + 3.0*sqrt(p1/(p2*p3))
       ) 
    +8.0*M_PI*G*gamma*sqrt(p2*p3/p1)*potential()
    +4.0*M_PI*G*gamma*P_phi*P_phi*pow(h_V/V_crit,5)/V_crit
    *( p2*p3*h_V + 6.0*vol*vol*dev_h1 );

  ///////////////////////////////////C2
  sc2 = -1.0/(gamma*Delta)
    *( 0.5*sqrt(p1*p3/p2)*A_V + dev_A2*vol )
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu1*c1)*sin(mu3*c3) + sin(mu2*c2)*sin(mu3*c3) )
    -2.0*L0/(gamma*lambda)
    *( ( A_V*p1*p3/(p2*p2) - dev_A2*p1*p3/p2 )*sin(mu2*c2)
      -( A_V*p3/p1 + dev_A2*p2*p3/p1 )*sin(mu1*c1)
      -( A_V*p1/p3 + dev_A2*p1*p2/p3 )*sin(mu3*c3)
       )
    -A_V*c2*cos(mu2*c2)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p1*p3/p2)/p2 + sin(mu1*c1) +sin(mu3*c3) )
    +A_V*c1*cos(mu1*c1)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p3/(p1*p2)) + p1/p2*( sin(mu2*c2)+sin(mu3*c3) ) )
    +A_V*c3*cos(mu3*c3)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p1/(p2*p3)) + p3/p2*( sin(mu1*c1)+sin(mu2*c2) ) )
    -A_V*L0*L0*onegamma2/gamma
    *( 5.0/2.0*pow(p1*p3/p2,1.5)/(p2*p2) - (p1*p1+p3*p3)/(p2*vol) 
       -1.5*sqrt(p2)*( pow(p1/p3,1.5)/p3+pow(p3/p1,1.5)/p1 ) + 3.0*sqrt(p2/(p1*p3))
       ) 
    +8.0*M_PI*G*gamma*sqrt(p1*p3/p2)*potential()
    +4.0*M_PI*G*gamma*P_phi*P_phi*pow(h_V/V_crit,5)/V_crit
    *( p1*p3*h_V + 6.0*vol*vol*dev_h2 );

  ///////////////////////////////////C3
  sc3 = -1.0/(gamma*Delta)
    *( 0.5*sqrt(p2*p1/p3)*A_V + dev_A3*vol )
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu1*c1)*sin(mu3*c3) + sin(mu2*c2)*sin(mu3*c3) )
    -2.0*L0/(gamma*lambda)
    *( ( A_V*p2*p1/(p3*p3) - dev_A3*p2*p1/p3 )*sin(mu3*c3)
      -( A_V*p1/p2 + dev_A3*p1*p3/p2 )*sin(mu2*c2)
      -( A_V*p2/p1 + dev_A3*p3*p2/p1 )*sin(mu1*c1)
       )
    -A_V*c3*cos(mu3*c3)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p2*p1/p3)/p3 + sin(mu2*c2) +sin(mu1*c1) )
    +A_V*c2*cos(mu2*c2)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p1/(p3*p2)) + p2/p3*( sin(mu1*c1)+sin(mu3*c3) ) )
    +A_V*c1*cos(mu1*c1)/(2.0*gamma*lambda)
    *( -2.0*L0*lambda*sqrt(p2/(p1*p3)) + p1/p3*( sin(mu3*c3)+sin(mu2*c2) ) )
    -A_V*L0*L0*onegamma2/gamma
    *( 5.0/2.0*pow(p2*p1/p3,1.5)/(p3*p3) - (p2*p2+p1*p1)/(p3*vol) 
       -1.5*sqrt(p3)*( pow(p2/p1,1.5)/p1+pow(p1/p2,1.5)/p2 ) + 3.0*sqrt(p3/(p2*p1))
       ) 
    +8.0*M_PI*G*gamma*sqrt(p2*p1/p3)*potential()
    +4.0*M_PI*G*gamma*P_phi*P_phi*pow(h_V/V_crit,5)/V_crit
    *( p2*p1*h_V + 6.0*vol*vol*dev_h3 );

  if(pot_switch==1)
    {
      sphi   = P_phi*vol*vol*pow(h_V/V_crit,6);
      sP_phi =-vol*dev_pot();
    }

  return 0;
}
