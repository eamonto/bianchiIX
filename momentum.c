/*
=========================================================================
momentum.c
=========================================================================
Here are implemented the expressions for the field momentum


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


////// SCALAR FIELD MOMENTUM
double field_momentum(void)
{
  double P_phi_aux = 0.0;

  if(class_effec_switch == 0) //Classical
    {
      P_phi_aux = class_bianchi_P_phi();
    }
  else
    { //Effective
      if(Ed_A_BS_switch == 0) //Edward's Equations
	{
	  P_phi_aux = Edward_P_phi();
	}
      else if(Ed_A_BS_switch == 1) //Asieh's Equations
	{
	  if(effec_effec_switch == 0) //Does not reduce to FRW
	    {
	      P_phi_aux = Asieh_notFRW_P_phi();
	    }
	  else //Reduces to FRW
	    {
	      P_phi_aux = Asieh_FRW_P_phi();
	    }
	}
      else //Exit
	{
	  printf("\n \t Not implemented\n");
	  exit(0);
	}
    }
  
  return P_phi_aux;
}


///// CLASSICAL FIELD MOMENTUM
double class_bianchi_P_phi(void)
{
  double aux = 0.0;
  static _Bool print_warning = 0;

  aux = 2.0*( 1.0/(8.0*M_PI*G*gamma*gamma)*
		    ( p1*p2*c1*c2 + p2*p3*c2*c3 + p1*p3*c1*c3
		      + 2.0*L0*( p1*p2*c3+p2*p3*c1+p3*p1*c2 )/r0
		      + onegamma2*L0*L0/(r0*r0)*
		      ( 2.0*(p1*p1 + p2*p2 + p3*p3)
			-(p1*p2/p3)*(p1*p2/p3)
			-(p2*p3/p1)*(p2*p3/p1)
			-(p3*p1/p2)*(p3*p1/p2)
			)
		      ) -p1*p2*p3*potential()
		    );
  
  if (print_warning==0 && aux < 0.0){
    printf("\n \t WARNING, the sqrt argument is negative. \n \n");
    print_warning=1;
  }
  
  return sqrt(fabs(aux));
}


///// MOMENTUM FOR EDWARD'S QUANTIZATION
double Edward_P_phi(void)
{
  double aux = 0.0;
  static _Bool print_warning = 0;

  aux = 2.0*( p1*p2*p3/(8.0*M_PI*G*gamma*gamma*Delta)*
		    ( sin(mu1*c1)*sin(mu2*c2)
		      + sin(mu2*c2)*sin(mu3*c3)
		      + sin(mu3*c3)*sin(mu1*c1)
		      )
		    + L0/(8.0*M_PI*G*gamma*gamma*lambda)*
		    ( pow(p1*p2,1.5)*sin(mu3*c3)/sqrt(p3)
		      + pow(p2*p3,1.5)*sin(mu1*c1)/sqrt(p1)
		      + pow(p3*p1,1.5)*sin(mu2*c2)/sqrt(p2)
		      )
		    + onegamma2*L0*L0/(32.0*M_PI*G*gamma*gamma)*
		    ( 2.0*(p1*p1 + p2*p2 + p3*p3)
		      - (p1*p2/p3)*(p1*p2/p3)
		      - (p2*p3/p1)*(p2*p3/p1)
		      - (p3*p1/p2)*(p3*p1/p2)
		      )
		    - p1*p2*p3*potential()
		    );
  
  if (print_warning==0 && aux < 0.0){
    printf("\n \t WARNING, the sqrt argument is negative. \n \n");
    print_warning=1;
  }
  
  return sqrt(fabs(aux));
}


///// MOMENTUM FOR ASIEH'S QUANTIZATION THAT NOT REDUCES TO FRW
double Asieh_notFRW_P_phi(void)
{
  double aux = 0.0;
  double A_V,h_V;
  static _Bool print_warning = 0;

  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  A_V = func_A(vol);
  h_V = func_h(vol);

  aux = pow(vol,4)*A_V*pow(h_V/V_crit,6)/(8.0*M_PI*G*gamma*gamma*Delta)
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu2*c2)*sin(mu3*c3) + sin(mu3*c3)*sin(mu1*c1) )
    -L0*A_V*pow(h_V/V_crit,4.0)/(4.0*M_PI*G*gamma*gamma*lambda)
    *( p1*p1*p2*p2*sin(mu3*c3) + p2*p2*p3*p3*sin(mu1*c1) + p1*p1*p3*p3*sin(mu2*c2) )
    +L0*L0*onegamma2*A_V*pow(h_V/V_crit,4)/(8.0*M_PI*G*gamma*gamma)
    *( 2.0*vol*(p1*p1 +p2*p2 +p3*p3) 
       -( pow(p1*p2,4.0) + pow(p1*p3,4.0) + pow(p2*p3,4.0) )
       *pow(h_V/V_crit,6.0) 
       )
    -vol*potential();

  aux = aux*2.0*pow(V_crit/h_V,6)/(vol*vol);

  if (print_warning==0 && aux < 0.0){
    printf("\n \t WARNING, the sqrt argument is negative. \n \n");
    print_warning=1;
  }

  return sqrt(fabs(aux));
}


///// MOMENTUM FOR ASIEH'S QUANTIZATION THAT REDUCES TO FRW
double Asieh_FRW_P_phi(void)
{
  double aux = 0.0;
  double A_V,h_V;
  static _Bool print_warning = 0;

  vol = sqrt(p1*p2*p3);

  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  A_V = func_A(vol);
  h_V = func_h(vol);

  aux = vol*A_V/(8.0*M_PI*G*gamma*gamma*Delta)
    *( sin(mu1*c1)*sin(mu2*c2) + sin(mu2*c2)*sin(mu3*c3) + sin(mu3*c3)*sin(mu1*c1) )
    -L0*A_V/(4.0*M_PI*G*gamma*gamma*lambda)
    *( p1*p2/p3*sin(mu3*c3) + p2*p3/p1*sin(mu1*c1) + p1*p3/p2*sin(mu2*c2) )
    +L0*L0*onegamma2*A_V/(8.0*M_PI*G*gamma*gamma)
    *( 2.0*( p1*sqrt(p1/(p2*p3)) +p2*sqrt(p2/(p1*p3)) +p3*sqrt(p3/(p1*p2)) ) 
       -( p1*p2/(p3*p3)*sqrt(p1*p2/p3) 
	 +p1*p3/(p2*p2)*sqrt(p1*p3/p2) 
	 +p2*p3/(p1*p1)*sqrt(p2*p3/p1) 
	  )
       )
    -vol*potential();

  aux = aux*2.0*pow(V_crit/h_V,6)/(vol*vol);

  if (print_warning==0 && aux < 0.0){
    printf("\n \t WARNING, the sqrt argument is negative. \n \n");
    print_warning=1;
  }
  
  return sqrt(fabs(aux));
}
