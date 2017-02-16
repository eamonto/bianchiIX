/*
=========================================================================
aux functions.c
=========================================================================
Here are implemented the auxiliar functions.


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

////// SCALAR FIELD POTENTIAL
double potential(void)
{
  double pot = 0.0;

  return pot;
}


////// DERIVATE OF THE SCALAR FIELD POTENTIAL
double dev_pot(void)
{
  double pot = 0.0;

  return pot;
}


// FUNCTION A, ASIEH
double func_A(double v)
{ 
  double val=1.0;

  if(v < V_crit) val = v*inv_V_crit;

  return val;
}


// FUNCTION f, ASIEH
double func_f(double v)
{ 
  double val=0.0,aux;

  aux = 1.0/3.0;

  val = 1.5*inv_V_crit*pow(v,aux)*( pow(v+V_crit,aux) - pow(fabs(v-V_crit),aux) );

  return val;
}


// FUNCTION g, ASIEH
double func_g(double v)
{ 
  double val=0.0;

  val = inv_V_crit*pow(v,1.0/3.0)*( sqrt(v+V_crit) - sqrt(fabs(v-V_crit)) );

  return val;
}


// FUNCTION h, ASIEH
double func_h(double v)
{ 
  double val=0.0;

  val = sqrt(v+V_crit) - sqrt(fabs(v-V_crit));

  return val;
}


// DERIVATIVE A FUNCTION
double dev_A(double v,int i)
{ 
  double val=0.0;

  if(v < V_crit) 
  {
    if     (i==0) val = inv_V_crit;
    else if(i==1) val = inv_V_crit*sqrt(p2*p3/p1); 
    else if(i==2) val = inv_V_crit*sqrt(p1*p3/p2); 
    else if(i==3) val = inv_V_crit*sqrt(p1*p2/p3); 
  }
    
  return val;
}


// DERIVATIVE f FUNCTION
double dev_f(double v)
{ 
  double val=0.0;

  val = 0.5*inv_V_crit*pow(v,-2.0/3.0)
    *( pow(v+V_crit,1.0/3.0) - pow(fabs(v-V_crit),1.0/3.0) )
    +0.5*inv_V_crit*pow(v,1.0/3.0)
    *( pow(v+V_crit,-2.0/3.0) - (v-V_crit)*pow(fabs(v-V_crit),-5.0/3.0) );

  return val;
}


// DERIVATIVE g FUNCTION
double dev_g(double v)
{ 
  double val=0.0;

  val = inv_V_crit*pow(v,-2.0/3.0)*( sqrt(v+V_crit) - sqrt(fabs(v-V_crit)) )/3.0
    + 0.5*inv_V_crit*pow(v,1.0/3.0)*( 1.0/sqrt(v+V_crit) - sqrt(fabs(v-V_crit))/(v-V_crit) );

  return val;
}


// DERIVATIVE h FUNCTION
double dev_h(double v,int i)
{ 
  double val=0.0;

  if(i==0)      val = 0.5*( 1.0/sqrt(v+V_crit) - sqrt(fabs(v-V_crit))/(v-V_crit) );
  else if(i==1) val = 0.25*sqrt(p2*p3/p1)*( 1.0/sqrt(v+V_crit) - sqrt(fabs(v-V_crit))/(v-V_crit) );
  else if(i==2) val = 0.25*sqrt(p1*p3/p2)*( 1.0/sqrt(v+V_crit) - sqrt(fabs(v-V_crit))/(v-V_crit) );
  else if(i==3) val = 0.25*sqrt(p1*p2/p3)*( 1.0/sqrt(v+V_crit) - sqrt(fabs(v-V_crit))/(v-V_crit) );

  return val;
}


// Auxiliar Sine function to find p3 from the constraint
double fsin1(double p3_f) //Initial Conditions at Maximal Volume
{ 
  double val=0.0, aux_cons=0.0, aux_vol=0.0;

  aux_cons = pow(2.0*M_PI*M_PI,1.0/3.0);
  aux_vol = sqrt(p1*p2*p3_f);
  
  //Edward's Equations
  if(Ed_A_BS_switch == 0)
    {
      val = aux_cons*lambda*sqrt(p1*p2*p3_f)*
	( 1.0/(p2*p2) + 1.0/(p3_f*p3_f) - 1.0/(p1*p1) );
    }
  //Asieh's Equations that Reduces to FRW
  else if (Ed_A_BS_switch == 1 && effec_effec_switch == 1)
    {
      val = aux_cons*lambda*sqrt(p1*p2*p3_f)*
	    ( 1.0/(p2*p2) + 1.0/(p3_f*p3_f) - 1.0/(p1*p1) );
    }
  //Asieh's Equations that NOT Reduces to FRW
  else if (Ed_A_BS_switch == 1 && effec_effec_switch == 0)
    {
      val = aux_cons*lambda*V_crit*V_crit/
	pow( sqrt(aux_vol+V_crit) - sqrt(fabs(aux_vol-V_crit)) , 2 )*
	( 1.0/(p2*p2) + 1.0/(p3_f*p3_f) - 1.0/(p1*p1) );
    }
  
  return val;
}


// Auxiliar Sine function to find p3 from the constraint
double fsin2(double p3_f) //Initial Conditions at Maximal Volume
{ 
  double val=0.0, aux_cons=0.0, aux_vol=0.0;

  aux_cons = pow(2.0*M_PI*M_PI,1.0/3.0);
  aux_vol = sqrt(p1*p2*p3_f);
  
  //Edward's Equations
  if(Ed_A_BS_switch == 0)
    {
      val = aux_cons*lambda*sqrt(p1*p2*p3_f)*
	( 1.0/(p1*p1) + 1.0/(p3_f*p3_f) - 1.0/(p2*p2) );
    }
  //Asieh's Equations that Reduces to FRW
  else if (Ed_A_BS_switch == 1 && effec_effec_switch == 1)
    {
      val = aux_cons*lambda*sqrt(p1*p2*p3_f)*
	( 1.0/(p1*p1) + 1.0/(p3_f*p3_f) - 1.0/(p2*p2) );
    }
  //Asieh's Equations that NOT Reduces to FRW
  else if (Ed_A_BS_switch == 1 && effec_effec_switch == 0)
    {
      val = aux_cons*lambda*V_crit*V_crit/
	pow( sqrt(aux_vol+V_crit) - sqrt(fabs(aux_vol-V_crit)) , 2 )*
	( 1.0/(p1*p1) + 1.0/(p3_f*p3_f) - 1.0/(p2*p2) );
    }
  
  return val;
}


// Auxiliar Sine function to find p3 from the constraint
double fsin3(double p3_f) //Initial Conditions at Maximal Volume
{ 
  double val=0.0, aux_cons=0.0, aux_vol=0.0;

  aux_cons = pow(2.0*M_PI*M_PI,1.0/3.0);
  aux_vol = sqrt(p1*p2*p3_f);

  //Edward's Equations
  if(Ed_A_BS_switch == 0)
    {
      val = aux_cons*lambda*sqrt(p1*p2*p3_f)*
	( 1.0/(p1*p1) + 1.0/(p2*p2) - 1.0/(p3_f*p3_f) );
    }
  //Asieh's Equations that Reduces to FRW
  else if (Ed_A_BS_switch == 1 && effec_effec_switch == 1)
    {
      val = aux_cons*lambda*sqrt(p1*p2*p3_f)*
	( 1.0/(p1*p1) + 1.0/(p2*p2) - 1.0/(p3_f*p3_f) );
    }
  //Asieh's Equations that NOT Reduces to FRW
  else if (Ed_A_BS_switch == 1 && effec_effec_switch == 0)
    {
      val = aux_cons*lambda*V_crit*V_crit/
	pow( sqrt(aux_vol+V_crit) - sqrt(fabs(aux_vol-V_crit)) , 2 )*
	( 1.0/(p1*p1) + 1.0/(p2*p2) - 1.0/(p3_f*p3_f) );
    }

  return val;
}


//Gravitational Constraint as function of p3
double func_grav_constraint(double p3_f)
{
  double val=0.0, aux_cons=0.0, aux_vol=0.0;

  aux_cons = pow(2.0*M_PI*M_PI,1.0/3.0);
  aux_vol = sqrt(p1*p2*p3_f);

  //Edward's Equations and Asieh's Equations that Reduces to FRW  
  if(Ed_A_BS_switch == 0 || (Ed_A_BS_switch == 1 && effec_effec_switch == 1))
    {
      val = p1*p2*p3_f/Delta*
	(fsin1(p3_f)*fsin2(p3_f) +fsin1(p3_f)*fsin3(p3_f) +fsin2(p3_f)*fsin3(p3_f))
	-2.0*aux_cons/lambda*
	( p1*p2*sqrt(p1*p2/p3_f)*fsin3(p3_f)
	 +p1*p3_f*sqrt(p1*p3_f/p2)*fsin2(p3_f)
         +p2*p3_f*sqrt(p2*p3_f/p1)*fsin1(p3_f)
	)
	+aux_cons*aux_cons*onegamma2*
	( 2.0*(p1*p1 + p2*p2 + p3_f*p3_f)
	  -( pow(p1*p2,13.0/4.0)/sqrt(p3_f)
	    +pow(p1*p3_f,13.0/4.0)/sqrt(p2)
	    +pow(p2*p3_f,13.0/4.0)/sqrt(p1)
	   )
	);
    }

  //Asieh's Equations that NOT Reduces to FRW  
  if(Ed_A_BS_switch == 1 && effec_effec_switch == 0)
    {
      val = pow(aux_vol,3)/(V_crit*V_crit*Delta)*
	pow( sqrt(aux_vol+V_crit) - sqrt(fabs(aux_vol-V_crit)) , 2 )*
	(fsin1(p3_f)*fsin2(p3_f) +fsin1(p3_f)*fsin3(p3_f) +fsin2(p3_f)*fsin3(p3_f))
	-2.0*aux_cons/lambda*
	( p1*p2*sqrt(p1*p2/p3_f)*fsin3(p3_f)
	 +p1*p3_f*sqrt(p1*p3_f/p2)*fsin2(p3_f)
         +p2*p3_f*sqrt(p2*p3_f/p1)*fsin1(p3_f)
	)
	+aux_cons*aux_cons*onegamma2*
	( 2.0*(p1*p1 + p2*p2 + p3_f*p3_f)
	  -( pow(p1*p2,13.0/4.0)/sqrt(p3_f)
	    +pow(p1*p3_f,13.0/4.0)/sqrt(p2)
	    +pow(p2*p3_f,13.0/4.0)/sqrt(p1)
	   )
	);
     }

  return val;
}


double find_max_vol(void)
{
  double error_aux,zero,step,center;
  int i;

  error_aux = 0.0;
  center = 0.0;
  step = p1*p2;
  i = 0;

  do
    {
      zero = center + step;
	      
      error_aux = -func_grav_constraint(zero);
	      
      if(error_aux<0.0) center = zero; 

      if(error_aux>0.0)  step = step/2.0;

      i++;
    }while( (fabs(error_aux) > 1.0e-13) && (i<100000) && step > 1.0e-15);

  printf("\n \t Error =%d , zero(p3) = %e, iterations =%d\n",(int)(error_aux) ,zero,i);

  return zero;
}
