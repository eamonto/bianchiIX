/*
=================================================================
header.h
=================================================================
This is the head of the program, here the global variables 
are declared. The program was developed with three different 
integrators: RK4, RK-Felberg and RK-Cash-Karp. These methods 
are implemented in the rk4.c and rk45.c files.


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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define G        1.0 //6.6732e-11       //Newton's Constant
#define hbar     1.0 //1.05459e-34      //Planck's Constant
#define C        1.0                    //Speed of Light 
#define gamma    0.23753295796592       //Barbero-Immirzi Parameter

#define epsilon  1.0                    //Orientation Parameter
#define L1       1.0                    //Fiducial Lenght 1
#define L2       1.0                    //Fiducial Lenght 2
#define L3       1.0                    //Fiducial Lenght 3
#define r0       2.0                    //Fiducial Radius


///////////////// GLOBAL VARIABLES /////////////////

////// DYNAMICAL VARIABLES
double c1,c2,c3,p1,p2,p3;
double phi,P_phi;

////// IN-OUT VARIABLES
char *initfile,*outputfile;
int std_out, time_output;

////// TIME PARAMETERS
double initial_time, final_time, run_time, dt;

////// PARAMETERS FOR THE DYNAMICAL EQUATIONS
//int Integrator;
int bianchi_switch;
int class_effec_switch;
int Ed_A_BS_switch;
int effec_effec_switch;
double L0;
int pot_switch;
int pot_select;

////// OBSERVABLES  
double vol;
double H1, H2, H3;
double a1, a2, a3;
double a_prom;
double omega;
double sigma2;
double density;
double constraint;
double expansion;
double dev_expansion;
double shear;
double Ricci;
double k1,k2,k3;
double curvature_param;
double x1,x2,x3;

//Auxiliar Constants
double Lp;
double V0;
double Delta;
double lambda;
double density_crit;
double V_crit;
double inv_V_crit;
double onegamma2;
double invgamma; 

//AUXILIARY DYNAMICAL VARIABLES
double mu1, mu2, mu3;
double mu1c1, mu2c2, mu3c3;

////// AUXILIARS FOR INTEGRATION 
double c1_p, c2_p, c3_p;
double p1_p, p2_p, p3_p;

double c1_a, c2_a, c3_a;
double p1_a, p2_a, p3_a;

double sc1, sc2, sc3;
double sp1, sp2, sp3;

double phi_p, phi_a, sphi;
double P_phi_p, P_phi_a, sP_phi;


///////////////LYBRARIES////////////////////////
#include <io_lib.h>
#include <initialize.h>
#include <observables.h>
#include <rk4.h>
#include <sources.h>
#include <momentum.h>
#include <aux_functions.h>
#include <thetadot.h>

