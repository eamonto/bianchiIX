/*
=========================================================================
io_lib.c
=========================================================================

Contains:
-> input-output routines.
-> create and remove files routines.
-> Usage message routine.

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


/////REMOVE ARCHIVE outputfile AND CREATE IT AGAIN
int create_remove_dir(void)
{
  char remove[200],copy[200],create[200];
  char *rm="rm -rf ",*cp="cp ",*mkdir="mkdir -p ";
  int out_return;

  strcpy(remove,rm);
  strcat(remove,outputfile);
  out_return = system(remove);

  strcpy(create,mkdir);
  strcat(create,outputfile);
  out_return = system(create);

  strcpy(copy,cp);
  strcat(copy,initfile);
  strcat(copy," ");
  strcat(copy,outputfile);
  out_return = system(copy);

  return out_return;
}


///////////VALIDATION OF INPUT FILES/////////////
int usage(void)
{
  printf("\n");
  printf("USAGE: ./exec <initfile> <outputfile> \n");
  printf("\n");
  exit(1);
}

///////////VALIDATE IF THE INITIAL CONDITIONS ARE CORRECT//////
int validate_initial_data(void)
{
  if (p1<=0 || p2<=0 || p3<=0){
    
    printf("\n");
    printf("\t p1,p2 and p3 need to be positive. \n");
    printf("\n");
    exit(1);
  }
  
  if ( mu1c1 < -M_PI/2.0 || mu2c2 < -M_PI/2.0 || mu3c3 < -M_PI/2.0 ){
    
    printf("\n");
    printf("\t Warning!! mu1c1, mu2c2, mu3c3 are not greater than -pi/2 \n");
    printf("\n");
    //    exit(1);
  }

  if ( mu1c1 > 3.0*M_PI/2.0 || mu2c2 > 3.0*M_PI/2.0 || mu3c3 > 3.0*M_PI/2.0 ){
    
    printf("\n");
    printf("\t Warning!! mu1c1, mu2c2, mu3c3 are not smaller than 3pi/2 \n");
    printf("\n");
    //    exit(1);
  }

  if (p2 == p3){
    if(p2 < 0.44789*pow(p1,1.34258) || p2 > 3.48989*pow(p1,0.58310)){
      
      printf("\n");
      printf("\t Warning!! You can be outside the physical region. \n");
      printf("\n");
    }
  }

  return 0;
}

/// IO ROUTINES

////////READ THE PARAMETERS/////
int read_param(void)
{
  int log_pf;
  FILE *pf;

  pf=fopen(initfile,"r");
  
  log_pf=fscanf(pf,"%lf%lf%lf%lf%lf%lf %lf%lf%lf %d %d%d %d%d %d",
		&mu1c1,&mu2c2,&mu3c3,&p1,&p2,&p3,
		&initial_time,&final_time,&dt,
		&time_output,
		&bianchi_switch,&class_effec_switch,
		&Ed_A_BS_switch,&effec_effec_switch,
		&std_out
		);

  fclose(pf);

  if (log_pf != 15) {
    printf("\n\t Error reading the parameters file.");
    printf("\n\t The program read %d initial values and they need to be 15. \n",log_pf);
    exit(1);  
  }

  return 0;
}


//// WRITE THE OBSERVABLES TO A FILE ////
int write_output(void)
{
  double spatial_Ricci;
  int log_pf;
  char newfile[200],*aux;  
  FILE *pf;

  if(density   < 1.0e-15)  density   = 0.0;
  if(shear     < 1.0e-15)  shear     = 0.0;
  if(fabs(Ricci)      < 1.0e-15)  Ricci     = 0.0;
  if(fabs(expansion)  < 1.0e-15)  expansion = 0.0;
  if(fabs(constraint) < 1.0e-15) constraint = 0.0;

  aux="/density.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,density/density_crit);
  fclose(pf);

  aux="/shear.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,shear);
  fclose(pf);

  aux="/volume.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,vol);
  fclose(pf);

  aux="/constraint.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,constraint);
  fclose(pf);

  aux="/expansion.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,expansion);
  fclose(pf);

  aux="/a123.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e \t %.15e \n",run_time,a_prom,a1,a2,a3);
  fclose(pf);

  aux="/H123.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e\t %.15e\t %.15e \n",run_time,H1,H2,H3);
  fclose(pf);


  spatial_Ricci = -0.5*(x1*x1 + x2*x2 + x3*x3 - 2.0*(x1*x2 + x1*x3+ x2*x3));

  aux="/curvature.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e  \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e\n"
		 ,run_time,Ricci,spatial_Ricci,curvature_param,x1,x2,x3);
  fclose(pf);

  aux="/cp.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t%.15e\t%.15e\t%.15e \t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e \n",
		 run_time,c1,c2,c3,p1,p2,p3,mu1*c1,mu2*c2,mu3*c3);
  fclose(pf);

  aux="/kasner.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e \t %.15e \n",run_time,k1,k2,k3,k1+k2+k3);
  fclose(pf);

  aux="/constants.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t%.15e \t%.15e \t%.15e \t%.15e\n"
		 ,run_time,field_momentum(),omega,sigma2,curvature_param);
  fclose(pf);

  return 0;
}
