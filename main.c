/*
========================================================================
main.c
========================================================================
Principal program. 
This program solves the Bianchi IX equations, classicals and effectives.


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


int main (int argc,char **argv)
{
  double Ttotal;
  int l;

  if(argc != 3) usage();        //VALIDATION OF INPUT FILES

  initfile   = argv[1];         //ASIGNATION OF FILE'S 
  outputfile = argv[2];         //NAMES
  
  read_param();                 //READ INITIAL CONDITIONS

  create_remove_dir();          //REMOVE OLD outputfile AND CREATE A NEW ONE

  initialize_all();             //INITIALIZATION OF ALL QUANTITIES

  validate_initial_data();      //VALIDATED IF THE INITIAL CONDITIONS ARE CORRECT

  if(std_out==1){
    printf("\n \t Start Simulation \n");
    printf("\t Time -> %4.6lf\n",run_time);
  }

  Ttotal = final_time-initial_time;

  l=0;

  do
    {

      if(l%time_output==0)
	{
	  compute_obs();          //COMPUTE OBSERVABLES
	  write_output();         //WRITE THE OBSERVABLES TO A FILE
	  if(std_out==1)
	    printf("\t Time -> %4.6lf\n",run_time);
	}

      //RK4 INTEGRATOR
      store_levels_rk4();
      
      rhs();
      evolution_rk4(1);
      rhs();
      evolution_rk4(2);
      rhs();
      evolution_rk4(3);
      rhs();
      evolution_rk4(4);
      
      run_time = run_time + dt;

      l++;
      
    }while(fabs(Ttotal)>fabs(initial_time-run_time));

  if(std_out==1)  printf("\n \t Finish Simulation \n");
  
  return 0;
}
