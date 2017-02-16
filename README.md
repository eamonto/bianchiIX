README
=================================================================
Program to solve the Effective Loop Quantum Cosmology Equations 
for Bianchi IX.

** Written by Edison Montoya (c) 2017**

** Up to date: Feb 16 2017**

In this file we specify how to use the program and we explain
the file's content and the functions.

To make the program works it is needed:

1. C compiler, example: gcc
2. make -GNU make utility to maintain groups of programs
3. Python and the package "Tkinter"


**The easier way to work with the program it is to execute the command**

    $make

This compiles the program and open a window where can be changed the
program parameters. Next you press the button "calculate" and the
program runs. You can see the results in the "output" file.
Note: To do this you need to have installed python and the package "Tkinter"


The **COMPILATION** is made with gcc, execute the command

    $make exec.out

It generates the archive exec.out


The **EXECUTION** is like this

    $./exec.out <parameterfile> <outputfile>


| Intput files           | Description                    |
| ---------------------- | ------------------------------ |
| `exec.out`             | Executable                     |
| `parametersfile`       | File with the parameters of the program (default param.txt) |
| `outputfile`           | Name of the directory where the results will be save        |


An easy way to execute the program is using the command

    $make run

It runs the program "exec.out". The user can
specify some variables into the Makefile:

    PARAMFILE = param.txt
    OUTPUTFILE = output

These are the default values.


The **OUTPUT** of the program is save into the output directory 
specified in the Makefile. The files are:


## Output Files  
--------------------

**a123.t**, Contains the evolution of the scale factors, it is written as:

    time, (a1a2a3)^1/3, a1, a2, a3

**cp.t**, Contains the evolution of the connections, it is written as:

    time, c1, c2, c3, p1, p2, p3 , mu1c1, mu2c2, mu3c3

**H123.t**, Contains the evolution of the Hubble parameters, it is written as:

    time, H1, H2, H3

**kasner.t**, Contains the evolution of the Kasner parameters, it is written as:

    time, k1, k2, k3, k1+k2+k3

**expansion.t**, Evolution of the expansion, it is written as:

    time, expansion

**constraint.t**, Evolution of the constraint relative error, it is written as:

    time, relative error of the constraint
 
**density.t**, Evolution of the density, it is written as:

    time, density
  
**volume.t**, Contains the evolution of the volume, it is written as:

    time, volume

**shear.t**, Contains the evolution of the shear, it is written as:

    time, shear

**curvature.t**, Contains the evolution of the Ricci scalar, it is written as:

    time, Ricci scalar, spatial Ricci, curvature parameter,x1 ,x2 ,x3

**constants.t**, Contains the evolution of the constants of motion, it is written as:

    time, P_phi, density parameter, shear parameter, curvature parameter

**param.f90**, It is a copy of the file with the parameters.


------------------------------------------------------------------------

** Now it is describe the files contained in this directory and the role
that they play in the solution of the problem. **

param.txt
------------------------------------------------------------------------
This file contains the parameters of the program:

Parameter     | Description
------------- | --------------
1.0           | mu1c1
1.0           | mu2c2
1.0           | mu3c3
100.0         | p1
100.0         | p2
100.0         | p3
0.0           | Initial Time
1.0           | Final Time
1.0e-4        | Time Step
10            | Number of time steps for write the output
0             | 0 = Bianchi I,  1 = Bianchi IX
1             | 0 = Classical,  1 = Effective
0             | Effective Equations; 0 = Edward,  1 = Asieh
0             | Asieh Equations; 0 = Reduce to FRW,  1 = Not Reduce to FRW
0             | Standart output; 0 = Off,  1 = On


header.h
------------------------------------------------------------------------
This is the head of the program, here the global variables 
are declared. The program was developed with a RK4 integrator.


------------------------------------------------------------------------
main.c
------------------------------------------------------------------------
Principal routine, it coordinates all the program in order to solve
the problem.


------------------------------------------------------------------------
initialize.c
------------------------------------------------------------------------
Initialization of all quantities. The routines are:

```C
// INITIALIZATION OF ALL QUANTITIES
int initialize_all(void);
```

------------------------------------------------------------------------
io_lib.c
------------------------------------------------------------------------
Contains:

* The input-output routines.
* The create and remove files routines.
* Usage message routine.

The routines are:

```C
///// REMOVE ARCHIVE outputfile AND CREATE IT AGAIN
int create_remove_dir(void);

///// VALIDATION OF INPUT FILES
int usage(void);

//// VALIDATED IF THE INITIAL CONDITIONS ARE CORRECT
int validate_intial_data(void);

///// READ THE PARAMETERS
int read_param(void);

///// WRITE THE OBSERVABLES TO A FILE
int write_output(void);
```



------------------------------------------------------------------------
rk4.c
------------------------------------------------------------------------
Implementation of the Runge-Kutta 4 method. The routines are:

```C
///// STORE THE CONDITIONS BEFORE THE INTEGRATION
int store_levels_rk4(void);

///// EVOLUTION WITH THE RK4 METHOD
int evolution_rk4(int k);
```


------------------------------------------------------------------------
observables.c
------------------------------------------------------------------------
Here are compute all the quantities that can give relevant 
information about the system dynamics. The routines are:

```C
///// COMPUTE OBSERVABLES
int compute_obs(void);
```


------------------------------------------------------------------------
momentum.c
------------------------------------------------------------------------
Here is compute the momentum for each theory. The routines are:

```C
////// SCALAR FIELD MOMENTUM
double field_momentum(void);

///// CLASSICAL FIELD MOMENTUM
double class_bianchi_P_phi(void);

///// MOMENTUM FOR EDWARD'S QUANTIZATION
double Edward_P_phi(void);

///// MOMENTUM FOR ASIEH'S QUANTIZATION THAT NOT REDUCES TO FRW
double Asieh_notFRW_P_phi(void);

///// MOMENTUM FOR ASIEH'S QUANTIZATION THAT REDUCES TO FRW
double Asieh_FRW_P_phi(void);
```


------------------------------------------------------------------------
aux_functions.c
------------------------------------------------------------------------
Here are compute the auxiliar functions for the theories with inverse
triad corrections. The routines are:

```C
double potential(void);           // SCALAR FIELD POTENTIAL
double dev_pot(void);             // DERIVATIVE OF THE SCALAR FIELD POTENTIAL
double func_A(double v);          // FUNCTION A, ASIEH
double func_f(double v);          // FUNCTION f, ASIEH
double func_g(double v);          // FUNCTION g, ASIEH
double func_h(double v);          // FUNCTION h, ASIEH
double dev_A(double v,int i);     // DERIVATIVE A FUNCTION
double dev_f(double v);           // DERIVATIVE f FUNCTION
double dev_g(double v);           // DERIVATIVE g FUNCTION
double dev_h(double v,int i);     // DERIVATIVE h FUNCTION
```


------------------------------------------------------------------------
sources.c
------------------------------------------------------------------------
Here are implemented the right hand side of the equations 
that we want to integrate. The routines are:

```C
int rhs(void);                // RIGHT HAND SIDE
int class_bianchi_rhs(void);  // CLASSICAL BIANCHI IX
int Edward_rhs(void);         // EDWARD EQUATIONS
int Asieh_notFRW_rhs(void);   // ASIEH EQUATIONS, NOT REDUCE TO FRW
int Asieh_FRW_rhs(void);      // ASIEH EQUATIONS, REDUCE TO FRW
```

------------------------------------------------------------------------
thetadot.c
------------------------------------------------------------------------
Here are implemented the routines to calculated the derivative of the 
expansion. The routines are:

```C
double class_thetadot(void);         // DERIVATIVE OF THE CLASSICAL EXPANSION
double Edward_thetadot(void);        // DERIVATIVE OF THE EXPANSION (EDWARD)
double Asieh_notFRW_thetadot(void);  // DERIVATIVE OF THE EXPANSION (ASIEH, NOT FRW)
double Asieh_FRW_thetadot(void);     // DERIVATIVE OF THE EXPANSION (ASIEH,FRW)
```