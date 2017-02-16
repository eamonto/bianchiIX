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
double fsin1(double p3_f);//Auxiliar Sine function to find p3 from the constraint
double fsin2(double p3_f);//Auxiliar Sine function to find p3 from the constraint
double fsin3(double p3_f);//Auxiliar Sine function to find p3 from the constraint
double func_grav_constraint(double p3_f);//Gravitational Constraint as function of p3
double find_max_vol(void);
