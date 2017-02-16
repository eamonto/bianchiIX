/////REMOVE ARCHIVE outputfile AND CREATE IT AGAIN
int create_remove_dir(void);

///////////VALIDATION OF INPUT FILES/////////////
int usage(void);

///////////VALIDATED IF THE INITIAL CONDITIONS ARE CORRECT//////
int validate_initial_data(void);

////////READ THE PARAMETERS/////
int read_param(void);

//// WRITE THE OBSERVABLES TO A FILE ////
int write_output(void);
