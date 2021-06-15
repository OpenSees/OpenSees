//A simple demo for the MUMPS interface 
//to run it, You just have to execute  the instruction within Scilab
//		exec double_example.sce; 


//*********************** MATRIX INITIALISATION ***********************//
// This matrix has to be a SciSparse, otherwise it won't work.
	exec('ex.sci');
	mat=sparse(a);

// Right Hand side setting
	RHS = ones(size(mat,1),1);
	
//****************** Initialisation of the Scilab MUMPS structure ******************//
timer();
[id]=initmumps();

//Here Job=-1, the next call will only initialise the C and Fortran structure
[id]=dmumps(id);
id.RHS=RHS;

//******************** CALL TO MUMPS FOR RESOLUTION ********************//
job=6;
id.JOB=job;

[id]=dmumps(id,mat);

// verification of the solution
solution=id.SOL;
norm_res=norm(mat*solution-RHS,'inf');
write(%io(2),norm_res);

//****************** DESTRUCTION OF THE MUMPS INSTANCE ******************//
job=-2;
id.JOB=job;
[id]=dmumps(id);
t=timer()
