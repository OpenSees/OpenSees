//A simple demo for the MUMPS interface, with the return of the schur complement
//to run it, You just have to execute  the instruction within Scilab
//		exec sparse_example.sce; 


//*********************** MATRIX INITIALISATION ***********************//
	n=10;
	mat=sprand(n,n,.5)+speye(n,n);
        size_schur=3;

// Right Hand side setting
	RHS = ones(n,1);

	
//****************** Initialisation of the Scilab MUMPS structure ******************//
timer();
[id]=initmumps();

//Here Job=-1, the next call will only initialise the C and Fortran structure
[id]=dmumps(id);

id.RHS=RHS;
id.VAR_SCHUR = [n-size_schur+1:n];

//******************** CALL TO MUMPS FOR RESOLUTION ON INTERNAL PROBLEM ************//
job=6;
id.JOB=job;

[id]=dmumps(id,mat);

// verification of the solution
solution=id.SOL;
norm1=norm(mat(1:n-size_schur,1:n-size_schur)*solution(1:n-size_schur) - ones(n-size_schur,1),'inf');
if norm1> 10^(-9) then
	write(%io(2),'WARNING: solution on internal problem may not be OK');
else
	write(%io(2),'SOLUTION on internal problem ok');
end


//******************* TRY REDUCED RHS FUNCTIONALITY **************//
id.JOB=3;
id.ICNTL(26)=1;

// Forward
[id]=dmumps(id,mat);

// Solve the problem on the Schur complement
id.REDRHS=id.SCHUR \ id.REDRHS;

// and reinject it to MUMPS
id.ICNTL(26)=2;
[id]=dmumps(id,mat);
solution=id.SOL;
norm1=norm(mat*solution-RHS,'inf')
if norm1> 10^(-9) then
	write(%io(2),'WARNING: solution on complete problem may not be OK');
else
	write(%io(2),'SOLUTION on complete problem ok');
end



//****************** DESTRUCTION OF THE MUMPS INSTANCE ******************//
job=-2;
id.JOB=job;
[id]=dmumps(id);
t=timer()
