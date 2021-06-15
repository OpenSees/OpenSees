% Simple example of using MUMPS in matlab

% initialization of a matlab MUMPS structure
id = initmumps;
% here JOB = -1, the call to MUMPS will initialize C 
% and fortran MUMPS structure
id = zmumps(id);
% load a sparse matrix
load lhr01;
mat = Problem.A;
n = size(mat,1);
mat = mat + sparse(1:n,1:n,i*ones(n,1));
% JOB = 6 means analysis+facto+solve
id.JOB = 6;
id.ICNTL(6) = 0;
% we set the right hand side
id.RHS = ones(size(mat,1),1);
%call to mumps
id = zmumps(id,mat);
% we see that there is a memory problem in INFOG(1) and INFOG(2)
id.INFOG(1)
id.INFOG(2)
% we activate the numerical maximun transversal 
id.ICNTL(6) = 6;
id = zmumps(id,mat);
if(norm(mat*id.SOL - ones(size(mat,1),1),'inf') > sqrt(eps))
	disp('WARNING : precision may not be OK');
else
	disp('SOLUTION OK');
end
norm(mat*id.SOL - ones(size(mat,1),1),'inf')
% destroy mumps instance
id.JOB = -2;
id = zmumps(id)

