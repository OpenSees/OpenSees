% Simple example of using MUMPS in matlab

% initialization of a matlab MUMPS structure
id = initmumps;
id.SYM = 0;

% here JOB = -1, the call to MUMPS will initialize C 
% and fortran MUMPS structure
id = dmumps(id);
% load a sparse matrix
load lhr01;
mat = Problem.A;
% JOB = 6 means analysis+facto+solve

%prob = UFget(373);
%mat = prob.A;
id.JOB = 6;
%%%%%%% BEGIN OPTIONAL PART TO ILLUSTRATE THE USE OF MAXIMUM TRANSVERSAL
id.ICNTL(7) = 5;
id.ICNTL(6) = 1;
id.ICNTL(8) = 7;
id.ICNTL(14) = 80;
% we set the rigth hand side
id.RHS = ones(size(mat,1),1);
%call to mumps
id = dmumps(id,mat);
% we see that there is a memory problem in INFOG(1) and INFOG(2)
id.INFOG(1)
id.INFOG(2)
% we activate the numerical maximun transversal 
fprintf('total number of nonzeros in factors %d\n', id.INFOG(10));

%%%%%%% END OPTIONAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(norm(mat*id.SOL - ones(size(mat,1),1),'inf') > sqrt(eps))
	disp('WARNING : precision may not be OK');
else
	disp('SOLUTION OK');
end
norm(mat*id.SOL - ones(size(mat,1),1),'inf')
% destroy mumps instance
SOL = id.SOL;
id.JOB = -2;
id = dmumps(id)
