%Example of using MUMPS in matlab with sparse right hansd side

% initialization of a matlab MUMPS structure
id = initmumps;
id = dmumps(id);
load lhr01;
mat = Problem.A;
% JOB = 6 means analysis+facto+solve
id.JOB = 6;
% we set the rigth hand side
id.RHS = ones(size(mat,1),2);
id.RHS(:,2) = 2*id.RHS(:,2);
id.RHS = sparse(id.RHS);
%call to mumps
id = dmumps(id,mat);
if(norm(mat*id.SOL - id.RHS,'inf') > sqrt(eps))
	disp('WARNING : precision may not be OK');
else
	disp('SOLUTION OK');
end
norm(mat*id.SOL - id.RHS,'inf')
% solution OK
% destroy mumps instance
id.JOB = -2;
id = dmumps(id)



