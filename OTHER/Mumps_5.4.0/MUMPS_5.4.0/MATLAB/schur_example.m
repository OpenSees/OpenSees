%Example of using MUMPS in matlab with schur option

% initialization of a matlab MUMPS structure
id = initmumps;
id = dmumps(id);
load lhr01;
mat = Problem.A;
themax = max(max(abs(mat)));
n = size(mat,1);
mat = mat+sparse(1:n,1:n,3*themax*ones(n,1));

% initialization of Schur option
id.VAR_SCHUR = [n-9:n];

% JOB = 6 means analysis+facto+solve
id.JOB = 6;
id.RHS = ones(size(mat,1),1);
%call to mumps
id = dmumps(id,mat);
disp('*** check solution restricted to mat(1:n-10,1:n-10)');
if(norm(mat(1:n-10,1:n-10)*id.SOL(1:n-10) - ones(n-10,1),'inf') > sqrt(eps))
	disp('WARNING : precision may not be OK');
else
	disp('SCHUR SOLUTION CHECK1 OK');
end
norm(mat(1:n-10,1:n-10)*id.SOL(1:n-10) - ones(n-10,1),'inf')


% we want to use Schur complement to solve 
% A * sol = rhs
% with sol = x   and rhs = rhs1
%            y             rhs2
%
% check that the complete solution verify
% y = S^(-1) * (rhs2 - A_{2,1} * A_{1,1}^(-1) * rhs1)
% and 
% x = A_{1,1}^(-1) * rhs1) - A_{1,2} * y
%
sol1 = id.SOL(1:n-10);
rhsy = ones(10,1)-mat(n-9:n,1:n-10)*sol1;

%%%%%%%%%%%%%%%%%%%
% TO CHANGE :
% usually the resolution below is replaced by an iterative scheme
y = id.SCHUR \ rhsy;
%%%%%%%%%%%%%%%%%%%%

rhsx = mat(1:n-10,n-9:n)*y;
id.JOB = 3;
id.RHS(1:n-10) = rhsx;
id = dmumps(id,mat);
rhsx = id.SOL(1:n-10);
x = sol1-rhsx;
sol = [x;y];
r = mat*sol - ones(n,1);
disp('*** check complete solution');
if( norm(r,'inf') > sqrt(eps))
	disp('WARNING : precision may not be OK');
else
	disp('SCHUR SOLUTION CHECK2 OK');
end
norm(r,'inf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NOW TRY REDUCED RHS FUNCTIONALITY 
%  (easier to use than previous
%   computations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id.JOB=3;
% Do forward solution step to obtain a reduced RHS
id.ICNTL(26)=1;
RHS=mat*ones(n,1);
id.RHS=RHS;
id = dmumps(id,mat);
% Solve the problem on the interface
id.REDRHS = id.SCHUR \ id.REDRHS;

% Do backward solution stage to expand the solution
id.ICNTL(26)=2;
id = dmumps(id,mat);
r = mat*id.SOL-RHS;
disp('*** check solution when REDRHS is used');
if( norm(r,'inf') > sqrt(eps))
	disp('WARNING : precision may not be OK');
else
	disp('SCHUR SOLUTION CHECK3 OK');
end
norm(r,'inf')

