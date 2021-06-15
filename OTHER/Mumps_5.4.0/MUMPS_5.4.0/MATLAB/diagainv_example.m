%Example of using MUMPS in matlab to compute diagonal of inverse of A

% Change to true to test example complex arithmetic
complex_arithmetic = false;

% initialization of a matlab MUMPS structure
id = initmumps;
if (complex_arithmetic)
  id = zmumps(id);
else
  id = dmumps(id);
end
load lhr01;
mat = Problem.A;

if (complex_arithmetic)
  % To test complex version
  mat = mat + i * speye(size(mat,1),size(mat,1)); 
end

% JOB = 4 means analysis+factorization
id.JOB = 4;
if (complex_arithmetic)
  id = zmumps(id,mat);
else
  id = dmumps(id,mat);
end

% Set the right hand side structure to requested entries of A-1
id.RHS = speye(size(mat,1),size(mat,1)); % Sparse format required

%call MUMPS solution phase to compute diagonal entries of A-1
id.ICNTL(30)=1; % Ask for A-1 entries
id.JOB=3;


if (complex_arithmetic)
  id = zmumps(id,mat);
else
  id = dmumps(id,mat);
end

% diagonal values have been computed in
% the (sparse) matrix id.SOL, which has
% the same structure as id.RHS

% Compare diagonal of inverse computed by Mumps and by matlab
disp(' ');
disp('Computing 2-norm of error on diagonal of inverse:');
norm(diag( diag(diag(inv(mat)))-id.SOL ),2)

% destroy mumps instance
id.JOB = -2;
if (complex_arithmetic)
  id = zmumps(id)
else
  id = dmumps(id)
end
