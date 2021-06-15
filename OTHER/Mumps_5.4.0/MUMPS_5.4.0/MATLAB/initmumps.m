function id = initmumps()
%
% id = initmumps
% it returns a default matlab MUMPS structure
%
% Use help mumps_help for detailed information
%
errmsg = nargoutchk(1,1,nargout);
if(~isempty(errmsg))
     disp(errmsg);
     return;
end
id = struct('SYM',0,'JOB',-1,'ICNTL',zeros(1,60)-9998,'CNTL',zeros(1,15)-9998,'PERM_IN',-9999,'COLSCA',-9999,'ROWSCA',-9999,'RHS',-9999,'INFOG',zeros(1,80)-9998,'RINFOG',zeros(1,40)-9998,'VAR_SCHUR',-9999,'SCHUR',-9999,'INST',-9999,'SOL',-9999,'REDRHS',-9999,'PIVNUL_LIST',-9999,'MAPPING',-9999,'SYM_PERM',-9999,'UNS_PERM',-9999,'TYPE',0,'KEEP',zeros(1,500)-9998,'DKEEP',zeros(1,230)-9998);
