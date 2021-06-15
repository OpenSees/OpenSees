function id=zmumps(id,mat)

//************************************************************************************************************** 
// [id] = zmumps(id,mat)
// id is a structure (see details in initmumps.m and MUMPS documentation)
// mat is an optional parameter if the job id.job = -1 or -2
// mat is a square sparse matrix
// informations are return in id fields
// 
// *************************************************************************************************************


if (typeof(id) ~= "StructMumps") then
  disp("Error. Please call initmumps first.");
  return;
end
arithtype=1;

if id.JOB == -2 then   
     if id.INST==-9999 then
         disp('Error. Uninitialized instance. MUMPS should be called with JOB=-1 first.');
         return;
     end
     if id.TYPE ~= arithtype then
	disp('Error. You are trying to call z/d version on a d/z instance');
        return;
     end
     // call the C routine zmumpsc
 
     zmumpsc(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS);
     id = [];
     return;
end


if id.JOB == -1 then
	if id.INST~=-9999 then
		disp('Error. Already initialized instance.');
	return;
	end
	// call the C routine zmumpsc
	[inform,rinform,sol,inst,schu,redrhs,pivnul_list,sym_perm,uns_perm,icntl,cntl] = zmumpsc(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS);
        id.INFOG = inform;
	id.RINFOG = rinform;
	id.SOL = sol;
	id.INST = inst;
	id.SCHUR = schu;
	id.REDRHS = redrhs;
	id.PIVNUL_LIST = pivnul_list;
	id.SYM_PERM = sym_perm;
	id.UNS_PERM = uns_perm;
        id.TYPE=arithtype;
        id.ICNTL=icntl;
	id.CNTL=cntl;
 	clear inform rinform sol inst schu redrhs pivnul_list sym_perm uns_perm icntl cntl
	return;
	
end

if id.INST ==-9999 then
	disp('Uninitialized instance');
	return;
end 
// call the C routine zmumpsc

if id.TYPE ~= arithtype then
	disp('You are trying to call z/d version on a d/z instance');
end

[inform,rinform,sol,inst,schu,redrhs,pivnul_list,sym_perm,uns_perm,icntl,cntl] = zmumpsc(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS, mat);
id.INFOG = inform;
id.RINFOG = rinform;
id.SOL = sol;
id.INST = inst;
if (id.JOB == 2|id.JOB==4|id.JOB==6) then
	if id.SYM == 0 then
                id.SCHUR=schu';
        else 
        	id.SCHUR=triu(schu)+tril(schu',-1);
        end
end 
id.REDRHS = redrhs;
id.PIVNUL_LIST = pivnul_list;
id.SYM_PERM(sym_perm) = [1:size(mat,1)];
id.UNS_PERM = uns_perm;
id.ICNTL=icntl;
id.CNTL=cntl;
clear inform rinform sol inst schu redrhs pivnul_list sym_perm uns_perm icntl cntl

endfunction
