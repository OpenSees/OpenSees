function id = initmumps()
//
// id = initmumps
// it returns a default Scilab MUMPS mlist (structure)
//

id = mlist(["StructMumps";"SYM";"JOB";"ICNTL";"CNTL";"PERM_IN";"COLSCA";"ROWSCA";"RHS";"INFOG";"RINFOG";"VAR_SCHUR";"SCHUR";"INST";"SOL";"REDRHS";"PIVNUL_LIST";"SYM_PERM";"UNS_PERM";"TYPE"],0,-1,zeros(1,60)-9998,zeros(1,15)-9998,-9999,-9999,-9999,-9999,zeros(1,80)-9998,zeros(1,40)-9998,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,0);

endfunction

