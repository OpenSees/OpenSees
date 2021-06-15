function printmumpsstat(id)
%
% printmumpsstat(id)
% print mumps info
%

disp(['After analysis : Estimated operations                ' num2str(id.RINFOG(1))]);
disp(['After analysis : Estimated space for factors         ' int2str(id.INFOG(3))]);
disp(['After analysis : Estimated integer space             ' int2str(id.INFOG(4))]);
disp(['After analysis : Estimated max front size            ' int2str(id.INFOG(5))]);
disp(['After analysis : Number of node in the tree          ' int2str(id.INFOG(6))]);
disp(['After analysis : Estimated total size (Mbytes)       ' int2str(id.INFOG(17))]);

disp(['After factorization : Assembly operations            ' num2str(id.RINFOG(2))]);
disp(['After factorization : Elimination operations         ' num2str(id.RINFOG(3))]);
disp(['After factorization : Real/Complex space to store LU         ' int2str(id.INFOG(9))]);
disp(['After factorization : Integer space to store LU                 ' int2str(id.INFOG(10))]);
disp(['After factorization : Largest front size             ' int2str(id.INFOG(11))]);
disp(['After factorization : Number of off-diagonal pivots  ' int2str(id.INFOG(12))]);
disp(['After factorization : Number of delayed pivots       ' int2str(id.INFOG(13))]);
disp(['After factorization : Number of memory compresses    ' int2str(id.INFOG(14))]);
disp(['After factorization : Total size needed (Mbytes)     ' int2str(id.INFOG(19))]);
