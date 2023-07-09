      SUBROUTINE trussf(eleObj,modl,tang,resid,isw,error) 
      
!DEC$ IF DEFINED (_DLL)
!DEC$ ATTRIBUTES DLLEXPORT :: TRUSSF
!DEC$ END IF

      use elementTypes
      use elementAPI
      implicit none
      
      type(eleObject)::eleObj
      type(modelState)::modl
      double precision tang(4, *)
      double precision resid(1)
      integer::isw;
      integer::error;

      integer :: tag, nd1, nd2, matTag, numCrd, i, j, numDOF
      real *8, pointer::theParam(:)
      integer, pointer::theNodes(:)

      double precision A, dx, dy, L, cs, sn
      double precision dLength, force, k

      integer :: iData(3);
      integer :: matTags(2);
c      integer (c_int), target :: matTags(2);
      
      type(c_ptr) :: theCMatPtr
      type(c_ptr), pointer :: theCMatPtrPtr(:)
      type(matObject), pointer :: theMat

      double precision dData(1), nd1Crd(2), nd2Crd(2)
      double precision d1(2), d2(2), tran(4)
      double precision strs(1), strn(1), tng(1)
      
      integer numData, err, matType

c     outside functions called
c      integer OPS_GetIntInput, OPS_GetDoubleInput, OPS_InvokeMaterial
c      integer OPS_GetNodeCrd, OPS_AllocateElement, OPS_GetNodeDisp

      IF (isw.eq.ISW_INIT) THEN
         
c     get the input data  - tag? nd1? nd2? A? matTag?

         numData = 3
         err = OPS_GetIntInput(numData, iData)
         tag = iData(1);
         nd1 = iData(2);
         nd2 = iData(3);

         numData = 1
         err = OPS_GetDoubleInput(numData, dData)
         A = dData(1);

         numData = 1
         err = OPS_GetIntInput(numData, iData)
         matTag = iData(1);

c     Allocate the element state 

         eleObj%tag = tag
         eleObj%nnode = 2
         eleObj%ndof = 4
         eleObj%nparam = 4
         eleObj%nstate = 0  
         eleObj%nmat = 1

         matTags(1) = matTag;
         matType = OPS_UNIAXIAL_MATERIAL_TYPE;
         err = OPS_AllocateElement(eleObj, matTags, matType)

c         theCMatPtr = theCMatPtrPtr(2); 
c         j=OPS_InvokeMaterialDirectly(theCMatPtr, modl, strn, strs,
c     +        tng, isw)

c         element sets material functions
c         call c_f_pointer(eleObj%mats, theCMatPtrPtr, [1]);
c         theCMatPtrPtr(1) = theCMatPtr; 
         
c     Initialize the element properties

         call c_f_pointer(eleObj%param, theParam, [4]);
         call c_f_pointer(eleObj%node, theNodes, [2]);

         numCrd = 2;
         err = OPS_GetNodeCrd(nd1, numCrd, nd1Crd);
         err = OPS_GetNodeCrd(nd2, numCrd, nd2Crd);

         dx = nd2Crd(1)-nd1Crd(1);
         dy = nd2Crd(2)-nd1Crd(2);

         L = sqrt(dx*dx + dy*dy);

         if (L == 0.0) then
c            OPS_Error("Warning - truss element has zero length\n", 1);
            return;
         end if

         cs = dx/L;
         sn = dy/L;

         theParam(1) = A;
         theParam(2) = L;
         theParam(3) = cs;
         theParam(4) = sn;

         theNodes(1) = nd1;
         theNodes(2) = nd2;

      ELSE

         IF (isw == ISW_COMMIT) THEN

            call c_f_pointer(eleObj%mats, theCMatPtrPtr, [1]);
            theCMatPtr = theCMatPtrPtr(1); 

            j=OPS_InvokeMaterialDirectly(theCMatPtr, modl, strn, strs,
     +       tng, isw)
            
         ELSE IF (isw == ISW_REVERT_TO_START) THEN

            call c_f_pointer(eleObj%mats, theCMatPtrPtr, [1]);
            theCMatPtr = theCMatPtrPtr(1); 

            j=OPS_InvokeMaterialDirectly(theCMatPtr, modl, strn, strs,
     +       tng, isw)

         ELSE IF (isw == ISW_FORM_TANG_AND_RESID) THEN

            call c_f_pointer(eleObj%param, theParam, [4]);
            call c_f_pointer(eleObj%node, theNodes, [2]);
            call c_f_pointer(eleObj%mats, theCMatPtrPtr, [1]);
            theCMatPtr = theCMatPtrPtr(1); 

            A = theParam(1);
            L = theParam(2);
            cs = theParam(3);
            sn = theParam(4);
            nd1 = theNodes(1);
            nd2 = theNodes(2);

            numDOF = 2;
            err = OPS_GetNodeDisp(nd1, numDOF, d1);
            err = OPS_GetNodeDisp(nd2, numDOF, d2);    

            tran(1) = -cs;
            tran(2) = -sn;
            tran(3) = cs;
            tran(4) = sn;
            
            dLength = 0.0;
            do 10 i = 1,2
               dLength = dLength - (d2(i)-d1(i)) * tran(i);
 10         continue

            strn(1) = dLength/L;

c            i = 0
c            i=OPS_InvokeMaterial(eleObj, i, modl, strn, strs, tng, isw)
            j=OPS_InvokeMaterialDirectly(theCMatPtr, modl, strn, strs,
     +       tng, isw)

            force = A*strs(1);
            k = A*tng(1)/L;

            do 20 i =1,4
               resid(i) = force * tran(i);
               do 30 j = 1,4
                  tang(i,j) = k * tran(i)*tran(j);
 30            continue
 20         continue

         END IF

      END IF

c     return error code
      error = 0

      END SUBROUTINE trussf


      SUBROUTINE localinit() 
      
!DEC$ IF DEFINED (_DLL)
!DEC$ ATTRIBUTES DLLEXPORT :: LOCALINIT
!DEC$ END IF

      use elementAPI
      implicit none
      integer::error;
      character *60:: msg;
      msg = 'trussf - Written by fmk UC Berkeley Copyright 2008'
      error = OPS_Error(msg);
      END SUBROUTINE localInit

