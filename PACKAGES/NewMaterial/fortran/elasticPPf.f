      SUBROUTINE ELASTICPPF(matObj,model,strain,tang,stress,isw,error) 
      
!DEC$ IF DEFINED (_DLL)
!DEC$ ATTRIBUTES DLLEXPORT :: ELASTICPPF
!DEC$ END IF
      use materialTypes
      use materialAPI
      implicit none
      
      type(matObject)::matObj
      type(modelState)::model
      real *8::strain 
      real *8::tang
      real *8::stress
      integer::isw
      integer::error;

      real *8, pointer::theParam(:)
      real *8, pointer::cState(:)
      real *8, pointer::tState(:)

      real *8 E, eyp, ep, fyp, fyn, f, fYieldSurface
      real *8 sigtrial, trialStrain, trialStress, trialTangent
      
      integer, target :: iData(2)
      integer, pointer ::iPtr(:)
      integer numData, err
      real *8, target :: dData(2)
      real *8, pointer:: dPtr(:)

c     outside functions called
c      integer OPS_GetIntInput, OPS_GetDoubleInput, OPS_AllocateMaterial

      IF (isw.eq.ISW_INIT) THEN

c     get the input data  - tag? E? eyp? 

         numData = 1
         iPtr=>iData;
         err = OPS_GetIntInput(numData, iPtr)
         numData = 2
         dPtr=>dData;
         err = OPS_GetDoubleInput(numData, dPtr)
         
c     Allocate the element state 
         matObj%tag = idata(1)
         matObj%nparam = 2
         matObj%nstate = 2  
         
         err = OPS_AllocateMaterial(matObj)
         call c_f_pointer(matObj%theParam, theParam, [2]);
         
c     Initialize the element properties
         theParam(1) = dData(1);
         theParam(2) = dData(2);

      ELSE
         
         call c_f_pointer(matObj%theParam, theParam, [2]);
         call c_f_pointer(matObj%cState, cState, [2]);         
         call c_f_pointer(matObj%tState, tState, [2]);         

         IF (isw == ISW_COMMIT) THEN
            
            trialStrain = tState(1)
            
            E = theParam(1)
            eyp = theParam(2)
            ep = cState(2)
            
            fyp = E*eyp
            fyn = -fyp
            
c     compute trial stress
            sigtrial = E * ( trialStrain - ep );
            
c     evaluate yield function
            IF ( sigtrial >= 0.0 ) THEN
               f =  sigtrial - fyp
            ELSE
               f = -sigtrial + fyn
            END IF
            
            fYieldSurface = - E * DBL_EPSILON;
            
            IF ( f > fYieldSurface ) THEN
c     plastic
               IF ( sigtrial > 0.0 ) THEN
                  ep = ep + f / E
               ELSE 
                  ep = ep - f / E
               END IF
            END IF
            
            cState(1) = trialStrain;    
            cState(2) = ep;
            
         ELSE IF (isw == ISW_REVERT_TO_START) THEN
            cState(1) = 0.0;
            cState(2) = 0.0;
            tState(1) = 0.0;
            tState(2) = 0.0;

         ELSE IF (isw == ISW_FORM_TANG_AND_RESID) THEN
            
            trialStrain = strain;
            
            E = theParam(1);
            eyp = theParam(2);
            ep = cState(2);

            fyp = E*eyp;
            fyn = -fyp;
            
c     compute trial stress
            sigtrial = E * ( trialStrain - ep );

c     evaluate yield function
            IF ( sigtrial >= 0.0 ) THEN
               f =  sigtrial - fyp;
            ELSE
               f = -sigtrial + fyn;
            END IF
            
            fYieldSurface = - E * DBL_EPSILON;
            
            IF ( f <= fYieldSurface ) THEN
               
               trialStress = sigtrial;
               trialTangent = E;
               
            ELSE
               
c     plastic
               IF ( sigtrial > 0.0 ) THEN
                  trialStress = fyp;
               ELSE
                  trialStress = fyn;
               END IF
               
               trialTangent = 0.0;
            END IF
            
            tState(1) = trialStrain;
            stress = trialStress;
            tang = trialTangent;
         END IF
      END IF

c     return error code
      error = 0

      END SUBROUTINE elasticPPf


      SUBROUTINE LOCALINIT() 
      
!DEC$ IF DEFINED (_DLL)
!DEC$ ATTRIBUTES DLLEXPORT :: LOCALINIT
!DEC$ END IF

      use materialAPI
      implicit none
      integer::error;
      character *60:: msg;
      msg = 'elasticPPf - Written by fmk UC Berkeley Copyright 2008'
      error = OPS_Error(msg);
      END SUBROUTINE localInit
