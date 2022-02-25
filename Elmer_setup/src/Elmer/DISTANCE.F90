!########################################################################      
!### A user function to compute new SMB
!#######################################################################
      FUNCTION DISTANCE_ANOM(Model,nodenumber,VarIn)  RESULT(VarOut)
       USE Types
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(5) !Time,distance1,distance2
       REAL(kind=dp) :: VarOut
       !-----------------
       !REAL :: norm
       !norm=MAX(1.0,((VarIn(3)+VarIn(2))/2))
       
       VarOut=(VarIn(2)+(VarIn(1)-floor(VarIn(1)))*(VarIn(3)-VarIn(2)))
       
      END FUNCTION DISTANCE_ANOM

      FUNCTION DISTANCE_DIFF(Model,nodenumber,VarIn)  RESULT(VarOut)
            USE Types
            implicit none
            !-----------------
            TYPE(Model_t) :: Model
            INTEGER :: nodenumber
            REAL(kind=dp) :: VarIn(2) !distancep,distancen
            REAL(kind=dp) :: VarOut
            !-----------------
            !REAL :: norm
            !norm=MAX(1.0,((VarIn(3)+VarIn(2))/2))
            
            VarOut= VarIn(1)-VarIn(2)
            
      END FUNCTION DISTANCE_DIFF
      
      FUNCTION DistanceCond(Model,nodenumber,VarIn) RESULT(VarOut)
            USE DefUtils
            implicit none
            !-----------------
            TYPE(Model_t) :: Model
            INTEGER :: nodenumber
            REAL(kind=dp) :: VarIn,VarOut
      
            IF (VarIn.LT.0.5) THEN
              VarOut = +1.0
            ELSE
              VarOut = -1.0
            END IF
      End FUNCTION DistanceCond
!#####################################################################
