
       FUNCTION passive_xc(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: VarIn 
       REAL(kind=dp) :: VarOut 
       LOGICAL,SAVE :: FirstTime=.TRUE.
       REAL(kind=dp),SAVE :: xc

       IF (FirstTime) THEN
         xc = ListGetConstReal( Model % Constants, 'x_crit', UnFoundFatal=.TRUE. )
         FirstTime=.FALSE.
       ENDIF

       IF (VarIn.GT.xc) THEN
        VarOut=+1.0
       ELSE
        VarOut=-1.0_dp
      END IF

       End FUNCTION passive_xc

       FUNCTION passive_xct(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: VarIn 
       REAL(kind=dp) :: VarOut 
       LOGICAL,SAVE :: FirstTime=.TRUE.
       REAL(kind=dp),SAVE :: xc0,rate
       REAL(kind=dp) :: xc,Time

       IF (FirstTime) THEN
         xc0 = ListGetConstReal( Model % Constants, 'x_crit', UnFoundFatal=.TRUE. )
         rate = ListGetConstReal( Model % Constants, 'x_crit rate', UnFoundFatal=.TRUE. )
         FirstTime=.FALSE.
       ENDIF

       Time = GetTime()
       xc=xc0+rate*Time
       IF (VarIn.GT.xc) THEN
        VarOut=+1.0
       ELSE
        VarOut=-1.0_dp
      END IF

       End FUNCTION passive_xct


       FUNCTION passive_mask(Model,nodenumber,VarIn) RESULT(VarOut)
        USE DefUtils
        implicit none
        !-----------------
        TYPE(Model_t) :: Model
        INTEGER :: nodenumber
        REAL(kind=dp),INTENT(IN) :: VarIn(5) ! Time,distance1p,distance2p,distance1n,distance2n
        REAL(kind=dp) :: VarOut, norm_ad, norm_re, norm_t
 
        

        IF (VarIn(2).GT.0.01.AND.VarIn(3).GT.0.01) THEN
         VarOut=-1.0_dp
        ELSE IF (VarIn(4).GT.0.01.AND.VarIn(5).GT.0.01) THEN
         VarOut=+1.0
        ELSE
          norm_re = VarIn(2)/(VarIn(2)+VarIn(5))
          norm_ad = VarIn(4)/(VarIn(4)+VarIn(3))
          norm_t = VarIn(1)-floor(VarIn(1))
          IF (norm_re.GT.norm_t.OR.norm_re.LT.norm_t.AND.norm_re.GT.0.01) THEN
            VarOut=-1.0_dp
          ELSE
            VarOut=1.0
          END IF
        END IF
 
       End FUNCTION passive_mask

       FUNCTION passive_mask_old(Model,nodenumber,VarIn) RESULT(VarOut)
        USE DefUtils
        implicit none
        !-----------------
        TYPE(Model_t) :: Model
        INTEGER :: nodenumber
        REAL(kind=dp),INTENT(IN) :: VarIn(2) ! mask_nodes, zs
        REAL(kind=dp) :: VarOut 
        LOGICAL,SAVE :: FirstTime=.TRUE.
        REAL(kind=dp),SAVE :: xc0,rate
        REAL(kind=dp) :: xc,Time
 

        IF (VarIn(1).GT.0.5) THEN
         VarOut=-1.0_dp
        ELSE
         VarOut=+1.0
       END IF
 
        End FUNCTION passive_mask_old

        FUNCTION passive_mask_simple(Model,nodenumber,VarIn) RESULT(VarOut)
          USE DefUtils
          implicit none
          !-----------------
          TYPE(Model_t) :: Model
          INTEGER :: nodenumber
          REAL(kind=dp),INTENT(IN) :: VarIn ! mask_nodes, zs
          REAL(kind=dp) :: VarOut 
          LOGICAL,SAVE :: FirstTime=.TRUE.
          REAL(kind=dp),SAVE :: xc0,rate
          REAL(kind=dp) :: xc,Time
   
  
          IF (VarIn.GT.0.5) THEN
           VarOut=-1.0_dp
          ELSE
           VarOut=+1.0
         END IF
   
          End FUNCTION passive_mask_simple
  
        FUNCTION passive_mask_norm(Model,nodenumber,VarIn) RESULT(VarOut)
          USE DefUtils
          implicit none
          !-----------------
          TYPE(Model_t) :: Model
          INTEGER :: nodenumber
          REAL(kind=dp),INTENT(IN) :: VarIn(5) ! Time,distance1p,distance2p,distance1n,distance2n
          REAL(kind=dp) :: VarOut, norm_ad, norm_re, norm_t
   
          
  
          IF (VarIn(2).GT.0.01.AND.VarIn(3).GT.0.01) THEN
           VarOut=-1.0_dp
          ELSE IF (VarIn(4).GT.0.01.AND.VarIn(5).GT.0.01) THEN
           VarOut=+1.0
          ELSE
            norm_re = VarIn(2)/(VarIn(2)+VarIn(5))
            norm_ad = VarIn(4)/(VarIn(4)+VarIn(3))
            norm_t = VarIn(1)-floor(VarIn(1))
            VarOut= norm_re/norm_t
          END IF
   
         End FUNCTION passive_mask_norm

        FUNCTION INV_MASK(Model,nodenumber,VarIn) RESULT(VarOut)
          USE DefUtils
          implicit none
          !-----------------
          TYPE(Model_t) :: Model
          INTEGER :: nodenumber
          REAL(kind=dp),INTENT(IN) :: VarIn ! Time,distance1p,distance2p,distance1n,distance2n
          REAL(kind=dp) :: VarOut
   
          VarOut = -VarIn+1
         
   
          End FUNCTION INV_MASK
  