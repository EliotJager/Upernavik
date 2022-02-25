
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
        REAL(kind=dp),INTENT(IN) :: VarIn(2) ! mask_nodes, zs
        REAL(kind=dp) :: VarOut 
        LOGICAL,SAVE :: FirstTime=.TRUE.
        REAL(kind=dp),SAVE :: xc0,rate
        REAL(kind=dp) :: xc,Time
 

        IF (VarIn(1).GT.0.5.OR.VarIn(2).GT.10.0) THEN
         VarOut=-1.0_dp
        ELSE
         VarOut=+1.0
       END IF
 
        End FUNCTION passive_mask
