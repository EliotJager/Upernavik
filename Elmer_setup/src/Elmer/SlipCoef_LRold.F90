FUNCTION SlipCoef(Model,nodenumber,VarIn)  RESULT(VarOut)
  USE DefUtils
  implicit none

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(3)  !slc, SSAVelocity1, SSAVelocity2
  REAL(kind=dp) :: VarOut

  TYPE(ValueList_t), POINTER :: material
  REAL(kind=dp) :: un
  REAL(kind=dp) :: m
  REAL(kind=dp),Parameter :: umin=100._dp,TauMin=0.1_dp

  ! inquire SSA friction exponent from Material properties
  material => GetMaterial()
  IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('SlipCoef_USF', "No material found?")
  ENDIF
  m = ListGetConstReal( Material, 'SSA Friction Exponent',UnFoundFatal=.TRUE.)

  ! compute slipcoef for require exponenet from input values
  ! coresponding to linear friction
  un=sqrt(VarIn(2)*VarIn(2)+VarIn(3)*VarIn(3))

  IF (un.GT.1._dp) THEN
   VarOut=VarIn(1)*(un**(1._dp-m))
  ELSE
   ! velocity was smaller than 1 => no ice area
   ! assume we should have TauMin for umin
   VarOut=TauMin/(umin**m)
  ENDIF

 End FUNCTION SlipCoef
 
 FUNCTION SlipCoefv2(Model,nodenumber,VarIn)  RESULT(VarOut)
  USE DefUtils
  implicit none

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(4)  !slc,SSAVelocity1, SSAVelocity2, m
  REAL(kind=dp) :: VarOut

  REAL(kind=dp) :: un
  REAL(kind=dp),Parameter ::umin=100._dp,TauMin=0.1_dp
  
  ! compute slipcoef for require exponent from input values
  ! coresponding to linear friction
  un=sqrt(VarIn(2)*VarIn(2)+VarIn(3)*VarIn(3))

  IF (un.GT.1._dp) THEN
   VarOut=VarIn(1)*(un**(1._dp-VarIn(4)))
  ELSE
   ! velocity was smaller than 1 => no ice area
   ! assume we should have TauMin for umin
   VarOut=TauMin/(umin**VarIn(4))
  ENDIF
  
 End FUNCTION SlipCoefv2


 FUNCTION CalculTaub (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(3)    ! slc, SSAvelocity 1, SSAvelocity 2
   REAL (KIND=dp) :: VarOut      ! taub
   INTEGER :: nodenumber
   
   TYPE(ValueList_t), POINTER :: material
   REAL(kind=dp) :: un
   REAL(kind=dp) :: m

   ! inquire SSA friction exponent from Material properties
   material => GetMaterial()
   IF (.NOT. ASSOCIATED(material)) THEN
     CALL Fatal('SlipCoef_USF', "No material found?")
   ENDIF
   m = ListGetConstReal( Material, 'SSA Friction Exponent',UnFoundFatal=.TRUE.)  
   
   un=sqrt(VarIn(2)*VarIn(2)+VarIn(3)*VarIn(3))
   
   VarOut=VarIn(1)*(un**(m))
    
 END FUNCTION CalculTaub

 FUNCTION CalculTaubRC (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(3)    ! slc, SSAvelocity 1, SSAvelocity 2
   REAL (KIND=dp) :: VarOut      ! taub
   INTEGER :: nodenumber
   
   TYPE(ValueList_t), POINTER :: material
   REAL(kind=dp) :: un
   REAL(kind=dp) :: m
   REAL(kind=dp) :: u0

   ! inquire SSA friction exponent from Material properties
   material => GetMaterial()
   IF (.NOT. ASSOCIATED(material)) THEN
     CALL Fatal('SlipCoef_USF', "No material found?")
   ENDIF
   m = ListGetConstReal( Material, 'SSA Friction Exponent',UnFoundFatal=.TRUE.)  
   u0 = ListGetConstReal( Material, 'SSA Friction Threshold Velocity',UnFoundFatal=.TRUE.)  
   
   un=sqrt(VarIn(2)*VarIn(2)+VarIn(3)*VarIn(3))
   
   VarOut=VarIn(1)*(un/(un+u0))**(m)
    
 END FUNCTION CalculTaubRC


 FUNCTION CalculSlc_ls (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(6)    ! alpha, a_ls, a_lsbed, ls, ls2015, bed
   REAL (KIND=dp) :: VarOut      ! slc
   INTEGER :: nodenumber
   
   REAL(kind=dp) :: alpha
         
   IF (VarIn(6).LE.-200.0.AND.VarIn(5).LE.30000.0) THEN
     alpha=VarIn(1)+(VarIn(2)+VarIn(3)*VarIn(6))*(VarIn(4)-VarIn(5))
   ELSE
     alpha=VarIn(1)
   ENDIF
   ! new alpha = old_alpha + (a_ls + a_lsbed*bed)*(ls_new-ls_old)

   VarOut=10**alpha
   
 END FUNCTION CalculSlc_ls

 FUNCTION CalculSlc_haf_u (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(8)    ! slipcoef, zs, zs_inv, zsl, bedrock, rhoi, rhow, multiple
   REAL (KIND=dp) :: VarOut      ! slc
   INTEGER :: nodenumber
   
   TYPE(ValueList_t), POINTER :: material
   REAL(kind=dp) :: lda
   REAL(kind=dp) :: h_af
   REAL(kind=dp) :: h_af_init
   REAL(kind=dp) :: hth

   ! inquire SSA friction exponent from Material properties
   material => GetMaterial()
   IF (.NOT. ASSOCIATED(material)) THEN
     CALL Fatal('SlipCoef_USF', "No material found?")
   ENDIF
   hth = ListGetConstReal( Material, 'SSA Friction Threshold Height',UnFoundFatal=.TRUE.)  

   ! h_af = zs -((rhow/rhoi)*(zsl - bedrock)+bedrock)
   h_af = VarIn(2) - ((VarIn(7)/VarIn(6))*(VarIn(4) - VarIn(5))+VarIn(5))
   h_af_init = MAX(hth, VarIn(3) - ((VarIn(7)/VarIn(6))*(VarIn(4) - VarIn(5))+VarIn(5)))
   !IF (h_af.GT.75.) THEN
   !  lda=+1.0
   !ELSE
   !  VarOut=h_af/MIN(75.0,h_af_init)
   !ENDIF
    
   lda=h_af/MIN(h_af_init,hth)

   IF (lda.LE.0.) THEN
     lda=+0.0
   ENDIF

   ! haf superior to h_threshold 
   IF (h_af.GE.hth) THEN
     lda=+1.0
   ENDIF

   VarOut=VarIn(1)*lda*VarIn(8)
    
 END FUNCTION CalculSlc_haf_u


 FUNCTION CalculMultiple (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(3)    ! SSAVelocity1, SSAVelocity2, multiple
   REAL (KIND=dp) :: VarOut      ! multiple
   INTEGER :: nodenumber
   
   REAL(kind=dp) :: un

   un=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
   
   IF (VarIn(3).LE.0.1) THEN
     VarOut=1.0
   ELSE IF (un.GE.20000) THEN
     VarOut=VarIn(3)*10
   ELSE
     VarOut=VarIn(3)
   ENDIF
    
 END FUNCTION CalculMultiple

 FUNCTION CalculSlc_haf (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(7)    ! slipcoef, zs, zs_inv, zsl, bedrock, rhoi, rhow
   REAL (KIND=dp) :: VarOut      ! slc
   INTEGER :: nodenumber
   
   REAL(kind=dp) :: lda
   REAL(kind=dp) :: h_af
   REAL(kind=dp) :: h_af_init

   ! h_af = zs -((rhow/rhoi)*(zsl - bedrock)+bedrock)
   h_af = VarIn(2) -((VarIn(7)/VarIn(6))*(VarIn(4) - VarIn(5))+VarIn(5))
   h_af_init = MAX(0.9,VarIn(3) -((VarIn(7)/VarIn(6))*(VarIn(4) - VarIn(5))+VarIn(5)))
   !IF (h_af.GT.75.) THEN
   !  lda=+1.0
   !ELSE
   !  VarOut=h_af/MIN(75.0,h_af_init)
   !ENDIF
    
   lda=h_af/h_af_init

   IF (lda.LE.0.) THEN
     lda=+0.0
   ENDIF

   ! haf superior to h_threshold 
   IF (h_af.GE.75) THEN
     lda=+1.0
   ENDIF

   VarOut=VarIn(1)*lda
    
 END FUNCTION CalculSlc_haf


 FUNCTION CalculSlc_haf_old (model, nodenumber, VarIn) RESULT(VarOut)
   USE DefUtils
   IMPLICIT NONE 
   TYPE(Model_t) :: model
   REAL (KIND=dp) :: VarIn(7)    ! slipcoef, zs, zs_inv, zsl, bedrock, rhoi, rhow
   REAL (KIND=dp) :: VarOut      ! slc
   INTEGER :: nodenumber
   
   REAL(kind=dp) :: lda
   REAL(kind=dp) :: h_af
   REAL(kind=dp) :: h_af_init

   ! h_af = zs -((rhow/rhoi)*(zsl - bedrock)+bedrock)
   h_af = VarIn(2) -((VarIn(7)/VarIn(6))*(VarIn(4) - VarIn(5))+VarIn(5))
   h_af_init = MAX(0.9,VarIn(3) -((VarIn(7)/VarIn(6))*(VarIn(4) - VarIn(5))+VarIn(5)))
   !IF (h_af.GT.75.) THEN
   !  lda=+1.0
   !ELSE
   !  VarOut=h_af/MIN(75.0,h_af_init)
   !ENDIF
    
   lda=h_af/h_af_init

   IF (lda.LE.0.) THEN
     lda=+0.0
   ENDIF

   ! bedrock superior to zsl
   IF (VarIn(5).GE.VarIn(4)) THEN
     lda=+1.0
   ENDIF

   VarOut=VarIn(1)*lda
    
 END FUNCTION CalculSlc_haf_old

 FUNCTION CalculSlc_LS_LinearRegression (model, nodenumber, VarIn) RESULT(VarOut)
  USE DefUtils
  IMPLICIT NONE 
  TYPE(Model_t) :: model
  REAL (KIND=dp) :: VarIn(9)    ! alpha_anom, c1, c2, c3, c4, b_lim, bed, distance, multiple
  REAL (KIND=dp) :: VarOut      ! slc
  INTEGER :: nodenumber
  
  REAL(kind=dp) :: alpha

  ! alpha = alpha_anom + (c2*bed+b_lim)*10**(c3*bed)*(10**(c4*ls))
  alpha=VarIn(1)+VarIn(2)+(VarIn(3)*VarIn(7)+VarIn(6))*10**(VarIn(4)*VarIn(7)+VarIn(5)*VarIn(8))

  VarOut=VarIn(9)*10**alpha
   
END FUNCTION CalculSlc_LS_LinearRegression