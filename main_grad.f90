  module precisions
  implicit none
  INTEGER, parameter :: dp=kind(1.0d0)
  Include 'mpif.h'
  logical     :: logik, ouvert 
  end module precisions
!************************************************************************************************************************************************************

  MODULE GlobalParam
  USE precisions
  IMPLICIT NONE
  INTEGER ::  ImpreE, ImpreF, argunit = 6
  INTEGER :: isave, physical_model, MyUnit = 30
  INTEGER::  Nr, Ntheta, iterfinal, cond_lim
  INTEGER::  ORDRE, TS_GEOM, NUMB=1
  REAL (KIND = DP)::  H_0, amplitude,  frottcoeff, disscoeff
  REAL (KIND = DP) ::  omega, alpha, beta, P_0
  REAL (KIND = DP) :: period_time, lambda, GAMMA
  REAL (KIND = DP) ::  CFL, TIMEOUT, angle, g, phi2
  REAL (KIND = DP) ::  Rplus, Rminus, dev, q0_jump
  REAL (KIND = DP)  ::  L1, L2
  REAL (KIND = DP), PARAMETER ::  EPS = 1.d-8
  REAL (KIND = DP), PARAMETER::  pi = 4.0d0*ATAN(1.0d0)

 INTEGER , PARAMETER  :: Nv_Prim = 12, H_pv= 1
 INTEGER , PARAMETER  :: U_pv= 2, V_pv= 3, ein_pv = 7
 INTEGER , PARAMETER  :: P11_pv= 4, P12_pv= 5, P22_pv= 6, press_pv= 8
 INTEGER , PARAMETER  :: sound_ar_pv= 9, sound_br_pv= 10
  INTEGER , PARAMETER  :: sound_ath_pv= 11, sound_bth_pv= 12 
END MODULE GlobalParam
!************************************************************************************************************************************************************

MODULE ModeleInterface
  USE GlobalParam
  USE precisions
  IMPLICIT NONE

CONTAINS
!--------------------------------------------------------
  
 FUNCTION InternalEn(h, p11, p22) RESULT(InternalE)
 REAL (KIND = DP)      :: h, p11, p22, InternalE
  InternalE =  0.5d0*(g*h + p11 + p22)
 END FUNCTION InternalEn
!--------------------------------------------------------

 FUNCTION Pression(h, p11) RESULT(press)
 REAL (KIND = DP)      :: h, p11, press
  press =  0.5d0*g*h**2.d0 + h*p11
 END FUNCTION Pression
!--------------------------------------------------------

FUNCTION Sound_ar(h, p11) RESULT(Sound_a)
 REAL (KIND = DP)      :: h, p11, Sound_a
  Sound_a =  dsqrt(g*h + 3.d0*p11)
 END FUNCTION Sound_ar
!--------------------------------------------------------

FUNCTION Sound_atheta(h, p22) RESULT(sound_ath)
 REAL (KIND = DP)      :: h, p22, sound_ath
  sound_ath =  dsqrt(g*h + 3.d0*p22)
 END FUNCTION Sound_atheta
!--------------------------------------------------------

FUNCTION Sound_br( p11) RESULT(sound_b)
 REAL (KIND = DP)      :: h, p11, sound_b
  sound_b =  dsqrt(p11)
 END FUNCTION Sound_br
!--------------------------------------------------------

FUNCTION Sound_btheta( p22) RESULT(sound_bth)
 REAL (KIND = DP)      :: h, p22, sound_bth
  sound_bth =  dsqrt(p22)
 END FUNCTION Sound_btheta
!--------------------------------------------------------

END MODULE ModeleInterface

!************************************************************************************************************************************************************
  PROGRAM code2D
  USE precisions
  use GlobalParam
  USE ModeleInterface
    IMPLICIT NONE  
    INTEGER  :: iv, I, IT, ir, il, Nl
    INTEGER  :: nbr, reste, ncpu, rang, nfile
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:):: CONS, FLUX, Prim
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:):: Umaxtampon,Vmaxtampon,  h_new
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:):: MaxVP, MinVp, R, L, bottom
    REAL (KIND = DP)  :: T1_CPU, T2_CPU, TIME,  TIME2
    REAL (KIND = DP)  :: Lr, Ltheta, DR, Dtheta, Dh, DT, dt2, dl
    REAL (KIND = DP) :: Cthetamax, crmax
    REAL(KIND=DP) :: Rmax, RMIN,  Lmin, Lmax, UMAX, vmax, dlmin
    REAL(KIND=DP)  :: TAMPONRMIN, TAMPONRMAX, tampondt
    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:):: NamesOfPrimVar
    CHARACTER*3                   :: numrg
  
  !//UNIQUEMENT POUR LE PARALLELE
  !//Variables Calcul Parallele
    INTEGER :: code
  !//Initialisation du calcul parallele
    CALL MPI_INIT(code)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  !----------------------------------------------------------------------------------------------------
    CALL CPU_TIME(T1_CPU)
  !----------------------------------------------------------------------------------------------------
    CALL LECTURE_DONNEES()
  !----------------------------------------------------------------------------------------------------
    Lr = Rplus - Rminus
    Ltheta = 2.d0*pi
 
    DR = Lr/DFLOAT(Nr)

    nbr = int( nr/ncpu )
  reste = nr-(ncpu)*nbr

  if (rang .lt. (reste)) then
      nbr = nbr +1
  endif
  nr = nbr 

    Dtheta = Ltheta/DFLOAT(Ntheta)
    Dl = Dtheta 
    Nl = Ntheta

   ALLOCATE( Prim(Nv_Prim,0:Nr+1,0:Ntheta+1),&
   & CONS(7,1:Nr,1:Ntheta), FLUX(7,0:Nr,0:Ntheta))  
   ALLOCATE( NamesOfPrimVar(Nv_Prim) )
   ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim), R(0:Nr+1), L(1:Ntheta))
   ALLOCATE( Umaxtampon(1:Nr,1:Ntheta),&
    &Vmaxtampon(1:Nr,1:Ntheta), h_new(1:Nr,1:Ntheta) )
   ALLOCATE( bottom(1:Nr) )
   bottom = 0.d0
   flux=0.d0; cons = 0.d0; prim = 0.d0; umax = 0.d0; h_new(:,:) = 0.d0
   vmax = 0.d0; Crmax = 0.d0; cthetamax = 0.d0; tampondt = 0.d0

   isave = -1
   if (rang == 0) then 
   print*, 'Nr =', Nr,'Ntheta =', Ntheta, 'Nl =', Nl
   print*, 'Dr =', Dr, 'Dtheta =', Dtheta, 'rang =', rang
   !----------------------------------------------------------------------------------------------------  
    NamesOfPrimVar(H_pv)               = "Depht Lenght"
    NamesOfPrimVar(U_pv)               = "Velocity r"
    NamesOfPrimVar(V_pv)               = "Velocity th"
    NamesOfPrimVar(P11_pv)             = "tensor P11"
    NamesOfPrimVar(P12_pv)             = "tensor P12"
    NamesOfPrimVar(P22_pv)             = "tensor P22"
    NamesOfPrimVar(ein_pv)             = "energy int"
    NamesOfPrimVar(press_pv)           = "pression r"
    NamesOfPrimVar(sound_ar_pv)        = "a-waves r"
    NamesOfPrimVar(sound_br_pv)        = "b-waves r"
    NamesOfPrimVar(sound_ath_pv)       = "a-waves th"
    NamesOfPrimVar(sound_bth_pv)       = "b-waves th"
  !----------------------------------------------------------------------------------------------------
    DO iv = 1, Nv_Prim
       WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ",&
       & iv, Nr, Ntheta, 'rang =', rang
    END DO
    WRITE(6,*) " >>> End of LECTURE_DONNEES"
  endif
  !----------------------------------------------------------------------------------------------------
  !INITIALISATION  
   CALL INITIALISATION(Prim, CONS, DR, Dtheta, R, L, Lr, Ltheta, bottom, rang, reste)
   CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
   CALL COMMUNICATION(Prim, ncpu, rang, R)
  !----------------------------------------------------------------------------------------------------
   ouvert = .true.
   Nfile = int(TIMEOUT/period_time) + 2
   TIME = 0.D0
   TIME2  = period_time ;
   IT = 1
   CALL PutonScreen(rang)
  !----------------------------------------------------------------------------------------------------
  TAMPONRMIN = MINVAL(R(:))
  TAMPONRMAX = MaxVal(R(:))

  CALL MPI_ALLREDUCE(TAMPONRMIN, RMIN, NUMB, MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,code)
  CALL MPI_ALLREDUCE(TAMPONRMAX, RMAX, NUMB,MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)

   Lmin = MINVAL(L(:))
   Lmax = MaxVal(L(:))
  IF (rang == 0) THEN   
  WRITE(6,'( 2(A10,E15.6))')  " RMIN = ", RMIN, &
       &                      " Rmax = ", Rmax
  WRITE(6,'( 2(A10,E15.6))') " Lmin = ", Lmin, &
       &                     " Lmax = ", Lmax
   ENDIF     
!----------------------------------------------------------------------------------------------------
   
 !  logik=.true.
   CALL Ecriture_donnees(R,L, Prim, time, DR, Dtheta, bottom, rang, Nfile, reste)
   if (cond_lim == 8) then
    call newton(R,L, ncpu, rang, dr)
   endif
 ! logik=.false.

! BOUCLE SUR LE TEMPS
!----------------------------------------------------------------------------------------------------
Time_loop: do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT))
!----------------------------------------------------------------------------------------------------
  DO ir = 1, Nr    
  DO il = 1, Ntheta
   Umaxtampon(ir,il) = dmax1(dabs(Prim(U_pv,ir,il)  - &
    & Prim(sound_ar_pv, ir,il) ), dabs(Prim(U_pv,ir,il)  + Prim(sound_ar_pv, ir,il) ))
  ENDDO
  ENDDO

   Umax = MAXVAL(Umaxtampon(1:Nr,1:Ntheta))

  DO ir = 1, Nr    
  DO il = 1, Ntheta
   vmaxtampon(ir,il) = dmax1( (dabs(Prim(v_pv,ir,il) &
   & - Prim(sound_ath_pv, ir,il)) )/r(ir), &
   & (dabs(Prim(v_pv,ir,il)  + Prim(sound_ath_pv, ir,il)) )/r(ir))
   
  ENDDO
  ENDDO
   vmax = MAXVAL(vmaxtampon(1:Nr,1:Ntheta))

  if (ntheta.gt.1) then
      tamponDt = dMIN1( Dr/(umax), dtheta/(vmax) ) ! dMIN1( Dr/(umax + crmax), rminus*dtheta/(vmax+ cthetamax)) 
  else
     tampondt = Dr/(umax) 
  endif

    tamponDT   = CFL*tampondt
    dt2 = 0.5d0*dt
CALL MPI_ALLREDUCE(TAMPONdt, dt, NUMB,MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,code)
if (rang == 0 .and. it == 1) then 
    WRITE(6,*) " Dt = ", Dt, 'dtheta =', Dtheta, 'dr =', dr, rang
endif

! IF (COND_LIM == 7 .AND. IT == 1) THEN

!     WRITE(numrg,'(i3.3)') rang
!     OPEN(MyUnit+10,FILE = './resu/BOTTOM_'//numrg//'.out')
!     DO IR = 1, NR
!           BOTTOM(IR) = (5.6d0)**2.d0*1.d-4*(1.d0/(Rplus) - 1.d0/r(ir)) ! (Rminus + 2.D0*L1)

!         WRITE(MyUnit+10,'(2(E20.13,1X))') R(IR), BOTTOM(IR)
!     ENDDO
!     CLOSE(MyUnit+10)

! ENDIF

IF (COND_LIM == 7 .AND. IT == 1) THEN
   WRITE(numrg,'(i3.3)') rang
    OPEN(MyUnit+10,FILE = './resu/BOTTOM2_'//numrg//'.out')
    DO IR = 1, NR
         IF (R(IR) - RMINUS <= 2.D0*L1) THEN 
          BOTTOM(IR) = DEV*((R(IR) - RMINUS - L1)**2.D0 - L1**2.D0 )**2.D0/L1**4.D0 
         ELSE
          BOTTOM(IR) =   (r(ir) - 2.D0*L1 -  RMINUS)*Dtan(angle) !
         ENDIF
       
        WRITE(MyUnit+10,'(2(E20.13,1X))') R(IR), BOTTOM(IR)
    ENDDO
   CLOSE(MyUnit+10)
 ENDIF
    

    TIME = TIME + DT
    !if (time > 28.d0) period_time = 0.1d0

!   CALL euler_method( CONS, Prim,R, dr, DT2, it)

!  !----------------------------------------------------------------------------------------------------
!   CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
!   CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
!   CALL COMMUNICATION(Prim, ncpu, rang, R)
  !----------------------------------------------------------------------------------------------------
     CALL HLLC_r_sub1(prim,flux, dt, R, L, DR, dtheta, IT, time)
    !call ts_nonconserv_rk2( CONS, Prim, alpha*DT, R, DR, it )
    CALL godunov_r_sub1(cons,flux,dt,R, DR, Prim)
   ! call ts_nonconserv_rk2(CONS, Prim, (1.d0 -alpha)*DT, R, DR, it)
    CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
    CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
    CALL COMMUNICATION(Prim, ncpu, rang, R)
!   !---------------------------------------------------------------------------------------------------- 
      CALL HLLC_r_sub2(prim,flux, dt, R, L, DR, dtheta, IT, time)
      CALL godunov_r_sub2(cons,flux,dt, R, DR)
      CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
      CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
      CALL COMMUNICATION(Prim, ncpu, rang, R)
 !---------------------------------------------------------------------------------------------------- 
   IF (ntheta .gt. 1) then
       if (IT == 1) print*, 'je calcule suivant theta'

        CALL HLLC_theta_sub1(prim,flux, dt, R,L, DR,dtheta, IT, time, rang, reste) 
        CALL godunov_theta_sub1(cons,flux,dt,Dtheta, R, dr)
        CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
        CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
        CALL COMMUNICATION(Prim, ncpu, rang, R)
       !----------------------------------------------------------------------------------------------------

        CALL HLLC_theta_sub2(prim,flux, dt, R, L, DR, dtheta, IT, time, rang, reste) 
        CALL godunov_theta_sub2(cons,flux,dt,Dtheta, R, dr)
        CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
        CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
        CALL COMMUNICATION(Prim, ncpu, rang, R)
  ENDIF  
 !----------------------------------------------------------------------------------------------------
 ! CALL rk2( CONS, Prim,R, dr, DT, it)

  CALL euler_method( CONS, Prim,R, dr, DT, it)
  CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
  CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
  CALL COMMUNICATION(Prim, ncpu, rang, R)

 !----------------------------------------------------------------------------------------------------



!   CALL Coriolis( CONS, Prim,R, dr, DT, it)
!   CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
!   CALL cond_lim_jump(Prim, TIME, R, L, dr, dtheta, rang, ncpu)
!   CALL COMMUNICATION(Prim, ncpu, rang, R)
  !----------------------------------------------------------------------------------------------------
  IT = IT + 1
  !----------------------------------------------------------------------------------------------------
   IF (TIME2.LE.TIME) THEN
    if (cond_lim == 8) then 
     call newton(R,L, ncpu, rang, dr)
    endif

  	CALL Ecriture_donnees(R,L, Prim, time, DR, Dtheta, bottom, rang, nfile, reste)
    CALL PutonScreen(rang)
    
   END IF
  !----------------------------------------------------------------------------------------------------
    IF (TIME2.LE.TIME) THEN
      PRINT*, 'EN', IT, 'ITERATIONS, ', ' TIME:', TIME
      TIME2 = TIME + period_time
    END IF
  !----------------------------------------------------------------------------------------------------
  ENDDO TIME_LOOP
   !! FIN BOUCLE SUR LE TEMPS

   CALL Ecriture_donnees(R,L, Prim, time, DR, Dtheta, bottom, rang, nfile, reste)
  !----------------------------------------------------------------------------------------------------
   CALL CPU_TIME(T2_CPU)
   DEALLOCATE(Prim, R, L, CONS, FLUX, MinVP, MaxVp, NamesOfPrimVar)  

   DEALLOCATE( Umaxtampon,Vmaxtampon, h_new, bottom )

  !----------------------------------------------------------------------------------------------------
  PRINT*, 'L EXECUTION DU PROGRAMME A PRIS', T2_CPU - T1_CPU
  PRINT*, 'EN', IT-1, 'ITERATIONS, ', ' TIME:', TIME
  CALL MPI_FINALIZE(code)
  CONTAINS
  !************************************************************************************************************************************************************

  SUBROUTINE LECTURE_DONNEES()
  USE precisions
  use GlobalParam
  IMPLICIT NONE
    OPEN(UNIT=21, FILE = 'data.inp', STATUS = 'OLD')
    READ(21,*) physical_model   ! 1 sv/ 2 shear 
    READ(21,*) cond_lim !COND LIM :1 box/2 absorbtion/3 batteur/4 jump;
    READ(21,*) ORDRE, TS_GEOM  ! ORDRE DU SHEMA:1 GODUNOV,2 MUSCEL; TS_GEOM:0 DANS SHEMA,1 SPLITTING AVEC ALPHA
    READ(21,*) angle           ! inclination angle
    READ(21,*) Nr, Ntheta      ! NUMBER OF CELLS
    READ(21,*) Rplus, Rminus   ! RADIUS
    READ(21,*) TIMEOUT         ! OUTPUT TIME
    READ(21,*) iterfinal       ! Iteration final
    READ(21,*) g, CFL          ! acceleration due to gravity
    READ(21,*) H_0, dev             ! INITIAL DEPTH; hydraulic jump hauteur de deversoir (obstacle)
    READ(21,*) phi2, P_0       ! patit enstrophy
    READ(21,*) frottcoeff, disscoeff  ! cf, cr
    READ(21,*) amplitude, q0_jump     ! amplitude des perturbation
    READ(21,*) ImpreE, ImpreF, period_time ! POUR IMPRIMER LES FICHIER SUR ECRAN ET DANS LE FICHIER
    READ(21,*) omega, lambda, GAMMA                 ! FOR WATER GLASS; lambda for nonstationar solution valeur de  p11 initialement
    READ(21,*) alpha, beta                 ! JEU DE TS GEOMETRIQUE ALPHA*TS(U^N)+(1-ALPHA)*TS(U')
    READ(21,*) L1, L2 
    close(21)
  return
  END SUBROUTINE
!************************************************************************************************************************************************************

SUBROUTINE INITIALISATION(Prim, CONS, DR, Dtheta,&
& R, L, Lr, Ltheta, bottom, rang, reste)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ir, il, rang, reste
REAL (KIND = DP) :: Dl,Dtheta, DR, Lr, Ltheta, U_0, fr
REAL (KIND = DP) :: ein, pres, u, v, h, p11, p12, p22, r1
REAL (KIND = DP) :: R(0:Nr+1), L(1:Ntheta)
REAL (KIND = DP) :: Prim(Nv_Prim, 0:Nr+1,0:Ntheta+1)
REAL (KIND = DP) :: CONS(7,1:Nr,1:Ntheta), bottom(1:Nr)

if (rang .eq. 0) R(0) = Rminus - dr*0.5D0
if (rang .eq. ncpu - 1) R(nr+1) = Rplus + dr*0.5d0
DO ir = 1, Nr
  if (rang .lt. (reste)) then
    R(ir) =Rminus + 0.5D0*DR + (ir-1)*DR + nr*rang*dr 
  else 
    R(ir) =Rminus + 0.5D0*DR + (ir-1)*DR+ &
    &(reste)*(nr+1)*dr + nr*(rang - reste)*dr
  endif


DO il = 1, Ntheta 

  Dl = Dtheta !R(ir)*
  L(il) = 0.5D0*Dl + (il-1)*Dl
!---------------------------------------------------------------------
 IF (COND_LIM == 2) THEN  ! barrage sur h
  
  if ( r(ir) .le. Rminus + 0.5d0*Lr) then  
   Prim(H_pv, ir, il) = 0.02d0
  else
   Prim(H_pv, ir, il) = 0.01d0 
  endif
  !Prim(U_pv,ir, il)  = 0.d0
  Prim(V_pv,ir, il)  = 0.d0
  Prim(P11_pv,ir,il) = 1.d-6
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) = 1.d-6
ELSEIF (cond_lim == 3) THEN !  verre d'eau qui tourne sur le disque
  Prim(H_pv, ir, il) = H_0  + omega*omega*r(ir)*r(ir)/(2.d0*g) ! 2.5d0*1.d-5 ! 
  Prim(U_pv,ir, il)  = 0.d0
  Prim(V_pv,ir, il)  = omega*R(ir)
  Prim(P11_pv,ir,il) = 0.d0 
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) = 0.d0  
ELSEIF (cond_lim == 4) THEN !  nonstationar exacte
  Prim(H_pv, ir, il) = H_0  
  Prim(U_pv,ir, il)  = r(ir)
  Prim(V_pv,ir, il)  = 0.d0 
  Prim(P11_pv,ir,il) = P_0 
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) = P_0
ELSEIF (COND_LIM == 5 ) THEN  ! rien qui se passe
  Prim(H_pv, ir, il) = H_0
  Prim(U_pv,ir, il)  = 0.d0
  Prim(V_pv,ir, il)  = 0.d0
  Prim(P11_pv,ir,il) = 1.d-7 
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) = 1.d-7 
ELSEIF (COND_LIM == 7 ) THEN 

!          IF (R(IR) - RMINUS <= L1  ) THEN 
!           BOTTOM(IR) = DEV*((R(IR) - RMINUS - L1)**2.D0 - L1**2.D0 )**2.D0/L1**4.D0 
!          ELSEIF (R(IR) - RMINUS <= L1+ L2 .and.  R(IR) - RMINUS > L1) THEN 
!           BOTTOM(IR) = DEV*(R(IR) - RMINUS - L1 - L2)**4.D0/L2**4.D0 
!          ELSE
!           BOTTOM(IR) =   (r(ir) - L1 - L2 -  RMINUS)*Dtan(angle) !
!          ENDIF

   
    !   BOTTOM(IR) =   (5.6d0*0.01d0)**2.d0*(1.d0/((Rplus)) - 1.d0/r(ir) ) ! (Rminus + 2.D0*L1)
      
 IF (R(IR) - RMINUS <= 2.D0*L1) THEN 
  BOTTOM(IR) = DEV*((R(IR) - RMINUS - L1)**2.D0 - L1**2.D0 )**2.D0/L1**4.D0 
 ELSE
  BOTTOM(IR) =   (r(ir) - 2.D0*L1 -  RMINUS)*Dtan(angle) !
 ENDIF

if ( R(IR) - RMINUS <= 0.3d0*lr) then ! *(1.d0+ 1.d-3*(dsin(10.d0*L(il) ) ))
Prim(H_pv, ir, il) = 1.5d0*H_0 !0.01d0 !- BOTTOM(IR) ! 1.5d0*H_0 !
else
Prim(H_pv, ir, il) = H_0 !- BOTTOM(IR) !H_0 
endif

  Prim(U_pv,ir, il)  = -q0_jump/(2.d0*pi*r(ir)*Prim(H_pv, ir, il)  )   
  Prim(V_pv,ir, il)  = 0.d0 ! 1.d-3*r(ir)
  Prim(P11_pv,ir,il) = PHI2*Prim(H_pv, ir, il)**2.d0*(1.d0 + eps) !*Prim(H_pv, ir, il)**2.d0 !
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) =  eps*PHI2*Prim(H_pv, ir, il)**2.d0*(1.d0 + eps) !0.5d0*PHI2*Prim(H_pv, ir, il)**2.d0*(1.d0 + eps)  !*Prim(H_pv, ir, il)**2.d0  !*PHI2*Prim(H_pv, ir, il)**2.d0 



ELSEIF (COND_LIM == 8 ) THEN 

  Prim(H_pv, ir, il) = H_0
  Prim(U_pv,ir, il)  = -q0_jump/( r(ir)*H_0  )  ! -2.d-1 !
  Prim(V_pv,ir, il)  = 0.d0
  Prim(P11_pv,ir,il) = eps
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) = eps
ELSEIF (COND_LIM == 9 ) THEN ! barrage sur u 

  Prim(H_pv, ir, il) = H_0
  if (r(ir) .le. Rminus + 0.5d0*lr) then
   Prim(U_pv,ir, il)  = 1.d0
  else
   Prim(U_pv,ir, il)  = -1.d0
  endif
  Prim(V_pv,ir, il)  = 0.d0
  Prim(P11_pv,ir,il) = EPS
  Prim(P12_pv,ir,il) = 0.d0
  Prim(P22_pv,ir,il) = EPS

ELSEIF (COND_LIM == 10 ) THEN ! nonstationar nonspherical turbulent solution 

  Prim(H_pv, ir, il) = 1.d0
  Prim(U_pv,ir, il)  = r(ir)*beta*DSIN(L(il))
  Prim(V_pv,ir, il)  = -r(ir)*beta*DCOS(L(il))
  Prim(P11_pv,ir,il) = lambda*dcos(L(il))*dcos(L(il)) + GAMMA*dsin(L(il))*dsin(L(il))
  Prim(P12_pv,ir,il) = (gamma - lambda)/2.d0*dsin(2.d0*L(il))
  Prim(P22_pv,ir,il) = gamma*dcos(L(il))*dcos(L(il)) + lambda*dsin(L(il))*dsin(L(il))
ENDIF
  h = Prim(H_pv, ir, il); p11 = Prim(P11_pv, ir, il); p22 = Prim(P22_pv, ir, il)

  Prim(sound_ar_pv,ir,il)= Sound_ar(h, p11)  !DSQRT( g*Prim(H_pv, ir, il) + 3.d0*Prim(P11_pv, ir, il) )
  Prim(sound_br_pv,ir,il)= Sound_br( p11) !DSQRT(Prim(P11_pv, ir, il) )
  Prim(sound_ath_pv,ir,il)= Sound_atheta(h, p22) !DSQRT( g*Prim(H_pv, ir, il) + 3.d0*Prim(P22_pv, ir, il) )
  Prim(sound_bth_pv,ir,il)=  Sound_btheta(p22) !DSQRT(Prim(P22_pv, ir, il) )
  Prim(press_pv,ir,il) = Pression(h, p11) !g*Prim(H_pv, ir, il)*Prim(H_pv, ir, il)/2.d0 + Prim(P11_pv, ir, il)*Prim(H_pv, ir, il)
  Prim(ein_pv,ir,il) = InternalEn(h, p11, p22)!Prim(ein_pv,ir,il) = 0.5d0*( Prim(H_pv, ir, il)*g + Prim(P11_pv, ir, il) + Prim(P22_pv, ir, il) )

ENDDO
ENDDO

  DO ir = 1, Nr
  DO il = 1, Ntheta 
  h = Prim(H_pv, ir, il); u = Prim(U_pv,ir, il); v = Prim(V_pv,ir, il); 
  p11 = Prim(P11_pv,ir,il); p12 = Prim(P12_pv,ir,il); p22 = Prim(P22_pv,ir,il)
  ein = Prim(ein_pv,ir,il)

  CONS(1,ir,il) = R(ir)*h
  CONS(2,ir,il) = R(ir)*h*u
  CONS(3,ir,il) = R(ir)*h*v
  CONS(4,ir,il) = R(ir)*h*p11
  CONS(5,ir,il) = R(ir)*p12
  CONS(6,ir,il) = R(ir)*h*p22
  CONS(7,ir,il) = R(ir)*h*(ein+ 0.5d0*(u**2.d0 + v**2.d0) ) 
  ENDDO
  ENDDO
  return
  END SUBROUTINE
  !************************************************************************************************************************************************************

  SUBROUTINE NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
  USE precisions
  use GlobalParam
   IMPLICIT NONE
   INTEGER :: ir,il, it
   REAL (KIND = DP) :: p11, p12, p22, h, u, v, dr, r1
   REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), CONS(7,1:Nr,1:Ntheta), R(0:Nr+1)
  ! CALCUL NOUVELS VARIABLES PRIMITIVES
  DO ir = 1, Nr    
  DO il = 1, Ntheta
   Prim(H_pv,ir, il) = CONS(1,ir, il)/R(ir)

!    IF (R(IR) - RMINUS <= 2.D0*L1) THEN 
!     BOTTOM(IR) = DEV*((R(IR) - RMINUS - L1)**2.D0 - L1**2.D0 )**2.D0/L1**4.D0
!     IF (BOTTOM(IR) .GT. Prim(H_pv,ir, il) ) THEN
!      Prim(H_pv,ir, il) = BOTTOM(IR)
!     ENDIF
!    ENDIF

!    if (cons(1,ir, il) .le. 1.d-8) then
!     print*, 'Prim(H_pv,ir, il) =', Prim(H_pv,ir, il), 'cons(1) =', cons(1,ir, il), 'it = ', it, 'ir=', ir, 'il = ', il
!    ! stop
!    endif

  Prim(U_pv,ir, il) = CONS(2,ir, il)/CONS(1,ir, il)
  Prim(V_pv,ir, il) = CONS(3,ir, il)/CONS(1,ir, il)
!------------------------------------------------------------------------------------------------------------------------
  if (dabs(Prim(V_pv,ir, il)).le.1.d-8) then
   Prim(V_pv,ir, il)=0.d0; cons(3, ir, il) = 0.d0
  endif
  if (dabs(Prim(U_pv,ir, il)).le.1.d-8) then 
   Prim(U_pv,ir, il)=0.d0; cons(2, ir, il) = 0.d0
  endif
!------------------------------------------------------------------------------------------------------------------------
     Prim(P11_pv,ir, il) = CONS(4,ir, il)/CONS(1,ir, il)
     Prim(P11_pv,ir, il) = dmax1(Prim(P11_pv,ir, il), 1.d-8)
     !Prim(P11_pv,ir, il) = eps
    
     Prim(P12_pv,ir, il) = CONS(5,ir, il)/R(ir)
     !if (Prim(P12_pv,ir, il) >= 0.d0 ) then
      Prim(P12_pv,ir, il) = dmin1(Prim(P12_pv,ir, il), &
        &dsqrt(Prim(P11_pv,ir, il)*Prim(P22_pv,ir, il)) ) 
     !else
      Prim(P12_pv,ir, il) = dmax1(Prim(P12_pv,ir, il), &
        &-dsqrt(Prim(P11_pv,ir, il)*Prim(P22_pv,ir, il) ) ) 
     !endif
    ! Prim(P12_pv,ir, il) = 0.d0

     Prim(P22_pv,ir, il) = CONS(6,ir, il)/CONS(1,ir, il)
     Prim(P22_pv,ir, il) = dmax1(Prim(P22_pv,ir, il), 1.d-8)
     !Prim(P22_pv,ir, il) = eps
!------------------------------------------------------------------------------------------------------------------------
     p11 = Prim(P11_pv,ir, il); p22 = Prim(P22_pv,ir, il); 
     p12 = Prim(P12_pv,ir, il); h = Prim(H_pv,ir, il)
     u = Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
!------------------------------------------------------------------------------------------------------------------------
   
     Prim(ein_pv,ir, il) = InternalEn(h, p11, p22) 
     Prim(ein_pv,ir, il) = dmax1(Prim(ein_pv,ir, il), 0.d0)

     Prim(sound_ar_pv,ir, il) = DSQRT(g*h + 3.d0*p11)
     Prim(sound_br_pv,ir, il) = DSQRT(p11)  
     Prim(sound_ath_pv,ir, il) = DSQRT(g*h + 3.d0*p22)
     Prim(sound_bth_pv,ir, il) = DSQRT(p22) 
     Prim(press_pv,ir, il) = g*h**2.d0/2.d0 + p11*h

  END DO  
  END DO
  return
  END SUBROUTINE
!**********************************************************************************************************************************************************

SUBROUTINE cond_lim_jump(Prim, TIME, R,L, dr, dtheta, rang, ncpu)
USE precisions
USE GlobalParam
IMPLICIT NONE
INTEGER ::  ir, il, iv, rang, ncpu 
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), R(0:Nr+1), L(1:Ntheta)
REAL (KIND = DP) :: TIME, dr, dtheta, dl, cd
REAL (KIND = DP) :: h, p11, p22,u,v, p12, r1

DO il = 1, Ntheta
 
IF (COND_LIM == 2 .or. COND_LIM == 9) THEN
    if (rang .eq. 0) then  
    DO iv = 1, 6
     Prim(iv,0, il)        = Prim(iv,1, il) 
    ENDDO
    Prim(u_pv,0, il)      = -Prim(u_pv,1, il)
    Prim(p12_pv,0, il)    = -Prim(p12_pv,1, il)
  elseif (rang .eq. ncpu - 1) then 
    DO iv = 1, 6
     Prim(iv,Nr+1, il)     = Prim(iv,Nr, il)
    ENDDO
    Prim(u_pv,nr+1, il)   = -Prim(u_pv,nr, il)
    Prim(p12_pv,nr+1, il)   = -Prim(p12_pv,nr, il)
  endif

elseif (cond_lim == 3) then ! water glass
  if (ncpu > 1) then 

if (rang .eq. 0) then 
  Prim(h_pv, 0, il) =     H_0 + omega*omega*R(0)*R(0)/(2.d0*g) ! Prim(h_pv, 1, il)!
  Prim(u_pv, 0, il) =    0.d0 !
  Prim(v_pv, 0, il) =    omega*R(0) !  Prim(v_pv, 1, il) ! 
  Prim(p11_pv, 0, il) =  0.d0 !Prim(p11_pv, 1, il) !
  Prim(p12_pv, 0, il) =  0.d0 !Prim(p12_pv, 1, il)
  Prim(p22_pv, 0, il) =  0.d0 !Prim(p22_pv, 1, il) !1.d-8 !

elseif (rang .eq. ncpu - 1) then 
  Prim(h_pv, Nr+1, il) =   H_0 + omega*omega*R(Nr + 1)*R(Nr + 1)/(2.d0*g) ! Prim(h_pv, Nr, il) !
  Prim(u_pv, Nr+1, il) =   0.d0 !Prim(u_pv, Nr, il) !Prim(u_pv, Nr, il) !
  Prim(v_pv, Nr+1, il) =   omega*R(nr+1) ! Prim(v_pv, Nr, il) ! 
  Prim(p11_pv,Nr+1, il)  = 0.d0 !Prim(p11_pv,Nr, il)
  Prim(p12_pv, Nr+1, il) = 0.d0 !Prim(p12_pv, Nr, il)
  Prim(p22_pv, Nr+1, il) = 0.d0 !Prim(p22_pv, Nr, il)

endif
else
  Prim(h_pv, 0, il) =     H_0 + omega*omega*R(0)*R(0)/(2.d0*g) ! Prim(h_pv, 1, il)!
  Prim(u_pv, 0, il) =    0.d0 !
  Prim(v_pv, 0, il) =    omega*R(0) !  Prim(v_pv, 1, il) ! 
  Prim(p11_pv, 0, il) =  0.d0 !Prim(p11_pv, 1, il) !
  Prim(p12_pv, 0, il) =  0.d0 !Prim(p12_pv, 1, il)
  Prim(p22_pv, 0, il) =  0.d0 !Prim(p22_pv, 1, il) !1.d-8 !


  Prim(h_pv, Nr+1, il) =   H_0 + omega*omega*R(Nr + 1)*R(Nr + 1)/(2.d0*g) ! Prim(h_pv, Nr, il) !
  Prim(u_pv, Nr+1, il) =   0.d0 !Prim(u_pv, Nr, il) !Prim(u_pv, Nr, il) !
  Prim(v_pv, Nr+1, il) =   omega*R(nr+1) ! Prim(v_pv, Nr, il) ! 
  Prim(p11_pv,Nr+1, il)  = 0.d0 !Prim(p11_pv,Nr, il)
  Prim(p12_pv, Nr+1, il) = 0.d0 !Prim(p12_pv, Nr, il)
  Prim(p22_pv, Nr+1, il) = 0.d0 !Prim(p22_pv, Nr, il)

  endif
  ELSEIF (COND_LIM == 4) THEN ! nonstationar exacte

if (ncpu > 1) then 
if (rang .eq. 0) then 
  Prim(h_pv, 0, il) =     H_0/(1.d0 + time)**2.d0
  Prim(u_pv, 0, il) =    R(0)/(1.d0 + time)
  Prim(v_pv, 0, il) =    0.d0
  Prim(p11_pv, 0, il) =  P_0/(1.d0 + time)**2.d0
  Prim(p12_pv, 0, il) =  0.d0 
  Prim(p22_pv, 0, il) =  P_0/(1.d0 + time)**2.d0

elseif  (rang .eq. ncpu -1) then 

  Prim(h_pv, Nr+1, il) =    H_0/(1.d0 + time)**2.d0
  Prim(u_pv, Nr+1, il) =   R(nr+1)/(1.d0 + time)
  Prim(v_pv, Nr+1, il) =   0.d0
  Prim(p11_pv,Nr+1, il)  = P_0/(1.d0 + time)**2.d0
  Prim(p12_pv, Nr+1, il) = 0.d0 
  Prim(p22_pv, Nr+1, il) = P_0/(1.d0 + time)**2.d0
endif

else
  Prim(h_pv, 0, il) =     H_0/(1.d0 + time)**2.d0
  Prim(u_pv, 0, il) =    R(0)/(1.d0 + time)
  Prim(v_pv, 0, il) =    0.d0
  Prim(p11_pv, 0, il) =  P_0/(1.d0 + time)**2.d0
  Prim(p12_pv, 0, il) =  0.d0 
  Prim(p22_pv, 0, il) =  P_0/(1.d0 + time)**2.d0


  Prim(h_pv, Nr+1, il) =    H_0/(1.d0 + time)**2.d0
  Prim(u_pv, Nr+1, il) =   R(nr+1)/(1.d0 + time)
  Prim(v_pv, Nr+1, il) =   0.d0
  Prim(p11_pv,Nr+1, il)  = P_0/(1.d0 + time)**2.d0
  Prim(p12_pv, Nr+1, il) = 0.d0 
  Prim(p22_pv, Nr+1, il) = P_0/(1.d0 + time)**2.d0
endif

ELSEIF (COND_LIM == 5) THEN ! rien qui se passe

if (rang .eq. 0) then 
  DO iv = 1, 6
   Prim(iv,0, il)        = Prim(iv,1, il) 
  ENDDO
  Prim(u_pv,0, il)      = -Prim(u_pv,1, il)
  Prim(p12_pv,0, il)    = -Prim(p12_pv,1, il)

elseif  (rang .eq. ncpu -1) then 
  DO iv = 1, 6
   Prim(iv,Nr+1, il)     = Prim(iv,Nr, il)
  ENDDO
  Prim(u_pv,nr+1, il)   = -Prim(u_pv,nr, il)
  Prim(p12_pv,nr+1, il)   = -Prim(p12_pv,nr, il)
endif


! hydraulic jump

ELSEIF (COND_LIM == 7) THEN 
    if (rang .eq. 0) then  
       
!         Prim(h_pv,0, il)     = Prim(h_pv,1, il)
!           Prim(v_pv,0, il)     =  Prim(v_pv,1, il) !0.01d0*r(0) !-Prim(v_pv,1, il)
!           Prim(p11_pv,0, il)     = Prim(p11_pv,1, il)
!           Prim(p22_pv,0, il)     = Prim(p22_pv,1, il)
        
!        if ( (prim(h_pv, 0, il)) .le. dev) then 
!           Prim(u_pv,0, il)     = - Prim(u_pv,1, il)
!             Prim(p12_pv,0, il)     = - Prim(p12_pv,1, il)

!        else
!             Cd = pi/(pi + 2.d0) + 0.08d0*(Prim(h_pv,0, il) - dev)/dev
!             Prim(u_pv,0, il) = - (2.d0*(1.d0/(r(0)*Prim(h_pv,0, il))*(2.d0/3.d0*r(1)*Cd*dsqrt(2.d0*G*(Prim(h_pv,1, il) - dev)**3.d0))) &
!               &  - Prim(u_pv,1, il) )
            
!           !Prim(u_pv,0, il) =  - dsqrt(g*Prim(h_pv,1, il) + 3.d0*Prim(p11_pv,1, il) )
!           !Prim(u_pv,0, il) =  -(2.d0*dsqrt(g*Prim(h_pv,1, il) + 3.d0*Prim(p11_pv,1, il)  ) &
!             !  &  - Prim(u_pv,1, il) )
!         !  Prim(u_pv,0, il) =  Prim(u_pv,1, il)
!              Prim(p12_pv,0, il)     = - Prim(p12_pv,1, il)
!        endif
          

  do iv = 1,6
    Prim(iv,0, il)     = Prim(iv,1, il)
  enddo

      p11 = Prim(P11_pv,0, il); p22 = Prim(P22_pv,0, il);
       p12 = Prim(P12_pv,0, il); h = Prim(H_pv,0, il)
      u = Prim(U_pv,0, il); v = Prim(V_pv,0, il)
! !------------------------------------------------------------------------------------------------------------------------
   
      Prim(ein_pv,0, il) =InternalEn(h, p11, p22) !CONS(7,ir, il)/CONS(1,ir, il) - (u*u+v*v )/2.d0

      Prim(sound_ar_pv,0, il) = DSQRT(g*h + 3.d0*p11)
      Prim(sound_br_pv,0, il) = DSQRT(p11)  
      Prim(sound_ath_pv,0, il) = DSQRT(g*h + 3.d0*p22)
      Prim(sound_bth_pv,0, il) = DSQRT(p22)
      Prim(press_pv,0, il) = g*h**2.d0/2.d0 + p11*h

 

 elseif  (rang .eq. ncpu -1) then 

    PRIM(H_PV, nr+1, IL)= H_0
    PRIM(U_PV, nr+1, IL)= -q0_jump*(1.d0 + 1.d-2*dsin(16.d0*L(il)))&
    &/(2.d0*pi*R(nr+1)*PRIM(H_PV, nr+1, IL) ) !*(1.d0 + 1.d-2*dsin(16.d0*L(il)))
    PRIM(V_PV, nr+1, IL) = 0.d0 !1.d-3*r(nr+1) !  PRIM(V_PV, nr, IL)
    PRIM(P11_PV, nr+1, IL)=0.5d0*PHI2*Prim(H_pv, nr+1, il)**2.d0  !0.5d0*PHI2*Prim(H_pv, ir, il)**2.D0 
     !PRIM(P11_PV, nr+1, IL)= 0.d0

    PRIM(P12_PV, nr+1, IL)= PRIM(P12_PV, nr, IL) !0.d0

    PRIM(P22_PV, nr+1, IL)= 0.5d0*PHI2*Prim(H_pv, nr+1, il)**2.d0  !0.5d0*PHI2*Prim(H_pv, ir, il)**2.D0 
    ! PRIM(P22_PV, nr+1, IL)= 0.d0

endif
!---------------------cond lim 8 !!!!!!!!!!!!!!

ELSEIF (COND_LIM == 8) THEN 
    if (rang .eq. 0) then  
     PRIM(H_PV, 0, IL)= PRIM(H_PV, 1, IL)
     PRIM(U_PV, 0, IL)= -q0_jump/( R(0)*PRIM(H_PV, 0, IL) )
     PRIM(V_PV, 0, IL)= 0.d0 !PRIM(V_PV, 1, IL) 
     PRIM(P11_PV, 0, IL)= eps !PRIM(P11_PV, 1, IL)
     PRIM(P12_PV, 0, IL)= 0.d0 !PRIM(P12_PV, 1, IL)
     PRIM(P22_PV, 0, IL)=  eps !PRIM(P22_PV, 1, IL)



elseif  (rang .eq. ncpu -1) then 

     PRIM(H_PV, NR+1, IL)= PRIM(H_PV, NR, IL)
     PRIM(U_PV, NR+1, IL)= -q0_jump/( R(NR+1)*PRIM(H_PV, NR+1, IL) )
     PRIM(V_PV, NR+1, IL)= 0.d0 !PRIM(V_PV, NR, IL) 
     PRIM(P11_PV, NR+1, IL)= 0.d0 !PRIM(P11_PV, NR, IL)
     PRIM(P12_PV, NR+1, IL)= 0.d0 !PRIM(P12_PV, NR, IL)
     PRIM(P22_PV, NR+1, IL)=  0.d0 !PRIM(P22_PV, NR, IL)

endif

elseif (cond_lim == 10) then ! nonstationar turbulent nonspherical solution 

if (rang .eq. 0) then 
  Prim(h_pv, 0, il)=1.d0/(1.d0 + beta*time**2.d0) 
  Prim(u_pv, 0, il)=r(0)*beta*(beta*time*dcos(l(il))&
  & + dsin(l(il)))/(1.d0 + beta**2.d0*time**2.d0)
  Prim(v_pv, 0, il)=r(0)*beta*(beta*time*dsin(l(il))&
  & - dcos(l(il)))/(1.d0 + beta**2.d0*time**2.d0)

  Prim(p11_pv, 0, il)=( (lambda + gamma*beta**2.d0*time**2.d0)*dcos(l(il))*dcos(l(il)) + &
  & (gamma + lambda*beta**2.d0*time**2.d0)*dsin(l(il))*dsin(l(il)) + &
  & (lambda - gamma)*beta*time*dsin(2.d0*l(il)) )/(1.d0 + beta**2.d0*time**2.d0)**2.d0

  Prim(p12_pv, 0, il)=(lambda - gamma)*(2.d0*beta*time*dcos(2.d0*l(il)) &
    &+ (time**2.d0*beta**2.d0 - 1.d0)*dsin(2.d0*l(il)))/&
  & (2.d0*(1.d0 + beta**2.d0*time**2.d0)**2.d0)

  Prim(p22_pv, 0, il)=( (lambda + gamma*beta**2.d0*time**2.d0)*dsin(l(il))*dsin(l(il)) + &
  & (gamma + lambda*beta**2.d0*time**2.d0)*dcos(l(il))*dcos(l(il)) + &
  & (gamma - lambda)*beta*time*dsin(2.d0*l(il)) )/(1.d0 + beta**2.d0*time**2.d0)**2.d0
elseif  (rang .eq. ncpu -1) then

  Prim(h_pv, Nr+1, il)=1.d0/(1.d0 + beta*time**2.d0)
  Prim(u_pv, Nr+1, il)=r(nr+1)*beta*(beta*time*dcos(l(il))&
  & + dsin(l(il)))/(1.d0 + beta**2.d0*time**2.d0)
  Prim(v_pv, Nr+1, il)=r(nr+1)*beta*(beta*time*dsin(l(il))&
  & - dcos(l(il)))/(1.d0 + beta**2.d0*time**2.d0)

  Prim(p11_pv,Nr+1, il)=( (lambda + gamma*beta**2.d0*time**2.d0)*dcos(l(il))*dcos(l(il)) + &
  & (gamma + lambda*beta**2.d0*time**2.d0)*dsin(l(il))*dsin(l(il)) + &
  & (lambda - gamma)*beta*time*dsin(2.d0*l(il)) )/(1.d0 + beta**2.d0*time**2.d0)**2.d0

  Prim(p12_pv, Nr+1, il)= (lambda - gamma)*(2.d0*beta*time*dcos(2.d0*l(il))&
  & + (time**2.d0*beta**2.d0 - 1.d0)*dsin(2.d0*l(il)))/&
  & (2.d0*(1.d0 + beta**2.d0*time**2.d0)**2.d0)

  Prim(p22_pv, Nr+1, il)=( (lambda + gamma*beta**2.d0*time**2.d0)*dsin(l(il))*dsin(l(il)) + &
  & (gamma + lambda*beta**2.d0*time**2.d0)*dcos(l(il))*dcos(l(il)) + &
  & (gamma - lambda)*beta*time*dsin(2.d0*l(il)) )/(1.d0 + beta**2.d0*time**2.d0)**2.d0
endif
ENDIF

h = Prim(H_pv, 0, il); p11 = Prim(P11_pv, 0, il); p22 = Prim(P22_pv, 0, il)

  Prim(sound_ar_pv,0,il) = Sound_ar(h, p11) !DSQRT( g*Prim(H_pv, 0, il) + 3.d0*Prim(P11_pv, 0, il) )
  Prim(sound_br_pv,0,il) = Sound_br(p11) !DSQRT(Prim(P11_pv, 0, il) )
  Prim(sound_ath_pv,0,il) = Sound_atheta(h,p22) !DSQRT( g*Prim(H_pv, 0, il) + 3.d0*Prim(P22_pv, 0, il) )
  Prim(sound_bth_pv,0,il) = Sound_btheta(p22) !DSQRT(Prim(P22_pv, 0, il) )
  Prim(press_pv,0,il) = Pression(h,p11) !g*Prim(H_pv, 0, il)*Prim(H_pv, 0, il)/2.d0 + Prim(P11_pv, 0, il)*Prim(H_pv, 0, il)
  Prim(ein_pv,0,il) = InternalEn(h, p11, p22) ! ( Prim(H_pv, 0, il)*g + Prim(P11_pv, 0, il) + Prim(P22_pv, 0, il) )/2.d0


if (rang .eq. 0) then  
  h = Prim(H_pv, 0, il); 
  p11 = Prim(P11_pv, 0, il); p22 = Prim(P22_pv, 0, il)

  Prim(sound_ar_pv,0,il) = Sound_ar(h, p11) 
  Prim(sound_br_pv,0,il) = Sound_br(p11) 
  Prim(sound_ath_pv,0,il) = Sound_atheta(h,p22) 
  Prim(sound_bth_pv,0,il) = Sound_btheta(p22) 
  Prim(press_pv,0,il) = Pression(h,p11) 
  Prim(ein_pv,0,il) = InternalEn(h, p11, p22) 

elseif  (rang .eq. ncpu -1) then 
  h = Prim(H_pv,Nr+1, il); 
  p11 = Prim(P11_pv, Nr+1, il); p22 = Prim(P22_pv, Nr+1, il)

  Prim(sound_ar_pv,Nr+1,il) = Sound_ar(h, p11)
  Prim(sound_br_pv,Nr+1,il) = Sound_br(p11) 
  Prim(sound_ath_pv,Nr+1,il) =Sound_atheta(h,p22) 
  Prim(sound_bth_pv,Nr+1,il) = Sound_btheta(p22) 
  Prim(press_pv,Nr+1,il) = Pression(h,p11) 
  Prim(ein_pv,Nr+1,il) = InternalEn(h, p11, p22) 
endif

 enddo

do iv = 1, Nv_Prim
do ir = 1, Nr
  Prim(iv,ir, 0)    = Prim(iv,ir, Ntheta)
  Prim(iv,ir, Ntheta +1) = Prim(iv,ir, 1)  
enddo
end do

return
END SUBROUTINE
!************************************************************************************************************************************************************

Subroutine HLLC_r_sub1(prim,flux, dt, R, L, DR, dtheta, IT, time)
USE precisions
USE GlobalParam
IMPLICIT NONE
INTEGER :: ir, il, it
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), FLUX(7,0:Nr,0:Ntheta)
REAL (KIND = DP) :: R(0:Nr+1), L(1:Ntheta)
real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
real(kind=dp) :: EL, ER,p22l,p22r,vr,vl
real(kind=dp) :: pstar,ustar,Estar,hstar, p12star, p22star
real(kind=dp) :: dt, dt2, time, dr, dtheta
real(kind=dp) :: u0, p0, v0, e0, p110, p120, p220, h0, cd
logical ::  machin, bidule
dt2 = 0.5d0*dt


DO il=1,ntheta
DO ir=0,nr

    !! etat gauche et droite
    ul=Prim(U_pv,ir,il);  ur=Prim(U_pv,ir+1, il)
    hl=Prim(h_pv,ir,il);  hr=Prim(h_pv,ir+1, il)
    vl=prim(v_pv,ir,il);  vr=Prim(v_pv,ir+1, il)

    p11l=prim(p11_pv,ir,il);  p11r=Prim(p11_pv,ir+1, il)
    p12l=prim(p12_pv,ir,il);  p12r=Prim(p12_pv,ir+1, il)
    p22l=prim(p22_pv,ir,il);  p22r=Prim(p22_pv,ir+1, il)

    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !! calcul des pression
    pl=g*hl*hl*0.5d0+hl*p11l
    pr=g*hr*hr*0.5d0+hr*p11r
    !! calcul des vitesses du son
    cl=dsqrt(g*hl+3.d0*p11l); cr=dsqrt(g*hr+3.d0*p11r)
    ! davis
    sr=ur+cr; if(ul+cl> sr) sr = ul+cl
    sr = dmax1(ul+cl, ur+cr)
    !if (dabs(sr) .le. 1.d-8) sr = 0.d0
    Sl=ul-cl; if (ur-cr < sl) sl = ur - cr
    !sr = - sl
    sl = dmin1(ul-cl, ur - cr)
   ! if (dabs(sl) .le. 1.d-8) sl = 0.d0
    ! etat star 
    ml=hl*(ul-sl)
    mr=hr*(ur-sr)

    ustar = (pl - pr + ml*ul - mr*ur)/(ml - mr)
    pstar= ( (ur-ul)*mr*ml-mr*pl+ml*pr )/(ml-mr)
!     if (dabs(ustar) .le. 1.d-8) then
!       ustar = 0.d0
!     endif

!     IF (cond_lim == 7) THEN
!       IF (ir == 0 .and. rang == 0) THEN
!         if (Prim(h_pv,1,il) .le. dev ) then
!           ustar = 0.d0
!         else
!          ustar = - dsqrt(Prim(U_pv,1, il) + 3.d0*Prim(p11_pv,1, il))
!        endif
!       ENDIF
!      ENDIF
!----------------------------------------------------------------------------------------------------
 !CALCUL LES FLUX

IF (ustar.ge.0.d0) THEN
   IF (sl.ge.0.d0) THEN
  !!etat gauche
   flux(1,ir,il)=hl*ul
   flux(2,ir,il)=hl*ul*ul+pl
   flux(3,ir,il)=hl*ul*vl
   flux(4,ir,il)=0.d0 !hl*p11l*ul !! equation inutile pour le cas r
   flux(5,ir,il)=p12l*ul
   flux(6,ir,il)=hl*p22l*ul
   flux(7,ir,il)=hl*El*ul+pl*ul

   ELSE !if (sl .lt. 0.d0 .or. bidule) then 

   !! etat star gauche
   hstar=ml/(ustar-sl)
   Estar=El+(pl*ul- pstar*ustar)/ml
   p12star = p12l*hstar/hl
  ! p22star = hl*p22l*(ul - sl)/(hstar*(ustar - sl))
   !! remplissage des flux
   flux(1,ir,il)=hstar*ustar
   flux(2,ir,il)=hstar*ustar*ustar+pstar
   flux(3,ir,il)=hstar*ustar*vl  
   flux(4,ir,il)=0.d0 !hstar*p11l*ustar ! equation inutile pour le cas r
   flux(5,ir,il)=p12star*ustar
   flux(6,ir,il)=hstar*p22l*ustar
   flux(7,ir,il)=hstar*Estar*ustar+pstar*ustar
  ENDIF
ELSE
   IF (sr.ge.0.d0) THEN
   
  !!etat droit etoile
   hstar=mr/(ustar-sr)
   Estar=Er+(pr*ur-pstar*ustar)/mr
   p12star = p12r*hstar/hr
   !p22star = hr*p22r*(ur - sr)/(hstar*(ustar - sr))
   !remplissage des flux
   flux(1,ir,il)=hstar*ustar
   flux(2,ir,il)=hstar*ustar*ustar+pstar   
   flux(3,ir,il)=hstar*ustar*vr
   flux(4,ir,il)=0.d0 !hstar*p11r*ustar ! equation inutile pour le cas r
   flux(5,ir,il)=p12star*ustar
   flux(6,ir,il)=hstar*p22r*ustar
   flux(7,ir,il)=hstar*Estar*ustar+pstar*ustar 
ELSE
  !!etat droit 
   flux(1,ir,il)=hr*ur
   flux(2,ir,il)=hr*ur*ur+pr
   flux(3,ir,il)=hr*ur*vr
   flux(4,ir,il)=0.d0 !hr*p11r*ur !! equation inutile pour le cas x
   flux(5,ir,il)=p12r*ur
   flux(6,ir,il)=hr*p22r*ur
   flux(7,ir,il)=hr*Er*ur+pr*ur
ENDIF


ENDIF



ENDDO
! IF (cond_lim == 7 .and. rang ==  0) THEN

!   Cd = pi/(pi + 2.d0) + 0.08d0*(Prim(h_pv,1, il) - dev)/dev
 
!   !Prim(u_pv,0, il) =  -( 2.d0*(1.d0/(r(0)*Prim(h_pv,0, il))*(2.d0/3.d0*Cd*dsqrt(2.d0*G*(Prim(h_pv,0, il) - dev)**3.d0))) &
            
!   ! u0 = - dsqrt(Prim(U_pv,1, il) + 3.d0*Prim(p11_pv,1, il))
!  h0 = Prim(h_pv,1, il) !-  ((2.d0/3.d0*Cd*r(1)*dsqrt(2.d0*G*(Prim(h_pv,0, il) - dev)**3.d0)) - r(1)*Prim(h_pv,1, il)*&
!    ! &                                         Prim(u_pv,1, il))/(r(1)*dsqrt(g*Prim(h_pv,1, il) + 3.d0*Prim(p11_pv,1, il)))
!   !h0 = Prim(h_pv,1, il)


! if (h0 <= dev) then
!   u0 = 0.d0
! else 
!   u0 = - dsqrt(g*Prim(h_pv,1, il) + 3.d0*Prim(p11_pv,1, il) ) !1.d0/(r(1)*Prim(h_pv,1, il))*(2.d0/3.d0*Cd*r(1)*dsqrt(2.d0*G*(Prim(h_pv,1, il) - dev)**3.d0))

! endif

!  ! print*, '(2.d0/3.d0*Cd*dsqrt(2.d0*G*(Prim(h_pv,0, il) - dev)**3.d0))=', (2.d0/3.d0*Cd*dsqrt(2.d0*G*(Prim(h_pv,0, il) - dev)**3.d0))
!   !h0 = Prim(h_pv,1, il) -  ((2.d0/3.d0*Cd*dsqrt(2.d0*G*(Prim(h_pv,0, il) - dev)**3.d0)) - r(1)*Prim(h_pv,1, il)*&
!   !  &                                         Prim(u_pv,1, il))/(r(1)*dsqrt(Prim(U_pv,1, il) + 3.d0*Prim(p11_pv,1, il)))
   
!    v0 = Prim(v_pv,1, il)
!    p120 = Prim(p12_pv,1, il)
!    p110 = Prim(p11_pv,1, il)
!    p220 = Prim(p22_pv,1, il)
!    p0 = g*h_0*h_0*0.5d0 + h_0*p110
!    E0 = (u0*u0+v0*v0+g*h0+p110+p220)*0.5D0

!    flux(1,0,il)=h0*u0
!    flux(2,0,il)=h_0*u0*u0+p0
!    flux(3,0,il)=h0*u0*v0
!    flux(4,0,il)=0.d0 !hr*p11r*ur !! equation inutile pour ce cas
!    flux(5,0,il)=p120*u0
!    flux(6,0,il)=h0*p220*u0
!    flux(7,0,il)=h0*E0*u0+p0*u0
!   ENDIF

ENDDO
RETURN
END subroutine HLLC_r_sub1
 !************************************************************************************************************************************************************


!************************************************************************************************************************************************************

Subroutine HLLC_r_sub2(prim,flux, dt, R,L, DR, dtheta, IT, time)
USE precisions
use GlobalParam
IMPLICIT NONE
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), FLUX(7,0:Nr,0:Ntheta)
REAL (KIND = DP) :: R(0:nr+1), L(1:ntheta)
real(kind=dp) :: ul, ur, hl, hr, p11l, p11r, p12r, p12l, ml, mr, sr, sl, p22l, p22r
real(kind=dp) :: EL, ER,vr,vl
real(kind=dp) :: vstar,Estar,hp12star
real(kind=dp) :: dt, dr, time, dtheta
INTEGER :: ir,il, it

do ir=0,nr
do il=1,ntheta

ul=prim(U_pv,ir,il);     ur=prim(U_pv,ir+1, il)
hl=prim(h_pv,ir,il);     hr=prim(h_pv,ir+1, il)
vl=prim(v_pv,ir,il);     vr=prim(v_pv,ir+1, il)
p11l=prim(p11_pv,ir,il); p11r=prim(p11_pv,ir+1, il)
p12l=prim(p12_pv,ir,il); p12r=prim(p12_pv,ir+1, il)
p22l=prim(p22_pv,ir,il); p22r=prim(p22_pv,ir+1, il)

!! calcul des énergie
El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0

! sr=dmax1(dsqrt(p11l),dsqrt(p11r)) !
!sl=dmin1(-dsqrt(p11l),-dsqrt(p11r))!

sl= - dsqrt(p11l)
 sr=dsqrt(p11r)

vstar         = (hl*(p12l- sl*vl) - hr*(p12r-sr*vr))/(sr*hr - hl*sl)
hp12star      = (hr*hl)/(hr*sr-hl*sl)*(sr*p12l-sl*p12r+sl*sr*(vr-vl))

Flux(1,ir,il) = 0.d0
Flux(2,ir,il) = 0.d0
flux(3,ir,il) = hp12star
Flux(4,ir,il) = 0.d0
Flux(5,ir,il) = vstar !! non conservatif
flux(6,ir,il) = 2.d0*hp12star*vstar !!inutile
flux(7,ir,il) = hp12star*vstar

ENDDO
ENDDO
return
end subroutine HLLC_r_sub2

!************************************************************************************************************************************************************

SUBROUTINE PutonScreen(rang)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  REAL (KIND = DP) :: MinVp(Nv_Prim), MaxVp(Nv_Prim)
   REAL (KIND = DP) :: tamponminvp(nv_prim), tamponmaxvp(nv_prim)
  integer:: iv, rang

  DO iv = 1, Nv_Prim
    TAMPONMinVp(iv)     = MINVAL(Prim(iv,1:Nr,1:Ntheta))
    TAMPONMaxVp(iv)     = MAXVAL(Prim(iv,1:Nr,1:Ntheta))
  ENDDO
 
  DO iv = 1, Nv_Prim
    CALL MPI_ALLREDUCE(TAMPONMinVp(iv), MinVp(iv),&
    & numb, MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(TAMPONMaxVp(iv), MaxVp(iv), numb,&
    & MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
  END DO

  if (rang == 0) then 
      WRITE(argunit,*)
      WRITE(argunit,'(65(":"))') 
      WRITE(argunit,'("::  ",47("-"),12(" "),"::")')
      WRITE(argunit,'(":: | Time  = ", E10.3, 1x, &
           &   "  ||     kt =  ",I9,3x,"| ",10x,"::")' ) TIME, it
      WRITE(argunit,'(":: |    Dt = ", E10.3,1x, &
           &   "  ||    Cfl =  ",E10.3,"  | ",10x,"::")' ) Dt , CFL
      WRITE(argunit,'("::  ",47("-"),12(" "),"::")')
      WRITE(argunit,'("::",61(" "),"::")') 
      WRITE(argunit,'("::   ",12x, "    ",  2(3x,A12), 11x," ::" )') " Minimum ", " Maximum "

      DO iv = 1, Nv_Prim
         WRITE(argunit,'("::   ",A12, " => ",  2(3x,E12.5), 11x," ::" )') &
              & TRIM(NamesOfPrimVar(iv)), MinVp(iv) , MaxVp(iv)
      END DO

      WRITE(argunit,'("::   ",A12, " => ",  2(3x,E12.5), 11x," ::" )')  
      WRITE(argunit,'(65(":"))') 
  endif
  END SUBROUTINE PutonScreen
!****************************************************************************************************************************************************************

SUBROUTINE Ecriture_donnees(R, L, Prim, time, DR, Dtheta, bottom,rang, Nfile, reste)
USE precisions
USE GlobalParam
IMPLICIT NONE
INTEGER :: ir, il, im, rang, Nfile, reste
REAL (KIND = DP) :: R(0:Nr+1), L(1:Ntheta),XX(0:Nr), YY(0:Ntheta), h_new(1:Nr,1:Ntheta)
REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), bottom(1:Nr)
REAL (KIND = DP) :: time, DR, Dl, dtheta, h, p0, U, V
REAL (KIND = DP) :: P11, P12, P22, alphaTB(1:Nr,1:Ntheta), traceP
REAL (KIND = DP) :: grad(1:2,1:Nr,1:Ntheta), norm_grad(1:Nr,1:Ntheta), shliren(1:Nr,1:Ntheta)
REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:):: err_h, err_u,err_p11
REAL (KIND = DP) :: tampon_grad_min, tampon_grad_max, grad_min, grad_max 
REAL (KIND = DP) :: error_h,error_u,error_p11
REAL (KIND = DP) ::  TamponError_h, TamponError_u, TamponError_p11
CHARACTER(LEN=4) :: NB, Zero="0000"
CHARACTER*4 :: numrg

allocate(err_h(1:Nr,1:Ntheta), err_u(1:Nr,1:Ntheta),err_p11(1:Nr,1:Ntheta))
error_h=0.d0; error_u =0.d0;error_p11=0.d0
err_h(1:Nr,1:Ntheta)=0.d0; err_u(1:Nr,1:Ntheta)=0.d0;err_p11(1:Nr,1:Ntheta)=0.d0
grad= 0.d0; norm_grad = 0.d0; shliren = 0.d0
alphaTB(:,:) = 0.d0
 WRITE(numrg,'(i4.4)') rang


 !print*, './resu/test_'//numrg//'.out'
 OPEN(MyUnit+rang,FILE = './resu/test_'//numrg//'.out')
 !if (logik) WRITE(MyUnit+rang,*) Nr, Ntheta, Nfile
!OPEN(MyUnit+rang+2,FILE = './resu/alpha_trP_'//numrg//'.out')

Do ir = 1, Nr
DO il = 1, Ntheta
    
  if (ntheta <=  1) then
    dl = 0.d0
    L(il) = 0.d0
  else
    Dl = Dtheta !R(ir)*
    L(il) = 0.5D0*Dl + (il-1)*Dl
  endif

  IF (COND_LIM == 7 ) THEN


         ! BOTTOM(IR) =   (5.6d0*0.01d0)**2.d0*(1.d0/((Rplus)) - 1.d0/r(ir))
          
          IF (R(IR) - RMINUS <= 2.D0*L1) THEN 
          BOTTOM(IR) = DEV*((R(IR) - RMINUS - L1)**2.D0 - L1**2.D0 )**2.D0/L1**4.D0 
         ELSE
          BOTTOM(IR) =   (r(ir) - 2.D0*L1 -  RMINUS)*Dtan(angle) !
         ENDIF
  else
    BOTTOM(IR) =   0.d0
  ENDIF

   ! bottom = 0.d0
    h = Prim(H_pv,ir, il)
    p0 = (0.5d0*PHI2*Prim(H_pv, ir, il)**2.d0)**2.d0

      P11 = Prim(P11_pv,ir, il);
      p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)

      traceP = Prim(P11_pv,ir, il)+ Prim(P22_pv,ir, il) 
      u = Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
       alphaTB(ir,il) = disscoeff*h**2.d0*dabs(traceP/H**2.d0 - phi2)/(traceP**2.d0 + 1.d-6)
    !  WRITE(MyUnit+rang+2,'(2(E20.13,1X))') Prim(P11_pv,ir, il)+ Prim(P22_pv,ir, il) ,  alphaTB(ir,il)


    WRITE(MyUnit+rang,'(8(E20.13,1X))') R(ir)*DCOS(L(1) ), R(ir)*DSIN( L(1) ), &
    &Prim(H_pv,ir, 1) + BOTTOM(ir), &
    & Prim(v_pv,ir, 1), Prim(u_pv, ir, 1), & 
    &  (prim(p11_pv,ir,1)*prim(p22_pv,ir,1) - prim(p12_pv,ir,1)**2.d0)/(Prim(H_pv, ir, 1)**2.d0),&
    & dabs(Prim(U_pv,ir, 1))/dsqrt(g*Prim(H_pv,ir, 1) + 3.d0*prim(p11_pv,ir,1)) , &
    & r(ir)*(h)*Prim(U_pv,ir, 1)!(prim(p11_pv,ir,il)*Prim(p22_pv,ir, il) -Prim(p12_pv,ir, il)**2.d0)/Prim(H_pv,ir, il)**2.d0
    END DO
  END DO
  WRITE(MyUnit+rang,*)


 IF (COND_LIM == 3) THEN 
   OPEN(MyUnit+1+rang, FILE = './resu/ANALITYCAL_SOLUTION_WATER_GLASS_'//numrg//'.out')   
   Do ir = 1, Nr
   DO il = 1, Ntheta

     WRITE(MyUnit+1+rang,'(8(E20.13,1x))') R(ir)*DCOS( L(il) ), R(ir)*DSIN( L(il) ),&
     & H_0 + omega*omega*R(ir)*R(ir)/(2.d0*g), 0.d0, &
     &  omega*R(ir), 0.d0,0.d0, 0.d0
   END DO
   END DO
   WRITE(MyUnit+1+rang,*)

elseif (COND_LIM == 4) THEN 
   OPEN(MyUnit+1+rang, FILE = './resu/err_NonStatSol_'//numrg//'.out')   


   Do ir = 1, Nr
   DO il = 1, Ntheta
     Dl = Dtheta 
    L(il) = 0.5D0*Dl + (il-1)*Dl
        h   = H_0/(1.d0 + time)**2.d0
        u = r(ir)/(1.d0 + time)
        p11 = P_0/(1.d0 + time)**2.d0
   
        
          if (dabs(u).le.1.d-8) then 
           u=0.d0; 
          endif
  
           err_h(ir,il) = dabs( (h - Prim(H_pv,ir,il))/h)
           err_u(ir,il) = dabs( (u - Prim(u_pv,ir,il) )/u)
           err_p11(ir,il) = dabs( (p11 - Prim(p11_pv,ir,il))/p11)
              END DO
   END DO

       TamponError_h = maxval(err_h(:, :)) ! SUM(err_h(:, :))/(nr*ntheta) 
       TamponError_u =maxval(err_u(:, :))  !SUM(err_u(:, :))/(nr*ntheta) !
       TamponError_p11 = maxval(err_p11(:, :)) ! SUM(err_p11(:, :))/(nr*ntheta) !

      CALL MPI_ALLREDUCE(TamponError_h, error_h, NUMB,MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
      CALL MPI_ALLREDUCE(TamponError_u, error_u, NUMB,MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
      CALL MPI_ALLREDUCE(TamponError_p11, error_p11, NUMB,MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
       WRITE(MyUnit+1+rang,'(4(E20.13,1X))') DLOG(Dr*Dtheta/(Lr*Ltheta) ),&
       &  DLOG(error_h),DLOG(error_u), DLOG(error_p11)




 !ELSEIF (COND_LIM == 10) THEN 

!   OPEN(MyUnit+1+rang,FILE = './resu/ANALITYCAL_TURB_SOLUTION'//numrg//'.out')   
!   Do ir = 1, Nr
!   DO il = 1, Ntheta
!   H=1.d0/(1.d0 + beta*time**2.d0)
!   U=r(IR)*beta*(beta*time*dcos(l(il)) + dsin(l(il)))/(1.d0 + beta**2.d0*time**2.d0)
!   V=r(IR)*beta*(beta*time*dsin(l(il)) - dcos(l(il)))/(1.d0 + beta**2.d0*time**2.d0)

!   P11=( (lambda + gamma*beta**2.d0*time**2.d0)*dcos(l(il))*dcos(l(il)) + &
!   & (gamma + lambda*beta**2.d0*time**2.d0)*dsin(l(il))*dsin(l(il)) + &
!   & (lambda - gamma)*beta*time*dsin(2.d0*l(il)) )/(1.d0 + beta**2.d0*time**2.d0)**2.d0

!   P12= (lambda - gamma)*(2.d0*beta*time*dcos(2.d0*l(il)) + (time**2.d0*beta**2.d0 - 1.d0)*dsin(2.d0*l(il)))/&
!   & (2.d0*(1.d0 + beta**2.d0*time**2.d0)**2.d0)

!   P22=( (lambda + gamma*beta**2.d0*time**2.d0)*dsin(l(il))*dsin(l(il)) + &
!   & (gamma + lambda*beta**2.d0*time**2.d0)*dcos(l(il))*dcos(l(il)) + &
!   & (gamma - lambda)*beta*time*dsin(2.d0*l(il)) )/(1.d0 + beta**2.d0*time**2.d0)**2.d0


!     WRITE(MyUnit+1+rang,'(8(E20.13,1x))') R(ir)*DCOS( L(il)/R(ir) ), R(ir)*DSIN( L(il)/R(ir) ), H, U, V, P11, P12, P22
!   END DO
!   END DO
!   WRITE(MyUnit+1+rang,*)

 ENDIF
 
if (rang .eq. 0) R(0) = Rminus - dr*0.5D0
if (rang .eq. ncpu - 1) R(nr+1) = Rplus + dr*0.5d0

do ir = 1, Nr
  if (rang .lt. (reste)) then
    R(ir) =Rminus + 0.5D0*DR + (ir-1)*DR + nr*rang*dr 
  else 
    R(ir) =Rminus + 0.5D0*DR + (ir-1)*DR+ (reste)*(nr+1)*dr + nr*(rang - reste)*dr
  endif
  do il = 1, Ntheta
  Dl = Dtheta ! R(ir)*
  L(il) = 0.5D0*Dl + (il-1)*Dl
  grad(1,ir, il) = (Prim(H_pv,ir+1, il) - Prim(H_pv,ir-1, il) )/(2.d0*dr)
  grad(2,ir, il) = (Prim(H_pv,ir, il+1) - Prim(H_pv,ir, il-1) )/(2.d0*dl)
  norm_grad(ir, il) = dsqrt(grad(1,ir, il)**2.d0 + grad(2,ir, il)**2.d0)

enddo
enddo
  numb = 1
  tampon_grad_max= maxval(norm_grad(:, :))
  tampon_grad_min = minval(norm_grad(:, :))
  CALL MPI_ALLREDUCE(tampon_grad_max, grad_max,&
    &NUMB, MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
  CALL MPI_ALLREDUCE(tampon_grad_min, grad_min,&
  & NUMB, MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,code)
  do ir = 1, Nr; do il = 1, Ntheta
  shliren(ir, il) = dexp(- (norm_grad(ir, il) )/(grad_max))
  !print*,  ' shliren(ir, il) = ', shliren(ir, il), 'grad_max', grad_max
  enddo;enddo

isave = isave + 1

! *****************************************************************************************************************************************************************
 if (ntheta > 1 ) then 

      do ir = 0, Nr

         if (rang .lt. (reste)) then
          XX(ir) = Rminus +  (ir)*DR + nr*rang*dr 
         else 
           XX(ir) = Rminus + (ir)*DR+ (reste)*(nr+1)*dr + nr*(rang - reste)*dr
         endif

       !XX(ir) = rminus + dfloat(ir)*dr + nr*rang*dr 

      ENDDO
     do il = 0, Ntheta
       YY(il) = 0.d0 + dfloat(il)*dtheta
     ENDDO 


     WRITE(unit=NB, fmt="(I4)") isave
     NB    = ADJUSTL(NB)
     im    = LEN_TRIM(NB) 
     print*, "./resu/test_"//numrg//"_"//Zero(1:4-im)//TRIM(NB)//".vtk" 
     OPEN(UNIT=MyUnit+rang+1, FILE="./resu/test_"//numrg//"_"//Zero(1:4-im)//TRIM(NB)//".vtk" )
     WRITE(MyUnit+rang+1,'(''# vtk DataFile Version 3.0'')')
     WRITE(MyUnit+rang+1,'(''vtk output'')')
     WRITE(MyUnit+rang+1,'(''ASCII'')')
     WRITE(MyUnit+rang+1,'(''DATASET STRUCTURED_GRID'')')
     WRITE(MyUnit+rang+1,FMT='(''DIMENSIONS'',I8,I8,I8)') Nr+1, Ntheta+1, 1
  !   WRITE(MyUnit,FMT='(''ORIGIN'',3(E11.4,1x))') rminus, 0.d0, 0.d0
  !   WRITE(MyUnit,FMT='(''SPACING'',3(E11.4,1x))') dr,dtheta, 0.0001d0
      WRITE(MyUnit+rang+1,FMT='(''POINTS '',I8, '' float '')') (Nr+1)*(Ntheta+1) 
      do il = 0, Ntheta
      do ir = 0, Nr
       WRITE(MyUnit+rang+1,'(3(E11.4,1x))')  XX(ir)*dCos(YY(il)),XX(ir)*dsin(YY(il)), 0.d0
     ENDDO ; ENDDO
       
     WRITE(MyUnit+rang+1,*) ' '
     WRITE(MyUnit+rang+1,FMT='(''CELL_DATA '',I9)') (Nr)*(Ntheta)*1
     WRITE(MyUnit+rang+1,*) ' '

     WRITE(MyUnit+rang+1,FMT='(''SCALARS '',A6, '' float 1'')') 'depth'
     WRITE(MyUnit+rang+1,'(''LOOKUP_TABLE default'')')     
    DO il=1,Ntheta ; DO ir=1,Nr
       WRITE(MyUnit+rang+1,'(G11.4)')  Prim(H_pv,ir, il) + BOTTOM(ir)        
     ENDDO ; ENDDO 

     WRITE(MyUnit+rang+1,FMT='(''VECTORS '',A12, '' float'')') 'Vitesse(m/s)'
     DO il=1,Ntheta ; DO ir=1,Nr
       WRITE(MyUnit+rang+1,'(3(E11.4,1x))')  Prim(U_pv,ir, il), Prim(v_pv,ir, il), 0.d0
     ENDDO ; ENDDO 
  
     WRITE(MyUnit+rang+1,FMT='(''SCALARS '',A6, '' float 1'')') 'Entropy'
     WRITE(MyUnit+rang+1,'(''LOOKUP_TABLE default'')')     
    DO il=1,Ntheta ; DO ir=1,Nr
      P11 = Prim(P11_pv,ir, il); h = Prim(h_pv,ir, il)
      p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
       WRITE(MyUnit+rang+1,'(G11.4)') (p11*p22 - p12**2.d0)/((h)**2.d0)
     ENDDO ; ENDDO 

     WRITE(MyUnit+rang+1,FMT='(''SCALARS '',A6, '' float 1'')') 'V(m/s)'
     WRITE(MyUnit+rang+1,'(''LOOKUP_TABLE default'')')     
    DO il=1,Ntheta ; DO ir=1,Nr
       WRITE(MyUnit+rang+1,'(G11.4)') Prim(v_pv,ir, il)   
     ENDDO ; ENDDO 
 

     WRITE(MyUnit+rang+1,FMT='(''SCALARS '',A6, '' float 1'')') 'shliren'
     WRITE(MyUnit+rang+1,'(''LOOKUP_TABLE default'')')     
    DO il=1,Ntheta ; DO ir=1,Nr
       WRITE(MyUnit+rang+1,'(G11.4)')  shliren(ir, il)  
     ENDDO ; ENDDO 

     WRITE(MyUnit+rang+1,*) ' '
  close(MyUnit+rang+1)


    
!     WRITE(6,*) " FILE = ", "./resu/test1_test3_pert5%_q2L_16pomp_"//numrg//"_"//Zero(1:4-im)//TRIM(NB)//".tp"
!     OPEN(UNIT=MyUnit+2+rang, FILE="./resu/test3_pert5%_q2L_16pomp_"//numrg//"_"//Zero(1:4-im)//TRIM(NB)//".tp" )
!     WRITE(MyUnit+2+rang,'(A)') 'TITLE="This is a title"'  
!     WRITE(MyUnit+2+rang,'(A)') 'VARIABLES= "X", "Y" , "H","U", "V", "P11","P12", "P22" '
!     WRITE(MyUnit+2+rang,*) 'ZONE I=', Nr,', J=', Ntheta, 'DATAPACKING=POINT' 
   
!     Do ir = 1, Nr  
!     DO il = 1, Ntheta
!           WRITE (MyUnit+2+rang,'(8(E16.8,1x))') R(ir)*DCOS( L(il)  ), R(ir)*DSIN( L(il)  ), Prim(H_pv,ir, il)+ BOTTOM(ir) , &
!           & Prim(U_pv,ir, il), Prim(V_pv,ir, il), prim(p11_pv,ir,il), Prim(p12_pv,ir, il), prim(p22_pv,ir,il)
!     END DO
!     END DO
!  ! *****************************************************************************************************************************************************************

!    WRITE(6,*) " FILE = ", "./resu/Radial_DamBreak_WithFriction_"//numrg//"_"//Zero(1:4-im)//TRIM(NB)//".csv"  
!    OPEN(UNIT=MyUnit+3+rang, FILE="./resu/Radial_DamBreak_WithFriction_"//numrg//"_"//Zero(1:4-im)//TRIM(NB)//".csv"  )
!    WRITE(MyUnit+3+rang,'(A)') '"X", "Y",  "H"'  
!    Do ir = 1, Nr 
!    DO il = 1, Ntheta  
    
!              WRITE (MyUnit+3+rang ,*)  R(ir)*DCOS(L(il) ), ',' , R(ir)*DSIN( L(il)), ',',  Prim(H_pv,ir, il)!+ BOTTOM(ir)
!    END DO
!    END DO

endif

deallocate(err_h, err_u,err_p11)
return
END SUBROUTINE
!************************************************************************************************************************************************************

SUBROUTINE euler_method( CONS, Prim,R, dr, DT, it)
  USE precisions
  USE GlobalParam
  implicit none 
  INTEGER :: ir, il, IT,k, j
  REAL (KIND = DP) :: DT, traceP, H, p11, p22, p12, u, v,Q, dr, modU, detP
  REAL (KIND = DP) :: CONS(7,1:Nr, 1:Ntheta), R(0:Nr+1)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1)
  REAL (KIND = DP) :: fracp11,fracp22, erreur, nu
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:) :: TS_phys, D
  allocate(TS_phys(7), D(3))
  TS_phys = 0.d0; D = 0.d0; nu = 1.d-6
  traceP = 0.d0; H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0
!------------------------------------------------------------------------------------------------------------------------

  DO  ir = 1, Nr
  DO  il  = 1, Ntheta  
!------------------------------------------------------------------------------------------------------------------------ 
 !TERME SOURCE
H = Prim(H_pv,ir, il)
P11 = Prim(P11_pv,ir, il)
p22 = Prim(P22_pv,ir, il)
p12 = Prim(P12_pv,ir, il)

traceP = Prim(P11_pv,ir, il)+ Prim(P22_pv,ir, il) 

u = Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22)

alpha =  disscoeff*((traceP/(H**2.d0) - phi2))/( traceP/(H**2.d0)  ) !disscoeff*((traceP/(H**2.d0) - phi2))/( traceP**2.d0/(H**2.d0)  ) !
alpha = dmax1(alpha, 0.d0)

modU = dsqrt(u**2.d0 + v**2.d0)

Q =  alpha*traceP**(1.5d0) !alpha*traceP*(modU) !alpha*traceP*(modU)**3.d0 !*(H_0**3.d0/h**3.d0)

D(1) =   -(2.d0)*alpha*(traceP**(0.5d0))*p11 !-2.d0*alpha*(modU)*p11 !-(2.d0)*alpha*(traceP**(0.5d0))*p11 !-
D(2) =  -(2.d0/h)*alpha*(traceP**(0.5d0))*p12!-(2.d0/h)*alpha*(modU)*(p12) !
D(3) =    -(2.d0)*alpha*(traceP**(0.5d0))*p22 !-2.d0*alpha*(modU)*(p22) !

!TS_phys(2) =  -(5.6d0*0.01d0)**2.d0*g/R(ir)*H-&
! & frottcoeff*R(ir)*dsqrt( u**2.d0 + v**2.d0 )*u

TS_phys(2) =  -g*R(ir)*Dtan(angle)*H  -frottcoeff*R(ir)*dsqrt(u**2.d0 + v**2.d0)*u  !-nu*R(ir)*u/h!
IF ( (R(IR) - RMINUS) <= 2.D0*L1) THEN
TS_phys(2) =  -(4.D0*DEV/L1**4.D0*( (R(IR) -&
& RMINUS - L1)**2.D0 - L1**2.D0)*(R(IR) - RMINUS - L1))*g*R(ir)*H-&
& frottcoeff*R(ir)*dsqrt( u**2.d0 + v**2.d0 )*u  !nu*R(ir)*u/h !
ENDIF

TS_phys(3) = -frottcoeff*R(ir)*dsqrt(u**2.d0 + v**2.d0 )*v  ! -nu*R(ir)*v/h !
TS_phys(4) =  R(ir)*D(1) 
TS_phys(5) =  R(ir)*D(2)
TS_phys(6) =  R(ir)*D(3) 

! TS_phys(7) = -(5.6d0*0.01d0)**2.d0*g/R(ir)*H*U   &
!  & -frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - R(ir)*Q   

TS_phys(7) = -g*R(ir)*Dtan(angle)*H*U - R(ir)*Q  &
& - frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**3.d0 !-nu*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**2.d0/h !

IF ( (R(IR) - RMINUS) <= 2.D0*L1) THEN
 TS_phys(7) = -(4.D0*DEV/L1**4.D0*( (R(IR) - RMINUS - L1)**2.D0 &
  &- L1**2.D0)*(R(IR) - RMINUS - L1))*g*R(ir)*H*U   &
 &- R(ir)*Q  -frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) !-nu*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**2.d0/h !
ENDIF  
!------------------------------------------------------------------------------------------------------------------------
DO k=1,7
  CONS(k,ir, il) = CONS(k,ir, il) + DT*TS_phys(k) 
ENDDO 

  ENDDO
  ENDDO
   CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
  return
  deallocate(TS_phys, D)
  END SUBROUTINE
!************************************************************************************************************************************************************

SUBROUTINE rk2( CONS, Prim,R, dr, DT, it)
  USE precisions
  USE GlobalParam
  implicit none 
  INTEGER :: ir, il, IT,k
  REAL (KIND = DP) :: DT, traceP, H, p11, p22, p12, u, v,Q, dr, modU, detP
  REAL (KIND = DP) :: CONS(7,1:Nr, 1:Ntheta), R(0:Nr+1)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1)
  REAL (KIND = DP) :: fracp11,fracp22, erreur
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:) :: TS_phys, D, TS
  REAL (KIND = DP), ALLOCATABLE:: cons1(:, :, :)
  allocate(TS_phys(7, 1:Nr, 1:Ntheta ), D(3, 1:Nr, 1:Ntheta), &
    &TS(7, 1:Nr, 1:Ntheta), cons1(7, 1:Nr, 1:Ntheta))
  TS_phys = 0.d0; D = 0.d0; TS = 0.d0; CONS1 = 0.d0
  traceP = 0.d0; H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0
!------------------------------------------------------------------------------------------------------------------------
  DO  ir = 1, Nr 
  DO il  = 1, Ntheta  
!------------------------------------------------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------------------------------------------------
   !TERME SOURCE
  H = Prim(H_pv,ir, il)
  P11 = Prim(P11_pv,ir, il);
  p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)

  traceP = Prim(P11_pv,ir, il)+ Prim(P22_pv,ir, il) 

  u = Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
  fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);

  alpha = disscoeff*h**2.d0*dabs(traceP/H**2.d0 - phi2)&
  &/(traceP**2.d0 +  phi2*g*H_0*h**2.d0*eps)

  modU = dsqrt(u**2.d0 + v**2.d0)
  Q = disscoeff*dabs(traceP/h**2.d0 - phi2)/(traceP/h**2.d0 + eps)*(modU**3.d0)

  D(1, ir, il) = -2.d0*alpha*p11*(modU**3.d0)
  D(2, ir, il) = -2.d0*alpha/H*p12*(modU**3.d0)
  D(3, ir, il) = -2.d0*alpha*p22*(modU**3.d0)

    TS_phys(2, ir, il) =  -g*R(ir)*Dtan(angle)*H &
    &-frottcoeff*R(ir)*dsqrt(u**2.d0 + v**2.d0)*u  
   IF ( (R(IR) - RMINUS) <= 2.D0*L1) THEN
     TS_phys(2, ir, il) =  -(4.D0*DEV/L1**4.D0*( (R(IR) - RMINUS - L1)**2.D0 &
      &- L1**2.D0)*(R(IR) - RMINUS - L1))*g*R(ir)*H-&
     & frottcoeff*R(ir)*dsqrt( u**2.d0 + v**2.d0 )*u  
   ENDIF
  TS_phys(3, ir, il) =  -frottcoeff*R(ir)*dsqrt(u**2.d0 + v**2.d0 )*v  
  TS_phys(4, ir, il) =  R(ir)*D(1, ir, il) 
  TS_phys(5, ir, il) =  R(ir)*D(2, ir, il)
  TS_phys(6, ir, il) =  R(ir)*D(3, ir, il) 
  TS_phys(7, ir, il) = -g*R(ir)*Dtan(angle)*H*U&
  & -frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**3.d0 - R(ir)*Q 
  IF ( (R(IR) - RMINUS) <= 2.D0*L1) THEN
   TS_phys(7, ir, il) = -(4.D0*DEV/L1**4.D0*( (R(IR) &
    &- RMINUS - L1)**2.D0 - L1**2.D0)*(R(IR) - RMINUS - L1))*g*R(ir)*H*U   &
   &-frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - R(ir)*Q  
  ENDIF
!------------------------------------------------------------------------------------------------------------------------
  DO k=1,7
    CONS1(k,ir, il) = CONS(k,ir, il) + DT*TS_phys(k, ir, il) 
  ENDDO 
ENDDO
ENDDO
CALL NOUVELL_VARIABLE_PRIM(Prim,CONS1, R, dr, it)


DO  ir = 1, Nr
  DO il  = 1, Ntheta  

      H = Prim(H_pv,ir, il)
      P11 = Prim(P11_pv,ir, il);
      p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)

      traceP = Prim(P11_pv,ir, il)+ Prim(P22_pv,ir, il) 

      u = Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
      fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);

      alpha = disscoeff*h**2.d0*dabs(traceP/H**2.d0 -&
      & phi2)/(traceP**2.d0 +  phi2*g*H_0*h**2.d0*eps)

      modU = dsqrt(u**2.d0 + v**2.d0)
      Q = disscoeff*dabs(traceP/h**2.d0 - phi2)/(traceP/h**2.d0 + eps)*(modU**3.d0)
    
      D(1, ir, il) = -2.d0*alpha*p11*(modU**3.d0)
      D(2, ir, il) = -2.d0*alpha/H*p12*(modU**3.d0)
      D(3, ir, il) = -2.d0*alpha*p22*(modU**3.d0)

      TS(2, ir, il) =  -g*R(ir)*Dtan(angle)*H -&
      &frottcoeff*R(ir)*dsqrt(u**2.d0 + v**2.d0)*u  
      IF ( (R(IR) - RMINUS) <= 2.D0*L1) THEN
       TS(2, ir, il) =  -(4.D0*DEV/L1**4.D0*( (R(IR) - &
        &RMINUS - L1)**2.D0 - L1**2.D0)*(R(IR) - RMINUS - L1))*g*R(ir)*H-&
       & frottcoeff*R(ir)*dsqrt( u**2.d0 + v**2.d0 )*u  
      ENDIF
      TS(3, ir, il) =  -frottcoeff*R(ir)*dsqrt(u**2.d0 + v**2.d0 )*v  
      TS(4, ir, il) =  R(ir)*D(1, ir, il) 
      TS(5, ir, il) =  R(ir)*D(2, ir, il)
      TS(6, ir, il) =  R(ir)*D(3, ir, il) 
      TS(7, ir, il) = -g*R(ir)*Dtan(angle)*H*U -&
      &frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**3.d0 - R(ir)*Q 
      IF ( (R(IR) - RMINUS) <= 2.D0*L1) THEN
       TS(7, ir, il) = -(4.D0*DEV/L1**4.D0*( (R(IR) - &
        &RMINUS - L1)**2.D0 - L1**2.D0)*(R(IR) - RMINUS - L1))*g*R(ir)*H*U   &
       &-frottcoeff*R(ir)*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - R(ir)*Q  
      ENDIF
!------------------------------------------------------------------------------------------------------------------------
      DO k=1,7
        CONS(k,ir, il) = CONS(k,ir, il) + 0.5d0*DT*(TS_phys(k,ir, il) + TS(k,ir, il))
      ENDDO 
     !erreur=2.d0*CONS(7,ir, il)-( CONS(2,ir, il)**2.d0+CONS(3,ir, il)**2.d0 )/CONS(1,ir, il)-&
    ! & g*CONS(1,ir, il)**2.d0/R(ir)-( CONS(4,ir, il) + CONS(6,ir, il) )
    ENDDO
  ENDDO
    CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
  
  return
  deallocate(TS_phys, D, ts, cons1)
  END SUBROUTINE
!***************************************************************************************************************************

SUBROUTINE ts_nonconserv_rk2(CONS, Prim, DT, R, DR, it)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ir, il, k, it
  REAL (KIND = DP) :: DT, dt2, dr, H, p11, p22, p12, u, v
  real(Kind=dp)    :: rhu, rh, rhv, rhe, rhp22
  REAL (KIND = DP) :: CONS(7,1:Nr, 1:Ntheta), R(0:Nr+1)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1)
  REAL (KIND = DP), ALLOCATABLE :: TS(:,:,:), TS1(:,:,:), cons1(:,:,:)
  allocate( TS(7, 1:Nr, 1:Ntheta), TS1(7, 1:Nr, 1:Ntheta), CONS1(7,1:Nr, 1:Ntheta)  )
  H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0
  TS = 0.d0; ts1 = 0.d0; CONS1 = 0.D0
  dt2 = 0.5d0*dt

  DO  ir = 1, Nr
  DO il  = 1, Ntheta  
    !TERME SOURCE
    H = Prim(H_pv,ir, il)
    P11 = Prim(P11_pv,ir, il); p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
    u =  Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
   
    TS1(2, ir, il) = h*(v*v + p22) + 0.5d0*g*h*h
    TS1(3, ir, il) =  - h*(u*v + p12) 
    TS1(5, ir, il) = v*(2.d0*p22 - p11)
    TS1(6, ir, il) = -2.d0*h*(p12*v+ p22*u)
  ENDDO
  ENDDO

  DO  ir = 1, Nr 
  DO il  = 1, Ntheta  
    DO k=1,7
      CONS1(k,ir, il) = CONS(k,ir, il) + DT*TS1(k,ir, il) 
    ENDDO 
    rhu=cons1(2,ir,il);rhv=cons1(3,ir,il);
    rh=cons1(1,ir,il);rhp22=cons1(6,ir,il); rhE=cons1(7,ir,il)
    cons1(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir))-rhp22-(rhu*rhu+rhv*rhv)/rh
  ENDDO
  ENDDO 
  
  CALL NOUVELL_VARIABLE_PRIM(Prim,CONS1, R, dr, it)

  DO ir = 1, Nr
  DO il = 1,Ntheta  
    !TERME SOURCE
    H = Prim(H_pv,ir, il)
    P11 = Prim(P11_pv,ir, il); p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
    u =  Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
   
    TS(2, ir, il) = h*(v*v + p22) + 0.5d0*g*h*h
    TS(3, ir, il) =  - h*(u*v + p12) 
    TS(5, ir, il) = v*(2.d0*p22 - p11)
    TS(6, ir, il) = -2.d0*h*(p12*v+ p22*u)
  ENDDO
  ENDDO

  DO  ir = 1, Nr
  DO il  = 1, Ntheta  
    DO k=1,7
      CONS(k,ir, il) = CONS(k,ir, il) + DT2*(TS1(k,ir, il) + TS(k,ir, il)) ! + DT*(TS(k,ir, il))  !
    ENDDO
     rhu=cons(2,ir,il);rhv=cons(3,ir,il);
     rh=cons(1,ir,il);rhp22=cons(6,ir,il); rhE=cons(7,ir,il)
     cons(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir) )-rhp22-(rhu*rhu+rhv*rhv)/rh
  ENDDO
  ENDDO 

  CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)

deallocate(ts, ts1, cons1)
return
END SUBROUTINE
!****************************************************************************************************************************************

Subroutine HLLC_theta_sub1(prim,flux, dt, R, L, DR, dtheta, IT, time, rang, reste)
USE precisions
use GlobalParam
IMPLICIT NONE
REAL (KIND = DP):: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), FLUX(7,0:Nr,0:Ntheta)
REAL(KIND=DP):: R(0:Nr+1), L(1:Ntheta)
real(kind=dp):: ul,ur,hl,hr,pr,pl,p11l,p11r
real(kind=dp)::p12r,p12l,cl,cr,ml,mr,sl,sr, rr, rl
real(kind=dp):: EL, ER,p22l,p22r,vr,vl,temp
real(kind=dp):: pstar,ustar,Estar,hstar,hp12star, p12star
real(kind=dp):: dt, dr, time, dtheta
INTEGER :: ir,il, it, rang, reste

DO il=0,ntheta
DO ir=1,nr

    !! etat gauche et droite
    ul=Prim(v_pv,ir,il);  ur=Prim(v_pv,ir, il+1)
    hl=Prim(h_pv,ir,il);  hr=Prim(h_pv,ir, il+1)
    vl=prim(u_pv,ir,il);  vr=Prim(u_pv,ir, il+1)

    p11l=prim(p22_pv,ir,il);  p11r=Prim(p22_pv,ir, il+1)
    p12l=prim(p12_pv,ir,il);  p12r=Prim(p12_pv,ir, il+1)
    p22l=prim(p11_pv,ir,il);  p22r=Prim(p11_pv,ir, il+1)

    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !! calcul des pression
    pl=g*hl*hl*0.5d0+hl*p11l
    pr=g*hr*hr*0.5d0+hr*p11r
    !! calcul des vitesses du son
    cl=dsqrt(g*hl+3.d0*p11l); cr=dsqrt(g*hr+3.d0*p11r)
    ! davis
    sr= (ur+cr); if( (ul+cl)> sr) sr = (ul+cl)
    !sr = dmax1(ul+cl, ur+cr)
    !if (dabs(sr) .le. 1.d-8) sr = 0.d0
    Sl=(ul-cl); if ( (ur-cr) < sl) sl = (ur - cr)
    !sr = - sl
   ! sl = dmin1(ul-cl, ur - cr)
   ! if (dabs(sl) .le. 1.d-8) sl = 0.d0
    ! etat star 
    ml=hl*(ul-sl)
    mr=hr*(ur-sr)

    ustar = (pl - pr + ml*ul - mr*ur)/(ml - mr)
    pstar= ( (ur-ul)*mr*ml-mr*pl+ml*pr )/(ml-mr)
!     if (dabs(ustar) .le. 1.d-8) then
!       ustar = 0.d0
!     endif

!     IF (cond_lim == 7) THEN
!       IF (ir == 0 .and. rang == 0) THEN
!         if (Prim(h_pv,1,il) .le. dev ) then
!           ustar = 0.d0
!         else
!          ustar = - dsqrt(Prim(U_pv,1, il) + 3.d0*Prim(p11_pv,1, il))
!        endif
!       ENDIF
!      ENDIF
!----------------------------------------------------------------------------------------------------
 !CALCUL LES FLUX

IF (ustar.ge.0.d0) THEN
   IF (sl.ge.0.d0) THEN
  !!etat gauche
   flux(1,ir,il)=hl*ul
   flux(2,ir,il)=hl*ul*ul+pl
   flux(3,ir,il)=hl*ul*vl
   flux(4,ir,il)=0.d0 !hl*p11l*ul !! equation inutile pour le cas r
   flux(5,ir,il)=p12l*ul
   flux(6,ir,il)=hl*p22l*ul
   flux(7,ir,il)=hl*El*ul+pl*ul

   ELSE !if (sl .lt. 0.d0 .or. bidule) then 

   !! etat star gauche
   hstar=ml/(ustar-sl)
   Estar=El+(pl*ul- pstar*ustar)/ml
   p12star = p12l*hstar/hl
  ! p22star = hl*p22l*(ul - sl)/(hstar*(ustar - sl))
   !! remplissage des flux
   flux(1,ir,il)=hstar*ustar
   flux(2,ir,il)=hstar*ustar*ustar+pstar
   flux(3,ir,il)=hstar*ustar*vl  
   flux(4,ir,il)=0.d0 !hstar*p11l*ustar ! equation inutile pour le cas r
   flux(5,ir,il)=p12star*ustar
   flux(6,ir,il)=hstar*p22l*ustar
   flux(7,ir,il)=hstar*Estar*ustar+pstar*ustar
  ENDIF
ELSE
   IF (sr.ge.0.d0) THEN
   
  !!etat droit etoile
   hstar=mr/(ustar-sr)
   Estar=Er+(pr*ur-pstar*ustar)/mr
   p12star = p12r*hstar/hr
   !p22star = hr*p22r*(ur - sr)/(hstar*(ustar - sr))
   !remplissage des flux
   flux(1,ir,il)=hstar*ustar
   flux(2,ir,il)=hstar*ustar*ustar+pstar   
   flux(3,ir,il)=hstar*ustar*vr
   flux(4,ir,il)=0.d0 !hstar*p11r*ustar ! equation inutile pour le cas r
   flux(5,ir,il)=p12star*ustar
   flux(6,ir,il)=hstar*p22r*ustar
   flux(7,ir,il)=hstar*Estar*ustar+pstar*ustar 
ELSE
  !!etat droit 
   flux(1,ir,il)=hr*ur
   flux(2,ir,il)=hr*ur*ur+pr
   flux(3,ir,il)=hr*ur*vr
   flux(4,ir,il)=0.d0 !hr*p11r*ur !! equation inutile pour le cas x
   flux(5,ir,il)=p12r*ur
   flux(6,ir,il)=hr*p22r*ur
   flux(7,ir,il)=hr*Er*ur+pr*ur
ENDIF
ENDIF

  !! inversion des flux :
  temp=flux(2,ir,il)
  flux(2,ir,il)=flux(3,ir,il)  
  flux(3,ir,il)=temp
  temp=flux(6,ir,il)
  flux(6,ir,il)=flux(4,ir,il)  
  flux(4,ir,il)=temp
ENDDO
ENDDO
return
end subroutine HLLC_theta_sub1
!************************************************************************************************************************************************************

Subroutine HLLC_theta_sub2(prim, flux, dt, R,L,DR,dtheta, IT, time, rang, reste)
USE precisions
USE GlobalParam
IMPLICIT NONE
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1), FLUX(7,0:Nr,0:Ntheta)
  REAL(KIND=DP) :: R(0:Nr+1), L(1:Ntheta)
  real(kind=dp) :: ul,ur,hl,hr,p11l,p11r,p12r
   real(kind=dp) ::p12l,ml,mr,sr,p22l,p22r, sl, rl, rr
  real(kind=dp) :: EL, ER,vr,vl
  real(kind=dp) :: vstar,Estar,hp12star
  real(kind=dp):: dt, dr, time, dtheta,temp
  INTEGER :: ir,il, it, rang, reste
  
do ir=1,nr
do il=0,ntheta

    ul=prim(v_pv,ir,il);     ur=prim(v_pv,ir, il+1)
    hl=prim(h_pv,ir,il);     hr=prim(h_pv,ir, il+1)
    vl=prim(u_pv,ir,il);     vr=prim(u_pv,ir, il+1)
    p11l=prim(p22_pv,ir,il); p11r=prim(p22_pv,ir, il+1)
    p12l=prim(p12_pv,ir,il); p12r=prim(p12_pv,ir, il+1)
    p22l=prim(p11_pv,ir,il); p22r=prim(p11_pv,ir, il+1)
    
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0

   !sr=dmax1(dsqrt(p11l),dsqrt(p11r)) !
  ! sl=dmin1(-dsqrt(p11l),-dsqrt(p11r))!

   sr=dsqrt(p11r)
     sl = -dsqrt(p11l)

  vstar         = (hl*(p12l- sl*vl) - hr*(p12r-sr*vr))/(sr*hr - hl*sl)
  hp12star      = (hr*hl)/(hr*sr-hl*sl)*(sr*p12l-sl*p12r+sl*sr*(vr-vl))

  Flux(1,ir,il) = 0.d0
  Flux(2,ir,il) = 0.d0
  flux(3,ir,il) = hp12star
  Flux(4,ir,il) = 0.d0
  Flux(5,ir,il) = vstar !! non conservatif
  flux(6,ir,il) = 2.d0*hp12star*vstar !!inutile
  flux(7,ir,il) = hp12star*vstar

  !! inversion des flux :
  temp=flux(2,ir,il)
  flux(2,ir,il)=flux(3,ir,il)  
  flux(3,ir,il)=temp
  temp=flux(6,ir,il)
  flux(6,ir,il)=flux(4,ir,il)  
  flux(4,ir,il)=temp  
ENDDO
ENDDO
return
end subroutine HLLC_theta_sub2
!****************************************************************************************************

subroutine godunov_r_sub1(cons, flux, dt, R, DR, Prim)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  real(KIND=dp) :: cons(7,1:Nr,1:Ntheta), FLUX(7,0:Nr,0:Ntheta)
  real(KIND=dp) :: R(0:Nr+1), Prim(Nv_Prim,0:Nr+1,0:Ntheta+1)
  real(Kind=dp) :: rhu, rh, rhv, rhe, rhp22, dt, DR
  real(Kind=dp) :: h, u, v, p11, p12, p22
  INTEGER :: k, ir, il

  do ir=1,nr 
  do il=1,ntheta
     h = Prim(h_pv,ir, il); u = Prim(u_pv,ir, il); 
     v = Prim(v_pv,ir, il)
    p11 = Prim(p11_pv,ir, il);
     p12 = Prim(p12_pv,ir, il); p22 = Prim(p22_pv,ir, il)

    cons(1,ir,il)=cons(1,ir,il)+dt/DR*( ( r(ir) &
      &- 0.5d0*dr)*(Flux(1,ir-1,il) )- ( r(ir) + 0.5d0*dr)*( flux(1,ir,il) ) ) 

!****************************************************************************************************
!     cons(2,ir,il)=cons(2,ir,il)+dt/DR*( r(ir-1)*( Flux(2,ir-1,il))- r(ir)*(flux(2,ir,il)) )

!     cons(3,ir,il) = cons(3,ir,il)+dt/DR*( r(ir-1)*(Flux(3,ir-1,il) )- r(ir)*( flux(3,ir,il) ) )

!     cons(5,ir,il) = cons(5,ir,il)+dt/DR*( r(ir-1)*(Flux(5,ir-1,il) )- r(ir)*( flux(5,ir,il) ) )

!     cons(6,ir,il) = cons(6,ir,il)+dt/DR*( r(ir-1)*(Flux(6,ir-1,il) )- r(ir)*( flux(6,ir,il) ) )

!     cons(7,ir,il) = cons(7,ir,il)+dt/DR*( r(ir-1)*(Flux(7,ir-1,il) )- r(ir)*( flux(7,ir,il)) )
 !*****************************************************************************************************

    cons(2,ir,il)=cons(2,ir,il)+dt/DR*( ( r(ir) - 0.5d0*dr)*( Flux(2,ir-1,il) - &
    &      0.5d0*g*h*h - h*v*v-h*p22)-&
    & ( r(ir) + 0.5d0*dr)*(flux(2,ir,il)- 0.5d0*g*h*h - h*v*v -h*p22) )

    cons(3,ir,il) = cons(3,ir,il)+dt/DR*( ( r(ir) - 0.5d0*dr)*(Flux(3,ir-1,il) + &
    &       h*u*v + h*p12)- ( r(ir) + 0.5d0*dr)*( flux(3,ir,il) + h*u*v + h*p12) )

    cons(5,ir,il) = cons(5,ir,il)+dt/DR*( ( r(ir)- 0.5d0*dr)*(Flux(5,ir-1,il) - &
      & 2.d0*v*p22 +v*p11)- ( r(ir) + 0.5d0*dr)*( flux(5,ir,il) -2.d0*v*p22 +v*p11) )

    cons(6,ir,il) = cons(6,ir,il)+dt/DR*( ( r(ir) -0.5d0*dr)*(Flux(6,ir-1,il) + &
      & 2.d0*h*(p12*v + p22*u) )- ( r(ir) + 0.5d0*dr)*( flux(6,ir,il) +2.d0*h*(p12*v + p22*u)) )

    cons(7,ir,il) = cons(7,ir,il)+dt/DR*( ( r(ir) &
      &- 0.5d0*dr)*(Flux(7,ir-1,il) )- ( r(ir) + 0.5d0*dr)*( flux(7,ir,il)) )

 !! corection de la variable p11
 rhu=cons(2,ir,il);rhv=cons(3,ir,il);rh=cons(1,ir,il);rhp22=cons(6,ir,il); rhE=cons(7,ir,il)
 cons(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir))-rhp22-(rhu*rhu+rhv*rhv)/rh
end do
end do
end subroutine godunov_r_sub1

 !************************************************************************************************************************************************************
 
 subroutine godunov_theta_sub1(cons,flux,dt,Dtheta, R, dr)
  USE precisions
  USE GlobalParam
  IMPLICIT NONE
  real(KIND=dp) :: cons(7,1:Nr,1:Ntheta), flux(7,0:Nr,0:Ntheta), R(0:Nr+1)
  real(Kind=dp) :: rhu,rh,rhv,rhe,rhp11,dt,Dtheta, dr
  INTEGER :: k,ir,il
  
DO ir=1,nr 
DO il=1,Ntheta
  DO k=1,7
   cons(k,ir,il)=cons(k,ir,il)+dt/dtheta*(Flux(k,ir,il-1)-flux(k,ir,il))
  ENDDO
   !! corection de la variable p22
   rhu=cons(2,ir,il); rhv=cons(3,ir,il); rh=cons(1,ir,il); rhp11=cons(4,ir,il); rhE=cons(7,ir,il)
   cons(6,ir,il)=2.d0*rhE-g*rh*rh/r(ir)-rhp11-(rhu*rhu+rhv*rhv)/rh
ENDDO
ENDDO
end subroutine godunov_theta_sub1
!************************************************************************************************************************************************************

subroutine godunov_r_sub2(cons,flux,dt, R, DR)
USE precisions
use GlobalParam
implicit none
real(KIND=dp) :: cons(7,1:Nr,1:Ntheta), FLUX(7,0:Nr,0:Ntheta), R(0:Nr+1)
real(Kind=dp) :: rhu,rh,rhv,rhe,p11,dt,DR, rhp11
INTEGER :: k,ir,il

 do ir=1,nr
 do il=1,ntheta
  p11=cons(4,ir,il)/cons(1,ir,il)
  !------------------------------------------------------------------------------------------------------------------------
  do k=1,4
   cons(k,ir,il)=cons(k,ir,il)+dt/DR*( ( r(ir) -&
    &0.5d0*dr)*Flux(k,ir-1,il)-( r(ir) + 0.5d0*dr)*flux(k,ir,il) )
  end do
  !------------------------------------------------------------------------------------------------------------------------
    cons(5,ir,il)=cons(5,ir,il)+dt/DR*( ( r(ir) - &
      &0.5d0*dr)*Flux(5,ir-1,il)-( r(ir) + 0.5d0*dr)*flux(5,ir,il) )*p11 !! ici dans le flux on a stocké vstar
  !------------------------------------------------------------------------------------------------------------------------
   do k=6,7
    cons(k,ir,il)=cons(k,ir,il)+dt/DR*(( r(ir) &
      &- 0.5d0*dr)*Flux(k,ir-1,il)-( r(ir) + 0.5d0*dr)*flux(k,ir,il))
   end do
  !------------------------------------------------------------------------------------------------------------------------
  !! corection de la variable p22
   rhu=cons(2,ir,il); rhv=cons(3,ir,il); 
   rh=cons(1,ir,il); rhp11=cons(4,ir,il); rhE=cons(7,ir,il)
   cons(6,ir,il)=2.d0*rhE-(rhu*rhu+rhv*rhv)/rh-g*rh*rh/(R(ir)) -rhp11
end do
end do
end subroutine godunov_r_sub2
 !************************************************************************************************************************************************************

subroutine godunov_theta_sub2(cons,flux,dt,Dtheta, R, dr)
USE precisions
USE GlobalParam
implicit none
real(KIND=dp) :: cons(7,1:Nr,1:Ntheta), flux(7,0:Nr,0:Ntheta), R(0:Nr+1)
real(Kind=dp) :: rhu, rh, rhv, rhe, p22, dt, Dtheta, dr, dl, rhp22
INTEGER :: k, ir, il

DO ir=1,nr
DO il=1,Ntheta
  p22=cons(6,ir,il)/cons(1,ir,il)
   DO k=1,4
    cons(k,ir,il)=cons(k,ir,il) + dt/dtheta*(Flux(k,ir,il-1)-flux(k,ir,il))
   ENDDO
    cons(5,ir,il)=cons(5,ir,il) + dt/dtheta*(Flux(5,ir,il-1)-flux(5,ir,il))*p22
   DO k=6,7
    cons(k,ir,il)=cons(k,ir,il) + dt/dtheta*(Flux(k,ir,il-1)-flux(k,ir,il))
   ENDDO
  !! corection de la variable p11
  rhu=cons(2,ir,il); rhv=cons(3,ir,il); rh=cons(1,ir,il); rhp22=cons(6,ir,il); rhE=cons(7,ir,il)
  cons(4,ir,il)=2.d0*rhE-(rhu*rhu+rhv*rhv)/rh-g*rh*rh/r(ir)-rhp22
ENDDO
ENDDO
END subroutine godunov_theta_sub2
!****************************************************************************************************************************

SUBROUTINE ts_nonconserv_rk4( CONS, Prim, DT, R, DR, it )
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ir, il, k, it
  REAL (KIND = DP) :: DT, dt2, dr, H, p11, p22, p12, u, v
  real(Kind=dp)    :: rhu, rh, rhv, rhe, rhp22
  REAL (KIND = DP) :: CONS(7,1:Nr, 1:Ntheta), R(0:Nr+1)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1)
  REAL (KIND = DP), ALLOCATABLE :: TS1(:,:,:), TS2(:,:,:), TS3(:,:,:), TS4(:,:,:)
  REAL (KIND = DP), ALLOCATABLE ::cons1(:,:,:), cons2(:,:,:), cons3(:,:,:)
  allocate( TS1(7, 1:Nr, 1:Ntheta), TS2(7, 1:Nr, 1:Ntheta), &
    &TS3(7, 1:Nr, 1:Ntheta),TS4(7, 1:Nr, 1:Ntheta) )
  allocate(CONS1(7,1:Nr, 1:Ntheta),&
    &CONS2(7,1:Nr, 1:Ntheta),CONS3(7,1:Nr, 1:Ntheta) )
  H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0
  TS1(:,:,:) = 0.d0; ts2(:,:,:) = 0.d0; 
   ts3(:,:,:) = 0.d0; ts4(:,:,:) = 0.d0;
  CONS1(:,:,:) = 0.D0; 
  CONS2(:,:,:) = 0.D0; CONS3(:,:,:) = 0.D0
  dt2 = 0.5d0*dt

  DO  ir = 1, Nr
  DO il  = 1, Ntheta  
    !TERME SOURCE
    H = Prim(H_pv,ir, il)
    P11 = Prim(P11_pv,ir, il); p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
    u =  Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
   
    TS1(2, ir, il) = h*(v*v + p22) + 0.5d0*g*h*h
    TS1(3, ir, il) =  - h*(u*v + p12) 
    TS1(5, ir, il) = h*v*(2.d0*p22 - p11)
    TS1(6, ir, il) = -2.d0*h*(p12*v+ p22*u)
  ENDDO
  ENDDO

  DO  ir = 1, Nr
  DO il  = 1, Ntheta  
    DO k=1,7
      CONS1(k,ir, il) = CONS(k,ir, il) + DT*TS1(k,ir, il) 
    ENDDO 
    rhu=cons1(2,ir,il);rhv=cons1(3,ir,il);rh=cons1(1,ir,il);
    rhp22=cons1(6,ir,il); rhE=cons1(7,ir,il)
    cons1(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir))-rhp22-(rhu*rhu+rhv*rhv)/rh
  ENDDO
  ENDDO 
  
  CALL NOUVELL_VARIABLE_PRIM(Prim,CONS1, R, dr, it)
!-----------------------------------------------------------------------------------------------------------
  DO ir = 1, Nr
  DO il = 1,Ntheta  
    !TERME SOURCE
    H = Prim(H_pv,ir, il)
    P11 = Prim(P11_pv,ir, il); 
    p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
    u =  Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
   
    TS2(2, ir, il) = h*(v*v + p22) + 0.5d0*g*h*h
    TS2(3, ir, il) =  - h*(u*v + p12) 
    TS2(5, ir, il) = h*v*(2.d0*p22 - p11)
    TS2(6, ir, il) = -2.d0*h*(p12*v+ p22*u)
  ENDDO
  ENDDO

  DO  ir = 1, Nr  
  DO il  = 1, Ntheta  
    DO k=1,7
      CONS2(k,ir, il) = CONS(k,ir, il) + DT2*TS2(k,ir, il) 
    ENDDO
     rhu=cons2(2,ir,il);rhv=cons2(3,ir,il);
     rh=cons2(1,ir,il);rhp22=cons2(6,ir,il); rhE=cons2(7,ir,il)
     cons2(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir))-rhp22-(rhu*rhu+rhv*rhv)/rh
  ENDDO
  ENDDO 

  CALL NOUVELL_VARIABLE_PRIM(Prim,CONS2, R, dr, it)
  !-----------------------------------------------------------------------------------------------------------

  DO ir = 1, Nr
  DO il = 1,Ntheta  
    !TERME SOURCE
    H = Prim(H_pv,ir, il)
    P11 = Prim(P11_pv,ir, il);
     p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
    u =  Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
   
    TS3(2, ir, il) = h*(v*v + p22) + 0.5d0*g*h*h
    TS3(3, ir, il) =  - h*(u*v + p12) 
    TS3(5, ir, il) = h*v*(2.d0*p22 - p11)
    TS3(6, ir, il) = -2.d0*h*(p12*v+ p22*u)
  ENDDO
  ENDDO

  DO  ir = 1, Nr 
  DO il  = 1, Ntheta  
    DO k=1,7
      CONS3(k,ir, il) = CONS(k,ir, il) + DT2*TS3(k,ir, il) 
    ENDDO
     rhu=cons3(2,ir,il);rhv=cons3(3,ir,il);
     rh=cons3(1,ir,il);rhp22=cons3(6,ir,il); rhE=cons3(7,ir,il)
     cons3(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir) )-rhp22-(rhu*rhu+rhv*rhv)/rh
  ENDDO
  ENDDO 

  CALL NOUVELL_VARIABLE_PRIM(Prim,CONS3, R, dr, it)

!-----------------------------------------------------------------------------------------------------------

  DO ir = 1, Nr
  DO il = 1,Ntheta  
    !TERME SOURCE
    H = Prim(H_pv,ir, il)
    P11 = Prim(P11_pv,ir, il);
     p22 = Prim(P22_pv,ir, il); p12 = Prim(P12_pv,ir, il)
    u =  Prim(U_pv,ir, il); v = Prim(V_pv,ir, il)
   
    TS4(2, ir, il) = h*(v*v + p22) + 0.5d0*g*h*h
    TS4(3, ir, il) =  - h*(u*v + p12) 
    TS4(5, ir, il) = h*v*(2.d0*p22 - p11)
    TS4(6, ir, il) = -2.d0*h*(p12*v+ p22*u)
  ENDDO
  ENDDO

  DO  ir = 1, Nr

  DO il  = 1, Ntheta  
    DO k=1,7
      CONS(k,ir, il) = CONS(k,ir, il) + DT*(TS1(k,ir, il) + &
        &2.d0*TS2(k,ir, il)+ 2.d0*TS3(k,ir, il)+ TS4(k,ir, il))/6.d0
    ENDDO
     rhu=cons(2,ir,il);rhv=cons(3,ir,il);
     rh=cons(1,ir,il);rhp22=cons(6,ir,il); rhE=cons(7,ir,il)
     cons(4,ir,il)=2.d0*rhE-g*rh*rh/(R(ir) )-rhp22-(rhu*rhu+rhv*rhv)/rh
  ENDDO
  ENDDO 
CALL NOUVELL_VARIABLE_PRIM(Prim,CONS, R, dr, it)
deallocate(ts1, ts2,ts3,ts4, cons1, cons2,cons3)
return
END SUBROUTINE


SUBROUTINE COMMUNICATION(Prim, ncpu, rang, R)
!Subroutine permettant l echange des donnees communes entre processeurs
  USE precisions
  use GlobalParam
 IMPLICIT NONE
!Variables entree
 INTEGER,INTENT(IN) :: ncpu,rang
 !Variables entree-sortie
 real(kind=dp) :: Prim(Nv_prim, 0:Nr+1,0:Ntheta+1 ), R(0:Nr+1)
 !Variables locales
 INTEGER :: dest, source, tag, code
 INTEGER :: status(MPI_STATUS_SIZE)
 REAL(KIND=DP) :: h, u, v, p11, p12, p22, tampon_env_r, tampon_rec_r, r1
 REAL(KIND=DP),ALLOCATABLE ::  tampon_env_h(:), tampon_env_u(:), &
 &tampon_env_v(:), tampon_env_p11(:), tampon_env_p12(:), tampon_env_p22(:)  
 REAL(KIND=DP),ALLOCATABLE ::  tampon_rec_h(:), tampon_rec_u(:),&
 & tampon_rec_v(:), tampon_rec_p11(:), tampon_rec_p12(:), tampon_rec_p22(:)
 INTEGER :: i, taille_max, nbpaq, nbel, el_deb, el_fin

ALLOCATE(tampon_env_h(1:Ntheta), tampon_env_u(1:Ntheta), &
  &tampon_env_v(1:Ntheta), tampon_env_p11(1:Ntheta),&
&  tampon_env_p12(1:Ntheta), tampon_env_p22(1:Ntheta))  

ALLOCATE(tampon_rec_h(1:Ntheta), tampon_rec_u(1:Ntheta),&
& tampon_rec_v(1:Ntheta), tampon_rec_p11(1:Ntheta),&
& tampon_rec_p12(1:Ntheta), tampon_rec_p22(1:Ntheta)) 
tag = 0
tampon_env_h= 0.d0; tampon_env_u= 0.d0; tampon_env_v= 0.d0;
 tampon_env_p11= 0.d0; tampon_env_p12= 0.d0; tampon_env_p22= 0.d0 
tampon_rec_h= 0.d0; tampon_rec_u= 0.d0; tampon_rec_v= 0.d0; 
tampon_rec_p11= 0.d0; tampon_rec_p12= 0.d0; tampon_rec_p22= 0.d0 
tampon_env_r = 0.d0;  tampon_rec_r= 0.d0
!Communication avec le proc de rang inferieur
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

IF(rang.NE.0) THEN

tampon_env_h(1:Ntheta)=Prim(h_pv,1, 1:Ntheta)  
tampon_env_u(1:Ntheta)=Prim(u_pv,1, 1:Ntheta)  
tampon_env_v(1:Ntheta)=Prim(v_pv,1, 1:Ntheta)  
tampon_env_p11(1:Ntheta)= Prim(p11_pv,1, 1:Ntheta) 
tampon_env_p12(1:Ntheta)= Prim(p12_pv,1, 1:Ntheta) 
tampon_env_p22(1:Ntheta)= Prim(p22_pv,1, 1:Ntheta)  
tampon_env_r= R(1)

!Envoi vers proc de rang inferieur
dest   = rang-1

CALL MPI_SEND(tampon_env_h,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_u,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_v,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_p11,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_p12,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_p22,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_r,numb,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)

!Reception du proc de rang inferieur
source = rang-1

CALL MPI_RECV(tampon_rec_h,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_u,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_v,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_p11,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_p12,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_p22,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_r,numb,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)

Prim(h_pv, 0, 1:Ntheta)=tampon_rec_h(1:Ntheta)
Prim(u_pv, 0, 1:Ntheta)=tampon_rec_u(1:Ntheta)
Prim(v_pv, 0, 1:Ntheta)=tampon_rec_v(1:Ntheta)
Prim(p11_pv, 0, 1:Ntheta)=tampon_rec_p11(1:Ntheta)
Prim(p12_pv, 0, 1:Ntheta)=tampon_rec_p12(1:Ntheta)
Prim(p22_pv, 0, 1:Ntheta)=tampon_rec_p22(1:Ntheta)
r(0) = tampon_rec_r
do il = 1, ntheta
  h = Prim(h_pv, 0, il); p11 = Prim(p11_pv, 0, il); p22 = Prim(p22_pv, 0, il)

  Prim(sound_ar_pv,0,il)= Sound_ar(h, p11)  
  Prim(sound_br_pv,0,il)= Sound_br( p11) 
  Prim(sound_ath_pv,0,il)= Sound_atheta(h, p22) 
  Prim(sound_bth_pv,0,il)=  Sound_btheta(p22) 
  Prim(press_pv,0,il) = Pression(h, p11) 
  Prim(ein_pv,0,il) = InternalEn(h, p11, p22)
enddo   

ENDIF

!Communication avec le proc de rang superieur
!--------------------------------------------

IF (rang.NE.ncpu-1) THEN

tampon_env_h(1:Ntheta)=Prim(h_pv,Nr, 1:Ntheta)  
tampon_env_u(1:Ntheta)=Prim(u_pv,Nr, 1:Ntheta)  
tampon_env_v(1:Ntheta)=Prim(v_pv,Nr, 1:Ntheta)  
tampon_env_p11(1:Ntheta)= Prim(p11_pv,Nr, 1:Ntheta) 
tampon_env_p12(1:Ntheta)= Prim(p12_pv,Nr, 1:Ntheta) 
tampon_env_p22(1:Ntheta)= Prim(p22_pv,Nr, 1:Ntheta)     
tampon_env_r= r(nr) 

!Envoi vers proc de rang superieur
dest   = rang+1

CALL MPI_SEND(tampon_env_h,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_u,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_v,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_p11,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_p12,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_p22,Ntheta,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
CALL MPI_SEND(tampon_env_r,numb,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)

!Reception du proc de rang superieur
source = rang+1

CALL MPI_RECV(tampon_rec_h,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_u,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_v,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_p11,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_p12,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_p22,Ntheta,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
CALL MPI_RECV(tampon_rec_r,numb,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)

!DEALLOCATE(tampon)

Prim(h_pv,nr+1, 1:Ntheta)=tampon_rec_h(1:Ntheta)
Prim(u_pv,nr+1, 1:Ntheta)=tampon_rec_u(1:Ntheta)
Prim(v_pv,nr+1, 1:Ntheta)=tampon_rec_v(1:Ntheta)
Prim(p11_pv,nr+1, 1:Ntheta)=tampon_rec_p11(1:Ntheta)
Prim(p12_pv,nr+1, 1:Ntheta)=tampon_rec_p12(1:Ntheta)
Prim(p22_pv,nr+1, 1:Ntheta)=tampon_rec_p22(1:Ntheta)
r(nr+1) = tampon_rec_r
do il = 1, ntheta
  h = Prim(h_pv, nr+1, il); p11 = Prim(p11_pv, nr+1, il); p22 = Prim(p22_pv, nr+1, il)

  Prim(sound_ar_pv,nr+1,il)=Sound_ar(h, p11)  
  Prim(sound_br_pv,nr+1,il)=Sound_br( p11) 
  Prim(sound_ath_pv,nr+1,il)=Sound_atheta(h, p22) 
  Prim(sound_bth_pv,nr+1,il)=Sound_btheta(p22) 
  Prim(press_pv,nr+1,il)=Pression(h, p11) 
  Prim(ein_pv,nr+1,il)=InternalEn(h, p11, p22)
    enddo   
ENDIF

DEALLOCATE(tampon_env_h, tampon_env_u, tampon_env_v, tampon_env_p11, tampon_env_p12,tampon_env_p22)
DEALLOCATE(tampon_rec_h, tampon_rec_u, tampon_rec_v, tampon_rec_p11, tampon_rec_p12,tampon_rec_p22)
RETURN
END SUBROUTINE COMMUNICATION
!-------------------------------------------------------------------------------------------------------------------------

subroutine newton(R,L, ncpu, rang, dr)
  USE precisions
  use GlobalParam
  INTEGER :: ir, il, k, iter, iter_max = 1000000
  INTEGER,INTENT(IN) :: ncpu,rang
  REAL (KIND = DP) :: DT, error, err, R(0:Nr+1), L(1:Ntheta), dr
  REAL (KIND = DP) :: fh(1:Nr,1:Ntheta), fprimh(1:Nr,1:Ntheta)
  REAL (KIND = DP) :: h1(1:Nr,1:Ntheta),h2(1:Nr,1:Ntheta),hprim(1:Nr,1:Ntheta), h_new(1:Nr,1:Ntheta)
  INTEGER :: dest, source, tag, code
  INTEGER :: status(MPI_STATUS_SIZE)
  CHARACTER(LEN=3) :: NB, Zero="000"
  CHARACTER*3 :: numrg

  err = 1.d0
  iter = 0; h1 = 0.d0; h2 = 0.d0; 
  h_new(:,:) = h_0; h1(:,:) = h_0;fh = 0.d0; fprimh = 0.d0; hprim(:,:) = 0.d0
if (ntheta .le. 1) then
  L(:) = 0.d0
endif
WRITE(numrg,'(i3.3)') rang

do while (iter .le. iter_max .and. err .gt. 1.d-16)

do ir = 1, Nr; do il = 1, ntheta
fh(ir,il) = q0_jump**2.d0/(2.d0*R(ir)**2.d0) + g*h_new(ir,il)**3.d0 -&
&      (q0_jump**2.d0/(2.d0*(rplus)**2.d0*h_0**2.d0) + g*h_0)*h_new(ir,il)**2.d0

!hprim = 0.5d0*(h_new(ir+1,il) - h_new(ir,il) )

fprimh(ir,il) = 3.d0*g*h_new(ir,il)**2.d0 - &
&                    2.d0*(q0_jump**2.d0/(2.d0*(rplus)**2.d0*h_0**2.d0) + g*h_0)*h_new(ir,il)
h1(ir,il)=h_new(ir,il)

if ( dabs(fh(ir,il)/fprimh(ir,il)) .gt. 1.d-9 ) then 
h_new(ir,il) = h_new(ir,il)  - fh(ir,il)/fprimh(ir,il)
else
h_new(ir,il) = h1(ir,il) 
endif
h2(ir,il) = dabs( h_new(ir,il) - h1(ir,il) )  

if (h_new(ir,il).le. 1.d-10 .and. cond_lim == 8) then 
print*, 'fh(ir,il)', fh(ir,il), 'fprimh(ir,il)', fprimh(ir,il)
print*,'h_new(ir,il)', h_new(ir,il), 'fh(ir,il)/fprimh(ir,il)', fh(ir,il)/fprimh(ir,il)
print*, 'h is negative'
!pause
endif

 !print*, 'fh(ir,il)', fh(ir,il), 'fprimh(ir,il)', fprimh(ir,il), 'h_new(ir,il)', h_new(ir,il), 'fh(ir,il)/fprimh(ir,il)', fh(ir,il)/fprimh(ir,il)
  enddo; enddo
  
  !call MPI_ALLREDUCE(err, error, NUMB, MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,code)
  !print*, 'error =', error, rang
  err = sum(h2(1:Nr,1:Ntheta))/(Nr*Ntheta)
  iter = iter + 1
enddo

OPEN(MyUnit+rang+1,FILE = './resu/newton_'//numrg//'.out')
do ir = 1, Nr; do il = 1, ntheta
WRITE(MyUnit+rang+1,'(3(E20.13,1X))') R(ir)*DCOS(L(il) ), R(ir)*DSIN( L(il) ), h_new(ir, il)
enddo; enddo
WRITE(MyUnit+rang+1,*)   

END SUBROUTINE newton




!************************************************************************************************************************************************


SUBROUTINE Coriolis( CONS, Prim,R, dr, DT, it)
  USE precisions
  USE GlobalParam
  implicit none 
  INTEGER :: ir, il, IT,k
  REAL (KIND = DP) :: DT, H, u, v,dr, f
  REAL (KIND = DP) :: CONS(7,1:Nr, 1:Ntheta), R(0:Nr+1)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nr+1,0:Ntheta+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:) :: TS_phys, D
  allocate(TS_phys(7), D(3))
  TS_phys = 0.d0; D = 0.d0
  H= 0.d0;  u= 0.d0; v= 0.d0
  f = 0.01d0
!------------------------------------------------------------------------------------------------------------------------
  DO  ir = 1, Nr
  DO  il  = 1, Ntheta  
!------------------------------------------------------------------------------------------------------------------------ 

       !TERME SOURCE
      H = Prim(H_pv,ir, il)
      u = Prim(u_pv,ir, il)
      v = Prim(v_pv,ir, il)

      TS_phys(2) = r(ir)*h*f*v + f**2.d0*r(ir)**2.d0*h/4.d0
      TS_phys(3) = -r(ir)*f*u*h
      TS_phys(7) = f**2.d0*r(ir)**2.d0*h*u/4.d0

!------------------------------------------------------------------------------------------------------------------------
      DO k=1,7
        CONS(k,ir, il) = CONS(k,ir, il) + DT*TS_phys(k) 
      ENDDO 
  ENDDO
  ENDDO
  return

  deallocate(TS_phys, D)
  END SUBROUTINE Coriolis

END PROGRAM code2D


   
