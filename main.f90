
module precisions
implicit none
intEinGer, parameter :: dp=kind(1.0d0)
end module precisions
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
      
PROGRAM code2D
USE precisions
  IMPLICIT NONE  
  INTEGER                                         ::  ImpreE, ImpreF, Nv_Prim, argunit = 6
  INTEGER                                         ::  Nx, Ny, I, IT,iterfinal, cond_lim, ix, iy
  INTEGER                                         ::  iv, H_pv, U_pv, V_pv, Ein_pv, P11_pv,  P12_pv, P22_pv, Sound_pv, Pres_pv

  REAL (KIND = DP)                                :: DX, Dy, Dh, DT, Lx, Ly, X0, Y0, T1_CPU, T2_CPU, H_0, amplitude
  REAL (KIND = DP)                                :: CFL, UMAX, TIME, TIMEOUT, period_time, TIME2
  REAL (KIND = DP)                                :: FLUXi(6), Gi(6)
  REAL (KIND = DP), PARAMETER                     :: EPS = 1.d-8
  REAL (KIND = DP), ALLOCATABLE                   :: Prim(:,:,:), Ein(:,:), Sound(:,:), Pression(:,:), X(:), Y(:), PrimG(:), PrimD(:)
  REAL (KIND = DP), ALLOCATABLE                   :: CONS(:,:,:), FLUX(:,:,:), GFLUX(:,:,:)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:)     :: MaxVP, MinVp
  REAL (KIND = DP)                                :: Hstar, Pstar, Ustar, Estar, Cmax
  REAL(KIND=DP)                                   :: xmax, XMIN, XMAX1, YMIN, YMAX 

  CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:)    :: NamesOfPrimVar

!----------------------------------------------------------------------------------------------------
pi = 4.0d0*ATAN(1.0d0)
CALL CPU_TIME(T1_CPU)
CALL LECTURE_DONNEES( angle, Nx,Ny, Lx, Ly, TIMEOUT, iterfinal, g, CFL, H_0, amplitude,  ImpreE, ImpreF)
!----------------------------------------------------------------------------------------------------
  Nv_Prim  = 9
  H_pv     = 1
  U_pv     = 2
  V_pv     = 3
  P11_pv   = 4
  P12_pv   = 5
  P22_pv   = 6
  Ein_pv   = 7
  Sound_pv = 8
  Pres_pv  = 9
  !---------------------------------------------------------------------------------------

DX = Lx/DFLOAT(Nx)
Dy = Ly/DFLOAT(Ny)
Dh =  MIN(Dx,Dy)
X0 = 0.D0
Y0 = 0.d0

 print*, 'Nx =', Nx,'Ny =', Ny, 'DX =', DX, 'DY =', DY
 !----------------------------------------------------------------------------------------------------
 ALLOCATE( Ein(0:Nx+1,0:Ny+1), Sound(0:Nx+1,0:Ny+1), X(1:Nx), Y(1:Ny), Pression(0:Nx+1,0:Ny+1) )  
 ALLOCATE( Prim(6,0:Nx+1,0:Ny+1), CONS(6,1:Nx,1:Ny), FLUX(6,0:Nx+1,0:Ny+1), GFLUX(6,0:Nx+1,0:Ny+1) )  
 ALLOCATE( NamesOfPrimVar(Nv_Prim) )
 ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim) )
!----------------------------------------------------------------------------------------------------
  NamesOfPrimVar(H_pv)         = "Depht Lenght"
  NamesOfPrimVar(U_pv)         = "Velocity (x)"
  NamesOfPrimVar(V_pv)         = "Velocity (y)"
  NamesOfPrimVar(V_pv)         = "Velocity (y)"
  NamesOfPrimVar(P11_pv)       = "component of tensor P11"
  NamesOfPrimVar(P12_pv)       = "component of tensor P12"
  NamesOfPrimVar(P22_pv)       = "component of tensor P22"
  NamesOfPrimVar(Sound_pv)     = "Sound Speed"
  NamesOfPrimVar(Pres_pv)      = "Pressure"
  !----------------------------------------------------------------------------------------------------

  DO iv = 1, Nv_Prim
     WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ", iv, Nx, Ny
  END DO
  WRITE(6,*) " >>> End of LECTURE_DONNEES"
  !----------------------------------------------------------------------------------------------------
  !INITIALISATION
  CALL INITIALISATION(X0, Y0, DX, DY, X, Y, Prim, SOUND, Ein, Nx, Ny,   CONS, angle, g, Pression, H_0, cond_lim, Lx, Ly, amplitude, Nv_Prim)

 ouvert = .true.
 !----------------------------------------------------------------------------------------------------

 TIME = 0.D0
 period_time = 1.0d0;
 TIME2  = period_time ;
 IT = 1

 CALL PutonScreen()

 !----------------------------------------------------------------------------------------------------
 XMIN = MINVAL(X(:))
 XMAX1 = MaxVal(X(:))

 YMIN = MINVAL(Y(:))
 YMAX = MaxVal(Y(:))
  
 WRITE(6,'( 2(A10,E15.6))')  " Xmin = ", XMIN, &
       &                     " Xmax = ", XMAX1
  WRITE(6,'( 2(A10,E15.6))') " Ymin = ", YMIN, &
       &                     " Ymax = ", YMAX
!----------------------------------------------------------------------------------------------------

 CALL Ecriture_donnees(X,Y, Prim, Ein, Pression, Nx,Ny, time, DX, DY )
!----------------------------------------------------------------------------------------------------
!-------------------------do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT)!-------------------
 ! BOUCLE SUR LE TEMPS

Time_loop: do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT))

  FLUX  = 0.D0
  GFLUX = 0.d0

  Umax = MAXVAL(SQRT((Prim(U_pv,1:Nx,1:Ny))**2.d0 + (Prim(V_pv,1:Nx,1:Ny))**2.d0))
  Cmax = MAXVAL(  Sound(1:Nx,1:Ny) )
    
  DT   = Dh/( Umax + Cmax)
  DT   = CFL*DT

IF( It == 1 ) WRITE(6,*) " Dt = ", Dt
 
  TIME = TIME + DT
!----------------------------------------------------------------------------------------------------
 		CALL Condition_lim_box(Prim, Ein, Nx, Ny, Pression, SOUND)
 
!----------------------------------------------------------------------------------------------------
 ! ETATS GAUCHE ET DROITE

DO ix    = 0, Nx
 DO iy   = 1, Ny

	PrimG(:) = Prim(:, ix,    iy)
  SOUNDG   = SOUND(ix, iy)
  PresG    = Pression(ix, iy)
  EinG     = Ein(ix, iy)

	PrimD(:)   = Prim(:,ix+1,  iy)
  SOUNDD     = SOUND(ix+1, iy)
  PresD      = Pression(ix+1, iy)
  EinD       = Ein(ix+1, iy)
	
!----------------------------------------------------------------------------------------------------
 ! APPEL SOLVEUR

  CALL HLLC( PrimG, PrimD, SOUNDG, SOUNDD, PresG, PresD, EinG, EinD, FLUXi, Gi, EPS)
  FLUX(:, ix, iy) = FLUXi(:)

ENDDO
ENDDO
!----------------------------------------------------------------------------------------------------
 !SHEMA DE GODUNOV
 CALL SHEMA_DE_GODUNOV_x(Nx, Ny, DX, DY, DT, CONS, FLUX, GFLUX)
!----------------------------------------------------------------------------------------------------
 
 IF( Ny > 1 ) THEN

  DO iy = 0, Ny
  DO ix = 1, Nx
  
    PrimG(:) = Prim(:, ix, iy)
    SOUNDG   = SOUND(ix, iy)
    PresG    = Pression(ix, iy)
    EinG     = Ein(ix, iy)

    PrimD(:)    = Prim(:, ix, iy+1)
    SOUNDD = SOUND(ix, iy+1)
    PresD  = Pression(ix, iy+1)
    EinD     = Ein(ix, iy+1)
  
!----------------------------------------------------------------------------------------------------
 ! APPEL SOLVEUR

CALL HLLC( PrimG, PrimD, SOUNDG, SOUNDD, PresG, PresD, EinG, EinD, FLUXi, Gi, EPS)
     GFLUX(:,ix, iy) = Gi(:)
   
  ENDDO
  ENDDO

ENDIF
!----------------------------------------------------------------------------------------------------
 !SHEMA DE GODUNOV

CALL SHEMA_DE_GODUNOV_y(Nx, Ny, DX, DY, DT, CONS, FLUX, GFLUX)

  
!----------------------------------------------------------------------------------------------------

CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND, EPS, Nx, Ny, CONS, g, Pression, Nv_Prim, angle, angle)


!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

IT = IT + 1
    IF (TIME2.LE.TIME) THEN
  	CALL Ecriture_donnees(X,Y, Prim, Ein, Pression, Nx,Ny, time, DX, DY )
    END IF

 
    IF (TIME2.LE.TIME) THEN
      PRINT*, 'EN', IT, 'ITERATIONS, ', ' TIME:', TIME
      TIME2 = TIME + period_time
    END IF

    CALL PutonScreen()


  ENDDO TIME_LOOP
!----------------------------------------------------------------------------------------------------
 ! FIN BOUCLE SUR LE TEMPS

 CALL Ecriture_donnees(X,Y, Prim, Ein, Pression, Nx,Ny, time, DX, DY )

 DEALLOCATE(Prim, Ein, Sound, X, Y, Pression, CONS, FLUX, GFLUX, MinVPression, MaxVp, NamesOfPrimVar)  

 CALL CPU_TIME(T2_CPU)

  PRINT*, 'L EXECUTION DU PROGRAMME A PRIS', T2_CPU - T1_CPU
  PRINT*, 'EN', IT-1, 'ITERATIONS, ', ' TIME:', TIME
 
STOP

CONTAINS
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE Calcule_vitesse_du_son()

END SUBROUTINE Calcule_vitesse_du_son
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

Subroutine HLLC_x_sub1(prim,flux,g,Nx,Ny)
implicit none
REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), flux(7,0:nx+1,0:ny+1)
real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
real(kind=dp) :: EL, ER
real(kind=dp) :: pstar,ustar,Estar,hstar
INTEGER :: ix,iy,K


do ix=1,nx+1
 do iy=1,ny
  !! etat gauche et droite
  ul=prim(U_pv,ix-1,iy);ur=prim(U_pv,ix,iy)
  hl=prim(h_pv,ix-1,iy);hr=prim(h_pv,ix,iy)
  vl=prim(v_pv,ix-1,iy);vr=prim(v_pv,ix,iy)
  p11l=prim(p11_pv,ix-1,iy);p11r=prim(p11_pv,ix,iy)
  p12l=prim(p11_pv,ix-1,iy);p12r=prim(p12_pv,ix,iy)
  p22l=prim(p22_pv,ix-1,iy);p22r=prim(p22_pv,ix,iy)
  !! calcul des énergie
  El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
  Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
  !! calcul des pression
  pl=g*hl*hl*0.5d0+hl*p11l
  pr=g*hr*hr*0.5d0+hr*p11r
  !! calcul des vitesses du son
  cl=dsqrt(g*hl+3.d0*p11l);cr=dsqrt(g*hr+3.d0*p11r)
  ! davis
  sr=dmax1(ul+cl,ur+cr)
  Sl=DMIN1(ul-cl,ur-cr)
  ! etat star 
  ml=hl*(ul-sl)
  mr=hr*(ur-sr)
  ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
  pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)


if (ustar.ge.0.d0)
if (sl.ge.0.d0) THEN
!!etat gauche
flux(1,ix,iy)=hl*ul
flux(2,ix,iy)=hl*ul*ul+pl
flux(3,ix,iy)=hl*ul*vl
flux(4,ix,iy)=p11l*ul !! equation inutile pour le cas x
flux(5,ix,iy)=p12l*ul
flux(6,ix,iy)=hl*p22l*ul
flux(7,ix,iy)=hl+El*ul+pl*ul
ELSE
!! etat star gauche
hstar=ml/(ustar-sl)
p12star=p12l*(ul-sl)/(ustar-sl)
Estar=El+(pl*ul-pstar*ustar)/ml
!! remplissage des flux
flux(1,ix,iy)=hstar*ustar
flux(2,ix,iy)=hstar*ustar*ustar+pstar
flux(3,ix,iy)=hstar*ustar*vl
flux(4,ix,iy)=p11l*ustar !! equation inutile pour le cas x
flux(5,ix,iy)=p12star*ustar
flux(6,ix,iy)=hstar*p22l*ustar
flux(7,ix,iy)=hstar+Estar*ustar+pstar*ustar
!!etat gauche etoile
end if
ELSE
if (sr.ge.0.d0)
!!etat droit etoile
!!etat star droit
hstar=mr/(ustar-sr)
p12star=p12r*(ul-sr)/(ustar-sr)
Estar=Er+(pr*ur-pstar*ustar)/mr
!remplissage des flux
flux(1,ix,iy)=hstar*ustar
flux(2,ix,iy)=hstar*ustar*ustar+pstar
flux(3,ix,iy)=hstar*ustar*vr
flux(4,ix,iy)=p11r*ustar !! equation inutile pour le cas x
flux(5,ix,iy)=p12star*ustar
flux(6,ix,iy)=hstar*p22r*ustar
flux(7,ix,iy)=hstar+Estar*ustar+pstar*ustar
ELSE
!!etat droit
flux(1,ix,iy)=hr*ur
flux(2,ix,iy)=hr*ur*ur+pr
flux(3,ix,iy)=hr*ur*vr
flux(4,ix,iy)=p11r*ur !! equation inutile pour le cas x
flux(5,ix,iy)=p12r*ur
flux(6,ix,iy)=hr*p22r*ur
flux(7,ix,iy)=hr+Er*ur+pr*ur
end if
end if


end do
end do
return
end subroutine HLLC_x_sub1

Subroutine HLLC_x_sub2(prim,flux,g,Nx,Ny)
implicit none
REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), flux(7,0:nx+1,0:ny+1)
real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
real(kind=dp) :: EL, ER
real(kind=dp) :: pstar,ustar,Estar,hstar
INTEGER :: ix,iy,K
do ix=1,nx+1
do iy=1,ny

end do
end do

return
end subroutine HLLC_x_sub2


  SUBROUTINE PutonScreen()

      REAL (KIND = DP)     :: MinVp(Nv_Prim), MaxVp(Nv_Prim)

      IF ( MOD(IT,ImpreE) == 0 ) THEN
              do iv = 1, Nv_Prim
              MinVp(iv)     = MINVAL(Prim(iv,1:Nx,1:Ny))
              MaxVp(iv)     = MAXVAL(Prim(iv,1:Nx,1:Ny))
              enddo
          
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
    END IF

  END SUBROUTINE PutonScreen

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE Ecriture_donnees(X,Y, Prim, Ein, Pression, Nx,Ny, time, DX, DY )
USE precisions
IMPLICIT NONE
INTEGER :: ix, iy, Nx, Ny , MyUnit = 30, il
REAL (KIND = DP) :: X(1:Nx), Y(1:Ny), Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
REAL (KIND = DP) ::   time, DX, DY
CHARACTER(LEN=3) :: NB, Zero="000"
! ECRITURE DES RESULTATS

    OPEN(19,FILE = 'test.out')
  
Do ix = 1, Nx
  DO iy = 1, Ny
     
      WRITE(19,10) X(ix), Y(iy),Prim(1,ix, iy),Prim(2,ix, iy),Prim(3,ix, iy),Prim(5,ix, iy),Prim(6,ix, iy), Ein(ix, iy), Pression(ix, iy)  
     END DO
  END DO
    WRITE(19,*)

  10 FORMAT(8(E20.13,1X))

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE LECTURE_DONNEES( angle, Nx,Ny, Lx, Ly, TIMEOUT, iterfinal, g, CFL, H_0, amplitude,  ImpreE, ImpreF)
USE precisions
IMPLICIT NONE
INTEGER :: Nx, Ny, iterfinal,   cond_lim, ImpreF, ImpreE
 REAL (KIND = DP) :: Lx, Ly, TIMEOUT,   H_0, H_1 , amplitude
 REAL (KIND = DP) :: g, CFL, angle , frottcoeff
 REAL (KIND = DP) :: dev, q0_jump
OPEN(UNIT=21, FILE = 'data_2D.inp', STATUS = 'OLD')
    READ(21,*) angle    ! inclination angle
  READ(21,*) Nx, Ny   ! NUMBER OF CELLS
  READ(21,*) Lx, Ly     ! DOMAIN LENGTH 
  READ(21,*) TIMEOUT  ! OUTPUT TIME
  READ(21,*) iterfinal ! Iteration final
    READ(21,*) g, CFL        ! acceleration due to gravity
    READ(21,*) H_0      ! stationary unstable solution (H_0, U_0, Phi_0 = 0)
    READ(21,*) amplitude ! amplitude des perturbation
    READ(21,*) ImpreE, ImpreF
 close(21)

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE INITIALISATION(X0, Y0, DX, DY, X, Y, Prim, SOUND, Ein, Nx, Ny,   CONS, angle, g, Pression, H_0, cond_lim, Lx, Ly, amplitude, Nv_Prim)
USE precisions
IMPLICIT NONE
INTEGER :: ix, iy, Nx, Ny,   Nv_Prim
INTEGER :: H_pv, U_pv, V_pv, Ein_pv, P11_pv,P12_pv,P22_pv, Sound_pv, Pres_pv
REAL (KIND = DP) :: X0, Y0, DY, DX,   angle, g, H_0 , Lx, Ly, amplitude
REAL (KIND = DP) :: X(1:Nx), Y(1:Ny), Prim(6, 0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1)
REAL (KIND = DP) :: Pression(0:Nx+1,0:Ny+1), CONS(6,1:Nx,1:Ny), SOUND(0:Nx+1,0:Ny+1)
REAL (KIND = DP) :: Pi

pi = 4.0d0*ATAN(1.0d0)
    
 DO ix = 1, Nx 
    X(ix) = X0 + 0.5D0*DX + (ix-1)*DX  

 DO iy = 1, Ny 
    Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY

        Prim(U_pv,(ix, iy) = 3.d0
        Prim(V_pv,ix, iy) = 0.d0
        Prim(H_pv, ix, iy) = H_0*(1.d0 +  amplitude*dsin(2.d0*Pi*x(ix)/Lx) ) 
        Prim(P11_pv,ix,iy) = 0.d0
        Prim(P12_pv,ix,iy) = 0.d0
        Prim(P22_pv,ix,iy) = 0.d0
        SOUND(ix,iy) = DSQRT(g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy))
        Pression(ix,iy) = g*H(ix,iy)*H(ix,iy)/2.d0 + Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
        Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + (Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy)) )/2.d0

enddo 
END DO

! VARIABLE CONSERVATIVES 
 DO ix = 1, Nx
  DO iy = 1, Ny 
   
  CONS(1,ix,iy) = Prim(H_pv, ix, iy)
  CONS(2,ix,iy) = Prim(H_pv, ix, iy)*Prim(U_pv, ix, iy)
  CONS(3,ix,iy) = Prim(H_pv, ix, iy)*Prim(V_pv, ix, iy)
  CONS(4,ix,iy) = Prim(H_pv, ix, iy)*(Ein(ix,iy)+ ((Prim(U_pv, ix, iy))**2.d0 + (Prim(V_pv, ix, iy))**2.d0 )/2.d0)
  CONS(5,ix,iy) = Prim(P12_pv, ix, iy)
  CONS(6,ix,iy) = Prim(H_pv, ix, iy)*Prim(P22_pv, ix, iy)

END DO
END DO
return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND, EPS, Nx, Ny, CONS, g, Pression, Nv_Prim, angle, angle)
USE precisions
 IMPLICIT NONE
 INTEGER :: ix,iy, Nx, Ny, Nv_Prim
 INTEGER :: H_pv, U_pv, V_pv, Ein_pv, P11_pv,P12_pv,P22_pv, Sound_pv, Pres_pv
 REAL (KIND = DP) :: EPS, phi2
 REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
 REAL (KIND = DP) :: SOUND(0:Nx+1,0:Ny+1) , CONS(6,1:Nx,1:Ny) 

 REAL (KIND = DP) :: angle, g 

! CALCUL NOUVELS VARIABLES PRIMITIVES
 DO ix = 1, Nx
   DO iy = 1, Ny
        Prim(H_pv,ix, iy) = CONS(1,ix, iy)
        Prim(U_pv,ix, iy) = CONS(2,ix, iy)/CONS(1,ix, iy)
        Prim(V_pv,ix, iy) = CONS(3,ix, iy)/CONS(1,ix, iy)
        Prim(P12_pv,ix, iy) = CONS(5,ix, iy)
        Prim(P22_pv,ix, iy) = CONS(6,ix, iy)/CONS(1,ix, iy)
        Prim(P11_pv,ix, iy) = 2.d0*CONS(4,ix, iy)/CONS(1,ix, iy) - ((Prim(U_pv,ix, iy))**2.d0+(Prim(V_pv,ix, iy))**2.d0 ) - g*Prim(H_pv,ix, iy) - Prim(P22_pv,ix, iy)
        Ein(ix, iy) = CONS(4,ix, iy)/CONS(1,ix, iy) - ((Prim(U_pv,ix, iy))**2.d0+(Prim(V_pv,ix, iy))**2.d0 )/2.d0
        SOUND(ix, iy) = DSQRT( g*Prim(H_pv,ix, iy) + 3.d0*P11(ix, iy)) 
        Pression(ix, iy) = g*(Prim(H_pv,ix, iy))**2.d0/2.d0 + Prim(P11_pv,ix, iy)*Prim(H_pv,ix, iy)
END DO  
END DO
return
END SUBROUTINE

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:



SUBROUTINE Condition_lim_box(Prim, Ein, Nx, Ny, Pression, SOUND )
USE precisions
IMPLICIT NONE
INTEGER :: Nx, Ny, ix, iy, nombre
REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), SOUND(0:Nx+1,0:Ny+1) , Pression(0:Nx+1,0:Ny+1)

  do iy = 1, Ny
   Prim(0, iy) = Prim(Nx, iy)
   Ein(0, iy) = Ein(Nx, iy)
   SOUND(0, iy) = SOUND(Nx, iy)
   Pression(0, iy) = Pression(Nx, iy)
  enddo

  do iy = 1, Ny
   Prim(Nx+1, iy) = Prim(1, iy)
   Ein(Nx+1, iy) = Ein(1, iy)
   SOUND(Nx+1, iy) = SOUND(1, iy)
   Pression(Nx+1, iy) = Pression(1, iy)
  enddo
  

    do ix = 1, Nx
   Prim(ix, 0) = Prim(ix, 1)
   Ein(ix, 0)      = Ein(ix, 1)
   SOUND(ix, 0)    = SOUND(ix, 1)
   Pression(ix, 0) = Pression(ix, 1)
   enddo


  do ix = 1, Nx
   Prim(ix, Ny+1)  = Prim(ix, Ny)
   Ein(ix, Ny+1)   = Ein(ix, Ny)
   SOUND(ix, Ny+1) = SOUND(ix, Ny)
   Pression(ix, Ny+1) = Pression(ix, Ny)
  enddo
return
END SUBROUTINE





!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE SHEMA_DE_GODUNOV_x(Nx, Ny, DX, DY, DT, CONS, FLUX, GFLUX)
USE precisions
IMPLICIT NONE
INTEGER :: ix, iy, Nx, Ny
REAL (KIND = DP) :: DX, DY, DT
REAL (KIND = DP) :: CONS(6,1:Nx,1:Ny), FLUX(6,0:Nx+1,0:Ny+1), GFLUX(6,0:Nx+1,0:Ny+1)

  DO ix =1,Nx
    do iy = 1, Ny
    CONS(:,ix, iy) = CONS(:,ix, iy) - DT/DX*(FLUX(:,ix, iy) - FLUX(:,ix-1, iy))     

  END DO
  ENDDO

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE SHEMA_DE_GODUNOV_y(Nx, Ny, DX, DY, DT, CONS, FLUX, GFLUX)
USE precisions
IMPLICIT NONE
INTEGER :: ix, iy, Nx, Ny
REAL (KIND = DP) :: DX, DY, DT
REAL (KIND = DP) :: CONS(6,1:Nx,1:Ny), FLUX(6,0:Nx+1,0:Ny+1), GFLUX(6,0:Nx+1,0:Ny+1)

  DO ix =1,Nx
    do iy = 1, Ny
    CONS(:,ix, iy) = CONS(:,ix, iy) - DT/DY*(GFLUX(:,ix, iy) - GFLUX(:,ix, iy-1))    

  END DO
  ENDDO

return
END SUBROUTINE


!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

SUBROUTINE HLLC( PrimL, PrimR, CL, CR, PL, PR, EL, ER, FHLLC, GHLLC, EPS, Phi1L, Phi1R)
USE precisions
IMPLICIT NONE
REAL (KIND = DP) :: HL, HR, UL, UR, Vl, VR, EL, ER, CL, CR, PR, PL, SRx, SLx, Sstarx, SRy, SLy, Sstary, UMAX, EPS, Phi1R, Phi1L
REAL (KIND = DP) :: CONSl(6), CONSr(6), FL(6), FR(6), GL(6), GR(6), FHLLC(6), GHLLC(6), CONSstarlx(6), CONSstarrx(6), CONSstarly(6), CONSstarry(6)

!---------------------------------------------------------------------------------------------------
 ! CALCUL SL ET SR (THE SLOWEST AND FASTEST SIGNALS) 
  SLX = DMIN1 ( UL - CL, UR - CR )
  SRX = DMAX1 ( UR + CR, UL + CL )
  SLY = DMIN1 ( Vl - CL, VR - CR )
  SRY = DMAX1 ( VR + CR, VL + CL )
! UMAX = DMAX1( UMAX, DABS(SL), DABS(SR))

  Sstarx = (PR - PL + HL*UL*(Slx - UL) - HR*UR*(SRx - UR))/(HL*(SLx - UL) - HR*(SRx - UR))
  Sstary = (PR - PL + HL*VL*(Sly - VL) - HR*VR*(SRy - VR))/(HL*(SLy - VL) - HR*(SRy - VR))
!----------------------------------------------------------------------------------------------------
 !CALCUL LES FLUX
 Consl(1) = HL
 Consl(2) = HL*UL
 Consl(3) = HL*VL
 Consl(4) = HL*( (UL*UL + VL*VL)/2.D0 + EL)
 

 Consr(1) = HR
 Consr(2) = HR*UR
 Consr(3) = HR*VR
 Consr(4) = HR*( (UR*UR + VR*VR)/2.D0 + ER)
 

 CONSstarlx(1) = HL*(SLx - UL)/(SLx - Sstarx)
 CONSstarlx(2) = HL*(SLx - UL)/(SLx - Sstarx)*Sstarx
 CONSstarlx(3) = HL*(SLx - UL)/(SLx - Sstarx)*VL
 CONSstarlx(4) = HL*(SLx - UL)/(SLx - Sstarx)*(EL + (UL*UL + VL*VL)/2.d0 + (Sstarx - UL)*(Sstarx + PL/(HL*(SLx-UL))))

 CONSstarrx(1) = HR*(SRx - UR)/(SRx - Sstarx)
 CONSstarrx(2) = HR*(SRx - UR)/(SRx - Sstarx)*Sstarx
 CONSstarrx(3) = HR*(SRx - UR)/(SRx - Sstarx)*VR
 CONSstarrx(4) = HR*(SRx - UR)/(SRx - Sstarx)*(ER + (UR*UR + VR*VR)/2.d0 +  (Sstarx - UR)*(Sstarx + PR/(HR*(SRx-UR))))

 CONSstarly(1) = HL*(SLy - VL)/(SLy - Sstary)
 CONSstarly(2) = HL*(SLy - VL)/(SLy - Sstary)*UL
 CONSstarly(3) = HL*(SLy - VL)/(SLy - Sstary)*Sstary
 CONSstarly(4) = HL*(SLy - VL)/(SLy - Sstary)*(EL + (UL*UL + VL*VL)/2.d0 + (Sstary - VL)*(Sstary + PL/(HL*(SLy-VL))))

 CONSstarry(1) = HR*(SRy - VR)/(SRy - Sstary)
 CONSstarry(2) = HR*(SRy - VR)/(SRy - Sstary)*UR
 CONSstarry(3) = HR*(SRy - VR)/(SRy - Sstary)*Sstary
 CONSstarry(4) = HR*(SRy - VR)/(SRy - Sstary)*(ER + (UR*UR + VR*VR)/2.d0 +  (Sstary - VR)*(Sstary + PR/(HR*(SRy-VR))))

  FL(1) = HL*UL
  FL(2) = HL*UL*UL + PL
  FL(3) = HL*UL*VL
  FL(4) = HL*UL*( (UL*UL + VL*VL)/2.D0 + EL) + PL*UL

  FR(1) = HR*UR
  FR(2) = HR*UR*UR + PR
  FR(3) = HR*UR*VR
  FR(4) = HR*UR*((UR*UR + VR*VR)/2.D0 + ER) + PR*UR


  GL(1) = HL*VL
  GL(2) = HL*UL*VL
  GL(3) = HL*VL*VL + PL
  GL(4) = HL*VL*((UL*UL + VL*VL)/2.D0 + EL) + PL*VL

  GR(1) = HR*VR
  GR(2) = HR*UR*VR 
  GR(3) = HR*VR*VR + PR
  GR(4) = HR*VR*((UR*UR + VR*VR)/2.D0 + ER) + PR*VR

  IF (DABS(Sstarx) .LT. EPS) THEN
    Sstarx = 0.D00
  END IF

  IF (DABS(Sstary) .LT. EPS) THEN
    Sstary = 0.D0
  END IF

   IF (SLx .GE. 0.d0 ) THEN
    FHLLC = FL
   ELSE IF ( SLx .LE. 0.d0 .AND. Sstarx .GE. 0.d0) THEN
    FHLLC = FL + SLx*(CONSstarlx - Consl)
   ELSE IF (Sstarx .LE. 0.d0 .AND. SRx .GE. 0.d0)  THEN
    FHLLC = FR + SRx*(CONSstarrx - Consr)
   ELSE IF (SRx .LE. 0.d0) THEN 
    FHLLC = FR    
   END IF

   IF (SLy .GE. 0.d0 ) THEN
    GHLLC = GL
   ELSE IF ( SLy .LE. 0.d0 .AND. Sstary .GE. 0.d0) THEN
    GHLLC = GL + SLy*(CONSstarly - Consl)
   ELSE IF (Sstary .LE. 0.d0 .AND. SRy .GE. 0.d0)  THEN
    GHLLC = GR + SRy*(CONSstarry - Consr)
   ELSE IF (SRy .LE. 0.d0) THEN 
    GHLLC = GR    
   END IF

RETURN
END SUBROUTINE 

END PROGRAM code2D


 
