  module precisions
  implicit none
  INTEGER, parameter :: dp=kind(1.0d0)
  Include 'mpif.h'
  end module precisions
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  MODULE GlobalParam
    USE precisions
    IMPLICIT NONE
    INTEGER                                         ::  ncpu, RANG, NBX, NBX0, numb, nombre
    INTEGER                                         ::  dest, source, tag, code
    INTEGER                                         ::  status(MPI_STATUS_SIZE)
    INTEGER                                         ::  ImpreE, ImpreF, Nv_Prim, argunit = 6, isave, penteX, penteY, test_shear  
    INTEGER                                         ::  Nx, Ny, iterfinal, cond_lim, oneDX, oneDY, twoD, method_source_term, WALLS
    INTEGER                                         ::  H_pv, U_pv, V_pv, P11_pv,  P12_pv, P22_pv
    REAL (KIND = DP)                                ::  X0, Y0, H_0, amplitude,  frottcoeff, disscoeff, Lx, Ly
    REAL (KIND = DP)                                ::  CFL, TIMEOUT, period_time, pi, angle, g, phi2
    REAL (KIND = DP), PARAMETER                     ::  EPS = 1.d-8
     REAL (KIND = DP)                                ::  lambda, gamma, beta
  END MODULE GlobalParam
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      

  PROGRAM code2D
  USE precisions
  use GlobalParam
    IMPLICIT NONE  
    INTEGER                                         :: iv, I, IT, ix, iy
    REAL (KIND = DP), ALLOCATABLE                   :: Ein(:,:), Sound_ax(:,:),Sound_bx(:,:), Sound_ay(:,:), Sound_by(:,:),  Pression(:,:)
    REAL (KIND = DP), ALLOCATABLE                   :: CONS(:,:,:), FLUX(:,:,:), Prim(:,:,:), pente(:,:,:), UmaxTampon(:,:), VmaxTampon(:,:)
    REAL (KIND = DP), ALLOCATABLE, DIMENSION(:)     :: MaxVP, MinVp, X, Y
    REAL (KIND = DP)                                :: T1_CPU, T2_CPU, TIME,  TIME2
    REAL (KIND = DP)                                :: DX, Dy, Dh, DT, dt2, tampondt
    REAL (KIND = DP)                                :: Hstar, Pstar, Ustar, Estar, Cmax
    REAL(KIND=DP)                                   :: xmax, XMIN, XMAX1, YMIN, YMAX, UMAX, Vmax, tampoUMAX, tampoVmax
    REAL(KIND=DP)                                   :: TAMPONXMIN, TAMPONXMAX, TAMPONYMIN, TAMPONYMAX

    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:)    :: NamesOfPrimVar
!////////////////////////UNIQUEMENT POUR LE PARALLELE//////////////////////////
!//Variables Calcul Parallele
!//Initialisation du calcul parallele
     CALL MPI_INIT(code)
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,code)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  !----------------------------------------------------------------------------------------------------
    pi = 4.0d0*ATAN(1.0d0)
    Nv_Prim  = 6
    H_pv     = 1
    U_pv     = 2
    V_pv     = 3
    P11_pv   = 4
    P12_pv   = 5
    P22_pv   = 6
  !---------------------------------------------------------------------------------------------------- 
    CALL CPU_TIME(T1_CPU)
  !----------------------------------------------------------------------------------------------------
    CALL LECTURE_DONNEES()
  !----------------------------------------------------------------------------------------------------
    DX = Lx/DFLOAT(Nx)
    Dy = Ly/DFLOAT(Ny)
    Dh =  dMIN1(Dx,Dy)
    X0 = 0.D0
    Y0 = 0.d0
    isave = -1
     IF(RANG .LT. NCPU-1) THEN
        NBX = int( Nx/NCPU )
        !NBY = int(Ny)
     ENDIF

      IF(RANG .EQ. NCPU-1) THEN
        NBX = Nx-(NCPU-1)*int(Nx/NCPU )
        !NBY = NY-(NCPU-1)*int(Ny/NCPU )
      ENDIF

      NBX0 = int( Nx/NCPU )
      !NBY0 = int( NY/NCPU )
      Nx = NBX 
      !Ny = NBY

  print*, 'Nx =', Nx,'Ny =', Ny, 'NBX =', NBX, 'NBX0', NBX0, 'RANG =', RANG, 'DX =', DX, 'DY =', DY
  !----------------------------------------------------------------------------------------------------
   !----------------------------------------------------------------------------------------------------
   ALLOCATE( Ein(0:Nx+1,0:Ny+1), Sound_ax(0:Nx+1,0:Ny+1),Sound_bx(0:Nx+1,0:Ny+1), Sound_ay(0:Nx+1,0:Ny+1),Sound_by(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1) )  
   ALLOCATE( Prim(6,0:Nx+1,0:Ny+1), CONS(7,1:Nx,1:Ny), FLUX(7,0:Nx,0:Ny), pente(6,0:Nx+1,0:Ny+1))  
   ALLOCATE( NamesOfPrimVar(Nv_Prim) , UmaxTampon(1:Nx,1:Ny), VmaxTampon(1:Nx,1:Ny))
   ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim), X(1:Nx), Y(1:Ny))
   flux(:,:,:)=0.d0; cons(:,:,:) = 0.d0; prim(:,:,:) = 0.d0; Sound_ay(:,:) = 0.d0; Sound_by(:,:) = 0.d0; Pression(:,:) = 0.d0
  !----------------------------------------------------------------------------------------------------
    NamesOfPrimVar(H_pv)         = "Depht Lenght"
    NamesOfPrimVar(U_pv)         = "Velocity (x)"
    NamesOfPrimVar(V_pv)         = "Velocity (y)"
    NamesOfPrimVar(P11_pv)       = "tensor P11"
    NamesOfPrimVar(P12_pv)       = "tensor P12"
    NamesOfPrimVar(P22_pv)       = "tensor P22"
  !---------------------------------------------------------------------------------------------------
  IF (RANG == 0) THEN 
    DO iv = 1, Nv_Prim
       WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ", iv, Nx, Ny, 'RANG = ', RANG
    END DO
    WRITE(6,*) " >>> End of LECTURE_DONNEES"
  ENDIF
  !----------------------------------------------------------------------------------------------------
  !INITIALISATION  
   CALL INITIALISATION(DX, DY, X, Y, Prim, SOUND_ax, SOUND_bx,SOUND_ay, SOUND_by,Ein,CONS, Pression)
  !----------------------------------------------------------------------------------------------------
   TIME = 0.D0; numb = 1
   period_time = 10.0d0;
   TIME2  = period_time ;
   IT = 1
  CALL PutonScreen()
  !----------------------------------------------------------------------------------------------------
   TAMPONXMIN = MINVAL(X(:))
   TAMPONXMAX = MAXVAL(X(:))

   TAMPONYMIN = MINVAL(Y(:))
   TAMPONYMAX = MAXVAL(Y(:))
  !----------------------------------------------------------------------------------------------------
   CALL MPI_ALLREDUCE(TAMPONXMIN, XMIN, numb, MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,code)
   CALL MPI_ALLREDUCE(TAMPONXMAX, XMAX1,numb, MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
 !----------------------------------------------------------------------------------------------------
   CALL MPI_ALLREDUCE(TAMPONYMIN, YMIN, numb, MPI_REAL8,MPI_MIN, MPI_COMM_WORLD,code)
   CALL MPI_ALLREDUCE(TAMPONYMAX, YMAX, numb, MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,code)
  !----------------------------------------------------------------------------------------------------
     IF (rang == 0) THEN
      WRITE(6,'( 2(A10,E15.6))')  " Xmin = ", XMIN, &
           &                      " Xmax = ", XMAX1
      WRITE(6,'( 2(A10,E15.6))') " Ymin = ",  YMIN, &
           &                     " Ymax = ",  YMAX
     ENDIF
  !----------------------------------------------------------------------------------------------------
   !----------------------------------------------------------------------------------------------------

  !  CALL Ecriture_donnees(X,Y, Prim, Ein, Pression, time, DX, DY )

   !-------------------------do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT)!-------------------
   ! BOUCLE SUR LE TEMPS

    Time_loop: do while ((IT .le. iterfinal) .and. (TIME .le. TIMEOUT))
    !----------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     do iy = 1, Ny
      UmaxTampon(ix, iy) = dmax1( Dabs(Prim(U_pv,ix,iy)+ Sound_ax(ix,iy)), Dabs(Prim(U_pv,ix,iy)- Sound_ax(ix,iy)) )
      VmaxTampon(ix, iy) = dmax1( Dabs(Prim(v_pv,ix,iy)+ Sound_ay(ix,iy)), Dabs(Prim(v_pv,ix,iy)- Sound_ay(ix,iy)) )
     enddo
    enddo
   tampoUmax = MAXVAL( UmaxTampon(1:Nx,1:Ny))
     CALL MPI_ALLREDUCE(tampoUmax,Umax, numb,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,code)
    tampovmax = MAXVAL(VmaxTampon(1:Nx,1:Ny))
    CALL MPI_ALLREDUCE(tampovmax,vmax, numb,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,code)
   
   !----------------------------------------------------------------------------------------------------    
   if (ny > 1) then 
      DT   = dmin1(Dx/DABS(Umax), Dy/DABS(Vmax) )
  else
    DT   = Dx/DABS(Umax)
  endif
      tampondt=dt
      CALL MPI_ALLREDUCE(tampondt,dt, numb,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)
      DT   = CFL*DT; dt2 = 0.5d0*dt
   !----------------------------------------------------------------------------------------------------
       IF (RANG == 0) THEN 
        IF( It == 1 ) WRITE(6,*) " Dt = ", Dt, 'RANG =',  rang
       ENDIF

       TIME = TIME + DT
   !----------------------------------------------------------------------------------------------------
 !                    if (method_source_term == 1) then 
!                    CALL euler_method(DT2, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
!                 else if (method_source_term == 2) then 
!                    CALL rk2(DT2, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
!                 else if (method_source_term == 3) then 
!                    CALL rk4(DT2, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
!                 endif
!               !----------------------------------------------------------------------------------------------------
!                  CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
!              !----------------------------------------------------------------------------------------------------
!                call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------
               IF (cond_lim == 1) THEN 
                  CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
               ELSE IF (cond_lim == 2) THEN 
                  CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
               ELSE IF (cond_lim == 3) THEN 
                  CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
              ELSEIF (cond_lim == 5) THEN 
                  CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
               ENDIF


     if (nx.gt.1) then
            
             CALL PENTE1_x(Prim, PENTE)
            !PENTE = 0.D0
            !----------------------------------------------------------------------------------------------------
             call HLLC_x_sub1(prim,flux, pente, cons, it, dt2) 

            !----------------------------------------------------------------------------------------------------
             call godunov_x_sub1(cons,flux,dt2,dx)

            !----------------------------------------------------------------------------------------------------
             CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
             call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------

              IF (cond_lim == 1) THEN 
                CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
              ELSE IF (cond_lim == 2) THEN 
               CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
              ELSE IF (cond_lim == 3) THEN 
               CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
              ELSEIF (cond_lim == 5) THEN 
                CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
              ENDIF
             
             CALL PENTE2_x(Prim, PENTE)

             call HLLC_x_sub2(prim,flux, cons, pente, dt2, it) 

             call godunov_x_sub2(cons,flux,dt2,dx)
            

             CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
             !----------------------------------------------------------------------------------------------------
              call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------


             IF (cond_lim == 1) THEN 
              CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 2) THEN 
              CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 3) THEN 
              CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
             ELSEIF (cond_lim == 5) THEN 
                CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
             ENDIF
    endif
  !----------------------------------------------------------------------------------------------------
  
  if (ny.gt.1) then
              CALL PENTE1_y( Prim, PENTE)
              call HLLC_y_sub1(prim,flux, cons, pente, it, dt2, dy)
              call godunov_y_sub1(cons,flux,dt2,dy)
              CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
            !----------------------------------------------------------------------------------------------------
              call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
            !----------------------------------------------------------------------------------------------------

             IF (cond_lim == 1) THEN 
              CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 2) THEN 
              CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 3) THEN 
              CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
            ELSEIF (cond_lim == 5) THEN 
                CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
             ENDIF

             CALL PENTE2_y( Prim, PENTE)
             call HLLC_y_sub2(prim,flux, cons, pente, it, dt2, dy) 
             call godunov_y_sub2(cons,flux,dt2,dy)
             CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
             call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------

             IF (cond_lim == 1) THEN 
              CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 2) THEN 
              CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 3) THEN 
              CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
             ELSEIF (cond_lim == 5) THEN 
              CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
             ENDIF
  end if

  if (ny.gt.1) then
              CALL PENTE1_y( Prim, PENTE)
              call HLLC_y_sub1(prim,flux, cons, pente, it, dt2, dy)
              call godunov_y_sub1(cons,flux,dt2,dy)
              CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
            !----------------------------------------------------------------------------------------------------
              call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
            !----------------------------------------------------------------------------------------------------

             IF (cond_lim == 1) THEN 
              CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 2) THEN 
              CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 3) THEN 
              CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
            ELSEIF (cond_lim == 5) THEN 
                CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
             ENDIF

             CALL PENTE2_y( Prim, PENTE)
             call HLLC_y_sub2(prim,flux, cons, pente, it, dt2, dy) 
             call godunov_y_sub2(cons,flux,dt2,dy)
             CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
             call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------

             IF (cond_lim == 1) THEN 
              CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 2) THEN 
              CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 3) THEN 
              CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
             ELSEIF (cond_lim == 5) THEN 
              CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
             ENDIF
  end if


     if (nx.gt.1) then
            
             CALL PENTE1_x(Prim, PENTE)
            !PENTE = 0.D0
            !----------------------------------------------------------------------------------------------------
             call HLLC_x_sub1(prim,flux, pente, cons, it, dt2) 

            !----------------------------------------------------------------------------------------------------
             call godunov_x_sub1(cons,flux,dt2,dx)

            !----------------------------------------------------------------------------------------------------
             CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
             call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------

              IF (cond_lim == 1) THEN 
                CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
              ELSE IF (cond_lim == 2) THEN 
               CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
              ELSE IF (cond_lim == 3) THEN 
               CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
              ELSEIF (cond_lim == 5) THEN 
                CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
              ENDIF
             
             CALL PENTE2_x(Prim, PENTE)
             !print*, pente(4,:,:),pente(5,:,:),pente(6,:,:)

             call HLLC_x_sub2(prim,flux, cons, pente, dt2, it) 

             call godunov_x_sub2(cons,flux,dt2,dx)
            

             CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
             !----------------------------------------------------------------------------------------------------
              call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
             !----------------------------------------------------------------------------------------------------


             IF (cond_lim == 1) THEN 
              CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 2) THEN 
              CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
             ELSE IF (cond_lim == 3) THEN 
              CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
             ELSEIF (cond_lim == 5) THEN 
                CALL cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
             ENDIF
    endif
  !----------------------------------------------------------------------------------------------------

!    if (method_source_term == 1) then 
!     CALL euler_method(DT, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
!    else if (method_source_term == 2) then 
!     CALL rk2(DT, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
!     else if (method_source_term == 3) then 
!     CALL rk4(DT, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
!    endif


  !  CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  !----------------------------------------------------------------------------------------------------
 !   call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )

!             IF (cond_lim == 1) THEN 
!               CALL Condition_lim_box(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
!              ELSE IF (cond_lim == 2) THEN 
!               CALL CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
!              ELSE IF (cond_lim == 3) THEN 
!               CALL Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
!              ENDIF
    IT = IT + 1
  !----------------------------------------------------------------------------------------------------
     IF (TIME2.LE.TIME) THEN
     CALL PutonScreen()
   	!CALL Ecriture_donnees(X,Y, Prim, Ein, Pression, time, DX, DY )
     PRINT*, 'EN', IT, 'ITERATIONS, ', ' TIME:', TIME
       TIME2 = TIME + period_time
    END IF

    ENDDO TIME_LOOP
  !----------------------------------------------------------------------------------------------------
   ! FIN BOUCLE SUR LE TEMPS

   CALL Ecriture_donnees(X,Y, Prim, Ein, Pression, time, DX, DY )
  !----------------------------------------------------------------------------------------------------
   DEALLOCATE(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, X, Y, Pression, CONS, FLUX, MinVP, MaxVp, NamesOfPrimVar)  
  !----------------------------------------------------------------------------------------------------
   CALL CPU_TIME(T2_CPU)
  !----------------------------------------------------------------------------------------------------
    PRINT*, 'L EXECUTION DU PROGRAMME A PRIS', T2_CPU - T1_CPU
    PRINT*, 'EN', IT-1, 'ITERATIONS, ', ' TIME:', TIME
  !----------------------------------------------------------------------------------------------------
    STOP
    CALL MPI_FINALIZE(CODE)
  !----------------------------------------------------------------------------------------------------
  CONTAINS
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  !!!suivant x

  Subroutine HLLC_x_sub1(prim,flux, pente, cons, it, dt)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix,iy, it
  real(kind=dp) :: dt, dt2
  REAL (KIND = DP) :: cons(7,1:Nx,1:Ny)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED
  REAL (KIND = DP) :: FLUX(7,0:Nx,0:Ny)
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), pente(6,0:Nx+1,0:Ny+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:, :) :: EinG, EinD, PressionG, PressionD
  real(kind=dp)    :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
  real(kind=dp)    :: EL, ER,p22l,p22r,vr,vl
  real(kind=dp)    :: pstar,ustar,Estar,hstar,p12star

  ALLOCATE(consG(7,1:Nx,1:Ny),consD(7,1:Nx,1:Ny), FLUXG(7,0:Nx,0:Ny),FLUXD(7,0:Nx,0:Ny),PrimG(Nv_Prim,0:Nx+1,0:Ny+1), PrimD(Nv_Prim,0:Nx+1,0:Ny+1), CONS_PRED(7,1:Nx,1:Ny))
  ALLOCATE(EinG(0:Nx+1,0:Ny+1), EinD(0:Nx+1,0:Ny+1), PressionG(0:Nx+1,0:Ny+1), PressionD(0:Nx+1,0:Ny+1))
  dt2 = 0.5d0*dt; PrimG(:,:,:) = 0.D0; PrimD(:,:,:)= 0.D0; FLUXG(:,:,:) = 0.D0; FLUXD(:,:,:)= 0.D0; CONSG(:,:,:) = 0.D0; CONSD(:,:,:)= 0.D0
  EinG= 0.D0; EinD= 0.D0; PressionG= 0.D0; PressionD= 0.D0; CONS_PRED = 0.D0

    do ix = 1, Nx
      do iy = 1, Ny
          do iv = 1, Nv_Prim

              PrimG(iv, ix, iy) = Prim(iv, ix, iy) - 0.5d0*Pente(iv, ix,iy)
              PrimD(iv, ix, iy) = Prim(iv, ix, iy) + 0.5d0*Pente(iv, ix,iy)
          enddo
             EinG(ix, iy) = 0.5d0*(g*PrimG(H_pv, ix, iy) + PrimG(p11_pv, ix, iy)+PrimG(p22_pv, ix, iy))
             EinD(ix, iy) = 0.5d0*(g*PrimD(H_pv, ix, iy) + PrimD(p11_pv, ix, iy)+PrimD(p22_pv, ix, iy))

             PressionG(ix, iy) = 0.5d0*g*PrimG(H_pv, ix, iy)**2.d0 + PrimG(H_pv, ix, iy)*PrimG(p11_pv, ix, iy)
             PressionD(ix, iy) = 0.5d0*g*PrimD(H_pv, ix, iy)**2.d0 + PrimD(H_pv, ix, iy)*PrimD(p11_pv, ix, iy)    
      enddo
    enddo


   IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
        CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
   ENDIF
 
CALL PRIM_TO_CONS_FLUX_sub1x(PrimG, EinG,PressionG, CONSG, FLUXG)
CALL PRIM_TO_CONS_FLUX_sub1x(PrimD, EinD,PressionD, CONSD, FLUXD)

CALL CONS_PREDICTION_sub1x( DX, DT, CONSG, CONS_PRED, FLUXG, FLUXD)

!  if (method_source_term == 1) then 
!   CALL euler_method(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  else if (method_source_term == 2) then 
!   CALL rk2(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!   else if (method_source_term == 3) then 
!   CALL rk4(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  endif

   
   IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
   ENDIF

  CALL PRIM_TO_CONS_FLUX_sub1x(PrimG, EinG,PressionG, CONSG, FLUXG)
  CALL CONS_PREDICTION_sub1x(DX, DT, CONSD, CONS_PRED, FLUXG, FLUXD)

!  if (method_source_term == 1) then 
!   CALL euler_method(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  else if (method_source_term == 2) then 
!   CALL rk2(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  else if (method_source_term == 3) then 
!   CALL rk4(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  endif

 CALL NOUVELL_VARIABLE_PRIM(PrimG, EinG, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionG, it)
!----------------------------------------------------------------------------------------------------
 call COMMUNICATION(PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!----------------------------------------------------------------------------------------------------

  CALL NOUVELL_VARIABLE_PRIM(PrimD, EinD, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionD, it)
!----------------------------------------------------------------------------------------------------
  call COMMUNICATION(PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it )
!----------------------------------------------------------------------------------------------------


  IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF

   do ix=0,nx
   do iy=1,ny
    !! etat gauche et droite
    ul=PrimD(U_pv,ix,iy);  ur=primG(U_pv,ix+1,iy)
    hl=PrimD(h_pv,ix,iy);  hr=primG(h_pv,ix+1,iy)
    vl=PrimD(v_pv,ix,iy);  vr=primG(v_pv,ix+1,iy)

    p11l=PrimD(p11_pv,ix,iy);  p11r=primG(p11_pv,ix+1,iy)
    p12l=PrimD(p12_pv,ix,iy);  p12r=primG(p12_pv,ix+1,iy)
    p22l=PrimD(p22_pv,ix,iy);  p22r=primG(p22_pv,ix+1,iy)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !! calcul des pression
    pl=g*hl*hl*0.5d0+hl*p11l
    pr=g*hr*hr*0.5d0+hr*p11r
    !! calcul des vitesses du son
    cl=dsqrt(g*hl+3.d0*p11l);cr=dsqrt(g*hr+3.d0*p11r)
    ! davis
    !sr=dmax1(ul+cl,ur+cr)
    !Sl=DMIN1(ul-cl,ur-cr)

    
    Sl=ul-cl; if (ur-cr < sl) sl = ur - cr
    sr=ur+cr; if(ul+cl> sr) sr = ul+cl

    ml=hl*(ul-sl)
    mr=hr*(ur-sr)
    ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
    pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)

  if (ustar.ge.0.d0) then
   if (sl.ge.0.d0) THEN
  !!etat gauche
   flux(1,ix,iy)=hl*ul
   flux(2,ix,iy)=hl*ul*ul+pl
   flux(3,ix,iy)=hl*ul*vl
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12l*ul
   flux(6,ix,iy)=hl*p22l*ul
   flux(7,ix,iy)=hl*El*ul+pl*ul
   ELSE
   !! etat star gauche
   hstar=ml/(ustar-sl)
  ! p12star=p12l*(ul-sl)/(ustar-sl)
   p12star = p12l*hstar/hl
   Estar=El+(pl*ul-pstar*ustar)/ml
   !! remplissage des flux

   flux(1,ix,iy)=hstar*ustar!hl*ul + sl*(hstar - hl)
   flux(2,ix,iy)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
   flux(3,ix,iy)=hstar*ustar*vl!hl*ul*vl + sl*vl*(hstar - hl)   
   flux(4,ix,iy)=0.d0    !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar! + sl*(p12star - p12l)
   flux(6,ix,iy)=hstar*p22l*ustar!hl*p22l*ul + sl*p22l*(hstar - hl)
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar!hl*El*ul + pl*ul + sl*( hstar*Estar - hl*El ) 
   
   endif
  ELSE
   if (sr.ge.0.d0) then
  !!etat droit etoile
  !!etat star droit
   hstar=mr/(ustar-sr)
  ! p12star=p12r*(ul-sr)/(ustar-sr)
   p12star = p12r*hstar/hr
   Estar=Er+(pr*ur-pstar*ustar)/mr
   !remplissage des flux

   flux(1,ix,iy)=hstar*ustar!hl*ul + sl*(hstar - hl)
   flux(2,ix,iy)=hstar*ustar*ustar+pstar! hl*ul*ul + pl + sl*(hstar*ustar - hl*ul)   
   flux(3,ix,iy)=hstar*ustar*vr!hl*ul*vl + sl*vl*(hstar - hl)   
   flux(4,ix,iy)=0.d0    !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar! + sl*(p12star - p12l)
   flux(6,ix,iy)=hstar*p22r*ustar!hl*p22l*ul + sl*p22l*(hstar - hl)
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar!hl*El*ul + pl*ul + sl*( hstar*Estar - hl*El ) 
   
   ELSE
  !!etat droit
   flux(1,ix,iy)=hr*ur
   flux(2,ix,iy)=hr*ur*ur+pr
   flux(3,ix,iy)=hr*ur*vr
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12r*ur
   flux(6,ix,iy)=hr*p22r*ur
   flux(7,ix,iy)=hr*Er*ur+pr*ur
   end if
  end if

  end do
  end do

  DEALLOCATE(consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED, EinG, EinD, PressionG, PressionD )
  return
  end subroutine HLLC_x_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Subroutine HLLC_x_sub2(prim,flux, cons, pente, dt, it)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
    INTEGER :: ix,iy, it
  real(kind=dp) :: dt, dt2
  REAL (KIND = DP) :: cons(7,1:Nx,1:Ny)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED
  REAL (KIND = DP) :: FLUX(7,0:Nx,0:Ny)
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), pente(6,0:Nx+1,0:Ny+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:, :) :: EinG, EinD, PressionG, PressionD
  real(kind=dp) :: ul, ur, hl, hr, p11l, p11r, p12r, p12l, ml, mr, sr, sl, p22l, p22r
  real(kind=dp) :: EL, ER,vr,vl
  real(kind=dp) :: vstar,Estar,hp12star
  ALLOCATE(consG(7,1:Nx,1:Ny),consD(7,1:Nx,1:Ny), FLUXG(7,0:Nx,0:Ny),FLUXD(7,0:Nx,0:Ny),PrimG(6,0:Nx+1,0:Ny+1), PrimD(6,0:Nx+1,0:Ny+1), CONS_PRED(7,1:Nx,1:Ny))
  ALLOCATE(EinG(0:Nx+1,0:Ny+1), EinD(0:Nx+1,0:Ny+1), PressionG(0:Nx+1,0:Ny+1), PressionD(0:Nx+1,0:Ny+1))
  dt2 = 0.5d0*dt; PrimG(:,:,:) = 0.D0; PrimD(:,:,:)= 0.D0; FLUXG(:,:,:) = 0.D0; FLUXD(:,:,:)= 0.D0; CONSG(:,:,:) = 0.D0; CONSD(:,:,:)= 0.D0
  EinG= 0.D0; EinD= 0.D0; PressionG= 0.D0; PressionD= 0.D0; CONS_PRED = 0.D0

do ix = 1, Nx
  do iy = 1, Ny
     do iv = 1, Nv_Prim
        PrimG(iv, ix, iy) = Prim(iv, ix, iy) - 0.5d0*Pente(iv, ix,iy)
        PrimD(iv, ix, iy) = Prim(iv, ix, iy) + 0.5d0*Pente(iv, ix,iy)
     enddo
       EinG(ix, iy) = 0.5d0*(g*PrimG(H_pv, ix, iy) + PrimG(p11_pv, ix, iy)+PrimG(p22_pv, ix, iy))
       EinD(ix, iy) = 0.5d0*(g*PrimD(H_pv, ix, iy) + PrimD(p11_pv, ix, iy)+PrimD(p22_pv, ix, iy))

       PressionG(ix, iy) = 0.5d0*g*PrimG(H_pv, ix, iy)**2.d0 + PrimG(H_pv, ix, iy)*PrimG(p11_pv, ix, iy)
       PressionD(ix, iy) = 0.5d0*g*PrimD(H_pv, ix, iy)**2.d0 + PrimD(H_pv, ix, iy)*PrimD(p11_pv, ix, iy)
  enddo
enddo
    
   IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
        CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF
!----------------------------------------------------------------------------------------------------

 CALL PRIM_TO_CONS_FLUX_sub2x(PrimG, EinG, CONSG, FLUXG)
 CALL PRIM_TO_CONS_FLUX_sub2x(PrimD, EinD, CONSD, FLUXD)
!----------------------------------------------------------------------------------------------------
 CALL CONS_PREDICTION_sub2x( DX, DT, CONSG, CONS_PRED, FLUXG, FLUXD)
!----------------------------------------------------------------------------------------------------
 !  if (method_source_term == 1) then 
!   CALL euler_method(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  else if (method_source_term == 2) then 
!   CALL rk2(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  else if (method_source_term == 3) then 
!   CALL rk4(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  endif


  IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
  ELSEIF (cond_lim == 5) THEN 
    CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF
!----------------------------------------------------------------------------------------------------

 CALL PRIM_TO_CONS_FLUX_sub2x(PrimG, EinG, CONSG, FLUXG)
 CALL CONS_PREDICTION_sub2x(DX, DT, CONSD, CONS_PRED, FLUXG, FLUXD)

!---------------------------------------------------------------------------------------------------- 
 !  if (method_source_term == 1) then 
!   CALL euler_method(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  else if (method_source_term == 2) then 
!   CALL rk2(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  else if (method_source_term == 3) then 
!   CALL rk4(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  endif

!----------------------------------------------------------------------------------------------------
  CALL NOUVELL_VARIABLE_PRIM(PrimG, EinG, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionG, it)
!----------------------------------------------------------------------------------------------------
  CALL COMMUNICATION(PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!----------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------- 
  CALL NOUVELL_VARIABLE_PRIM(PrimD, EinD, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionD, it)
  call COMMUNICATION(PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it )
!----------------------------------------------------------------------------------------------------

  IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
  ELSEIF (cond_lim == 5) THEN 
    CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF

  do ix=0,nx
  do iy=1,ny
    ul=PrimD(U_pv,ix,iy);     ur=primG(U_pv,ix+1,iy)
    hl=PrimD(h_pv,ix,iy);     hr=primG(h_pv,ix+1,iy)
    vl=PrimD(v_pv,ix,iy);     vr=primG(v_pv,ix+1,iy)
    p11l=PrimD(p11_pv,ix,iy); p11r=primG(p11_pv,ix+1,iy)
    p12l=PrimD(p12_pv,ix,iy); p12r=primG(p12_pv,ix+1,iy)
    p22l=PrimD(p22_pv,ix,iy); p22r=primG(p22_pv,ix+1,iy)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !sr=dsqrt(dmax1(p11l,p11r)+1.d-8)
     !sr=dsqrt(p11r)
     sl= - dsqrt(p11l)
     sr=dsqrt(p11r)
    ! vstar=(hl*p12l-hr*p12r+sr*(hl*vl+hr*vr))/(sr*(hr+hl))
     !hp12star=(hr*hl)/(hr+hl)*(p12l+p12r+sr*(vl-vr))
   !   sr=dmax1(dsqrt(p11l),dsqrt(p11r)) 
   !   sl=dmin1(-dsqrt(p11l),-dsqrt(p11r))

  vstar         = (hl*(p12l- sl*vl) - hr*(p12r-sr*vr))/(sr*hr - hl*sl)
  hp12star      = (hr*hl)/(hr*sr-hl*sl)*(sr*p12l-sl*p12r+sl*sr*(vr-vl))
  Flux(1,ix,iy) = 0.d0
  Flux(2,ix,iy) = 0.d0
  flux(3,ix,iy) = hp12star
  Flux(4,ix,iy) = 0.d0
  Flux(5,ix,iy) = vstar !! non conservatif
  flux(6,ix,iy) = 2.d0*hp12star*vstar !!inutile
  flux(7,ix,iy) = hp12star*vstar

  end do
  end do

  DEALLOCATE(EinG, EinD, PressionG, PressionD, consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED )
  return
  end subroutine HLLC_x_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_x_sub1(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  real(KIND=dp) :: cons(7,1:nx,1:ny), FLUX(7,0:Nx,0:Ny)
  real(Kind=dp) :: hu, h, hv, he, hp22, dt, dx
  INTEGER :: k, ix, iy

  do ix=1,nx
   do iy=1,ny
!------------------------------------------------------------------------------------------------------------------------
    do k=1,7
    cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
   end do
!------------------------------------------------------------------------------------------------------------------------
 !! corection de la variable p11
  hu=cons(2,ix,iy);hv=cons(3,ix,iy);h=cons(1,ix,iy);hp22=cons(6,ix,iy);hE=cons(7,ix,iy)
  cons(4,ix,iy)=2.d0*hE-g*h*h-hp22-(hu*hu+hv*hv)/h
  end do
  end do
  end subroutine godunov_x_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  subroutine godunov_x_sub2(cons,flux,dt,dx)
  USE precisions
  use GlobalParam
  implicit none
  real(KIND=dp) :: cons(7,1:nx,1:ny), FLUX(7,0:Nx,0:Ny)
  real(Kind=dp) :: hu,h,hv,he,hp11,p11,dt,dx
  INTEGER :: k,ix,iy

  do ix=1,nx
   do iy=1,ny
        p11=cons(4,ix,iy)/cons(1,ix,iy)
       !------------------------------------------------------------------------------------------------------------------------
        do k=1,4
          cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
        end do
       !------------------------------------------------------------------------------------------------------------------------
          cons(5,ix,iy)=cons(5,ix,iy)+dt/dx*(Flux(5,ix-1,iy)-flux(5,ix,iy))*p11 !! ici dans le flux on a stocké vstar
       !------------------------------------------------------------------------------------------------------------------------
         do k=6,7
          cons(k,ix,iy)=cons(k,ix,iy)+dt/dx*(Flux(k,ix-1,iy)-flux(k,ix,iy))
         end do
       !------------------------------------------------------------------------------------------------------------------------
        !! corection de la variable p22
        hu=cons(2,ix,iy); hv=cons(3,ix,iy); h=cons(1,ix,iy); hp11=cons(4,ix,iy); hE=cons(7,ix,iy)
        cons(6,ix,iy)=2.d0*hE-(hu*hu+hv*hv)/h-g*h*h-hp11
   end do
  end do
  end subroutine godunov_x_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  !!! suivant y

  Subroutine HLLC_y_sub1(prim,flux, cons, pente, it, dt, dy)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER :: ix,iy, it
  real(kind=dp) :: dt, dt2, dy
  REAL (KIND = DP) :: cons(7,1:Nx,1:Ny)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED
  REAL (KIND = DP) :: FLUX(7,0:Nx,0:Ny)
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), pente(6,0:Nx+1,0:Ny+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:, :) :: EinG, EinD, PressionG, PressionD
  real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
  real(kind=dp) :: EL, ER,p22l,p22r,vr,vl,temp
  real(kind=dp) :: pstar,ustar,Estar,hstar,p12star
  ALLOCATE(consG(7,1:Nx,1:Ny),consD(7,1:Nx,1:Ny), FLUXG(7,0:Nx,0:Ny),FLUXD(7,0:Nx,0:Ny),PrimG(6,0:Nx+1,0:Ny+1), PrimD(6,0:Nx+1,0:Ny+1), CONS_PRED(7,1:Nx,1:Ny))
  ALLOCATE(EinG(0:Nx+1,0:Ny+1), EinD(0:Nx+1,0:Ny+1), PressionG(0:Nx+1,0:Ny+1), PressionD(0:Nx+1,0:Ny+1))
  dt2 = 0.5d0*dt; PrimG(:,:,:) = 0.D0; PrimD(:,:,:)= 0.D0; FLUXG(:,:,:) = 0.D0; FLUXD(:,:,:)= 0.D0; CONSG(:,:,:) = 0.D0; CONSD(:,:,:)= 0.D0
  EinG(:, :) = 0.D0; EinD(:, :)= 0.D0; PressionG(:, :) = 0.D0; PressionD(:, :) = 0.D0; CONS_PRED(:,:,:) = 0.D0

do ix = 1, Nx
  do iy = 1, Ny
     do iv = 1, Nv_Prim
        PrimG(iv, ix, iy) = Prim(iv, ix, iy) - 0.5d0*Pente(iv, ix,iy)
        PrimD(iv, ix, iy) = Prim(iv, ix, iy) + 0.5d0*Pente(iv, ix,iy)
     enddo
       EinG(ix, iy) = 0.5d0*(g*PrimG(H_pv, ix, iy) + PrimG(p11_pv, ix, iy)+PrimG(p22_pv, ix, iy))
       EinD(ix, iy) = 0.5d0*(g*PrimD(H_pv, ix, iy) + PrimD(p11_pv, ix, iy)+PrimD(p22_pv, ix, iy))

       PressionG(ix, iy) = 0.5d0*g*PrimG(H_pv, ix, iy)**2.d0 + PrimG(H_pv, ix, iy)*PrimG(p11_pv, ix, iy)
       PressionD(ix, iy) = 0.5d0*g*PrimD(H_pv, ix, iy)**2.d0 + PrimD(H_pv, ix, iy)*PrimD(p11_pv, ix, iy)
  enddo
enddo
    
   IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
        CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF
 
CALL PRIM_TO_CONS_FLUX_sub1y(PrimG, EinG, CONSG, FLUXG)
CALL PRIM_TO_CONS_FLUX_sub1y(PrimD, EinD, CONSD, FLUXD)

CALL CONS_PREDICTION_sub1y( Dy, DT, CONSG, CONS_PRED, FLUXG, FLUXD)

!  if (method_source_term == 1) then 
!   CALL euler_method(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  else if (method_source_term == 2) then 
!   CALL rk2(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  else if (method_source_term == 3) then 
!   CALL rk4(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!  endif


 IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
 ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
 ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
 ELSEIF (cond_lim == 5) THEN 
    CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
 ENDIF

CALL PRIM_TO_CONS_FLUX_sub1y(PrimG, EinG,CONSG, FLUXG)

CALL CONS_PREDICTION_sub1y(Dy, DT, CONSD, CONS_PRED, FLUXG, FLUXD)

!  if (method_source_term == 1) then 
!   CALL euler_method(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  else if (method_source_term == 2) then 
!   CALL rk2(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  else if (method_source_term == 3) then 
!   CALL rk4(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!  endif


CALL NOUVELL_VARIABLE_PRIM(PrimG, EinG, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionG, it)
!----------------------------------------------------------------------------------------------------
 call COMMUNICATION(PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it )
!----------------------------------------------------------------------------------------------------
 CALL NOUVELL_VARIABLE_PRIM(PrimD, EinD, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionD, it)
!----------------------------------------------------------------------------------------------------
 call COMMUNICATION(PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it )
!----------------------------------------------------------------------------------------------------

  IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
  ELSEIF (cond_lim == 5) THEN 
    CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF

  do ix=1,nx
   do iy=0,ny
    !! etat gauche et droite
    ul=PrimD(v_pv,ix,iy);   ur=PrimG(v_pv,ix,iy+1)
    hl=PrimD(h_pv,ix,iy);   hr=PrimG(h_pv,ix,iy+1)
    vl=PrimD(U_pv,ix,iy);   vr=primG(U_pv,ix,iy+1)

    p11l=PrimD(p22_pv,ix,iy); p11r=primG(p22_pv,ix,iy+1)
    p12l=PrimD(p12_pv,ix,iy); p12r=primG(p12_pv,ix,iy+1)
    p22l=PrimD(p11_pv,ix,iy); p22r=primG(p11_pv,ix,iy+1)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !! calcul des pression
    pl=g*hl*hl*0.5d0+hl*p11l
    pr=g*hr*hr*0.5d0+hr*p11r
    !! calcul des vitesses du son
    cl=dsqrt(g*hl+3.d0*p11l);cr=dsqrt(g*hr+3.d0*p11r)
    ! davis
!     sr=dmax1(ul+cl,ur+cr)
!     Sl=DMIN1(ul-cl,ur-cr)

    Sl=ul-cl; if (ur-cr < sl) sl = ur - cr
    sr=ur+cr; if(ul+cl> sr) sr = ul+cl
    ! etat star 
    ml=hl*(ul-sl)
    mr   =hr*(ur-sr)
    ustar=(ul*ml-ur*mr+pl-pr)/(ml-mr)
    pstar=((ul-ur)*mr*ml+mr*pl-ml*pr)/(mr-ml)


  if (ustar.ge.0.d0) then
   if (sl.ge.0.d0) THEN
   !!etat gauche
   flux(1,ix,iy)=hl*ul
   flux(2,ix,iy)=hl*ul*ul+pl
   flux(3,ix,iy)=hl*ul*vl
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12l*ul
   flux(6,ix,iy)=hl*p22l*ul
   flux(7,ix,iy)=hl*El*ul+pl*ul
   ELSE

   !! etat star gauche
   hstar=ml/(ustar-sl)
   !p12star=p12l*(ul-sl)/(ustar-sl)
    p12star = p12l*hstar/hl
   Estar=El+(pl*ul-pstar*ustar)/ml
   !! remplissage des flux
   flux(1,ix,iy)=hstar*ustar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar
   flux(3,ix,iy)=hstar*ustar*vl
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar
   flux(6,ix,iy)=hstar*p22l*ustar
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar
   endif
  ELSE
   if (sr.ge.0.d0) then
  !!etat droit etoile
  !!etat star droit
   hstar=mr/(ustar-sr)
   !p12star=p12r*(ul-sr)/(ustar-sr)
    p12star = p12r*hstar/hr
   Estar=Er+(pr*ur-pstar*ustar)/mr
   !remplissage des flux
   flux(1,ix,iy)=hstar*ustar
   flux(2,ix,iy)=hstar*ustar*ustar+pstar
   flux(3,ix,iy)=hstar*ustar*vr
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12star*ustar
   flux(6,ix,iy)=hstar*p22r*ustar
   flux(7,ix,iy)=hstar*Estar*ustar+pstar*ustar
   ELSE
  !!etat droit
   flux(1,ix,iy)=hr*ur
   flux(2,ix,iy)=hr*ur*ur+pr
   flux(3,ix,iy)=hr*ur*vr
   flux(4,ix,iy)=0.d0 !! equation inutile pour le cas x
   flux(5,ix,iy)=p12r*ur
   flux(6,ix,iy)=hr*p22r*ur
   flux(7,ix,iy)=hr*Er*ur+pr*ur
   end if
  end if
  !! inversion des flux :
  temp=flux(2,ix,iy)
  flux(2,ix,iy)=flux(3,ix,iy)  
  flux(3,ix,iy)=temp
  temp=flux(6,ix,iy)
  flux(6,ix,iy)=flux(4,ix,iy)  
  flux(4,ix,iy)=temp
  end do
  end do
  DEALLOCATE(consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED, EinG, EinD, PressionG, PressionD )
  return
  end subroutine HLLC_y_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Subroutine HLLC_y_sub2(prim,flux, cons, pente, it, dt, dy)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER :: ix,iy, it, iv
  real(kind=dp) :: dt, dt2, dy
  REAL (KIND = DP) :: cons(7,1:Nx,1:Ny)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED
  REAL (KIND = DP) :: FLUX(7,0:Nx,0:Ny)
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), pente(6,0:Nx+1,0:Ny+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:, :) :: EinG, EinD, PressionG, PressionD
  real(kind=dp) :: ul,ur,hl,hr,pr,pl,p11l,p11r,p12r,p12l,cl,cr,ml,mr,sl,sr
  real(kind=dp) :: EL, ER,p22l,p22r,vr,vl,temp
  real(kind=dp) :: pstar,vstar,Estar,hstar,hp12star
  ALLOCATE(consG(7,1:Nx,1:Ny),consD(7,1:Nx,1:Ny), FLUXG(7,0:Nx,0:Ny),FLUXD(7,0:Nx,0:Ny),PrimG(6,0:Nx+1,0:Ny+1), PrimD(6,0:Nx+1,0:Ny+1), CONS_PRED(7,1:Nx,1:Ny))
  ALLOCATE(EinG(0:Nx+1,0:Ny+1), EinD(0:Nx+1,0:Ny+1), PressionG(0:Nx+1,0:Ny+1), PressionD(0:Nx+1,0:Ny+1))
  dt2 = 0.5d0*dt; PrimG(:,:,:) = 0.D0; PrimD(:,:,:)= 0.D0; FLUXG(:,:,:) = 0.D0; FLUXD(:,:,:)= 0.D0; CONSG(:,:,:) = 0.D0; CONSD(:,:,:)= 0.D0
  EinG(:, :)= 0.D0; EinD(:, :) = 0.D0; PressionG(:, :)= 0.D0; PressionD(:, :)= 0.D0; CONS_PRED(:,:,:) = 0.D0

  do ix = 1, Nx
    do iy = 1, Ny
        do iv = 1, Nv_Prim
          PrimG(iv, ix, iy) = Prim(iv, ix, iy) - 0.5d0*Pente(iv, ix,iy)
          PrimD(iv, ix, iy) = Prim(iv, ix, iy) + 0.5d0*Pente(iv, ix,iy)
        enddo
           EinG(ix, iy) = 0.5d0*(g*PrimG(H_pv, ix, iy) + PrimG(p11_pv, ix, iy)+PrimG(p22_pv, ix, iy))
           EinD(ix, iy) = 0.5d0*(g*PrimD(H_pv, ix, iy) + PrimD(p11_pv, ix, iy)+PrimD(p22_pv, ix, iy))

           PressionG(ix, iy) = 0.5d0*g*PrimG(H_pv, ix, iy)**2.d0 + PrimG(H_pv, ix, iy)*PrimG(p11_pv, ix, iy)
           PressionD(ix, iy) = 0.5d0*g*PrimD(H_pv, ix, iy)**2.d0 + PrimD(H_pv, ix, iy)*PrimD(p11_pv, ix, iy)
    enddo
  enddo
    
   IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
  ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
        CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF
 
  CALL PRIM_TO_CONS_FLUX_sub2y(PrimG, EinG, CONSG, FLUXG)
  CALL PRIM_TO_CONS_FLUX_sub2y(PrimD, EinD, CONSD, FLUXD)

  CALL CONS_PREDICTION_sub2y( Dy, DT, CONSG, CONS_PRED, FLUXG, FLUXD)

!     if (method_source_term == 1) then 
!       CALL euler_method(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!     else if (method_source_term == 2) then 
!       CALL rk2(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!     else if (method_source_term == 3) then 
!       CALL rk4(DT2, CONS_PRED, PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it)
!     endif


  IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
  ELSEIF (cond_lim == 5) THEN 
        CALL cond_lim_turb_sol(PrimG, EinG, PressionG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF

  CALL PRIM_TO_CONS_FLUX_sub2y(PrimG, EinG, CONSG, FLUXG)

  CALL CONS_PREDICTION_sub2y(Dy, DT, CONSD, CONS_PRED, FLUXG, FLUXD)
    
!     if (method_source_term == 1) then 
!       CALL euler_method(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!     else if (method_source_term == 2) then 
!       CALL rk2(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!     else if (method_source_term == 3) then 
!       CALL rk4(DT2, CONS_PRED, PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it)
!     endif

   CALL NOUVELL_VARIABLE_PRIM(PrimG, EinG, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionG, it)
  !----------------------------------------------------------------------------------------------------
   call COMMUNICATION(PrimG, EinG, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionG, it )
  !----------------------------------------------------------------------------------------------------

   CALL NOUVELL_VARIABLE_PRIM(PrimD, EinD, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_PRED, PressionD, it)
  !----------------------------------------------------------------------------------------------------
   call COMMUNICATION(PrimD, EinD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,PressionD, it )
  !----------------------------------------------------------------------------------------------------

  IF (cond_lim == 1) THEN 
      CALL Condition_lim_box(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
   ELSE IF (cond_lim == 2) THEN 
     CALL CondLimABSORPTION(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  ELSE IF (cond_lim == 3) THEN 
     CALL Condition_lim_batteur(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, dy)
  ELSEIF (cond_lim == 5) THEN 
    CALL cond_lim_turb_sol(PrimD, EinD, PressionD, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by, TIME, X,Y)
  ENDIF
 
  do ix=1,nx
  do iy=0,ny
    ul=PrimD(v_pv,ix,iy);     ur=primG(v_pv,ix,iy+1)
    hl=PrimD(h_pv,ix,iy);     hr=primG(h_pv,ix,iy+1)
    vl=PrimD(u_pv,ix,iy);     vr=primG(u_pv,ix,iy+1)
    p22l=PrimD(p11_pv,ix,iy); p22r=primG(p11_pv,ix,iy+1)
    p12l=PrimD(p12_pv,ix,iy); p12r=primG(p12_pv,ix,iy+1)
    p11l=PrimD(p22_pv,ix,iy); p11r=primG(p22_pv,ix,iy+1)
    !! calcul des énergie
    El=(ul*ul+vl*vl+g*hl+p11l+p22l)*0.5D0
    Er=(ur*ur+vr*vr+g*hr+p11r+p22r)*0.5D0
    !sr=dqrt(dmax1(p11l,p11r))
    !sr=dmax1(dsqrt(p11l),dsqrt(p11r)) 
    !sl=dmin1(-dsqrt(p11l),-dsqrt(p11r))
     sr=dsqrt(p11r)
     sl = -dsqrt(p11l)
     !vstar=(hl*p12l-hr*p12r+sr*(hl*vl+hr*vr))/(sr*(hr+hl))
     !hp12star=(hr*hl)/(hr+hl)*(p12l+p12r+sr*(vl-vr))
     
      vstar         = (hl*(p12l- sl*vl) - hr*(p12r-sr*vr))/(sr*hr - hl*sl)
     hp12star      = (hr*hl)/(hr*sr-hl*sl)*(sr*p12l-sl*p12r+sl*sr*(vr-vl))
    Flux(1,ix,iy)=0.d0
    Flux(3,ix,iy)=0.d0
    flux(2,ix,iy)=hp12star
    Flux(6,ix,iy)=0.d0
    Flux(5,ix,iy)=vstar !! non conservatif
    flux(4,ix,iy)=2.d0*hp12star*vstar !!inutile
    flux(7,ix,iy)=hp12star*vstar
  end do
  end do

  DEALLOCATE(consG,consD, FLUXG,FLUXD, PrimG, PrimD, CONS_PRED, EinG, EinD, PressionG, PressionD )
  return
  end subroutine HLLC_y_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_y_sub1(cons,flux,dt,dy)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  real(KIND=dp) :: cons(7,1:nx,1:ny),flux(7,0:nx,0:ny)
  real(Kind=dp) :: hu,h,hv,he,hp11,dt,dy
  INTEGER :: k,ix,iy

  do ix=1,nx
   do iy=1,ny

    do k=1,7
    cons(k,ix,iy)=cons(k,ix,iy)+dt/dy*(Flux(k,ix,iy-1)-flux(k,ix,iy))
   end do

  !! corection de la variable p22
  hu=cons(2,ix,iy); hv=cons(3,ix,iy); h=cons(1,ix,iy); hp11=cons(4,ix,iy); hE=cons(7,ix,iy)
  cons(6,ix,iy)=2.d0*hE-g*h*h-hp11-(hu*hu+hv*hv)/h
  end do
  end do
  end subroutine godunov_y_sub1
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  subroutine godunov_y_sub2(cons,flux,dt,dy)
  USE precisions
  USE GlobalParam
  implicit none
  real(KIND=dp) :: cons(7,1:nx,1:ny), flux(7,0:nx,0:ny)
  real(Kind=dp) :: hu, h, hv, he, hp22, p22, dt, dy
  INTEGER :: k, ix, iy

  do ix=1,nx
   do iy=1,ny
       p22=cons(6,ix,iy)/cons(1,ix,iy)
       do k=1,4
        cons(k,ix,iy)=cons(k,ix,iy)+dt/dy*(Flux(k,ix,iy-1)-flux(k,ix,iy))
       end do
        cons(5,ix,iy)=cons(5,ix,iy)+dt/dy*(Flux(5,ix,iy-1)-flux(5,ix,iy))*p22
       do k=6,7
        cons(k,ix,iy)=cons(k,ix,iy)+dt/dy*(Flux(k,ix,iy-1)-flux(k,ix,iy))
       end do
        !! corection de la variable p11
        hu=cons(2,ix,iy);hv=cons(3,ix,iy);h=cons(1,ix,iy);hp22=cons(6,ix,iy);hE=cons(7,ix,iy)
        cons(4,ix,iy)=2.d0*hE-(hu*hu+hv*hv)/h-g*h*h-hp22
   end do
  end do

  end subroutine godunov_y_sub2
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE PutonScreen()
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER              :: nbre
  REAL (KIND = DP)     :: MinVp(Nv_Prim), MaxVp(Nv_Prim)
  REAL (KIND = DP)     :: tamponMinVp_h, tamponMinVp_u, tamponMinVp_v, tamponMinVp_p11, tamponMinVp_p12, tamponMinVp_p22
  REAL (KIND = DP)     :: tamponMaxVp_h, tamponMaxVp_u, tamponMaxVp_v, tamponMaxVp_p11, tamponMaxVp_p12, tamponMaxVp_p22
  nbre = 1
    tamponMinVp_h       = MINVAL(Prim(H_pv,1:Nx,1:Ny))
    tamponMinVp_u       = MINVAL(Prim(u_pv,1:Nx,1:Ny))
    tamponMinVp_v       = MINVAL(Prim(v_pv,1:Nx,1:Ny))
    tamponMinVp_p11     = MINVAL(Prim(p11_pv,1:Nx,1:Ny))
    tamponMinVp_p12     = MINVAL(Prim(p12_pv,1:Nx,1:Ny))
    tamponMinVp_p22     = MINVAL(Prim(p22_pv,1:Nx,1:Ny))


    tamponMaxVp_h       = MAXVAL(Prim(H_pv,1:Nx,1:Ny)) 
    tamponMaxVp_u       = MAXVAL(Prim(u_pv,1:Nx,1:Ny)) 
    tamponMaxVp_v       = MAXVAL(Prim(v_pv,1:Nx,1:Ny)) 
    tamponMaxVp_p11     = MAXVAL(Prim(p11_pv,1:Nx,1:Ny)) 
    tamponMaxVp_p12     = MAXVAL(Prim(p12_pv,1:Nx,1:Ny)) 
    tamponMaxVp_p22     = MAXVAL(Prim(p22_pv,1:Nx,1:Ny)) 


    CALL MPI_ALLREDUCE(tamponMinVp_h,MinVp(H_pv),      nbre,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMinVp_u,MinVp(u_pv),      nbre,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMinVp_v,MinVp(v_pv),      nbre,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMinVp_p11,MinVp(p11_pv),  nbre,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMinVp_p12,MinVp(p12_pv),  nbre,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMinVp_p22,MinVp(p22_pv),  nbre,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,code)


    CALL MPI_ALLREDUCE(tamponMaxVp_h,MaxVp(H_pv),      nbre,MPI_REAL8,MPI_Max,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMaxVp_u,MaxVp(u_pv),      nbre,MPI_REAL8,MPI_Max,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMaxVp_v,MaxVp(v_pv),      nbre,MPI_REAL8,MPI_Max,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMaxVp_p11,MaxVp(p11_pv),  nbre,MPI_REAL8,MPI_Max,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMaxVp_p12,MaxVp(p12_pv),  nbre,MPI_REAL8,MPI_Max,MPI_COMM_WORLD,code)
    CALL MPI_ALLREDUCE(tamponMaxVp_p22,MaxVp(p22_pv),  nbre,MPI_REAL8,MPI_Max,MPI_COMM_WORLD,code)


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
  ENDSUBROUTINE PutonScreen
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE Ecriture_donnees(X,Y, Prim, Ein, Pression, time, DX, DY )
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy , MyUnit = 30, il
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:):: err_p11,err_p12,err_p22,  err_h,err_u,err_v
  REAL (KIND = DP) ::error_p11,error_p12,error_p22, error_h,error_u,error_v
  REAL (KIND = DP) :: p11,p12,p22,h,u,v 
  REAL (KIND = DP) :: time, DX, DY
  CHARACTER(LEN=3) :: NB, Zero="000"
  CHARACTER*3 :: numrg

allocate(err_p11(1:Nx,1:Ny),err_p12(1:Nx,1:Ny),err_p22(1:Nx,1:Ny), err_h(1:Nx,1:Ny),err_u(1:Nx,1:Ny),err_v(1:Nx,1:Ny))
  error_p11=0.d0;error_p12=0.d0;error_p22=0.d0; error_h=0.d0;error_u=0.d0;error_v=0.d0
  err_p11(1:Nx,1:Ny)=0.d0;err_p12(1:Nx,1:Ny)=0.d0;err_p22(1:Nx,1:Ny)=0.d0
  err_h(1:Nx,1:Ny)=0.d0;err_u(1:Nx,1:Ny)=0.d0;err_v(1:Nx,1:Ny)=0.d0
  p11=0.d0;p12=0.d0;p22=0.d0; h=0.d0; u=0.d0; v=0.d0

  WRITE(numrg,'(i3.3)') rang
if (cond_lim == 5) then
       OPEN(MyUnit+3,FILE = './resu/err_turb_test.out')
    
   Do ix = 1, Nx
     DO iy = 1, Ny

      h   = 1.d0/( 1.d0 + beta**2.d0*time**2.d0 )
        u = beta*(beta*time*x(ix) + y(iy))/(1.d0 + beta**2.d0*time**2.d0)
        v = beta*(-x(ix) + beta*time*y(iy))/(1.d0 + beta**2.d0*time**2.d0)
        
        p11 = (lambda + gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
        p12 = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
        p22 = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0

        
         if (dabs(p12).le.1.d-8) then 
           p12=0.d0
          endif
           p22 = dmax1(p22, 1.d-8)
           p11 = dmax1(p11, 1.d-8)
          err_h(ix,iy) = dabs(h - Prim(H_pv,ix, iy) )/h
     
          err_u(ix,iy) = dabs( (u - Prim(u_pv,ix, iy) )/u)
        
          err_v(ix,iy) = dabs( (v- Prim(v_pv,ix, iy) )/v )

         err_p11(ix,iy) = dabs(p11 - Prim(p11_pv,ix, iy) )
 
         err_p12(ix,iy) = dabs(p12 - Prim(p12_pv,ix, iy) )
  
        err_p22(ix,iy) = dabs(p22 - Prim(p22_pv,ix, iy) )

        END DO
    END DO
 
       error_h = SUM(err_h(:, :))/(nx*ny) !maxval(err_h(:, :))
       error_u = SUM(err_u(:, :))/(nx*ny) !maxval(err_u(:, :))
       error_v = SUM(err_v(:, :))/(nx*ny) !maxval(err_v(:, :))

       error_p11 =  SUM(err_p11(:, :))/(nx*ny) !  maxval(err_p11(:, :)) !
       error_p12 =  SUM(err_p12(:, :))/(nx*ny) ! maxval(err_p12(:, :)) !
       error_p22 =  SUM(err_p22(:, :))/(nx*ny) ! maxval(err_p22(:, :)) !

       !WRITE(MyUnit+3,'(4(E20.13,1X))') DLOG(DX/(Lx) ), DLOG(error_p11/lambda),DLOG(error_p12/(gamma+lambda)), DLOG(error_p22/gamma)
       WRITE(MyUnit+3,'(4(E20.13,1X))')  DLOG(DX/(Lx) ), DLOG(error_h),DLOG(error_u), DLOG(error_v)
endif

! if (twoD == 0) then
!    OPEN(MyUnit+1,FILE = "./resu/case1_300p_"//numrg//".out")
!     Do ix = 1, Nx
!       DO iy = 1, Ny
!           WRITE(MyUnit+1,11) X(ix), Y(iy), Prim(H_pv,ix, iy), Prim(U_pv,ix, iy), Prim(V_pv,ix, iy),prim(p11_pv,ix,iy), Prim(p12_pv,ix, iy), Prim(p22_pv,ix, iy) 
!       END DO
!     END DO
!    WRITE(MyUnit+1,*)
!    11 FORMAT(8(E20.13,1X))
! endif

! IF (twoD == 1) then 
!                      isave = isave + 1
!                     WRITE(unit=NB, fmt="(I3)") isave
!                     NB    = ADJUSTL(NB)
!                     il    = LEN_TRIM(NB) 

!                    WRITE(6,*) " FILE = ", "./resu/OR2case1_box_500p300p_cfl05_1m3_05mNONSTRANGEULER_"//numrg//"_"//Zero(1:3-il)//TRIM(NB)//".vtk" 
!                    OPEN(UNIT=MyUnit, FILE="./resu/OR2case1_box_500p300p_cfl05_1m3_05mNONSTRANGEULER_"//numrg//"_"//Zero(1:3-il)//TRIM(NB)//".vtk" )
!                    WRITE(MyUnit,'(''# vtk DataFile Version 2.0'')')
!                    WRITE(MyUnit,'(''Rectilinear 3D Dataset'')')
!                    WRITE(MyUnit,'(''ASCII'')')
!                    WRITE(MyUnit,'(''           '')')
!                    WRITE(MyUnit,'(''DATASET STRUCTURED_POINTS'')')
!                    WRITE(MyUnit,FMT='(''DIMENSIONS'',I8,I8,I8)') Nx+1, Ny+1, 2
!                    WRITE(MyUnit,FMT='(''ORIGIN'',3(E11.4,1x))') 0.d0, 0.d0, 0.d0
!                    WRITE(MyUnit,FMT='(''SPACING'',3(E11.4,1x))') dx,dy, 0.0001d0
!                    WRITE(MyUnit,*) ' '
!                    WRITE(MyUnit,FMT='(''CELL_DATA '',I9)') Nx*Ny*1
!                    WRITE(MyUnit,*) ' '

!                    WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'depth'
!                    WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
!                   DO iy=1,Ny ; DO ix=1,Nx
!                      WRITE(MyUnit,'(G11.4)') Prim(H_pv,ix, iy)      
!                    ENDDO ; ENDDO 

!                    WRITE(MyUnit,FMT='(''VECTORS '',A12, '' float'')') 'Vitesse(m/s)'
!                    DO iy=1,Ny ; DO ix=1,Nx
!                      WRITE(MyUnit,'(3(E11.4,1x))') Prim(U_pv,ix, iy), Prim(V_pv,ix, iy), 0.d0
!                    ENDDO ; ENDDO 

!                    WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P11'
!                    WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
!                   DO iy=1,Ny ; DO ix=1,Nx
!                      WRITE(MyUnit,'(G11.4)') Prim(P11_pv,ix, iy)       
!                    ENDDO ; ENDDO 

!                    WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P22'
!                    WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
!                   DO iy=1,Ny ; DO ix=1,Nx
!                      WRITE(MyUnit,'(G11.4)') Prim(P22_pv,ix, iy)         
!                    ENDDO ; ENDDO 

!                    WRITE(MyUnit,FMT='(''SCALARS '',A6, '' float 1'')') 'P12'
!                    WRITE(MyUnit,'(''LOOKUP_TABLE default'')')     
!                   DO iy=1,Ny ; DO ix=1,Nx
!                      WRITE(MyUnit,'(G11.4)') Prim(P12_pv,ix, iy)     
!                    ENDDO ; ENDDO 

!                    WRITE(MyUnit,*) ' '
!                 close(MyUnit)

!               !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
               
!                   WRITE(unit=NB, fmt="(I3)") isave
!                   NB    = ADJUSTL(NB)
!                   il    = LEN_TRIM(NB) 
!                   WRITE(6,*) " FILE = ", "./resu/OR2case1_box_500p300p_cfl05_1m3_05mNONSTRANGEULER_"//numrg//"_"//Zero(1:3-il)//TRIM(NB)//".tp"
!                   OPEN(UNIT=MyUnit+2, FILE="./resu/OR2case1_box_500p300p_cfl05_1m3_05mNONSTRANGEULER_"//numrg//"_"//Zero(1:3-il)//TRIM(NB)//".tp" )
!                   WRITE(MyUnit+2,'(A)') 'TITLE="This is a title"'  
!                   WRITE(MyUnit+2,'(A)') 'VARIABLES= "X", "Y" , "H","U", "V", "P11","P12", "P22" '
!                   WRITE(MyUnit+2,*) 'ZONE I=', NX,', J=', Ny,'DATAPACKING=POINT' 
                 
!                  DO  ix=1, Nx
!                   DO  iy=1,Ny
                     
!                         WRITE (MyUnit+2,'(9(E16.8,1x))') X(ix), Y(iy), Prim(H_pv,ix, iy), Prim(U_pv,ix, iy), Prim(V_pv,ix, iy),prim(p11_pv,ix,iy), &
!                         &                                                                           Prim(p12_pv,ix, iy), prim(p22_pv,ix,iy), &
!                         &                                                        (prim(p11_pv,ix,iy)*prim(p22_pv,ix,iy) - Prim(p12_pv,ix, iy)**2.d0)/Prim(H_pv,ix, iy)**2.d0
!                      END DO
!                   END DO


!               !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

!                    WRITE(6,*) " FILE = ", "./resu/OR2case1_box_500p300p_cfl05_1m3_05mNONSTRANGEULER_"//numrg//"_"//Zero(1:3-il)//TRIM(NB)//".csv"
!                    OPEN(UNIT=MyUnit+3, FILE="./resu/OR2case1_box_500p300p_cfl05_1m3_05mNONSTRANGEULER_"//numrg//"_"//Zero(1:3-il)//TRIM(NB)//".csv"  )
                   
!                     WRITE(MyUnit+3,'(A)') '"X", "Y",  "H"'  
                
!                 DO  ix=1, Nx 
!                   DO  iy=1,Ny
                    
!                           WRITE (MyUnit+3 ,*)  X(ix), ',' , Y(iy), ',' , &   
!                               &         Prim(H_pv,ix, iy)*100.d0
                      
!                      END DO
!                    END DO
!   !               !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

!   !                    WRITE(6,*) " FILE = ", "./resu/P11surH2_case2_300p_100p_1m045m"
!   !                    OPEN(UNIT=MyUnit+4, FILE="./resu/P11surH2_case2_300p_100p_1m045m.csv" )
                     
!   !                     WRITE(MyUnit+4,'(A)') '"X", "Y",  "P11/H^2"'  
                  
!   !                 DO  ix=1, Nx 
!   !                   DO  iy=1,Ny
                      
!   !                           WRITE (MyUnit+4 ,*)  X(ix), ',' , Y(iy), ',' , &   
!   !                               &         prim(p11_pv,ix,iy)/Prim(H_pv,ix, iy)**2.d0
                        
!   !                      END DO
!   !                    END DO
!   !               !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
!   !                  WRITE(6,*) " FILE = ", "./resu/entropie_case2_300p_100p_1m045m"
!   !                    OPEN(UNIT=MyUnit+5, FILE="./resu/entropie_case2_300p_100p_1m045m.csv" )
                     
!   !                     WRITE(MyUnit+5,'(A)') '"X", "Y",  "P11/H^2"'  
                  
!   !                 DO  ix=1, Nx 
!   !                   DO  iy=1,Ny
                      
!   !                           WRITE (MyUnit+5 ,*)  X(ix), ',' , Y(iy), ',' , &   
!   !                               &         (prim(p11_pv,ix,iy)*prim(p11_pv,ix,iy) - Prim(p12_pv,ix, iy)**2.d0)/Prim(H_pv,ix, iy)**2.d0
                        
!   !                      END DO
!   !                    END DO
!   !               !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
!   ENDIF 
  DEALLOCATE(err_p11,err_p12,err_p22)
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE LECTURE_DONNEES()
  USE precisions
  use GlobalParam
    IMPLICIT NONE
   !REAL (KIND = DP) :: Lx, Ly
   OPEN(UNIT=21, FILE = 'data_2D.inp', STATUS = 'OLD')
    READ(21,*) penteX, penteY, WALLS   
    READ(21,*) twoD, oneDX, oneDY      
    READ(21,*) cond_lim, method_source_term, test_shear ! 1 box/ 2 absorbtion/ 3 batteur/ 4 jump; cond_lim ;   test_shear
    READ(21,*) angle                 ! inclination angle
    READ(21,*) Nx, Ny                 ! NUMBER OF CELLS
    READ(21,*) Lx, Ly                 ! DOMAIN LENGTH 
    READ(21,*) TIMEOUT                ! OUTPUT TIME
    READ(21,*) iterfinal              ! Iteration final
    READ(21,*) g, CFL                 ! acceleration due to gravity
    READ(21,*) H_0                    ! stationary unstable solution (H_0, U_0, Phi_0 = 0)
    READ(21,*) phi2                   ! patit enstrophy
    READ(21,*) frottcoeff, disscoeff  ! cf, cr
    READ(21,*) amplitude              ! amplitude des perturbation
    READ(21,*) ImpreE, ImpreF
    READ(21,*) lambda, gamma, beta
    close(21)
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE INITIALISATION(DX, DY, X, Y, Prim, SOUND_ax, SOUND_bx,SOUND_ay, SOUND_by,Ein,CONS, Pression)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER :: ix, iy
  REAL (KIND = DP) :: DY, DX, U_0 !Lx, Ly,
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: Prim(6, 0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Pression(0:Nx+1,0:Ny+1), CONS(7,1:Nx,1:Ny), Ein(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1), SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1), SOUND_by(0:Nx+1,0:Ny+1)
  REAL(KIND=DP) :: X_centre,Y_centre,a,b, theta, x1, y1
  X_centre = 0.5D0; Y_centre = 0.5D0; a = 0.2D0; b = 0.2d0; theta = 0.d0
 
  pi = 4.0d0*ATAN(1.0d0)
  U_0 = DSQRT(g*Dtan(angle)*H_0/frottcoeff) 
  
  DO ix = 1, Nx
      X(ix) = X0 + 0.5D0*DX + (ix-1)*DX  + nbx0*rang*dx  
   DO iy = 1, Ny
      Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY

          IF (cond_lim == 1 ) THEN ! ! cond_lim: 1 box/ 2 absorption / 3 batteur/ 4 hydraulic jump
             
          IF (twoD == 1) THEN 
            Prim(H_pv, ix, iy) = H_0*(1.d0 + amplitude*dsin(2.d0*Pi*y(iy)/Ly)+  amplitude*dsin(2.d0*Pi*x(ix)/Lx) )
          ELSE IF (oneDX == 1)  THEN 
            Prim(H_pv, ix, iy) = H_0*(1.d0 + amplitude*dsin(2.d0*Pi*x(ix)/Lx) )
          ELSE IF (oneDY == 1) THEN 
            Prim(H_pv, ix, iy) = H_0*(1.d0 +  amplitude*dsin(2.d0*Pi*y(iy)/Ly))
          ENDIF

          ELSE IF (cond_lim == 2 ) THEN 
             Prim(H_pv, ix, iy) = H_0
          ELSE IF (cond_lim == 3 ) THEN 
             Prim(H_pv, ix, iy)= H_0
          ELSEIF (cond_lim == 5 ) THEN 
             Prim(H_pv, ix, iy)= 1.d0
             Prim(U_pv,ix, iy) = beta*y(iy)
             Prim(V_pv,ix, iy) = -beta*x(ix)
             Prim(P11_pv,ix,iy) = lambda 
             Prim(P12_pv,ix,iy) =  0.d0
             Prim(P22_pv,ix,iy) = gamma
          ENDIF
      
          
!         IF (twoD == 1) THEN 
!               IF ( penteX == 1 ) THEN 
!                  Prim(U_pv,ix, iy) = U_0
!                  Prim(V_pv,ix, iy) = 0.d0
!               ELSE IF (penteY == 1) THEN
!                 Prim(U_pv,ix, iy) = 0.d0
!                 Prim(V_pv,ix, iy) =  U_0 
!               ENDIF
!               Prim(P11_pv,ix,iy) = 0.5d0*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)*phi2 + eps
!               Prim(P12_pv,ix,iy) = 0.d0
!               Prim(P22_pv,ix,iy) = 0.5d0*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)*phi2   + eps

!           ELSE IF (oneDX == 1) THEN 
!               Prim(U_pv,ix, iy) = U_0
!               Prim(V_pv,ix, iy) =  0.d0
!               Prim(P11_pv,ix,iy) = Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)*phi2  
!               Prim(P12_pv,ix,iy) =  0.d0
!               Prim(P22_pv,ix,iy) = EPS
!           ELSE IF (oneDY == 1) THEN 
!              Prim(U_pv,ix, iy) = 0.d0
!              Prim(V_pv,ix, iy) =  U_0 
!              Prim(P11_pv,ix,iy) = eps
!              Prim(P12_pv,ix,iy) =  0.d0
!              Prim(P22_pv,ix,iy) = Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)*phi2
!           ENDIF

          SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
          SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
          SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
          SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
          Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 + Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
          Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0

   ENDDO
  ENDDO

! IF (test_shear == 1) then ! faire pas parallele ici
! !           DO ix = 1, Nx
! !                 X(ix) = X0 + 0.5D0*DX + (ix-1)*DX
! !           DO iy = 1, Ny
! !                 Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY
       
! !                       Prim(H_pv, ix, iy) = 0.125d0*H_0
! !                       Prim(U_pv,ix, iy) = 0.D0
! !                       Prim(V_pv,ix, iy) = 0.d0

! !                     Prim(P11_pv,ix,iy) = 1.d-7
! !                     Prim(P12_pv,ix,iy) = 0.0d0
! !                     Prim(P22_pv,ix,iy) = 1.d-7

! !                     SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
! !                     SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
! !                     SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
! !                     SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
! !                     Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 + Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
! !                     Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0
! !           ENDDO
! !           ENDDO

! !           theta = theta*dacos(-1.d0)/180.d0

! !      DO ix = 1, Nx
! !           X(ix) = 0.5D0*DX + (ix-1)*DX  
! !       DO iy = 1, Ny
! !         Y(iy) = 0.5D0*DY + (iy-1)*DY
! !         x1 = (X(ix)-X_centre)*dcos(theta)+(Y(iy)-Y_centre)*dsin(theta)
! !         y1 = (Y(iy)-Y_centre)*dcos(theta)-(X(ix)-X_centre)*dsin(theta)
! !         !Selection
! !         IF((x1/a)**2.d0+(y1/b)**2.d0.LE. 1.d0) THEN
! !                ! IF(obstacle .Gt.0) THEN
! !                   PRIM(H_pv,ix,iy)   = H_0
! !                   Prim(U_pv,ix, iy) = 0.D0
! !                   Prim(V_pv,ix, iy) = 0.d0

! !                     Prim(P11_pv,ix,iy) = 1.d-7
! !                     Prim(P12_pv,ix,iy) = 0.0d0
! !                     Prim(P22_pv,ix,iy) = 1.d-7

! !                     SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
! !                     SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
! !                     SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
! !                     SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
! !                     Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 + Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
! !                     Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0
! !               ENDIF
! !       ENDDO
! !       ENDDO

!           DO ix = 1, Nx/3
!                 X(ix) = X0 + 0.5D0*DX + (ix-1)*DX  
!           DO iy = 1, Ny
!                 Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY
     
!                       Prim(H_pv, ix, iy) = H_0
!                       Prim(U_pv,ix, iy) = 0.0D0
!                       Prim(V_pv,ix, iy) = 0.0d0
    
!                     Prim(P11_pv,ix,iy) = 1.d-4
!                     Prim(P12_pv,ix,iy) = 0.0d0
!                     Prim(P22_pv,ix,iy) = 1.d-4
!                     SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
!                     SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
!                     SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
!                     SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
!                     Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 !+ Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
!                     Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0
!           ENDDO
!           ENDDO

!           DO ix = Nx/3, 2*Nx/3
!                 X(ix) = X0 + 0.5D0*DX + (ix-1)*DX  
!           DO iy = 1, Ny
!                 Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY
     
!                       Prim(H_pv, ix, iy) = H_0
!                       Prim(U_pv,ix, iy) = 0.0D0
!                       Prim(V_pv,ix, iy) = 0.2d0
    
!                     Prim(P11_pv,ix,iy) = 1.d-4
!                     Prim(P12_pv,ix,iy) = 0.0d0
!                     Prim(P22_pv,ix,iy) = 1.d-4
!                     SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
!                     SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
!                     SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
!                     SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
!                     Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 + Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
!                     Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0
!           ENDDO
!           ENDDO

!            DO ix = 2*Nx/3, Nx
!                 X(ix) = X0 + 0.5D0*DX + (ix-1)*DX  
!           DO iy = 1, Ny
!                 Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY
     
!                       Prim(H_pv, ix, iy) = H_0
!                       Prim(U_pv,ix, iy) = 0.0D0
!                       Prim(V_pv,ix, iy) = 0.0d0
    
!                     Prim(P11_pv,ix,iy) = 1.d-4
!                     Prim(P12_pv,ix,iy) = 0.0d0
!                     Prim(P22_pv,ix,iy) = 1.d-4
!                     SOUND_ax(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P11_pv, ix, iy) )
!                     SOUND_bx(ix,iy) = DSQRT(Prim(P11_pv, ix, iy) )
!                     SOUND_ay(ix,iy) = DSQRT( g*Prim(H_pv, ix, iy) + 3.d0*Prim(P22_pv, ix, iy) )
!                     SOUND_by(ix,iy) = DSQRT(Prim(P22_pv, ix, iy) )
!                     Pression(ix,iy) = g*Prim(H_pv, ix, iy)*Prim(H_pv, ix, iy)/2.d0 + Prim(P11_pv, ix, iy)*Prim(H_pv, ix, iy)
!                     Ein(ix,iy) = ( Prim(H_pv, ix, iy)*g + Prim(P11_pv, ix, iy) + Prim(P22_pv, ix, iy) )/2.d0
!           ENDDO
!           ENDDO


! ENDIF

  ! VARIABLE CONSERVATIVES 
  DO ix = 1, Nx
    DO iy = 1, Ny 
      CONS(1,ix,iy) = Prim(H_pv, ix, iy)
      CONS(2,ix,iy) = Prim(H_pv, ix, iy)*Prim(U_pv, ix, iy)
      CONS(3,ix,iy) = Prim(H_pv, ix, iy)*Prim(V_pv, ix, iy)
      CONS(4,ix,iy) = Prim(H_pv, ix, iy)*Prim(P11_pv, ix, iy)
      CONS(5,ix,iy) = Prim(P12_pv, ix, iy)
      CONS(6,ix,iy) = Prim(H_pv, ix, iy)*Prim(P22_pv, ix, iy)
      CONS(7,ix,iy) = Prim(H_pv, ix, iy)*(Ein(ix,iy)+ (Prim(U_pv, ix, iy)**2.d0 + Prim(V_pv, ix, iy)**2.d0 )/2.d0)
    END DO
  END DO
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

  SUBROUTINE NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx, SOUND_ay, SOUND_by,CONS, Pression, it)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
   INTEGER :: ix,iy, it
  REAL (KIND = DP) :: p11, p12, p22
   REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), CONS(7,1:Nx,1:Ny), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
   REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1), SOUND_bx(0:Nx+1,0:Ny+1),SOUND_ay(0:Nx+1,0:Ny+1), SOUND_by(0:Nx+1,0:Ny+1)

  ! CALCUL NOUVELS VARIABLES PRIMITIVES
   DO ix = 1, Nx
     DO iy = 1, Ny

          Prim(H_pv,ix, iy) = CONS(1,ix, iy)
          Prim(U_pv,ix, iy) = CONS(2,ix, iy)/CONS(1,ix, iy)
          Prim(V_pv,ix, iy) = CONS(3,ix, iy)/CONS(1,ix, iy)
!------------------------------------------------------------------------------------------------------------------------
          if (dabs(Prim(V_pv,ix, iy)).le.1.d-8) then
           Prim(V_pv,ix, iy)=0.d0; cons(3, ix, iy) = 0.d0
          endif

          if (dabs(Prim(U_pv,ix, iy)).le.1.d-8) then 
           Prim(U_pv,ix, iy)=0.d0; cons(2, ix, iy) = 0.d0
          endif
!------------------------------------------------------------------------------------------------------------------------
          Prim(P11_pv,ix, iy) = CONS(4,ix, iy)/CONS(1,ix, iy)
          Prim(P11_pv,ix, iy) = dmax1(Prim(P11_pv,ix, iy), 1.d-8)

          Prim(P12_pv,ix, iy) = CONS(5,ix, iy)
          Prim(P12_pv,ix, iy) = dmin1(Prim(P12_pv,ix, iy), dsqrt(Prim(P11_pv,ix, iy)*Prim(P22_pv,ix, iy)))
          Prim(P12_pv,ix, iy) = dmax1(Prim(P12_pv,ix, iy), -dsqrt(Prim(P11_pv,ix, iy)*Prim(P22_pv,ix, iy) ))
          if (dabs(Prim(p12_pv,ix, iy)).le.1.d-8) then 
           Prim(p12_pv,ix, iy)=0.d0; cons(5, ix, iy) = 0.d0
          endif

          Prim(P22_pv,ix, iy) = CONS(6,ix, iy)/CONS(1,ix, iy)
          Prim(P22_pv,ix, iy) = dmax1(Prim(P22_pv,ix, iy), 1.d-8)
!------------------------------------------------------------------------------------------------------------------------
          p11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)

!------------------------------------------------------------------------------------------------------------------------
          Ein(ix, iy) = CONS(7,ix, iy)/CONS(1,ix, iy) - ((Prim(U_pv,ix, iy))**2.d0+(Prim(V_pv,ix, iy))**2.d0 )/2.d0
          SOUND_ax(ix, iy) = DSQRT( g*Prim(H_pv,ix, iy) + 3.d0*p11)
          SOUND_bx(ix, iy) = DSQRT(p11)  
          SOUND_ay(ix, iy) = DSQRT( g*Prim(H_pv,ix, iy) + 3.d0*p22)
          SOUND_by(ix, iy) = DSQRT(p22)  
          Pression(ix, iy) = g*(Prim(H_pv,ix, iy))**2.d0/2.d0 + p11*Prim(H_pv,ix, iy)
  END DO  
  END DO
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE Condition_lim_box(Prim, Ein, Pression, SOUND_ax,  SOUND_bx , SOUND_ay,  SOUND_by)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER ::  ix, iy, iv
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1),Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1) , SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)
  REAL(KIND=DP),ALLOCATABLE ::  tampon_env_h( :),tampon_env_u( :),tampon_env_v( :),tampon_env_p11( :),tampon_env_p12( :),tampon_env_p22( :)
  REAL(KIND=DP),ALLOCATABLE ::  tampon_rec_h( :),tampon_rec_u( :),tampon_rec_v( :),tampon_rec_p11( :),tampon_rec_p12( :),tampon_rec_p22( :)
   
  ALLOCATE(tampon_env_h(1:Ny), tampon_env_u(1:Ny),tampon_env_v(1:Ny),tampon_env_p11(1:Ny),tampon_env_p12(1:Ny),tampon_env_p22(1:Ny))
  ALLOCATE(tampon_rec_h(1:Ny), tampon_rec_u(1:Ny),tampon_rec_v(1:Ny),tampon_rec_p11(1:Ny),tampon_rec_p12(1:Ny),tampon_rec_p22(1:Ny))

  tag = 0;  tampon_env_h = 0; tampon_env_u= 0; tampon_env_v = 0; tampon_env_p11 = 0; tampon_env_p12 = 0; tampon_env_p22 = 0; 
  tampon_rec_h = 0; tampon_rec_u= 0; tampon_rec_v = 0; tampon_rec_p11 = 0; tampon_rec_p12 = 0; tampon_rec_p22 = 0; 
  nombre = Ny

  IF (penteX == 1) then 
    IF (NCPU == 1) THEN 
        !------------------------------------------------------------------------------------------------------------------------
             do iy = 1, Ny
             Prim(:,0, iy)    = Prim(:,Nx, iy)

             Ein(0, iy)       = Ein(Nx, iy)
             SOUND_ax(0, iy)  = SOUND_ax(Nx, iy)
             SOUND_bx(0, iy)  = SOUND_bx(Nx, iy)
             SOUND_ay(0, iy)  = SOUND_ay(Nx, iy)
             SOUND_by(0, iy)  = SOUND_by(Nx, iy)
             Pression(0, iy)  = Pression(Nx, iy)
            enddo
        !------------------------------------------------------------------------------------------------------------------------
            do iy = 1, Ny
             Prim(:,Nx+1, iy)   = Prim(:,1, iy)

             Ein(Nx+1, iy)       = Ein(1, iy)
             SOUND_ax(Nx+1, iy)  = SOUND_ax(1, iy)
             SOUND_bx(Nx+1, iy)  = SOUND_bx(1, iy)
             SOUND_ay(Nx+1, iy)  = SOUND_ay(1, iy)
             SOUND_by(Nx+1, iy)  = SOUND_by(1, iy)
             Pression(Nx+1, iy)  = Pression(1, iy)
            enddo   
 ELSE
          IF (rang == 0) THEN
              
                tampon_env_h(1:Ny) = Prim(H_pv, 1, 1:Ny) 
                tampon_env_u(1:Ny) = Prim(u_pv, 1, 1:Ny)
                tampon_env_v(1:Ny) = Prim(v_pv, 1, 1:Ny)
                tampon_env_p11(1:Ny) = Prim(p11_pv, 1, 1:Ny)
                tampon_env_p12(1:Ny) = Prim(p12_pv, 1, 1:Ny)
                tampon_env_p22(1:Ny) = Prim(p22_pv, 1, 1:Ny) 

               !Envoi vers proc de rang inferieur
                dest   = ncpu-1
             
                CALL MPI_SEND(tampon_env_h,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_u,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_v,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_p11,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_p12,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_p22,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                

                 !Reception du proc de rang ncpu -1
                source = ncpu-1
                
                CALL MPI_RECV(tampon_rec_h,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_u,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_v,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_p11,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_p12,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_p22,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
             
                 Prim(H_pv, 0, 1:Ny) = tampon_rec_h( 1:Ny) 
                 Prim(u_pv, 0, 1:Ny) = tampon_rec_u( 1:Ny) 
                 Prim(v_pv, 0, 1:Ny) = tampon_rec_v( 1:Ny) 
                 Prim(p11_pv, 0, 1:Ny) = tampon_rec_p11( 1:Ny) 
                 Prim(p12_pv, 0, 1:Ny) = tampon_rec_p12( 1:Ny) 
                 Prim(p22_pv, 0, 1:Ny) = tampon_rec_p22( 1:Ny) 

                  SOUND_ax(0, 1:Ny)  = DSQRT( g*Prim(H_pv,0, 1:Ny) + 3.d0*Prim(P11_pv,0, 1:Ny))
                  SOUND_bx(0, 1:Ny)  = DSQRT(Prim(P11_pv,0, 1:Ny))  
                  SOUND_ay(0, 1:Ny)  = DSQRT( g*Prim(H_pv,0, 1:Ny) + 3.d0*Prim(P22_pv,0, 1:Ny))
                  SOUND_by(0, 1:Ny)  = DSQRT(Prim(P22_pv,0, 1:Ny))  
                  Pression(0, 1:Ny)  = g*(Prim(H_pv,0, 1:Ny))**2.d0/2.d0 + Prim(P11_pv,0, 1:Ny)*Prim(H_pv,0, 1:Ny)
            ENDIF

          IF (rang == ncpu-1) THEN
              
                tampon_env_h(1:Ny) = Prim(H_pv, Nx, 1:Ny) 
                tampon_env_u(1:Ny) = Prim(u_pv, Nx, 1:Ny)
                tampon_env_v(1:Ny) = Prim(v_pv, Nx, 1:Ny)
                tampon_env_p11(1:Ny) = Prim(p11_pv, Nx, 1:Ny)
                tampon_env_p12(1:Ny) = Prim(p12_pv, Nx, 1:Ny)
                tampon_env_p22(1:Ny) = Prim(p22_pv, Nx, 1:Ny) 
              
                !Envoi vers proc de rang 0
                dest   = 0

                CALL MPI_SEND(tampon_env_h,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_u,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_v,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_p11,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_p12,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
                CALL MPI_SEND(tampon_env_p22,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
               

                !Reception du proc de rang 0
                source = 0
                
                CALL MPI_RECV(tampon_rec_h,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_u,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_v,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_p11,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_p12,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
                CALL MPI_RECV(tampon_rec_p22,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
               

                 Prim(H_pv, Nx+1, 1:Ny) = tampon_rec_h( 1:Ny) 
                 Prim(u_pv, Nx+1, 1:Ny) = tampon_rec_u( 1:Ny) 
                 Prim(v_pv, Nx+1, 1:Ny) = tampon_rec_v( 1:Ny) 
                 Prim(p11_pv, Nx+1, 1:Ny) = tampon_rec_p11( 1:Ny) 
                 Prim(p12_pv, Nx+1, 1:Ny) = tampon_rec_p12( 1:Ny) 
                 Prim(p22_pv, Nx+1, 1:Ny) = tampon_rec_p22( 1:Ny) 

                  SOUND_ax( Nx+1, 1:Ny)  = DSQRT( g*Prim(H_pv, Nx+1, 1:Ny) + 3.d0*Prim(P11_pv, Nx+1, 1:Ny))
                  SOUND_bx( Nx+1, 1:Ny)  = DSQRT(Prim(P11_pv, Nx+1, 1:Ny))  
                  SOUND_ay( Nx+1, 1:Ny)  = DSQRT( g*Prim(H_pv, Nx+1, 1:Ny) + 3.d0*Prim(P22_pv, Nx+1, 1:Ny))
                  SOUND_by( Nx+1, 1:Ny)  = DSQRT(Prim(P22_pv, Nx+1, 1:Ny))  
                  Pression( Nx+1, 1:Ny)  = g*(Prim(H_pv, Nx+1, 1:Ny))**2.d0/2.d0 + Prim(P11_pv, Nx+1, 1:Ny)*Prim(H_pv, Nx+1, 1:Ny)

            ENDIF
  
ENDIF
        !------------------------------------------------------------------------------------------------------------------------
            do ix = 1, Nx
             Prim(:,ix, 0)       = Prim(:,ix, 1)
             Prim(V_pv,ix, 0)    = - Prim(V_pv,ix, 1)
             Prim(p12_pv,ix, 0)  = - Prim(p12_pv,ix, 1)

            Ein(ix, 0)          = (Prim(H_pv, ix, 0)*g + Prim(P11_pv, ix, 0) + Prim(P22_pv, ix, 0) )/2.d0
            SOUND_ax(ix, 0)     = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P11_pv,ix, 0)) 
            SOUND_bx(ix, 0)     = DSQRT(Prim(P11_pv,ix, 0))
            SOUND_ay(ix, 0)     = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P22_pv,ix, 0)) 
            SOUND_by(ix, 0)     = DSQRT(Prim(P22_pv,ix, 0))  
            Pression(ix, 0)     = g*(Prim(H_pv,ix, 0))**2.d0/2.d0 + Prim(P11_pv,ix, 0)*Prim(H_pv,ix, 0)


!              Prim(:,ix, 0)    = Prim(:,ix, Ny)
           
!              Ein(ix, 0)       = Ein(ix, Ny)
!              SOUND_ax(ix, 0)  = SOUND_ax(ix, Ny)
!              SOUND_bx(ix, 0)  = SOUND_bx(ix, Ny)
!              SOUND_ay(ix, 0)  = SOUND_ay(ix, Ny)
!              SOUND_by(ix, 0)  = SOUND_by(ix, Ny)
!              Pression(ix, 0)  = Pression(ix, Ny)
            enddo
        !------------------------------------------------------------------------------------------------------------------------
            do ix = 1, Nx
             Prim(:,ix, Ny+1)   = Prim(:,ix, Ny)
             Prim(V_pv,ix, Ny+1)    = - Prim(V_pv,ix, Ny)
             Prim(p12_pv,ix, Ny+1)    = - Prim(p12_pv,ix, Ny)

             Ein(ix, Ny +1)         = (Prim(H_pv, ix, Ny +1)*g + Prim(P11_pv, ix, Ny +1) + Prim(P22_pv, ix, Ny +1) )/2.d0
             SOUND_ax(ix, Ny +1)    = DSQRT( g*Prim(H_pv,ix, Ny +1) + 3.d0*Prim(P11_pv,ix, Ny +1))
             SOUND_bx(ix, Ny +1)    = DSQRT(Prim(P11_pv,ix,Ny +1))
             SOUND_ay(ix, Ny +1)    = DSQRT( g*Prim(H_pv,ix, Ny +1) + 3.d0*Prim(P22_pv,ix, Ny +1))
             SOUND_by(ix, Ny +1)    = DSQRT(Prim(P22_pv,ix, Ny +1)) 
             Pression(ix, Ny +1)    = g*(Prim(H_pv,ix, Ny +1))**2.d0/2.d0 + Prim(P11_pv,ix,Ny +1)*Prim(H_pv,ix,Ny +1)

!              Prim(:,ix, Ny+1)   = Prim(:,ix, 1)
    
!              Ein(ix, Ny+1)       = Ein(ix, 1)
!              SOUND_ax(ix, Ny+1)  = SOUND_ax(ix, 1)
!              SOUND_bx(ix, Ny+1)  = SOUND_bx(ix, 1)
!              SOUND_ay(ix, Ny+1)  = SOUND_ay(ix, 1)
!              SOUND_by(ix, Ny+1)  = SOUND_by(ix, 1)
!              Pression(ix, Ny+1)  = Pression(ix, 1)
            enddo

ELSE IF (penteY == 1) then 
            do iy = 1, Ny

                 Prim(:,0, iy)    = Prim(:,1, iy)
                 Prim(u_pv,0, iy)    = - Prim(u_pv,1, iy)
                 Prim(p12_pv,0, iy)    = - Prim(p12_pv,1, iy)

                Ein(0, iy)          = (Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
                SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
                SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
                SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
                SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
                Pression(0, iy)     = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
                enddo
            !------------------------------------------------------------------------------------------------------------------------
                do iy = 1, Ny

                 Prim(:,nx+1, iy)   = Prim(:,Nx, iy)
                 Prim(u_pv,nx+1, iy)    = - Prim(u_pv,Nx, iy)
                 Prim(p12_pv,nx+1, iy)    = - Prim(p12_pv,Nx, iy)

                  Ein(nx+1, iy)          = (Prim(H_pv, nx+1, iy)*g + Prim(P11_pv, nx+1, iy) + Prim(P22_pv, nx+1, iy) )/2.d0
                  SOUND_ax(nx+1, iy)     = DSQRT( g*Prim(H_pv,nx+1, iy) + 3.d0*Prim(P11_pv,nx+1, iy)) 
                  SOUND_bx(nx+1, iy)     = DSQRT(Prim(P11_pv,nx+1, iy))
                  SOUND_ay(nx+1, iy)     = DSQRT( g*Prim(H_pv,nx+1, iy) + 3.d0*Prim(P22_pv,nx+1, iy)) 
                  SOUND_by(nx+1, iy)     = DSQRT(Prim(P22_pv,nx+1, iy))  
                  Pression(nx+1, iy)     = g*(Prim(H_pv,nx+1, iy))**2.d0/2.d0 + Prim(P11_pv,nx+1, iy)*Prim(H_pv,nx+1, iy)
                enddo   
            !------------------------------------------------------------------------------------------------------------------------
                do ix = 1, Nx
                Prim(:,ix, 0)    = Prim(:,ix, Ny)

                 Ein(ix, 0)       = Ein(ix, Ny)
                 SOUND_ax(ix, 0)  = SOUND_ax(ix, Ny)
                 SOUND_bx(ix, 0)  = SOUND_bx(ix, Ny)
                 SOUND_ay(ix, 0)  = SOUND_ay(ix, Ny)
                 SOUND_by(ix, 0)  = SOUND_by(ix, Ny)
                 Pression(ix, 0)  = Pression(ix, Ny)

                enddo
            !------------------------------------------------------------------------------------------------------------------------
                do ix = 1, Nx
                    Prim(:,ix, Ny+1)   = Prim(:,ix, 1)

                 Ein(ix, Ny+1)      = Ein(ix, 1)
                 SOUND_ax(ix, Ny+1)  = SOUND_ax(ix, 1)
                 SOUND_bx(ix, Ny+1)  = SOUND_bx(ix, 1)
                 SOUND_ay(ix, Ny+1)  = SOUND_ay(ix, 1)
                 SOUND_by(ix, Ny+1)  = SOUND_by(ix, 1)
                 Pression(ix, Ny+1)  = Pression(ix, 1)
                enddo

  ENDIF 
  return
  END SUBROUTINE
  !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
  
  SUBROUTINE CondLimABSORPTION(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by)
  USE precisions
  use GlobalParam
    IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) ::  SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)  

!------------------------------------------------------------------------------------------------------------------------

IF (RANG == 0) THEN 
    do iy = 1, Ny
     Prim(:,0,  iy)  = Prim(:,1, iy)

     Ein(0,     iy)  = Ein(1, iy)
     SOUND_ax(0, iy) = SOUND_ax(1,  iy)
     SOUND_bx(0, iy) = SOUND_bx(1,  iy)
     SOUND_ay(0, iy) = SOUND_ay(1,  iy)
     SOUND_by(0, iy) = SOUND_by(1,  iy)
     Pression(0,iy)  = Pression(1, iy)
    enddo
ENDIF
!------------------------------------------------------------------------------------------------------------------------

IF (RANG == NCPU -1) THEN 
    do iy = 1, Ny
     Prim(:,Nx+1,  iy) = Prim(:,Nx, iy)

     Ein(Nx+1,     iy)  = Ein(Nx, iy)
     SOUND_ax(Nx+1, iy) = SOUND_ax(Nx, iy)
     SOUND_bx(Nx+1, iy) = SOUND_bx(Nx, iy)
     SOUND_ay(Nx+1, iy) = SOUND_ay(Nx, iy)
     SOUND_by(Nx+1, iy) = SOUND_by(Nx, iy)
     Pression(Nx+1,iy)  = Pression(Nx,iy)
    enddo   
ENDIF
!------------------------------------------------------------------------------------------------------------------------

IF (ncpu == 1) THEN 
    do iy = 1, Ny
     Prim(:,0,  iy)  = Prim(:,1, iy)

     Ein(0,     iy)  = Ein(1, iy)
     SOUND_ax(0, iy) = SOUND_ax(1,  iy)
     SOUND_bx(0, iy) = SOUND_bx(1,  iy)
     SOUND_ay(0, iy) = SOUND_ay(1,  iy)
     SOUND_by(0, iy) = SOUND_by(1,  iy)
     Pression(0,iy)  = Pression(1, iy)

     Prim(:,Nx+1,  iy) = Prim(:,Nx, iy)

     Ein(Nx+1,     iy)  = Ein(Nx, iy)
     SOUND_ax(Nx+1, iy) = SOUND_ax(Nx, iy)
     SOUND_bx(Nx+1, iy) = SOUND_bx(Nx, iy)
     SOUND_ay(Nx+1, iy) = SOUND_ay(Nx, iy)
     SOUND_by(Nx+1, iy) = SOUND_by(Nx, iy)
     Pression(Nx+1,iy)  = Pression(Nx,iy)
    enddo   
ENDIF



!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix,  0)  = Prim(:,ix, 1)

     Ein(ix,     0)  = Ein(ix, 1)
     SOUND_ax(ix, 0) = SOUND_ax(ix,  1)
     SOUND_bx(ix, 0) = SOUND_bx(ix,  1)
     SOUND_ay(ix, 0) = SOUND_ay(ix,  1)
     SOUND_by(ix, 0) = SOUND_by(ix,  1)
     Pression(ix,0)  = Pression(ix, 1)
    enddo
!------------------------------------------------------------------------------------------------------------------------
    do ix = 1, Nx
     Prim(:,ix,  Ny+1) = Prim(:,ix, Ny)

     Ein(ix,     Ny+1) = Ein(ix,     Ny)
     SOUND_ax(ix, Ny+1) = SOUND_ax(ix, Ny)
     SOUND_bx(ix, Ny+1) = SOUND_bx(ix, Ny)
     SOUND_ay(ix, Ny+1) = SOUND_ay(ix, Ny)
     SOUND_by(ix, Ny+1) = SOUND_by(ix, Ny)
     Pression(ix,Ny+1)  = Pression(ix,Ny)
    enddo
  return
  END SUBROUTINE
  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE Condition_lim_batteur(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,TIME, dy)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: q_0,u_0, yi,xi, omega, TIME, dy
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)  
    pi = 4.0d0*ATAN(1.0d0)
    omega = 6.19012d0
!------------------------------------------------------------------------------------------------------------------------
   if (angle .le. 1.d-8 .and. frottcoeff .le. 1.d-8 .and. disscoeff .le. 1.d-8) then 
      u_0 = 3.d0
   else 
      u_0 = DSQRT(H_0*g*Dtan(angle)/frottcoeff)
   end if 

   q_0 = H_0*u_0  
!------------------------------------------------------------------------------------------------------------------------


IF (penteX == 1) then 
  IF (RANG == 0) THEN  
          do iy = 1, Ny
            Prim(H_pv,0, iy) = H_0*(1.d0 + amplitude*DSIN(omega*TIME))

            IF(Ny> 1 ) THEN
!               yi =  REAL(iy-0.5)/REAL(Ny) 
!               Prim(H_pv,0, iy) = H_0*(1.d0 + amplitude*DSIN(omega*TIME)+ amplitude*SIN(2.0*pi*yi)*SIN(omega*Time/2.0))

              Y(iy) =0.5D0*DY + (iy-1)*DY
              !yi =  REAL(iy-0.5)/REAL(Ny) 
              Prim(H_pv,0, iy) = H_0*(1.d0 + amplitude*DCOS(4.0*pi*y(iy)/Ly)*DSIN(omega*Time/2.0)*DSIN(omega*Time))
            ENDIF
       
             Prim(U_pv,0, iy)    = q_0/Prim(H_pv,0, iy) 
             Prim(V_pv,0, iy)    = 0.d0 
         
           IF (twoD == 1) THEN 
             Prim(P11_pv,0, iy)  = 0.5D0*Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2
             Prim(P12_pv,0, iy)  = 0.d0
             Prim(P22_pv,0, iy)  = 0.5D0*Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2
           ELSE IF (oneDX == 1) THEN 
             Prim(P11_pv,0, iy)  = Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2
             Prim(P12_pv,0, iy)  = 0.d0
             Prim(P22_pv,0, iy)  = EPS 
           ENDIF 

           Ein(0, iy)          = (Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
           SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
           SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
           SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
           SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
           Pression(0, iy)     = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
          enddo
  ENDIF
        !-------------------------------------------------------------------------------------------------------------------
  IF (RANG == NCPU -1) THEN       
          do iy = 1, Ny
           Prim(H_pv,Nx+1, iy)   = Prim(H_pv,Nx, iy)
           Prim(U_pv,Nx+1, iy)   = Prim(U_pv,Nx, iy)
           Prim(V_pv,Nx+1, iy)   = Prim(V_pv,Nx, iy) 
           Prim(P11_pv,Nx+1, iy) = Prim(P11_pv,Nx, iy)
           Prim(P12_pv,Nx+1, iy) = Prim(P12_pv,Nx, iy)
           Prim(P22_pv,Nx+1, iy) = Prim(P22_pv,Nx, iy)
           Ein(Nx+1, iy)         = Ein(Nx, iy)
           SOUND_ax(Nx+1, iy)    = SOUND_ax(Nx, iy)
           SOUND_bx(Nx+1, iy)    = SOUND_bx(Nx, iy)
           SOUND_ay(Nx+1, iy)    = SOUND_ay(Nx, iy)
           SOUND_by(Nx+1, iy)    = SOUND_by(Nx, iy)
           Pression(Nx+1, iy)    = Pression(Nx, iy)
          enddo
  ENDIF
        !--------------------------------------------------------------------------------------------------------------------
        IF (WALLS == 1) THEN 
              do ix = 1, Nx
             
                     DO iv = 1 , Nv_Prim
                      Prim(iv,ix, 0) = Prim(iv,ix, 1)
                     ENDDO
                     Prim(V_pv, ix, 0)   = - Prim(V_pv, ix,1)
                     Prim(p12_pv, ix, 0) = - Prim(p12_pv, ix, 1)
                   
                      Ein(ix, 0)          = (Prim(H_pv, ix, 0)*g + Prim(P11_pv, ix, 0) + Prim(P22_pv, ix, 0) )/2.d0
                      SOUND_ax(ix, 0)     = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P11_pv,ix, 0)) 
                      SOUND_bx(ix, 0)     = DSQRT(Prim(P11_pv,ix, 0))
                      SOUND_ay(ix, 0)     = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P22_pv,ix, 0)) 
                      SOUND_by(ix, 0)     = DSQRT(Prim(P22_pv,ix, 0))  
                      Pression(ix, 0)     = g*(Prim(H_pv,ix, 0))**2.d0/2.d0 + Prim(P11_pv,ix, 0)*Prim(H_pv,ix, 0)
          
                
                     DO iv = 1 , Nv_Prim
                       Prim(iv,ix, Ny +1)   = Prim(iv,ix, Ny)
                     ENDDO

                     Prim(V_pv,ix, Ny +1)   = -Prim(V_pv,ix, Ny) 
                     Prim(P12_pv,ix, Ny +1) = -Prim(P12_pv,ix, Ny)
                     
                     Ein(ix, Ny +1)         = (Prim(H_pv, ix, Ny +1)*g + Prim(P11_pv, ix, Ny +1) + Prim(P22_pv, ix, Ny +1) )/2.d0
                     SOUND_ax(ix, Ny +1)    = DSQRT( g*Prim(H_pv,ix, Ny +1) + 3.d0*Prim(P11_pv,ix, Ny +1))
                     SOUND_bx(ix, Ny +1)    = DSQRT(Prim(P11_pv,ix,Ny +1))
                     SOUND_ay(ix, Ny +1)    = DSQRT( g*Prim(H_pv,ix, Ny +1) + 3.d0*Prim(P22_pv,ix, Ny +1))
                     SOUND_by(ix, Ny +1)    = DSQRT(Prim(P22_pv,ix, Ny +1)) 
                     Pression(ix, Ny +1)    = g*(Prim(H_pv,ix, Ny +1))**2.d0/2.d0 + Prim(P11_pv,ix,Ny +1)*Prim(H_pv,ix,Ny +1)
           
              end do
        ELSE
            do ix = 1, Nx
               DO iv = 1 , Nv_Prim
                Prim(iv,ix, 0)   = Prim(iv,ix, Ny)
               ENDDO

               Ein(ix, 0)          = Ein(ix, Ny)
               SOUND_ax(ix, 0)     = SOUND_ax(ix, Ny)
               SOUND_bx(ix, 0)     = SOUND_bx(ix, Ny)
               SOUND_ay(ix, 0)     = SOUND_ay(ix, Ny)
               SOUND_by(ix, 0)     = SOUND_by(ix, Ny)
               Pression(ix, 0)     = Pression(ix, Ny)

               DO iv = 1 , Nv_Prim
                Prim(iv,ix, Ny +1)   = Prim(iv,ix, 1)
               ENDDO
             
               Ein(ix, Ny +1)         = Ein(ix, 1)
               SOUND_ax(ix, Ny +1)    = SOUND_ax(ix, 1)
               SOUND_bx(ix, Ny +1)    = SOUND_bx(ix, 1)
               SOUND_ay(ix, Ny +1)   = SOUND_ay(ix, 1)
               SOUND_by(ix, Ny +1)    = SOUND_by(ix, 1)
               Pression(ix, Ny +1)    = Pression(ix, 1)
            end do
        ENDIF
       
   ELSE IF (penteY == 1) then 
     do ix = 1, Nx
             Prim(H_pv,ix, 0) = H_0*(1.d0 + amplitude*DSIN(omega*TIME))

            IF(Nx> 1 ) THEN
              xi =  REAL(ix-0.5)/REAL(Nx) 
              Prim(H_pv,ix, 0) = H_0*(1.d0 + amplitude*DSIN(omega*TIME)+ amplitude*SIN(2.0*pi*xi)*SIN(omega*Time/2.0))
            ENDIF
         
           Prim(U_pv,ix, 0)    = 0.d0 
           Prim(V_pv,ix, 0)    = q_0/Prim(H_pv,ix, 0)
          
           IF (twoD == 1) THEN 
             Prim(P11_pv,ix, 0)  = 0.5D0*Prim(H_pv,ix, 0)*Prim(H_pv,ix, 0)*phi2
             Prim(P12_pv,ix, 0)  = 0.d0
             Prim(P22_pv,ix, 0)  = 0.5D0*Prim(H_pv,ix, 0)*Prim(H_pv,ix, 0)*phi2  
           ENDIF
           IF (oneDY == 1) THEN
             Prim(P11_pv,ix, 0)  = EPS 
             Prim(P12_pv,ix, 0)  = 0.d0
             Prim(P22_pv,ix, 0)  =  Prim(H_pv,ix, 0)*Prim(H_pv,ix, 0)*phi2 
           ENDIF

           Ein(ix, 0)          = (Prim(H_pv, ix, 0)*g + Prim(P11_pv, ix, 0) + Prim(P22_pv, ix, 0) )/2.d0
           SOUND_ax(ix, 0)     = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P11_pv,ix, 0)) 
           SOUND_bx(ix, 0)     = DSQRT(Prim(P11_pv,ix, 0))
           SOUND_ay(ix, 0)     = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P22_pv,ix, 0)) 
           SOUND_by(ix, 0)     = DSQRT(Prim(P22_pv,ix, 0))  
           Pression(ix, 0)     = g*(Prim(H_pv,ix, 0))**2.d0/2.d0 + Prim(P11_pv,ix, 0)*Prim(H_pv,ix, 0)
          enddo
        !------------------------------------------------------------------------------------------------------------------------
          IF (WALLS == 1) THEN 
            do iy = 1, Ny
          
                   DO iv = 1 , Nv_Prim
                    Prim(iv,0, iy)   = Prim(iv,1, iy)
                   ENDDO
                   Prim(u_pv, 0, iy) = - Prim(u_pv, 1, iy)
                   Prim(p12_pv,0, iy) = - Prim(p12_pv, 1, iy)
                 
                      Ein(0, iy)          = (Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
                      SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
                      SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
                      SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
                      SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
                      Pression(0, iy)     = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
              
              

                
                   DO iv = 1 , Nv_Prim
                     Prim(iv,Nx+1, iy)   = Prim(iv,Nx, iy)
                   ENDDO
                   Prim(u_pv,Nx+1, iy)   = -Prim(u_pv,Nx, iy) 
                   Prim(P12_pv,Nx+1, iy) = -Prim(P12_pv,Nx, iy)
                   
                  Ein(nx+1, iy)          = (Prim(H_pv, nx+1, iy)*g + Prim(P11_pv, nx+1, iy) + Prim(P22_pv, nx+1, iy) )/2.d0
                  SOUND_ax(nx+1, iy)     = DSQRT( g*Prim(H_pv,nx+1, iy) + 3.d0*Prim(P11_pv,nx+1, iy)) 
                  SOUND_bx(nx+1, iy)     = DSQRT(Prim(P11_pv,nx+1, iy))
                  SOUND_ay(nx+1, iy)     = DSQRT( g*Prim(H_pv,nx+1, iy) + 3.d0*Prim(P22_pv,nx+1, iy)) 
                  SOUND_by(nx+1, iy)     = DSQRT(Prim(P22_pv,nx+1, iy))  
                  Pression(nx+1, iy)     = g*(Prim(H_pv,nx+1, iy))**2.d0/2.d0 + Prim(P11_pv,nx+1, iy)*Prim(H_pv,nx+1, iy)
            
            end do
          ELSE
             do iy = 1, Ny
             
                   DO iv = 1 , Nv_Prim
                    Prim(iv,0, iy)   = Prim(iv,Nx, iy)
                   ENDDO

                   Ein(0, iy)          = Ein(Nx, iy)
                   SOUND_ax(0, iy)     = SOUND_ax(Nx, iy)
                   SOUND_bx(0, iy)     = SOUND_bx(Nx, iy)
                   SOUND_ay(0, iy)     = SOUND_ay(Nx, iy)
                   SOUND_by(0, iy)     = SOUND_by(Nx, iy)
                   Pression(0, iy)     = Pression(Nx, iy)
             
              
       
                   DO iv = 1 , Nv_Prim
                    Prim(iv,Nx+1, iy)   = Prim(iv,1, iy)
                   ENDDO
                 
                   Ein(Nx+1, iy)         = Ein(1, iy)
                   SOUND_ax(Nx+1, iy)    = SOUND_ax(1, iy)
                   SOUND_bx(Nx+1, iy)    = SOUND_bx(1, iy)
                   SOUND_ay(Nx+1, iy)    = SOUND_ay(1, iy)
                   SOUND_by(Nx+1, iy)    = SOUND_by(1, iy)
                   Pression(Nx+1, iy)    = Pression(1, iy)
         
             end do
          ENDIF

   ENDIF
   return
   END SUBROUTINE
  !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

! SUBROUTINE Cond_lim_jump(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,TIME, dy)
!   USE precisions
!   use GlobalParam
!   IMPLICIT NONE
!   INTEGER ::  ix, iy
!   REAL (KIND = DP) :: q_0,u_0, yi,xi, omega, TIME, dy
!   REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
!   REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)  
!     pi = 4.0d0*ATAN(1.0d0)
!     omega = 6.19012d0
! !------------------------------------------------------------------------------------------------------------------------
!    if (angle .le. 1.d-8 .and. frottcoeff .le. 1.d-8 .and. disscoeff .le. 1.d-8) then 
!       u_0 = 3.d0
!    else 
!       u_0 = DSQRT(H_0*g*Dtan(angle)/frottcoeff)
!    end if 

!    q_0 = H_0*u_0  
! !------------------------------------------------------------------------------------------------------------------------

!   IF (RANG == 0) THEN  
!           do iy = 1, Ny
!             Prim(H_pv,0, iy) = H_0
!             Prim(U_pv,0, iy)    = q_0/Prim(H_pv,0, iy) 
!             Prim(V_pv,0, iy)    = 0.d0 
         
!            IF (twoD == 1) THEN 
!              Prim(P11_pv,0, iy)  = 0.5D0*Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2
!              Prim(P12_pv,0, iy)  = 0.d0
!              Prim(P22_pv,0, iy)  = 0.5D0*Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2
!            ELSE IF (oneDX == 1) THEN 
!              Prim(P11_pv,0, iy)  = Prim(H_pv,0, iy)*Prim(H_pv,0, iy)*phi2
!              Prim(P12_pv,0, iy)  = 0.d0
!              Prim(P22_pv,0, iy)  = EPS 
!            ENDIF 

!            Ein(0, iy)          = (Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
!            SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
!            SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
!            SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
!            SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
!            Pression(0, iy)     = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
!           enddo
!   ENDIF
!         !-------------------------------------------------------------------------------------------------------------------
!   IF (RANG == NCPU -1) THEN       
!           do iy = 1, Ny
!            Prim(H_pv,Nx+1, iy)   = Prim(H_pv,Nx, iy)
!            if (Prim(H_pv,Nx+1, iy) .le. dev) then 
!             Prim(U_pv,Nx+1, iy)   = - Prim(U_pv,Nx, iy)
!            else
!             Cd = pi/(pi + 2.d0) + 0.08d0*(PRIM(H_PV, Nx+1, iy) - dev)/dev
!             PRIM(U_PV, Nx+1, iy)=1.d0/(PRIM(H_PV, Nx+1,iy))*(2.d0/3.d0*Cd*dsqrt(2.d0*G*(PRIM(H_PV, Nx+1, iy) - dev)**3.d0))
!            endif
!            Prim(V_pv,Nx+1, iy)   = Prim(V_pv,Nx, iy) 
!            Prim(P11_pv,Nx+1, iy) = Prim(P11_pv,Nx, iy)
!            Prim(P12_pv,Nx+1, iy) = Prim(P12_pv,Nx, iy)
!            Prim(P22_pv,Nx+1, iy) = Prim(P22_pv,Nx, iy)
          
!           Ein(nx+1, iy)          = (Prim(H_pv, nx+1, iy)*g + Prim(P11_pv, nx+1, iy) + Prim(P22_pv, nx+1, iy) )/2.d0
!           SOUND_ax(nx+1, iy)     = DSQRT( g*Prim(H_pv,nx+1, iy) + 3.d0*Prim(P11_pv,nx+1, iy)) 
!           SOUND_bx(nx+1, iy)     = DSQRT(Prim(P11_pv,nx+1, iy))
!           SOUND_ay(nx+1, iy)     = DSQRT( g*Prim(H_pv,nx+1, iy) + 3.d0*Prim(P22_pv,nx+1, iy)) 
!           SOUND_by(nx+1, iy)     = DSQRT(Prim(P22_pv,nx+1, iy))  
!           Pression(nx+1, iy)     = g*(Prim(H_pv,nx+1, iy))**2.d0/2.d0 + Prim(P11_pv,nx+1, iy)*Prim(H_pv,nx+1, iy)
!           enddo
!   ENDIF
!         !--------------------------------------------------------------------------------------------------------------------
        
!    return
!    END SUBROUTINE
  !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE cond_lim_turb_sol(Prim, Ein, Pression, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,TIME, X,Y)
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  INTEGER ::  ix, iy
  REAL (KIND = DP) :: TIME, h, p11, p22
  REAL (KIND = DP) :: X(1:Nx), Y(1:Ny)
  REAL (KIND = DP) :: Prim(6,0:Nx+1,0:Ny+1), Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1) ,SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1) ,SOUND_by(0:Nx+1,0:Ny+1)  
 
!  IF (RANG == 0) THEN  
!   do iy = 1, Ny
!      Prim(H_pv,0, iy) = 1.d0/(1.d0 + beta**2.d0*time**2.d0)
!      Prim(U_pv,0, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(beta*time*(-0.5D0*DX) + y(iy)) 
!      Prim(V_pv,0, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(-0.5D0*DX + beta*time*y(iy))  
!      Prim(P11_pv,0, iy) = (lambda +gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
!      Prim(P12_pv,0, iy) = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
!      Prim(P22_pv,0, iy) = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0

!      Ein(0, iy)         =  ( Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
!      SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
!      SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
!      SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
!      SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
!      Pression(0, iy)    = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)
!   enddo
! ENDIF
! !------------------------------------------------------------------------------------------------------------------------
!   IF (RANG == NCPU -1) THEN   
!   do iy = 1, Ny
!    Prim(H_pv,Nx+1, iy)   = 1.d0/(1.d0 + beta**2.d0*time**2.d0)
!    Prim(U_pv,Nx+1, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(beta*time*(Lx + 0.5D0*DX) + y(iy))
!    Prim(V_pv,Nx+1, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(-(Lx + 0.5D0*DX) + beta*time*y(iy))  
!    Prim(P11_pv,Nx+1, iy) = (lambda +gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
!    Prim(P12_pv,Nx+1, iy) = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
!    Prim(P22_pv,Nx+1, iy) = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0


!    Ein(Nx+1, iy)         = ( Prim(H_pv, Nx+1, iy)*g + Prim(P11_pv, Nx+1, iy) + Prim(P22_pv, Nx+1, iy) )/2.d0
!    SOUND_ax(Nx+1, iy)    = DSQRT( g*Prim(H_pv,Nx+1, iy) + 3.d0*Prim(P11_pv,Nx+1, iy)) 
!    SOUND_bx(Nx+1, iy)    = DSQRT(Prim(P11_pv,Nx+1, iy))
!    SOUND_ay(Nx+1, iy)    = DSQRT( g*Prim(H_pv,Nx+1, iy) + 3.d0*Prim(P22_pv,Nx+1, iy)) 
!    SOUND_by(Nx+1, iy)    = DSQRT(Prim(P22_pv,Nx+1, iy)) 
!    Pression(Nx+1, iy)    = g*(Prim(H_pv,Nx+1, iy))**2.d0/2.d0 + Prim(P11_pv,Nx+1, iy)*Prim(H_pv,Nx+1, iy)
!   enddo
! ENDIF

IF (ncpu == 1) then 
   do iy = 1, Ny
      Y(iy) = Y0 + 0.5D0*DY + (iy-1)*DY

     Prim(H_pv,0, iy) = 1.d0/(1.d0 + beta**2.d0*time**2.d0)
     Prim(U_pv,0, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(beta*time*(-0.5D0*DX) + y(iy)) 
     Prim(V_pv,0, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(-0.5D0*DX + beta*time*y(iy))  
     Prim(P11_pv,0, iy) = (lambda +gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
     Prim(P12_pv,0, iy) = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
     Prim(P22_pv,0, iy) = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0

     Ein(0, iy)         =  ( Prim(H_pv, 0, iy)*g + Prim(P11_pv, 0, iy) + Prim(P22_pv, 0, iy) )/2.d0
     SOUND_ax(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P11_pv,0, iy)) 
     SOUND_bx(0, iy)     = DSQRT(Prim(P11_pv,0, iy))
     SOUND_ay(0, iy)     = DSQRT( g*Prim(H_pv,0, iy) + 3.d0*Prim(P22_pv,0, iy)) 
     SOUND_by(0, iy)     = DSQRT(Prim(P22_pv,0, iy))  
     Pression(0, iy)    = g*(Prim(H_pv,0, iy))**2.d0/2.d0 + Prim(P11_pv,0, iy)*Prim(H_pv,0, iy)



   Prim(H_pv,Nx+1, iy)   = 1.d0/(1.d0 + beta**2.d0*time**2.d0)
   Prim(U_pv,Nx+1, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(beta*time*(Lx + 0.5D0*DX) + y(iy))
   Prim(V_pv,Nx+1, iy)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(-(Lx + 0.5D0*DX) + beta*time*y(iy))  
   Prim(P11_pv,Nx+1, iy) = (lambda +gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
   Prim(P12_pv,Nx+1, iy) = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
   Prim(P22_pv,Nx+1, iy) = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0


   Ein(Nx+1, iy)         = ( Prim(H_pv, Nx+1, iy)*g + Prim(P11_pv, Nx+1, iy) + Prim(P22_pv, Nx+1, iy) )/2.d0
   SOUND_ax(Nx+1, iy)    = DSQRT( g*Prim(H_pv,Nx+1, iy) + 3.d0*Prim(P11_pv,Nx+1, iy)) 
   SOUND_bx(Nx+1, iy)    = DSQRT(Prim(P11_pv,Nx+1, iy))
   SOUND_ay(Nx+1, iy)    = DSQRT( g*Prim(H_pv,Nx+1, iy) + 3.d0*Prim(P22_pv,Nx+1, iy)) 
   SOUND_by(Nx+1, iy)    = DSQRT(Prim(P22_pv,Nx+1, iy)) 
   Pression(Nx+1, iy)    = g*(Prim(H_pv,Nx+1, iy))**2.d0/2.d0 + Prim(P11_pv,Nx+1, iy)*Prim(H_pv,Nx+1, iy)
   enddo
endif
!------------------------------------------------------------------------------------------------------------------------
  do ix = 1, Nx
      X(ix) = X0 + 0.5D0*DX + (ix-1)*DX  + nbx0*rang*dx  

   Prim(H_pv,ix, 0)   = 1.d0/(1.d0 + beta**2.d0*time**2.d0)
   Prim(U_pv,ix, 0)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(beta*time*x(ix) -0.5d0*dy)
   Prim(V_pv,ix, 0)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(-x(ix) + beta*time*(-0.5d0*dy))  
   Prim(P11_pv,ix, 0) = (lambda +gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
   Prim(P12_pv,ix, 0) = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
   Prim(P22_pv,ix, 0) = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0

   Ein(ix, 0)         = ( Prim(H_pv, ix, 0)*g + Prim(P11_pv, ix, 0) + Prim(P22_pv,ix, 0) )/2.d0
   SOUND_ax(ix, 0)    = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P11_pv,ix, 0)) 
   SOUND_bx(ix, 0)    = DSQRT(Prim(P11_pv,ix, 0))
   SOUND_ay(ix, 0)    = DSQRT( g*Prim(H_pv,ix, 0) + 3.d0*Prim(P22_pv,ix, 0)) 
   SOUND_by(ix, 0)    = DSQRT(Prim(P22_pv,ix, 0)) 
   Pression(ix, 0)    = g*(Prim(H_pv,ix, 0))**2.d0/2.d0 + Prim(P11_pv,ix, 0)*Prim(H_pv,ix, 0)


   Prim(H_pv,ix, Ny +1)   = 1.d0/(1.d0 + beta**2.d0*time**2.d0)
   Prim(U_pv,ix, Ny +1)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(beta*time*x(ix) + Ly + 0.5d0*dy)
   Prim(V_pv,ix, Ny +1)   = beta/(1.d0 + beta**2.d0*time**2.d0)*(-x(ix) + beta*time*(Ly + 0.5d0*dy))  
   Prim(P11_pv,ix, Ny +1) = (lambda +gamma*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0
   Prim(P12_pv,ix, Ny +1) = (lambda - gamma)*beta*time/(1.d0 + beta**2.d0*time**2.d0)**2.d0
   Prim(P22_pv,ix, Ny +1) = (gamma + lambda*beta**2.d0*time**2.d0)/(1.d0 + beta**2.d0*time**2.d0)**2.d0

   Ein(ix, Ny +1)         = ( Prim(H_pv, ix, Ny +1)*g + Prim(P11_pv, ix, Ny +1) + Prim(P22_pv,ix, Ny +1) )/2.d0
   SOUND_ax(ix, Ny +1)    = DSQRT( g*Prim(H_pv,ix, Ny +1) + 3.d0*Prim(P11_pv,ix, Ny +1)) 
   SOUND_bx(ix, Ny +1)    = DSQRT(Prim(P11_pv,ix, Ny +1))
   SOUND_ay(ix, Ny +1)    = DSQRT( g*Prim(H_pv,ix, Ny +1) + 3.d0*Prim(P22_pv,ix, Ny +1)) 
   SOUND_by(ix, Ny +1)    = DSQRT(Prim(P22_pv,ix, Ny +1)) 
   Pression(ix, Ny +1)    = g*(Prim(H_pv,ix, Ny +1))**2.d0/2.d0 + Prim(P11_pv,ix, Ny +1)*Prim(H_pv,ix, Ny +1)
   end do

  return
  END SUBROUTINE
  !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  SUBROUTINE euler_method(DT, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, iy, IT,k
  REAL (KIND = DP) :: DT, traceP, H, p11, p22, p12, u, v,Q, alpha
  REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1),SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1),SOUND_by(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1),fracp11,fracp22,erreur
  REAL (KIND = DP), ALLOCATABLE :: TS(:)
  ALLOCATE(TS(7))
  traceP = 0.d0; H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0; fracp11 = 0.d0; fracp22 = 0.d0; Q = 0.d0; alpha = 0.d0; TS(:) = 0.D0
!------------------------------------------------------------------------------------------------------------------------
  DO  ix = 1, Nx
    DO iy  = 1, Ny  
!------------------------------------------------------------------------------------------------------------------------ 
        if (cons(1,ix, iy).le.1.d-8) then
            print*, 'pas de l eau', cons(1,ix,iy), ix, iy, it
           stop
        end if
!------------------------------------------------------------------------------------------------------------------------
       !TERME SOURCE

            H = Prim(H_pv,ix, iy)
            P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
            traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
            fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
            alpha = disscoeff*(traceP - (H**2.d0)*phi2 + eps)/( traceP**2.d0) ! + eps*phi2*g*h_0
            alpha = dmax1(alpha, 0.d0)
            Q = alpha*traceP*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
     
            IF (penteY  == 1) THEN       
                TS(2) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                TS(3) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                TS(4) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                TS(5) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                TS(6) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                TS(7) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
           ELSE IF (penteX  == 1) THEN 
                TS(2) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                TS(3) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                TS(4) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                TS(5) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                TS(6) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                TS(7) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
            ENDIF

        !------------------------------------------------------------------------------------------------------------------------

            do k=1,7
              CONS(k,ix, iy) = CONS(k,ix, iy) + DT*TS(k) 
            end do 


          ENDDO
          ENDDO
        
 DEALLOCATE(TS)
  return
  END SUBROUTINE
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 SUBROUTINE rk2(DT, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, iy, IT,k
  REAL (KIND = DP) :: DT, dt2, traceP, H, p11, p22, p12, u, v,Q, alpha
  REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1),SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1),SOUND_by(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1), fracp11, fracp22, erreur
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: TS1, TS2, CONS_k1

  ALLOCATE(TS1(7,1:Nx, 1:Ny ), TS2(7,1:Nx, 1:Ny), CONS_k1(7,1:Nx, 1:Ny))
  traceP = 0.d0; H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0; fracp11 = 0.d0; fracp22 = 0.d0; Q = 0.d0; alpha = 0.d0
  cons_k1 = 0.d0; TS1 = 0.d0; TS2 = 0.d0
  dt2 = 0.5d0*dt
!------------------------------------------------------------------------------------------------------------------------
  DO  ix = 1, Nx
    DO iy  = 1, Ny  
!------------------------------------------------------------------------------------------------------------------------ 
        if (cons(1,ix, iy).le.1.d-8) then
           print*, 'pas de l eau', cons(1,ix,iy), ix, iy, it
           stop
        end if

       !TERME SOURCE
            H = Prim(H_pv,ix, iy)
            P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
            traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
            fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
            alpha = disscoeff*(traceP - (H**2.d0)*phi2)/( traceP**2.d0)
             alpha = dmax1(alpha, 0.d0)
            Q = alpha*traceP**((dsqrt(u**2.d0 + v**2.d0))**3.d0)

        IF (penteY  == 1) THEN       
            TS1(2,ix, iy) =  - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
            TS1(3,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
            TS1(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(7,ix, iy) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
       ELSE IF (penteX  == 1) THEN 
            TS1(2,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
            TS1(3,ix, iy) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
            TS1(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(7,ix, iy) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
        ENDIF
        !------------------------------------------------------------------------------------------------------------------------
            do k=1,7
              CONS_k1(k,ix, iy) = CONS(k,ix, iy) + DT*TS1(k,ix, iy) 
            end do 
              

        ENDDO
        ENDDO  
        !------------------------------------------------------------------------------------------------------------------------
          CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_k1, Pression, it)
          call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
        !------------------------------------------------------------------------------------------------------------------------
          DO  ix = 1, Nx
            DO iy  = 1, Ny  
                H = Prim(H_pv,ix, iy)
                P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
                traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
                fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
                alpha = disscoeff*(traceP - (H**2.d0)*phi2 )/( traceP**2.d0)
                 alpha = dmax1(alpha, 0.d0)
                Q = alpha*traceP**((dsqrt(u**2.d0 + v**2.d0))**3.d0)
     
                IF (penteY  == 1) THEN       
                    TS2(2,ix, iy) =  - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS2(3,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS2(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(7,ix, iy) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
               ELSE IF (penteX  == 1) THEN 
                    TS2(2,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS2(3,ix, iy) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS2(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(7,ix, iy) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
                ENDIF
            !------------------------------------------------------------------------------------------------------------------------

                do k=1,7
                  CONS(k,ix, iy) = CONS(k,ix, iy) + DT2*(TS1(k,ix, iy) + TS2(k,ix, iy) )
                end do 
          ENDDO
      ENDDO
        !------------------------------------------------------------------------------------------------------------------------
      CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS, Pression, it)
      call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
     
        !------------------------------------------------------------------------------------------------------------------------
  DEALLOCATE(TS1, TS2, cons_k1)
  return
  END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 SUBROUTINE rk4(DT, CONS, Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it)
  USE precisions
  use GlobalParam
  implicit none 
  INTEGER :: ix, iy, IT,k
  REAL (KIND = DP) :: DT, dt2, traceP, H, p11, p22, p12, u, v,Q, alpha
  REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny)
  REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1),SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1),SOUND_by(0:Nx+1,0:Ny+1)
  REAL (KIND = DP) :: Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1), fracp11, fracp22, erreur
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: TS1, TS2, TS3, TS4
  REAL (KIND = DP), ALLOCATABLE, DIMENSION(:,:,:)  :: CONS_k1, CONS_k2,CONS_k3

  ALLOCATE(TS1(7,1:Nx, 1:Ny ), TS2(7,1:Nx, 1:Ny), TS3(7,1:Nx, 1:Ny), TS4(7,1:Nx, 1:Ny))
  ALLOCATE(CONS_k1(7,1:Nx, 1:Ny), CONS_k2(7,1:Nx, 1:Ny),CONS_k3(7,1:Nx, 1:Ny))
  traceP = 0.d0; H= 0.d0; p11= 0.d0; p22= 0.d0; p12= 0.d0; u= 0.d0; v= 0.d0; fracp11 = 0.d0; fracp22 = 0.d0; Q = 0.d0; alpha = 0.d0
  cons_k1 = 0.d0;cons_k2 = 0.d0;cons_k3 = 0.d0; TS1 = 0.d0; TS2 = 0.d0; TS3 = 0.d0; TS4 = 0.d0
  dt2 = 0.5d0*dt
!------------------------------------------------------------------------------------------------------------------------
  DO  ix = 1, Nx
    DO iy  = 1, Ny  
!------------------------------------------------------------------------------------------------------------------------ 
        if (cons(1,ix, iy).le.1.d-8) then
           print*, 'pas de l eau', cons(1,ix,iy), ix, iy, it
           stop
        end if
!------------------------------------------------------------------------------------------------------------------------
       !TERME SOURCE
            H = Prim(H_pv,ix, iy)
            P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
            traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
            fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
            alpha = disscoeff*(traceP - (H**2.d0)*phi2 + eps)/( traceP**2.d0)
             alpha = dmax1(alpha, 0.d0)
            Q = alpha*traceP*((dsqrt(u**2.d0 + v**2.d0))**3.d0)

        IF (penteY  == 1) THEN       
            TS1(2,ix, iy) =  - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u 
            TS1(3,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
            TS1(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(7,ix, iy) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q 
       ELSE IF (penteX  == 1) THEN 
            TS1(2,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u 
            TS1(3,ix, iy) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
            TS1(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
            TS1(7,ix, iy) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q 
        ENDIF
        !------------------------------------------------------------------------------------------------------------------------
            do k=1,7
              CONS_k1(k,ix, iy) = CONS(k,ix, iy) + DT*TS1(k,ix, iy) 
            end do 
         
              ! ts1(4,ix, iy)= - dt*fracp11*(TS1(2,ix, iy)**2.d0 +TS1(3,ix, iy)**2.d0 )/h                               
              ! ts1(6,ix, iy)= - dt*fracp22*(TS1(2,ix, iy)**2.d0 +TS1(3,ix, iy)**2.d0 )/h                          
              ! CONS_k1(4,ix, iy) = CONS_k1(4,ix, iy) + dt*TS1(4,ix, iy) 
               !CONS_k1(6,ix, iy) = CONS_k1(6,ix, iy) + dt*TS1(6,ix, iy) 

        ENDDO
        ENDDO  
        !------------------------------------------------------------------------------------------------------------------------
          CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_k1, Pression, it)
          call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
        !------------------------------------------------------------------------------------------------------------------------
          DO  ix = 1, Nx
            DO iy  = 1, Ny  
                H = Prim(H_pv,ix, iy)
                P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
                traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
                fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
                alpha = disscoeff*(traceP - (H**2.d0)*phi2 + eps)/( traceP**2.d0)
                 alpha = dmax1(alpha, 0.d0)
                Q = alpha*traceP*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
     
                IF (penteY  == 1) THEN       
                    TS2(2,ix, iy) =  - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS2(3,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS2(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(7,ix, iy) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
               ELSE IF (penteX  == 1) THEN 
                    TS2(2,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS2(3,ix, iy) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS2(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS2(7,ix, iy) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
                ENDIF
            !------------------------------------------------------------------------------------------------------------------------

                do k=1,7
                  CONS_k2(k,ix, iy) = CONS(k,ix, iy) + DT2*(TS2(k,ix, iy) )
                end do 
              
                 !  ts2(4,ix, iy)= - 0.25d0*dt2*fracp11*(TS2(2,ix, iy)**2.d0 +TS2(3,ix, iy)**2.d0 )/h                               !fracp11*erreur
                 !  ts2(6,ix, iy)= - 0.25d0*dt2*fracp22*(TS2(2,ix, iy)**2.d0 +TS2(3,ix, iy)**2.d0 )/h                               !fracp22*erreur

                  ! CONS_k2(4,ix, iy) = CONS_k2(4,ix, iy) + dt2*TS2(4,ix, iy) !!verifier le signe
                  ! CONS_k2(6,ix, iy) = CONS_k2(6,ix, iy) + dt2*TS2(6,ix, iy)
          ENDDO
      ENDDO
        !------------------------------------------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------------------------------------------
          CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_k2, Pression, it)
          call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
        !------------------------------------------------------------------------------------------------------------------------
          DO  ix = 1, Nx
            DO iy  = 1, Ny  
                H = Prim(H_pv,ix, iy)
                P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
                traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
                fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
                alpha = disscoeff*(traceP - (H**2.d0)*phi2 + eps)/( traceP**2.d0)
                alpha = dmax1(alpha, 0.d0)
                Q = alpha*traceP*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
     
                IF (penteY  == 1) THEN       
                    TS3(2,ix, iy) =  - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS3(3,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS3(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS3(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS3(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS3(7,ix, iy) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
               ELSE IF (penteX  == 1) THEN 
                    TS3(2,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS3(3,ix, iy) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS3(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS3(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS3(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS3(7,ix, iy) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
                ENDIF
            !------------------------------------------------------------------------------------------------------------------------

                do k=1,7
                  CONS_k3(k,ix, iy) = CONS(k,ix, iy) + DT2*(TS3(k,ix, iy) )
                end do 
              
                 !  ts3(4,ix, iy)= - 0.25d0*dt2*fracp11*(TS3(2,ix, iy)**2.d0+ TS3(3,ix, iy)**2.d0 )/h                               !fracp11*erreur
                 !  ts3(6,ix, iy)= - 0.25d0*dt2*fracp22*(TS2(2,ix, iy)**2.d0 +TS3(3,ix, iy)**2.d0 )/h                               !fracp22*erreur

                  ! CONS_k3(4,ix, iy) = CONS_k3(4,ix, iy) + dt2*TS3(4,ix, iy) !!verifier le signe
                  ! CONS_k3(6,ix, iy) = CONS_k3(6,ix, iy) + dt2*TS3(6,ix, iy)
          ENDDO
      ENDDO
        !------------------------------------------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------------------------------------------
          CALL NOUVELL_VARIABLE_PRIM(Prim, Ein, SOUND_ax, SOUND_bx,  SOUND_ay, SOUND_by, CONS_k3, Pression, it)
          call COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
        !------------------------------------------------------------------------------------------------------------------------
          DO  ix = 1, Nx
            DO iy  = 1, Ny  
                H = Prim(H_pv,ix, iy)
                P11 = Prim(P11_pv,ix, iy); p22 = Prim(P22_pv,ix, iy); p12 = Prim(P12_pv,ix, iy)
                traceP = P11+ P22 ; u =  Prim(U_pv,ix, iy); v = Prim(V_pv,ix, iy)
                fracp11=p11/(P11+ P22); fracp22=p22/(P11+ P22);
                alpha = disscoeff*(traceP - (H**2.d0)*phi2 + eps)/( traceP**2.d0)
                alpha = dmax1(alpha, 0.d0)
                Q = alpha*traceP*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
     
                IF (penteY  == 1) THEN       
                    TS4(2,ix, iy) =  - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS4(3,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS4(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS4(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS4(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS4(7,ix, iy) =   g*Dtan(angle)*H*v - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
               ELSE IF (penteX  == 1) THEN 
                    TS4(2,ix, iy) = g*Dtan(angle)*H - frottcoeff*dsqrt( u**2.d0 + v**2.d0 )*u ! g*Dtan(angle)*H  
                    TS4(3,ix, iy) = - frottcoeff*dsqrt(u**2.d0 + v**2.d0 )*v  
                    TS4(4,ix, iy) = - 2.d0*alpha*p11*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS4(5,ix, iy) = - 2.d0*alpha/H*p12*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS4(6,ix, iy) = - 2.d0*alpha*p22*((dsqrt(u**2.d0 + v**2.d0))**3.d0)
                    TS4(7,ix, iy) =   g*Dtan(angle)*H*U - frottcoeff*(dsqrt(u**2.d0 + v**2.d0))**(3.d0) - Q ! g*Dtan(angle)*H*U 
                ENDIF
            !------------------------------------------------------------------------------------------------------------------------

                do k=1,7
                  CONS(k,ix, iy) = CONS(k,ix, iy) + DT*(TS1(k,ix, iy) + 2.d0*TS2(k,ix, iy)+ 2.d0*TS3(k,ix, iy)+ TS4(k,ix, iy) )/6.d0
                end do 
              
                  ! ts4(4,ix, iy)= - dt*fracp11*((TS1(k,ix, iy) + 2.d0*TS2(k,ix, iy)+ 2.d0*TS3(k,ix, iy)+ TS4(k,ix, iy) )**2.d0 +(TS1(k,ix, iy) + 2.d0*TS2(k,ix, iy)+ 2.d0*TS3(k,ix, iy)+ TS4(k,ix, iy) )**2.d0 )/(36.d0*h )                              !fracp11*erreur
                  ! ts4(6,ix, iy)= - dt*fracp22*((TS1(k,ix, iy) + 2.d0*TS2(k,ix, iy)+ 2.d0*TS3(k,ix, iy)+ TS4(k,ix, iy) )**2.d0 +(TS1(k,ix, iy) + 2.d0*TS2(k,ix, iy)+ 2.d0*TS3(k,ix, iy)+ TS4(k,ix, iy) )**2.d0 )/(36.d0*h )                             !fracp22*erreur

                  ! CONS(4,ix, iy) = CONS(4,ix, iy) + dt*TS4(4,ix, iy) !!verifier le signe
                  ! CONS(6,ix, iy) = CONS(6,ix, iy) + dt*TS4(6,ix, iy)
          ENDDO
      ENDDO
        !------------------------------------------------------------------------------------------------------------------------
     
  DEALLOCATE(TS1, TS2, cons_k1, cons_k2, cons_k3)
  return
  END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



SUBROUTINE PENTE1_x( Prim, PENTE)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy
REAL (KIND = DP) :: pent, pente11, pente12
REAL (KIND = DP) :: PENTE(6,0:Nx+1,0:Ny+1), Prim(6,0:Nx+1,0:Ny+1)
pente = 0.d0

 DO ix= 1,Nx
  DO iy = 1, Ny
    DO iv = 1, Nv_Prim
      PENTE11 = (Prim(iv, ix,iy) - Prim(iv, ix- 1, iy)) 
      PENTE12=   (Prim(iv, ix + 1,iy) - Prim(iv, ix, iy)) 
      call minmod(pente11,pente12,pent)    
      PENTE(iv, ix, iy)=pent
    END DO
  END DO
END DO

! DO iy = 0, Ny+1
!   DO iv = 1, Nv_Prim
!      PENTE(iv, 0, iy) =  0.d0 
!      PENTE(iv, Nx + 1, iy) =  0.d0 
!   END DO
! END DO

! DO ix = 0, Nx+1
!   DO iv = 1, Nv_Prim
!      PENTE(iv, ix, 0) =  0.d0 
!      PENTE(iv, ix, Ny + 1) =  0.d0 
!   END DO
! END DO

return
END SUBROUTINE

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PENTE2_x( Prim, PENTE)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy
REAL (KIND = DP) :: pent, pente11, pente12
REAL (KIND = DP) :: PENTE(6,0:Nx+1,0:Ny+1), Prim(6,0:Nx+1,0:Ny+1)
pente = 0.d0 

 DO ix= 1,Nx
  DO iy = 1, Ny
     DO iv = 1, Nv_Prim
  PENTE11 = (Prim(iv, ix,iy) - Prim(iv, ix- 1, iy)) 
  PENTE12=   (Prim(iv, ix + 1,iy) - Prim(iv, ix, iy)) 
  call vleer(pente11,pente12,pent)    
  PENTE(iv, ix, iy)=pent

END DO
END DO
END DO

! DO iy = 0, Ny+1
!   DO iv = 1, Nv_Prim
!      PENTE(iv, 0, iy) =  0.d0 
!      PENTE(iv, Nx + 1, iy) =  0.d0 
!   END DO
! END DO

! DO ix = 0, Nx+1
!   DO iv = 1, Nv_Prim
!      PENTE(iv, ix, 0) =  0.d0 
!      PENTE(iv, ix, Ny + 1) =  0.d0 
!   END DO
! END DO

return
END SUBROUTINE

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PENTE1_y( Prim, PENTE)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy
REAL (KIND = DP) :: pent, pente11, pente12
REAL (KIND = DP) :: Prim(Nv_prim, 0:Nx+1, 0:Ny+1), PENTE(6,0:Nx+1,0:Ny+1)
 
pente = 0.d0 
 DO ix= 1,Nx
  DO iy = 1, Ny
    DO iv = 1, Nv_Prim
      PENTE11 = (Prim(iv, ix,iy) - Prim(iv, ix, iy-1)) 
      PENTE12=   (Prim(iv, ix,iy +1) - Prim(iv, ix, iy)) 
      call minmod(pente11,pente12,pent)    
      PENTE(iv, ix, iy)=pent

    END DO
  END DO
END DO

! DO iy = 1, Ny
!   DO iv = 1, Nv_Prim
!      PENTE(iv,  0,   iy) =  0.d0 
!      PENTE(iv, Nx+1, iy) =  0.d0 
!   END DO
! END DO

! DO ix = 1, Nx
!   DO iv = 1, Nv_Prim
!      PENTE(iv, ix, 0) =  0.d0 
!      PENTE(iv, ix, Ny + 1) =  0.d0 
!   END DO
! END DO

return
END SUBROUTINE

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
SUBROUTINE PENTE2_y( Prim, PENTE)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy, iv
REAL (KIND = DP) :: pent, pente11, pente12
REAL (KIND = DP) :: Prim(Nv_prim, 0:Nx+1, 0:Ny+1), PENTE(6,0:Nx+1,0:Ny+1)
 
pente = 0.d0 
 
 DO ix= 1,Nx
  DO iy = 1, Ny
    DO iv = 1, Nv_Prim
      PENTE11 = (Prim(iv, ix,iy) - Prim(iv, ix, iy-1)) 
      PENTE12=   (Prim(iv, ix,iy +1) - Prim(iv, ix, iy)) 
      call vleer(pente11,pente12,pent)    
      PENTE(iv, ix, iy)=pent

    END DO
  END DO
END DO
!print*, 'sub pente', pente(u_pv,:,:)
!pause
! DO iy = 1, Ny
!   DO iv = 1, Nv_Prim
!      PENTE(iv,  0,   iy) =  0.d0 
!      PENTE(iv, Nx+1, iy) =  0.d0 
!   END DO
! END DO

! DO ix = 1, Nx
!   DO iv = 1, Nv_Prim
!      PENTE(iv, ix,    0) =  0.d0 
!      PENTE(iv, ix, Ny+1) =  0.d0 
!   END DO
! END DO

return
END SUBROUTINE

 !limiteur de pente:

subroutine vleer(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,prod   
if ((DABS(s1) .lt. 1.d-8) .and. (DABS(s1) .lt. 1.d-8)) then
  slim = 0.d0
  return
endif

prod=s1*s2
if(prod>0.d0.and.dabs(s1+s2)>1.d-6) then
  slim=2.d0*prod/(s1+s2)
else
  slim=0.d0
endif
return
end subroutine vleer

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

subroutine minmod(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,ss1,ss2

if(s1*s2 > 0.d0) then
  slim=dabs(s1)
  if(dabs(s2)<slim) slim=dabs(s2)
  if(s1 < 0.d0)slim=-slim
else
  slim=0.d0
endif
return
end subroutine minmod
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!    module limiteur de pentes superbee

subroutine superbee(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,ss1,ss2      

if(s1*s2.gt.1.d-6) then
  ss1=dabs(s1)
  ss2=dabs(s2)
  slim=dmax1(dmin1(2.d0*ss1,ss2),dmin1(ss1,2.d0*ss2))
  if(s1.lt.0.d0)slim=-slim
else
  slim=0.d0
endif
return
end subroutine superbee
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !    module limiteur de pentes Van Albada

subroutine valbada(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim  

if(s1.ne.0.d0.or.s2.ne.0.d0) then
  slim=s1*s2*(s1+s2)/(s1*s1+s2*s2)
else
  slim=0.d0
endif
return
end subroutine valbada
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !    module limiteur de pentes MC

subroutine vmc(s1,s2,slim)
use precisions
implicit none
real(kind=dp) :: s1,s2,slim,ss1,ss2   

if(s1*s2.gt.0.d0) then
  ss1=dabs(s1)
  ss2=dabs(s2)
  slim=dmin1(2.d0*ss1,2.d0*ss2,0.5d0*(ss1+ss2))
  if(s1.lt.0.d0) slim=-slim
else
  slim=0.d0
endif
return
end subroutine vmc

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PRIM_TO_CONS_FLUX_sub1x(Prim, Ein,Pression, CONS, FLUX)
USE precisions
  use GlobalParam
IMPLICIT NONE
REAL (KIND = DP), DIMENSION(0:Nx+1,0:Ny+1 ):: Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression
REAL (KIND = DP) :: CONS(7,1:Nx,1:Ny),  FLUX(7,0:Nx,0:Ny), Prim(Nv_Prim,0:Nx+1,0:Ny+1)
REAL (KIND = DP) :: h, u, v, p11, p12, p22
INTEGER :: ix, iy

! VARIABLE CONSERVATIVES 
   DO ix = 1, Nx
    DO iy = 1, Ny
      h = Prim(H_pv, ix, iy); u = Prim(U_pv, ix, iy); v = Prim(V_pv, ix, iy)
      p11 = Prim(P11_pv, ix, iy);  p12 = Prim(P12_pv, ix, iy); p22 = Prim(P22_pv, ix, iy)

        CONS(1,ix,iy) = h
        CONS(2,ix,iy) = h*u
        CONS(3,ix,iy) = h*v
        CONS(4,ix,iy) = h*p11
        CONS(5,ix,iy) = p12
        CONS(6,ix,iy) = h*p22
        CONS(7,ix,iy) = h*(Ein(ix,iy)+ (u**2.d0 + v**2.d0 )/2.d0)

       flux(1,ix,iy)=h*u
       flux(2,ix,iy)=h*u*u+Pression(ix, iy)
       flux(3,ix,iy)=h*u*v
       flux(4,ix,iy)=p11*u !! equation inutile pour le cas x
       flux(5,ix,iy)=p12*u
       flux(6,ix,iy)=h*p22*u
       flux(7,ix,iy)=h*u*(Ein(ix, iy) + 0.5d0*u**2.d0 +0.5d0*v**2.d0)+Pression(ix, iy)*u
  
  END DO
  END DO


DO iy = 1, Ny
   h = Prim(H_pv, 0, iy); u = Prim(U_pv, 0, iy); v = Prim(V_pv, 0, iy)
  p11 = Prim(P11_pv, 0, iy);  p12 = Prim(P12_pv, 0, iy); p22 = Prim(P22_pv, 0, iy)
   flux(1,0,iy)=h*u
   flux(2,0,iy)=h*u*u+Pression(ix, iy)
   flux(3,0,iy)=h*u*v
   flux(4,0,iy)=p11*u !! equation inutile pour le cas x
   flux(5,0,iy)=p12*u
   flux(6,0,iy)=h*p22*u
   flux(7,0,iy)=h*u*(Ein(ix, iy) + 0.5d0*u**2.d0 +0.5d0*v**2.d0)+Pression(ix, iy)*u
END DO

DO ix = 1, Nx
  h = Prim(H_pv, ix, 0); u = Prim(U_pv, ix, 0); v = Prim(V_pv, ix, 0)
  p11 = Prim(P11_pv, ix, 0);  p12 = Prim(P12_pv, ix, 0); p22 = Prim(P22_pv, ix, 0)

   flux(1,ix,0)=h*u
   flux(2,ix,0)=h*u*u+Pression(ix, iy)
   flux(3,ix,0)=h*u*v
   flux(4,ix,0)=p11*u !! equation inutile pour le cas x
   flux(5,ix,0)=p12*u
   flux(6,ix,0)=h*p22*u
   flux(7,ix,0)=h*u*(Ein(ix, iy) + 0.5d0*u**2.d0 +0.5d0*v**2.d0)+Pression(ix, iy)*u
END DO
return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PRIM_TO_CONS_FLUX_sub2x(Prim, Ein, CONS, FLUX)
USE precisions
use GlobalParam
IMPLICIT NONE
REAL (KIND = DP), DIMENSION(0:Nx+1,0:Ny+1 ):: Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression
REAL (KIND = DP) :: CONS(7,1:Nx,1:Ny),  FLUX(7,0:Nx,0:Ny), Prim(Nv_Prim,0:Nx+1,0:Ny+1)
REAL (KIND = DP) :: h, u, v, p11, p12, p22
INTEGER :: ix, iy
! VARIABLE CONSERVATIVES 
   DO ix = 1, Nx
    DO iy = 1, Ny
      h = Prim(H_pv, ix, iy); u = Prim(U_pv, ix, iy); v = Prim(V_pv, ix, iy)
      p11 = Prim(P11_pv, ix, iy);  p12 = Prim(P12_pv, ix, iy); p22 = Prim(P22_pv, ix, iy)

      CONS(1,ix,iy) = h
      CONS(2,ix,iy) = h*u
      CONS(3,ix,iy) = h*v
      CONS(4,ix,iy) = h*p11
      CONS(5,ix,iy) = p12
      CONS(6,ix,iy) = h*p22
      CONS(7,ix,iy) = h*(Ein(ix,iy)+ (u**2.d0 + v**2.d0 )/2.d0)

      Flux(1,ix,iy) = 0.d0
      Flux(2,ix,iy) = 0.d0
      flux(3,ix,iy) = h*p12
      Flux(4,ix,iy) = 0.d0
      Flux(5,ix,iy) = v ! non conservatif
      flux(6,ix,iy) = 0.d0   !2.d0*hp12star*vstar !!inutile
      flux(7,ix,iy) = h*p12*v
  
  END DO
  END DO

  DO ix = 1, Nx
      h = Prim(H_pv, ix, 0); u = Prim(U_pv, ix, 0); v = Prim(V_pv, ix, 0)
      p11 = Prim(P11_pv, ix, 0);  p12 = Prim(P12_pv, ix, 0); p22 = Prim(P22_pv, ix, 0)

      Flux(1,ix,0) = 0.d0
      Flux(2,ix,0) = 0.d0
      flux(3,ix,0) = h*p12
      Flux(4,ix,0) = 0.d0
      Flux(5,ix,0) = v ! non conservatif
      flux(6,ix,0) = 0.d0   !2.d0*hp12star*vstar !!inutile
      flux(7,ix,0) = h*p12*v
  END DO

  DO iy = 1, Ny
    h = Prim(H_pv, 0, iy); u = Prim(U_pv, 0, iy); v = Prim(V_pv, 0, iy)
    p11 = Prim(P11_pv, 0, iy);  p12 = Prim(P12_pv, 0, iy); p22 = Prim(P22_pv, 0, iy)

  Flux(1,0,iy) = 0.d0
  Flux(2,0,iy) = 0.d0
  flux(3,0,iy) = h*p12
  Flux(4,0,iy) = 0.d0
  Flux(5,0,iy) = v ! non conservatif
  flux(6,0,iy) = 0.d0   !2.d0*hp12star*vstar !!inutile
  flux(7,0,iy) = h*p12*v
  END DO

return
END SUBROUTINE
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PRIM_TO_CONS_FLUX_sub1y(Prim, Ein, CONS, FLUX)
USE precisions
use GlobalParam
IMPLICIT NONE
REAL (KIND = DP), DIMENSION(0:Nx+1,0:Ny+1 ):: Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by
REAL (KIND = DP) :: CONS(7,1:Nx,1:Ny),  FLUX(7,0:Nx,0:Ny), Prim(Nv_Prim,0:Nx+1,0:Ny+1)
REAL (KIND = DP) :: h, u, v, p11, p12, p22, temp
INTEGER :: ix, iy

! VARIABLE CONSERVATIVES 
   DO ix = 1, Nx
    DO iy = 1, Ny
      h = Prim(H_pv, ix, iy); v = Prim(v_pv, ix, iy); u = Prim(u_pv, ix, iy)
      p22 = Prim(P22_pv, ix, iy);  p12 = Prim(P12_pv, ix, iy); p11 = Prim(P11_pv, ix, iy)

      CONS(1,ix,iy) = h
      CONS(2,ix,iy) = h*u
      CONS(3,ix,iy) = h*v
      CONS(4,ix,iy) = h*p11
      CONS(5,ix,iy) = p12
      CONS(6,ix,iy) = h*p22
      CONS(7,ix,iy) = h*(Ein(ix,iy)+ (u**2.d0 + v**2.d0 )/2.d0)

     flux(1,ix,iy)=h*v
     flux(2,ix,iy)=h*u*v 
     flux(3,ix,iy)=h*v*v + 0.5d0*g*h**2.d0 + h*p22
     flux(4,ix,iy)=h*v*p11 
     flux(5,ix,iy)=p12*v
     flux(6,ix,iy)=h*p22*v !! equation inutile pour le cas y
     flux(7,ix,iy)=h*v*(Ein(ix, iy) + 0.5d0*u**2.d0 +0.5d0*v**2.d0)+(0.5d0*g*h**2.d0 + h*p22)*v

  END DO
  END DO

  DO ix = 1, Nx
      h = Prim(H_pv, ix, 0); v = Prim(U_pv, ix, 0); u = Prim(V_pv, ix, 0)
      p22 = Prim(P11_pv, ix, 0);  p12 = Prim(P12_pv, ix, 0); p11 = Prim(P22_pv, ix, 0)

     flux(1,ix,0)=h*v
     flux(2,ix,0)=h*u*v 
     flux(3,ix,0)=h*v*v+0.5d0*g*h**2.d0+h*p22
     flux(4,ix,0)=h*v*p11 
     flux(5,ix,0)=p12*v
     flux(6,ix,0)=h*p22*v !! equation inutile pour le cas y
     flux(7,ix,0)=h*v*(Ein(ix, 0) + 0.5d0*u**2.d0 +0.5d0*v**2.d0)+(0.5d0*g*h**2.d0 + h*p22)*v
  END DO

    DO iy = 1, Ny
      h = Prim(H_pv, 0, iy); v = Prim(U_pv, 0, iy); u = Prim(V_pv, 0, iy)
      p22 = Prim(P11_pv, 0, iy);  p12 = Prim(P12_pv, 0, iy); p11 = Prim(P22_pv, 0, iy)

     flux(1,0,iy)=h*v
     flux(2,0,iy)=h*u*v 
     flux(3,0,iy)=h*v*v+0.5d0*g*h**2.d0+h*p22
     flux(4,0,iy)=h*v*p11 
     flux(5,0,iy)=p12*v
     flux(6,0,iy)=h*p22*v !! equation inutile pour le cas y
     flux(7,0,iy)=h*v*(Ein(0, iy) + 0.5d0*u**2.d0 +0.5d0*v**2.d0)+(0.5d0*g*h**2.d0 + h*p22)*v   
  END DO

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE PRIM_TO_CONS_FLUX_sub2y(Prim, Ein, CONS, FLUX)
USE precisions
  use GlobalParam
IMPLICIT NONE
REAL (KIND = DP), DIMENSION(0:Nx+1,0:Ny+1 ):: Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by
REAL (KIND = DP) :: CONS(7,1:Nx,1:Ny),  FLUX(7,0:Nx,0:Ny), Prim(Nv_Prim,0:Nx+1,0:Ny+1)
REAL (KIND = DP) :: h, u, v, p11, p12, p22, temp
INTEGER :: ix, iy

! VARIABLE CONSERVATIVES 
   DO ix = 1, Nx
    DO iy = 1, Ny
      h = Prim(H_pv, ix, iy); v = Prim(v_pv, ix, iy); u = Prim(u_pv, ix, iy)
      p22 = Prim(P22_pv, ix, iy);  p12 = Prim(P12_pv, ix, iy); p11 = Prim(P11_pv, ix, iy)

      CONS(1,ix,iy) = h
      CONS(2,ix,iy) = h*u
      CONS(3,ix,iy) = h*v
      CONS(4,ix,iy) = h*p11
      CONS(5,ix,iy) = p12
      CONS(6,ix,iy) = h*p22
      CONS(7,ix,iy) = h*(Ein(ix,iy)+ (u**2.d0 + v**2.d0 )/2.d0)

     flux(1,ix,iy)=0.d0
     flux(2,ix,iy)=h*p12 
     flux(3,ix,iy)=0.d0
     flux(4,ix,iy)=0.d0 
     flux(5,ix,iy)=u
     flux(6,ix,iy)=0.d0 !! equation inutile pour le cas y
     flux(7,ix,iy)=h*p12*u

  END DO
  END DO

  DO ix = 0, Nx
      h = Prim(H_pv, ix, 0); v = Prim(U_pv, ix, 0); u = Prim(V_pv, ix, 0)
      p22 = Prim(P11_pv, ix, 0);  p12 = Prim(P12_pv, ix, 0); p11 = Prim(P22_pv, ix, 0)

     flux(1,ix,0)=0.d0
     flux(2,ix,0)=h*p12 
     flux(3,ix,0)=0.d0
     flux(4,ix,0)=0.d0
     flux(5,ix,0)=u
     flux(6,ix,0)=0.d0!! equation inutile pour le cas y
     flux(7,ix,0)=h*p12*u
  END DO

    DO iy = 1, Ny
      h = Prim(H_pv, 0, iy); v = Prim(U_pv, 0, iy); u = Prim(V_pv, 0, iy)
      p22 = Prim(P11_pv, 0, iy);  p12 = Prim(P12_pv, 0, iy); p11 = Prim(P22_pv, 0, iy)

     flux(1,0,iy)=0.d0
     flux(2,0,iy)=h*u*v 
     flux(3,0,iy)=0.d0
     flux(4,0,iy)=0.d0
     flux(5,0,iy)=u
     flux(6,0,iy)=0.d0 !! equation inutile pour le cas y
     flux(7,0,iy)=h*p12*u   
  END DO


return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE CONS_PREDICTION_sub1x( DX, DT, CONS, CONS_PRED, FLUXGG, FLUXDD)
USE precisions
  use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy
REAL (KIND = DP) :: DX, DT, hu, hv, h, hp22, hE
REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny), FLUXGG(7,0:Nx,0:Ny), FLUXDD(7,0:Nx,0:Ny), CONS_PRED(7,1:Nx, 1:Ny)

CONS_PRED = 0.d0

DO ix=1,Nx
  DO iy=1,Ny
    CONS_PRED(:,ix, iy) = CONS(:,ix, iy) - 0.5D0*DT/DX*(FLUXDD(:,ix, iy) - FLUXGG(:,ix, iy))  

hu=CONS_PRED(2,ix,iy); hv=CONS_PRED(3,ix,iy); h=CONS_PRED(1,ix,iy); hp22=CONS_PRED(6,ix,iy); hE=CONS_PRED(7,ix,iy)
CONS_PRED(4,ix,iy)=2.d0*hE-g*h*h-hp22-(hu*hu+hv*hv)/h

  END DO
END DO


return
END SUBROUTINE

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE CONS_PREDICTION_sub2x( DX, DT, CONS, CONS_PRED, FLUXGG, FLUXDD)
USE precisions
  use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy, k
REAL (KIND = DP) :: DX, DT, hu, hv, h, hp22, hE, p11, hp11
REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny), FLUXGG(7,0:Nx,0:Ny), FLUXDD(7,0:Nx,0:Ny), CONS_PRED(7,1:Nx, 1:Ny)

CONS_PRED = 0.d0

 do ix=1,nx
   do iy=1,ny
    p11=cons(4,ix,iy)/cons(1,ix,iy)
!------------------------------------------------------------------------------------------------------------------------
    do k=1,4
      CONS_PRED(K,ix, iy) = CONS(k,ix, iy) - 0.5D0*DT/DX*(FLUXDD(k,ix, iy) - FLUXGG(k,ix, iy)) 
   end do
!------------------------------------------------------------------------------------------------------------------------
      CONS_PRED(5,ix,iy)=CONS(5,ix,iy)-0.5d0*dt/dx*(FluxDD(5,ix,iy)-fluxGG(5,ix,iy))*p11 !! ici dans le flux on a stocké vstar
!------------------------------------------------------------------------------------------------------------------------
   do k=6,7
    CONS_PRED(k,ix,iy)=CONS(k,ix,iy)-0.5d0*dt/dx*(FluxDD(k,ix,iy)-fluxGG(k,ix,iy))
   end do
  
    !! corection de la variable p22
    hu=CONS_PRED(2,ix,iy); hv=CONS_PRED(3,ix,iy); h=CONS_PRED(1,ix,iy); hp11=CONS_PRED(4,ix,iy); hE=CONS_PRED(7,ix,iy)
    CONS_PRED(6,ix,iy)=2.d0*hE-(hu*hu+hv*hv)/h-g*h*h-hp11
   end do
  end do

return
END SUBROUTINE

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE CONS_PREDICTION_sub1y( Dy, DT, CONS, CONS_PRED, FLUXGG, FLUXDD)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: iy, ix
REAL (KIND = DP) :: Dy, DT, hu, hv, h, hp22, hE, hp11
REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny), FLUXGG(7,0:Nx,0:Ny), FLUXDD(7,0:Nx,0:Ny), CONS_PRED(7,1:Nx, 1:Ny)

CONS_PRED = 0.d0

  DO ix=1,Nx
    DO iy=1,Ny
      CONS_PRED(:,ix, iy) = CONS(:,ix, iy) - 0.5D0*DT/Dy*(FLUXDD(:,ix, iy) - FLUXGG(:,ix, iy))   

      !! corection de la variable p22
      hu=CONS_PRED(2,ix,iy); hv=CONS_PRED(3,ix,iy); h=CONS_PRED(1,ix,iy); hp11=CONS_PRED(4,ix,iy); hE=CONS_PRED(7,ix,iy)
      CONS_PRED(6,ix,iy)=2.d0*hE-g*h*h-hp11-(hu*hu+hv*hv)/h
    END DO
  END DO

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE CONS_PREDICTION_sub2y(Dy, DT, CONS, CONS_PRED, FLUXGG, FLUXDD)
USE precisions
use GlobalParam
IMPLICIT NONE
INTEGER :: ix, iy, k
REAL (KIND = DP) :: Dy, DT, hu, hv, h, hp22, hE, p22
REAL (KIND = DP) :: CONS(7,1:Nx, 1:Ny), FLUXGG(7,0:Nx,0:Ny), FLUXDD(7,0:Nx,0:Ny), CONS_PRED(7,1:Nx, 1:Ny)

CONS_PRED = 0.d0

do ix=1,nx
   do iy=1,ny
       p22=cons(6,ix,iy)/cons(1,ix,iy)
     
       do k=1,4
        CONS_PRED(k,ix,iy)=CONS(k,ix,iy)-0.5d0*dt/dy*(FluxDD(k,ix,iy)-fluxGG(k,ix,iy))
       end do
        CONS_PRED(5,ix,iy)=CONS(5,ix,iy)-0.5d0*dt/dy*(FluxDD(5,ix,iy)-fluxGG(5,ix,iy))*p22
       do k=6,7
        CONS_PRED(k,ix,iy)=CONS(k,ix,iy)-0.5d0*dt/dy*(FluxDD(k,ix,iy)-fluxGG(k,ix,iy))
       end do

        !! corection de la variable p11
        hu=CONS_PRED(2,ix,iy);hv=CONS_PRED(3,ix,iy);h=CONS_PRED(1,ix,iy);hp22=CONS_PRED(6,ix,iy);hE=CONS_PRED(7,ix,iy)
        CONS_PRED(4,ix,iy)=2.d0*hE-(hu*hu+hv*hv)/h-g*h*h-hp22
   end do
  end do
! print*, 'cons_pred', CONS_PRED(1,:,:)
! pause

return
END SUBROUTINE
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SUBROUTINE COMMUNICATION(Prim, Ein, SOUND_ax,SOUND_bx, SOUND_ay,SOUND_by,Pression, it )
  !Subroutine permettant l echange des donnees communes entre processeurs
   use precisions
   use GlobalParam
   IMPLICIT NONE
   !Variables entree-sortie
    REAL (KIND = DP) :: Prim(Nv_Prim,0:Nx+1,0:Ny+1)
    REAL (KIND = DP) :: SOUND_ax(0:Nx+1,0:Ny+1), SOUND_bx(0:Nx+1,0:Ny+1), SOUND_ay(0:Nx+1,0:Ny+1),SOUND_by(0:Nx+1,0:Ny+1)
    REAL (KIND = DP) :: Ein(0:Nx+1,0:Ny+1), Pression(0:Nx+1,0:Ny+1)
   !Variables locales
   REAL(KIND=DP),ALLOCATABLE ::  tampon_env_h( :),tampon_env_u( :),tampon_env_v( :),tampon_env_p11( :),tampon_env_p12( :),tampon_env_p22( :)
   REAL(KIND=DP),ALLOCATABLE ::  tampon_rec_h( :),tampon_rec_u( :),tampon_rec_v( :),tampon_rec_p11( :),tampon_rec_p12( :),tampon_rec_p22( :)
   INTEGER :: i, it

 ALLOCATE(tampon_env_h(1:Ny), tampon_env_u(1:Ny),tampon_env_v(1:Ny),tampon_env_p11(1:Ny),tampon_env_p12(1:Ny),tampon_env_p22(1:Ny))
 ALLOCATE(tampon_rec_h(1:Ny), tampon_rec_u(1:Ny),tampon_rec_v(1:Ny),tampon_rec_p11(1:Ny),tampon_rec_p12(1:Ny),tampon_rec_p22(1:Ny))

  tag = 0;  tampon_env_h = 0; tampon_env_u= 0; tampon_env_v = 0; tampon_env_p11 = 0; tampon_env_p12 = 0; tampon_env_p22 = 0; 
  tampon_rec_h = 0; tampon_rec_u= 0; tampon_rec_v = 0; tampon_rec_p11 = 0; tampon_rec_p12 = 0; tampon_rec_p22 = 0; 
  nombre = Ny
  !Communication avec le proc de rang inferieur

  IF(rang.NE.0) THEN
    
      tampon_env_h(1:Ny) = Prim(H_pv, 1, 1:Ny) 
      tampon_env_u(1:Ny) = Prim(u_pv, 1, 1:Ny)
      tampon_env_v(1:Ny) = Prim(v_pv, 1, 1:Ny)
      tampon_env_p11(1:Ny) = Prim(p11_pv, 1, 1:Ny)
      tampon_env_p12(1:Ny) = Prim(p12_pv, 1, 1:Ny)
      tampon_env_p22(1:Ny) = Prim(p22_pv, 1, 1:Ny) 

     !Envoi vers proc de rang inferieur
      dest   = rang-1
   
      CALL MPI_SEND(tampon_env_h,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_u,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_v,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_p11,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_p12,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_p22,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      

       !Reception du proc de rang inferieur
      source = rang-1
      
      CALL MPI_RECV(tampon_rec_h,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_u,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_v,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_p11,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_p12,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_p22,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
   
       Prim(H_pv, 0, 1:Ny) = tampon_rec_h( 1:Ny) 
       Prim(u_pv, 0, 1:Ny) = tampon_rec_u( 1:Ny) 
       Prim(v_pv, 0, 1:Ny) = tampon_rec_v( 1:Ny) 
       Prim(p11_pv, 0, 1:Ny) = tampon_rec_p11( 1:Ny) 
       Prim(p12_pv, 0, 1:Ny) = tampon_rec_p12( 1:Ny) 
       Prim(p22_pv, 0, 1:Ny) = tampon_rec_p22( 1:Ny) 

        SOUND_ax(0, 1:Ny)  = DSQRT( g*Prim(H_pv,0, 1:Ny) + 3.d0*Prim(P11_pv,0, 1:Ny))
        SOUND_bx(0, 1:Ny)  = DSQRT(Prim(P11_pv,0, 1:Ny))  
        SOUND_ay(0, 1:Ny)  = DSQRT( g*Prim(H_pv,0, 1:Ny) + 3.d0*Prim(P22_pv,0, 1:Ny))
        SOUND_by(0, 1:Ny)  = DSQRT(Prim(P22_pv,0, 1:Ny))  
        Pression(0, 1:Ny)  = g*(Prim(H_pv,0, 1:Ny))**2.d0/2.d0 + Prim(P11_pv,0, 1:Ny)*Prim(H_pv,0, 1:Ny)
  ENDIF

  !Communication avec le proc de rang superieur
  !--------------------------------------------

  IF (rang.NE.ncpu-1) THEN
    
      tampon_env_h(1:Ny) = Prim(H_pv, Nx, 1:Ny) 
      tampon_env_u(1:Ny) = Prim(u_pv, Nx, 1:Ny)
      tampon_env_v(1:Ny) = Prim(v_pv, Nx, 1:Ny)
      tampon_env_p11(1:Ny) = Prim(p11_pv, Nx, 1:Ny)
      tampon_env_p12(1:Ny) = Prim(p12_pv, Nx, 1:Ny)
      tampon_env_p22(1:Ny) = Prim(p22_pv, Nx, 1:Ny) 
    
      !Envoi vers proc de rang superieur
      dest   = rang+1

      CALL MPI_SEND(tampon_env_h,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_u,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_v,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_p11,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_p12,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
      CALL MPI_SEND(tampon_env_p22,Ny,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,code)
     

      !Reception du proc de rang superieur
      source = rang+1
      
      CALL MPI_RECV(tampon_rec_h,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_u,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_v,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_p11,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_p12,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
      CALL MPI_RECV(tampon_rec_p22,Ny,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,code)
     

       Prim(H_pv, Nx+1, 1:Ny) = tampon_rec_h( 1:Ny) 
       Prim(u_pv, Nx+1, 1:Ny) = tampon_rec_u( 1:Ny) 
       Prim(v_pv, Nx+1, 1:Ny) = tampon_rec_v( 1:Ny) 
       Prim(p11_pv, Nx+1, 1:Ny) = tampon_rec_p11( 1:Ny) 
       Prim(p12_pv, Nx+1, 1:Ny) = tampon_rec_p12( 1:Ny) 
       Prim(p22_pv, Nx+1, 1:Ny) = tampon_rec_p22( 1:Ny) 


        SOUND_ax( Nx+1, 1:Ny)  = DSQRT( g*Prim(H_pv, Nx+1, 1:Ny) + 3.d0*Prim(P11_pv, Nx+1, 1:Ny))
        SOUND_bx( Nx+1, 1:Ny)  = DSQRT(Prim(P11_pv, Nx+1, 1:Ny))  
        SOUND_ay( Nx+1, 1:Ny)  = DSQRT( g*Prim(H_pv, Nx+1, 1:Ny) + 3.d0*Prim(P22_pv, Nx+1, 1:Ny))
        SOUND_by( Nx+1, 1:Ny)  = DSQRT(Prim(P22_pv, Nx+1, 1:Ny))  
        Pression( Nx+1, 1:Ny)  = g*(Prim(H_pv, Nx+1, 1:Ny))**2.d0/2.d0 + Prim(P11_pv, Nx+1, 1:Ny)*Prim(H_pv, Nx+1, 1:Ny)

  ENDIF


  DEALLOCATE(tampon_env_h, tampon_env_u,tampon_env_v,tampon_env_p11,tampon_env_p12,tampon_env_p22)
  DEALLOCATE(tampon_rec_h, tampon_rec_u,tampon_rec_v,tampon_rec_p11,tampon_rec_p12,tampon_rec_p22)

  RETURN
  END SUBROUTINE COMMUNICATION


  

  END PROGRAM code2D


   
