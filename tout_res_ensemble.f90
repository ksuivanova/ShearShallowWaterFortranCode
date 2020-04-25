module precisions
implicit none
integer, parameter :: dp=kind(1.0d0)    
end module precisions

PROGRAM reecrir_res
USE precisions
IMPLICIT NONE
INTEGER :: I, Nx_max, iu, Nx, Ny, entier,nfile, j, k, entierglob
INTEGER,parameter :: ncpu=3
REAL (KIND = DP), ALLOCATABLE :: U(:,:,:), V(:,:,:), H(:,:,:), E(:,:,:), Phi1(:,:,:), P(:,:,:), X(:,:), Y(:,:) !,Q(:,:,:)
REAL (KIND = DP) :: phi2
character*3      :: numrg

OPEN(19,FILE = 'data_2D.inp', status = 'old', action = 'read')
    do i=1,3
      READ(19, *)  
    enddo   
    READ(19,*) Nx_max, Ny

   print*, 'Nx_max =', Nx_max, Ny
    
    do i=1,3
      READ(19, *)  
    enddo   

    READ(19,*) phi2
    close(19)

      OPEN(19, FILE = 'RES000.out')
        READ(19,*) Nx, Ny, nfile
      close(19)
     print*, 'Nx, Ny, nfile =', Nx, Ny, nfile
    
    ALLOCATE( U(1:Nx_max, 1:Ny, nfile), V(1:Nx_max, 1:Ny, nfile), H(1:Nx_max, 1:Ny, nfile), E(1:Nx_max, 1:Ny, nfile),&
    &         Phi1(1:Nx_max, 1:Ny, nfile), X(1:Nx_max,  nfile), Y(1:Ny, nfile), P(1:Nx_max, 1:Ny, nfile)  )  ! , Q(1:Nx_max, 1:Ny, nfile)
    
     entierglob = 0
    do iu = 0, Ncpu-1
      WRITE(numrg,'(i3.3)') iu
      OPEN(19+iu, FILE = 'RES'//numrg//'.out')
      READ(19+iu,*) Nx, Ny, nfile
      print*, 'Nx, Ny, nfile =', Nx, Ny, nfile
      DO j=1,nfile
        
        entier = entierglob
          
          DO i=1,Nx
       
           entier=entier+1
           print*, 'entier 1 =', entier
          Do k = 1, Ny  
          read(19+IU,'(8(E20.13,1X))') X(entier, j), Y(k, j), U(entier, k, j), V(entier, k, j), Phi1(entier, k, j), H(entier,k, j), E(entier, k, j), P(entier, k, j) !, Q(entier, k, j)
     
          enddo
          enddo
         print*, 'entier 3=', entier, 'entierglob = ', entierglob
        READ(19+iu,*)
     
      enddo
       
       entierglob = entier
       print*, 'entier 4 =', entier, 'entierglob =', entierglob
      close(19+iu)
    enddo   
   

      OPEN(50,FILE = 'RES_complet.out')


 DO J = 1,nfile
   DO I = 1,Nx_max
    DO k  = 1, Ny
  

      WRITE(50,24) X(I, j), Y(k, j), U(I, k, j), V(I, k, j), Phi1(I, k,j), H(I, k, j), E(I, k, j), P(I, k, j) !, Q(I, k, j)
      
      ENDDO
   ENDDO
   write(50,*)
enddo

    24 FORMAT(8(E16.8,1X))


END PROGRAM reecrir_res
