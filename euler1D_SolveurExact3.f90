PROGRAM euler1D_SolveurExact
IMPLICIT NONE

INTEGER :: N, i, iter
INTEGER, PARAMETER :: DP=kind(1.d0)
REAL (kind=DP) :: t1_cpu, t2_cpu, L, xd, tfin, dt, pg, pd, rog, rod, ug, ud, x0, dx
REAL (kind=DP) :: temps, gama, eps
REAL (kind=DP) :: ul, ur, Pl, Pr, rol, ror, cl, cr, umax, cfl
REAL (kind=DP), ALLOCATABLE :: u(:), P(:), ro(:), e(:), c(:), x(:)
REAL (kind=DP), ALLOCATABLE :: Cons(:,:), Flux(:,:)

CALL cpu_time(t1_cpu)

! ----------------------------------------------------------------------------------
! Lecture donnees
OPEN(unit=10,file='data.inp',status='old')
	READ(10,*) N 	! Nombre de mailles
	READ(10,*) L 	! Longeur du tube
	READ(10,*) xd 	! Position de la separation 
	READ(10,*) tfin ! Temps final
	READ(10,*) pg 	! Pression gauche membrane
	READ(10,*) pd 	! Pression droite membrane
	READ(10,*) rog 	! Masse volumique gauche membrane
	READ(10,*) rod 	! Masse volumique droite membrane
	READ(10,*) ug 	! Vitesse gauche membrane
	READ(10,*) ud 	! Vitesse droite membrane
	READ(10,*) eps 	! Critere de convergence
CLOSE(10)

! ----------------------------------------------------------------------------------
! Initialisation
ALLOCATE(u(0:N+1), P(0:N+1), ro(0:N+1), e(0:N+1), c(0:N+1), x(1:N))
ALLOCATE(Cons(1:N,3), Flux(0:N,3))

x0=0.d0
dx=L/dfloat(N)
temps=0.d0

! Air
gama = 1.4d0

! ----------------------------------------------------------------------------------
! Initialisations : Variables primitives
DO i=1,N
	x(i)=x0+dfloat(i-1)*dx+0.5d0*dx
	IF (x(i)<=xd) THEN
		P(i)=pg
		u(i)=ug
		ro(i)=rog
		e(i)=pg/((gama-1.d0)*rog)
		c(i)=dsqrt(gama*pg/rog)
	ELSE
		P(i)=pd
		u(i)=ud
		ro(i)=rod
		e(i)=pd/((gama-1.d0)*rod)
		c(i)=dsqrt(gama*pd/rod)
	END IF
END DO

! ----------------------------------------------------------------------------------
! Initialisations : Variables conservatives
DO i=1,N
	Cons(i,1)=ro(i)
	Cons(i,2)=ro(i)*u(i)
	Cons(i,3)=ro(i)*(e(i)+u(i)*u(i)/2.d0)
END DO

Flux=0.d0

! ----------------------------------------------------------------------------------
! ----------------------------------------
! Boucle sur le temps
iter=0
DO WHILE (temps<tfin)

! ----------------------------------------------------------------------------------
! CL : Absorption
u(0)=u(1)
P(0)=P(1)
ro(0)=ro(1)
 c(0)=c(1)

u(N+1)=u(N)
P(N+1)=P(N)
ro(N+1)=ro(N)
 c(N+1)=c(N)

! ----------------------------------------------------------------------------------
! Etats gauche et droite
umax=0.d0
DO i=0,N
	ul=u(i)
	Pl=P(i)
	rol=ro(i)
	cl=c(i)

	ur=u(i+1)
	Pr=P(i+1)
	ror=ro(i+1)
	cr=c(i+1)

! ----------------------------------------------------------------------------------
! Appel Solveur Exact
	CALL SolveurExact(eps, gama, umax, ul, ur, Pl, Pr, rol, ror, cl, cr, Flux(i,1), Flux(i,2), Flux(i,3))
END DO

! ----------------------------------------------------------------------------------
! Mise a jour du pas de temps
cfl=0.9d0
dt=cfl*dx/umax

! ----------------------------------------------------------------------------------
! Schema de Godunov : Calcul nouvelles variables conservatives
DO i=1,N
	Cons(i,:)=Cons(i,:)-dt/dx*(Flux(i,:)-Flux(i-1,:))
END DO

! ----------------------------------------------------------------------------------
! Calcul nouvelles variables primitives
DO i=1,N
	ro(i)=Cons(i,1)
	u(i)=Cons(i,2)/Cons(i,1)
	IF (dabs(u(i))<1.d-9) THEN
		u(i)=0.d0
		Cons(i,2)=0.d0
	END IF
	e(i)=Cons(i,3)/Cons(i,1)-u(i)*u(i)/2.d0
	P(i)=(gama-1)*ro(i)*e(i)
	c(i)=dsqrt(gama*P(i)/ro(i))
END DO

temps=temps+dt
iter=iter+1
! PRINT*, 'iter =', iter
END DO
! Fin boucle sur le temps
! ----------------------------------------
! ----------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
! Ecriture des donnees
OPEN(10,file='res')
DO i=1,N
	WRITE(10,11) x(i), u(i), P(i)/1.d5, ro(i)
END DO
CLOSE(10)
11 format(5(e11.4,1x))

DEALLOCATE(u, P, ro, e, c, x, Cons, Flux)

CALL cpu_time(t2_cpu)
PRINT*, 'L execution du programme a pris : ', t2_cpu-t1_cpu
PRINT*, 'En ', iter, 'iterations.'

END PROGRAM euler1D_SolveurExact

! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

SUBROUTINE SolveurExact(eps, gama, umax, ul, ur, Pl, Pr, rol, ror, cl, cr, F1, F2, F3)
IMPLICIT NONE

INTEGER :: iteration
INTEGER, PARAMETER :: DP=kind(1.d0)
REAL (kind=DP) :: ul, ur, Pl, Pr, rol, ror, cl, cr
REAL (kind=DP) :: gama, z, ulmoy, urmoy, u0star
REAL (kind=DP) :: ustar, eps, Plstar, Prstar, Plstarprim, Prstarprim, Constl, Constr
REAL (kind=DP) :: Machl, Machr, clstar, crstar, Pstar
REAL (kind=DP) :: rolstar, rorstar, shockspeedl, shockspeedr
REAL (kind=DP) :: rarefactionheadl, rarefactionheadr, rarefactiontaill, rarefactiontailr
REAL (kind=DP) :: u, P, ro
REAL (kind=DP) :: F1, F2, F3, Sl, Sr, umax

! ----------------------------------------------------------------------------------
! Initialisation
z=cr/cl*(Pl/Pr)**((gama-1.d0)/2.d0/gama)
ulmoy=ul+2.d0*cl/(gama-1.d0)
urmoy=ur-2.d0*cr/(gama-1.d0)

u0star=(ulmoy*z+urmoy)/(1.d0+z)

! ----------------------------------------------------------------------------------
! Iteration pour trouver ustar
Plstar=0.d0
Prstar=1.d0
Constl=gama*Pl/cl
Constr=gama*Pr/cr
iteration=1
DO WHILE (abs(1.d0-Plstar/Prstar)>eps)

! Calcul de ustar
	IF (iteration==1) THEN
		ustar=u0star
	ELSE
		ustar=ustar-(Plstar-Prstar)/(Plstarprim-Prstarprim)
	END IF

! Choc face à gauche
	IF (ustar<=ul) THEN
		Machl=(gama+1.d0)/4.d0*(ustar-ul)/cl- &
		& dsqrt(1.d0+((gama+1.d0)/4.d0*(ustar-ul)/cl)*((gama+1.d0)/4.d0*(ustar-ul)/cl))
		Plstar=Pl+Constl*(ustar-ul)*Machl
		Plstarprim=2.d0*Constl*Machl*Machl*Machl/(1.d0+Machl*Machl)
	END IF

! Choc face à droite
	IF (ustar>=ur) THEN
		Machr=(gama+1.d0)/4.d0*(ustar-ur)/cr+ &
		& dsqrt(1.d0+((gama+1.d0)/4.d0*(ustar-ur)/cr)*((gama+1.d0)/4.d0*(ustar-ur)/cr))
		Prstar=Pr+Constr*(ustar-ur)*Machr
		Prstarprim=2.d0*Constr*Machr*Machr*Machr/(1.d0+Machr*Machr)
	END IF

! Detentes face à gauche
	IF (ustar>ul) THEN
		clstar=cl-(gama-1.d0)/2.d0*(ustar-ul)
		Plstar=Pl*(clstar/cl)**(2.d0*gama/(gama-1.d0))
		Plstarprim=-gama*Plstar/clstar
	END IF

! Detentes face à droite
	IF (ustar<ur) THEN
		crstar=cr+(gama-1.d0)/2.d0*(ustar-ur)
		Prstar=Pr*(crstar/cr)**(2.d0*gama/(gama-1.d0))
		Prstarprim=gama*Prstar/crstar
	END IF

iteration=iteration+1
END DO

! ----------------------------------------------------------------------------------
! Termes manquant à calculer
Pstar=Plstar

! Choc face à gauche
IF (ustar<=ul) THEN
	shockspeedl=ul+cl*Machl
	rolstar=rol*(ul-shockspeedl)/(ustar-shockspeedl)
END IF

! Choc face à droite
IF (ustar>=ur) THEN
	shockspeedr=ur+cr*Machr
	rorstar=ror*(ur-shockspeedr)/(ustar-shockspeedr)
END IF

! Detentes face à gauche
IF (ustar>ul) THEN
	rarefactionheadl=ul-cl
	rarefactiontaill=ustar-clstar
	rolstar=rol*(PStar/Pl)**(1.d0/gama)
END IF

! Detentes face à droite
IF (ustar<ur) THEN
	rarefactionheadr=ur+cr
	rarefactiontailr=ustar+crstar
	rorstar=ror*(PStar/Pr)**(1.d0/gama)
END IF


! ----------------------------------------------------------------------------------
! Determination des variables primitives

IF (ustar>=0.d0) THEN
	IF (Pstar>=Pl) THEN
		IF (shockspeedl>=0.d0) THEN
			u=ul
			P=Pl
			ro=rol
		ELSE
			u=ustar
			P=Pstar
			ro=rolstar
		END IF
	ELSE
		IF (rarefactionheadl>=0.d0) THEN
			u=ul
			P=Pl
			ro=rol
		ELSE
			IF (rarefactiontaill>=0.d0) THEN
				u=2.d0/(gama+1)*(cl+(gama-1)/2.d0*ul)
				P=Pl*(2.d0/(gama+1)+(gama-1)/(gama+1)/cl*ul)**(2.d0*gama/(gama-1))
				ro=rol*(2.d0/(gama+1)+(gama-1)/(gama+1)/cl*ul)**(2.d0/(gama-1))
			ELSE
				u=ustar
				P=Pstar
				ro=rolstar
			END IF
		END IF
	END IF
ELSE
	IF (Pstar>=Pr) THEN
		IF (shockspeedr>=0.d0) THEN
			u=ustar
			P=Pstar
			ro=rorstar
		ELSE
			u=ur
			P=Pr
			ro=ror
		END IF
	ELSE
		IF (rarefactionheadr>=0.d0) THEN
			IF (rarefactiontailr>=0.d0) THEN
				u=ustar
				P=Pstar
				ro=rorstar
			ELSE
				u=2.d0/(gama+1)*(-cr+(gama-1)/2.d0*ur)
				P=Pr*(2.d0/(gama+1)-(gama-1)/(gama+1)/cr*ur)**(2.d0*gama/(gama-1))
				ro=ror*(2.d0/(gama+1)-(gama-1)/(gama+1)/cr*ur)**(2.d0/(gama-1))
			END IF
		ELSE
			u=ur
			P=Pr
			ro=ror
		END IF
	END IF
END IF

! ----------------------------------------------------------------------------------
! Calcul et retour des Flux et umax
F1=ro*u
F2=ro*u*u+P
F3=u*P*gama/(gama-1.d0)+ro*u*u*u/2.d0

Sl=dmin1(ul-cl, ur-cr)
Sr=dmax1(ul+cl, ur+cr)
umax=dmax1(umax,dabs(Sl),dabs(Sr))

RETURN
END SUBROUTINE