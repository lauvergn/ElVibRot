!======================================================================
!
!      Calculation of Tdef Tcor Trot at Q
!      Tdef = Tdef2 * d2./dQ1dQ2 + Tdef1 * d./dQ1 + vep
!      Tcor = Tcor2 * d./dQ1*Jx  + Tcor1 * Jx
!      Trot = Trot  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!======================================================================
      SUBROUTINE calc_f2_f1Q_ana(Qdyn,                                  &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind) ::  Qdyn(mole%nb_var)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!--- variables for Tnum --------------------------------------------------
!
!      Tdef = Tdef2(1,2) * d2./dQ1dQ2 + Tdef1(1) * d./dQ1 + vep
!      Tcor = Tcor2(1,x) * d./dQ1*Jx  + Tcor1(x) * Jx
!      Trot = Trot(x,y)  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!      nrho: type of normalization
!            2 => clever choice ( Q=R =>1.; Q=val => sin(val); Q=cos(val) =>1. Q=dih =>1.)
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
      real (kind=Rkind) :: Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
      real (kind=Rkind) :: Trot(3,3)

      real (kind=Rkind) :: rho

!-------------------------------------------------------------------------

      real (kind=Rkind), parameter :: mH   = 1837.1526464003414_Rkind ! mH ! Tnum
      real (kind=Rkind), parameter :: mH2  = TWO*mH ! mH2

      integer       :: i,iQdyn

      real (kind=Rkind) :: R,th,phi,c,s
      real (kind=Rkind) :: BH2

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING calc_f2_f1Q_ana'
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qdyn',Qdyn
       END IF
!-----------------------------------------------------------

      Tdef2(:,:) = ZERO
      Tdef1(:)   = ZERO
      vep        = ZERO
      rho        = ZERO
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO

      R   = Qdyn(7)+Qdyn(13)

      th  = Qdyn(8)
      phi = Qdyn(9)

      BH2 = ONE/(mH*R**2)

      s = sin(th)
      c = cos(th)

      IF (debug) THEN
        write(out_unitp,*) 'R,th,phi',R,th,phi
        write(out_unitp,*) 'BH2',BH2
      END IF

      rho = s

      Tdef2(1,1)   = -HALF/(mH*TWO)
      Tdef2(2,2)   = -HALF/(mH*TWO)
      Tdef2(3,3)   = -HALF/(mH*TWO)

      Tdef2(4,4)   = -BH2
      Tdef2(5,5)   = -BH2/(s*s)

      Tdef1(4)     = -BH2 * c/s

      DO i=1,3
        Trot(i,i) = -HALF
      END DO


!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN
        write(out_unitp,*) 'END calc_f2_f1Q_ana'
      END IF
!-----------------------------------------------------------


      end subroutine calc_f2_f1Q_ana
