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
      SUBROUTINE calc_f2_f1Q_ana(Qsym0,                                 &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
      USE mod_system
      USE mod_Tnum
      USE mod_Constant
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind) :: Qsym0(mole%nb_var)
      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
      real (kind=Rkind) :: Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
      real (kind=Rkind) :: Trot(3,3)

      real (kind=Rkind) :: rho


      integer :: nb_Func,ndimFunc
      real (kind=Rkind), allocatable :: d0Func(:),d1Func(:,:)
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING calc_f2_f1Q_ana'
         write(out_unitp,*) 'mole%nb_act',mole%nb_act
         write(out_unitp,*) 'mole%nb_var',mole%nb_var
         write(out_unitp,*) 'Qsym0',Qsym0
         IF (debug) THEN
           write(out_unitp,*)
           CALL Write_CoordType(mole)
           write(out_unitp,*)
         END IF
         write(out_unitp,*) 'JJ',para_Tnum%JJ
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

      CALL get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      IF (ndimFunc /= 1 .OR. nb_Func /= 5) STOP 'pb with QML'

      CALL get_Qmodel_d0d1Func(d0Func,d1Func,Qsym0(1:1),nb_Func,ndimFunc)

      Tdef2(:,:) = -HALF * d0Func(5)
      Tdef1(:)   = -HALF * d1Func(1,5)
      vep        = ZERO
      rho        = ONE
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO

!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN
        write(out_unitp,*) 'END calc_f2_f1Q_ana'
      END IF
!-----------------------------------------------------------

      end subroutine calc_f2_f1Q_ana
