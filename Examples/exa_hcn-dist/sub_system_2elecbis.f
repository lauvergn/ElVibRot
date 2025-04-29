c
C================================================================
C    calc_Op : calculation of the potential and scalar operator matrices
c    mat_V(nb_be,nb_be) and Mat_Scal(nb_be,nb_be,nb_ScalOp)
c    nb_be : nb of electronic surfaces
c    Q are the coordinates in active order or dyn order
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Qop,nb_QOp,mole,calc_ScalOp,pot_cplx)

      USE EVR_system_m
      USE mod_Tnum
      IMPLICIT NONE

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

      integer           :: nb_be,nb_ScalOp,nb_QOp
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qop(nb_QOp)

      integer           :: nb_QOp_loc=2

      IF (nb_be == 1 ) THEN
        write(out_unit,*) ' ERROR sub_system for 2 PES'
        STOP
      END IF

      CALL Mat_pot0(mat_V,Qop(1:nb_QOp_loc),nb_be,nb_QOp_loc) 

      IF (pot_cplx) THEN
        CALL Mat_im_pot0(mat_imV,Qop(1:nb_QOp_loc),nb_be,nb_QOp_loc)
      END IF

      IF (calc_ScalOp) THEN
        CALL Mat_ScalarOp(mat_ScalOp,Qop(1:nb_QOp_loc),mole,
     *                    nb_be,nb_ScalOp,nb_QOp_loc)
      END IF

      END SUBROUTINE calcN_op
C================================================================
C     pot0(x) 1 D 2 surfaces
C================================================================
      SUBROUTINE Mat_pot0(mat_V,Qop,nb_be,nb_QOp)
      USE EVR_system_m
      IMPLICIT NONE

      integer           :: nb_be,nb_QOp
      real (kind=Rkind) :: Qop(nb_QOp)
      real (kind=Rkind) :: mat_V(nb_be,nb_be)

      real (kind=Rkind) :: kdiag(nb_QOp)
      real (kind=Rkind) :: Qeq1(nb_QOp)
      real (kind=Rkind) :: Qeq2(nb_QOp)
      real (kind=Rkind) :: DQ(nb_QOp)
      real (kind=Rkind), parameter :: e11  = 0.00_Rkind
      real (kind=Rkind), parameter :: e22  = 0.03_Rkind
      real (kind=Rkind), parameter :: e12  = 0.005_Rkind

       kdiag(:) = [ 0.1_Rkind,0.5_Rkind]
       Qeq1(:)  = [-0.7_Rkind,0.5_Rkind]
       Qeq2(:)  = [-0.6_Rkind,0.4_Rkind]

       !Write(6,*) 'Qop',Qop
       !Write(6,*) 'kdiag',kdiag

       DQ(:) = Qop-Qeq1
       !Write(6,*) 'DQ1',DQ
       mat_V(1,1) = e11 + HALF * dot_product(kdiag*DQ,DQ)

       DQ(:) = Qop-Qeq2
       !Write(6,*) 'DQ2',DQ
       mat_V(2,2) = e22 + HALF * dot_product(kdiag*DQ,DQ)

       mat_V(1,2) = e12
       mat_V(2,1) = e12

       !write(6,*) 'mat_V',mat_V

      END SUBROUTINE Mat_pot0
C================================================================
C    fonction im_pot0(x) imaginary part of pot0
C================================================================
      SUBROUTINE Mat_im_pot0(mat_imV,Qop,nb_be,nb_QOp)
      USE EVR_system_m
      IMPLICIT NONE

      integer           :: nb_be,nb_QOp
      real (kind=Rkind) :: mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: Qop(nb_QOp)

      integer           :: i
      real (kind=Rkind) :: im_pot0 ! imaginary function

      mat_imV = ZERO
      mat_imV(1,1) = -0.001_Rkind

      END SUBROUTINE Mat_im_pot0
c
C================================================================
C    Scalar Opertor (here dipole moment, 3 components (x, y, z)
C================================================================
      SUBROUTINE Mat_ScalarOp(mat_ScalOp,Qop,mole,
     *                        nb_be,nb_ScalOp,nb_QOp)
      USE EVR_system_m
      USE mod_Tnum
      IMPLICIT NONE

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole


      integer           :: nb_be,nb_ScalOp,nb_QOp
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qop(nb_QOp)


       mat_ScalOp(:,:,:) = ZERO

       mat_ScalOp(1,2,1) = ONE
       mat_ScalOp(2,1,1) = ONE

       END SUBROUTINE Mat_ScalarOp
C
C=======================================================================
C    fonction pot_rest(x)
C=======================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE EVR_system_m
      IMPLICIT NONE
      real (kind=Rkind) :: pot_rest


      real (kind=Rkind) :: Qact(1)
      integer           :: nb_inact2n
      real (kind=Rkind) :: Delta_Qact(nb_inact2n)

      !-----------------------------------------------------------------
      pot_rest = ZERO
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      STOP 'The function pot_rest MUST be written'
      !-----------------------------------------------------------------


      END FUNCTION pot_rest

C=======================================================================
C    sub hessian
C=======================================================================
      SUBROUTINE sub_hessian(h)
      USE EVR_system_m
      IMPLICIT NONE

       real (kind=Rkind) :: h

      !-----------------------------------------------------------------
       h = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine sub_hessian MUST be written'
      !-----------------------------------------------------------------


      END SUBROUTINE sub_hessian
C=======================================================================
C     analytical gradient along the path
C=======================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qdyn,mole,deriv,num,step)
      USE EVR_system_m
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the CoordType and Tnum ---------------------------------
      TYPE (CoordType)    :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: 
     *                    d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: step
      logical           :: deriv,num

      real (kind=Rkind) :: Qact1(mole%nb_act1)

      !----- for debuging ----------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='d0d1d2_g'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_var',mole%nb_var
        write(out_unit,*) 'nb_act1',mole%nb_act1
        write(out_unit,*) 'nb_inact22,nb_inact21',
     *                   mole%nb_inact22,mole%nb_inact21
        write(out_unit,*) 'nb_inact2n',mole%nb_inact2n
        write(out_unit,*) 'deriv',deriv
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      Qact1(:) = Qdyn(mole%liste_QactTOQdyn(1:mole%nb_act1))
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      d0g(:)     = ZERO
      d1g(:,:)   = ZERO
      d2g(:,:,:) = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine d0d1d2_g MUST be written'
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'd0g at Qact:',Qact1
        write(out_unit,*) d0g(:)
        write(out_unit,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE d0d1d2_g
C=======================================================================
C     analytical hessian along the path (only d0h is used!!)
C=======================================================================
      SUBROUTINE d0d1d2_h(d0h,d1h,d2h,Qdyn,mole,deriv,num,step)
      USE EVR_system_m
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the CoordType and Tnum ----------------------------------
      TYPE (CoordType)    :: mole

      real (kind=Rkind) :: Qdyn(mole%nb_var)


      real (kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
      real (kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
      real (kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)

      real (kind=Rkind) :: step
      logical           :: deriv,num


      real (kind=Rkind) :: Qact1(mole%nb_act1)

      !----- for debuging ----------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='d0d1d2_h'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_var',mole%nb_var
        write(out_unit,*) 'nb_act1',mole%nb_act1
        write(out_unit,*) 'nb_inact22,nb_inact21',
     *            mole%nb_inact22,mole%nb_inact21
        write(out_unit,*) 'nb_inact2n',mole%nb_inact2n
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      Qact1(:) = Qdyn(mole%liste_QactTOQdyn(1:mole%nb_act1))
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      d0h(:,:) = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine d0d1d2_h MUST be written'
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Qact1',Qact1
        write(out_unit,*) 'd0h at Qact1'
        CALL Write_Mat_MPI(d0h,6,4)
        write(out_unit,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE d0d1d2_h
C=======================================================================
C     analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
C     for the variable iq
C=======================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE EVR_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      integer           :: iq,nb_act
      real (kind=Rkind) :: Qact(nb_act)
      integer           :: nderiv,it
      TYPE (Type_dnS)   :: dnQflex


      ! for debuging ---------------------------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      ! for debuging ---------------------------------------------------


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_act',nb_act
        write(out_unit,*) 'iq',iq
      END IF
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      CALL sub_ZERO_TO_dnS(dnQflex)

      ! Zero order dervivative
      dnQflex%d0 = ZERO
      ! First order dervivatives
      dnQflex%d1(:) = ZERO
      ! Second order dervivatives
      dnQflex%d2(:,:) = ZERO
      ! Third order dervivatives
      dnQflex%d3(:,:,:) = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine calc_dnQflex MUST be written'
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(out_unit,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

       END SUBROUTINE calc_dnQflex
C================================================================
C     analytical derivative (Qopt Qopt' Qopt" Qopt'") calculation
c     for the variable i_Qdyn
C     derivative with respect to Qact1(:) coordinates
C================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_Qdyn,
     *                        d0Qopt,d1Qopt,d2Qopt,d3Qopt,
     *                        Qdyn,mole,nderiv)
      USE EVR_system_m
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the CoordType and Tnum ----------------------------------
      TYPE (CoordType)    :: mole

      integer           :: i_Qdyn

      real (kind=Rkind) :: Qdyn(mole%nb_var)

      integer           :: nderiv

      real (kind=Rkind) :: d0Qopt
      real (kind=Rkind) :: d1Qopt(mole%nb_act1)
      real (kind=Rkind) :: d2Qopt(mole%nb_act1,mole%nb_act1)
      real (kind=Rkind) ::
     *                   d3Qopt(mole%nb_act1,mole%nb_act1,mole%nb_act1)


      !local variables
      real (kind=Rkind) :: Qact1(mole%nb_act1)


      !----- for debuging ----------------------------------------------
      character (len=*), parameter :: name_sub='d0d1d2d3_Qeq'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !----- for debuging ----------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_inact20,nb_act',
     *                     mole%nb_inact20,mole%nb_act
        write(out_unit,*) 'nb_var',mole%nb_var
        write(out_unit,*) 'i_Qdyn',i_Qdyn
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      Qact1(:) = Qdyn(mole%liste_QactTOQdyn(1:mole%nb_act1))
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      d0Qopt        = ZERO
      d1Qopt(:)     = ZERO
      d2Qopt(:,:)   = ZERO
      d3Qopt(:,:,:) = ZERO
      !-----------------------------------------------------------------

      STOP 'The subroutine d0d1d2d3_Qeq MUST be written'


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Qact1 : ',Qact1
        write(out_unit,*) 'd0Qopt : ',d0Qopt
        write(out_unit,*) 'd1Qopt : ',d1Qopt
        write(out_unit,*) 'd2Qopt : ',d2Qopt
        write(out_unit,*) 'd3Qopt : ',d3Qopt
        write(out_unit,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE d0d1d2d3_Qeq