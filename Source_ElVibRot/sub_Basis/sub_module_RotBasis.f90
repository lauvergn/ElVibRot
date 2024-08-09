!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
! MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!    Copyright 2015 David Lauvergnat [1]
!      with contributions of
!        Josep Maria Luis (optimization) [2]
!        Ahai Chen (MPI) [1,4]
!        Lucien Dupuy (CRP) [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Institut de Química Computacional and Departament de Química,
!        Universitat de Girona, Catalonia, Spain
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================
      MODULE mod_RotBasis_Param
      USE EVR_system_m
      use mod_nDindex
      IMPLICIT NONE

      PRIVATE

        TYPE RotBasis_Param
          integer                       :: Jrot                        = -1 !  J value
          integer                       :: nb_Rot                      = 0  !  size of the operator (2*Jrot+1)
          integer                       :: nb_term                     = 0
          integer                       :: tab_der_TO_iterm(-3:0,-3:0) = 0
                                                    ! i1 or i2 =-3,-2,-1   => Jz, Jy, Jx
                                                    ! ex: -2,-1            => JyJx+JxJy operator
                                                    ! ex: -2, 0 or 0,-2    => Jy operator

          integer,          allocatable :: tab_iterm_TO_der(:,:)            !  ...(2,nb_term)
          real(kind=Rkind), allocatable :: tab_RotOp(:,:,:)                 ! tab_RotOp(nb_Rot,nb_Rot,0:nb_term)
                                                                            ! tab_RotOp(:,:,0) is not used but it needs when tab_der_TO_iterm(0,0)=0
        CONTAINS
          PROCEDURE, PRIVATE, PASS(RotBasis_Para1) :: RotBasis_Param2TORotBasis_Param1
          GENERIC,   PUBLIC  :: assignment(=) => RotBasis_Param2TORotBasis_Param1
        END TYPE RotBasis_Param

      PUBLIC :: RotBasis_Param, alloc_RotBasis_Param,                   &
                dealloc_RotBasis_Param, Init_RotBasis_Param,            &
                Write_RotBasis_Param,RotBasis_Param2TORotBasis_Param1

      CONTAINS

      ! symmetrized version => JxJy+JyJx ...
      SUBROUTINE alloc_RotBasis_Param(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2,iterm
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_RotBasis_Param'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot < 0 .OR. RotBasis_Para%nb_Rot > 0) RETURN
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      RotBasis_Para%Jrot    = Jrot
      RotBasis_Para%nb_Rot  = Jrot+Jrot+1

      RotBasis_Para%nb_term = 0
      DO J1=-3,-1
      DO J2=J1,-1
         RotBasis_Para%nb_term   = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,J2) = RotBasis_Para%nb_term
         RotBasis_Para%tab_der_TO_iterm(J2,J1) = RotBasis_Para%nb_term
      END DO
      END DO
      DO J1=-3,-1
         RotBasis_Para%nb_term  = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,0) = RotBasis_Para%nb_term
         RotBasis_Para%tab_der_TO_iterm(0,J1) = RotBasis_Para%nb_term
      END DO
      RotBasis_Para%tab_der_TO_iterm(0,0) = 0


      CALL alloc_NParray(RotBasis_Para%tab_iterm_TO_der,                &
                                         [2,RotBasis_Para%nb_term], &
                        "RotBasis_Para%tab_iterm_TO_der",name_sub, [1,0])

      DO J1=-3,0
      DO J2=J1,0
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
         RotBasis_Para%tab_iterm_TO_der(:,iterm) = [J1,J2]
      END DO
      END DO

      CALL alloc_NParray(RotBasis_Para%tab_RotOp,                       &
 [RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot,RotBasis_Para%nb_term], &
                        "RotBasis_Para%tab_RotOp",name_sub, [1,1,0])
      RotBasis_Para%tab_RotOp(:,:,:) = ZERO



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE alloc_RotBasis_Param

      SUBROUTINE alloc_RotBasis_Param_old(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_RotBasis_Param_old'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot < 0 .OR. RotBasis_Para%nb_Rot > 0) RETURN
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      RotBasis_Para%Jrot    = Jrot
      RotBasis_Para%nb_Rot  = Jrot+Jrot+1

      RotBasis_Para%nb_term = 0
      DO J1=-3,-1
      DO J2=-3,-1
         RotBasis_Para%nb_term   = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,J2) = RotBasis_Para%nb_term
      END DO
      END DO
      DO J1=-3,-1
         RotBasis_Para%nb_term  = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,0) = RotBasis_Para%nb_term
         RotBasis_Para%tab_der_TO_iterm(0,J1) = RotBasis_Para%nb_term
      END DO
      RotBasis_Para%tab_der_TO_iterm(0,0) = 0


      CALL alloc_NParray(RotBasis_Para%tab_RotOp,                       &
 [RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot,RotBasis_Para%nb_term], &
                        "RotBasis_Para%tab_RotOp",name_sub, [1,1,0])
      RotBasis_Para%tab_RotOp(:,:,:) = ZERO



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE alloc_RotBasis_Param_old


      SUBROUTINE dealloc_RotBasis_Param(RotBasis_Para)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='dealloc_RotBasis_Param'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      RotBasis_Para%Jrot    = -1
      RotBasis_Para%nb_Rot  = 0
      RotBasis_Para%nb_term = 0

      RotBasis_Para%tab_der_TO_iterm(:,:) = 0

      IF (allocated(RotBasis_Para%tab_RotOp)) THEN
        CALL dealloc_NParray(RotBasis_Para%tab_iterm_TO_der,            &
                            "RotBasis_Para%tab_iterm_TO_der",name_sub)
      END IF

      IF (allocated(RotBasis_Para%tab_RotOp)) THEN
        CALL dealloc_NParray(RotBasis_Para%tab_RotOp,                   &
                            "RotBasis_Para%tab_RotOp",name_sub)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE dealloc_RotBasis_Param

      SUBROUTINE Init_RotBasis_Param(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2,K,iterm,iterm1,iterm2
      complex (kind=Rkind) :: s2
      complex (kind=Rkind) :: is2

      complex(kind=Rkind), allocatable :: PWang(:,:)     ! PWang(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: PWangDag(:,:)  ! PWangDag(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: test_inv(:,:)  ! PWangDag(nb_Rot,nb_Rot)

      complex(kind=Rkind), allocatable :: Jx(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jy(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jz(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Ji(:,:,:)  ! Ji(-J:J,-J:J,-3:-1)

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Init_RotBasis_Param'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot < 0) RETURN
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      s2  = cmplx( sqrt(HALF), ZERO,kind=Rkind )
      is2 = cmplx( ZERO, sqrt(HALF),kind=Rkind )


      CALL alloc_RotBasis_Param(RotBasis_Para,Jrot)

!---------------------------------------------------------------------
      ! -1- First the Wang basis
      CALL alloc_NParray(PWang,[Jrot,RotBasis_Para%nb_Rot],         &
                        'PWang',name_sub,[-Jrot,1] )

      CALL alloc_NParray(PWangDag,[RotBasis_Para%nb_Rot,Jrot],      &
                        'PWangDag',name_sub,[1,-Jrot] )


      PWangDag(:,:) = ZERO
      PWangDag(1,0) = ONE
      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) =  s2
        ELSE
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) = -s2
        END IF
      END DO

      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) = -is2
        ELSE
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) =  is2
        END IF
      END DO
      PWang(:,:) = Transpose(conjg(PWangDag))

      IF (debug) THEN
        write(out_unit,*) 'PWangDag (in line)'
        CALL Write_Mat(PWangDag,out_unit,5)

        write(out_unit,*) 'PWang (in column)'
        CALL Write_Mat(PWang,out_unit,5)
      END IF

      !CALL alloc_NParray(test_inv,[RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot],         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(PWang,PWangDag)
      !write(out_unit,*) 'PWang . PWangDag'
      !CALL Write_Mat(test_inv,out_unit,5)
      !test_inv = matmul(PWangDag,PWang)
      !write(out_unit,*) 'PWangDag . PWang'
      !CALL Write_Mat(test_inv,out_unit,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)


!---------------------------------------------------------------------
      ! -2- Jx,Jy,Jz in the I JKM > basis
      CALL alloc_NParray(Ji,[Jrot,Jrot,-1],'Ji',name_sub,[-Jrot,-Jrot,-3] )

      CALL alloc_NParray(Jx,[Jrot,Jrot],'Jx',name_sub,[-Jrot,-Jrot] )
      Jx(:,:) = ZERO
      CALL alloc_NParray(Jy,[Jrot,Jrot],'Jy',name_sub,[-Jrot,-Jrot] )
      Jy(:,:) = ZERO
      CALL alloc_NParray(Jz,[Jrot,Jrot],'Jz',name_sub,[-Jrot,-Jrot] )
      Jz(:,:) = ZERO


      DO K=-Jrot,Jrot
        Jz(K,K) = cmplx(K,0,kind=Rkind)

        IF (K < Jrot) THEN
          Jx(K,K+1) =      cmplx(funcP_JK(Jrot,K),kind=Rkind)
          Jy(K,K+1) = -EYE*cmplx(funcP_JK(Jrot,K),kind=Rkind)
        END IF

        IF (K > -Jrot) THEN
          Jx(K,K-1) =      cmplx(funcM_JK(Jrot,K),kind=Rkind)
          Jy(K,K-1) =  EYE*cmplx(funcM_JK(Jrot,K),kind=Rkind)
        END IF

      END DO

      !write(out_unit,*)
      !write(out_unit,*) 'Re Jx'
      !CALL Write_Mat(real(Jx),out_unit,5)
      !write(out_unit,*) 'Jy'
      !CALL Write_Mat(aimag(Jy),out_unit,5)
      !write(out_unit,*) 'Jz'
      !CALL Write_Mat(real(Jz),out_unit,5)

!---------------------------------------------------------------------
      ! -3- Jx,Jy,Jz in the Wang basis
      Ji(:,:,-1) = matmul(PWangDag,matmul(Jx,PWang))
      Ji(:,:,-2) = matmul(PWangDag,matmul(Jy,PWang))
      Ji(:,:,-3) = matmul(PWangDag,matmul(Jz,PWang))

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'Im(Jx), SumRe(Jx)',sum(abs(real(Ji(:,:,-1))))
        CALL Write_Mat(aimag(Ji(:,:,-1)),out_unit,5)
        write(out_unit,*) 'Im(Jy), SumRe(Jy)',sum(abs(real(Ji(:,:,-2))))
        CALL Write_Mat(aimag(Ji(:,:,-2)),out_unit,5)
        write(out_unit,*) 'Im(Jz), SumRe(Jz)',sum(abs(real(Ji(:,:,-3))))
        CALL Write_Mat(aimag(Ji(:,:,-3)),out_unit,5)
      END IF
      !CALL alloc_NParray(test_inv,[RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot],         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(Jx,Jx)+matmul(Jy,Jy)+matmul(Jz,Jz)
      !write(out_unit,*) 'J^2'
      !CALL Write_Mat(test_inv,out_unit,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)

!---------------------------------------------------------------------
      ! -4- Transfert Jx,Jy .... (Jx*Jy+JyJx) ... in RotBasis_Para%tab_RotOp
      DO J1=-3,-1
        iterm = RotBasis_Para%tab_der_TO_iterm(J1,J1)
        RotBasis_Para%tab_RotOp(:,:,iterm) =                            &
                                     real(matmul(Ji(:,:,J1),Ji(:,:,J1)))
        DO J2=J1+1,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
         RotBasis_Para%tab_RotOp(:,:,iterm) =                           &
                                  real(matmul(Ji(:,:,J1),Ji(:,:,J2))) + &
                                  real(matmul(Ji(:,:,J2),Ji(:,:,J1)))
        END DO
      END DO
      ! A minus sign is added because the Hamitonian has i*Ji and the Wang basis are imaginary
      DO J1=-3,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,0)
         RotBasis_Para%tab_RotOp(:,:,iterm) = -aimag(Ji(:,:,J1))
      END DO

      ! Id
      iterm = RotBasis_Para%tab_der_TO_iterm(0,0)
      RotBasis_Para%tab_RotOp(:,:,iterm) = ZERO
      DO K=1,RotBasis_Para%nb_Rot
        RotBasis_Para%tab_RotOp(K,K,iterm) = ONE
      END DO



!---------------------------------------------------------------------
      CALL dealloc_NParray(PWangDag,'PWangDag',name_sub)
      CALL dealloc_NParray(PWang,'PWang',name_sub)
      CALL dealloc_NParray(Jx,'Jx',name_sub)
      CALL dealloc_NParray(Jy,'Jy',name_sub)
      CALL dealloc_NParray(Jz,'Jz',name_sub)
      CALL dealloc_NParray(Ji,'Ji',name_sub)


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

!STOP 'Init_RotBasis_Param'
!---------------------------------------------------------------------

      END SUBROUTINE Init_RotBasis_Param

      SUBROUTINE Init_RotBasis_Param_old(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2,K,iterm,iterm1,iterm2
      complex (kind=Rkind) :: s2
      complex (kind=Rkind) :: is2

      complex(kind=Rkind), allocatable :: PWang(:,:)     ! PWang(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: PWangDag(:,:)  ! PWangDag(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: test_inv(:,:)  ! PWangDag(nb_Rot,nb_Rot)

      complex(kind=Rkind), allocatable :: Jx(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jy(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jz(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Ji(:,:,:)  ! Ji(-J:J,-J:J,-3:-1)

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Init_RotBasis_Param_old'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot <= 0 .OR. RotBasis_Para%nb_Rot > 0) RETURN
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      s2  = cmplx( sqrt(HALF), ZERO,kind=Rkind )
      is2 = cmplx( ZERO, sqrt(HALF),kind=Rkind )


      CALL alloc_RotBasis_Param(RotBasis_Para,Jrot)

!---------------------------------------------------------------------
      ! -1- First the Wang basis
      CALL alloc_NParray(PWang,[Jrot,RotBasis_Para%nb_Rot],         &
                        'PWang',name_sub,[-Jrot,1] )

      CALL alloc_NParray(PWangDag,[RotBasis_Para%nb_Rot,Jrot],      &
                        'PWangDag',name_sub,[1,-Jrot] )


      PWangDag(:,:) = CZERO
      PWangDag(1,0) = CONE
      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) =  s2
        ELSE
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) = -s2
        END IF
      END DO

      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) = -is2
        ELSE
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) =  is2
        END IF
      END DO
      PWang(:,:) = Transpose(conjg(PWangDag))

      IF (debug) THEN
        write(out_unit,*) 'PWangDag (in line)'
        CALL Write_Mat(PWangDag,out_unit,5)

        write(out_unit,*) 'PWang (in column)'
        CALL Write_Mat(PWang,out_unit,5)
      END IF

      !CALL alloc_NParray(test_inv,[RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot],         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(PWang,PWangDag)
      !write(out_unit,*) 'PWang . PWangDag'
      !CALL Write_Mat(test_inv,out_unit,5)
      !test_inv = matmul(PWangDag,PWang)
      !write(out_unit,*) 'PWangDag . PWang'
      !CALL Write_Mat(test_inv,out_unit,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)


!---------------------------------------------------------------------
      ! -2- Jx,Jy,Jz in the I JKM > basis
      CALL alloc_NParray(Ji,[Jrot,Jrot,-1],'Ji',name_sub,[-Jrot,-Jrot,-3] )

      CALL alloc_NParray(Jx,[Jrot,Jrot],'Jx',name_sub,[-Jrot,-Jrot] )
      Jx(:,:) = CZERO
      CALL alloc_NParray(Jy,[Jrot,Jrot],'Jy',name_sub,[-Jrot,-Jrot] )
      Jy(:,:) = CZERO
      CALL alloc_NParray(Jz,[Jrot,Jrot],'Jz',name_sub,[-Jrot,-Jrot] )
      Jz(:,:) = CZERO


      DO K=-Jrot,Jrot
        Jz(K,K) = cmplx(K,0,kind=Rkind)

        IF (K < Jrot) THEN
          Jx(K,K+1) =      cmplx(funcP_JK(Jrot,K),kind=Rkind)
          Jy(K,K+1) = -EYE*cmplx(funcP_JK(Jrot,K),kind=Rkind)
        END IF

        IF (K > -Jrot) THEN
          Jx(K,K-1) =      cmplx(funcM_JK(Jrot,K),kind=Rkind)
          Jy(K,K-1) =  EYE*cmplx(funcM_JK(Jrot,K),kind=Rkind)
        END IF

      END DO

      !write(out_unit,*)
      !write(out_unit,*) 'Re Jx'
      !CALL Write_Mat(real(Jx),out_unit,5)
      !write(out_unit,*) 'Jy'
      !CALL Write_Mat(aimag(Jy),out_unit,5)
      !write(out_unit,*) 'Jz'
      !CALL Write_Mat(real(Jz),out_unit,5)

!---------------------------------------------------------------------
      ! -3- Jx,Jy,Jz in the Wang basis
      Ji(:,:,-1) = matmul(PWangDag,matmul(Jx,PWang))
      Ji(:,:,-2) = matmul(PWangDag,matmul(Jy,PWang))
      Ji(:,:,-3) = matmul(PWangDag,matmul(Jz,PWang))

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'Im(Jx), SumRe(Jx)',sum(abs(real(Ji(:,:,-1))))
        CALL Write_Mat(aimag(Ji(:,:,-1)),out_unit,5)
        write(out_unit,*) 'Im(Jy), SumRe(Jy)',sum(abs(real(Ji(:,:,-2))))
        CALL Write_Mat(aimag(Ji(:,:,-2)),out_unit,5)
        write(out_unit,*) 'Im(Jz), SumRe(Jz)',sum(abs(real(Ji(:,:,-3))))
        CALL Write_Mat(aimag(Ji(:,:,-3)),out_unit,5)
      END IF
      !CALL alloc_NParray(test_inv,[RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot],         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(Jx,Jx)+matmul(Jy,Jy)+matmul(Jz,Jz)
      !write(out_unit,*) 'J^2'
      !CALL Write_Mat(test_inv,out_unit,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)

!---------------------------------------------------------------------
      ! -4- Transfert Jx,Jy .... Jx*Jy ... in RotBasis_Para%tab_RotOp
      DO J1=-3,-1
      DO J2=-3,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
         RotBasis_Para%tab_RotOp(:,:,iterm) = real(matmul(Ji(:,:,J1),Ji(:,:,J2)))
      END DO
      END DO
      DO J1=-3,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,0)
         RotBasis_Para%tab_RotOp(:,:,iterm) = aimag(Ji(:,:,J1))
      END DO

      ! J^2
      DO J1=-3,-1
        iterm = RotBasis_Para%tab_der_TO_iterm(J1,J1)
        RotBasis_Para%tab_RotOp(:,:,0) = RotBasis_Para%tab_RotOp(:,:,0) + &
            RotBasis_Para%tab_RotOp(:,:,iterm)
      END DO


      test_inv=matmul(Ji(:,:,-1),Ji(:,:,-2))-matmul(Ji(:,:,-2),Ji(:,:,-1))
      !write(out_unit,*) '[Jx,Jy]'
      !CALL Write_Mat(test_inv,out_unit,5)
      test_inv = test_inv + EYE*Ji(:,:,-3)
      write(out_unit,*) '[Jx,Jy] + i Jz = 0 ?',sum(abs(test_inv))
      !CALL Write_Mat(test_inv,out_unit,5)

      test_inv=matmul(Ji(:,:,-3),Ji(:,:,-1))-matmul(Ji(:,:,-1),Ji(:,:,-3))
      !write(out_unit,*) '[Jz,Jx]'
      !CALL Write_Mat(test_inv,out_unit,5)
      test_inv = test_inv+EYE*Ji(:,:,-2)
      write(out_unit,*) '[Jz,Jx] + i Jy = 0 ?',sum(abs(test_inv))
      !CALL Write_Mat(test_inv,out_unit,5)

      test_inv=matmul(Ji(:,:,-2),Ji(:,:,-3))-matmul(Ji(:,:,-3),Ji(:,:,-2))
      !write(out_unit,*) '[Jy,Jz]'
      !CALL Write_Mat(test_inv,out_unit,5)
      test_inv = test_inv+EYE*Ji(:,:,-1)
      write(out_unit,*) '[Jy,Jz] + i Jx = 0 ?',sum(abs(test_inv))
      !CALL Write_Mat(test_inv,out_unit,5)



!---------------------------------------------------------------------
      CALL dealloc_NParray(PWangDag,'PWangDag',name_sub)
      CALL dealloc_NParray(PWang,'PWang',name_sub)
      CALL dealloc_NParray(Jx,'Jx',name_sub)
      CALL dealloc_NParray(Jy,'Jy',name_sub)
      CALL dealloc_NParray(Jz,'Jz',name_sub)
      CALL dealloc_NParray(Ji,'Ji',name_sub)


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

!STOP 'Init_RotBasis_Param'
!---------------------------------------------------------------------

      END SUBROUTINE Init_RotBasis_Param_old

      FUNCTION funcP_JK(J,K)
        real (kind=Rkind) :: funcP_JK
        integer, intent(in) :: J,K

        funcP_JK = HALF*sqrt(real(J*J+J-K*K-K,kind=Rkind))

      END  FUNCTION funcP_JK
      FUNCTION funcM_JK(J,K)
        real (kind=Rkind) :: funcM_JK
        integer, intent(in) :: J,K

        funcM_JK = HALF*sqrt(real(J*J+J-K*K+K,kind=Rkind))

      END  FUNCTION funcM_JK
      SUBROUTINE Write_RotBasis_Param(RotBasis_Para)

      TYPE (RotBasis_Param), intent(in) :: RotBasis_Para


      integer :: J1,J2,iterm
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_RotBasis_Param'
      !logical,parameter :: debug=.FALSE.
      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      write(out_unit,*) 'RotBasis_Para%Jrot   ',RotBasis_Para%Jrot
      write(out_unit,*) 'RotBasis_Para%nb_Rot ',RotBasis_Para%nb_Rot
      write(out_unit,*) 'RotBasis_Para%nb_term',RotBasis_Para%nb_term

      IF (allocated(RotBasis_Para%tab_RotOp)) THEN
        write(out_unit,*) 'J(-1) => -i*Jx, J(-2) => -i*Jy, J(-3) => -i*Jz'
        DO J1=-3,-1
        DO J2=-3,-1
          iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
          write(out_unit,*) 'JiJj op.',J1,J2,' iterm: ',iterm
          CALL Write_Mat(RotBasis_Para%tab_RotOp(:,:,iterm),out_unit,5)
        END DO
        END DO
        DO J1=-3,-1
          iterm = RotBasis_Para%tab_der_TO_iterm(J1,0)
          write(out_unit,*) 'Im(Ji op.)',J1,' iterm: ',iterm
          CALL Write_Mat(RotBasis_Para%tab_RotOp(:,:,iterm),out_unit,5)
        END DO

        write(out_unit,*) 'Id'
        CALL Write_Mat(RotBasis_Para%tab_RotOp(:,:,0),out_unit,5)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Write_RotBasis_Param

      SUBROUTINE RotBasis_Param2TORotBasis_Param1(RotBasis_Para1,       &
                                                         RotBasis_Para2)

      CLASS (RotBasis_Param), intent(inout) :: RotBasis_Para1
      TYPE (RotBasis_Param),  intent(in)    :: RotBasis_Para2

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RotBasis_Param2TORotBasis_Param1'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      CALL alloc_RotBasis_Param(RotBasis_Para1,RotBasis_Para2%Jrot)

      IF (allocated(RotBasis_Para2%tab_RotOp)) THEN
        RotBasis_Para1%tab_RotOp(:,:,:) = RotBasis_Para2%tab_RotOp(:,:,:)
      END IF

      IF (allocated(RotBasis_Para2%tab_iterm_TO_der)) THEN
        RotBasis_Para1%tab_iterm_TO_der(:,:) = RotBasis_Para2%tab_iterm_TO_der(:,:)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_RotBasis_Param(RotBasis_Para2)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE RotBasis_Param2TORotBasis_Param1



      END MODULE mod_RotBasis_Param
