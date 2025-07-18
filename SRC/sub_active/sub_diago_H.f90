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
      SUBROUTINE sub_diago_H(H,E,Vec,n,sym)
      USE EVR_system_m
      USE mod_Constant, ONLY: get_Conv_au_TO_unit
      IMPLICIT NONE

!------ active Matrix H Vec E ------------------------------------

      integer           :: n
      real (kind=Rkind) :: H(n,n)
      real (kind=Rkind) :: Vec(n,n)
      real (kind=Rkind) :: E(n)
      logical           :: sym

!----- divers --------------------------------------------------------
      real (kind=Rkind) :: auTOcm_inv
      real (kind=Rkind) :: A,B

      real (kind=Rkind), allocatable :: Hsave(:,:),r(:)
      logical           :: residual = .FALSE.

      integer :: i


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_diago_H'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'n',n
         flush(out_unit)
       END IF
!---------------------------------------------------------------------------------------
!       H matrix diagonalisation
!---------------------------------------------------------------------------------------
  auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

  IF(keep_MPI) THEN
    IF (residual) THEN
      CALL alloc_NParray(Hsave,[n,n],'Hsave',name_sub)
      CALL alloc_NParray(r,[n],'r',name_sub)
      Hsave(:,:) = H(:,:)
    END IF
  ENDIF

  IF (.NOT. sym) THEN
    IF(keep_MPI) CALL diagonalization(H,E,Vec,n,4,1,.TRUE.)
  ELSE
    IF(keep_MPI) CALL diagonalization(H,E,Vec,n,3,1,.TRUE.)
  END IF

  IF (debug) THEN
    write(out_unit,*) 'OpMin,OpMax (ua)  : ',minval(E),maxval(E)
    write(out_unit,*) 'OpMin,OpMax (cm-1): ',minval(E)*auTOcm_inv,maxval(E)*auTOcm_inv
  END IF

  IF(keep_MPI) THEN
    DO i=1,n
      A = maxval(Vec(:,i))
      B = minval(Vec(:,i))
      IF (abs(A) < abs(B)) Vec(:,i) = -Vec(:,i)
    END DO
  ENDIF
  
  IF(keep_MPI) THEN
    IF (residual) THEN
      DO i=1,n
        r(:) = matmul(Hsave,Vec(:,i))-E(i)*Vec(:,i)
        write(out_unit,*) 'residual (cm-1)',i,sqrt(dot_product(r,r))*&
                                                         auTOcm_inv
        write(out_unit,*) 'residual (au)',i,sqrt(dot_product(r,r))
        write(out_unit,*) 'err (cm-1)',i,dot_product(Vec(:,i),r)*    &
                                                         auTOcm_inv
        write(out_unit,*) 'err (au)',i,dot_product(Vec(:,i),r)
      END DO
      CALL dealloc_NParray(Hsave,'Hsave',name_sub)
      CALL dealloc_NParray(r,'r',name_sub)
    END IF
  END IF !for keep_MPI!

!------------------------------------------------------
  IF (debug) THEN

    write(out_unit,*) ' level energy :'
    DO i=1,n
      write(out_unit,*) i,E(i)*auTOcm_inv,(E(i)-E(1))*auTOcm_inv
    END DO

    write(out_unit,*) ' Vec:'
    CALL Write_Mat_MPI(Vec,out_unit,5)

    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
!------------------------------------------------------

END SUBROUTINE sub_diago_H
!=======================================================================================

!
!=====================================================================
!
!     diagonalisation complexe
!
!=====================================================================
      SUBROUTINE sub_diago_CH(CH,CE,CVec,n)
      USE EVR_system_m
      USE mod_Constant, ONLY: get_Conv_au_TO_unit
      IMPLICIT NONE
      !
      !------ active Matrix H Vec E ------------------------------------
      integer              :: n
      complex (kind=Rkind) :: CH(n,n)
      complex (kind=Rkind) :: CVec(n,n)
      complex (kind=Rkind) :: CE(n)
      complex (kind=Rkind) :: S(n,n)

      real (kind=Rkind)    :: auTOcm_inv


!----- divers --------------------------------------------------------
      integer              :: i


!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_diago_CH'
      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

!=====================================================================
!
!     H matrix diagonalisation
!
!=====================================================================

      CALL diagonalization(CH,CE,CVec,diago_type=4)

!=====================================================================
!
!
!=====================================================================

!-----------------------------------------------------------
      IF (debug) THEN

        write(out_unit,*) ' level energy :'
        DO i=1,n
          write(out_unit,*) i,CE(i)*cmplx(auTOcm_inv,kind=Rkind),     &
                             (CE(i)-CE(1))*cmplx(auTOcm_inv,kind=Rkind)
        END DO

        write(out_unit,*) ' Vec:'
        CALL Write_Mat_MPI(CVec,out_unit,5)

        IF (debug) THEN
          S = matmul(transpose(CVec),CVec)
          CALL Write_Mat_MPI(S,out_unit,5,info='eigenvector Overlap1')
        END IF

stop
        write(out_unit,*) 'END sub_diago_CH'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_diago_CH
