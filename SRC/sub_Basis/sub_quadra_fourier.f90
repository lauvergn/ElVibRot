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
!=============================================================
!
!      determination des tous les Ln(xi)=serie_fourier(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
  SUBROUTINE sub_quadra_fourier(base,nosym,nstep,nb_shift,tab_shift)
      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base
      logical       :: num,nosym
      real (kind=Rkind) :: step

      integer  :: nstep,nb_shift
      integer  :: tab_shift(nb_shift)
      integer  :: ib,ishift,nb_i
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      logical :: deriv
      real (kind=Rkind) :: d0,d1,d2,d3,dx
      integer       :: i,ii,k,l,nq
!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING sub_quadra_fourier'
         write(out_unit,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------

       IF (base%nb <= 0) STOP 'ERROR nb<=0'
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_fourier et nb_quadra
!      "cos 1 0" or "cos": all terms
!      "cos 2 0" : only sine
!      "cos 2 -1": only cosine
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis) THEN
        IF (print_level > -1) write(out_unit,*) '    Basis: Fourier series'
        IF (print_level > -1) write(out_unit,*) '      old nb_fourier',base%nb
        IF (mod(base%nb,nb_shift) == 0) THEN
          nb_i = base%nb/nb_shift
        ELSE
          nb_i = base%nb/nb_shift + 1
        END IF
        base%nb = nb_i * nb_shift
        IF (print_level > -1) write(out_unit,*) '      nb_fourier',base%nb
        IF (print_level > -1) write(out_unit,*) '      nb_shift:',nb_shift,'tab',tab_shift
      ELSE
        nstep    = 1
        nb_i     = base%nb
        nb_shift = 1
      END IF

      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unit,*) '      old nb_quadra',nq
          IF ( nq < nb_i*nstep ) nq = nb_i*nstep
          IF (nosym .AND. nq == nb_i*nstep .AND. mod(nq,2) == 0) nq = nq + 1
          IF (print_level > -1) write(out_unit,*) '      new nb_quadra',nq
        END IF

        IF (print_level > -1) write(out_unit,*) 'fourier: nb,nq',base%nb,nq
        CALL Set_nq_OF_basis(base,nq)
        CALL alloc_xw_OF_basis(base)
        CALL gauss_fourier(base%x(1,:),base%w,nq)
        IF (nosym) THEN
          dx = base%x(1,2)-base%x(1,1)
          base%x = base%x + dx*HALF
          !base%x = base%x - dx*HALF
        END IF
        base%wrho(:) = base%w(:)
      END IF
      base%rho(:)  = ONE

!     calcul des valeurs de la serie de fourier et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.


      CALL alloc_dnb_OF_basis(base)
      ib = 0
      DO ishift=1,nb_shift

        DO i=1,nb_i
          ib = ib + 1
          ii = nstep*i+tab_shift(ishift)
          IF (print_level > -1 .AND. base%print_info_OF_basisDP)  THEN
            IF (mod(ii,2) == 0 .AND. ii < 11) THEN
              write(out_unit,*) 'basis function: fourrier',ii,'or sin',int(ii/2)
            ELSE IF (mod(ii,2) == 1 .AND. ii < 11) THEN
              write(out_unit,*) 'basis function: fourrier',ii,'or cos',int(ii/2)
            END IF
          END IF
          base%tab_ndim_index(1,ib) = ii
          DO k=1,nq
            CALL d0d1d2d3fourier(base%x(1,k),d0,d1,d2,d3,ii)
            base%dnRGB%d0(k,ib)     = d0
            base%dnRGB%d1(k,ib,1)   = d1
            base%dnRGB%d2(k,ib,1,1) = d2
!           write(out_unit,*) ib,k,ii,d0,d1,d2
          END DO
        END DO

      END DO

      IF (base%nb == nq .AND. mod(base%nb,2) == 0) THEN
        base%dnRGB%d0(:,nq) = base%dnRGB%d0(:,nq) / sqrt(TWO)
        base%dnRGB%d1(:,nq,:) = base%dnRGB%d1(:,nq,:) / sqrt(TWO)
        base%dnRGB%d2(:,nq,:,:) = base%dnRGB%d2(:,nq,:,:) / sqrt(TWO)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unit,*) 'END sub_quadra_fourier'
      END IF
!-----------------------------------------------------------

  end subroutine sub_quadra_fourier
  SUBROUTINE RB_TO_RG_fourier(RB,RG,base)
  USE EVR_system_m
  USE mod_basis
  IMPLICIT NONE

  !---------------------------------------------------------------------
  !---------- variables passees en argument ----------------------------
  real (kind=Rkind), intent(in)    :: RB(:)
  real (kind=Rkind), intent(inout) :: RG(:)

  TYPE (basis), intent(in)         :: base

  integer :: nb,nq
  !----- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------
  nq = get_nq_FROM_basis(base)
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING RB_TO_RG_fourier'
    write(out_unit,*) 'nb,nq',base%nb,nq
    write(out_unit,*) 'RB',RB(:)
  END IF
  !-----------------------------------------------------------

  IF (base%nb < size(RB) .OR. nq /= size(RG)) THEN
    write(out_unit,*) ' ERROR in RB_TO_RG_fourier'
    write(out_unit,*) ' nb is inconsistent with size(RB)',base%nb,size(RB)
    write(out_unit,*) ' nq is inconsistent with size(RG)',nq,size(RG)
    STOP ' ERROR in RB_TO_RG_fourier: inconsistent sizes'
  END IF

  RG = matmul(base%dnRGB%d0(:,1:nb),RB(1:nb))

!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'RG',RG(:)
    write(out_unit,*) 'END RB_TO_RG_fourier'
  END IF
!-----------------------------------------------------------

end subroutine RB_TO_RG_fourier

