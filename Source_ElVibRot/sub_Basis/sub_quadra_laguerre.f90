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
!  x.L"(x,n)+(1-x)L'(x,n)+n.L(x,n)=0
!  x.L'(x,n)-n.L(x,n)+n.L(x,n-1)=0
!
!  (n+1)L(x,n+1)+(x-2n-1)L(x,n)+n.L(x,n-1) = 0
!     with L(x,0)=1 and L(x,1)=1-x
! Normalization: Int[0,+inf] L(x,n)L(x,m).Exp(-x)dx = delta(n,m)
!
!=============================================================
      SUBROUTINE sub_quadra_laguerre(base)

      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE

!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base

!---------------------------------------------------------------------
!---------------------------------------------------------------------
      integer  :: iq,ib,nq
      real (kind=Rkind) :: B
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_laguerre'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)

       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------


       IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis .AND. print_level > -1) THEN
        write(out_unit,*) '    Basis: Laguerre polynomia'
        write(out_unit,*) '      nb_Laguerre',base%nb
      END IF
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unit,*) '      old nb_quadra',nq
          IF ( nq < base%nb ) nq = base%nb + 1
          IF (print_level > -1) write(out_unit,*) '      new nb_quadra',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)

        CALL alloc_xw_OF_basis(base)

        CALL laguerre_compute (nq,base%x(1,:),base%w,ZERO)
        !write(6,*) 'x_Laguerre',base%x(1,:)
        base%w(:) = base%w * exp(base%x(1,:))
        !B=30._Rkind
        !DO iq=1,nq
        !  base%x(1,iq) = B/nq*(iq-HALF)
        !  base%w(iq)   = B/nq
        !END DO
        base%wrho(:) = base%w(:)
        !base%xPOGridRep_done = .TRUE.
      END IF
      base%rho(:)  = ONE

      CALL alloc_dnb_OF_basis(base)

      IF (base%print_info_OF_basisDP .AND. print_level > -1)                    &
                        write(out_unit,*) '      All Laguerre polynomials'
        base%tab_ndim_index(1,:) = [(ib,ib=1,base%nb)]
        DO ib=1,base%nb
        DO iq=1,nq
          CALL d0d1d2poly_laguerre_weight(base%x(1,iq),                         &
                                          base%dnRGB%d0(iq,ib),                 &
                                          base%dnRGB%d1(iq,ib,1),               &
                                          base%dnRGB%d2(iq,ib,1,1),ib-1,0)
        END DO
        END DO
!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_laguerre
