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
!      Particule-in-a-box basis set
!
!=============================================================
      SUBROUTINE sub_quadra_dfst(base,nosym)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)      :: base
      logical           :: nosym


      integer           :: ib
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      !---------------------------------------------------------------------
      logical            :: deriv
      real (kind=Rkind)  :: d0,d1,d2,d3,dx
      integer            :: i,ii,k,nq,m
      logical :: nb_OK
      integer                        :: err_allo
      character (len=:), allocatable :: err_msg
      !---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_dfst'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------

       IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur sub_quadra_box et nb_quadra
!----------------------------------------------------------------------------
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

      !find if nb+1=2**m
      m = log(real(base%nb+1,kind=Rkind))/log(TWO)
      nb_OK = (base%nb+1 == 2**m .OR. base%nb+1 == 2**(m+1) .OR. base%nb+1 == 2**(m-1))
      IF (.NOT. nb_OK) THEN
        write(out_unitp,*) '    Basis: ',name_sub
        write(out_unitp,*) '      nb_dfst',base%nb
        write(out_unitp,*) '      nb_dfst+1 IS NOT 2**m'
        STOP ' ERROR in sub_quadra_dfst: nb_dfst+1 IS NOT 2**m'
      END IF
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: ',name_sub
        write(out_unitp,*) '      nb_dfst',base%nb
      END IF
      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          write(out_unitp,*) '      old nb_quadra',nq
          IF ( nq /= base%nb ) nq = base%nb
          write(out_unitp,*) '      new nb_quadra',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)
        CALL alloc_xw_OF_basis(base)

        dx = pi/(base%nb+1)
        base%x(1,:)  = [(dx*i,i=1,base%nb)]
        base%w(:)    = dx
        base%wrho(:) = base%w(:)
        base%rho(:)  = ONE
      END IF
      base%rho(:)  = ONE

      CALL base%dsft%alloc(nq+1, err_allo,err_msg)
      IF (err_allo /= 0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Allocation error in QDUtil_alloc_FFT_OOURA'
        write(out_unitp,*) 'Error message: ',err_msg
        STOP 'ERROR in sub_quadra_dfst: Allocation error in QDUtil_alloc_FFT_OOURA'
      END IF
  

      CALL alloc_dnb_OF_basis(base)


      DO ib=1,base%nb
        base%tab_ndim_index(1,ib) = ib
        IF (debug) write(out_unitp,*) 'basis, particle in a dfst[0,Pi]:',ib
        DO k=1,nq
          CALL d0d1d2d3box(base%x(1,k),d0,d1,d2,d3,ib)
          base%dnRGB%d0(k,ib)     = d0
          base%dnRGB%d1(k,ib,1)   = d1
          base%dnRGB%d2(k,ib,1,1) = d2
         !write(out_unitp,*) ib,k,d0,d1,d2
        END DO
    END DO

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base,write_all=.TRUE.)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_dfst
