!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
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
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
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
      SUBROUTINE sub_quadra_FT(base,nosym,nstep,nb_shift,tab_shift)

      USE mod_system
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
      integer       :: i,ii,k,nq
!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical :: debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_FT'
         write(out_unitp,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------

      IF (base%nb .EQ. 0) RETURN
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_fourier et nb_quadra
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: FT'
        write(out_unitp,*) '      old nb_FT',base%nb
        IF (mod(base%nb,nb_shift) == 0) THEN
          nb_i = base%nb/nb_shift
        ELSE
          nb_i = base%nb/nb_shift + 1
        END IF
        base%nb = nb_i * nb_shift
        write(out_unitp,*) '      nb_FT',base%nb
        write(out_unitp,*) '      nb_shift:',nb_shift,'tab',tab_shift
      END IF

      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          write(out_unitp,*) '      old nb_quadra',nq
          IF ( nq < nb_i*nstep ) nq = nb_i*nstep
          IF (nq == nb_i*nstep .AND. mod(nq,2) == 0) nq = nq + 1
          write(out_unitp,*) '      new nb_quadra',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)
        CALL alloc_xw_OF_basis(base)

        CALL gauss_fourier(base%x(1,:),base%w,nq)
        dx = base%x(1,2)-base%x(1,1)
        base%x = base%x + dx*HALF
        base%wrho(:) = base%w(:)
      END IF
      base%rho(:)  = ONE


!     calcul des valeurs de la serie de fourier et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.


      IF (.NOT. base%cplx) THEN
        write(out_unitp,*) ' ERROR in sub_quadra_FT'
        write(out_unitp,*) ' the basis must be complex !!'
        write(out_unitp,*) ' base%cplx',base%cplx
        STOP
      END IF


      CALL alloc_dnb_OF_basis(base)
      ib = 0
      DO ishift=1,nb_shift

        DO i=1,nb_i
          ib = ib + 1
          ii = nstep*i+tab_shift(ishift)
          IF (base%print_info_OF_basisDP) write(out_unitp,*) 'base cos/sin:',ii
          base%tab_ndim_index(1,ib) = ii
          DO k=1,nq
            CALL d0d1d2d3fourier(base%x(1,k),d0,d1,d2,d3,ii)
            base%dnCGB%d0(k,ib)     = cmplx(d0,kind=Rkind)
            base%dnCGB%d1(k,ib,1)   = cmplx(d1,kind=Rkind)
            base%dnCGB%d2(k,ib,1,1) = cmplx(d2,kind=Rkind)
!           write(out_unitp,*) ib,k,ii,d0,d1,d2
          END DO
        END DO

      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_fourier'
      END IF
!-----------------------------------------------------------

      RETURN
      end subroutine sub_quadra_FT

