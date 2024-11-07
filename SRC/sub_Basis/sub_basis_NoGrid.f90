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
!      Basis set for several diabatic PES
!        For this basis nq must be equal to nb
!
!=============================================================
  SUBROUTINE Init_NoGrid_basis(base)
    USE EVR_system_m
    USE mod_basis
    IMPLICIT NONE

    !---------------------------------------------------------------------
    !---------- variables passees en argument ----------------------------
     TYPE (basis), intent(inout) :: base


      integer           :: i,Read_symab

!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='Init_NoGrid_basis'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'nb',base%nb
       END IF
!-----------------------------------------------------------

     IF (base%nb <= 0) STOP 'ERROR nb<=0'
     base%with_grid = .FALSE.

     base%ndim = 1
     ! nq is not defined
!----------------------------------------------------------------------------

      base%primitive_done = .TRUE.
      base%packed_done    = .TRUE.

      CALL dealloc_nDindex(base%nDindB)
      base%nDindB%packed = .TRUE.
      CALL init_nDindexPrim(base%nDindB,1,[base%nb])
      base%nDindB%With_L      = .TRUE.
      base%nDindB%Tab_L(:)    = 0
      base%nDindB%Tab_Norm(:) = ZERO


      CALL alloc_SymAbelian(base%P_SymAbelian,base%nb)
      Read_symab = Get_Read_symabOFSymAbelian(base%P_SymAbelian)
      CALL Set_ReadsymabOFSymAbelian(base%P_SymAbelian,Read_symab)
      DO i=1,base%nb
        CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,i,-1)
      END DO
      CALL Set_nbPERsym_FROM_SymAbelian(base%P_SymAbelian)
      !CALL Write_SymAbelian(base%P_SymAbelian)

      base%ndim = 1


!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE Init_NoGrid_basis

