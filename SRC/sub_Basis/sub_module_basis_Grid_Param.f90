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
      MODULE mod_Basis_Grid_Param
      USE EVR_system_m
      IMPLICIT NONE

        PRIVATE

        TYPE Basis_Grid_Param
          integer :: LGrid_max = -1
          integer :: nq        =  0
          integer :: nq_init   =  0
        CONTAINS
          PROCEDURE, PRIVATE, PASS(Basis_Grid_Para1) :: Basis_Grid_Param2TOBasis_Grid_Param1
          GENERIC,   PUBLIC  :: assignment(=) => Basis_Grid_Param2TOBasis_Grid_Param1
        END TYPE Basis_Grid_Param

        PUBLIC Basis_Grid_Param, Write_Basis_Grid_Param,                &
               Basis_Grid_ParamTOBasis_Grid_Param_init,                 &
               Basis_Grid_Param_initTOBasis_Grid_Param

      CONTAINS

      SUBROUTINE Write_Basis_Grid_Param(Basis_Grid_Para,Rec_line)

      TYPE (Basis_Grid_Param), intent(in) :: Basis_Grid_Para
      character (len=*), optional         :: Rec_line

      integer :: L
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_Basis_Grid_Param'
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

      IF (present(Rec_line)) THEN
       write(out_unit,*) trim(Rec_line),'LGridmax',Basis_Grid_Para%LGrid_max
       write(out_unit,*) trim(Rec_line),'nq',Basis_Grid_Para%nq
       write(out_unit,*) trim(Rec_line),'nq_ini',Basis_Grid_Para%nq_init

      ELSE
       write(out_unit,*) 'LGridmax',Basis_Grid_Para%LGrid_max
       write(out_unit,*) 'nq',Basis_Grid_Para%nq
       write(out_unit,*) 'nq_ini',Basis_Grid_Para%nq_init
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Write_Basis_Grid_Param

      SUBROUTINE Basis_Grid_Param2TOBasis_Grid_Param1(Basis_Grid_Para1, &
                                                      Basis_Grid_Para2)


      CLASS (Basis_Grid_Param), intent(inout) :: Basis_Grid_Para1
      TYPE (Basis_Grid_Param),  intent(in) :: Basis_Grid_Para2


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Basis_Grid_Param2TOBasis_Grid_Param1'
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


      Basis_Grid_Para1%LGrid_max = Basis_Grid_Para2%LGrid_max
      Basis_Grid_Para1%nq        = Basis_Grid_Para2%nq
      Basis_Grid_Para1%nq_init   = Basis_Grid_Para2%nq_init

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_Grid_Param(Basis_Grid_Para1)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Basis_Grid_Param2TOBasis_Grid_Param1

      SUBROUTINE Basis_Grid_ParamTOBasis_Grid_Param_init(Basis_Grid_Para)


      TYPE (Basis_Grid_Param), intent(inout) :: Basis_Grid_Para


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Basis_Grid_ParamTOBasis_Grid_Param_init'
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


      Basis_Grid_Para%nq_init   = Basis_Grid_Para%nq

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_Grid_Param(Basis_Grid_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Basis_Grid_ParamTOBasis_Grid_Param_init

      SUBROUTINE Basis_Grid_Param_initTOBasis_Grid_Param(Basis_Grid_Para)


      TYPE (Basis_Grid_Param), intent(inout) :: Basis_Grid_Para


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Basis_Grid_Param_initTOBasis_Grid_Param'
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


      Basis_Grid_Para%nq   = Basis_Grid_Para%nq_init

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_Grid_Param(Basis_Grid_Para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Basis_Grid_Param_initTOBasis_Grid_Param

      END MODULE mod_Basis_Grid_Param
