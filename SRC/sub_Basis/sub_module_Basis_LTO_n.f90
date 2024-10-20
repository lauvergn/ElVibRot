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
      MODULE mod_Basis_L_TO_n
      USE EVR_system_m
      IMPLICIT NONE

        PRIVATE

        TYPE Basis_L_TO_n
          integer                     :: L_TO_n_type     = 0        ! for the initalization (after we use the Tab_L_TO_n)

          integer                     :: Lmax            = -1      ! parameter for the number of points of SparseGrid
          !!! relation between L and n: n(L) = A + B * L**expo or n(L) = Tab_L_TO_n(L)

          integer                     :: A                = 1
          integer                     :: B                = 1
          integer                     :: expo             = 1
          integer                     :: C                = 0       !  use nq(L) with basis parameters. ???? with Lexpo_TO_nq > 1
                                                                    ! Then when L>LB add (L-LB)*L_TO_nq_C
          integer                     :: max_n            = huge(1) ! value such n(L)<= max_n
          integer, allocatable        :: Tab_L_TO_n(:)
          integer, allocatable        :: Tab_n_TO_L(:)
          logical, allocatable        :: skip_deltaL(:)

        CONTAINS
          PROCEDURE, PRIVATE, PASS(L_TO_n_para1) :: L_TO_n_para2_TO_L_TO_n_para1
          GENERIC,   PUBLIC  :: assignment(=) => L_TO_n_para2_TO_L_TO_n_para1
        END TYPE Basis_L_TO_n

        PUBLIC :: Basis_L_TO_n, init_Basis_L_TO_n,                      &
                  Write_Basis_L_TO_n, Set_Basis_L_TO_n,                 &
                  dealloc_Basis_L_TO_n,                                 &
                  Get_n_FROM_Basis_L_TO_n, Get_L_FROM_Basis_L_TO_n,     &
                  check_Basis_L_TO_n

      CONTAINS

      SUBROUTINE alloc_Basis_L_TO_n(L_TO_n_para,Lmax)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para
      integer,             intent (in)   :: Lmax


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_Basis_L_TO_n'
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

      CALL dealloc_Basis_L_TO_n(L_TO_n_para)

      IF (Lmax < 0 ) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Wrong Lmax value',Lmax
        STOP
      END IF

      L_TO_n_para%Lmax           = Lmax

      CALL alloc_NParray(L_TO_n_para%tab_L_TO_n,[Lmax],               &
                        "L_TO_n_para%tab_L_TO_n",name_sub,[0] )

      CALL alloc_NParray(L_TO_n_para%skip_deltaL,[Lmax],               &
                        "L_TO_n_para%skip_deltaL",name_sub,[0] )

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE alloc_Basis_L_TO_n

      SUBROUTINE dealloc_Basis_L_TO_n(L_TO_n_para)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='dealloc_Basis_L_TO_n'
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

      L_TO_n_para%L_TO_n_type = 0
      L_TO_n_para%Lmax        = -1
      L_TO_n_para%A           = 1
      L_TO_n_para%B           = 1
      L_TO_n_para%C           = 0
      L_TO_n_para%expo        = 1


      IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
        CALL dealloc_NParray(L_TO_n_para%tab_L_TO_n,                    &
                            "L_TO_n_para%tab_L_TO_n",name_sub)
      END IF

      IF (allocated(L_TO_n_para%tab_n_TO_L)) THEN
        CALL dealloc_NParray(L_TO_n_para%tab_n_TO_L,                    &
                            "L_TO_n_para%tab_n_TO_L",name_sub)
      END IF

      IF (allocated(L_TO_n_para%skip_deltaL)) THEN
        CALL dealloc_NParray(L_TO_n_para%skip_deltaL,                    &
                            "L_TO_n_para%skip_deltaL",name_sub)
      END IF
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE dealloc_Basis_L_TO_n

      SUBROUTINE Write_Basis_L_TO_n(L_TO_n_para,Rec_line)
      TYPE (Basis_L_TO_n), intent(in) :: L_TO_n_para
      character (len=*), optional     :: Rec_line

      integer :: L
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_Basis_L_TO_n'
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
        write(out_unit,*) trim(Rec_line),'L_TO_n_type     ',L_TO_n_para%L_TO_n_type
        write(out_unit,*) trim(Rec_line),'Lmax            ',L_TO_n_para%Lmax
        write(out_unit,*) trim(Rec_line),'A               ',L_TO_n_para%A
        write(out_unit,*) trim(Rec_line),'B               ',L_TO_n_para%B
        write(out_unit,*) trim(Rec_line),'C               ',L_TO_n_para%C
        write(out_unit,*) trim(Rec_line),'expo            ',L_TO_n_para%expo
        write(out_unit,*) trim(Rec_line),'max_n           ',L_TO_n_para%max_n

        write(out_unit,*) trim(Rec_line),'alloc tab_L_TO_n',allocated(L_TO_n_para%tab_L_TO_n)
        IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
          write(out_unit,*) trim(Rec_line),'tab_L_TO_n'
          write(out_unit,*) trim(Rec_line),':  ',L_TO_n_para%tab_L_TO_n
        END IF

        write(out_unit,*) trim(Rec_line),'alloc tab_n_TO_L',allocated(L_TO_n_para%tab_n_TO_L)
        IF (allocated(L_TO_n_para%tab_n_TO_L)) THEN
          write(out_unit,*) trim(Rec_line),'tab_n_TO_L'
          write(out_unit,*) trim(Rec_line),':  ',L_TO_n_para%tab_n_TO_L
        END IF

        write(out_unit,*) trim(Rec_line),'alloc skip_deltaL',allocated(L_TO_n_para%skip_deltaL)
        IF (allocated(L_TO_n_para%skip_deltaL)) THEN
          write(out_unit,*) trim(Rec_line),'skip_deltaL',L_TO_n_para%skip_deltaL
        END IF

      ELSE
        write(out_unit,*) 'L_TO_n_type     ',L_TO_n_para%L_TO_n_type
        write(out_unit,*) 'Lmax            ',L_TO_n_para%Lmax
        write(out_unit,*) 'A               ',L_TO_n_para%A
        write(out_unit,*) 'B               ',L_TO_n_para%B
        write(out_unit,*) 'C               ',L_TO_n_para%C
        write(out_unit,*) 'expo            ',L_TO_n_para%expo
        write(out_unit,*) 'max_n           ',L_TO_n_para%max_n

        write(out_unit,*) 'alloc tab_L_TO_n',allocated(L_TO_n_para%tab_L_TO_n)
        IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
          write(out_unit,*) 'tab_L_TO_n'
          write(out_unit,*) '   ',L_TO_n_para%tab_L_TO_n
        END IF

        write(out_unit,*) 'alloc tab_n_TO_L',allocated(L_TO_n_para%tab_n_TO_L)
        IF (allocated(L_TO_n_para%tab_n_TO_L)) THEN
          write(out_unit,*) 'tab_n_TO_L'
          write(out_unit,*) '   ',L_TO_n_para%tab_n_TO_L
        END IF

        write(out_unit,*) 'alloc skip_deltaL',allocated(L_TO_n_para%skip_deltaL)
        IF (allocated(L_TO_n_para%skip_deltaL)) THEN
          write(out_unit,*) 'skip_deltaL',L_TO_n_para%skip_deltaL
        END IF

      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Write_Basis_L_TO_n

      SUBROUTINE init_Basis_L_TO_n(L_TO_n_para,Lmax)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para
      integer,             intent (in)   :: Lmax

      integer  :: L,errBasis_L_TO_n,nmin,nmax,n


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='init_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)

      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

    L_TO_n_para%Lmax = Lmax
    IF (L_TO_n_para%L_TO_n_type /=2) THEN
      IF (allocated(L_TO_n_para%tab_L_TO_n))                            &
           CALL dealloc_NParray(L_TO_n_para%tab_L_TO_n,"L_TO_n_para%tab_L_TO_n",name_sub)

      IF (Lmax < 0 ) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Wrong Lmax value',Lmax
        STOP
      END IF

      L_TO_n_para%Lmax = Lmax
      CALL alloc_NParray(L_TO_n_para%tab_L_TO_n,[Lmax],               &
                        "L_TO_n_para%tab_L_TO_n",name_sub,[0] )
    END IF

      SELECT CASE (L_TO_n_para%L_TO_n_type)
      CASE (0) ! n(L) = A + B * L**expo
        DO L=0,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%A + L_TO_n_para%B * L**L_TO_n_para%expo
        END DO

      CASE (1) ! Delta_n(L) = B * L**(expo-1)  => n(L) = A + B * L*(L+1)/2 (with expo=2)
        L_TO_n_para%tab_L_TO_n(0) = L_TO_n_para%A
        DO L=1,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%tab_L_TO_n(L-1) + L_TO_n_para%B * L**(L_TO_n_para%expo-1)
        END DO

      CASE (2) ! The tab_L_TO_n is already present
        CONTINUE !

      CASE (3) ! Delta_n(L) = B * L**(expo-1)  => n(L) = A + B * L*(L-1)/2 (with expo=2)
        L_TO_n_para%tab_L_TO_n(0) = L_TO_n_para%A
        DO L=1,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%tab_L_TO_n(L-1) + L_TO_n_para%B * int(ONE+log(real(L,kind=Rkind)))
        END DO

      CASE (4) ! n(L) = A+ B * tanh(L/4)
        DO L=0,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%A +       &
            int( real(L_TO_n_para%B,kind=Rkind) * tanh(real(L,kind=Rkind)/FOUR) )
        END DO

      CASE Default
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) '  WRONG L_TO_n_type',L_TO_n_para%L_TO_n_type
        STOP
      END SELECT

      WHERE (L_TO_n_para%tab_L_TO_n > L_TO_n_para%max_n)
        L_TO_n_para%tab_L_TO_n = L_TO_n_para%max_n
      END WHERE

      CALL check_Basis_L_TO_n(L_TO_n_para%tab_L_TO_n,errBasis_L_TO_n)

      IF (errBasis_L_TO_n /= 0) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) '  Problem with initialization'
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
        STOP
      END IF

      nmax = get_n_FROM_Basis_L_TO_n(L_TO_n_para,Lmax)
      nmin = get_L_FROM_Basis_L_TO_n(L_TO_n_para,nmax) ! to set up the table

      ! set skip_deltaL(:)
      !write(out_unit,*) 'init skip_deltaL'
      CALL alloc_NParray(L_TO_n_para%skip_deltaL,[Lmax],                      &
                        "L_TO_n_para%skip_deltaL",name_sub,[0] )
      L_TO_n_para%skip_deltaL(:) = .FALSE.
      DO L=1,Lmax ! start at L=1, because L-1 is used
        L_TO_n_para%skip_deltaL(L)= (L_TO_n_para%tab_L_TO_n(L) == L_TO_n_para%tab_L_TO_n(L-1))
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE init_Basis_L_TO_n

      FUNCTION Get_n_FROM_Basis_L_TO_n(L_TO_n_para,L,L2) RESULT(n)
      integer  :: n
      TYPE (Basis_L_TO_n), intent(in)    :: L_TO_n_para
      integer,             intent (in)   :: L
      integer, optional,   intent (in)   :: L2

      integer :: Ltmp,u,L2_loc,n2B,n2C

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Get_n_FROM_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'L',L
        write(out_unit,*) 'L2?',present(L2)
        IF (present(L2)) write(out_unit,*) 'L2',L2
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------



    IF (L < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) 'L',L
      write(out_unit,*) 'L2?',present(L2)
      IF (present(L2)) write(out_unit,*) 'L2',L2
      CALL Write_Basis_L_TO_n(L_TO_n_para)
      write(out_unit,*) ' ERROR:  L < 0',L
      STOP 'ERROR in Get_n_FROM_Basis_L_TO_n: L < 0'
    END IF

    L2_loc = L
    IF (present(L2)) L2_loc = L2

    IF (L2_loc < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) 'L',L
      write(out_unit,*) 'L2?',present(L2)
      IF (present(L2)) write(out_unit,*) 'L2',L2
      CALL Write_Basis_L_TO_n(L_TO_n_para)
      write(out_unit,*) ' ERROR:  L2 < 0',L2_loc
      STOP 'ERROR in Get_n_FROM_Basis_L_TO_n: L2 < 0'
    END IF

    IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
      u = ubound(L_TO_n_para%Tab_L_TO_n,dim=1)

      IF (L <= u) THEN
        n = L_TO_n_para%tab_L_TO_n(L)
      ELSE
        n = L_TO_n_para%Tab_L_TO_n(u)
        IF (u > 0) n = n +                                                      &
                   (l-u)*(L_TO_n_para%Tab_L_TO_n(u)-L_TO_n_para%Tab_L_TO_n(u-1))
      END IF
    ELSE

      SELECT CASE (L_TO_n_para%L_TO_n_type)
      CASE (0) ! A + B * L**expo
        IF (L <= L2_loc) THEN
          n   = L_TO_n_para%A + L_TO_n_para%B * L     **L_TO_n_para%expo
        ELSE
          n2B = L_TO_n_para%A + L_TO_n_para%B * L2_loc**L_TO_n_para%expo
          n2C = L_TO_n_para%A + L_TO_n_para%C * L2_loc**L_TO_n_para%expo
          n   = L_TO_n_para%A + (n2B-n2C)     + L_TO_n_para%C * L     **L_TO_n_para%expo
        END IF

      CASE (1) ! Delta_n(L) = B * L**(expo-1)  => n(L) = A + B * L*(L-1)/2 (with expo=2)
        n = L_TO_n_para%A
        DO Ltmp=1,L
          n = n + L_TO_n_para%B * Ltmp**(L_TO_n_para%expo-1)
        END DO

      CASE (2) ! read tab_L_TO_n
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) '   L_TO_n_type=2 and tab_L_TO_n is not allocated !'
        STOP 'ERROR in Get_n_FROM_Basis_L_TO_n: L_TO_n_type=2 and tab_L_TO_n is not allocated'

      CASE (3)
        n = L_TO_n_para%A
        DO Ltmp=1,L
          n = n + L_TO_n_para%B * int(ONE+log(real(Ltmp,kind=Rkind)))
        END DO

      CASE (4)
        n = L_TO_n_para%A +       &
          int( real(L_TO_n_para%B,kind=Rkind) * tanh(real(L,kind=Rkind)/FOUR) )

      CASE Default
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) '  WRONG L_TO_n_type',L_TO_n_para%L_TO_n_type
        STOP 'ERROR in Get_n_FROM_Basis_L_TO_n: WRONG L_TO_n_type'
      END SELECT

    END IF

    n = min(n,L_TO_n_para%max_n)

    IF (n < 0) THEN
     write(out_unit,*) 'ERROR in ',name_sub
     write(out_unit,*) '  n < 0, n:',n
     write(out_unit,*) 'L',L
     write(out_unit,*) 'L2?',present(L2)
     IF (present(L2)) write(out_unit,*) 'L2',L2
     CALL Write_Basis_L_TO_n(L_TO_n_para)
     write(out_unit,*)
     STOP 'ERROR in Get_n_FROM_Basis_L_TO_n: n<0'
    END IF

    IF (n == 0) THEN
     write(out_unit,*) 'WARNING in ',name_sub, '  n == 0'
     flush(out_unit)
    END IF
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'L,max_n,n',L,L_TO_n_para%max_n,n
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END FUNCTION Get_n_FROM_Basis_L_TO_n

      FUNCTION Get_L_FROM_Basis_L_TO_n(L_TO_n_para,n) RESULT(L)
      integer  :: L
      TYPE (Basis_L_TO_n), intent(inout)    :: L_TO_n_para
      integer,             intent (in)   :: n


      integer :: nni,nn0,nn,nmin,nmax

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Get_L_FROM_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'n',n
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

    nmax = get_n_FROM_Basis_L_TO_n(L_TO_n_para,L_TO_n_para%Lmax)
    nmin = get_n_FROM_Basis_L_TO_n(L_TO_n_para,0)

    IF (n < nmin .OR. n > nmax) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      CALL Write_Basis_L_TO_n(L_TO_n_para)
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) 'n is out of range. n:',n
      write(out_unit,*) '   range:',nmin,nmax
      STOP
    END IF

    IF (.NOT. allocated(L_TO_n_para%tab_n_TO_L)) THEN
      CALL alloc_NParray(L_TO_n_para%tab_n_TO_L,[nmax],'tab_n_TO_L',name_sub,[nmin])

      DO L=L_TO_n_para%Lmax,0,-1
        nn = get_n_FROM_Basis_L_TO_n(L_TO_n_para,L)
        L_TO_n_para%tab_n_TO_L(nmin:nn) = L
      END DO
    END IF

    L = L_TO_n_para%tab_n_TO_L(n)



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'L',L
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END FUNCTION Get_L_FROM_Basis_L_TO_n

      SUBROUTINE Set_Basis_L_TO_n(L_TO_n_para,A,B,C,expo,Tab_L_TO_n,max_n,L_TO_n_type)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para
      integer, optional,   intent (in)   :: A,B,C,expo,L_TO_n_type,max_n
      integer, optional,   intent (in), allocatable   :: Tab_L_TO_n(:)


      integer  :: L,Aloc,Bloc,Cloc,expoloc,L_TO_n_type_loc,max_n_loc
      integer  :: errBasis_L_TO_n


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_Basis_L_TO_n'
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

      Aloc = 1
      IF (present(A)) Aloc=A

      Bloc = 1
      IF (present(B)) Bloc=B

      Cloc = 0
      IF (present(C)) Cloc=C

      expoloc = 1
      IF (present(expo)) expoloc=expo

      L_TO_n_type_loc = 0
      IF (present(L_TO_n_type)) L_TO_n_type_loc=L_TO_n_type

      max_n_loc = huge(1)
      IF (present(max_n)) max_n_loc=max_n

      CALL dealloc_Basis_L_TO_n(L_TO_n_para)

      L_TO_n_para%A           = Aloc
      L_TO_n_para%B           = Bloc
      L_TO_n_para%C           = Cloc
      L_TO_n_para%expo        = expoloc
      L_TO_n_para%max_n       = max_n_loc
      L_TO_n_para%L_TO_n_type = L_TO_n_type_loc

      IF (present(Tab_L_TO_n)) THEN
        IF (.NOT. allocated(Tab_L_TO_n)) STOP ' in Set_Basis_L_TO_n: Tab_L_TO_n is not allocated'
        L_TO_n_para%L_TO_n_type = 2
        L_TO_n_para%Tab_L_TO_n  = Tab_L_TO_n
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Set_Basis_L_TO_n

      SUBROUTINE check_Basis_L_TO_n(tab_L_TO_n,errBasis_L_TO_n)

      integer, allocatable,     intent (in)    :: tab_L_TO_n(:)
      integer,                  intent (inout) :: errBasis_L_TO_n


      integer :: L
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='check_Basis_L_TO_n'
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
      errBasis_L_TO_n = 0

      IF (.NOT. allocated(tab_L_TO_n)) errBasis_L_TO_n = 1
      IF (errBasis_L_TO_n /= 0) RETURN

      IF (tab_L_TO_n( lbound(tab_L_TO_n,dim=1) ) < 0)  errBasis_L_TO_n = -1 ! tab_L_TO_n(0) < 1
      IF (errBasis_L_TO_n /= 0) RETURN


      DO L=lbound(tab_L_TO_n,dim=1)+1,ubound(tab_L_TO_n,dim=1) ! tab_L_TO_n(L) < tab_L_TO_n(L-1)
        IF ( tab_L_TO_n( L) < tab_L_TO_n( L-1) ) THEN
          errBasis_L_TO_n = -1
          RETURN
        END IF
      END DO



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE check_Basis_L_TO_n

      SUBROUTINE L_TO_n_para2_TO_L_TO_n_para1(L_TO_n_para1,L_TO_n_para2)
      CLASS (Basis_L_TO_n), intent(inout) :: L_TO_n_para1
      TYPE (Basis_L_TO_n),  intent(in)    :: L_TO_n_para2


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='L_TO_n_para2_TO_L_TO_n_para1'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'L_TO_n_para2'
        CALL Write_Basis_L_TO_n(L_TO_n_para2)
        flush(out_unit)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      L_TO_n_para1%L_TO_n_type = L_TO_n_para2%L_TO_n_type
      L_TO_n_para1%Lmax        = L_TO_n_para2%Lmax
      L_TO_n_para1%A           = L_TO_n_para2%A
      L_TO_n_para1%B           = L_TO_n_para2%B
      L_TO_n_para1%C           = L_TO_n_para2%C
      L_TO_n_para1%expo        = L_TO_n_para2%expo
      L_TO_n_para1%max_n       = L_TO_n_para2%max_n

      IF (allocated(L_TO_n_para2%tab_L_TO_n))  L_TO_n_para1%tab_L_TO_n   = L_TO_n_para2%tab_L_TO_n
      IF (allocated(L_TO_n_para2%tab_n_TO_L))  L_TO_n_para1%tab_n_TO_L   = L_TO_n_para2%tab_n_TO_L
      IF (allocated(L_TO_n_para2%skip_deltaL)) L_TO_n_para1%skip_deltaL  = L_TO_n_para2%skip_deltaL

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'L_TO_n_para1'
        CALL Write_Basis_L_TO_n(L_TO_n_para1)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE L_TO_n_para2_TO_L_TO_n_para1

      END MODULE mod_Basis_L_TO_n
