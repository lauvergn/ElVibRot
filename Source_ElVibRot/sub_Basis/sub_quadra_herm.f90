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
      SUBROUTINE sub_quadra_hermite(base,paire)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      integer          :: paire
      logical          :: num
      real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical                       :: deriv
      integer                       :: i,iq,nq,nb_nosym
      logical                       :: Print_basis
      TYPE (File_t)                 :: herm_file
      integer                       :: idum,nio,err_io
      real(kind=Rkind)              :: wdum
      logical, parameter            :: WithRead = .TRUE.

      real(kind=Rkind), allocatable :: x_loc(:,:),w_loc(:)


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
         write(out_unitp,*) 'nb',base%nb
         write(out_unitp,*) 'nq',nq
       END IF
!-----------------------------------------------------------
       deriv = .TRUE.
       num   = .FALSE.

       Print_basis = base%print_info_OF_basisDP .AND. print_level > -1

       IF (base%nb <= 0) THEN
         base%nb = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)
       END IF

       IF (base%nb <= 0) THEN
         write(out_unitp,*) 'ERROR in sub_quadra_hermite'
         write(out_unitp,*) 'nb<=0',base%nb
         STOP 'ERROR nb<=0'
       END IF

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (Print_basis) THEN
        write(out_unitp,*) '    Basis: Hermite polynomials'
        write(out_unitp,*) '      nb_hermite',base%nb
      END IF

      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (Print_basis) write(out_unitp,*) '      nb_hermite',base%nb
          IF (Print_basis) write(out_unitp,*) '      old nb_quadra',nq
          IF (paire == -1) THEN
            nb_nosym = base%nb
          ELSE IF (paire == 1) THEN ! odd
            nb_nosym = 2*base%nb - 1
          ELSE ! even
            nb_nosym = 2*base%nb
          END IF
        END IF
        flush(out_unitp)
        SELECT CASE (base%Nested)
        CASE(1)

          IF (mod(base%nq_max_Nested,2) == 0) base%nq_max_Nested = base%nq_max_Nested + 1

          IF (base%check_nq_OF_basis .AND. mod(nq,2) == 0)    nq = nq + 1 ! here nq must be odd : 2n-1
          IF (base%check_nq_OF_basis .AND. nq < 2*nb_nosym-1) nq = 2*nb_nosym-1
          IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq

          CALL alloc_NParray(x_loc,[1,nq],'x_loc',name_sub)
          CALL alloc_NParray(w_loc,[nq],  'w_loc',name_sub)

          CALL grid_HermiteNested1(x_loc,w_loc,nq,base%nq_max_Nested)

        CASE(2) ! not yet
          IF (mod(base%nq_max_Nested,2) == 0) base%nq_max_Nested = base%nq_max_Nested + 1

          IF (base%check_nq_OF_basis .AND. mod(nq,2) == 0)    nq = nq + 1 ! here nq must be odd : 2n-1
          IF (base%check_nq_OF_basis .AND. nq < 2*nb_nosym-1) nq = 2*nb_nosym-1
          IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq

          CALL alloc_NParray(x_loc,[1,nq],'x_loc',name_sub)
          CALL alloc_NParray(w_loc,[nq],  'w_loc',name_sub)


          STOP 'not yet nested2'
        CASE Default

          IF (base%check_nq_OF_basis) THEN
            IF (nq < nb_nosym) nq = nb_nosym + 1
          END IF
          IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq
          flush(out_unitp)

          CALL alloc_NParray(x_loc,[1,nq],'x_loc',name_sub)
          CALL alloc_NParray(w_loc,[nq],  'w_loc',name_sub)

          CALL hercom(nq,x_loc,w_loc)

          IF (WithRead) THEN
            herm_file%name = trim(EVRT_path) //                                 &
              '/Internal_data/HermQuadra/herm' // TO_string(nq) // '.txt'

            write(out_unitp,*) 'herm_file%name: ',herm_file%name
            flush(out_unitp)
            CALL file_open(herm_file,nio,old=.TRUE.,err_file=err_io)

            IF (err_io /= 0) THEN
              write(out_unitp,*) 'WARNNING in ',name_sub
              write(out_unitp,*) 'cannot find quadrature file: ',herm_file%name
            ELSE

              DO iq=1,nq
                read(nio,*,iostat=err_io) idum,x_loc(:,iq),wdum,w_loc(iq)
              END DO

            END IF
            CALL file_close(herm_file)

          END IF
        END SELECT
        !base%xPOGridRep_done = .TRUE.

        nq = nq + base%nq_extra
        CALL Set_nq_OF_basis(base,nq)
        CALL alloc_xw_OF_basis(base)
        base%w(:) = ZERO

        base%x(:,1:size(w_loc)) = x_loc
        base%w(1:size(w_loc))   = w_loc

        IF (allocated(base%x_extra)) THEN
          base%x(:,size(w_loc)+1:nq) = base%x_extra(:,:)
        END IF

        base%wrho(:) = base%w(:)
        base%rho(:)  = ONE

        IF (allocated(x_loc)) CALL dealloc_NParray(x_loc,'x_loc',name_sub)
        IF (allocated(w_loc)) CALL dealloc_NParray(w_loc,'w_loc',name_sub)
      END IF

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN
        IF (Print_basis) write(out_unitp,*) '      even Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i-1,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_0_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nq,deriv,num,step)
      ELSE IF (paire == 1) THEN
        IF (Print_basis) write(out_unitp,*) '      odd Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_1_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nq,deriv,num,step)
      ELSE
        IF (Print_basis) write(out_unitp,*) '      All Hermite polynomials'
        base%tab_ndim_index(1,:) = [(i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_exp_grille(                                     &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nq,deriv,num,step)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base,write_all=.TRUE.)
        write(out_unitp,*) 'END sub_quadra_hermite'
      END IF
!-----------------------------------------------------------

      end subroutine sub_quadra_hermite

  SUBROUTINE sub_quadra_hermiteAB(base,paire)
    USE mod_system
    USE mod_basis
    USE ADdnSVM_m
    IMPLICIT NONE
  
    !---------- variables passees en argument ----------------------------
    TYPE (basis), intent(inout) :: base
    integer,      intent(in)    :: paire
  
    integer           :: nq,iq,ib
    real (kind=Rkind) :: d0u,d0x,d1u(1),d2u(1,1)
    TYPE (dnS_t)      :: x,u,dnf,dnbu
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='sub_quadra_hermiteAB'
    logical, parameter :: debug = .TRUE.
    !logical, parameter :: debug = .FALSE.
    !-----------------------------------------------------------
    nq = get_nq_FROM_basis(base)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nb',base%nb
      write(out_unitp,*) 'nq',nq
    END IF
    !-----------------------------------------------------------

    ! first call the usual HO basis
    CALL sub_quadra_hermite(base,paire)
    
    DO iq=1,nq
      ! first the new grid point
      d0u = base%x(1,iq) ! the grid point in in the ]-inf,+inf[ range
      base%x(1,iq) = (base%B(1)-base%A(1))*HALF * tanh(d0u) + (base%B(1)+base%A(1))*HALF 
      d0x = base%x(1,iq) ! the grid point in the ]A,B[ range

      ! 2d the new weight x rho ! Warning the new rho is added in wrho
      x = variable(d0x,nVar=1,nderiv=3)
      u = atanh((x+x-base%A(1)-base%B(1))/(base%B(1)-base%A(1)))
      IF (abs(get_d0(u) - d0u) > ONETENTH**9) STOP 'ERROR get_d0(u) /= d0u'

      !base%wrho(iq) = base%wrho(iq) * d1u(1)

      !basis derivatives
      !base%dnRGB%d2(iq,:,1,1) = (d2u(1,1)* base%dnRGB%d1(iq,:,1) + d1u(1)**2 * base%dnRGB%d2(iq,:,1,1))
      !base%dnRGB%d1(iq,:,1)   = d1u(1)  * base%dnRGB%d1(iq,:,1)

      DO ib=1,base%nb
        CALL set_dnS(dnf,d0=base%dnRGB%d0(iq,ib),       &
          d1=[base%dnRGB%d1(iq,ib,1)], &
          d2=reshape([base%dnRGB%d2(iq,ib,1,1)],shape=[1,1]))

        dnbu = dnF_OF_dnS(dnf,u)*sqrt(deriv(u,ider=1))
        base%dnRGB%d0(iq,ib)     = get_d0(dnbu)
        base%dnRGB%d1(iq,ib,:)   = get_d1(dnbu)
        base%dnRGB%d2(iq,ib,:,:) = get_d2(dnbu)

        !base%dnRGB%d2(iq,ib,1,1) = d2u(1,1)* base%dnRGB%d1(iq,ib,1) + d1u(1)**2 * base%dnRGB%d2(iq,ib,1,1)
        !base%dnRGB%d1(iq,ib,1)   = d1u(1)  * base%dnRGB%d1(iq,ib,1)
      END DO
    END DO

    !-----------------------------------------------------------
    IF (debug) THEN
      CALL RecWrite_basis(base,write_all=.TRUE.)
      write(out_unitp,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------

  end subroutine sub_quadra_hermiteAB

  SUBROUTINE sub_quadra_hermite_half(base,paire)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      integer          :: paire
      logical          :: num
      real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical                       :: deriv
      integer                       :: i,iq,nq
      logical                       :: Print_basis
      TYPE (File_t)             :: herm_file
      integer                       :: idum,nio,err_io
      real(kind=Rkind)              :: wdum
      logical, parameter            :: WithRead = .TRUE.

      real(kind=Rkind), allocatable :: x_loc(:,:),w_loc(:)


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite_half'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
         write(out_unitp,*) 'nb',base%nb
         write(out_unitp,*) 'nq',nq
       END IF
!-----------------------------------------------------------
       deriv = .TRUE.
       num   = .FALSE.

       Print_basis = base%print_info_OF_basisDP .AND. print_level > -1

       IF (base%nb <= 0) THEN
         base%nb = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)
       END IF

       IF (base%nb <= 0) THEN
         write(out_unitp,*) 'ERROR in ',name_sub
         write(out_unitp,*) 'nb<=0',base%nb
         STOP 'ERROR in sub_quadra_hermite_half: nb<=0'
       END IF

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (Print_basis) THEN
        write(out_unitp,*) '    Basis: Hermite polynomials in [0,+inf]'
        write(out_unitp,*) '      nb_hermite',base%nb
      END IF

      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN


        IF (Print_basis) write(out_unitp,*) '      nb_hermite',base%nb
        IF (Print_basis) write(out_unitp,*) '      old nb_quadra',nq
        flush(out_unitp)

        IF (base%check_nq_OF_basis) THEN
          IF (nq < base%nb) nq = base%nb+1
        END IF
        IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq
        flush(out_unitp)

        !nb_nosym = 2*base%nb
        !nq_nosym = 2*nq




        CALL alloc_NParray(x_loc,[1,2*nq],'x_loc',name_sub)
        CALL alloc_NParray(w_loc,[2*nq],  'w_loc',name_sub)

        CALL hercom(2*nq,x_loc,w_loc)

        IF (WithRead) THEN
          herm_file%name = trim(EVRT_path) //                                   &
            '/Internal_data/HermQuadra/herm' // TO_string(2*nq) // '.txt'

          write(out_unitp,*) 'herm_file%name: ',herm_file%name
          CALL file_open(herm_file,nio,old=.TRUE.,err_file=err_io)

          IF (err_io /= 0) THEN
            write(out_unitp,*) 'WARNNING in ',name_sub
            write(out_unitp,*) 'cannot find quadrature file: ',herm_file%name
          ELSE

            DO iq=1,2*nq
              read(nio,*,iostat=err_io) idum,x_loc(:,iq),wdum,w_loc(iq)
            END DO

          END IF
          CALL file_close(herm_file)

        END IF
        !base%xPOGridRep_done = .TRUE.

        CALL Set_nq_OF_basis(base,nq+base%nq_extra)
        CALL alloc_xw_OF_basis(base)
        base%w(:) = ZERO

        base%x(:,1:nq) = x_loc(:,nq+1:2*nq) ! the positive part of the grid
        base%w(  1:nq) = w_loc(  nq+1:2*nq)

        IF (allocated(base%x_extra)) THEN
          base%x(:,nq+1:nq+base%nq_extra) = base%x_extra(:,:)
        END IF

        base%wrho(:) = base%w(:)
        base%rho(:)  = ONE

        IF (allocated(x_loc)) CALL dealloc_NParray(x_loc,'x_loc',name_sub)
        IF (allocated(w_loc)) CALL dealloc_NParray(w_loc,'w_loc',name_sub)
      END IF

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      nq = nq + base%nq_extra

      IF (paire == 0) THEN
        IF (Print_basis) write(out_unitp,*) '      even Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i-1,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_0_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nq,deriv,num,step)
      ELSE
        IF (Print_basis) write(out_unitp,*) '      odd Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_1_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nq,deriv,num,step)
      END IF
     
     ! because the integration is done on the positive part of the grid,
     !   the basis function must be multiply by sqrt(2)
     base%dnRGB%d0 = base%dnRGB%d0 * sqrt(TWO)
     base%dnRGB%d1 = base%dnRGB%d1 * sqrt(TWO)
     base%dnRGB%d2 = base%dnRGB%d2 * sqrt(TWO)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base,write_all=.TRUE.)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

end subroutine sub_quadra_hermite_half
SUBROUTINE sub_quadra_hermite_halfRight(base,paire)

  USE mod_system
  USE mod_basis
  IMPLICIT NONE

!---------- variables passees en argument ----------------------------
  TYPE (basis)     :: base
  integer          :: paire
  logical          :: num
  real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  logical                       :: deriv
  integer                       :: i,iq,nq
  logical                       :: Print_basis
  TYPE (File_t)             :: herm_file
  integer                       :: idum,nio,err_io
  real(kind=Rkind)              :: wdum
  logical, parameter            :: WithRead = .TRUE.

  real(kind=Rkind), allocatable :: x_loc(:,:),w_loc(:)


!----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_quadra_hermite_halfRight'
  !logical, parameter :: debug = .TRUE.
  logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
   nq = get_nq_FROM_basis(base)

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',name_sub
     write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
     write(out_unitp,*) 'nb',base%nb
     write(out_unitp,*) 'nq',nq
   END IF
!-----------------------------------------------------------
   deriv = .TRUE.
   num   = .FALSE.

   Print_basis = base%print_info_OF_basisDP .AND. print_level > -1

   IF (base%nb <= 0) THEN
     base%nb = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)
   END IF

   IF (base%nb <= 0) THEN
     write(out_unitp,*) 'ERROR in ',name_sub
     write(out_unitp,*) 'nb<=0',base%nb
     STOP 'ERROR in sub_quadra_hermite_halfRight: nb<=0'
   END IF

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
  IF (Print_basis) THEN
    write(out_unitp,*) '    Basis: Hermite polynomials in [0,+inf['
    write(out_unitp,*) '      nb_hermite',base%nb
  END IF

  base%packed            = .TRUE.
  base%packed_done       = .TRUE.


  IF (.NOT. base%xPOGridRep_done) THEN


    IF (Print_basis) write(out_unitp,*) '      nb_hermite',base%nb
    IF (Print_basis) write(out_unitp,*) '      old nb_quadra',nq
    flush(out_unitp)

    IF (base%check_nq_OF_basis) THEN
      IF (nq < base%nb) nq = base%nb+1
    END IF
    IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq
    flush(out_unitp)

    !nb_nosym = 2*base%nb
    !nq_nosym = 2*nq




    CALL alloc_NParray(x_loc,[1,2*nq],'x_loc',name_sub)
    CALL alloc_NParray(w_loc,[2*nq],  'w_loc',name_sub)

    CALL hercom(2*nq,x_loc,w_loc)

    IF (WithRead) THEN
      herm_file%name = trim(EVRT_path) //                                   &
        '/Internal_data/HermQuadra/herm' // TO_string(2*nq) // '.txt'

      write(out_unitp,*) 'herm_file%name: ',herm_file%name
      CALL file_open(herm_file,nio,old=.TRUE.,err_file=err_io)

      IF (err_io /= 0) THEN
        write(out_unitp,*) 'WARNNING in ',name_sub
        write(out_unitp,*) 'cannot find quadrature file: ',herm_file%name
      ELSE

        DO iq=1,2*nq
          read(nio,*,iostat=err_io) idum,x_loc(:,iq),wdum,w_loc(iq)
        END DO

      END IF
      CALL file_close(herm_file)

    END IF
    !base%xPOGridRep_done = .TRUE.

    CALL Set_nq_OF_basis(base,nq+base%nq_extra)
    CALL alloc_xw_OF_basis(base)
    base%w(:) = ZERO

    base%x(:,1:nq) = x_loc(:,nq+1:2*nq) ! the positive part of the grid
    base%w(  1:nq) = w_loc(  nq+1:2*nq)

    IF (allocated(base%x_extra)) THEN
      base%x(:,nq+1:nq+base%nq_extra) = base%x_extra(:,:)
    END IF

    base%wrho(:) = base%w(:)
    base%rho(:)  = ONE

    IF (allocated(x_loc)) CALL dealloc_NParray(x_loc,'x_loc',name_sub)
    IF (allocated(w_loc)) CALL dealloc_NParray(w_loc,'w_loc',name_sub)
  END IF

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
  deriv = .TRUE.
  CALL alloc_dnb_OF_basis(base)

  nq = nq + base%nq_extra

  IF (Print_basis) write(out_unitp,*) '      odd Hermite polynomials'
  base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
  CALL d0d1d2poly_Hermite_1_exp_grille(                                     &
                         base%x(1,:),                                       &
                         base%dnRGB%d0(:,:),                                &
                         base%dnRGB%d1(:,:,1),                              &
                         base%dnRGB%d2(:,:,1,1),                            &
                         base%nb,nq,deriv,num,step)

 ! because the integration is done on the positive part of the grid,
 !   the basis function must be multiply by sqrt(2)
 base%dnRGB%d0 = base%dnRGB%d0 * sqrt(TWO)
 base%dnRGB%d1 = base%dnRGB%d1 * sqrt(TWO)
 base%dnRGB%d2 = base%dnRGB%d2 * sqrt(TWO)

!-----------------------------------------------------------
  IF (debug) THEN
    CALL RecWrite_basis(base,write_all=.TRUE.)
    write(out_unitp,*) 'END ',name_sub
  END IF
!-----------------------------------------------------------

end subroutine sub_quadra_hermite_halfRight
SUBROUTINE sub_quadra_hermite_halfLeft_new(base,paire)

  USE mod_system
  USE mod_basis
  IMPLICIT NONE

!---------- variables passees en argument ----------------------------
  TYPE (basis)     :: base
  integer          :: paire
  logical          :: num
  real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  logical                       :: deriv
  integer                       :: i,iq,nq
  logical                       :: Print_basis
  TYPE (File_t)             :: herm_file
  integer                       :: idum,nio,err_io
  real(kind=Rkind)              :: wdum
  logical, parameter            :: WithRead = .TRUE.

  real(kind=Rkind), allocatable :: x_loc(:,:),w_loc(:)


!----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_quadra_hermite_halfLeft_new'
  !logical, parameter :: debug = .TRUE.
  logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
   nq = get_nq_FROM_basis(base)

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',name_sub
     write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
     write(out_unitp,*) 'nb',base%nb
     write(out_unitp,*) 'nq',nq
   END IF
!-----------------------------------------------------------
   deriv = .TRUE.
   num   = .FALSE.

   Print_basis = base%print_info_OF_basisDP .AND. print_level > -1

   IF (base%nb <= 0) THEN
     base%nb = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)
   END IF

   IF (base%nb <= 0) THEN
     write(out_unitp,*) 'ERROR in ',name_sub
     write(out_unitp,*) 'nb<=0',base%nb
     STOP 'ERROR in sub_quadra_hermite_halfLeft_new: nb<=0'
   END IF

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
  IF (Print_basis) THEN
    write(out_unitp,*) '    Basis: Hermite polynomials in ]-inf,0]'
    write(out_unitp,*) '      nb_hermite',base%nb
  END IF

  base%packed            = .TRUE.
  base%packed_done       = .TRUE.


  IF (.NOT. base%xPOGridRep_done) THEN


    IF (Print_basis) write(out_unitp,*) '      nb_hermite',base%nb
    IF (Print_basis) write(out_unitp,*) '      old nb_quadra',nq
    flush(out_unitp)

    IF (base%check_nq_OF_basis) THEN
      IF (nq < base%nb) nq = base%nb+1
    END IF
    IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq
    flush(out_unitp)

    !nb_nosym = 2*base%nb
    !nq_nosym = 2*nq




    CALL alloc_NParray(x_loc,[1,2*nq+1],'x_loc',name_sub)
    CALL alloc_NParray(w_loc,[2*nq+1],  'w_loc',name_sub)

    CALL hercom(2*nq+1,x_loc,w_loc)

    IF (WithRead) THEN
      herm_file%name = trim(EVRT_path) //                                   &
        '/Internal_data/HermQuadra/herm' // TO_string(2*nq+1) // '.txt'

      write(out_unitp,*) 'herm_file%name: ',herm_file%name
      CALL file_open(herm_file,nio,old=.TRUE.,err_file=err_io)

      IF (err_io /= 0) THEN
        write(out_unitp,*) 'WARNNING in ',name_sub
        write(out_unitp,*) 'cannot find quadrature file: ',herm_file%name
      ELSE

        DO iq=1,2*nq+1
          read(nio,*,iostat=err_io) idum,x_loc(:,iq),wdum,w_loc(iq)
        END DO

      END IF
      CALL file_close(herm_file)

    END IF
    !base%xPOGridRep_done = .TRUE.

    CALL Set_nq_OF_basis(base,nq+1+base%nq_extra)
    CALL alloc_xw_OF_basis(base)
    base%w(:) = ZERO

    base%x(:,1:nq+1) = x_loc(:,1:nq+1) ! the negative part of the grid
    base%w(  1:nq+1) = w_loc(  1:nq+1)
    write(6,*) 'x=0?',base%x(:,nq+1)

    IF (allocated(base%x_extra)) THEN
      base%x(:,nq+1:nq+base%nq_extra) = base%x_extra(:,:)
    END IF

    base%wrho(:) = base%w(:)
    base%rho(:)  = ONE

    IF (allocated(x_loc)) CALL dealloc_NParray(x_loc,'x_loc',name_sub)
    IF (allocated(w_loc)) CALL dealloc_NParray(w_loc,'w_loc',name_sub)
  END IF

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
  deriv = .TRUE.
  CALL alloc_dnb_OF_basis(base)

  nq = nq+1 + base%nq_extra

  IF (Print_basis) write(out_unitp,*) '      odd Hermite polynomials'
  base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
  CALL d0d1d2poly_Hermite_1_exp_grille(                                     &
                         base%x(1,:),                                       &
                         base%dnRGB%d0(:,:),                                &
                         base%dnRGB%d1(:,:,1),                              &
                         base%dnRGB%d2(:,:,1,1),                            &
                         base%nb,nq,deriv,num,step)

 ! because the integration is done on the positive part of the grid,
 !   the basis function must be multiply by sqrt(2)
 base%dnRGB%d0 = base%dnRGB%d0 * sqrt(TWO)
 base%dnRGB%d1 = base%dnRGB%d1 * sqrt(TWO)
 base%dnRGB%d2 = base%dnRGB%d2 * sqrt(TWO)

!-----------------------------------------------------------
  IF (debug) THEN
    CALL RecWrite_basis(base,write_all=.TRUE.)
    write(out_unitp,*) 'END ',name_sub
  END IF
!-----------------------------------------------------------

end subroutine sub_quadra_hermite_halfLeft_new
SUBROUTINE sub_quadra_hermite_halfLeft(base,paire)

  USE mod_system
  USE mod_basis
  IMPLICIT NONE

!---------- variables passees en argument ----------------------------
  TYPE (basis)     :: base
  integer          :: paire
  logical          :: num
  real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  logical                       :: deriv
  integer                       :: i,iq,nq
  logical                       :: Print_basis
  TYPE (File_t)             :: herm_file
  integer                       :: idum,nio,err_io
  real(kind=Rkind)              :: wdum
  logical, parameter            :: WithRead = .TRUE.

  real(kind=Rkind), allocatable :: x_loc(:,:),w_loc(:)


!----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_quadra_hermite_halfLeft'
  !logical, parameter :: debug = .TRUE.
  logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
   nq = get_nq_FROM_basis(base)

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',name_sub
     write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
     write(out_unitp,*) 'nb',base%nb
     write(out_unitp,*) 'nq',nq
   END IF
!-----------------------------------------------------------
   deriv = .TRUE.
   num   = .FALSE.

   Print_basis = base%print_info_OF_basisDP .AND. print_level > -1

   IF (base%nb <= 0) THEN
     base%nb = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)
   END IF

   IF (base%nb <= 0) THEN
     write(out_unitp,*) 'ERROR in ',name_sub
     write(out_unitp,*) 'nb<=0',base%nb
     STOP 'ERROR in sub_quadra_hermite_halfLeft: nb<=0'
   END IF

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
  IF (Print_basis) THEN
    write(out_unitp,*) '    Basis: Hermite polynomials in ]-inf,0]'
    write(out_unitp,*) '      nb_hermite',base%nb
  END IF

  base%packed            = .TRUE.
  base%packed_done       = .TRUE.


  IF (.NOT. base%xPOGridRep_done) THEN


    IF (Print_basis) write(out_unitp,*) '      nb_hermite',base%nb
    IF (Print_basis) write(out_unitp,*) '      old nb_quadra',nq
    flush(out_unitp)

    IF (base%check_nq_OF_basis) THEN
      IF (nq < base%nb) nq = base%nb+1
    END IF
    IF (Print_basis) write(out_unitp,*) '      new nb_quadra',nq
    flush(out_unitp)

    !nb_nosym = 2*base%nb
    !nq_nosym = 2*nq




    CALL alloc_NParray(x_loc,[1,2*nq],'x_loc',name_sub)
    CALL alloc_NParray(w_loc,[2*nq],  'w_loc',name_sub)

    CALL hercom(2*nq,x_loc,w_loc)

    IF (WithRead) THEN
      herm_file%name = trim(EVRT_path) //                                   &
        '/Internal_data/HermQuadra/herm' // TO_string(2*nq) // '.txt'

      write(out_unitp,*) 'herm_file%name: ',herm_file%name
      CALL file_open(herm_file,nio,old=.TRUE.,err_file=err_io)

      IF (err_io /= 0) THEN
        write(out_unitp,*) 'WARNNING in ',name_sub
        write(out_unitp,*) 'cannot find quadrature file: ',herm_file%name
      ELSE

        DO iq=1,2*nq
          read(nio,*,iostat=err_io) idum,x_loc(:,iq),wdum,w_loc(iq)
        END DO

      END IF
      CALL file_close(herm_file)

    END IF
    !base%xPOGridRep_done = .TRUE.

    CALL Set_nq_OF_basis(base,nq+base%nq_extra)
    CALL alloc_xw_OF_basis(base)
    base%w(:) = ZERO

    base%x(:,1:nq) = x_loc(:,1:nq) ! the negative part of the grid
    base%w(  1:nq) = w_loc(  1:nq)

    IF (allocated(base%x_extra)) THEN
      base%x(:,nq+1:nq+base%nq_extra) = base%x_extra(:,:)
    END IF

    base%wrho(:) = base%w(:)
    base%rho(:)  = ONE

    IF (allocated(x_loc)) CALL dealloc_NParray(x_loc,'x_loc',name_sub)
    IF (allocated(w_loc)) CALL dealloc_NParray(w_loc,'w_loc',name_sub)
  END IF

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
  deriv = .TRUE.
  CALL alloc_dnb_OF_basis(base)

  nq = nq + base%nq_extra

  IF (Print_basis) write(out_unitp,*) '      odd Hermite polynomials'
  base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
  CALL d0d1d2poly_Hermite_1_exp_grille(                                     &
                         base%x(1,:),                                       &
                         base%dnRGB%d0(:,:),                                &
                         base%dnRGB%d1(:,:,1),                              &
                         base%dnRGB%d2(:,:,1,1),                            &
                         base%nb,nq,deriv,num,step)

 ! because the integration is done on the positive part of the grid,
 !   the basis function must be multiply by sqrt(2)
 base%dnRGB%d0 = base%dnRGB%d0 * sqrt(TWO)
 base%dnRGB%d1 = base%dnRGB%d1 * sqrt(TWO)
 base%dnRGB%d2 = base%dnRGB%d2 * sqrt(TWO)

!-----------------------------------------------------------
  IF (debug) THEN
    CALL RecWrite_basis(base,write_all=.TRUE.)
    write(out_unitp,*) 'END ',name_sub
  END IF
!-----------------------------------------------------------

end subroutine sub_quadra_hermite_halfLeft
!=============================================================
!
!      determination des tous les Hn(xi)=poly_h(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_hermitebox(base,paire)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base
      integer paire
      logical num
      real(kind=Rkind)  step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical  :: deriv
      real(kind=Rkind) :: A,B,xeq,scale,dx
      integer  :: nq,i,j,nqo

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermitebox'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nqo = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_hermitebox'
         write(out_unitp,*) 'nb,nq',base%nb,nqo
         write(out_unitp,*) 'num,step',num,step
       END IF
!-----------------------------------------------------------

       deriv = .TRUE.
       num   = .FALSE.

       IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%print_info_OF_basisDP .AND. print_level > -1)            &
                 write(out_unitp,*) '    Basis: Hermite polynomia+BoxAB'
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in sub_quadra_hermitebox'
        write(out_unitp,*) 'xPOGridRep_done=t is not possible for this basis!'
        write(out_unitp,*) 'CHECK your data'
        STOP
      END IF
      IF (base%check_nq_OF_basis) THEN
        IF (print_level > -1) write(out_unitp,*) '      nb_hermite',base%nb
        IF (print_level > -1) write(out_unitp,*) '      old nb_quadra',nqo
        IF (paire .EQ. -1) THEN
          IF ( nqo < base%nb ) nqo = base%nb + 1
        ELSE
          IF ( nqo < 2*base%nb ) nqo = 2*base%nb+1
        END IF
        IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nqo
      END IF
      CALL Set_nq_OF_basis(base,nqo)

      CALL alloc_xw_OF_basis(base)
!     - The true Gauss-Hermite grid -----------------------------------------
!       with a small number of grid points (nb)
      nq = base%nb
      CALL hercom(nq,base%x(1,:),base%w)
!     - The true Gauss-Hermite grid -----------------------------------------


!     - We add grid points between A and B (regularly spaced) ---------------
!       without weight (0.)
!     ------ for the domain [A,B] -------------------------------------------
      A = base%x(1,1)
      B = base%x(1,nq)
      IF (base%print_info_OF_basisDP .AND. print_level > -1)            &
                                   write(out_unitp,*) 'HermBox: A,B',A,B
      dx = (B-A)/real(nqo-nq,kind=Rkind)
      base%x(1,nq+1:nqo) =                                              &
         [(A+dx*(real(i,kind=Rkind)-HALF),i=1,nqo-nq)]
      base%w(nq+1:nqo) = ZERO


!     - sort the grid points --------------------------------------------------
      DO i=1,nqo
      DO j=i+1,nqo
        IF (base%x(1,i) > base%x(1,j) ) THEN
          A           = base%x(1,j)
          base%x(1,j) = base%x(1,i)
          base%x(1,i) = A
          A         = base%w(j)
          base%w(j) = base%w(i)
          base%w(i) = A
        END IF
      END DO
      END DO
!     - We add grid points between A and B (regularly spaced) ---------------


      base%wrho(:) = base%w(:)
      base%rho(:)  = ONE

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN ! even
        IF (base%print_info_OF_basisDP .AND. print_level > -1)                  &
                       write(out_unitp,*) '      even Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i-1,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_0_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)
      ELSE IF (paire == 1) THEN ! odd
        IF (base%print_info_OF_basisDP .AND. print_level > -1)                  &
                        write(out_unitp,*) '      odd Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_1_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)
      ELSE
        IF (base%print_info_OF_basisDP .AND. print_level > -1)                  &
                        write(out_unitp,*) '      All Hermite polynomials'
        base%tab_ndim_index(1,:) = [(i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_exp_grille(                                     &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)
      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_hermite'
      END IF
!-----------------------------------------------------------

      end subroutine sub_quadra_hermitebox


      SUBROUTINE grid_HermiteNested1(x,w,nq,nqmax)

      USE mod_system
      USE mod_basis
      USE BasisMakeGrid
      IMPLICIT NONE

      integer, intent(in) :: nq,nqmax
      real(kind=Rkind), intent(inout) :: x(nq)
      real(kind=Rkind), intent(inout) :: w(nq)

!---------------------------------------------------------------------
!---------- local variables ----------------------------
      TYPE (basis)     :: base
      logical          :: num,deriv
      real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      !real(kind=Rkind) :: A,B,xeq,scale,dx
      !integer  :: nq,iq,i,j,nb_nosym,nq0,dnq,nqo
      real(kind=Rkind), allocatable :: x_loc(:)
      real(kind=Rkind), allocatable :: w_loc(:)
      integer  :: dnb


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='grid_HermiteNested1'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nq,nqmax',nq,nqmax
      END IF
!-----------------------------------------------------------

      deriv = .TRUE.
      num   = .FALSE.

      IF (nq > nqmax .OR. mod(nqmax-nq,2) /= 0) THEN
        write(out_unitp,*) 'nq,nqmax',nq,nqmax
        STOP 'ERROR in grid_HermiteNested1'
      END IF


      CALL alloc_NParray(x_loc,[nqmax],"x_loc",name_sub)
      CALL alloc_NParray(w_loc,[nqmax],"w_loc",name_sub)
      CALL hercom(nqmax,x_loc(:),w_loc(:))
      !write(out_unitp,*) 'old w(:)',w_loc(:)
      !write(out_unitp,*) 'old x(:)',x_loc(:)

      ! new grid
      dnb = (nqmax-nq)/2
      x(:) = x_loc(1+dnb:nqmax-dnb)


      !weight
      CALL Set_nq_OF_basis(base,nq)
      base%nb   = (nq+1)/2
      base%ndim = 1
      CALL alloc_dnb_OF_basis(base)
      CALL d0d1d2poly_Hermite_exp_grille(                                       &
                             x(:),                                              &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nq,deriv,num,step)

      CALL Weight_OF_grid(w,base%dnRGB%d0,base%nb,nq)

      CALL dealloc_dnb_OF_basis(base)

      CALL dealloc_NParray(x_loc,"x_loc",name_sub)
      CALL dealloc_NParray(w_loc,"w_loc",name_sub)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE grid_HermiteNested1

      SUBROUTINE sub_quadra_HermiteNested2(base,paire)

      USE mod_system
      USE mod_basis
      USE BasisMakeGrid
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base
      integer paire
      logical num
      real(kind=Rkind)  step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical  :: deriv
      real(kind=Rkind) :: A,B,xeq,scale,dx
      integer  :: nb,nq,iq,i,j,nb_nosym,nq0,dnq,nq1,nq2,nqo
      real(kind=Rkind), allocatable :: x_loc(:)
      real(kind=Rkind), allocatable :: w_loc(:)


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_quadra_HermiteNested2'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nqo = get_nq_FROM_basis(base)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb,nq',base%nb,nqo
        write(out_unitp,*) 'num,step',num,step
      END IF
!-----------------------------------------------------------


      deriv = .TRUE.
      num   = .FALSE.

      IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%print_info_OF_basisDP) write(out_unitp,*) '    Basis: Hermite polynomials+Nested'
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'xPOGridRep_done=t is not possible for this basis!'
        write(out_unitp,*) 'CHECK your data'
        STOP
      END IF
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '      nb_hermite',base%nb
        write(out_unitp,*) '      old nb_quadra',nqo
        IF (paire == -1) THEN
          nb_nosym = base%nb
        ELSE IF (paire == 1) THEN ! even
          nb_nosym = 2*base%nb - 1
        ELSE ! no sym
          nb_nosym = 2*base%nb
        END IF
        IF (mod(nqo,2) == 0) nqo = nqo + 1
        IF (nqo < 2*nb_nosym-1) THEN
          nqo = 2*nb_nosym-1
        END IF
        write(out_unitp,*) '      new nb_quadra',nqo
      END IF

      IF (base%L_SparseGrid < 0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'base%L_SparseGrid MUST larger than -1',base%L_SparseGrid
        write(out_unitp,*) 'CHECK the fortran !!'
        STOP
      END IF
      IF (base%nq_max_Nested < 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'nq_max_Nested < 1',base%nq_max_Nested
        write(out_unitp,*) 'CHECK your data or the fortran!!'
        STOP
      END IF
      IF (mod(base%nq_max_Nested,2) == 0) base%nq_max_Nested = base%nq_max_Nested + 1

      ! nested Hermite Grid (here without weight)
      IF (base%L_SparseGrid <= base%nq_max_Nested/2) THEN
        nqo      = base%L_SparseGrid * 2 + 1
        nb_nosym = base%L_SparseGrid + 1
        nq0      = base%nq_max_Nested/2 + 1
        dnq      = nqo/2
        nq1      = nq0-dnq
        nq2      = nq0+dnq
      ELSE
        nqo      = base%nq_max_Nested
        nb_nosym = base%nq_max_Nested
        nq1      = 1
        nq2      = base%nq_max_Nested
      END IF
      write(out_unitp,*) 'base%nq_max_Nested',base%nq_max_Nested
      write(out_unitp,*) 'base%L_SparseGrid,nq,nb_nosym,nq1,nq2',base%L_SparseGrid,nqo,nb_nosym,nq1,nq2
        CALL Set_nq_OF_basis(base,nqo)

      CALL alloc_xw_OF_basis(base)

      CALL alloc_NParray(x_loc,[base%nq_max_Nested],                        &
                        "x_loc","sub_quadra_HermiteNested2")
      CALL alloc_NParray(w_loc,[base%nq_max_Nested],                        &
                        "w_loc","sub_quadra_HermiteNested2")

      CALL hercom(base%nq_max_Nested,x_loc(:),w_loc(:))
      !write(out_unitp,*) 'old w(:)',w_loc(:)
      !write(out_unitp,*) 'old x(:)',x_loc(:)
      base%x(1,:) = x_loc(nq1:nq2)
      write(out_unitp,*) 'new x(:)',base%x(1,:)



      base%rho(:)  = ONE

      ! weights calculation
      nb = base%nb ! save nb
      base%nb = nb_nosym
      deriv = .FALSE.

      CALL alloc_dnb_OF_basis(base)
      CALL d0d1d2poly_Hermite_exp_grille(                                       &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)

      CALL Weight_OF_grid_basis(base)
      CALL dealloc_dnb_OF_basis(base)
      base%nb = nb
      write(out_unitp,*) 'new w',base%w

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN
        IF (base%print_info_OF_basisDP) write(out_unitp,*) '      even Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i-1,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_0_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)
      ELSE IF (paire == 1) THEN
        IF (base%print_info_OF_basisDP) write(out_unitp,*) '      odd Hermite polynomials'
        base%tab_ndim_index(1,:) = [(2*i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_1_exp_grille(                                   &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)
      ELSE
        IF (base%print_info_OF_basisDP) write(out_unitp,*) '      All Hermite polynomials'
        base%tab_ndim_index(1,:) = [(i,i=1,base%nb)]
        CALL d0d1d2poly_Hermite_exp_grille(                                     &
                             base%x(1,:),                                       &
                             base%dnRGB%d0(:,:),                                &
                             base%dnRGB%d1(:,:,1),                              &
                             base%dnRGB%d2(:,:,1,1),                            &
                             base%nb,nqo,deriv,num,step)
      END IF

      CALL dealloc_NParray(x_loc,"x_loc","sub_quadra_HermiteNested2")
      CALL dealloc_NParray(w_loc,"w_loc","sub_quadra_HermiteNested2")

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_HermiteNested2

     SUBROUTINE sub_quadra_hermite_cuba(base)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      logical          :: num = .FALSE.
      real(kind=Rkind) :: step = ZERO
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical          :: deriv = .TRUE.
      integer          :: p,d,option
      integer          :: tab_nb(base%ndim)
      real(kind=Rkind) :: x(base%ndim)
      real(kind=Rkind) :: d0(base%ndim)
      real(kind=Rkind) :: d1(base%ndim)
      real(kind=Rkind) :: d2(base%ndim)
      integer          :: i,iB,iQ,id,jd,nq,LB,LG


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite_cuba'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (base%L_SparseBasis < 0) THEN
         LB = base%nb-1
       ELSE
         LB = base%L_SparseBasis
       END IF

       IF (base%L_SparseGrid < 0) THEN
         LG = nq-1
       ELSE
         LG = base%L_SparseGrid
       END IF

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_hermite_cuba'
         write(out_unitp,*) 'LB,LG',LB,LG
         write(out_unitp,*) 'ndim',base%ndim
       END IF
!-----------------------------------------------------------
       IF (base%ndim < 1) STOP ' Problem with ndim'

       IF (base%ndim == 1) STOP ' You MUST use "HO" basis set'

       ! ndim is the "n" parameter of Burkardt subroutine
       ! d is the largest degree of the polynomials: d=2*(nq-1) = 2*LG
       IF (LB < 0) STOP 'ERROR LB<=0'

       ! test if nq >= nb
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: Hermite polynomials, cubature'
        write(out_unitp,*) '      LB,LG',LB,LG
      END IF
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old LG',LG
          IF ( LG < LB ) LG = LB + 1
          IF (print_level > -1) write(out_unitp,*) '      new LG',LG
          IF (print_level > -1) write(out_unitp,*) '           d',2*LG
        END IF

        d = 2*LG ! degree
        ! calculation of the grid point number (with the EN_02... size subroutine)
        SELECT CASE (d)
        CASE(0)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_01_1_size'
          CALL en_r2_01_1_size(base%ndim,nq)
        CASE(2)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_02_xiu_size'
          CALL en_r2_02_xiu_size(base%ndim,nq)
        CASE(4)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_05_1_size'
          option = 1
          CALL en_r2_05_1_size(base%ndim,option,nq)
        CASE(6)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_07_3_size'
          option = 1
          CALL en_r2_07_3_size(base%ndim,option,nq)
        CASE(8)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_09_1_size'
          option = 1
          CALL en_r2_09_1_size(base%ndim,option,nq)
        CASE(10)
         IF (print_level > -1)  write(out_unitp,*) ' en_r2_11_1_size'
          option = 1
          CALL en_r2_11_1_size(base%ndim,option,nq)
        CASE DEFAULT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' we cannot use cubature with d > 10'
          write(out_unitp,*) 'd',d
          write(out_unitp,*) 'nb,nq',base%nb,nq
          write(out_unitp,*) 'ndim',base%ndim
          nq = -1
          !RETURN
        END SELECT
        IF (base%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) 'maximum degree:',d
          write(out_unitp,*) 'Cubature point number:',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)


        IF (nq < 1) RETURN
        CALL alloc_xw_OF_basis(base)


        ! the grid points
        SELECT CASE (d)
        CASE(0)
          CALL en_r2_01_1(base%ndim,nq,base%x,base%w)
        CASE(2)
          CALL en_r2_02_xiu(base%ndim,nq,base%x,base%w)
        CASE(4)
          option = 1
          CALL en_r2_05_1(base%ndim,option,nq,base%x,base%w)
        CASE(6)
          option = 1
          CALL en_r2_07_3(base%ndim,option,nq,base%x,base%w)
        CASE(8)
          option = 1
          CALL en_r2_09_1(base%ndim,option,nq,base%x,base%w)
        CASE(10)
          option = 1
          CALL en_r2_11_1(base%ndim,option,nq,base%x,base%w)
        CASE DEFAULT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' we cannot use cubature with d > 10'
          write(out_unitp,*) 'nb,nq',base%nb,nq
          write(out_unitp,*) 'ndim',base%ndim
          nq = -1
          STOP
        END SELECT

        DO iQ=1,nq
          x(:) = base%x(:,iQ)
          base%w(iQ) = base%w(iQ) / exp(-dot_product(x,x))
        END DO

        IF (base%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) 'cubature',nq
          DO iQ=1,nq
            write(out_unitp,*) iQ,base%x(:,iQ),base%w(iQ)
          END DO
        END IF

        base%wrho(:) = base%w(:)
        !base%xPOGridRep_done = .TRUE.
      END IF
      base%rho(:)  = ONE



      ! number of basis functions:

      CALL dealloc_nDindex(base%nDindB)
      CALL alloc_nDindex(base%nDindB,ndim=base%ndim)

      base%nDindB%packed = .TRUE.
      tab_nb(:) = LB+1
      CALL init_nDindexPrim(base%nDindB,ndim=base%ndim,                 &
                            Type_OF_nDindex=0,                          &
                            MaxNorm=real(LB,kind=Rkind),                &
                            nDsize=tab_nb)
      IF (debug) CALL Write_nDindex(base%nDindB,name_sub)
      base%nb = base%nDindB%max_nDI
      IF (base%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) 'Basis function number:',base%nb
      END IF


!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)


      DO iQ=1,nq
        x(:) = base%x(:,iQ)
        DO iB=1,base%nb
          CALL calc_nDindex(base%nDindB,iB,tab_nb)

          base%tab_ndim_index(:,iB) = tab_nb(:)


          DO i=1,base%ndim
            CALL d0d1d2poly_Hermite_exp(x(i),tab_nb(i)-1,                       &
                                       d0(i),d1(i),d2(i),deriv,num,step)
          END DO
          base%dnRGB%d0(iQ,iB) = product(d0)

          DO id=1,base%ndim
            base%dnRGB%d1(iQ,iB,id)    = d1(id)
            base%dnRGB%d2(iQ,iB,id,id) = d2(id)
            DO i=1,base%ndim
              IF (i == id) CYCLE
              base%dnRGB%d1(iQ,iB,id)    = base%dnRGB%d1(iQ,iB,id)    * d0(i)
              base%dnRGB%d2(iQ,iB,id,id) = base%dnRGB%d2(iQ,iB,id,id) * d0(i)
            END DO
          END DO

          DO id=1,base%ndim
          DO jd=id+1,base%ndim

            base%dnRGB%d2(iQ,iB,id,jd) = d1(id) * d1(jd)
            DO i=1,base%ndim
              IF (i == id .OR. i ==jd) CYCLE
              base%dnRGB%d2(iQ,iB,id,jd) = base%dnRGB%d2(iQ,iB,id,jd) * d0(i)
            END DO
            base%dnRGB%d2(iQ,iB,jd,id) = base%dnRGB%d2(iQ,iB,id,jd)
          END DO
          END DO

        END DO
      END DO

!      CALL sort_basis(base)
!      base%check_basis = .TRUE.
!      write(out_unitp,*) 'coucou ',name_sub
!      CALL check_ortho_basis(base,.FALSE.)
!      base%check_basis = .FALSE.
!      write(out_unitp,*) 'old w',base%w
!      CALL Weight_OF_grid_basis(base)
!      write(out_unitp,*) 'new w',base%w

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_hermite_cuba'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_hermite_cuba
      SUBROUTINE sub_quadra_hermite_cuba_DML(base,err_grid)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      logical          :: num = .FALSE.
      real(kind=Rkind) :: step = ZERO
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical          :: deriv = .TRUE.
      integer          :: p,d,option
      integer          :: tab_nb(base%ndim)
      real(kind=Rkind) :: x(base%ndim)
      real(kind=Rkind) :: d0(base%ndim)
      real(kind=Rkind) :: d1(base%ndim)
      real(kind=Rkind) :: d2(base%ndim)
      integer          :: i,iB,iQ,id,jd,idum,nq


      integer                  :: LB,LG
      integer                  :: nio,err_io
      TYPE (File_t)            :: cubature_file
      logical                  :: err_grid

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite_cuba_DML'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (base%L_SparseBasis < 0) THEN
         LB = base%nb-1
       ELSE
         LB = base%L_SparseBasis
       END IF

       IF (base%L_SparseGrid < 0) THEN
         LG = nq-1
       ELSE
         LG = base%L_SparseGrid
       END IF

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) '    Basis: Hermite polynomials, cubture'
         write(out_unitp,*) '      LB,LG',LB,LG
       END IF
!-----------------------------------------------------------
       err_Grid = .FALSE.

       IF (base%ndim < 1) STOP ' Problem with ndim'

       IF (base%ndim == 1) STOP ' You MUST use "HO" basis set'

       ! ndim is the "n" parameter of Burkardt subroutine
       ! d is the largest degree of the polynomials: d=2*(nq-1)
       IF (LB < 0) STOP 'ERROR LB<0'

       ! test if nq >= nb
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: Hermite polynomia, cubature'
        write(out_unitp,*) '      LB,LG',LB,LG
      END IF
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old LG',LG
          IF (LG < LB) LG = LB + 1
          IF (print_level > -1) write(out_unitp,*) '      new LG',LG
        END IF

        IF (LG /= 0) THEN
          cubature_file%name = trim(EVRT_path) // '/Internal_data/HermCuba/' //  &
             TO_String(base%ndim) // 'D_deg' // TO_String(LG)


          write(out_unitp,*) 'cubature_file%name: ',cubature_file%name
          CALL file_open(cubature_file,nio,old=.TRUE.,err_file=err_io)

          IF (err_io == 0) THEN
            read(nio,*,iostat=err_io) nq
          ELSE
            write(out_unitp,*) 'WARNNING in ',name_sub
            write(out_unitp,*) 'cannot find cubature file: ',cubature_file%name
            err_Grid = .TRUE.
          END IF
          write(out_unitp,*) '      cubature nq',nq
          CALL file_close(cubature_file)

          CALL Set_nq_OF_basis(base,nq)


          err_Grid = err_Grid .OR. (nq < 1)
          IF (err_Grid) RETURN

          CALL alloc_xw_OF_basis(base)

          CALL file_open(cubature_file,nio,old=.TRUE.,err_file=err_io)
          read(nio,*,iostat=err_io)

          DO iq=1,nq
            read(nio,*,iostat=err_io) idum,base%x(:,iq),base%w(iq)
          END DO
          CALL file_close(cubature_file)

        ELSE ! LG = 0 (degree zero => only one point (0., 0., ...0.)
          nq = 1
          CALL Set_nq_OF_basis(base,nq)

          CALL alloc_xw_OF_basis(base)
          base%x(:,1) = ZERO
          base%w(1)   = sqrt(pi)**base%ndim
        END IF
        base%wrho(:) = base%w(:)
        !base%xPOGridRep_done = .TRUE.
      END IF
      base%rho(:)  = ONE


      ! number of basis functions:

      CALL dealloc_nDindex(base%nDindB)
      CALL alloc_nDindex(base%nDindB,ndim=base%ndim)

      base%nDindB%packed = .TRUE.
      tab_nb(:) = LB+1
      CALL init_nDindexPrim(base%nDindB,ndim=base%ndim,                         &
                            Type_OF_nDindex=0,                                  &
                            MaxNorm=real(LB,kind=Rkind),                        &
                            nDsize=tab_nb)
      IF (debug) CALL Write_nDindex(base%nDindB,name_sub)
      base%nb = base%nDindB%max_nDI
      IF (base%print_info_OF_basisDP) write(out_unitp,*) 'Basis function number:',base%nb


!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)


      DO iQ=1,nq
        x(:) = base%x(:,iQ)
        DO iB=1,base%nb
          CALL calc_nDindex(base%nDindB,iB,tab_nb)

          base%tab_ndim_index(:,iB) = tab_nb(:)


          DO i=1,base%ndim
            CALL d0d1d2poly_Hermite_exp(x(i),tab_nb(i)-1,                       &
                                       d0(i),d1(i),d2(i),deriv,num,step)
          END DO
          base%dnRGB%d0(iQ,iB) = product(d0)

          DO id=1,base%ndim
            base%dnRGB%d1(iQ,iB,id)    = d1(id)
            base%dnRGB%d2(iQ,iB,id,id) = d2(id)
            DO i=1,base%ndim
              IF (i == id) CYCLE
              base%dnRGB%d1(iQ,iB,id)    = base%dnRGB%d1(iQ,iB,id)    * d0(i)
              base%dnRGB%d2(iQ,iB,id,id) = base%dnRGB%d2(iQ,iB,id,id) * d0(i)
            END DO
          END DO

          DO id=1,base%ndim
          DO jd=id+1,base%ndim

            base%dnRGB%d2(iQ,iB,id,jd) = d1(id) * d1(jd)
            DO i=1,base%ndim
              IF (i == id .OR. i ==jd) CYCLE
              base%dnRGB%d2(iQ,iB,id,jd) = base%dnRGB%d2(iQ,iB,id,jd) * d0(i)
            END DO
            base%dnRGB%d2(iQ,iB,jd,id) = base%dnRGB%d2(iQ,iB,id,jd)
          END DO
          END DO

        END DO
      END DO

!      CALL sort_basis(base)
!      base%check_basis = .TRUE.
!      write(out_unitp,*) 'coucou ',name_sub
!      CALL check_ortho_basis(base,.FALSE.)
!      base%check_basis = .FALSE.
!      write(out_unitp,*) 'old w',base%w
!      CALL Weight_OF_grid_basis(base)
!      write(out_unitp,*) 'new w',base%w

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_hermite_cuba_DML'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_hermite_cuba_DML

      SUBROUTINE transfo_Q_TO_tQ(base)
      USE mod_system
      USE mod_basis
      USE mod_dnSVM
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base


      TYPE (Type_dnS)    :: dntQ,dntQ_inv
      integer            :: i,dnErr
      real(kind=Rkind)   :: q,d1,d2
      integer, parameter :: itype = 2

      CALL alloc_dnS(dntQ,nb_var_deriv=1,nderiv=3)
      CALL alloc_dnS(dntQ_inv,nb_var_deriv=1,nderiv=3)

      IF (.NOT. allocated(base%cte_Transfo) ) THEN
        CALL alloc_NParray(base%cte_Transfo, [20,1],'base%cte_Transfo','transfo_Q_TO_tQ')
      END IF

      DO i=1,get_nq_FROM_basis(base)
        q  = base%x(1,i)
        CALL sub_dntf(-itype,dntQ_inv,q,base%cte_Transfo(:,1),dnErr)
        IF (dnErr /= 0) THEN
          write(out_unitp,*) ' ERROR in transfo_Q_TO_tQ'
          write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, iQdyn:',base%iQdyn(1)
          STOP 'ERROR in sub_dntf called from transfo_Q_TO_tQ'
        END IF

        base%x(1,i) = dntQ_inv%d0

        CALL sub_dntf(itype,dntQ,dntQ_inv%d0,base%cte_Transfo(:,1),dnErr)
        IF (dnErr /= 0) THEN
          write(out_unitp,*) ' ERROR in transfo_Q_TO_tQ'
          write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, iQdyn:',base%iQdyn(1)
          STOP 'ERROR in sub_dntf called from transfo_Q_TO_tQ'
        END IF

        d1 = dntQ%d1(1)
        d2 = dntQ%d2(1,1)


        base%w(i)   = base%w(i)*base%rho(i)/d1
        base%rho(i) = d1
!       the wrho(:) is unchanged  !!
!       but nrho(1) has to be changed
        base%nrho(1) = itype  ! it means, the variable, x, is substituted by cos(theta)

        base%dnRGB%d2(i,:,1,1) = d2 * base%dnRGB%d1(i,:,1) + d1*d1 * base%dnRGB%d2(i,:,1,1)
        base%dnRGB%d1(i,:,1)   = d1 * base%dnRGB%d1(i,:,1)
      END DO

      CALL dealloc_dnS(dntQ)
      CALL dealloc_dnS(dntQ_inv)


      END SUBROUTINE transfo_Q_TO_tQ
