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
      MODULE mod_BFGS
      USE EVR_system_m

      IMPLICIT NONE

      TYPE param_BFGS

      integer           :: max_iteration        = 10

      real (kind=Rkind) :: max_grad             = 0.000015_Rkind
      real (kind=Rkind) :: RMS_grad             = 0.000010_Rkind

      real (kind=Rkind) :: max_step             = 0.000060_Rkind
      real (kind=Rkind) :: RMS_step             = 0.000040_Rkind
      real (kind=Rkind) :: largest_step         = 0.5_Rkind

      logical           :: calc_hessian         = .FALSE.
      logical           :: read_hessian         = .FALSE.
      logical           :: calc_hessian_always  = .FALSE.
      integer           :: calc_hessian_every_n = huge(1)
      real (kind=Rkind), allocatable :: hessian_inv_init(:,:)


      integer           :: nb_neg               = 0 ! number of negative hessian eigenvalues (calc_hessian_always=t)
      integer           :: max_Q_ForPrinting    = 11

      real (kind=Rkind), allocatable :: TS_vector(:)


      END TYPE param_BFGS

      CONTAINS

      SUBROUTINE Read_param_BFGS(para_BFGS,nb_Opt)
      TYPE (param_BFGS), intent(inout) :: para_BFGS
      integer,           intent(in)    :: nb_Opt

      integer           :: max_iteration,max_Q_ForPrinting

      real (kind=Rkind) :: max_grad
      real (kind=Rkind) :: RMS_grad

      real (kind=Rkind) :: max_step
      real (kind=Rkind) :: RMS_step
      real (kind=Rkind) :: largest_step
      logical           :: calc_hessian,read_hessian,calc_hessian_always
      integer           :: calc_hessian_every_n
      integer           :: nb_neg
      real (kind=Rkind) :: Norm2,TS_vector(nb_Opt)


      integer :: err_io
      NAMELIST /BFGS/ max_iteration,max_grad,RMS_grad,max_step,RMS_step,largest_step, &
                      calc_hessian,read_hessian,calc_hessian_every_n,calc_hessian_always, &
                      nb_neg,TS_vector,max_Q_ForPrinting

        max_iteration       = 10
        nb_neg              = 0
        TS_vector           = ZERO

        max_grad            = 0.000015_Rkind
        RMS_grad            = 0.000010_Rkind

        max_step            = 0.000060_Rkind
        RMS_step            = 0.000040_Rkind
        largest_step        = 0.5_Rkind
        calc_hessian        = .FALSE.
        read_hessian        = .FALSE.
        calc_hessian_always = .FALSE.
        calc_hessian_every_n = -1
        max_Q_ForPrinting   = 11

        read(in_unit,BFGS,IOSTAT=err_io)
        IF (err_io < 0) THEN
           write(out_unit,*) ' ERROR in Read_param_BFGS'
           write(out_unit,*) '  while reading the "BFGS" namelist'
           write(out_unit,*) ' end of file or end of record'
           write(out_unit,*) ' Check your data !!'
           STOP ' ERROR in Read_param_BFGS while reading the "BFGS" namelist'
        END IF
        IF (err_io > 0) THEN
          write(out_unit,*) ' ERROR in Read_param_BFGS'
          write(out_unit,*) '  while reading the "BFGS" namelist'
          write(out_unit,*) ' Probably you are using a wrong keyword'
          write(out_unit,*) ' Check your data !!'
          STOP ' ERROR in Read_param_BFGS while reading the "BFGS" namelist'
       END IF
        IF (print_level > 1) write(out_unit,BFGS)

        IF (calc_hessian_every_n == 0) THEN
          write(out_unit,*) ' ERROR in Read_param_BFGS'
          write(out_unit,*) '  calc_hessian_every_n MUST be > 0'
          write(out_unit,*) '  calc_hessian_every_n ',calc_hessian_every_n
          write(out_unit,*) ' Check your data !!'
          STOP 'STOP in Read_param_BFGS: calc_hessian_every_n MUST be > 0'
        END IF
        IF (calc_hessian_always) THEN
          IF (calc_hessian_every_n > 1) THEN
            write(out_unit,*) ' ERROR in Read_param_BFGS'
            write(out_unit,*) '  Inconsistent parameters values:'
            write(out_unit,*) '  calc_hessian_always',calc_hessian_always
            write(out_unit,*) '  calc_hessian_every_n > 1',calc_hessian_every_n
            write(out_unit,*) ' Check your data !!'
            STOP 'STOP in Read_param_BFGS: Inconsistent parameters values, calc_hessian_always and calc_hessian_every_n'
          ELSE ! IF (calc_hessian_every_n < 0) THEN
            calc_hessian_every_n = 1
          END IF
        ELSE
          IF (calc_hessian_every_n < 0) calc_hessian_every_n = huge(1)
        END IF

        IF (nb_neg > 0 .AND. calc_hessian_every_n /= 1) THEN
          write(out_unit,*) ' ERROR in Read_param_BFGS'
          write(out_unit,*) '  For nb_neg > 0, calc_hessian_always MUST be .TRUE. or calc_hessian_every_n=1'
          write(out_unit,*) '  nb_neg             ',nb_neg
          write(out_unit,*) '  calc_hessian_always',calc_hessian_always
          write(out_unit,*) '  calc_hessian_every_n',calc_hessian_every_n
          write(out_unit,*) ' Check your data !!'
          STOP 'STOP in Read_param_BFGS: inconsistent calc_hessian_every_n and nb_neg values.'
        END IF
        IF (calc_hessian .AND. read_hessian) THEN
          write(out_unit,*) ' ERROR in Read_param_BFGS'
          write(out_unit,*) '  Both calc_hessian and read_hessian are .TRUE., it not possible'
          write(out_unit,*) ' Check your data !!'
          STOP 'STOP in Read_param_BFGS:  Both calc_hessian and read_hessian are .TRUE., it not possible.'
        END IF
        para_BFGS%max_iteration       = max_iteration

        para_BFGS%max_grad            = max_grad
        para_BFGS%RMS_grad            = RMS_grad

        para_BFGS%max_step            = max_step
        para_BFGS%RMS_step            = RMS_step

        para_BFGS%largest_step        = largest_step


        para_BFGS%calc_hessian        = calc_hessian
        para_BFGS%read_hessian        = read_hessian
        para_BFGS%calc_hessian_every_n = calc_hessian_every_n
        !para_BFGS%calc_hessian_always = calc_hessian_always

        para_BFGS%nb_neg              = nb_neg

        Norm2 = dot_product(TS_Vector,TS_Vector)
        IF (Norm2 > ONETENTH**6) THEN
          para_BFGS%TS_Vector = TS_Vector / sqrt(Norm2)
        ELSE
          para_BFGS%TS_Vector = TS_Vector
          para_BFGS%TS_Vector = ZERO
        END IF

      END SUBROUTINE Read_param_BFGS

      SUBROUTINE Write_param_BFGS(para_BFGS)
      TYPE (param_BFGS), intent(in)   :: para_BFGS

      write(out_unit,*) 'WRITE param_BFGS'
      write(out_unit,*)
      write(out_unit,*) '  max_iteration  ',para_BFGS%max_iteration
      write(out_unit,*)
      write(out_unit,*) '  max_grad       ',para_BFGS%max_grad
      write(out_unit,*) '  RMS_grad       ',para_BFGS%RMS_grad
      write(out_unit,*) '  max_step       ',para_BFGS%max_step
      write(out_unit,*) '  RMS_step       ',para_BFGS%RMS_step
      write(out_unit,*) '  largest_step   ',para_BFGS%largest_step

      write(out_unit,*)
      write(out_unit,*) '  calc_hessian (initialization)',para_BFGS%calc_hessian
      write(out_unit,*) '  read_hessian (initialization)',para_BFGS%read_hessian
      write(out_unit,*) '  calc_hessian_every_n',para_BFGS%calc_hessian_every_n
      write(out_unit,*) '  calc_hessian_always ',(para_BFGS%calc_hessian_every_n == 1)
      IF (para_BFGS%calc_hessian_every_n == 1)  write(out_unit,*) '  BFGS procedure is forced to Newton-Raphson one'

      write(out_unit,*) '  nb_neg (=1 => TS)',para_BFGS%nb_neg
      write(out_unit,*) 'END WRITE param_BFGS'

  END SUBROUTINE Write_param_BFGS
  SUBROUTINE dealloc_param_BFGS(para_BFGS)
        TYPE (param_BFGS), intent(inout)   :: para_BFGS
  
        para_BFGS%max_iteration        = 10

        para_BFGS%max_grad             = 0.000015_Rkind
        para_BFGS%RMS_grad             = 0.000010_Rkind
        para_BFGS%max_step             = 0.000060_Rkind
        para_BFGS%RMS_step             = 0.000040_Rkind
        para_BFGS%largest_step         = 0.5_Rkind
  
        para_BFGS%calc_hessian         = .FALSE.
        para_BFGS%read_hessian         = .FALSE.
        para_BFGS%calc_hessian_always  = .FALSE.
        para_BFGS%calc_hessian_every_n = huge(1)

        para_BFGS%nb_neg               = 0 ! number of negative hessian eigenvalues (calc_hessian_always=t)
        para_BFGS%max_Q_ForPrinting    = 11

        IF (allocated(para_BFGS%TS_vector))        deallocate(para_BFGS%TS_vector)
        IF (allocated(para_BFGS%hessian_inv_init)) deallocate(para_BFGS%hessian_inv_init)

  END SUBROUTINE dealloc_param_BFGS
  SUBROUTINE Sub_BFGS(BasisnD,xOpt_min,SQ,nb_Opt,                   &
                          para_Tnum,mole,PrimOp,Qact,para_BFGS)

      USE EVR_system_m
      USE mod_dnSVM
      use mod_Coord_KEO, only: CoordType, tnum, alloc_array, dealloc_array
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save
      real (kind=Rkind), intent(inout) :: Qact(:)


!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: BasisnD

!----- variables pour la namelist minimum ----------------------------
      TYPE (PrimOp_t)  :: PrimOp
      integer          :: nb_scalar_Op
      logical          :: calc_scalar_Op

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t) :: para_ReadOp
      logical         :: Save_FileGrid,Save_MemGrid


!----- local variables -----------------------------------------------
!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)          :: para_AllOp_loc

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis)       :: para_AllBasis_loc
      TYPE (basis)                :: basis_temp

!----- for the optimization -------------------------------------------
      TYPE (param_BFGS) :: para_BFGS
      integer, intent(in) :: nb_Opt
      real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt)
      real (kind=Rkind), allocatable :: grad(:),hessian(:,:)

!---------- working variables ----------------------------------------
  TYPE (param_dnMatOp) :: dnMatOp(1)
  integer              :: nderiv_alloc,nbcol,err

!---------------------------------------------------------------------
!      logical,parameter :: debug= .FALSE.
      logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Sub_BFGS'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'xopt_min',xopt_min(:)
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------

        Qact(1:nb_Opt) = xopt_min(:)


        write(out_unit,*) 'Qact',Qact
        !---------JML-------------------------------------
        write(out_unit,*) 'mole%nb_act',mole%nb_act
        write(out_unit,*) 'RMS_grad',para_BFGS%RMS_grad
        write(out_unit,*) 'RMS_step',para_BFGS%RMS_step
        write(out_unit,*) 'max_iteration',para_BFGS%max_iteration

        IF (para_BFGS%calc_hessian) THEN
          write(out_unit,*) ' The initial hessian is calculated'
          !-------- allocation -----------------------------------------------
          CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_Opt,PrimOp%nb_elec,nderiv=2)
          CALL alloc_NParray(para_BFGS%hessian_inv_init,[nb_Opt,nb_Opt],  &
                            'para_BFGS%hessian_inv_init',name_sub)

          IF (allocated(hessian)) THEN
            CALL dealloc_NParray(hessian,'hessian',name_sub)
          END IF
          CALL alloc_NParray(hessian,[nb_Opt,nb_Opt],'hessian',name_sub)

          IF (allocated(grad)) THEN
            CALL dealloc_NParray(grad,'grad',name_sub)
          END IF
          CALL alloc_NParray(grad,[nb_Opt],'grad',name_sub)

          !-------- end allocation --------------------------------------------

          !----- Hessian and gradient ------------------------------------
          CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp)

          write(out_unit,*) 'Energy',Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)

          CALL Get_Grad_FROM_Tab_OF_dnMatOp(grad,dnMatOp)
          write(out_unit,*) 'grad'
          CALL Write_Vec_MPI(grad,out_unit,5)
          CALL dealloc_NParray(grad,'grad',name_sub)

          CALL Get_Hess_FROM_Tab_OF_dnMatOp(hessian,dnMatOp) ! for the ground state
          !----- End Hessian ------------------------------------

          !-------- deallocation ---------------------------------------------
          CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
          !-------- end deallocation -----------------------------------------

        ELSE IF (para_BFGS%read_hessian) THEN
          write(out_unit,*) ' The initial hessian is read in internal coordinates'
          !-------- allocation -----------------------------------------------
          CALL alloc_NParray(para_BFGS%hessian_inv_init,[nb_Opt,nb_Opt],  &
                            'para_BFGS%hessian_inv_init',name_sub)
          IF (allocated(hessian)) THEN
            CALL dealloc_NParray(hessian,'hessian',name_sub)
          END IF
          CALL alloc_NParray(hessian,[nb_Opt,nb_Opt],'hessian',name_sub)
          !-------- end allocation --------------------------------------------

          !----- Hessian ------------------------------------
          read(in_unit,*,IOSTAT=err) nbcol
          IF (err /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' "End of file", while reading nbcol'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF (print_level > 1) write(out_unit,*)'nbcol=',nbcol

          CALL Read_Mat(hessian,in_unit,nbcol,err)
          IF (err /= 0) THEN
            write(out_unit,*) 'ERROR ',name_sub
            write(out_unit,*) ' while reading the matrix "hessian"'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
          !----- End Hessian ------------------------------------

        END IF

        IF (allocated(hessian)) THEN
          IF (print_level > 1 .OR. nb_Opt <= para_BFGS%max_Q_ForPrinting) THEN
            write(out_unit,*) 'hessian'
            CALL Write_Mat_MPI(hessian,out_unit,5)
          ELSE
            write(out_unit,*) 'The hessian is too large. => No printing'
          END IF

          para_BFGS%hessian_inv_init = inv_OF_Mat_TO(hessian)

          IF (print_level > 1 .OR. nb_Opt <= para_BFGS%max_Q_ForPrinting) THEN
            write(out_unit,*) 'hessian inverse'
            CALL Write_Mat_MPI(para_BFGS%hessian_inv_init,out_unit,5)
          ELSE
            write(out_unit,*) 'The hessian inverse is too large. => No printing'
          END IF

          !-------- deallocation ---------------------------------------------
          CALL dealloc_NParray(hessian,'hessian',name_sub)
          !-------- end deallocation -----------------------------------------

        END IF

        IF (para_BFGS%max_iteration < 1) STOP 'STOP in Sub_BFGS: max_iteration < 1'
        !-------- allocation -----------------------------------------------
        CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_Opt,PrimOp%nb_elec,nderiv=1)
        !-------- end allocation --------------------------------------------

        CALL dfpmin_new(Qact,dnMatOp, mole,PrimOp,para_Tnum,para_BFGS,&
          para_BFGS%RMS_grad,para_BFGS%RMS_step,para_BFGS%max_iteration)

        xopt_min(:) = Qact(1:nb_Opt)

        !---------JMLend-------------------------------------
        !--------------------------------------------------

        !--------------------------------------------------
        ! this subroutine print the  matrix of derived type.
        CALL Write_MatOFdnS(dnMatOp(1)%tab_dnMatOp(:,:,1))
        !--------------------------------------------------


        !-------- deallocation ---------------------------------------------
        CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
        IF (allocated(para_BFGS%TS_vector)) THEN
          CALL dealloc_NParray(para_BFGS%TS_vector,'para_BFGS%TS_vector',name_sub)
        END IF
        IF (allocated(para_BFGS%hessian_inv_init)) THEN
          CALL dealloc_NParray(para_BFGS%hessian_inv_init,'para_BFGS%hessian_inv_init',name_sub)
        END IF
        !-------- end deallocation -----------------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Sub_BFGS

!!!!!!!!!!!!!!!!!!!JML!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------------
SUBROUTINE dfpmin_new(Qact,dnMatOp,mole,PrimOp,para_Tnum,para_BFGS,    &
                      gtol,tolx,itmax)
!---------------------------------------------------------------------------
!
 USE EVR_system_m
 use mod_Coord_KEO, only: CoordType, tnum, alloc_array, dealloc_array
 USE mod_PrimOp
 USE mod_basis
 USE mod_Op
 USE mod_Auto_Basis
!
!
!  Given a starting point p(1:n) that is a vector of length n, Broyden-Fletxer-
! Goldfrab-Shano varian of Davidon-Fletcher-Powell minimization is performed
! on a subroutine dfunc(...0), using this gradient as calculated by a routine dfunc(...1).
! The convergence requirement on zeroing the gradient by is input as gtol.
! Returned quantities are p(1:n) (the location of the minimum, iter (the number
! of iterations that were performed), and fret (the minimum value of the
! function). The routine lnsrch is called to perform approximate line
! minimizations, Parameters ITMAX
! is the maximum allowed number of iterations; STPMX is the scaled maximun step
! lengh allowed in the line search; TOLX is the convergence criterium in the
! x values
!
 implicit none

 real (kind=Rkind), pointer    :: p(:) => null()
 real (kind=Rkind), intent(in) ::  gtol, tolx
 integer,           intent(in) :: itmax
 integer :: n
 logical :: check
 real (kind=Rkind), parameter :: EPS=epsilon(p), STPMX=100
 real (kind=Rkind) :: xxxg, den, fp, temp, sum, stpmax, test, fret, fac, fae, fad, sumdg, sumxi
 integer :: its, iun,i,j
 real (kind=Rkind), allocatable :: g(:),dg(:),pnew(:),hdg(:),hessin(:,:),xi(:)

 real (kind=Rkind), target :: Qact(:)
 TYPE (param_dnMatOp)      :: dnMatOp(1)
 TYPE (CoordType)          :: mole
 TYPE (PrimOp_t)           :: PrimOp
 TYPE (Tnum)               :: para_Tnum
 TYPE (param_BFGS)         :: para_BFGS


 p => Qact(1:mole%nb_act)
 n =  mole%nb_act
 !write(out_unit,*) 'n',n
 !write(out_unit,*) 'gtol,tolx,itmax',gtol,tolx,itmax
 allocate (g(n),hdg(n),pnew(n),dg(n),xi(n),hessin(n,n))

 call dfunc(p,g,fp,dnMatOp,mole,PrimOp,para_Tnum,1) ! Calculate starting function value and gradient.
 write(out_unit,*)
 write(out_unit,*) 'Iteration=    0 '
 write(out_unit,*) ' Geometry:'
 write(out_unit,*) ' p = ', p
 write(out_unit,*) ' Energy = ',fp

 xxxg=sqrt(dot_product(g,g)/n)
!!!!!!!!!!!!!!
  test=ZERO            ! Test of convergence on zero gradient
  den=max(fp,ONE)
  do i=1,n
   temp=abs(g(i))*max(abs(p(i)),ONE)/den
   if (temp > test) test=temp
  end do
!!!!!!!!!!!!!!
 write(out_unit,*) ' RMS Gradient = ',xxxg
 write(out_unit,*) ' Test on gradient convergence = ',test
 flush(out_unit)

 IF (allocated(para_BFGS%hessian_inv_init)) THEN
   write(out_unit,*) ' The initial hessian is transferred'
   hessin(:,:) = para_BFGS%hessian_inv_init(:,:)
 ELSE
   hessin(:,:) = ZERO
   do i=1,n
     hessin(i,i)=ONE    !!!!!!!!!! Here Eckart has 1.0d0/a20(i)
   end do
 END IF

 xi=-g
 sum = dot_product(p,p)
 stpmax=STPMX*max(sqrt(sum),real(n,kind=Rkind))
 do its=1,ITMAX                 ! Main loop over the iterations.
  call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,dnMatOp,mole,PrimOp,para_Tnum)
  ! The new function evaluation occurs in lnsrch; save the value in fp for the
  ! next line search. It is usually safe to ignore the value of check.
  fp   = fret
  xi   = pnew-p    ! update the line direction
  IF (print_level > 1) write(out_unit,*) 'DeltaQ',xi(:)
  p    = pnew      ! and the current point
  test = ZERO
  do i=1,n     ! Test the convergence in Delta(x)
   temp=abs(xi(i))/max(abs(p(i)),ONE)
   if (temp > test) test=temp
  end do
  dg=g                  ! Save the old gradient
  call dfunc(p,g,fret,dnMatOp,mole,PrimOp,para_Tnum,1)  ! and get a new gradient

  test=ZERO            ! Test of convergence on zero gradient
  den=max(fret,ONE)
  do i=1,n
   temp=abs(g(i))*max(abs(p(i)),ONE)/den
   if (temp > test) test=temp
  end do
  xxxg=sqrt(dot_product(g,g)/n)

  write(out_unit,*)
  write(out_unit,'(a,i4)') ' Iteration= ',its
  IF (print_level > 0) write(out_unit,*) ' Geometry: '
  IF (print_level > 0) write(out_unit,*) ' p = ', p
  write(out_unit,*) ' Energy',fret
  write(out_unit,*) ' RMS Gradient',xxxg
  write(out_unit,*) ' Test on gradient convergence = ', test
  flush(out_unit)

   if (test < tolx) then  !!! Testing what happen if this part is removed
   write(out_unit,*) ' Geometry coordinates converged !! RMS step criteria'
   flush(out_unit)
   deallocate (g,hdg,pnew,dg,xi,hessin)
   return
  end if
!
  if ( test < gtol ) then
   write(out_unit,*) ' Gradient converged !!'
   write(out_unit,*) ' Optimized energy: ',fret
   write(out_unit,*) ' Optimized RMS gradient: ', xxxg, test
   write(out_unit,*) ' gtol: ', gtol
   write(out_unit,*) ' Optimized geometry: '
   write(out_unit,*) ' p = ', p
   deallocate (g,hdg,pnew,dg,xi,hessin)
   return
  end if
!
  dg=g-dg          ! Compute difference of gradient
  do i=1,n         ! and difference times current matrix.
    hdg(i) = dot_product(hessin(:,i),dg)
  end do
  fac   = dot_product(dg,xi) ! Calculate dot products for the denominators.
  fae   = dot_product(dg,hdg)
  sumdg = dot_product(dg,dg)
  sumxi = dot_product(xi,xi)

  if (fac**2 > EPS*sumdg*sumxi) then ! Skip update if fac is not
   fac=ONE/fac                        ! sufficiently positive.
   fad=ONE/fae
   do i=1,n         ! The vector that makes BFGS different from DFP.
    dg(i)=fac*xi(i)-fad*hdg(i)
   end do
   do i=1,n         ! The BFGS updating formula:
    do j=1,n
     hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
    end do
   end do
  end if
  do i=1,n          ! Now calculate the next direction to go,
    xi(i) = dot_product(hessin(:,i),g)
  end do
  xi=-xi
 end do             ! and go back for another iteration.
 deallocate (g,hdg,pnew,dg,xi,hessin)
 write(out_unit,*) 'Too many iterations in dfpmin'
 write(out_unit,*) ' energy: ',fret
 write(out_unit,*) ' Geometry: '
 write(out_unit,*) ' p = ', p
 stop
 return
 end subroutine dfpmin_new

!----------------------------------------------------------------------
  subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,dnMatOp,mole,PrimOp,para_Tnum)
!---------------------------------------------------------------------
!
 USE EVR_system_m
 USE mod_Coord_KEO, only: CoordType, tnum, alloc_array, dealloc_array
 USE mod_PrimOp
 USE mod_basis
 USE mod_Op
 USE mod_Auto_Basis
!
 implicit none
 LOGICAL :: check
 integer :: n,i
 real (kind=Rkind) :: g(n),p(n),x(n),xold(n)
 real (kind=Rkind), PARAMETER :: ALF=ONETENTH**4,TOLX=ONETENTH**9
 real (kind=Rkind) :: fold,f,stpmax,sum,slope,test,temp,alamin,alam,tmplam
 real (kind=Rkind) :: rhs1, rhs2, a, b, alam2, disc, f2
!
 TYPE (param_dnMatOp) :: dnMatOp(1)
 TYPE (CoordType)     :: mole
 TYPE (PrimOp_t)      :: PrimOp
 TYPE (Tnum)          :: para_Tnum
!
 check=.false.
 sum=sqrt(dot_product(p,p))
!
! write(out_unit,*) 'sum=', sum, 'stpmax=', stpmax
!
 if (sum > stpmax) then
  do i=1,n
   p(i)=p(i)*(stpmax/sum)
  end do
 endif
 slope = dot_product(g,p)
 test=ZERO
 do i=1,n
  temp=abs(p(i))/max(abs(xold(i)),ONE)
  if(temp.gt.test)test=temp
 end do
 alamin=TOLX/test
 alam=ONE
1 continue
 do  i=1,n
  x(i)=xold(i)+alam*p(i)
 end do
 call dfunc(x,g,f,dnMatOp,mole,PrimOp,para_Tnum,0)
 if(alam < alamin)then
  x=xold
  check=.true.
  return
 else if(f <= fold+ALF*alam*slope) then
  return
 else
  if(alam == ONE)then
    tmplam=-slope/(TWO*(f-fold-slope))
   else
    rhs1=f-fold-alam*slope
    rhs2=f2-fold-alam2*slope
    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
    if(a == ZERO)then
     tmplam=-slope/(TWO*b)
    else
     disc=b*b-THREE*a*slope
     if (disc < ZERO) then
      tmplam=HALF*alam
     else if (b < ZERO) then
      tmplam=(-b+sqrt(disc))/(THREE*a)
     else
      tmplam=-slope/(b+sqrt(disc))
     end if
    endif
    if(tmplam > HALF*alam) tmplam=HALF*alam
   endif
  endif
  alam2=alam
  f2=f
  alam=max(tmplam,ONETENTH*alam)
  goto 1
  end subroutine
!  (C) Copr. 1986-92 Numerical Recipes Software Bc21D#,#5,15!".
!
!---------------------------------------------------------------------------
  SUBROUTINE dfunc(xt,df,f,dnMatOp,mole,PrimOp,para_Tnum,nderiv_dnE)
!---------------------------------------------------------------------------
 USE EVR_system_m
 USE mod_dnSVM
 USE mod_Coord_KEO, only: CoordType, tnum, alloc_array, dealloc_array,get_Qact0,sub_QactTOdnx
 USE mod_PrimOp
 USE mod_basis
 USE mod_Op
 USE mod_Auto_Basis
!
! calculation of the gradient at xt
!
 IMPLICIT none
 integer, intent(in)  :: nderiv_dnE
 TYPE (param_dnMatOp) :: dnMatOp(1)
 TYPE (CoordType)     :: mole
 TYPE (PrimOp_t)      :: PrimOp
 TYPE (Tnum)          :: para_Tnum

 real(kind=Rkind), intent(in)    :: xt(mole%nb_act)
 real(kind=Rkind), intent(inout) :: df(mole%nb_act), f

 integer :: i
 real (kind=Rkind)    :: Qact(mole%nb_var)
 TYPE (Type_dnVec)    :: dnx


 !write(out_unit,*) 'dfunc subroutine',mole%nb_act,nderiv_dnE
 !write(out_unit,*) 'xt = ', xt

 Qact(:) = ZERO
 Qact(1:mole%nb_act)=xt

 IF (print_level > 1) THEN
   write(out_unit,*) '=== Current geometry (not recenter) ==========='
   CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv=0)

   CALL get_Qact0(Qact,mole%tab_Qtransfo(mole%itActive)%ActiveTransfo)
   CALL sub_QactTOdnx(Qact,dnx,mole,nderiv=0,Gcenter=.FALSE.,WriteCC=.TRUE.)

   CALL dealloc_dnSVM(dnx)
 END IF

! write(out_unit,*) 'Qact = ',Qact
! flush(out_unit)
!
! The subroutine below enables to calculate the energy, the gradient and/or the hessian
! nderiv_dnE = 0 : => energy only
! nderiv_dnE = 1 : => energy and gradient (here nderiv_dnE = 1)
! nderiv_dnE = 2 : => energy, gradient and hessian
! the derivatives are done with respect to the number of active coordinates

 CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)

 CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv_dnE)

! the results are in the matrix of derived type MatdnE(:,:)
!   Remark: we have a matrix because ElVibRot can deal with several diabatic electronic states
!
!  For one electronic state (which is the only possibility with on-the-fly calculation)
!    the energy   is in : MatdnE(1,1)%d0
!    the gradient is in : MatdnE(1,1)%d1(:)     (if nderiv_dnE >= 1)
!    the hessian  is in : MatdnE(1,1)%d2(:,:)   (if nderiv_dnE >= 2)

!  Remark, the allocation of MatdnE have to be done with "nderiv_alloc" larger than "nderiv_dnE"
!     With nderiv_alloc=2 and nderiv_dnE=1, you can calculate the energy and the gradient
!     and update the hessian in MatdnE(1,1)%d2(:,:), if you want.

!  write(out_unit,*) ' MatdnE(1,1)%d0 = ', MatdnE(1,1)%d0
!  write(out_unit,*) ' MatdnE(1,1)%d1 = ', MatdnE(1,1)%d1
!
  f = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp,1)
  if (nderiv_dnE >= 1) CALL Get_Grad_FROM_Tab_OF_dnMatOp(df,dnMatOp,1)

  write(out_unit,*) ' Energy = ',f
  IF (nderiv_dnE >= 1) write(out_unit,*) ' Active modes gradient norm = ',sqrt(dot_product(df,df))
  IF (print_level > 1) write(out_unit,*) ' Active modes gradient = ', df
  !write(out_unit,*) 'end dfunc subroutine'
  !flush(out_unit)

 END SUBROUTINE dfunc

  SUBROUTINE Sub_Newton(BasisnD,xOpt_min,SQ,nb_Opt,para_Tnum,mole,PrimOp,Qopt,para_BFGS)

        USE EVR_system_m
        USE mod_dnSVM
        use mod_Coord_KEO, only: CoordType, tnum, alloc_array, dealloc_array
        USE mod_PrimOp
        USE mod_basis
        USE mod_Op
        USE mod_Auto_Basis
        IMPLICIT NONE
        
        !----- for the CoordType and Tnum --------------------------------------
        TYPE (CoordType) :: mole
        TYPE (Tnum)    :: para_Tnum
        logical        :: Cart_Transfo_save
        real (kind=Rkind) :: Qact(mole%nb_var)
        
        
        !----- for the basis set ----------------------------------------------
        TYPE (basis)          :: BasisnD
        
        !----- variables pour la namelist minimum ----------------------------
        TYPE (PrimOp_t)  :: PrimOp
        integer          :: nb_scalar_Op
        logical          :: calc_scalar_Op
        
        !----- variables for the construction of H ---------------------------
        TYPE (ReadOp_t)         :: para_ReadOp
        logical                     :: Save_FileGrid,Save_MemGrid
        
        
        !----- local variables -----------------------------------------------
        !----- variables for the construction of H ---------------------------
        TYPE (param_AllOp)          :: para_AllOp_loc
        
        !----- for the basis set ----------------------------------------------
        TYPE (param_AllBasis)       :: para_AllBasis_loc
        TYPE (basis)                :: basis_temp
        
        !----- for the optimization -------------------------------------------
        TYPE (param_BFGS) :: para_BFGS
        integer, intent(in) :: nb_Opt
        real (kind=Rkind), intent(inout) :: Qopt(:)
        logical :: TS_Vector

        real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt)
        real (kind=Rkind), allocatable :: grad(:),hess(:,:),hess0(:,:),diag(:),Vec(:,:),Vec_ext(:,:),tVec(:,:),mDQit(:)
        real (kind=Rkind) :: Ene,Ene0,sc_mDQit
        real (kind=Rkind) :: max_grad,RMS_grad,max_step,RMS_step,norm_disp,Norm2
        logical :: conv
        integer :: it,iq,iqq,its
        !---------- working variables ----------------------------------------
        TYPE (param_dnMatOp) :: dnMatOp(1)
        integer              :: nderiv_alloc,nbcol,err
        
        !---------------------------------------------------------------------
        logical,parameter :: debug= .FALSE.
        !logical,parameter :: debug= .TRUE.
        character (len=*), parameter :: name_sub='Sub_Newton'
        !---------------------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*)
          write(out_unit,*) 'BEGINNING ',name_sub
          write(out_unit,*) 'xopt_min',xopt_min(:)
          write(out_unit,*)
          write(out_unit,*) 'mole%nb_act',mole%nb_act
          write(out_unit,*) 'RMS_grad',para_BFGS%RMS_grad
          write(out_unit,*) 'RMS_step',para_BFGS%RMS_step
          write(out_unit,*) 'max_iteration',para_BFGS%max_iteration
          flush(out_unit)
        END IF
        !---------------------------------------------------------------------
        !Ene0 = huge(ONE)
        Qopt(1:nb_Opt) = xopt_min(:)
        !write(out_unit,*) 'Qopt',Qopt

        CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_Opt,PrimOp%nb_elec,nderiv=2)
        CALL alloc_NParray(hess,[nb_Opt,nb_Opt],'hess',name_sub)
        CALL alloc_NParray(grad,[nb_Opt],       'grad',name_sub)
        CALL alloc_NParray(vec, [nb_Opt,nb_Opt],'vec',name_sub)
        CALL alloc_NParray(diag,[nb_Opt],       'diag',name_sub)
        write(out_unit,*) '=================================================='
        DO it=0,para_BFGS%max_iteration

          CALL get_Qact0(Qact,mole%tab_Qtransfo(mole%itActive)%ActiveTransfo)
          Qact(1:nb_Opt) = Qopt(:)
          IF (debug) write(out_unit,*) 'Qact',Qact
          IF (mod(it,para_BFGS%calc_hessian_every_n) == 0) THEN
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=2)
            CALL Get_Hess_FROM_Tab_OF_dnMatOp(hess,dnMatOp) ! for the ground state
            IF (debug) CALL Write_Mat_MPI(hess, out_unit, 5, info='hess')
            hess0 = hess
          ELSE
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=1)
            hess = hess0
          END IF
        
          !IF (debug) CALL Write_Tab_OF_dnMatOp(dnMatOp)
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(grad,dnMatOp)  
          IF (debug) CALL Write_Vec_MPI(grad, out_unit, 5, info='grad')
          Ene = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
          IF (it ==0) Ene0 = Ene

          IF (allocated(para_BFGS%TS_Vector)) THEN
            TS_Vector = ( para_BFGS%nb_neg == 1 .AND. dot_product(para_BFGS%TS_Vector,para_BFGS%TS_Vector) > HALF )
          ELSE
            TS_Vector = .FALSE.
          END IF

          IF (TS_Vector) THEN
            ! make the new vectors
            CALL alloc_NParray(vec_ext, [nb_Opt,nb_Opt+1],'vec_ext',name_sub)
            CALL diagonalization(hess,diag,vec_ext(:,2:nb_Opt+1),nb_Opt)
            vec_ext(:,1) = para_BFGS%TS_Vector
            vec_ext(:,:) = Ortho_GramSchmidt(vec_ext)
            IF (debug) CALL Write_Mat_MPI(vec_ext, out_unit, 5, info='vec_ortho')
            iqq = 0
            DO iq=1,nb_Opt+1
              IF (dot_product(vec_ext(:,iq),vec_ext(:,iq)) > HALF) THEN
                iqq = iqq+1
                vec(:,iqq) = vec_ext(:,iq)
              END IF
            END DO
            CALL dealloc_NParray(vec_ext,'vec_ext',name_sub)

            tvec = transpose(vec)

            !new hessian on the vec
            hess = matmul(tvec,matmul(hess,Vec))
            !remove some coupling and change the sign of hess(its,its)
            its = 1
            Norm2 = -abs(hess(its,its))
            hess(:,its)   = ZERO
            hess(its,:)   = ZERO
            hess(its,its) = Norm2
            !transform the new hessian in the original basis
            hess = matmul(vec,matmul(hess,tvec))
          ELSE
            CALL diagonalization(hess,diag,Vec,nb_Opt)
            IF (debug) write(out_unit,*) 'diag',diag ; flush(out_unit)
            write(out_unit,*) ' Number hessian negative eigenvalue(s): ',count(diag < ZERO)

            IF (para_BFGS%nb_neg == 0) THEN
              IF (count(diag < ZERO) /= 0) THEN 
                write(out_unit,*) 'WARNING: some hessian eigenvalues are negative!'
                DO iq=1,nb_Opt
                  IF (diag(iq) < ZERO) write(out_unit,*) 'negative eigenvalue',iq,diag(iq)
                END DO
              END IF
              tvec = transpose(vec)
              DO iq=1,nb_Opt
                Vec(:,iq) = Vec(:,iq) * abs(diag(iq))
              END DO
              hess = matmul(Vec,tVec)
            ELSE IF (count(diag < ZERO) /= para_BFGS%nb_neg) THEN
              write(out_unit,*) 'ERROR in ',name_sub
              write(out_unit,*) '    Wrong number of negative hessian eigenvalues!'
              write(out_unit,*) '    Expected: ',para_BFGS%nb_neg
              write(out_unit,*) '    it has: ',count(diag < ZERO)
              STOP 'ERROR in Sub_Newton: Wrong number of negative hessian eigenvalues'
            END IF
          END IF
          !write(out_unit,*) 'hess',hess
          !write(out_unit,*) 'hess?',matmul(Vec,tVec)
        
          mDQit = LinearSys_Solve(hess,grad)
          IF (debug) CALL Write_Vec_MPI(mDQit, out_unit, 5, info='mDQit')

        
          max_grad = maxval(abs(grad))
          RMS_grad = sqrt(dot_product(grad,grad)/nb_Opt)
          max_step = maxval(abs(mDQit))
          RMS_step = sqrt(dot_product(mDQit,mDQit)/nb_Opt)
        
        
          write(out_unit,*) '--------------------------------------------------'
          write(out_unit,*) 'it,E',it,Ene

          conv = (max_grad <= para_BFGS%max_grad)
          write(out_unit,*) 'it,max_grad,treshold',it,max_grad,para_BFGS%max_grad,conv
          conv = (RMS_grad <= para_BFGS%RMS_grad)
          write(out_unit,*) 'it,RMS_grad,treshold',it,RMS_grad,para_BFGS%RMS_grad,conv
          conv = (max_step <= para_BFGS%max_step)
          write(out_unit,*) 'it,max_step,treshold',it,max_step,para_BFGS%max_step,conv
          conv = (RMS_step <= para_BFGS%RMS_step)
          write(out_unit,*) 'it,RMS_step,treshold',it,RMS_step,para_BFGS%RMS_step,conv
        
          norm_disp = sqrt(dot_product(mDQit,mDQit))
          IF (norm_disp > para_BFGS%Largest_step) THEN
            write(out_unit,*) ' The displacements are too large.'
            IF (print_level > 1) write(out_unit,*) ' The displacements:',-mDQit
            write(out_unit,*) ' Its norm:',norm_disp
            write(out_unit,*) '  => They are scaled by ',para_BFGS%Largest_step/norm_disp
            mDQit = mDQit * para_BFGS%Largest_step/norm_disp
          END IF
        
          DO iq=1,nb_Opt
            write(out_unit,*) 'iq,Q(iq),grad(iq),DelatQ(iq)',iq,Qopt(iq),grad(iq),-mDQit(iq)
          END DO
        
          conv = (max_grad <= para_BFGS%max_grad) .AND.                               &
                 (RMS_grad <= para_BFGS%RMS_grad) .AND.                               &
                 (max_step <= para_BFGS%max_step) .AND.                               &
                 (RMS_step <= para_BFGS%RMS_step)
        
        
          IF (conv) THEN 
            Qopt(1:nb_Opt) = Qopt(1:nb_Opt)-mDQit
            EXIT
          END IF

          IF (para_BFGS%nb_neg == 0) THEN
            sc_mDQit = ONE
            DO
              CALL get_Qact0(Qact,mole%tab_Qtransfo(mole%itActive)%ActiveTransfo)
              Qact(1:nb_Opt) = Qopt(:)-sc_mDQit*mDQit
              IF (debug) write(out_unit,*) 'Qact',Qact
              CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)
              Ene = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
              write(out_unit,*) 'sc_mDQit,Ene,Ene0',sc_mDQit,Ene,Ene0
              IF (Ene < Ene0) THEN
                Qopt(:) = Qact(1:nb_Opt)
                Ene0 = Ene
                EXIT
              END IF
              sc_mDQit = sc_mDQit * HALF
            END DO
          ELSE
            Qopt(1:nb_Opt) = Qopt(1:nb_Opt)-mDQit
          END IF



        END DO
        IF (para_BFGS%max_iteration > 0) THEN
          write(out_unit,*) 'Geometry optimization is converged?',conv
          write(out_unit,*) 'Optimized geometry:'
          DO iq=1,nb_Opt
            write(out_unit,*) 'iq,Qopt(iq),',iq,Qopt(iq)
          END DO
        ELSE
          write(out_unit,*) 'No optimization (Max_it=0)'
        END IF
        write(out_unit,*) '=================================================='
        flush(out_unit)
        
        xopt_min(:) = Qopt(1:nb_Opt)
        
        !--------------------------------------------------
        ! this subroutine print the  matrix of derived type.
        IF (debug) CALL Write_MatOFdnS(dnMatOp(1)%tab_dnMatOp(:,:,1))
        !--------------------------------------------------
        
        
        !-------- deallocation ---------------------------------------------
        CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
        !-------- end deallocation -----------------------------------------
        
        
        !---------------------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'END ',name_sub
          flush(out_unit)
        END IF
        !---------------------------------------------------------------------

  END SUBROUTINE Sub_Newton

END MODULE mod_BFGS
