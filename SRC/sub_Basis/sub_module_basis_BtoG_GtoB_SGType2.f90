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
MODULE mod_basis_BtoG_GtoB_SGType2
USE EVR_system_m
USE mod_basis_set_alloc
IMPLICIT NONE


TYPE TypeRDP
  integer :: n2 = 0
  integer :: n3 = 0
  real(kind=Rkind), allocatable :: RDP(:,:)
END TYPE TypeRDP

CONTAINS

SUBROUTINE alloc_TabRDP_at_iG(TabRDP,iG)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)
integer                    :: iG

IF (allocated(TabRDP)) THEN
  !write(out_unit,*) 'BEGINNING alloc_TabRDP_at_iG'

  IF (.NOT. allocated(TabRDP(iG)%RDP)) THEN
    CALL alloc_NParray(TabRDP(iG)%RDP,                                  &
                                    [TabRDP(iG)%n2, TabRDP(iG)%n3], &
                      'TabRDP(iG)%RDP','alloc_TabRDP_at_iG')
    TabRDP(iG)%RDP = ZERO
  END IF

!write(out_unit,*) 'END alloc_TabRDP_at_iG'
END IF



END SUBROUTINE alloc_TabRDP_at_iG

SUBROUTINE dealloc_TabRDP_at_iG(TabRDP,iG)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)
integer               :: iG

IF (allocated(TabRDP)) THEN
  !write(out_unit,*) 'BEGINNING dealloc_TabRDP_at_iG'

    !write(out_unit,*) 'iG, Shape RDP',iG,shape(TabRDP(iG)%RDP)

    IF (allocated(TabRDP(iG)%RDP))                                      &
      CALL dealloc_NParray(TabRDP(iG)%RDP,'TabRDP(iG)%RDP','dealloc_TabRDP_at_iG')

!write(out_unit,*) 'END dealloc_TabRDP_at_iG'
END IF



END SUBROUTINE dealloc_TabRDP_at_iG

SUBROUTINE dealloc_TabRDP(TabRDP)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)


integer               :: iG


IF (allocated(TabRDP)) THEN
  !write(out_unit,*) 'BEGINNING dealloc_TabRDP'
  !write(out_unit,*) 'dealloc TabRDP',size(TabRDP)

  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    !write(out_unit,*) 'iG, Shape RDP',iG,shape(TabRDP(iG)%RDP)

    TabRDP(iG)%n2 = 0
    TabRDP(iG)%n3 = 0

    IF (allocated(TabRDP(iG)%RDP))                                      &
      CALL dealloc_NParray(TabRDP(iG)%RDP,'TabRDP(iG)%RDP','dealloc_TabRDP')

  END DO

  deallocate(TabRDP)
  !write(out_unit,*) 'alloc TabRDP ?',allocated(TabRDP)
  !write(out_unit,*) 'END dealloc_TabRDP'
END IF



END SUBROUTINE dealloc_TabRDP


SUBROUTINE Write_TabRDP(TabRDP)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)


integer               :: iG,i2,i3


IF (allocated(TabRDP)) THEN
  write(out_unit,*) '======== TabRDP ============================',shape(TabRDP)
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    IF (allocated(TabRDP(iG)%RDP)) THEN
      write(out_unit,*) '    iG',iG
      DO i3=1,TabRDP(iG)%n3
      DO i2=1,TabRDP(iG)%n2
        write(out_unit,*) 'i2,i3',i2,i3,TabRDP(iG)%RDP(i2,i3)
      END DO
      END DO
    ELSE
      write(out_unit,*) '    iG',iG,' RDP is not allocated'
    END IF

  END DO

END IF


END SUBROUTINE Write_TabRDP

SUBROUTINE SumSq_TabRDP(TabRDP)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)

integer               :: iG,i1,i2,i3
real (kind=Rkind)     :: SS,S


SS = ZERO
IF (allocated(TabRDP)) THEN
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    IF (allocated(TabRDP(iG)%RDP)) S  = sum(TabRDP(iG)%RDP**2)
    SS = SS + S
  END DO
  write(out_unit,*) 'SumSq TabRDP',SS

END IF


END SUBROUTINE SumSq_TabRDP

SUBROUTINE Size_TabRDP(TabRDP,nb_BG)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)
integer  :: nb_BG

integer               :: iG
real(kind=Rkind)  :: Rnb_BG


nb_BG  = 0
Rnb_BG = ZERO
IF (allocated(TabRDP)) THEN
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    IF (allocated(TabRDP(iG)%RDP)) THEN
      Rnb_BG = Rnb_BG + real(size(TabRDP(iG)%RDP),kind=Rkind)
    ELSE
      Rnb_BG = Rnb_BG + real(TabRDP(iG)%n2*TabRDP(iG)%n3,kind=Rkind)
    END IF
  END DO
END IF

IF (Rnb_BG < 1024._Rkind) THEN
  write(out_unit,*) 'size WPG',int(Rnb_BG),' Words'
ELSE IF (Rnb_BG < 1024._Rkind**2) THEN
  write(out_unit,*) 'size WPG',int(Rnb_BG/(1024._Rkind**1))*Rkind,' kW (Words)'
ELSE
  write(out_unit,*) 'size WPG',int(Rnb_BG/(1024._Rkind**2))*Rkind,' MW (Words)'
END IF
flush(out_unit)

END SUBROUTINE Size_TabRDP

SUBROUTINE Norm_OFF_Diff_TabRDP(TabRDP1,TabRDP2)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP)      :: TabRDP1(:),TabRDP2(:)


integer            :: iG
real(kind=Rkind)   :: Norm


IF (TabRDP1(1)%n2 /= TabRDP2(1)%n2 .OR.                                 &
    TabRDP1(1)%n3 /= TabRDP2(1)%n3 .OR.                                 &
    size(TabRDP1) /= size(TabRDP2) ) THEN
  write(out_unit,*) ' ERROR in Norm_OFF_Diff_TabRDP'
  write(out_unit,*) ' incompatible TabRDP1 and TabRDP2'
  STOP
END IF

Norm = ZERO
DO iG=1,ubound(TabRDP1,dim=1)
  Norm = Norm + sum( (TabRDP1(iG)%RDP(:,:) - TabRDP2(iG)%RDP(:,:))**2 )
END DO
write(out_unit,*) 'Norm_OFF_Diff_TabRDP',Norm


END SUBROUTINE Norm_OFF_Diff_TabRDP
SUBROUTINE TabRDP2_TO_TabRDP1(TabRDP1,TabRDP2)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP1(:),TabRDP2(:)


integer            :: iG
real(kind=Rkind)   :: Norm

IF (allocated(TabRDP1)) CALL dealloc_TabRDP(TabRDP1)

IF (allocated(TabRDP2)) THEN
  allocate(TabRDP1( lbound(TabRDP2,dim=1) : ubound(TabRDP2,dim=1) ))
  !write(out_unit,*) 'alloc TabRDP',size(TabRDP1)


  DO iG=lbound(TabRDP1,dim=1),ubound(TabRDP1,dim=1)

    IF (allocated(TabRDP2(iG)%RDP)) THEN
      TabRDP1(iG)%n2 = TabRDP2(iG)%n2
      TabRDP1(iG)%n3 = TabRDP2(iG)%n3

      CALL alloc_NParray(TabRDP1(iG)%RDP,                               &
                                   [TabRDP1(iG)%n2,TabRDP1(iG)%n3], &
                        'TabRDP1(iG)%RDP','TabRDP2_TO_TabRDP1')
      TabRDP1(iG)%RDP(:,:) = TabRDP2(iG)%RDP(:,:)
    END IF

  END DO
END IF

END SUBROUTINE TabRDP2_TO_TabRDP1

SUBROUTINE Set_BgG_FOR_id(BgG,ind_Grid,ind_Basis,tab_ba,D,LG,id,WithAlloc_RDP)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)

integer                       :: D,LG,id
logical                       :: WithAlloc_RDP
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:LG,D)

integer          :: iG,i,li,nb,nbb,nq,nqq,iq,ll,nb_BG
logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

CALL dealloc_TabRDP(BgG)

IF (debug) THEN
  write(out_unit,*) 'BEGINNING in Set_BgG_FOR_id',id
  write(out_unit,*) 'tab_q, for B',ind_Basis(id+1)%tab_q
  write(out_unit,*) 'tab_q, for G',ind_Grid(id)%tab_q
  write(out_unit,*) 'tab_q, for g',0
END IF


allocate(BgG(ind_Grid(id)%MaxnD))
!write(out_unit,*) 'alloc TabRDP',size(BgG)

DO iG=1,ind_Grid(id)%MaxnD

  nqq = 1
  !write(out_unit,*) 'iG,l(:)',ind_Grid(id)%tab_ind(:,iG)
  DO i=1,ind_Grid(id)%ndim
    li = ind_Grid(id)%tab_ind(i,iG)
    iq = ind_Grid(id)%tab_q(i)
    IF (iq < 1 .OR. iq > D) CYCLE
    nqq = nqq * get_nq_FROM_basis(tab_ba(li,iq))
  END DO

  nbb = ind_Basis(id+1)%MaxnD


  BgG(iG)%n2 = nbb
  BgG(iG)%n3 = nqq

  IF (WithAlloc_RDP) THEN
    CALL alloc_NParray(BgG(iG)%RDP,[BgG(iG)%n2,BgG(iG)%n3],          &
                      'BgG(iG)%RDP','Set_BgG_FOR_id')
    BgG(iG)%RDP(:,:) = ZERO
  END IF

  !IF (debug) write(out_unit,*) 'iG, Shape RDP',iG,':',shape(BgG(iG)%RDP)

END DO

IF (debug) THEN
  !CALL Write_TabRDP(BgG)
  CALL Size_TabRDP(BgG,nb_BG)
  write(out_unit,*) 'END in Set_BgG_FOR_id',id
  flush(out_unit)
END IF
END SUBROUTINE Set_BgG_FOR_id

SUBROUTINE Transfer_WP_TO_BgG(WPG,BgG)
USE EVR_system_m
IMPLICIT NONE

real(kind=Rkind)   :: WPG(:)
TYPE(TypeRDP)      :: BgG(:)


integer               :: i_WP,iG,nqq


IF (BgG(1)%n2 /= 1) THEN
  write(out_unit,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unit,*) ' BgG is not completely on the G grid'
  STOP
END IF
IF (size(WPG) /= sum(BgG(:)%n3)) THEN
  write(out_unit,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unit,*) ' The sizes of WPG and BgG are not compatible'
  write(out_unit,*) 'size(WPG)',size(WPG)
  write(out_unit,*) 'sum(BgG(:)%n3)',sum(BgG(:)%n3)
  STOP
END IF


i_WP = 0
DO iG=1,ubound(BgG,dim=1)
  nqq = BgG(iG)%n3
  BgG(iG)%RDP(1,:) = WPG(i_WP+1:i_WP+nqq)
  i_WP = i_WP+nqq
END DO
!write(out_unit,*) 'i_WP',i_WP


END SUBROUTINE Transfer_WP_TO_BgG

SUBROUTINE Transfer_BgG_TO_WP(BgG,WPG)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP)      :: BgG(:)
real(kind=Rkind)   :: WPG(:)


integer               :: i_WP,iG,nqq


IF (BgG(1)%n2 /= 1) THEN
  write(out_unit,*) ' ERROR in Transfer_BgG_TO_WP'
  write(out_unit,*) ' BgG is not completely on the G grid'
  STOP
END IF
IF (size(WPG) /= sum(BgG(:)%n3)) THEN
  write(out_unit,*) ' ERROR in Transfer_BgG_TO_WP'
  write(out_unit,*) ' The sizes of WPG and BgG are not compatible'
  write(out_unit,*) 'size(WPG)',size(WPG)
  write(out_unit,*) 'sum(BgG(:)%n3)',sum(BgG(:)%n3)
  STOP
END IF


i_WP = 0
DO iG=1,ubound(BgG,dim=1)
  nqq = BgG(iG)%n3
  WPG(i_WP+1:i_WP+nqq) = BgG(iG)%RDP(1,:)
  i_WP = i_WP+nqq
END DO
!write(out_unit,*) 'i_WP',i_WP


END SUBROUTINE Transfer_BgG_TO_WP

SUBROUTINE Norm_OFF_Diff_WP_BgG(WPG,BgG)
USE EVR_system_m
IMPLICIT NONE

TYPE(TypeRDP)      :: BgG(:)
real(kind=Rkind)   :: WPG(:)


integer            :: i_WP,iG,nqq
real(kind=Rkind)   :: Norm


IF (BgG(1)%n2 /= 1) THEN
  write(out_unit,*) ' ERROR in Norm_OFF_Diff_WP_BgG'
  write(out_unit,*) ' BgG is not completely on the G grid'
  STOP
END IF
IF (size(WPG) /= sum(BgG(:)%n3)) THEN
  write(out_unit,*) ' ERROR in Norm_OFF_Diff_WP_BgG'
  write(out_unit,*) ' The sizes of WPG and BgG are not compatible'
  write(out_unit,*) 'size(WPG)',size(WPG)
  write(out_unit,*) 'sum(BgG(:)%n3)',sum(BgG(:)%n3)
  STOP
END IF


i_WP = 0
Norm = ZERO
DO iG=1,ubound(BgG,dim=1)
  nqq = BgG(iG)%n3
  Norm = Norm + sum( (BgG(iG)%RDP(:,1) - WPG(i_WP+1:i_WP+nqq))**2 )
  i_WP = i_WP+nqq
END DO
write(out_unit,*) 'Norm_OFF_Diff_WP_BgG',Norm


END SUBROUTINE Norm_OFF_Diff_WP_BgG

SUBROUTINE BgG_TO_BbG(BgG,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)
real(kind=Rkind)   :: WSG(:)



integer            :: D,LG,LB,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:LG,D)

TYPE(TypeRDP), allocatable    :: BgG_new(:)


integer               :: i_WP,iG,iq,iqq,nq,nqq,ibb,nbb,iqqi,iqqf,ib,nb,li,i,iGm1,lbb,lqq,ibbNew
integer               :: tabl(D)
integer               :: tabnq(D)


integer :: max_nq,max_nb,Lmax
real(kind=Rkind), allocatable  :: b(:)
integer, allocatable  :: tab_ibbnew_AT_ibb(:)


logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

IF (debug) THEN
  write(out_unit,*)
  write(out_unit,*) 'BEGINNING BgG_TO_BbG',id
  CALL Write_TabRDP(BgG)
  write(out_unit,*) 'WSG(:)',WSG(:)
END IF

CALL Set_BgG_FOR_id(BgG_new,ind_Grid,ind_Basis,tab_ba,D,LG,id,.FALSE.)

IF (debug) THEN
  write(out_unit,*) 'size(BgG),size(BgG_new)',size(BgG),size(BgG_new)
  flush(out_unit)
END IF

! id index of "b" of BbG
IF (size(BgG) /= ind_Grid(id-1)%MaxnD) THEN
  write(out_unit,*) ' ERROR in BgG_TO_BbG'
  write(out_unit,*) ' size(BgG) /= ind_Grid(id-1)%MaxnD',size(BgG),ind_Grid(id-1)%MaxnD
  STOP
END IF

IF (size(BgG_new) /= ind_Grid(id)%MaxnD) THEN
  write(out_unit,*) ' ERROR in BgG_TO_BbG'
  write(out_unit,*) ' size(BgG_new) /= ind_Grid(id)%MaxnD',size(BgG_new),ind_Grid(id)%MaxnD
  STOP
END IF

iGm1 = maxval(ind_Grid(id-1)%indD_OF_Dm1)
IF (iGm1 /= size(BgG_new)) THEN
  write(out_unit,*) ' ERROR in BgG_TO_BbG'
  write(out_unit,*) ' size(BgG_new) /= iGm1',size(BgG_new),iGm1
  STOP
END IF


max_nq = 0
Lmax   = ubound(tab_ba,dim=1)
DO iq=1,ubound(tab_ba,dim=2)
  nq = get_nq_FROM_basis(tab_ba(Lmax,iq))
  IF (nq > max_nq) max_nq = nq
END DO
max_nb = maxval(tab_ba(:,:)%nb)

!write(out_unit,*) 'nb threads (BasisTOGrid_maxth):',BasisTOGrid_maxth

DO iG=1,size(BgG)

  DO i=1,ind_Grid(id-1)%ndim
    tabl(i) = ind_Grid(id-1)%tab_ind(i,iG)
    iq = ind_Grid(id-1)%tab_q(i)
    IF (iq < 1 .OR. iq > D) CYCLE
    tabnq(i) = get_nq_FROM_basis(tab_ba(tabl(i),iq))
  END DO
  !write(out_unit,*) 'iG,tabl ',iG,tabl(1:ind_Grid(id-1)%ndim)
  !write(out_unit,*) 'iG,tabnq',iG,tabnq(1:ind_Grid(id-1)%ndim)

  nbb  = BgG(iG)%n2
  nqq  = product(tabnq(2:ind_Grid(id-1)%ndim))
  iGm1 = ind_Grid(id-1)%indD_OF_Dm1(iG)

  CALL alloc_TabRDP_at_iG(BgG_new,iGm1)

  !write(out_unit,*) '---------------------' ; flush(out_unit)
  !write(out_unit,*) 'iG,nbb,nqq',iG,nbb,nqq ; flush(out_unit)
  allocate(tab_ibbnew_AT_ibb(nbb))
  ibbNew = 0
  DO ibb=1,nbb
    tab_ibbnew_AT_ibb(ibb) = ibbNew
    lbb = ind_Basis(id)%i_TO_l(ibb)
    ibbNew = ibbNew + tab_ba(LB-lbb,id)%nb
  END DO

  nq = tabnq(1)


  !$OMP   PARALLEL &
  !$OMP   DEFAULT(NONE) &
  !$OMP   SHARED(tab_ibbnew_AT_ibb,nbb,nqq,LB,tabl,tabnq,nq) &
  !$OMP   SHARED(BgG,BgG_new,iG,iGm1,max_nq,max_nb) &
  !$OMP   SHARED(id,nb_mult_GTOB,tab_ba,ind_Basis,ind_Grid,WSG) &
  !$OMP   PRIVATE(ibb,ibbNew,lbb,li,nb,iqq,iqqi,iqqf,b,iq) &
  !$OMP   NUM_THREADS(BasisTOGrid_maxth)
   CALL alloc_NParray(b,  [max_nb],'b',  'BgG_TO_BbG')
  !$OMP   DO SCHEDULE(DYNAMIC)
  DO iqq=1,nqq
    iqqi = (iqq-1)*nq
    iqqf = iqqi + nq



    DO ibb=1,nbb
      ibbNew = tab_ibbnew_AT_ibb(ibb)
      lbb    = ind_Basis(id)%i_TO_l(ibb)

      li     = min(tabl(1),LB-lbb)
      iq     = ind_Grid(id-1)%tab_q(1)
      nb     = min( tab_ba(li,iq)%nb , tab_ba(LB-lbb,id)%nb )

      li = tabl(1)
      CALL RG_TO_RB_basis(BgG(iG)%RDP(ibb,iqqi+1:iqqf),b(1:nb),tab_ba(li,iq))

      IF (id == 1) THEN
        b(1:nb) = b(1:nb) *  WSG(iG)
      END IF


      DO ib=1,nb
        !$OMP ATOMIC
        BgG_new(iGm1)%RDP(ibbNew+ib,iqq) = BgG_new(iGm1)%RDP(ibbNew+ib,iqq) + b(ib)
      END DO

    END DO

  END DO
  !$OMP   END DO
   CALL dealloc_NParray(b,  'b',  'BgG_TO_BbG')
  !$OMP   END PARALLEL

  CALL dealloc_TabRDP_at_iG(BgG,iG)
  deallocate(tab_ibbnew_AT_ibb)


END DO


CALL TabRDP2_TO_TabRDP1(BgG,BgG_new)
CALL dealloc_TabRDP(BgG_new)


IF (debug) THEN
  CALL Write_TabRDP(BgG)
  CALL SumSq_TabRDP(BgG)
  write(out_unit,*) 'END BgG_TO_BbG',id
  flush(out_unit)
END IF

END SUBROUTINE BgG_TO_BbG


SUBROUTINE BbG_TO_BgG(BgG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)

integer            :: D,LG,LB,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:LG,D)

TYPE(TypeRDP), allocatable    :: BgG_new(:)

integer               :: i_WP,iq,iqq,nq,nqq,ibb,nbb,iqqi,iqqf,ib,nb,li,liq,i,lbb,lqq,ibbNew
integer               :: iG,iGm1,iG1,iG2,max_nb

integer, allocatable  :: tab_ibbnew_AT_ibb(:)

logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

IF (debug) THEN
  write(out_unit,*)
  write(out_unit,*) 'BEGINNING BbG_TO_BgG',id
  CALL Write_TabRDP(BgG)
END IF


CALL Set_BgG_FOR_id(BgG_new,ind_Grid,ind_Basis,tab_ba,D,LG,id-1,.FALSE.)

IF (debug) THEN
  write(out_unit,*) 'size(BgG_new),size(BgG)',size(BgG_new),size(BgG)
  flush(out_unit)
END IF

! id index of "b" of BgG_new
IF (size(BgG_new) /= ind_Grid(id-1)%MaxnD) THEN
  write(out_unit,*) ' ERROR in BbG_TO_BgG'
  write(out_unit,*) ' size(BgG_new) /= ind_Grid(id-1)%MaxnD',size(BgG_new),ind_Grid(id-1)%MaxnD
  STOP
END IF

iGm1 = maxval(ind_Grid(id-1)%indD_OF_Dm1)
IF (iGm1 /= size(BgG)) THEN
  write(out_unit,*) ' ERROR in BbG_TO_BgG'
  write(out_unit,*) ' size(BgG_new) /= iGm1',size(BgG),iGm1
  STOP
END IF


iG2 = 0
DO iGm1=1,size(BgG)
iG1 = iG2
iG2 = iG1 + count(ind_Grid(id-1)%indD_OF_Dm1 == iGm1)
DO iG=iG1+1,iG2

  CALL alloc_TabRDP_at_iG(BgG_new,iG)


  nbb  = BgG_new(iG)%n2
  nqq  = BgG(iGm1)%n3



  !write(out_unit,*) '---------------------' ; flush(out_unit)
  !write(out_unit,*) 'iG,iGm1,nbb,nqq',iG,iGm1,nbb,nqq ; flush(out_unit)
  allocate(tab_ibbnew_AT_ibb(nbb))
  ibbNew = 0
  DO ibb=1,nbb
    tab_ibbnew_AT_ibb(ibb) = ibbNew
    lbb = ind_Basis(id)%i_TO_l(ibb)
    ibbNew = ibbNew + tab_ba(LB-lbb,id)%nb
    !write(out_unit,*) 'ibb,lbb,nb,ibbNew',ibb,lbb,tab_ba(LB-lbb,id)%nb,ibbNew
  END DO

  liq = ind_Grid(id-1)%tab_ind(1,iG)
  nq  = get_nq_FROM_basis(tab_ba(liq,id))
  max_nb = maxval(tab_ba(:,:)%nb)

  !$OMP   PARALLEL &
  !$OMP   DEFAULT(NONE) &
  !$OMP   SHARED(tab_ibbnew_AT_ibb,nbb,nqq,nq,liq,LB,max_nb) &
  !$OMP   SHARED(BgG,BgG_new,iG,iGm1) &
  !$OMP   SHARED(id,nb_mult_BTOG,tab_ba,ind_Basis,ind_Grid) &
  !$OMP   PRIVATE(ibb,ibbNew,lbb,li,nb,iqq,iqqi,iqqf) &
  !$OMP   NUM_THREADS(BasisTOGrid_maxth)
  !$OMP   DO SCHEDULE(DYNAMIC)
  DO iqq=1,nqq
    iqqi = (iqq-1)*nq
    iqqf = iqqi + nq

    DO ibb=1,nbb
      ibbNew = tab_ibbnew_AT_ibb(ibb) ! local variable

      lbb = ind_Basis(id)%i_TO_l(ibb)
      li = min( liq , LB-lbb )
      nb = tab_ba(li,id)%nb

      CALL RB_TO_RG_basis(BgG(iGm1)%RDP(ibbNew+1:ibbNew+nb,iqq),&
                          BgG_new(iG)%RDP(ibb,iqqi+1:iqqi+nq),tab_ba(liq,id))
    END DO


  END DO
  !$OMP   END DO
  !$OMP   END PARALLEL

  deallocate(tab_ibbnew_AT_ibb)


END DO
CALL dealloc_TabRDP_at_iG(BgG,iGm1)
END DO


CALL TabRDP2_TO_TabRDP1(BgG,BgG_new)
CALL dealloc_TabRDP(BgG_new)

IF (debug) THEN
  CALL Write_TabRDP(BgG)
  CALL SumSq_TabRDP(BgG)
  write(out_unit,*) 'END BbG_TO_BgG',id
  flush(out_unit)
END IF

END SUBROUTINE BbG_TO_BgG

SUBROUTINE Norm_OF_BgG(BgG,WSG,indGrid,tab_ba,D,LG)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE


TYPE(TypeRDP), allocatable    :: BgG(:)
real(kind=Rkind), allocatable :: WSG(:)

integer              :: D,LG
TYPE(TypeDInd)       :: indGrid
TYPE(basis)          :: tab_ba(0:LG,D)



integer               :: iG,i,nqq,iqq
real(kind=Rkind)      :: w,Norm
integer  :: tabiq(D),tabnq(D),tabl(D)

IF (BgG(1)%n2 /= 1) THEN
  write(out_unit,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unit,*) ' BgG is not completely on the G grid'
  STOP
END IF

IF (size(BgG) /= size(WSG) ) THEN
  write(out_unit,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unit,*) ' The number of SG are different'
  write(out_unit,*) ' size(BgG) /= size(WSG)',size(BgG),size(WSG)
  STOP
END IF

Norm = ZERO
DO iG=1,indGrid%MaxnD
  DO i=1,D
    tabl(i)  = indGrid%tab_ind(i,iG)
    tabnq(i) = get_nq_FROM_basis(tab_ba(tabl(i),i))
  END DO
  nqq = product(tabnq)

  DO iqq=1,nqq

    ! the weight with Smolyak coeficient
    CALL InD_TO_tabi(iqq,D,tabnq,tabiq)
    w = WSG(iG)
    DO i=1,D
      w = w * tab_ba(tabl(i),i)%w(tabiq(i))
    END DO
    Norm = Norm + BgG(iG)%RDP(1,iqq)**2 * w
  END DO
END DO

write(out_unit,*) 'Norm of BgB',Norm

!-------------------------------------------
!-------------------------------------------
END SUBROUTINE Norm_OF_BgG

SUBROUTINE sub_G_TO_B_new(RVecG,RVecB,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

real(kind=Rkind)   :: RVecG(:)
real(kind=Rkind)   :: RVecB(:)
integer                     :: D

real(kind=Rkind)            :: WSG(:)
integer                     :: LG(D),LB(D)
TYPE(TypeDInd), allocatable :: ind_Grid(:)
TYPE(TypeDInd), allocatable :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:maxval(LG),D)

TYPE(TypeRDP), allocatable  :: WPG(:)
TYPE(TypeRDP), allocatable  :: WPB(:)

integer :: id

integer :: nb_BG

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (debug) THEN
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  CALL time_perso('sub_G_TO_B')
  write(out_unit,*) '====================================='
  flush(out_unit)
  !write(out_unit,*) 'RVecG',RVecG
END IF



CALL Set_BgG_FOR_id(WPG,ind_Grid,ind_Basis,tab_ba,D,LG(1),0,.TRUE.)
IF (debug) CALL Size_TabRDP(WPG,nb_BG)
CALL Transfer_WP_TO_BgG(RVecG,WPG)

DO id=1,D

  CALL BgG_TO_BbG(WPG,WSG,ind_Grid,ind_Basis,tab_ba,D,LG(id),LB(id),id)
  IF (debug) CALL Size_TabRDP(WPG,nb_BG)

END DO


RVecB(:) = WPG(1)%RDP(:,1)

CALL dealloc_TabRDP(WPG)
CALL dealloc_TabRDP(WPB)

IF (debug) THEN
  !write(out_unit,*) 'RVecB',RVecB
  write(out_unit,*) '====================================='
  CALL time_perso('sub_G_TO_B')
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  flush(out_unit)
END IF

END SUBROUTINE sub_G_TO_B_new
SUBROUTINE sub_B_TO_G_new(RVecB,RVecG,ind_Grid,ind_Basis,tab_ba,D,LG,LB)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

real(kind=Rkind)   :: RVecG(:)
real(kind=Rkind)   :: RVecB(:)

integer            :: D
integer            :: LG(D),LB(D)
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:maxval(LG),D)

TYPE(TypeRDP), allocatable    :: WPG(:)

integer :: id
integer :: nb_BG

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

!CALL Check_mem()
IF (debug) THEN
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  CALL time_perso('sub_B_TO_G')
  write(out_unit,*) '====================================='
  !write(out_unit,*) 'RVecB',RVecB
  flush(out_unit)
END IF

CALL Set_BgG_FOR_id(WPG,ind_Grid,ind_Basis,tab_ba,D,LG(D),D,.TRUE.)
IF (debug) CALL Size_TabRDP(WPG,nb_BG)
WPG(1)%RDP(:,1) = RVecB(:)

DO id=D,1,-1

  CALL BbG_TO_BgG(WPG,ind_Grid,ind_Basis,tab_ba,D,LG(id),LB(id),id)
  IF (debug) CALL Size_TabRDP(WPG,nb_BG)

END DO
CALL Transfer_BgG_TO_WP(WPG,RVecG)
CALL dealloc_TabRDP(WPG)

IF (debug) THEN
  !write(out_unit,*) 'RVecG',RVecG
  write(out_unit,*) '====================================='
  CALL time_perso('sub_B_TO_G')
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  flush(out_unit)
END  IF

!CALL UnCheck_mem()

END SUBROUTINE sub_B_TO_G_new
SUBROUTINE sub_G_TO_B(RVecG,RVecB,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

real(kind=Rkind)   :: RVecG(:)
real(kind=Rkind)   :: RVecB(:)

real(kind=Rkind)            :: WSG(:)
integer                     :: D,LG,LB
TYPE(TypeDInd), allocatable :: ind_Grid(:)
TYPE(TypeDInd), allocatable :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:LG,D)

TYPE(TypeRDP), allocatable  :: WPG(:)
TYPE(TypeRDP), allocatable  :: WPB(:)

integer :: id

integer :: nb_BG

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (debug) THEN
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  CALL time_perso('sub_G_TO_B')
  write(out_unit,*) '====================================='
  flush(out_unit)
  write(out_unit,*) 'RVecG',RVecG
END IF



CALL Set_BgG_FOR_id(WPG,ind_Grid,ind_Basis,tab_ba,D,LG,0,.TRUE.)
IF (debug) CALL Size_TabRDP(WPG,nb_BG)
CALL Transfer_WP_TO_BgG(RVecG,WPG)

DO id=1,D

  CALL BgG_TO_BbG(WPG,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)
  IF (debug) CALL Size_TabRDP(WPG,nb_BG)

END DO


RVecB(:) = WPG(1)%RDP(:,1)

CALL dealloc_TabRDP(WPG)
CALL dealloc_TabRDP(WPB)

IF (debug) THEN
  write(out_unit,*) 'RVecB',RVecB
  write(out_unit,*) '====================================='
  CALL time_perso('sub_G_TO_B')
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  flush(out_unit)
END IF

END SUBROUTINE sub_G_TO_B
SUBROUTINE sub_B_TO_G(RVecB,RVecG,ind_Grid,ind_Basis,tab_ba,D,LG,LB)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

real(kind=Rkind)   :: RVecG(:)
real(kind=Rkind)   :: RVecB(:)

integer            :: D,LG,LB
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(basis)   :: tab_ba(0:LG,D)

TYPE(TypeRDP), allocatable    :: WPG(:)

integer :: id
integer :: nb_BG

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

!CALL Check_mem()
IF (debug) THEN
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  CALL time_perso('sub_B_TO_G_old')
  write(out_unit,*) '====================================='
  !write(out_unit,*) 'RVecB',RVecB
  flush(out_unit)
END IF

CALL Set_BgG_FOR_id(WPG,ind_Grid,ind_Basis,tab_ba,D,LG,D,.TRUE.)
IF (debug) CALL Size_TabRDP(WPG,nb_BG)
WPG(1)%RDP(:,1) = RVecB(:)

DO id=D,1,-1


  !CALL dealloc_TabRDP(WPG)
  !CALL Set_BgG_FOR_id(WPG,ind_Grid,ind_Basis,tab_ba,D,LG,id,.FALSE.)

  CALL BbG_TO_BgG(WPG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)
  IF (debug) write(out_unit,*) 'id',id

  IF (debug) CALL Size_TabRDP(WPG,nb_BG)

END DO

CALL Transfer_BgG_TO_WP(WPG,RVecG)
CALL dealloc_TabRDP(WPG)

IF (debug) THEN
  !write(out_unit,*) 'RVecG',RVecG
  write(out_unit,*) '====================================='
  CALL time_perso('sub_B_TO_G_old')
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  write(out_unit,*) '====================================='
  flush(out_unit)
END  IF
!stop
!CALL UnCheck_mem()

END SUBROUTINE sub_B_TO_G
RECURSIVE SUBROUTINE DerivOp_TO_RVecG_SGType2(RVecG,nq,ind_Grid,ind_Basis,tab_ba,D,LG,LB,tab_der)
USE EVR_system_m
USE mod_module_DInd
IMPLICIT NONE

  integer, intent(in)              :: nq
  real (kind=Rkind), intent(inout) :: RVecG(:)
  integer, optional                :: tab_der(2)

integer                     :: D,LG,LB
TYPE(TypeDInd), allocatable :: ind_Grid(:)
TYPE(TypeDInd), allocatable :: ind_Basis(:)
TYPE(basis)                 :: tab_ba(0:LG,D)

integer          :: tab_der_loc(2),dnba_ind(2)
integer          :: iq1_d,dnba_ind1(2),iq2_d,dnba_ind2(2)

integer            :: iG,i,nqq,li,iq1,nq1,iq,iq2,nq2,iq3,nq3,l2,iqd
integer            :: tabiq(D),tabl(D)

TYPE(TypeRDP),    allocatable :: WPG(:)


integer,          allocatable :: tabnq(:)
real(kind=Rkind), allocatable :: RGgG(:,:,:),Rg(:)
real(kind=Rkind), pointer     :: B3GG(:,:)

!-- for debuging --------------------------------------------------
integer :: err_mem,memory
character (len=*), parameter :: name_sub='DerivOp_TO_RVecG_SGType2'
logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
!------------------------------------------------------------------

IF (debug) THEN
  write(out_unit,*) 'BEGINNING ',name_sub
  write(out_unit,*) 'RvecG(:)',RvecG(:)
  flush(out_unit)
END IF

IF (present(tab_der)) THEN
  tab_der_loc(:) = tab_der(:)
ELSE
  tab_der_loc(:) = 0
END IF
WHERE (tab_der_loc < 0) tab_der_loc = 0

IF (debug) write(out_unit,*) 'tab_der_loc',tab_der_loc
!analysis of tab_der_loc
iq1_d = 0
iq2_d = 0
DO i=1,D
  dnba_ind(:) = tab_ba(0,i)%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))
  IF (count(dnba_ind(:) > 0) > 0) THEN
    IF (iq1_d == 0) THEN
      iq1_d = i
      dnba_ind1(:) = dnba_ind(:)
    ELSE
      iq2_d = i
      dnba_ind2(:) = dnba_ind(:)
    END IF
  END IF
END DO

IF (debug) write(out_unit,*) 'dnba_ind2(:)',dnba_ind2(:)
IF (debug) write(out_unit,*) 'iq1_d,iq2_d',iq1_d,iq2_d


IF (iq1_d == 0) THEN
  write(out_unit,*) ' ERROR in ',name_sub
  write(out_unit,*) ' iq1_d == 0: no derivative'
  write(out_unit,*) ' it should never appended!'
  write(out_unit,*) ' CHECK the fortran!!!'
  STOP
END IF

CALL Set_BgG_FOR_id(WPG,ind_Grid,ind_Basis,tab_ba,D,LG,0,.TRUE.)
CALL Transfer_WP_TO_BgG(RvecG,WPG)

!CALL Write_TabRDP(WPG)

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(WPG,D,ind_Grid,tab_ba,iq1_d,iq2_d,dnba_ind1,dnba_ind2,nb_mult_OpPsi) &
!$OMP   PRIVATE(iG,i,li,nqq,nq1,nq2,nq3,iq1,iq2,iq3,iqd,l2) &
!$OMP   PRIVATE(tabnq,dnba_ind,RGgG,B3GG) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)
allocate(tabnq(D))
!$OMP   DO SCHEDULE(DYNAMIC)
DO iG=1,size(WPG)
  DO i=1,D
    li       = ind_Grid(0)%tab_ind(i,iG)
    tabnq(i) = get_nq_FROM_basis(tab_ba(li,i))
  END DO
  nqq = product(tabnq)

  ! first derivative: iq1_d1
  iqd         = iq1_d
  dnba_ind(:) = dnba_ind1(:)

  IF (iqd == 1) THEN
    nq1 = 1
  ELSE
    nq1 = product(tabnq(1:iqd-1))
  END IF
  nq2 = tabnq(iqd)
  l2  = ind_Grid(0)%tab_ind(iqd,iG)
  IF (iqd == D) THEN
    nq3 = 1
  ELSE
    nq3 = product(tabnq(iqd+1:D))
  END IF

  IF (nq1*nq2*nq3 /= nqq) THEN
    STOP 'ERROR with nqi'
  END IF

  allocate(RGgG(nq1,nq2,nq3))
  RGgG(:,:,:) = reshape(WPG(iG)%RDP,[nq1,nq2,nq3])

  CALL Get3_MatdnRGG(tab_ba(l2,iqd),B3GG,dnba_ind)

  DO iq3=1,nq3
  DO iq1=1,nq1
    RGgG(iq1,:,iq3) = matmul(B3GG,RGgG(iq1,:,iq3))
  END DO
  END DO

  WPG(iG)%RDP = reshape(RGgG,[1,nqq])
  deallocate(RGgG)

  IF (iq2_d == 0) CYCLE ! second derivative: iq2_d, if needed
  iqd         = iq2_d
  dnba_ind(:) = dnba_ind2(:)

  IF (iqd == 1) THEN
    nq1 = 1
  ELSE
    nq1 = product(tabnq(1:iqd-1))
  END IF
  nq2 = tabnq(iqd)
  l2  = ind_Grid(0)%tab_ind(iqd,iG)
  IF (iqd == D) THEN
    nq3 = 1
  ELSE
    nq3 = product(tabnq(iqd+1:D))
  END IF

  IF (nq1*nq2*nq3 /= nqq) THEN
    STOP 'ERROR with nqi'
  END IF

  allocate(RGgG(nq1,nq2,nq3))
  RGgG(:,:,:) = reshape(WPG(iG)%RDP,[nq1,nq2,nq3])

  CALL Get3_MatdnRGG(tab_ba(l2,iqd),B3GG,dnba_ind)

  DO iq3=1,nq3
  DO iq1=1,nq1
    RGgG(iq1,:,iq3) = matmul(B3GG,RGgG(iq1,:,iq3))
  END DO
  END DO

  WPG(iG)%RDP = reshape(RGgG,[1,nqq])
  deallocate(RGgG)

END DO
!$OMP   END DO
deallocate(tabnq)
!$OMP   END PARALLEL

CALL Transfer_BgG_TO_WP(WPG,RvecG)
CALL dealloc_TabRDP(WPG)

 IF (debug) THEN
   write(out_unit,*) 'RvecG(:)',RvecG(:)
   write(out_unit,*) 'END ',name_sub
 END IF

END SUBROUTINE DerivOp_TO_RVecG_SGType2


END MODULE mod_basis_BtoG_GtoB_SGType2
