Program test
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real32,real64,real128,int32,int64
  IMPLICIT NONE

integer, parameter :: Rkind=real64

integer, parameter :: nb_Q0=6
real(kind=Rkind)   :: th,Q0(nb_Q0)
integer            :: i,nb_pts

integer                        :: nb,nq           ! numbers of basis functions and grid points
integer                        :: nb_vec          ! number of eigenvectors/eigenvalues
real (kind=Rkind), allocatable :: EigenVal(:)     ! eigenvalues.               EigenVal(i)
real (kind=Rkind), allocatable :: EigenVecB(:,:)  ! eigenvectors on the basis. EigenVecB(:,i)
real (kind=Rkind), allocatable :: EigenVecG(:,:)  ! eigenvectors on the grid.  EigenVecG(:,i)
real (kind=Rkind), allocatable :: RhoWeight(:)    ! rho(Q).Weight(Q), on the grid points.


Q0(:) = [0.900000,3.187000_Rkind, 2.179000_Rkind,0.000000_Rkind,3.141593_Rkind,0.000000_Rkind]

CALL init_EVR()

CALL get_nb_nq(nb,nq)
allocate(EigenVal(nb))
allocate(EigenVecB(nb,nb))
allocate(EigenVecG(nq,nb))
allocate(RhoWeight(nq))
write(OUTPUT_UNIT,*) 'END init_EVR'

CALL levels_EVR(EigenVal,EigenVecB,EigenVecG,RhoWeight,nb,nq,nb_vec)
write(OUTPUT_UNIT,*) 'nb_vec',nb_vec
write(OUTPUT_UNIT,*) 'EigenVal(:)',EigenVal(1:nb_vec)

write(OUTPUT_UNIT,*) 'modify Q0'
nb_pts = 10
DO i=1,2*nb_pts-1
  th = -1._Rkind+real(i,kind=8)/real(nb_pts,kind=8)
  Q0(:) = [th,3.187_Rkind, 2.179_Rkind,0._Rkind,3.141593_Rkind,0._Rkind]
  CALL Modify_TnumRefGeom_Q0(Q0,nb_Q0,.TRUE.)
  CALL levels_EVR(EigenVal,EigenVecB,EigenVecG,RhoWeight,nb,nq,nb_vec)
  write(OUTPUT_UNIT,*) 'EigenVal(:)',EigenVal(1:nb_vec)
END DO
CALL finalyze_EVR()
END
