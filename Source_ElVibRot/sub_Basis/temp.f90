SUBROUTINE OpRV_TO_RV_basis(Vin,Vout,lBin,lBout,base,der)
  USE mod_system
  IMPLICIT NONE
    
  real (kind=Rkind), intent(in)           :: Vin(:)
  real (kind=Rkind), intent(inout)        :: Vout(:)
  logical,           intent(in)           :: lBin,lBout

  TYPE (basis),      intent(in), target   :: base
  integer,           intent(in), optional :: der(2)

  integer :: nb,nq,n_Vin,n_Vout,der_basis(2)
  real (kind=Rkind), pointer       :: Mat(:,:)

  !----- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------
  nq     = get_nq_FROM_basis(base)
  nb     = get_nb_FROM_basis(base)
  n_Vin  = size(Vin)
  n_Vout = size(Vout)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING OpRV_TO_RV_basis'
    write(out_unitp,*) 'nb,nq',nb,nq
    write(out_unitp,*) 'n_Vin',n_Vin
    write(out_unitp,*) 'n_Vout',n_Vout
    write(out_unitp,*) 'Vin',Vin(:)
  END IF
  !-----------------------------------------------------------

  IF (NewBasisEl .AND. base%ndim == 0) THEN
    G(:) = B
  ELSE
    IF (nb_basis < nb_B .OR. nq /= size(G)) THEN
      write(out_unitp,*) ' ERROR in RB_TO_RG_basis'
      write(out_unitp,*) ' nb_basis is inconsistent with nb_B',nb_basis,nb_B
      write(out_unitp,*) ' nq is inconsistent with size(G)',nq,size(G)
      STOP ' ERROR in RB_TO_RG_basis: inconsistent sizes'
    END IF
    nb_mult_BTOG = nb_mult_BTOG + int(nb_B*nq,kind=ILkind)

    IF (present(der)) THEN
      der_basis = base%Tabder_Qdyn_TO_Qbasis(der(:))
    ELSE
      der_basis = 0
    END IF

    IF (base%dnBBRep_done) THEN
      IF (der_basis(1) == 0 .AND. der_basis(2) == 0) THEN
        Binter = B
      ELSE IF (der_basis(1) > 0 .AND. der_basis(2) == 0) THEN
        Mat => base%dnRBB%d1(:,:,der_basis(1))
        Binter = matmul(Mat,B)
      ELSE IF (der_basis(1) == 0 .AND. der_basis(2) > 0) THEN
        Mat => base%dnRBB%d1(:,:,der_basis(2))
        Binter = matmul(Mat,B)
      ELSE ! der_basis(1) > 0 and der_basis(2) > 0
        Mat => base%dnRBB%d2(:,:,der_basis(1),der_basis(2))
        Binter = matmul(Mat,B)
      END IF

      Mat => base%dnRGB%d0(:,1:nb_B)
      G = matmul(Mat,Binter)

      deallocate(Binter)

    ELSE
      IF (der_basis(1) == 0 .AND. der_basis(2) == 0) THEN
        Mat => base%dnRGB%d0(:,1:nb_B)
      ELSE IF (der_basis(1) > 0 .AND. der_basis(2) == 0) THEN
        Mat => base%dnRGB%d1(:,1:nb_B,der_basis(1))
      ELSE IF (der_basis(1) == 0 .AND. der_basis(2) > 0) THEN
        Mat => base%dnRGB%d1(:,1:nb_B,der_basis(2))
      ELSE ! der_basis(1) > 0 and der_basis(2) > 0
        Mat => base%dnRGB%d2(:,1:nb_B,der_basis(1),der_basis(2))
      END IF

      G = matmul(Mat,B)
    END IF

  END IF

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'Vout',Vout(:)
    write(out_unitp,*) 'END OpRV_TO_RV_basis'
  END IF
  !-----------------------------------------------------------
END SUBROUTINE OpRV_TO_RV_basis