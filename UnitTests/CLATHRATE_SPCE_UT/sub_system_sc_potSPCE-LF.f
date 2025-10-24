c
c================================================================
c    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
c================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Qcart,nb_cart,mole,
     *                    calc_ScalOp,pot_cplx)

      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE

c----- for the Tnum --------------------------------------
      TYPE (CoordType) :: mole


      integer           :: nb_be,nb_ScalOp,nb_cart
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qcart(nb_cart)

      real (kind=Rkind) :: pot0,im_pot0


      !write(6,*) 'size Qcart',size(Qcart),nb_cart


      IF (nb_be == 1 ) THEN
        call sc_sp(reshape(Qcart(1:6),(/3,2/)),mat_V(1,1))
        mat_V(1,1) = mat_V(1,1) / 219474.63d0 ! the conversion factor comes from the sc_sp subroutine
        mat_V(1,1) = min(mat_V(1,1),ONE)
        !mat_V(1,1) = ZERO
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Qcart)
        IF (calc_ScalOp) THEN
          CALL sub_dipole(mat_ScalOp(1,1,:),Qcart,mole)
        END IF
      ELSE
        write(6,*) ' ERROR in calc_op'
        write(6,*) ' more than ONE electronic surface is impossible'
        write(6,*) ' Rq: nb_be',nb_be
        STOP
      END IF

      RETURN
      END
C=======================================================================
C    sub hessian
C=======================================================================
      SUBROUTINE sub_hessian(h)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: h

      !-----------------------------------------------------------------
       h = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine sub_hessian MUST be make'
      !-----------------------------------------------------------------


      END SUBROUTINE sub_hessian
C================================================================
C    subroutine calculant le gradient
C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)
      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE

c----- for the Tnum --------------------------------------
      TYPE (CoordType) :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: 
     *                d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qsym0(mole%nb_var)
      real (kind=Rkind) :: step
      logical       :: deriv,num


      d0g = 0.d0


      END
C================================================================
C    fonction pot_rest(x)
C================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: pot_rest
       real (kind=Rkind) :: Qact(1)
       integer       :: nb_inact2n
       real (kind=Rkind) :: Delta_Qact(nb_inact2n)

       pot_rest = 0.d0

       END
C================================================================
C    fonction im_pot0(x)
C================================================================
      FUNCTION im_pot0(Qsym0)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: im_pot0
       real (kind=Rkind) :: Qsym0(1)
       real (kind=Rkind) :: z

       z = 0.d0

       im_pot0 = z

       RETURN
       END
C================================================================
C    subroutine calculant la matrice hessienne
C    en fonction de x=cos(theta)
C================================================================
       SUBROUTINE d0d1d2_h(d0h,d1h,d2h,
     *                     Qsym0,mole,deriv,num,step)

      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE

c----- for the Tnum --------------------------------------
      TYPE (CoordType) :: mole

      
      real (kind=Rkind) ::  Qsym0(mole%nb_var)


      real (kind=Rkind) :: step
      logical deriv,num

      real (kind=Rkind) :: d0h
      real (kind=Rkind) :: d1h
      real (kind=Rkind) :: d2h

c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      IF (debug) THEN
      write(6,*)
      write(6,*) 'BEGINNING d0d1d2_h'
      END IF
c---------------------------------------------------------------------


      STOP 'd0d1d2_h'

      END
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
C================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_qsym,
     *                        d0req,d1req,d2req,d3req,
     *                        Qsym0,mole,nderiv)

      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE

c----- for the Tnum --------------------------------------
      TYPE (CoordType) :: mole

       integer i_qsym
       real (kind=Rkind) ::  Qsym0(mole%nb_var)

       integer nderiv

       real (kind=Rkind) ::  d0req
       real (kind=Rkind) ::  d1req
       real (kind=Rkind) ::  d2req
       real (kind=Rkind) ::  d3req


c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING d0d1d2d3_Qeq'
        write(6,*) 'i_qsym',i_qsym
      END IF
c---------------------------------------------------------------------

      STOP 'd0d1d2d3_Qeq'

      RETURN
      END

C================================================================
C    analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
c    for the variable iq
C================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

       integer :: iq,nb_act
       real (kind=Rkind) ::  Qact(nb_act)
       integer :: nderiv,it
       TYPE (Type_dnS)   :: dnQflex
       STOP 'dnQflex'
       END
c
C================================================================
c    dipole read
C================================================================
      SUBROUTINE sub_dipole(dip,Q,mole)
      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE

c----- for the Tnum --------------------------------------
      TYPE (CoordType) :: mole

      real (kind=Rkind) :: Q(mole%nb_var)
      real (kind=Rkind) :: DQ(mole%nb_var)
      real (kind=Rkind) :: dip(3)

      real (kind=Rkind), parameter :: Qref(3) = 
     *  (/ 1.869713_Rkind,1.869713_Rkind,1.745803_Rkind /)

      real (kind=Rkind), parameter :: MuX = 0.51523633071858288d0
      real (kind=Rkind), parameter :: MuY = ZERO
      real (kind=Rkind), parameter :: MuZ = 0.43212653883837798d0
      real (kind=Rkind), parameter :: dMuX(3) = (/
     x   -2.89121043244283428d-002,-0.19662504130644246d0,
     x    0.11731395143422289d0/)
      real (kind=Rkind), parameter :: dMuY(3) = ZERO
      real (kind=Rkind), parameter :: dMuZ(3) = (/
     z -0.19456328715229212d0,    5.40468620138319714d-003,
     z -0.34043381267231121d0/)

      DQ(:) = Q(:)-Qref(:)


      dip(1) = MuX + dot_product(DQ,dMuX)
      dip(2) = MuY + dot_product(DQ,dMuY)
      dip(3) = MuZ + dot_product(DQ,dMuZ)

      !dip = dip * 0.5291772084d0

      END
c***********************************************************************
      subroutine sc_sp(h2coor,vpot)
c***********************************************************************
c     routine to evaluate the interaction potential for one H2 molecule
c     inside the small cage cavity (20 water molecules)
c     of the sII structure of clathrate hydrate
c
c     input :  array h2coor(3,2)  [bohr] (H atoms cartesian coordinates
c                                          of the H2 molecule)
c     output:  vpot H2-large_cage [cm-1]
c-------------------
c     bigcag(3,ntot): coordinates
c                     for the water molecules in the large cage [bohr]
c***********************************************************************
      implicit real*8 (a-h,o-z)
      parameter (two=2.0d0,pi=3.14159265358979323846d0)
c     nwater = number of water molecules for the cage
c     ntot   = total number of atoms of the cage
      parameter (nwater=20,ntot=nwater*3,au2cm=219474.63d0) 
      dimension coord(3,5),smallc(3,ntot)
      dimension h2coor(3,2)
cccccc
      smallc(1, 1)=   -0.31649d0
      smallc(2, 1)=    0.38003d0
      smallc(3, 1)=   -7.07291d0
      smallc(1, 2)=    3.70808d0
      smallc(2, 2)=   -2.65146d0
      smallc(3, 2)=   -5.67704d0
      smallc(1, 3)=    6.94833d0
      smallc(2, 3)=    0.28393d0
      smallc(3, 3)=   -2.77947d0
      smallc(1, 4)=    4.85546d0
      smallc(2, 4)=    5.15550d0
      smallc(3, 4)=   -2.41820d0
      smallc(1, 5)=    0.34577d0
      smallc(2, 5)=    5.17496d0
      smallc(3, 5)=   -5.09665d0
      smallc(1, 6)=    6.94765d0
      smallc(2, 6)=   -1.90366d0
      smallc(3, 6)=    2.02974d0
      smallc(1, 7)=    4.78465d0
      smallc(2, 7)=    1.66045d0
      smallc(3, 7)=    5.21254d0
      smallc(1, 8)=    3.53305d0
      smallc(2, 8)=    6.04447d0
      smallc(3, 8)=    2.61917d0
      smallc(1, 9)=    0.34023d0
      smallc(2, 9)=   -0.38558d0
      smallc(3, 9)=    7.05546d0
      smallc(1,10)=   -3.68434d0
      smallc(2,10)=    2.64591d0
      smallc(3,10)=    5.65960d0
      smallc(1,11)=   -6.92460d0
      smallc(2,11)=   -0.28948d0
      smallc(3,11)=    2.76202d0
      smallc(1,12)=   -4.83172d0
      smallc(2,12)=   -5.16106d0
      smallc(3,12)=    2.40075d0
      smallc(1,13)=   -0.32204d0
      smallc(2,13)=   -5.18051d0
      smallc(3,13)=    5.07920d0
      smallc(1,14)=   -4.76091d0
      smallc(2,14)=   -1.66601d0
      smallc(3,14)=   -5.22999d0
      smallc(1,15)=   -6.92391d0
      smallc(2,15)=    1.89811d0
      smallc(3,15)=   -2.04719d0
      smallc(1,16)=   -3.50932d0
      smallc(2,16)=   -6.05003d0
      smallc(3,16)=   -2.63661d0
      smallc(1,17)=    1.76223d0
      smallc(2,17)=   -6.66343d0
      smallc(3,17)=   -2.91488d0
      smallc(1,18)=   -1.73849d0
      smallc(2,18)=    6.65787d0
      smallc(3,18)=    2.89744d0
      smallc(1,19)=    3.76898d0
      smallc(2,19)=   -6.16184d0
      smallc(3,19)=    1.94674d0
      smallc(1,20)=   -3.74524d0
      smallc(2,20)=    6.15629d0
      smallc(3,20)=   -1.96419d0
      smallc(1,21)=   -0.08391d0
      smallc(2,21)=    2.06393d0
      smallc(3,21)=   -6.37888d0
      smallc(1,22)=   -0.38327d0
      smallc(2,22)=    0.63072d0
      smallc(3,22)=   -8.89059d0
      smallc(1,23)=    2.26809d0
      smallc(2,23)=   -1.61376d0
      smallc(3,23)=   -6.14702d0
      smallc(1,24)=    2.98751d0
      smallc(2,24)=   -4.03968d0
      smallc(3,24)=   -4.71534d0
      smallc(1,25)=    5.81406d0
      smallc(2,25)=   -0.74362d0
      smallc(3,25)=   -3.79378d0
      smallc(1,26)=    8.55826d0
      smallc(2,26)=    0.20693d0
      smallc(3,26)=   -3.65894d0
      smallc(1,27)=    5.57003d0
      smallc(2,27)=    3.46539d0
      smallc(3,27)=   -2.48266d0
      smallc(1,28)=    4.41648d0
      smallc(2,28)=    5.40842d0
      smallc(3,28)=   -0.65338d0
      smallc(1,29)=    1.86723d0
      smallc(2,29)=    5.18989d0
      smallc(3,29)=   -4.06892d0
      smallc(1,30)=   -1.01662d0
      smallc(2,30)=    5.52546d0
      smallc(3,30)=   -3.91669d0
      smallc(1,31)=    6.94789d0
      smallc(2,31)=   -1.14343d0
      smallc(3,31)=    0.35843d0
      smallc(1,32)=    8.71272d0
      smallc(2,32)=   -2.27583d0
      smallc(3,32)=    2.37222d0
      smallc(1,33)=    5.50892d0
      smallc(2,33)=    0.38643d0
      smallc(3,33)=    4.10640d0
      smallc(1,34)=    3.23140d0
      smallc(2,34)=    0.90630d0
      smallc(3,34)=    5.83701d0
      smallc(1,35)=    3.97118d0
      smallc(2,35)=    4.50982d0
      smallc(3,35)=    3.52700d0
      smallc(1,36)=    4.50123d0
      smallc(2,36)=    7.36725d0
      smallc(3,36)=    3.44631d0
      smallc(1,37)=   -1.07313d0
      smallc(2,37)=    0.67903d0
      smallc(3,37)=    6.56527d0
      smallc(1,38)=    0.30262d0
      smallc(2,38)=   -0.39329d0
      smallc(3,38)=    8.89117d0
      smallc(1,39)=   -4.79969d0
      smallc(2,39)=    1.68568d0
      smallc(3,39)=    4.56178d0
      smallc(1,40)=   -3.06076d0
      smallc(2,40)=    4.01515d0
      smallc(3,40)=    4.60718d0
      smallc(1,41)=   -6.20152d0
      smallc(2,41)=   -1.97259d0
      smallc(3,41)=    2.63721d0
      smallc(1,42)=   -8.59075d0
      smallc(2,42)=   -0.55116d0
      smallc(3,42)=    3.48779d0
      smallc(1,43)=   -3.24481d0
      smallc(2,43)=   -5.17079d0
      smallc(3,43)=    3.32426d0
      smallc(1,44)=   -4.35616d0
      smallc(2,44)=   -5.47084d0
      smallc(3,44)=    0.65457d0
      smallc(1,45)=   -0.08946d0
      smallc(2,45)=   -3.49661d0
      smallc(3,45)=    5.77323d0
      smallc(1,46)=   -0.37321d0
      smallc(2,46)=   -6.28524d0
      smallc(3,46)=    6.54489d0
      smallc(1,47)=   -3.20767d0
      smallc(2,47)=   -0.91185d0
      smallc(3,47)=   -5.85446d0
      smallc(1,48)=   -5.48518d0
      smallc(2,48)=   -0.39198d0
      smallc(3,48)=   -4.12384d0
      smallc(1,49)=   -6.92415d0
      smallc(2,49)=    1.13788d0
      smallc(3,49)=   -0.37587d0
      smallc(1,50)=   -8.68898d0
      smallc(2,50)=    2.27028d0
      smallc(3,50)=   -2.38967d0
      smallc(1,51)=   -3.94744d0
      smallc(2,51)=   -4.51537d0
      smallc(3,51)=   -3.54444d0
      smallc(1,52)=   -4.47749d0
      smallc(2,52)=   -7.37280d0
      smallc(3,52)=   -3.46376d0
      smallc(1,53)=   -0.05385d0
      smallc(2,53)=   -6.44310d0
      smallc(3,53)=   -2.75819d0
      smallc(1,54)=    2.40286d0
      smallc(2,54)=   -6.47979d0
      smallc(3,54)=   -1.20400d0
      smallc(1,55)=    0.07759d0
      smallc(2,55)=    6.43755d0
      smallc(3,55)=    2.74075d0
      smallc(1,56)=   -2.37912d0
      smallc(2,56)=    6.47424d0
      smallc(3,56)=    1.18655d0
      smallc(1,57)=    4.83460d0
      smallc(2,57)=   -4.66786d0
      smallc(3,57)=    2.00789d0
      smallc(1,58)=    2.35681d0
      smallc(2,58)=   -5.77211d0
      smallc(3,58)=    3.05359d0
      smallc(1,59)=   -4.84346d0
      smallc(2,59)=    4.68511d0
      smallc(3,59)=   -1.99286d0
      smallc(1,60)=   -4.81635d0
      smallc(2,60)=    7.54743d0
      smallc(3,60)=   -2.50152d0
cccccc
c     calculate total potential inside the cage + h2
      vpot=0.0d0
      do i=1,nwater
c        oxygen
         coord(1,1)=smallc(1,i)
         coord(2,1)=smallc(2,i)
         coord(3,1)=smallc(3,i)
c        first hydrogen of water molecules
         coord(1,2)=smallc(1,nwater+2*i-1)
         coord(2,2)=smallc(2,nwater+2*i-1)
         coord(3,2)=smallc(3,nwater+2*i-1)
c        second hydrogen of water molecules
         coord(1,3)=smallc(1,nwater+2*i)
         coord(2,3)=smallc(2,nwater+2*i)
         coord(3,3)=smallc(3,nwater+2*i)
c     H2 molecule at the center of the cage
         coord(1,4)=h2coor(1,1)
         coord(2,4)=h2coor(2,1)
         coord(3,4)=h2coor(3,1)
         coord(1,5)=h2coor(1,2)
         coord(2,5)=h2coor(2,2)
         coord(3,5)=h2coor(3,2)
c        now call SPC/E subroutine for WATER-H2:
c        subroutine oh2_h2(zcoord,Eint) 
         call oh2_h2(coord,vfit)
         vpot=vpot+vfit
      enddo
      return 
      end
c***********************************************************************
      subroutine oh2_h2(zcoord,Eint)
c***********************************************************************
c     routine to evaluate to interaction potential between
c     one H2 molecule and one H2O molecule within the SPC/E model,
c     see Alavi et al. J. Chem. Phys. 123, 024507 (2005)
c     and Berendsen et al. J. Phys. Chem. 91, 6269 (1987)
c
c     distances in [a.u.]
c     angles    in [radians]
c     Eint      in [cm-1]
c***********************************************************************
      implicit real*8 (a-h,o-z)
      dimension zcoord(3,5),x(3,6),d(9)
c     zcoord(i,1) = Oxygen coords
c     zcoord(i,2) = H atom coords
c     zcoord(i,3) = H atom coords
c     zcoord(i,4) = H2 first atom coords
c     zcoord(i,5) = H2 second atom coords
      parameter(sigma=3.10134d0,epsi=0.430624d0,
     &       calkcm=83.592d0,ang2au=0.529177d0,au2cm=219474.63d0,
     &       hwater=0.4238d0,oxygen=-0.8476d0,hcm=-0.9864d0,h2=0.4932d0)
      pi=dacos(-1.0d0)
      degree=pi/180.0d0
cccccc
      do i=1,9
         d(i)=0.0d0
      enddo
c     x(i,1) are the coordinates of the Oxygen atom
      x(1,1)=zcoord(1,1)
      x(2,1)=zcoord(2,1)
      x(3,1)=zcoord(3,1)
c     H atoms in WATER molecule
      x(1,2)=zcoord(1,2)
      x(2,2)=zcoord(2,2)
      x(3,2)=zcoord(3,2)
      x(1,3)=zcoord(1,3)
      x(2,3)=zcoord(2,3)
      x(3,3)=zcoord(3,3)
c
      x(1,4)=zcoord(1,4)
      x(2,4)=zcoord(2,4)
      x(3,4)=zcoord(3,4)
      x(1,5)=zcoord(1,5)
      x(2,5)=zcoord(2,5)
      x(3,5)=zcoord(3,5)
c     x(i,6) are the coordinates of the H_2 center of mass
      x(1,6)=(zcoord(1,4)+zcoord(1,5))/2.0d0
      x(2,6)=(zcoord(2,4)+zcoord(2,5))/2.0d0
      x(3,6)=(zcoord(3,4)+zcoord(3,5))/2.0d0
c
      R=0.0d0
      do j=1,3
         R=R+(x(j,1)-x(j,6))**2
      enddo
      R=dsqrt(R)
      do j=1,3
         d(1)=d(1)+(x(j,1)-x(j,4))**2
         d(2)=d(2)+(x(j,1)-x(j,5))**2
         d(3)=d(3)+(x(j,1)-x(j,6))**2
         d(4)=d(4)+(x(j,2)-x(j,4))**2
         d(5)=d(5)+(x(j,2)-x(j,5))**2
         d(6)=d(6)+(x(j,2)-x(j,6))**2
         d(7)=d(7)+(x(j,3)-x(j,4))**2
         d(8)=d(8)+(x(j,3)-x(j,5))**2
         d(9)=d(9)+(x(j,3)-x(j,6))**2
      enddo
      do i=1,9
         d(i)=dsqrt(d(i))
      enddo
      v_lj=(sigma/(R*ang2au))**12-(sigma/(R*ang2au))**6
c     v_lj now is kcalmol-1: below
      v_lj=4.0d0*epsi*v_lj
c     v_lj now is cm-1: below
      v_lj=v_lj*calkcm
c     v_coul now is in a.u.
      v_coul=oxygen*h2/d(1)+h2*oxygen/d(2)+hcm*oxygen/d(3)+
     &       hwater*h2/d(4)+hwater*h2/d(5)+hwater*hcm/d(6)+
     &       hwater*h2/d(7)+hwater*h2/d(8)+hwater*hcm/d(9)
c     v_coul now is in cm-1: below
      v_coul=v_coul*au2cm
c     Vpot in cm-1: below
      Eint= v_lj + v_coul
      return
      end
