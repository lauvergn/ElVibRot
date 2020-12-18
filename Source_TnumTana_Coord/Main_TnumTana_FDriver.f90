!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
 PROGRAM Main_TnumTana_FDriver
 IMPLICIT NONE

  integer, parameter :: nt=10**4
  integer, parameter :: nat=5
  real (kind=8) :: Qact(3*nat-6),Qcart(3*nat)
  integer :: InputUnit,OutputUnit


  integer :: i
  character (len=*), parameter :: name_sub='Main_TnumTana'

  open(newunit=OutputUnit,file='res_driver')
  open(newunit=InputUnit,file='dat_driver')

  CALL Init_InputUnit_Driver(InputUnit)
  CALL Init_OutputUnit_Driver(OutputUnit)
  write(6,*) 'InputUnit',InputUnit
  write(6,*) 'OutputUnit',OutputUnit


  Qact(:) = 0.5d0
  write(6,*) 'Qact (initial values)',Qact
  CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
  CALL cart_TO_Qact(Qact,size(Qact),Qcart,size(Qcart))
  write(6,*) 'Qact (from cart_TO_Qact)',Qact

  write(6,*) 'Beginning loop:',nt
 !$OMP   PARALLEL &
 !$OMP   DEFAULT(NONE) &
 !$OMP   PRIVATE(i,Qact,Qcart)

 !$OMP   DO SCHEDULE(STATIC)
  DO i=1,nt
    IF (mod(i,100) == 0) write(6,'(".")',advance='no')
    Qact(1) = 0.5d0 + real(i,kind=8)*0.001d0
    CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
  END DO
 !$OMP   END DO
 !$OMP   END PARALLEL
  write(6,*)
  write(6,*) 'END loop'

  DO i=1,3*nat,3
    write(6,*) (i-1)/3+1,Qcart(i:i+2)
  END DO



 END PROGRAM Main_TnumTana_FDriver
