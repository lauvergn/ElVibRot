      SUBROUTINE SECOND( T )
*
      REAL       T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the system function ETIME.
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..
c     REAL               ETIME
c     EXTERNAL           ETIME
*     ..
*     .. Executable Statements ..
*

      INTEGER, save :: count(2)=(/0,0/),count_rate,count_max

c     T1 = ETIME( TARRAY )
c     T  = TARRAY( 1 )

      CALL SYSTEM_CLOCK(count(1),count_rate,count_max)
      T  = real(count(1)-count(2))/real(count_rate)
      count(2) = count(1)
      

      RETURN
*
*     End of SECOND
*
      END
