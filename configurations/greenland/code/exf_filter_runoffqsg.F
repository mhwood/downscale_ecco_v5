#include "EXF_OPTIONS.h"

      SUBROUTINE EXF_FILTER_RUNOFFQSG(
     U                          arr, myThid )

C     ==================================================================
C     SUBROUTINE EXF_FILTER_RUNOFFQSG
C     ==================================================================
C
C     o this routine is desgined to apply a wet mask to the 
C     o runoff qsg field when using BOTH iceplume and shelfice
C     o the original implementation (EXF_FILTER_RL) was masking plumes 
C     o below the ice shelves
C
C     ==================================================================
C     SUBROUTINE EXF_FILTER_RUNOFFQSG
C     ==================================================================

      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
c#include "EXF_CONSTANTS.h"

C     == routine arguments ==
      _RL arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wateroricemask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid

C     == local variables ==
      INTEGER bi,bj
      INTEGER i,j,k

C     == end of interface ==

C     filter forcing field array

        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)

C         Make a mask which is 0 in the water
C         The original mask is 0 where there is an ice shelf

C       Initialize it with zeros first
        DO j = 1,sNy
          DO i = 1,sNx
            wateroricemask = 0. _d 0
          ENDDO
        ENDDO

C       Then, make it 1 where there is water at any depth
          DO j = 1,sNy
          DO i = 1,sNx
          DO k = 1,Nr
            IF ( maskC(i,j,k,bi,bj) .NE. 0. ) THEN
              wateroricemask = 1. _d 0
            ENDIF
          ENDDO
          ENDDO
          ENDDO

C         Now, use the new mask to mask out the Qsg
            DO j = 1,sNy
             DO i = 1,sNx
               IF ( wateroricemask(i,j,bi,bj) .EQ. 0. ) THEN
                 arr(i,j,bi,bj) = 0. _d 0
               ENDIF
             ENDDO
            ENDDO

         ENDDO
        ENDDO

      RETURN
      END
