*
*----------------------------------------------------------------------*
*                                                                      *
C     Exact Riemann Solver for the Time-Dependent                      *
C     One Dimensional Euler Equations                                  *
*                                                                      *
C     Name of program: HE-E1RPEX                                       *
*                                                                      *
C     Purpose: to solve the Riemann problem exactly,                   *
C              for the time dependent one dimensional                  *
C              Euler equations for an ideal gas                        *
*                                                                      *
C     Input  file: e1rpex.ini                                          *
C     Output file: e1rpex.out (exact solution)                         *
*                                                                      *
C     Programer: E. F. Toro                                            *
*                                                                      *
C     Last revision: 31st May 1999                                     *
*                                                                      *
C     Theory is found in Ref. 1, Chapt. 4 and in original              *
C     references therein                                               *
*                                                                      *
C     1. Toro, E. F., "Riemann Solvers and Numerical                   *
C                      Methods for Fluid Dynamics"                     *
C                      Springer-Verlag, 1997                           *
C                      Second Edition, 1999                            *
*                                                                      *
C     This program is part of                                          *
*                                                                      *
C     NUMERICA                                                         *
C     A Library of Source Codes for Teaching,                          *
C     Research and Applications,                                       *
C     by E. F. Toro                                                    *
C     Published by NUMERITEK LTD, 1999                                 *
C     Website: www.numeritek.com                                       *
*                                                                      *
*----------------------------------------------------------------------*
*
      IMPLICIT NONE
*
C     Declaration of variables:
*
      INTEGER I, CELLS
*
      REAL    GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,
     & DL, UL, PL, CL, DR, UR, PR, CR,
     & DIAPH1, DOMLEN, DS, DX, PM, PSCALE, PS, S,
     & TIMEOU, UM, US, XPOS
*
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
*
C     Input variables
*
C     DOMLEN   : Domain length
C     DIAPH1   : Position of diaphragm 1
C     CELLS    : Number of computing cells
C     GAMMA    : Ratio of specific heats
C     TIMEOU   : Output time
C     DL       : Initial density  on left state
C     UL       : Initial velocity on left state
C     PL       : Initial pressure on left state
C     DR       : Initial density  on right state
C     UR       : Initial velocity on right state
C     PR       : Initial pressure on right state
C     PSCALE   : Normalising constant
*
C      OPEN(UNIT = 1,FILE = 'EulerClas.ini',STATUS = 'UNKNOWN')
*
C     Initial data and parameters are read in
*

      WRITE(*,*) ' ADIABATIC INDEX OF THE GAS: ' 
      READ (*,*)   GAMMA

      WRITE(*,*) ' TIME FOR THE SOLUTION: ' 
      READ (*,*)   TIMEOU

C     -----
C     LEFT STATE
C     -----

      WRITE(*,*) ' -- LEFT STATE -- ' 
      WRITE(*,*) '      PRESSURE     : ' 
      READ (*,*) PL 
      WRITE(*,*) '      DENSITY      : ' 
      READ (*,*) DL 
      WRITE(*,*) '      FLOW VELOCITY: ' 
      READ (*,*) UL

C     ------ 
C     RIGHT STATE 
C     ------

      WRITE(*,*) ' -- RIGHT STATE -- ' 
      WRITE(*,*) '      PRESSURE     : ' 
      READ (*,*) PR 
      WRITE(*,*) '      DENSITY      : ' 
      READ (*,*) DR 
      WRITE(*,*) '      FLOW VELOCITY: ' 
      READ (*,*) UR

C     ------ 
C     RIGHT STATE 
C     ------

      WRITE(*,*) '      DOMAIN LENGTH: ' 
      READ (*,*) DOMLEN
      WRITE(*,*) '      INTERFACE: ' 
      READ (*,*) DIAPH1
      WRITE(*,*) '      NUMBER OF CELLS: ' 
      READ (*,*) CELLS
      WRITE(*,*) '      FLOW VELOCITY: ' 
      READ (*,*) UR
      WRITE(*,*) '      FLOW VELOCITY: ' 
      READ (*,*) UR
      WRITE(*,*) '      FLOW VELOCITY: ' 
      READ (*,*) UR

      PSCALE = 1.0
*
      CLOSE(1)
*
C     Compute gamma related constants
*
      G1 = (GAMMA - 1.0)/(2.0*GAMMA)
      G2 = (GAMMA + 1.0)/(2.0*GAMMA)
      G3 = 2.0*GAMMA/(GAMMA - 1.0)
      G4 = 2.0/(GAMMA - 1.0)
      G5 = 2.0/(GAMMA + 1.0)
      G6 = (GAMMA - 1.0)/(GAMMA + 1.0)
      G7 = (GAMMA - 1.0)/2.0
      G8 = GAMMA - 1.0
*
C     Compute sound speeds
*
      CL = SQRT(GAMMA*PL/DL)
      CR = SQRT(GAMMA*PR/DR)
*
C     The pressure positivity condition is tested for
*
      IF(G4*(CL+CR).LE.(UR-UL))THEN
*
C        The initial data is such that vacuum is generated.
C        Program stopped.
*
         WRITE(6,*)
         WRITE(6,*)'***Vacuum is generated by data***'
         WRITE(6,*)'***Program stopped***'
         WRITE(6,*)
*
         STOP
      ENDIF
*
C     Exact solution for pressure and velocity in star
C     region is found
*
      CALL STARPU(PM, UM, PSCALE)
*
      DX = DOMLEN/REAL(CELLS)
*
C     Complete solution at time TIMEOU is found
*
      OPEN(UNIT = 1,FILE = 'sol-RiemannClas.dat',STATUS = 'UNKNOWN')
*
      DO 10 I = 1, CELLS
*
         XPOS = (REAL(I) - 0.5)*DX
         S    = (XPOS - DIAPH1)/TIMEOU
*
C        Solution at point (X,T) = ( XPOS - DIAPH1,TIMEOU)
C        is found
*
         CALL SAMPLE(PM, UM, S, DS, US, PS)
*
C        Exact solution profiles are written to e1rpex.out.
*
         WRITE(1, 20)XPOS, DS, PS/PSCALE, US, PS/DS/G8/PSCALE
*
 10   CONTINUE
*
      CLOSE(1)
*
 20   FORMAT(5(F14.6, 2X))
*
      END
*
*----------------------------------------------------------------------*
*
      SUBROUTINE STARPU(P, U, PSCALE)
*
      IMPLICIT NONE
*
C     Purpose: to compute the solution for pressure and
C              velocity in the Star Region
*
C     Declaration of variables
*
      INTEGER I, NRITER
*
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,
     & CHANGE, FL, FLD, FR, FRD, P, POLD, PSTART,
     & TOLPRE, U, UDIFF, PSCALE
*
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      DATA TOLPRE, NRITER/1.0E-06, 20/
*
C     Guessed value PSTART is computed
*
      CALL GUESSP(PSTART)
*
      POLD  = PSTART
      UDIFF = UR - UL
*
      WRITE(6,*)'----------------------------------------'
      WRITE(6,*)'   Iteration number      Change  '
      WRITE(6,*)'----------------------------------------'
*
      DO 10 I = 1, NRITER
*
         CALL PREFUN(FL, FLD, POLD, DL, PL, CL)
         CALL PREFUN(FR, FRD, POLD, DR, PR, CR)
         P      = POLD - (FL + FR + UDIFF)/(FLD + FRD)
         CHANGE = 2.0*ABS((P - POLD)/(P + POLD))
         WRITE(6, 30)I, CHANGE
         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(P.LT.0.0)P = TOLPRE
         POLD  = P
*
 10   CONTINUE
*
      WRITE(6,*)'Divergence in Newton-Raphson iteration'
*
 20   CONTINUE
*
C     Compute velocity in Star Region
*
      U = 0.5*(UL + UR + FR - FL)
*
      WRITE(6,*)'---------------------------------------'
      WRITE(6,*)'   Pressure        Velocity'
      WRITE(6,*)'---------------------------------------'
      WRITE(6,40)P/PSCALE, U
      WRITE(6,*)'---------------------------------------'
*
 30   FORMAT(5X, I5,15X, F12.7)
 40   FORMAT(2(F14.6, 5X))
*
      END
*
*----------------------------------------------------------------------*
*
      SUBROUTINE GUESSP(PM)
*
C     Purpose: to provide a guessed value for pressure
C              PM in the Star Region. The choice is made
C              according to adaptive Riemann solver using
C              the PVRS, TRRS and TSRS approximate
C              Riemann solvers. See Sect. 9.5 of Chapt. 9
C              of Ref. 1
*
      IMPLICIT NONE
*
C     Declaration of variables
*
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,
     & GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,
     & CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,
     & PTL, PTR, QMAX, QUSER, UM
*
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
*
      QUSER = 2.0
*
C     Compute guess pressure from PVRS Riemann solver
*
      CUP  = 0.25*(DL + DR)*(CL + CR)
      PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
      PPV  = MAX(0.0, PPV)
      PMIN = MIN(PL,  PR)
      PMAX = MAX(PL,  PR)
      QMAX = PMAX/PMIN
*
      IF(QMAX.LE.QUSER.AND.
     & (PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
*
C        Select PVRS Riemann solver
*
         PM = PPV
      ELSE
         IF(PPV.LT.PMIN)THEN
*
C           Select Two-Rarefaction Riemann solver
*
            PQ  = (PL/PR)**G1
            UM  = (PQ*UL/CL + UR/CR +
     & G4*(PQ - 1.0))/(PQ/CL + 1.0/CR)
            PTL = 1.0 + G7*(UL - UM)/CL
            PTR = 1.0 + G7*(UM - UR)/CR
            PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
         ELSE
*
C           Select Two-Shock Riemann solver with
C           PVRS as estimate
*
            GEL = SQRT((G5/DL)/(G6*PL + PPV))
            GER = SQRT((G5/DR)/(G6*PR + PPV))
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
         ENDIF
      ENDIF
*
      END
*
*----------------------------------------------------------------------*
*
      SUBROUTINE PREFUN(F,FD,P,DK,PK,CK)
*
C     Purpose: to evaluate the pressure functions
C              FL and FR in exact Riemann solver
C              and their first derivatives
*
      IMPLICIT NONE
*
C     Declaration of variables
*
      REAL    AK, BK, CK, DK, F, FD, P, PK, PRATIO, QRT,
     & GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
*
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
*
      IF(P.LE.PK)THEN
*
C        Rarefaction wave
*
         PRATIO = P/PK
         F    = G4*CK*(PRATIO**G1 - 1.0)
         FD   = (1.0/(DK*CK))*PRATIO**(-G2)
      ELSE
*
C        Shock wave
*
         AK  = G5/DK
         BK  = G6*PK
         QRT = SQRT(AK/(BK + P))
         F   = (P - PK)*QRT
         FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
      ENDIF
*
      END
*
*----------------------------------------------------------------------*
*
      SUBROUTINE SAMPLE(PM, UM, S, D, U, P)
*
C     Purpose: to sample the solution throughout the wave
C              pattern. Pressure PM and velocity UM in the
C              Star Region are known. Sampling is performed
C              in terms of the 'speed' S = X/T. Sampled
C              values are D, U, P
*
C     Input variables : PM, UM, S, /GAMMAS/, /STATES/
C     Output variables: D, U, P
*
      IMPLICIT NONE
*
C     Declaration of variables
*
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,
     & GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,
     & C, CML, CMR, D, P, PM, PML, PMR,  S,
     & SHL, SHR, SL, SR, STL, STR, U, UM
*
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
*
      IF(S.LE.UM)THEN
*
C        Sampling point lies to the left of the contact
C        discontinuity
*
         IF(PM.LE.PL)THEN
*
C           Left rarefaction
*
            SHL = UL - CL
*
            IF(S.LE.SHL)THEN
*
C              Sampled point is left data state
*
               D = DL
               U = UL
               P = PL
            ELSE
               CML = CL*(PM/PL)**G1
               STL = UM - CML
*
               IF(S.GT.STL)THEN
*
C                 Sampled point is Star Left state
*
                  D = DL*(PM/PL)**(1.0/GAMMA)
                  U = UM
                  P = PM
               ELSE
*
C                 Sampled point is inside left fan
*
                  U = G5*(CL + G7*UL + S)
                  C = G5*(CL + G7*(UL - S))
                  D = DL*(C/CL)**G4
                  P = PL*(C/CL)**G3
               ENDIF
            ENDIF
         ELSE
*
C           Left shock
*
            PML = PM/PL
            SL  = UL - CL*SQRT(G2*PML + G1)
*
            IF(S.LE.SL)THEN
*
C              Sampled point is left data state
*
               D = DL
               U = UL
               P = PL
*
            ELSE
*
C              Sampled point is Star Left state
*
               D = DL*(PML + G6)/(PML*G6 + 1.0)
               U = UM
               P = PM
            ENDIF
         ENDIF
      ELSE
*
C        Sampling point lies to the right of the contact
C        discontinuity
*
         IF(PM.GT.PR)THEN
*
C           Right shock
*
            PMR = PM/PR
            SR  = UR + CR*SQRT(G2*PMR + G1)
*
            IF(S.GE.SR)THEN
*
C              Sampled point is right data state
*
               D = DR
               U = UR
               P = PR
            ELSE
*
C              Sampled point is Star Right state
*
               D = DR*(PMR + G6)/(PMR*G6 + 1.0)
               U = UM
               P = PM
            ENDIF
         ELSE
*
C           Right rarefaction
*
            SHR = UR + CR
*
            IF(S.GE.SHR)THEN
*
C              Sampled point is right data state
*
               D = DR
               U = UR
               P = PR
            ELSE
               CMR = CR*(PM/PR)**G1
               STR = UM + CMR
*
               IF(S.LE.STR)THEN
*
C                 Sampled point is Star Right state
*
                  D = DR*(PM/PR)**(1.0/GAMMA)
                  U = UM
                  P = PM
               ELSE
*
C                 Sampled point is inside left fan
*
                  U = G5*(-CR + G7*UR + S)
                  C = G5*(CR - G7*(UR - S))
                  D = DR*(C/CR)**G4
                  P = PR*(C/CR)**G3
               ENDIF
            ENDIF
         ENDIF
      ENDIF
*
      END
*
*----------------------------------------------------------------------*
*

