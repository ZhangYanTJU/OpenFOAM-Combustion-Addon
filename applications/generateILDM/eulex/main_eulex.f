C  Test Example for the Call of Subroutine EULEX
C
C  Explicit Extrapolation Integrator
C  for Non-Stiff Systems of Ordinary Differential Equations
C  (Based on the Explicit Euler Discretization)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(2)
      EXTERNAL FCN
C
C  Print-Parameter KFLAG
C  KFLAG=0 : No Output
C  KFLAG=1 : Integration Monitor
C  KFLAG=2 : Intermediate Solution Points T,Y(1),...
C  KFLAG=3 : Integration Monitor and Intermediate Solution Points
C
      KFLAG=1
C  Number of Equations
      N=2
C  Starting  Point of Integration
      T=0.D0
C  Initial Values Y(T)
      Y(1)=0.D0
      Y(2)=0.D0
C  Final Point of Integration
      TEND=20.D0
C  Desired Accuracy
      TOL=1.D-5
C  Maximum Permitted Stepsize
      HMAX=TEND-T
C  Initial Stepsize Guess (Here=TOL)
      H=TOL
C  Call of EULEX with these Parameters
      CALL EULEX (N,FCN,T,Y,TEND,TOL,HMAX,H,KFLAG)
C
      WRITE(6,101) T,(Y(I),I=1,N)
101   FORMAT('  SOLUTION AT T=  ',E25.16,/,10X,4E25.16)
      STOP
      END
C
C
C  Example E5 from Testset Detest (Due to Hull et al.)
      SUBROUTINE FCN(N,T,Y,DY)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y(2),DY(2)
C
      DY(1)=Y(2)
      DY(2)=DSQRT(1.D0+Y(2)*Y(2))/(25.D0-T)
C
      RETURN
      END
