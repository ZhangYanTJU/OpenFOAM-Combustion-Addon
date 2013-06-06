      program LIMEX_Demo
c
c-----------------------------------------------------------------------
c
c     Main driver routine as an example for the integration
c     of DAE-systems with LIMEX 4.2A.
c
c     Small artificial example with an analytical Jacobian.
c
c-----------------------------------------------------------------------
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
      dimension  IFail ( 3 ), Iopt ( 30 ), IPos ( 3 ) 
c
      dimension  Ropt ( 5 ) 
c
      dimension  y ( 3 ), ys ( 3 ), aTol ( 3 ), rTol ( 3 ) 
c
      external   Fcn, Jacobian
c
c-----------------------------------------------------------------------
c
c     Size of system
c
c-----------------------------------------------------------------------
c
      n = 3
c
c-----------------------------------------------------------------------
c
c     Time interval 
c
c-----------------------------------------------------------------------
c
      t_Begin = 0.0d0
      t_End   = 0.1108d0
c
c-----------------------------------------------------------------------
c
c     Initial values and derivatives
c
c-----------------------------------------------------------------------
c
      y(1) = - 1.0d0
      y(2) =   5.0d0
      y(3) =   4.0d0 * y(1) / ( 3.0d0 * y(1) - 8.0d0 )
c
      ys(1) = 0.0d0
      ys(2) = 0.0d0
      ys(3) = 0.0d0
c
c-----------------------------------------------------------------------
c
c     Required accuracies of integration
c
c-----------------------------------------------------------------------
c
      rTol(1) = 1.0d-8
      aTol(1) = 1.0d-8
c
c-----------------------------------------------------------------------
c
c     Initial stepsize
c
c-----------------------------------------------------------------------
c
      h = 1.0d-3
c
c-----------------------------------------------------------------------
c
c     Control parameters 
c
c     Iopt(1)  Integration monitoring 
c              = 0 : no output
c              = 1 : standard output
c              = 2 : additional integration monitor
c
c     Iopt(2)  Unit number for monitor output
c
c     Iopt(3)  Solution output
c              = 0 : no output
c              = 1 : initial and solution values
c              = 2 : additional solution values at intermediate
c                    points
c
c     Iopt(4)  Unit number for solution output
c
c     Iopt(5)  Singular on nonsingular matrix B
c              = 0 : matrix B may be singular
c              = 1 : matrix B is nonsingular
c
c     Iopt(6)  Consistent initial value determination
c              = 0 : no determination of CIVs
c              = 1 : determination of CIVs
c
c     Iopt(7)  Numerical or analytical computation of the Jacobian
c              = 0 : Numerical approximation of the Jacobian
c              = 1 : Analytical computation of the Jacobian
c
c     Iopt(8)  Lower bandwith of the Jacobian if numerical
c              evaluated
c
c     Iopt(9)  Upper bandwith of the Jacobian if numerical
c              evaluated
c
c     Iopt(10) Control of reuse of the Jacobian
c              = 0 : no Jacobian reuse
c              = 1 : reuse of the Jacobian 
c
c     Iopt(11) Switch for error tolerances
c              = 0 : both rTol and aTol are scalars
c              = 1 : both rTol and aTol are vectors
c
c     Iopt(12) Sw1tch for one step mode
c              = 0 : one step mode off  
c              = 1 : one step mode on
c              = 2 : one step mode on, return on dense output
c                    points
c
c     Iopt(13) Dense output option
c              = 0 : no dense output
c              = 1 : dense output on equidistant points in the
c                    whole integration interval
c              = 2 : dense output on equidistant points in
c                    every integration interval
c              = 3 : dense output at additional points, that
c                    the distance between two dense output  
c                    points is at least Ropt(2)
c
c     Iopt(14) The number of equidistant output points if
c              Iopt(13) = 1 or 2
c
c     Iopt(15) Unit number for dense solution output
c
c     Iopt(16) Type of call
c              = 0 : first LIMEX call
c              = 1 : continuation call
c
c     Iopt(17) Switch for behavior of LIMEX on t_End
c              = 0 : integrate exactly up to t_End
c              = 1 : LIMEX may internally integrate to t > t_End
c
c     Iopt(18) Generation of of a PostScript plot of the Jacobian
c              = 0 : no PostScript plot produced
c              - 1 : in the initial step the file 'Jacobian.ps'
c                    will be produced with the initial Jacobian
c              > 0 : in step Iopt(18) the file 'Jacobian.ps'
c                    will be produced with the current Jacobian
c 
c     Iopt(19) - Iopt(23) 
c              not used in LIMEX versions 4.2A1 and 4.2A2
c
c     Ropt(1)  Maximum allowed stepsize
c
c     Ropt(2)  Maximum distance between two dense output points,
c              if Iopt(13) = 3
c
c     Ropt(3)  Upper bound for t, used only if Iopt(17) = 1
c
c     Ropt(4), Ropt(5) 
c              not used in LIMEX versions 4.2A1 and 4.2A2
c
c-----------------------------------------------------------------------
c
      Iopt(1)  = 2
      Iopt(2)  = 6
      Iopt(3)  = 1
      Iopt(4)  = 6
      Iopt(5)  = 0
      Iopt(6)  = 0
      Iopt(7)  = 0
      Iopt(8)  = 3
      Iopt(9)  = 3
      Iopt(10) = 0
      Iopt(11) = 0
      Iopt(12) = 0
      Iopt(13) = 0
      Iopt(14) = 0
      Iopt(15) = 0
      Iopt(16) = 0
      Iopt(17) = 0
      Iopt(18) = 0
c
      Ropt(1) = t_End - t_Begin
      Ropt(2) = 0.0005d0
      Ropt(3) = 0.1108d0
c
c-----------------------------------------------------------------------
c
c     Any component may be negative.
c
c-----------------------------------------------------------------------
c
      IPos(1) = 0
      IPos(2) = 0
      IPos(3) = 0
c
c-----------------------------------------------------------------------
c
c     Select call typ
c
c-----------------------------------------------------------------------
c
      write ( *, 9000 ) 
c
      read ( *, * ) iType
c
c-----------------------------------------------------------------------
c
c     Call of integration routine
c
c     a) simple call (one step mode off, no continuation)
c 
c-----------------------------------------------------------------------
c
      if ( iType .eq. 1 ) then
c
         call LIMEX ( n, Fcn, Jacobian, t_Begin, t_End, y, ys, 
     2                rTol, aTol, h, Iopt, Ropt, IPos, IFail )
c
c-----------------------------------------------------------------------
c
c     Call of integration routine
c
c     b) one step mode call with dense output
c 
c-----------------------------------------------------------------------
c
      else if ( iType .eq. 2 ) then
c
         IFail(1) = 0
c
         Iopt(12) = 2 
         Iopt(13) = 2
         Iopt(14) = 3
         Iopt(15) = 10
c
         do while ( IFail(1) .eq. 0 .and. t_Begin .lt. t_End ) 
c
            call LIMEX ( n, Fcn, Jacobian, t_Begin, t_End, y, ys, 
     2                   rTol, aTol, h, Iopt, Ropt, IPos, IFail )
c
         end do
c
c-----------------------------------------------------------------------
c
c     Call of integration routine
c
c     c) continuation calls (may be combined with one step mode) 
c
c-----------------------------------------------------------------------
c
      else
c
         IFail(1) = 0
c
         Iopt(17) = 1
c
         Delta_t = 0.1d0 * ( t_End - t_Begin )
         t       = t_Begin
c
         do i = 1, 10 
c
            t = t + Delta_t
c
            do while ( IFail(1) .eq. 0 .and. t_Begin .lt. t ) 
c
               call LIMEX ( n, Fcn, Jacobian, t_Begin, t, y, ys,
     2                      rTol, aTol, h, Iopt, Ropt, IPos, IFail )
c
            end do
c
            if ( IFail(1) .lt. 0 ) go to 110
c
         end do
c
  110    continue
c
      end if
c
c-----------------------------------------------------------------------
c
c     Write statistic of run
c
c-----------------------------------------------------------------------
c
      write ( *, 9010 ) ' Nr. of F-Evaluations (Int) : ', Iopt(24)
      write ( *, 9020 ) ' Nr. of F-Evaluations (Jac) : ', Iopt(25)
      write ( *, 9020 ) '             Nr. of LU-Dec. : ', Iopt(26)
      write ( *, 9020 ) '        Nr. of solver calls : ', Iopt(27)
      write ( *, 9020 ) '               Nr. of steps : ', Iopt(28)
      write ( *, 9020 ) '     Nr. of Jac. evaluation : ', Iopt(29)
c
c-----------------------------------------------------------------------
c
c     Format statements
c
c-----------------------------------------------------------------------
c
 9000 format ( ' Enter type of call : 1 => simple call', /, 
     2         27x, '(one step mode off, no continuation)', /, 22x,
     3         '2 => one step mode', /, 22x, '3 => continuation calls' )
 9010 format ( /, a, i8 )
 9020 format ( a, i8 )
c
c-----------------------------------------------------------------------
c
      end
      subroutine Fcn ( n, nz, t, y, f, b, ir, ic, Info )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     Example for the call of subroutine LIMEX
c     integrator of differential-algebraic systems. 
c
c     Defines the right-hand side f and the left-hand 
c     side matrix B(t,y).
c
c-----------------------------------------------------------------------
c
      dimension y ( n ), f ( n ), b ( * )
c
      dimension ir ( * ), ic ( * )
c
c-----------------------------------------------------------------------
c
c     Error indicator
c
c-----------------------------------------------------------------------
c
      Info = 0
c
c-----------------------------------------------------------------------
c
c     Right-hand side of the DAE
c
c-----------------------------------------------------------------------
c
      f(1) = - 0.8d0 * y(1) + 10.0d0 * y(2) - 0.6d0 * y(1) * y(3)
      f(2) = - 10.0d0 + 1.6d0 * y(3) / y(2)
      f(3) =   0.8d0 * y(1) + 1.6d0 * y(3) - 0.6d0 * y(1) * y(3)
c
c-----------------------------------------------------------------------
c
c     Left-hand side of the DAE and the non-zero entries of matrix B 
c     in coordinate format.
c
c-----------------------------------------------------------------------
c
      nz = 2
c
      ir(1) = 1
      ic(1) = 1
      b(1)  = 1.0d0
c
      ir(2) = 2
      ic(2) = 2
      b(2)  = 1.0d0 / y(2)
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Jacobian ( n, t, y, ys, Jac, LDJac, ml, mu, 
     2                      Full_or_Band, JacInfo )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     Example for the call of subroutine LIMEX
c     integrator of differential-algebraic systems. 
c
c     Defines the Jacobian of the residual of the DAE
c
c-----------------------------------------------------------------------
c
      integer          n, LDJac, ml, mu, Full_or_Band, JacInfo
c
      double precision Jac
c
      dimension        Jac ( LDJac, * ) 
c
      dimension        y ( * ), ys ( * )
c
c-----------------------------------------------------------------------
c
c     Define the Jacobian matrix.
c
c-----------------------------------------------------------------------
c
      if ( Full_or_Band .eq. 0 ) then
         k = 0
      else
         k = mu + 1 - 1
      end if
c
      Jac(1+k,1) = - 0.8d0 - 0.6d0 * y(3)
      Jac(2+k,1) = 0.0d0
      Jac(3+k,1) = 0.8d0 - 0.6d0 * y(3)
c
      if ( Full_or_Band .eq. 0 ) then
         k = 0
      else
         k = mu + 1 - 2
      end if
c
      Jac(1+k,2) = 10.0d0
      Jac(2+k,2) = ( ys(2) - 1.6d0 * y(3) ) / ( y(2) * y(2) )
      Jac(3+k,2) = 0.0d0
c
      if ( Full_or_Band .eq. 0 ) then
         k = 0
      else
         k = mu + 1 - 3
      end if
c
      Jac(1+k,3) = - 0.6d0 * y(1)
      Jac(2+k,3) = 1.6d0 / y(2)
      Jac(3+k,3) = 1.6d0 - 0.6d0 * y(1)
c
      JacInfo = 0
c
c-----------------------------------------------------------------------
c
      return
      end
