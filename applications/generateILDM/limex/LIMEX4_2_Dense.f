c-----------------------------------------------------------------------
c
c                           LIMEX4_2_Dense
c
c-----------------------------------------------------------------------
c
c     This is a collection of subroutines for LIMEX, the extrapolation
c     integrator for  the solution  of linearly-implicit differential-
c     algebraic systems of the form
c
c          B (t,y) * y' (t) = f (t,y)
c
c     with B a (n,n)-matrix of rank less or equal n.
c
c     This  collection contains  routines for  the computation and the
c     evaluation of right derivatives and Hermite polynomials.
c
c-----------------------------------------------------------------------
c
c     Copyright (C) 2000, Konrad-Zuse-Zentrum fuer Informationstechnik
c     Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c            **************************************************
c            **                                              **
c            **    This is version 4.2 of March, 17, 2000    **
c            **                                              **
c            **************************************************
c
c-----------------------------------------------------------------------
c
c     Overview of current versions:
c
c     4.2A1  Non-sparse dense  or  banded  Jacobians. Direct  solvers:
c            LAPACK routines. Based on BLAS and LAPACK routines.
c
c     4.2A2  Non-sparse dense or banded Jacobians. Direct solvers: NAG
c            routines. Based on NAG routines (Mark 16).
c
c     4.2B1  Sparse Jacobians. Direct  solvers: MA28 routines from the
c            Harwell subroutine library, iterative  solvers: GMRES and
c            BICGSTAB  with a  variable ILU  preconditioner. Based  on
c            BLAS and LAPACK routines.
c
c     4.2B2  Sparse Jacobians. Direct solvers: NAG routines, iterative
c            solvers: GMRES  and  BICGSTAB  with a  variable ILU  pre-
c            conditioner. Based on NAG routines (Mark 16).
c
c     Versions with other solver routines are available on request.
c
c-----------------------------------------------------------------------
c
c     NOTICE: "The LIMEX  program may be used SOLELY  for educational,
c     research, and benchmarking  purposes by no-profit organizations.
c     Commercial and other organizations  may make use of LIMEX SOLELY
c     for benchmarking  purposes only. LIMEX may  be modified by or on
c     behalf of the user for  such use but  at no time  shall LIMEX or
c     any such  modified version  of LIMEX become the  property of the
c     user. LIMEX  is provided  without warranty  of any  kind, either
c     expressed  or implied.  Neither the  authors nor their employers
c     shall be liable  for any direct or consequential  loss or damage
c     whatsoever  arising out  of the  use or  misuse of  LIMEX by the
c     user. LIMEX  must not  be sold. You  may make  copies  of LIMEX,
c     but  this NOTICE  and the  Copyright notice  must appear  in all
c     copies. Any  other  use  of LIMEX  requires written  permission.
c     Your use of LIMEX is an implicit agreement to these conditions."
c
c     LIMEX is even  available for commercial use. Contact the authors
c     who will provide details of price and condition of use.
c
c-----------------------------------------------------------------------
c
c     Contents:
c
c     Comp_Herm : computes right derivatives and the coefficient of an
c                 Hermite interpolation polynomial.
c
c     Eval_Herm : evaluates an Hermite polynomial.
c
c-----------------------------------------------------------------------
c
      subroutine Comp_Herm ( n, Max_Size, Dense, k, ipt, nj, Work )
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     This  subroutine computes  approximations  for right derivatives
c     the coefficients of a Hermite interpolation polynomial.
c
c-----------------------------------------------------------------------
c
c     Copyright (C) 2000, Konrad-Zuse-Zentrum fuer Informationstechnik
c     Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c     Arguments:
c
c-----------------------------------------------------------------------
c
c     n          Integer variable, must be set by caller on input, not
c                modified. Size of the system.
c
c     Max_Size   Integer variable, must be set by caller on input, not
c                modified. Leading  dimension of  the  array  Dense as 
c                declared in LIMEX.
c
c     Dense      Real array of size at least Max_Size * ipt(k+1), must 
c                be set by caller on input. The values of intermediate 
c                solutions resp. their differences  as computed in the
c                current integration step. On exit the coefficients of
c                the Hermite interpolation poynomial.
c
c     k          Integer variable, must be set by caller on input. The 
c                currrent order of the extrapolation scheme.
c
c     ipt        Integer  array of size at least k + 1, must be set by  
c                caller on  input. Pointer array  to the  intermediate 
c                solutions. 
c
c     nj         Integer array  of size at least k + 1, must be set by      
c                caller on input. The stepsize dividing sequence.
c
c     Work       Real  array of  size at  least n * (k+1), need not be
c                ste by caller on input. Work array.
c
c-----------------------------------------------------------------------
c
c     Subroutines and functions called
c
c-----------------------------------------------------------------------
c
c     From the BLAS library:
c     ----------------------
c
c     daxpy      Adds a scalar multiple of a vector to another vector
c
c     dcopy      Copies real vectors
c
c-----------------------------------------------------------------------
c
c     Declaration of the arguments
c
c-----------------------------------------------------------------------
c
      integer            n, Max_Size, k, ipt ( * ), nj ( * )
c
      double precision   Dense ( Max_Size, * ), Work ( n, * )
c
c-----------------------------------------------------------------------
c
c     Local scalars
c
c-----------------------------------------------------------------------
c
      integer            i, ifac, j, j1, j2
c
      double precision   fac, tmp
c
c-----------------------------------------------------------------------
c
c     Internal parameters
c
c-----------------------------------------------------------------------
c
      double precision   one
c
      parameter ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
c     Main loop until the current order k.
c
c-----------------------------------------------------------------------
c
      ifac = 1
c
      do j = 1, k
c
         ifac = ifac * j
         fac  = one / dble ( ifac )
c
c-----------------------------------------------------------------------
c
c     Store the  normalized right  differences in the  first column of
c     the extrapolation tableau.
c
c-----------------------------------------------------------------------
c
         do j1 = j, k + 1 
c
            tmp = fac * dble ( nj(j1) ** j )
            j2  = ipt(j1)
c
            do i = 1, n
               Work(i,j1) = tmp * Dense(i,j2)
            end do
c
         end do
c
c-----------------------------------------------------------------------
c
c     Extrapolation of the right differences.
c
c-----------------------------------------------------------------------
c
         do j1 = j, k
c
            do j2 = j1, j, - 1
c
               tmp = one / ( dble ( nj(j1+1) ) / dble ( nj(j2) ) - one )
c
               do i = 1, n
                  Work(i,j2) =   Work(i,j2+1)
     2                         + tmp * ( Work(i,j2+1) - Work(i,j2) )
               end do   
c
            end do  
c
         end do
c
c-----------------------------------------------------------------------
c
c     Compute higher order differences.
c
c-----------------------------------------------------------------------
c
         if ( j .lt. k ) then
c
            do j1 = j, k
               do i = ipt(j1+1), ipt(j1+1) - j1 + j, - 1
                  call daxpy ( n, -one, Dense(1,i-1), 1, Dense(1,i), 
     2                         1 ) 
               end do
            end do
c
         end if
c
c-----------------------------------------------------------------------
c
c     Restore estimates for the normalized derivatives.
c
c-----------------------------------------------------------------------
c
         call dcopy ( n, Work(1,j), 1, Dense(1,j+2), 1 )
c
      end do 
c
c-----------------------------------------------------------------------
c
c     Compute coefficients of the Hermite interpolation polynomial.
c
c-----------------------------------------------------------------------
c
      do j = 1, k + 1
         call daxpy ( n, -one, Dense(1,j), 1, Dense(1,j+1), 1 )
      end do   
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Eval_Herm ( n, Max_Size, Dense, k, tFac, y_Interp )
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     This subroutine evaluates a Hermite interpolation polynomial.
c
c-----------------------------------------------------------------------
c
c     Copyright (C) 2000, Konrad-Zuse-Zentrum fuer Informationstechnik
c     Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c     Arguments:
c
c-----------------------------------------------------------------------
c
c     n          Integer variable, must be set by caller on input, not
c                modified. Size of the system.
c
c     Max_Size   Integer variable, must be set by caller on input, not
c                modified. Leading  dimension of  the  array  Dense as
c                declared in LIMEX.
c
c     Dense      Real array  of size at least  Max_Size * k+1, must be
c                set by caller on input. The coefficients of a Hermite
c                interpolation poynomial.
c
c     k          Integer variable, must be set by caller on input. The
c                currrent order of the extrapolation scheme.
c
c     tFac       Real  variable, must  be set  by caller on input. The
c                ratio between t - t0 and the actual stepsize, if t is
c                the current dense output point within this step.
c
c     y_Interp   Real  array of  size n, need not be  set by caller on
c                input. On  exit the estimated  solution values at the
c                dense output  point, i.e. the values  of the  Hermite 
c                interpolation polynomials.
c
c-----------------------------------------------------------------------
c
c     Subroutines and functions called
c
c-----------------------------------------------------------------------
c
c     From the BLAS library:
c     ----------------------
c
c     dcopy      Copies real vectors
c
c-----------------------------------------------------------------------
c
c     Declaration of the arguments
c
c-----------------------------------------------------------------------
c
      integer            n, Max_Size, k
c
      double precision   Dense ( Max_Size, * ), tFac, y_Interp ( * )
c
c-----------------------------------------------------------------------
c
c     Local scalars
c
c-----------------------------------------------------------------------
c
      integer            i, j
c
      double precision   tmp
c
c-----------------------------------------------------------------------
c
c     Internal parameters
c
c-----------------------------------------------------------------------
c
      double precision   one
c
      parameter ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
c     Evaluate  the Hermite  interpolation polynomials  by the  Horner
c     scheme.
c
c-----------------------------------------------------------------------
c
      tmp = tFac - one
c
      call dcopy ( n, Dense(1,k+2), 1, y_Interp, 1 )
c
      do j = 1, k
         do i = 1, n
            y_Interp(i) = Dense(i,k+2-j) + tmp * y_Interp(i)
         end do
      end do
c
      do i = 1, n
         y_Interp(i) = Dense(i,1) + tFac * y_Interp(i)
      end do
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine JacPlot ( n, LenJac, Jac, LDJac, ml, mu, ir, ic, nStep,
     2                     mode )
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     This subroutine generates a Postscript plot of the Jacobian.
c
c-----------------------------------------------------------------------
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c     Arguments:
c
c-----------------------------------------------------------------------
c
c     n          Integer variable, must be set by caller on input, not
c                modified. Size of the system.
c
c     LenJac     Integer variable, must be set by caller on input, not
c                modified. Number of  non-zero entries in the Jacobian
c                if specified by parameter mode as a sparse matrix.  
c
c     Jac        Real array of size (LDJac,*), must  be set by  caller
c                on  input only  if  mode = 0 or 1, not  modified. The 
c                values of the Jacobian matrix in full or banded form. 
c
c     LDJac      Integer variable, must be set by caller on input only
c                if mode = 0 or 1, not  modified. Leading dimension of 
c                array Jac.
c
c     ml         Integer variable, must be set by caller on input only
c                if mode = 1, not  modified. For  banded Jacobians the
c                number of lower diagonals.
c
c     mu         Integer variable, must be set by caller on input only
c                if mode = 1, not  modified. For  banded Jacobians the
c                number of upper diagonals.
c
c     ir         Integer array, must be set by caller on input only if
c                mode = 2, not modified. Row  indices for the non-zero
c                entries of the Jacobian in sparse form.
c
c     ic         Integer array, must be set by caller on input only if
c                mode = 2, not  modified. Column indices for the  non-
c                zero entries of the Jacobian is sparse form.
c
c     nStep      Integer variable, must be set by caller on input, not
c                modified. Current step number, which  is displayed in
c                the title.
c
c     mode       Integer variable, must be set by caller on input, not
c                modified. Determines the Jacobian matrix format.     
c                = 0 : dense Jacobian
c                = 1 : dense banded Jacobian
c                = 2 : sparse Jacobian in coordinate format
c
c-----------------------------------------------------------------------
c
c     Subroutines and functions called:
c
c     none
c
c-----------------------------------------------------------------------
c
c     Declaration of the arguments
c
c-----------------------------------------------------------------------
c
      integer            n, LenJac, LDJac, ml, mu, ir ( * ), ic ( * ),
     2                   nStep, mode
c
      double precision   Jac ( LDJac, * )
c
c-----------------------------------------------------------------------
c
c     Local scalars
c
c-----------------------------------------------------------------------
c
      integer            i, j
c
      real               Bottom, cm_dot, Left, Margin, Right, Scaling, 
     2                   Top, xTitle, yTitle 
c
c-----------------------------------------------------------------------
c
c     Internal parameters, sizes may be changed.
c
c-----------------------------------------------------------------------
c
      real               PaperSize, PlotSize
c
      parameter        ( PaperSize = 21.0, PlotSize = 15.0 )
c
c-----------------------------------------------------------------------
c
c     Define conversion  factor from  cm to  dot, margin  and scaling 
c     factor.
c 
c-----------------------------------------------------------------------
c
      cm_dot  = 72.0 / 2.54
      Margin  = 0.5 * ( PaperSize - PlotSize )
      Scaling = 0.125 * PlotSize * cm_dot / float ( n + 1 ) + 10.0
c
c-----------------------------------------------------------------------
c
c     Define title coordinates and bounding box.
c
c-----------------------------------------------------------------------
c
      xTitle = 0.5 * PaperSize
      yTitle = 3.0 + PlotSize
c
      Left   = Margin * cm_dot - Scaling
      Right  = ( Margin + PlotSize ) * cm_dot + Scaling
      Bottom = 2.0 * cm_dot - Scaling
      Top    = ( 3.35 + PlotSize ) * cm_dot + Scaling
c
c-----------------------------------------------------------------------
c
c     Open file and write postscript header information.
c
c-----------------------------------------------------------------------
c
      open ( 99, file = 'Jacobian.ps' )
c
      write ( 99, 9000 ) '%!'
      write ( 99, 9000 ) '%%Creator: JacPlot in LIMEX4.2'
      write ( 99, 9010 ) '%%BoundingBox:', Left, Bottom, Right, Top
      write ( 99, 9000 ) '%%EndComments'
      write ( 99, 9000 ) '/cm {72 mul 2.54 div} def'
      write ( 99, 9000 ) '/mc {72 div 2.54 mul} def'
      write ( 99, 9000 ) '/pnum { 72 div 2.54 mul 20 string'
      write ( 99, 9000 ) 'cvs print ( ) print} def'
      write ( 99, 9000 )
     2   '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
      write ( 99, 9000 ) 'gsave'
      write ( 99, 9000 ) 
     2   '/Helvetica findfont 0.5 cm scalefont setfont'
      write ( 99, 9020 ) xTitle, ' cm ', yTitle, ' cm moveto'
      write ( 99, 9030 ) '(Jacobian matrix from LIMEX4.2, step = ', 
     2   nStep, ') Cshow'
      write ( 99, 9040 ) Margin, 'cm 2.0 cm translate'
      write ( 99, 9050 ) PlotSize, 'cm', n + 1, 'div dup scale'
c
c-----------------------------------------------------------------------
c
c     Draw frame.
c
c-----------------------------------------------------------------------
c
      write ( 99, 9000 ) '0.25 setlinewidth'
      write ( 99, 9000 ) 'newpath'
      write ( 99, 9000 ) '0 0 moveto'
      write ( 99, 9060 ) n + 1, 0, 'lineto'
      write ( 99, 9060 ) n + 1, n + 1,'lineto'
      write ( 99, 9060 ) 0, n + 1, 'lineto'
      write ( 99, 9000 ) 'closepath stroke'
c
c-----------------------------------------------------------------------
c
c     Prepare matrix plot.
c
c-----------------------------------------------------------------------
c
      write ( 99, 9000 ) '1 1 translate'
      write ( 99, 9000 ) '0.8 setlinewidth'
      write ( 99, 9000 ) '/p {moveto 0 -.40 rmoveto'
      write ( 99, 9000 ) '           0  .80 rlineto stroke} def'
c     
c-----------------------------------------------------------------------
c
c     If mode = 0 the matrix is dense and full.
c
c-----------------------------------------------------------------------
c
      if ( mode .eq. 0 ) then
c
         do i = 1, n
            do j = 1, n
               if ( Jac(j,i) .ne. 0.0d0 ) 
     2            write ( 99, 9060 ) i - 1, n - j, 'p'
            end do
         end do
c
c-----------------------------------------------------------------------
c
c     If mode = 1 the matrix is dense and banded.
c
c-----------------------------------------------------------------------
c
      else if ( mode .eq. 1 ) then
c
         do i = 1, n
            do j = max ( 1, i - mu ), min ( n, i + ml )
               if ( Jac(mu+1+j-i,i) .ne. 0.0d0 )
     2            write ( 99, 9060 ) i - 1, n - j, 'p'
            end do
         end do
c
c-----------------------------------------------------------------------
c
c     If mode = 2 the matrix is sparse.
c
c-----------------------------------------------------------------------
c
      else if ( mode .eq. 2 ) then
c
         do i = 1, LenJac
            write ( 99, 9060 ) ic(i) - 1, n - ir(i), 'p'
         end do
c
      end if
c
c-----------------------------------------------------------------------
c
c     Finalize and close postscript file.
c
c-----------------------------------------------------------------------
c
      write ( 99, 9000 ) 'showpage'
c
      close ( 99 )
c
c-----------------------------------------------------------------------
c
c     Format-statements
c
c-----------------------------------------------------------------------
c
 9000 format ( a )
 9010 format ( a, 4 ( 1x, f9.2 ) )
 9020 format ( 2 ( 1x, f9.2, a ) )
 9030 format ( a, i6, a )
 9040 format ( f9.2, 1x, a )
 9050 format ( f9.2, 1x, a, 1x, i6, 1x, a )
 9060 format ( 2 ( i6, 1x ), a )
c
c-----------------------------------------------------------------------
c
      return
      end
