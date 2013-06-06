      Program LIMEX_Demo
c
c-----------------------------------------------------------------------
c
c     Main driver routine as an example for the integration
c     of DAE-systems with LIMEX version 4.2A.
c
c     This example models the startup phase of an automobile 
c     catalytic converter. May be used solely for benchmarking
c     and test purposes.
c
c     This example computes consistent initial values and
c     evaluates the (banded) Jacobian numerically.
c
c-----------------------------------------------------------------------
c
      implicit double precision ( a-h, o-z ) 
c
c-----------------------------------------------------------------------
c
c     Parameters
c
c     Max_Nr_of_PDEs      Maximal number of PDEs
c
c     Max_Gridsize        Maximal number of gridpoints
c
c     Max_Size            Maximal system size
c
c     Nr_Int_Par          Integer model parameters
c
c     Nr_Real_Par         Real model parameters
c
c-----------------------------------------------------------------------
c
c     A similar parameter statement occurs in most subroutines.
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      dimension IFail ( 3 ), Iopt ( 30 ), IPos ( Max_Size )
c
      dimension Ropt  ( 5 )
c
      dimension y     ( Max_Size ), ys   ( Max_Size ),
     2          aTol  ( Max_Size ), rTol ( Max_Size )
c
      external  Fcn, Jacobian
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPAr ( Nr_Int_Par ), rPAr ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
c     Define problem 
c
c-----------------------------------------------------------------------
c
      call Define_Problem ( )
c
c-----------------------------------------------------------------------
c
c     Define grid
c
c-----------------------------------------------------------------------
c 
      call Define_Grid ( )
c
c-----------------------------------------------------------------------
c
c     Define initial values
c
c-----------------------------------------------------------------------
c
      call Define_Initial_Values ( y, ys ) 
c
c-----------------------------------------------------------------------
c
c     Define integration specifications
c
c-----------------------------------------------------------------------
c
c     Size of system
c
c-----------------------------------------------------------------------
c
      n = Nr_of_PDEs * Gridsize
c
c-----------------------------------------------------------------------
c
c     Time interval
c
c-----------------------------------------------------------------------
c
      t_Begin =  0.0d0
      t_End   = 28.0d0
c
c-----------------------------------------------------------------------
c
c     Required accuracies of integration
c
c-----------------------------------------------------------------------
c
      rTol(1) = 1.0d-3
      aTol(1) = 1.0d-3
c
c-----------------------------------------------------------------------
c
c     Initial stepsize
c
c-----------------------------------------------------------------------
c
      h = 1.0d-5
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
      Iopt(4)  = 10
      Iopt(5)  = 0
      Iopt(6)  = 1
      Iopt(7)  = 0
      Iopt(8)  = 11 
      Iopt(9)  = 11
      Iopt(10) = 1
      Iopt(11) = 0
      Iopt(12) = 0
      Iopt(13) = 0
      Iopt(14) = 0
      Iopt(15) = 0
      Iopt(16) = 0
      Iopt(17) = 0
      Iopt(18) = 12
c
      Ropt(1) = t_End - t_Begin
      Ropt(2) = 0.0d0
      Ropt(3) = 0.0d0
c
c-----------------------------------------------------------------------
c
c     Any component may be negative.
c
c-----------------------------------------------------------------------
c
      do i = 1, n
         IPos(i) = 0
      end do
c
c-----------------------------------------------------------------------
c
c     Call of integration routine
c
c     simple call (one step mode off, no continuation)
c
c-----------------------------------------------------------------------
c
      call LIMEX ( n, Fcn, Jacobian, t_Begin, t_End, y, ys,
     2             rTol, aTol, h, Iopt, Ropt, IPos, IFail )
c
c-----------------------------------------------------------------------
c
      write ( *, '(/,a,i8)' ) ' Nr. of F-Evaluations (Int) : ', Iopt(24)
      write ( *, '(a,i8)' )   ' Nr. of F-Evaluations (Jac) : ', Iopt(25)
      write ( *, '(a,i8)' )   '             Nr. of LU-Dec. : ', Iopt(26)
      write ( *, '(a,i8)' )   '        Nr. of solver calls : ', Iopt(27)
      write ( *, '(a,i8)' )   '               Nr. of steps : ', Iopt(28)
      write ( *, '(a,i8)' )   '     Nr. of Jac. evaluation : ', Iopt(29)
c
c-----------------------------------------------------------------------
c
      end
      subroutine Define_Problem ( )
C
      implicit double precision ( a-h, o-z )
c 
c-----------------------------------------------------------------------
c
c     This subroutine is an user routine for the basic definition 
c     of the problem. 
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c                                                     
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation   
c     Gridsize     size of grid for space discretization     
c     x_Left       left bound of domain                       
c     x_Right      right bound of domain                     
c     iPar         array for real user parameters          
c     rPar         array for integer user parameters          
c                                                             
c-----------------------------------------------------------------------
c 
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar 
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
c     Number of PDEs
c
c-----------------------------------------------------------------------
c
      Nr_of_PDEs = 11
c
c-----------------------------------------------------------------------
c
c     Number of grid points
c
c-----------------------------------------------------------------------
c
      Gridsize = 321
c
c-----------------------------------------------------------------------
c
c     Domain
c
c-----------------------------------------------------------------------
c
      x_Left  = 0.0
      x_Right = 0.16 
c
c-----------------------------------------------------------------------
c
c     User model parameters
c
c-----------------------------------------------------------------------
c
      ipar(1)  = 1
c
      rpar(1)  =     60.0
      rpar(2)  =     62.0
      rpar(3)  =      0.0
      rpar(4)  =      0.0
      rpar(5)  =      0.0
      rpar(6)  =      0.656
      rpar(7)  =   2622.0
      rpar(8)  =  55050.0
      rpar(9)  =      0.0
      rpar(10) =      0.0
      rpar(11) =      0.0
      rpar(12) =      1.0E-03
      rpar(13) =      1.08966
      rpar(14) =      0.0
      rpar(15) =      0.0
      rpar(16) =      0.0
      rpar(17) =      1.1
      rpar(18) =    659.0
      rpar(19) =      8.5E-04
      rpar(20) =      0.0
      rpar(21) =      0.0
      rpar(22) =      0.0
      rpar(23) =      0.45
      rpar(24) =   7150.0
      rpar(25) =      1.3E-02
      rpar(26) =      1.0
      rpar(27) =      1.0
      rpar(28) =      0.0
      rpar(29) =     42.0
      rpar(30) =      1.392E+12
      rpar(31) =  14556.0
      rpar(32) =   2080.0
      rpar(33) =    361.0
      rpar(34) =      1.928E+06
      rpar(35) =      0.0
      rpar(36) =      0.0
      rpar(37) =      0.0
      rpar(38) =     28.0
      rpar(39) =      6.699E+10
      rpar(40) =  12556.0
      rpar(41) =     65.5
      rpar(42) =    961.0
      rpar(43) = 283200.0
      rpar(44) =      0.0
      rpar(45) =      0.0
      rpar(46) =      0.0
      rpar(47) =      2.0
      rpar(48) = 242000.0
      rpar(49) =      0.0
      rpar(50) =      3.98
      rpar(51) =  11611.0
      rpar(52) =     32.0
      rpar(53) =      0.0
      rpar(54) =      0.0
      rpar(55) =      0.0
      rpar(56) =      0.156
      rpar(57) =      0.2
      rpar(58) =      0.9
      rpar(59) =      0.265
      rpar(60) =      0.0
      rpar(61) =      0.0
      rpar(62) =      0.0
      rpar(63) =      5.0E-05
      rpar(64) =      1.0E-04
      rpar(65) =      3.0E-04
      rpar(66) =      1.0E-04
      rpar(67) =      0.0
      rpar(68) =      0.0
      rpar(69) =      0.0
      rpar(70) =      1.0
      rpar(71) =     20.0
      rpar(72) =      0.0
      rpar(73) =      0.0
      rpar(74) =      0.0
      rpar(75) =      5.0E-03
      rpar(76) =      3.66
      rpar(77) =      0.0
      rpar(78) =      0.0
      rpar(79) =      0.0
      rpar(80) =      1.0E-03
      rpar(81) =      0.0
      rpar(82) =      0.0
      rpar(83) =      0.0
      rpar(84) =      2.0E-02
      rpar(85) =      5.0E-02
      rpar(86) =      1.0E-02
      rpar(87) =      1.0E-03
      rpar(88) =      1.0E-01
      rpar(89) =      0.0
      rpar(90) =      0.0
      rpar(91) =      0.0
      rpar(92) =     29.0
      rpar(93) =      0.0
      rpar(94) =      0.0
      rpar(95) =    100.0
      rpar(96) =      1.0E-04
      rpar(97) =      6.7E-03
      rpar(98) =      3.5E-04
      rpar(99) =     20.0
      rpar(100)=    400.0
c
      return
      end
      subroutine Define_Initial_Values ( Solution, Derivatives )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine is an user routine for the local definition 
c     of initial values u(x,t0).                                 
c                                                            
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize    size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension  Solution ( * ), Derivatives ( * )
c
c-----------------------------------------------------------------------
c
c     BeginOfInitialValues
c
c     KAT-1
c
c-----------------------------------------------------------------------
c
      do i = 1, Gridsize
c
         j = ( i - 1 ) * Nr_of_PDEs 
c
         Solution(j+1)  = 20.0d0 
         Solution(j+2)  = 20.0d0 
         Solution(j+3)  = 20.0d0 
         Solution(j+4)  =  0.0d0 
         Solution(j+5)  =  0.0d0 
         Solution(j+6)  =  0.0d0 
         Solution(j+7)  =  0.0d0 
         Solution(j+8)  =  0.0d0 
         Solution(j+9)  =  0.0d0 
         Solution(j+10) =  0.2d0 
         Solution(j+11) =  0.2d0 
c
      end do
c
      do i = 1, Nr_of_PDEs * Gridsize
         Derivatives(i) = 0.0d0
      end do
c
c-----------------------------------------------------------------------
c 
c EndOfInitialValues
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Diff_Coeff ( ix, x, t, u, dval )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine is an user routine for the local definition   
c     of the problem describing diffusion coefficient matrix dval    
c                                                               
c-----------------------------------------------------------------------
c
c     Input parameter 
c
c     ix           actual position in discretization process 
c     x            actual position in space (=xgrid(ix+1/2) 
c     t            actual time                             
c     u            array of length npde, containing actual solution 
c                  for actual (x,t)-values (linear interpolation   
c                    of u(x(ix),t) and u(x(ix+1),t) )             
c                                                          
c
c     Output parameter                                       
c                                                       
c     dval         actual values for diffusion coefficients as       
c                  square matrix dval(npde,npde) 
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension    u ( * ), dval ( Max_Nr_of_PDEs, * )
c
c-----------------------------------------------------------------------
c
c     Initialization
c
      npde = Nr_of_PDEs
c
      do k = 1, npde
         do j = 1, npde
            dval(j,k) = 0.d0
         end do
      end do
c
c-----------------------------------------------------------------------
c
c     Begin of diffusion coefficients
c
c-----------------------------------------------------------------------
c
      r_i = rpar(1)
      r_a = rpar(2)
c
      eps_w = (r_a*r_a - r_i*r_i)/r_i/r_i
      eps_g = rpar(6)
      eps_k = 1. - eps_g
c
      xlam_g = rpar(12)
      xlam_k = rpar(19)
      xlam_w = rpar(25)
c
      d_eff_co   = rpar(64)
      d_eff_c3h6 = rpar(63)
      d_eff_h2   = rpar(65)
      d_eff_o2   = rpar(66)
c
c-----------------------------------------------------------------------
c
c     Diffusion coefficients fluid  nd boundary 
c
c-----------------------------------------------------------------------
c
      dval(1,1) = eps_w*xlam_w
      dval(2,2) = eps_g*xlam_g
      dval(3,3) = eps_k*xlam_k
      dval(5,5) = eps_g*D_eff_co
      dval(4,4) = dval(5,5)
      dval(7,7) = eps_g*D_eff_c3h6
      dval(6,6) = dval(7,7)
      dval(9,9) = eps_g*D_eff_h2
      dval(8,8) = dval(9,9)
      dval(11,11) = eps_g*D_eff_o2
      dval(10,10) = dval(11,11)
c
c-----------------------------------------------------------------------
c
c     End of diffusion coefficients
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine PDE ( ix, x, t, u, ux, duxx, ffun, bmat )
c 
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine is an user routine for the local definition 
c     of the right hand side of the PDE.                         
c                                                               
c-----------------------------------------------------------------------
c                                              
c     Input parameter
c                                              
c     ix           actual position in actual grid           
c                  (number of actual node)                 
c     x            actual position in space               
c     t            actual time                           
c     u            array containing actual solution
c                  for actual (x,t)-values                        
c     ux           array containing actual space 
c                  derivative of solution for actual (x,t)-values
c     duxx         matrix containing actual diffusion           
c                  terms for actual (x,t)-values               
c                                                       
c     Output parameter                                    
c                                                    
c     ffun         actual rhs of PDE                
c     bmat         actual lhs of the DAE
c                                                  
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right,
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension    u ( * ), ux ( * ), duxx ( Max_Nr_of_PDEs, * )
c
      dimension    ffun ( * ), bmat ( Max_Nr_of_PDEs, * ) 
c
c-----------------------------------------------------------------------
c
      dimension    conv ( Max_Nr_of_PDEs ), q ( Max_Nr_of_PDEs )
c
      save
c
      data         icount / 0 /
c
c-----------------------------------------------------------------------
c
c     Begin of equations
c
c-----------------------------------------------------------------------
c
      if ( icount .eq. 0 ) then
c
         r_a = rpar(2)/1000.
         r_i = rpar(1)/1000.
c
         eps_g = rpar(6)
         eps_k = 1. - eps_g
         eps_w = 0.1
c
         cp_w = rpar(23)
         cp_k = rpar(17)
         cp_g = rpar(13)
c
         rho_w = rpar(24)
         rho_k = rpar(18)
c
         pges = rpar(70)*1.0d5
         rr   = 8314.
         xm_luft = rpar(92)
c
      end if
c
      TI    = u(2) + 273.15
      RHO_G = PGES * XM_LUFT/(RR * TI)
c
      bmat(1,1) = eps_w * rho_w * cp_w
      bmat(2,2) = eps_g * rho_g * cp_g 
      bmat(3,3) = eps_k * rho_k * cp_k
      bmat(5,5) = eps_g * rho_g
      bmat(4,4) = bmat(5,5)
      bmat(7,7) = eps_g * rho_g
      bmat(6,6) = bmat(7,7)
      bmat(9,9) = eps_g * rho_g
      bmat(8,8) = bmat(9,9)
      bmat(11,11) = eps_g * rho_g
      bmat(10,10) =  bmat(11,11)
c        
c-----------------------------------------------------------------------
c
      if( icount .eq. 0 ) then
c
         r_i  = rpar(1)/1000.
         pi   = 3.1415927
         CP_G = rpar(13) 
         xmp  = rpar(84)
         GZ   = xmp/pi/r_i/r_i
c
      end if
c
      conv(1) = 0.0
      conv(2) = gz*cp_g*ux(2) 
      conv(3) = 0.0
      conv(4) = 0.0 
      conv(5) = gz *ux(5)
      conv(6) = 0.0 
      conv(7) = gz *ux(7)
      conv(8) = 0.0 
      conv(9) = gz *ux(9)
      conv(10) = 0.0 
      conv(11) = gz *ux(11)
c
c-----------------------------------------------------------------------
c
      if( icount .eq. 0 ) then
c
         pi   = 3.1415927
         pges = rpar(70)*1.0d5
c
         xm_luft = rpar(92)
         xm_co   =  rpar(38)
         xm_c3h6 = rpar(29)
         xm_h2   = rpar(47)
         xm_o2   = rpar(52)
c
         rr = 8314.0 
c
         xkz_co = rpar(39)
         xkn_co = rpar(41)
         ez_co  = rpar(40)
         en_co  = rpar(42)
c
         xkz_c3h6 = rpar(30)
         xkn_c3h6 = rpar(32)
         ez_c3h6  =  rpar(31)
         en_c3h6  = rpar(33)
c
         xk_coc3h6 = rpar(50)
         e_coc3h6 =  rpar(51)
c
         xhr_co   = rpar(43)
         xhr_c3h6 = rpar(34)
         xhr_h2   = rpar(48)
c
         beta_co   = rpar(57)
         beta_c3h6 = rpar(56)
         beta_h2   = rpar(58)
         beta_o2   = rpar(59)
c
         av = rpar(7)
         ax = rpar(8)
c
         r_i = rpar(1)/1000.
         r_a = rpar(2)/1000.
c
         alpha_a = rpar(75)
         tu      = rpar(71)
         xnu     = rpar(76)
         xlam_g  = rpar(12)
         d_h     = rpar(80)
         xlam_w  = rpar(25)
c
         alpha_i = xnu*xlam_g/d_h/20.
c
         D_BM    = rpar(97)
         XLAM_BM = rpar(96)
         A_H     = pi * r_i * r_i 
         XL      = 0.16d0 
         R_B     = 2./3. * r_i
         A_B     = PI * r_b * r_b
         XK_S    = RPAR(98)
c
         HILFA = log(r_i/r_b)       /(2.*PI*XK_S   *XL)
     2     + log((r_i + d_bm)/r_i) / (2.*PI*XLAM_BM*XL)
     3     + log((r_a + r_i)/2/r_i) /(2.*PI*XLAM_W *XL)
c
         HILFB =  1./(4. * PI * alpha_a * r_a * r_a)
     2     + log(2*r_a/(r_a + r_i)) /(2.*PI*XLAM_W *XL)
c
         xk_u = 1./HILFB / A_B
         xk_i = 1./HILFA / A_B 
c
         xfakhemm = rpar(28)
         xhoch    = rpar(27)
         xfako2   = rpar(26)
c
         icount = 1
c
      end if
c
c-----------------------------------------------------------------------
c
      TI = u(3) + 273.15
c
      T_WAND = U(1)
      T_GAS  = U(2)
      T_KAT  = U(3)
      W_COS  = U(4)
      W_COG  = U(5)
      W_C3H6S= U(6)
      W_C3H6G= U(7)
      W_H2S  = U(8)
      W_H2G  = U(9)
      W_O2S  = U(10)
      W_O2G  = U(11)

      RHO_G = PGES * XM_LUFT/(RR * TI)
      Y_COS   = XM_LUFT/XM_CO   * W_COS
      Y_C3H6S = XM_LUFT/XM_C3H6 * W_C3H6S
      Y_H2S   = XM_LUFT/XM_H2   * W_H2S
      Y_O2S   = XM_LUFT/XM_O2   * W_O2S
c
      XXKZ_CO = XKZ_CO * DEXP(-EZ_CO/ TI)
      XXKN_CO = XKN_CO * DEXP( EN_CO/ TI)
      XXKZ_C3H6 = XKZ_C3H6 * DEXP(-EZ_C3H6/ TI)
      XXKN_C3H6 = XKN_C3H6 * DEXP( EN_C3H6/ TI)
      XXK_COC3H6 = XK_COC3H6 * DEXP( E_COC3H6/ TI)
c
c-----------------------------------------------------------------------
c
      G1CO = 1. + xfakhemm * XXKN_CO * abs(Y_COS)
      G1HC = 1. + xfakhemm * XXKN_C3H6 * abs(Y_C3H6S)
c
      R_CO   = XXKZ_CO  * Y_COS * Y_O2S/ TI/G1CO**xhoch 
      R_C3H6 = XXKZ_C3H6 * Y_C3H6S * Y_O2S/ TI/G1HC**xhoch
      R_H2   = XXKZ_CO * Y_H2S*Y_O2S/TI/G1CO**xhoch
      R_O2   = 0.5*R_CO + 4.5*R_C3H6 + 0.5*R_H2
      R_O2   = xfako2*R_O2
c
      q(1) = xk_i/r_i  * (T_KAT - T_WAND)
     2             - xk_u/r_i/r_i*r_a  * (T_WAND - tu) 
      q(2) = alpha_i * av * (T_KAT - T_GAS)
      q(3) = - xk_i/r_i * (T_KAT - T_WAND)
     2               - alpha_i * av * (T_KAT - T_GAS)
     3               + xhr_co * ax * 10. * R_CO
     4               + xhr_c3h6 * ax * 10. * R_C3H6
     5               + xhr_h2 * ax * 10. * R_H2
      q(4) = rho_g * beta_co * av * (W_COG - W_COS)
     2               - 1.0 * ax * 10. * xM_CO * R_CO
      q(5) = - rho_g * beta_co * av * (W_COG - W_COS)
c
      q(6) = rho_g * beta_C3H6 * av * (W_C3H6G - W_C3H6S)
     2               - 1.0 * ax * 10. * xM_C3H6 * R_C3H6
      q(7) = - rho_g * beta_C3H6 * av * (W_C3H6G - W_C3H6S)
c
      q(8) = rho_g * beta_H2 * av * (W_H2G - W_H2S)
     2               - 1.0 * ax * 10. * xM_H2 * R_H2
      q(9) = - rho_g * beta_H2 * av * (W_H2G - W_H2S)
c
      q(10) = rho_g * beta_O2 * av * (W_O2G - W_O2S)
     2               - 1.0 * ax * 10. * xM_O2 * R_O2
      q(11)= - rho_g * beta_O2 * av * (W_O2G - W_O2S)
c
c-----------------------------------------------------------------------
c
      do i = 1, Nr_of_PDEs
         ffun(i) = duxx(i,i) - conv(i) + q(i)
      end do
c
c-----------------------------------------------------------------------
c
c     End of equations
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Boundary_Conditions ( ix, x, t, u, alpha, beta, gamma, 
     2                                 delta ) 
c 
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine is an user routine for the local definition  
c     of the boundary describing functions alpha, beta, gamma.
c                                                             
c-----------------------------------------------------------------------
c                                             
c     Input parameter
c                                               
c     ix           actual position in actual grid           
c                  ix=1 indicates left boundary            
c                  ix=nx indicates right boundary         
c     x            actual position in space              
c     t            actual time                          
c     u            array of length npde, containing actual solution  
c                  for actual (x,t)-values             
c                                                   
c
c  Output parameter
c                 
c     alpha        actual values for boundary condition function alpha
c     beta         actual values for boundary condition function beta
c     gamma        actual values for boundary condition function gamma   
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c 
      dimension    u ( * )
c
      dimension    alpha ( * ), beta ( * ), gamma ( * ), delta ( * ) 
c
c-----------------------------------------------------------------------
c
      dimension    alphal ( 11, 1 ), betal ( 11, 1 ), gammal ( 11, 1 )
c
      dimension    alphar ( 11, 1 ), betar ( 11, 1 ), gammar ( 11, 1 )
c
c-----------------------------------------------------------------------
c
c     Begin of boundary conditions
c
c-----------------------------------------------------------------------
c
c     Initialization
c
      npde = Nr_of_PDEs
      nx   = Gridsize
c
c-----------------------------------------------------------------------
c
      istroem = ipar(1)
c
      r_i = rpar(1)/1000.
      pi  = 3.1415927
      a   = pi * r_i * r_i
      xmp = rpar(84)
      gz0 = xmp/a
c
      D_eff_co   = rpar(64)
      D_eff_c3h6 = rpar(63)
      D_eff_h2   = rpar(65)
      D_eff_o2   = rpar(66)
c
      xlam_w = rpar(25)
      xlam_k = rpar(19)
      xlam_g = rpar(12)
c
      cp_g = rpar(13)
c
      eps_g = rpar(6)
      eps_k = 1. - eps_g
      eps_w = 0.1
c
      gz = float(istroem) * gz0
c
      w_co_zu   = rpar(85)
      w_c3h6_zu = rpar(86)
      w_h2_zu   = rpar(87)
      w_o2_zu   = rpar(88)
c
      temp_anfang = rpar(99)
      temp_ende   = rpar(100)
c
      t_steig = rpar(95)
      t_zu = t/t_steig*(temp_ende - temp_anfang) + temp_anfang
c
      if(t.gt.100.d0) t_zu=temp_ende
c
c-----------------------------------------------------------------------
c
      alphal(1,1) = 0.0 
      betal(1,1)  = -XLAM_w* EPS_w
      gammal(1,1) = 0.0
c
      alphar(1,1) = 0.0 
      betar(1,1)  = -XLAM_w* EPS_w
      gammar(1,1) = 0.0
c
      alphal(3,1) = 0.0 
      betal(3,1)  = -XLAM_k*EPS_k
      gammal(3,1) = 0.0
c
      alphar(3,1) = 0.0 
      betar(3,1)  = -XLAM_k*EPS_k
      gammar(3,1) = 0.0
c
      alphal(4,1) = 0.0 
      betal(4,1)  = 1.0
      gammal(4,1) = 0.0
c
      alphar(4,1) = 0.0 
      betar(4,1)  = 1.0
      gammar(4,1) = 0.0
c
      alphal(6,1) = 0.0 
      betal(6,1)  = 1.0
      gammal(6,1) = 0.0
c
      alphar(6,1) = 0.0 
      betar(6,1)  = 1.0
      gammar(6,1) = 0.0
c
      alphal(8,1) = 0.0 
      betal(8,1)  = 1.0
      gammal(8,1) = 0.0
c
      alphar(8,1) = 0.0 
      betar(8,1)  = 1.0
      gammar(8,1) = 0.0
c
      alphal(10,1) = 0.0 
      betal(10,1)  = 1.0
      gammal(10,1) = 0.0
c
      alphar(10,1) = 0.0 
      betar(10,1)  = 1.0
      gammar(10,1) = 0.0
c
c-----------------------------------------------------------------------
c
      alphal(2,1) = gz*cp_g 
      betal(2,1)  = -eps_g*xlam_g
      gammal(2,1) = gz*cp_g*t_zu
c
      alphar(2,1) = 0.0 
      betar(2,1)  = -xlam_g*eps_g
      gammar(2,1) = 0.0
c
      alphal(5,1) = gz 
      betal(5,1)  = -D_eff_co*eps_g
      gammal(5,1) = gz*w_co_zu
c
      alphar(5,1) = 0.0 
      betar(5,1)  = -D_eff_co*eps_g
      gammar(5,1) = 0.0
c
      alphal(7,1) = gz 
      betal(7,1)  = -D_eff_c3h6*eps_g
      gammal(7,1) = gz*w_c3h6_zu
c
      alphar(7,1) = 0.0 
      betar(7,1)  = -D_eff_c3h6*eps_g
      gammar(7,1) = 0.0
c
      alphal(9,1) = gz 
      betal(9,1)  = -D_eff_h2*eps_g
      gammal(9,1) = gz*w_h2_zu
c
      alphar(9,1) = 0.0 
      betar(9,1)  = -D_eff_h2*eps_g
      gammar(9,1) = 0.0
c
      alphal(11,1) = gz 
      betal(11,1)  = -D_eff_o2*eps_g
      gammal(11,1) = gz*w_o2_zu
c
      alphar(11,1) = 0.0 
      betar(11,1)  = -D_eff_o2*eps_g
      gammar(11,1) = 0.0
c
c-----------------------------------------------------------------------
c
c     For left boundary
c
c-----------------------------------------------------------------------
c
      if ( ix .eq. 1 ) then
c
         do i = 1, npde
c
            alpha(i) = alphal(i,1)
            beta(i)  = betal(i,1)
            gamma(i) = gammal(i,1)
            delta(i) = 0.d0
c
         end do
c
      end if
c
c-----------------------------------------------------------------------
c
c     For right boundary
c
c-----------------------------------------------------------------------
c
      if ( ix .eq. nx ) then
c
         do i = 1, npde
c
            alpha(i) = alphar(i,1)
            beta(i)  = betar(i,1)
            gamma(i) = gammar(i,1) 
            delta(i) = 0.d0
c
         end do
c
      end if
c
c-----------------------------------------------------------------------
c
c     End of boundary conditions
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Define_Grid ( )  
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This is an user routine which creates a grid in the domain
c     [ x_Left, x_Right ]              
c                                                                   
c-----------------------------------------------------------------------
c
c     All grid parameters are transferred to other subroutines via 
c     the common block Grid_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Grid_Info
c
c     Grid         the gridpoints
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer     Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension   iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      common      / Grid_Info   / Grid
c
      dimension   Grid ( Max_Gridsize )
c
c-----------------------------------------------------------------------
c
c     Begin of initial grid
c
c-----------------------------------------------------------------------
c
      delta_x = ( x_Right - x_Left ) / dble ( Gridsize - 1 )
c
      Grid(1) = x_Left
c
      do i = 2, Gridsize
         Grid(i) = Grid(i-1) + delta_x
      end do   
c
c-----------------------------------------------------------------------
c
c     End of initial grid
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Fcn ( n, nzv, t, u, urhs, bv, irv, icv, Info )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     Computes the left-hand and right-hand side of the DA-system
c
c-----------------------------------------------------------------------
c
c     Computes the lhs and rhs between of the DAE.
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c     All grid parameters are transferred to other subroutines via
c     the common block Grid_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Grid_Info
c
c     Grid         the gridpoints
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      common     / Grid_Info / Grid
c
      dimension  Grid ( Max_Gridsize )
c
c-----------------------------------------------------------------------
c
      dimension u     ( n ), urhs ( n ) 
c
      dimension bv    ( * )
c
      dimension irv   ( * ), icv  ( * )
c
      dimension umid  ( Max_Size ), xmid ( Max_Size )
c
      dimension uhl   ( Max_Nr_of_PDEs ), uhm ( Max_Nr_of_PDEs ), 
     2          uhr   ( Max_Nr_of_PDEs ) 
c
      dimension ux    ( Max_Nr_of_PDEs ) 
c 
      dimension duxx  ( Max_Nr_of_PDEs, Max_Nr_of_PDEs) 
c
      dimension bmat  ( Max_Nr_of_PDEs, Max_Nr_of_PDEs )
c
      dimension ffun  ( Max_Nr_of_PDEs ) 
c
      dimension dvall ( Max_Nr_of_PDEs, Max_Nr_of_PDEs ), 
     2          dvalr ( Max_Nr_of_PDEs, Max_Nr_of_PDEs )
c
      dimension alpha ( Max_Nr_of_PDEs ), beta  ( Max_Nr_of_PDEs ), 
     2          gamma ( Max_Nr_of_PDEs ), delta ( Max_Nr_of_PDEs ) 
c
c-----------------------------------------------------------------------
c
c     Initialization
c
c-----------------------------------------------------------------------
c 
      Info = 0
c
      npde = Nr_of_PDEs
      nx   = Gridsize
c
      zero = 0.0
      half = 0.5
      one  = 1.0
      two  = 2.0
c
      ih = 0
c
c-----------------------------------------------------------------------
c
c     Set midpoints of intervals and associated solution values
c
c-----------------------------------------------------------------------
c
      do i = 1, nx - 1
c
         xmid(i) = half * ( Grid(i) + Grid(i+1) )
c
         k = ( i - 1 ) * npde
c
         do j = 1 + k, npde + k 
            umid(j) = half * ( u(j) + u(j+npde) )
         end do
c
      end do
c
c-----------------------------------------------------------------------
c
c     Left boundary
c
c-----------------------------------------------------------------------
c
      do j = 1, npde 
         uhm(j) = u(j)
         uhr(j) = u(npde+j)
      end do 
c
      dxr = Grid(2) - Grid(1)
c
c-----------------------------------------------------------------------
c
c     Get functions for left boundary conditions
c
c-----------------------------------------------------------------------
c
      call Boundary_Conditions ( 1, Grid(1), t, uhm, alpha, beta,
     2                           gamma, delta ) 
c
c-----------------------------------------------------------------------
c
c     Get diffusion coefficient at x(1)    
c
c-----------------------------------------------------------------------
c
      call Diff_Coeff ( 1, Grid(1), t, uhm(1), dvall ) 
c
c-----------------------------------------------------------------------
c
c     Get diffusion coefficient at x(1+1/2)    
c
c-----------------------------------------------------------------------
c
      call Diff_Coeff ( 1, xmid(1), t, umid(1), dvalr )
c
c-----------------------------------------------------------------------
c
c     Approximate ux at left boundary
c
c-----------------------------------------------------------------------
c
      do j = 1, npde
c
         if ( beta(j) .ne. zero ) then
c
            ux(j) = ( gamma(j) - alpha(j) * uhm(j) ) / beta(j)
c
         else
c
            ux(j) = ( uhr(j) - uhm(j) ) / dxr
c
         endif 
c
      end do
c
c-----------------------------------------------------------------------
c
c     Approximate duxx at left boundary
c
c-----------------------------------------------------------------------
c
      call duu_xx_Left ( Grid(1), Grid(2), uhm, uhr, xmid(1),
     2                   dvall, dvalr, ux, duxx )
c
c-----------------------------------------------------------------------
c
c     Get function for right hand side of PDE  
c
c-----------------------------------------------------------------------
c
      call PDE ( 1, Grid(1), t, uhm, ux, duxx, ffun, bmat ) 
c
c-----------------------------------------------------------------------
c
c     Form rhs for limex
c
c-----------------------------------------------------------------------
c
      do j = 1, npde 
c
         if ( beta(j) .ne. zero ) then
c
            urhs(j) = ffun(j)
c
         else
c
            urhs(j) = alpha(j) * uhm(j) - gamma(j)
c
         endif
c
      end do  
c
c-----------------------------------------------------------------------
c
c     Form lhs for limex
c
c-----------------------------------------------------------------------
c
      do j = 1, npde
c
         ih      = ih + 1
         irv(ih) = j
         icv(ih) = j
c
         if ( beta(j) .ne. zero ) then
c
            bv(ih) = bmat(j,j)
c
         else
c
            bv(ih) = zero 
c
         end if
c
      end do   
c
c-----------------------------------------------------------------------
c
c     Interior nodes
c
c-----------------------------------------------------------------------
c
c     Get diffusion coefficient at x(il-1/2)
c
c-----------------------------------------------------------------------
c
      do lx = 2, nx - 1
c
         lxn  = lx * npde
         lxn1 = lxn - npde
c
         dxl = dxr
         dxr = Grid(lx+1) - Grid(lx)
c
         do j = 1, npde
c
            do k = 1, npde
               dvall(j,k) = dvalr(j,k)
            end do
c
            uhl(j)  = uhm(j)
            uhm(j)  = uhr(j)
            uhr(j)  = u(lxn+j)
c
         end do   
c
c-----------------------------------------------------------------------
c
c     Get diffusion coefficient at x(lx+1/2)    
c
c-----------------------------------------------------------------------
c
         call Diff_Coeff ( lx, xmid(lx), t, umid(lxn1+1), dvalr ) 
c
c-----------------------------------------------------------------------
c
c     Approximate du/dx
c
c-----------------------------------------------------------------------
c
         call du_dx ( Grid(lx-1), Grid(lx), Grid(lx+1), uhl, uhm, uhr,
     2                ux )
c
c-----------------------------------------------------------------------
c
c                  2    2
c     Approximate d u/dx 
c
c-----------------------------------------------------------------------
c
         call duu_xx ( Grid(lx-1), Grid(lx), Grid(lx+1), uhl, uhm, uhr,
     2                 xmid(lx-1), xmid(lx), dvall, dvalr, duxx )
c
c-----------------------------------------------------------------------
c
c     Get function for right hand side of PDE  
c
c-----------------------------------------------------------------------
c
         call PDE ( lx, Grid(lx), t, uhm, ux, duxx, ffun, bmat ) 
c
c-----------------------------------------------------------------------
c
c     Form rhs and lhs for LIMEX
c
c-----------------------------------------------------------------------
c
         do j = 1, npde
            urhs(lxn1+j) = ffun(j)
         end do
c
         do j = 1, npde
c
            ih      = ih + 1
            irv(ih) = lxn1 + j
            icv(ih) = irv(ih) 
            bv(ih)  = bmat(j,j)
c
         end do  
c
      end do  
c
c-----------------------------------------------------------------------
c
c     Right boundary
c
c-----------------------------------------------------------------------
c
      lxn  = nx * npde
      lxn1 = lxn - npde
c
      dxl = dxr
c
      do j = 1, npde
c
         do k = 1, npde
            dvall(j,k) = dvalr(j,k)
         end do
c
         uhl(j) = uhm(j)
         uhm(j) = uhr(j)
c
      end do    
c
c-----------------------------------------------------------------------
c
c     Get functions for right boundary conditions
c
c-----------------------------------------------------------------------
c
      call Boundary_Conditions ( nx, Grid(nx), t, uhm, alpha, beta,
     2                           gamma, delta )
c
c-----------------------------------------------------------------------
c
c     Get diffusion coefficient at x(nx)    
c
c-----------------------------------------------------------------------
c
      call Diff_Coeff ( nx, Grid(nx), t, uhm(1), dvalr ) 
c
c-----------------------------------------------------------------------
c
c     Approximate ux at right boundary
c
c-----------------------------------------------------------------------
c
      do j = 1, npde
c
         if ( beta(j) .ne. zero ) then
c
            ux(j) = ( gamma(j) - alpha(j) * uhm(j) ) / beta(j)
c
         else
c
            ux(j) = ( uhm(j) - uhl(j) ) / dxl
c
         end if 
c
      end do   
c
c-----------------------------------------------------------------------
c
c     Approximate duxx at right boundary
c
c-----------------------------------------------------------------------
c
      call duu_xx_Right ( Grid(nx-1), Grid(nx), uhl, uhm,
     2                    xmid(nx-1), dvall, dvalr, ux, duxx )
c
c-----------------------------------------------------------------------
c
c     Get functions for right hand side of PDE  
c
c-----------------------------------------------------------------------
c
      call PDE ( nx, Grid(nx), t, uhm, ux, duxx, ffun, bmat )
c
c-----------------------------------------------------------------------
c
c     Form rhs for LIMEX
c
c-----------------------------------------------------------------------
c
      do j = 1, npde 
c
         if ( beta(j) .ne. zero ) then
c
            urhs(lxn1+j) = ffun(j)
c
         else
c        
            urhs(lxn1+j) = alpha(j) * uhm(j) - gamma(j)
c
         endif
c
      end do   
c
      do j = 1, npde
c
         ih      = ih + 1
         irv(ih) = lxn1 + j
         icv(ih) = irv(ih)
c
         if( beta(j) .ne. zero ) then
c
            bv(ih) = bmat(j,j)
c
         else
c
            bv(ih) = zero
c
         end if
c
      end do
c 
c-----------------------------------------------------------------------
c
c     Count the non-zero variables in matrix B
c
c-----------------------------------------------------------------------
c
      nzv = ih
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
      dimension        Jac ( LDJAC, * )
c
      dimension        y ( * ), ys ( * )
c
c-----------------------------------------------------------------------
c
c     Dummy routine for the catalysator model. 
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine du_dx ( xl, xm, xr, ul, um, ur, ux )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine approximates the gradients of u by
c     
c     d u(x,t)                1
c     -------- = ---------------------------  *
c       dx       Delta(x(i-1))+Delta(x(i+1))
c
c              Delta(x(i-1))                 Delta(x(i))
c            ( -------------(u(i+1)-u(i)) + -------------(u(i)-u(i-1)) )
c               Delta(x(i))                 Delta(x(i-1))
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension ux ( * ), ul ( * ), um ( * ), ur ( * )
c 
c-----------------------------------------------------------------------
c
      dxl  = xm - xl
      dxr  = xr - xm
c
      dxl_dxr = dxl / dxr
      dxr_dxl = dxr / dxl
c
      sum_dx_rec = 1.0 / ( xr - xl )
c 
      do j = 1, Nr_of_PDEs
         ux(j) = sum_dx_rec * 
     2    ( dxl_dxr * ( ur(j) - um(j) ) + dxr_dxl * ( um(j) - ul(j) ) )
      end do  
c 
c-----------------------------------------------------------------------
c
      return
      end
      subroutine duu_xx ( xl, xm, xr, ul, um, ur, xlh, xrh, dlh, drh, 
     2                    duxx )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine approximates the second derivatives by
c
c      2
c     d u(x,t)                2
c     -------- = -------------------------  *
c         2      Delta(x(i-1))+Delta(x(i))
c       dX
c
c                   u(i+1)-u(i)   u(i)-u(i-1)
c                (  ----------- - -------------  )
c                   Delta(x(i))   Delta(x(i-1))
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension  ul ( * ), um ( * ), ur ( * ) 
c
      dimension  dlh ( Max_Nr_of_PDEs, * ), drh ( Max_Nr_of_PDEs, * ) 
c
      dimension  duxx ( Max_Nr_of_PDEs, * )
c 
c-----------------------------------------------------------------------
c
      npde = Nr_of_PDEs
c
      dxh = 2.0 / ( xr - xl )
      dxl = dxh / ( xm - xl )
      dxr = dxh / ( xr - xm )
c
      do j = 1, npde
c
         uxr = ( ur(j) - um(j) ) * dxr  
         uxl = ( um(j) - ul(j) ) * dxl
c
         do k = 1, npde
            duxx(k,j) = drh(k,j) * uxr - dlh(k,j) * uxl
         end do  
c
      end do   
c 
c-----------------------------------------------------------------------
c
      return
      end
      subroutine duu_xx_Left ( xm, xr, um, ur, xrh, dlh, drh, ux, duxx )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine approximates the second derivatives on the
c     left boundary
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension  um ( * ), ur ( * )
c
      dimension  dlh ( Max_Nr_of_PDEs, * ), drh ( Max_Nr_of_PDEs, * ) 
c
      dimension  ux ( * )
c
      dimension  duxx ( Max_Nr_of_PDEs, * )
c
      dimension  uxl ( Max_Nr_Of_PDEs ), uxr ( Max_Nr_of_PDEs )
c
c-----------------------------------------------------------------------
c 
      npde = Nr_of_PDEs
c
      dxr  = xr - xm
      dxri = 1.0 / dxr
      dxh  = 2.0 * dxri
c
      do j = 1, npde
c
         uxr(j) = ( ur(j) - um(j) ) * dxri
         uxl(j) = ux(j)
c
      end do
c
      do j = 1, npde
         do k = 1, npde
            duxx(k,j) = ( drh(k,j) * uxr(j) - dlh(k,j) * uxl(j) ) * dxh
         end do   
      end do  
c 
c-----------------------------------------------------------------------
c
      return
      end
      subroutine duu_xx_Right ( xl, xm, ul, um, xlh, dlh, drh, ux, 
     2                          duxx )
c
      implicit double precision ( a-h, o-z )
c
c-----------------------------------------------------------------------
c
c     This subroutine approximates the second derivatives on the
c     right boundary
c
c-----------------------------------------------------------------------
c
c     All model parameters are transferred to other subroutines via
c     the common block Problem_Info.
c
c-----------------------------------------------------------------------
c
c     Parameters in common block Problem_Info
c
c     Nr_of_PDEs   dimension of partial differential equation
c     Gridsize     size of grid for space discretization
c     x_Left       left bound of domain
c     x_Right      right bound of domain
c     iPar         array for real user parameters
c     rPar         array for integer user parameters
c
c-----------------------------------------------------------------------
c
      parameter ( Max_Nr_of_PDEs  =                            11,
     2            Max_Gridsize    =                           650,
     3            Max_Size        = Max_Nr_of_PDEs * Max_Gridsize,
     4            Nr_Int_Par      =                           100,
     5            Nr_Real_Par     =                           100 )
c
c-----------------------------------------------------------------------
c
      common     / Problem_Info / Nr_of_PDEs, Gridsize, x_Left, x_Right, 
     2                            iPar, rPar
c
      integer    Nr_of_PDEs, Gridsize
c
      double precision x_Left, x_Right
c
      dimension  iPar ( Nr_Int_Par ), rPar ( Nr_Real_Par )
c
c-----------------------------------------------------------------------
c
      dimension  ul ( * ), um ( * )
c
      dimension  dlh ( Max_Nr_of_PDEs, * ), drh ( Max_Nr_of_PDEs, * )
c
      dimension  ux ( * )
c
      dimension  duxx ( Max_Nr_of_PDEs, * )
c
      dimension  uxl ( Max_Nr_of_PDEs ), uxr ( Max_Nr_of_PDEs )
c
c-----------------------------------------------------------------------
c 
      npde = Nr_of_PDEs
c
      dxl  = xm - xl
      dxli = 1.0 / dxl
      dxh  = 2.0 * dxli
c
      do j = 1, npde
c
         uxr(j) = ux(j)
         uxl(j) = ( um(j) - ul(j) ) * dxli
c
      end do   
c
      do j = 1, npde
         do k = 1, npde
            duxx(k,j) = ( drh(k,j) * uxr(j) - dlh(k,j) * uxl(j) ) * dxh
         end do   
      end do   
c 
c-----------------------------------------------------------------------
c
      return
      end
