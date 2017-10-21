************************************************************************
*
*     File  lusoltest1.f
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     This program solves a sparse linear system Ax = b using LUSOL.
*
*     09 Oct 2003: First version.
*                  Michael Saunders, SOL, Stanford University.
*                  Tests lu1fac and lu6sol on example data A, b.
*                  Data files A1764.txt, b1764.txt provided by
*                  Tao Gang, tao-gang@utulsa.edu.
*     07 Feb 2004: A1764.txt contains data as  Aij  i  j.
*                  Tao Gang sent new files
*                     A6805.txt   b6805.txt
*                    A10009.txt  b10009.txt
*                  with entries  i  j  Aij.
*                  Now read generic files Afile.txt bfile.txt in that form.
*-----------------------------------------------------------------------

      program          lusoltest1
      implicit	       none

      ! Externals      

      integer          idamax
      double precision dnormi
      external         dnormi, idamax

      ! Storage for A, b

      integer          maxm, maxn, maxnz
      parameter	      (maxn   =  25000,
     $                 maxm   =   maxn,
     $                 maxnz  = 500000)
      double precision Aij(maxnz), b( maxm ), xexact( maxn )
      integer          iA(maxnz), jA(maxnz)

      ! Storage for LUSOL.

      integer          lena
      parameter	      (lena = 20000000)
      double precision parmlu(30),  a(lena), w(maxn)
      integer          luparm(30),
     &                 indc(lena),  indr(lena),
     &                 lenc(maxn),  lenr(maxm),
     &                 ip(maxm),    iq(maxn),
     &                 iploc(maxn), iqloc(maxm),
     &                 ipinv(maxm), iqinv(maxn),
     &                 locc(maxn),  locr(maxm)

      ! Local storage.

      integer          i     , inform, j     , k     ,
     &                 m     , mode  , matfil, msgfil,
     &                 n     , nelem , nnzero, rhsfil
      double precision Ak    , Amax  , factol, test  ,
     &                 bnorm , rnorm , xnorm 
      double precision rhs(maxm)
      double precision r(maxn), x(maxn)


      !-----------------------------------------------------------------
      ! Define files.
      !-----------------------------------------------------------------
      msgfil = 6
      matfil = 12
      rhsfil = 13
      write(msgfil, '(a)') ' ', ' ==========',
     &                          ' lusoltest1',
     &                          ' =========='
      open( unit=matfil, file='Afile.txt', status='OLD' )
      open( unit=rhsfil, file='bfile.txt', status='OLD' )

      !-----------------------------------------------------------------
      ! Read data for A
      !-----------------------------------------------------------------
      nnzero = 0
      m      = 0
      n      = 0

      do k = 1, maxnz
         read( matfil, *, end=300 ) i, j, Ak
         if (k .gt. maxnz) go to 250
         nnzero = k
         iA(k)  = i
         jA(k)  = j
         Aij(k) = Ak
         m      = max( m, i )
         n      = max( n, j )
      end do
      close( matfil )
      go to 300

  250 write(msgfil, *) 'Too much data in A file.  Increase maxnz'
      write(msgfil, *) 'Current maxnz =', maxnz
      go to 900

  300 write(msgfil, *) 'A  read successfully'
      write(msgfil, '(a, i8)')
     &   ' m      =', m,
     &   ' n      =', n,
     &   ' nnzero =', nnzero

      if (m .gt. maxm  .or.  n .gt. maxn) then
         write(msgfil, '(a)') ' ', ' These exceed maxm or maxn'
         go to 900
      end if

      !-----------------------------------------------------------------
      ! Read data for b
      !-----------------------------------------------------------------
      do i = 1, m
         read( rhsfil, *, end=350 ) b(i)
      end do
      close( rhsfil )
      write(msgfil, '(a)') ' ', ' b  read successfully', ' '
      go to 400

  350 write(msgfil, *) 'Not enough data in b file.'
      go to 900

      !------------------------------------------------------------------
      ! Set parameters for LUSOL's lu1fac.
      !------------------------------------------------------------------
  400 factol    = 2.0d+0   ! Stability tolerance

      luparm(1) = msgfil   ! File number for printed messages
      luparm(2) = 10       ! Print level. >= 0 to get singularity info.
                           !              >=10 to get more LU statistics.
                           !              >=50 to get info on each pivot.
      luparm(3) = 5        ! maxcol
      luparm(6) = 1        ! Threshold Pivoting: 0 = TPP, 1 = TRP, 2 = TCP
      luparm(8) = 1        ! keepLU
      parmlu(1) = factol   ! Ltol1:  max |Lij| during Factor
      parmlu(2) = factol   ! Ltol2:  max |Lij| during Update 
      parmlu(3) = 3.0d-13  ! small:  drop tolerance
      parmlu(4) = 3.7d-11  ! Utol1:  absolute tol for small Uii
      parmlu(5) = 3.7d-11  ! Utol2:  relative tol for small Uii
      parmlu(6) = 3.0d+0   ! Uspace: 
      parmlu(7) = 0.3d+0   ! dens1
      parmlu(8) = 0.5d+0   ! dens2

      !-----------------------------------------------------------------
      ! Load A into (a, indc, indr).
      !-----------------------------------------------------------------
      do k = 1, nnzero
         a(k)    = Aij(k)
         indc(k) = iA(k)
         indr(k) = jA(k)
      end do

      !------------------------------------------------------------------
      ! Factor  A = L U.
      !------------------------------------------------------------------
      nelem = nnzero

      call lu1fac( m    , n    , nelem, lena , luparm, parmlu,
     $             a    , indc , indr , ip   , iq    ,
     $             lenc , lenr , locc , locr ,
     $             iploc, iqloc, ipinv, iqinv, w     , inform )

      if (inform .gt. 1) then
         write(msgfil, '(/ a, i4)') ' lu1fac error.   inform =', inform
         go to 900
      end if

      Amax   = parmlu(10) ! This is the largest element in A.
                          ! We use it below as an estimate of ||A||_inf,
                          ! even though it isn't a proper norm.

      !------------------------------------------------------------------
      ! SOLVE  A x = b.
      ! Save b first because
      ! call lu6sol( 5, ... ) overwrites rhs when computing x.
      !------------------------------------------------------------------
      do i = 1, m
         rhs(i) = b(i)
      end do

      mode   = 5
      call lu6sol( mode, m, n, rhs, x,
     $             lena, luparm, parmlu,
     $             a, indc, indr, ip, iq,
     $             lenc, lenr, locc, locr,
     $             inform )

*     ------------------------------------------------------------------
*     Set r = b - Ax.
*     Find norm of r and x.
*     ------------------------------------------------------------------
      do i = 1, m
         r(i) = b(i)
      end do

      do k = 1, nnzero
         i    = iA(k)
         j    = jA(k)
         r(i) = r(i) - Aij(k)*x(j)
      end do

      bnorm  = dnormi( m, b )
      rnorm  = dnormi( m, r )
      xnorm  = dnormi( n, x )
      test   = rnorm / (Amax*xnorm)

      write(msgfil, *) ' '
      write(msgfil, '(a, 1p, e10.1)')
     &   ' bnorm =', bnorm,
     &   ' rnorm =', rnorm,
     &   ' xnorm =', xnorm
      if     (test .le. 1d-8) then
         write(msgfil, *) 'rnorm is small.  The test seems successful.'
      elseif (test .le. 1d-4) then
         write(msgfil, *) 'rnorm is not very small.  '
         write(msgfil, *) 'The LU factors may not be good enough.'
         write(msgfil, *) 'Try a smaller factol and/or Rook Pivoting'
      else
         write(msgfil, *) 'rnorm is rather large.  '
         write(msgfil, *) 'The LU factorization has probably failed.'
         write(msgfil, *) 'Try a smaller factol and/or Rook Pivoting'
      end if

  900 stop

      end ! program lusoltest1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnormi( n, x )

      implicit         none
      integer          n
      double precision x(*)

*     ===============================================================
*     dnormi  returns the infinity-norm of the vector x.
*     ===============================================================
      integer          j

      dnormi = 0.0d+0

      do j = 1, n
         dnormi = max( dnormi, abs(x(j)) )
      end do

      end ! double precision function dnormi
