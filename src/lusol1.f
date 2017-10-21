!***********************************************************************
!
!     File lusol1.f
!
!     lu1fac   lu1fad   lu1gau   lu1mar   lu1mRP   lu1mCP   lu1mSP
!     lu1pen   lu1mxc   lu1mxr   lu1or1   lu1or2   lu1or3   lu1or4
!     lu1pq1   lu1pq2   lu1pq3   lu1rec   lu1slk
!     lu1ful   lu1DPP   lu1DCP
!
! 26 Apr 2002: TCP implemented using heap data structure.
! 01 May 2002: lu1DCP implemented.
! 07 May 2002: lu1mxc must put 0.0 at top of empty columns.
! 09 May 2002: lu1mCP implements Markowitz with cols searched
!              in heap order.
!              Often faster (searching 20 or 40 cols) but more dense.
! 11 Jun 2002: TRP implemented.
!              lu1mRP implements Markowitz with Threshold Rook Pivoting.
!              lu1mxc maintains max col elements.  (Previously lu1max.)
!              lu1mxr maintains max row elements.
! 12 Jun 2002: lu1mCP seems too slow on big problems (e.g. memplus).
!              Disabled it for the moment.  (Use lu1mar + TCP.)
! 14 Dec 2002: TSP implemented.
!              lu1mSP implements Markowitz with
!              Threshold Symmetric Pivoting.
! 07 Mar 2003: character*1, character*2 changed to f90 form.
!              Comments changed from * in column to ! in column 1.
!              Comments kept within column 72 to avoid compiler warning.
! 19 Dec 2004: Hdelete(...) has new input argument Hlenin.
! 21 Dec 2004: Print Ltol and Lmax with e10.2 instead of e10.1.
! 26 Mar 2006: lu1fad: Ignore nsing from lu1ful.
!              lu1DPP: nsing redefined (but not used by lu1fad).
!              lu1DCP: nsing redefined (but not used by lu1fad).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1fac( m    , n    , nelem, lena , luparm, parmlu,
     &                   a    , indc , indr , ip   , iq    ,
     &                   lenc , lenr , locc , locr ,
     &                   iploc, iqloc, ipinv, iqinv, w     , inform )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena)   , w(n)
      integer            indc(lena), indr(lena), ip(m)   , iq(n),
     &                   lenc(n)   , lenr(m)   ,
     &                   iploc(n)  , iqloc(m)  , ipinv(m), iqinv(n)
      integer            locc(n)   , locr(m)
!     ------------------------------------------------------------------
!     lu1fac computes a factorization A = L*U, where A is a sparse
!     matrix with m rows and n columns, P*L*P' is lower triangular
!     and P*U*Q is upper triangular for certain permutations P, Q
!     (which are returned in the arrays ip, iq).
!     Stability is ensured by limiting the size of the elements of L.
!
!     The nonzeros of A are input via the parallel arrays a, indc, indr,
!     which should contain nelem entries of the form    aij,    i,    j
!     in any order.  There should be no duplicate pairs         i,    j.
!
!     ******************************************************************
!     *        Beware !!!   The row indices i must be in indc,         *
!     *              and the column indices j must be in indr.         *
!     *              (Not the other way round!)                        *
!     ******************************************************************
!
!     It does not matter if some of the entries in a(*) are zero.
!     Entries satisfying  abs( a(i) ) .le. parmlu(3)  are ignored.
!     Other parameters in luparm and parmlu are described below.
!
!     The matrix A may be singular.  On exit, nsing = luparm(11) gives
!     the number of apparent singularities.  This is the number of
!     "small" diagonals of the permuted factor U, as judged by
!     the input tolerances Utol1 = parmlu(4) and  Utol2 = parmlu(5).
!     The diagonal element diagj associated with column j of A is
!     "small" if
!                 abs( diagj ) .le. Utol1
!     or
!                 abs( diagj ) .le. Utol2 * max( uj ),
!
!     where max( uj ) is the maximum element in the j-th column of U.
!     The position of such elements is returned in w(*).  In general,
!     w(j) = + max( uj ),  but if column j is a singularity,
!     w(j) = - max( uj ).  Thus, w(j) .le. 0 if column j appears to be
!     dependent on the other columns of A.
!
!     NOTE: lu1fac (like certain other sparse LU packages) does not
!     treat dense columns efficiently.  This means it will be slow
!     on "arrow matrices" of the form
!                  A = (x       a)
!                      (  x     b)
!                      (    x   c)
!                      (      x d)
!                      (x x x x e)
!     if the numerical values in the dense column allow it to be
!     chosen LATE in the pivot order.
!
!     With TPP (Threshold Partial Pivoting), the dense column is
!     likely to be chosen late.
!
!     With TCP (Threshold Complete Pivoting), if any of a,b,c,d
!     is significantly larger than other elements of A, it will
!     be chosen as the first pivot and the dense column will be
!     eliminated, giving reasonably sparse factors.
!     However, if element e is so big that TCP chooses it, the factors
!     will become dense.  (It's hard to win on these examples!)
!     ==================================================================
!
!
!     Notes on the array names
!     ------------------------
!
!     During the LU factorization, the sparsity pattern of the matrix
!     being factored is stored twice: in a column list and a row list.
!
!     The column list is ( a, indc, locc, lenc )
!     where
!           a(*)    holds the nonzeros,
!           indc(*) holds the indices for the column list,
!           locc(j) points to the start of column j in a(*) and indc(*),
!           lenc(j) is the number of nonzeros in column j.
!
!     The row list is    (    indr, locr, lenr )
!     where
!           indr(*) holds the indices for the row list,
!           locr(i) points to the start of row i in indr(*),
!           lenr(i) is the number of nonzeros in row i.
!
!
!     At all stages of the LU factorization, ip contains a complete
!     row permutation.  At the start of stage k,  ip(1), ..., ip(k-1)
!     are the first k-1 rows of the final row permutation P.
!     The remaining rows are stored in an ordered list
!                          ( ip, iploc, ipinv )
!     where
!           iploc(nz) points to the start in ip(*) of the set of rows
!                     that currently contain nz nonzeros,
!           ipinv(i)  points to the position of row i in ip(*).
!
!     For example,
!           iploc(1) = k   (and this is where rows of length 1 begin),
!           iploc(2) = k+p  if there are p rows of length 1
!                          (and this is where rows of length 2 begin).
!
!     Similarly for iq, iqloc, iqinv.
!     ==================================================================
!
!
!     00 Jun 1983  Original version.
!     00 Jul 1987  nrank  saved in luparm(16).
!     12 Apr 1989  ipinv, iqinv added as workspace.
!     26 Apr 1989  maxtie replaced by maxcol in Markowitz search.
!     16 Mar 1992  jumin  saved in luparm(19).
!     10 Jun 1992  lu1fad has to move empty rows and cols to the bottom
!                  (via lu1pq3) before doing the dense LU.
!     12 Jun 1992  Deleted dense LU (lu1ful, lu1vlu).
!     25 Oct 1993  keepLU implemented.
!     07 Feb 1994  Added new dense LU (lu1ful, lu1den).
!     21 Dec 1994  Bugs fixed in lu1fad (nrank) and lu1ful (ipvt).
!     08 Aug 1995  Use ip instead of w as parameter to lu1or3 (for F90).
!     13 Sep 2000  TPP and TCP options implemented.
!     17 Oct 2000  Fixed troubles due to A = empty matrix (Todd Munson).
!     01 Dec 2000  Save Lmax, Umax, etc. after both lu1fad and lu6chk.
!                  lu1fad sets them when keepLU = false.
!                  lu6chk sets them otherwise, and includes items
!                  from the dense LU.
!     11 Mar 2001  lu6chk now looks at diag(U) when keepLU = false.
!     26 Apr 2002  New TCP implementation using heap routines to
!                  store largest element in each column.
!                  New workspace arrays Ha, Hj, Hk required.
!                  For compatibility, borrow space from a, indc, indr
!                  rather than adding new input parameters.
!     01 May 2002  lu1den changed to lu1DPP (dense partial  pivoting).
!                  lu1DCP implemented       (dense complete pivoting).
!                  Both TPP and TCP now switch to dense mode and end.
!
!     Systems Optimization Laboratory, Stanford University.
!  ---------------------------------------------------------------------
!
!
!  INPUT PARAMETERS
!
!  m      (not altered) is the number of rows in A.
!  n      (not altered) is the number of columns in A.
!  nelem  (not altered) is the number of matrix entries given in
!         the arrays a, indc, indr.
!  lena   (not altered) is the dimension of  a, indc, indr.
!         This should be significantly larger than nelem.
!         Typically one should have
!            lena > max( 2*nelem, 10*m, 10*n, 10000 )
!         but some applications may need more.
!         On machines with virtual memory it is safe to have
!         lena "far bigger than necessary", since not all of the
!         arrays will be used.
!  a      (overwritten) contains entries   Aij  in   a(1:nelem).
!  indc   (overwritten) contains the indices i in indc(1:nelem).
!  indr   (overwritten) contains the indices j in indr(1:nelem).
!
!  luparm input parameters:                                Typical value
!
!  luparm( 1) = nout     File number for printed messages.         6
!
!  luparm( 2) = lprint   Print level.                              0
!                   <  0 suppresses output.
!                   =  0 gives error messages.
!                  >= 10 gives statistics about the LU factors.
!                  >= 50 gives debug output from lu1fac
!                        (the pivot row and column and the
!                        no. of rows and columns involved at
!                        each elimination step).
!
!  luparm( 3) = maxcol   lu1fac: maximum number of columns         5
!                        searched allowed in a Markowitz-type
!                        search for the next pivot element.
!                        For some of the factorization, the
!                        number of rows searched is
!                        maxrow = maxcol - 1.
!
!  luparm( 6) = 0    =>  TPP: Threshold Partial   Pivoting.        0
!             = 1    =>  TRP: Threshold Rook      Pivoting.
!             = 2    =>  TCP: Threshold Complete  Pivoting.
!             = 3    =>  TSP: Threshold Symmetric Pivoting.
!             = 4    =>  TDP: Threshold Diagonal  Pivoting.
!                             (TDP not yet implemented).
!                        TRP and TCP are more expensive than TPP but
!                        more stable and better at revealing rank.
!                        Take care with setting parmlu(1), especially
!                        with TCP.
!                        NOTE: TSP and TDP are for symmetric matrices
!                        that are either definite or quasi-definite.
!                        TSP is effectively TRP for symmetric matrices.
!                        TDP is effectively TCP for symmetric matrices.
!
!  luparm( 8) = keepLU   lu1fac: keepLU = 1 means the numerical    1
!                        factors will be computed if possible.
!                        keepLU = 0 means L and U will be discarded
!                        but other information such as the row and
!                        column permutations will be returned.
!                        The latter option requires less storage.
!
!  parmlu input parameters:                                Typical value
!
!  parmlu( 1) = Ltol1    Max Lij allowed during Factor.
!                                                  TPP     10.0 or 100.0
!                                                  TRP      4.0 or  10.0
!                                                  TCP      5.0 or  10.0
!                                                  TSP      4.0 or  10.0
!                        With TRP and TCP (Rook and Complete Pivoting),
!                        values less than 25.0 may be expensive
!                        on badly scaled data.  However,
!                        values less than 10.0 may be needed
!                        to obtain a reliable rank-revealing
!                        factorization.
!  parmlu( 2) = Ltol2    Max Lij allowed during Updates.            10.0
!                        during updates.
!  parmlu( 3) = small    Absolute tolerance for       eps**0.8 = 3.0d-13
!                        treating reals as zero.
!  parmlu( 4) = Utol1    Absolute tol for flagging    eps**0.67= 3.7d-11
!                        small diagonals of U.
!  parmlu( 5) = Utol2    Relative tol for flagging    eps**0.67= 3.7d-11
!                        small diagonals of U.
!                        (eps = machine precision)
!  parmlu( 6) = Uspace   Factor limiting waste space in  U.      3.0
!                        In lu1fac, the row or column lists
!                        are compressed if their length
!                        exceeds Uspace times the length of
!                        either file after the last compression.
!  parmlu( 7) = dens1    The density at which the Markowitz      0.3
!                        pivot strategy should search maxcol
!                        columns and no rows.
!                        (Use 0.3 unless you are experimenting
!                        with the pivot strategy.)
!  parmlu( 8) = dens2    the density at which the Markowitz      0.5
!                        strategy should search only 1 column,
!                        or (if storage is available)
!                        the density at which all remaining
!                        rows and columns will be processed
!                        by a dense LU code.
!                        For example, if dens2 = 0.1 and lena is
!                        large enough, a dense LU will be used
!                        once more than 10 per cent of the
!                        remaining matrix is nonzero.
!
!
!  OUTPUT PARAMETERS
!
!  a, indc, indr     contain the nonzero entries in the LU factors of A.
!         If keepLU = 1, they are in a form suitable for use
!         by other parts of the LUSOL package, such as lu6sol.
!         U is stored by rows at the start of a, indr.
!         L is stored by cols at the end   of a, indc.
!         If keepLU = 0, only the diagonals of U are stored, at the
!         end of a.
!  ip, iq    are the row and column permutations defining the
!         pivot order.  For example, row ip(1) and column iq(1)
!         defines the first diagonal of U.
!  lenc(1:numl0) contains the number of entries in nontrivial
!         columns of L (in pivot order).
!  lenr(1:m) contains the number of entries in each row of U
!         (in original order).
!  locc(1:n) = 0 (ready for the LU update routines).
!  locr(1:m) points to the beginning of the rows of U in a, indr.
!  iploc, iqloc, ipinv, iqinv  are undefined.
!  w      indicates singularity as described above.
!  inform = 0 if the LU factors were obtained successfully.
!         = 1 if U appears to be singular, as judged by lu6chk.
!         = 3 if some index pair indc(l), indr(l) lies outside
!             the matrix dimensions 1:m , 1:n.
!         = 4 if some index pair indc(l), indr(l) duplicates
!             another such pair.
!         = 7 if the arrays a, indc, indr were not large enough.
!             Their length "lena" should be increase to at least
!             the value "minlen" given in luparm(13).
!         = 8 if there was some other fatal error.  (Shouldn't happen!)
!         = 9 if no diagonal pivot could be found with TSP or TDP.
!             The matrix must not be sufficiently definite
!             or quasi-definite.
!
!  luparm output parameters:
!
!  luparm(10) = inform   Return code from last call to any LU routine.
!  luparm(11) = nsing    No. of singularities marked in the
!                        output array w(*).
!  luparm(12) = jsing    Column index of last singularity.
!  luparm(13) = minlen   Minimum recommended value for  lena.
!  luparm(14) = maxlen   ?
!  luparm(15) = nupdat   No. of updates performed by the lu8 routines.
!  luparm(16) = nrank    No. of nonempty rows of U.
!  luparm(17) = ndens1   No. of columns remaining when the density of
!                        the matrix being factorized reached dens1.
!  luparm(18) = ndens2   No. of columns remaining when the density of
!                        the matrix being factorized reached dens2.
!  luparm(19) = jumin    The column index associated with DUmin.
!  luparm(20) = numL0    No. of columns in initial  L.
!  luparm(21) = lenL0    Size of initial  L  (no. of nonzeros).
!  luparm(22) = lenU0    Size of initial  U.
!  luparm(23) = lenL     Size of current  L.
!  luparm(24) = lenU     Size of current  U.
!  luparm(25) = lrow     Length of row file.
!  luparm(26) = ncp      No. of compressions of LU data structures.
!  luparm(27) = mersum   lu1fac: sum of Markowitz merit counts.
!  luparm(28) = nUtri    lu1fac: triangular rows in U.
!  luparm(29) = nLtri    lu1fac: triangular rows in L.
!  luparm(30) =
!
!
!
!  parmlu output parameters:
!
!  parmlu(10) = Amax     Maximum element in  A.
!  parmlu(11) = Lmax     Maximum multiplier in current  L.
!  parmlu(12) = Umax     Maximum element in current  U.
!  parmlu(13) = DUmax    Maximum diagonal in  U.
!  parmlu(14) = DUmin    Minimum diagonal in  U.
!  parmlu(15) = Akmax    Maximum element generated at any stage
!                        during TCP factorization.
!  parmlu(16) = growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax
!  parmlu(17) =
!  parmlu(18) =
!  parmlu(19) =
!  parmlu(20) = resid    lu6sol: residual after solve with U or U'.
!  ...
!  parmlu(30) =
!  ---------------------------------------------------------------------

      character          mnkey*1, kPiv(0:3)*2
      integer            lPiv
      logical            keepLU, TCP, TPP, TRP, TSP
      double precision   Lmax, Ltol

      double precision   zero         ,  one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

!     Grab relevant input parameters.

      nelem0 = nelem
      nout   = luparm(1)
      lprint = luparm(2)
      lPiv   = luparm(6)
      keepLU = luparm(8) .ne. 0

      Ltol   = parmlu(1)    ! Limit on size of Lij
      small  = parmlu(3)    ! Drop tolerance

      TPP    = lPiv .eq. 0  ! Threshold Partial   Pivoting (normal).
      TRP    = lPiv .eq. 1  ! Threshold Rook      Pivoting
      TCP    = lPiv .eq. 2  ! Threshold Complete  Pivoting.
      TSP    = lPiv .eq. 3  ! Threshold Symmetric Pivoting.
      kPiv(0)= 'PP'
      kPiv(1)= 'RP'
      kPiv(2)= 'CP'
      kPiv(3)= 'SP'

!     Initialize output parameters.

      inform = 0
      minlen = nelem + 2*(m + n)
      numl0  = 0
      lenL   = 0
      lenU   = 0
      lrow   = 0
      mersum = 0
      nUtri  = m
      nLtri  = 0
      ndens1 = 0
      ndens2 = 0
      nrank  = 0
      nsing  = 0
      jsing  = 0
      jumin  = 0

      Amax   = zero
      Lmax   = zero
      Umax   = zero
      DUmax  = zero
      DUmin  = zero
      Akmax  = zero

      if (m .gt. n) then
         mnkey  = '>'
      else if (m .eq. n) then
         mnkey  = '='
      else
         mnkey  = '<'
      end if

!     Float version of dimensions.

      dm     = m
      dn     = n
      delem  = nelem

!     Initialize workspace parameters.

      luparm(26) = 0             ! ncp
      if (lena .lt. minlen) go to 970

!     ------------------------------------------------------------------
!     Organize the  aij's  in  a, indc, indr.
!     lu1or1  deletes small entries, tests for illegal  i,j's,
!             and counts the nonzeros in each row and column.
!     lu1or2  reorders the elements of  A  by columns.
!     lu1or3  uses the column list to test for duplicate entries
!             (same indices  i,j).
!     lu1or4  constructs a row list from the column list.
!     ------------------------------------------------------------------
      call lu1or1( m   , n    , nelem, lena , small,
     &             a   , indc , indr , lenc , lenr,
     &             Amax, numnz, lerr , inform )

      if (nout. gt. 0  .and.  lprint .ge. 10) then
         densty = 100.0d+0 * delem / (dm * dn)
         write(nout, 1000) m, mnkey, n, nelem, Amax, densty
      end if
      if (inform .ne. 0) go to 930

      nelem  = numnz

      call lu1or2( n, nelem, lena, a, indc, indr, lenc, locc )
      call lu1or3( m, n, lena, indc, lenc, locc, ip,
     &             lerr, inform )

      if (inform .ne. 0) go to 940

      call lu1or4( m, n, nelem, lena,
     &             indc, indr, lenc, lenr, locc, locr )

!     ------------------------------------------------------------------
!     Set up lists of rows and columns with equal numbers of nonzeros,
!     using  indc(*)  as workspace.
!     ------------------------------------------------------------------
      call lu1pq1( m, n, lenr, ip, iploc, ipinv, indc(nelem + 1) )
      call lu1pq1( n, m, lenc, iq, iqloc, iqinv, indc(nelem + 1) )

!     ------------------------------------------------------------------
!     For TCP, allocate Ha, Hj, Hk at the end of a, indc, indr.
!     Then compute the factorization  A = L*U.
!     ------------------------------------------------------------------
      if (TPP .or. TSP) then
         lenH   = 1
         lena2  = lena
         locH   = lena
         lmaxr  = 1
      else if (TRP) then
         lenH   = 1             ! Dummy
         lena2  = lena  - m     ! Reduced length of      a
         locH   = lena          ! Dummy
         lmaxr  = lena2 + 1     ! Start of Amaxr      in a
      else if (TCP) then
         lenH   = n             ! Length of heap
         lena2  = lena  - lenH  ! Reduced length of      a, indc, indr
         locH   = lena2 + 1     ! Start of Ha, Hj, Hk in a, indc, indr
         lmaxr  = 1             ! Dummy
      end if

      call lu1fad( m     , n     , nelem , lena2 , luparm, parmlu,
     &             a     , indc  , indr  , ip    , iq    ,
     &             lenc  , lenr  , locc  , locr  ,
     &             iploc , iqloc , ipinv , iqinv , w     ,
     &             lenH  ,a(locH), indc(locH), indr(locH), a(lmaxr),
     &             inform, lenL  , lenU  , minlen, mersum,
     &             nUtri , nLtri , ndens1, ndens2, nrank ,
     &             Lmax  , Umax  , DUmax , DUmin , Akmax )

      luparm(16) = nrank
      luparm(23) = lenL
      if (inform .eq. 7) go to 970
      if (inform .eq. 9) go to 985
      if (inform .gt. 0) go to 980

      if ( keepLU ) then
!        ---------------------------------------------------------------
!        The LU factors are at the top of  a, indc, indr,
!        with the columns of  L  and the rows of  U  in the order
!
!        ( free )   ... ( u3 ) ( l3 ) ( u2 ) ( l2 ) ( u1 ) ( l1 ).
!
!        Starting with ( l1 ) and ( u1 ), move the rows of  U  to the
!        left and the columns of  L  to the right, giving
!
!        ( u1 ) ( u2 ) ( u3 ) ...   ( free )   ... ( l3 ) ( l2 ) ( l1 ).
!
!        Also, set  numl0 = the number of nonempty columns of L.
!        ---------------------------------------------------------------
         lu     = 0
         ll     = lena  + 1
         lm     = lena2 + 1
         ltopl  = ll - lenL - lenU
         lrow   = lenU

         do k = 1, nrank
            i       =   ip(k)
            lenUk   = - lenr(i)
            lenr(i) =   lenUk
            j       =   iq(k)
            lenLk   = - lenc(j) - 1
            if (lenLk .gt. 0) then
                numl0        = numl0 + 1
                iqloc(numl0) = lenLk
            end if

            if (lu + lenUk .lt. ltopl) then
!              =========================================================
!              There is room to move ( uk ).  Just right-shift ( lk ).
!              =========================================================
               do idummy = 1, lenLk
                  ll       = ll - 1
                  lm       = lm - 1
                  a(ll)    = a(lm)
                  indc(ll) = indc(lm)
                  indr(ll) = indr(lm)
               end do
            else
!              =========================================================
!              There is no room for ( uk ) yet.  We have to
!              right-shift the whole of the remaining LU file.
!              Note that ( lk ) ends up in the correct place.
!              =========================================================
               llsave = ll - lenLk
               nmove  = lm - ltopl

               do idummy = 1, nmove
                  ll       = ll - 1
                  lm       = lm - 1
                  a(ll)    = a(lm)
                  indc(ll) = indc(lm)
                  indr(ll) = indr(lm)
               end do

               ltopl  = ll
               ll     = llsave
               lm     = ll
            end if

!           ======================================================
!           Left-shift ( uk ).
!           ======================================================
            locr(i) = lu + 1
            l2      = lm - 1
            lm      = lm - lenUk

            do l = lm, l2
               lu       = lu + 1
               a(lu)    = a(l)
               indr(lu) = indr(l)
            end do
         end do

!        ---------------------------------------------------------------
!        Save the lengths of the nonempty columns of  L,
!        and initialize  locc(j)  for the LU update routines.
!        ---------------------------------------------------------------
         do k = 1, numl0
            lenc(k) = iqloc(k)
         end do

         do j = 1, n
            locc(j) = 0
         end do

!        ---------------------------------------------------------------
!        Test for singularity.
!        lu6chk  sets  nsing, jsing, jumin, Lmax, Umax, DUmax, DUmin
!        (including entries from the dense LU).
!        inform = 1  if there are singularities (nsing gt 0).
!        ---------------------------------------------------------------
         call lu6chk( 1, m, n, w, lena, luparm, parmlu,
     &                a, indc, indr, ip, iq,
     &                lenc, lenr, locc, locr, inform )
         nsing  = luparm(11)
         jsing  = luparm(12)
         jumin  = luparm(19)
         Lmax   = parmlu(11)
         Umax   = parmlu(12)
         DUmax  = parmlu(13)
         DUmin  = parmlu(14)

      else
!        ---------------------------------------------------------------
!        keepLU = 0.  L and U were not kept, just the diagonals of U.
!        lu1fac will probably be called again soon with keepLU = .true.
!        11 Mar 2001: lu6chk revised.  We can call it with keepLU = 0,
!                     but we want to keep Lmax, Umax from lu1fad.
!        05 May 2002: Allow for TCP with new lu1DCP.  Diag(U) starts
!                     below lena2, not lena.  Need lena2 in next line.
!        ---------------------------------------------------------------
         call lu6chk( 1, m, n, w, lena2, luparm, parmlu,
     &                a, indc, indr, ip, iq,
     &                lenc, lenr, locc, locr, inform )
         nsing  = luparm(11)
         jsing  = luparm(12)
         jumin  = luparm(19)
         DUmax  = parmlu(13)
         DUmin  = parmlu(14)
      end if

      go to 990

!     ------------
!     Error exits.
!     ------------
  930 inform = 3
      if (lprint .ge. 0) write(nout, 1300) lerr, indc(lerr), indr(lerr)
      go to 990

  940 inform = 4
      if (lprint .ge. 0) write(nout, 1400) lerr, indc(lerr), indr(lerr)
      go to 990

  970 inform = 7
      if (lprint .ge. 0) write(nout, 1700) lena, minlen
      go to 990

  980 inform = 8
      if (lprint .ge. 0) write(nout, 1800)
      go to 990

  985 inform = 9
      if (lprint .ge. 0) write(nout, 1900)

!     Store output parameters.

  990 nelem      = nelem0
      luparm(10) = inform
      luparm(11) = nsing
      luparm(12) = jsing
      luparm(13) = minlen
      luparm(15) = 0
      luparm(16) = nrank
      luparm(17) = ndens1
      luparm(18) = ndens2
      luparm(19) = jumin
      luparm(20) = numl0
      luparm(21) = lenL
      luparm(22) = lenU
      luparm(23) = lenL
      luparm(24) = lenU
      luparm(25) = lrow
      luparm(27) = mersum
      luparm(28) = nUtri
      luparm(29) = nLtri

      parmlu(10) = Amax
      parmlu(11) = Lmax
      parmlu(12) = Umax
      parmlu(13) = DUmax
      parmlu(14) = DUmin
      parmlu(15) = Akmax

      Agrwth = Akmax  / (Amax + 1.0d-20)
      Ugrwth = Umax   / (Amax + 1.0d-20)
      if ( TPP ) then
         growth = Ugrwth
      else ! TRP or TCP or TSP
         growth = Agrwth
      end if
      parmlu(16) = growth

!     ------------------------------------------------------------------
!     Print statistics for the LU factors.
!     ------------------------------------------------------------------
      ncp    = luparm(26)
      condU  = DUmax / max( DUmin, 1.0d-20 )
      dincr  = lenL + lenU - nelem
      dincr  = dincr * 100.0d+0 / max( delem, one )
      avgmer = mersum
      avgmer = avgmer / dm
      nbump  = m - nUtri - nLtri

      if (nout. gt. 0  .and.  lprint .ge. 10) then
         if ( TPP ) then
            write(nout, 1100) avgmer, lenL, lenL+lenU, ncp, dincr,
     &                        nUtri, lenU, Ltol, Umax, Ugrwth,
     &                        nLtri, ndens1, Lmax

         else
            write(nout, 1120) kPiv(lPiv), avgmer,
     &                        lenL, lenL+lenU, ncp, dincr,
     &                        nUtri, lenU, Ltol, Umax, Ugrwth,
     &                        nLtri, ndens1, Lmax, Akmax, Agrwth
         end if

         write(nout, 1200) nbump, ndens2, DUmax, DUmin, condU
      end if

      return

 1000 format(' m', i12, ' ', a, 'n', i12, '  Elems', i9,
     &       '  Amax', 1p, e10.1, '  Density', 0p, f7.2)
 1100 format(' Merit', 0p, f8.1, '  lenL', i9, '  L+U', i11,
     &       '  Cmpressns', i5, '  Incres', 0p, f8.2
     & /     ' Utri', i9, '  lenU', i9, '  Ltol', 1p, e10.2,
     &       '  Umax', e10.1, '  Ugrwth', e8.1
     & /     ' Ltri', i9, '  dense1', i7, '  Lmax', e10.2)
 1120 format(' Mer', a2, 0p, f8.1, '  lenL', i9, '  L+U', i11,
     &       '  Cmpressns', i5, '  Incres', 0p, f8.2
     & /     ' Utri', i9, '  lenU', i9, '  Ltol', 1p, e10.2,
     &       '  Umax', e10.1, '  Ugrwth', e8.1
     & /     ' Ltri', i9, '  dense1', i7, '  Lmax', e10.2,
     &       '  Akmax', e9.1, '  Agrwth', e8.1)
 1200 format(' bump', i9, '  dense2', i7, '  DUmax', 1p, e9.1,
     &       '  DUmin', e9.1, '  condU', e9.1)
 1300 format(/ ' lu1fac  error...  entry  a(', i8, ')  has an illegal',
     &         ' row or column index'
     &       //' indc, indr =', 2i8)
 1400 format(/ ' lu1fac  error...  entry  a(', i8, ')  has the same',
     &         ' indices as an earlier entry'
     &       //' indc, indr =', 2i8)
 1700 format(/ ' lu1fac  error...  insufficient storage'
     &       //' Increase  lena  from', i10, '  to at least', i10)
 1800 format(/ ' lu1fac  error...  fatal bug',
     &         '   (sorry --- this should never happen)')
 1900 format(/ ' lu1fac  error...  TSP used but',
     &         ' diagonal pivot could not be found')

      end ! subroutine lu1fac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1fad( m     , n     , nelem , lena  , luparm, parmlu,
     &                   a     , indc  , indr  , ip    , iq    ,
     &                   lenc  , lenr  , locc  , locr  ,
     &                   iploc , iqloc , ipinv , iqinv , w     ,
     &                   lenH  , Ha    , Hj    , Hk    , Amaxr ,
     &                   inform, lenL  , lenU  , minlen, mersum,
     &                   nUtri , nLtri , ndens1, ndens2, nrank ,
     &                   Lmax  , Umax  , DUmax , DUmin , Akmax )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), Amaxr(m), w(n)
      double precision   Ha(lenH)
      integer            indc(lena), indr(lena), ip(m), iq(n)
      integer            lenc(n)   , lenr(m)
      integer            locc(n)   , locr(m)
      integer            iploc(n)  , iqloc(m), ipinv(m), iqinv(n)
      integer            Hj(lenH)  , Hk(lenH)
      double precision   Lmax

!     ------------------------------------------------------------------
!     lu1fad  is a driver for the numerical phase of lu1fac.
!     At each stage it computes a column of  L  and a row of  U,
!     using a Markowitz criterion to select the pivot element,
!     subject to a stability criterion that bounds the elements of  L.
!
!     00 Jan 1986  Version documented in LUSOL paper:
!                  Gill, Murray, Saunders and Wright (1987),
!                  Maintaining LU factors of a general sparse matrix,
!                  Linear algebra and its applications 88/89, 239-270.
!
!     02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest
!                  element in each column is now kept at the start of
!                  the column, i.e. in position locc(j) of a and indc.
!                  This should speed up the Markowitz searches.
!                  To save time on highly triangular matrices, we wait
!                  until there are no further columns of length 1
!                  before setting and maintaining that property.
!
!     12 Apr 1989  ipinv and iqinv added (inverses of ip and iq)
!                  to save searching ip and iq for rows and columns
!                  altered in each elimination step.  (Used in lu1pq2)
!
!     19 Apr 1989  Code segmented to reduce its size.
!                  lu1gau does most of the Gaussian elimination work.
!                  lu1mar does just the Markowitz search.
!                  lu1mxc moves biggest elements to top of columns.
!                  lu1pen deals with pending fill-in in the row list.
!                  lu1pq2 updates the row and column permutations.
!
!     26 Apr 1989  maxtie replaced by maxcol, maxrow in the Markowitz
!                  search.  maxcol, maxrow change as density increases.
!
!     25 Oct 1993  keepLU implemented.
!
!     07 Feb 1994  Exit main loop early to finish off with a dense LU.
!                  densLU tells lu1fad whether to do it.
!     21 Dec 1994  Bug fixed.  nrank was wrong after the call to lu1ful.
!     12 Nov 1999  A parallel version of dcopy gave trouble in lu1ful
!                  during left-shift of dense matrix D within a(*).
!                  Fixed this unexpected problem here in lu1fad
!                  by making sure the first and second D don't overlap.
!
!     13 Sep 2000  TCP (Threshold Complete Pivoting) implemented.
!                  lu2max added
!                  (finds aijmax from biggest elems in each col).
!                  Utri, Ltri and Spars1 phases apply.
!                  No switch to Dense CP yet.  (Only TPP switches.)
!     14 Sep 2000  imax needed to remember row containing aijmax.
!     22 Sep 2000  For simplicity, lu1mxc always fixes
!                  all modified cols.
!                  (TPP spars2 used to fix just the first maxcol cols.)
!     08 Nov 2000: Speed up search for aijmax.
!                  Don't need to search all columns if the elimination
!                  didn't alter the col containing the current aijmax.
!     21 Nov 2000: lu1slk implemented for Utri phase with TCP
!                  to guard against deceptive triangular matrices.
!                  (Utri used to have aijtol >= 0.9999 to include
!                  slacks, but this allows other 1s to be accepted.)
!                  Utri now accepts slacks, but applies normal aijtol
!                  test to other pivots.
!     28 Nov 2000: TCP with empty cols must call lu1mxc and lu2max
!                  with ( lq1, n, ... ), not just ( 1, n, ... ).
!     23 Mar 2001: lu1fad bug with TCP.
!                  A col of length 1 might not be accepted as a pivot.
!                  Later it appears in a pivot row and temporarily
!                  has length 0 (when pivot row is removed
!                  but before the column is filled in).  If it is the
!                  last column in storage, the preceding col also thinks
!                  it is "last".  Trouble arises when the preceding col
!                  needs fill-in -- it overlaps the real "last" column.
!                  (Very rarely, same trouble might have happened if
!                  the drop tolerance caused columns to have length 0.)
!
!                  Introduced ilast to record the last row in row file,
!                             jlast to record the last col in col file.
!                  lu1rec returns ilast = indr(lrow + 1)
!                              or jlast = indc(lcol + 1).
!                  (Should be an output parameter, but didn't want to
!                  alter lu1rec's parameter list.)
!                  lu1rec also treats empty rows or cols safely.
!                  (Doesn't eliminate them!)
!
!     26 Apr 2002: Heap routines added for TCP.
!                  lu2max no longer needed.
!                  imax, jmax used only for printing.
!     01 May 2002: lu1DCP implemented (dense complete pivoting).
!                  Both TPP and TCP now switch to dense LU
!                  when density exceeds dens2.
!     06 May 2002: In dense mode, store diag(U) in natural order.
!     09 May 2002: lu1mCP implemented (Markowitz TCP via heap).
!     11 Jun 2002: lu1mRP implemented (Markowitz TRP).
!     28 Jun 2002: Fixed call to lu1mxr.
!     14 Dec 2002: lu1mSP implemented (Markowitz TSP).
!     15 Dec 2002: Both TPP and TSP can grab cols of length 1
!                  during Utri.
!     19 Dec 2004: Hdelete(...) has new input argument Hlenin.
!     26 Mar 2006: lu1fad returns nrank  = min( mrank, nrank )
!                  and ignores nsing from lu1ful
!     26 Mar 2006: Allow for empty columns before calling Hbuild.
!
!     Systems Optimization Laboratory, Stanford University.
!     ------------------------------------------------------------------

      logical            Utri, Ltri, spars1, spars2, dense
      logical            densLU, keepLU
      logical            TCP, TPP, TRP, TSP
      double precision   Lij, Ltol, small
      integer            Hlen, Hlenin, hops, h, lPiv

      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one  = 1.0d+0 )

!     ------------------------------------------------------------------
!     Local variables
!     ---------------
!
!     lcol   is the length of the column file.  It points to the last
!            nonzero in the column list.
!     lrow   is the analogous quantity for the row file.
!     lfile  is the file length (lcol or lrow) after the most recent
!            compression of the column list or row list.
!     nrowd  and  ncold  are the number of rows and columns in the
!            matrix defined by the pivot column and row.  They are the
!            dimensions of the submatrix D being altered at this stage.
!     melim  and  nelim  are the number of rows and columns in the
!            same matrix D, excluding the pivot column and row.
!     mleft  and  nleft  are the number of rows and columns
!            still left to be factored.
!     nzchng is the increase in nonzeros in the matrix that remains
!            to be factored after the current elimination
!            (usually negative).
!     nzleft is the number of nonzeros still left to be factored.
!     nspare is the space we leave at the end of the last row or
!            column whenever a row or column is being moved to the end
!            of its file.  nspare = 1 or 2 might help reduce the
!            number of file compressions when storage is tight.
!
!     The row and column ordering permutes A into the form
!
!                        ------------------------
!                         \                     |
!                          \         U1         |
!                           \                   |
!                            --------------------
!                            |\
!                            | \
!                            |  \
!            P A Q   =       |   \
!                            |    \
!                            |     --------------
!                            |     |            |
!                            |     |            |
!                            | L1  |     A2     |
!                            |     |            |
!                            |     |            |
!                            --------------------
!
!     where the block A2 is factored as  A2 = L2 U2.
!     The phases of the factorization are as follows.
!
!     Utri   is true when U1 is being determined.
!            Any column of length 1 is accepted immediately (if TPP).
!
!     Ltri   is true when L1 is being determined.
!            lu1mar exits as soon as an acceptable pivot is found
!            in a row of length 1.
!
!     spars1 is true while the density of the (modified) A2 is less
!            than the parameter dens1 = parmlu(7) = 0.3 say.
!            lu1mar searches maxcol columns and maxrow rows,
!            where  maxcol = luparm(3),  maxrow = maxcol - 1.
!            lu1mxc is used to keep the biggest element at the top
!            of all remaining columns.
!
!     spars2 is true while the density of the modified A2 is less
!            than the parameter dens2 = parmlu(8) = 0.6 say.
!            lu1mar searches maxcol columns and no rows.
!            lu1mxc could fix up only the first maxcol cols (with TPP).
!            22 Sep 2000:  For simplicity, lu1mxc fixes all
!                          modified cols.
!
!     dense  is true once the density of A2 reaches dens2.
!            lu1mar searches only 1 column (the shortest).
!            lu1mxc could fix up only the first column (with TPP).
!            22 Sep 2000:  For simplicity, lu1mxc fixes all
!                          modified cols.
!
!     ------------------------------------------------------------------
!     To eliminate timings, comment out all lines containing "time".
!     ------------------------------------------------------------------

!     integer            eltime, mktime

!     call timer ( 'start ', 3 )
!     ntime  = (n / 4.0)


      nout   = luparm(1)
      lprint = luparm(2)
      maxcol = luparm(3)
      lPiv   = luparm(6)
      keepLU = luparm(8) .ne. 0

      TPP    = lPiv .eq. 0  ! Threshold Partial   Pivoting (normal).
      TRP    = lPiv .eq. 1  ! Threshold Rook      Pivoting
      TCP    = lPiv .eq. 2  ! Threshold Complete  Pivoting.
      TSP    = lPiv .eq. 3  ! Threshold Symmetric Pivoting.

      densLU = .false.
      maxrow = maxcol - 1
      ilast  = m                 ! Assume row m is last in the row file.
      jlast  = n                 ! Assume col n is last in the col file.
      lfile  = nelem
      lrow   = nelem
      lcol   = nelem
      minmn  = min( m, n )
      maxmn  = max( m, n )
      nzleft = nelem
      nspare = 1

      if ( keepLU ) then
         lu1    = lena   + 1
      else
!        Store only the diagonals of U in the top of memory.
         ldiagU = lena   - n
         lu1    = ldiagU + 1
      end if

      Ltol   = parmlu(1)
      small  = parmlu(3)
      uspace = parmlu(6)
      dens1  = parmlu(7)
      dens2  = parmlu(8)
      Utri   = .true.
      Ltri   = .false.
      spars1 = .false.
      spars2 = .false.
      dense  = .false.

!     Check parameters.

      Ltol   = max( Ltol, 1.0001d+0 )
      dens1  = min( dens1, dens2 )

!     Initialize output parameters.
!     lenL, lenU, minlen, mersum, nUtri, nLtri, ndens1, ndens2, nrank
!     are already initialized by lu1fac.

      Lmax   = zero
      Umax   = zero
      DUmax  = zero
      DUmin  = 1.0d+20
      if (nelem .eq. 0) Dumin = zero
      Akmax  = zero
      hops   = 0

      ! More initialization.

      if (TPP .or. TSP) then ! Don't worry yet about lu1mxc.
         aijmax = zero
         aijtol = zero
         Hlen   = 1

      else ! TRP or TCP
         ! Move biggest element to top of each column.
         ! Set w(*) to mark slack columns (unit vectors).

         call lu1mxc( 1, n, iq, a, indc, lenc, locc )
         call lu1slk( m, n, lena, iq, iqloc, a, locc, w )
      end if

      if (TRP) then
         ! Find biggest element in each row.

         call lu1mxr( 1, m, ip, Amaxr,
     &                a, indc, lenc, locc, indr, lenr, locr )
      end if

      if (TCP) then
         ! Set Ha(1:Hlen) = biggest element in each column,
         !     Hj(1:Hlen) = corresponding column indices.

         Hlen  = 0
         do kk = 1, n
            Hlen     = Hlen + 1
            j        = iq(kk)
            if (lenc(j) .gt. 0) then
               lc    = locc(j)
               amax  = abs( a(lc) )
            else
               amax  = zero
            end if
            Ha(Hlen) = amax
            Hj(Hlen) = j
            Hk(j)    = Hlen
         end do

         ! Build the heap, creating new Ha, Hj and setting Hk(1:Hlen).

         call Hbuild( Ha, Hj, Hk, Hlen, Hlen, hops )
      end if

!     ------------------------------------------------------------------
!     Start of main loop.
!     ------------------------------------------------------------------
      mleft  = m + 1
      nleft  = n + 1

      do 800 nrowu = 1, minmn

!        mktime = (nrowu / ntime) + 4
!        eltime = (nrowu / ntime) + 9
         mleft  = mleft - 1
         nleft  = nleft - 1

!        Bail out if there are no nonzero rows left.

         if (iploc(1) .gt. m) go to 900

         ! For TCP, the largest Aij is at the top of the heap.

         if ( TCP ) then
            aijmax = Ha(1)      ! Marvelously easy !
            Akmax  = max( Akmax, aijmax )
            aijtol = aijmax / Ltol
         end if

!        ===============================================================
!        Find a suitable pivot element.
!        ===============================================================

         if ( Utri ) then
!           ------------------------------------------------------------
!           So far all columns have had length 1.
!           We are still looking for the (backward) triangular part of A
!           that forms the first rows and columns of U.
!           ------------------------------------------------------------

            lq1    = iqloc(1)
            lq2    = n
            if (m   .gt.   1) lq2 = iqloc(2) - 1

            if (lq1 .le. lq2) then  ! There are more cols of length 1.
               if (TPP .or. TSP) then
                  jbest  = iq(lq1)  ! Grab the first one.

               else ! TRP or TCP    ! Scan all columns of length 1.
                  jbest  = 0

                  do lq = lq1, lq2
                     j      = iq(lq)
                     if (w(j) .gt. zero) then ! Accept a slack
                        jbest  = j
                        go to 250
                     end if

                     lc     = locc(j)
                     amax   = abs( a(lc) )
                     if (TRP) then
                        i      = indc(lc)
                        aijtol = Amaxr(i) / Ltol
                     end if

                     if (amax .ge. aijtol) then
                        jbest  = j
                        go to 250
                     end if
                  end do
               end if

  250          if (jbest .gt. 0) then
                  lc     = locc(jbest)
                  ibest  = indc(lc)
                  mbest  = 0
                  go to 300
               end if
            end if

!           This is the end of the U triangle.
!           We will not return to this part of the code.
!           TPP and TSP call lu1mxc for the first time
!           (to move biggest element to top of each column).

            if (lprint .ge. 50) then
               write(nout, 1100) 'Utri ended.  spars1 = true'
            end if
            Utri   = .false.
            Ltri   = .true.
            spars1 = .true.
            nUtri  =  nrowu - 1
            if (TPP .or. TSP) then
               call lu1mxc( lq1, n, iq, a, indc, lenc, locc )
            end if
         end if

         if ( spars1 ) then
!           ------------------------------------------------------------
!           Perform a Markowitz search.
!           Search cols of length 1, then rows of length 1,
!           then   cols of length 2, then rows of length 2, etc.
!           ------------------------------------------------------------
!           call timer ( 'start ', mktime )

          ! if (TPP) then ! 12 Jun 2002: Next line disables lu1mCP below
            if (TPP .or. TCP) then
               call lu1mar( m    , n     , lena  , maxmn,
     &                      TCP  , aijtol, Ltol  , maxcol, maxrow,
     &                      ibest, jbest , mbest ,
     &                      a    , indc  , indr  , ip   , iq,
     &                      lenc , lenr  , locc  , locr ,
     &                      iploc, iqloc )

            else if (TRP) then
               call lu1mRP( m    , n     , lena  , maxmn,
     &                      Ltol , maxcol, maxrow,
     &                      ibest, jbest , mbest ,
     &                      a    , indc  , indr  , ip   , iq,
     &                      lenc , lenr  , locc  , locr ,
     &                      iploc, iqloc , Amaxr )

!           else if (TCP) then ! Disabled by test above
!              call lu1mCP( m    , n     , lena  , aijtol,
!    &              ibest, jbest , mbest ,
!    &              a    , indc  , indr  ,
!    &              lenc , lenr  , locc  ,
!    &              Hlen , Ha    , Hj    )

            else if (TSP) then
               call lu1mSP( m    , n     , lena  , maxmn,
     &                      Ltol , maxcol,
     &                      ibest, jbest , mbest ,
     &                      a    , indc  , iq    , locc , iqloc )
               if (ibest .eq. 0) go to 990
            end if

!           call timer ( 'finish', mktime )

            if ( Ltri ) then

!              So far all rows have had length 1.
!              We are still looking for the (forward) triangle of A
!              that forms the first rows and columns of L.

               if (mbest .gt. 0) then
                   Ltri   = .false.
                   nLtri  =  nrowu - 1 - nUtri
                   if (lprint .ge. 50) then
                      write(nout, 1100) 'Ltri ended.'
                   end if
               end if

            else

!              See if what's left is as dense as dens1.

               if (nzleft  .ge.  (dens1 * mleft) * nleft) then
                   spars1 = .false.
                   spars2 = .true.
                   ndens1 =  nleft
                   maxrow =  0
                   if (lprint .ge. 50) then
                      write(nout, 1100) 'spars1 ended.  spars2 = true'
                   end if
               end if
            end if

         else if ( spars2 .or. dense ) then
!           ------------------------------------------------------------
!           Perform a restricted Markowitz search,
!           looking at only the first maxcol columns.  (maxrow = 0.)
!           ------------------------------------------------------------
!           call timer ( 'start ', mktime )

          ! if (TPP) then ! 12 Jun 2002: Next line disables lu1mCP below
            if (TPP .or. TCP) then
               call lu1mar( m    , n     , lena  , maxmn,
     &              TCP  , aijtol, Ltol  , maxcol, maxrow,
     &              ibest, jbest , mbest ,
     &              a    , indc  , indr  , ip   , iq,
     &              lenc , lenr  , locc  , locr ,
     &              iploc, iqloc )

            else if (TRP) then
               call lu1mRP( m    , n     , lena  , maxmn,
     &              Ltol , maxcol, maxrow,
     &              ibest, jbest , mbest ,
     &              a    , indc  , indr  , ip   , iq,
     &              lenc , lenr  , locc  , locr ,
     &              iploc, iqloc , Amaxr )

!           else if (TCP) then ! Disabled by test above
!              call lu1mCP( m    , n     , lena  , aijtol,
!    &              ibest, jbest , mbest ,
!    &              a    , indc  , indr  ,
!    &              lenc , lenr  , locc  ,
!    &              Hlen , Ha    , Hj    )

            else if (TSP) then
               call lu1mSP( m    , n     , lena  , maxmn,
     &                      Ltol , maxcol,
     &                      ibest, jbest , mbest ,
     &                      a    , indc  , iq    , locc , iqloc )
               if (ibest .eq. 0) go to 985
            end if

!           call timer ( 'finish', mktime )

!           See if what's left is as dense as dens2.

            if ( spars2 ) then
               if (nzleft  .ge.  (dens2 * mleft) * nleft) then
                   spars2 = .false.
                   dense  = .true.
                   ndens2 =  nleft
                   maxcol =  1
                   if (lprint .ge. 50) then
                      write(nout, 1100) 'spars2 ended.  dense = true'
                   end if
               end if
            end if
         end if

!        ---------------------------------------------------------------
!        See if we can finish quickly.
!        ---------------------------------------------------------------
         if ( dense  ) then
            lenD   = mleft * nleft
            nfree  = lu1 - 1

            if (nfree .ge. 2 * lenD) then

!              There is room to treat the remaining matrix as
!              a dense matrix D.
!              We may have to compress the column file first.
!              12 Nov 1999: D used to be put at the
!                           beginning of free storage (lD = lcol + 1).
!                           Now put it at the end     (lD = lu1 - lenD)
!                           so the left-shift in lu1ful will not
!                           involve overlapping storage
!                           (fatal with parallel dcopy).
!
               densLU = .true.
               ndens2 = nleft
               lD     = lu1 - lenD
               if (lcol .ge. lD) then
                  call lu1rec( n, .true., luparm, lcol,
     &                         lena, a, indc, lenc, locc )
                  lfile  = lcol
                  jlast  = indc(lcol + 1)
               end if

               go to 900
            end if
         end if

!        ===============================================================
!        The best  aij  has been found.
!        The pivot row  ibest  and the pivot column  jbest
!        Define a dense matrix  D  of size  nrowd  by  ncold.
!        ===============================================================
  300    ncold  = lenr(ibest)
         nrowd  = lenc(jbest)
         melim  = nrowd  - 1
         nelim  = ncold  - 1
         mersum = mersum + mbest
         lenL   = lenL   + melim
         lenU   = lenU   + ncold
         if (lprint .ge. 50) then
            if (nrowu .eq. 1) then
               write(nout, 1100) 'lu1fad debug:'
            end if
            if ( TPP .or. TRP .or. TSP ) then
               write(nout, 1200) nrowu, ibest, jbest, nrowd, ncold
            else ! TCP
               jmax   = Hj(1)
               imax   = indc(locc(jmax))
               write(nout, 1200) nrowu, ibest, jbest, nrowd, ncold,
     &                           imax , jmax , aijmax
            end if
         end if

!        ===============================================================
!        Allocate storage for the next column of  L  and next row of  U.
!        Initially the top of a, indc, indr are used as follows:
!
!                   ncold       melim       ncold        melim
!
!        a      |...........|...........|ujbest..ujn|li1......lim|
!
!        indc   |...........|  lenr(i)  |  lenc(j)  |  markl(i)  |
!
!        indr   |...........| iqloc(i)  |  jfill(j) |  ifill(i)  |
!
!              ^           ^             ^           ^            ^
!              lfree   lsave             lu1         ll1          oldlu1
!
!        Later the correct indices are inserted:
!
!        indc   |           |           |           |i1........im|
!
!        indr   |           |           |jbest....jn|ibest..ibest|
!
!        ===============================================================
         if ( keepLU ) then
!           relax
         else
!           Always point to the top spot.
!           Only the current column of L and row of U will
!           take up space, overwriting the previous ones.
            lu1    = ldiagU + 1
         end if
         ll1    = lu1   - melim
         lu1    = ll1   - ncold
         lsave  = lu1   - nrowd
         lfree  = lsave - ncold

!        Make sure the column file has room.
!        Also force a compression if its length exceeds a certain limit.

         limit  = uspace * lfile  +  m  +  n  +  1000
         minfre = ncold  + melim
         nfree  = lfree  - lcol
         if (nfree .lt. minfre  .or.  lcol .gt. limit) then
            call lu1rec( n, .true., luparm, lcol,
     &                   lena, a, indc, lenc, locc )
            lfile  = lcol
            jlast  = indc(lcol + 1)
            nfree  = lfree - lcol
            if (nfree .lt. minfre) go to 970
         end if

!        Make sure the row file has room.

         minfre = melim + ncold
         nfree  = lfree - lrow
         if (nfree .lt. minfre  .or.  lrow .gt. limit) then
            call lu1rec( m, .false., luparm, lrow,
     &                   lena, a, indr, lenr, locr )
            lfile  = lrow
            ilast  = indr(lrow + 1)
            nfree  = lfree - lrow
            if (nfree .lt. minfre) go to 970
         end if

!        ===============================================================
!        Move the pivot element to the front of its row
!        and to the top of its column.
!        ===============================================================
         lpivr  = locr(ibest)
         lpivr1 = lpivr + 1
         lpivr2 = lpivr + nelim

         do 330 l = lpivr, lpivr2
            if (indr(l) .eq. jbest) go to 335
  330    continue

  335    indr(l)     = indr(lpivr)
         indr(lpivr) = jbest

         lpivc  = locc(jbest)
         lpivc1 = lpivc + 1
         lpivc2 = lpivc + melim

         do 340 l = lpivc, lpivc2
            if (indc(l) .eq. ibest) go to 345
  340    continue

  345    indc(l)     = indc(lpivc)
         indc(lpivc) = ibest
         abest       = a(l)
         a(l)        = a(lpivc)
         a(lpivc)    = abest

         if ( keepLU ) then
!           relax
         else
!           Store just the diagonal of U, in natural order.
!!!         a(ldiagU + nrowu) = abest ! This was in pivot order.
            a(ldiagU + jbest) = abest
         end if

         !==============================================================
         ! Delete pivot col from heap.
         ! Hk tells us where it is in the heap.
         !==============================================================
         if ( TCP ) then
            kbest  = Hk(jbest)
            Hlenin = Hlen
            call Hdelete( Ha, Hj, Hk, Hlenin, Hlen, n, kbest, h )
            hops   = hops + h
         end if

!        ===============================================================
!        Delete the pivot row from the column file
!        and store it as the next row of  U.
!        set  indr(lu) = 0     to initialize jfill ptrs on columns of D,
!             indc(lu) = lenj  to save the original column lengths.
!        ===============================================================
         a(lu1)    = abest
         indr(lu1) = jbest
         indc(lu1) = nrowd
         lu        = lu1

         diag      = abs( abest )
         Umax      = max(  Umax, diag )
         DUmax     = max( DUmax, diag )
         DUmin     = min( DUmin, diag )

         do 360 lr   = lpivr1, lpivr2
            lu       = lu + 1
            j        = indr(lr)
            lenj     = lenc(j)
            lenc(j)  = lenj - 1
            lc1      = locc(j)
            last     = lc1 + lenc(j)

            do 350 l = lc1, last
               if (indc(l) .eq. ibest) go to 355
  350       continue

  355       a(lu)      = a(l)
            indr(lu)   = 0
            indc(lu)   = lenj
            Umax       = max( Umax, abs( a(lu) ) )
            a(l)       = a(last)
            indc(l)    = indc(last)
            indc(last) = 0       ! Free entry
!???        if (j .eq. jlast) lcol = lcol - 1
  360    continue

!        ===============================================================
!        Delete the pivot column from the row file
!        and store the nonzeros of the next column of  L.
!        Set  indc(ll) = 0     to initialize markl(*) markers,
!             indr(ll) = 0     to initialize ifill(*) row fill-in cntrs,
!             indc(ls) = leni  to save the original row lengths,
!             indr(ls) = iqloc(i)    to save parts of  iqloc(*),
!             iqloc(i) = lsave - ls  to point to the nonzeros of  L
!                      = -1, -2, -3, ... in mark(*).
!        ===============================================================
         indc(lsave) = ncold
         if (melim .eq. 0) go to 700

         ll     = ll1 - 1
         ls     = lsave
         abest  = one / abest

         do 390 lc   = lpivc1, lpivc2
            ll       = ll + 1
            ls       = ls + 1
            i        = indc(lc)
            leni     = lenr(i)
            lenr(i)  = leni - 1
            lr1      = locr(i)
            last     = lr1 + lenr(i)

            do 380 l = lr1, last
               if (indr(l) .eq. jbest) go to 385
  380       continue

  385       indr(l)    = indr(last)
            indr(last) = 0       ! Free entry
!???        if (i .eq. ilast) lrow = lrow - 1

            a(ll)      = - a(lc) * abest
            Lij        = abs( a(ll) )
            Lmax       = max( Lmax, Lij )
!!!!! DEBUG
!           if (Lij .gt. Ltol) then
!              write( *  ,*) ' Big Lij!!!', nrowu
!              write(nout,*) ' Big Lij!!!', nrowu
!           end if

            indc(ll)   = 0
            indr(ll)   = 0
            indc(ls)   = leni
            indr(ls)   = iqloc(i)
            iqloc(i)   = lsave - ls
  390    continue

!        ===============================================================
!        Do the Gaussian elimination.
!        This involves adding a multiple of the pivot column
!        to all other columns in the pivot row.
!
!        Sometimes more than one call to lu1gau is needed to allow
!        compression of the column file.
!        lfirst  says which column the elimination should start with.
!        minfre  is a bound on the storage needed for any one column.
!        lu      points to off-diagonals of u.
!        nfill   keeps track of pending fill-in in the row file.
!        ===============================================================
         if (nelim .eq. 0) go to 700
         lfirst = lpivr1
         minfre = mleft + nspare
         lu     = 1
         nfill  = 0

! 400    call timer ( 'start ', eltime )
  400    call lu1gau( m     , melim , ncold , nspare, small ,
     &                lpivc1, lpivc2, lfirst, lpivr2, lfree , minfre,
     &                ilast , jlast , lrow  , lcol  , lu    , nfill ,
     &                a     , indc  , indr  ,
     &                lenc  , lenr  , locc  , locr  ,
     &                iqloc , a(ll1), indc(ll1),
     &                        a(lu1), indr(ll1), indr(lu1) )
!        call timer ( 'finish', eltime )

         if (lfirst .gt. 0) then

!           The elimination was interrupted.
!           Compress the column file and try again.
!           lfirst, lu and nfill have appropriate new values.

            call lu1rec( n, .true., luparm, lcol,
     &                   lena, a, indc, lenc, locc )
            lfile  = lcol
            jlast  = indc(lcol + 1)
            lpivc  = locc(jbest)
            lpivc1 = lpivc + 1
            lpivc2 = lpivc + melim
            nfree  = lfree - lcol
            if (nfree .lt. minfre) go to 970
            go to 400
         end if

!        ===============================================================
!        The column file has been fully updated.
!        Deal with any pending fill-in in the row file.
!        ===============================================================
         if (nfill .gt. 0) then

!           Compress the row file if necessary.
!           lu1gau has set nfill to be the number of pending fill-ins
!           plus the current length of any rows that need to be moved.

            minfre = nfill
            nfree  = lfree - lrow
            if (nfree .lt. minfre) then
               call lu1rec( m, .false., luparm, lrow,
     &                      lena, a, indr, lenr, locr )
               lfile  = lrow
               ilast  = indr(lrow + 1)
               lpivr  = locr(ibest)
               lpivr1 = lpivr + 1
               lpivr2 = lpivr + nelim
               nfree  = lfree - lrow
               if (nfree .lt. minfre) go to 970
            end if

!           Move rows that have pending fill-in to end of the row file.
!           Then insert the fill-in.

            call lu1pen( m     , melim , ncold , nspare, ilast,
     &                   lpivc1, lpivc2, lpivr1, lpivr2, lrow ,
     &                   lenc  , lenr  , locc  , locr  ,
     &                   indc  , indr  , indr(ll1), indr(lu1) )
         end if

!        ===============================================================
!        Restore the saved values of  iqloc.
!        Insert the correct indices for the col of L and the row of U.
!        ===============================================================
  700    lenr(ibest) = 0
         lenc(jbest) = 0

         ll          = ll1 - 1
         ls          = lsave

         do 710  lc  = lpivc1, lpivc2
            ll       = ll + 1
            ls       = ls + 1
            i        = indc(lc)
            iqloc(i) = indr(ls)
            indc(ll) = i
            indr(ll) = ibest
  710    continue

         lu          = lu1 - 1

         do 720  lr  = lpivr, lpivr2
            lu       = lu + 1
            indr(lu) = indr(lr)
  720    continue

!        ===============================================================
!        Free the space occupied by the pivot row
!        and update the column permutation.
!        Then free the space occupied by the pivot column
!        and update the row permutation.
!
!        nzchng is found in both calls to lu1pq2, but we use it only
!        after the second.
!        ===============================================================
         call lu1pq2( ncold, nzchng,
     &                indr(lpivr), indc( lu1 ), lenc, iqloc, iq, iqinv )

         call lu1pq2( nrowd, nzchng,
     &                indc(lpivc), indc(lsave), lenr, iploc, ip, ipinv )

         nzleft = nzleft + nzchng

!        ===============================================================
!        lu1mxr resets Amaxr(i) in each modified row i.
!        lu1mxc moves the largest aij to the top of each modified col j.
!        28 Jun 2002: Note that cols of L have an implicit diag of 1.0,
!                     so lu1mxr is called with ll1, not ll1+1, whereas
!                        lu1mxc is called with          lu1+1.
!        ===============================================================
         if (Utri .and. TPP) then
            ! Relax -- we're not keeping big elements at the top yet.

         else
            if (TRP  .and.  melim .gt. 0) then
               call lu1mxr( ll1, ll, indc, Amaxr,
     &                      a, indc, lenc, locc, indr, lenr, locr )
            end if

            if (nelim .gt. 0) then
               call lu1mxc( lu1+1, lu, indr, a, indc, lenc, locc )

               if (TCP) then ! Update modified columns in heap
                  do kk = lu1+1, lu
                     j    = indr(kk)
                     k    = Hk(j)
                     v    = abs( a(locc(j)) ) ! Biggest aij in column j
                     call Hchange( Ha, Hj, Hk, Hlen, n, k, v, j, h )
                     hops = hops + h
                  end do
               end if
            end if
         end if

!        ===============================================================
!        Negate lengths of pivot row and column so they will be
!        eliminated during compressions.
!        ===============================================================
         lenr(ibest) = - ncold
         lenc(jbest) = - nrowd

!        Test for fatal bug: row or column lists overwriting L and U.

         if (lrow .gt. lsave) go to 980
         if (lcol .gt. lsave) go to 980

!        Reset the file lengths if pivot row or col was at the end.

         if (ibest .eq. ilast) then
            lrow = locr(ibest)
         end if

         if (jbest .eq. jlast) then
            lcol = locc(jbest)
         end if
  800 continue

!     ------------------------------------------------------------------
!     End of main loop.
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     Normal exit.
!     Move empty rows and cols to the end of ip, iq.
!     Then finish with a dense LU if necessary.
!     ------------------------------------------------------------------
  900 inform = 0
      call lu1pq3( m, lenr, ip, ipinv, mrank )
      call lu1pq3( n, lenc, iq, iqinv, nrank )
      nrank  = min( mrank, nrank )

      if ( densLU ) then
!        call timer ( 'start ', 17 )
         call lu1ful( m     , n    , lena , lenD , lu1 , TPP,
     &                mleft , nleft, nrank, nrowu,
     &                lenL  , lenU , nsing,
     &                keepLU, small,
     &                a     , a(lD), indc , indr , ip  , iq,
     &                lenc  , lenr , locc , ipinv, locr )
!***     21 Dec 1994: Bug in next line.
!***     nrank  = nrank - nsing.  Changed to next line:
!***     nrank  = minmn - nsing

!***     26 Mar 2006: Previous line caused bug with m<n and nsing>0.
!        Don't mess with nrank any more.  Let end of lu1fac handle it.
!        call timer ( 'finish', 17 )
      end if

      minlen = lenL  +  lenU  +  2*(m + n)
      go to 990

!     Not enough space free after a compress.
!     Set  minlen  to an estimate of the necessary value of  lena.

  970 inform = 7
      minlen = lena  +  lfile  +  2*(m + n)
      go to 990

!     Fatal error.  This will never happen!
!    (Famous last words.)

  980 inform = 8
      go to 990

!     Fatal error with TSP.  Diagonal pivot not found.

  985 inform = 9

!     Exit.

  990 continue
!     call timer ( 'finish', 3 )
      return

 1100 format(/ 1x, a)
 1200 format(  ' nrowu', i7,
     &       '   i,jbest', 2i7,
     &       '   nrowd,ncold', 2i6,
     &       '   i,jmax', 2i7,
     &       '   aijmax', 1p, e10.2)

      end ! subroutine lu1fad

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1gau( m     , melim , ncold , nspare, small ,
     &                   lpivc1, lpivc2, lfirst, lpivr2, lfree , minfre,
     &                   ilast , jlast , lrow  , lcol  , lu    , nfill ,
     &                   a     , indc  , indr  ,
     &                   lenc  , lenr  , locc  , locr  ,
     &                   mark  , al    , markl ,
     &                           au    , ifill , jfill )

      implicit           double precision (a-h,o-z)
      double precision   a(*)        , al(melim)   , au(ncold)
      integer            indc(*)     , indr(*)     , lenc(*)  , lenr(*),
     &                   mark(*)     , markl(melim),
     &                   ifill(melim), jfill(ncold)
      integer            locc(*)     , locr(*)

!     ------------------------------------------------------------------
!     lu1gau does most of the work for each step of
!     Gaussian elimination.
!     A multiple of the pivot column is added to each other column j
!     in the pivot row.  The column list is fully updated.
!     The row list is updated if there is room, but some fill-ins may
!     remain, as indicated by ifill and jfill.
!
!
!  Input:
!     ilast    is the row    at the end of the row    list.
!     jlast    is the column at the end of the column list.
!     lfirst   is the first column to be processed.
!     lu + 1   is the corresponding element of U in au(*).
!     nfill    keeps track of pending fill-in.
!     a(*)     contains the nonzeros for each column j.
!     indc(*)  contains the row indices for each column j.
!     al(*)    contains the new column of L.  A multiple of it is
!              used to modify each column.
!     mark(*)  has been set to -1, -2, -3, ... in the rows
!              corresponding to nonzero 1, 2, 3, ... of the col of L.
!     au(*)    contains the new row of U.  Each nonzero gives the
!              required multiple of the column of L.
!
!  Workspace:
!     markl(*) marks the nonzeros of L actually used.
!              (A different mark, namely j, is used for each column.)
!
!  Output:
!     ilast     New last row    in the row    list.
!     jlast     New last column in the column list.
!     lfirst    = 0 if all columns were completed,
!               > 0 otherwise.
!     lu        returns the position of the last nonzero of U
!               actually used, in case we come back in again.
!     nfill     keeps track of the total extra space needed in the
!               row file.
!     ifill(ll) counts pending fill-in for rows involved in the new
!               column of L.
!     jfill(lu) marks the first pending fill-in stored in columns
!               involved in the new row of U.
!
!     16 Apr 1989: First version of lu1gau.
!     23 Apr 1989: lfirst, lu, nfill are now input and output
!                  to allow re-entry if elimination is interrupted.
!     23 Mar 2001: Introduced ilast, jlast.
!     27 Mar 2001: Allow fill-in "in situ" if there is already room
!                  up to but NOT INCLUDING the end of the
!                  row or column file.
!                  Seems safe way to avoid overwriting empty rows/cols
!                  at the end.  (May not be needed though, now that we
!                  have ilast and jlast.)
!     ------------------------------------------------------------------

      logical            atend

      do 600 lr = lfirst, lpivr2
         j      = indr(lr)
         lenj   = lenc(j)
         nfree  = lfree - lcol
         if (nfree .lt. minfre) go to 900

!        ---------------------------------------------------------------
!        Inner loop to modify existing nonzeros in column  j.
!        Loop 440 performs most of the arithmetic involved in the
!        whole LU factorization.
!        ndone counts how many multipliers were used.
!        ndrop counts how many modified nonzeros are negligibly small.
!        ---------------------------------------------------------------
         lu     = lu + 1
         uj     = au(lu)
         lc1    = locc(j)
         lc2    = lc1 + lenj - 1
         atend  = j .eq. jlast
         ndone  = 0
         if (lenj .eq. 0) go to 500

         ndrop  = 0

         do 440 l = lc1, lc2
            i        =   indc(l)
            ll       = - mark(i)
            if (ll .gt. 0) then
               ndone     = ndone + 1
               markl(ll) = j
               a(l)      = a(l)  +  al(ll) * uj
               if (abs( a(l) ) .le. small) then
                  ndrop  = ndrop + 1
               end if
            end if
  440    continue

!        ---------------------------------------------------------------
!        Remove any negligible modified nonzeros from both
!        the column file and the row file.
!        ---------------------------------------------------------------
         if (ndrop .eq. 0) go to 500
         k      = lc1

         do 480 l = lc1, lc2
            i        = indc(l)
            if (abs( a(l) ) .le. small) go to 460
            a(k)     = a(l)
            indc(k)  = i
            k        = k + 1
            go to 480

!           Delete the nonzero from the row file.

  460       lenj     = lenj    - 1
            lenr(i)  = lenr(i) - 1
            lr1      = locr(i)
            last     = lr1 + lenr(i)

            do 470 lrep = lr1, last
               if (indr(lrep) .eq. j) go to 475
  470       continue

  475       indr(lrep) = indr(last)
            indr(last) = 0
            if (i .eq. ilast) lrow = lrow - 1
  480    continue

!        Free the deleted elements from the column file.

         do 490  l  = k, lc2
            indc(l) = 0
  490    continue
         if (atend) lcol = k - 1

!        ---------------------------------------------------------------
!        Deal with the fill-in in column j.
!        ---------------------------------------------------------------
  500    if (ndone .eq. melim) go to 590

!        See if column j already has room for the fill-in.

         if (atend) go to 540
         last   = lc1  + lenj - 1
         l1     = last + 1
         l2     = last + (melim - ndone)
         ! 27 Mar 2001: Be sure it's not at or past end of the col file.
         if (l2 .ge. lcol) go to 520

         do 510 l = l1, l2
            if (indc(l) .ne. 0) go to 520
  510    continue
         go to 540

!        We must move column j to the end of the column file.
!        First, leave some spare room at the end of the
!        current last column.

  520    do 522  l  = lcol + 1, lcol + nspare
            lcol    = l
            indc(l) = 0     ! Spare space is free.
  522    continue

         atend   = .true.
         jlast   = j
         l1      = lc1
         lc1     = lcol + 1
         locc(j) = lc1

         do 525  l     = l1, last
            lcol       = lcol + 1
            a(lcol)    = a(l)
            indc(lcol) = indc(l)
            indc(l)    = 0      ! Free space.
  525    continue

!        ---------------------------------------------------------------
!        Inner loop for the fill-in in column j.
!        This is usually not very expensive.
!        ---------------------------------------------------------------
  540    last          = lc1 + lenj - 1
         ll            = 0

         do 560   lc   = lpivc1, lpivc2
            ll         = ll + 1
            if (markl(ll)  .eq. j    ) go to 560
            aij        = al(ll) * uj
            if (abs( aij ) .le. small) go to 560
            lenj       = lenj + 1
            last       = last + 1
            a(last)    = aij
            i          = indc(lc)
            indc(last) = i
            leni       = lenr(i)

!           Add 1 fill-in to row i if there is already room.
!           27 Mar 2001: Be sure it's not at or past the end
!                        of the row file.

            l          = locr(i) + leni
            if (     l  .ge. lrow) go to 550
            if (indr(l) .gt.    0) go to 550

            indr(l)    = j
            lenr(i)    = leni + 1
            go to 560

!           Row i does not have room for the fill-in.
!           Increment ifill(ll) to count how often this has
!           happened to row i.  Also, add m to the row index
!           indc(last) in column j to mark it as a fill-in that is
!           still pending.
!
!           If this is the first pending fill-in for row i,
!           nfill includes the current length of row i
!           (since the whole row has to be moved later).
!
!           If this is the first pending fill-in for column j,
!           jfill(lu) records the current length of column j
!           (to shorten the search for pending fill-ins later).

  550       if (ifill(ll) .eq. 0) nfill     = nfill + leni + nspare
            if (jfill(lu) .eq. 0) jfill(lu) = lenj
            nfill      = nfill     + 1
            ifill(ll)  = ifill(ll) + 1
            indc(last) = m + i
  560    continue

         if ( atend ) lcol = last

!        End loop for column  j.  Store its final length.

  590    lenc(j) = lenj
  600 continue

!     Successful completion.

      lfirst = 0
      return

!     Interruption.  We have to come back in after the
!     column file is compressed.  Give lfirst a new value.
!     lu and nfill will retain their current values.

  900 lfirst = lr
      return

      end ! subroutine lu1gau

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1mar( m    , n     , lena  , maxmn,
     &                   TCP  , aijtol, Ltol  , maxcol, maxrow,
     &                   ibest, jbest , mbest ,
     &                   a    , indc  , indr  , ip   , iq,
     &                   lenc , lenr  , locc  , locr ,
     &                   iploc, iqloc )

      implicit           double precision (a-h,o-z)
      logical            TCP
      double precision   Ltol      , a(lena)
      integer            indc(lena), indr(lena), ip(m)   , iq(n)   ,
     &                   lenc(n)   , lenr(m)   , iploc(n), iqloc(m)
      integer            locc(n)   , locr(m)

!     ------------------------------------------------------------------
!     lu1mar  uses a Markowitz criterion to select a pivot element
!     for the next stage of a sparse LU factorization,
!     subject to a Threshold Partial Pivoting stability criterion (TPP)
!     that bounds the elements of L.
!
!     00 Jan 1986  Version documented in LUSOL paper:
!                  Gill, Murray, Saunders and Wright (1987),
!                  Maintaining LU factors of a general sparse matrix,
!                  Linear algebra and its applications 88/89, 239-270.
!
!     02 Feb 1989  Following Suhl and Aittoniemi (1987), the largest
!                  element in each column is now kept at the start of
!                  the column, i.e. in position locc(j) of a and indc.
!                  This should speed up the Markowitz searches.
!
!     26 Apr 1989  Both columns and rows searched during spars1 phase.
!                  Only columns searched during spars2 phase.
!                  maxtie replaced by maxcol and maxrow.
!     05 Nov 1993  Initializing  "mbest = m * n"  wasn't big enough when
!                  m = 10, n = 3, and last column had 7 nonzeros.
!     09 Feb 1994  Realised that "mbest = maxmn * maxmn" might overflow.
!                  Changed to    "mbest = maxmn * 1000".
!     27 Apr 2000  On large example from Todd Munson,
!                  that allowed  "if (mbest .le. nz1**2) go to 900"
!                  to exit before any pivot had been found.
!                  Introduced kbest = mbest / nz1.
!                  Most pivots can be rejected with no integer multiply.
!                  True merit is evaluated only if it's as good as the
!                  best so far (or better).  There should be no danger
!                  of integer overflow unless A is incredibly
!                  large and dense.
!
!     10 Sep 2000  TCP, aijtol added for Threshold Complete Pivoting.
!
!     Systems Optimization Laboratory, Stanford University.
!     ------------------------------------------------------------------

      double precision   abest, lbest
      double precision   zero,           one,          gamma
      parameter        ( zero = 0.0d+0,  one = 1.0d+0, gamma  = 2.0d+0 )

!     gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper.

!     ------------------------------------------------------------------
!     Search cols of length nz = 1, then rows of length nz = 1,
!     then   cols of length nz = 2, then rows of length nz = 2, etc.
!     ------------------------------------------------------------------
      abest  = zero
      lbest  = zero
      ibest  = 0
      kbest  = maxmn + 1
      mbest  = -1
      ncol   = 0
      nrow   = 0
      nz1    = 0

      do nz = 1, maxmn
!        nz1    = nz - 1
!        if (mbest .le. nz1**2) go to 900
         if (kbest .le. nz1   ) go to 900
         if (ibest .gt. 0) then
            if (ncol  .ge. maxcol) go to 200
         end if
         if (nz    .gt. m     ) go to 200

!        ---------------------------------------------------------------
!        Search the set of columns of length  nz.
!        ---------------------------------------------------------------
         lq1    = iqloc(nz)
         lq2    = n
         if (nz .lt. m) lq2 = iqloc(nz + 1) - 1

         do 180 lq = lq1, lq2
            ncol   = ncol + 1
            j      = iq(lq)
            lc1    = locc(j)
            lc2    = lc1 + nz1
            amax   = abs( a(lc1) )

!           Test all aijs in this column.
!           amax is the largest element (the first in the column).
!           cmax is the largest multiplier if aij becomes pivot.

            if ( TCP ) then
               if (amax .lt. aijtol) go to 180 ! Nothing in whole column
            end if

            do 160 lc = lc1, lc2
               i      = indc(lc)
               len1   = lenr(i) - 1
!              merit  = nz1 * len1
!              if (merit .gt. mbest) go to 160
               if (len1  .gt. kbest) go to 160

!              aij  has a promising merit.
!              Apply the stability test.
!              We require  aij  to be sufficiently large compared to
!              all other nonzeros in column  j.  This is equivalent
!              to requiring cmax to be bounded by Ltol.

               if (lc .eq. lc1) then

!                 This is the maximum element, amax.
!                 Find the biggest element in the rest of the column
!                 and hence get cmax.  We know cmax .le. 1, but
!                 we still want it exactly in order to break ties.
!                 27 Apr 2002: Settle for cmax = 1.

                  aij    = amax
                  cmax   = one

!                 cmax   = zero
!                 do 140 l = lc1 + 1, lc2
!                    cmax  = max( cmax, abs( a(l) ) )
!  140            continue
!                 cmax   = cmax / amax
               else

!                 aij is not the biggest element, so cmax .ge. 1.
!                 Bail out if cmax will be too big.

                  aij    = abs( a(lc) )
                  if ( TCP ) then ! Absolute test for Complete Pivoting
                     if (aij         .lt.  aijtol) go to 160
                  else !!! TPP
                     if (aij * Ltol  .lt.  amax  ) go to 160
                  end if
                  cmax   = amax / aij
               end if

!              aij  is big enough.  Its maximum multiplier is cmax.

               merit  = nz1 * len1
               if (merit .eq. mbest) then

!                 Break ties.
!                 (Initializing mbest < 0 prevents getting here if
!                 nothing has been found yet.)
!                 In this version we minimize cmax
!                 but if it is already small we maximize the pivot.

                  if (lbest .le. gamma  .and.  cmax .le. gamma) then
                     if (abest .ge. aij ) go to 160
                  else
                     if (lbest .le. cmax) go to 160
                  end if
               end if

!              aij  is the best pivot so far.

               ibest  = i
               jbest  = j
               kbest  = len1
               mbest  = merit
               abest  = aij
               lbest  = cmax
               if (nz .eq. 1) go to 900
  160       continue

!           Finished with that column.

            if (ibest .gt. 0) then
               if (ncol .ge. maxcol) go to 200
            end if
  180    continue

!        ---------------------------------------------------------------
!        Search the set of rows of length  nz.
!        ---------------------------------------------------------------
! 200    if (mbest .le. nz*nz1) go to 900
  200    if (kbest .le. nz    ) go to 900
         if (ibest .gt. 0) then
            if (nrow  .ge. maxrow) go to 290
         end if
         if (nz    .gt. n     ) go to 290

         lp1    = iploc(nz)
         lp2    = m
         if (nz .lt. n) lp2 = iploc(nz + 1) - 1

         do 280 lp = lp1, lp2
            nrow   = nrow + 1
            i      = ip(lp)
            lr1    = locr(i)
            lr2    = lr1 + nz1

            do 260 lr = lr1, lr2
               j      = indr(lr)
               len1   = lenc(j) - 1
!              merit  = nz1 * len1
!              if (merit .gt. mbest) go to 260
               if (len1  .gt. kbest) go to 260

!              aij  has a promising merit.
!              Find where  aij  is in column  j.

               lc1    = locc(j)
               lc2    = lc1 + len1
               amax   = abs( a(lc1) )
               do 220 lc = lc1, lc2
                  if (indc(lc) .eq. i) go to 230
  220          continue

!              Apply the same stability test as above.

  230          aij    = abs( a(lc) )
               if ( TCP ) then   !!! Absolute test for Complete Pivoting
                  if (aij .lt. aijtol) go to 260
               end if

               if (lc .eq. lc1) then

!                 This is the maximum element, amax.
!                 Find the biggest element in the rest of the column
!                 and hence get cmax.  We know cmax .le. 1, but
!                 we still want it exactly in order to break ties.
!                 27 Apr 2002: Settle for cmax = 1.

                  cmax   = one

!                 cmax   = zero
!                 do 240 l = lc1 + 1, lc2
!                    cmax  = max( cmax, abs( a(l) ) )
! 240             continue
!                 cmax   = cmax / amax
               else

!                 aij is not the biggest element, so cmax .ge. 1.
!                 Bail out if cmax will be too big.

                  if ( TCP ) then
                     ! relax
                  else
                     if (aij * Ltol  .lt.  amax) go to 260
                  end if
                  cmax   = amax / aij
               end if

!              aij  is big enough.  Its maximum multiplier is cmax.

               merit  = nz1 * len1
               if (merit .eq. mbest) then

!                 Break ties as before.
!                 (Initializing mbest < 0 prevents getting here if
!                 nothing has been found yet.)

                  if (lbest .le. gamma  .and.  cmax .le. gamma) then
                     if (abest .ge. aij ) go to 260
                  else
                     if (lbest .le. cmax) go to 260
                  end if
               end if

!              aij  is the best pivot so far.

               ibest  = i
               jbest  = j
               kbest  = len1
               mbest  = merit
               abest  = aij
               lbest  = cmax
               if (nz .eq. 1) go to 900
  260       continue

!           Finished with that row.

            if (ibest .gt. 0) then
               if (nrow .ge. maxrow) go to 290
            end if
  280    continue

!        See if it's time to quit.

  290    if (ibest .gt. 0) then
            if (nrow .ge. maxrow  .and.  ncol .ge. maxcol) go to 900
         end if

!        Press on with next nz.

         nz1    = nz
         if (ibest .gt. 0) kbest  = mbest / nz1
      end do

  900 return

      end ! subroutine lu1mar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1mRP( m    , n     , lena  , maxmn,
     &                   Ltol , maxcol, maxrow,
     &                   ibest, jbest , mbest ,
     &                   a    , indc  , indr  , ip   , iq,
     &                   lenc , lenr  , locc  , locr ,
     &                   iploc, iqloc , Amaxr )

      implicit
     &     none
      integer
     &     m, n, lena, maxmn, maxcol, maxrow, ibest, jbest, mbest
      double precision
     &     Ltol
      integer
     &     indc(lena), indr(lena), ip(m)   , iq(n)   ,
     &     lenc(n)   , lenr(m)   , iploc(n), iqloc(m),
     &     locc(n)   , locr(m)
      double precision
     &     a(lena)   , Amaxr(m)

!     ------------------------------------------------------------------
!     lu1mRP  uses a Markowitz criterion to select a pivot element
!     for the next stage of a sparse LU factorization,
!     subject to a Threshold Rook Pivoting stability criterion (TRP)
!     that bounds the elements of L and U.
!
!     11 Jun 2002: First version of lu1mRP derived from lu1mar.
!     11 Jun 2002: Current version of lu1mRP.
!     ------------------------------------------------------------------

      integer
     &     i, j, kbest, lc, lc1, lc2, len1,
     &     lp, lp1, lp2, lq, lq1, lq2, lr, lr1, lr2,
     &     merit, ncol, nrow, nz, nz1
      double precision
     &     abest, aij, amax, atoli, atolj

      double precision   zero
      parameter        ( zero = 0.0d+0 )

!     ------------------------------------------------------------------
!     Search cols of length nz = 1, then rows of length nz = 1,
!     then   cols of length nz = 2, then rows of length nz = 2, etc.
!     ------------------------------------------------------------------
      abest  = zero
      ibest  = 0
      kbest  = maxmn + 1
      mbest  = -1
      ncol   = 0
      nrow   = 0
      nz1    = 0

      do nz = 1, maxmn
!        nz1    = nz - 1
!        if (mbest .le. nz1**2) go to 900
         if (kbest .le. nz1   ) go to 900
         if (ibest .gt. 0) then
            if (ncol  .ge. maxcol) go to 200
         end if
         if (nz    .gt. m     ) go to 200

!        ---------------------------------------------------------------
!        Search the set of columns of length  nz.
!        ---------------------------------------------------------------
         lq1    = iqloc(nz)
         lq2    = n
         if (nz .lt. m) lq2 = iqloc(nz + 1) - 1

         do 180 lq = lq1, lq2
            ncol   = ncol + 1
            j      = iq(lq)
            lc1    = locc(j)
            lc2    = lc1 + nz1
            amax   = abs( a(lc1) )
            atolj  = amax / Ltol    ! Min size of pivots in col j

!           Test all aijs in this column.

            do 160 lc = lc1, lc2
               i      = indc(lc)
               len1   = lenr(i) - 1
!              merit  = nz1 * len1
!              if (merit .gt. mbest) go to 160
               if (len1  .gt. kbest) go to 160

!              aij  has a promising merit.
!              Apply the Threshold Rook Pivoting stability test.
!              First we require aij to be sufficiently large
!              compared to other nonzeros in column j.
!              Then  we require aij to be sufficiently large
!              compared to other nonzeros in row    i.

               aij    = abs( a(lc) )
               if (aij         .lt.  atolj   ) go to 160
               if (aij * Ltol  .lt.  Amaxr(i)) go to 160

!              aij  is big enough.

               merit  = nz1 * len1
               if (merit .eq. mbest) then

!                 Break ties.
!                 (Initializing mbest < 0 prevents getting here if
!                 nothing has been found yet.)

                  if (abest .ge. aij) go to 160
               end if

!              aij  is the best pivot so far.

               ibest  = i
               jbest  = j
               kbest  = len1
               mbest  = merit
               abest  = aij
               if (nz .eq. 1) go to 900
  160       continue

!           Finished with that column.

            if (ibest .gt. 0) then
               if (ncol .ge. maxcol) go to 200
            end if
  180    continue

!        ---------------------------------------------------------------
!        Search the set of rows of length  nz.
!        ---------------------------------------------------------------
! 200    if (mbest .le. nz*nz1) go to 900
  200    if (kbest .le. nz    ) go to 900
         if (ibest .gt. 0) then
            if (nrow  .ge. maxrow) go to 290
         end if
         if (nz    .gt. n     ) go to 290

         lp1    = iploc(nz)
         lp2    = m
         if (nz .lt. n) lp2 = iploc(nz + 1) - 1

         do 280 lp = lp1, lp2
            nrow   = nrow + 1
            i      = ip(lp)
            lr1    = locr(i)
            lr2    = lr1 + nz1
            atoli  = Amaxr(i) / Ltol   ! Min size of pivots in row i

            do 260 lr = lr1, lr2
               j      = indr(lr)
               len1   = lenc(j) - 1
!              merit  = nz1 * len1
!              if (merit .gt. mbest) go to 260
               if (len1  .gt. kbest) go to 260

!              aij  has a promising merit.
!              Find where  aij  is in column j.

               lc1    = locc(j)
               lc2    = lc1 + len1
               amax   = abs( a(lc1) )
               do lc = lc1, lc2
                  if (indc(lc) .eq. i) go to 230
               end do

!              Apply the Threshold Rook Pivoting stability test.
!              First we require aij to be sufficiently large
!              compared to other nonzeros in row    i.
!              Then  we require aij to be sufficiently large
!              compared to other nonzeros in column j.

  230          aij    = abs( a(lc) )
               if (aij         .lt.  atoli) go to 260
               if (aij * Ltol  .lt.  amax ) go to 260

!              aij  is big enough.

               merit  = nz1 * len1
               if (merit .eq. mbest) then

!                 Break ties as before.
!                 (Initializing mbest < 0 prevents getting here if
!                 nothing has been found yet.)

                  if (abest .ge. aij ) go to 260
               end if

!              aij  is the best pivot so far.

               ibest  = i
               jbest  = j
               kbest  = len1
               mbest  = merit
               abest  = aij
               if (nz .eq. 1) go to 900
  260       continue

!           Finished with that row.

            if (ibest .gt. 0) then
               if (nrow .ge. maxrow) go to 290
            end if
  280    continue

!        See if it's time to quit.

  290    if (ibest .gt. 0) then
            if (nrow .ge. maxrow  .and.  ncol .ge. maxcol) go to 900
         end if

!        Press on with next nz.

         nz1    = nz
         if (ibest .gt. 0) kbest  = mbest / nz1
      end do

  900 return

      end ! subroutine lu1mRP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1mCP( m     , n     , lena  , aijtol,
     &                   ibest , jbest , mbest ,
     &                   a     , indc  , indr  ,
     &                   lenc  , lenr  , locc  ,
     &                   Hlen  , Ha    , Hj    )

      integer            m, n, lena, ibest, jbest, mbest, Hlen
      integer            indc(lena), indr(lena),
     &                   lenc(n)   , lenr(m)   , locc(n), Hj(Hlen)
      double precision   aijtol
      double precision   a(lena)   , Ha(Hlen)

!     ------------------------------------------------------------------
!     lu1mCP  uses a Markowitz criterion to select a pivot element
!     for the next stage of a sparse LU factorization,
!     subject to a Threshold Complete Pivoting stability criterion (TCP)
!     that bounds the elements of L and U.
!
!     09 May 2002: First version of lu1mCP.
!                  It searches columns only, using the heap that
!                  holds the largest element in each column.
!     09 May 2002: Current version of lu1mCP.
!     ------------------------------------------------------------------

      integer
     &     j, kheap, lc, lc1, lc2, lenj, maxcol, ncol, nz1
      double precision
     &     abest, aij, amax, cmax, lbest

      double precision   zero,           one,          gamma
      parameter        ( zero = 0.0d+0,  one = 1.0d+0, gamma  = 2.0d+0 )

!     gamma  is "gamma" in the tie-breaking rule TB4 in the LUSOL paper.

!     ------------------------------------------------------------------
!     Search up to maxcol columns stored at the top of the heap.
!     The very top column helps initialize mbest.
!     ------------------------------------------------------------------
      abest  = zero
      lbest  = zero
      ibest  = 0
      jbest  = Hj(1)               ! Column at the top of the heap
      lenj   = lenc(jbest)
      mbest  = lenj * Hlen         ! Bigger than any possible merit
      maxcol = 40                  ! ??? Big question
      ncol   = 0                   ! No. of columns searched

      do 500 kheap = 1, Hlen

         amax   = Ha(kheap)
         if (amax .lt. aijtol) go to 500

         ncol   = ncol + 1
         j      = Hj(kheap)

!        ---------------------------------------------------------------
!        This column has at least one entry big enough (the top one).
!        Search the column for other possibilities.
!        ---------------------------------------------------------------
         lenj   = lenc(j)
         nz1    = lenj - 1
         lc1    = locc(j)
         lc2    = lc1 + nz1
!--      amax   = abs( a(lc1) )

!        Test all aijs in this column.
!        amax is the largest element (the first in the column).
!        cmax is the largest multiplier if aij becomes pivot.

         do 160 lc = lc1, lc2
            i      = indc(lc)
            len1   = lenr(i) - 1
            merit  = nz1 * len1
            if (merit .gt. mbest) go to 160

!           aij  has a promising merit.

            if (lc .eq. lc1) then

!              This is the maximum element, amax.
!              Find the biggest element in the rest of the column
!              and hence get cmax.  We know cmax .le. 1, but
!              we still want it exactly in order to break ties.
!              27 Apr 2002: Settle for cmax = 1.

               aij    = amax
               cmax   = one

!              cmax   = zero
!              do 140 l = lc1 + 1, lc2
!                 cmax  = max( cmax, abs( a(l) ) )
!  140         continue
!              cmax   = cmax / amax
            else

!              aij is not the biggest element, so cmax .ge. 1.
!              Bail out if cmax will be too big.

               aij    = abs( a(lc) )
               if (aij   .lt.  aijtol    ) go to 160
               cmax   = amax / aij
            end if

!           aij  is big enough.  Its maximum multiplier is cmax.

            if (merit .eq. mbest) then

!              Break ties.
!              (Initializing mbest "too big" prevents getting here if
!              nothing has been found yet.)
!              In this version we minimize cmax
!              but if it is already small we maximize the pivot.

               if (lbest .le. gamma  .and.  cmax .le. gamma) then
                  if (abest .ge. aij ) go to 160
               else
                  if (lbest .le. cmax) go to 160
               end if
            end if

!           aij  is the best pivot so far.

            ibest  = i
            jbest  = j
            mbest  = merit
            abest  = aij
            lbest  = cmax
            if (merit .eq. 0) go to 900 ! Col or row of length 1
  160    continue

         if (ncol .ge. maxcol) go to 900
  500 continue

  900 return

      end ! subroutine lu1mCP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1mSP( m    , n     , lena  , maxmn,
     &                   Ltol , maxcol,
     &                   ibest, jbest , mbest ,
     &                   a    , indc  , iq    , locc , iqloc )

      implicit
     &     none
      integer
     &     m, n, lena, maxmn, maxcol, ibest, jbest, mbest
      double precision
     &     Ltol
      integer
     &     indc(lena), iq(n), iqloc(m), locc(n)
      double precision
     &     a(lena)

!     ------------------------------------------------------------------
!     lu1mSP  is intended for symmetric matrices that are either
!     definite or quasi-definite.
!     lu1mSP  uses a Markowitz criterion to select a pivot element for
!     the next stage of a sparse LU factorization of a symmetric matrix,
!     subject to a Threshold Symmetric Pivoting stability criterion
!     (TSP) restricted to diagonal elements to preserve symmetry.
!     This bounds the elements of L and U and should have rank-revealing
!     properties analogous to Threshold Rook Pivoting for unsymmetric
!     matrices.
!
!     14 Dec 2002: First version of lu1mSP derived from lu1mRP.
!                  There is no safeguard to ensure that A is symmetric.
!     14 Dec 2002: Current version of lu1mSP.
!     ------------------------------------------------------------------

      integer
     &     i, j, kbest, lc, lc1, lc2,
     &     lq, lq1, lq2, merit, ncol, nz, nz1
      double precision
     &     abest, aij, amax, atolj

      double precision   zero
      parameter        ( zero = 0.0d+0 )

!     ------------------------------------------------------------------
!     Search cols of length nz = 1, then cols of length nz = 2, etc.
!     ------------------------------------------------------------------
      abest  = zero
      ibest  = 0
      kbest  = maxmn + 1
      mbest  = -1
      ncol   = 0
      nz1    = 0

      do nz = 1, maxmn
!        nz1    = nz - 1
!        if (mbest .le. nz1**2) go to 900
         if (kbest .le. nz1   ) go to 900
         if (ibest .gt. 0) then
            if (ncol  .ge. maxcol) go to 200
         end if
         if (nz    .gt. m     ) go to 200

!        ---------------------------------------------------------------
!        Search the set of columns of length  nz.
!        ---------------------------------------------------------------
         lq1    = iqloc(nz)
         lq2    = n
         if (nz .lt. m) lq2 = iqloc(nz + 1) - 1

         do 180 lq = lq1, lq2
            ncol   = ncol + 1
            j      = iq(lq)
            lc1    = locc(j)
            lc2    = lc1 + nz1
            amax   = abs( a(lc1) )
            atolj  = amax / Ltol    ! Min size of pivots in col j

!           Test all aijs in this column.
!           Ignore everything except the diagonal.

            do 160 lc = lc1, lc2
               i      = indc(lc)
               if (i .ne. j) go to 160     ! Skip off-diagonals.
!              merit  = nz1 * nz1
!              if (merit .gt. mbest) go to 160
               if (nz1   .gt. kbest) go to 160

!              aij  has a promising merit.
!              Apply the Threshold Partial Pivoting stability test
!              (which is equivalent to Threshold Rook Pivoting for
!              symmetric matrices).
!              We require aij to be sufficiently large
!              compared to other nonzeros in column j.

               aij    = abs( a(lc) )
               if (aij .lt. atolj  ) go to 160

!              aij  is big enough.

               merit  = nz1 * nz1
               if (merit .eq. mbest) then

!                 Break ties.
!                 (Initializing mbest < 0 prevents getting here if
!                 nothing has been found yet.)

                  if (abest .ge. aij) go to 160
               end if

!              aij  is the best pivot so far.

               ibest  = i
               jbest  = j
               kbest  = nz1
               mbest  = merit
               abest  = aij
               if (nz .eq. 1) go to 900
  160       continue

!           Finished with that column.

            if (ibest .gt. 0) then
               if (ncol .ge. maxcol) go to 200
            end if
  180    continue

!        See if it's time to quit.

  200    if (ibest .gt. 0) then
            if (ncol .ge. maxcol) go to 900
         end if

!        Press on with next nz.

         nz1    = nz
         if (ibest .gt. 0) kbest  = mbest / nz1
      end do

  900 return

      end ! subroutine lu1mSP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1pen( m     , melim , ncold , nspare, ilast,
     &                   lpivc1, lpivc2, lpivr1, lpivr2, lrow ,
     &                   lenc  , lenr  , locc  , locr  ,
     &                   indc  , indr  , ifill , jfill )

      integer            indc(*)     , indr(*)     , lenc(*), lenr(*),
     &                   ifill(melim), jfill(ncold)
      integer            locc(*)     , locr(*)

!     ------------------------------------------------------------------
!     lu1pen deals with pending fill-in in the row file.
!     ifill(ll) says if a row involved in the new column of L
!               has to be updated.  If positive, it is the total
!               length of the final updated row.
!     jfill(lu) says if a column involved in the new row of U
!               contains any pending fill-ins.  If positive, it points
!               to the first fill-in in the column that has yet to be
!               added to the row file.
!
!     16 Apr 1989: First version of lu1pen.
!     23 Mar 2001: ilast used and updated.
!     ------------------------------------------------------------------

         ll     = 0

         do 650 lc  = lpivc1, lpivc2
            ll      = ll + 1
            if (ifill(ll) .eq. 0) go to 650

            ! Another row has pending fill.
            ! First, add some spare space at the end
            ! of the current last row.

            do 620  l  = lrow + 1, lrow + nspare
               lrow    = l
               indr(l) = 0
  620       continue

            ! Now move row i to the end of the row file.

            i       = indc(lc)
            ilast   = i
            lr1     = locr(i)
            lr2     = lr1 + lenr(i) - 1
            locr(i) = lrow + 1

            do 630 lr = lr1, lr2
               lrow       = lrow + 1
               indr(lrow) = indr(lr)
               indr(lr)   = 0
  630       continue

            lrow    = lrow + ifill(ll)
  650    continue

!        Scan all columns of  D  and insert the pending fill-in
!        into the row file.

         lu     = 1

         do 680 lr = lpivr1, lpivr2
            lu     = lu + 1
            if (jfill(lu) .eq. 0) go to 680
            j      = indr(lr)
            lc1    = locc(j) + jfill(lu) - 1
            lc2    = locc(j) + lenc(j)   - 1

            do 670 lc = lc1, lc2
               i      = indc(lc) - m
               if (i .gt. 0) then
                  indc(lc)   = i
                  last       = locr(i) + lenr(i)
                  indr(last) = j
                  lenr(i)    = lenr(i) + 1
               end if
  670       continue
  680    continue

      end ! subroutine lu1pen

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1mxc( k1, k2, iq,
     &                   a, indc, lenc, locc )

      implicit           none
      integer            k1, k2
      integer            iq(k2), indc(*), lenc(*), locc(*)
      double precision   a(*)

!     ------------------------------------------------------------------
!     lu1mxc  moves the largest element in each of columns iq(k1:k2)
!     to the top of its column.
!     If k1 > k2, nothing happens.
!
!     06 May 2002: (and earlier)
!                  All columns k1:k2 must have one or more elements.
!     07 May 2002: Allow for empty columns.  The heap routines need to
!                  find 0.0 as the "largest element".
!     ------------------------------------------------------------------

      integer             i, j, k, l, lc, lc1, lc2, lenj
      double precision    amax
      double precision    zero
      parameter         ( zero = 0.0d+0 )

      do k = k1, k2
         j      = iq(k)
         lc1    = locc(j)
         lenj   = lenc(j)

         if (lenj .eq. 0) then
            a(lc1) = zero
         else

            ! The next 10 lines are equivalent to
            ! l      = idamax( lenc(j), a(lc1), 1 )  +  lc1 - 1
            ! >>>>>>>>
            lc2    = lc1 + lenc(j) - 1
            amax   = abs( a(lc1) )
            l      = lc1

            do lc = lc1+1, lc2
               if (amax .lt. abs( a(lc) )) then
                  amax   =  abs( a(lc) )
                  l      =  lc
               end if
            end do
            ! >>>>>>>>

            if (l .gt. lc1) then
               amax      = a(l)
               a(l)      = a(lc1)
               a(lc1)    = amax
               i         = indc(l)
               indc(l)   = indc(lc1)
               indc(lc1) = i
            end if
         end if
      end do

      end ! subroutine lu1mxc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1mxr( k1, k2, ip, Amaxr,
     &                   a, indc, lenc, locc, indr, lenr, locr )

      implicit           none
      integer            k1, k2
      integer            ip(k2),
     &                   indc(*), lenc(*), locc(*),
     &                   indr(*), lenr(*), locr(*)
      double precision   a(*), Amaxr(*)

!     ------------------------------------------------------------------
!     lu1mxr  finds the largest element in each of row ip(k1:k2)
!     and stores it in Amaxr(*).  The nonzeros are stored column-wise
!     in (a,indc,lenc,locc) and their structure is row-wise
!     in (  indr,lenr,locr).
!     If k1 > k2, nothing happens.
!
!     11 Jun 2002: First version of lu1mxr.
!                  Allow for empty columns.
!     19 Dec 2004: Declare Amaxr(*) instead of Amaxr(k2)
!                  to stop grumbles from the NAG compiler.
!                  (Mentioned by Mick Pont.)
!     ------------------------------------------------------------------

      integer             i, j, k, lc, lc1, lc2, lr, lr1, lr2
      double precision    amax
      double precision    zero
      parameter         ( zero = 0.0d+0 )

      do k = k1, k2
         amax   = zero
         i      = ip(k)

         ! Find largest element in row i.

         lr1    = locr(i)
         lr2    = lr1 + lenr(i) - 1

         do lr = lr1, lr2
            j      = indr(lr)

            ! Find where  aij  is in column  j.

            lc1    = locc(j)
            lc2    = lc1 + lenc(j) - 1
            do lc = lc1, lc2
               if (indc(lc) .eq. i) go to 230
            end do

  230       amax   = max( amax, abs( a(lc) ) )
         end do

         Amaxr(i) = amax
      end do

      end ! subroutine lu1mxr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1or1( m, n, nelem, lena, small,
     &                   a, indc, indr, lenc, lenr,
     &                   Amax, numnz, lerr, inform )

      implicit           double precision (a-h,o-z)
      double precision   a(lena)
      integer            indc(lena), indr(lena)
      integer            lenc(n), lenr(m)

!     ------------------------------------------------------------------
!     lu1or1  organizes the elements of an  m by n  matrix  A  as
!     follows.  On entry, the parallel arrays   a, indc, indr,
!     contain  nelem  entries of the form     aij,    i,    j,
!     in any order.  nelem  must be positive.
!
!     Entries not larger than the input parameter  small  are treated as
!     zero and removed from   a, indc, indr.  The remaining entries are
!     defined to be nonzero.  numnz  returns the number of such nonzeros
!     and  Amax  returns the magnitude of the largest nonzero.
!     The arrays  lenc, lenr  return the number of nonzeros in each
!     column and row of  A.
!
!     inform = 0  on exit, except  inform = 1  if any of the indices in
!     indc, indr  imply that the element  aij  lies outside the  m by n
!     dimensions of  A.
!
!     xx Feb 1985: Original version.
!     17 Oct 2000: a, indc, indr now have size lena to allow nelem = 0.
!     ------------------------------------------------------------------
      do 10 i = 1, m
         lenr(i) = 0
   10 continue

      do 20 j = 1, n
         lenc(j) = 0
   20 continue

      Amax   = 0.0d+0
      numnz  = nelem
      l      = nelem + 1

      do 100 ldummy = 1, nelem
         l      = l - 1
         if (abs( a(l) ) .gt. small) then
            i      = indc(l)
            j      = indr(l)
            Amax   = max( Amax, abs( a(l) ) )
            if (i .lt. 1  .or.  i .gt. m) go to 910
            if (j .lt. 1  .or.  j .gt. n) go to 910
            lenr(i) = lenr(i) + 1
            lenc(j) = lenc(j) + 1
         else

!           Replace a negligible element by last element.  Since
!           we are going backwards, we know the last element is ok.

            a(l)    = a(numnz)
            indc(l) = indc(numnz)
            indr(l) = indr(numnz)
            numnz   = numnz - 1
         end if
  100 continue

      inform = 0
      return

  910 lerr   = l
      inform = 1
      return

      end ! subroutine lu1or1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1or2( n, numa, lena, a, inum, jnum, len, loc )
      integer            inum(lena), jnum(lena), len(n)
      integer            loc(n)
      double precision   a(lena), ace, acep

!     ------------------------------------------------------------------
!     lu1or2  sorts a list of matrix elements  a(i,j)  into column
!     order, given  numa  entries  a(i,j),  i,  j  in the parallel
!     arrays  a, inum, jnum  respectively.  The matrix is assumed
!     to have  n  columns and an arbitrary number of rows.
!
!     On entry,  len(*)  must contain the length of each column.
!
!     On exit,  a(*) and inum(*)  are sorted,  jnum(*) = 0,  and
!     loc(j)  points to the start of column j.
!
!     lu1or2  is derived from  mc20ad,  a routine in the Harwell
!     Subroutine Library, author J. K. Reid.
!
!     xx Feb 1985: Original version.
!     17 Oct 2000: a, inum, jnum now have size lena to allow nelem = 0.
!     ------------------------------------------------------------------

!     Set  loc(j)  to point to the beginning of column  j.

      l = 1
      do 150 j  = 1, n
         loc(j) = l
         l      = l + len(j)
  150 continue

!     Sort the elements into column order.
!     The algorithm is an in-place sort and is of order  numa.

      do 230 i = 1, numa
!        Establish the current entry.
         jce     = jnum(i)
         if (jce .eq. 0) go to 230
         ace     = a(i)
         ice     = inum(i)
         jnum(i) = 0

!        Chain from current entry.

         do 200 j = 1, numa

!           The current entry is not in the correct position.
!           Determine where to store it.

            l        = loc(jce)
            loc(jce) = loc(jce) + 1

!           Save the contents of that location.

            acep = a(l)
            icep = inum(l)
            jcep = jnum(l)

!           Store current entry.

            a(l)    = ace
            inum(l) = ice
            jnum(l) = 0

!           If next current entry needs to be processed,
!           copy it into current entry.

            if (jcep .eq. 0) go to 230
            ace = acep
            ice = icep
            jce = jcep
  200    continue
  230 continue

!     Reset loc(j) to point to the start of column j.

      ja = 1
      do 250 j  = 1, n
         jb     = loc(j)
         loc(j) = ja
         ja     = jb
  250 continue

      end ! subroutine lu1or2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1or3( m, n, lena,
     &                   indc, lenc, locc, iw,
     &                   lerr, inform )

      integer            indc(lena), lenc(n), iw(m)
      integer            locc(n)

!     ------------------------------------------------------------------
!     lu1or3  looks for duplicate elements in an  m by n  matrix  A
!     defined by the column list  indc, lenc, locc.
!     iw  is used as a work vector of length  m.
!
!     xx Feb 1985: Original version.
!     17 Oct 2000: indc, indr now have size lena to allow nelem = 0.
!     ------------------------------------------------------------------

      do 100 i = 1, m
         iw(i) = 0
  100 continue

      do 200 j = 1, n
         if (lenc(j) .gt. 0) then
            l1    = locc(j)
            l2    = l1 + lenc(j) - 1

            do 150 l = l1, l2
               i     = indc(l)
               if (iw(i) .eq. j) go to 910
               iw(i) = j
  150       continue
         end if
  200 continue

      inform = 0
      return

  910 lerr   = l
      inform = 1
      return

      end ! subroutine lu1or3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1or4( m, n, nelem, lena,
     &                   indc, indr, lenc, lenr, locc, locr )

      integer            indc(lena), indr(lena), lenc(n), lenr(m)
      integer            locc(n), locr(m)

!     ------------------------------------------------------------------
!     lu1or4     constructs a row list  indr, locr
!     from a corresponding column list  indc, locc,
!     given the lengths of both columns and rows in  lenc, lenr.
!
!     xx Feb 1985: Original version.
!     17 Oct 2000: indc, indr now have size lena to allow nelem = 0.
!     ------------------------------------------------------------------

!     Initialize  locr(i)  to point just beyond where the
!     last component of row  i  will be stored.

      l      = 1
      do 10 i = 1, m
         l       = l + lenr(i)
         locr(i) = l
   10 continue

!     By processing the columns backwards and decreasing  locr(i)
!     each time it is accessed, it will end up pointing to the
!     beginning of row  i  as required.

      l2     = nelem
      j      = n + 1

      do 40 jdummy = 1, n
         j  = j - 1
         if (lenc(j) .gt. 0) then
            l1 = locc(j)

            do 30 l = l1, l2
               i        = indc(l)
               lr       = locr(i) - 1
               locr(i)  = lr
               indr(lr) = j
   30       continue

            l2     = l1 - 1
         end if
   40 continue

      end ! subroutine lu1or4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1pq1( m, n, len, iperm, loc, inv, num )

      integer            len(m), iperm(m), loc(n), inv(m), num(n)

!     ------------------------------------------------------------------
!     lu1pq1  constructs a permutation  iperm  from the array  len.
!
!     On entry:
!     len(i)  holds the number of nonzeros in the i-th row (say)
!             of an m by n matrix.
!     num(*)  can be anything (workspace).
!
!     On exit:
!     iperm   contains a list of row numbers in the order
!             rows of length 0,  rows of length 1,..., rows of length n.
!     loc(nz) points to the first row containing  nz  nonzeros,
!             nz = 1, n.
!     inv(i)  points to the position of row i within iperm(*).
!     ------------------------------------------------------------------

!     Count the number of rows of each length.

      nzero  = 0
      do 10  nz  = 1, n
         num(nz) = 0
         loc(nz) = 0
   10 continue

      do 20  i   = 1, m
         nz      = len(i)
         if (nz .eq. 0) then
            nzero   = nzero   + 1
         else
            num(nz) = num(nz) + 1
         end if
   20 continue

!     Set starting locations for each length.

      l      = nzero + 1
      do 60  nz  = 1, n
         loc(nz) = l
         l       = l + num(nz)
         num(nz) = 0
   60 continue

!     Form the list.

      nzero  = 0
      do 100  i   = 1, m
         nz       = len(i)
         if (nz .eq. 0) then
            nzero    = nzero + 1
            iperm(nzero) = i
         else
            l        = loc(nz) + num(nz)
            iperm(l) = i
            num(nz)  = num(nz) + 1
         end if
  100 continue

!     Define the inverse of iperm.

      do 120 l  = 1, m
         i      = iperm(l)
         inv(i) = l
  120 continue

      end ! subroutine lu1pq1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1pq2( nzpiv, nzchng,
     &                   indr , lenold, lennew, iqloc, iq, iqinv )

      integer            indr(nzpiv), lenold(nzpiv), lennew(*),
     &                   iqloc(*)   , iq(*)        , iqinv(*)

!     ===============================================================
!     lu1pq2 frees the space occupied by the pivot row,
!     and updates the column permutation iq.
!
!     Also used to free the pivot column and update the row perm ip.
!
!     nzpiv   (input)    is the length of the pivot row (or column).
!     nzchng  (output)   is the net change in total nonzeros.
!
!     14 Apr 1989  First version.
!     ===============================================================

      nzchng = 0

      do 200  lr  = 1, nzpiv
         j        = indr(lr)
         indr(lr) = 0
         nz       = lenold(lr)
         nznew    = lennew(j)

         if (nz .ne. nznew) then
            l        = iqinv(j)
            nzchng   = nzchng + (nznew - nz)

!           l above is the position of column j in iq  (so j = iq(l)).

            if (nz .lt. nznew) then

!              Column  j  has to move towards the end of  iq.

  110          next        = nz + 1
               lnew        = iqloc(next) - 1
               if (lnew .ne. l) then
                  jnew        = iq(lnew)
                  iq(l)       = jnew
                  iqinv(jnew) = l
               end if
               l           = lnew
               iqloc(next) = lnew
               nz          = next
               if (nz .lt. nznew) go to 110
            else

!              Column  j  has to move towards the front of  iq.

  120          lnew        = iqloc(nz)
               if (lnew .ne. l) then
                  jnew        = iq(lnew)
                  iq(l)       = jnew
                  iqinv(jnew) = l
               end if
               l           = lnew
               iqloc(nz)   = lnew + 1
               nz          = nz   - 1
               if (nz .gt. nznew) go to 120
            end if

            iq(lnew) = j
            iqinv(j) = lnew
         end if
  200 continue

      end ! subroutine lu1pq2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1pq3( n, len, iperm, iw, nrank )

      integer            len(n), iperm(n)
      integer            iw(n)

!     ------------------------------------------------------------------
!     lu1pq3  looks at the permutation  iperm(*)  and moves any entries
!     to the end whose corresponding length  len(*)  is zero.
!
!     09 Feb 1994: Added work array iw(*) to improve efficiency.
!     ------------------------------------------------------------------

      nrank  = 0
      nzero  = 0

      do 10 k = 1, n
         i    = iperm(k)

         if (len(i) .eq. 0) then
            nzero        = nzero + 1
            iw(nzero)    = i
         else
            nrank        = nrank + 1
            iperm(nrank) = i
         end if
   10 continue

      do 20 k = 1, nzero
         iperm(nrank + k) = iw(k)
   20 continue

      end ! subroutine lu1pq3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1rec( n, reals, luparm, ltop,
     &                   lena, a, ind, len, loc )

      logical            reals
      integer            luparm(30), ltop
      double precision   a(lena)
      integer            ind(lena), len(n)
      integer            loc(n)

!     ------------------------------------------------------------------
!     00 Jun 1983: Original version of lu1rec followed John Reid's
!                  compression routine in LA05.  It recovered
!                  space in ind(*) and optionally a(*)
!                  by eliminating entries with ind(l) = 0.
!                  The elements of ind(*) could not be negative.
!                  If len(i) was positive, entry i contained
!                  that many elements, starting at  loc(i).
!                  Otherwise, entry i was eliminated.
!
!     23 Mar 2001: Realised we could have len(i) = 0 in rare cases!
!                  (Mostly during TCP when the pivot row contains
!                  a column of length 1 that couldn't be a pivot.)
!                  Revised storage scheme to
!                     keep        entries with       ind(l) >  0,
!                     squeeze out entries with -n <= ind(l) <= 0,
!                  and to allow len(i) = 0.
!                  Empty items are moved to the end of the compressed
!                  ind(*) and/or a(*) arrays are given one empty space.
!                  Items with len(i) < 0 are still eliminated.
!
!     27 Mar 2001: Decided to use only ind(l) > 0 and = 0 in lu1fad.
!                  Still have to keep entries with len(i) = 0.
!
!     On exit:
!     ltop         is the length of useful entries in ind(*), a(*).
!     ind(ltop+1)  is "i" such that len(i), loc(i) belong to the last
!                  item in ind(*), a(*).
!     ------------------------------------------------------------------

      nempty = 0

      do 10 i = 1, n
         leni = len(i)
         if (leni .gt. 0) then
            l      = loc(i) + leni - 1
            len(i) = ind(l)
            ind(l) = - (n + i)
         else if (leni .eq. 0) then
            nempty = nempty + 1
         end if
   10 continue

      k      = 0
      klast  = 0    ! Previous k
      ilast  = 0    ! Last entry moved.

      do 20 l = 1, ltop
         i    = ind(l)
         if (i .gt. 0) then
            k      = k + 1
            ind(k) = i
            if (reals) a(k) = a(l)

         else if (i .lt. -n) then

!           This is the end of entry  i.

            i      = - (i + n)
            ilast  = i
            k      = k + 1
            ind(k) = len(i)
            if (reals) a(k) = a(l)
            loc(i) = klast + 1
            len(i) = k     - klast
            klast  = k
         end if
   20 continue

!     Move any empty items to the end, adding 1 free entry for each.

      if (nempty .gt. 0) then
         do i = 1, n
            if (len(i) .eq. 0) then
               k      = k + 1
               loc(i) = k
               ind(k) = 0
               ilast  = i
            end if
         end do
      end if

      nout   = luparm(1)
      lprint = luparm(2)
      if (lprint .ge. 50) write(nout, 1000) ltop, k, reals, nempty
      luparm(26) = luparm(26) + 1  ! ncp

!     Return ilast in ind(ltop + 1).

      ltop        = k
      ind(ltop+1) = ilast
      return

 1000 format(' lu1rec.  File compressed from', i10, '   to', i10, l3,
     &       '  nempty =', i8)

      end ! subroutine lu1rec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1slk( m, n, lena, iq, iqloc, a, locc, w )

      implicit           none
      integer            m, n, lena
      integer            iq(n), iqloc(m), locc(n)
      double precision   a(lena), w(n)

!     ------------------------------------------------------------------
!     lu1slk  sets w(j) > 0 if column j is a unit vector.
!
!     21 Nov 2000: First version.  lu1fad needs it for TCP.
!                  Note that w(*) is nominally an integer array,
!                  but the only spare space is the double array w(*).
!     ------------------------------------------------------------------

      integer            j, lc1, lq, lq1, lq2

      do j = 1, n
         w(j) = 0.0d+0
      end do

      lq1    = iqloc(1)
      lq2    = n
      if (m .gt. 1) lq2 = iqloc(2) - 1

      do lq = lq1, lq2
         j      = iq(lq)
         lc1    = locc(j)
         if (abs( a(lc1) ) .eq. 1.0d+0) then
            w(j) = 1.0d+0
         end if
      end do

      end ! subroutine lu1slk

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1ful( m     , n    , lena , lenD , lu1 , TPP,
     &                   mleft , nleft, nrank, nrowu,
     &                   lenL  , lenU , nsing,
     &                   keepLU, small,
     &                   a     , d    , indc , indr , ip  , iq,
     &                   lenc  , lenr , locc , ipinv, ipvt )

      implicit           double precision (a-h,o-z)
      logical            TPP       , keepLU
      double precision   a(lena)   , d(lenD)
      integer            indc(lena), indr(lena), ip(m)   , iq(n),
     &                   lenc(n)   , lenr(m)   , ipinv(m)
      integer            locc(n)   , ipvt(m)

!     ------------------------------------------------------------------
!     lu1ful computes a dense (full) LU factorization of the
!     mleft by nleft matrix that remains to be factored at the
!     beginning of the nrowu-th pass through the main loop of lu1fad.
!
!     02 May 1989: First version.
!     05 Feb 1994: Column interchanges added to lu1DPP.
!     08 Feb 1994: ipinv reconstructed, since lu1pq3 may alter ip.
!     ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )


!     ------------------------------------------------------------------
!     If lu1pq3 moved any empty rows, reset ipinv = inverse of ip.
!     ------------------------------------------------------------------
      if (nrank .lt. m) then
         do 100 l    = 1, m
            i        = ip(l)
            ipinv(i) = l
  100    continue
      end if

!     ------------------------------------------------------------------
!     Copy the remaining matrix into the dense matrix D.
!     ------------------------------------------------------------------
!     call dload ( lenD, zero, d, 1 )
      do j = 1, lenD
         d(j) = zero
      end do
      ipbase = nrowu - 1
      ldbase = 1 - nrowu

      do 200 lq = nrowu, n
         j      = iq(lq)
         lc1    = locc(j)
         lc2    = lc1 + lenc(j) - 1

         do 150 lc = lc1, lc2
            i      = indc(lc)
            ld     = ldbase + ipinv(i)
            d(ld)  = a(lc)
  150    continue

         ldbase = ldbase + mleft
  200 continue

!     ------------------------------------------------------------------
!     Call our favorite dense LU factorizer.
!     ------------------------------------------------------------------
      if ( TPP ) then
         call lu1DPP( d, mleft, mleft, nleft, small, nsing,
     &                ipvt, iq(nrowu) )
      else
         call lu1DCP( d, mleft, mleft, nleft, small, nsing,
     &                ipvt, iq(nrowu) )
      end if

!     ------------------------------------------------------------------
!     Move D to the beginning of A,
!     and pack L and U at the top of a, indc, indr.
!     In the process, apply the row permutation to ip.
!     lkk points to the diagonal of U.
!     ------------------------------------------------------------------
      call dcopy ( lenD, d, 1, a, 1 )

      ldiagU = lena   - n
      lkk    = 1
      lkn    = lenD  - mleft + 1
      lu     = lu1

      do 450  k = 1, min( mleft, nleft )
         l1     = ipbase + k
         l2     = ipbase + ipvt(k)
         if (l1 .ne. l2) then
            i      = ip(l1)
            ip(l1) = ip(l2)
            ip(l2) = i
         end if
         ibest  = ip( l1 )
         jbest  = iq( l1 )

         if ( keepLU ) then
!           ===========================================================
!           Pack the next column of L.
!           ===========================================================
            la     = lkk
            ll     = lu
            nrowd  = 1

            do 410  i = k + 1, mleft
               la     = la + 1
               ai     = a(la)
               if (abs( ai ) .gt. small) then
                  nrowd    = nrowd + 1
                  ll       = ll    - 1
                  a(ll)    = ai
                  indc(ll) = ip( ipbase + i )
                  indr(ll) = ibest
               end if
  410       continue

!           ===========================================================
!           Pack the next row of U.
!           We go backwards through the row of D
!           so the diagonal ends up at the front of the row of  U.
!           Beware -- the diagonal may be zero.
!           ===========================================================
            la     = lkn + mleft
            lu     = ll
            ncold  = 0

            do 420  j = nleft, k, -1
               la     = la - mleft
               aj     = a(la)
               if (abs( aj ) .gt. small  .or.  j .eq. k) then
                  ncold    = ncold + 1
                  lu       = lu    - 1
                  a(lu)    = aj
                  indr(lu) = iq( ipbase + j )
               end if
  420       continue

            lenr(ibest) = - ncold
            lenc(jbest) = - nrowd
            lenL        =   lenL + nrowd - 1
            lenU        =   lenU + ncold
            lkn         =   lkn  + 1

         else
!           ===========================================================
!           Store just the diagonal of U, in natural order.
!           ===========================================================
            a(ldiagU + jbest) = a(lkk)
         end if

         lkk    = lkk  + mleft + 1
  450 continue

      end ! subroutine lu1ful

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1DPP( a, lda, m, n, small, nsing,
     &                   ipvt, iq )

      implicit           none
      integer            lda, m, n, nsing
      integer            ipvt(m), iq(n)
      double precision   a(lda,n), small

!     ------------------------------------------------------------------
!     lu1DPP factors a dense m x n matrix A by Gaussian elimination,
!     using row interchanges for stability, as in dgefa from LINPACK.
!     This version also uses column interchanges if all elements in a
!     pivot column are smaller than (or equal to) "small".  Such columns
!     are changed to zero and permuted to the right-hand end.
!
!     As in LINPACK, ipvt(*) keeps track of pivot rows.
!     Rows of U are interchanged, but we don't have to physically
!     permute rows of L.  In contrast, column interchanges are applied
!     directly to the columns of both L and U, and to the column
!     permutation vector iq(*).
!
!     02 May 1989: First version derived from dgefa
!                  in LINPACK (version dated 08/14/78).
!     05 Feb 1994: Generalized to treat rectangular matrices
!                  and use column interchanges when necessary.
!                  ipvt is retained, but column permutations are applied
!                  directly to iq(*).
!     21 Dec 1994: Bug found via example from Steve Dirkse.
!                  Loop 100 added to set ipvt(*) for singular rows.
!     26 Mar 2006: nsing redefined (see below).
!                  Changed to implicit none.
!     ------------------------------------------------------------------
!
!     On entry:
!
!        a       Array holding the matrix A to be factored.
!
!        lda     The leading dimension of the array  a.
!
!        m       The number of rows    in  A.
!
!        n       The number of columns in  A.
!
!        small   A drop tolerance.  Must be zero or positive.
!
!     On exit:
!
!        a       An upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        nsing   Number of singularities detected.
!                26 Mar 2006: nsing redefined to be more meaningful.
!                Users may define rankU = n - nsing and regard
!                U as upper-trapezoidal, with the first rankU columns
!                being triangular and the rest trapezoidal.
!                It would be better to return rankU, but we still
!                return nsing for compatibility (even though lu1fad
!                no longer uses it).
!
!        ipvt    Records the pivot rows.
!
!        iq      A vector to which column interchanges are applied.
!     ------------------------------------------------------------------

      double precision    t
      integer             idamax, i, j, k, kp1, l, last, lencol, rankU
      double precision    zero         ,  one
      parameter         ( zero = 0.0d+0,  one = 1.0d+0 )


      rankU  = 0
      k      = 1
      last   = n

!     ------------------------------------------------------------------
!     Start of elimination loop.
!     ------------------------------------------------------------------
   10 kp1    = k + 1
      lencol = m - k + 1

      ! Find l, the pivot row.

      l       = idamax( lencol, a(k,k), 1 ) + k - 1
      ipvt(k) = l

      if (abs( a(l,k) ) .le. small) then
         !==============================================================
         ! Do column interchange, changing old pivot column to zero.
         ! Reduce "last" and try again with same k.
         !==============================================================
         j        = iq(last)
         iq(last) = iq(k)
         iq(k)    = j

         do i = 1, k - 1
            t         = a(i,last)
            a(i,last) = a(i,k)
            a(i,k)    = t
         end do

         do i = k, m
            t         = a(i,last)
            a(i,last) = zero
            a(i,k)    = t
         end do

         last     = last - 1
         if (k .le. last) go to 10

      else
         rankU  = rankU + 1
         if (k .lt. m) then
            !===========================================================
            !Do row interchange if necessary.
            !===========================================================
            if (l .ne. k) then
               t      = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
            end if

            !===========================================================
            ! Compute multipliers.
            ! Do row elimination with column indexing.
            !===========================================================
            t = - one / a(k,k)
            call dscal ( m-k, t, a(kp1,k), 1 )

            do j = kp1, last
               t    = a(l,j)
               if (l .ne. k) then
                  a(l,j) = a(k,j)
                  a(k,j) = t
               end if
               call daxpy ( m-k, t, a(kp1,k), 1, a(kp1,j), 1 )
            end do

            k = k + 1
            if (k .le. last) go to 10
         end if
      end if

      ! Set ipvt(*) for singular rows.

      do 100 k = last + 1, m
         ipvt(k) = k
  100 continue

      nsing  = n - rankU

      end ! subroutine lu1DPP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu1DCP( a, lda, m, n, small, nsing,
     &                   ipvt, iq )

      implicit           none
      integer            lda, m, n, nsing
      integer            ipvt(m), iq(n)
      double precision   a(lda,n), small

!     ------------------------------------------------------------------
!     lu1DCP factors a dense m x n matrix A by Gaussian elimination,
!     using Complete Pivoting (row and column interchanges) for
!     stability.
!     This version also uses column interchanges if all elements in a
!     pivot column are smaller than (or equal to) "small".  Such columns
!     are changed to zero and permuted to the right-hand end.
!
!     As in LINPACK's dgefa, ipvt(*) keeps track of pivot rows.
!     Rows of U are interchanged, but we don't have to physically
!     permute rows of L.  In contrast, column interchanges are applied
!     directly to the columns of both L and U, and to the column
!     permutation vector iq(*).
!
!     01 May 2002: First dense Complete Pivoting, derived from lu1DPP.
!     07 May 2002: Another break needed at end of first loop.
!     26 Mar 2006: Cosmetic mods while looking for "nsing" bug when m<n.
!                  nsing redefined (see below).
!                  Changed to implicit none.
!     ------------------------------------------------------------------
!
!     On entry:
!
!        a       Array holding the matrix A to be factored.
!
!        lda     The leading dimension of the array  a.
!
!        m       The number of rows    in  A.
!
!        n       The number of columns in  A.
!
!        small   A drop tolerance.  Must be zero or positive.
!
!     On exit:
!
!        a       An upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        nsing   Number of singularities detected.
!                26 Mar 2006: nsing redefined to be more meaningful.
!                Users may define rankU = n - nsing and regard
!                U as upper-trapezoidal, with the first rankU columns
!                being triangular and the rest trapezoidal.
!                It would be better to return rankU, but we still
!                return nsing for compatibility (even though lu1fad
!                no longer uses it).
!
!        ipvt    Records the pivot rows.
!
!        iq      A vector to which column interchanges are applied.
!     ------------------------------------------------------------------

      double precision    aijmax, ajmax, t
      integer             idamax, i, imax, j, jlast, jmax, jnew,
     &                    k, kp1, l, last, lencol, rankU
      double precision    zero         ,  one
      parameter         ( zero = 0.0d+0,  one = 1.0d+0 )


      rankU  = 0
      lencol = m + 1
      last   = n

      !-----------------------------------------------------------------
      ! Start of elimination loop.
      !-----------------------------------------------------------------
      do k = 1, n
         kp1    = k + 1
         lencol = lencol - 1

         ! Find the biggest aij in row imax and column jmax.

         aijmax = zero
         imax   = k
         jmax   = k
         jlast  = last

         do j = k, jlast
   10       l      = idamax( lencol, a(k,j), 1 ) + k - 1
            ajmax  = abs( a(l,j) )

            if (ajmax .le. small) then
               !========================================================
               ! Do column interchange, changing old column to zero.
               ! Reduce  "last"  and try again with same j.
               !========================================================
               jnew     = iq(last)
               iq(last) = iq(j)
               iq(j)    = jnew

               do i = 1, k - 1
                  t         = a(i,last)
                  a(i,last) = a(i,j)
                  a(i,j)    = t
               end do

               do i = k, m
                  t         = a(i,last)
                  a(i,last) = zero
                  a(i,j)    = t
               end do

               last   = last - 1
               if (j .le. last) go to 10 ! repeat
               go to 200                 ! break
            end if

            ! Check if this column has biggest aij so far.

            if (aijmax .lt. ajmax) then
                aijmax  =   ajmax
                imax    =   l
                jmax    =   j
            end if

            if (j .ge. last) go to 200   ! break
         end do

  200    ipvt(k) = imax

         if (jmax .ne. k) then
            !==========================================================
            ! Do column interchange (k and jmax).
            !==========================================================
            jnew     = iq(jmax)
            iq(jmax) = iq(k)
            iq(k)    = jnew

            do i = 1, m
               t         = a(i,jmax)
               a(i,jmax) = a(i,k)
               a(i,k)    = t
            end do
         end if

         if (k .lt. m) then
            !===========================================================
            ! Do row interchange if necessary.
            !===========================================================
            t         = a(imax,k)
            if (imax .ne. k) then
               a(imax,k) = a(k,k)
               a(k,k)    = t
            end if

            !===========================================================
            ! Compute multipliers.
            ! Do row elimination with column indexing.
            !===========================================================
            t      = - one / t
            call dscal ( m-k, t, a(kp1,k), 1 )

            do j = kp1, last
               t         = a(imax,j)
               if (imax .ne. k) then
                  a(imax,j) = a(k,j)
                  a(k,j)    = t
               end if
               call daxpy ( m-k, t, a(kp1,k), 1, a(kp1,j), 1 )
            end do

         else
            go to 500               ! break
         end if

         if (k .ge. last) go to 500 ! break
      end do

      ! Set ipvt(*) for singular rows.

  500 do k = last + 1, m
         ipvt(k) = k
      end do

      nsing  = n - rankU

      end ! subroutine lu1DCP
