************************************************************************
*
*     File  lusol8a.f
*
*     lu8rpc
*
*     Sparse LU update: Replace Column
*     LUSOL's sparse implementation of the Bartels-Golub update.
*
* 01 May 2002: Derived from LUSOL's original lu8a.f file.
* 01 May 2002: Current version of lusol8a.f.
* 15 Sep 2004: Test nout. gt. 0 to protect write statements.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lu8rpc( mode1, mode2, m, n, jrep, v, w,
     $                   lena, luparm, parmlu,
     $                   a, indc, indr, ip, iq,
     $                   lenc, lenr, locc, locr,
     $                   inform, diag, vnorm )

      implicit           double precision(a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), v(m), w(n)
      integer            indc(lena), indr(lena), ip(m), iq(n)
      integer            lenc(n), lenr(m)
      integer            locc(n), locr(m)

*     ------------------------------------------------------------------
*     lu8rpc  updates the LU factorization  A = L*U  when column  jrep
*     is replaced by some vector  a(new).
*
*     lu8rpc  is an implementation of the Bartels-Golub update,
*     designed for the case where A is rectangular and/or singular.
*     L is a product of stabilized eliminations (m x m, nonsingular).
*     P U Q is upper trapezoidal (m x n, rank nrank).
*
*     If  mode1 = 0,  the old column is taken to be zero
*                     (so it does not have to be removed from  U).
*
*     If  mode1 = 1,  the old column need not have been zero.
*
*     If  mode2 = 0,  the new column is taken to be zero.
*                     v(*)  is not used or altered.
*
*     If  mode2 = 1,  v(*)  must contain the new column  a(new).
*                     On exit,  v(*)  will satisfy  L*v = a(new).
*
*     If  mode2 = 2,  v(*)  must satisfy  L*v = a(new).
*
*     The array  w(*)  is not used or altered.
*
*     On entry, all elements of  locc  are assumed to be zero.
*     On a successful exit (inform ne 7), this will again be true.
*
*     On exit:
*     inform = -1  if the rank of U decreased by 1.
*     inform =  0  if the rank of U stayed the same.
*     inform =  1  if the rank of U increased by 1.
*     inform =  2  if the update seemed to be unstable
*                  (diag much bigger than vnorm).
*     inform =  7  if the update was not completed (lack of storage).
*     inform =  8  if jrep is not between 1 and n.
*
*     -- Jan 1985: Original F66 version.
*     -- Jul 1987: Modified to maintain U in trapezoidal form.
*     10 May 1988: First f77 version.
*     16 Oct 2000: Added test for instability (inform = 2).
*     ------------------------------------------------------------------

      logical            singlr
      parameter        ( zero = 0.0d+0 )

      nout   = luparm(1)
      lprint = luparm(2)
      nrank  = luparm(16)
      lenL   = luparm(23)
      lenU   = luparm(24)
      lrow   = luparm(25)
      Utol1  = parmlu(4)
      Utol2  = parmlu(5)
      nrank0 = nrank
      diag   = zero
      vnorm  = zero
      if (jrep .lt. 1) go to 980
      if (jrep .gt. n) go to 980

*     ------------------------------------------------------------------
*     If mode1 = 0, there are no elements to be removed from  U
*     but we still have to set  krep  (using a backward loop).
*     Otherwise, use lu7zap to remove column  jrep  from  U
*     and set  krep  at the same time.
*     ------------------------------------------------------------------
      if (mode1 .eq. 0) then
         krep   = n + 1

   10    krep   = krep - 1
         if (iq(krep) .ne. jrep) go to 10
      else
         call lu7zap( m, n, jrep, krep,
     $                lena, lenU, lrow, nrank,
     $                a, indr, ip, iq, lenr, locr )
      end if

*     ------------------------------------------------------------------
*     Insert a new column of u and find klast.
*     ------------------------------------------------------------------

      if (mode2 .eq. 0) then
         klast  = 0
      else
         if (mode2 .eq. 1) then

*           Transform v = a(new) to satisfy  L*v = a(new).

            call lu6sol( 1, m, n, v, w, lena, luparm, parmlu,
     $                   a, indc, indr, ip, iq,
     $                   lenc, lenr, locc, locr, inform )
         end if

*        Insert into  U  any nonzeros in the top of  v.
*        row  ip(klast)  will contain the last nonzero in pivotal order.
*        Note that  klast  will be in the range  ( 0, nrank ).

         call lu7add( m, n, jrep, v,
     $                lena, luparm, parmlu,
     $                lenL, lenU, lrow, nrank,
     $                a, indr, ip, lenr, locr,
     $                inform, klast, vnorm )
         if (inform .eq. 7) go to 970
      end if

*     ------------------------------------------------------------------
*     In general, the new column causes U to look like this:
*
*                 krep        n                 krep  n
*
*                ....a.........          ..........a...
*                 .  a        .           .        a  .
*                  . a        .            .       a  .
*                   .a        .             .      a  .
*        P U Q =     a        .    or        .     a  .
*                    b.       .               .    a  .
*                    b .      .                .   a  .
*                    b  .     .                 .  a  .
*                    b   ......                  ..a...  nrank
*                    c                             c
*                    c                             c
*                    c                             c     m
*
*     klast points to the last nonzero "a" or "b".
*     klast = 0 means all "a" and "b" entries are zero.
*     ------------------------------------------------------------------

      if (mode2 .eq. 0) then
         if (krep .gt. nrank) go to 900
      else if (nrank .lt. m) then

*        Eliminate any "c"s (in either case).
*        Row nrank + 1 may end up containing one nonzero.

         call lu7elm( m, n, jrep, v,
     $                lena, luparm, parmlu,
     $                lenL, lenU, lrow, nrank,
     $                a, indc, indr, ip, iq, lenr, locc, locr,
     $                inform, diag )
         if (inform .eq. 7) go to 970

         if (inform .eq. 1) then

*           The nonzero is apparently significant.
*           Increase nrank by 1 and make klast point to the bottom.

            nrank = nrank + 1
            klast = nrank
         end if
      end if

      if (nrank .lt. n) then

*        The column rank is low.
*
*        In the first case, we want the new column to end up in
*        position nrank, so the trapezoidal columns will have a chance
*        later on (in lu7rnk) to pivot in that position.
*
*        Otherwise the new column is not part of the triangle.  We
*        swap it into position nrank so we can judge it for singularity.
*        lu7rnk might choose some other trapezoidal column later.

         if (krep .lt. nrank) then
            klast     = nrank
         else
            iq(krep ) = iq(nrank)
            iq(nrank) = jrep
            krep      = nrank
         end if
      end if

*     ------------------------------------------------------------------
*     If krep .lt. klast, there are some "b"s to eliminate:
*
*                  krep
*
*                ....a.........
*                 .  a        .
*                  . a        .
*                   .a        .
*        P U Q =     a        .  krep
*                    b.       .
*                    b .      .
*                    b  .     .
*                    b   ......  nrank
*
*     If krep .eq. klast, there are no "b"s, but the last "a" still
*     has to be moved to the front of row krep (by lu7for).
*     ------------------------------------------------------------------

      if (krep .le. klast) then

*        Perform a cyclic permutation on the current pivotal order,
*        and eliminate the resulting row spike.  krep becomes klast.
*        The final diagonal (if any) will be correctly positioned at
*        the front of the new krep-th row.  nrank stays the same.

         call lu7cyc( krep, klast, ip )
         call lu7cyc( krep, klast, iq )

         call lu7for( m, n, krep, klast,
     $                lena, luparm, parmlu,
     $                lenL, lenU, lrow,
     $                a, indc, indr, ip, iq, lenr, locc, locr,
     $                inform, diag )
         if (inform .eq. 7) go to 970
         krep   = klast

*        Test for instability (diag much bigger than vnorm).

         singlr = vnorm .lt. Utol2 * abs( diag )
         if ( singlr ) go to 920
      end if

*     ------------------------------------------------------------------
*     Test for singularity in column krep (where krep .le. nrank).
*     ------------------------------------------------------------------

      diag   = zero
      iw     = ip(krep)
      singlr = lenr(iw) .eq. 0

      if (.not. singlr) then
         l1     = locr(iw)
         j1     = indr(l1)
         singlr = j1 .ne. jrep

         if (.not. singlr) then
            diag   = a(l1)
            singlr = abs( diag ) .le. Utol1          .or.
     $               abs( diag ) .le. Utol2 * vnorm
         end if
      end if

      if ( singlr  .and.  krep .lt. nrank ) then

*        Perform cyclic permutations to move column jrep to the end.
*        Move the corresponding row to position nrank
*        then eliminate the resulting row spike.

         call lu7cyc( krep, nrank, ip )
         call lu7cyc( krep, n    , iq )

         call lu7for( m, n, krep, nrank,
     $                lena, luparm, parmlu,
     $                lenL, lenU, lrow,
     $                a, indc, indr, ip, iq, lenr, locc, locr,
     $                inform, diag )
         if (inform .eq. 7) go to 970
      end if

*     Find the best column to be in position nrank.
*     If singlr, it can't be the new column, jrep.
*     If nothing satisfactory exists, nrank will be decreased.

      if ( singlr  .or.  nrank .lt. n ) then
         jsing  = 0
         if ( singlr ) jsing = jrep

         call lu7rnk( m, n, jsing,
     $                lena, luparm, parmlu,
     $                lenL, lenU, lrow, nrank,
     $                a, indc, indr, ip, iq, lenr, locc, locr,
     $                inform, diag )
      end if

*     ------------------------------------------------------------------
*     Set inform for exit.
*     ------------------------------------------------------------------

  900 if (nrank .eq. nrank0) then
         inform =  0
      else if (nrank .lt. nrank0) then
         inform = -1
         if (nrank0 .eq. n) then
            if (nout. gt. 0  .and.  lprint .ge. 0)
     &           write(nout, 1100) jrep, diag
         end if
      else
         inform =  1
      end if
      go to 990

*     Instability.

  920 inform = 2
      if (nout. gt. 0  .and.  lprint .ge. 0)
     &     write(nout, 1200) jrep, diag
      go to 990

*     Not enough storage.

  970 inform = 7
      if (nout. gt. 0  .and.  lprint .ge. 0)
     &     write(nout, 1700) lena
      go to 990

*     jrep  is out of range.

  980 inform = 8
      if (nout. gt. 0  .and.  lprint .ge. 0)
     &     write(nout, 1800) m, n, jrep

*     Exit.

  990 luparm(10) = inform
      luparm(15) = luparm(15) + 1
      luparm(16) = nrank
      luparm(23) = lenL
      luparm(24) = lenU
      luparm(25) = lrow
      return

 1100 format(/ ' lu8rpc  warning.  Singularity after replacing column.',
     $       '    jrep =', i8, '    diag =', 1p, e12.2 )
 1200 format(/ ' lu8rpc  warning.  Instability after replacing column.',
     $       '    jrep =', i8, '    diag =', 1p, e12.2 )
 1700 format(/ ' lu8rpc  error...  Insufficient storage.',
     $         '    lena =', i8)
 1800 format(/ ' lu8rpc  error...  jrep  is out of range.',
     $         '    m =', i8, '    n =', i8, '    jrep =', i8)

      end ! subroutine lu8rpc
