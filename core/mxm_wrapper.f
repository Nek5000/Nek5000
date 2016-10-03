      subroutine mxm(a,n1,b,n2,c,n3)
c
c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.
c
      real a(n1,n2),b(n2,n3),c(n1,n3)
c
      include 'SIZE'
      include 'CTIMER'
      include 'OPCTR'
      include 'TOTAL'
c
      integer*8 tt

#ifdef TIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'mxm   '
      endif
      isbcnt = n1*n3*(2*n2-1)
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount = dcount + (isbcnt)
      etime1 = dnekclock()
#endif


#ifdef BGQ
      if (n2 == 8 .and. mod(n1,4) == 0 &
c        .and. MOD(LOC(a),tt)==0 & 
c        .and. MOD(LOC(b),tt)==0 & 
c        .and. MOD(LOC(c),tt)==0 & 
        ) then
        call mxm_bgq_8(a, n1, b, n2, c, n3)  
        goto 111
      endif
      if (n2 == 16 .and. mod(n1,4) == 0 &
c        .and. MOD(LOC(a),tt)==0 & 
c        .and. MOD(LOC(b),tt)==0 & 
c        .and. MOD(LOC(c),tt)==0 & 
         ) then
        call mxm_bgq_16(a, n1, b, n2, c, n3)  
        goto 111
      endif
      tt = 32
      if (n2 == 10 .and. mod(n1,4) == 0 .and. mod(n3,2) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
          ) then
        call mxm_bgq_10(a, n1, b, n2, c, n3)  
        goto 111
      endif
      if (n2 == 6 .and. mod(n1,4) == 0 .and. mod(n3,2) == 0 &
        .and. MOD(LOC(a),tt)==0 & 
        .and. MOD(LOC(b),tt)==0 & 
        .and. MOD(LOC(c),tt)==0 & 
       ) then
        call mxm_bgq_6(a, n1, b, n2, c, n3)  
        goto 111
      endif
#endif

#ifdef XSMM
      call libxsmm_dgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)
      goto 111
#endif

#ifdef BLAS_MXM
      call dgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)
      goto 111
#endif

      call mxmf2(a,n1,b,n2,c,n3)

 111  continue
#ifdef TIMER
      tmxmf = tmxmf + dnekclock() - etime1  
#endif
      return
      end
