c-----------------------------------------------------------------------
      subroutine makeq

C     Generate forcing function for the solution of a passive scalar.
C     !! NOTE: Do not change the content of the array BQ until the current
C              time step is completed 

      include 'SIZE'
      include 'TOTAL'

      logical  ifturb,if_conv_std
      common /SCRUZ/ wrk(lx1,ly1,lz1,lelt)

      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f

      call whatfld (ifturb)

      if (ifturb)                                          call maketq
      if (.not.ifturb                   .and.if_conv_std)  call makeuq
      if (ifadvc(ifield).and..not.ifchar.and.if_conv_std)  call convab

      if (iftran) then
         if(ifcvode) then
           n = nx1*ny1*nz1*nelfld(ifield)
           call wlaplacian (wrk,t(1,1,1,1,ifield-1),
     &                      vdiff(1,1,1,1,ifield),ifield)
           call add2    (bq(1,1,1,1,ifield-1),wrk,ntot)
           call invcol2 (bq(1,1,1,1,ifield-1),vtrans(1,1,1,1,ifield),n)
         else
           if (ifmvbd) then       ! ifchar is false
              call admesht
              call makeabq
              call makebdq
           elseif (ifchar.and.ifadvc(ifield)) then
              call makeabq
              call convch
           else
              call makeabq
              call makebdq
           endif
         endif
      endif

      return
      end

