c-----------------------------------------------------------------------
      subroutine makeq

C     Generate forcing function for the solution of a passive scalar.
C     !! NOTE: Do not change the content of the array BQ until the current
C              time step is completed 

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      logical  ifturb,if_conv_std

      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f

      call whatfld (ifturb)

      if (ifturb)                                          call maketq
      if (.not.ifturb                   .and.if_conv_std)  call makeuq
      if (ifadvc(ifield).and..not.ifchar.and.if_conv_std)  call convab

      if (iftran.and..not.ifcvode) then

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

      return
      end
c-----------------------------------------------------------------------
