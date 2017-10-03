c----------------------------------------------------------------------
      subroutine usr_particles_init
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTDATA'
      include 'CMTPART'

      logical ifreadpart

c     zero before anything happens
      call rzero(ptdum,iptlen)
      call rzero(pttime,iptlen)

      ! begin timer
      ptdum(1) = dnekclock()

      nr   = lr     ! Mandatory for proper striding
      ni   = li     ! Mandatory
      nrgp = lrgp
      nigp = ligp
      nrf  = lrf
      nif  = lif

      call rzero(rpart,lr*llpart)
      call izero(ipart,li*llpart)

      call read_particle_input(ifreadpart)
      call set_bounds_box
      call set_part_pointers
      if(.not.ifreadpart) call place_particles     ! n initialized here
      call output_particle_options
      call move_particles_inproc          ! initialize fp & cr comm handles
         ntmp  = iglsum(n,1)
         if (nid.eq.0) write(6,*) 'Passed move_particles_inproc'
      call init_interpolation ! barycentric weights for interpolation
         ntmp  = iglsum(n,1)
         if (nid.eq.0) write(6,*) 'Passed init_interpolation'
      if (two_way.eq.1) then
         call compute_neighbor_el_proc    ! compute list of neigh. el. ranks 
            ntmp  = iglsum(n,1)
            if (nid.eq.0) write(6,*) 'Passed compute_neighbor_el_proc'
         call particles_solver_nearest_neighbor ! nearest neigh
            ntmp  = iglsum(n,1)
            if (nid.eq.0) write(6,*) 'Passed particles_solver_nearest_n'
         call point_to_grid_corr_init    ! for gamma correction integrat
            ntmp  = iglsum(n,1)
            if (nid.eq.0) write(6,*) 'Passed point_to_grid_corr_init'
         call spread_props_grid           ! put particle props on grid
            ntmp  = iglsum(n,1)
            if (nid.eq.0) write(6,*) 'Passed spread_props_grid'

         do i = 2,nitspl
            call interp_props_part_location ! interpolate
            call correct_spl
            call particles_solver_nearest_neighbor ! nearest neigh
            call spread_props_grid           ! put particle props on grid
            if (nid.eq.0) write(6,*) i,'Pre-SPL iteration'
         enddo
      endif

c     rxv = rxbo(2,1) - 2.*(0.039/5.)
c     rxv = -0.4

c     rdumv = 0.
c     do i=1,n
c        rdumv = rdumv + rpart(jvol,i)
c     enddo
c     rdumt = glsum(rdumv,1)
c        do i=1,n
c           if (rpart(jx,i) .gt. rxv) then
c              rpart(jspl,i) = phi_desire*vol_distrib/rdumt
c              rpart(jspl,i) = 1. ! decrease with distance!!
c           endif
c        enddo
c           call spread_props_grid           ! put particle props on grid

c     red_interp = (lx1+1)/2 ! 0         = full spectral interpolation 
c     call init_interpolation ! barycentric weights for interpolation
c        ntmp  = iglsum(n,1)
c        if (nid.eq.0) write(6,*) 'Passed reinit_interpolation'

      ! end timer 
      pttime(1) = pttime(1) + dnekclock() - ptdum(1)

      ntmp  = iglsum(n,1)
      if (nid.eq.0) write(6,*) 'Passed usr_particles_init'

      return
      end
c----------------------------------------------------------------------
      subroutine set_bounds_box
c
c     set domain and element bounds for a box geometry. Notice that
c     this ONLY works with non curved elements.
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTTIMERS'
      include 'CMTPART'

      real   xdrange(2,3)
      common /domainrange/ xdrange
      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      ! begin timer
      ptdum(2) = dnekclock()

      if(istep.eq.0.or.istep.eq.1)then
        call domain_size(xdrange(1,1),xdrange(2,1),xdrange(1,2)
     $                  ,xdrange(2,2),xdrange(1,3),xdrange(2,3))
        ntot = lx1*ly1*lz1*nelt
        nxyz = lx1*ly1*lz1
        do ie = 1,nelt
           xerange(1,1,ie) = vlmin(xm1(1,1,1,ie),nxyz)
           xerange(2,1,ie) = vlmax(xm1(1,1,1,ie),nxyz)
           xerange(1,2,ie) = vlmin(ym1(1,1,1,ie),nxyz)
           xerange(2,2,ie) = vlmax(ym1(1,1,1,ie),nxyz)
           xerange(1,3,ie) = vlmin(zm1(1,1,1,ie),nxyz)
           xerange(2,3,ie) = vlmax(zm1(1,1,1,ie),nxyz)
        enddo  
      endif

      ! end timer
      pttime(2) = pttime(2) + dnekclock() - ptdum(2)

      return
      end
c----------------------------------------------------------------------
      subroutine place_particles
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange
      real   xdrange(2,3)
      common /domainrange/ xdrange

      integer icalld
      save    icalld
      data    icalld  /-1/

      integer nwe
      real    dum

      ! begin timer
      ptdum(3) = dnekclock()

c     setup items
      pi    = 4.0*atan(1.) ! pi
      nlxyze = nx1*ny1*nz1*nelt
      call rzero(ptw,nlxyze*8)

c     region to distribute particles
c     now done in particles.ini file

c     some computed parameters
      nwe         = int(nw/np)                ! num. part per proc
      deltax      = (xerange(2,1,1) - xerange(1,1,1))/nx1
      deltaf      = df_dx*deltax                 ! gaussian filter half
c     dp          = deltaf/df_dp                 ! particle diameter
         df_dp       = deltaf/((dp(1)+dp(2))/2.) ! average
      rsig        = deltaf/(2.*sqrt(2.*log(2.))) ! gaussian filter std.

      vol_distrib = (rxbo(2,1)-rxbo(1,1))*(rxbo(2,2)-rxbo(1,2))*
     >              (rxbo(2,3)-rxbo(1,3))

c     correct nwe if discrepancy
      nw_tmp      = iglsum(nwe,1)
      if ((nw_tmp .ne. nw) .and. (nid.eq.0)) nwe = nwe + (nw - nw_tmp)

c     main loop to distribute particles
      if (icalld .lt. 0) then
         rdum   = ran2(-nrandseed*np-nid-1) ! initialize random number generator
         n      = 0
         icalld = icalld + 1
      endif

      do i = 1,nwe
         n = n + 1
         if (n.gt.llpart)then 
            write(6,*)'Not enough space to store more particles'
            call exitt
         endif

         do j=0,2
            rval = unif_random(rxbo(1,j+1),rxbo(2,j+1))
            rpart(jx+j,n)  = rval
            rpart(jx1+j,n) = rval
            rpart(jx2+j,n) = rval
            rpart(jx3+j,n) = rval
         enddo

c        set some rpart values for later use
         rpart(jdp,n)   = unif_random(dp(1),dp(2)) ! particle diameter
         tau_p          = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0  ! part. time scale stokes


         tau_p = 10.*dt
         rho_p = 18.*mu_0*tau_p/rpart(jdp,n)**2


         rpart(jtaup,n) = tau_p      ! particle time scale
         rpart(jrhop,n) = rho_p      ! material density of particle
         rpart(jvol,n)  = pi*rpart(jdp,n)**3/6.! particle volume
         rpart(jspl,n)  = 1.         ! super particle loading
         if (nitspl.gt.0) then
         rdumv = pi*((dp(1)+dp(2))/2.)**3/6.
         if (two_way.eq.1) then
c           rpart(jspl,n) =  phi_desire*vol_distrib/(nw*rpart(jvol,n))
c           rpart(jspl,n) = phi_desire*vol_distrib/(nw*rdumv)
         endif
         endif
         rpart(jgam,n)  = 1.          ! initial integration correction

         rpart(jtemp,n)  = 294. ! intial temp as fluid air
         rpart(jtempf,n) = 294. ! intial temp as fluid air

c        set global particle id (3 part tag)
         ipart(jpid1,n) = nid 
         ipart(jpid2,n) = n 
         ipart(jpid3,n) = istep
      enddo


      if (nitspl.gt.0) then
      if (istep.eq. 0 .or.istep .eq.1) then ! be careful doing this in this way with injecting
      rdumv = 0.
      do i=1,n
         rdumv = rdumv + rpart(jvol,i)
      enddo
      rdumt = glsum(rdumv,1)
      if (two_way.eq.1) then
      do i=1,n
            rpart(jspl,i) =  phi_desire*vol_distrib/(rdumt)
      enddo
      endif
      endif
      endif

c     check if zstart and zlen is alright for a 2d case
      if (.not. if3d) then
          if (abs(zstart-1.0) .gt. 1E-16) then
             write(6,*)'***particle zstart is not right for 2d case'
             call exitt
          elseif(abs(zlen) .gt. 1E-16) then
             write(6,*)'***particle zlen is not right for 2d case'
             call exitt
         endif
      endif

      ! end timer
      pttime(3) = pttime(3) + dnekclock() - ptdum(3)
      return
      end
c-----------------------------------------------------------------------
      function unif_random(rxl,rxr)
c     must initialize ran2 first
      real xl,xr,unif_random

      rdum       = ran2(2)
      rlen       = rxr - rxl
      unif_random= rxl + rdum*rlen

      return
      end
c-----------------------------------------------------------------------
      function unif_random_norm(rxl,rxr,rstd)
c     must initialize ran2 first
      real xl,xr,unif_random_norm,rstd,rxfne(1000),rcdf(1000)

      rmu = (rxr + rxl)/2.
      nxfn  = 1000
      rxlf  = rmu - 5.*rstd
      rxrf  = rmu + 5.*rstd
      rdxf  = (rxrf-rxlf)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = rxlf + (i-1.)*rdxf
         rcdf(i)  = 0.5*(1. + erf((rxfne(i)-rmu)/(rstd*sqrt(2.))))
      enddo

      rdum = unif_random(0.,1.)

!     find lower min value for inverse sampling
      idum = 0
      rmin = 100.
      do i=1,nxfn
         if (abs(rdum - rcdf(i)) .lt. rmin) then
            rmin = abs(rdum -rcdf(i))
            idum = i
         endif
      enddo
      ml = idum
      if (rdum .lt. rcdf(idum)) ml = ml + 1

      if (rdum .gt. rcdf(nxfn)) then
         unif_random_norm = rxrf 
      elseif (rdum .lt. rcdf(1)) then
         unif_random_norm = rxlf
      else
         rm = (rxfne(ml+1) - rxfne(ml))/(rcdf(ml+1) - rcdf(ml))
         unif_random_norm = rxfne(ml) + rm*(rdum - rcdf(ml))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_part_pointers
      include 'SIZE'
      include 'CMTPART'

      ! begin timer
      ptdum(4) = dnekclock()

c     ipart pointers ------------------------------------------------
      jrc   = 1 ! Pointer to findpts return code
      jpt   = 2 ! Pointer to findpts return processor id
      je0   = 3 ! Pointer to findpts return element id
      je00  = 4 ! Pointer to findpts return element id
      jps   = 5 ! Pointer to proc id for data swap
      jpid1 = 6 ! initial proc number
      jpid2 = 7 ! initial local particle id
      jpid3 = 8 ! initial time step introduced
      jpnn  = 9 ! initial time step introduced
      jpid  = 10 ! initial time step introduced
      jrco  = 11 ! initial time step introduced
      jai   = 12 ! Pointer to auxiliary integers

      nai = ni - (jai-1)  ! Number of auxiliary integers
      if (nai.le.0) call exitti('Error in nai:$',ni)

c     rpart pointers ------------------------------------------------
      jr  = 1         ! Pointer to findpts return rst variables
      jd  = jr + 3    ! Pointer to findpts return distance
      jx  = jd + 1    ! Pointer to findpts input x value
      jy  = jx + 1    ! Pointer to findpts input y value
      jz  = jy + 1    ! Pointer to findpts input z value
      jv0 = jz + 1    ! particle velocity at this timestep
      ju0 = jv0 + 3   ! fluid velocity at this time step
      jf0 = ju0 + 3   ! particle total force at this timestep
      jq0 = jf0 + 3   ! temperature forcing
      jg0 = jq0 + 1   ! work done by forces
      jquu= jg0 + 1   ! undisturbed unsteady temp. forcing
      jqqs= jquu+ 1   ! quasi-steady temp. forcing

c     forcing
      ii  = jqqs + 1
c     if (part_force(1).ne.0) then ! user specified force
         jfusr = ii
         ii    = ii + 3
c     endif
c     if (part_force(2).ne.0) then ! quasi-steady force
         jfqs  = ii
         ii    = ii + 3
c     endif
c     if (part_force(3).ne.0) then ! undisturbed force
         jfun  = ii
         ii    = ii + 3
c     endif
c     if (part_force(4).ne.0) then ! inviscid unsteady force
         jfiu  = ii
         ii    = ii + 3
c     endif

c     other parameters (some may not be used; all at part. location)
      jtaup   = ii          ! particle time scale
      jcd     = jtaup   + 1 ! drag coeff
      jdrhodt = jcd     + 3 ! density material time derivative
      jre     = jdrhodt + 1 ! Relative Reynolds number
      jDuDt   = jre     + 1 ! fluid velocity time derivative
      jtemp   = jDuDt   + 3 ! part. temperature 
      jtempf  = jtemp   + 1 ! fluid temperature at part. loc.
      jrho    = jtempf  + 1 ! fluid denisty 
      jrhop   = jrho    + 1 ! particle material density
      ja      = jrhop   + 1 ! fluid mach number
      jvol    = ja      + 1 ! particle volume 
      jvol1   = jvol    + 1 ! particle volume fraction at part. loc.
      jdp     = jvol1   + 1 ! particle diameter
      jgam    = jdp     + 1 ! spread to grid correction
      jspl    = jgam    + 1 ! super particle loading
      jcmiu   = jspl    + 1 ! added mass coefficient

c     bdf/ext integration
      jx1 = jcmiu+1 ! Pointer to xyz at t^{n-1}
      jx2 = jx1 +3 ! Pointer to xyz at t^{n-1}
      jx3 = jx2 +3 ! Pointer to xyz at t^{n-1}

      jv1 = jx3+ 3 ! Pointer to particle velocity at t^{n-1}
      jv2 = jv1+ 3 ! Pointer to particle velocity at t^{n-2}
      jv3 = jv2+ 3 ! Pointer to particle velocity at t^{n-3}

      ju1 = jv3+ 3 ! Pointer to fluid velocity at t^{n-1}
      ju2 = ju1+ 3 ! Pointer to fluid velocity at t^{n-2}
      ju3 = ju2+ 3 ! Pointer to fluid velocity at t^{n-3}

      jar = ju3+ 3 ! Pointer to auxiliary reals

      nar = nr - (jar-1)  ! Number of auxiliary reals
      if (nar.le.0) call exitti('Error in nar:$',nr)

c     ghost particle integer pointers -------------------------------
      jgppid1 = 1 ! initial proc number
      jgppid2 = 2 ! initial local particle id
      jgppid3 = 3 ! initial time step introduced
      jgpps   = 4 ! Pointer to proc id for data swap
      jgppt   = 5 ! findpts return processor id
      jgpes   = 6 ! Destination element to be sent to

c     ghost particle real pointers ----------------------------------
      jgpx    = 1 ! ghost particle xloc
      jgpy    = 2 ! ghost particle yloc
      jgpz    = 3 ! ghost particle zloc
      jgpfh   = 4 ! ghost particle hydrodynamic xforce (i+1 > y, i+2 > z)
      jgpvol  = jgpfh+3  ! ghost particle volume
      jgpgam  = jgpvol+1 ! spreading correction (if used)
      jgpspl  = jgpgam+1 ! super particle loading
      jgpg0   = jgpspl+1 ! 
      jgpq0   = jgpg0 +1 ! 
      jgpv0   = jgpq0 +1 ! velocity (3 components)

      ! end timer
      pttime(4) = pttime(4) + dnekclock() - ptdum(4)

      return
      end
c----------------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
c between 0.0 and 1.0 (exclusive of the endpoint values). 
c Call with idum a negative integer to initialize; thereafter, do not alter 
c idum between successive deviates in a sequence. RNMX should approximate the 
c largest floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then 
         idum1=max(-idum,1) 
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1 
            if (idum1.lt.0) idum1=idum1+IM1 
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1) 
      endif
      k=idum1/IQ1 
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1 
      k=idum2/IQ2 
      idum2=IA2*(idum2-k*IQ2)-k*IR2 
      if (idum2.lt.0) idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1 
      if(iy.lt.1)iy=iy+IMM1 
      ran2=min(AM*iy,RNMX)
      return
      END
c----------------------------------------------------------------------
