c----------------------------------------------------------------------
      subroutine lpm_stokes_particles_solver
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer             stage,nstage
      common /tstepstage/ stage,nstage

      nstage_part = 3
      if (abs(time_integ) .eq. 2) nstage_part = 1

      if (istep.eq.0) then
         call lpm_init(0)
      else

         call set_tstep_coef_part(dt) ! in nek5000 with rk3
         do stage=1,nstage_part
            call lpm_usr_particles_solver
         enddo
         call compute_phig_qtl(dt,usrdiv) ! nek5000 (see Zwick 2018)
      endif

      if(mod(istep,iostep).eq.0.or. istep.eq.1) then
         call lpm_usr_particles_io
      endif

      return
      end
c----------------------------------------------------------------------
c     setup routines
c----------------------------------------------------------------------
      subroutine lpm_init(idum)
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LPM'

c     common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
c     integer gfirst, inoassignd, resetFindpts, pload

      icmtp = idum

      ! for 2d cases
      zlpart   = 0.0
      zrpart   = dp(1)

      call set_part_pointers
      call set_bounds_box
      call set_part_params 
      call place_particles
      call move_particles_inproc

c     call move_particles_inproc          
      if (two_way.gt.1) then
c        call compute_neighbor_el_proc    
         call create_extra_particles
         call send_ghost_particles
         call spread_props_grid           
      endif
      call init_interpolation
      call interp_props_part_location 

      if (time_integ .lt. 0) call pre_sim_collisions ! e.g., settling p
   
c     resetFindpts = 0
c     call computeRatio
c     call reinitialize
c     call printVerify

      return
      end
c----------------------------------------------------------------------
      subroutine place_particles
c
c     Place particles in this routine, also called for injection
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer icalld
      save    icalld
      data    icalld  /-1/

      real     unif_random,unif_random_norm,unif_random_cyl
     >        ,unif_random_sphere
      external unif_random,unif_random_norm,unif_random_cyl,
     >         unif_random_sphere

      icalld = icalld + 1

      if (ipart_restartr .eq. 0) then

c        correct nwe if discrepancy on rank 0
         nwe         = int(nw/np)                ! num. part per proc
         nw_tmp      = iglsum(nwe,1)
         ndef        = nw - nw_tmp
         if (nid .lt. ndef) nwe = nwe + 1

c        main loop to distribute particles
         do i_pt_part = 1,nwe
            n = n + 1
  754 continue
            ! Error checking
            if (n.gt.llpart)then 
               if (nid.eq.0)
     >            write(6,*)'Not enough space to store more particles'
               call exitt
            endif

            ! sphere or cylinder
            if (rxco(1) .gt. -1E7) then

               ! distribute in cylinder 
               if (rxco(6) .gt. -1E7) then
               rnx = rxco(1)
               rny = rxco(2)
               rnz = rxco(3)

               rx0 = rxco(4)
               ry0 = rxco(5)
               rz0 = 1.0
               if (if3d) rz0 = rxco(6)

               rmag = sqrt(rnx**2 + rny**2 + rnz**2)
               rnx = rnx/rmag
               rny = rny/rmag
               rnz = rnz/rmag

               rin  = rxco(7)
               rout = rxco(8)
               rht  = 0.
               if (if3d) rht  = rxco(9)

               rrad = unif_random_cyl(rin,rout)
               rthet= unif_random(0.,2.*pi)
               rxtr = unif_random(-rht/2.,rht/2.)

               do j=0,2
               ! x cylinder
               if ( abs(rnx).gt.abs(rny).and.abs(rnx).gt.abs(rnz)) then
                  if (j.eq. 0) rdum = rx0 + rxtr
                  if (j.eq. 1) rdum = ry0 + rrad*cos(rthet)
                  if (j.eq. 2) rdum = rz0 + rrad*sin(rthet)
               endif
               ! y cylinder
               if ( abs(rny).gt.abs(rnx).and.abs(rny).gt.abs(rnz)) then
                  if (j.eq. 0) rdum = rx0 + rrad*sin(rthet)
                  if (j.eq. 1) rdum = ry0 + rxtr
                  if (j.eq. 2) rdum = rz0 + rrad*cos(rthet)
               endif
               ! z cylinder
               if ( abs(rnz).gt.abs(rnx).and.abs(rnz).gt.abs(rny)) then
                  if (j.eq. 0) rdum = rx0 + rrad*cos(rthet)
                  if (j.eq. 1) rdum = ry0 + rrad*sin(rthet)
                  if (j.eq. 2) rdum = rz0 + rxtr
               endif
                  rpart(jx +j,n) = rdum
                  rpart(jx1+j,n) = rdum
                  rpart(jx2+j,n) = rdum
                  rpart(jx3+j,n) = rdum
               enddo

               ! distribute sphere
               else
               rx0 = rxco(1)
               ry0 = rxco(2)
               rz0 = 1.0
               if (if3d) rz0 = rxco(3)
               rin  = rxco(4)
               rout = rxco(5)

               rrad  = unif_random_sphere(rin,rout)
               rthet1= unif_random(0.,1.)
               rthet1= acos(2.*rthet1 - 1) !tricky tricky
               rthet2= unif_random(0.,2.*pi)

               do j=0,2
                  if (j.eq. 0) rdum = rx0 + rrad*sin(rthet1)*cos(rthet2)
                  if (j.eq. 1) rdum = ry0 + rrad*sin(rthet1)*sin(rthet2)
                  if (j.eq. 2) rdum = rz0 + rrad*cos(rthet1)
                  rpart(jx +j,n) = rdum
                  rpart(jx1+j,n) = rdum
                  rpart(jx2+j,n) = rdum
                  rpart(jx3+j,n) = rdum
               enddo
               endif

               ! check if box defined and if it is, particles must fit 
               ! in that box even if spere or cylinder distribution.
               ! If not, try placing again
               if (rxbo(1,1) .gt. -1E7) then
                  do j=0,ldim-1
                     if ((rpart(jx+j,n) .gt. rxbo(2,j+1)) .or.
     >                   (rpart(jx+j,n) .lt. rxbo(1,j+1))) then
                        goto 754
                     endif
                  enddo
               endif
                  
            ! distribute in box
            else
               do j=0,2
                  rpart(jx +j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
                  rpart(jx1+j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
                  rpart(jx2+j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
                  rpart(jx3+j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
               enddo
            endif

c           set some rpart values for later use
            if (dp_std .gt. 0) then
               rpart(jdp,n) = unif_random_norm(dp(1),dp_std)
            else
               rpart(jdp,n) = unif_random(dp(1),dp(2))
            endif
            rpart(jtaup,n) = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0  ! particle time scale
            rpart(jrhop,n) = rho_p                               ! particle density 
            rpart(jvol,n) = pi*rpart(jdp,n)**3/6.      ! particle volume
            rpart(jspl,n)  = rspl                                ! super particle loading
            rpart(jrpe,n) = rpart(jspl,n)**(1./3.)*rpart(jdp,n)/2.
            rpart(jvol,n) = rpart(jspl,n)*rpart(jvol,n)
         
            rpart(jtemp,n)  = tp_0                               ! intial particle temp
            rpart(jtempf,n) = tp_0                               ! intial fluid temp (overwritten)
            rpart(jrho,n)   = param(1)                           ! initial fluid density (overwritten interp)
         
c           set global particle id (3 part tag)
            ipart(jpid1,n) = nid 
            ipart(jpid2,n) = i_pt_part
            ipart(jpid3,n) = icalld
         enddo

      elseif (ipart_restartr .gt. 0) then
         ! read in data
         nread_part = 1
         do j=1,nread_part
            call read_parallel_restart_part
         enddo

      elseif (ipart_restartr .lt. 0) then

         call read_ascii_restart_part
      endif


      ! force 2d z to be 1
      if (.not. if3d) then
         do i=1,n
            rpart(jx +2,i) = 1.
            rpart(jx1+2,i) = 1.
            rpart(jx2+2,i) = 1.
            rpart(jx3+2,i) = 1.

            rpart(ju0+2,i) = 0.0
            rpart(ju1+2,i) = 0.0
            rpart(ju2+2,i) = 0.0
            rpart(ju3+2,i) = 0.0
         enddo
      endif

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
      include 'LPM'

c     common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
c     integer gfirst, inoassignd, resetFindpts, pload


      if(istep.eq.0.or.istep.eq.1)then
c     if((istep.eq.0) .or. (istep.eq.1).or.(resetFindpts.eq.1)) then 
        call domain_size(xdrange(1,1),xdrange(2,1),xdrange(1,2)
     $                  ,xdrange(2,2),xdrange(1,3),xdrange(2,3))

       if (.not. if3d) then
          xdrange(1,3) = zlpart
          xdrange(2,3) = zrpart
       endif
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine set_part_params
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*132 deathmessage

      if (icalld .lt. 0) then
         rdum   = ran2(-nrandseed*np-nid-1) ! initialize
         icalld = icalld + 1
      endif

c     setup items
      nlxyze = nx1*ny1*nz1*nelt
      call rzero(rpart,lr*llpart)
      call izero(ipart,li*llpart)
      call rzero(ptw,nlxyze*7)
      call rzero(ptdum,iptlen)
      call rzero(pttime,iptlen)

      ! do nothing, no spreading
      if (npro_method .eq. 0) then
         d2chk(2) = 0.0
      ! gaussian set by user input parameters
      elseif (npro_method .eq. 1) then
         rsig = dfilt/(2.*sqrt(2.*log(2.))) ! gaussian filter std. * DP
         d2chk(2) = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))
      endif

      mu_0   = abs(param(2))

      return
      end
c----------------------------------------------------------------------
      subroutine set_check_spl_params
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      character*132 deathmessage

      rdeff_max = dp(2)/2.
      do i = 1,n
         if (rpart(jrpe,i) .gt. rdeff_max) rdeff_max=rpart(jrpe,i)
      enddo
      rdeff_max = glmax(rdeff_max,1)
      rdeff_max = rdeff_max*2. ! to get diameter

      d2chk(2) = d2chk(2)*rdeff_max
      d2chk(3) = 0.6     *rdeff_max

      d2chk(1) = max(d2chk(2),d2chk(3))

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_set_dt(rdt_part)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real dt_dum,dt_col,cflp,cflt

      common /save_dt_part/ rdpe_max, rdpe_min

      integer icalld
      save    icalld
      data    icalld  /-1/

      if (istep .gt. 1 .or. time_integ .lt. 0) then
         icalld = icalld + 1
c        if (icalld .le. 20) then
            rdpe_max = 0.
            rdpe_min = 100.
            do i=1,n
               if (rpart(jrpe,i).lt.rdpe_min) rdpe_min = rpart(jrpe,i)
               if (rpart(jrpe,i).gt.rdpe_max) rdpe_max = rpart(jrpe,i)
            enddo
            rdpe_min = glmin(rdpe_min,1)
            rdpe_max = glmax(rdpe_max,1)
            rdpe_min = 2.*rdpe_min ! to get diameter
            rdpe_max = 2.*rdpe_max ! to get diameter
c        endif

c       ! particle cfl, particles cant move due to velocity
        dt_dum = rdt_part
        dt_part = 1000.
        cflp = 0.10
        do i=1,n
           rvmag  = sqrt(rpart(jv0,i)**2 + rpart(jv0+1,i)**2 
     >                   + rpart(jv0+2,i)**2)
        
           cflt = dt_dum*rvmag/rdpe_min
           if (cflt .lt. cflp) dt_part = dt_dum
           if (cflt .ge. cflp) dt_part = cflp*rdpe_min/rvmag ! make sure smallest small overlap
        enddo
        dt_part  = glmin(dt_part,1)
        rdt_part = min(dt_dum,dt_part)
        if (dt_part .le. 1E-14) rdt_part = dt_dum

   
         if (two_way .gt. 2) then
            ! resolving collisions
            rm1      = rho_p*pi/6.*rdpe_min**3 ! max
            rm2      = rho_p*pi/6.*rdpe_max**3 ! max
            rm12     = 1./(1./rm1 + 1./rm2)
            n_resolve= 10
            dt_col   = sqrt(rm12/ksp*(log(e_rest)**2+pi**2))/n_resolve
            rdt_part = min(rdt_part,dt_col)
         endif
      endif

      

      return
      end
c----------------------------------------------------------------------
      subroutine set_part_pointers
      include 'SIZE'
      include 'LPM'

      nr   = lr     ! Mandatory for proper striding
      ni   = li     ! Mandatory
      nrgp = lrgp
      nigp = ligp
      nrf  = lrf
      nif  = lif
      n    = 0

c     ipart pointers ------------------------------------------------
      jrc   = 1 ! Pointer to findpts return code
      jpt   = 2 ! Pointer to findpts return processor id
      je0   = 3 ! Pointer to findpts return element id
      jps   = 4 ! Pointer to proc id for data swap
      jpid1 = 5 ! initial proc number
      jpid2 = 6 ! initial local particle id
      jpid3 = 7 ! initial time step introduced
      jicx  = 8 ! initial time step introduced
      jicy  = 9 ! initial time step introduced
      jicz  = 10 ! initial time step introduced
      jicn  = 11 ! initial time step introduced
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

      jfcol  = ii  ! collisional force

c     other parameters (some may not be used; all at part. location)
      jtaup   = jfcol   + 3 ! particle time scale
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
      jrpe    = jdp     + 1 ! particle effective diameter spl
      jgam    = jrpe    + 1 ! spread to grid correction
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
      jgpiic  = 1 ! if gp used in collisions or not
      jgpps   = 2 ! Pointer to proc id for data swap
      jgppt   = 3 ! findpts return processor id
      jgpes   = 4 ! Destination element to be sent to
      jgpicx  = 5
      jgpicy  = 6
      jgpicz  = 7
      jgppid1 = 8 
      jgppid2 = 9
      jgppid3 = 10

c     ghost particle real pointers ----------------------------------
      jgpx    = 1 ! ghost particle xloc
      jgpy    = 2 ! ghost particle yloc
      jgpz    = 3 ! ghost particle zloc
      jgpfh   = 4 ! ghost particle hydrodynamic xforce (i+1 > y, i+2 > z)
      jgpvol  = jgpfh+3  ! ghost particle volume
      jgprpe  = jgpvol+1  ! ghost particle effective radius
      jgpspl  = jgprpe+1 ! spreading correction (if used)
      jgpg0   = jgpspl+1 ! 
      jgpq0   = jgpg0 +1 ! 
      jgpv0   = jgpq0 +1 ! velocity (3 components)

      return
      end
c----------------------------------------------------------------------
c     particle force routines
c----------------------------------------------------------------------
      subroutine lpm_usr_particles_solver
c
c     call routines in ordered way - main solver structure
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage

      if (istep .gt. time_delay) then

      if (stage.eq.1) then

         if (inject_rate .gt. 0) then
            if ( modulo(istep,inject_rate) .eq. 0) call place_particles
         endif

         ! Update where particle is stored at
         ptdum(3) = dnekclock()
            call move_particles_inproc
         pttime(3) = pttime(3) + dnekclock() - ptdum(3)

         if (two_way.gt.1) then
            ! Create ghost/wall particles
            ptdum(4) = dnekclock()
               call create_extra_particles
            pttime(4) = pttime(4) + dnekclock() - ptdum(4)

            ! Send ghost particles
            ptdum(5) = dnekclock()
               call send_ghost_particles
            pttime(5) = pttime(5) + dnekclock() - ptdum(5)
   
            ! Projection to Eulerian grid
            ptdum(6) = dnekclock()
               call spread_props_grid
            pttime(6) = pttime(6) + dnekclock() - ptdum(6)
   
         endif
      endif
            ptdum(4) = dnekclock()
            call create_ghost_particles_col
            pttime(4) = pttime(4) + dnekclock() - ptdum(4)

            ptdum(5) = dnekclock()
            call send_ghost_particles
            pttime(5) = pttime(5) + dnekclock() - ptdum(5)

      ! Interpolate Eulerian properties to particle location
      ptdum(7) = dnekclock()
         call interp_props_part_location
      pttime(7) = pttime(7) + dnekclock() - ptdum(7)

      ! Evaluate particle force models
      ptdum(8) = dnekclock()
         call usr_particles_forcing  
      pttime(8) = pttime(8) + dnekclock() - ptdum(8)

      ! Integrate in time
      ptdum(9) = dnekclock()
         if (abs(time_integ) .eq. 1) call rk3_integrate
         if (abs(time_integ) .eq. 2) call bdf_integrate
      pttime(9) = pttime(9) + dnekclock() - ptdum(9)

      ! Update forces
      ptdum(10) = dnekclock()
         call compute_forcing_post_part
      pttime(10) = pttime(10) + dnekclock() - ptdum(10)

      endif ! time_delay

      return
      end
c----------------------------------------------------------------------
      subroutine spread_props_grid
c
c     spread particle properties at fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'LPM'

      real               phig(lx1,ly1,lz1,lelt)
      common /otherpvar/ phig

      common /myparth/ i_fp_hndl, i_cr_hndl

      real    xx,yy,zz,vol,pfx,pfy,pfz,pmass,pmassf,vcell,multfc,multfci
     >       ,qgqf,rvx,rvy,rvz,rcountv(8,nelt),rx2(3)
     >       ,rxyzp(n_walls*2,3)
      integer e

      common /lpm_fix/ phigdum,phigvdum
      real phigdum(lx1,ly1,lz1,lelt,3),phigvdum(lx1,ly1,lz1,lelt)

      real xrun(lx1),yrun(ly1),zrun(lz1)

      integer xrune(lx1),yrune(lz1),zrune(lz1)


      integer el_save(2,lbmax),nsave

      real rproj(11,llpart_gp+llpart)
      integer iproj(4,llpart_gp+llpart)

      logical partl


      nxyz = lx1*ly1*lz1

      nxyzdum = nxyz*nfldpm*neltb
      call rzero(ptwb,nxyzdum)

      d2chk2_sq = d2chk(2)**2

      do ip=1,n
         rsigp   = rsig*rpart(jrpe,ip)*2.
         multfci = 1./(sqrt(2.*pi)**2 * rsigp**2) ! exponential
         if (if3d) multfci = multfci**(1.5d+0)
         ralph   = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))*
     >           rpart(jrpe,ip)*2.
         ralph2   = ralph**2
         rbexpi   = 1./(-2.*rsigp**2)

         multfc = multfci
         if (.not. if3d) 
     >          multfc = multfc/(xdrange(2,3)-xdrange(1,3))

         rproj(1 ,ip) = rbexpi
                    
         rproj(2 ,ip) = rpart(jx,ip)
         rproj(3 ,ip) = rpart(jy,ip)
         rproj(4 ,ip) = rpart(jz,ip)
         rproj(5 ,ip) = -rpart(jf0,ip)*multfc
         rproj(6 ,ip) = -rpart(jf0+1,ip)*multfc
         rproj(7 ,ip) = -rpart(jf0+2,ip)*multfc
         rproj(8 ,ip) = rpart(jvol,ip)*multfc
c        rproj(9 ,ip) = -(rpart(jg0,ip) + rpart(jq0,ip))*multfc
         rproj(9 ,ip) = rpart(jv0  ,ip)*rproj(8,ip)
         rproj(10,ip) = rpart(jv0+1,ip)*rproj(8,ip)
         rproj(11,ip) = rpart(jv0+2,ip)*rproj(8,ip)

         ii = ipart(jicx,ip)
         jj = ipart(jicy,ip)
         kk = ipart(jicz,ip)
         nn = ipart(jicn,ip)
         iproj(1,ip)  = ii
         iproj(2,ip)  = jj
         iproj(3,ip)  = kk
         iproj(4,ip)  = nn
      enddo

      do ip=1,nfptsgp
         rsigp   = rsig*rptsgp(jgprpe,ip)*2.
         multfci = 1./(sqrt(2.*pi)**2 * rsigp**2) ! exponential
         if (if3d) multfci = multfci**(1.5d+0)
         ralph   = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))*
     >           rptsgp(jgprpe,ip)*2.
         ralph2   = ralph**2
         rbexpi   = 1./(-2.*rsigp**2)

         multfc = multfci
         if (.not. if3d) 
     >          multfc = multfc/(xdrange(2,3)-xdrange(1,3))

         rproj(1 ,ip+n) = rbexpi
         rproj(2 ,ip+n) = rptsgp(jgpx,ip)
         rproj(3 ,ip+n) = rptsgp(jgpy,ip)
         rproj(4 ,ip+n) = rptsgp(jgpz,ip)
         rproj(5 ,ip+n) = -rptsgp(jgpfh,ip)*multfc
         rproj(6 ,ip+n) = -rptsgp(jgpfh+1,ip)*multfc
         rproj(7 ,ip+n) = -rptsgp(jgpfh+2,ip)*multfc
         rproj(8 ,ip+n) = rptsgp(jgpvol,ip)*multfc
c        rproj(9 ,ip+n) = -(rptsgp(jgpg0,ip) + rptsgp(jgpq0,ip))*multfc
         rproj(9 ,ip+n) = rptsgp(jgpv0  ,ip)*rproj(8,ip+n)
         rproj(10,ip+n) = rptsgp(jgpv0+1,ip)*rproj(8,ip+n)
         rproj(11,ip+n) = rptsgp(jgpv0+2,ip)*rproj(8,ip+n)

         ii    = iptsgp(jgpicx,ip)
         jj    = iptsgp(jgpicy,ip)
         kk    = iptsgp(jgpicz,ip)
         iproj(1,ip+n)  = ii
         iproj(2,ip+n)  = jj
         iproj(3,ip+n)  = kk
         iproj(4,ip+n)  = ii + ndxgp*jj + ndxgp*ndygp*kk
      enddo

      ndum = n+nfptsgp

      if (ifrectp .eq. 1) then
         do ip=1,ndum
            ii    = iproj(1,ip)
            jj    = iproj(2,ip)
            kk    = iproj(3,ip)
            ndum  = iproj(4,ip)
          
            ilow  = ii-1
            ihigh = ii+1
            jlow  = jj-1
            jhigh = jj+1
            klow  = kk-1
            khigh = kk+1
         
            d2chk2_e = exp(d2chk_sq*rproj(1,ip))
          
            do ie=1,neltb
         
               if (lpm_el_map(1,ie) .gt. ndum) cycle
               if (lpm_el_map(2,ie) .lt. ndum) cycle 
         
               if (lpm_el_map(3,ie) .gt. ihigh) cycle
               if (lpm_el_map(4,ie) .lt. ilow)  cycle
               if (lpm_el_map(5,ie) .gt. jhigh) cycle
               if (lpm_el_map(6,ie) .lt. jlow)  cycle
               if (lpm_el_map(7,ie) .gt. khigh) cycle
               if (lpm_el_map(8,ie) .lt. klow)  cycle
         
               do i=1,nx1
                  xrun(i) = (xm1b(i,1,1,1,ie) - rproj(2,ip))**2
                  xrune(i) = 1
                  if (xrun(i) .gt. d2chk2_sq) xrune(i) = -1
               enddo
               do j=1,ny1
                  yrun(j) = (xm1b(1,j,1,2,ie) - rproj(3,ip))**2
                  yrune(j) = 1
                  if (yrun(j) .gt. d2chk2_sq) yrune(j) = -1
               enddo
               do k=1,nz1
                  zrun(k) = (xm1b(1,1,k,3,ie) - rproj(4,ip))**2
                  if (.not. if3d) zrun(k) = 0.0
                  zrune(k) = 1
                  if (zrun(k) .gt. d2chk2_sq) zrune(k) = -1
               enddo
         
               do k=1,nz1
                  if (zrune(k) .lt. 0) cycle
               do j=1,ny1
                  if (yrune(j) .lt. 0) cycle
                  rdum1 = yrun(j) + zrun(k)
               do i=1,nx1
               if (mod_gp_grid(i,j,k,ie,4).ne.iproj(4,ip)) cycle
               rdist2 = rdum1 + xrun(i)
                  rexp = exp(rdist2*rproj(1,ip))
                  
                  do jj=1,7
                     j1 = jj+4
                     ptwb(i,j,k,jj,ie) = ptwb(i,j,k,jj,ie) 
     >                                    + rproj(j1,ip)*rexp
                  enddo
               enddo
               enddo
               enddo
            enddo
         enddo
      elseif (ifrectp .eq. 0) then
         do ip=1,ndum
            ii    = iproj(1,ip)
            jj    = iproj(2,ip)
            kk    = iproj(3,ip)
            ndum  = iproj(4,ip)
          
            ilow  = ii-1
            ihigh = ii+1
            jlow  = jj-1
            jhigh = jj+1
            klow  = kk-1
            khigh = kk+1
         
            d2chk2_e = exp(d2chk_sq*rproj(1,ip))
          
            do ie=1,neltb
         
               if (lpm_el_map(1,ie) .gt. ndum) cycle 
               if (lpm_el_map(2,ie) .lt. ndum) cycle 
         
               if (lpm_el_map(3,ie) .gt. ihigh) cycle
               if (lpm_el_map(4,ie) .lt. ilow)  cycle
               if (lpm_el_map(5,ie) .gt. jhigh) cycle
               if (lpm_el_map(6,ie) .lt. jlow)  cycle
               if (lpm_el_map(7,ie) .gt. khigh) cycle
               if (lpm_el_map(8,ie) .lt. klow)  cycle
         
               do i=1,nxyz
                  if (mod_gp_grid(i,1,1,ie,4).ne.iproj(4,ip)) cycle
                  rdist2  = (xm1b(i,1,1,1,ie) - rproj(2,ip))**2 +
     >                      (xm1b(i,1,1,2,ie) - rproj(3,ip))**2
                  if(if3d) rdist2 = rdist2                      +
     >                      (xm1b(i,1,1,3,ie) - rproj(4,ip))**2
                  if (rdist2 .gt. d2chk2_sq) cycle
               
                  rexp = exp(rdist2*rproj(1,ip))
                  
                  do j=1,7
                     j1 = j+4
                     ptwb(i,1,1,j,ie) = ptwb(i,1,1,j,ie) 
     >                                    + rproj(j1,ip)*rexp
                  enddo
                enddo
             enddo
          enddo
      endif

      ! now send xm1b to the processors in nek that hold xm1

      neltbc = neltb
      ndum = 6*neltbc
      call icopy(er_mapc,er_map,ndum)
      do ie=1,neltbc
         er_mapc(5,ie) = er_mapc(2,ie)
         er_mapc(6,ie) = er_mapc(2,ie)

      enddo
      nl = 0
      nii = 6
      njj = 6
      nrr = nxyz*nfldpm
      call fgslib_crystal_tuple_transfer(i_cr_hndl,neltbc,lbmax
     >                  , er_mapc,nii,partl,nl,ptwb,nrr,njj)

      ! add the fields from the bins to ptw array
      nlxyze = lx1*ly1*lz1*lelt
      call rzero(ptw,nlxyze*7)
      do ie=1,neltbc
         iee = er_mapc(1,ie)
         do ip=1,7
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            ptw(i,j,k,iee,ip) = ptw(i,j,k,iee,ip) + ptwb(i,j,k,ip,ie)
         enddo
         enddo
         enddo
         enddo
      enddo

#ifdef CMTNEK
      ! filtering makes velocity field smoother for p*grad(phig) in CMT

c     call dsavg(ptw(1,1,1,1,1))
c     call dsavg(ptw(1,1,1,1,2))
c     call dsavg(ptw(1,1,1,1,3))
c     call dsavg(ptw(1,1,1,1,4))
c     call dsavg(ptw(1,1,1,1,5))
c     call dsavg(ptw(1,1,1,1,6))
c     call dsavg(ptw(1,1,1,1,7))

      wght = 1
      ncut = 1
      call filter_s0(ptw(1,1,1,1,1),wght,ncut,'phip') 
      call filter_s0(ptw(1,1,1,1,2),wght,ncut,'phip') 
      call filter_s0(ptw(1,1,1,1,3),wght,ncut,'phip') 
      call filter_s0(ptw(1,1,1,1,4),wght,ncut,'phip') 
      call filter_s0(ptw(1,1,1,1,5),wght,ncut,'phip') 
      call filter_s0(ptw(1,1,1,1,6),wght,ncut,'phip') 
      call filter_s0(ptw(1,1,1,1,7),wght,ncut,'phip') 

c     call dsavg(ptw(1,1,1,1,1))
c     call dsavg(ptw(1,1,1,1,2))
c     call dsavg(ptw(1,1,1,1,3))
c     call dsavg(ptw(1,1,1,1,4))
c     call dsavg(ptw(1,1,1,1,5))
c     call dsavg(ptw(1,1,1,1,6))
c     call dsavg(ptw(1,1,1,1,7))
#endif

      rvfmax = 0.74
      rvfmin = 1E-15
      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         if (ptw(i,j,k,ie,4) .gt. rvfmax) ptw(i,j,k,ie,4) = rvfmax
         if (ptw(i,j,k,ie,4) .lt. rvfmin) ptw(i,j,k,ie,4) = rvfmin
         phig(i,j,k,ie) = 1. - ptw(i,j,k,ie,4)
      enddo
      enddo
      enddo
      enddo

c     wght = 1
c     ncut = 1
c     call filter_s0(phig,wght,ncut,'phip') 
c     call filter_s0(ptw(1,1,1,1,5),wght,ncut,'phip') 
c     call filter_s0(ptw(1,1,1,1,6),wght,ncut,'phip') 
c     call filter_s0(ptw(1,1,1,1,7),wght,ncut,'phip') 

      nxyze = nxyz*nelt

      call gradm1(phigdum(1,1,1,1,1),
     >            phigdum(1,1,1,1,2),
     >            phigdum(1,1,1,1,3),
     >            phig)


      call col2(phigdum(1,1,1,1,1),pr,nxyze)
      call col2(phigdum(1,1,1,1,2),pr,nxyze)
      if(if3d) call col2(phigdum(1,1,1,1,3),pr,nxyze)

c     call dsavg(phigdum(1,1,1,1,1))
c     call dsavg(phigdum(1,1,1,1,2))
c     call dsavg(phigdum(1,1,1,1,3))

c     nspread = 6
c     rle = 5.94E-4
c     rgrad_max = 101325/(nspread*rle)

c     do i=1,nxyze
c        if (abs(phigdum(i,1,1,1,2)) .lt. rgrad_max) then
c           if (phigdum(i,1,1,1,2) .lt. 0)
c    >         phigdum(i,1,1,1,2) = - rgrad_max
c           if (phigdum(i,1,1,1,2) .gt. 0)
c    >         phigdum(i,1,1,1,2) = rgrad_max
c        endif
c     enddo


      call opdiv(phigvdum,ptw(1,1,1,1,5),
     >                    ptw(1,1,1,1,6),
     >                    ptw(1,1,1,1,7))


      call col2(phigvdum,pr,nxyze)
      call chsign(phigvdum,nxyze)

c     call dsavg(phigvdum)

      return
      end
c----------------------------------------------------------------------
      subroutine rk3_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage
      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      common /PARTRK3/ kv_stage_p
      real   kv_stage_p(llpart,7)

      integer fdim
      real    pmass

      jx0 = jx

c     rk3 stage one items ---------------------------------------------
      if (stage.eq.1) then
c        used for time derivative of v in iu force
         call get_bdf_ext_coefs(beta,alpha,times)

c        move data to previous positions
         do j=0,ndim-1
         do i=1,n
            rpart(ju3+j,i)=rpart(ju2+j,i)
            rpart(ju2+j,i)=rpart(ju1+j,i)
            rpart(ju1+j,i)=rpart(ju0+j,i)
            rpart(jv3+j,i)=rpart(jv2+j,i)
            rpart(jv2+j,i)=rpart(jv1+j,i)
            rpart(jv1+j,i)=rpart(jv0+j,i)
            rpart(jx3+j,i)=rpart(jx2+j,i)
            rpart(jx2+j,i)=rpart(jx1+j,i)
            rpart(jx1+j,i)=rpart(jx0+j,i)
         enddo
         enddo

         do i=1,n
            kv_stage_p(i,1) = rpart(jx0  ,i)
            kv_stage_p(i,2) = rpart(jx0+1,i)
            kv_stage_p(i,3) = rpart(jx0+2,i)
            kv_stage_p(i,4) = rpart(jv0  ,i)
            kv_stage_p(i,5) = rpart(jv0+1,i)
            kv_stage_p(i,6) = rpart(jv0+2,i)
            kv_stage_p(i,7) = rpart(jtemp,i)
         enddo
      endif

c     all rk3 stages items --------------------------------------------
      do i=1,n
         rpart(jv0  ,i) = tcoef(1,stage)*kv_stage_p(i,4) +
     >                    tcoef(2,stage)*rpart(jv0  ,i)  +
     >                    tcoef(3,stage)*rpart(jf0  ,i)
         rpart(jv0+1,i) = tcoef(1,stage)*kv_stage_p(i,5) +
     >                    tcoef(2,stage)*rpart(jv0+1,i)  +
     >                    tcoef(3,stage)*rpart(jf0+1,i)
         if (if3d)
     >   rpart(jv0+2,i) = tcoef(1,stage)*kv_stage_p(i,6) +
     >                    tcoef(2,stage)*rpart(jv0+2,i)  +
     >                    tcoef(3,stage)*rpart(jf0+2,i)

         rpart(jx0  ,i) = tcoef(1,stage)*kv_stage_p(i,1) +
     >                    tcoef(2,stage)*rpart(jx0  ,i)  +
     >                    tcoef(3,stage)*rpart(jv0  ,i)
         rpart(jx0+1,i) = tcoef(1,stage)*kv_stage_p(i,2) +
     >                    tcoef(2,stage)*rpart(jx0+1,i)  +
     >                    tcoef(3,stage)*rpart(jv0+1,i)
         if (if3d)
     >   rpart(jx0+2,i) = tcoef(1,stage)*kv_stage_p(i,3) +
     >                    tcoef(2,stage)*rpart(jx0+2,i)  +
     >                    tcoef(3,stage)*rpart(jv0+2,i)

         rpart(jtemp,i) = tcoef(1,stage)*kv_stage_p(i,7) +
     >                    tcoef(2,stage)*rpart(jtemp,i)  +
     >                    tcoef(3,stage)*rpart(jq0  ,i)
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_forcing
c
c     calculate the rhs of particle equation
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real vel_diff,pmass,pmassf

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,7)

      if ((part_force(2) .ne. 0) .or. (part_force(3) .ne. 0))
     >   call calc_substantial_derivative

      do i=1,n
         rpart(jfcol  ,i) = 0.
         rpart(jfcol+1,i) = 0.
         rpart(jfcol+2,i) = 0.
      enddo

      do i=1,n
c        setup values ------------------------------------------------
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if(part_force(3).ne.0) pmass = pmass + rpart(jcmiu,i)*pmassf ! am

         call get_part_use_block(i,0)

c        momentum rhs ------------------------------------------------
         lpmforce(1) = 0
         lpmforce(2) = 0
         lpmforce(3) = 0
         lpmforcec(1) = 0 ! coupled to fluid
         lpmforcec(2) = 0 ! coupled to fluid
         lpmforcec(3) = 0 ! coupled to fluid
         call lpm_usr_f  ! user/body force
            rpart(jfusr+0,i) = lpmforce(1) + lpmforcec(1)
            rpart(jfusr+1,i) = lpmforce(2) + lpmforcec(2)
            rpart(jfusr+2,i) = lpmforce(3) + lpmforcec(3)

         if (time_integ .gt. 0) then
            call lpm_f_qs
               rpart(jfqs+0,i) = lpmforce(1)
               rpart(jfqs+1,i) = lpmforce(2)
               rpart(jfqs+2,i) = lpmforce(3)
            call lpm_f_un
               rpart(jfun+0,i) = lpmforce(1)
               rpart(jfun+1,i) = lpmforce(2)
               rpart(jfun+2,i) = lpmforce(3)
            call lpm_f_iu
               rpart(jfiu+0,i) = lpmforce(1)
               rpart(jfiu+1,i) = lpmforce(2)
               rpart(jfiu+2,i) = lpmforce(3)
         endif

         call usr_particles_f_col(i)     ! colision force 

         rfx = rpart(jfusr+0,i) + 
     >         rpart(jfqs +0,i) +
     >         rpart(jfun +0,i) +
     >         rpart(jfiu +0,i) +
     >         rpart(jfcol+0,i)
         rfy = rpart(jfusr+1,i) + 
     >         rpart(jfqs +1,i) +
     >         rpart(jfun +1,i) +
     >         rpart(jfiu +1,i) +
     >         rpart(jfcol+1,i)
         if (if3d)
     >   rfz = rpart(jfusr+2,i) + 
     >         rpart(jfqs +2,i) +
     >         rpart(jfun +2,i) +
     >         rpart(jfiu +2,i) +
     >         rpart(jfcol+2,i)

         ! mass weighted total force
         rpart(jf0+0,i) = rfx/pmass
         rpart(jf0+1,i) = rfy/pmass
         rpart(jf0+2,i) = rfz/pmass

c        energy rhs --------------------------------------------------
         lpmheat = 0
         if (time_integ .gt. 0) then
            call lpm_q_uu
               rpart(jquu,i) = lpmheat
            call lpm_q_qs
               rpart(jqqs,i) = lpmheat

            rqq = rpart(jquu,i) +
     >            rpart(jqqs,i)

            pmass = rpart(jvol,i)*rpart(jrhop,i)
            rpart(jq0,i) = rqq/(pmass*cp_p)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_part_use_block(i,ipre)
c
c     set common block for particle-user use
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
#ifdef CMTNEK
      include 'CMTDATA'
#endif

      real MixtPerf_C_GRT_part
      external MixtPerf_C_GRT_part

      if (ipre .eq. 0) then
      lpmx_p(1)    = rpart(jx,i)
      lpmx_p(2)    = rpart(jy,i)
      lpmx_p(3)    = rpart(jz,i)
      lpmvol_p     = rpart(jvol,i)
      lpmvolfrac_p = rpart(jvol1,i)
      lpmvolfrac_f = 1.0 - lpmvolfrac_p
      lpmtau_p     = rpart(jtaup,i)
      lpmv_p(1)    = rpart(jv0+0,i)
      lpmv_p(2)    = rpart(jv0+1,i)
      lpmv_p(3)    = rpart(jv0+2,i)
      lpmv_f(1)    = rpart(ju0+0,i)
      lpmv_f(2)    = rpart(ju0+1,i)
      lpmv_f(3)    = 0.0
      if(if3d) lpmv_f(3)    = rpart(ju0+2,i)
      lpmvisc_f    = mu_0
      lpmdiam_p    = rpart(jdp,i)
      lpmdens_p    = rpart(jrhop,i)
      lpmdens_f    = rpart(jrho,i)
      lpmtemp_p    = rpart(jtemp,i)
      lpmtemp_f    = rpart(jtempf,i)
      lpmDuDt(1)   = rpart(jDuDt+0,i)
      lpmDuDt(2)   = rpart(jDuDt+1,i)
      lpmDuDt(3)   = 0.0
      if(if3d) lpmDuDt(3)   = rpart(jDuDt+2,i)
      lpmkappa_f   = abs(param(8))
      lpmvdiff_pf  = sqrt((lpmv_p(1) - lpmv_f(1))**2+
     >                    (lpmv_p(2) - lpmv_f(2))**2+
     >                    (lpmv_p(3) - lpmv_f(3))**2)
      if (abs(lpmvdiff_pf) .lt. 1E-6) lpmvdiff_pf=1E-6

      lpmmach_p    = MixtPerf_C_GRT_part(gmaref,rgasref,
     >                                   rpart(jtempf,i),icmtp)
      if (icmtp .eq. 0) then
         lpmmach_p  = 0.0
      elseif (icmtp .eq. 1) then
         lpmmach_p = lpmvdiff_pf/lpmmach_p
      endif
      lpmre_p      = rpart(jrho,i)*rpart(jdp,i)*lpmvdiff_pf/mu_0 ! Re
      lpmi         = i
      elseif (ipre .eq. 1) then

      rpart(jx,i)        =  lpmx_p(1)    
      rpart(jy,i)        =  lpmx_p(2)    
      rpart(jz,i)        =  lpmx_p(3)    
      rpart(jvol,i)      =  lpmvol_p     
      rpart(jvol1,i)     =  lpmvolfrac_p 
      lpmvolfrac_f       =  1.0 - lpmvolfrac_p
      rpart(jtaup,i)     =  lpmtau_p     
      rpart(jv0+0,i)     =  lpmv_p(1)    
      rpart(jv0+1,i)     =  lpmv_p(2)    
      rpart(jv0+2,i)     =  lpmv_p(3)    
      rpart(ju0+0,i)     =  lpmv_f(1)    
      rpart(ju0+1,i)     =  lpmv_f(2)    
      rpart(ju0+2,i)     =  lpmv_f(3)    
      rpart(jdp,i)       =  lpmdiam_p    
      rpart(jrhop,i)     =  lpmdens_p    
      rpart(jrho,i)      =  lpmdens_f    
      rpart(jtemp,i)     =  lpmtemp_p    
      rpart(jtempf,i)    =  lpmtemp_f    
      rpart(jDuDt+0,i)   =  lpmDuDt(1)   
      rpart(jDuDt+1,i)   =  lpmDuDt(2)   
      rpart(jDuDt+2,i)   =  lpmDuDt(3)   

      lpmvdiff_pf  = sqrt((rpart(ju0  ,i)-rpart(jv0  ,i))**2+
     >                    (rpart(ju0+1,i)-rpart(jv0+1,i))**2+
     >                    (rpart(ju0+2,i)-rpart(jv0+2,i))**2)
      if (abs(lpmvdiff_pf) .lt. 1E-6) lpmvdiff_pf=1E-6

      lpmmach_p    = MixtPerf_C_GRT_part(gmaref,rgasref,
     >                                   rpart(jtempf,i),icmtp)
      lpmmach_p    = lpmvdiff_pf/lpmmach_p
      lpmre_p      = rpart(jrho,i)*rpart(jdp,i)*lpmvdiff_pf/mu_0 ! Re

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_forcing_post_part
c
c     post calculate forces due to factoring of equations
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real uvel(0:2), vvel(0:2), pmass, pmassf,S_qs

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,7)

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      pi  = 4.0d+0*atan(1.0d+0)

      do i=1,n
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if (part_force(3).ne.0) pmass =pmass + rpart(jcmiu,i)*pmassf

         call get_part_use_block(i,0)

         ! get user coupled to fluid forcing
         lpmforcec(1) = 0 ! coupled to fluid
         lpmforcec(2) = 0 ! coupled to fluid
         lpmforcec(3) = 0 ! coupled to fluid
         call lpm_usr_f   ! user/body force

c        momentum forcing to fluid
         do j=0,ndim-1
            rdum = 0.

            rdvdt = rpart(jf0+j,i) ! note already divided by Mp + am
            ram_s = rdvdt*rpart(jcmiu,i)*pmassf
            rpart(jfiu+j,i) = rpart(jfiu+j,i) - ram_s

c           note that no coupled f_un in this formulation
            rdum = rdum + lpmforcec(j+1)
            rdum = rdum + rpart(jfqs+j,i)
            rdum = rdum + rpart(jfiu+j,i)

            rpart(jf0+j,i) = rdum ! now the force to couple with gas
         enddo

c        energy forcing to fluid (quasi-steady)
         rpart(jg0,i) = rpart(jv0  ,i)*rpart(jfqs  ,i) + !force work
     >                  rpart(jv0+1,i)*rpart(jfqs+1,i) +
     >                  rpart(jv0+2,i)*rpart(jfqs+2,i)
         rpart(jg0,i) = rpart(jg0,i) + 
     >                  rpart(ju0  ,i)*rpart(jfiu  ,i) + !iu
     >                  rpart(ju0+1,i)*rpart(jfiu+1,i) +
     >                  rpart(ju0+2,i)*rpart(jfiu+2,i)

         rdum = 0.
         rdum = rdum + rpart(jqqs,i)

         rpart(jq0,i) = rdum

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_substantial_derivative
c 
c     calculate rhs of NS, which is just pressure gradient. used for
c     undisturbed and invisicid unsteady force
c
c     no forcing included...should it be?
c
c     rhs_fluidp(i,j,k,e,l)
c   
c        l = 1     >    dP/dx
c        l = 2     >    dP/dy
c        l = 3     >    dP/dz
c        l = 4     >    -P* div(phi_p v), v is eulerian particle vel
c        l = 5     >    P * d phi_g/dx
c        l = 6     >    P * d phi_g/dy
c        l = 7     >    P * d phi_g/dz
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'LPM'

      real               phig(lx1,ly1,lz1,lelt)
      common /otherpvar/ phig

      integer e
      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1),
     >         pm1(lx1,ly1,lz1,lelt,3),udum(lx1,ly1,lz1)

      nxyz=nx1*ny1*nz1
      nlxyze = lx1*ly1*lz1*lelt

      call rzero(rhs_fluidp,nlxyze*7)

      ! if pressure solved on different mesh, map to vel mesh
      if (lx2 .ne. lx1) then
         call mappr(pm1(1,1,1,1,1),pr,pm1(1,1,1,1,2),pm1(1,1,1,1,3))
      else
         call copy(pm1(1,1,1,1,1),pr(1,1,1,1),nlxyze)
      endif

c     compute grad pr
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        pm1(1,1,1,e,1),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,1) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,2) = 1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,3) = 1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
      enddo

      ! div (phi_p * v)
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! x dir
     >                                        ptw(1,1,1,e,5),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! y dir
     >                                        ptw(1,1,1,e,6),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! z dir
     >                                        ptw(1,1,1,e,7),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,4) = -pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,4)
      enddo
      enddo
      enddo
      enddo

c     compute grad phi_g
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        phig(1,1,1,e),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,5) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,6) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,7) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,5) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,5)
         rhs_fluidp(i,j,k,e,6) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,6)
         rhs_fluidp(i,j,k,e,7) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,7)
      enddo
      enddo
      enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_col(i)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real mcfac,mcfac_wall, rdum, rdum3(3)

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      real rpx1(3), rpx2(3), rpx0(3),r1(3),r2(3),r3(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      integer col_list(llpart,3),ncol

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)


      ptdum(11) = dnekclock()

      mcfac       = 2.*sqrt(ksp)*log(e_rest)/
     >              sqrt(log(e_rest)**2+pi**2)
      mcfac_wall  = 2.*sqrt(ksp_wall)*log(e_rest_wall)/
     >              sqrt(log(e_rest_wall)**2+pi**2)

      if (two_way .gt. 2) then

      rrp1   = rpart(jrpe ,i)
      rvol1  = rpart(jvol ,i)
      rrho1  = rpart(jrhop,i)
      rx1(1) = rpart(jx   ,i)
      rx1(2) = rpart(jx+1 ,i)
      rx1(3) = rpart(jx+2 ,i)
      rv1(1) = rpart(jv0  ,i)
      rv1(2) = rpart(jv0+1,i)
      rv1(3) = rpart(jv0+2,i)
c     rx1(1) = alpha(1)*rpart(jx1  ,i) +alpha(2)*rpart(jx2  ,i)+
c    >         alpha(3)*rpart(jx3  ,i)
c     rx1(2) = alpha(1)*rpart(jx1+1,i)+alpha(2)*rpart(jx2+1,i)+
c    >         alpha(3)*rpart(jx3+1,i)
c     rx1(3) = alpha(1)*rpart(jx1+2,i)+alpha(2)*rpart(jx2+2,i)+
c    >         alpha(3)*rpart(jx3+2,i)
c     rv1(1) = alpha(1)*rpart(jv1  ,i) +alpha(2)*rpart(jv2  ,i)+
c    >         alpha(3)*rpart(jv3  ,i)
c     rv1(2) = alpha(1)*rpart(jv1+1,i)+alpha(2)*rpart(jv2+1,i)+
c    >         alpha(3)*rpart(jv3+1,i)
c     rv1(3) = alpha(1)*rpart(jv1+2,i)+alpha(2)*rpart(jv2+2,i)+
c    >         alpha(3)*rpart(jv3+2,i)

      icx1 = ipart(jicx,i)
      icy1 = ipart(jicy,i)
      
      icz1 = 0
      if (if3d) icz1 = ipart(jicz,i)

      icxm = icx1 -1
      icxp = icx1 +1
      icym = icy1 -1
      icyp = icy1 +1

      iczm = 0
      iczp = 0
      if (if3d) then
         iczm = icz1 -1
         iczp = icz1 +1
      endif

      ndum1 = ipart(jicn,i)

      ncol = 0

c     let every particle search for itself
c        particles in local elements
         do j = i+1,n
c           if( i .eq. j) cycle

            icx2  = ipart(jicx,j)
            icy2  = ipart(jicy,j)
            icz2  = ipart(jicz,j)
            ndum2 = ipart(jicn,j)
c           if (ndum1 .ne. ndum2) cycle
            if ((icx2.ge.icxm).and.(icx2.le.icxp)) then
            if ((icy2.ge.icym).and.(icy2.le.icyp)) then
            if ((icz2.ge.iczm).and.(icz2.le.iczp)) then

c           if (abs(rpart(jx,i) -rpart(jx,j)) .gt. d2chk(3)) cycle
c           if (abs(rpart(jy,i) -rpart(jy,j)) .gt. d2chk(3)) cycle
c           if (abs(rpart(jz,i) -rpart(jz,j)) .gt. d2chk(3)) cycle

            rrp2   = rpart(jrpe ,j)
            rvol2  = rpart(jvol ,j)
            rrho2  = rpart(jrhop,j)
            rx2(1) = rpart(jx   ,j)
            rx2(2) = rpart(jx+1 ,j)
            rx2(3) = rpart(jx+2 ,j)
            rv2(1) = rpart(jv0  ,j)
            rv2(2) = rpart(jv0+1,j)
            rv2(3) = rpart(jv0+2,j)

            idum = 1
            call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rpart(jfcol,j),idum)

            endif
            endif
            endif

 1234 continue
         enddo




c        search list of ghost particles
         do j = 1,nfptsgp
c        do j = 1,0
            icx2 = iptsgp(jgpicx,j)
            icy2 = iptsgp(jgpicy,j)
            icz2 = iptsgp(jgpicz,j)
            if ((icx2.ge.icxm).and.(icx2.le.icxp)) then
            if ((icy2.ge.icym).and.(icy2.le.icyp)) then
            if ((icz2.ge.iczm).and.(icz2.le.iczp)) then

            rrp2   = rptsgp(jgprpe ,j)
            rvol2  = rptsgp(jgpvol ,j)
            rrho2  = rho_p            ! assume same density. Need2fix
            rx2(1) = rptsgp(jgpx   ,j)
            rx2(2) = rptsgp(jgpx+1 ,j)
            rx2(3) = rptsgp(jgpx+2 ,j)
            rv2(1) = rptsgp(jgpv0  ,j)
            rv2(2) = rptsgp(jgpv0+1,j)
            rv2(3) = rptsgp(jgpv0+2,j)

c           if (ncol .gt. 0) then           
c              do k=1,ncol
c                 if (iptsgp(jgppid1,j) .eq. col_list(k,1)) then
c                 if (iptsgp(jgppid2,j) .eq. col_list(k,2)) then
c                 if (iptsgp(jgppid3,j) .eq. col_list(k,3)) then
c                    goto 1235
c                 endif
c                 endif
c                 endif
c              enddo
c           endif

            idum = 0
            call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rdum3,idum)

c           ncol = ncol +1
c           col_list(ncol,1) = iptsgp(jgppid1,j)
c           col_list(ncol,2) = iptsgp(jgppid2,j)
c           col_list(ncol,3) = iptsgp(jgppid3,j)

            endif
            endif
            endif
 1235 continue
         enddo

         ! plane wall collisions
         do j = 1,np_walls
            rnx = plane_wall_coords(1,j)
            rny = plane_wall_coords(2,j)
            rnz = plane_wall_coords(3,j)
            rpx = plane_wall_coords(4,j)
            rpy = plane_wall_coords(5,j)
            rpz = rx1(3)
            if (if3d) rpz = plane_wall_coords(6,j)


            rd    = -(rnx*rpx + rny*rpy + rnz*rpz)

            rdist = abs(rnx*rpart(jx,i)+rny*rpart(jy,i)+rnz*rpart(jz,i)
     >                    +rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            rrp2   = 0.
            rvol2  = 1.
            rrho2  = 1E8
            rx2(1) = rpart(jx  ,i) - rdist*rnx
            rx2(2) = rpart(jx+1,i) - rdist*rny
            rx2(3) = rpart(jx+2,i) - rdist*rnz
            rv2(1) = 0.
            rv2(2) = 0.
            rv2(3) = 0.

            idum = 0
            call compute_collide(mcfac_wall,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rdum3,idum)
         enddo

         ! cylinder wall collisions
         do j = 1,nc_walls
            rnx = cyl_wall_coords(1,j)
            rny = cyl_wall_coords(2,j)
            rnz = cyl_wall_coords(3,j)
            rpx = cyl_wall_coords(4,j)
            rpy = cyl_wall_coords(5,j)
            rpz = rx1(3)
            if (if3d) rpz = cyl_wall_coords(6,j)


            rrad = cyl_wall_coords(7,j)

            rx2(1) = rpart(jx,i)
            rx2(2) = rpart(jy,i)
            rx2(3) = rpart(jz,i)
            ! for now only works with cylinders aligned with axes at
            ! origin
            if (rnz .gt. 0.5) then
               rtheta = atan2(rpart(jy,i)-rpy,rpart(jx,i)-rpx)
               rx2(1) = rpx+rrad*cos(rtheta)
               rx2(2) = rpy+rrad*sin(rtheta)
            elseif (rnx .gt. 0.5) then
               rtheta = atan2(rpart(jz,i)-rpz,rpart(jy,i)-rpy)
               rx2(2) = rpy+rrad*cos(rtheta)
               rx2(3) = rpz+rrad*sin(rtheta)
            elseif (rny .gt. 0.5) then
               rtheta = atan2(rpart(jx,i)-rpx,rpart(jz,i)-rpz)
               rx2(3) = rpz+rrad*cos(rtheta)
               rx2(1) = rpx+rrad*sin(rtheta)
            endif

            rrp2   = 0.
            rvol2  = 1.
            rrho2  = 1E8
            rv2(1) = 0.
            rv2(2) = 0.
            rv2(3) = 0.

            idum = 0
            call compute_collide(mcfac_wall,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rdum3,idum)
         enddo

      endif


      pttime(11) = pttime(11) + dnekclock() - ptdum(11)

      return
      end
c----------------------------------------------------------------------
      subroutine compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,fcf1,
     >                                                      fcf2,iflg)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real fcf1(3),fcf2(3),er,eta,ere,mcfac
      integer iflg

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      rthresh  = rrp1 + rrp2
      rthresh2 = rthresh**2

      rxdiff  = rx2(1) - rx1(1)
      rsum2 = rxdiff**2
      if (rsum2 .gt. rthresh2) goto 1511

      rydiff = rx2(2) - rx1(2)
      rydiff2 = rydiff**2
      rsum2 = rsum2 + rydiff2
      if (rsum2 .gt. rthresh2) goto 1511

      if (if3d) then
         rzdiff = rx2(3) - rx1(3)
         rzdiff2 = rzdiff**2
         rsum2 = rsum2 + rzdiff2
         if (rsum2 .gt. rthresh2) goto 1511
      endif

      rdiff = sqrt(rsum2)
      rm1   = rrho1*rvol1
      rm2   = rrho2*rvol2

      ! wall
      if (rm2 .gt. 1E7) then
         rksp_use = ksp_wall
         re_rest_use = e_rest_wall
      ! other particles
      else
         rksp_use = ksp
         re_rest_use = e_rest
      endif

      rmult = 1./sqrt(1./rm1 + 1./rm2)
      eta  = mcfac*rmult

      ! first, handle normal collision part
      rbot     = 1./rdiff
      rn_12x   = rxdiff*rbot
      rn_12y   = rydiff*rbot
      rn_12z   = rzdiff*rbot

      rdelta12 = rthresh - rdiff

      rv12_mag = (rv2(1) - rv1(1))*rn_12x +
     >           (rv2(2) - rv1(2))*rn_12y +
     >           (rv2(3) - rv1(3))*rn_12z

      rv12_mage = rv12_mag*eta

      rksp_max = rksp_use*rdelta12

      rnmag = -rksp_max - rv12_mage

      rfn1 = rnmag*rn_12x
      rfn2 = rnmag*rn_12y
      rfn3 = 0.0
      if (if3d) rfn3 = rnmag*rn_12z

      fcf1(1) = fcf1(1) + rfn1
      fcf1(2) = fcf1(2) + rfn2
      fcf1(3) = fcf1(3) + rfn3

      if (iflg .eq. 1) then
          fcf2(1) = fcf2(1) - rfn1
          fcf2(2) = fcf2(2) - rfn2
          fcf2(3) = fcf2(3) - rfn3
      endif


 1511 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine send_ghost_particles
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      logical partl         

c     ndum = iglsum(nfptsgp,1)
c     if (nid .eq. 0) write(6,*) 'Before nfptsgp',ndum

      ! send ghost particles
      call fgslib_crystal_tuple_transfer(i_cr_hndl,nfptsgp,llpart_gp
     $           , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpps) ! jgpps is overwri

c     ndum = iglsum(nfptsgp,1)
c     if (nid .eq. 0) write(6,*) 'After nfptsgp',ndum

      ! sort by element for ghost particles for projection performance
c     call fgslib_crystal_tuple_sort    (i_cr_hndl,nfptsgp
c    $              , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpes,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_extra_particles
c
c     create ghost and wall particles
c
c     bc_part = -1,1  => non-periodic search
c     bc_part = 0  => periodic search
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      character*132 deathmessage

      call create_ghost_particles_full 

      nmax   = iglmax(n,1)
      ngpmax = iglmax(nfptsgp,1)

      if (nmax .gt. llpart) then
         if (nid.eq.0) write(6,1) nmax, llpart, nid
         call exitt
      elseif (ngpmax .gt. llpart_gp) then
         if (nid.eq.0) write(6,2) ngpmax, llpart_gp, nid
         call exitt
      endif
    1 format('Max number of real particles:',
     >   i9,'. Not moving because llpart =',i9, ' on nid = ', i9)
    2 format('Max number of ghost particles:',
     >   i9,'. Not moving because llpart_gp =',i9, ' on nid = ', i9)

      return
      end
c----------------------------------------------------------------------
      subroutine create_ghost_particles_full
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in rptsgp(jgpx,j),rptsgp(jgpy,j),
c     rptsgp(jgpz,j), while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),ntypesl(7),iadd2(3),gpsave(27)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      xdlen = binb(2) - binb(1)
      ydlen = binb(4) - binb(3)
      zdlen = -1.
      if (if3d) zdlen = binb(6) - binb(5)
      if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) xdlen = -1
      if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) ydlen = -1
      if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      nfptsgp = 0

      if (npro_method .eq. 0) return ! only for projection

      do ip=1,n
         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip  = ipart(jicx,ip)
         jjp  = ipart(jicy,ip)
         kkp  = ipart(jicz,ip)

         rxl = binb(1) + rdxgp*iip
         rxr = rxl + rdxgp
         ryl = binb(3) + rdygp*jjp
         ryr = ryl + rdygp
         rzl = 0.0
         rzr = 0.0
         if (if3d) then
            rzl = binb(5) + rdzgp*kkp
            rzr = rzl + rdzgp
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 2 .or. iflgsum .eq. 3) cycle

            ndum  = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndum/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)

            iadd2(1) = iig
            iadd2(2) = jjg
            iadd2(3) = kkg
            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 1 .or. iflgsum .eq. 3) cycle

            ndum  = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndum/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)
            
            iadd2(1) = iig
            iadd2(2) = jjg
            iadd2(3) = kkg
            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 1 .or. iflgsum .eq. 2) cycle

            ndum  = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndum/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)
            
            iadd2(1) = iig
            iadd2(2) = jjg
            iadd2(3) = kkg
            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  333 continue
         enddo
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine create_ghost_particles_full_new
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in rptsgp(jgpx,j),rptsgp(jgpy,j),
c     rptsgp(jgpz,j), while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),ntypesl(7),iadd2(3),gpsave(27)


! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      xdlen = binb(2) - binb(1)
      ydlen = binb(4) - binb(3)
      zdlen = -1.
      if (if3d) zdlen = binb(6) - binb(5)
      if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) xdlen = -1
      if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) ydlen = -1
      if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      nfptsgp = 0

      do ip=1,n
         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip = ipart(jicx,ip) 
         jjp = ipart(jicy,ip) 
         kkp = ipart(jicz,ip) 

         rxl = binb(1) + rdxgp*iip
         rxr = rxl + rdxgp
         ryl = binb(3) + rdygp*jjp
         ryr = ryl + rdygp
         rzl = 0.0
         rzr = 0.0
         if (if3d) then
            rzl = binb(5) + rdzgp*kkp
            rzr = rzl + rdzgp
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle
               

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 2 .or. iflgsum .eq. 3) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndumn/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank


            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)

            idum = iip
            jdum = jjp
            kdum = kkp

            if (iflgx .eq. 1) then
            if (iip .eq. 0)       idum = ndxgp
            if (iip .eq. ndxgp-1) idum = -1
            endif
            if (iflgy .eq. 1) then
            if (jjp .eq. 0)       jdum = ndygp
            if (jjp .eq. ndygp-1) jdum = -1
            endif
            if (iflgz .eq. 1) then
            if (kkp .eq. 0)       kdum = ndzgp
            if (kkp .eq. ndzgp-1) kdum = -1
            endif

            iadd2(1) = idum
            iadd2(2) = jdum
            iadd2(3) = kdum
            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 1 .or. iflgsum .eq. 3) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndumn/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)
            
            idum = iip
            jdum = jjp
            kdum = kkp

            if (iflgx .eq. 1) then
            if (iip .eq. 0)       idum = ndxgp
            if (iip .eq. ndxgp-1) idum = -1
            endif
            if (iflgy .eq. 1) then
            if (jjp .eq. 0)       jdum = ndygp
            if (jjp .eq. ndygp-1) jdum = -1
            endif
            if (iflgz .eq. 1) then
            if (kkp .eq. 0)       kdum = ndzgp
            if (kkp .eq. ndzgp-1) kdum = -1
            endif

            iadd2(1) = idum
            iadd2(2) = jdum
            iadd2(3) = kdum

            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(2))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 1 .or. iflgsum .eq. 2) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndumn/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)
            
            idum = iip
            jdum = jjp
            kdum = kkp

            if (iflgx .eq. 1) then
            if (iip .eq. 0)       idum = ndxgp
            if (iip .eq. ndxgp-1) idum = -1
            endif
            if (iflgy .eq. 1) then
            if (jjp .eq. 0)       jdum = ndygp
            if (jjp .eq. ndygp-1) jdum = -1
            endif
            if (iflgz .eq. 1) then
            if (kkp .eq. 0)       kdum = ndzgp
            if (kkp .eq. ndzgp-1) kdum = -1
            endif

            iadd2(1) = idum
            iadd2(2) = jdum
            iadd2(3) = kdum

            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  333 continue
         enddo
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine create_ghost_particles_col
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in rptsgp(jgpx,j),rptsgp(jgpy,j),
c     rptsgp(jgpz,j), while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),ntypesl(7),iadd2(3),gpsave(27)


! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      xdlen = binb(2) - binb(1)
      ydlen = binb(4) - binb(3)
      zdlen = -1.
      if (if3d) zdlen = binb(6) - binb(5)
      if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) xdlen = -1
      if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) ydlen = -1
      if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      nfptsgp = 0

      do ip=1,n
         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip = ipart(jicx,ip) 
         jjp = ipart(jicy,ip) 
         kkp = ipart(jicz,ip) 

         rxl = binb(1) + rdxgp*iip
         rxr = rxl + rdxgp
         ryl = binb(3) + rdygp*jjp
         ryr = ryl + rdygp
         rzl = 0.0
         rzr = 0.0
         if (if3d) then
            rzl = binb(5) + rdzgp*kkp
            rzr = rzl + rdzgp
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle
               

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 2 .or. iflgsum .eq. 3) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndumn/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank


            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)

            idum = iip
            jdum = jjp
            kdum = kkp

            if (iflgx .eq. 1) then
            if (iip .eq. 0)       idum = ndxgp
            if (iip .eq. ndxgp-1) idum = -1
            endif
            if (iflgy .eq. 1) then
            if (jjp .eq. 0)       jdum = ndygp
            if (jjp .eq. ndygp-1) jdum = -1
            endif
            if (iflgz .eq. 1) then
            if (kkp .eq. 0)       kdum = ndzgp
            if (kkp .eq. ndzgp-1) kdum = -1
            endif

            iadd2(1) = idum
            iadd2(2) = jdum
            iadd2(3) = kdum
            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 1 .or. iflgsum .eq. 3) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndumn/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)
            
            idum = iip
            jdum = jjp
            kdum = kkp

            if (iflgx .eq. 1) then
            if (iip .eq. 0)       idum = ndxgp
            if (iip .eq. ndxgp-1) idum = -1
            endif
            if (iflgy .eq. 1) then
            if (jjp .eq. 0)       jdum = ndygp
            if (jjp .eq. ndygp-1) jdum = -1
            endif
            if (iflgz .eq. 1) then
            if (kkp .eq. 0)       kdum = ndzgp
            if (kkp .eq. ndzgp-1) kdum = -1
            endif

            iadd2(1) = idum
            iadd2(2) = jdum
            iadd2(3) = kdum

            
            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            distchk = 0.0
            dist = 0.0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (if3d) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (d2chk(3))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
c           if (iflgsum .eq. 1 .or. iflgsum .eq. 2) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg
            ndumn = ndumn/nreach
            nrank = modulo(ndumn,np)

            if (nrank .eq. nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
       
            call check_periodic_gp(rxnew,rxdrng,iadd)
       
            iadd(1)  = 0
            iadd(2)  = nrank
            iadd(3)  = ipart(ip,je0)
            
            idum = iip
            jdum = jjp
            kdum = kkp

            if (iflgx .eq. 1) then
            if (iip .eq. 0)       idum = ndxgp
            if (iip .eq. ndxgp-1) idum = -1
            endif
            if (iflgy .eq. 1) then
            if (jjp .eq. 0)       jdum = ndygp
            if (jjp .eq. ndygp-1) jdum = -1
            endif
            if (iflgz .eq. 1) then
            if (kkp .eq. 0)       kdum = ndzgp
            if (kkp .eq. ndzgp-1) kdum = -1
            endif

            iadd2(1) = idum
            iadd2(2) = jdum
            iadd2(3) = kdum

            call add_a_ghost_particle(rxnew,iadd,iadd2,ip)
                        
  333 continue
         enddo
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine check_periodic_gp(rxnew,rxdrng,iadd)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
c
      real rxnew(3), rxdrng(3)
      integer iadd(3), irett(3), ntype, ntypel(7)

      xloc = rxnew(1)
      yloc = rxnew(2)
      zloc = rxnew(3)

      xdlen = rxdrng(1)
      ydlen = rxdrng(2)
      zdlen = rxdrng(3)

      ii = iadd(1)
      jj = iadd(2)
      kk = iadd(3)

      irett(1) = 0
      irett(2) = 0
      irett(3) = 0

      if (xdlen .gt. 0 ) then
      if (ii .ge. ndxgp) then
         xloc = xloc - xdlen
         irett(1) = 1
         goto 123
      endif
      endif
      if (xdlen .gt. 0 ) then
      if (ii .lt. 0) then
         xloc = xloc + xdlen
         irett(1) = 1
         goto 123
      endif
      endif

  123 continue    
      if (ydlen .gt. 0 ) then
      if (jj .ge. ndygp) then
         yloc = yloc - ydlen
         irett(2) = 1
         goto 124
      endif
      endif
      if (ydlen .gt. 0 ) then
      if (jj .lt. 0) then
         yloc = yloc + ydlen
         irett(2) = 1
         goto 124
      endif
      endif
  124 continue

      if (if3d) then
         if (zdlen .gt. 0 ) then
         if (kk .ge. ndzgp) then
            zloc = zloc - zdlen
            irett(3) = 1
            goto 125
         endif
         endif
         if (zdlen .gt. 0 ) then
         if (kk .lt. 0) then
            zloc = zloc + zdlen
            irett(3) = 1
            goto 125
         endif
         endif
      endif
  125 continue

      rxnew(1) = xloc
      rxnew(2) = yloc
      rxnew(3) = zloc

      return
      end
c----------------------------------------------------------------------
      subroutine add_a_ghost_particle(rxnew,iadd,iadd2,i)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real    rxnew(3)
      integer iadd(3),iadd2(3)

      nfptsgp = nfptsgp + 1

      rptsgp(jgpx,nfptsgp)    = rxnew(1)       ! x loc
      rptsgp(jgpy,nfptsgp)    = rxnew(2)       ! y loc
      rptsgp(jgpz,nfptsgp)    = rxnew(3)       ! z loc
      rptsgp(jgpfh,nfptsgp)   = rpart(jf0,i)   ! hyd. force x
      rptsgp(jgpfh+1,nfptsgp) = rpart(jf0+1,i) ! hyd. force y
      rptsgp(jgpfh+2,nfptsgp) = rpart(jf0+2,i) ! hyd. force z
      rptsgp(jgpvol,nfptsgp)  = rpart(jvol,i)  ! particle volum
      rptsgp(jgprpe,nfptsgp)  = rpart(jrpe,i)  ! particle rp eff
      rptsgp(jgpspl,nfptsgp)  = rpart(jspl,i)  ! spl
      rptsgp(jgpg0,nfptsgp)   = rpart(jg0,i)   ! work done by forc
      rptsgp(jgpq0,nfptsgp)   = rpart(jq0,i)   ! heating from part 
      rptsgp(jgpv0,nfptsgp)   = rpart(jv0,i)   ! particle velocity
      rptsgp(jgpv0+1,nfptsgp) = rpart(jv0+1,i) ! particle velocity
      rptsgp(jgpv0+2,nfptsgp) = rpart(jv0+2,i) ! particle velocity

      iptsgp(jgpiic,nfptsgp)  = iadd(1)        ! use in collisions
      iptsgp(jgpps,nfptsgp)   = iadd(2)        ! overwritten mpi
      iptsgp(jgppt,nfptsgp)   = iadd(2)        ! dest. mpi rank
      iptsgp(jgpes,nfptsgp)   = iadd(3)        ! dest. elment
      iptsgp(jgpicx,nfptsgp)  = iadd2(1)
      iptsgp(jgpicy,nfptsgp)  = iadd2(2)
      iptsgp(jgpicz,nfptsgp)  = iadd2(3)
      iptsgp(jgppid1,nfptsgp) = ipart(jpid1,i)
      iptsgp(jgppid2,nfptsgp) = ipart(jpid2,i)
      iptsgp(jgppid3,nfptsgp) = ipart(jpid3,i)

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_bins
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl

      real pfx,pfy,pfz,vol,qgqf,rvx,rvy,rvz,rexp,multfc,multfci,rx2(3)
     >     ,rxyzp(6,3),d2new(3)
      integer ntypesl(7), ngp_trim(lbmax), ifac(3), icount(3)

      real    rxnew(3), rxdrng(3)
      integer iadd(3)


      logical partl

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (if3d) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ! compute binb
      xmin = 1E8
      ymin = 1E8
      zmin = 1E8
      xmax = 0.
      ymax = 0.
      zmax = 0.
      do i=1,n
         rduml = rpart(jx,i) - d2chk(1)
         rdumr = rpart(jx,i) + d2chk(1)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = rpart(jy,i) - d2chk(1)
         rdumr = rpart(jy,i) + d2chk(1)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         rduml = rpart(jz,i) - d2chk(1)
         rdumr = rpart(jz,i) + d2chk(1)
         if (rduml .lt. zmin) zmin = rduml
         if (rdumr .gt. zmax) zmax = rdumr
      enddo

      binb(1) = glmin(xmin,1)
      binb(2) = glmax(xmax,1)
      binb(3) = glmin(ymin,1)
      binb(4) = glmax(ymax,1)
      binb(5) = 0.0
      binb(6) = 0.0
      if(if3d) binb(5) = glmin(zmin,1)
      if(if3d) binb(6) = glmax(zmax,1)

      binb(1) = max(binb(1),xdrange(1,1))
      binb(2) = min(binb(2),xdrange(2,1))
      binb(3) = max(binb(3),xdrange(1,2))
      binb(4) = min(binb(4),xdrange(2,2))
      if(if3d)binb(5) = max(binb(5),xdrange(1,3))
      if(if3d)binb(6) = min(binb(6),xdrange(2,3))

      if (bc_part(1) .eq. 0) then
         binb(1) = xdrange(1,1)
         binb(2) = xdrange(2,1)
      endif
      if (bc_part(3) .eq. 0) then
         binb(3) = xdrange(1,2)
         binb(4) = xdrange(2,2)
      endif
      if (bc_part(5) .eq. 0 .and. if3d) then
         binb(5) = xdrange(1,3)
         binb(6) = xdrange(2,3)
      endif

      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      icount(1) = 0
      icount(2) = 0
      icount(3) = 0
      d2new(1) = d2chk(1)
      d2new(2) = d2chk(1)
      d2new(3) = d2chk(1)

      ndxgp = floor( (binb(2) - binb(1))/d2new(1))
      ndygp = floor( (binb(4) - binb(3))/d2new(2))
      ndzgp = 1
      if (if3d) ndzgp = floor( (binb(6) - binb(5))/d2new(3))

      if (ndxgp*ndygp*ndzgp .gt. np) then
         nmax = 1000
         d2chk_save = d2chk(1)
         
         do i=1,nmax
         do j=0,ndim-1
            ifac(j+1) = 1 + i
            d2new(j+1) = (binb(2+2*j) - binb(1+2*j))/ifac(j+1)
            nbb = ifac(1)*ifac(2)*ifac(3)
            if(d2new(j+1) .lt. d2chk_save .or. nbb .gt. np) then
               icount(j+1) = icount(j+1) + 1
               ifac(j+1) = ifac(j+1) - icount(j+1)
               d2new(j+1) = (binb(2+2*j) - binb(1+2*j))/ifac(j+1)
            endif
         enddo
            if (icount(1) .gt. 0) then
            if (icount(2) .gt. 0) then
            if (icount(3) .gt. 0) then
               exit
            endif
            endif
            endif
         enddo
      endif


! -------------------------------------------------------
c SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! how many spacings in each direction
      ndxgp = floor( (binb(2) - binb(1))/d2new(1))
      ndygp = floor( (binb(4) - binb(3))/d2new(2))
      ndzgp = 1
      if (if3d) ndzgp = floor( (binb(6) - binb(5))/d2new(3))

      if (nid .eq. 0) write(6,*) ndxgp, ndygp, ndzgp, ifac, d2new

      ! grid spacing for that many spacings
      rdxgp = (binb(2) - binb(1))/real(ndxgp)
      rdygp = (binb(4) - binb(3))/real(ndygp)
      rdzgp = 1.
      if (if3d) rdzgp = (binb(6) - binb(5))/real(ndzgp)

      ninc = 2
      if (bc_part(1) .ne. 0) then
         binb(1) = binb(1) - ninc/2*rdxgp
         binb(2) = binb(2) + ninc/2*rdxgp
         ndxgp = ndxgp + ninc
         ic = 0
         if (binb(1) .lt. xdrange(1,1)) then
            ic = ic + 1
            binb(1) = xdrange(1,1)
         endif
         if (binb(2) .gt. xdrange(2,1)) then
            ic = ic + 1
            binb(2) = xdrange(2,1)
         endif
         ndxgp = ndxgp - ic
      endif
      if (bc_part(3) .ne. 0) then
         binb(3) = binb(3) - ninc/2*rdygp
         binb(4) = binb(4) + ninc/2*rdygp
         ndygp = ndygp + ninc
         ic = 0
         if (binb(3) .lt. xdrange(1,2)) then
            ic = ic + 1
            binb(3) = xdrange(1,2)
         endif
         if (binb(4) .gt. xdrange(2,2)) then
            ic = ic + 1
            binb(4) = xdrange(2,2)
         endif
         ndygp = ndygp - ic
      endif
      if (bc_part(5) .ne. 0) then
      if (if3d) then
         binb(5) = binb(5) - ninc/2*rdzgp
         binb(6) = binb(6) + ninc/2*rdzgp
         ndzgp = ndzgp + ninc
         ic = 0
         if (binb(5) .lt. xdrange(1,3)) then
            ic = ic + 1
            binb(5) = xdrange(1,3)
         endif
         if (binb(6) .gt. xdrange(2,3)) then
            ic = ic + 1
            binb(6) = xdrange(2,3)
         endif
         ndzgp = ndzgp - ic
      endif
      endif

      nbin_now = ndxgp*ndygp*ndzgp
c     if (nbin_now .eq. nbin_save) return ! except not all the time..

      nreach = nbin_new/np
      nreach = max(1,nreach)

c     if (nid.eq.0) write(6,*) 'Setting up bins'
c     if (nid.eq.0) write(6,*) ndxgp,ndygp,ndzgp,d2chk(1), np,nbin_now

      if (npro_method .eq. 0) return ! only for projection

! see which bins are in which elements
      neltb = 0
      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rxval = xm1(i,j,k,ie)
         ryval = ym1(i,j,k,ie)
         rzval = 0.
         if(if3d) rzval = zm1(i,j,k,ie)

         if (rxval .gt. binb(2)) goto 1233
         if (rxval .lt. binb(1)) goto 1233
         if (ryval .gt. binb(4)) goto 1233
         if (ryval .lt. binb(3)) goto 1233
         if (if3d .and. rzval .gt. binb(6)) goto 1233
         if (if3d .and. rzval .lt. binb(5)) goto 1233

         ii    = floor((rxval-binb(1))/rdxgp) 
         jj    = floor((ryval-binb(3))/rdygp) 
         kk    = floor((rzval-binb(5))/rdzgp) 
         if (.not. if3d) kk = 0
         if (ii .eq. ndxgp) ii = ndxgp - 1
         if (jj .eq. ndygp) jj = ndygp - 1
         if (kk .eq. ndzgp) kk = ndzgp - 1
         ndum  = ii + ndxgp*jj + ndxgp*ndygp*kk
         ndum = ndum/nreach
         nrank = modulo(ndum,np)

         neltb = neltb + 1
         if(neltb .gt. lbmax) then
           write(6,*) 'increase lbmax',nid,neltb,lbmax
           call exitt
         endif

c        ptw(i,j,k,ie,1) = ii
c        ptw(i,j,k,ie,2) = jj
c        ptw(i,j,k,ie,3) = kk
c        ptw(i,j,k,ie,4) = ndum

         er_map(1,neltb) = ie
         er_map(2,neltb) = nid
         er_map(3,neltb) = ndum
         er_map(4,neltb) = nrank
         er_map(5,neltb) = nrank
         er_map(6,neltb) = nrank

         if (neltb .gt. 1) then
         do il=1,neltb-1
            if (er_map(1,il) .eq. ie) then
            if (er_map(4,il) .eq. nrank) then
               neltb = neltb - 1
               goto 1233
            endif
            endif
         enddo
         endif


 1233 continue
      enddo
      enddo
      enddo
      enddo

c     call output_two_way_io
c     call exitt


      nxyz = lx1*ly1*lz1
      do ie=1,neltb
         iee = er_map(1,ie)
         call copy(xm1b(1,1,1,1,ie), xm1(1,1,1,iee),nxyz)
         call copy(xm1b(1,1,1,2,ie), ym1(1,1,1,iee),nxyz)
         call copy(xm1b(1,1,1,3,ie), zm1(1,1,1,iee),nxyz)
      enddo

      neltbb = neltb

      do ie=1,neltbb
         call icopy(er_maps(1,ie),er_map(1,ie),6)
      enddo


      nl = 0
      nii = 6
      njj = 6
      nrr = nxyz*3
      call fgslib_crystal_tuple_transfer(i_cr_hndl,neltb,lbmax
     >                  , er_map,nii,partl,nl,xm1b,nrr,njj)


      call partition_bins

c     call output_two_way_io
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine partition_bins
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl

      logical partl

      real xerange(2,3,lbmax)
      common /lpm_elementrange/ xerange

      integer irmv(lbmax)

         nxyz = lx1*ly1*lz1

c     if (nbin_now .lt. np) then

c        do i=1,np
c           lpm_binsort(1,i) = i-1
c           lpm_binsort(2,i) = n
c           lpm_binsort(3,i) = i-1
c        enddo
c        
c        nps   = 1 ! index of new proc for doing stuff
c        nglob = 2 ! unique key to sort by
c        nkey  = 1 ! number of keys (just 1 here)
c        ndum = 3
c        call fgslib_crystal_ituple_transfer(i_cr_hndl,lpm_binsort,
c    >                 ndum,np,lpmax,nps)
c        call fgslib_crystal_ituple_sort(i_cr_hndl,lpm_binsort,
c    >                 ndum,np,nglob,nkey)

c        call bcast(lpm_binsort, 3*np*isize)

c        np_high = np
c        np_low  = np - nbin_now +1
c        np_low  = max(np_low,nbin_now)

c        if (n .gt. 0) then
c        do i=np_high,np_low,-1
c           if (lpm_binsort(1,i) .eq. nid) then
c              ii = np_high + 1 - i
c              neltb_save = neltb
c              do ie=1,neltb_save
c                 neltb = neltb + 1
c                 call icopy(er_map(1,neltb),er_map(1,ie),4)
c                 er_map(5,ie)    = er_map(4,ie)
c                 er_map(6,ie)    = er_map(4,ie)
c                 er_map(5,neltb) = lpm_binsort(1,ii)
c                 er_map(6,neltb) = lpm_binsort(1,ii)

c                 call copy(xm1b(1,1,1,1,neltb), xm1b(1,1,1,1,ie),nxyz)
c                 call copy(xm1b(1,1,1,2,neltb), xm1b(1,1,1,2,ie),nxyz)
c                 call copy(xm1b(1,1,1,3,neltb), xm1b(1,1,1,3,ie),nxyz)
c              enddo
c           endif
c        enddo
c        endif

c     endif

c     nl = 0
c     nii = 6
c     njj = 6
c     nrr = nxyz*3
c     nkey = 3
c     call fgslib_crystal_tuple_transfer(i_cr_hndl,neltb,lbmax
c    >                  , er_map,nii,partl,nl,xm1b,nrr,njj)


      nl = 0
      nii = 6
      njj = 6
      nrr = nxyz*3
      nkey = 3
      call fgslib_crystal_tuple_sort    (i_cr_hndl,neltb
     $              , er_map,nii,partl,nl,xm1b,nrr,nkey,1)



c     nd_min = iglmin(neltb,1)
c     nd_max = iglmax(neltb,1)
c     nd_sum = iglsum(neltb,1)
c     if (nid.eq.0) write(6,*) 'Finished setting up bins'
c     if (nid.eq.0) write(6,*) ndxgp,ndygp,ndzgp,nd_sum,nd_min,nd_max
      
      do ie=1,neltb
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rxval = xm1b(i,j,k,1,ie)
         ryval = xm1b(i,j,k,2,ie)
         rzval = 0.
         if(if3d) rzval = xm1b(i,j,k,3,ie)
         
         ii    = floor((rxval-binb(1))/rdxgp) 
         jj    = floor((ryval-binb(3))/rdygp) 
         kk    = floor((rzval-binb(5))/rdzgp) 
         if (.not. if3d) kk = 0
         if (ii .eq. ndxgp) ii = ndxgp - 1
         if (jj .eq. ndygp) jj = ndygp - 1
         if (kk .eq. ndzgp) kk = ndzgp - 1
         ndum  = ii + ndxgp*jj + ndxgp*ndygp*kk

c        ptw(i,j,k,er_map(1,ie),1) = ii
c        ptw(i,j,k,er_map(1,ie),2) = jj
c        ptw(i,j,k,er_map(1,ie),3) = kk
c        ptw(i,j,k,er_map(1,ie),4) = ndum
      
         mod_gp_grid(i,j,k,ie,1) = ii
         mod_gp_grid(i,j,k,ie,2) = jj
         mod_gp_grid(i,j,k,ie,3) = kk
         mod_gp_grid(i,j,k,ie,4) = ndum
   
      enddo
      enddo
      enddo
      enddo

c     call output_two_way_io
c     call exitt

      do ie=1,neltb
         xerange(1,1,ie) = vlmin(xm1b(1,1,1,1,ie),nxyz)
         xerange(2,1,ie) = vlmax(xm1b(1,1,1,1,ie),nxyz)
         xerange(1,2,ie) = vlmin(xm1b(1,1,1,2,ie),nxyz)
         xerange(2,2,ie) = vlmax(xm1b(1,1,1,2,ie),nxyz)
         xerange(1,3,ie) = vlmin(xm1b(1,1,1,3,ie),nxyz)
         xerange(2,3,ie) = vlmax(xm1b(1,1,1,3,ie),nxyz)

         ilow  = floor((xerange(1,1,ie) - binb(1))/rdxgp)
         ihigh = floor((xerange(2,1,ie) - binb(1))/rdxgp)
         jlow  = floor((xerange(1,2,ie) - binb(3))/rdygp)
         jhigh = floor((xerange(2,2,ie) - binb(3))/rdygp)
         klow  = floor((xerange(1,3,ie) - binb(5))/rdzgp)
         khigh = floor((xerange(2,3,ie) - binb(5))/rdzgp)
         if (.not. if3d) then
            klow = 0
            khigh = 0
         endif

         lpm_el_map(1,ie) = ilow  + ndxgp*jlow  + ndxgp*ndygp*klow
         lpm_el_map(2,ie) = ihigh + ndxgp*jhigh + ndxgp*ndygp*khigh
         lpm_el_map(3,ie) = ilow
         lpm_el_map(4,ie) = ihigh
         lpm_el_map(5,ie) = jlow
         lpm_el_map(6,ie) = jhigh
         lpm_el_map(7,ie) = klow
         lpm_el_map(8,ie) = khigh
      enddo



      nbin_save = nbin_now


      return
      end
c-----------------------------------------------------------------------
      subroutine extra_wall_particle_exp(rxyzp,rx2,ic)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
c
      real rxyzp(n_walls*2,3), rx2(3), rx22(3)
      integer ic, ip

      ! plane wall collisions
      ic = 0
      do j = 1,np_walls
         rnx = plane_wall_coords(1,j)
         rny = plane_wall_coords(2,j)
         rnz = plane_wall_coords(3,j)
         rpx = plane_wall_coords(4,j)
         rpy = plane_wall_coords(5,j)
         rpz = 1.0
         if (if3d) rpz = plane_wall_coords(6,j)

         rd    = -(rnx*rpx + rny*rpy + rnz*rpz)

         rdist = abs(rnx*rx2(1)+rny*rx2(2)+rnz*rx2(3)+rd)
         rdist = rdist*2.
         rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

         if (rdist .gt. d2chk(1)) cycle
         ic = ic + 1

         rxyzp(ic,1) = rx2(1) - rdist*rnx
         rxyzp(ic,2) = rx2(2) - rdist*rny
         rxyzp(ic,3) = rx2(3) - rdist*rnz
      enddo

      ! cylinder wall collisions
      do j = 1,nc_walls
         rnx = cyl_wall_coords(1,j)
         rny = cyl_wall_coords(2,j)
         rnz = cyl_wall_coords(3,j)
         rpx = cyl_wall_coords(4,j)
         rpy = cyl_wall_coords(5,j)
         rpz = 1.0
         if (if3d) rpz = cyl_wall_coords(6,j)
         rrad = cyl_wall_coords(7,j)

         rx22(1) = rx2(1)
         rx22(2) = rx2(2)
         rx22(3) = rx2(3)
         ! for now only works with cylinders aligned with axes at
         ! origin
         if (rnz .gt. 0.5) then
            rtheta = atan2(rx22(2),rx22(1))
            rx22(1) = rpx + rrad*cos(rtheta)
            rx22(2) = rpy + rrad*sin(rtheta)
         elseif (rnx .gt. 0.5) then
            rtheta = atan2(rx22(3),rx22(2))
            rx22(2) = rpy + rrad*cos(rtheta)
            rx22(3) = rpz + rrad*sin(rtheta)
         elseif (rny .gt. 0.5) then
            rtheta = atan2(rx22(1),rx22(3))
            rx22(3) = rpz + rrad*cos(rtheta)
            rx22(1) = rpx + rrad*sin(rtheta)
         endif

         rx2d = rx22(1) - rx2(1)
         ry2d = rx22(2) - rx2(2)
         rz2d = rx22(3) - rx2(3)

         rdist = sqrt(rx2d**2 + ry2d**2 + rz2d**2)
         rx2d = rx2d/rdist
         ry2d = ry2d/rdist
         rz2d = rz2d/rdist

         rdist = rdist*2.

         if (rdist .gt. d2chk(1)) cycle
         ic = ic + 1

         rxyzp(ic,1) = rx2(1) + rx2d*rdist
         rxyzp(ic,2) = rx2(2) + ry2d*rdist
         rxyzp(ic,3) = rx2(3) + rz2d*rdist

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_bdf_ext_coefs(beta,alpha,times)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real beta(0:3),alpha(0:3),times(0:3)
      real c(0:8)

      integer ilast,ncoef
      save    ilast,ncoef
      data    ilast,ncoef / -9 , 0 /

      do i=3,1,-1
         times(i)=times(i-1)
      enddo
      times(0) = time

      call rzero(beta ,4)
      call rzero(alpha,4)
      if (istep.ne.ilast) then
         ilast = istep
         ncoef = ncoef + 1
         ncoef = min(ncoef,3) ! Maximum 3rd order in time
      endif
      ncoefm1 = ncoef - 1

      call fd_weights_full(times(0),times(1),ncoefm1,0,alpha(1))
      call fd_weights_full(times(0),times(0),ncoef,1,c)
      do j=0,ncoef
         beta(j) = c(ncoef+1+j)
      enddo
      do j=1,ncoef
         beta(j) = -beta(j)  ! Change sign, for convenience
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine update_particle_location
c     check if particles are outside domain
c     > if bc_part = 0 then it is periodic
c     > if bc_part = -1,1 then particles are killed (outflow)
      include 'SIZE'
      include 'LPM'

      integer in_part(llpart)

      jx0 = jx

      do i=1,n
         in_part(i) = 0
         do j=0,ndim-1
            if (rpart(jx0+j,i).lt.xdrange(1,j+1))then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx0+j,i) = xdrange(2,j+1) - 
     &                             abs(xdrange(1,j+1) - rpart(jx0+j,i))
                  rpart(jx1+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx1+j,i))
                  rpart(jx2+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx2+j,i))
                  rpart(jx3+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx3+j,i))
                  goto 1512
                endif
            endif
            if (rpart(jx0+j,i).gt.xdrange(2,j+1))then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx0+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx0+j,i) - xdrange(2,j+1))
                  rpart(jx1+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx1+j,i) - xdrange(2,j+1))
                  rpart(jx2+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx2+j,i) - xdrange(2,j+1))
                  rpart(jx3+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx3+j,i) - xdrange(2,j+1))
                  goto 1512
                endif
            endif
            if (ipart(jrc,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get here
            endif
 1512 continue
         enddo
      enddo

      ic = 0
      do i=1,n
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               call copy(rpart(1,ic),rpart(1,i),nr)
               call icopy(ipart(1,ic),ipart(1,i),ni)
            endif
         endif
      enddo
      n = ic

      return
      end
c-----------------------------------------------------------------------
c     interpolation routines
c-----------------------------------------------------------------------
      subroutine interp_props_part_location
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'LPM'
      include 'GEOM'

      common /myparth/ i_fp_hndl, i_cr_hndl

      if (red_interp .eq. 0) then

      elseif (red_interp .eq. 1) then

      if (ifrectp .eq. 0) then
      ! set location but don't move...
      call fgslib_findpts(i_fp_hndl !  stride     !   call fgslib_findpts( ihndl,
     $        , ipart(jrc,1),li        !   $             rcode,1,
     $        , ipart(jpt,1),li        !   &             proc,1,
     $        , ipart(je0,1),li        !   &             elid,1,
     $        , rpart(jr ,1),lr        !   &             rst,ndim,
     $        , rpart(jd ,1),lr        !   &             dist,1,
     $        , rpart(jx ,1),lr        !   &             pts(    1),1,
     $        , rpart(jy ,1),lr        !   &             pts(  n+1),1,
     $        , rpart(jz ,1),lr ,n)    !   &             pts(2*n+1),1,n)

      ! when data is non-local
      call fgslib_findpts_eval(i_fp_hndl,rpart(ju0+0,1),nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)  ,ni,
     &                             rpart(jr,1)   ,nr,n,
     &                             vx)
      call fgslib_findpts_eval(i_fp_hndl,rpart(ju0+1,1),nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)  ,ni,
     &                             rpart(jr,1)   ,nr,n,
     &                             vy)
      call fgslib_findpts_eval(i_fp_hndl,rpart(ju0+2,1),nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)  ,ni,
     &                             rpart(jr,1)   ,nr,n,
     &                             vz)
      call fgslib_findpts_eval(i_fp_hndl,rpart(jrho,1)  ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)   ,ni,
     &                             rpart(jr,1)    ,nr,n,
     &                             vtrans)
      call fgslib_findpts_eval(i_fp_hndl,rpart(jtempf,1)  ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)   ,ni,
     &                             rpart(jr,1)    ,nr,n,
     &                             t)
      call fgslib_findpts_eval(i_fp_hndl,rpart(jDuDt+0,1) ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)     ,ni,
     &                             rpart(jr,1)      ,nr,n,
     &                             rhs_fluidp(1,1,1,1,1))
      call fgslib_findpts_eval(i_fp_hndl,rpart(jDuDt+1,1) ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)     ,ni,
     &                             rpart(jr,1)      ,nr,n,
     &                             rhs_fluidp(1,1,1,1,2))
      call fgslib_findpts_eval(i_fp_hndl,rpart(jDuDt+2,1) ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)     ,ni,
     &                             rpart(jr,1)      ,nr,n,
     &                             rhs_fluidp(1,1,1,1,3))
      call fgslib_findpts_eval(i_fp_hndl,rpart(jvol1,1)   ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)     ,ni,
     &                             rpart(jr,1)      ,nr,n,
     &                             ptw(1,1,1,1,4))


c        call dzlib_interp

      else
         call dzlib_findpts
         call dzlib_nflist


c        jps = jpid1-1     ! Pointer to temporary proc id for swapping
c        do i=1,n        ! Can't use jpt because it messes up particle info
c           ipart(jps,i) = ipart(jpt,i)
c        enddo
c        call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
c    $           , ipart,ni,partl,nl,rpart,nr,jps)

c        nplow = iglmin(n,1)
c        nptot = iglsum(n,1)
c        nphigh = iglmax(n,1)
c        if(nid.eq.0)write(6,*)'interp diag      :',nptot,nplow,nphigh

c     call fgslib_crystal_tuple_sort    (i_cr_hndl,n
c    $              , ipart,ni,partl,0,rpart,nr,je0,1)

c     ! note that all data is local now
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+0,1),nr,
c    &                                   ipart(je0,1)  ,ni,
c    &                                   rpart(jr,1)   ,nr,n,
c    &                                   vx)
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+1,1),nr,
c    &                                   ipart(je0,1)  ,ni,
c    &                                   rpart(jr,1)   ,nr,n,
c    &                                   vy)
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+2,1),nr,
c    &                                   ipart(je0,1)  ,ni,
c    &                                   rpart(jr,1)   ,nr,n,
c    &                                   vz)
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(jrho,1)  ,nr,
c    &                                   ipart(je0,1)   ,ni,
c    &                                   rpart(jr,1)    ,nr,n,
c    &                                   vtrans)
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(jDuDt+0,1) ,nr,
c    &                                   ipart(je0,1)     ,ni,
c    &                                   rpart(jr,1)      ,nr,n,
c    &                                   rhs_fluidp(1,1,1,1,1))
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(jDuDt+1,1) ,nr,
c    &                                   ipart(je0,1)     ,ni,
c    &                                   rpart(jr,1)      ,nr,n,
c    &                                   rhs_fluidp(1,1,1,1,2))
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(jDuDt+2,1) ,nr,
c    &                                   ipart(je0,1)     ,ni,
c    &                                   rpart(jr,1)      ,nr,n,
c    &                                   rhs_fluidp(1,1,1,1,3))
c     call fgslib_findpts_eval_local(i_fp_hndl,rpart(jvol1,1)   ,nr,
c    &                                   ipart(je0,1)     ,ni,
c    &                                   rpart(jr,1)      ,nr,n,
c    &                                   ptw(1,1,1,1,4))


c     call move_particles_inproc

         call dzlib_interp
      endif

      ! when data is non-local
c     call fgslib_findpts_eval(i_fp_hndl,rpart(ju0+0,1),nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)  ,ni,
c    &                             rpart(jr,1)   ,nr,n,
c    &                             vx)
c     call fgslib_findpts_eval(i_fp_hndl,rpart(ju0+1,1),nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)  ,ni,
c    &                             rpart(jr,1)   ,nr,n,
c    &                             vy)
c     call fgslib_findpts_eval(i_fp_hndl,rpart(ju0+2,1),nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)  ,ni,
c    &                             rpart(jr,1)   ,nr,n,
c    &                             vz)
c     call fgslib_findpts_eval(i_fp_hndl,rpart(jrho,1)  ,nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)   ,ni,
c    &                             rpart(jr,1)    ,nr,n,
c    &                             vtrans)
c     call fgslib_findpts_eval(i_fp_hndl,rpart(jDuDt+0,1) ,nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)     ,ni,
c    &                             rpart(jr,1)      ,nr,n,
c    &                             rhs_fluidp(1,1,1,1,1))
c     call fgslib_findpts_eval(i_fp_hndl,rpart(jDuDt+1,1) ,nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)     ,ni,
c    &                             rpart(jr,1)      ,nr,n,
c    &                             rhs_fluidp(1,1,1,1,2))
c     call fgslib_findpts_eval(i_fp_hndl,rpart(jDuDt+2,1) ,nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)     ,ni,
c    &                             rpart(jr,1)      ,nr,n,
c    &                             rhs_fluidp(1,1,1,1,3))
c     call fgslib_findpts_eval(i_fp_hndl,rpart(jvol1,1)   ,nr,
c    &                             ipart(jrc,1)  ,ni,
c    &                             ipart(jpt,1)  ,ni,
c    &                             ipart(je0,1)     ,ni,
c    &                             rpart(jr,1)      ,nr,n,
c    &                             ptw(1,1,1,1,4))


      ! check if values outside reasonable bounds
      rvfmax = 0.74
      rvfmin = 1E-4
      do i=1,n
         if (rpart(jvol1,i) .lt. rvfmin) rpart(jvol1,i) = rvfmin
         if (rpart(jvol1,i) .gt. rvfmax) rpart(jvol1,i) = rvfmax
      enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dzlib_interp
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xerange(2,3,lbmax)
      common /lpm_elementrange/ xerange

      common /myparth/ i_fp_hndl, i_cr_hndl

      logical partl

      integer infpts_map(llpart),nfpts1
      common /infc_list/ nfpts1,infpts_map

      ! first send grid data to bins

      nxyz = lx1*ly1*lz1
      do ie=1,neltbb
         iee = er_maps(1,ie)
         call copy(ptwb(1,1,1,1,ie),vx(1,1,1,iee)          ,nxyz)
         call copy(ptwb(1,1,1,2,ie),vy(1,1,1,iee)          ,nxyz)
         call copy(ptwb(1,1,1,3,ie),vz(1,1,1,iee)          ,nxyz)
         call copy(ptwb(1,1,1,4,ie),vtrans(1,1,1,iee,1)    ,nxyz)
         call copy(ptwb(1,1,1,5,ie),t(1,1,1,iee,1)         ,nxyz)
         call copy(ptwb(1,1,1,6,ie),rhs_fluidp(1,1,1,iee,1),nxyz)
         call copy(ptwb(1,1,1,7,ie),rhs_fluidp(1,1,1,iee,2),nxyz)
         call copy(ptwb(1,1,1,8,ie),rhs_fluidp(1,1,1,iee,3),nxyz)
         call copy(ptwb(1,1,1,9,ie),ptw(1,1,1,iee,4)       ,nxyz)
      enddo

      neltbc = neltbb
      ndum = 6*neltbc
      call icopy(er_mapc,er_maps,ndum)

      nl = 0
      nii = 6
      njj = 6
      nrr = nxyz*nfldpm
      call fgslib_crystal_tuple_transfer(i_cr_hndl,neltbc,lbmax
     >                  , er_mapc,nii,partl,nl,ptwb,nrr,njj)


      do i=1,n

         if (infpts_map(i) .gt. 0) cycle
         iee   = ipart(je0,i) + 1
         nrank = ipart(jpt,i)

         rrdum = 1.0
         if(if3d) rrdum = rpart(jr+2,i)

         do ie=1,neltbc
            if (er_mapc(1,ie) .eq. iee  ) then
            if (er_mapc(2,ie) .eq. nrank) then

            call init_baryinterp(rpart(jr,i),rpart(jr+1,i),rrdum,nxyz)

            call baryinterp(ptwb(1,1,1,1,ie),rpart(ju0+0,i)  ,nxyz)
            call baryinterp(ptwb(1,1,1,2,ie),rpart(ju0+1,i)  ,nxyz)
            call baryinterp(ptwb(1,1,1,3,ie),rpart(ju0+2,i)  ,nxyz)
            call baryinterp(ptwb(1,1,1,4,ie),rpart(jrho,i)   ,nxyz)
            call baryinterp(ptwb(1,1,1,5,ie),rpart(jtempf,i) ,nxyz)
            call baryinterp(ptwb(1,1,1,6,ie),rpart(jDuDt+0,i),nxyz)
            call baryinterp(ptwb(1,1,1,7,ie),rpart(jDuDt+1,i),nxyz)
            call baryinterp(ptwb(1,1,1,8,ie),rpart(jDuDt+2,i),nxyz)
            call baryinterp(ptwb(1,1,1,9,ie),rpart(jvol1,i)  ,nxyz)

            endif
            endif

         enddo
      enddo



      




      return
      end
c-----------------------------------------------------------------------
      subroutine baryinterp(field,pofx,nxyz)
      ! used for 3d interpolation only
      include 'SIZE'
      common /BARYPARAMS/ xgll, ygll, zgll, wxgll, wygll, wzgll
      real xgll(lx1), ygll(ly1), zgll(lz1),
     >     wxgll(lx1), wygll(ly1), wzgll(lz1)
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot
      real field(1),pofx,top

      pofx = 0.00
      do i=1,nxyz
         pofx =  pofx + rep(i,1,1)*field(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine init_baryinterp(x,y,z,nxyz)
      ! used for 3d interpolation only
      include 'SIZE'
      common /BARYPARAMS/ xgll, ygll, zgll, wxgll, wygll, wzgll
      real xgll(lx1), ygll(ly1), zgll(lz1),
     >     wxgll(lx1), wygll(ly1), wzgll(lz1)
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot

      real x, y, z, repy, repz,repx,diff
      real bwgtx(lx1),bwgty(ly1),bwgtz(lz1)

      bot= 0.00
      do k=1,nz1
         diff = z - zgll(k)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgtz(k) = wzgll(k)/diff
      enddo
      do i=1,nx1
         diff = x - xgll(i)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgtx(i) = wxgll(i)/diff
      enddo 
      do j=1,nx1
         diff = y-ygll(j)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgty(j) = wygll(j)/diff
      enddo

      do k=1,nz1
      do j=1,ny1
         repdum = bwgty(j)*bwgtz(k)
      do i=1,nx1
         rep(i,j,k) =  repdum* bwgtx(i)
         bot        =  bot + rep(i,j,k)
      enddo
      enddo
      enddo 

      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rep(i,j,k) =  rep(i,j,k)/bot
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine init_interpolation
      include 'SIZE' 
      include 'INPUT' 
 
      ! calculates the barycentric lagrange weights
 
      common /BARYPARAMS/ xgll, ygll, zgll, wxgll, wygll, wzgll
      real xgll(lx1), ygll(ly1), zgll(lz1),
     >     wxgll(lx1), wygll(ly1), wzgll(lz1)

c     get gll points in all directions
      call zwgll(xgll,wxgll,lx1)
      call zwgll(ygll,wygll,ly1)
      call rone(zgll,lz1)
      if(if3d) call zwgll(zgll,wzgll,lz1)
c     set all weights to ones first
      call rone(wxgll,lx1)
      call rone(wygll,ly1)
      call rone(wzgll,lz1)
 
      ! calc x bary weights
      do j=1,nx1
         do k=1,nx1
            if (j .NE. k) then
               wxgll(j) = wxgll(j)/(xgll(j) - xgll(k))
            endif
         enddo
      enddo
      ! calc y bary weights
      do j=1,nx1
         do k=1,nx1
            if (j .NE. k) then
               wygll(j) = wygll(j)/(ygll(j) - ygll(k))
            endif
         enddo
      enddo
      ! calc z bary weights
      if (if3d) then
      do j=1,nx1
         do k=1,nx1
            if (j .NE. k) then
               wzgll(j) = wzgll(j)/(zgll(j) - zgll(k))
            endif
         enddo
      enddo
      endif

      return 
      end
c-----------------------------------------------------------------------
      subroutine dzlib_nflist
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xerange(2,3,lbmax)
      common /lpm_elementrange/ xerange

      real rpart_nf(lr,llpart)
      common /rpart_nfc/ rpart_nf
      integer ipart_nf(li,llpart)
      common /ipart_nfc/ ipart_nf
      integer infpts_map(llpart),nfpts1
      common /infc_list/ nfpts1,infpts_map

      common /myparth/ i_fp_hndl, i_cr_hndl

      nfpts1 = 0
      do i=1,n
         infpts_map(i) = -1
         if (ipart(jrc,i) .eq. 2) then
            nfpts1 = nfpts1 + 1
            call copy(rpart_nf(1,nfpts1),rpart(1,i),nr)
            call icopy(ipart_nf(1,nfpts1),ipart(1,i),ni)
            infpts_map(nfpts1) = i
         endif
      enddo

      call fgslib_findpts(i_fp_hndl !  stride     !   call fgslib_findpts( ihndl,
     $        , ipart_nf(jrc,1),li        !   $             rcode,1,
     $        , ipart_nf(jpt,1),li        !   &             proc,1,
     $        , ipart_nf(je0,1),li        !   &             elid,1,
     $        , rpart_nf(jr ,1),lr        !   &             rst,ndim,
     $        , rpart_nf(jd ,1),lr        !   &             dist,1,
     $        , rpart_nf(jx ,1),lr        !   &             pts(    1),1,
     $        , rpart_nf(jy ,1),lr        !   &             pts(  n+1),1,
     $        , rpart_nf(jz ,1),lr ,nfpts1)    !   &             pts(2*n+1),1,n)

      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(ju0+0,1),nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)  ,ni,
     &                             rpart_nf(jr,1)   ,nr,nfpts1,
     &                             vx)
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(ju0+1,1),nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)  ,ni,
     &                             rpart_nf(jr,1)   ,nr,nfpts1,
     &                             vy)
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(ju0+2,1),nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)  ,ni,
     &                             rpart_nf(jr,1)   ,nr,nfpts1,
     &                             vz)
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(jrho,1)  ,nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)   ,ni,
     &                             rpart_nf(jr,1)    ,nr,nfpts1,
     &                             vtrans)
      call fgslib_findpts_eval(i_fp_hndl,rpart(jtempf,1)  ,nr,
     &                             ipart(jrc,1)  ,ni,
     &                             ipart(jpt,1)  ,ni,
     &                             ipart(je0,1)   ,ni,
     &                             rpart(jr,1)    ,nr,n,
     &                             t)
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(jDuDt+0,1) ,nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)     ,ni,
     &                             rpart_nf(jr,1)      ,nr,nfpts1,
     &                             rhs_fluidp(1,1,1,1,1))
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(jDuDt+1,1) ,nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)     ,ni,
     &                             rpart_nf(jr,1)      ,nr,nfpts1,
     &                             rhs_fluidp(1,1,1,1,2))
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(jDuDt+2,1) ,nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)     ,ni,
     &                             rpart_nf(jr,1)      ,nr,nfpts1,
     &                             rhs_fluidp(1,1,1,1,3))
      call fgslib_findpts_eval(i_fp_hndl,rpart_nf(jvol1,1)   ,nr,
     &                             ipart_nf(jrc,1)  ,ni,
     &                             ipart_nf(jpt,1)  ,ni,
     &                             ipart_nf(je0,1)     ,ni,
     &                             rpart_nf(jr,1)      ,nr,nfpts1,
     &                             ptw(1,1,1,1,4))

      do i=1,nfpts1
         j = infpts_map(i) 
         call copy(rpart(1,j),rpart_nf(1,i),nr)
         call icopy(ipart(1,j),ipart_nf(1,i),ni)
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine dzlib_findpts
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xerange(2,3,lbmax)
      common /lpm_elementrange/ xerange

      logical iztrue

      ! must set ipart of jpt,je0 and rpart of jr,js,jt

      ifnd = 0
      do ip = 1,n
         xloc = rpart(jx,ip)
         yloc = rpart(jy,ip)
         zloc = rpart(jz,ip)
         ipart(jrc ,ip) = 2
         do ie=1,neltb
            iztrue = .false.
            if (zloc.ge.xerange(1,3,ie).and.zloc.le.xerange(2,3,ie))
     >         iztrue = .true.
            if (.not. if3d) iztrue = .true.

            if (iztrue) then
            if (xloc.ge.xerange(1,1,ie).and.xloc.le.xerange(2,1,ie))then
            if (yloc.ge.xerange(1,2,ie).and.yloc.le.xerange(2,2,ie))then
                ipart(je0 ,ip) = er_map(1,ie)-1
                ipart(jrc ,ip) = 0
                ipart(jpt ,ip) = er_map(2,ie)
                rpart(jd  ,ip) = 1.0 
                rloc = -1.0 + 2.0*(xloc - xerange(1,1,ie))/
     $                 (xerange(2,1,ie)-xerange(1,1,ie))
                sloc = -1.0 + 2.0*(yloc - xerange(1,2,ie))/
     $                 (xerange(2,2,ie)-xerange(1,2,ie))
                rpart(jr  ,ip) = rloc
                rpart(jr+1,ip) = sloc
                if (if3d) then
                tloc = -1.0 + 2.0*(zloc - xerange(1,3,ie))/
     $                 (xerange(2,3,ie)-xerange(1,3,ie))
                rpart(jr+2,ip) = tloc
                endif
                exit
            endif
            endif
            endif
            if (ie .eq. neltb) then
c              ifnd = ifnd + 1
c              rpart(jtempf,ip) = -100
c              write(6,*) nid,xloc,yloc,zloc
            endif
         enddo
      enddo

c     ifnd = iglsum(ifnd,1)
c     if (nid .eq. 0) write(6,*) 'not found',ifnd

      return
      end
c----------------------------------------------------------------------
c     particle input/output/ restart routines
c----------------------------------------------------------------------
      subroutine lpm_usr_particles_io
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'LPM'

      rdumt = dnekclock()

c     always output fields if 2 or 4 way coupled 
      if (two_way .gt. 1) then
         call output_two_way_io
      endif

c     output diagnostics to logfile
      if (npio_method .lt. 0) then
         call output_particle_timers 
         call output_particle_diagnostics
      endif

      call output_parallel_part
      call output_parallel_restart_part

      pttime(1) = pttime(1) - (dnekclock() - rdumt)

      return
      end
!-----------------------------------------------------------------------
      subroutine output_parallel_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      common /myparth/ i_fp_hndl, i_cr_hndl

      character*20 vtufile
      character*14 vtufile1
      character*50 dumstr

      integer icalld
      save    icalld
      data    icalld  /0/

      logical partl         

      integer vtu,vtu1,prevs(2,np)
      integer*4 iint
      real stride_len

      icalld = icalld+1

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,np
         prevs(1,i) = i-1
         prevs(2,i) = n
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call fgslib_crystal_ituple_transfer(i_cr_hndl,prevs,
     >                 ndum,np,np,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,prevs,
     >                 ndum,np,nglob,nkey)

      stride_len = 0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! -----------------------------------------------------
! SEND PARTICLES TO A SMALLER NUMBER OF RANKS TO OUTPUT
! -----------------------------------------------------
      do i = 1,n
         idum = int((stride_len + i)/llpart)
         ipart(jps,i) = idum
      enddo
      nl = 0
      call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     >                  , ipart,ni,partl,nl,rpart,nr,jps)

! -------------------
! GET TOTAL PARTICLES
! -------------------
      nptot = iglsum(n,1)
      npmax_set = int(nptot/llpart) + 1 ! use as few procs as possible
      npmax = np
      if (npmax .gt. npmax_set) npmax = npmax_set

! -----------------------
! WRITE PARALLEL VTU FILE
! -----------------------
      if (nid.eq.0) then
         vtu1 = 107
         write(vtufile1,'(A4,I5.5,A5)')'part',icalld,'.pvtu'
         open(unit=vtu1,file=vtufile1,status='replace',
     >             access='stream',form="unformatted")

         write(vtu1) '<VTKFile '
         write(vtu1) 'type="PUnstructuredGrid" '
         write(vtu1) 'version="1.0" '
         if (lpm_endian .eq. 0) then
            write(vtu1) 'byte_order="LittleEndian"> '
         elseif (lpm_endian .eq. 1) then
            write(vtu1) 'byte_order="BigEndian"> '
         endif
         
         write(vtu1) '<PUnstructuredGrid GhostLevel="0">'

         write(vtu1) '<PPoints> '
         call vtu_write_dataarray(vtu1,"Position    ",3,0,0)
         write(vtu1) '</PPoints> '
         write(vtu1) '<PPointData> '
         call vtu_write_dataarray(vtu1,"VelocityP   ",3,0,0)
         call vtu_write_dataarray(vtu1,"VelocityF   ",3,0,0)
         call vtu_write_dataarray(vtu1,"ForceUSR    ",3,0,0)
         if (part_force(1) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceQS     ",3,0,0)
         if (part_force(2) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceUN     ",3,0,0)
         if (part_force(3) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceIU     ",3,0,0)
         if (two_way .gt. 2)
     >      call vtu_write_dataarray(vtu1,"ForceC      ",3,0,0)
         if (part_force(1) .ne. 0 .or. part_force(3) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceWork   ",1,0,0)
         if (part_force(4) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"HeatQS      ",1,0,0)
         if (part_force(5) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"HeatUN      ",1,0,0)
         if (part_force(4) .ne. 0 .and. part_force(5) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"TemperatureP",1,0,0)
         call vtu_write_dataarray(vtu1,"TemperatureF",1,0,0)
         call vtu_write_dataarray(vtu1,"DensityF    ",1,0,0)
         if (two_way .ge. 2)
     >      call vtu_write_dataarray(vtu1,"ParticleVF  ",1,0,0)
         call vtu_write_dataarray(vtu1,"RadiusCG    ",1,0,0)
         call vtu_write_dataarray(vtu1,"Diameter    ",1,0,0)
         call vtu_write_dataarray(vtu1,"ID_1        ",1,0,0)
         call vtu_write_dataarray(vtu1,"ID_2        ",1,0,0)
         call vtu_write_dataarray(vtu1,"ID_3        ",1,0,0)
         write(vtu1) '</PPointData> '

         write(vtu1) '<PCells> '
         write(vtu1) '<PDataArray '
         write(vtu1) 'type="Int32" '
         write(vtu1) 'Name="connectivity" '
         write(vtu1) 'format="ascii"> '
         write(vtu1) '</PDataArray> '
         write(vtu1) '<PDataArray '
         write(vtu1) 'type="Int32" '
         write(vtu1) 'Name="offsets" '
         write(vtu1) 'format="ascii"> '
         write(vtu1) '</PDataArray> '
         write(vtu1) '<PDataArray '
         write(vtu1) 'type="Int32" '
         write(vtu1) 'Name="types" '
         write(vtu1) 'format="ascii"> '
         write(vtu1) '</PDataArray> '
         write(vtu1) '</PCells> '
         
         do j=0,npmax-1
            write(vtufile,'(A4,I5.5,A2,I5.5,A4)')
     >              'part',icalld,'_p',j,'.vtu'
            write(vtu1) '<Piece Source="',trim(vtufile),'"/>'
         enddo

         write(vtu1) '</PUnstructuredGrid> '
         write(vtu1) '</VTKFile> '

         close(vtu1)

      endif


! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      if (n .gt. 0) then
      write(vtufile,'(A4,I5.5,A2,I5.5,A4)')'part',icalld,'_p',nid,'.vtu'


      vtu=867+nid


      open(unit=vtu,file=vtufile,status='replace',
     >             access='stream',form="unformatted")

! ------------
! FRONT MATTER
! ------------
      write(vtu) '<VTKFile '
      write(vtu) 'type="UnstructuredGrid" '
      write(vtu) 'version="1.0" '
      if (lpm_endian .eq. 0) then
         write(vtu) 'byte_order="LittleEndian"> '
      elseif (lpm_endian .eq. 1) then
         write(vtu) 'byte_order="BigEndian"> '
      endif

      write(vtu) '<UnstructuredGrid> '

      write(vtu) '<FieldData> ' 
      write(vtu) '<DataArray '  ! time
      write(vtu) 'type="Float32" '
      write(vtu) 'Name="TIME" '
      write(vtu) 'NumberOfTuples="1" '
      write(vtu) 'format="ascii"> '
#ifdef CMTNEK
      write(dumstr,'(ES30.16)') time_cmt
#else
      write(dumstr,'(ES30.16)') time
#endif
      write(vtu) trim(dumstr), ' '
      write(vtu) '</DataArray> '

      write(vtu) '<DataArray '  ! cycle
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="CYCLE" '
      write(vtu) 'NumberOfTuples="1" '
      write(vtu) 'format="ascii"> '
      write(dumstr,'(I0)') istep
      write(vtu) trim(dumstr), ' '
      write(vtu) '</DataArray> '
      write(vtu) '</FieldData>'

      write(vtu) '<Piece '
      write(dumstr,'(A16,I0,A2)') 'NumberOfPoints="',n,'" '
      write(vtu) trim(dumstr), ' '
      write(vtu) 'NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      write(vtu) '<Points> '
      iint = 0
      call vtu_write_dataarray(vtu,"Position    ",3,iint,1)
      write(vtu) '</Points> '

! ----
! DATA 
! ----
      idum = 0
      write(vtu) '<PointData> '
      iint = iint + 3*wdsize*n + isize
      call vtu_write_dataarray(vtu,"VelocityP   ",3,iint,1)
      iint = iint + 3*wdsize*n + isize
      call vtu_write_dataarray(vtu,"VelocityF   ",3,iint,1)
      iint = iint + 3*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ForceUSR    ",3,iint,1)
      if (part_force(1) .ne. 0) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceQS     ",3,iint,1)
      endif
      if (part_force(2) .ne. 0) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceUN     ",3,iint,1)
      endif
      if (part_force(3) .ne. 0) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceIU     ",3,iint,1)
      endif
      if (two_way .gt. 2) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceC      ",3,iint,1)
      endif
      idum = 3
      if (part_force(1) .ne. 0 .or. part_force(3) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceWork   ",1,iint,1)
         idum = 1
      endif
      if (part_force(4) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"HeatQS      ",1,iint,1)
         idum = 1
      endif
      if (part_force(5) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"HeatUN      ",1,iint,1)
         idum = 1
      endif
      if (part_force(4) .ne. 0 .and. part_force(5) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"TemperatureP",1,iint,1)
         idum = 1
      endif
      iint = iint + idum*wdsize*n + isize
      call vtu_write_dataarray(vtu,"TemperatureF",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"DensityF    ",1,iint,1)
      if (two_way .ge. 2) then
         iint = iint + 1*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ParticleVF  ",1,iint,1)
      endif

      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"RadiusCG    ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"Diameter    ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ID_1        ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ID_2        ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ID_3        ",1,iint,1)
      write(vtu) '</PointData> '

! ----------
! END MATTER
! ----------
      write(vtu) '<Cells> '
      write(vtu) '<DataArray '
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="connectivity" '
      write(vtu) 'format="ascii"> '
      write(vtu) '</DataArray> '
      write(vtu) '<DataArray '
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="offsets" '
      write(vtu) 'format="ascii"> '
      write(vtu) '</DataArray> '
      write(vtu) '<DataArray '
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="types" '
      write(vtu) 'format="ascii"> '
      write(vtu) '</DataArray> '
      write(vtu) '</Cells> '
      write(vtu) '</Piece> '
      write(vtu) '</UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu) '<AppendedData encoding="raw"> _'

      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jx,i)
         write(vtu) rpart(jy,i)
         if (.not. if3d) then
            write(vtu) 0.0
         else
            write(vtu) rpart(jz,i)
         endif
      enddo
      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jv0,i)
         write(vtu) rpart(jv0+1,i)
         write(vtu) rpart(jv0+2,i)
      enddo
      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(ju0,i)
         write(vtu) rpart(ju0+1,i)
         write(vtu) rpart(ju0+2,i)
      enddo
      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jfusr,i)
         write(vtu) rpart(jfusr+1,i)
         write(vtu) rpart(jfusr+2,i)
      enddo
      if (part_force(1) .ne. 0) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfqs,i)
            write(vtu) rpart(jfqs+1,i)
            write(vtu) rpart(jfqs+2,i)
         enddo
      endif
      if (part_force(2) .ne. 0) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfun,i)
            write(vtu) rpart(jfun+1,i)
            write(vtu) rpart(jfun+2,i)
         enddo
      endif
      if (part_force(3) .ne. 0) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfiu,i)
            write(vtu) rpart(jfiu+1,i)
            write(vtu) rpart(jfiu+2,i)
         enddo
      endif
      if (two_way .gt. 2) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfcol,i)
            write(vtu) rpart(jfcol+1,i)
            write(vtu) rpart(jfcol+2,i)
         enddo
      endif
      if (part_force(1) .ne. 0 .or. part_force(3) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jg0,i)
         enddo
      endif
      if (part_force(4) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jqqs,i)
         enddo
      endif
      if (part_force(5) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jquu,i)
         enddo
      endif
      if (part_force(4) .ne. 0 .and. part_force(5) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jtemp,i)
         enddo
      endif
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jtempf,i)
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jrho,i)
      enddo
      if (two_way .ge. 2) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jvol1,i)
         enddo
      endif
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jrpe,i)
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jdp,i)
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         rdum = real(int(ipart(jpid1,i)))
         write(vtu) rdum
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         rdum = real(int(ipart(jpid2,i)))
         write(vtu) rdum
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         rdum = real(int(ipart(jpid3,i)))
         write(vtu) rdum
      enddo

      write(vtu) ' </AppendedData> '
      write(vtu) '</VTKFile> '

      close(vtu)
      endif

! ---------------------------------------------------
! SEND PARTICLES BACK TO ORIGINAL PROCESSOR WHEN DONE
! ---------------------------------------------------
      call move_particles_inproc

      return
      end
c----------------------------------------------------------------------
      subroutine vtu_write_dataarray(vtu,dataname,ncomp,idist,ip)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      integer vtu,ncomp
      integer*4 idist
      character*12 dataname
      character*50 dumstr

      if (ip .eq. 0) then
         write(vtu) '<PDataArray '
         write(vtu) 'type="Float64" '
         write(dumstr,*)'Name="',trim(dataname),'" '
         write(vtu) trim(dumstr), ' '
         write(dumstr,'(A20,I0,A2)')'NumberOfComponents="',ncomp,'" '
         write(vtu) trim(dumstr), ' '
         write(vtu) 'format="append"/> '
      else
         write(vtu) '<DataArray '
         write(vtu) 'type="Float64" '
         write(dumstr,*)'Name="',trim(dataname),'" '
         write(vtu) trim(dumstr), ' '
         write(dumstr,'(A20,I0,A2)')'NumberOfComponents="',ncomp,'" '
         write(vtu) trim(dumstr), ' '
         write(vtu) 'format="append" '
         write(dumstr,'(A8,I0,A3)') 'offset="',idist,'"/>'
         write(vtu) trim(dumstr), ' '
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_parallel_restart_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /0/

      character*10 filename

      logical partl         ! This is a dummy placeholder, used in cr()

      integer vtu,vtu1,prevs(3,mp)
      integer*4 iint
      integer stride_len

      integer*8 disp
      integer*4 nptot
      integer pth
      real*4 rnptot, rnptot_new

! ----------------------------------------
! Setup file names to write to mpi
! ----------------------------------------
      icalld = icalld+1
      write(filename,'(A5,I5.5)') 'rpart', icalld
      pth = 174

      if (nid.eq.0) then
         open(unit=pth,iostat=istat,file=filename,status='old')
         if (stat .eq. 0) close(pth, status='delete')
      endif
      nptot = iglsum(n,1)
      rnptot = real(nptot)

      disp   = 0
      icount = 0
      if (nid .eq. 0) icount = 1
      iorank = -1
      call byte_open_mpi(filename,pth,.false.,ierr)
      call byte_set_view(disp,pth)
      call byte_write_mpi(rnptot,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,mp
         prevs(1,i) = i-1
         prevs(2,i) = nid
         prevs(3,i) = n
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 2 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum  = 3
      call fgslib_crystal_ituple_transfer(i_cr_hndl,prevs,
     >                 ndum,mp,mp,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,prevs,
     >                 ndum,mp,nglob,nkey)

      stride_len = 0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(3,i)
      enddo
      endif

! -----------------------------------------------------
! SEND PARTICLES TO A SMALLER NUMBER OF RANKS TO OUTPUT
! -----------------------------------------------------
      do i = 0,n-1
         idum = int((stride_len + i)/llpart)
         ipart(jps,i+1) = idum
      enddo
      nl = 0
      call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     >                  , ipart,ni,partl,nl,rpart,nr,jps)

! ----------------------------------------
! Calculate how many processors to be used
! ----------------------------------------
      npmax = min(nptot/llpart+1,mp)
      stride_len = 0
      if (nid .le. npmax-1 .and. nid. ne. 0) stride_len = nid*llpart

! ----------------------
! COPY DATA FOR EASY I/O
! ----------------------
      do i = 1,n
         rfpts(1,i) = rpart(jx,i)
         rfpts(2,i) = rpart(jy,i)
         rfpts(3,i) = rpart(jz,i)
         rfpts(4,i) = rpart(jx1+0,i)
         rfpts(5,i) = rpart(jx1+1,i)
         rfpts(6,i) = rpart(jx1+2,i)
         rfpts(7,i) = rpart(jx2+0,i)
         rfpts(8,i) = rpart(jx2+1,i)
         rfpts(9,i) = rpart(jx2+2,i)
         rfpts(10,i) = rpart(jx3+0,i)
         rfpts(11,i) = rpart(jx3+1,i)
         rfpts(12,i) = rpart(jx3+2,i)
      
         rfpts(13,i) = rpart(jv0,i)
         rfpts(14,i) = rpart(jv0+1,i)
         rfpts(15,i) = rpart(jv0+2,i)
         rfpts(16,i) = rpart(jv1+0,i)
         rfpts(17,i) = rpart(jv1+1,i)
         rfpts(18,i) = rpart(jv1+2,i)
         rfpts(19,i) = rpart(jv2+0,i)
         rfpts(20,i) = rpart(jv2+1,i)
         rfpts(21,i) = rpart(jv2+2,i)
         rfpts(22,i) = rpart(jv3+0,i)
         rfpts(23,i) = rpart(jv3+1,i)
         rfpts(24,i) = rpart(jv3+2,i)
      
         rfpts(25,i) = rpart(ju0,i)
         rfpts(26,i) = rpart(ju0+1,i)
         rfpts(27,i) = rpart(ju0+2,i)
         rfpts(28,i) = rpart(ju1+0,i)
         rfpts(29,i) = rpart(ju1+1,i)
         rfpts(30,i) = rpart(ju1+2,i)
         rfpts(31,i) = rpart(ju2+0,i)
         rfpts(32,i) = rpart(ju2+1,i)
         rfpts(33,i) = rpart(ju2+2,i)
         rfpts(34,i) = rpart(ju3+0,i)
         rfpts(35,i) = rpart(ju3+1,i)
         rfpts(36,i) = rpart(ju3+2,i)
      
         rfpts(37,i) = rpart(jdp,i)
         rfpts(38,i) = rpart(jspl,i)
         rfpts(39,i) = rpart(jtemp,i)
         rfpts(40,i) = real(ipart(jpid1,i))
         rfpts(41,i) = real(ipart(jpid2,i))
         rfpts(42,i) = real(ipart(jpid3,i))
      enddo

! -------------------------------
! WRITE PARTICLES TO RESTART FILE
! -------------------------------
      disp   = isize + stride_len*lrf*isize ! is there real*4 var?
      icount = n*lrf
      iorank = -1
      call byte_open_mpi(filename,pth,.false.,ierr)
      call byte_set_view(disp,pth)
      call byte_write_mpi(rfpts,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

! ----------------------------------------
! Move particles back to where they belong
! ----------------------------------------
      call move_particles_inproc


      return
      end
c----------------------------------------------------------------------
      subroutine read_parallel_restart_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*10 filename

      logical partl         ! This is a dummy placeholder, used in cr()

      integer vtu,vtu1,prevs(2,mp)
      integer*4 iint
      integer stride_len

      integer*8 disp
      integer*4 nptot
      integer pth, nread
      real*4 rnptot

      rpi    = 4.0d+0*atan(1.d+0) ! pi

! ------------------------------
! LOAD TOTAL NUMBER OF PARTICLES
! ------------------------------
      write(filename,'(A5,I5.5)') 'rpart', ipart_restartr

      pth = 185

      disp = 0
      icount = 0
      if (nid .eq. 0) icount = 1
      iorank = -1
      call byte_open_mpi(filename,pth,.true.,ierr)
      call byte_set_view(disp,pth)
      call byte_read_mpi(rnptot,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

      nptot = int(rnptot)
      call bcast(nptot, isize)

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      npmax = min(nptot/llpart+1,mp)
      stride_len = 0
      if (nid .le. npmax-1 .and. nid. ne. 0) stride_len = nid*llpart

      nread = llpart
      if (nid .gt. npmax-1) nread = 0

      ndiff = nptot - (npmax-1)*llpart
      if (nid .eq. npmax-1) nread = ndiff

! -------------------------
! Parallel MPI file read in
! -------------------------
      disp   = isize + stride_len*lrf*isize ! is there real*4 var?
      icount = nread*lrf
      iorank = -1
      call byte_open_mpi(filename,pth,.true.,ierr)
      call byte_set_view(disp,pth)
      call byte_read_mpi(rfpts,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

! ----------------------------------------
! Assign values read in to rpart and ipart
! ----------------------------------------
      i = n ! if there are previous particles
      do ii = 1,nread
         i = n + ii
         rpart(jx,i)    =  rfpts(1,ii) 
         rpart(jy,i)    =  rfpts(2,ii)
         rpart(jz,i)    =  rfpts(3,ii)  
         rpart(jx1+0,i) =  rfpts(4,ii)  
         rpart(jx1+1,i) =  rfpts(5,ii)
         rpart(jx1+2,i) =  rfpts(6,ii)  
         rpart(jx2+0,i) =  rfpts(7,ii)  
         rpart(jx2+1,i) =  rfpts(8,ii)
         rpart(jx2+2,i) =  rfpts(9,ii)  
         rpart(jx3+0,i) =  rfpts(10,ii)
         rpart(jx3+1,i) =  rfpts(11,ii)
         rpart(jx3+2,i) =  rfpts(12,ii)
                                        
         rpart(jv0,i)   =  rfpts(13,ii)
         rpart(jv0+1,i) =  rfpts(14,ii)
         rpart(jv0+2,i) =  rfpts(15,ii)
         rpart(jv1+0,i) =  rfpts(16,ii)
         rpart(jv1+1,i) =  rfpts(17,ii)
         rpart(jv1+2,i) =  rfpts(18,ii)
         rpart(jv2+0,i) =  rfpts(19,ii)
         rpart(jv2+1,i) =  rfpts(20,ii)
         rpart(jv2+2,i) =  rfpts(21,ii)
         rpart(jv3+0,i) =  rfpts(22,ii)
         rpart(jv3+1,i) =  rfpts(23,ii)
         rpart(jv3+2,i) =  rfpts(24,ii)
                                        
         rpart(ju0,i)   =  rfpts(25,ii)
         rpart(ju0+1,i) =  rfpts(26,ii)
         rpart(ju0+2,i) =  rfpts(27,ii)
         rpart(ju1+0,i) =  rfpts(28,ii)
         rpart(ju1+1,i) =  rfpts(29,ii)
         rpart(ju1+2,i) =  rfpts(30,ii)
         rpart(ju2+0,i) =  rfpts(31,ii)
         rpart(ju2+1,i) =  rfpts(32,ii)
         rpart(ju2+2,i) =  rfpts(33,ii)
         rpart(ju3+0,i) =  rfpts(34,ii)
         rpart(ju3+1,i) =  rfpts(35,ii)
         rpart(ju3+2,i) =  rfpts(36,ii)
                                        
         rpart(jdp,i)   =  rfpts(37,ii)
         rpart(jspl,i)  =  rfpts(38,ii)
         rpart(jtemp,i) =  rfpts(39,ii)
         ipart(jpid1,i) =  nint(rfpts(40,ii))
         ipart(jpid2,i) =  nint(rfpts(41,ii))
         ipart(jpid3,i) =  nint(rfpts(42,ii))
      
         ! extra stuff
         rpart(jtaup,i) = rpart(jdp,i)**2*rho_p/18.0d+0/mu_0
         rpart(jrhop,i) = rho_p      ! material density of particle
         rpart(jvol,i) = rpart(jspl,i)*rpi*rpart(jdp,i)**3/6.d+0! particle volume
         rpart(jgam,i)  = 1.          ! initial integration correction
         rpart(jrpe,i) = rpart(jspl,i)**(1./3.)*rpart(jdp,i)/2.

      enddo
      n = i

      return
      end
c----------------------------------------------------------------------
      subroutine read_ascii_restart_part
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real unif_random
      external unif_random

      character*10 filename

      ntotal = iglsum(n,1)
      if (ntotal .gt. nw) return

      write(filename,'(A5,I5.5)') 'rpart', abs(ipart_restartr)
      if (nid .eq. 0) then
         iread = 753
         open(unit=iread,file=filename,action='read',status='old')
         read(iread,*) ndum
         ndum1 = abs(ndum)
         do i=1,ndum1
            n = n + 1
            if (ndum .gt. 0) then
               read(iread,*) rpart(jx,n),rpart(jy,n),rpart(jz,n)
     >                   ,rpart(jv0+0,n), rpart(jv0+1,n), rpart(jv0+2,n)

               if (dp_std .gt. 0) then
                  rpart(jdp,n) = unif_random_norm(dp(1),dp_std)
               else
                  rpart(jdp,n) = unif_random(dp(1),dp(2))
               endif
               rpart(jtaup,n) = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0
               rpart(jrhop,n) = rho_p
               rpart(jvol,n)  = pi*rpart(jdp,n)**3/6.
               rpart(jspl,n)  = rspl
               rpart(jrpe,n)  = rpart(jspl,n)**(1./3.)*rpart(jdp,n)/2.
               rpart(jvol,n)  = rpart(jspl,n)*rpart(jvol,n)

               rpart(jtemp,n)  = tp_0
               rpart(jtempf,n) = tp_0
               rpart(jrho,n)   = param(1)

               ipart(jpid1,n) = nid
               ipart(jpid2,n) = i
               ipart(jpid3,n) = istep

            endif
            if (ndum .lt. 0) then
               read(iread,*) rpart(jx,n),rpart(jy,n),rpart(jz,n)
     >                   ,rpart(jv0+0,n), rpart(jv0+1,n), rpart(jv0+2,n)
     >                   ,rpart(jdp,n), rpart(jspl,n), rpart(jtemp,n)
     >                   ,ipart(jpid1,n), ipart(jpid2,n), ipart(jpid3,n)

               rpart(jtaup,n) = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0
               rpart(jrhop,n) = rho_p
               rpart(jvol,n)  = rpart(jspl,n)*pi*rpart(jdp,n)**3/6.d+0
               rpart(jgam,n)  = 1.
               rpart(jrpe,n)  = rpart(jspl,n)**(1./3.)*rpart(jdp,n)/2.
            endif

            rfluc = 1E-2
            rfluc = rfluc*rpart(jrpe,n)
            rpart(jx,n) = rpart(jx,n) + unif_random(-rfluc,rfluc)
            rpart(jy,n) = rpart(jy,n) + unif_random(-rfluc,rfluc)
            rpart(jz,n) = rpart(jz,n) + unif_random(-rfluc,rfluc)

            do j=0,2
               rpart(jx1+j,n) = rpart(jx+j,n)
               rpart(jx2+j,n) = rpart(jx+j,n)
               rpart(jx3+j,n) = rpart(jx+j,n)
               rpart(jv1+j,n) = rpart(jv0+j,n)
               rpart(jv2+j,n) = rpart(jv0+j,n)
               rpart(jv3+j,n) = rpart(jv0+j,n)
            enddo

         enddo
         close(iread)
      endif

      call update_particle_location   
      call move_particles_inproc

      return
      end
c----------------------------------------------------------------------
      subroutine output_two_way_io
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'LPM'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real rfpfluid(3),rfpfluidl(3),msum,msum_tot(3,2)

      itmp = 1
      call outpost2(ptw(1,1,1,1,1),         ! fhyd_x
     >              ptw(1,1,1,1,2),         ! fhyd_y
     >              ptw(1,1,1,1,3),         ! fhyd_z
     >              ptw(1,1,1,1,4),         ! phi_p (but not if lx1!=lx2
     >              ptw(1,1,1,1,4),         ! phi_p
     >              itmp          ,        
     >              'ptw')

c     eulerian integrations -----------------------------------------
c     fluid momentum 
      msum_tot(1,1) = glsc3(bm1,vtrans,vx,nx1*ny1*nz1*nelv)
      msum_tot(2,1) = glsc3(bm1,vtrans,vy,nx1*ny1*nz1*nelv)
      msum_tot(3,1) = glsc3(bm1,vtrans,vz,nx1*ny1*nz1*nelv)
c     particle volume fraction
      vf_part_e     = glsc2(bm1,ptw(1,1,1,1,4),nx1*ny1*nz1*nelt)
c     particle forces on fluid
      rfpfluid(1)   = glsc2(bm1,ptw(1,1,1,1,1),nx1*ny1*nz1*nelt)
      rfpfluid(2)   = glsc2(bm1,ptw(1,1,1,1,2),nx1*ny1*nz1*nelt)
      rfpfluid(3)   = glsc2(bm1,ptw(1,1,1,1,3),nx1*ny1*nz1*nelt)

      rdum = xdrange(2,3) - xdrange(1,3)
      if (.not.if3d) vf_part_e   = vf_part_e*rdum   ! Here:
      if (.not.if3d) rfpfluid(1) = rfpfluid(1)*rdum ! for 2d, assume
      if (.not.if3d) rfpfluid(2) = rfpfluid(2)*rdum ! z thicknes 
      if (.not.if3d) rfpfluid(3) = rfpfluid(3)*rdum 


c     lagrangian integrations ---------------------------------------
c     particle momentum
      do ieq=0,2
         msum = 0.0
         rsum = 0.0
         do i=1,n
           msum = msum + 
     >       rpart(jv0+ieq,i)*rpart(jrhop,i)*rpart(jvol,i)
           rsum = rsum + rpart(jf0+ieq,i)
         enddo
         msum_tot(ieq+1,2) = glsum(msum,1)
         rfpfluidl(1+ieq)  = glsum(rsum,1)
      enddo
c     particle volume fraction
      msum = 0.0
      do i=1,n
         msum = msum + rpart(jvol,i)
      enddo
      vf_part_l = glsum(msum,1)

      vf_rel_error = abs(vf_part_l - vf_part_e)/vf_part_l*100.0

c     print to files ------------------------------------------------
c     print properties to logfile
      if (nid.eq.0) write(6,500) "E Fluid Momentum :              ", 
     >                  istep, msum_tot(1,1),msum_tot(2,1),msum_tot(3,1)
      if (nid.eq.0) write(6,500) "E Particle forces:              ", 
     >                  istep, rfpfluid(1),rfpfluid(2),rfpfluid(3)
      if (nid.eq.0) write(6,500) "E Particle Volume:              ", 
     >                  istep, vf_part_e
      if (nid.eq.0) write(6,500) "L Particle Momentum :           ", 
     >                  istep, msum_tot(1,2),msum_tot(2,2),msum_tot(3,2)
      if (nid.eq.0) write(6,500) "L Particle forces:              ", 
     >                  istep, rfpfluidl(1),rfpfluidl(2),rfpfluidl(3)
      if (nid.eq.0) write(6,500) "L Particle Volume:              ", 
     >                  istep, vf_part_l
      if (nid.eq.0) write(6,500) "VF Relative Error %              ", 
     >                  istep, vf_rel_error

  500 FORMAT(A30,I20,3ES20.10)

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_input_defaults
      include 'SIZE'
      include 'LPM'

      ! set some defaults
      nw = 0
      rxbo(1,1) = -1E8
      rxbo(2,1) = -1E8
      rxbo(1,2) = -1E8
      rxbo(2,2) = -1E8
      rxbo(1,3) = -1E8
      rxbo(2,3) = -1E8
      rxco(1) = -1E8
      rxco(2) = -1E8
      rxco(3) = -1E8
      rxco(4) = -1E8
      rxco(5) = -1E8
      rxco(6) = -1E8
      rxco(7) = -1E8
      rxco(8) = -1E8
      rxco(9) = -1E8
      dp(1) = 0.
      dp(2) = 0.
      dp_std = -1.
      tp_0  = 273.
      rho_p = 2500.
      cp_p  = 840.
      do i=1,5
         part_force(i) = 0
      enddo
      time_integ  = 1
      two_way     = 1
      red_interp  = 1
      npio_method = 1
      inject_rate = 0
      time_delay  = 0
      nrandseed   = 1
      lpm_endian  = 0
      npro_method = 1
      rspl        = 1.
      dfilt       = 2.
      ralphdecay  = 1E-3
      do i=1,6
         bc_part(i) = 1
      enddo
      ipart_restartr = 0
      ksp            = 10.
      e_rest         = 0.9
      ksp_wall       = ksp
      e_rest_wall    = e_rest

      np_walls = 0
      nc_walls = 0

      nbin_save = 0

      ifrectp = 0

c     zlpart = 0.0
c     zrpart = 1.0
c     if3d = if3d


      return
      end
c----------------------------------------------------------------------
c     subroutine read_particle_input_par
c     include 'SIZE'
c     include 'INPUT'
c     include 'TSTEP'
c     include 'PARALLEL'
c     include 'LPM'
c
c     only used for load balanced code since there is no par file in 
c     old version of the code. To convert from standard nek5000 input:
c         periodic* = only needs to be included for it to be periodic
c         so yes or no doesn't matter
c
c         remove commas seperating parameter inputs
c
c         must be lowercase key words with one space then an equal sign
c         followed by another space

c     character*72 dum_str

c     character*200 buffer, label
c     integer pos, fh, ios, line, dum
c     parameter(fh = 15)

c     call lpm_input_defaults

c     ios  = 0
c     line = 0

c     open(fh, file='particles.par')
c
c     do while (ios == 0)
c        read(fh, '(A)', iostat=ios) buffer
c        if (ios == 0) then
c           line = line + 1
c     
c           ! Find the first instance of whitespace.  Split label and data.
c           pos = scan(buffer, '=')
c           label = buffer(1:pos)
c           buffer = buffer(pos+1:)
c     
c           select case (label)
c           case ('npart =')
c              read(buffer, *, iostat=ios) nw
c              if(nid.eq.0)write(6,*) 'Read npart: ', nw
c           case ('distributebox =')
c              read(buffer, *, iostat=ios) rxbo(1,1),rxbo(2,1),
c    >                                     rxbo(1,2),rxbo(2,2),
c    >                                     rxbo(1,3),rxbo(2,3)
c              if(nid.eq.0)write(6,*) 'Read distributebox '
c           case ('distributecylinder =')
c              read(buffer, *, iostat=ios) rxco(1),rxco(2),rxco(3),
c    >                                     rxco(4),rxco(5),rxco(6), 
c    >                                     rxco(7),rxco(8),rxco(9)
c              if(nid.eq.0)write(6,*) 'Read distributecylinder '
c           case ('distributesphere =')
c              read(buffer, *, iostat=ios) rxco(1),rxco(2),rxco(3),
c    >                                     rxco(4),rxco(5)
c              if(nid.eq.0)write(6,*) 'Read distributesphere '
c           case ('diameter =')
c              read(buffer, *, iostat=ios) dp(1)
c              if(nid.eq.0)write(6,*) 'Read diameter: ', dp(1)
c              dp(2) = dp(1)
c           case ('diameteruniform =')
c              read(buffer, *, iostat=ios) dp(1), dp(2)
c              if(nid.eq.0)write(6,*) 'Read diameteruniform: ', 
c    >                         dp(1), dp(2)
c           case ('diametergaussian =')
c              read(buffer, *, iostat=ios) dp(1), dp_std
c              if(nid.eq.0)write(6,*) 'Read diametergaussian: ', 
c    >                         dp(1), dp_std
c              dp(2) = dp(1)
c           case ('temperature =')
c              read(buffer, *, iostat=ios) tp_0
c              if(nid.eq.0)write(6,*) 'Read temperature: ', tp_0
c           case ('density =')
c              read(buffer, *, iostat=ios) rho_p
c              if(nid.eq.0)write(6,*) 'Read density: ', rho_p
c           case ('specificheat =')
c              read(buffer, *, iostat=ios) cp_p
c              if(nid.eq.0)write(6,*) 'Read specificheat: ', cp_p
c           case ('forceqs =')
c              read(buffer, *, iostat=ios) part_force(1)
c              if(nid.eq.0)write(6,*) 'Read forceqs: ', part_force(1)
c           case ('forceun =')
c              read(buffer, *, iostat=ios) part_force(2)
c              if(nid.eq.0)write(6,*) 'Read forceun: ', part_force(2)
c           case ('forceiu =')
c              read(buffer, *, iostat=ios) part_force(3)
c              if(nid.eq.0)write(6,*) 'Read forceiu: ', part_force(3)
c           case ('heatqs =')
c              read(buffer, *, iostat=ios) part_force(4)
c              if(nid.eq.0)write(6,*) 'Read heatqs: ', part_force(4)
c           case ('heatun =')
c              read(buffer, *, iostat=ios) part_force(5)
c              if(nid.eq.0)write(6,*) 'Read heatun: ', part_force(5)
c           case ('timestepper =') 
c              read(buffer, *, iostat=ios) time_integ
c              if(nid.eq.0)write(6,*) 'Read timestepper: ', time_integ
c           case ('coupling =') 
c              read(buffer, *, iostat=ios) two_way
c              if(nid.eq.0)write(6,*) 'Read coupling: ', two_way
c           case ('interpolation') 
c              read(buffer, *, iostat=ios) red_interp
c              if(nid.eq.0)write(6,*) 'Read interpolation: ', red_interp
c           case ('io =')
c              read(buffer, *, iostat=ios) npio_method
c              if(nid.eq.0)write(6,*) 'Read io: ', npio_method
c           case ('injectionstep =')
c              read(buffer, *, iostat=ios) inject_rate
c              if(nid.eq.0)write(6,*) 'Read injectionstep: ',inject_rate
c           case ('delaystep =')
c              read(buffer, *, iostat=ios) time_delay
c              if(nid.eq.0)write(6,*) 'Read delaystep: ', time_delay
c           case ('seed =')
c              read(buffer, *, iostat=ios) nrandseed
c              if(nid.eq.0)write(6,*) 'Read seed: ', nrandseed
c           case ('projection =')
c              read(buffer, *, iostat=ios) npro_method
c              if(nid.eq.0)write(6,*) 'Read projection: ', npro_method
c           case ('coarsegrain =')
c              read(buffer, *, iostat=ios) rspl
c              if(nid.eq.0)write(6,*) 'Read coarsegrain: ', rspl
c           case ('filter =')
c              read(buffer, *, iostat=ios) dfilt
c              if(nid.eq.0)write(6,*) 'Read filter: ', dfilt
c           case ('alpha =')
c              read(buffer, *, iostat=ios) ralphdecay
c              if(nid.eq.0)write(6,*) 'Read alpha: ', ralphdecay

c           case ('wallp01 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp02 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp03 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp04 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp05 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp06 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp07 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp08 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp09 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallc01 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc02 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc03 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc04 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc05 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc06 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc07 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc08 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc09 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('periodicx =')
c              read(buffer, *, iostat=ios)
c              if(nid.eq.0)write(6,*) 'Read periodicx '
c              bc_part(1) = 0
c              bc_part(2) = 0
c           case ('periodicy =')
c              read(buffer, *, iostat=ios)
c              if(nid.eq.0)write(6,*) 'Read periodicy '
c              bc_part(3) = 0
c              bc_part(4) = 0
c           case ('periodicz =')
c              read(buffer, *, iostat=ios)
c              if(nid.eq.0)write(6,*) 'Read periodicz '
c              bc_part(5) = 0
c              bc_part(6) = 0
c           case ('restartstep =')
c              read(buffer, *, iostat=ios) ipart_restartr
c              if(nid.eq.0)write(6,*)'Read restartstep: ',ipart_restartr
c           case ('spring =')
c              read(buffer, *, iostat=ios) ksp
c              if(nid.eq.0)write(6,*) 'Read spring: ', ksp
c           case ('restitution =')
c              read(buffer, *, iostat=ios) e_rest
c              if(nid.eq.0)write(6,*) 'Read restitution: ', e_rest
c           case default
c              if(nid.eq.0)write(6,*) 'Skipping label at line', line
c           end select
c        end if
c     enddo

c     close(fh)


c     return
c     end
c----------------------------------------------------------------------
c     subroutine add_wall_list(ninc,buf,list,ncomp)
c     include 'SIZE'
c     include 'TOTAL'
c     include 'LPM'

c     integer ninc,ncom
c     real list(ncomp,n_walls)
c     character*200 buf

c     ninc = ninc + 1
c     if (ninc .gt. n_walls) then
c        if (nid.eq.0) then
c            write(6,*) 'Increase max number particle wall',ncomp
c            call exitt
c         endif
c     endif
c     if (ncomp .eq. 7) then
c        read(buf, *, iostat=ios) cyl_wall_coords(1,nc_walls)
c    >                           ,cyl_wall_coords(2,nc_walls)
c    >                           ,cyl_wall_coords(3,nc_walls)
c    >                           ,cyl_wall_coords(4,nc_walls)
c    >                           ,cyl_wall_coords(5,nc_walls)
c    >                           ,cyl_wall_coords(6,nc_walls)
c    >                           ,cyl_wall_coords(7,nc_walls)
c        if(nid.eq.0)write(6,*) 'Read wall_cyl number ', ninc
c     elseif(ncomp .eq. 6) then
c        read(buf, *, iostat=ios) plane_wall_coords(1,np_walls)
c    >                           ,plane_wall_coords(2,np_walls)
c    >                           ,plane_wall_coords(3,np_walls)
c    >                           ,plane_wall_coords(4,np_walls)
c    >                           ,plane_wall_coords(5,np_walls)
c    >                           ,plane_wall_coords(6,np_walls)
c         if(nid.eq.0)write(6,*) 'Read wall_plane number ',ninc
c      endif

c     return
c     end
c----------------------------------------------------------------------
      subroutine output_particle_diagnostics
      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      integer icalld
      save icalld
      data icalld /0/

      real rtpart
      save rtpart

      real rdiags_part(3,2), rvels(4)

      if (icalld .eq. 0) then
         rtpart = time_cmt
         if (icmtp .eq. 0) rtpart = time
         icalld = icalld + 1
      endif
      rtdum  = time_cmt - rtpart
      if (icmtp .eq. 0) rtdum = time - rtpart

      do i=1,3
         rdiags_part(i,1) =  1E8
         rdiags_part(i,2) = -1E8
      enddo
      do i=1,4
         rvels(i) = 0.0
      enddo
      
      do i=1,n
         rvx = rpart(jv0,i)
         rvy = rpart(jv0+1,i)
         rvz = rpart(jv0+2,i)
         rvmag = sqrt(rvx**2 + rvy**2 + rvz**2)

         rxx = rpart(jx,i)
         ryy = rpart(jy,i)
         rzz = rpart(jz,i)
         rrp = rpart(jrpe,i)

         ! maxes
         if (rxx + rrp .gt. rdiags_part(1,2)) rdiags_part(1,2) = rxx+rrp
         if (ryy + rrp .gt. rdiags_part(2,2)) rdiags_part(2,2) = ryy+rrp
         if (rzz + rrp .gt. rdiags_part(3,2)) rdiags_part(3,2) = rzz+rrp

         ! mins
         if (rxx - rrp .lt. rdiags_part(1,1)) rdiags_part(1,1) = rxx-rrp
         if (ryy - rrp .lt. rdiags_part(2,1)) rdiags_part(2,1) = ryy-rrp
         if (rzz - rrp .lt. rdiags_part(3,1)) rdiags_part(3,1) = rzz-rrp

         ! velocities
         if ( abs(rvx)   .gt. abs(rvels(1)) )  rvels(1) = rvx
         if ( abs(rvy)   .gt. abs(rvels(2)) )  rvels(2) = rvy
         if ( abs(rvz)   .gt. abs(rvels(3)) )  rvels(3) = rvz
         if ( abs(rvmag) .gt. abs(rvels(4)) )  rvels(4) = rvmag
      enddo

      ! compute globally now
      do i=1,3
         rdum = rdiags_part(i,1)
         rdiags_part(i,1) =  glmin(rdum,1)
         rdum = rdiags_part(i,2)
         rdiags_part(i,2) =  glmax(rdum,1)
      enddo
      do i=1,4
         rdum  = rvels(i)
         rdum1 = glmin(rdum,1)
         rdum  = rvels(i)
         rdum2 = glmax(rdum,1)
         rvels(i) = rdum1
         if (abs(rdum1) .lt. abs(rdum2)) rvels(i) = rdum2
      enddo

      nplow = iglmin(n,1)
      nptot = iglsum(n,1)
      nphigh = iglmax(n,1)

      nplowg = iglmin(nfptsgp,1)
      nptotg = iglsum(nfptsgp,1)
      nphighg = iglmax(nfptsgp,1)

      if (nid .eq. 0) then
      write(6,*)'----- START PARTICLE DIAGNOSTICS: -----'
      write(6,*)'NPART TOTAL      :',nptot,nplow,nphigh
      write(6,*)'NGHOST TOTAL     :',nptotg,nplowg,nphighg
      write(6,*)'XMIN,XMAX        :',rdiags_part(1,1),rdiags_part(1,2)
      write(6,*)'YMIN,YMAX        :',rdiags_part(2,1),rdiags_part(2,2)
      write(6,*)'ZMIN,ZMAX        :',rdiags_part(3,1),rdiags_part(3,2)
      write(6,*)'MAX(VX,VY,VZ,|V|):',rvels(1),rvels(2),rvels(3),rvels(4)
      write(6,*)'PTIME            :',rtdum
      write(6,*)'----- END PARTICLE DIAGNOSTICS: -----'
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_particle_timers
      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      rftime_t = 0.
      rptime_t = 0.

      if(nid.eq.0) then
         write(6,*) 'TIMER H: ', istep
      endif

      do i = 1,iptlen
         rdum  = pttime(i)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/mp
         if(nid.eq.0)  write(6,*) 'TIMER #:',i,rtime

         ! fluid and particle total time: note i == iptlen is f_col
         if (i .eq. 1) rftime_t = rtime
         if ((i .gt. 1) .and. (i.ne.iptlen)) rptime_t = rptime_t +
     >                                                  rtime
      enddo

      
      if (nid.eq.0) then
         write (6,*) 'TOTAL F:', rftime_t
         write (6,*) 'TOTAL P:', rptime_t
      endif


      return
      end
c----------------------------------------------------------------------
c     effeciently move particles between processors routines
c----------------------------------------------------------------------
      subroutine move_particles_inproc
c     Interpolate fluid velocity at current xyz points and move
c     data to the processor that owns the points.
c     Input:    n = number of points on this processor
c     Output:   n = number of points on this processor after the move
c     Code checks for n > llpart and will not move data if there
c     is insufficient room.
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LPM'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

c     common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
c     integer gfirst, inoassignd, resetFindpts, pload

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         ! This is a dummy placeholder, used in cr()

      nl = 0                ! No logicals exchanged

      if (icalld1.eq.0) then
         call set_check_spl_params 

         tolin = 1.e-12
         if (wdsize.eq.4) tolin = 1.e-6
         call intpts_setup  (tolin,i_fp_hndl)
         call fgslib_crystal_setup (i_cr_hndl,nekcomm,np)
         icalld1 = icalld1 + 1
      endif

      call update_particle_location  
      call setup_bins

      ! instead get which bin it is in
      do i=1,n
         ! check if particles are greater or less than binb bounds....
         ! add below....

         ii    = floor((rpart(jx,i)-binb(1))/rdxgp) 
         jj    = floor((rpart(jy,i)-binb(3))/rdygp) 
         kk    = floor((rpart(jz,i)-binb(5))/rdzgp) 
         if (.not. if3d) kk = 0
         if (ii .eq. ndxgp) ii = ndxgp - 1
         if (jj .eq. ndygp) jj = ndygp - 1
         if (kk .eq. ndzgp) kk = ndzgp - 1
         ndum  = ii + ndxgp*jj + ndxgp*ndygp*kk
         ndumn = ndum/nreach
         nrank = modulo(ndumn, np)

         ipart(jicx,i) = ii
         ipart(jicy,i) = jj
         ipart(jicz,i) = kk
         ipart(jicn,i) = ndum

         ipart(jpt,i) = nrank
      enddo

      jps = jpid1-1     ! Pointer to temporary proc id for swapping
      do i=1,n        ! Can't use jpt because it messes up particle info
         ipart(jps,i) = ipart(jpt,i)
      enddo
      call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     $           , ipart,ni,partl,nl,rpart,nr,jps)


c     if (nbin_now .ne. nbin_save) call partition_bins

c     if (nbin_now .lt. np) then
c        np_high = np
c        np_low  = np - nbin_now +1
c        np_low  = max(np_low,nbin_now)

c        if (n .gt. 0) then
c        do i=np_high,np_low,-1
c           if (lpm_binsort(1,i) .eq. nid) then
c              ii = np_high + 1 - i

c              write(6,*) 'yoo',ii

c              nhalf = n/2 + 1
c              nhalf = max(nhalf,0)

c              do ip =1,nhalf
c                 ipart(jpt,ip) = lpm_binsort(1,ii)
c              enddo
c           endif
c        enddo
c        endif

c        jps = jpid1-1     ! Pointer to temporary proc id for swapping
c        do i=1,n        ! Can't use jpt because it messes up particle info
c           ipart(jps,i) = ipart(jpt,i)
c        enddo

c        call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
c    $           , ipart,ni,partl,nl,rpart,nr,jps)

c     endif


      return
      end
c-----------------------------------------------------------------------
      subroutine intpts_setup(tolin,ih)
c
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      tol = tolin
      if (tolin.lt.0) tol = 1e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by
c
      if(nid.eq.0) write(6,*) 'initializing intpts(), tol=', tol
      call fgslib_findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
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
      function unif_random(rxl,rxr)
c
c     must initialize ran2 first
c
      real xl,xr,unif_random

      rdum       = ran2(2)
      rlen       = rxr - rxl
      unif_random= rxl + rdum*rlen

      return
      end
c-----------------------------------------------------------------------
      function unif_random_norm(rmu,rstd)
c
c     must initialize ran2 first
c
      real xl,xr,unif_random_norm,rstd,rxfne(1000),rcdf(1000)

      nxfn  = 1000
      rxlf  = max(0.0,rmu-5.*rstd)
      rxrf  = rmu+5.*rstd
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
      function unif_random_sphere(rxl,rxr)
c
c     must initialize ran2 first
c
      parameter (nxfn = 10000)
      real xl,xr,unif_random_sphere,rxfne(nxfn)

      real    unif_random
      external unif_random

      rxlf  = rxl
      rxrf  = rxr
      rdxf  = (1. - (rxlf/rxrf)**3)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = (rxlf/rxrf)**3 + (i-1.)*rdxf
      enddo

      rdum = unif_random((rxl/rxr)**3,rxr/rxr)

      unif_random_sphere = rxr*(rdum)**(1./3.)

      return
      end
c-----------------------------------------------------------------------
      function unif_random_cyl(rxl,rxr)
c
c     must initialize ran2 first
c
      parameter (nxfn = 10000)
      real xl,xr,unif_random_cyl,rxfne(nxfn)

      real    unif_random
      external unif_random

      rxlf  = rxl
      rxrf  = rxr
      rdxf  = (1. - (rxlf/rxrf)**2)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = (rxlf/rxrf)**2 + (i-1.)*rdxf
      enddo

      rdum = unif_random((rxl/rxr)**2,rxr/rxr)

      unif_random_cyl = rxr*sqrt(rdum)

      return
      end
c----------------------------------------------------------------------
C> Compute coefficients for Runge-Kutta stages \cite{TVDRK}
      subroutine set_tstep_coef_part(dt_in)

      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      real dt_in

      dt_cmt = dt_in

      tcoef(1,1) = 0.0
      tcoef(2,1) = 1.0 
      tcoef(3,1) = dt_cmt
      tcoef(1,2) = 3.0/4.0
      tcoef(2,2) = 1.0/4.0 
      tcoef(3,2) = dt_cmt/4.0 
      tcoef(1,3) = 1.0/3.0
      tcoef(2,3) = 2.0/3.0 
      tcoef(3,3) = dt_cmt*2.0/3.0 

      return
      end
c----------------------------------------------------------------------
      subroutine compute_phig_qtl(rdt_in,div)
c
c     Computes modified divergence constraint for multiphase dense
c     compressible flow
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real div(lx2,ly2,lz2,lelv)

      common /phig_qtl_blk/ phig_last,phig_qtl
      real phig_last(lx1,ly1,lz1,lelt,4),phig_qtl(lx1,ly1,lz1,lelt)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)

      integer icalld
      save    icalld
      data    icalld  /-1/

      icalld = icalld + 1

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         phig_last(ix,iy,iz,ie,1) = 1. - ptw(ix,iy,iz,ie,4)
         phig_qtl(ix,iy,iz,ie) = 0.
      enddo
      enddo
      enddo
      enddo

      if (icalld .eq. 0) then
         do ie=1,nelt
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            phig_last(ix,iy,iz,ie,2) = 1. - ptw(ix,iy,iz,ie,4)
            phig_last(ix,iy,iz,ie,3) = 1. - ptw(ix,iy,iz,ie,4)
            phig_last(ix,iy,iz,ie,4) = 1. - ptw(ix,iy,iz,ie,4)
         enddo
         enddo
         enddo
         enddo

         goto 123
      endif
      
      if (icalld .lt. 5) goto 123

      do ie=1,nelt
         call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                   phig_last(1,1,1,ie,1),lx1,if3d)
         
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            phig_qtl(ix,iy,iz,ie) = phig_last(ix,iy,iz,ie,1) -
     >                              phig_last(ix,iy,iz,ie,2)
            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)/rdt_in

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vx(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dx
     >             (ur(ix,iy,iz)*RXM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SXM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TXM1(ix,iy,iz,ie))

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vy(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dy
     >             (ur(ix,iy,iz)*RYM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SYM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TYM1(ix,iy,iz,ie))

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vz(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dz
     >             (ur(ix,iy,iz)*RZM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SZM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TZM1(ix,iy,iz,ie))


            phig_qtl(ix,iy,iz,ie) = -1.0*phig_qtl(ix,iy,iz,ie)/
     >                              phig_last(ix,iy,iz,ie,1)
         enddo
         enddo
         enddo
      enddo

  123 continue

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         phig_last(ix,iy,iz,ie,4) = phig_last(ix,iy,iz,ie,3)
         phig_last(ix,iy,iz,ie,3) = phig_last(ix,iy,iz,ie,2)
         phig_last(ix,iy,iz,ie,2) = phig_last(ix,iy,iz,ie,1)
      enddo
      enddo
      enddo
      enddo

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         div(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)
      enddo
      enddo
      enddo
      enddo

c     wght = 1.0
c     ncut = 1
c     call filter_s0(div(1,1,1,1),wght,ncut,'div') 

      return
      end
c-----------------------------------------------------------------------
      subroutine pre_sim_collisions
c
c     time stepping routine for pre-simulation collisions/settling
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage
      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      nstage_part = 3
      if (abs(time_integ) .eq. 2) nstage_part = 1

      ! pre simulation iteration for packed bed
      do istep=0,nsteps

         if(modulo(istep,iostep).eq.0) then
            call lpm_usr_particles_io
         endif

         do stage=1,nstage_part

            if (stage .eq. 1) then
               rdt_part = abs(param(12))
               call lpm_set_dt(rdt_part)
               dt_cmt     = rdt_part
               dt         = dt_cmt
               time       = time + dt_cmt
               time_cmt   = time_cmt + dt_cmt
               if (abs(time_integ).eq.1)
     >            call set_tstep_coef_part(rdt_part)
               if (nid.eq. 0) 
     >            write(6,*) 'pre-sim_io time',istep,time,dt_cmt
            endif

            call lpm_usr_particles_solver

         enddo

      enddo

      if (nid.eq.0) write(6,*) 'FINISHED PRE COLLISIONS - EXITING NOW'
      call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_userf(ix,iy,iz,e,ffxp,ffyp,ffzp,qvolp)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'LPM'

      integer ix,iy,iz,e

      real ffxp,ffyp,ffzp,qvolp

      ! particle forcing
      if (two_way .ge.2) then
         if (istep .gt. time_delay) then
            ffxp =  ptw(ix,iy,iz,e,1)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ffyp =  ptw(ix,iy,iz,e,2)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ffzp =  ptw(ix,iy,iz,e,3)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))

            qvolp=0.0
         else
            ffxp = 0.0
            ffyp = 0.0
            ffzp = 0.0
         endif
      else
         ffxp = 0.0
         ffyp = 0.0
         ffzp = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bdf_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)
      
      real s,pmass

      call get_bdf_ext_coefs(beta,alpha,times)

      jx0 = jx
c     move data to previous positions
      do j=0,ndim-1
      do i=1,n
         rpart(ju3+j,i)=rpart(ju2+j,i)
         rpart(ju2+j,i)=rpart(ju1+j,i)
         rpart(ju1+j,i)=rpart(ju0+j,i)
         rpart(jv3+j,i)=rpart(jv2+j,i)
         rpart(jv2+j,i)=rpart(jv1+j,i)
         rpart(jv1+j,i)=rpart(jv0+j,i)
         rpart(jx3+j,i)=rpart(jx2+j,i)
         rpart(jx2+j,i)=rpart(jx1+j,i)
         rpart(jx1+j,i)=rpart(jx0+j,i)
      enddo
      enddo

c     Solve for velocity at time t^n
      do i=1,n
        do j=0,ndim-1
          rhs = rpart(jf0+j,i)
     $        +     beta (1)*rpart(jv1+j,i)
     $        +     beta (2)*rpart(jv2+j,i)
     $        +     beta (3)*rpart(jv3+j,i)
          rpart(jv0+j,i) = rhs / beta(0)
          rhx = beta (1)*rpart(jx1+j,i)
     $        + beta (2)*rpart(jx2+j,i)
     $        + beta (3)*rpart(jx3+j,i) + rpart(jv0+j,i)
          rpart(jx0+j,i) = rhx / beta(0)     ! Implicit solve for x
        enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      FUNCTION MixtPerf_C_GRT_part(G,R,T,icmt)
      REAL G,R,T! INTENT(IN) 
      REAL MixtPerf_C_GRT_part
      integer icmt
      if (icmt .eq. 0) then
         MixtPerf_C_GRT_part = 1.
      else
         MixtPerf_C_GRT_part = SQRT(G*R*T)
      endif

      END
c-----------------------------------------------------------------------
c     Force models below
c-----------------------------------------------------------------------
      subroutine lpm_f_qs
c
c     quasi-steady force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (abs(part_force(1)).eq.1) then

         rdum = lpmvol_p*lpmdens_p/lpmtau_p

         if (part_force(1) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum  = rdum*(1. - 2.*rphip)/(1. - rphip)**3
         endif

         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

      elseif (abs(part_force(1)).eq.2) then
         rdum = lpmvol_p*lpmdens_p/lpmtau_p

         rmacr= 0.6
         rcd_std = 1.+0.15*lpmre_p**(0.687) + 
     >               0.42*(1.+42500./lpmre_p**(1.16))**(-1)
         rcd_mcr = 1.+0.15*lpmre_p**(0.684) + 
     >               0.513*(1. + 483./lpmre_p**(0.669))**(-1)
         rcd1 = rcd_std + (rcd_mcr - rcd_std)*lpmmach_p/rmacr
         
         rdum = rdum*rcd1

         if (part_force(1) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum  = rdum*(1. - 2.*rphip)/(1. - rphip)**3
         endif

         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

      elseif (abs(part_force(1)).eq.3) then
         rdum = lpmvol_p*lpmdens_p/lpmtau_p

         rcd_std = 1.+0.15*lpmre_p**(0.687) + 
     >               0.42*(1.+42500./lpmre_p**(1.16))**(-1)
         rdum = rdum*rcd_std

         if (part_force(1) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum  = rdum*(1. - 2.*rphip)/(1. - rphip)**3
         endif

         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

      elseif (abs(part_force(1)).eq.4) then
         ! ergun
         rbeta1 = 150.*lpmvolfrac_p*lpmvolfrac_p*lpmvisc_f/lpmvolfrac_f/
     >            lpmdiam_p**2 + 1.75*lpmvolfrac_p*lpmdens_f*
     >            lpmvdiff_pf/lpmdiam_p

         ! wen-yu
         rrep = lpmvolfrac_f*lpmre_p
         if (rrep .lt. 1000) then
            rcd = 24./rrep*(1. + 0.15*rrep**(0.687))
         else
            rcd = 0.44
         endif

         rbeta2 = 0.75*rcd*lpmvolfrac_f**(-2.65)*
     >        lpmvolfrac_p*lpmvolfrac_f*lpmdens_f*lpmvdiff_pf/lpmdiam_p

         ! stiching
         rs = 0.2
         rpp = atan(150 * 1.75*(lpmvolfrac_p - rs))/pi + 0.5
         rbeta = rpp*rbeta1 + (1.-rpp)*rbeta2


         rdum = lpmvol_p*rbeta/lpmvolfrac_p
         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

         if (abs(lpmvolfrac_p) .lt. 1E-12) then
             write(6,*) 'Use different drag model w/o volume frac /0'
             call exitt
         endif

      elseif (part_force(1).eq.0) then
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_f_un
c
c     undisturbed force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (abs(part_force(2)) .ne. 0) then
         lpmforce(1) = -lpmvol_p*lpmDuDt(1)
         lpmforce(2) = -lpmvol_p*lpmDuDt(2)
         lpmforce(3) = 0.0
         if (if3d) lpmforce(3) = -lpmvol_p*lpmDuDt(3)
      elseif (part_force(2) .eq. 0) then
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_f_iu
c
c     inviscid unsteady force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (abs(part_force(3)) .eq. 1) then
         ! shape
         rdum = 0.5

         ! volume fraction
         if (part_force(3) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum = rdum*(1.+2.*rphip)/(1.-rphip)
         endif

         ! set coeff for implicit dv/dt solve
         rpart(jcmiu,lpmi) = rdum
         
         ! but don't add a fluid contribution for now, but should later
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0

      elseif (abs(part_force(3)) .eq. 2) then
         ! shape
         rdum = 0.5
         ! mach
         rdum = rdum*(1. + 1.8*lpmmach_p**2 + 7.6*lpmmach_p**4)

         ! volume fraction
         if (part_force(3) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum = rdum*(1.+2.*rphip)/(1.-rphip)
         endif

         ! set coeff for implicit dv/dt solve
         rpart(jcmiu,lpmi) = rdum
         
         ! but don't add a fluid contribution for now, but should later
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0

      elseif (part_force(3) .eq. 0) then

         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_q_qs
c
c     quasi-steady heat transfer
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real nu_g,rdum,rpra,kappa_g

      nu_g    = lpmvisc_f/lpmdens_f ! kinematic viscosity
      rpra    = nu_g/lpmkappa_f
      rdum    = 2.*pi*lpmkappa_f*lpmdiam_p

      if (abs(part_force(4)) .gt. 0) then
         rdum    = rdum*(1. + 0.3*sqrt(lpmre_p)*rpra**(2./3.)) !re < 500
         lpmheat = rdum*(lpmtemp_f - lpmtemp_p)
      elseif (part_force(4) .eq. 0) then
         lpmheat = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_q_uu
c
c     undisturbed heat transfer
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (part_force(5) .eq. 0) then
         lpmheat = 0.0 ! nothing for now
      endif

      return
      end
c-----------------------------------------------------------------------
