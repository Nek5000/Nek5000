c----------------------------------------------------------------------
c routine called in case of particle calls only in .usr file (i.e.,
c  with nek5000 particles and not cmt-nek particles. must use 
c  bdf/ext time integration. Otherwise, cmt-nek will not call this fxn.
      subroutine stokes_particles
      include 'SIZE'
      include 'TOTAL'

      if (istep.eq.0) then
         call usr_particles_init
      else
         call usr_particles_solver
      endif

      if(mod(istep,iostep).eq.0.or. istep.eq.1) then
         call usr_particles_io(istep)
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine usr_particles_solver
c
c     call routines in ordered way - main solver structure
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      logical ifinject
      integer icalld
      save    icalld
      data    icalld  /-1/

      ! begin timer
      ptdum(5) = dnekclock()

c     should we inject particles
      ifinject = .false.
      if (inject_rate .gt. 0) then
      if ((mod(istep,inject_rate).eq.0)) then 
         ifinject = .true. 
      endif
      endif

      if (istep .gt. time_delay) then

c     scheme 1 --------------------------------------------------------
      if (time_integ .eq. 0) then           

c     rk3 integration -------------------------------------------------
      elseif (time_integ .eq. 1) then       
         if (stage.eq.1) then
            call update_particle_location   ! move outlier particles
            if (ifinject) call place_particles ! inject particles
            call move_particles_inproc      ! update mpi rank
         endif
         call interp_props_part_location    ! interpolate
         call usr_particles_forcing         ! fluid to part. forcing
         call rk3_integrate                 ! time integration
         call compute_forcing_post_part     ! update forces
         if (two_way.eq.1) then             ! part. to fluid forcing
            call particles_solver_nearest_neighbor    ! nn
            call spread_props_grid          ! put particle props on grid
         endif

c     Other -----------------------------------------------------------
      elseif (time_integ .eq. 2) then

      endif ! particle scheme

      endif ! time_delay

      ! end timer
      pttime(5) = pttime(5) + dnekclock() - ptdum(5)

      return
      end
c-----------------------------------------------------------------------
      subroutine rk3_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      common /PARTRK3/ kv_stage_p
      real   kv_stage_p(llpart,7)

      integer fdim
      real    pmass

      ! begin timer
      ptdum(10) = dnekclock()

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
         rpart(jx0  ,i) = tcoef(1,stage)*kv_stage_p(i,1) +
     >                    tcoef(2,stage)*rpart(jx0  ,i)  +
     >                    tcoef(3,stage)*rpart(jv0  ,i)
         rpart(jx0+1,i) = tcoef(1,stage)*kv_stage_p(i,2) +
     >                    tcoef(2,stage)*rpart(jx0+1,i)  +
     >                    tcoef(3,stage)*rpart(jv0+1,i)
         rpart(jx0+2,i) = tcoef(1,stage)*kv_stage_p(i,3) +
     >                    tcoef(2,stage)*rpart(jx0+2,i)  +
     >                    tcoef(3,stage)*rpart(jv0+2,i)
         rpart(jv0  ,i) = tcoef(1,stage)*kv_stage_p(i,4) +
     >                    tcoef(2,stage)*rpart(jv0  ,i)  +
     >                    tcoef(3,stage)*rpart(jf0  ,i)
         rpart(jv0+1,i) = tcoef(1,stage)*kv_stage_p(i,5) +
     >                    tcoef(2,stage)*rpart(jv0+1,i)  +
     >                    tcoef(3,stage)*rpart(jf0+1,i)
         rpart(jv0+2,i) = tcoef(1,stage)*kv_stage_p(i,6) +
     >                    tcoef(2,stage)*rpart(jv0+2,i)  +
     >                    tcoef(3,stage)*rpart(jf0+2,i)
         rpart(jtemp,i) = tcoef(1,stage)*kv_stage_p(i,7) +
     >                    tcoef(2,stage)*rpart(jtemp,i)  +
     >                    tcoef(3,stage)*rpart(jq0  ,i)
      enddo

      ! end timer
      pttime(10) = pttime(10) + dnekclock() - ptdum(10)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_bdf_ext_coefs(beta,alpha,times)
      include 'SIZE'
      include 'TOTAL'
      include 'CMTPART'

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

