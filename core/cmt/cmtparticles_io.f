c----------------------------------------------------------------------
      subroutine usr_particles_io(nistep) ! nistep not used, remove
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*15 locstring, velstring, datastring
      integer oldfile, wdata_type, newfile
      integer*8 disp, stride_len 
      integer status_mpi(MPI_STATUS_SIZE)
      integer prevs(0:np-1),npt_total,e
      real realtmp(4,llpart), vmag, rfpfluid(3),rfpfluidl(3),
     >     dum_fld(lx1,ly1,lz1,lelt),msum,msum_tot(3,2)

      ! begin timer
      ptdum(21) = dnekclock()

      if (.true.) then ! skip/no skip particle mpi io

c     setup files to write to mpi -------------------------------------
      icalld = icalld+1
      write(locstring,'(A7,I5.5,A3)') 'partxyz', icalld, '.3D' 
      write(datastring,'(A8,I5.5)') 'partdata', icalld

      do i = 1,n
         realtmp(1,i) = rpart(jx,i)
         realtmp(2,i) = rpart(jy,i)
         realtmp(3,i) = rpart(jz,i)

         vmag = rpart(jv0,i)**2 + rpart(jv0+1,i)**2
         if (if3d) vmag = vmag + rpart(jv0+2,i)**2
         vmag = sqrt(vmag)

         realtmp(4,i) = real(ipart(jai,i))
      enddo

      call MPI_Send(n, 1, MPI_INTEGER, 0, 0, nekcomm, ierr)
      npt_total = iglsum(n,1)

c     write particle data to read into script
      if (nid.eq. 0) then
c         output data so files can be easily converted to binary
c     dz comment for no particle io (begin)
c         open(364, file=datastring, action="write")
c         write(364,*) npt_total
c         close(364)
c     dz comment for no particle io (end)
          prevs(0) = n
          do i=1,np-1
             call MPI_Recv(prevs(i),1,MPI_INTEGER,i,
     >                     0,nekcomm,status_mpi,ierr)
          enddo
      endif
      call MPI_BCAST(prevs,np, MPI_INTEGER,0,nekcomm,ierr) 

      stride_len = 0.0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(i-1)
      enddo
      endif


c     print out particle values -----------------------------------
c     dz comment for no particle io (begin)
c     call MPI_FILE_OPEN(nekcomm, locstring,
c    >                   MPI_MODE_CREATE + MPI_MODE_WRONLY, 
c    >                   MPI_INFO_NULL, oldfile, ierr) 

c     disp = stride_len*4*8
c     call MPI_FILE_SET_VIEW(oldfile, disp, MPI_DOUBLE_PRECISION,
c    >                       MPI_DOUBLE_PRECISION, "native", 
c    >                       MPI_INFO_NULL, ierr) 
c     call MPI_FILE_WRITE(oldfile, realtmp(1,1), n*4,
c    >                  MPI_DOUBLE_PRECISION,
c    >                  MPI_STATUS_IGNORE, ierr) 
c     call MPI_FILE_CLOSE(oldfile, ierr) 
c     dz comment for no particle io (end)

c     output grid data for particles two-way coupled ------------------
c        note in visit, xvel = fhydx, yvel = fhydy,
c        zvel = fhydz, pressure = phi_p 
      itmp = 1
      call rzero(dum_fld,lx1*ly1*lz1*lelt)
      call outpost2(ptw(1,1,1,1,1),         ! fhyd_x
     >              ptw(1,1,1,1,2),         ! fhyd_y
     >              ptw(1,1,1,1,3),         ! fhyd_z
     >              ptw(1,1,1,1,4),         ! phi_p
     >              ptw(1,1,1,1,5),         ! hydro energy coupling
     >              itmp          ,        
     >              'ptw')

c     eulerian integrations -----------------------------------------
c     fluid momentum cmtnek
c     do ieq=2,4 
c        msum = 0.0
c        do e=1,nelt
c        do k=1,nz1
c        do j=1,ny1
c        do i=1,nx1
c           msum = msum + (u(i,j,k,ieq,e)*bm1(i,j,k,e))
c        enddo
c        enddo
c        enddo
c        enddo
c        msum_tot(ieq-1,1) = glsum(msum,1)
c     enddo 
c     fluid momentum nek5000
      msum_tot(1,1) = glsc3(bm1,vtrans,vx,nx1*ny1*nz1*nelv)
      msum_tot(2,1) = glsc3(bm1,vtrans,vy,nx1*ny1*nz1*nelv)
      msum_tot(3,1) = glsc3(bm1,vtrans,vz,nx1*ny1*nz1*nelv)
c     particle volume fraction
      vf_part_e     = glsc2(bm1,ptw(1,1,1,1,4),nx1*ny1*nz1*nelt)
c     particle forces on fluid
      rfpfluid(1)   = glsc2(bm1,ptw(1,1,1,1,1),nx1*ny1*nz1*nelt)
      rfpfluid(2)   = glsc2(bm1,ptw(1,1,1,1,2),nx1*ny1*nz1*nelt)
      rfpfluid(3)   = glsc2(bm1,ptw(1,1,1,1,3),nx1*ny1*nz1*nelt)

c     lagrangian integrations ---------------------------------------
c     particle momentum
      do ieq=0,2
         msum = 0.0
         rsum = 0.0
         do i=1,n
           msum = msum + 
     >       rpart(jspl,i)*rpart(jv0+ieq,i)*rpart(jrhop,i)*rpart(jvol,i)
           rsum = rsum + rpart(jspl,i)*rpart(jf0+ieq,i)
         enddo
         msum_tot(ieq+1,2) = glsum(msum,1)
         rfpfluidl(1+ieq)  = glsum(rsum,1)
      enddo
c     particle volume fraction
      msum = 0.0
      do i=1,n
         msum = msum + rpart(jspl,i)*rpart(jvol,i)
      enddo
      vf_part_l = glsum(msum,1)

c     print to files ------------------------------------------------
c     print properties to logfile
      if (nid.eq.0) write(6,500) "--- Eulerian Properties ------"
      if (nid.eq.0) write(6,500) "Fluid Momentum :              ", 
     >                  msum_tot(1,1),msum_tot(2,1),msum_tot(3,1)
      if (nid.eq.0) write(6,500) "Particle forces:              ", 
     >                  rfpfluid(1),rfpfluid(2),rfpfluid(3)         
      if (nid.eq.0) write(6,500) "Particle Volume:              ", 
     >                  vf_part_e
      if (nid.eq.0) write(6,500) "--- Lagrangian Properties --- "
      if (nid.eq.0) write(6,500) "Particle Momentum :           ", 
     >                  msum_tot(1,2),msum_tot(2,2),msum_tot(3,2)
      if (nid.eq.0) write(6,500) "Particle forces:              ", 
     >                  rfpfluidl(1),rfpfluidl(2),rfpfluidl(3)         
      if (nid.eq.0) write(6,500) "Particle Volume:              ", 
     >                  vf_part_l
      if (nid.eq.0) then 
         open(1511, file="MOM", action="write",position="append")
                    write(1511,"(9ES20.10)") 1.*istep,
     >                          msum_tot(1,1),msum_tot(1,2),
     >                          msum_tot(2,1),msum_tot(2,2),
     >                          msum_tot(3,1),msum_tot(3,2),
     >                          vf_part_e    ,vf_part_l
         close(1511)
      endif

c     print values for one particle ----------------------------------
      do i =1,n
         if (ipart(jpid1,i) .eq. 0) then   ! started on proc 0
         if (ipart(jpid2,i) .eq. 1) then   ! id 1 at time 0
         if (ipart(jpid3,i) .eq. 0) then   ! time 0
            open(1511, file="PID", action="write",position="append")
            write(6,600) 'PID: ', 1.*istep !2
     >          ,rpart(jx,i)    !3
     >          ,rpart(jy,i)    !4
     >          ,rpart(jz,i)    !5
     >          ,rpart(jv0,i)   !6
     >          ,rpart(jv0+1,i) !7
     >          ,rpart(jv0+2,i) !8
     >          ,real(ipart(jpnn,i))!9
            write(1511,"(20ES20.10)") 1.*istep !1
     >          ,rpart(jx,i)    !2
     >          ,rpart(jy,i)    !3
     >          ,rpart(jz,i)    !4
     >          ,rpart(jv0,i)   !5
     >          ,rpart(jv0+1,i) !6
     >          ,rpart(jv0+2,i) !7
     >          ,rpart(jfusr,i) !8
     >          ,rpart(jfusr+1,i) !9
     >          ,rpart(jfusr+2,i) !10
     >          ,rpart(jfqs,i) !11
     >          ,rpart(jfqs+1,i) !12
     >          ,rpart(jfqs+2,i) !13
     >          ,rpart(jfun,i) !14
     >          ,rpart(jfun+1,i) !15
     >          ,rpart(jfun+2,i) !16
     >          ,rpart(jfiu,i) !17
     >          ,rpart(jfiu+1,i) !18
     >          ,rpart(jfiu+2,i) !19
c    >          ,real(ipart(jpnn,i))!9
            close(1511)
         endif
         endif
         endif
      enddo      

  500 FORMAT(A30,9ES20.10)
  600 FORMAT(A9,8ES20.10)

      ! end timer
      pttime(21) = pttime(21) + dnekclock() - ptdum(21)

      endif ! end skip/no skip particle mpi io

      ! begin timer
      ptdum(22) = dnekclock()
      call output_along_line_avg
      ! end timer
      pttime(22) = pttime(22) + dnekclock() - ptdum(22)

c     output timings
      call output_particle_timers 

      return
      end
c----------------------------------------------------------------------
      subroutine output_along_line_avg
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*15 outstring
      integer*8 disp, stride_len 
      integer status_mpi(MPI_STATUS_SIZE)
      real send_vals(1+50)
      real rxgls(lx1),uf(lx1,ly1,lz1,lelt,50),rcount(50),rdum(50),
     >     rtmp(50)
      integer nlfl

      real   xdrange(2,3)
      common /domainrange/ xdrange
      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      common /running_avgs/ rec_vals
      real rec_vals(1+50,1000*15) !1000 elements, by nx1=15 max


      nlfl = 21 

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         uf(i,j,k,ie,1) = ptw(i,j,k,ie,1)
         uf(i,j,k,ie,2) = ptw(i,j,k,ie,2)
         uf(i,j,k,ie,3) = ptw(i,j,k,ie,3)
         uf(i,j,k,ie,4) = ptw(i,j,k,ie,4)
         uf(i,j,k,ie,5) = ptw(i,j,k,ie,5)
         uf(i,j,k,ie,6) = ptw(i,j,k,ie,6)
         uf(i,j,k,ie,7) = ptw(i,j,k,ie,7) 
         uf(i,j,k,ie,8) = ptw(i,j,k,ie,8) 
         uf(i,j,k,ie,9) = rhs_fluidp(i,j,k,ie,1)
         uf(i,j,k,ie,10)= rhs_fluidp(i,j,k,ie,2)
         uf(i,j,k,ie,11)= rhs_fluidp(i,j,k,ie,3)
         uf(i,j,k,ie,12)= rhs_fluidp(i,j,k,ie,4)
         uf(i,j,k,ie,13)= rhs_fluidp(i,j,k,ie,5)
         uf(i,j,k,ie,14)= rhs_fluidp(i,j,k,ie,6)
         uf(i,j,k,ie,15)= rhs_fluidp(i,j,k,ie,7)
         uf(i,j,k,ie,16)= vx(i,j,k,ie)
         uf(i,j,k,ie,17)= vy(i,j,k,ie)
         uf(i,j,k,ie,18)= vz(i,j,k,ie)
         uf(i,j,k,ie,19)= pr(i,j,k,ie)
         uf(i,j,k,ie,20)= vtrans(i,j,k,ie,1)
         uf(i,j,k,ie,21)= t(i,j,k,ie,1)
      enddo
      enddo
      enddo
      enddo

      icalld = icalld+1
      write(outstring,'(A8,I5.5)') 'avgsdata', icalld

      do i=1,nx1 
         rxgls(i) = (xgll(i) + 1.)*rleng/2.
      enddo
      
      rthresh = 1E-12
      rxs = xdrange(1,1)
      rxt = rxs
      icm = 1
      do while (abs(rxs - xdrange(2,1)) .ge. rthresh)

         do i=1,nx1
            rxs = rxt + rxgls(i)

            call rzero(rdum,nlfl)
            call rzero(rcount,nlfl)

            do ie=1,nelt
            do ik=1,nz1
            do ij=1,ny1
            do ii=1,nx1
               rxv = xm1(ii,ij,ik,ie)
               if (abs(rxv - rxs) .lt. rthresh) then
                  do j = 1,nlfl
                     rdum(j) = rdum(j) + uf(ii,ij,ik,ie,j)
                     rcount(j) = rcount(j) + 1.
                  enddo
               endif
            enddo
            enddo
            enddo
            enddo

            isz = 1
            rec_vals(1,icm) = rxs

            do j = 1,nlfl
               rec_vals(j+1,icm) =  glsum(rdum(j),1)
               rtmp(j) = glsum(rcount(j),1)
               rec_vals(j+1,icm) = rec_vals(j+1,icm)/rtmp(j)
            enddo
c           rec_vals(2,icm) = rec_vals(2,icm)
c           rec_vals(2,icm) = rtmp
c           call mpi_allreduce(rdum,rec_vals(2,icm),isz,MPI_REAL,MPI_SUM
c    >                           ,nekcomm,ierr)
c           call mpi_allreduce(icount,isum,isz,MPI_INTEGER,MPI_SUM
c    >                           ,nekcomm,ierr)
c           rec_vals(2,icm) = rec_vals(2,icm)/isum
            icm = icm + 1
         enddo

         rxt = rxs
      enddo

      if (nid.eq. 0) then
          open(364, file=outstring, action="write",position="append")
          do i =1,icm-1 ! last point has issues
             write(364,600) rec_vals(1,i),
     >                      rec_vals(2,i),
     >                      rec_vals(3,i),
     >                      rec_vals(4,i),
     >                      rec_vals(5,i),
     >                      rec_vals(6,i),
     >                      rec_vals(7,i),
     >                      rec_vals(8,i),
     >                      rec_vals(9,i),
     >                      rec_vals(10,i),
     >                      rec_vals(11,i),
     >                      rec_vals(12,i),
     >                      rec_vals(13,i),
     >                      rec_vals(14,i),
     >                      rec_vals(15,i),
     >                      rec_vals(16,i),
     >                      rec_vals(17,i),
     >                      rec_vals(18,i),
     >                      rec_vals(19,i),
     >                      rec_vals(20,i),
     >                      rec_vals(21,i),
     >                      rec_vals(22,i) ! add one
          enddo
          close(364)
      endif
      
  600 FORMAT(22ES20.10)
      return
      end
c----------------------------------------------------------------------
      subroutine read_particle_input(ifreadpart)
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'CMTPART'

      logical ifreadpart

      open(unit=81,file="particles.inp",form="formatted")
      read (81,*) nw
      read (81,*) dp(1)
      read (81,*) dp(2)
      read (81,*) cp_p
      read (81,*) kappa_g
      read (81,*) df_dx
      read (81,*) rleng
      read (81,*) nitspl
      read (81,*) ralphdecay
      read (81,*) time_integ
      read (81,*) bc_part(1)
      read (81,*) bc_part(2)
      read (81,*) bc_part(3)
      read (81,*) bc_part(4)
      read (81,*) bc_part(5)
      read (81,*) bc_part(6)
      read (81,*) two_way
      read (81,*) red_interp
      read (81,*) nrect_assume
      read (81,*) part_force(1)
      read (81,*) part_force(2)
      read (81,*) part_force(3)
      read (81,*) part_force(4)
      read (81,*) part_force(5)
      read (81,*) part_force(6)
      read (81,*) time_delay
      read (81,*) inject_rate
      read (81,*) ipart_restart 
      read (81,*) rho_p
      read (81,*) mu_0
      read (81,*) phi_desire
      read (81,*) rxbo(1,1)
      read (81,*) rxbo(2,1)
      read (81,*) rxbo(1,2)
      read (81,*) rxbo(2,2)
      read (81,*) rxbo(1,3)
      read (81,*) rxbo(2,3)
      read (81,*) nrandseed
      close(81)

      ifreadpart = .false.
      if (ipart_restart.gt.0) then
         ifreadpart = .true.

      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_particle_timers
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'CMTPART'

      if(mod(istep,iostep).eq.0.or. istep.eq.1) then
c        first compute totals
         rdum  = ftime/istep
         dtime = glsum(rdum,1)
         rtime_total = dtime/np ! both pinit and psolve in here

         rdum  = pttime(1)/istep
         dtime = glsum(rdum,1)
         rtime_pinit = dtime/np

         rdum  = pttime(5)/istep
         dtime = glsum(rdum,1)
         rtime_psolve = dtime/np

         rtime_fsolve = rtime_total - rtime_psolve - rtime_pinit
         rtime_fpsolve = rtime_fsolve + rtime_psolve

         if(nio.eq.0) 
     >      write(6,500) istep,rtime_fpsolve,'   ',
     >                 'TOTAL TIME                          '
         if(nio.eq.0) 
     >      write(6,500) istep,rtime_fsolve,'   ',
     >                 'FLUID TIME                          '
         if(nio.eq.0) 
     >      write(6,500) istep,rtime_psolve,'   ',
     >                 'PARTICLE TIME                       '

c        now compute detailed particle times - next level
         if(nio.eq.0)
     >   write(6,500) istep,rtime_pinit, '   ',
     >              '0.0 usr_particles_init                 '

         rdum  = pttime(2)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '0.1 set_bounds_box                     '

         rdum  = pttime(4)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '0.2 set_part_pointers                  '

         rdum  = pttime(3)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '0.3 place_particles                    '

         if(nio.eq.0)
     >   write(6,500) istep,rtime_psolve, '   ',
     >              '1.0 usr_particles_solver               '

         rdum  = pttime(17)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.1 update_particle_location           '

         rdum  = pttime(23)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.2 move_particles_inproc              '

         rdum  = pttime(24)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.2.1 particles_in_nid                 '

         rdum  = pttime(25)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.2.2 update_findpts_info              '

         rdum  = pttime(20)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.3 interp_props_part_location         '

         rdum  = pttime(18)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.3.1 baryinterp                       '

         rdum  = pttime(19)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.3.2 triinterp                        '

         rdum  = pttime(11)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.4 usr_particles_forcing              '

         rdum  = pttime(13)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.4.1 calc_substantial_derivative      '

         rdum  = pttime(10)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.5 rk3_integrate                      '

         rdum  = pttime(12)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.6 compute_forcing_post_part          '

         rdum  = pttime(14)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.7 particles_solver_nearest_neighbors '

         rdum  = pttime(16)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.7.1 create_ghost_particles           '

         rdum  = (pttime(14)-pttime(16)-pttime(15))/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.7.2 send ghost particles             '

         rdum  = pttime(15)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.7.3 search_nearest_neighbors         '

         rdum  = pttime(6)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.8 spread_props_grid                  '

         rdum  = pttime(7)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.8.1 local_part_to_grid               '

         rdum  = pttime(8)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.8.2 remote_part_to_grid              '

         rdum  = pttime(9)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.8.12.1 point_to_grid                 '

         rdum  = pttime(21)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.9 usr_particles_io                   '

         rdum  = pttime(22)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nio.eq.0)
     >   write(6,500) istep,rtime, '   ',
     >              '1.10 output_along_line_avg             '


      endif

  500 FORMAT(I13,ES13.4,A3,A39)

      return
      end
c----------------------------------------------------------------------
      subroutine output_particle_options
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      npt_total = iglsum(n,1) ! recompute for sanity check

c     print particle options -------------------------------------------
      if (nid.eq. 0) then
c     headers here
      write(6,100)
     >'/-------------------------------------------------------------//'
      write(6,100)
     >'/                       Particle Options                      //'
      write(6,100)
     >'/-------------------------------------------------------------//'

c     start printing options
      write(6,200) 
     >     'Number of particles                               : ',
     >     npt_total 
      write(6,200) 
     >     'Solver (0=bdf/ext tracer,1=rk3,2=bdf)             : ',
     >     time_integ
      write(6,300)
     >     'Reduced interpolation                             : ',
     >     red_interp , '/',lx1
      write(6,200)
     >     'Delay particles until timestep                    : ',
     >     time_delay 
      write(6,200)
     >     'Two-way coupled (0=no,1=yes)                      : ',
     >     two_way
      write(6,400)
     >     'Forces (0=no,1=yes,-1=on w/o corr;usr,qs,un,iu)   : ',
     >     part_force(1),',',part_force(2),',',part_force(3),',',
     >     part_force(4)
      write(6,500)
     >     'Particle density [kg/m^3]                         : ',
     >     rho_p
      write(6,500)
     >     'Particle diameter low [m]                         : ',
     >     dp(1)
      write(6,500)
     >     'Particle diameter high [m]                        : ',
     >     dp(2)
      write(6,500)
     >     'Fluid viscosity for drag [Pa s]                   : ',
     >     mu_0 
      write(6,500)
     >     'Average particle volume fraction                  : ',
     >     phi_desire

c     footers here
      write(6,100)
     >'//-------------------------------------------------------------/'
      write(6,100)
     >'/                     End Particle Options                    //'
      write(6,100)
     >'//-------------------------------------------------------------/'
      endif

  100 FORMAT(A65) 
  200 FORMAT(A52,I13) 
  300 FORMAT(A52,I10,A1,I2) 
  400 FORMAT(A52,I4,A1,I2,A1,I2,A1,I2) 
  500 FORMAT(A52,ES13.4)

      return
      end
c----------------------------------------------------------------------
