c-----------------------------------------------------------------------
      subroutine particles_solver_nearest_neighbor
c
c     this routine will let particles search for their nearest neighbors
c     using the ghost particle approach.
c
c     bc_part = -1,1  => non-periodic search
c     bc_part = 0  => periodic search
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange
      real   xdrange(2,3)
      common /domainrange/ xdrange

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      logical partl         ! dummy used in c_t_t()

      ! begin timer
      ptdum(14) = dnekclock()

c     distance to check in x,y,z space
      if (istep.eq.0.or.istep.eq.1) then
         rtmp = df_dx/(2.*nx1)*sqrt(-log(ralphdecay)/log(2.))*rleng
         if (rtmp .gt. 0.5) then
         if (nid .eq. 0) 
     >      write (6,*)'error: r/L_e is too large,change part.',
     >                                                      rtmp/rleng
            rtmp = 0.5*rleng
         else
         endif

         d2chk(1) = rtmp
         d2chk(2) = rtmp
         d2chk(3) = rtmp
      endif

      if (istep.eq.0.or.istep.eq.1) then
         ntmp = iglsum(nfptsgp,1)
         if (nid.eq.0) write(6,*) 'Passed init ghost parts', rtmp/rleng
      endif

c     create ghost particles
      if (nrect_assume .eq. 1) call create_ghost_particles_rect
      if (nrect_assume .eq. 0) call create_ghost_particles_gen

      if (istep.eq.0.or.istep.eq.1) then
         ntmp = iglsum(nfptsgp,1)
         if (nid.eq.0) write(6,*) 'Passed create_ghost_particles',ntmp
      endif

c     send ghost particles
      call crystal_tuple_transfer(i_cr_hndl,nfptsgp,llpart
     $           , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpps) ! jgpps is overwri

      if (istep.eq.0.or.istep.eq.1) then
         ntmp = iglmax(nfptsgp,1)
         if (nid.eq.0) write(6,*) 'Passed send ghost particles',ntmp
      endif

c     search nearest neighbors from this proc particles and
c     remote proc nearby particles (ghost particles)
c     call search_nearest_neighbor

      ! end timer
      pttime(14) = pttime(14) + dnekclock() - ptdum(14)

      return
      end
c-----------------------------------------------------------------------
      subroutine search_nearest_neighbor
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'
c
c     this routine implements a naive particle search. This should be
c     updated for the future with some kind of more recent algorithm for
c     nearest neighbor searching. Particles will need to search local
c     particles (in rpart and ipart) and remote particles that are
c     nearby but on different MPI ranks (rptsgp and iptsgp of length nfptsgp)
c
c     particles will check if they are within d2chk of each other
c

      integer nneigh

      ! begin timer
      ptdum(15) = dnekclock()

      d3 = 0.5*d2chk(1) ! user can change, but d2chk is robust max value
                        ! Note: 1/2*d2chk seems to work even w/outflow

c     let every particle search for itself
      do i = 1,n
         ipart(jai,i) = ipart(jpnn,i) ! for testing
         nneigh = 0
c        particles in local elements
         do j = 1,n
            if (i .ne. j) then
               pdist = abs(rpart(jx,i)-rpart(jx,j))**2  
     >                          + abs(rpart(jy,i)-rpart(jy,j))**2
     >                          + abs(rpart(jz,i)-rpart(jz,j))**2
               pdist = sqrt(pdist)
               if (pdist .gt. d3) goto 1109
               nneigh = nneigh + 1
            endif
1109        continue
         enddo

c        search list of ghost particles
         do j = 1,nfptsgp
            if (iptsgp(jgpes,j).eq. ipart(je0,i)) then ! exclude ghosts not
                                                      ! meant for this eleme
            pdist = abs(rpart(jx,i)-rptsgp(jgpx,j))**2  
     >                    + abs(rpart(jy,i)-rptsgp(jgpy,j))**2
     >                    + abs(rpart(jz,i)-rptsgp(jgpz,j))**2
            pdist = sqrt(pdist)
            if (pdist .gt. d3) goto 11092
            nneigh = nneigh + 1
            endif
11092       continue
         enddo
         ipart(jpnn,i) = nneigh
         ipart(jai,i) = ipart(jai,i) - ipart(jpnn,i) ! comptued distance
                                                     ! for testing
      enddo

      ! end timer
      pttime(15) = pttime(15) + dnekclock() - ptdum(15)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_ghost_particles_gen
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
      include 'CMTDATA'
      include 'CMTPART'

      common /myparth/ i_fp_hndl, i_cr_hndl

      real rcoords_gp(3,7*n) ! possible seven gp
      real rdums(4,7*n),rdist
      integer idums(3,7*n),ieremove(nelt)

      ! begin timer
      ptdum(16) = dnekclock()

      ic = 0
      do i = 1,n

         ie = ipart(je0,i) + 1

         isgnx = 1
         isgny = 1
         isgnz = 1
         if (rpart(jr+0,i) .lt. 0.) isgnx = -1
         if (rpart(jr+1,i) .lt. 0.) isgny = -1
         if (rpart(jr+2,i) .lt. 0.) isgnz = -1

         ! difference in d2chk?
         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i) + isgnx*d2chk(1)
         rcoords_gp(2,ic) = rpart(jy,i)
         rcoords_gp(3,ic) = rpart(jz,i)

         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i)
         rcoords_gp(2,ic) = rpart(jy,i) + isgny*d2chk(1)
         rcoords_gp(3,ic) = rpart(jz,i)

         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i)
         rcoords_gp(2,ic) = rpart(jy,i)
         rcoords_gp(3,ic) = rpart(jz,i) + isgnz*d2chk(1)

         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i) + isgnx*d2chk(1)
         rcoords_gp(2,ic) = rpart(jy,i) + isgny*d2chk(1)
         rcoords_gp(3,ic) = rpart(jz,i)

         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i)
         rcoords_gp(2,ic) = rpart(jy,i) + isgny*d2chk(1)
         rcoords_gp(3,ic) = rpart(jz,i) + isgnz*d2chk(1)

         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i) + isgnx*d2chk(1)
         rcoords_gp(2,ic) = rpart(jy,i)
         rcoords_gp(3,ic) = rpart(jz,i) + isgnz*d2chk(1)

         ic = ic + 1
         rcoords_gp(1,ic) = rpart(jx,i) + isgnx*d2chk(1)
         rcoords_gp(2,ic) = rpart(jy,i) + isgny*d2chk(1)
         rcoords_gp(3,ic) = rpart(jz,i) + isgnz*d2chk(1)

c        do j=1,7
c           if (idums(2,j) .ne. ipart(jps,i)) then
c           if (idums(1,j) .eq. 0) then
c              nfptsgp = nfptsgp + 1

c              ! make sure if future changes here to change in rect only
c              ! function, or else we may have bug later on
c              rptsgp(jgpx,nfptsgp)    = rpart(jx,i)    ! x loc
c              rptsgp(jgpy,nfptsgp)    = rpart(jy,i)    ! y log
c              rptsgp(jgpz,nfptsgp)    = rpart(jz,i)    ! z log
c              rptsgp(jgpfh,nfptsgp)   = rpart(jf0,i)   ! hyd. force x
c              rptsgp(jgpfh+1,nfptsgp) = rpart(jf0+1,i) ! hyd. force y
c              rptsgp(jgpfh+2,nfptsgp) = rpart(jf0+2,i) ! hyd. force z
c              rptsgp(jgpvol,nfptsgp)  = rpart(jvol,i)  ! particle volum
c              rptsgp(jgpgam,nfptsgp)  = rpart(jgam,i)  ! spread correct
c              rptsgp(jgpspl,nfptsgp)  = rpart(jspl,i)  ! super particle
c              rptsgp(jgpg0,nfptsgp)   = rpart(jg0,i)   ! work done by forc
c              rptsgp(jgpq0,nfptsgp)   = rpart(jq0,i)   ! heating from part 
c              rptsgp(jgpv0,nfptsgp)   = rpart(jv0,i)   ! particle velocity
c              rptsgp(jgpv0+1,nfptsgp) = rpart(jv0+1,i) ! particle velocity
c              rptsgp(jgpv0+2,nfptsgp) = rpart(jv0+2,i) ! particle velocity
c
c              iptsgp(jgppid1,nfptsgp) = ipart(jpid1,i)          ! part id 1 tag
c              iptsgp(jgppid2,nfptsgp) = ipart(jpid2,i)          ! part id 2 tag
c              iptsgp(jgppid3,nfptsgp) = ipart(jpid3,i)          ! part id 3 tag
c              iptsgp(jgpps,nfptsgp)   = iprocmap  ! overwritten mpi
c              iptsgp(jgppt,nfptsgp)   = iprocmap  ! dest. mpi rank
c              iptsgp(jgpes,nfptsgp)   = ielmap    ! dest. elment
c           endif
c           endif
c        enddo
      enddo


      nigpl = 3
      nr1  = 4
      nr2  = 3
            call findpts(i_fp_hndl !  stride     !   call findpts( ihndl,
     $           , idums(1,1),nigpl       !   $             rcode,1,
     $           , idums(2,1),nigpl       !   &             proc,1,
     $           , idums(3,1),nigpl       !   &             elid,1,
     $           , rdums(1,1),nr1        !   &             rst,ndim,
     $           , rdums(4,1),nr1        !   &             dist,1,
     $           , rcoords_gp(1,1),nr2        !   &             pts(    1),1,
     $           , rcoords_gp(2,1),nr2        !   &             pts(  n+1),1,
     $           , rcoords_gp(3,1),nr2 ,ic)    !   &             pts(2*n+1),1,n)


      
      nelmax = 0
      npmax  = 0
      nfptsgp = 0
      icdum = ic
      ic = 0
      do i=1,n
         do j=1,7
            ic = ic + 1
            if (idums(1,ic) .ne. 2) then
            if (idums(2,ic) .ne. nid) then

               iprocmap = idums(2,ic)
               ielmap   = idums(3,ic)

               nfptsgp = nfptsgp + 1
c              ! make sure if future changes here to change in rect only
c              ! function, or else we may have bug later on
               rptsgp(jgpx,nfptsgp)    = rpart(jx,i)    ! x loc
               rptsgp(jgpy,nfptsgp)    = rpart(jy,i)    ! y log
               rptsgp(jgpz,nfptsgp)    = rpart(jz,i)    ! z log
               rptsgp(jgpfh,nfptsgp)   = rpart(jf0,i)   ! hyd. force x
               rptsgp(jgpfh+1,nfptsgp) = rpart(jf0+1,i) ! hyd. force y
               rptsgp(jgpfh+2,nfptsgp) = rpart(jf0+2,i) ! hyd. force z
               rptsgp(jgpvol,nfptsgp)  = rpart(jvol,i)  ! particle volum
               rptsgp(jgpgam,nfptsgp)  = rpart(jgam,i)  ! spread correct
               rptsgp(jgpspl,nfptsgp)  = rpart(jspl,i)  ! super particle
               rptsgp(jgpg0,nfptsgp)   = rpart(jg0,i)   ! work done by forc
               rptsgp(jgpq0,nfptsgp)   = rpart(jq0,i)   ! heating from part 
               rptsgp(jgpv0,nfptsgp)   = rpart(jv0,i)   ! particle velocity
               rptsgp(jgpv0+1,nfptsgp) = rpart(jv0+1,i) ! particle velocity
               rptsgp(jgpv0+2,nfptsgp) = rpart(jv0+2,i) ! particle velocity
 
               iptsgp(jgppid1,nfptsgp) = ipart(jpid1,i)          ! part id 1 tag
               iptsgp(jgppid2,nfptsgp) = ipart(jpid2,i)          ! part id 2 tag
               iptsgp(jgppid3,nfptsgp) = ipart(jpid3,i)          ! part id 3 tag
               iptsgp(jgpps,nfptsgp)   = iprocmap  ! overwritten mpi
               iptsgp(jgppt,nfptsgp)   = iprocmap  ! dest. mpi rank
               iptsgp(jgpes,nfptsgp)   = ielmap    ! dest. elment

c              double particles possibly created, filter out
               do ii = 1,7 ! look back at most 7 ghost particles
                  if (iptsgp(jgpps,nfptsgp-ii) .eq. 
     >                               iptsgp(jgpps,nfptsgp)) then
                  if (iptsgp(jgpes,nfptsgp-ii) .eq. 
     >                               iptsgp(jgpes,nfptsgp)) then
                  if (iptsgp(jgppid1,nfptsgp-ii) .eq.
     >                               iptsgp(jgppid1,nfptsgp)) then
                  if (iptsgp(jgppid2,nfptsgp-ii) .eq.
     >                               iptsgp(jgppid2,nfptsgp)) then
                  if (iptsgp(jgppid3,nfptsgp-ii) .eq.
     >                               iptsgp(jgppid3,nfptsgp)) then
                     nfptsgp = nfptsgp - 1
                  endif
                  endif
                  endif
                  endif
                  endif
               enddo

            endif
            endif
         enddo
      enddo

      ! end timer
      pttime(16) = pttime(16) + dnekclock() - ptdum(16)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_ghost_particles_rect
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
      include 'CMTDATA'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      ! begin timer
      ptdum(16) = dnekclock()

      nfptsgp = 0
      do i = 1,n
         ie = ipart(je0,i) + 1
c        vector coordinates of what faces a particle is next to
         ii = 0
         jj = 0
         kk = 0
         if (abs(rpart(jx,i) - xerange(1,1,ie)).lt.d2chk(1)) ii=-1
         if (abs(rpart(jx,i) - xerange(2,1,ie)).lt.d2chk(1)) ii=1
         if (abs(rpart(jy,i) - xerange(1,2,ie)).lt.d2chk(2)) jj=-1
         if (abs(rpart(jy,i) - xerange(2,2,ie)).lt.d2chk(2)) jj=1
         if (abs(rpart(jz,i) - xerange(1,3,ie)).lt.d2chk(3)) kk=-1
         if (abs(rpart(jz,i) - xerange(2,3,ie)).lt.d2chk(3)) kk=1

         itype = abs(ii)+abs(jj)+abs(kk) ! face (1), edge (2), or
                                         ! corner (3) particle

         if (itype.eq.1) then          ! face particle
            call gp_create(ii,jj,kk,i,
     >       nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
         elseif (itype.eq.2) then      ! edge particle
            call gp_create(ii,jj,kk,i,
     >       nedgegp,el_edge_num,el_edge_proc_map,el_edge_el_map)
            if (abs(ii) + abs(jj) .eq. 2) then
               call gp_create(0,jj,kk,i,
     >          nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
               call gp_create(ii,0,kk,i,
     >          nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
            elseif (abs(ii) + abs(kk) .eq. 2) then
               call gp_create(0,jj,kk,i,
     >          nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
               call gp_create(ii,jj,0,i,
     >          nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
            elseif (abs(jj) + abs(kk) .eq. 2) then
               call gp_create(ii,0,kk,i,
     >          nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
               call gp_create(ii,jj,0,i,
     >          nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
            endif
         elseif (itype.eq.3) then       ! corner particle
            call gp_create(ii,jj,kk,i,
     >       ncornergp,el_corner_num,el_corner_proc_map,
     >       el_corner_el_map)
            call gp_create(0,jj,kk,i,
     >       nedgegp,el_edge_num,el_edge_proc_map,el_edge_el_map)
            call gp_create(ii,0,kk,i,
     >       nedgegp,el_edge_num,el_edge_proc_map,el_edge_el_map)
            call gp_create(ii,jj,0,i,
     >       nedgegp,el_edge_num,el_edge_proc_map,el_edge_el_map)
            call gp_create(ii,0,0,i,
     >       nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
            call gp_create(0,jj,0,i,
     >       nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
            call gp_create(0,0,kk,i,
     >       nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
         endif
      enddo

      ! end timer
      pttime(16) = pttime(16) + dnekclock() - ptdum(16)

      return
      end
c-----------------------------------------------------------------------
      subroutine gp_create(ii,jj,kk,i,
     >             nnl,el_tmp_num,el_tmp_proc_map,el_tmp_el_map)
c
c     this routine will create a ghost particle and append its position
c     to rptsgp and its processor and element to iptsgp. nfptsgp will then
c     be incremented. Note that ghost particles will not be created if 
c     they are to be created on the same processor. In the near future, 
c     this might not be true if periodic conditions are needed.
c
c     el_tmp_num holds vector coordinates of tmp=face,edge, or corners
c     el_tmp_proc_map holds MPI rank of neighbor elements in el_tmp_num
c                     order
c     el_tmp_el_map holds local element number of neighbor elements
c
c     ii,jj,kk are vectors that tell what element a ghost particle
c     should be sent to
c
c     i is which particle is creating the ghost particle from rpart,etc
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange
      real   xdrange(2,3)
      common /domainrange/ xdrange

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      integer el_tmp_proc_map(lelt,12)  ,el_tmp_el_map(lelt,12),
     >        el_tmp_num(36)

      real rdumpos(3)

      xdlen = xdrange(2,1) - xdrange(1,1)
      ydlen = xdrange(2,2) - xdrange(1,2)
      zdlen = xdrange(2,3) - xdrange(1,3)

      ic = 0
      do j=1,3*nnl-2,3
         ic = ic + 1
         if (el_tmp_num(j)  .eq.ii) then
         if (el_tmp_num(j+1).eq.jj) then
         if (el_tmp_num(j+2).eq.kk) then
            nfptsgp = nfptsgp + 1
            ie = ipart(je0,i)+1
            iitmp1 = 0
            iitmp2 = 0
            iitmp3 = 0
            ! note that altering locs is for bc in periodic ..
            xloc = rpart(jx,i)
            if (xloc+d2chk(1)*ii .gt. xdrange(2,1)) then
                 xloc = rpart(jx,i) - xdlen
                 iitmp1 = 1
            endif
            if (xloc+d2chk(1)*ii .lt. xdrange(1,1))then
                 xloc = rpart(jx,i) + xdlen
                 iitmp1 = 1
            endif
            yloc = rpart(jy,i)
            if (yloc+d2chk(2)*jj .gt. xdrange(2,2))then
                 yloc = rpart(jy,i) - ydlen
                 iitmp2 = 1
            endif
            if (yloc+d2chk(2)*jj .lt. xdrange(1,2))then
                 yloc = rpart(jy,i) + ydlen
                 iitmp2 = 1
            endif
            zloc = rpart(jz,i)
            if (zloc+d2chk(3)*kk .gt. xdrange(2,3))then
                 zloc = rpart(jz,i) - zdlen
                 iitmp3 = 1
            endif
            if (zloc+d2chk(3)*kk .lt. xdrange(1,3))then
                 zloc = rpart(jz,i) + zdlen
                 iitmp3 = 1
            endif
            rptsgp(jgpx,nfptsgp)    = xloc           ! x loc
            rptsgp(jgpy,nfptsgp)    = yloc           ! y log
            rptsgp(jgpz,nfptsgp)    = zloc           ! z log
            rptsgp(jgpfh,nfptsgp)   = rpart(jf0,i)   ! hyd. force x
            rptsgp(jgpfh+1,nfptsgp) = rpart(jf0+1,i) ! hyd. force y
            rptsgp(jgpfh+2,nfptsgp) = rpart(jf0+2,i) ! hyd. force z
            rptsgp(jgpvol,nfptsgp)  = rpart(jvol,i)  ! particle volum
            rptsgp(jgpgam,nfptsgp)  = rpart(jgam,i)  ! spread correct
            rptsgp(jgpspl,nfptsgp)  = rpart(jspl,i)  ! super particle
            rptsgp(jgpg0,nfptsgp)   = rpart(jg0,i)   ! work done by forc
            rptsgp(jgpq0,nfptsgp)   = rpart(jq0,i)   ! heating from part 
            rptsgp(jgpv0,nfptsgp)   = rpart(jv0,i)   ! particle velocity
            rptsgp(jgpv0+1,nfptsgp) = rpart(jv0+1,i) ! particle velocity
            rptsgp(jgpv0+2,nfptsgp) = rpart(jv0+2,i) ! particle velocity

            iptsgp(jgppid1,nfptsgp) = ipart(jpid1,i)          ! part id 1 tag
            iptsgp(jgppid2,nfptsgp) = ipart(jpid2,i)          ! part id 2 tag
            iptsgp(jgppid3,nfptsgp) = ipart(jpid3,i)          ! part id 3 tag

               ipdum  = el_tmp_proc_map(ie,ic)
               iedum  = el_tmp_el_map(ie,ic)

            iptsgp(jgpps,nfptsgp)   = ipdum  ! overwritten mpi
            iptsgp(jgppt,nfptsgp)   = ipdum  ! dest. mpi rank
            iptsgp(jgpes,nfptsgp)   = iedum    ! dest. elment

c           check if extra particles have been created on the same mpi
c           rank and also take care of boundary particles
            ibctype = abs(bc_part(1))+abs(bc_part(3))+abs(bc_part(5))

c           take care of periodic stuff first
            if (nid.eq.iptsgp(jgppt,nfptsgp)) then ! dont create gp on own rank 
                                                   ! unless moved and periodic
            if (ibctype .eq. 0) then            ! all three sides periodic
               if (iitmp1+iitmp2+iitmp3 .eq.0) then
                  nfptsgp=nfptsgp-1
                  goto 1511
               endif
            elseif (ibctype .eq. 1) then        ! only two sides periodic
               if (abs(bc_part(1)) .eq. 1) then
                  if (iitmp2+iitmp3 .eq. 0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               elseif (abs(bc_part(3)) .eq. 1) then
                  if (iitmp1+iitmp3 .eq. 0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               elseif (abs(bc_part(5)) .eq. 1) then
                  if (iitmp1+iitmp2 .eq. 0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               endif
            elseif (ibctype .eq. 2) then        ! only one side periodic
               if (bc_part(1) .eq. 0) then
                  if (iitmp1 .eq. 0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               elseif (bc_part(3) .eq. 0) then
                  if (iitmp2 .eq. 0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               elseif (bc_part(5) .eq. 0) then
                  if (iitmp3 .eq. 0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               endif
            elseif (ibctype .eq. 3) then        ! no sides periodic 
               nfptsgp=nfptsgp-1
               goto 1511
            endif
            endif ! end if(nid.eq. ...)

c           take care of non-periodic stuff second
            if (ibctype .gt. 0) then
               if (ibctype .eq. 3) then         ! no sides periodic
                  if (iitmp1+iitmp2+iitmp3 .gt.0) then
                     nfptsgp=nfptsgp-1
                     goto 1511
                  endif
               elseif (ibctype .eq.1) then      ! two sides periodic
                  if (abs(bc_part(1)) .eq. 1) then
                     if (iitmp1 .gt. 0) then
                        nfptsgp=nfptsgp-1
                        goto 1511
                     endif
                  elseif (abs(bc_part(3)) .eq. 1) then
                     if (iitmp2 .gt. 0) then
                        nfptsgp=nfptsgp-1
                        goto 1511
                     endif
                  elseif (abs(bc_part(5)) .eq. 1) then
                     if (iitmp3 .gt. 0) then
                        nfptsgp=nfptsgp-1
                        goto 1511
                     endif
                  endif
               elseif (ibctype .eq.2) then      ! one side periodic
                  if (bc_part(1) .eq. 0) then
                     if (iitmp2+iitmp3.gt.0) then
                        nfptsgp=nfptsgp-1
                        goto 1511
                     endif
                  elseif (bc_part(3) .eq. 0) then
                     if (iitmp1+iitmp3.gt.0) then
                        nfptsgp=nfptsgp-1
                        goto 1511
                     endif
                  elseif (bc_part(5) .eq. 0) then
                     if (iitmp1+iitmp2.gt.0) then
                        nfptsgp=nfptsgp-1
                        goto 1511
                     endif
                  endif
               endif
            endif
 1511 continue
         endif
         endif
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_neighbor_el_proc
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'
c
c     This routine is called once at the beginning of the particle
c     simulation. At the end of this routine, the common blocks
c     /neighbor_proc/ & /neighbor_el_number/ are set. The idea behind
c     this routine is to know what processor owns neighboring spectral
c     elements and what local element number the neighboring element is.
c
c     el_*_proc_map holds: *(face,edge,corner) neighboring element 
c                           MPI rank number
c     el_*_el_map holds:   *(face,edge,corner) neighboring element
c                           local numbers
c
c     The ordering of faces, edges, and corners are given in el_*_num
c
c     el_*_proc_map(i,j) and el_*_el_map(i,j) are ordered by elements 
c     1 <= i <= nelt, and 1 <= j <= 26, where j=1,nfacegp are element
c     faces, j=nfacegp+1,nfacegp+nedgegp are element edges, and 
c     j = nfacegp+nedgegp+1,nfacegp+nedgegp+ncornergp are corners

      real   xerange(2,3,lelt)
      common /elementrange/ xerange
      real   xdrange(2,3)
      common /domainrange/ xdrange

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      real  rimp(7,lelt*26)
      integer iimp(3,lelt*26)

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/
     >                 0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                 0,-1,1,  1,0,1,  0,1,1,  -1,0,1,
     >                 -1,-1,0, 1,-1,0, 1,1,0,  -1,1,0
     >              /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1
     >              /)
      nfacegp   = 6  ! number of faces
      nedgegp   = 12 ! number of edges
      ncornergp = 8  ! number of corners
      idum      = 0  ! dummy arguement

      icount = 0 
      do i=1,nelt
         do j=1,3*nfacegp-2,3   ! faces
            icount = icount + 1

            xlen = xerange(2,1,i)-xerange(1,1,i)
            ylen = xerange(2,2,i)-xerange(1,2,i)
            zlen = xerange(2,3,i)-xerange(1,3,i)

            xmid = (xerange(2,1,i)+xerange(1,1,i))/2.0
            ymid = (xerange(2,2,i)+xerange(1,2,i))/2.0
            zmid = (xerange(2,3,i)+xerange(1,3,i))/2.0

            isignxx = el_face_num(j) 
            isignyy = el_face_num(j+1) 
            isignzz = el_face_num(j+2) 

            xloc = xmid + (xlen+1.e-3)*isignxx/2.0
            yloc = ymid + (ylen+1.e-3)*isignyy/2.0
            zloc = zmid + (zlen+1.e-3)*isignzz/2.0
            call bounds_p_check(xloc,xdrange(1,1),xdrange(2,1),idum)
            call bounds_p_check(yloc,xdrange(1,2),xdrange(2,2),idum)
            call bounds_p_check(zloc,xdrange(1,3),xdrange(2,3),idum)

            rimp(1,icount) = xloc
            rimp(2,icount) = yloc
            rimp(3,icount) = zloc

            iimp(1,icount) = nid
            iimp(2,icount) = 0
            iimp(3,icount) = i-1
         enddo
         do j=1,3*nedgegp-2,3    ! edges
            icount = icount + 1

            xlen = xerange(2,1,i)-xerange(1,1,i)
            ylen = xerange(2,2,i)-xerange(1,2,i)
            zlen = xerange(2,3,i)-xerange(1,3,i)

            xmid = (xerange(2,1,i)+xerange(1,1,i))/2.0
            ymid = (xerange(2,2,i)+xerange(1,2,i))/2.0
            zmid = (xerange(2,3,i)+xerange(1,3,i))/2.0

            isignxx = el_edge_num(j) 
            isignyy = el_edge_num(j+1) 
            isignzz = el_edge_num(j+2) 

            xloc = xmid + (xlen+1.e-6)*isignxx/2.0
            yloc = ymid + (ylen+1.e-6)*isignyy/2.0
            zloc = zmid + (zlen+1.e-6)*isignzz/2.0
            call bounds_p_check(xloc,xdrange(1,1),xdrange(2,1),idum)
            call bounds_p_check(yloc,xdrange(1,2),xdrange(2,2),idum)
            call bounds_p_check(zloc,xdrange(1,3),xdrange(2,3),idum)

            rimp(1,icount) = xloc
            rimp(2,icount) = yloc
            rimp(3,icount) = zloc

            iimp(1,icount) = nid
            iimp(2,icount) = 0
            iimp(3,icount) = i-1
         enddo
         do j=1,3*ncornergp-2,3   ! corners
            icount = icount + 1

            xlen = xerange(2,1,i)-xerange(1,1,i)
            ylen = xerange(2,2,i)-xerange(1,2,i)
            zlen = xerange(2,3,i)-xerange(1,3,i)

            xmid = (xerange(2,1,i)+xerange(1,1,i))/2.0
            ymid = (xerange(2,2,i)+xerange(1,2,i))/2.0
            zmid = (xerange(2,3,i)+xerange(1,3,i))/2.0

            isignxx = el_corner_num(j) 
            isignyy = el_corner_num(j+1) 
            isignzz = el_corner_num(j+2) 

            xloc = xmid + (xlen+1.e-6)*isignxx/2.0
            yloc = ymid + (ylen+1.e-6)*isignyy/2.0
            zloc = zmid + (zlen+1.e-6)*isignzz/2.0
            call bounds_p_check(xloc,xdrange(1,1),xdrange(2,1),idum)
            call bounds_p_check(yloc,xdrange(1,2),xdrange(2,2),idum)
            call bounds_p_check(zloc,xdrange(1,3),xdrange(2,3),idum)

            rimp(1,icount) = xloc
            rimp(2,icount) = yloc
            rimp(3,icount) = zloc

            iimp(1,icount) = nid
            iimp(2,icount) = 0
            iimp(3,icount) = i-1
         enddo
      enddo

c     get processor and local element number of neighboring elemetns
      call findpts(i_fp_hndl !  stride     !   call findpts( ihndl,
     $           , iimp(2,1),3        !   $             rcode,1,
     $           , iimp(1,1),3        !   &             proc,1,
     $           , iimp(3,1),3        !   &             elid,1,
     $           , rimp(5,1),7        !   &             rst,ndim,
     $           , rimp(4,1),7        !   &             dist,1,
     $           , rimp(1,1),7        !   &             pts(    1),1,
     $           , rimp(2,1),7        !   &             pts(  n+1),1,
     $           , rimp(3,1),7 ,icount)    !   &             pts(2*n+1),1,n)

c     set common block values to be used later
      do i = 1,nelt
         nstride = (i-1)*(nfacegp+nedgegp+ncornergp)
         do j = 1,nfacegp
            ijloc = nstride + j
            el_face_proc_map(i,j) = iimp(1,ijloc)
            el_face_el_map(i,j) = iimp(3,ijloc)
         enddo
         nstride = nstride + nfacegp 
         do j = 1,nedgegp
            ijloc = nstride + j
            el_edge_proc_map(i,j) = iimp(1,ijloc)
            el_edge_el_map(i,j) = iimp(3,ijloc)
         enddo
         nstride = nstride + nedgegp 
         do j = 1,ncornergp
            ijloc = nstride + j
            el_corner_proc_map(i,j) = iimp(1,ijloc)
            el_corner_el_map(i,j) = iimp(3,ijloc)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bounds_p_check(xx,xl,xr,ifmove)
      include 'SIZE'
      include 'CMTPART'
c     
c     check if xx is between domain bounds of left (xl) and right (xr)
c
c     if it is outside of bounds, move periodically to other domain side
c     and set ifmove to 1 so we know if xx has been changed
c

      ifmove = 0
      if (xx .gt. xr) then
         xx = abs(xx - xr) + xl
         ifmove = 1
      endif
      if (xx .lt. xl) then
         xx = xr - abs(xx - xl) 
         ifmove = 1
      endif

      return
      end
c----------------------------------------------------------------------
