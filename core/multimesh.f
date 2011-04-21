c-----------------------------------------------------------------------
      subroutine get_session_info(intracomm)
      include 'mpif.h'
      include 'SIZE'
      include 'GLOBALCOM' 
      include 'TSTEP' 
      include 'INPUT' 
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal
      integer nid_global_root(0:nsessmax-1)
      character*132 SESSION_MULT(0:nsessmax-1), PATH_MULT(0:nsessmax-1)
      logical high

C     nsessmax = upper limit for number of sessions
C     npmax = upper limit for number of processors in each session 


C     Read from a file: number of sessions (nsessions), 
C     session name (SESSION_MULT(n-1)),
C     number of pocessors in each session (npsess(n-1)) 
C     and path (PATH_MULT(n-1))


      call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
      if ( mpi_is_initialized .eq. 0 ) call mpi_init (ierr)

      call mpi_comm_size(mpi_comm_world,np_global,ierr)
      call mpi_comm_rank(mpi_comm_world,nid_global,ierr)

      nid=nid_global
      nekcomm=mpi_comm_world

      ierr = 0
      IF(NID_GLOBAL.EQ.0) THEN
      OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD',ERR=24)
      read(8,*) nsessions
      do n=0,nsessions-1
         call blank(session_mult(n),132)
         call blank(path_mult(n)   ,132)
         read(8,10) session_mult(n)
         read(8,10) path_mult(n)
         read(8,*) npsess(n)
      end do
 10   FORMAT(A132)
      CLOSE(UNIT=8)
      GOTO 23
 24   ierr = 1
 23   END IF
      
      call err_chk(ierr,' Cannot open SESSION.NAME!$')

      call bcast(nsessions,4)

      if (nsessions.gt.2) 
     &  call exitti('More than 2 sessions are not supported!$',1)

      do n=0,nsessions-1
      call bcast(npsess(n),4)
      call bcast(session_mult(n),132)
      call bcast(path_mult(n),132)
      end do

      npall=0
      do n=0,nsessions-1
         npall=npall+npsess(n)
      end do
     
C     Check if number of processors in each session is consistent 
C     with the total number of processors

      if (npall.ne.np_global) 
     &  call exitti('Wrong number of processors!$',1)

C     Assign key for splitting into multiple groups

      nid_global_root_next=0
      do n=0,nsessions-1
         nid_global_root(n)=nid_global_root_next
         nid_global_root_next=nid_global_root(n)+npsess(n)
         if (nid_global.ge.nid_global_root(n).and.
     &   nid_global.lt.nid_global_root_next) 
     &    idsess=n
      end do
         
      call mpi_comm_split(mpi_comm_world,idsess,nid,intracomm,ierr)
 
      SESSION = SESSION_MULT(idsess)
      PATH = PATH_MULT(idsess)

C     Intercommunications set up only for 2 sessions 

      if (nsessions.gt.1) then

      if (idsess.eq.0) idsess_neighbor=1
      if (idsess.eq.1) idsess_neighbor=0
 
       call mpi_intercomm_create(intracomm,0,mpi_comm_world, 
     &  nid_global_root(idsess_neighbor), 10,intercomm,ierr)

      np_neighbor=npsess(idsess_neighbor)
      
      call iniproc(intracomm)

      high=.true.
      call mpi_intercomm_merge(intercomm, high, iglobalcomm, ierr)

 11   format(3i20,a7)
      
      call mpi_allreduce(nid_global,irecv,1,mpi_integer,mpi_min,
     & iglobalcomm,ierr)

      IFNEKNEK  = .true.

c     Initialize order of interface extrapolation to 1. Only used with NEKNEK.
      NINTER = 1

      end if 

      return
      end
C-----------------------------------------------------------------------
      subroutine multimesh_create

      include 'SIZE'
      include 'TOTAL'
      integer iList_all(ldim+2,nmaxcom)
      real    pts(ldim,nmaxcom)

      integer icalld
      save    icalld
      data    icalld  /0/

C     Set interpolation flag: points with boundary condition = 'int' get intflag=1. 
C     Boundary conditions are changed back to 'v' or 't'.

      if (icalld.eq.0) then
         call set_intflag
         icalld = icalld + 1
      end if 

C     Every processor sends all the points with intflag=1 to the processors of remote session
C     (send_points)

C    Root processor receives all the boundary points from remote session
C    and redistributes it among its processors for efficient 
C    localization of points (recv_points)

      call exchange_points(pts,iList_all,npoints_all)
      
C     Find which boundary points of remote session are
C     within the local mesh (intpts_locate) 

      call intpts_locate (pts,iList_all,npoints_all)

C     Sort which pocessors of remote session owns each point
C     and communicate the point info(iList to this processor
C     Points which are inside the mesh are masked and
C     their identities are stored in array (point_sort)

      call points_sort


      return

      end
C-------------------------------------------------------------
      subroutine set_intflag 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      CHARACTER CB*3

C     Set interpolation flag: points with boundary condition = 'int' get intflag=1. 
C     Boundary conditions are changed back to 'v' or 't'.

      if (ifflow) then
         ifield = 1
      elseif (ifheat) then
         ifield   = 2
      endif   

      nfaces = 2*ndim
      nel    = nelfld(ifield)
      
      nflag=nel*nfaces
      call izero(intflag,nflag)

      do 2010 iel=1,nel
      do 2010 iface=1,nfaces
         cb=cbc(iface,iel,ifield)
         if (cb.eq.'int') then
            intflag(iface,iel)=1
         if (ifield.eq.1) cbc(iface,iel,ifield)='v'
         if (ifield.eq.2) cbc(iface,iel,ifield)='t'
         end if
 2010 continue

      return
      end

c-------------------------------------------------------------
      subroutine exchange_points(pts,iList_all,npoints_all)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      include 'mpif.h'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      real    pts(ldim,nmaxcom)
 
      real rsend(ldim*nmaxl), rrecv(ldim*nmaxl)     
      integer isend((ldim+2)*nmaxl),irecv((ldim+2)*nmaxl)
      integer nsend(0:lp-1), nrecv(0:lp-1), iList_all(ldim+2,nmaxcom)
    
C     Look for boundary points with Diriclet b.c. (candidates for
C     interpolation)

      
      if (ifflow) then
         ifield = 1
      elseif (ifheat) then
         ifield   = 2
      endif   

      nfaces = 2*ndim
      nel    = nelfld(ifield)

      ip = 0
      do 2010 iel=1,nel
      do 2010 iface=1,nfaces
         if (intflag(iface,iel).eq.1) then
            call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
            do 100 iz=kz1,kz2
            do 100 iy=ky1,ky2
            do 100 ix=kx1,kx2
               call nekasgn (ix,iy,iz,iel)
               ip=ip+1

               if (if3d) then
               isend((ldim+2)*(ip-1)+1)=ix
               isend((ldim+2)*(ip-1)+2)=iy
               isend((ldim+2)*(ip-1)+3)=iz
               isend((ldim+2)*(ip-1)+4)=iel
               isend((ldim+2)*(ip-1)+5)=nid

               rsend(ldim*(ip-1)+1)=x
               rsend(ldim*(ip-1)+2)=y
               rsend(ldim*(ip-1)+3)=z
               else
               isend((ldim+2)*(ip-1)+1)=ix
               isend((ldim+2)*(ip-1)+2)=iy
               isend((ldim+2)*(ip-1)+3)=iel
               isend((ldim+2)*(ip-1)+4)=nid

               rsend(ldim*(ip-1)+1)=x
               rsend(ldim*(ip-1)+2)=y
               end if


               if (ip.gt.nmaxl) then
                  write(6,*) nid,
     &            ' ABORT: nbp (current ip) too large',ip,nmaxl
                  call exitt
               endif

  100       continue
         endif
 2010 continue

      nbp = ip

      call izero(iList_all,(ldim+2)*nmaxcom)
      call rzero(pts,  ldim*nmaxcom)

C     Scatter boundary points, its coordinates and its local identities 
C     to the the processors of remote session (for load balancing)  

C     Determine amount of points to send to each processor

      ndistrib = nbp/np_neighbor
      iremainder = nbp - ndistrib*np_neighbor

      do id=0,np_neighbor-1
        if(id.le.iremainder-1) then
            nsend(id) = ndistrib + 1
            else 
            nsend(id) = ndistrib
         end if
         len=isize
         call mpi_irecv (nrecv(id),len,mpi_byte, id, id, 
     &                                    intercomm, msg,ierr)
         len=isize 
         call mpi_send  (nsend(id),len,mpi_byte, id, nid, 
     &                                    intercomm, ierr)

         call mpi_wait (msg,status,ierr)

      enddo

      call mpi_barrier(intercomm,ierr)


      iaddress=1

      ip=0

      do id=0,np_neighbor-1
         len=(ldim+2)*nrecv(id)*isize
         call mpi_irecv (irecv, len ,mpi_byte, id, 100+id, 
     &                                     intercomm, msg, ierr)
         len=(ldim+2)*nsend(id)*isize
         call mpi_send (isend(iaddress),len,mpi_byte, id, 100+nid, 
     &                                           intercomm, ierr)
         iaddress=iaddress+(ldim+2)*nsend(id)

         call mpi_wait (msg, status, ierr)

         do n=1,nrecv(id)
            ip = ip + 1

            if (ip.gt.nmaxcom) then
               write(6,*) nid,'ABORT: increase nmaxcom',ip,nmaxcom
               call exitt
            endif

            
           do j=1,ldim+2
              iList_all(j,ip)=irecv((ldim+2)*(n-1)+j)
           enddo

         enddo ! n

      enddo    ! id

      call mpi_barrier(intercomm,ierr)

      npoints_all = ip

      iaddress=1

      ip=0

      do id=0,np_neighbor-1
         len=ldim*nrecv(id)*wdsize
         call mpi_irecv (rrecv, len ,mpi_byte, id, 200+id, 
     &                                  intercomm, msg, ierr)
         len=ldim*nsend(id)*wdsize
         call mpi_send (rsend(iaddress),len,mpi_byte, id, 200+nid, 
     &                                        intercomm, ierr)
         iaddress=iaddress+ldim*nsend(id)

         call mpi_wait (msg, status, ierr)

         do n=1,nrecv(id)
            ip = ip + 1

            if (ip.gt.nmaxcom) then
               write(6,*) 'ABORT: increase nmaxcom' 
               call exitt
            endif
            
               do j=1,ldim
               pts(j,ip)=rrecv(ldim*(n-1)+j)
               end do

         enddo ! n

      enddo ! id

      call mpi_barrier(intercomm,ierr)

      npoints_all = ip    

      return
      end
C-------------------------------------------------------------------------
      subroutine intpts_locate
     $ (pts,iList_all,npoints_all)

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'NEKNEK'

      real    pts(ldim,nmaxcom)
      real    dist_all(nmaxcom)
      real    rst_all(nmaxcom*ldim)
      integer rcode_all(nmaxcom),elid_all(nmaxcom),proc_all(nmaxcom)
      real    dist(nmaxcom)
      integer iList_all(ldim+2,nmaxcom)
      character*3 CB

        call intpts_setup(-1.0,inth_multi) ! use default tolerance
        

        call findpts(inth_multi,rcode_all,1,
     &                 proc_all,1,
     &                 elid_all,1,
     &                 rst_all,ndim,
     &                 dist_all,1,
     &                 pts(1,1),ndim,
     &                 pts(2,1),ndim,
     &                 pts(3,1),ndim,npoints_all)


c     Rearrange arrays so that only the points which are found within 
c     the mesh are marked for interpolation: those are truly internal 
c     points (rcode_all=0) plus some of the "boundary points" 
c     (rcode_all=1) which fall on the boundaries of the elements 
c     of another domain (essentially also internal points) or 
c     wall or periodic boundaries 
c
c     Points with rcode_all(i)=1 which fall on the "interpolation" 
c     boundaries of another domain are treated with boundary 
c     conditions and not with interpolation routines, and are 
c     therefore excluded.
c
c     rcode_all, proc_all ... npoints_all represent all boundary points,
c     rcode, proc ... npoints - only interpolation points 
c     (used by interpolation routines)
c     Contained in a common block /multipts_r/, /multipts_i/ (in NEKNEK)

C

C     Exclude "true" boundary points 

      if (ifflow) then
         ifield = 1
      elseif (ifheat) then
         ifield   = 2
      endif   

      ip=0

      do 100 i=1,npoints_all
         
          if (rcode_all(i).lt.2) then

            iel  = elid_all(i) + 1

         if (rcode_all(i).eq.1.and.dist_all(i).gt.1e-12) then
            if (ndim.eq.2) 
     &      write(6,'(A,3E15.7)') 
     &   ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',
     &   (pts(k,i),k=1,ndim),dist_all(i)  
            if (ndim.eq.3) 
     &      write(6,'(A,4E15.7)') 
     &   ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',
     &   (pts(k,i),k=1,ndim),dist_all(i)  
            GOTO 100
         end if
            ip=ip+1
            rcode(ip) = rcode_all(i)
            elid(ip)  = elid_all(i)
            proc(ip)  = proc_all(i)
            do j=1,ndim
            rst(ndim*(ip-1)+j)   = rst_all(ndim*(i-1)+j)
            end do
            do n=1,ndim+2
            iList(n,ip) = iList_all(n,i)
            end do   
          end if
 100   continue

      npoints=ip

      write(6,'(a7,i7,1x,a10)') 'found', npoints, session
      
      return
      end
  
C---------------------------------------------------------------------
      subroutine points_sort
      include 'SIZE'
      include 'TOTAL'  
      include 'NEKNEK' 
      include 'mpif.h'    
      integer isend((ldim+2)*nmaxl), irecv((ldim+1)*nmaxl)
      integer status(mpi_status_size)

c     Sort all the points in the iList serviced by the local 
c     processor with rank nid according to the rank of remote 
c     processor (id) who needs (owns) this particular point

c     All necessary information is contained in common block /proclist/

c     nbp_send(id), id=0...lp-1 is the number of points sent by 
c     each local processor with rank nid to the remote processor 
c     with rank id. If nbp_send(id)=0, this processor is marked 
c     as non-sending 

c     nbp_recv(id), id=0...lp-1 is the number of points received 
c     by each local processor with rank nid from the remote processor 
c     with rank id. If nbp_recv(id)=0, this processor is marked as 
c     non-receiving 

c     idp(n,id)=1...nint is the index of each point within the array 
c     of all points serviced by the local processor nid (to be sent 
c     to different remote processors) Here id is the rank of remote 
c     processor who needs this point, n=1...nbp_send(id) is the local 
c     point index among the points sent by local processor nid to 
c     remote processor id

c     iden(4,n,id) is the identity information for the points 
c     received by the local processor nid from the remote 
c     processor id=0...lp-1 to attach the interpolated value 
c     to the corresponding point owned by the local processor 
c     nid (in array valint, see subroutine get_values called 
c     at every time step). Here n=1...nbp_recv(id) is the local 
c     point index among the points received by local processor 
c     nid from remote processor id 



C     iden(1,:,:) = ix
C     iden(2,:,:) = iy
C     iden(3,:,:) = iz
C     iden(4,:,:) = iel

      call izero(nbp_send,lp)
      call izero(nbp_recv,lp)
      call izero(idp,nmaxl*lp)
      call izero(iden,(ldim+1)*nmaxl*lp)
      
      if (IFFLOW) then
         IFIELD = 1
      else if (IFHEAT) then
         IFIELD   = 2
      end if   

      NEL    = NELFLD(IFIELD)
      NXYZ  =NX1*NY1*NZ1
      NTOT  =NXYZ*NEL

      call izero(imask,ntot)


c     Remote processor id which needs each interpolated point i 
c     (i=1...npoints) serviced by the local processor nid is contained 
c     in iList(ldim+2,:).  This allows to calculate nbp_send(id) 
c     and idp(n,id), where n=1...nbp_send(id) 
c     nbp_send(id) - how many points local processor nid will send to remote processor id,
c     and idp(n,id) - index of each sent point within the interpolation points serviced by the local processor nid 

      do i=1,npoints
         id=iList(ldim+2,i)
         nbp_send(id)=nbp_send(id)+1
         n=nbp_send(id)
         idp(n,id)=i
      end do



c     Local processor nid sends the information about the interpolated 
c     points to each processor id of the remote session  
c     telling how many points to expect (nbp_send(id))and their identity contained 
c     in iList:ix,iy,iz,iel. 

c     Local processor nid receives the information about the interpolated 
c     points from each processor id of remote session (nbp_recv(id)) 
c     telling how many points to expect and their identity contained in 
c     iList:ix,iy,iz,iel.
         


      call mpi_barrier(intercomm, ierr)

! Non-blocking communicaton to receive nbp_recb(id) from each processor
! and send nbp_send(id) of each processor of the remote session 

      do id=0,np_neighbor-1  
         len=isize
         call mpi_irecv (nbp_recv(id),len,mpi_byte, id, id, 
     $                                      intercomm, msg,ierr)
    
         call mpi_send  (nbp_send(id),len,mpi_byte, id, nid, 
     &                                          intercomm, ierr)

         call mpi_wait (msg,status,ierr)
      end do


      do id=0,np_neighbor-1
         ifsend(id)=.false.
         ifrecv(id)=.false.
         if (nbp_recv(id).gt.0) ifrecv(id)=.true.
         if (nbp_send(id).gt.0) ifsend(id)=.true.
      end do

      call mpi_barrier(intercomm, ierr)

      do id=0,np_neighbor-1 
         len=(ldim+1)*nbp_recv(id)*isize
      call mpi_irecv (irecv,len,mpi_byte, id, 
     & 100+id,intercomm,msg,ierr)
      
         do n=1,nbp_send(id)
            il=idp(n,id)
            do j=1,ldim+1
               isend((ldim+1)*(n-1)+j)=iList(j,il)
            end do
         end do
         len=(ldim+1)*nbp_send(id)*isize
      call mpi_send (isend,len,mpi_byte, id, 100+nid, 
     & intercomm, ierr)

      call mpi_wait(msg,status,ierr)

C     Receiving processor nid masks which points he receives from remote processor id as interpolation points
C     using the identity information contained in array 'irecv'. Those points are masked with imask=1 
C     (imask=0 for all other points). 

      do n=1,nbp_recv(id)
         ix=irecv((ldim+1)*(n-1)+1)
         iy=irecv((ldim+1)*(n-1)+2)
         iz = 1
         if (if3d) iz=irecv((ldim+1)*(n-1)+3)
         ie=irecv((ldim+1)*(n-1)+ldim+1)
         iden(1,n,id)=ix
         iden(2,n,id)=iy
         if (if3d) iden(3,n,id)=iz
         iden(ldim+1,n,id)=ie
         imask(ix,iy,iz,ie)=1
      end do

      end do

      call mpi_barrier(intercomm, ierr)



      return
      end


C----------------------------------------------------------------
      subroutine get_values(which_field)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'mpif.h' 

      character*3 which_field(nfld)
      real field(lx1*ly1*lz1*lelt,nfldmax)
      real rsend(nfldmax*nmaxl), rrecv(nfldmax*nmaxl)
      real fieldout(nfldmax,nmaxl)

      integer status(mpi_status_size)

C     Information about communication of points is contained in common block /proclist/

C     nbp_send(id), id=0...lp-1 is the number of points sent by each local processor with rank nid to 
C     the remote processor with rank id. If nbp_send(id)=0, this processor is marked as non-sending 

C     nbp_recv(id), id=0...lp-1 is the number of points received by each local processor with rank nid from 
C     the remote processor with rank id. If nbp_recv(id)=0, this processor is marked as non-receiving 

C     idp(n,id)=1...nint is the index of each point within the array of all points serviced by the local processor nid
C     (to be sent to different remote processors)
C     Here id is the rank of remote processor who needs this point,
C     n=1...nbp_send(id) is the local point index among the points sent by local processor nid to remote processor id

C     nprocsend - total number of processors in the remote session to whom local processor nid sends the points (all id's with nbp_send(id)>0,
C     id=0...lp-1)
  
C     idsend(n), n=1...nprocsend is the rank of n'th remote processor to whom local processor nid send the points
C     (note: n is just the local running number of the receiving processors >=1, explaining that array idsend is from 1 to lp,
C     not from 0 to lp-1)   

C     nprocrecv - total number of processors in the remote session from whom local processor nid receives the points (all id's with nbp_recv(id)>0,
C     id=0...lp-1)
  
C     idrecv(n), n=1...nprocsend is the rank of n'th remote processor from whom local processor nid receives the points 
C     (note: n is just the local running number of the sending processors >=1, explaining that array idrecv is from 1 to lp,
C     not from 0 to lp-1)   

C     iden(4,n,id) is the identity information for the points received by the local processor nid from the remote processor id=0...lp-1
C     to attach the interpolated value to the corresponding point owned by the local processor nid (in array valint, see subroutine get_values
C     called at every time step). Here n=1...nbp_recv(id) is the local point index among the points received by local processor nid from remote processor id
 
C     iden(1,:,:) = ix
C     iden(2,:,:) = iy
C     iden(3,:,:) = iz
C     iden(4,:,:) = iel


C     Put field values for the field ifld=1,nfld to the working array field (:,:,:,:,nfld) used byt the findpts_value routine
C     according to the field identificator which_field(ifld) 


      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt
      do ifld=1,nfld
        if (which_field(ifld).eq.'t' ) call copy(field(1,ifld),t ,nt)
        if (which_field(ifld).eq.'vx') call copy(field(1,ifld),vx,nt)
        if (which_field(ifld).eq.'vy') call copy(field(1,ifld),vy,nt)
        if (which_field(ifld).eq.'vz') call copy(field(1,ifld),vz,nt)

C     Find interpolation values      

         call findpts_eval(inth_multi,fieldout(ifld,1),nfldmax,
     &                     rcode,1,
     &                     proc,1,
     &                     elid,1,
     &                     rst,ndim,npoints,
     &                     field(1,ifld))

       end do  


C     Send interpolation values to the corresponding processors 
C     of remote session

      call mpi_barrier(intercomm, ierr)

      do id=0, np_neighbor-1
 
         if (ifrecv(id)) then
          len=nfld*nbp_recv(id)*wdsize
          call mpi_irecv (rrecv,len,mpi_byte,id,id,intercomm,msg,ierr)
       end if


       if (ifsend(id)) then     
          do n=1,nbp_send(id)
             il=idp(n,id)
             do ifld=1,nfld
                rsend(nfld*(n-1)+ifld)=fieldout(ifld,il)
             end do
           end do
           len=nfld*nbp_send(id)*wdsize
           call mpi_send (rsend,len,mpi_byte, id, nid, intercomm, ierr)
       end if
  


      if (ifrecv(id)) then
         call mpi_wait (msg,status,ierr)
    
         do n=1,nbp_recv(id)

C           extract point identity

            ix = iden(1,n,id)
            iy = iden(2,n,id)
            iz = 1
            if (if3d) iz=iden(3,n,id)
            ie=iden(ldim+1,n,id)      

               do ifld=1,nfld
                  valint(ix,iy,iz,ie,ifld)=rrecv(nfld*(n-1)+ifld)
               end do

         end do      

      end if

      end do

      call mpi_barrier(intercomm,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine value_exchange(value, oper)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'mpif.h' 
      character*4 oper
      integer status(mpi_status_size)

C     Time step synchronization between different sessions


      len=wdsize

      if (oper.eq.'recv') call mpi_irecv (value, len ,mpi_byte, 0, 0, 
     & intercomm, msg, ierr)

      if (oper.eq.'send') then
      if (nid.eq.0) then
         do id=0, np_neighbor-1
       call mpi_send (value,len,mpi_byte, id, 0, 
     & intercomm, ierr)
         end do   
      end if   
      end if

      if (oper.eq.'recv') call mpi_wait (msg, status, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk_set_xfer
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      real l2,linf
      character*3 which_field(nfldmax)

c     nfld is the number of fields to interpolate.
c     nfld = 3 for just veliocities, nfld = 4 for velocities + temperature

      which_field(1)='vx'
      which_field(2)='vy'
      which_field(3)='vz'
      if (nfld.gt.3) which_field(4)='t'

      if (nsessions.gt.1) call get_values(which_field)

      return
      end
c------------------------------------------------------------------------
      subroutine bcopy
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'

      n    = nx1*ny1*nz1*nelt

      do k=1,nfld
         call copy(bdrylg(1,k,2),bdrylg(1,k,1),n)
         call copy(bdrylg(1,k,1),bdrylg(1,k,0),n)
         call copy(bdrylg(1,k,0),valint(1,1,1,1,k),n)
      enddo

c     Order of extrpolation is contolled by the parameter NINTER contained in NEKNEK
c     First order interface extrapolation, NINTER=1 (time lagging) is activated. It is unconditionally stable.
c     If you want to use higher-order interface extrapolation schemes, you need to increase ngeom to ngeom=3-5 for scheme to be stable

      if (NINTER.eq.1.or.istep.eq.0) then
       c0=1
       c1=0
       c2=0
       else if (NINTER.eq.2.or.istep.eq.1) then
         c0=2
         c1=-1
         c2=0
       else 
         c0=3
         c1=-3
         c2=1
      endif
     
      do k=1,nfld
      do i=1,n
         ubc(i,1,1,1,k) = 
     $      c0*bdrylg(i,k,0)+c1*bdrylg(i,k,1)+c2*bdrylg(i,k,2)
      enddo
      enddo

      return
      end
C---------------------------------------------------------------------
      subroutine setintercomm(nekcommtrue,nptrue) 
      include 'SIZE' 
      INCLUDE 'GLOBALCOM'
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal 
      
      nekcommtrue=nekcomm
      nekcomm=iglobalcomm

      nptrue=np
      np=np_global

      return
      end

c-----------------------------------------------------------------------
      subroutine unsetintercomm(nekcommtrue,nptrue)
      include 'SIZE' 
      INCLUDE 'GLOBALCOM'
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal 

      nekcomm=nekcommtrue
      np=nptrue

      return
      end
c------------------------------------------------------------------------


