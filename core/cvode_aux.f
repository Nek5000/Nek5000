      subroutine add_fcvfun_usr(ydot,jstart)
c                                               
c     plug-in specific contributions to rhs     
c                                               
      include 'SIZE'                            
      include 'TOTAL'                          
                                               
      common /scrns/ sumQ(lx1,ly1,lz1,lelt)     
     $              ,sumV(lx1,ly1,lz1,lelt)    
     $              ,cpwt(lx1,ly1,lz1,lelt)    
     $              ,sumT(lx1,ly1,lz1,lelt)    
     $              ,dtmp(lx1,ly1,lz1,lelt)   
      real ydot(1)
                                                
      if(.not.ifvarp(1)) return                 
      ntotv = nx1*ny1*nz1*nelv                

      if (.not. ifvcor) then

         dp0thdt = 0.0

      else

         dd = gamma0   ! CVref/CPref ! Note CVref denotes the inverse CPref  
         dd = (dd - 1.)/dd                                           
         xfacr= dd * dp0thdt                        
         call rzero (dtmp,ntotv)                   
         call cadd(dtmp,xfacr,ntotv)                
         call invcol2(dtmp,vtrans(1,1,1,1,2),ntotv) 
         call add2(ydot,dtmp,ntotv)                

      endif

c      if(nid.eq.0) write(*,*) 'this is dpdt and jstart', jstart, dp0thdt
      ydot(jstart) = dp0thdt                          

      return
      end
c----------------------------------------------------------------------
      subroutine cv_unpack_sol(y)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         call copy(t(1,1,1,1,i-1),y(j),ntot)
         j = j + ntot
      enddo

      p0th = y(j)
c      if(nid.eq.0) write(*,*) 'this is p0th and j', j, p0th

      return
      end
c----------------------------------------------------------------------
      subroutine cv_pack_sol(y)
c
c     copy the internal nek array (t) to the cvode solution (y)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         call copy (y(j),t(1,1,1,1,i-1),ntot)
         j = j + ntot
      enddo

      y(j) = p0th
c      if(nid.eq.0) write(*,*) 'this is y(j) and j', j, y(j)

      return
      end
c----------------------------------------------------------------------
      subroutine make_p0th
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      common /scrns/ sumQ(lx1,ly1,lz1,lelt)
     $              ,sumV(lx1,ly1,lz1,lelt)
     $              ,cpwt(lx1,ly1,lz1,lelt)
     $              ,sumT(lx1,ly1,lz1,lelt)    
     $              ,dtmp(lx1,ly1,lz1,lelt)

      nxyz = nx1*ny1*nz1
      ntotv = nxyz*nelv

      if (.not. ifvcor) then
         dp0thdt = 0.0
         return
      endif

      termQ = 0.0
      termV = 0.0
      call rzero(sumQ,ntotv)
      call rzero(sumV,ntotv)
      call rzero(cpwt,ntotv)

      call rone(dtmp,ntotv)
      dd = gamma0 ! CVref/CPref ! Note CVref denotes the inverse CPref
      dd = -1.0*(dd - 1.)/dd
      call cmult(dtmp,dd,ntotv)
c  here we should add the variable molecular weight and cp factor
c      call invcol3(dtmp,vtrans,wave,ntotv)
c      call invcol2(dtmp,vtrans(1,1,1,1,2),ntotv)
      call cadd(dtmp,1.0,ntotv)
      call col2(dtmp,bm1,ntotv)
      p0alph1 = p0th / glsum(dtmp,ntotv)

      call copy   (sumQ,QTL,ntotv)

      call col2   (sumQ,   bm1,ntotv)
      termQ   = glsum(sumQ,ntotv)
      call volume_change(termV)
c      if(nid.eq.0) write(*,*) 'dV/dt term and div is', termV, termQ

      dp0thdt = p0alph1*(termQ - termV)

c      if(nid.eq.0) write(*,*) 'ydot for p0 is', dp0thdt

      return
      end

c----------------------------------------------------------------------
      subroutine volume_change(termV)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      character cb*3

      nxyz1= nx1*ny1*nz1
      ntot1= nxyz1*nelv
      nfaces = 2*ndim

      termA = 0.0
      termVL= 0.0

      do 100 iel=1,nelv
      do 100 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.'v  ' .or. cb.eq.'V  ' .or. cb.eq.'mv ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
            ia = 0
            do 10 iz=kz1,kz2
            do 10 iy=ky1,ky2
            do 10 ix=kx1,kx2
               ia =ia + 1
               termxyz = vx(ix,iy,iz,iel)*unx(ia,1,iface,iel)
     $                 + vy(ix,iy,iz,iel)*uny(ia,1,iface,iel)
     $                 + vz(ix,iy,iz,iel)*unz(ia,1,iface,iel)
               termA  = termA + area(ia,1,iface,iel)
               termVL = termVL+ termxyz * area(ia,1,iface,iel)
 10         continue
         endif
 100  continue

      termV= glsum(termVL,1)  ! sum across processors

      return
      end

      subroutine add_qthermal
c
c    add pressure term to thermal divergence for closed systems: 
c
c    qtl := qtl(standard) + (1-(g0-1)/g0/W/Cp) dp0th/dt

      include 'SIZE'
      include 'TOTAL'
      
      common /scrns/ w1(lx1,ly1,lz1,lelt)
     $              ,w2(lx1,ly1,lz1,lelt)
     $              ,tx(lx1,ly1,lz1,lelt)
     $              ,ty(lx1,ly1,lz1,lelt)
     $              ,tz(lx1,ly1,lz1,lelt)

      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelv

         call rone(w2,ntot)
         dd = gamma0 ! CVref/CPref ! Note CVref denotes the inverse CPref
         dd = -1.0*(dd - 1.)/dd
         call cmult(w2,dd,ntot)
c  here we should add the variable molecular weight and cp factor
c      call invcol3(dtmp,vtrans,wave,ntotv)
c      call invcol2(dtmp,vtrans(1,1,1,1,2),ntotv)

         call cadd(w2,1.0,ntot)
         dd =-dp0thdt/p0th
c         if(nid.eq.0) write(6,*) 'add_qthermal ifield, dd=', ifield, dd
         call cmult(w2,dd,ntot)
         call add2 (QTL,w2,ntot)

      return
      end
