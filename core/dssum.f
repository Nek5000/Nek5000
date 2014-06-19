      subroutine setupds(gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'

      integer gs_handle
      integer vertex(1)
      integer*8 glo_num(1),ngv

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      t0 = dnekclock()

c     Global-to-local mapping for gs
      call set_vert(glo_num,ngv,nx,nel,vertex,.false.)

c     Initialize gather-scatter code
      ntot      = nx*ny*nz*nel
      call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

c     call gs_chkr(glo_num)

      t1 = dnekclock() - t0
      if (nio.eq.0) then
         write(6,1) t1,gs_handle,nx,ngv,melg
    1    format('   setupds time',1pe11.4,' seconds ',2i3,2i12)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dssum(u,nx,ny,nz)
      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'NONCON'
      include 'PARALLEL'
      include 'TSTEP'
      real u(1)

      parameter (lface=lx1*ly1)
      common /nonctmp/ uin(lface,2*ldim),uout(lface)

      ifldt = ifield
c     if (ifldt.eq.0)       ifldt = 1
      if (ifldt.eq.ifldmhd) ifldt = 1
c     write(6,*) ifldt,ifield,gsh_fld(ifldt),imesh,' ifldt'

      if (ifsync) call nekgsync()

#ifndef NOTIMER
      if (icalld.eq.0) then
         tdsmx=0.
         tdsmn=0.
      endif
      icalld=icalld+1
      etime1=dnekclock()
#endif

c
c                 T         ~  ~T  T
c     Implement QQ   :=   J Q  Q  J
c 
c 
c                  T
c     This is the J  part,  translating child data
c
c      call apply_Jt(u,nx,ny,nz,nel)
c
c
c
c                 ~ ~T
c     This is the Q Q  part
c
      call gs_op(gsh_fld(ifldt),u,1,1,0)  ! 1 ==> +
c
c
c 
c     This is the J  part,  interpolating parent solution onto child
c
c      call apply_J(u,nx,ny,nz,nel)
c
c
#ifndef NOTIMER
      timee=(dnekclock()-etime1)
      tdsum=tdsum+timee
      ndsum=icalld
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dsop(u,op,nx,ny,nz)
      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TSTEP'
      include  'CTIMER'

      real u(1)
      character*3 op
      character*10 s1,s2
c
c     o gs recognized operations:
c     
c             o "+" ==> addition.
c             o "*" ==> multiplication.
c             o "M" ==> maximum.
c             o "m" ==> minimum.
c             o "A" ==> (fabs(x)>fabs(y)) ? (x) : (y), ident=0.0.
c             o "a" ==> (fabs(x)<fabs(y)) ? (x) : (y), ident=MAX_DBL
c             o "e" ==> ((x)==0.0) ? (y) : (x),        ident=0.0.
c     
c             o note: a binary function pointer flavor exists.
c     
c     
c     o gs level:
c     
c             o level=0 ==> pure tree
c             o level>=num_nodes-1 ==> pure pairwise
c             o level = 1,...num_nodes-2 ==> mix tree/pairwise.
c      
c
      ifldt = ifield
c     if (ifldt.eq.0)       ifldt = 1
      if (ifldt.eq.ifldmhd) ifldt = 1

c     if (nio.eq.0) 
c    $   write(6,*) istep,' dsop: ',op,ifield,ifldt,gsh_fld(ifldt)

      if(ifsync) call nekgsync()

      if (op.eq.'+  ') call gs_op(gsh_fld(ifldt),u,1,1,0)
      if (op.eq.'sum') call gs_op(gsh_fld(ifldt),u,1,1,0)
      if (op.eq.'SUM') call gs_op(gsh_fld(ifldt),u,1,1,0)

      if (op.eq.'*  ') call gs_op(gsh_fld(ifldt),u,1,2,0)
      if (op.eq.'mul') call gs_op(gsh_fld(ifldt),u,1,2,0)
      if (op.eq.'MUL') call gs_op(gsh_fld(ifldt),u,1,2,0)

      if (op.eq.'m  ') call gs_op(gsh_fld(ifldt),u,1,3,0)
      if (op.eq.'min') call gs_op(gsh_fld(ifldt),u,1,3,0)
      if (op.eq.'mna') call gs_op(gsh_fld(ifldt),u,1,3,0)
      if (op.eq.'MIN') call gs_op(gsh_fld(ifldt),u,1,3,0)
      if (op.eq.'MNA') call gs_op(gsh_fld(ifldt),u,1,3,0)

      if (op.eq.'M  ') call gs_op(gsh_fld(ifldt),u,1,4,0)
      if (op.eq.'max') call gs_op(gsh_fld(ifldt),u,1,4,0)
      if (op.eq.'mxa') call gs_op(gsh_fld(ifldt),u,1,4,0)
      if (op.eq.'MAX') call gs_op(gsh_fld(ifldt),u,1,4,0)
      if (op.eq.'MXA') call gs_op(gsh_fld(ifldt),u,1,4,0)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine vec_dssum(u,v,w,nx,ny,nz)
c
c     Direct stiffness summation of the face data, for field U.
c
c     Boundary condition data corresponds to component IFIELD of 
c     the CBC array.
c
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
      include 'CTIMER'

      REAL U(1),V(1),W(1)

      if(ifsync) call nekgsync()

#ifndef NOTIMER
      if (icalld.eq.0) tvdss=0.0d0
      if (icalld.eq.0) tgsum=0.0d0
      icalld=icalld+1
      nvdss=icalld
      etime1=dnekclock()
#endif

c
c============================================================================
c     execution phase
c============================================================================
c
      ifldt = ifield
c     if (ifldt.eq.0)       ifldt = 1
      if (ifldt.eq.ifldmhd) ifldt = 1

      call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,1,0)

#ifndef NOTIMER
      timee=(dnekclock()-etime1)
      tvdss=tvdss+timee
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif

      return
      end

c-----------------------------------------------------------------------
      subroutine vec_dsop(u,v,w,nx,ny,nz,op)
c
c     Direct stiffness summation of the face data, for field U.
c
c     Boundary condition data corresponds to component IFIELD of 
c     the CBC array.
c
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
      include 'CTIMER'
c
      real u(1),v(1),w(1)
      character*3 op

c============================================================================
c     execution phase
c============================================================================

      ifldt = ifield
c     if (ifldt.eq.0)       ifldt = 1
      if (ifldt.eq.ifldmhd) ifldt = 1

c     write(6,*) 'opdsop: ',op,ifldt,ifield
      if(ifsync) call nekgsync()

      if (op.eq.'+  ' .or. op.eq.'sum' .or. op.eq.'SUM')
     $   call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,1,0)


      if (op.eq.'*  ' .or. op.eq.'mul' .or. op.eq.'MUL')
     $   call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,2,0)


      if (op.eq.'m  ' .or. op.eq.'min' .or. op.eq.'mna'
     $                .or. op.eq.'MIN' .or. op.eq.'MNA')
     $   call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,3,0)


      if (op.eq.'M  ' .or. op.eq.'max' .or. op.eq.'mxa'
     $                .or. op.eq.'MAX' .or. op.eq.'MXA')
     $   call gs_op_many(gsh_fld(ifldt),u,v,w,u,u,u,ndim,1,4,0)


      return
      end
c-----------------------------------------------------------------------
      subroutine nvec_dssum(u,stride,n,gs_handle)

c     Direct stiffness summation of the array u for n fields
c
      include 'SIZE'
      include 'CTIMER'

      real u(1)
      integer n,stride,gs_handle

      if(ifsync) call nekgsync()

#ifndef NOTIMER
      icalld=icalld+1
      nvdss=icalld
      etime1=dnekclock()
#endif
      call gs_op_fields(gs_handle,u,stride,n,1,1,0)

#ifndef NOTIMER
      timee=(dnekclock()-etime1)
      tvdss=tvdss+timee
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif

      return
      end

c----------------------------------------------------------------------
      subroutine matvec3(uout,Jmat,uin,iftrsp,n1,n2)
c
      include 'SIZE'
c
      real Jmat (n1,n1,2)
      real uin   (1)
      real uout  (1)
      logical iftrsp
c
      common /matvtmp/ utmp(lx1,ly1)
c
      if (ndim.eq.2) then
         call mxm (Jmat(1,1,1),n1,uin,n1,uout,n2)
      else
         if (iftrsp) then
            call transpose(uout,n2,uin,n1)
         else
            call copy     (uout,uin,n1*n2)
         endif
         call mxm (Jmat(1,1,1),n1,uout,n1,utmp,n2)
         call mxm (utmp,n2,Jmat(1,1,2),n1,uout,n1)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine matvec3t(uout,Jmat,uin,iftrsp,n1,n2)
c
      include 'SIZE'
c
      real Jmat (n1,n1,2)
      real uin   (n1,n2)
      real uout  (n1,n2)
      logical iftrsp
c
      common /matvtmp/ utmp(lx1*ly1)
c
      call transpose(utmp,n2,uin,n1)
      call mxm (Jmat(1,1,2),n1,utmp,n1,uout,n2)
      call mxm (uout,n2,Jmat(1,1,1),n1,utmp,n1)
      if (iftrsp) then
         call copy     (uout,utmp,n1*n2)
      else
         call transpose(uout,n2,utmp,n1)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine matvect (out,d,vec,n1,n2)
      dimension d(n1,n2),out(1),vec(1)
c
c   handle non-square matrix in mat-vec mult -- TRANSPOSE
c    N1 is still the number of rows
c    N2 is still the number of cols
c
c
      call mxm(vec,1,d,n1,out,n2)
c
      return
      end
c-----------------------------------------------------------------------
c      subroutine opq_in_place(a,b,c)
c      include 'SIZE'
c      real a(1),b(1),c(1)
c
c      call q_in_place(a)
c      call q_in_place(b)
c      if (ndim .eq.3) call q_in_place(c)
c
c      return
c      end
c-----------------------------------------------------------------------
      subroutine vectof_add(b,a,ie,iface,nx,ny,nz)
C
C     Copy vector A to the face (IFACE) of B
C     IFACE is the input in the pre-processor ordering scheme.
C
      DIMENSION A(NX,NY)
      DIMENSION B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        B(IX,IY,IZ,IE) = B(IX,IY,IZ,IE) + A(k,1)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine zero_f(b,ie,iface,nx,ny,nz)
C
C     ZERO the face (IFACE) of B
C     IFACE is the input in the pre-processor ordering scheme.
C
      DIMENSION B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
c
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        B(IX,IY,IZ,IE) = 0.
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ftovec_0(a,b,ie,iface,nx,ny,nz)
C
C     Copy the face (IFACE) of B to vector A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      DIMENSION A(NX,NY)
      DIMENSION B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        A(k,1)=B(IX,IY,IZ,IE)
        B(IX,IY,IZ,IE)=0.0
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ftovec(a,b,ie,iface,nx,ny,nz)
C
C     Copy the face (IFACE) of B to vector A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      real A(NX,NY)
      real B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        A(k,1)=B(IX,IY,IZ,IE)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vectof(b,a,ie,iface,nx,ny,nz)
C
C     Copy vector A to the face (IFACE) of B
C     IFACE is the input in the pre-processor ordering scheme.
C
      real A(NX,NY)
      real B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        B(IX,IY,IZ,IE) = A(k,1)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ftoveci(a,b,ie,iface,nx,ny,nz)
C
C     Copy the face (IFACE) of B to vector A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      integer A(NX,NY)
      integer B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        A(k,1)=B(IX,IY,IZ,IE)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vectofi(b,a,ie,iface,nx,ny,nz)
C
C     Copy vector A to the face (IFACE) of B
C     IFACE is the input in the pre-processor ordering scheme.
C
      integer A(NX,NY)
      integer B(NX,NY,NZ,1)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        B(IX,IY,IZ,IE) = A(k,1)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine apply_Jt(u,nx,ny,nz,nel)
      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'NONCON'
      include 'PARALLEL'
      include 'TSTEP'
      real u(1)
c
      parameter (lface=lx1*ly1)
      common /nonctmp/ uin(lface,2*ldim),uout(lface)
c
c 
c                  T
c     This is the J  part,  translating child data
c
      do ie = 1 , nel
c        Note, we zero out u() on this face after extracting, for
c        consistency reasons discovered during Jerry's thesis. 
c        Thus,  "ftovec_0" rather than ftovec().   (iface -- Ed notation)
         do iface = 1 , 2*ndim
            im = mortar(iface,ie)
            if (im.ne.0) then
               call ftovec_0(uin(1,iface),u,ie,iface,nx,ny,nz)
            endif
         enddo
         do iface=1,2*ndim
            im = mortar(iface,ie)
            if (im.ne.0) then
               if (if3d) then
                 call matvec3t
     $               (uout,Jmat(1,1,1,im),uin(1,iface),ifJt(im),nx,nx)
               else
                 call matvect (uout,Jmat(1,1,1,im),uin(1,iface),nx,nx)
               endif
               call vectof_add(u,uout,ie,iface,nx,ny,nz)
            endif
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine apply_J(u,nx,ny,nz,nel)
      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'NONCON'
      include 'PARALLEL'
      include 'TSTEP'
      real u(1)
c
      parameter (lface=lx1*ly1)
      common /nonctmp/ uin(lface,2*ldim),uout(lface)
c
c     This is the J  part,  interpolating parent solution onto child
c
c
      do ie = 1 , nel
         do iface = 1 , 2*ndim
            im = mortar(iface,ie)
            if (im.ne.0) then
               call ftovec(uin(1,iface),u,ie,iface,nx,ny,nz)
            endif
         enddo
         do iface=1,2*ndim
            im = mortar(iface,ie)
            if (im.ne.0) then
               call matvec3
     $            (uout,Jmat(1,1,1,im),uin(1,iface),ifJt(im),nx,nz)
               call vectof (u,uout,ie,iface,nx,ny,nz)
            endif
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine h1_proj(u,nx,ny,nz)
      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'NONCON'
      include 'PARALLEL'
      include 'TSTEP'
      real u(1)
c
      parameter (lface=lx1*ly1)
      common /nonctmp/ uin(lface,2*ldim),uout(lface)

      if(ifsync) call nekgsync()

#ifndef NOTIMER
      if (icalld.eq.0) then
         tdsmx=0.
         tdsmn=0.
      endif
      icalld=icalld+1
      etime1=dnekclock()
#endif
c
      ifldt = ifield
c     if (ifldt.eq.0) ifldt = 1
      nel = nelv
      if (ifield.ge.2) nel=nelt
      ntot = nx*ny*nz*nel


c
c
c                        ~  ~T  
c     Implement   :=   J Q  Q  Mu
c 
c 
c                  T
c
      call col2  (u,umult,ntot)
c
c                 ~ ~T
c     This is the Q Q  part
c
      call gs_op(gsh_fld(ifldt),u,1,1,0) 
c
c 
c     This is the J  part,  interpolating parent solution onto child
c
      call apply_J(u,nx,ny,nz,nel)
c
c
#ifndef NOTIMER
      timee=(dnekclock()-etime1)
      tdsum=tdsum+timee
      ndsum=icalld
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dssum_msk(u,mask,nx,ny,nz)
      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'NONCON'
      include 'PARALLEL'
      include 'TSTEP'
      real u(1),mask(1)
c
      parameter (lface=lx1*ly1)
      common /nonctmp/ uin(lface,2*ldim),uout(lface)

      if(ifsync) call nekgsync()

#ifndef NOTIMER
      if (icalld.eq.0) then
         tdsmx=0.
         tdsmn=0.
      endif
      icalld=icalld+1
      etime1=dnekclock()
#endif
c
      ifldt = ifield
c     if (ifldt.eq.0) ifldt = 1
      nel = nelv
      if (ifield.ge.2) nel=nelt
      ntot = nx*ny*nz*nel


c
c                    T           ~  ~T  T
c     Implement Q M Q   :=   J M Q  Q  J
c 
c 
c                  T
c     This is the J  part,  translating child data
c
      call apply_Jt(u,nx,ny,nz,nel)
c
c
c
c                 ~ ~T
c     This is the Q Q  part
c
      call gs_op(gsh_fld(ifldt),u,1,1,0) 
      call col2  (u,mask,ntot)
c
c 
c     This is the J  part,  interpolating parent solution onto child
c
      call apply_J(u,nx,ny,nz,nel)
c
c
#ifndef NOTIMER
      timee=(dnekclock()-etime1)
      tdsum=tdsum+timee
      ndsum=icalld
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dssum_msk2(u,mask,binv,nx,ny,nz)
      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'NONCON'
      include 'PARALLEL'
      include 'TSTEP'
      real u(1),mask(1),binv(1)
c
      parameter (lface=lx1*ly1)
      common /nonctmp/ uin(lface,2*ldim),uout(lface)

      if(ifsync) call nekgsync()

#ifndef NOTIMER
      if (icalld.eq.0) then
         tdsmx=0.
         tdsmn=0.
      endif
      icalld=icalld+1
      etime1=dnekclock()
#endif

c
      ifldt = ifield
c     if (ifldt.eq.0) ifldt = 1
      nel = nelv
      if (ifield.ge.2) nel=nelt
      ntot = nx*ny*nz*nel


c
c
c                    T           ~  ~T  T
c     Implement Q M Q   :=   J M Q  Q  J
c 
c 
c                  T
c     This is the J  part,  translating child data
c
      call apply_Jt(u,nx,ny,nz,nel)
c
c
c
c                 ~ ~T
c     This is the Q Q  part
c
      call gs_op(gsh_fld(ifldt),u,1,1,0) 
      call col3  (u,mask,binv,ntot)
c
c 
c     This is the J  part,  interpolating parent solution onto child
c
      call apply_J(u,nx,ny,nz,nel)
c
c
#ifndef NOTIMER
      timee=(dnekclock()-etime1)
      tdsum=tdsum+timee
      ndsum=icalld
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif
c
      return
      end
c-----------------------------------------------------------------------
