C-----------------------------------------------------------------------
C
      subroutine set_rsb(nelin)
C
      include 'basics.inc'
      include 'basicsp.inc'
C
      integer rsb_dim
C
      integer n,nv
      integer ia(1),ja(1),vertex(1),color(1)
      real*8  vals(1),x_centroid(1),y_centroid(1),z_centroid(1)
C
      equivalence (ia,u)
      equivalence (ja,v)
      equivalence (color,w)
      equivalence (vals,p)
      equivalence (vertex,t)
      equivalence (x_centroid,wkv1)
      equivalence (y_centroid,wkv2)
      equivalence (z_centroid,wkv3)
C     
      if (if3d) then
         rsb_dim = 1
      else
         rsb_dim = 0
      endif
      call compute_adj(ia,ja,vals,x_centroid,y_centroid,z_centroid
     $     ,nrow,vertex,nv,nelin)
 
      call out_csr_matl(ia,ja,vals,nrow)
      n=nrow
c
c     write(21,*) (ia(k),k=1,n+1)
c     nc = ia(n+1)-ia(1)
c     write(21,*) (ja(k),k=1,nc)
c     write(21,*) (vals(k),k=1,nc)
c     write(6,*) 'quit in set_rsb',n,nc
c     call exit
c
c     call hmt_spy_csrmat(ia,ja,vals,nrow,'PER..Adj.')
C
      call izero(color,n)
C
      call begin_rsb(ia,ja,vals,x_centroid,y_centroid,z_centroid,color,n
     $     ,vertex,nv,rsb_dim)
C
      return
      end
C     
C-----------------------------------------------------------------------
C
      subroutine reset_rsb()
C
      call end_rsb()
C
      return
      end

C
C-----------------------------------------------------------------------
C
      function irsb()
C
      integer c_rsb
C
      irsb = c_rsb()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine rsb()
C
      integer c_rsb
C
      ierr = c_rsb()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine set_level()
C
      call c_set_level()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine set_seed()
C
      call c_set_seed()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine set_tolerance()
C
      call c_set_tolerance()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine set_max_repeats()
C
      call c_set_max_repeats()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine set_max_lanczos()
C
      call c_set_max_lanczos()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine prs_rsb_defaults()
C
      call c_prs_defaults()
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine rsb_drawel(i)
C
C
      integer I
C
      include 'basics.inc'
      include 'basicsp.inc'
C
      integer color(1)
C
      equivalence (color,w)
C
      call hmt_drawel_color(I,color(I))
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine proc_map(draw_flag,nelin)
C
      integer draw_flag
      integer color(1)
C
      include 'basics.inc'
      include 'basicsp.inc'
C
      equivalence (color,w)
C
      call c_extract_proc_map(color)
      CALL REFRESH
      CALL DRGRID
      if (draw_flag.eq.1) then
         call drmesh_color_nln(nelin)
      else
         call drmesh_color(nelin)
      endif
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine sep_map(draw_flag,nelin)
C  
      integer draw_flag
C
      if (draw_flag.eq.1) then
         call drmesh_color_nln(nelin)
      else
         call drmesh_color(nelin)
      endif
C
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine sep_proc_map(draw_flag,nelin)
C
      integer draw_flag
C     
      if (draw_flag.eq.1) then
         call drmesh_color_nln(nelin)
      else
         call drmesh_color(nelin)
      endif
C
      return
      end
C
c-----------------------------------------------------------------------
c
      subroutine cdrmesh_color()
c
c     Redraw preprocessor mesh with color given by integer array ic(1:nel)
c
      include 'basics.inc'
      include 'basicsp.inc'
      integer color(1)
c
      equivalence (color,w)
c
      do ie=1,nel
         call hmt_drawel_color(ie,color(ie))
      enddo
      return
      end
C
c-----------------------------------------------------------------------
c
      subroutine drmesh_color(nelin)
c
c     Redraw preprocessor mesh with color given by integer array ic(1:nel)
c
      include 'basics.inc'
      include 'basicsp.inc'
      integer color(1)
c
      equivalence (color,w)
c
      do ie=1,nelin
         call hmt_drawel_color(ie,color(ie))
      enddo
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine drmesh_color_nln(nelin)
c
c     Redraw preprocessor mesh with color given by integer array ic(1:nel)
c
      include 'basics.inc'
      include 'basicsp.inc'
      integer color(1)
c
      equivalence (color,w)
c
      do ie=1,nelin
         call hmt_drawel_color_nln(ie,color(ie))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_adj(ia,ja,vals,x_centroid,y_centroid
     $     ,z_centroid,n,vertex,nv,nelin)
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      integer n,nv
      integer ia(1),ja(1),vertex(1)
      real*8  vals(1),x_centroid(1),y_centroid(1),z_centroid(1)
c
      integer iseg,i_min,nc,k,ie,ic,kl  
      integer start_pt(1),start_el(1),ind(1),ninseg(1)
      real    dx(1),dy(1),dz(1),ww(1)
c
      logical ifseg(1)
      logical if_any_per
c
      equivalence (start_pt,rxm1)
      equivalence (start_el,rym1)
      equivalence (ind     ,rzm1)
      equivalence (dx      ,sxm1)
      equivalence (dy      ,sym1)
      equivalence (dz      ,szm1)
      equivalence (ninseg  ,txm1)
      equivalence (ifseg   ,tym1)
      equivalence (ww      ,tzm1)
c
      integer      va(1)
      equivalence (va      ,sxm1)
C
C
C     only works if each element is defined by exactly 8 pts!!!
C
      n = nelin
      nc = 2**ndim
      nv = n*nc
c
      kfld = 2
      if (ifflow) kfld=1
c
c
c     Lexicographical sort to establish adj. array:
c     Any elements sharing associated vertices share and edge.
c
      k  = 0
      do ie=1,nel
         xmin = x(ie,1)
         xmax = x(ie,1)
         ymin = y(ie,1)
         ymax = y(ie,1)
         zmin = z(ie,1)
         zmax = z(ie,1)
         x_centroid(ie) = 0.0
         y_centroid(ie) = 0.0
         z_centroid(ie) = 0.0
         kl   = k+1
         do ic=1,nc
            k=k+1
            xp(k)=x(ie,ic)
            xmin = min(x(ie,ic),xmin)
            xmax = max(x(ie,ic),xmax)
            x_centroid(ie)=x_centroid(ie)+x(ie,ic)
            yp(k)=y(ie,ic)
            ymin = min(y(ie,ic),ymin)
            ymax = max(y(ie,ic),ymax)
            y_centroid(ie)=y_centroid(ie)+y(ie,ic)
            if (if3d) then
               zp(k)=z(ie,ic)
               zmin = min(z(ie,ic),zmin)
               zmax = max(z(ie,ic),zmax)
               z_centroid(ie)=z_centroid(ie)+z(ie,ic)
            else
               zmin          =0.
               zmax          =0.
               z_centroid(ie)=0.
            endif
            vertex(k) = 9999999
            start_pt(k) = k
            start_el(k) = ie
         enddo
         x_centroid(ie)=x_centroid(ie)/nc
         y_centroid(ie)=y_centroid(ie)/nc
         z_centroid(ie)=z_centroid(ie)/nc
c
c     check for outflow bc's ... assume 4-faced/elements ... hard coded
c            
         if (cbc(1,ie,kfld).eq.'O  '.or.cbc(1,ie,kfld).eq.'o  '.or.
     $       cbc(1,ie,kfld).eq.'ON '.or.cbc(1,ie,kfld).eq.'on ') then
            vertex((ie-1)*nc+1) = -9999999
            vertex((ie-1)*nc+2) = -9999999
            vertex((ie-1)*nc+5) = -9999999
            vertex((ie-1)*nc+6) = -9999999
         endif
         if (cbc(2,ie,kfld).eq.'O  '.or.cbc(2,ie,kfld).eq.'o  '.or.
     $       cbc(2,ie,kfld).eq.'ON '.or.cbc(2,ie,kfld).eq.'on ') then
            vertex((ie-1)*nc+2) = -9999999
            vertex((ie-1)*nc+3) = -9999999
            vertex((ie-1)*nc+6) = -9999999
            vertex((ie-1)*nc+7) = -9999999
         endif
         if (cbc(3,ie,kfld).eq.'O  '.or.cbc(3,ie,kfld).eq.'o  '.or.
     $       cbc(3,ie,kfld).eq.'ON '.or.cbc(3,ie,kfld).eq.'on ') then
            vertex((ie-1)*nc+3) = -9999999
            vertex((ie-1)*nc+4) = -9999999
            vertex((ie-1)*nc+7) = -9999999
            vertex((ie-1)*nc+8) = -9999999
         endif
         if (cbc(4,ie,kfld).eq.'O  '.or.cbc(4,ie,kfld).eq.'o  '.or.
     $       cbc(4,ie,kfld).eq.'ON '.or.cbc(4,ie,kfld).eq.'on ') then
            vertex((ie-1)*nc+4) = -9999999
            vertex((ie-1)*nc+1) = -9999999
            vertex((ie-1)*nc+8) = -9999999
            vertex((ie-1)*nc+5) = -9999999
         endif
         if (if3d) then
           if (cbc(5,ie,kfld).eq.'O  '.or.cbc(5,ie,kfld).eq.'o  '.or.
     $         cbc(5,ie,kfld).eq.'ON '.or.cbc(5,ie,kfld).eq.'on ') then
            vertex((ie-1)*nc+1) = -9999999
            vertex((ie-1)*nc+2) = -9999999
            vertex((ie-1)*nc+3) = -9999999
            vertex((ie-1)*nc+4) = -9999999
           endif
           if (cbc(6,ie,kfld).eq.'O  '.or.cbc(6,ie,kfld).eq.'o  '.or.
     $         cbc(6,ie,kfld).eq.'ON '.or.cbc(6,ie,kfld).eq.'on ') then
            vertex((ie-1)*nc+5) = -9999999
            vertex((ie-1)*nc+6) = -9999999
            vertex((ie-1)*nc+7) = -9999999
            vertex((ie-1)*nc+8) = -9999999
           endif
         endif
         dxe = .001*(xmax-xmin)
         dye = .001*(ymax-ymin)
         dze = .001*(zmax-zmin)
         call cfill(dx(kl),dxe,nc)
         call cfill(dy(kl),dye,nc)
         call cfill(dz(kl),dze,nc)
      enddo
c
      if (k.ne.nv) then
         write(6,*) '# elements          = ',nel
         write(6,*) '# vertices          = ',nv
         write(6,*) 'k                   = ',k
      else
         write(6,*) '# elements          = ',nel
         write(6,*) '# vertices          = ',nv
      endif
c
c     Lexicographical sorting scheme... modifies x,y,z,dx,dy,dz
c
      call lfalse(ifseg   ,nv)
      nseg        = 1
      ifseg(1)    = .true.
      ninseg(1)   = nv
c
      do j=1,ndim
c
c       Sort within each segment
c
         ioff=1
         do iseg=1,nseg
            if  (j.eq.1) then
               call rank (xp(ioff),ind,ninseg(iseg))
            elseif  (j.eq.2) then
               call rank (yp(ioff),ind,ninseg(iseg))
            else
               call rank (zp(ioff),ind,ninseg(iseg))
            endif
            call iswap (start_pt(ioff),ww,ind,ninseg(iseg))
            call iswap (start_el(ioff),ww,ind,ninseg(iseg))
            call  swap (xp      (ioff),ww,ind,ninseg(iseg))
            call  swap (yp      (ioff),ww,ind,ninseg(iseg))
            call  swap (zp      (ioff),ww,ind,ninseg(iseg))
            call  swap (dx      (ioff),ww,ind,ninseg(iseg))
            call  swap (dy      (ioff),ww,ind,ninseg(iseg))
            call  swap (dz      (ioff),ww,ind,ninseg(iseg))
            ioff=ioff+ninseg(iseg)
         enddo
c
c       Check for jumps in current coordinate
c
         if (j.eq.1) then
            do i=2,nv
               if (abs(xp(i)-xp(i-1)).gt.min(dx(i),dx(i-1)))
     $              ifseg(i)=.true.
            enddo
         elseif (j.eq.2) then
            do i=2,nv
               if (abs(yp(i)-yp(i-1)).gt.min(dy(i),dy(i-1)))
     $              ifseg(i)=.true.
            enddo
         elseif (j.eq.3) then
            do i=2,nv
               if (abs(zp(i)-zp(i-1)).gt.min(dz(i),dz(i-1)))
     $              ifseg(i)=.true.
            enddo
         endif
c
c       Count up number of different segments
c
         nseg = 0
         do i=1,nv
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo
c
c     Data now sorted lexicographically!
c
c     To sort additional incoming data:
c
c     >
c     >      do k=1,n
c     >         t_sort(k) = t_in(init_loc(k))
c     >      enddo
c     >
c
c     Generate vertex numbering 1 ... nseg
c     Indicate vertex in outflow bc by making it negative
c
      write(6,*) '# vertices (unique) = ',nseg,' (before periodic)'
c
c---- pff change here ----------------------------------------------------
c
      ioff=0
      do iseg=1,nseg
         do j=1,ninseg(iseg)
            va(start_pt(ioff+j)) = iseg
         enddo
         ioff = ioff + ninseg(iseg)
      enddo
c
c     Now that unique vertex numbering established based upon
c     physical coincidence, collapse all vertex pairs joined by
c     *periodic* bc coincidence.
c
c     call period_coinc(va,start_pt,nv)
c
      call period_bc_check(if_any_per)
      if (if_any_per) then
         call prs('Check periodicity? (y/n)$')
         call res(ans,1)
         if (ans.eq.'y'.or.ans.eq.'Y') call period_coinc(va,start_pt,nv)
      endif
c
c
c     Re-Sort entire list by *unique* vertex numbers
c
      l = 0
      do ie=1,nel
      do ic=1,nc
         l = l+1
         start_pt(l) = l
         start_el(l) = ie
      enddo
      enddo
c
      write(6,*) 'call irank 3d:',nv
      call irank (va,ind,nv)
      call iswap (va      ,ww,ind,nv)
      call iswap (start_pt,ww,ind,nv)
      call iswap (start_el,ww,ind,nv)
c
c     re-Count number of segments
c
      call lfalse(ifseg   ,nv)
      ifseg(1) =.true.
      do i=2,nv
         if (va(i).ne.va(i-1)) ifseg(i)=.true.
      enddo
c
c     Count up number of different segments = num unique vertices
c
      nseg = 0
      do i=1,nv
         if (ifseg(i)) then
            nseg = nseg+1
            ninseg(nseg) = 1
         else
            ninseg(nseg) = ninseg(nseg) + 1
         endif
      enddo
c---- to here --------------------------------------------------------
c
      write(6,*) '# vertices (unique) = ',nseg,' (after  periodic)'
      ioff=0
      ioutflow=0
      do iseg=1,nseg
         i_min = iseg
c
c        this loop finds if any element vertex has an outflow bc
         do j=1,ninseg(iseg)
            i_min = min(i_min,vertex(start_pt(ioff+j)))
         enddo
         if ((i_min.ne.iseg).and.(i_min.ne.-9999999)) then
            write(6,*) 'compute_adj() :: vertex map error!'
            call exit(0)
         endif
         if (i_min.lt.0) then
            i_min = -iseg
            ioutflow=ioutflow+1
         endif
c
c        This loop sets all shared vertices to "min" vertex value
         do j=1,ninseg(iseg)
            vertex(start_pt(ioff+j)) = i_min
         enddo
         ioff = ioff + ninseg(iseg)
      enddo
      vertex(nv+1) = nseg
      vertex(nv+2) = ioutflow
c
c     Now, sort by segment, element #s associated with each vertex
c
      ioff=1
      do iseg=1,nseg
         call irank (start_el(ioff)   ,ind,ninseg(iseg))
         call iswap (start_pt(ioff),ww,ind,ninseg(iseg))
         call iswap (start_el(ioff),ww,ind,ninseg(iseg))
         ioff=ioff+ninseg(iseg)
      enddo
c
c     Now generate upper 1/2 of adjacency graph: (i,j) format.
c
      k   =0
      ioff=0
      do iseg=1,nseg
         do i=1,ninseg(iseg)
            ie=ioff+i
c           do j=i+1,ninseg(iseg)      (for upper 1/2 only)
            do j=1,ninseg(iseg)
               je = ioff+j
               k  = k+1
               ia(k) = start_el(ie)
               ja(k) = start_el(je)
            enddo
         enddo
         ioff=ioff+ninseg(iseg)
      enddo
      ne_tmp = k
c
c     Now sort by row-number (i)
c
      call irank(ia,ind,ne_tmp)
      call iswap(ia,ww,ind,ne_tmp)
      call iswap(ja,ww,ind,ne_tmp)
c
c     Now get unique/sorted list of column numbers c
c     and overwrite ia,ja into CSR format.
c
      ie_l     = 1
      nrow     = 0
      ncol_tot = 0
      ia(ne_tmp+1) = 0
      do ie=2,ne_tmp+1
         row_sum = 0.0
         index = 0
         if (ia(ie).ne.ia(ie_l)) then
            n_in_col = ie-ie_l
C           write(6,*) n_in_col
            j1 = ie_l+1
            j2 = ie_l+n_in_col-1
c
            call irank(ja(ie_l)   ,ind,n_in_col)
            call iswap(ja(ie_l),ww,ind,n_in_col)
c
c           Overwrite ia/ja
c
c           At least one nonzero...
            n_in_col = 1
            i_ct=1
            ncol_tot     = ncol_tot+1
            ja(ncol_tot) = ja(ie_l)
c
            do jj=j1,j2
               if (ja(jj).ne.ja(jj-1)) then
                  n_in_col     = n_in_col+1
                  if (i_ct.eq.8) then
                     vals(ncol_tot) = 0.0
                  else
                     vals(ncol_tot)   = -1.0
                     row_sum = row_sum - 1.0
                  endif
                  if (ja(jj-1).eq.(nrow+1)) then
                     index = ncol_tot
                  endif
                  ncol_tot     = ncol_tot+1
                  ja(ncol_tot) = ja(jj)
                  i_ct = 1
               else
                  i_ct = i_ct + 1
               endif
            enddo
            if (i_ct.eq.8) then
               vals(ncol_tot) = 0.0
            else
               vals(ncol_tot)   = -1.0
               row_sum = row_sum - 1.0
            endif
            if (ja(jj-1).eq.(nrow+1)) then
               index = ncol_tot
            endif
            vals(index) = row_sum
c
            nrow      = nrow     + 1
            ia (nrow) = ncol_tot + 1
            ie_l = ie
         endif
      enddo
c
c     Shift ia by 1.
      do i=nrow,1,-1
         ia(i+1) = ia(i)
      enddo
      ia(1) = 1
c
      ncol = ia(nrow+1)-ia(1)
      write(6,*) ncol,nrow,' nnz,nrow',nseg,nv,vertex(nv+1),vertex(nv+2)
c
c     write(6,*) 'quit in compute_adj_3d'
c     call exit
c     call hmt_spy_csrmat(ia,ja,vals,nrow,'Cntr.Adj.')
c
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine rank(A,IND,N)
C
C     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
C
      integer ind(1)
      real      a(1)
C
      if (n.le.1) return
      DO 10 J=1,N
         IND(j)=j
   10 continue
C
      if (n.eq.1) return
      L=n/2+1
      ir=n
  100 CONTINUE
         IF (l.gt.1) THEN
            l=l-1
            indx=ind(l)
            q=a(indx)
         ELSE
            indx=ind(ir)
            q=a(indx)
            ind(ir)=ind(1)
            ir=ir-1
            if (ir.eq.1) then
               ind(1)=indx
               return
            endif
         ENDIF
         i=l
         j=l+l
  200    CONTINUE
         IF (J.le.IR) THEN
            IF (J.lt.IR) THEN
               IF ( A(IND(j)).lt.A(IND(j+1)) ) j=j+1
            ENDIF
            IF (q.lt.A(IND(j))) THEN
               IND(I)=IND(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GOTO 200
         ENDIF
         IND(I)=INDX
      GOTO 100
      END
C
C-----------------------------------------------------------------------
C
      subroutine irank(a,ind,n)
C
C     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
C
      integer a(1),ind(1),q
C
      if (n.le.1) return
      do 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.eq.1) return
      L=n/2+1
      ir=n
  100 continue
         IF (l.gt.1) THEN
            l=l-1
            indx=ind(l)
            q=a(indx)
         ELSE
            indx=ind(ir)
            q=a(indx)
            ind(ir)=ind(1)
            ir=ir-1
            if (ir.eq.1) then
               ind(1)=indx
               return
            endif
         ENDIF
         i=l
         j=l+l
  200    CONTINUE
         IF (J.le.IR) THEN
            IF (J.lt.IR) THEN
               IF ( a(ind(j)).lt.a(ind(j+1)) ) j=j+1
            ENDIF
            IF (q.lt.a(ind(j))) THEN
               ind(i)=ind(j)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GOTO 200
         ENDIF
         IND(I)=INDX
      GOTO 100
      end
C
C-----------------------------------------------------------------------
C
      subroutine hmt_spy_csrmat(ia,ja,vals,n,name9)
      integer ia(1),ja(1)
      real*8 vals(1)
c
      character*9 name9
      character*1 s(100)
c
      write(6,1) name9,n
    1 format(/,'CSR Mat:',a9,3x,'n =',i3,/)
c
      n100 = min(n,80)
      do i=1,n
         call blank(s,100)
         n1 = ia(i)
         n2 = ia(i+1)-1
         do jj=n1,n2
            j = ja  (jj)
            write(6,*) '(',j,vals(jj),')'
         enddo
         write(6,100) (s(k),k=1,n100)
      enddo
    5 format('X')
  100 format(100a1)
c
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine spy_csrmat(ia,ja,n,name9)
      integer ia(1),ja(1)
c
      character*9 name9
      character*1 s(100)
c
      write(6,1) name9,n
    1 format(/,'CSR Mat:',a9,3x,'n =',i3,/)
c
      n100 = min(n,80)
      do i=1,n
         call blank(s,100)
         n1 = ia(i)
         n2 = ia(i+1)-1
         do jj=n1,n2
            j = ja  (jj)
            if (j.le.n100) write(s(j),5)
         enddo
         write(6,100) (s(k),k=1,n100)
      enddo
    5 format('X')
  100 format(100a1)
c
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine cfill(a,b,n)
      real  a(1)
      do i = 1, n
         a(i) = b
      enddo
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine lfalse(a,n)
      logical  a(1)
      do i = 1, n
         a(i) = .false.
      enddo
      return
      end
C-----------------------------------------------------------------------
      subroutine period_coinc(vertex,start_pt,nv)
c
c     Adjust adjacency arrays and vertex to acct for
c     periodic bcs.
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      integer n,nv
      integer vertex(1),start_pt(1)
c
      logical ifcons
c
      integer vlist(4,2)
      integer v1,v2
C
      kfld = 2
      if (ifflow) kfld=1
C
      nc       = 2**ndim
      nfaces   = 2*ndim
      ncf      = 2**(ndim-1)
c
c     March over all elements, looking for periodic face.
c     Also, check for consistency of P-P bcs.
c
      do ie=1,nel
         do iface = 1,nfaces
            if (cbc(iface,ie,kfld).eq.'P  ') then
               je = bc(1,iface,ie,kfld)
               jf = bc(2,iface,ie,kfld)
c
               ifcons = .true.
               if (cbc(jf,je,kfld).ne.'P  '  ) ifcons = .false.
               if ( bc(1,jf,je,kfld).ne.ie   ) ifcons = .false.
               if ( bc(2,jf,je,kfld).ne.iface) ifcons = .false.
               if (.not.ifcons) then
                  call blank(line,70)
                  write(line,1)
                  call prs(line)
                  write(line,2) ie,iface
                  call prs(line)
                  write(line,2) je,jf
                  call prs(line)
    1             format('WARNING: inconsistent periodic BCs.$')
    2             format('Reset el/face2:',i6,i3,' to "p"$')
               endif
c
c              Find pairings of vertices based upon face info & geom.
c
               call find_match(vlist,iface,ie,jf,je)
c
               do i=1,ncf
                  iv1   = vlist(i,1)
                  iv2   = vlist(i,2)
                  v1    = vertex(iv1)
                  v2    = vertex(iv2)
c
              if (v2.gt.v1) 
     $        write(6,11) ie,iface,je,jf,' Replace vertex'
     $        ,iv2,v2,' with ',iv1,v1
c
              if (v1.gt.v2) 
     $        write(6,11) ie,iface,je,jf,' Replace vertex'
     $        ,iv1,v1,' with ',iv2,v2 
  11          format(2(i6,i3),a15,2i10,a6,2i10)
c
c  
c                 Replace larger of 2 vertices.  If same, do nothing.
                  if (v2.gt.v1) call reset_v(vertex,v1,v2,nv)
                  if (v1.gt.v2) call reset_v(vertex,v2,v1,nv)
c
c
               enddo
            endif
         enddo
      enddo
c
      return
      end
C-----------------------------------------------------------------------
      subroutine find_match(vlist,iface,ie,jf,je)
c
c     Return list of corresponding vertices when (iface,ie) is
c     "periodic" with (jf,je).
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      integer vlist(4,2)
c
      real xl(4,2)
      real yl(4,2)
      real zl(4,2)
c
c     order3 is an array which gives a cyclic permutation of vertices in 
c     counter-clockwise order when looking at a face from *outside* the cube.
c     Using pre-nek ordering, which is miserably non-symmetric...
      integer order3(4,6)
      save    order3
      data    order3  / 1, 2, 6, 5
     $                , 2, 3, 7, 6
     $                , 3, 4, 8, 7
     $                , 4, 1, 5, 8
     $                , 1, 4, 3, 2
     $                , 5, 6, 7, 8 /
c
c     Pairing of face (jf,je) to (iface,ie) is done by considering
c     4 rotations of (jf,je) and taking that which minimizes the
c     L2 measure of the separation of the vertices.
c
      nc  = 2**ndim
      ncf = 2**(ndim-1)
c
c     Load (iface,ie) in reverse order
c
      do ic = 1,ncf
         xl(ic,1)    = x(ie,order3(ncf+1-ic,iface))
         yl(ic,1)    = y(ie,order3(ncf+1-ic,iface))
         zl(ic,1)    = z(ie,order3(ncf+1-ic,iface))
         vlist(ic,1) =      order3(ncf+1-ic,iface) + nc*(ie-1)
      enddo
c
      dmin = 1.e22
      jrmn = 0
      do jrot = 0,ncf-1
         do jc=1,ncf
            jcr         = mod1(jc+jrot,ncf)
            xl(jc,2)    = x(je,order3(jcr,jf))
            yl(jc,2)    = y(je,order3(jcr,jf))
            zl(jc,2)    = z(je,order3(jcr,jf))
         enddo
         if (ndim.eq.2) call rzero(zl,8)
         d = 0.
         do ic=1,ncf                         !major error found 11/19/01 pff
            d = d + (xl(ic,1)-xl(ic,2))**2
     $            + (yl(ic,1)-yl(ic,2))**2
     $            + (zl(ic,1)-zl(ic,2))**2
         enddo
         if (d.lt.dmin) then
            dmin = d
            jrmn = jrot
         endif
      enddo
c
c     OK... distance minimizing rotation now found.  Choose this
c     as list
c
      do jc = 1,ncf
         jcr         = mod1(jc+jrmn,ncf)
         vlist(jc,2) = order3(jcr,jf) + nc*(je-1)
      enddo
c     write(6,1) jrmn,'rot1',(vlist(k,1),k=1,4),ie,iface
c     write(6,2) 'rot2',(vlist(k,2),k=1,4),je,jf,dmin
c   1 format(i1,a4,2x,6i7)
c   2 format(1x,a4,2x,6i7,1pe12.4)
c
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine reset_v(vertex,v1,v2,nv)
c
      integer vertex(1),v1,v2
c
c     Update vertex list, replacing v2 ref's w/ v1 refs.
c
      do i=1,nv
         if (vertex(i).eq.v2) vertex(i) = v1
      enddo
c
      return
      end
C
C-----------------------------------------------------------------------
C
      subroutine out_csrmat(acsr,ia,ja,n,name)
      real    acsr(1)
      integer ia(1),ja(1)
c
      character*9 name
      character*5 s(33)
c
      nnz = ia(n+1)-ia(1)
      write(6,*) 'outcs: ',name,n,nnz
      write(6,9) (acsr(k),k=1,nnz)
      write(6,8) (ja(k),k=1,nnz)
    8 format(10i10)
    9 format(1p10e10.2)
c
      write(6,1) name,n
    1 format(/,'CSR Mat:',a9,3x,'n =',i3,/)
c
      n33 = min(n,26)
      do i=1,n
         call blank(s,132)
         n1 = ia(i)
         n2 = ia(i+1)-1
         do jj=n1,n2
            j = ja  (jj)
            a = acsr(jj)
            if (a.ne.0..and.j.le.n33) write(s(j),5) a
         enddo
         write(6,33) (s(k),k=1,n33)
      enddo
    5 format(f5.2)
   33 format(33a5)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine period_check(ifld)
c
c     Adjust adjacency arrays and vertex to acct for
c     periodic bcs.
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      nc       = 2**ndim
      nfaces   = 2*ndim
      ncf      = 2**(ndim-1)
c
c     March over all elements, looking for periodic face.
c     Also, check for consistency of P-P bcs.
c
      do ie=1,nel
         do iface = 1,nfaces
            if (cbc(iface,ie,ifld).eq.'P  ') then
               je = bc(1,iface,ie,ifld)
               jf = bc(2,iface,ie,ifld)
c
               icons = 0
               if (cbc(jf,je,ifld).ne.'P  '  ) icons = 1
               if ( bc(1,jf,je,ifld).ne.ie   ) icons = 2
               if ( bc(2,jf,je,ifld).ne.iface) icons = 3
               if (icons.ne.0) then
                  call blank(line,70)
                  write(line,1) icons,cbc(jf,je,ifld)
                  call prs(line)
                  write(line,2) ie,iface,ifld
                  call prs(line)
                  write(line,2) je,jf,ifld
                  call prs(line)
    1             format('WARNING: inconsistent per. BCs:',i2,1x,a3,'$')
    2             format('Reset el/face1:',i6,2i3,' to "p"$')
                  cbc(iface,ie,ifld) = 'p  '
                  call rzero(bc(1,iface,ie,ifld),5)
c
               endif
c
            endif
         enddo
      enddo
c
      return
      end
C-----------------------------------------------------------------------
      subroutine period_bc_check(if_any_per)
c
c     Adjust adjacency arrays and vertex to acct for
c     periodic bcs.
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      logical if_any_per
C
      nc       = 2**ndim
      nfaces   = 2*ndim
      ncf      = 2**(ndim-1)
c
c     March over all elements, looking for periodic face.
c     Also, check for consistency of P-P bcs.
c
      ifld0=2
      ifld1=1
      if (ifflow) ifld0=1
      if (ifheat) ifld1=2+npscal
c
      if_any_per = .false.
      do ifld = ifld0,ifld1
      do ie=1,nel
         do iface = 1,nfaces
            if (cbc(iface,ie,ifld).eq.'P  ') if_any_per = .true.
         enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_csr_matl(ia,ja,vals,n)
      integer ia(1),ja(1)
      real*8 vals(1)
c
      open(unit=39,file='csr.dat',status='unknown')
      do i=1,n
         j0=ia(i  )
         j1=ia(i+1)-1
         do j=j0,j1
            write(39,39) i,ja(j),vals(j)
         enddo
      enddo
   39 format(2i10,1p1e19.9)
      close(unit=39)
c
      return
      end
c-----------------------------------------------------------------------
