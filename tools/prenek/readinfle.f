      subroutine check_inputfile(filen,parfle,re2fle,ifpar,ifre2,ierr)

c     RETURNS LINE AND PARFLE AS THE .REA AND .PAR FILE NAMES
C     RETURNS IFPAR AND IFRE2 LOGICALS FOR EXISTANCE OF FILES

      logical ifpar,ifre2,ifrea,isthere

      character file1*80,file1s(80)
      equivalence(file1,file1s(1))

      character temp*80,temps(80)
      equivalence(temp,temps(1))

      character filen*80,parfle*80,re2fle*80
      integer len

      ifpar=.false.
      ifre2=.false.    !default rea only case
      ifrea=.false.
      isthere=.false.
   
      ierr=0

      call blank(file1,80)
      call blank(temp,80)
      call blank(parfle,80)
      call blank(re2fle,80)
      temp=filen

      if (indx1(filen,'.re2',4).ne.0) then ! .re2 present
         re2fle= filen
         inquire(file=re2fle,exist=isthere)
         if(isthere) then
            ifre2 =.true.
            isthere=.false.
         endif

         len = ltrunc(filen,80)
         
         call chcopy(file1s,temps,len-4)
         call chcopy(file1s(len-3),'.par',4)
         inquire(file=file1,exist=isthere)
         if(isthere) then
            ifpar  =.true.
            parfle = file1
            isthere=.false.
         endif

         call chcopy(file1s,temps,len-4)
         call chcopy(file1s(len-3),'.rea',4)
         inquire(file=file1,exist=isthere)
         if(isthere) then
            ifrea  =.true.
            isthere=.false.
            filen=file1   !set filen to be the .rea file name
         endif
         if(.not.ifre2.or.(.not.ifpar.and..not.ifrea))then
            call prs('Error finding re2 and supporting files$')
            ierr=1
            return
         endif
      elseif (indx1(filen,'.rea',4).ne.0) then ! .rea present
         inquire(file=filen,exist=isthere)
         if(isthere) then
            ifrea =.true.
            isthere=.false.
         endif

c       IF .rea SPECIFIED, USE REA AND NOT ANY OTHER FILE TYPE
c        len = ltrunc(filen,80)
c        call chcopy(file1s,temps,len-4)
c        call chcopy(file1s(len-3),'.par',4)
c        inquire(file=file1,exist=isthere)
c        if(isthere) then
c           ifpar  =.true.
c           parfle = file1
c           isthere=.false.
c        endif
c
c        call chcopy(file1s,temps,len-4)
c        call chcopy(file1s(len-3),'.re2',4)
c        inquire(file=file1,exist=isthere)
c        if(isthere) then
c          ifre2  =.true.
c          re2fle = file1
c          isthere=.false.
c        endif
         if(.not.ifrea) then
            ierr=2
            call prs('Error finding .rea file$')
            return
         endif
      elseif (indx1(filen,'.par',4).ne.0) then ! .rea present
         parfle= filen
         inquire(file=parfle,exist=isthere)
         if(isthere) then
            ifpar =.true.
            isthere=.false.
         endif

         len = ltrunc(filen,80)
         call chcopy(file1s,temps,len-4)
         call chcopy(file1s(len-3),'.re2',4)
         inquire(file=file1,exist=isthere)
         if(isthere) then
            ifre2  =.true.
            re2fle = file1
            isthere=.false.
         endif

c        call chcopy(file1s,temps,len-4)
c        call chcopy(file1s(len-3),'.rea',4)
c        inquire(file=file1,exist=isthere)
c        if(isthere) then
c           ifrea  =.true.
c           isthere=.false.
c           filen  = file1
c        endif
         if(.not.ifpar.or..not.ifre2) then
            ierr=3
            call prs('Error finding .par or supporting files$')
            return
         endif
      else
         len = ltrunc(filen,80)

         file1=temp
         call chcopy(file1s(len+1),'.re2',4)
         inquire(file=file1,exist=isthere)
         if(isthere) then
            ifre2  =.true.
            re2fle = file1
            isthere=.false.
         endif

         file1=temp
         call chcopy(file1s(len+1),'.rea',4)
         inquire(file=file1,exist=isthere)
         if(isthere) then
            ifrea  =.true.
            isthere=.false.
            filen  = file1
         endif
          
         file1=temp
         call chcopy(file1s(len+1),'.par',4)
         inquire(file=file1,exist=isthere)
         if(isthere) then
            ifpar  =.true.
            parfle = file1
            isthere=.false.
         endif
         if(.not.ifrea.and..not.ifre2) then
           call prs('Error finding .rea or .re2 file$')
           ierr=4
           return
         endif
         if(ifre2.and.(.not.ifrea.and..not.ifpar)) then
           call prs('Error finding supporting .re2 file$')
           ierr=5
           return
         endif
      endif 

c     write(6,*) ifre2,re2fle,ifpar,parfle,ifrea,filen

      return
      end
c---------------------------------------------------------------
      subroutine setDefaultParam
      include 'basics.inc'

      ifxyo = .true.
      ifvo  = .false.
      ifpo  = .false.
      ifto  = .false.

      ifadvc(1) = .true.
      do i=1,MPSCAL
         ifadvc(i+1) = .true.
      enddo

      ifflow    = .false.
      ifheat    = .false.
      iftran    = .true.
      ifaxis    = .false.
      ifstrs    = .false.
      ifmvbd    = .false.
      ifchar    = .false.
      ifmgrid   = .false.
      ifmodel   = .false.
      ifkeps    = .false.
 
      if3d=.false.   !!!DUMMY SET FOR NOW
      ndim=2         !!!TO BE READ IN RE2

      return
      end
c------------------------------------------------------------
      subroutine param_read_par(parfle)
      include 'basics.inc'


      character*132 c_out,txt
      character*80 parfle

      call finiparser_load(parfle,ierr)
      if(ierr .ne. 0) return

      call par_verify(ierr)
c     if(ierr .ne. 0) return


      call finiparser_find(i_out,'temperature',ifnd)
      if(ifnd .eq. 1) then
        ifheat = .true.
        ifto   = .true.
      endif

      j = 0
      do i = 1,MPSCAL-1
         write(txt,"('scalar',i2.2)") i
         call finiparser_find(i_out,txt,ifnd)
         if (ifnd .eq. 1) then
            j = j + 1
            ifpsco(i) = .true.
         endif
      enddo

      if (index(c_out,'CHAR') .gt. 0) then
         ifchar = .true.
      else if (index(c_out,'STEADY') .gt. 0) then
         iftran = .false.
      endif

      call finiparser_find(i_out,'velocity',ifnd)
      if(ifnd .eq. 1) then
        ifflow = .true.
        ifvo   = .true.
        ifpo   = .true.
      endif

      call finiparser_getBool(i_out,'mesh:motion',ifnd)
      if(ifnd .eq. 1 .and. i_out .eq. 1) then
        ifmvbd = .true.
      endif

      call finiparser_getBool(i_out,'problemType:axiSymmetry',ifnd)
      if(ifnd .eq. 1) then
        ifaxis = .false.
        if(i_out .eq. 1) ifaxis = .true.
      endif

      call finiparser_getBool(i_out,
     &                        'problemType:stressFormulation',ifnd)
      if(ifnd .eq. 1) then
        ifstrs = .false.
        if(i_out .eq. 1) ifstrs = .true.
      endif

      call finiparser_getBool(i_out,'velocity:advection',ifnd)
      if(ifnd .eq. 1) then
        ifadvc(1) = .false.
        if(i_out .eq. 1) ifadvc(1) = .true.
      endif

      call finiparser_getBool(i_out,'temperature:advection',ifnd)
      if(ifnd .eq. 1) then
        ifadvc(2) = .false.
        if(i_out .eq. 1) ifadvc(2) = .true.
      endif

      do i = 1,mpscal-1
         write(txt,"('scalar',i2.2,a)") i,':advection'
         call finiparser_getBool(i_out,txt,ifnd)
         if(ifnd .eq. 1) then
           ifadvc(i+2) = .false.
           if(i_out .eq. 1) ifadvc(i+2) = .true.
         endif
      enddo

      call finiparser_getBool(i_out,'mesh:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifxyo = .false.
        if(i_out .eq. 1) ifxyo = .true.
      endif

      call finiparser_getBool(i_out,'velocity:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifvo = .false.
        if(i_out .eq. 1) ifvo = .true.
      endif

      call finiparser_getBool(i_out,'pressure:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifpo = .false.
        if(i_out .eq. 1) ifpo = .true.
      endif

      call finiparser_getBool(i_out,'temperature:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifto = .false.
        if(i_out .eq. 1) ifto = .true.
      endif

      j=0
      do i = 1,mpscal-1
         write(txt,"('scalar',i2.2,a)") i,':writeToFieldFile'
         call finiparser_getBool(i_out,txt,ifnd)
         if(ifnd .eq. 1) then
           j=j+1
           ifpsco(i) = .false.
           if(i_out .eq. 1) ifpsco(i) = .true.
         endif
      enddo
      param(23)=j
      npscal=param(23)


      return
      end
c-----------------------------------------------------------------------
      subroutine par_verify(ierr)

      INCLUDE 'PARDICT'


      character*132  key
      character*1024 val

      character*132 txt
      character*1   tx1(132)
      equivalence   (tx1,txt)

      ierr = 0

      call finiparser_getDictEntries(n)
      do i = 1,n
         call finiparser_getPair(key,val,i,ifnd)
         call capit(key,132)

         is = index(key,'_') ! ignore user keys
         if (is.eq.1) goto 10

         do j = 1,PARDICT_NKEYS ! do we find the key in the par-dictionary
            if(index(pardictkey(j),key).eq.1) goto 10

            is = index(key,'SCALAR')
            if(is .eq. 1) then
              call chcopy(txt,key,132)
              call chcopy(tx1(is+6),'%%',2)
              if(index(pardictkey(j),txt).eq.1) goto 10
            endif
         enddo
         write(6,*) 'ERROR: Par file contains unknown key ', key
         ierr = ierr + 1
   10 enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine readat_re2
      logical ifbswap

      call open_bin_file(ifbswap,zh)
      call bin_rd1(ifbswap,zh)

      return
      end
c-----------------------------------------------------------------------
      subroutine open_bin_file(ifbswap,zhght) ! open file & chk for byteswap
      include 'basics.inc'

      common /newfle/re2fle,parfle
      character re2fle*80,parfle*80

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv

      common /lfiles/ ifre2,ifpar
      logical ifre2,ifpar

      logical ifbswap,if_byte_swap_test

      integer fnami (33)
      character*132 fname
      equivalence (fname,fnami)

      character*132 hdr
      character*5 version
      real*4      test

      ierr=0
      zhght=0.0

      call izero(fnami,33)
      m = indx2(re2fle,80,' ',1)-1
      call chcopy(fname,re2fle,m)
c     fname=re2fle

      call byte_open(fname,ierr)
      if(ierr.ne.0) goto 100
      call byte_read(hdr,20,ierr)
      if(ierr.ne.0) goto 100

      read (hdr,1) version,nelgt,ndim,nelgv
    1 format(a5,i9,i3,i9)

      nel=nelgt
      if (nel.ge.nelm) then
         call prsii('NELM too small in basics.inc.$',nel,nelm)
         call prs  ('Recompile prenek.  ABORT$')
      endif
 
      if (ndim.eq.3.and.ifpar) if3d=.true.

      IF(IF3D.AND.NDIM.EQ.2)THEN
         CALL PRS
     $        ('WARNING: PARAMETER SET TO 2D BUT DATA IS FOR 3D$')
         CALL PRS('PLEASE ENTER Z-HEIGHT: $')
         CALL RER( ZHGHT)
      ENDIF
      IF(.NOT.IF3D.AND.NDIM.EQ.3) THEN
         CALL PRS(
     $        'ERROR: PARAMETER SET TO 3-D BUT DATA IS FOR 2-D$')
      ENDIF

      wdsizi = 4
      if(version.eq.'#v002') wdsizi = 8
      if(version.eq.'#v003') then
         wdsizi = 8
         param(32) = 1
      endif

      call byte_read(test,1,ierr)
      if(ierr.ne.0) goto 100
      ifbswap = if_byte_swap_test(test,ierr)
      if(ierr.ne.0) goto 100

      return

 100  call prs('Error opening or reading .re2 header. Abort.$')

      return
      end

c-----------------------------------------------------------------------
      subroutine bin_rd1(ifbswap,zh)  ! read mesh, curve, and bc info

      include 'basics.inc'

      logical ifbswap
 

                  i_bc = 2
      if (ifflow) i_bc = 1

                  nfldt = 1
      if (ifheat) nfldt = 2+npscal

c     if (param(32).gt.0) nfldt = i_bc + param(32)-1
      lcbc=18*nelm*(mpscal + 1)
      call blank(cbc,lcbc)

      call bin_rd1_mesh  (ifbswap,zh)   ! version 1 of binary reader
      call bin_rd1_curve (ifbswap)

      do ifield = i_bc,nfldt
         call bin_rd1_bc (cbc(1,1,ifield),bc(1,1,1,ifield),ifbswap)
      enddo

      ierr=0
      call byte_close(ierr)
      if (ierr.ne.0) call prs('Error closing re2 file. Abort $')

      return
      end
c-----------------------------------------------------------------------
      subroutine bin_rd1_mesh(ifbswap,zh)    ! version 1 of binary reader

      include 'basics.inc'

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv
      logical ifbswap

      integer e,eg,buf(55)

      nwds = (1 + ndim*(2**ndim))*(wdsizi/4) ! group + 2x4 for 2d, 3x8 for 3d

      if (nwds.gt.55) then
         call prs(' Error in bin_rd1_mesh: buf size')
         call exitt
      endif

      ierr  = 0
      ierr2 = 0
      do eg=1,nelgt 
         numapt(eg)=1  !Std. defaults
         letapt(eg)='A'
         igroup    =0
         call byte_read  (buf,nwds,ierr)
         call buf_to_xyz (buf,eg,ifbswap,zh,ierr2)
         ierr = ierr + ierr2
         if(ierr.ne.0) goto 100
      enddo

      return

 100  if (ierr.ne.0) call prs('Error reading .re2 mesh. Abort. $')

      return
      end

c-----------------------------------------------------------------------
      subroutine buf_to_xyz(buf,e,ifbswap,zhght,ierr)! version 1 of binary reader
      include 'basics.inc'

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv

      logical ifbswap

      integer e,eg,buf(0:49)

      nwds = (1 + ndim*(2**ndim))*(wdsizi/4) ! group + 2x4 for 2d, 3x8 for 3d

      if     (ifbswap.and.ierr.eq.0.and.wdsizi.eq.8) then
          call byte_reverse8(buf,nwds,ierr)
      elseif (ifbswap.and.ierr.eq.0.and.wdsizi.eq.4) then
          call byte_reverse (buf,nwds,ierr)
      endif
      if(ierr.ne.0) return

      if(wdsizi.eq.8) then
         call copyi4(igroup(e),buf(0),1) !0-1
         if (ndim.eq.3) then
            call copy  (x(1,e),buf( 2),8) !2 --17
            call copy  (y(1,e),buf(18),8) !18--33
            call copy  (z(1,e),buf(34),8) !34--49
         else
            call copy  (x(1,e),buf( 2),4) !2 --9
            call copy  (y(1,e),buf(10),4) !10--17
            if(if3d) then
               zero=0.0
               call copy (x(5,e),x(1,e),4)
               call copy (y(5,e),y(1,e),4)
               call copy (z(1,e),zero,4)
               call copy (z(5,e),zhght,4)
            endif
          endif
      else
         igroup(e) = buf(0)
         if (ndim.eq.3) then
            call copy4r(x(1,e),buf( 1),8)
            call copy4r(y(1,e),buf( 9),8)
            call copy4r(z(1,e),buf(17),8)
         else
            call copy4r(x(1,e),buf( 1),4)
            call copy4r(y(1,e),buf( 5),4)
c          write(6,*) e, 'EEEEE' , x(1,e),igroup(e)
            if(if3d) then
               zero=0.0
               call copy4r (x(5,e),x(1,e),4)
               call copy4r (y(5,e),y(1,e),4)
               call copy4r (z(1,e),zero,4)
               call copy4r (z(5,e),zhght,4)
            endif
         endif
      endif

      return
      end

c-----------------------------------------------------------------------      
      subroutine bin_rd1_curve (ifbswap) ! v. 1 of curve side reader

      include 'basics.inc'

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv

      logical ifbswap

      integer e,eg,buf(55)
      real rcurve

      nwds = (2 + 1 + 5)*(wdsizi/4) 
      len  = 4*nwds      ! 4 bytes / wd

      if (nwds.gt.55) then
         call prs(' Error in bin_rd1_curve: buf size')
         call exitt
      endif

      ierr = 0
      len1 = 4
      if(wdsizi.eq.8) then
         call byte_read(rcurve,2,ierr)
         if (ifbswap) call byte_reverse8(rcurve,2,ierr)
         ncurve = rcurve
      else
         call byte_read(ncurve,1,ierr)
         if (ifbswap) call byte_reverse(ncurve,1,ierr)
      endif
     
      if (ierr.ne.0) goto 100

      do k=1,ncurve
         call byte_read(buf,nwds,ierr)
         if (ierr.ne.0) goto 100

         if(wdsizi.eq.8) then
           if(ifbswap) call byte_reverse8(buf,nwds-2,ierr)
c          call copyi4(eg,buf(1),1)  !1,2
         else
           if (ifbswap) call byte_reverse(buf,nwds-1,ierr) ! last is char
c          eg  = buf(1)
         endif

         call buf_to_curve(buf)
      enddo
      
 100  if (ierr.ne.0) call prs('Error reading .re2 curved data. Abort.$')

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_curve(buf)    ! version 1 of binary reader
      include 'basics.inc'

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv

      integer e,eg,f,buf(30)

      if(wdsizi.eq.8) then
        call copyi4(eg,buf(1),1) !1-2
c       e  = gllel(eg)

        call copyi4(f,buf(3),1) !3-4

        call copy  ( curve(1,f,eg),buf(5) ,5) !5--14
        call chcopy(ccurve(  f,eg),buf(15),1)!15
      else
        eg = buf(1)
c       e  = gllel(eg)
        f  = buf(2)

        call copy4r( curve(1,f,eg),buf(3),5)
        call chcopy(ccurve(f,eg)  ,buf(8),1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bin_rd1_bc (cbl,bl,ifbswap) ! v. 1 of bc reader

      include 'basics.inc'

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv

      logical ifbswap

      character*3 cbl(6,nelm)
      real         bl(5,6,nelm)

      integer e,eg,buf(55)
      real rbc_max

      nwds = (2 + 1 + 5)*(wdsizi/4)   ! eg + iside + cbc + bc(5,:,:)
      len  = 4*nwds      ! 4 bytes / wd

      if (nwds.gt.55) then
         call prs(' Error in bin_rd1_bc: buf size')
         call exitt
      endif

      do e=1,nel   ! fill up cbc w/ default
      do k=1,6
         cbl(k,e) = 'E  '
      enddo
      enddo

      ierr=0

      if(wdsizi.eq.8) then
         call byte_read(rbc_max,2,ierr)
         if (ifbswap) call byte_reverse8(rbc_max,2,ierr) ! last is char
         nbc_max = rbc_max
      else
         call byte_read(nbc_max,1,ierr)
         if (ifbswap) call byte_reverse(nbc_max,1,ierr) ! last is char
      endif

      do k=1,nbc_max
         call byte_read(buf,nwds,ierr)
         if (ierr.ne.0) goto 100
         if(wdsizi.eq.8) then
           if (ifbswap) call byte_reverse8(buf,nwds-2,ierr)
           if (ierr.ne.0) goto 100
           call copyi4(eg,buf(1),1) !1&2 of buf
         else
           if (ifbswap) call byte_reverse(buf,nwds-1,ierr) ! last is char
           if (ierr.ne.0) goto 100
           eg  = buf(1)
         endif

         call buf_to_bc(cbl,bl,buf)
      enddo

      return 

 100  call prs('Error reading boundary data for re2. Abort.$')

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_bc(cbl,bl,buf)    ! version 1 of binary reader

      include 'basics.inc'

      common /ibinfle/ wdsizi,nelgt,nelgv
      integer wdsizi,nelgt,nelgv

      character*3 cbl(6,nelm)
      real         bl(5,6,nelm)

      integer e,eg,f,buf(30)

      if(wdsizi.eq.8) then
        call copyi4(eg,buf(1),1) !1-2
c       e  = gllel(eg)

        call copyi4(f,buf(3),1) !3-4

        call copy  (bl(1,f,eg),buf(5),5) !5--14
        call chcopy(cbl( f,eg),buf(15),3)!15-16

        if(nel.ge.1000000.and.cbl(f,e).eq.'P  ')
     $   call copyi4(bl(1,f,eg),buf(5),1) !Integer assign connecting P element

      else
        eg = buf(1)
c       e  = gllel(eg)
        f  = buf(2)

        call copy4r ( bl(1,f,eg),buf(3),5)
        call chcopy (cbl(  f,eg),buf(8),3)
c       write(6,*) cbl(f,eg)

        if (nel.ge.1 000 000.and.cbl(f,e).eq.'P  ')
     $     bl(1,f,eg) = buf(3) ! Integer assign of connecting periodic element
      endif

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest,ierr)

      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern

      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.

      test2 = bytetest
      call byte_reverse(test2,1,ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine copy4r(a,b,n)
      real   a(1)
      real*4 b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copyi4(a,b,n)
      integer a(1)
      real    b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function indx2(s1,l1,s2,l2)
      character*80 s1,s2

      n1=l1-l2+1
      indx2=0
      if (n1.lt.1) return

      do i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx2=i
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
