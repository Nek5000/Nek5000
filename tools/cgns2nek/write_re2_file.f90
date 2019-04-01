!****************************************************************************
!
!        Write re2 file for Nek5000
!
!********************************************************************
!
SUBROUTINE write_re2_file ()
  !
  USE CGNS
  USE MODULE_global
  !
  IMPLICIT NONE
  !
  CHARACTER (len=132) :: re2file
  CHARACTER (len=80) ::  hdr
  CHARACTER(1) :: dir,var1
  !
  INTEGER(iprec) :: i,ii,NEL,NDIM,NELV,var2,n,c_node(12),l,j,nbc
  INTEGER(iprec) :: ie,i1,i2,i3,i4,i5,i6,i7,i8,nmax,re2len,ierr
  REAL(high)     :: var3,rgroup,buf2(30),rbc
  REAL*4 test,buf(10) 
  !
  DATA   test  / 6.54321 /
  DATA c_node  / 9,10,11,12,17,18,19,20,13,14,15,16 /
  !
  !
  ! Write mesh coordinates, curvature data and boundary data in re2 file
  !
  re2file = cgns_file(1:cgns_len) //'.re2'
  re2len = INDEX (re2file, ' ') - 1
  !
  WRITE(6,'(A,A)') ' Writing ', re2file 
  CALL byte_open(re2file,ierr)
  !
  !   Coordinates
  !
  CALL blank(hdr,80)
  !
  WRITE(hdr,1) max_elem, 3, max_elem
1 FORMAT('#v002',i9,i3,i9,' this is the hdr')

  CALL byte_write(hdr,20,ierr)
  CALL byte_write(test,1,ierr)     ! write the endian discriminator
  !
  rgroup = 0. ! for now skip rgroup (not used at the moment)
  ii = 0
  !
  DO ie = 1,max_elem
     !
     CALL byte_write(rgroup, 2,ierr)
     !
     i1= elements(sec,ii+1)
     i2= elements(sec,ii+2)
     i3= elements(sec,ii+3)
     i4= elements(sec,ii+4)
     i5= elements(sec,ii+5)
     i6= elements(sec,ii+6)
     i7= elements(sec,ii+7)
     i8= elements(sec,ii+8)
     !
     buf2(1)   = xv(i1)
     buf2(2)   = xv(i2)
     buf2(3)   = xv(i3)
     buf2(4)   = xv(i4)
     buf2(5)   = xv(i5)
     buf2(6)   = xv(i6)
     buf2(7)   = xv(i7)
     buf2(8)   = xv(i8)

     buf2(9)   = yv(i1)
     buf2(10)  = yv(i2)
     buf2(11)  = yv(i3)
     buf2(12)  = yv(i4)
     buf2(13)  = yv(i5)
     buf2(14)  = yv(i6)
     buf2(15)  = yv(i7)
     buf2(16)  = yv(i8)

     buf2(17)  = zv(i1)
     buf2(18)  = zv(i2)
     buf2(19)  = zv(i3)
     buf2(20)  = zv(i4)
     buf2(21)  = zv(i5)
     buf2(22)  = zv(i6)
     buf2(23)  = zv(i7)
     buf2(24)  = zv(i8)

     CALL byte_write  (buf2(1) ,16, ierr)
     CALL byte_write  (buf2(9) ,16, ierr)
     CALL byte_write  (buf2(17),16, ierr)
     !
     ii = ii+inc
     !
  ENDDO
  !
  ! CURVED SIDE data
  !
  var2 = 0
  var3 = 0.0        
  !
  !  ***** CURVED SIDE DATA ***** '

  IF (sectype(sec) == HEXA_20 .OR. sectype(sec) == HEXA_27) THEN
     !
     IF (sectype(sec) == HEXA_20) THEN
        var2 = max_elem*8
        nmax = 8
     ELSE
        var2 = max_elem*12
        nmax = 12
     ENDIF
     !
     rbc = var2
     !
     CALL byte_write(rbc,2, ierr)
     !
     ii = 0
     CALL blank(var1,4)
     var1= 'm'
     !
     DO ie = 1,max_elem
        DO n = 1,nmax
           i= elements(sec,ii+c_node(n))
           buf2(1) = DBLE(ie) 
           buf2(2) = DBLE(n)
           buf(1) = xv(i)
           buf(2) = yv(i)
           buf(3) = zv(i)
           buf(4) = var3 
           buf(5) = var3 
           CALL copy48(buf2(3),buf(1),5)
           !call blank(buf2(8),8)
           CALL chcopy (buf2(8),var1,1)
           CALL byte_write (buf2,16,ierr)
        ENDDO
        ii = ii+inc
     ENDDO
     !
  ELSE
     rbc = var2
     CALL byte_write(rbc,2, ierr)
  ENDIF
  !
  !  ***** FLUID BOUNDARY CONDITIONS *****  '
  !
  nbc = 0
  !
  DO i = 1,max_elem
     DO n = 1,6
        IF(vel_type(i,n)  /=  'E  ' ) nbc = nbc+1
     ENDDO
  ENDDO
  !
  rbc=nbc
  CALL byte_write(rbc,2,ierr)
  !
  DO i = 1,max_elem
     DO n = 1,6
        IF(vel_type(i,n) /= 'E  ' )  THEN       
           buf2(1) = DBLE(i)
           buf2(2) = DBLE(n)
           buf2(3) =vel_val(i,n,1)
           buf2(4) =vel_val(i,n,2)
           buf2(5) =vel_val(i,n,3)
           buf2(6) =vel_val(i,n,4)
           buf2(7) =vel_val(i,n,5)
           CALL blank(buf2(8),8)
           CALL chcopy (buf2(8),vel_type(i,n),3)
           CALL byte_write (buf2,16,ierr)
        ENDIF
     ENDDO
  ENDDO
  
  RETURN
  !
END SUBROUTINE write_re2_file

SUBROUTINE blank(s,n)
  CHARACTER(1) s(1)
  DO i=1,n
     s(i)=' '
  ENDDO
  RETURN
END SUBROUTINE blank

SUBROUTINE chcopy(a,b,n)
  CHARACTER(1) A(1), B(1)

  DO I = 1, N
     A(I) = B(I)
  ENDDO
  !
  RETURN
END SUBROUTINE chcopy

SUBROUTINE copy48(a,b,n)
  REAL*8 a(1)
  REAL*4 b(1)
  DO 100 i = 1, n
100  a(i) = b(i)
     RETURN
   END SUBROUTINE copy48
