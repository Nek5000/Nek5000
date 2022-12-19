!--------------------------------------------------------------------
      subroutine right_hand_check
! check if there is non-right hand elements (3D)
! because if mesh is from ICEM, and mirror operation is made in ICEM,
! the exported exo file will contain non-right hand elements.
! this subroutine will:
! 1. check right-hand
! 2. fix if not
      use SIZE
      integer iel
      logical ifnonrighthand
      character nek_check

      write(6,*) 'performing non-right-hand check'
	  
!      do iel=1,num_elem
!         !write(6,*) 'performing non-right-hand check on element ',iel
!         call check_if_non_right(ifnonrighthand,iel)
!         if (ifnonrighthand) call fix_if_non_right(iel)		 
!      enddo

!      write(6,*) 'done: non-right-hand check'
	  
      write(6,*) 'using nek-method to do non-right-hand check? (y/n)'
      read (5,*) nek_check
 
      if (nek_check.eq.'y') then
        do iel=1,num_elem
          call nek_check_non_right_hand(iel)
       enddo
      endif
	  
      return
      end
!--------------------------------------------------------------------
      subroutine nek_check_non_right_hand(iel)
      use SIZE
      logical ifnonrighthand
      integer iel
      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,7,9,19,21,25,27/ ! for nek non-right-hand element check
      real XYZ(3,8)
      real V1,V2,V3,V4,V5,V6,V7,V8

      do iver = 1,8
       XYZ(1,iver) = xm1(hex8_to_hex27_vertex(iver),1,1,iel)
       XYZ(2,iver) = ym1(hex8_to_hex27_vertex(iver),1,1,iel)
       XYZ(3,iver) = zm1(hex8_to_hex27_vertex(iver),1,1,iel)       
      enddo
	  
      V1= VOLUM0(XYZ(1,2),XYZ(1,3),XYZ(1,5),XYZ(1,1))
      V2= VOLUM0(XYZ(1,4),XYZ(1,1),XYZ(1,6),XYZ(1,2))
      V3= VOLUM0(XYZ(1,1),XYZ(1,4),XYZ(1,7),XYZ(1,3))
      V4= VOLUM0(XYZ(1,3),XYZ(1,2),XYZ(1,8),XYZ(1,4))
      V5=-VOLUM0(XYZ(1,6),XYZ(1,7),XYZ(1,1),XYZ(1,5))
      V6=-VOLUM0(XYZ(1,8),XYZ(1,5),XYZ(1,2),XYZ(1,6))
      V7=-VOLUM0(XYZ(1,5),XYZ(1,8),XYZ(1,3),XYZ(1,7))
      V8=-VOLUM0(XYZ(1,7),XYZ(1,6),XYZ(1,4),XYZ(1,8))

      if ((V1.LE.0.0).OR.(V2.LE.0.0).OR. &
       (V3.LE.0.0).OR.(V4.LE.0.0).OR. &
       (V5.LE.0.0).OR.(V6.LE.0.0).OR. &
       (V7.LE.0.0).OR.(V8.LE.0.0)) then
   
      write(6,*) 'WARNINGb: Detected non-right-handed element.'
      write(6,*) 'at location:',XYZ(1,1),',',XYZ(2,1),',',XYZ(3,1)
      endif

      return
      end
!------------------------------------------------------------------------
      subroutine nek_check_non_right_hand_per_element(XYZorg,ifnonrighthand)
      use SIZE
      logical ifnonrighthand
!      integer iel
!      integer hex8_to_hex27_vertex(8)
!      data hex8_to_hex27_vertex /1,3,7,9,19,21,25,27/ ! for nek non-right-hand element check
      real*8 XYZ(3,8),XYZorg(3,8)
      real V1,V2,V3,V4,V5,V6,V7,V8

      do iver = 1,8
       XYZ(1,iver) = XYZorg(1,iver)
       XYZ(2,iver) = XYZorg(2,iver)
       XYZ(3,iver) = XYZorg(3,iver)
      enddo
	  
      ! swap 3-4
       XYZ(1,3) = XYZorg(1,4)
       XYZ(2,3) = XYZorg(2,4)
       XYZ(3,3) = XYZorg(3,4)
       XYZ(1,4) = XYZorg(1,3)
       XYZ(2,4) = XYZorg(2,3)
       XYZ(3,4) = XYZorg(3,3)
	   
      ! swap 7-8
       XYZ(1,7) = XYZorg(1,8)
       XYZ(2,7) = XYZorg(2,8)
       XYZ(3,7) = XYZorg(3,8)
       XYZ(1,8) = XYZorg(1,7)
       XYZ(2,8) = XYZorg(2,7)
       XYZ(3,8) = XYZorg(3,7)

	  
      ifnonrighthand = .false.
	  
      V1= VOLUM0(XYZ(1,2),XYZ(1,3),XYZ(1,5),XYZ(1,1))
      V2= VOLUM0(XYZ(1,4),XYZ(1,1),XYZ(1,6),XYZ(1,2))
      V3= VOLUM0(XYZ(1,1),XYZ(1,4),XYZ(1,7),XYZ(1,3))
      V4= VOLUM0(XYZ(1,3),XYZ(1,2),XYZ(1,8),XYZ(1,4))
      V5=-VOLUM0(XYZ(1,6),XYZ(1,7),XYZ(1,1),XYZ(1,5))
      V6=-VOLUM0(XYZ(1,8),XYZ(1,5),XYZ(1,2),XYZ(1,6))
      V7=-VOLUM0(XYZ(1,5),XYZ(1,8),XYZ(1,3),XYZ(1,7))
      V8=-VOLUM0(XYZ(1,7),XYZ(1,6),XYZ(1,4),XYZ(1,8))

      if ((V1.LE.0.0).OR.(V2.LE.0.0).OR. &
       (V3.LE.0.0).OR.(V4.LE.0.0).OR. &
       (V5.LE.0.0).OR.(V6.LE.0.0).OR. &
       (V7.LE.0.0).OR.(V8.LE.0.0)) then
   
      !write(6,*) 'WARNINGb: Detected non-right-handed element.'
      !write(6,*) 'at location:',XYZ(1,1),',',XYZ(2,1),',',XYZ(3,1)
      ifnonrighthand = .true.
      endif

      return
      end
!--------------------------------------------------------------------
      subroutine check_if_non_right(ifnonrighthand,iel)
      use SIZE
      logical ifnonrighthand
      integer iel
      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,9,7,19,21,27,25/
      real*8 hex8_vertex(3,8),vec12(3),vec14(3),vec15(3)
      real*8 vec1(3),AA,dot_prod

      do iver = 1,8
       hex8_vertex(1,iver) = xm1(hex8_to_hex27_vertex(iver),1,1,iel)
       hex8_vertex(2,iver) = ym1(hex8_to_hex27_vertex(iver),1,1,iel)
       hex8_vertex(3,iver) = zm1(hex8_to_hex27_vertex(iver),1,1,iel)       
      enddo
	  
      vec12(1) = hex8_vertex(1,2) - hex8_vertex(1,1)
      vec12(2) = hex8_vertex(2,2) - hex8_vertex(2,1)
      vec12(3) = hex8_vertex(3,2) - hex8_vertex(3,1)

      vec14(1) = hex8_vertex(1,4) - hex8_vertex(1,1)
      vec14(2) = hex8_vertex(2,4) - hex8_vertex(2,1)
      vec14(3) = hex8_vertex(3,4) - hex8_vertex(3,1)

      vec15(1) = hex8_vertex(1,5) - hex8_vertex(1,1)
      vec15(2) = hex8_vertex(2,5) - hex8_vertex(2,1)
      vec15(3) = hex8_vertex(3,5) - hex8_vertex(3,1)

      call cross_product(vec12,vec14,vec1,AA) 
      dot_prod = vec1(1)*vec15(1) + vec1(2)*vec15(2) + vec1(3)*vec15(3)
  
      if(dot_prod.gt.0.0) then
       ifnonrighthand = .FALSE.
      else
       ifnonrighthand = .TRUE.
       !write(6,*) 'non-right hand element detected'
      endif  

      return
      end
!-----------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine fix_if_non_right(iel)
! fix non-right hand element
! 1 <->19, 3 <-> 21, 7<->25, 9 <->27 
      use SIZE
      integer iel,iver
      real*8 xm2(27),ym2(27),zm2(27)
      character(3) cbc5,cbc6
      real bc55,bc56

! swap vertex
      do iver = 1,27
       xm2(iver) = xm1(iver,1,1,iel)
       ym2(iver) = ym1(iver,1,1,iel)
       zm2(iver) = zm1(iver,1,1,iel)
      enddo

      do iver = 1,9
       xm1(iver,1,1,iel) = xm2(iver+18)
       ym1(iver,1,1,iel) = ym2(iver+18)
       zm1(iver,1,1,iel) = zm2(iver+18)
      enddo
     
      do iver = 19,27
       xm1(iver,1,1,iel) = xm2(iver-18)
       ym1(iver,1,1,iel) = ym2(iver-18)
       zm1(iver,1,1,iel) = zm2(iver-18)
      enddo

! swap face 5 <-> 6

      cbc5 = cbc(5,iel)
      bc55 = bc(5,5,iel)
 
      cbc6 = cbc(6,iel)
      bc56 = bc(5,6,iel)

      cbc(5,iel)   = cbc6
      bc (5,5,iel) = bc56

      cbc(6,iel)   = cbc5
      bc (5,6,iel) = bc55 

      return
      end
!--------------------------------------------------------------------
      FUNCTION VOLUM0(P1,P2,P3,P0)
!
!                           3
!     Given four points in R , (P1,P2,P3,P0), VOLUM0 returns
!     the volume enclosed by the parallelagram defined by the
!     vectors { (P1-P0),(P2-P0),(P3-P0) }.  This routine has
!     the nice feature that if the 3 vectors so defined are
!     not right-handed then the volume returned is negative.

      REAL*8 P1(3),P2(3),P3(3),P0(3)

         U1=P1(1)-P0(1)
         U2=P1(2)-P0(2)
         U3=P1(3)-P0(3)

         V1=P2(1)-P0(1)
         V2=P2(2)-P0(2)
         V3=P2(3)-P0(3)

         W1=P3(1)-P0(1)
         W2=P3(2)-P0(2)
         W3=P3(3)-P0(3)

         CROSS1 = U2*V3-U3*V2
         CROSS2 = U3*V1-U1*V3
         CROSS3 = U1*V2-U2*V1

         VOLUM0  = W1*CROSS1 + W2*CROSS2 + W3*CROSS3
         
      RETURN
      END