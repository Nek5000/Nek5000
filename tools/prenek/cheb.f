      program quick_cheb
c
      write(6,*) 'input nbox:'
      read (5,*) n
      write(6,*) 'nbox:',n
c
      one = 1.
      pi  = 4.*atan(one)
      dt  = pi/n
      do i=n,0,-1
         t = i*dt
         y = cos(t)
         write(6,10) y
      enddo
   10 format(f16.5)
      stop
      end
