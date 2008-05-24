      program quick
c
      x = 2.5
      y = 0.7
c
      write(6,*) 'input scale:'
      read (5,*) s
      write(6,*) 'scale:',s
c
      write(6,10) x
      do i=1,7
         x = x+y
         y = s*y
         y = min(4.0,y)
         write(6,10) x
      enddo
   10 format(f16.5)
      stop
      end
