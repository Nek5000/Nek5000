      program quick
c
      x = 2.4
      y = 0.8
c
      write(6,*) 'input scale:'
      read (5,*) s
      write(6,*) 'scale:',s
c
      write(6,10) x
      do i=1,11
         x = x+y
         y = s*y
         y = min(1.8,y)
         write(6,10) x
      enddo
   10 format(f16.5)
      stop
      end
