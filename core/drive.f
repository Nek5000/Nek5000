      program NEKTON
c
      call nek_init(intracomm)
      call nek_solve()
      call nek_end()

      call exitt0()

      end
