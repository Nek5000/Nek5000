C------------------------------------------------------------------------
C NEK5000: The Open Source Spectral Element Solver
C COPYRIGHT (c) 2008 UCHICAGO ARGONNE, LLC
C
C The UChicago Argonne, LLC as Operator of Argonne National
C Laboratory holds copyright in the Software. The copyright holder
C reserves all rights except those expressly granted to licensees,
C and U.S. Government license rights.
C
C License
C
C    NEK5000 is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    NEK5000 is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with NEK5000.  If not, see <http://www.gnu.org/licenses/>.
C
C--------------------------------------------------------------------------

      program NEKTON
c
c     top level driver
c
      call nek_init()
      call nek_solve()
      call nek_end()

      call exitt()

      end
