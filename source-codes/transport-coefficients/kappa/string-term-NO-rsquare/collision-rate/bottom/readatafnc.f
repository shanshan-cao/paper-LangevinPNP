      subroutine readdata

      implicit none 

      include 'coms.in'

      integer i,j
      integer dummy_i, dummy_j

      open(unit=90,file="../../MM.dat")

      do i=1,101
         do j=1,1000000
              read(90,*) dummy_i, dummy_j, A2q(i,j)
         enddo
      enddo

      write(*,*) "potential ", A2q(1,1000000)


      end subroutine
