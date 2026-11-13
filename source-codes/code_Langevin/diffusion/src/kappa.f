      subroutine read_kappa_table

      implicit none
      include 'df_coms.f'

      integer i,j
      character*10 dummy_c
      integer dummy_n
      double precision dummy_f
      
      open(unit=10,file='FPcoeff.dat',
     &        status='old',form='formatted')

      i=0
      do while(i.le.gamma_nT)
         j=0
         do while(j.le.gamma_np)
            read(unit=10,fmt=*,err=2098,end=2097)
     &      dummy_f,dummy_f,tabDragA(i,j),tabKappaT(i,j),tabKappaL(i,j)
c            write(6,*) i,j,qhat_over_T3(i,j)
            j=j+1
         enddo
         i=i+1
      enddo
      write(6,*) "kappa table: ",tabDragA(gamma_nT,gamma_np),
     &      tabKappaT(gamma_nT,gamma_np),tabKappaL(gamma_nT,gamma_np)

 101  format(I4,2X,I4,4(2X,D12.6))

      close(10)
      return

 2097 continue
      write(6,*) 'ERROR: EOF reached in qhat table'
      write(6,*) 'terminating ...'
      stop
      return

 2098 continue
      write(6,*) 'READ-ERROR in qhat table'
      write(6,*) 'terminating ...'
      stop
      return

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_kappa(dragA_value,kappaT_value,kappaL_value,
     &                    temperature,momentum)

      implicit none
      include 'df_coms.f'

      integer i,j
      double precision temperature,momentum
      double precision dragA_value,kappaT_value,kappaL_value

      i=int((temperature-0.16d0)/delta_Te+0.5d0)  ! for temperature
      j=int((momentum-0d0)/delta_p+0.5d0) ! for momentum

      if(i.lt.0) i=0
      if(i.gt.gamma_nT) i=gamma_nT
      if(j.lt.0) j=0
      if(j.gt.gamma_np) j=gamma_np

      dragA_value=tabDragA(i,j)
      kappaT_value=2d0*tabKappaT(i,j)
      kappaL_value=2d0*tabKappaL(i,j)

      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           
      subroutine read_mass_table

      implicit none
      include 'df_coms.f'

      integer i
      double precision dummy_f
      
      open(unit=10,file='massT.dat',
     &        status='old',form='formatted')

      i=0
      do while(i.le.mass_nT)
         read(unit=10,fmt=*,err=2088,end=2087) dummy_f,massT(i)
         i=i+1
      enddo

      close(10)
      return

 2087 continue
      write(6,*) 'ERROR: EOF reached in qhat table'
      write(6,*) 'terminating ...'
      stop
      return

 2088 continue
      write(6,*) 'READ-ERROR in qhat table'
      write(6,*) 'terminating ...'
      stop
      return

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getMass(temperature,newMass)

      implicit none
      include 'df_coms.f'

      integer i
      double precision temperature,newMass

      i=int((temperature-0.16d0)/delta_massT+0.5d0)  ! for temperature

      if(i.lt.0) i=0
      if(i.gt.mass_nT) i=mass_nT

      newMass=massT(i)

      end subroutine

           

