      subroutine read_lsq_kappa_table

      implicit none
      include 'df_coms.f'

      integer i,j
      character*10 dummy_c
      integer dummy_n
      double precision dummy_t,dummy_p
      
      open(unit=10,file='kappa_tp.dat',
     &        status='old',form='formatted')

      i=0
      do while(i.lt.gamma_nT)
         j=1
         do while(j.le.gamma_np)
            read(unit=10,fmt=*,err=2098,end=2097)
     &   tabkapT(i,j),dummy_p,tabkaptp(i,j),tabkapfn(i,j)

c      write(6,*) i,j,tabkaptp(i,j),tabkapfn(i,j)

            j=j+1
         enddo
         i=i+1
      enddo


 101  format(I4,2X,I4,2(2X,D12.6))

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

      subroutine get_lsq_kappa(dragA_value,kappaT_value,kappaL_value,
     &                    temperature,momentum,energhq)

      implicit none
      include 'df_coms.f'

      integer i,j
      double precision temperature,momentum,energhq
      double precision kappaT_value,kappaL_value,dragA_value

      i=int((temperature-0.12d0)/delta_Te+0.5d0) ! for temperature
      j=int((momentum-0d0)/delta_p+0.5d0)    ! for momentum

      if(i.lt.0) i=0
      if(i.gt.gamma_nT-1) i=gamma_nT-1
      if(j.lt.1) j=1
      if(j.gt.gamma_np) j=gamma_np

      kappaT_value = tabkaptp(i,j)
      kappaL_value = tabkaptp(i,j)
      dragA_value = (tabkaptp(i,j)-temperature*tabkapfn(i,j))/
     &  (2.0*temperature*energhq)

c      write(6,*) i,j,temperature,momentum,tabkaptp(i,j),tabkapfn(i,j)




      end subroutine





           

