
      program main

      implicit none

      integer nt,np,i,j
      parameter (nt=30,np=50)
      double precision te_st(nt,np),p1_st(nt,np),e1_st(nt,np),
     &  kap_st(nt,np)
      double precision te_yu(nt,np),p1_yu(nt,np),e1_yu(nt,np),
     &   kap_yu(nt,np)
      double precision te_tot(nt,np),p1_tot(nt,np),e1_tot(nt,np),
     &   kap_tot(nt,np) 
      double precision fnc_kap
      double precision pi, KPfactor,runKT,preKT,Kfactor

      pi = 3.1415926535d0
      KPfactor=1.0d0
      preKT=0.27/0.3
      do i=1,nt
        do j=1,np
          te_st(i,j)=0.0d0
          p1_st(i,j)=0.0d0          
          e1_st(i,j)=0.0d0
          kap_st(i,j)=0d0
          te_yu(i,j)=0.0d0
          p1_yu(i,j)=0.0d0
          e1_yu(i,j)=0.0d0
          kap_yu(i,j)=0.0d0
          te_tot(i,j)=0.0d0
          p1_tot(i,j)=0.0d0
          e1_tot(i,j)=0.0d0
          kap_tot(i,j)=0.0d0
        end do
      end do
      fnc_kap=0.0d0




ccccccccccccccccccccccccccccccccccccccccccccccccc
c read Te p1,e1,kappa of string and yukawa

      open(20,file='kappa_tp-y.dat',status='replace')
      open(21,file='kappa-yukawanew.dat',status='replace')

      OPEN(36,FILE='kappa-string.dat',STATUS='OLD') 
         DO I=1,nt
           Do J=1,np
             READ(36,*) te_st(I,J),p1_st(I,J),e1_st(I,J),kap_st(I,J)
           END DO
         END DO
      CLOSE(36)

      OPEN(36,FILE='kappa-yukawa.dat',STATUS='OLD') 
         DO I=1,nt
           Do J=1,np
             READ(36,*) te_yu(I,J),p1_yu(I,J),e1_yu(I,J),kap_yu(I,J)
           END DO
         END DO
      CLOSE(36)

c            KPfactor=1.0d0
c            runKT=4.0*pi/9.0/log(2.0*e1*Te/0.04)/0.3
c            preKT=0.27/0.3
c            Kfactor=KPfactor*runKT*preKT
c            qhat_over_T3=qhat_over_T3*Kfactor
c            kappa=kappa*Kfactor


c         DO I=1,nt
c           Do J=1,np
c            runKT=4.0*pi/9.0/log(2.0*e1_yu(I,J)*te_yu(I,J)/0.04)/0.3
c            Kfactor=KPfactor*runKT*preKT
c            kap_yu(I,J)=kap_yu(I,J)*Kfactor
c            write(21,104) te_yu(i,j),p1_yu(i,j),e1_yu(i,j),kap_yu(i,j)  
c           END DO
c         END DO


         DO I=1,nt
           Do J=1,np
             te_tot(I,J)=te_yu(I,J)
             p1_tot(I,J)=p1_yu(I,J)
             e1_tot(I,J)=e1_yu(I,J)
             kap_tot(I,J)=kap_st(I,J)+kap_yu(I,J)
c             kap_tot(I,J)=kap_yu(I,J)
           END DO
         END DO

       do i=1,nt
         do j=1,np
           if(j.eq.1) then
              fnc_kap=(kap_tot(i,j+1)-kap_tot(i,j))/
     &        (e1_tot(i,j+1)-e1_tot(i,j))
           else
             if(j.eq.np) then
              fnc_kap=(kap_tot(i,np)-kap_tot(i,np-1))/
     &        (e1_tot(i,np)-e1_tot(i,np-1))      
             else
              fnc_kap=((kap_tot(i,j)-kap_tot(i,j-1))/
     &        (e1_tot(i,j)-e1_tot(i,j-1))+(kap_tot(i,j+1)
     &       -kap_tot(i,j))/(e1_tot(i,j+1)-e1_tot(i,j)))/2d0
             end if
           end if
             write(20,103) te_tot(i,j),p1_tot(i,j),kap_tot(i,j),fnc_kap    

         end do        
       end do


      close(20)
      close(21)
 103  format(E12.6,3(2X,E12.6))
 104  format(E12.6,3(2X,E12.6))
      end
         
