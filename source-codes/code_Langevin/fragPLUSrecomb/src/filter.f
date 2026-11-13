c**********************************************************************c
c                                                                      c
      program filter
c                                                                      c
c     1st version for HQ fragmentation: 08/18/2009 David Meidinger     c
c     2nd version for HQ frag.+recomb.: 05/21/2013 Shanshan Cao        c
c                                                                      c
c**********************************************************************c

      include 'common.f'
      logical isheavymeson
      logical isparton                 
      character*5 fortranunit
      character*77 filename
      integer i,j,k,rand_seed,numJudge
      double precision rlu,plength,recomb_prob
      double precision recomb_prob_br(1:nChannel)
      double precision recomb_random,temp_cell
      integer ind_l,ind_p,ind_temp,ind_beta,iChannel

      double precision max_ud(0:diffTemp*2),max_s(0:diffTemp*2)

      max_ud = (/0.063d0, 0.064d0, 0.065d0, 0.067d0, 0.068d0, 0.069d0,
     &           0.071d0, 0.072d0, 0.073d0, 0.075d0, 0.076d0, 0.077d0,
     &           0.079d0, 0.080d0, 0.081d0, 0.083d0, 0.084d0, 0.086d0,
     &           0.087d0, 0.088d0, 0.090d0, 0.091d0, 0.093d0, 0.094d0,
     &           0.096d0, 0.097d0, 0.099d0, 0.100d0, 0.102d0, 0.103d0,
     &           0.105d0, 0.107d0, 0.108d0, 0.110d0, 0.111d0, 0.113d0,
     &           0.115d0, 0.116d0, 0.118d0, 0.120d0, 0.121d0, 0.123d0,
     &           0.125d0, 0.126d0, 0.128d0, 0.130d0, 0.132d0, 0.133d0,
     &           0.135d0, 0.137d0, 0.139d0, 0.141d0, 0.142d0, 0.144d0,
     &           0.146d0, 0.148d0, 0.150d0, 0.152d0, 0.154d0, 0.155d0,
     &           0.157d0, 0.159d0, 0.161d0, 0.163d0, 0.165d0, 0.167d0,
     &           0.169d0, 0.171d0, 0.173d0, 0.175d0, 0.177d0
     &          /)

      max_s = (/0.063d0, 0.064d0, 0.065d0, 0.067d0, 0.068d0, 0.069d0,
     &           0.071d0, 0.072d0, 0.073d0, 0.075d0, 0.076d0, 0.077d0,
     &           0.079d0, 0.080d0, 0.081d0, 0.083d0, 0.084d0, 0.086d0,
     &           0.087d0, 0.088d0, 0.090d0, 0.091d0, 0.093d0, 0.094d0,
     &           0.096d0, 0.097d0, 0.099d0, 0.100d0, 0.102d0, 0.103d0,
     &           0.105d0, 0.107d0, 0.108d0, 0.110d0, 0.111d0, 0.113d0,
     &           0.115d0, 0.116d0, 0.118d0, 0.120d0, 0.121d0, 0.123d0,
     &           0.125d0, 0.126d0, 0.128d0, 0.130d0, 0.132d0, 0.133d0,
     &           0.135d0, 0.137d0, 0.139d0, 0.141d0, 0.142d0, 0.144d0,
     &           0.146d0, 0.148d0, 0.150d0, 0.152d0, 0.154d0, 0.155d0,
     &           0.157d0, 0.159d0, 0.161d0, 0.163d0, 0.165d0, 0.167d0,
     &           0.169d0, 0.171d0, 0.173d0, 0.175d0, 0.177d0
     &          /)


c      max_s = (/0.049d0, 0.050d0, 0.051d0, 0.053d0, 0.054d0, 0.055d0,
c     &          0.056d0, 0.057d0, 0.058d0, 0.059d0, 0.061d0, 0.062d0,
c     &          0.063d0, 0.064d0, 0.065d0, 0.067d0, 0.068d0, 0.069d0,
c     &          0.070d0, 0.072d0, 0.073d0, 0.074d0, 0.076d0, 0.077d0,
c     &          0.078d0, 0.080d0, 0.081d0, 0.083d0, 0.084d0, 0.085d0,
c     &          0.087d0, 0.088d0, 0.090d0, 0.091d0, 0.093d0, 0.094d0,
c     &          0.096d0, 0.097d0, 0.099d0, 0.100d0, 0.102d0, 0.103d0,
c     &          0.105d0, 0.106d0, 0.108d0, 0.110d0, 0.111d0, 0.113d0,
c     &          0.115d0, 0.116d0, 0.118d0, 0.120d0, 0.121d0, 0.123d0,
c     &          0.125d0, 0.126d0, 0.128d0, 0.130d0, 0.132d0, 0.133d0,
c     &          0.135d0, 0.137d0, 0.139d0, 0.141d0, 0.142d0, 0.144d0,
c     &          0.146d0, 0.148d0, 0.150d0, 0.152d0, 0.154d0
c     &         /)

      rand_seed=-1
      call sseed(rand_seed)

      call parameterRead

      numOut = 0

c debug
      if (hadr_flag.ne.1.and.hadr_flag.ne.2.and.hadr_flag.ne.3) then
         write(6,*) "Wrong value for hadr_flag ..."
         write(6,*) "Terminating ..."
         stop
      endif

c if necessary, open files
      if (INU.NE.5) then
         write (fortranunit,fmt='("ftn"I2.2)') INU
         call getenv(fortranunit,filename)
         open (unit=INU,file=filename,status='OLD',err=601)
      endif
      if (OUTU.NE.6) then
         write (fortranunit,fmt='("ftn"I2.2)') OUTU
         call getenv(fortranunit,filename)
         open (unit=OUTU,file=filename,err=602)
      endif
      if (ERRU.NE.6) then
         write (fortranunit,fmt='("ftn"I2.2)') ERRU
         call getenv(fortranunit,filename)
         open (unit=ERRU,file=filename,err=603)
      endif
      if (ANALYSIS) then
         write (fortranunit,fmt='("ftn"I2.2)') ANAU
         call getenv(fortranunit,filename)
         open (unit=ANAU,file=filename,err=604)
         call a_init
      endif

      call initpythia 

      call checkformat
c      fformat=19

      if(hadr_flag.eq.3) call read_recomb_table

      if (fformat.EQ.19) goto 19
      if (fformat.EQ.13) goto 13



c**********************************************************************c
c**** MAIN LOOP OSCAR 1997A *******************************************c


  19  continue
c event input 
         call readeventoscar

         if (ANALYSIS) call a_addevent_in
         remsamp=NUMSAMP

c recombination/fragmentation/decay

191      continue
         ipout=1
         count_recomb=0
         frag=.false.

         if(hadr_flag.eq.1) then ! pure fragmentation
            recomb_prob=0d0
         elseif(hadr_flag.eq.2) then ! pure recombination to D0 or D+
            recomb_prob=1d0          
            recomb_prob_br(1)=1d0
c            recomb_BR2=0d0
c            recomb_BR1=1d0
         endif
        
         ipin = 1
         do while (ipin.le.eh_npart)
c check if particle has to be decayed
            !write(6,*) ipin
            if ((HEAVYMESONDEC.AND.isheavymeson()).OR.isparton()) then
               recombORfrag=1
               judge_HM=0

c rescale thermal mass for recombination
               if(rescale_m.eq.1) then
                  original_mass=p_mass(ipin)
                  p_mass(ipin)=m_HQ
                  p_p0(ipin)=sqrt(p_px(ipin)**2+p_py(ipin)**2
     &                       +p_pz(ipin)**2+p_mass(ipin)**2)
               end if

c boost to the local rest frame of cell for recombination
               call rotbos(0d0,0d0,-c_vx(ipin),-c_vy(ipin),-c_vz(ipin),
     &                    p_px(ipin),p_py(ipin),p_pz(ipin),p_p0(ipin))
               temp_cell = Thydro(ipin)
               ind_temp = int((temp_cell-0.160d0)/0.002d0+0.5d0)
               temp_common = 0.160d0+ind_temp*0.002d0
               if(ind_temp.gt.diffTemp) ind_temp=diffTemp
               if(ind_temp.lt.0) recomb_prob=0d0

c               write(6,*) beta_cell,temp_cell,ind_temp,temp_common

               if(hadr_flag.eq.3.and.ind_temp.ge.0) then ! frag. + recomb.
                  plength=sqrt(p_px(ipin)**2+p_py(ipin)**2+
     &                    p_pz(ipin)**2)
                  ind_p=int(plength/prob_int)
                  ind_l=ind_p+1
                  if(ind_l.ge.prob_num) then ! this is the correct version
                     recomb_prob=0d0
                  else
                     iChannel=1
                     do while(iChannel.le.nChannel)
                       recomb_prob_br(iChannel)=prob_c(ind_l,iChannel)-
     &                 (prob_c(ind_l,iChannel)-prob_c(ind_l+1,iChannel))
     &                 *(plength-ind_p*prob_int)/prob_int
                       iChannel=iChannel+1
                     end do
                     max_3body_s=wigner_max(ind_l,5)
                     max_3body_p=wigner_max(ind_l,6)
                     recomb_prob=recomb_prob_br(nChannel)
                  endif
               endif

c 1/6/2015 the code is modified so that every c goes to D
c pay attention to boost back and forth
c               numJudge = 0
c               do while(judge_HM.eq.0)      
               
c               numJudge = numJudge+1
c               if(numJudge.gt.1D6) then
c                  write(6,*) "Large numJudge for ipin = ",ipin,
c     &                   beta_cell,recomb_BR1,recomb_random
c               endif

               recomb_random=rlu(0)

               if(recomb_random.lt.recomb_prob) then ! do recomb
                 ! big modification on 09/18/2019
                 if(HQid.eq.4) then ! charm 
                  if(recomb_random.lt.recomb_prob_br(1)) then ! D0 and D+
                     mxvalue = max_ud(ind_temp*2)
                     if(rlu(0).lt.0.5d0) then ! D0
                         id_hadF = 421
                         mass_hadF = 1.86d0
                         call force_recomb2(recombORfrag,1,0,1.86d0,1d0,
     &                                   id_hadF,mass_hadF)
                     else ! D+
                         id_hadF = 411
                         mass_hadF = 1.87d0                     
                        call force_recomb2(recombORfrag,1,0,1.87d0,0d0,
     &                                   id_hadF,mass_hadF)
                     endif
                  else if(recomb_random.lt.recomb_prob_br(2)) then ! D^*0 and D^*+
                     mxvalue = max_ud(ind_temp*2)
                     if(rlu(0).lt.0.5d0) then ! D^*0
                         id_hadF = 421
                         mass_hadF = 1.86d0                       
                         call force_recomb2(recombORfrag,1,0,2.01d0,1d0,
     &                                   id_hadF,mass_hadF)                         
                     else ! D^*+
c                        if(rlu(0).lt.0.68d0) then 
c                          id_hadF = 421
c                          mass_hadF = 1.86d0                      
c                        else
                          id_hadF = 411
                          mass_hadF = 1.87d0  
c                         endif                      
                          call force_recomb2(recombORfrag,1,0,2.01d0,
     &                    0.68d0,id_hadF,mass_hadF)   
                     endif
               
                  else if(recomb_random.lt.recomb_prob_br(3)) then ! D_1(2420)
                     mxvalue = max_ud(ind_temp*2)
                     id_hadF = 421
                     mass_hadF = 1.86d0  
                     call force_recomb2(recombORfrag,1,1,2.42d0,0.5d0,
     &                                   id_hadF,mass_hadF)
                  else if(recomb_random.lt.recomb_prob_br(4)) then ! D_0^*(2300)
                     mxvalue = max_ud(ind_temp*2)
                     id_hadF = 421
                     mass_hadF = 1.86d0                      
                     call force_recomb2(recombORfrag,1,1,2.30d0,0.5d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(5)) then ! D_1^*(2440)
                     mxvalue = max_ud(ind_temp*2)
c                        if(rlu(0).lt.0.68d0) then 
c                          id_hadF = 421
c                          mass_hadF = 1.86d0                      
c                        else
                          id_hadF = 411
                          mass_hadF = 1.87d0  
c                        endif                   
                     call force_recomb2(recombORfrag,1,1,2.44d0,0.5d0,
     &                                   id_hadF,mass_hadF)                      
                  else if(recomb_random.lt.recomb_prob_br(6)) then ! D_2^*(2420)
                     mxvalue = max_ud(ind_temp*2)
c                        if(rlu(0).lt.0.68d0) then 
c                          id_hadF = 421
c                          mass_hadF = 1.86d0                      
c                        else
                          id_hadF = 411
                          mass_hadF = 1.87d0  
c                        endif                      
                     call force_recomb2(recombORfrag,1,1,2.46d0,0.5d0,
     &                                   id_hadF,mass_hadF)  
                  else if(recomb_random.lt.recomb_prob_br(7)) then ! Ds (1968)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 431
                     mass_hadF = 2.0d0                       
                     call force_recomb2(recombORfrag,3,0,1.97d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(8)) then ! Ds* (2122)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 431
                     mass_hadF = 2.0d0                      
                     call force_recomb2(recombORfrag,3,0,2.12d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(9)) then ! Ds* (2536)
                     mxvalue = max_s(ind_temp*2)
                     if(rlu(0).lt. 0.5d0) then !D0
                        id_hadF = 421
                        mass_hadF = 1.86d0                      
                     else !D+
                        id_hadF = 411
                        mass_hadF = 1.87d0  
                     endif                       
                     call force_recomb2(recombORfrag,3,1,2.54d0,1d0,
     &                                    id_hadF,mass_hadF) 
                  else if(recomb_random.lt.recomb_prob_br(10)) then ! Ds_0 (2317)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 431
                     mass_hadF = 2.0d0                        
                     call force_recomb2(recombORfrag,3,1,2.32d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(11)) then ! Ds_1 (2460)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 431
                     mass_hadF = 2.0d0                      
                     call force_recomb2(recombORfrag,3,1,2.46d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(12)) then ! Ds_2 (2573)
                     mxvalue = max_s(ind_temp*2)
                     if(rlu(0).lt. 0.5d0) then !D0
                        id_hadF = 421
                        mass_hadF = 1.86d0                      
                     else !D+
                        id_hadF = 411
                        mass_hadF = 1.87d0  
                     endif                      
                     call force_recomb2(recombORfrag,3,1,2.57d0,1d0,
     &                                   id_hadF,mass_hadF)                      
                  else if(recomb_random.lt.recomb_prob_br(13)) then 
                     ! Lambda s wave 1/2+
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_ud(ind_temp*2)
                     id_hadF =4122
                     mass_hadF =2.286d0
                     call force_recomb3(recombORfrag,1,1,0,2.29d0,
     &                                   id_hadF,mass_hadF)
                  else if(recomb_random.lt.recomb_prob_br(14)) then
                     ! Lambda s wave 3/2+
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_ud(ind_temp*2) 
                     id_hadF =4122
                     mass_hadF =2.286d0
                     call force_recomb3(recombORfrag,1,1,0,2.86d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(15)) then
                     ! Lambda p wave together
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_ud(ind_temp*2) 
                     id_hadF = 4122
                     mass_hadF = 2.286d0                     
                     call force_recomb3(recombORfrag,1,1,1,2.60d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(16)) then 
c                     ! Sigma s wave 1/2+   ! Sigma++, Sigma+, and Sigma0
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_ud(ind_temp*2)
                     id_hadF = 4122
                     mass_hadF = 2.286d0     
                     call force_recomb3(recombORfrag,1,1,0,2.46d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(17)) then
c                     ! Sigma s wave 3/2+   ! Sigma*++, Sigma*+, and Sigma*0
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_ud(ind_temp*2)
                     id_hadF = 4122
                     mass_hadF = 2.286d0 
                     call force_recomb3(recombORfrag,1,1,0,2.520d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(18)) then
c                     ! Sigma p wave together
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_ud(ind_temp*2) 
                     id_hadF = 4122
                     mass_hadF = 2.286d0
                     call force_recomb3(recombORfrag,1,1,1,2.80d0,
     &                                   id_hadF,mass_hadF)
cccccccccccccccccccc add Cascade prob_coal_sp   (lisq)ccccccccccccccccccccc 
                  else if(recomb_random.lt.recomb_prob_br(19)) then
c                     ! Cascade s wave 1/2+  ! Cascade+, Cascade0, Cascade'+, Cascade'0
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_s(ind_temp*2) 
                     id_hadF = 4132
                     mass_hadF = 2.470d0
                     call force_recomb3(recombORfrag,1,3,0,2.47d0,
     &                                   id_hadF,mass_hadF)
                  else if(recomb_random.lt.recomb_prob_br(20)) then
c                     ! Cascade s wave 3/2+  ! Cascade*+, Cascade*0, Cascade'*+, Cascade'*0
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_s(ind_temp*2)
                     id_hadF = 4132
                     mass_hadF = 2.470d0                     
                     call force_recomb3(recombORfrag,1,3,0,2.645d0,
     &                                   id_hadF,mass_hadF)
                  else if(recomb_random.lt.recomb_prob_br(21)) then
c                     ! Cascade p wave together
                     mxvalue1 = max_ud(ind_temp*2) 
                     mxvalue2 = max_s(ind_temp*2)
c                     ！Cascade_c^0
                     id_hadF = 4132    
                     mass_hadF =  2.470d0                      
                     call force_recomb3(recombORfrag,1,3,1,2.815d0,
     &                                   id_hadF,mass_hadF)
cccccccccccccccccccc add Omega prob_coal_sp  (lisq) ccccccccccccccccccccc 
                  else if(recomb_random.lt.recomb_prob_br(22)) then
c                     ! Omega s wave 1/2+   ! Omega
                     mxvalue1 = max_s(ind_temp*2) 
                     mxvalue2 = max_s(ind_temp*2) 
c                     ！Omega_c^0
                     id_hadF = 4332    
                     mass_hadF = 2.695d0                    
                     call force_recomb3(recombORfrag,3,3,0,2.6950d0,
     &                                   id_hadF,mass_hadF)
                  else if(recomb_random.lt.recomb_prob_br(23)) then
c                     ! Omega s wave 3/2+  ! Omega*
                     mxvalue1 = max_s(ind_temp*2) 
                     mxvalue2 = max_s(ind_temp*2) 
                     id_hadF = 4332
                     mass_hadF = 2.695d0                   
                     call force_recomb3(recombORfrag,3,3,0,2.77d0,
     &                                   id_hadF,mass_hadF)
                  else if(recomb_random.lt.recomb_prob_br(24)) then
c                     ! Omega p wave together
                     mxvalue1 = max_s(ind_temp*2) 
                     mxvalue2 = max_s(ind_temp*2) 
                     id_hadF = 4332
                     mass_hadF = 2.695d0                    
                     call force_recomb3(recombORfrag,3,3,1,3.1d0,
     &                                   id_hadF,mass_hadF)
                  endif   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
                  else if(HQid.eq.5) then ! beauty quarks
                  if(recomb_random.lt.recomb_prob_br(1)) then ! B+ and B0
                     mxvalue = max_ud(ind_temp*2)
                     if(rlu(0).lt.0.5d0) then ! B+
                       id_hadF = 521
                       mass_hadF = 5.30d0                       
                        call force_recomb2(recombORfrag,1,0,5.28d0,1d0,
     &                                   id_hadF,mass_hadF)                         
                     else ! B0
                       id_hadF = 511
                       mass_hadF = 5.30d0                      
                        call force_recomb2(recombORfrag,1,0,5.28d0,0d0,
     &                                   id_hadF,mass_hadF)                          
                     endif
                  else if(recomb_random.lt.recomb_prob_br(2)) then ! B^*+ and B^*0
                     mxvalue = max_ud(ind_temp*2)
                     if(rlu(0).lt.0.5d0) then ! B^*+
                       id_hadF = 521
                       mass_hadF = 5.30d0                       
                         call force_recomb2(recombORfrag,1,0,5.32d0,
     &                    0.5d0,id_hadF,mass_hadF)
                     else ! B^*0
                       id_hadF = 511
                       mass_hadF = 5.30d0                      
                         call force_recomb2(recombORfrag,1,0,5.32d0,
     &                  0.5d0,id_hadF,mass_hadF)
                     endif
                  else if(recomb_random.lt.recomb_prob_br(3)) then ! B_1(5721)
                     mxvalue = max_ud(ind_temp*2)
                     id_hadF = 521
                      mass_hadF = 5.30d0 
                     call force_recomb2(recombORfrag,1,1,5.72d0,0.5d0,
     &                                   id_hadF,mass_hadF)                      
                  else if(recomb_random.lt.recomb_prob_br(4)) then ! B_0^*(?)
                     mxvalue = max_ud(ind_temp*2)
                     id_hadF = 511
                     mass_hadF = 5.30d0                     
                     call force_recomb2(recombORfrag,1,1,5.72d0,0.5d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(5)) then ! B_1^*(?)
                     mxvalue = max_ud(ind_temp*2)
                     id_hadF = 521
                      mass_hadF = 5.30d0                      
                     call force_recomb2(recombORfrag,1,1,5.72d0,0.5d0,
     &                                   id_hadF,mass_hadF)                      
                  else if(recomb_random.lt.recomb_prob_br(6)) then ! B_2^*(5747)
                     mxvalue = max_ud(ind_temp*2)
                     id_hadF = 521
                      mass_hadF = 5.30d0                      
                     call force_recomb2(recombORfrag,1,1,5.75d0,0.5d0,
     &                                   id_hadF,mass_hadF)                      
                  else if(recomb_random.lt.recomb_prob_br(7)) then ! Bs
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 531
                     mass_hadF = 5.40d0                      
                     call force_recomb2(recombORfrag,3,0,5.37d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(8)) then ! Bs*
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 531
                     mass_hadF = 5.40d0  
                     call force_recomb2(recombORfrag,3,0,5.42d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(9)) then ! Bs_0* (?)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 531
                     mass_hadF = 5.40d0                       
                     call force_recomb2(recombORfrag,3,1,5.42d0,1d0,
     &                                   id_hadF,mass_hadF)                     
                  else if(recomb_random.lt.recomb_prob_br(10)) then ! Bs_0 (?)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 531
                     mass_hadF = 5.40d0                       
                     call force_recomb2(recombORfrag,3,1,5.42d0,1d0,
     &                                   id_hadF,mass_hadF)                       
                  else if(recomb_random.lt.recomb_prob_br(11)) then ! Bs_1 (5830)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 531
                     mass_hadF = 5.40d0                      
                     call force_recomb2(recombORfrag,3,1,5.83d0,1d0,
     &                                   id_hadF,mass_hadF) 
                  else if(recomb_random.lt.recomb_prob_br(12)) then ! Bs_2 (5840)
                     mxvalue = max_s(ind_temp*2)
                     id_hadF = 531
                     mass_hadF = 5.40d0                      
                     call force_recomb2(recombORfrag,3,1,5.84d0,1d0,
     &                                   id_hadF,mass_hadF)                      
                !  else if(recomb_random.lt.recomb_prob_br(13)) then 
                !     ! Lambda s wave 1/2+
                !     mxvalue1 = max_ud(ind_temp*2) 
                !     mxvalue2 = max_ud(ind_temp*2) 
                !     call force_recomb3(recombORfrag,1,1,0,2.29d0)
                !  else if(recomb_random.lt.recomb_prob_br(14)) then
                !     ! Lambda s wave 3/2+
                !     mxvalue1 = max_ud(ind_temp*2) 
                !     mxvalue2 = max_ud(ind_temp*2) 
                !     call force_recomb3(recombORfrag,1,1,0,2.86d0)
                !  else if(recomb_random.lt.recomb_prob_br(15)) then
                !     ! Lambda p wave together
                !     mxvalue1 = max_ud(ind_temp*2) 
                !     mxvalue2 = max_ud(ind_temp*2) 
                !     call force_recomb3(recombORfrag,1,1,1,2.60d0)
                !  else if(recomb_random.lt.recomb_prob_br(16)) then
                !     ! Sigma s wave 1/2+
                !     mxvalue1 = max_ud(ind_temp*2) 
                !     mxvalue2 = max_ud(ind_temp*2) 
                !     call force_recomb3(recombORfrag,1,1,0,2.46d0)
                !  else if(recomb_random.lt.recomb_prob_br(17)) then
                !     ! Sigma s wave 3/2+
                !     mxvalue1 = max_ud(ind_temp*2) 
                !     mxvalue2 = max_ud(ind_temp*2) 
                !     call force_recomb3(recombORfrag,1,1,0,2.52d0)
                !  else if(recomb_random.lt.recomb_prob_br(18)) then
                !     ! Sigma p wave together
                !     mxvalue1 = max_ud(ind_temp*2) 
                !     mxvalue2 = max_ud(ind_temp*2) 
                !     call force_recomb3(recombORfrag,1,1,1,2.80d0)
                  endif
                 endif ! c or b

               else

c if not recombined, boost back to global c.m. frame for fragmentation
                  call rotbos(0d0,0d0,c_vx(ipin),c_vy(ipin),c_vz(ipin),
     &                      p_px(ipin),p_py(ipin),p_pz(ipin),p_p0(ipin))

c scale thermal mass back for recombination
                  if(rescale_m.eq.1) then
                     p_mass(ipin)=original_mass
                     p_p0(ipin)=sqrt(p_px(ipin)**2+p_py(ipin)**2
     &                          +p_pz(ipin)**2+p_mass(ipin)**2)
                  end if

                  call fragdecay
               endif

c               enddo ! (judge_HM is now set to 1)

               frag=.true.
            else
c if ONLYFRAG=2, skip unprocessed particles
               if (ONLYFRAG.NE.2) call copypartoscar
            endif

            ipin = ipin+1
          enddo
c 190     continue
c final number of particles in this event
         oeh_npart=ipout-1
c event output (suppress for unchanged events, when ONLYFRAG is not 0)
         if ((ONLYFRAG.EQ.0).OR.frag) then
            call writeeventoscar
            if (ANALYSIS) call a_addevent_out
         endif

c process events with frag=.true. NUMSAMP times
         if (frag.AND.(remsamp.NE.1)) then
            remsamp=remsamp-1
            write(6,*) remsamp,'samples left'
            goto 191
         endif

      goto 19
      
c**********************************************************************c
c**** MAIN LOOP URQMD FILE 13 *****************************************c

  13  continue
c event input 
         call readevent13
         if (ANALYSIS) call a_addevent_in
         remsamp=NUMSAMP
c fragmentation/decay
 131     continue
         ipout=1
         frag=.false.
         do 130 ipin=1,eh_npart
c check if particle has to be decayed
c (no need to check for partons, since URQMD won't produce them)
            if (HEAVYMESONDEC.AND.isheavymeson()) then
              call fragdecay   
              frag=.true.
            else
c if ONLYFRAG=2, skip unprocessed particles
              if (ONLYFRAG.NE.2) call copypart13
            endif
 130     continue  
c final number of particles in this event
         oeh_npart=ipout-1
c event output (suppress for unchanged events, when ONLYFRAG is not 0)
         if ((ONLYFRAG.EQ.0).OR.frag) then
            call writeevent13
            if (ANALYSIS) call a_addevent_out
         endif 
c process events with frag=.true. NUMSAMP times
         if (frag.AND.(remsamp.NE.1)) then
            remsamp=remsamp-1
            goto 131
         endif
      goto 13
      
c**********************************************************************c
c**** ERRORS while opening files: always stdout ***********************c

 601  write (6,*) 'ERROR while opening input file!'
      stop
 602  write (6,*) 'ERROR while opening output file'
      stop
 603  write (6,*) 'ERROR while opening error file'
      stop
 604  write (6,*) 'ERROR while opening analysis file'
      stop

 605  continue
      write(6,*) "program ends successfully:)"      

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_recomb_table

      include 'common.f'

      integer i,iChannel
      double precision dummy_pc,dummy_prob

      if(HQid.eq.4) then
         open(unit=11,file='recomb_c_raw.dat',
     &        status='old',form='formatted')
         open(unit=12,file='max_wigner_c.dat',
     &        status='old',form='formatted')
         prob_int=prob_int_c
         omega_HM=omega_cM
         omega_HB=omega_cB
         m_HQ=m_c
      elseif(HQid.eq.5) then
         open(unit=11,file='recomb_b_raw.dat',
     &        status='old',form='formatted')
         open(unit=12,file='max_wigner_b.dat',
     &        status='old',form='formatted')
         prob_int=prob_int_b
         omega_HM=omega_bM
         omega_HB=omega_bB
         m_HQ=m_b
      else
         write(*,*) "un-recognized heavy quark ID ..."
         stop
      endif
 
c      i=1
c      do while(i.le.prob_num)
c         read(unit=11,fmt=2196,err=2198,end=2197) dummy_pc,
c     &       prob_c(i,0:diffTemp)
c         read(unit=12,fmt=2196,err=2198,end=2197) dummy_pc,
c     &       DBR1(i,0:diffTemp)
cc         write(6,*) dummy_pc,prob_c(i,0),prob_c(i,1),prob_c(i,21)
c         i=i+1
c      enddo

      i=1
      do while(i.le.prob_num)
         read(unit=11,fmt=2196,err=2198,end=2197) dummy_pc,
     &       dummy_prob,prob_c_raw(i,1:nChannel)
         read(unit=12,fmt=2196,err=2198,end=2197) dummy_pc,
     &       wigner_max(i,1:10)
         ! calculate accumulative probability, ignore p wave channels
         prob_c(i,1)=prob_c_raw(i,1)
         iChannel=2
         do while(iChannel.le.nChannel)
c            if((iChannel.ge.3.and.iChannel.le.6).or.
c     &         (iChannel.ge.9.and.iChannel.le.12).or.
c     &         iChannel.eq.15.or.iChannel.eq.18.or.
c     &         iChannel.eq.21.or.iChannel.eq.24) then
c               prob_c(i,iChannel)=prob_c(i,iChannel-1)
c            else
              prob_c(i,iChannel)=prob_c(i,iChannel-1)
     &                        +prob_c_raw(i,iChannel)
c            endif
            iChannel=iChannel+1
         end do 
         if(prob_c(i,nChannel).gt.1d0) then
            iChannel=1
            do while(iChannel.le.nChannel)
               prob_c(i,iChannel)=prob_c(i,iChannel)/prob_c(i,nChannel)
               iChannel=iChannel+1
            end do
         end if

         write(6,*) dummy_pc,prob_c(i,24)
         i=i+1
      enddo


      write(6,*) "read in prob table successfully:)"

      return

c 2196 format(E12.4,55(2X,E12.4))
 2196 format(E10.4,55(2X,E10.4))

 2197 continue
      write(6,*) "end"
      return

 2198 continue
      write(6,*) "error"
      return

      end
