ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine force_recomb2(hadrflag,lightFlavor,orbital,
     &            massH,prob0,id_hadFO,mass_hadFO)
c do recombination of heavy-light quarks
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       include 'common.f' 
      
       integer hadrflag,i,loopcount,lightFlavor,orbital
       double precision rlu,phi,utheta,temp,plight
       double precision kmin,kmin2,plmax,plmaxcm,gacell,gacm
       double precision becell,becm,Elmax,Elmaxcm,pHQ,dum_random
       double precision reducedM0
       double precision massH,prob0,mass_hadFO
       integer id_hadFO
c generate one light parton according to thermal distribution
       temp=Thydro(ipin)
       if(temp.lt.0.1d0) return

       if(lightFlavor.eq.1) then
          mlight=m_ud
       elseif(lightFlavor.eq.3) then
          mlight=m_s
       else
          write(6,*) "Flavor not available!"
          stop
       endif

       plmax=15d0*temp
       Elmax=sqrt(plmax**2+mlight**2)
       becell=sqrt(c_vx(ipin)**2+c_vy(ipin)**2+c_vz(ipin)**2)
       gacell=1d0/sqrt(1d0-becell**2)
       plmaxcm=gacell*(plmax+becell*Elmax)
       Elmaxcm=gacell*(Elmax+becell*plmax)
       pHQ=sqrt(p_p0(ipin)**2-p_mass(ipin)**2)
       becm=(pHQ+plmaxcm)/(p_p0(ipin)+Elmaxcm)
       gacm=1d0/sqrt(1d0-becm**2)
       kmin=gacm*(pHQ-becm*p_p0(ipin))
       kmin2=kmin**2
       if(pHQ.lt.plmaxcm) kmin2=-1d0

c calculate parameters for Wigner function
       reducedM0=mlight*p_mass(ipin)/(mlight+p_mass(ipin))
       sigma0=1d0/sqrt(reducedM0*omega_HM)

       loopcount=0
       do while(hadrflag.eq.1)
          call judgehadmech2(orbital,hadrflag,kmin2,massH)
          loopcount=loopcount+1
          if(loopcount.gt.1D8) then
             write(6,*) "too many loops, no recomb"
             return
          endif
       enddo

c do recombination if hadrflag=2
       op_num(ipout)=ipout
c       op_px(ipout)=p_px(ipin)+px_light
c       op_py(ipout)=p_py(ipin)+py_light
c       op_pz(ipout)=p_pz(ipin)+pz_light

      
c       if(abs(p_id(ipin)).eq.4) then ! c quark
c          if(lightFlavor.eq.1) then
c             op_id(ipout)=4210 ! D0 or D+
c             op_mass(ipout)=1.9d0
c          elseif(lightFlavor.eq.3) then
c             op_id(ipout)=4310 ! Ds+
c             op_mass(ipout)=2.1d0
c          endif
c       else
c          write(6,*) "HF not available!"
c          stop
c       endif

       dum_random=rlu(0)

       if(abs(p_id(ipin)).eq.4) then ! c quark
          if(lightFlavor.eq.1) then
             if(dum_random.le.prob0) then
                op_id(ipout)= 421 ! D0
                op_mass(ipout)=1.86d0
             else
                op_id(ipout)=411 ! D+
                op_mass(ipout)=1.87d0
             endif
          else if(lightFlavor.eq.3) then
c             op_id(ipout)=431
c             op_mass(ipout)=2.0d0
             op_id(ipout)=id_hadFO
             op_mass(ipout)=mass_hadFO               
          else
              write(6,*) "Flavor not available!"
              stop
          endif
          
          if(p_id(ipin).lt.0) op_id(ipout)=-op_id(ipout)
       elseif(abs(p_id(ipin)).eq.5) then ! b quark
          if(lightFlavor.eq.1) then
             if(dum_random.le.prob0) then
                op_id(ipout)=521 ! B+
                op_mass(ipout)=5.3d0
             else
                op_id(ipout)=511 ! B0
                op_mass(ipout)=5.3d0
             endif
          else if(lightFlavor.eq.3) then
             op_id(ipout)=531
             op_mass(ipout)=5.4d0
          else
              write(6,*) "Flavor not available!"
              stop
          endif
          if(p_id(ipin).gt.0) op_id(ipout)=-op_id(ipout)
       else
          write(6,*) "Heavy meson not listed"
          stop
       endif

c in this version, recombine in 2-body c.m. first and then boost out
       op_px(ipout)=0d0
       op_py(ipout)=0d0
       op_pz(ipout)=0d0
       op_p0(ipout)=directM
!       op_p0(ipout)=op_mass(ipout)
!       op_p0(ipout)=massH

       call rotbos(0d0,0d0,vx_HM,vy_HM,vz_HM,
     &      op_px(ipout),op_py(ipout),op_pz(ipout),op_p0(ipout))

c       write(6,*) "1: ", directM,vx_HM,vy_HM,vz_HM,op_px(ipout),
c     &                   op_py(ipout),op_pz(ipout),op_p0(ipout)
c       call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
c     &      op_px(ipout),op_py(ipout),op_pz(ipout),op_p0(ipout))

c       op_p0(ipout)=sqrt(op_px(ipout)**2+op_py(ipout)**2+
c     &                  op_pz(ipout)**2+op_mass(ipout)**2)
       op_r0(ipout)=p_r0(ipin)
       op_rx(ipout)=p_rx(ipin)
       op_ry(ipout)=p_ry(ipin)
       op_rz(ipout)=p_rz(ipin) 
       oThydro(ipout)=Thydro(ipin)
       osdensity(ipout)=sdensity(ipin)
       oc_vx(ipout)=c_vx(ipin)
       oc_vy(ipout)=c_vy(ipin)
       oc_vz(ipout)=c_vz(ipin)
       op_wt(ipout)=p_wt(ipin)
       op_ipx(ipout)=p_ipx(ipin)
       op_ipy(ipout)=p_ipy(ipin)
       op_ipz(ipout)=p_ipz(ipin)
       call rotbos(0d0,0d0,c_vx(ipin),c_vy(ipin),c_vz(ipin),
     &            op_px(ipout),op_py(ipout),op_pz(ipout),op_p0(ipout))

       call decay2body(op_px(ipout),op_py(ipout),op_pz(ipout),
     &      op_p0(ipout),op_mass(ipout))

       ipout=ipout+1
       count_recomb=count_recomb+1
       judge_HM=1
       frag=.true.

       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine judgehadmech2(orbital,hadrflag,kminsquare,massH)
c judge the hadronization mechanism, frag. or recomb.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       include 'common.f'

       integer hadrflag,orbital
       double precision cmbetax,cmbetay,cmbetaz
c       double precision mxpt,mxvalue
       double precision qx_pair,qy_pair,qz_pair,qpair2,plight
       double precision probWigner,temp
       double precision rlu,utheta,phi
       double precision probmax,kminsquare,distribution,Egluon
       double precision massH,p0_HQ,px_HM,py_HM,pz_HM,p0_HM

c generate one light parton according to thermal distribution
       temp=Thydro(ipin)
c       temp = temp_common
c       temp=0.165d0

       if(orbital.eq.0) then
          if(kminsquare.lt.0) then
             probmax=1d0
          else
             probmax=exp(-kminsquare*sigma0**2)
          endif
       else if(orbital.eq.1) then 
          if(kminsquare*sigma0**2.lt.1d0) then
             probmax=exp(-1d0)
          else
             probmax=kminsquare*sigma0**2*exp(-kminsquare*sigma0**2)
          endif
       else
          write(6,*) "only s and p waves are availabe now"
          stop
       endif
 
c       mxpt=sqrt(2d0*temp**2+2d0*temp*sqrt(temp**2+mlight**2))
c       mxvalue=mxpt**2*exp(-sqrt(mlight**2+mxpt**2)/temp)
c       write(6,*) mxpt,mxvalue

c       mxvalue=0.07d0 ! for sampling efficiency, non-trivial
       plight=rlu(0)*15d0*temp
       p0_light=sqrt(plight**2+mlight**2)
       Egluon=sqrt(plight**2+m_g**2)
       distribution=6d0*plight**2/(exp(p0_light/temp)+1d0)+
     &              plight**2/(exp(Egluon/temp)-1d0)*16d0/6d0
       do while(rlu(0)*mxvalue.gt.distribution)
          plight=rlu(0)*15d0*temp
          p0_light=sqrt(plight**2+mlight**2)
          Egluon=sqrt(plight**2+m_g**2)
          distribution=6d0*plight**2/(exp(p0_light/temp)+1d0)+
     &              plight**2/(exp(Egluon/temp)-1d0)*16d0/6d0
       enddo

       utheta=rlu(0)*2d0-1d0
       phi=rlu(0)*2d0*pi
       px_light=plight*sqrt(1-utheta*utheta)*cos(phi)
       py_light=plight*sqrt(1-utheta*utheta)*sin(phi)
       pz_light=plight*utheta        

c       call rotbos(0d0,0d0,c_vx(ipin),c_vy(ipin),c_vz(ipin),
c     &            px_light,py_light,pz_light,p0_light)

cc if use pure recombination mechanism, i.e. fW = 1
c       if (hadr_flag.eq.2) then
c          hadrflag = 2
c          return
c       endif

c boost heavy/light quark to the c.m. frame of the heavy-light system

       cmbetax=(p_px(ipin)+px_light)/(p_p0(ipin)+p0_light)
       cmbetay=(p_py(ipin)+py_light)/(p_p0(ipin)+p0_light)
       cmbetaz=(p_pz(ipin)+pz_light)/(p_p0(ipin)+p0_light)

       call rotbos(0d0,0d0,-cmbetax,-cmbetay,-cmbetaz,
     &            px_light,py_light,pz_light,p0_light)
       call rotbos(0d0,0d0,-cmbetax,-cmbetay,-cmbetaz,
     &            p_px(ipin),p_py(ipin),p_pz(ipin),p_p0(ipin))

c calculate the relative momentum
       qx_pair=(p0_light*p_px(ipin)-p_p0(ipin)*px_light)/
     &         (p0_light+p_p0(ipin))
       qy_pair=(p0_light*p_py(ipin)-p_p0(ipin)*py_light)/
     &         (p0_light+p_p0(ipin))
       qz_pair=(p0_light*p_pz(ipin)-p_p0(ipin)*pz_light)/
     &         (p0_light+p_p0(ipin))
       qpair2=qx_pair*qx_pair+qy_pair*qy_pair+qz_pair*qz_pair

       if(orbital.eq.0) then ! s wave
          probWigner=exp(-qpair2*sigma0**2)  ! only part of the Wigner      
       else if(orbital.eq.1) then ! p wave
          probWigner=exp(-qpair2*sigma0**2)*qpair2*sigma0**2
       endif

       p0_HQ=p_p0(ipin)
c       probWigner=probWigner*(Gamma_HM*massH)**2/
c     &     ((Gamma_HM*massH)**2+((p0_HQ+p0_light)**2-massH**2)**2)

c use Monte Carlo to decide whether to recombine
       if(hadr_flag.eq.2.or.rlu(0)*probmax.lt.probWigner) then
          hadrflag=2

          directM=p0_HQ+p0_light
          vx_HM=cmbetax
          vy_HM=cmbetay
          vz_HM=cmbetaz

c          ! add/subtract a photon to conserve energy momentum
c          px_HM=0d0
c          py_HM=0d0
c          pz_HM=0d0
c          p0_HM=p0_HQ+p0_light
c          call decay2body(px_HM,py_HM,pz_HM,p0_HM,massH)
c          call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
c     &               px_HM,py_HM,pz_HM,p0_HM)
c          vx_HM=px_HM/p0_HM
c          vy_HM=py_HM/p0_HM
c          vz_HM=pz_HM/p0_HM
       endif 

c boost heavy/light quark back to the c.m. frame
       call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
     &            px_light,py_light,pz_light,p0_light)
       call rotbos(0d0,0d0,cmbetax,cmbetay,cmbetaz,
     &            p_px(ipin),p_py(ipin),p_pz(ipin),p_p0(ipin))

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine decay2body(px,py,pz,p0,m_final)

       ! m_initial decay to photon + m_final if m_initial > m_final
       ! m_initial + photon -> m_final if m_initial < m_final

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       include 'common.f'

       double precision px,py,pz,p0,m_initial,m_final
       double precision vx,vy,vz
       double precision utheta,phi,p_companion
       double precision rlu
       double precision m_pi

       m_pi=0.13d0       

       vx=px/p0
       vy=py/p0
       vz=pz/p0

       call rotbos(0d0,0d0,-vx,-vy,-vz,px,py,pz,p0)

       m_initial=p0

       if(dabs(m_final-m_initial).gt.m_pi) then
          p_companion = sqrt(((m_initial**2+m_pi**2-m_final**2)
     &                  /2d0/m_initial)**2-m_pi**2)
       else
          p_companion = dabs(m_initial**2-m_final**2)/2d0/m_initial
       endif       

       utheta=rlu(0)*2d0-1d0
       phi=rlu(0)*2d0*pi
       px=-p_companion*sqrt(1-utheta*utheta)*cos(phi)
       py=-p_companion*sqrt(1-utheta*utheta)*sin(phi)
       pz=-p_companion*utheta      
       p0=sqrt(px*px+py*py+pz*pz+m_final*m_final)

       call rotbos(0d0,0d0,vx,vy,vz,px,py,pz,p0)

       return
       end

