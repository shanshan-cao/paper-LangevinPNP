c**********************************************************************c
c                                                                      c
      subroutine writeeventoscar
c                                                                      c
c    writes the processed event in OSCAR 1997A format                  c
c                                                                      c
c    last change: 08/06/2009 David Meidinger                           c
c                                                                      c
c**********************************************************************c

      include 'common.f'

      integer j

c**** write event-header **********************************************c

c      write (unit=OUTU,fmt=911,err=321) 
c     .     eh_event, oeh_npart, eh_imppar, eh_angle, eh_ntim, eh_itim,
c     .     eh_evsamp

      write (unit=OUTU,fmt=911,err=321) 
     .     eh_event, oeh_npart, eh_imppar, eh_angle, eh_ntim, eh_itim,
     .     NUMSAMP

c**** write particles *************************************************c

      do 211 j=1,oeh_npart
c         write (unit=OUTU,fmt=921,err=322)
c     .        op_num(j), op_id(j), op_px(j), op_py(j), op_pz(j), 
c     .        op_p0(j), op_mass(j), op_rx(j), op_ry(j), op_rz(j),
c     .        op_r0(j), oThydro(j), oc_vx(j), oc_vy(j), oc_vz(j),
c     .        op_fx(j), op_fy(j), op_fz(j), op_f0(j)

         write (unit=OUTU,fmt=921,err=322)
     .        op_num(j), op_id(j), op_px(j), op_py(j), op_pz(j), 
     .        op_p0(j), op_mass(j), op_rx(j), op_ry(j), op_rz(j),
     .        op_r0(j), oThydro(j), oc_vx(j), oc_vy(j), oc_vz(j),
     .        op_ipx(j), op_ipy(j), op_ipz(j), op_wt(j)
 211  continue

      return

c**** errors **********************************************************c

 321  write (ERRU,*) 'ERROR while writing event-header'
      stop    

 322  write (ERRU,*) 'ERROR while writing particle data'
      stop

      end

