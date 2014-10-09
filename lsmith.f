      program lsmith
      implicit none
c------------------------------------------------------------------------------
c
c This function returns dsigma^2/dq^2dE for QE scattering from Llewellyn-Smith
c 
c------------------------------------------------------------------------------
      integer ie,iq,i
      integer nwpawc
      parameter (nwpawc = 100000)
      common/pawc/h(nwpawc)
      real dsigma,dsigmaL,dsigmaM,dsigmaH,enu,qsq,h,ma,corr
      real binwid,emin
      real rlsdd
      external rlsdd
      logical first,nubar
      data first /.true./ 
      logical nuclcorr,morefluxbins     
      integer nq
      parameter (nq=500)
      real qlim
      parameter (qlim=5.0)
c..BNL flux parameters
      integer ne
      real eflux(15)
      data eflux/0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 
     +         6.0, 8.0, 10.0, 15.0, 20.0/

c
c..initialization
      dsigma=0.
      nubar=.false.
      nuclcorr=.false.      ! apply nuclear correction for D2?
      morefluxbins=.true.  ! calculate for a larger number of E bins?
      ma=1.10               ! value for axial mass (in GeV)
c
      if (morefluxbins) then
        ne=50
        binwid=0.1
        emin=0.1
      else
        ne=15
      endif
c
      if (first) then
         first=.false.
         print *,''
         print *,'----------------------------------------------------'
         print *,'Llewellyn-Smith free nucleon QE cross section calc'
         print *,'nubar=',nubar
         print *,'nuclear correction=',nuclcorr
         print *,'more BNL flux bins?',morefluxbins
         print *,'assuming DIPOLE form for F_A with mA=',ma
         print *,'----------------------------------------------------'
         print *,''
	 call hlimit(nwpawc)
         do i=1,ne
c           call hbprof(100+i,'dsig/dQ2',nq,0.,qlim,0.,10000.,' ')
           call hbook1(100+i,'dsig/dQ2,',nq,0.,qlim,0.)
         enddo
      endif
c
      enu=emin-binwid
      do ie = 1,ne
       if (morefluxbins) then
         enu=enu+binwid
       else
         enu=eflux(ie)
       endif
       do iq=1,nq
         qsq=qlim*(iq-0.5)/nq
           dsigma=rlsdd(enu,qsq,nubar,ma)*1.E15*1.E32
c..apply nuclear correction for D2 if selected, this is a function of Q2 only
           if (nuclcorr) then
            if (qsq.le.0.14) then
              corr=0.45246+27.873*qsq-841.93*qsq**2+14938*qsq**3
     +             -0.15328E6*qsq**4+0.89676E6*qsq**5-0.27769E7*qsq**6
     +             +0.35278E7*qsq**7 
            else
              corr=1.0
            endif
            dsigma=dsigma*corr
           endif
c           call hfill(100+ie,qsq,dsigma,1.)
           call hf1(100+ie,qsq,dsigma)
       enddo
      enddo
      call hrput(0,'lsmith.hbook','N')
c
      return
      end
