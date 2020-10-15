c 
      subroutine qout_xu(iseq,lw,w,mgwk,wk,nwork)                          
c***********************************************************************
c     Purpose:  Output data for plotting or printing.
c***********************************************************************
c
      include 'cfl1.h'
      include 'cfl3.h'
      include 'cfl4.h'
      parameter(nn=maxbl)
c
c     maxbl   - maximum number of blocks
c     maxgr   - maximum number of grids
c     nplots  - maximum number of data sets to output via PLOT3D or print
c     ncycmax - maximum number of time-steps/cycles
c
      dimension w(mgwk)
      dimension lw(47,nn)
      dimension wk(nwork)
      dimension nmap(maxbl)
c
      integer*4 ip3ddim(3,nplots),np3d
c
      common /bin/ ibin
      common /chmra/ iovrlp(nn)
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv
      common /grbl/ rkap0g(nn,3),nblock,levelg(nn),igridg(nn),
     .        iflimg(nn,3),ifdsg(nn,3),iviscg(nn,3),
     .        jdimg(nn),kdimg(nn),idimg(nn),idiagg(nn,3),
     .        nblcg(nn),idegg(nn,3),
     .        jsg(nn),ksg(nn),isg(nn),jeg(nn),keg(nn),ieg(nn),mit(5,nn)
     .        ,kte1g(nn),kte2g(nn),ilamlog(nn),ilamhig(nn),jlamlog(nn),
     .        jlamhig(nn),klamlog(nn),klamhig(nn),iwfg(nn,3)
      common /grids/ ngrid,ncgg(maxgr),nblg(maxgr),iemg(maxgr),
     .        inewgg(maxgr)
      common /hstry/ rms(ncycmax),rmsb,rmstb1,rmstb2,clw(ncycmax),
     .        cdw(ncycmax),cdpw(ncycmax),cdvw(ncycmax),
     .        cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .        cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .        chdw(ncycmax),swetw(ncycmax),
     .        fmdotw(ncycmax),cfttotw(ncycmax),
     .        cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .        rmstr1(ncycmax),rmstr2(ncycmax),
     .        nneg1(ncycmax),nneg2(ncycmax),ntr,ihstry
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),istrag,impl,iru,
     .        nju,ijac,iaf,nres,levelb(5),mgflag,iconsf,mseq,
     .        ncyc1(5),levelt(5),nitfo1(5),mmx,ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /pl3d/ nplot3d,inpl3d(nplots,11)
      common /prout/ nprint,inpr(nplots,11)
      common /reyue/ reue,tinf,ivisc(3)
      common /skips/ iadvance(nn),iforce(nn),lfgm
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /moov/movie,nframes,icall1,lhdr
      common /conversion/ radtodeg
c
      ibin = 1
c
c*******************************************************************************
c     write output file in PLOT3D format
c     iptype = 0....output grid and q files at grid points
c            = 1....output grid and q files at cell centers
c            = 2....output grid and turb files at cell centers
c            > 2....output grid and function file at cell centers
c*******************************************************************************
c
      if (nplot3d.le.0) go to 231
      ncount    = 0
      np3d      = nplot3d
      ifunc     = 0
      do 60 n=1,nplot3d
      if (n.eq.1) then
      if (ibin.eq.0) then
         if (icall1.eq.0) write(3,*) np3d
         write(4,*) np3d
      else
         if (icall1.eq.0) write(3) np3d
         write(4) np3d
      end if
      end if
c
      nbl = inpl3d(n,1)
c
      if (nbl.gt.nblock) then
      write(11,777)nbl
  777 format(6h Block,i3,43h does not exist.  No plot3d output printed.)
      go to 60
      end if
c
      ncount = ncount+1
      i1 = inpl3d(n,3)
      i2 = inpl3d(n,4)
      i3 = inpl3d(n,5)
      j1 = inpl3d(n,6)
      j2 = inpl3d(n,7)
      j3 = inpl3d(n,8)
      k1 = inpl3d(n,9)
      k2 = inpl3d(n,10)
      k3 = inpl3d(n,11)
c
c     if (inpl3d(n,2).eq.1) then
      if (inpl3d(n,2).gt.0) then
c        if zone has function file output, all do, and all have
c        the same function output
         if (inpl3d(n,2).gt.2) ifunc = inpl3d(n,2)
         call lead(nbl,lw)
         i2 = min(idim-1,i2)
         j2 = min(jdim-1,j2)
         k2 = min(kdim-1,k2)
         i1 = min(idim-1,i1)
         j1 = min(jdim-1,j1)
         k1 = min(kdim-1,k1)
      end if
      ip3ddim(1,ncount) = (i2-i1)/i3+1
      ip3ddim(2,ncount) = (j2-j1)/j3+1
      ip3ddim(3,ncount) = (k2-k1)/k3+1
   60 continue
c
      if (ibin.eq.0) then
         if (i2d.eq.0) then
c            if (icall1.eq.0) 
c     .      write(3,*) ((ip3ddim(i,n),i=1,3),n=1,ncount)
            if (ifunc.eq.0) then
c               write(4,*) ((ip3ddim(i,n),i=1,3),n=1,ncount)
            else
c               write(4,*) ((ip3ddim(i,n),i=1,3),1,n=1,ncount)
            end if
         else
c            if (icall1.eq.0) 
c     .      write(3,*) ((ip3ddim(i,n),i=2,3),n=1,ncount)
            if (ifunc.eq.0) then
c               write(4,*) ((ip3ddim(i,n),i=2,3),n=1,ncount)
            else
c               write(4,*) ((ip3ddim(i,n),i=2,3),1,n=1,ncount)
            end if
         end if
      else
         if (i2d.eq.0) then
c            if (icall1.eq.0) 
c     .      write(3) ((ip3ddim(i,n),i=1,3),n=1,ncount)
            if (ifunc.eq.0) then
c               write(4) ((ip3ddim(i,n),i=1,3),n=1,ncount)
            else
c               write(4) ((ip3ddim(i,n),i=1,3),1,n=1,ncount)
            end if
         else
c            if (icall1.eq.0) 
c     .      write(3) ((ip3ddim(i,n),i=2,3),n=1,ncount)
            if (ifunc.eq.0) then
c               write(4) ((ip3ddim(i,n),i=2,3),n=1,ncount)
            else
c               write(4) ((ip3ddim(i,n),i=2,3),1,n=1,ncount)
            end if
         end if
      end if
c
c     correspondence between global block number and plot3d block
c     number stored in nmap(n) for n=1,nblock
c
      do 1097 n = 1,nblock
      nmap(n) = 1
      do 1098 nnn = 1,nplot3d
      m = inpl3d(nnn,1)
      if (n.eq.m) nmap(n) = nnn
 1098 continue
 1097 continue             
c
      if (lhdr .gt. 0) write(11,1096)
 1096 format(1h )
c
      do 70 n=1,nplot3d
c
      nbl = inpl3d(n,1)
c
      if (nbl.gt.nblock .or. nbl.le. 0) go to 70
c
      i1 = inpl3d(n,3)
      i2 = inpl3d(n,4)
      i3 = inpl3d(n,5)
      j1 = inpl3d(n,6)
      j2 = inpl3d(n,7)
      j3 = inpl3d(n,8)
      k1 = inpl3d(n,9)
      k2 = inpl3d(n,10)
      k3 = inpl3d(n,11)
c
      call lead(nbl,lw)
      iflag = 1
c
      if (inpl3d(n,2).eq.0) then
c
c     check for need to call bc and qface
c
         if (n.ne.1) then
            if (inpl3d(n,1).eq.inpl3d(n-1,1) .and. 
     .          inpl3d(n-1,2).eq.0) go to 63
         end if
         if (iadvance(nbl).ge.0) then
c   Need to temporarily set level .ne. lglobal so won't do turb bcs in bc
            leveltemp=level
            level=lglobal-1
            nttuse=max(ntt-1,1)
            call bc(1,nbl,lw,w,mgwk,wk,nwork,clw(nttuse))
            level=leveltemp
         end if
c
         call qface(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .              w(lbcj),w(lbck),w(lbci))
   63    ibwk = jdim*kdim*idim*5+1
         call plot3d_xu(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
     .               w(lq),w(lqi0),w(lqj0),w(lqk0),w(lx),w(ly),w(lz),
     .               wk,wk(ibwk),w(lblk),iflag,w(lvis),iovrlp(nbl),
     .               nbl,nmap,w(lbcj),w(lbck),w(lbci),
     .               w(lsj),w(lsk),w(lsi),w(lvj0),w(lvk0),w(lvi0))
         call resetwk(jdim,kdim,idim,2,wk(ibwk))
      else if (inpl3d(n,2).eq.1 .or. inpl3d(n,2).gt.2) then
         i2 = min(idim-1,i2)
         j2 = min(jdim-1,j2)
         k2 = min(kdim-1,k2)
         i1 = min(idim-1,i1)
         j1 = min(jdim-1,j1)
         k1 = min(kdim-1,k1)
         ibwk = jdim*kdim*idim*5+1
c         call plot3c(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
c     .               w(lq),w(lx),w(ly),w(lz),wk,wk(ibwk),w(lblk),iflag,
c     .               w(lvis),iovrlp(nbl),nbl,nmap,w(lsj),w(lsk),w(lsi),
c     .               w(lsnk0),ifunc,qav)
      else
c   Write out additional turb info to plot3d-type files
         if (ivmx .gt. 1) then
         i2 = min(idim-1,i2)
         j2 = min(jdim-1,j2)
         k2 = min(kdim-1,k2)
         i1 = min(idim-1,i1)
         j1 = min(jdim-1,j1)
         k1 = min(kdim-1,k1)
c  if ivmx .eq. 8,9, or .ge. 11, already have storage for ux:
         if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.ge.11) then
         ibwk = jdim*kdim*idim*5+1
         ibwk2 = ibwk +jdim*kdim*idim
         nroom=nwork-(ibwk2+jdim*kdim*9)
         if(nroom .lt. 0.) then
          write(11,'('' not enough memory for plot3t in cbsem.'')')
          write(11,'('' nroom='',i12)') nroom
          write(11,'('' not writing out turb info'')')
         else
          call plot3t(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
     .               w(lq),w(lx),w(ly),w(lz),wk,wk(ibwk),w(lblk),
     .               w(lvis),iovrlp(nbl),nbl,nmap,w(lsj),w(lsk),w(lsi),
     .               w(lsnk0),w(lux),w(lxib),w(lvol),w(lqj0),w(lqk0),
     .               w(lqi0),w(lbcj),w(lbck),w(lbci),wk(ibwk2),w(lcmuv))
         end if
         else
         ibwk = jdim*kdim*idim*5+1
         ibwk2 = ibwk +jdim*kdim*idim
         ibwk3 = ibwk2+(jdim-1)*(kdim-1)*(idim-1)*9
         nroom=nwork-(ibwk3+jdim*kdim*9)
         if(nroom .lt. 0.) then
          write(11,'('' not enough memory for plot3t in cbsem.'')')
          write(11,'('' nroom='',i12)') nroom
          write(11,'('' not writing out turb info'')')
         else
          call plot3t(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
     .               w(lq),w(lx),w(ly),w(lz),wk,wk(ibwk),w(lblk),
     .               w(lvis),iovrlp(nbl),nbl,nmap,w(lsj),w(lsk),w(lsi),
     .               w(lsnk0),wk(ibwk2),w(lxib),w(lvol),w(lqj0),w(lqk0),
     .               w(lqi0),w(lbcj),w(lbck),w(lbci),wk(ibwk3),w(lcmuv))
         end if
         end if
         end if
         call resetwk(jdim,kdim,idim,1,wk(ibwk))
      end if
   70 continue
c
      if (lhdr .gt. 0) write(11,1096)
c
c     for stationary grid cases, set icall1 flag to prevent output
c     of grid to plot3d file on subsequent calls (which occur only if 
c     abs(movie) > 0).  For dynamic grid cases, grid is output every
c     time plot3d routine is is called.
c
      if (iunst.eq.0) then
         icall1 = 1
      else
         icall1 = 0
      end if
c
  231 continue
c
      return
      end
