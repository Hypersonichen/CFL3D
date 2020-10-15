      subroutine plot3d_xu(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q, 
     .                  qi0,qj0,qk0,x,y,z,xw,blank2,blank,iflag,
     .                  vist3d,iover,nblk,nmap,bcj,bck,bci,sj,sk,si,
     .                  vj0,vk0,vi0)
c***********************************************************************
c     Purpose:  Write the output file at the grid points in PLOT3D
c     format and print solution data.
c
c     outputs grid/solution in single precision for use with FAST/PLOT3D
c***********************************************************************
c
      include 'cfl1.h'
      include 'cfl2.h'
      include 'cfl3.h'
      parameter(nn=maxbl)
c
      real*4    xw(jdim,kdim,idim,5),xmachw,alphww,reuew,timew
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension vist3d(jdim,kdim,idim)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4)
c
      dimension blank2(jdim,kdim,idim,2),blank(jdim,kdim,idim)
      dimension nmap(maxbl)
c
      character*4 charmov
c
      common /bin/ ibin
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /grbl/ rkap0g(nn,3),nblock,levelg(nn),igridg(nn),
     .        iflimg(nn,3),ifdsg(nn,3),iviscg(nn,3),
     .        jdimg(nn),kdimg(nn),idimg(nn),idiagg(nn,3),
     .        nblcg(nn),idegg(nn,3),
     .        jsg(nn),ksg(nn),isg(nn),jeg(nn),keg(nn),ieg(nn),mit(5,nn)
     .        ,kte1g(nn),kte2g(nn),ilamlog(nn),ilamhig(nn),jlamlog(nn),
     .        jlamhig(nn),klamlog(nn),klamhig(nn),iwfg(nn,3)
      common /motion/ utrans(nn),vtrans(nn),wtrans(nn),omegax(nn),
     .        omegay(nn),omegaz(nn),xorig(nn),yorig(nn),zorig(nn),
     .        dxmx(nn),dymx(nn),dzmx(nn),dthxmx(nn),dthymx(nn),
     .        dthzmx(nn),thetax(nn),thetay(nn),thetaz(nn),rfreqt(nn),
     .        rfreqr(nn),xorig0(nn),yorig0(nn),zorig0(nn),time2(nn),
     .        thetaxl(nn),thetayl(nn),thetazl(nn),itrans(nn),irotat(nn),
     .        xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,omegaxmc,
     .        omegaymc,omegazmc,xorigmc,yorigmc,zorigmc,xorig0mc,
     .        yorig0mc,zorig0mc,thetaxmc,thetaymc,thetazmc,dxmxmc,
     .        dymxmc,dzmxmc,dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .        rfreqrmc,time2mc,itransmc,irotatmc,idef
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),istrag,impl,iru,
     .        nju,ijac,iaf,nres,levelb(5),mgflag,iconsf,mseq,
     .        ncyc1(5),levelt(5),nitfo1(5),mmx,ngam,nsm(5),iipv
      common /patch/ windex(maxxe,2),ninter,iindex(intmax,6*nsub1+9),
     .        nblkpt(maxxe)
      common /blkcom/ nblkk(2,mxbli),nbli,limblk(2,6,mxbli),
     .        isva(2,2,mxbli),nblon(mxbli)
      common /pl3d/ nplot3d,inpl3d(nplots,11)
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /wallfun/ iwf(3)
      common /igrdtyp/ ip3dgrd,ialph
      common /moov/movie,nframes,icall1,lhdr
      common /conversion/ radtodeg
c
	common /outfile/ nfile,movnum
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c*******************************************************************************
c
c     设置输出的文件编号
c
	nfile = nfile+1
c	write(*,*)'nfile=',nfile
	write(charmov(1:4),'(i4.4)')nfile
	open(unit=234,file='./dat/movie_'//charmov//'.bin',
     .     form='unformatted',status='unknown')
c*******************************************************************************
c
c
c
c*******************************************************************************
c     set iblank (blank2) array
c*******************************************************************************
c
c     assign default iblank (blank2) array 
c
      do 2049 i=1,idim
      do 2049 j=1,jdim
      do 2049 k=1,kdim
      blank2(j,k,i,1) = 1.
      blank2(j,k,i,2) = 1.
 2049 continue
c
c     zero out edges and corners if desired (izero > 0)
c
      izero = 0
c     
      if (izero.gt.0) then
         do 2050 i=1,idim
         do 2050 j=1,jdim,jdim1
         do 2050 k=1,kdim,kdim1
         blank2(j,k,i,1) = 0.
 2050    continue
c
         if (idim.gt.2) then
            do 2051 i=1,idim,idim1
            do 2151 j=1,jdim
            blank2(j,1,i,1)    = 0.
            blank2(j,kdim,i,1) = 0.
 2151       continue
            do 2251 k=1,kdim
            blank2(1,k,i,1)    = 0.
            blank2(jdim,k,i,1) = 0.
 2251       continue
 2051       continue
         else
            do 2052 j=1,jdim
            do 2052 k=1,kdim
            blank2(j,k,idim,1) = 0.
 2052       continue
         end if
c
         if (jdim.eq.2) then
            do 2054 i=1,idim
            do 2054 k=1,kdim
            blank2(jdim,k,i,1) = 0.
 2054       continue
         end if
c
         if (kdim.eq.2) then
            do 2056 i=1,idim
            do 2056 k=1,kdim
            blank2(j,kdim,i,1) = 0.
 2056       continue
         end if
      end if
c
c     solid surface iblank (blank2) values (iblank=2 for solid surface)
c
      j11 = 1
      j22 = jdim1
      if(jdim.eq.2) then
        j11 = 1
        j22 = 1
      end if
      k11 = 1
      k22 = kdim1
      if(kdim.eq.2) then
        k11 = 1
        k22 = 1
      end if
      i11 = 1
      i22 = idim1
      if(idim.eq.2) then
        i11 = 1
        i22 = 1
      end if
      i = 1
      do 1190 nnn=1,2
      do 1191 j=j11,j22
      do 1191 k=k11,k22
      kk = min(k+1,k22)
      jj = min(j+1,j22)
      blank2(j,k,i,1) = 1. + max( bci(j,k,nnn),  bci(jj,k,nnn),
     .                            bci(j,kk,nnn), bci(jj,kk,nnn) )
 1191 continue
      i = idim
 1190 continue
      j = 1
      do 2190 nnn=1,2
      do 2191 i=i11,i22
      do 2191 k=k11,k22
      kk = min(k+1,k22)
      ii = min(i+1,i22)
      blank2(j,k,i,1) = 1. + max( bcj(k,i,nnn),  bcj(kk,i,nnn),
     .                            bcj(k,ii,nnn), bcj(kk,ii,nnn) )
 2191 continue
      j = jdim
 2190 continue
      k = 1
      do 3190 nnn=1,2
      do 3191 j=j11,j22
      do 3191 i=i11,i22
      jj = min(j+1,j22)
      ii = min(i+1,i22)
      blank2(j,k,i,1) = 1. + max( bck(j,i,nnn),  bck(jj,i,nnn),
     .                            bck(j,ii,nnn), bck(jj,ii,nnn) )
 3191 continue
      k = kdim
 3190 continue
c
c     iblank (blank2) array for patch interface boundaries
c  
c     Set blank2 array on block edges from nblkpt array.  Nblkpt array
c     corresponds to cell center grid.  For ambiguous points which arise
c     when going from cell-center data to grid point data (e.g. along block
c     edges, corners, or for interior points on a patch interface for which 
c     surrounding cell-centers lie in different blocks), blank2 value 
c     remains equal to previously set values.
c
      if (abs(ninter).gt.0) then
         do 1600 inter=1,abs(ninter)
         lmax1  = iindex(inter,1)
         nbl    = iindex(inter,lmax1+2)
         if (nbl.ne.nblk) go to 1600
         lst    = iindex(inter,2*lmax1+5)
         lcoord = iindex(inter,2*lmax1+3)/10
         lend   = iindex(inter,2*lmax1+3)-lcoord*10
         j21    = iindex(inter,2*lmax1+6)
         j22    = iindex(inter,2*lmax1+7)
         k21    = iindex(inter,2*lmax1+8)
         k22    = iindex(inter,2*lmax1+9)
c
         if (lcoord.eq.1) then
            if (lend.eq.1) i = 1
            if (lend.eq.2) i = idim
            if (jdim.gt.2 .and. kdim.gt.2) then
               do 1610 j=j21+1,j22-1
               do 1610 k=k21+1,k22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(k-1-k21) + (j-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21-1) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1610          continue
            else if (jdim.gt.2) then
c              connecting block numbers of the two cell-centered points
c              surrounding the current grid point
               k = k21
               do 1611 j=j21+1,j22-1
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1611          continue
            else if (kdim.gt.2) then
               j = j21
               do 1612 k=k21+1,k22-1
               ll = lst + (j22-j21)*(k-k21-1) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1612          continue
            end if
         end if
c
         if (lcoord.eq.2) then
            if (lend.eq.1) j = 1
            if (lend.eq.2) j = jdim
            if (idim.gt.2 .and .kdim.gt.2) then
               do 1620 i=k21+1,k22-1
               do 1620 k=j21+1,j22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(i-1-k21) + (k-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21-1) + (k-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1620          continue
            else if (idim.gt.2) then
               k = j21
               do 1621 i=k21+1,k22-1
               ll = lst + (j22-j21)*(i-k21-1) + (k-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1621          continue
            else if (kdim.gt.2) then
               i = k21
               do 1622 k=j21+1,j22-1
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1622          continue
            end if
         end if
c
         if (lcoord.eq.3) then
            if (lend.eq.1) k = 1
            if (lend.eq.2) k = kdim
            if (idim.gt.2 .and. jdim.gt.2) then
               do 1630 i=k21+1,k22-1
               do 1630 j=j21+1,j22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(i-1-k21) + (j-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21-1) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1630          continue
            else if (idim.gt.2) then
               j = j21
               do 1631 i=k21+1,k22-1
               ll = lst + (j22-j21)*(i-k21-1) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1631          continue
            else if (jdim.gt.2) then
               i = k21
               do 1632 j=j21+1,j22-1
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1632          continue
            end if
         end if
c
 1600    continue
      end if
c
c     iblank (blank2) array for 1:1 interface boundaries
c
      if(nbli.gt.0) then
        do 100 n=1,abs(nbli)
        if(nblon(n).ge.0) then
          if(nblk.eq.nblkk(1,n) .or. nblk.eq.nblkk(2,n)) then
            it = 1
            ir = 2
            if(nblk.eq.nblkk(2,n)) then
              it = 2
              ir = 1
            end if
c
c           allow for 1-1 blocking in same grid
c
            itime = 1
            if (nblkk(1,n).eq.nblkk(2,n)) itime = 2
            do 101 iti = 1,itime
            if (iti.gt.1) then
               it = 1
               ir = 2
            end if
c
            is = limblk(it,1,n)
            ie = limblk(it,4,n)
            js = limblk(it,2,n)
            je = limblk(it,5,n)
            ks = limblk(it,3,n)
            ke = limblk(it,6,n)
c
c           cell center indicies ---> grid point indicies
c
            if(isva(it,1,n)+isva(it,2,n).eq.5) then
c             i = constant interface
              if(js.gt.je) js = js+1
              if(js.lt.je) je = je+1
              if(ks.gt.ke) ks = ks+1
              if(ks.lt.ke) ke = ke+1
c             2d cases
              if(jdim.eq.2) then
                js = 1
                je = 2
              end if
              if(kdim.eq.2) then
                ks = 1
                ke = 2
              end if
            end if
            if(isva(it,1,n)+isva(it,2,n).eq.4) then
c             j = constant interface
              if(is.gt.ie) is = is+1
              if(is.lt.ie) ie = ie+1
              if(ks.gt.ke) ks = ks+1
              if(ks.lt.ke) ke = ke+1
c             2d cases
              if(idim.eq.2) then
                is = 1
                ie = 2
              end if
              if(kdim.eq.2) then
                ks = 1
                ke = 2
              end if
            end if
            if(isva(it,1,n)+isva(it,2,n).eq.3) then
c             k = constant interface
              if(js.gt.je) js = js+1
              if(js.lt.je) je = je+1
              if(is.gt.ie) is = is+1
              if(is.lt.ie) ie = ie+1
c             2d cases
              if(jdim.eq.2) then
                js = 1
                je = 2
              end if
              if(idim.eq.2) then
                is = 1
                ie = 2
              end if
            end if
c
            is1 = min(is,ie)
            ie1 = max(is,ie)
            js1 = min(js,je)

            je1 = max(js,je)
            ks1 = min(ks,ke)
            ke1 = max(ks,ke)
c
            do 110 i = is1,ie1
            do 110 j = js1,je1
            do 110 k = ks1,ke1
            blank2(j,k,i,1) = -float(nmap(nblkk(ir,n)))
  110       continue
  101       continue
          end if
        end if
  100   continue
      end if
c
c     iblank (blank2) array for embedded grids - the underlying
c     coarse grid areas are blanked out if the parameter ibembed > 0
c
      ibembed = 1
c
      if (ibembed.gt.0) then
         do 7500 nblc=1,nblock
         if (nblk.eq.nblc) go to 7500
         nblcc    = nblcg(nblc)
         if (nblcc.eq.nblk) then
            js = jsg(nblc)
            if (js.lt.jdimg(nblcc) .and. js.gt.1) js = js + 1
            ks = ksg(nblc)
            if (ks.lt.kdimg(nblcc) .and. ks.gt.1) ks = ks + 1
            is = isg(nblc)
            if (is.lt.idimg(nblcc) .and. is.gt.1) is = is + 1
            je = jeg(nblc)
            if (je.gt.2 .and. je.lt.jdimg(nblcc)) je = je - 1
            ke = keg(nblc)
            if (ke.gt.2 .and. ke.lt.kdimg(nblcc)) ke = ke - 1
            ie = ieg(nblc)
            if (ie.gt.2 .and. ie.lt.idimg(nblcc)) ie = ie - 1
            do 7501 i=is,ie
            do 7501 j=js,je
            do 7501 k=ks,ke
            blank2(j,k,i,1) = 0.
 7501       continue
         end if
 7500    continue
      end if
c
c     iblank (blank2) array for overlapped grids
c
      if (iover.eq.1) then
c
c     interior of faces (i=1 and i=idim)
c
      do 413 j=2,jdim1
      do 413 k=2,kdim1
      blank2(j,k,1,2)    = min(blank(j,k,1),      blank(j-1,k,1),
     .                         blank(j,k-1,1),    blank(j-1,k-1,1))
      blank2(j,k,idim,2) = min(blank(j,k,idim1),  blank(j-1,k,idim1),
     .                         blank(j,k-1,idim1),blank(j-1,k-1,idim1))
  413 continue
c
c     edges and corners (i=1 and i=idim)
c
      do 414 m=1,2
      k  = 1
      kk = 1
      if (m.eq.2) then
         k  = kdim
         kk = kdim1
      end if
      do 414 j=2,jdim1
      blank2(j,k,1,2)    = min(blank(j,kk,1),    blank(j-1,kk,1))
      blank2(j,k,idim,2) = min(blank(j,kk,idim1),blank(j-1,kk,idim1))
  414 continue
      do 415 m=1,2
      j  = 1
      jj = 1
      if (m.eq.2) then
         j  = jdim
         jj = jdim1
      end if
      do 415 k=2,kdim1
      blank2(j,k,1,2)    = min(blank(jj,k,1),    blank(jj,k-1,1))
      blank2(j,k,idim,2) = min(blank(jj,k,idim1),blank(jj,k-1,idim1))
  415 continue
      blank2(jdim,kdim,1,2)    = blank2(jdim-1,kdim-1,1,2)
      blank2(1,1,1,2)          = blank2(2,2,1,2)
      blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim-1,idim,2)
      blank2(1,1,idim,2)       = blank2(2,2,idim,2)
c
c     interior of faces (j=1 and j=jdim)
c
      do 513 i=2,idim1
      do 513 k=2,kdim1
      blank2(1,k,i,2)    = min(blank(1,k,i),      blank(1,k,i-1),
     .                         blank(1,k-1,i),    blank(1,k-1,i-1))
      blank2(jdim,k,i,2) = min(blank(jdim1,k,i),  blank(jdim1,k,i-1),
     .                         blank(jdim1,k-1,i),blank(jdim1,k-1,i-1))
  513 continue
c
c     edges and corners (j=1 and j=jdim)
c  
      do 514 m=1,2
      k  = 1
      kk = 1
      if (m.eq.2) then
         k  = kdim
         kk = kdim1
      end if
      do 514 i=2,idim1
      blank2(1,k,i,2)    = min(blank(1,kk,i),    blank(1,kk,i-1))
      blank2(jdim,k,i,2) = min(blank(jdim1,kk,i),blank(jdim1,kk,i-1))
  514 continue
      do 515 m=1,2
      i  = 1
      ii = 1
      if (m.eq.2) then
         i  = idim
         ii = idim1
      end if
      do 515 k=2,kdim1
      blank2(1,k,i,2)    = min(blank(1,k,ii),    blank(1,k-1,ii))
      blank2(jdim,k,i,2) = min(blank(jdim1,k,ii),blank(jdim1,k-1,ii))
  515 continue
      blank2(1,kdim,idim,2)    = blank2(1,kdim-1,idim-1,2)
      blank2(1,1,1,2)          = blank2(1,2,2,2)
      blank2(jdim,kdim,idim,2) = blank2(jdim,kdim-1,idim-1,2)
      blank2(jdim,1,1,2)       = blank2(jdim,2,2,2)
c
c     interior of faces (k=1 and k=kdim)
c
      do 613 i=2,idim1
      do 613 j=2,jdim1
      blank2(j,1,i,2)    = min(blank(j,1,i),      blank(j,1,i-1),
     .                         blank(j-1,1,i),    blank(j-1,1,i-1))
      blank2(j,kdim,i,2) = min(blank(j,kdim1,i),  blank(j,kdim1,i-1),
     .                         blank(j-1,kdim1,i),blank(j-1,kdim1,i-1))
  613 continue
c
c     edges and corners (k=1 and k=kdim)
c
      do 614 m=1,2
      j  = 1
      jj = 1
      if (m.eq.2) then
         j  = jdim
         jj = jdim1
      end if
      do 614 i=2,idim1
      blank2(j,1,i,2)    = min(blank(jj,1,i),    blank(jj,1,i-1))
      blank2(j,kdim,i,2) = min(blank(jj,kdim1,i),blank(jj,kdim1,i-1))
  614 continue
      do 615 m=1,2
      i  = 1
      ii = 1
      if (m.eq.2) then
         i  = idim
         ii = idim1
      end if
      do 615 j=2,jdim1
      blank2(j,1,i,2)    = min(blank(j,1,ii),      blank(j-1,1,ii))
      blank2(j,kdim,i,2) = min(blank(j,kdim1,ii),  blank(j-1,kdim1,ii))
  615 continue
      blank2(jdim,1,idim,2)    = blank2(jdim-1,1,idim-1,2)
      blank2(1,1,1,2)          = blank2(2,1,2,2)
      blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim,idim-1,2)
      blank2(1,kdim,1,2)       = blank2(2,kdim,2,2)
c
c    interior cells
c
      do 713 i=2,idim1
      do 713 j=2,jdim1
      do 713 k=2,kdim1
      blank2(j,k,i,2) = min(blank(j,k,i),    blank(j-1,k,i),
     .                      blank(j,k-1,i),  blank(j-1,k-1,i),
     .                      blank(j,k,i-1),  blank(j-1,k,i-1),
     .                      blank(j,k-1,i-1),blank(j-1,k-1,i-1) )
  713 continue
c
c     combining topology and overset grid blankings
c
      do 813 i=1,idim
      do 813 j=1,jdim
      do 813 k=1,kdim
      blank2(j,k,i,1) = min(blank2(j,k,i,1),blank2(j,k,i,2))
  813 continue
c
      end if
c
      if (iflag.eq.1) then
c
c*******************************************************************************
c      plot3d data
c*******************************************************************************
c
c      if(lhdr.gt.0) then
c      if (i2d .eq. 1) then
c      write(11,'(''writing plot3d file for JDIM X KDIM ='',i5,
c     .'' x '',i5,'' grid'')') jdim,kdim
c      write(11,'(''   plot3d files to be read with 2d/mgrid/blank'',
c     .'' qualifiers'')')
c      else
c      write(11,93)idim,jdim,kdim
c   93 format(45h writing plot3d file for IDIM X JDIM X KDIM =,
c     .i5,3h x ,i5,3h x ,i5,5h grid)
c      write(11,2042)
c 2042 format(3x,52hplot3d files to be read with /mgrid/blank qualifiers)
c      end if
c      end if
c
c     load grid into single precision array
c
      do 9000 i=1,idim
      do 9000 j=1,jdim
      do 9000 k=1,kdim
      xw(j,k,i,1) = x(j,k,i)
      xw(j,k,i,2) = y(j,k,i)
      xw(j,k,i,3) = z(j,k,i)
 9000 continue
c
c     ialph > 0 for a grid that was read in plot3d format with alpha measured
c               in the xy plane (TLNS3D convention)
c
c     output grid
c
c      if(ibin.eq.0) then
c        if(i2d.eq.0) then
c        if(ialph.eq.0) then
c        if(icall1.eq.0)
c     .  write(3,*)(((xw(j,k,i,1),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((xw(j,k,i,2),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((xw(j,k,i,3),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((int(blank2(j,k,i,1)),i=i1,i2,i3),j=j1,j2,j3),
c     .               k=k1,k2,k3)
c        else
c        if(icall1.eq.0)
c     .  write(3,*)(((xw(j,k,i,1),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((xw(j,k,i,3),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((-xw(j,k,i,2),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((int(blank2(j,k,i,1)),i=i1,i2,i3),j=j1,j2,j3),
c     .              k=k1,k2,k3)
c        end if
c        else
c        if(icall1.eq.0)
c     .  write(3,*)(((xw(j,k,i,1),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((xw(j,k,i,3),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
c     .            (((int(blank2(j,k,i,1)),i=i1,i2,i3),j=j1,j2,j3),
c     .               k=k1,k2,k3)
c        end if
c      else
        if(i2d.eq.0) then
        if(ialph.eq.0) then
        if(icall1.eq.0)
     .  write(233)(((xw(j,k,i,1),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
     .          (((xw(j,k,i,2),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
     .          (((xw(j,k,i,3),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3)
        else
        if(icall1.eq.0)
     .  write(233)(((xw(j,k,i,1),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
     .          (((xw(j,k,i,3),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
     .          (((-xw(j,k,i,2),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3)
        end if
        else
        if(icall1.eq.0)
     .  write(233)(((xw(j,k,i,1),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3),
     .          (((xw(j,k,i,3),i=i1,i2,i3),j=j1,j2,j3),k=k1,k2,k3)
        end if
c      end if
c
c     determine q values at grid points; load into single precision array
c
      iw = 0
      do 2000 i=i1,i2,i3
      iw = iw+1
      kw = 0
      do 2001 k=k1,k2,k3
      kw = kw+1
      jw = 0
      do 2002 j=j1,j2,j3
      jw = jw+1
      do 3000 l=1,5
      id  = i
      id1 = id-1
      if (id1.le.0)   id1 = 1
      if (id.gt.idim1) id = idim1
      jd = j
      kd = k
c
      if (k.eq.1 .or. k.eq.kdim) then
c     k=1/k=kdim faces
         if (k.eq.kdim) kd = kdim1
         if (j.eq.1 .or. j.eq.jdim) then
c        edge points and corner points
            if (j.eq.jdim) jd = jdim1
            xw(jw,kw,iw,l) = 0.5*(q(jd,kd,id,l)+q(jd,kd,id1,l))
         else
c        interior points on  k=1/k=kdim faces
            m = 2
            if (k.eq.kdim) m = 4
            xw(jw,kw,iw,l) = .25*(qk0(jd,id,l,m)   +qk0(jd-1,id,l,m)
     .                          + qk0(jd,id1,l,m) + qk0(jd-1,id1,l,m))
         end if
      else if (j.eq.1 .or. j.eq.jdim) then
c     interior points on j=1/j=jdim faces
         if (j.eq.jdim) jd = jdim1
         m = 2
         if (j.eq.jdim) m = 4
         xw(jw,kw,iw,l) = .25*(qj0(kd,id,l,m)   +qj0(kd-1,id,l,m)
     .                       + qj0(kd,id1,l,m)  +qj0(kd-1,id1,l,m))
      else
c     interior points
         if (i.eq.1 .or. i.eq.idim) go to 3000
         xw(jw,kw,iw,l) = .125*(q(j,k,i,l)    +q(j,k,i-1,l)
     .                        + q(j-1,k,i,l)  +q(j-1,k,i-1,l)
     .                        + q(j,k-1,i,l)  +q(j,k-1,i-1,l)
     .                        + q(j-1,k-1,i,l)+q(j-1,k-1,i-1,l))
      end if
 3000 continue
 2002 continue
 2001 continue
 2000 continue
c
c     interior points on i=1/i=idim faces
      iw = 0
      do 2008 i=i1,i2,i3
      iw = iw+1
      if (i.ne.1 .and. i.ne.idim) go to 2008
      m  = 2
      if (i.eq.idim) m = 4
      kw = 0
      do 2007 k=k1,k2,k3
      kw = kw+1
      jw = 0
      do 2006 j=j1,j2,j3
      jw = jw+1
      if (j.eq.1 .or. j.eq.jdim .or. k.eq.1 .or. k.eq.kdim) go to 2006
c     interior points
      do 2005 l=1,5
      xw(jw,kw,iw,l) = .25*(qi0(j,k,l,m)   +qi0(j,k-1,l,m)
     .                    + qi0(j-1,k,l,m) +qi0(j-1,k-1,l,m))
 2005 continue
 2006 continue
 2007 continue
 2008 continue
c
c     output solution (q)
c
      do 4000 i=1,iw
      do 4000 k=1,kw
      do 4000 j=1,jw
      xw(j,k,i,5) = xw(j,k,i,5)/gm1+0.5*(xw(j,k,i,2)**2
     .             +xw(j,k,i,3)**2+xw(j,k,i,4)**2)*xw(j,k,i,1)
      xw(j,k,i,2) = xw(j,k,i,1)*xw(j,k,i,2)
      xw(j,k,i,3) = xw(j,k,i,1)*xw(j,k,i,3)
      xw(j,k,i,4) = xw(j,k,i,1)*xw(j,k,i,4)
 4000 continue
c
      alphaw=radtodeg*(alpha+thetay(nblk))
      xmachw = xmach
      alphww = alphaw
      reuew  = reue
      timew  = time
c
c      if (ibin.eq.0) then
c         write(4,*) xmachw,alphww,reuew,timew
c         if (i2d.eq.0) then
c            if (ialph.eq.0) then
c               write(4,*)((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
c            else
c               write(4,*)((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
c     .                    (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
c     .                    (((-xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
c     .                    (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw) 
c            end if
c         else
c            write(4,*)((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
c     .                ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
c         end if
c      else
c         write(4)  xmachw,alphww,reuew,timew
c*************************************************************************
c
c     输出流场信息
c
c*************************************************************************
         if (i2d.eq.0) then
            if (ialph.eq.0) then
               write(234)((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
            else
               write(234)((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .                  (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .                  (((-xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .                  (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw) 
            end if
         else
            write(234)  ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .                ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
         end if
c      end if
      end if
c
      close(234)
c
      return
      end
