c
c     Sampling flowfield ...
c
      subroutine sampling(jdim,kdim,idim,x,y,z,q,myid,ntime_accum)
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5)
      character*6 tchar2
      character*16 charmov
c
      common /twod/ i2d
c
cccc
c      write(charmov2,'(i0)')nt
c      open(1025, file='flow_'//trim(charmov2)//'.dat')
cccc
#ifdef SAMPLING_BIN

      !write(charmov(1:8),'(i8.8)')ntime_accum
      write(charmov,'(i0)')ntime_accum
c#ifndef DIST_MPI
      !write(tchar2(1:3),'(i3.3)')myid
      write(tchar2,'(i0)')myid
      ifile = 224+myid
      !open(ifile,file='sample_'//charmov//'_'//tchar2//'.bin',
      open(ifile,file=
     .'sample_'//trim(charmov)//'_'//trim(tchar2)//'.bin',
     .form='unformatted',status='unknown')
c#endif
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
      
      if(i2d.eq.1)then
      i = 1
      do 1001 k = 1,kdim
      do 1001 j = 1,jdim
      write(ifile)x(j,k,i),z(j,k,i),(q(j,k,i,ivar),ivar=1,2),
     .(q(j,k,i,ivar),ivar=4,5)
 1001 continue
c
      else
c
      do 1002 k = 1,kdim
      do 1002 j = 1,jdim
      do 1002 i = 1,idim
      write(ifile)x(j,k,i),y(j,k,i),z(j,k,i),(q(j,k,i,ivar),ivar=1,5)
 1002 continue
c
      endif
      close(ifile)


#else


      !write(charmov(1:8),'(i8.8)')ntime_accum
      write(charmov,'(i0)')ntime_accum
c#ifndef DIST_MPI
      !write(tchar2(1:3),'(i3.3)')myid
      write(tchar2,'(i0)')myid
      ifile = 224+myid
      !open(ifile,file='sample_'//charmov//'_'//tchar2//'.dat',
      open(ifile,file=
     .'sample_'//trim(charmov)//'_'//trim(tchar2)//'.dat',
     .status='unknown') 
c#endif
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
      
      if(i2d.eq.1)then
      i = 1
      write(ifile,*)'variables= "x" "z" "rho" "u" "w" "p"'
      write(ifile,*)'zone f= point,j=',jdim,'k=',kdim
      do 1001 k = 1,kdim
      do 1001 j = 1,jdim
      write(ifile,*)x(j,k,i),z(j,k,i),(q(j,k,i,ivar),ivar=1,2),
     .(q(j,k,i,ivar),ivar=4,5)
 1001 continue
c
      else
c
      write(ifile,*)'variables= "x" "y" "z" "rho" "u" "v" "w" "p"'
      write(ifile,*)'zone f= point,i=',idim,'j=',jdim,'k=',kdim
      do 1002 k = 1,kdim
      do 1002 j = 1,jdim
      do 1002 i = 1,idim
      write(ifile,*)x(j,k,i),y(j,k,i),z(j,k,i),(q(j,k,i,ivar),ivar=1,5)
 1002 continue
c
      endif
      close(ifile)


#endif




      return
      end
