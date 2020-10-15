      module module_dissipation 
#ifdef HYBRID
      real,allocatable:: tfai(:,:,:)
      type dissArr
c---------------------------L. Chen: legacy code... need to fix this--
      real,allocatable:: fai(:,:,:)
c      real,allocatable:: dqqx(:,:,:,:),
c     .                   dqqy(:,:,:,:),
c     .                   dqqz(:,:,:,:)
      real,allocatable:: dqx(:,:,:,:),
     .                   dqy(:,:,:,:),
     .                   dqz(:,:,:,:)
c      real,allocatable:: csj(:,:,:,:),csk(:,:,:,:),csi(:,:,:,:)
c      real,allocatable:: cr(:,:,:),ci(:,:,:),ct(:,:,:)
c---------------------------------------------------------------------
      end type

      type:: pDissArr
      type(dissArr),pointer:: pv
      end type
#endif
#ifdef HYBRID
      type(pDissArr),allocatable:: dissArr_blks(:)
#endif

      contains 
  
#ifdef HYBRID 
c      subroutine diss_alloc(maxbl)
c      
c      end subroutine

      subroutine artif4(jdim,kdim,idim,qq,qqj0,qqk0,qqi0,bci,
     .                  bcj,bck,sj,sk,si,dqx,dqy,dqz)
c***********************************************************************
c	L.Chern: For Hybrid upwind/central scheme, 4th a.v. only
c*********************************************************************** 
      dimension qq(jdim,kdim,idim,5),q(jdim,kdim,idim,5),
     .          qqj0(kdim,idim-1,5,4),qqk0(jdim,idim-1,5,4),
     .          qqi0(jdim,kdim,5,4),qj0(kdim,idim-1,5,4),
     .          qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension dqx(jdim,kdim,idim-1,5),dqy(jdim,kdim,idim-1,5),
     .          dqz(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension dz(-1:jdim,-1:kdim,-1:idim,5),
     .          d3z(0:jdim-1,0:kdim-1,0:idim-1,5)
      dimension emuj4(0:jdim-1,0:kdim-1,0:idim-1),
     .          psi(0:jdim,0:kdim,0:idim)
      dimension r1(0:jdim,0:kdim,0:idim),
     .          r2(0:jdim,0:kdim,0:idim),
     .          r3(0:jdim,0:kdim,0:idim)
      dimension r(jdim-1,kdim-1,idim-1)
c
	common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
	common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
	common /twod/ i2d
c
	eps2  = 1./4.
	eps4  = 1./128.  !1./256.
	eta   = 0.1      !  eta取值范围为: 0 ~ 1
	idim1 = idim-1
	jdim1 = jdim-1
	kdim1 = kdim-1
	idim2 = idim-2
	jdim2 = jdim-2
	kdim2 = kdim-2
c
c     原始量转换成守恒量
c
      do 5001 n = 1,4
      do 5001 k = 1,kdim
	do 5001 i = 1,idim1
	qj0(k,i,1,n) = qqj0(k,i,1,n)
	qj0(k,i,2,n) = qqj0(k,i,1,n)*qqj0(k,i,2,n)
	qj0(k,i,3,n) = qqj0(k,i,1,n)*qqj0(k,i,3,n) 
	qj0(k,i,4,n) = qqj0(k,i,1,n)*qqj0(k,i,4,n)
	qj0(k,i,5,n) = qqj0(k,i,5,n)/gm1+0.5*qqj0(k,i,1,n)*
     .			(qqj0(k,i,2,n)**2+qqj0(k,i,3,n)**2+qqj0(k,i,4,n)**2)  
 5001 continue
c
      do 5002 n = 1,4
      do 5002 j = 1,jdim
	do 5002 i = 1,idim1
	qk0(j,i,1,n) = qqk0(j,i,1,n)
	qk0(j,i,2,n) = qqk0(j,i,1,n)*qqk0(j,i,2,n)
	qk0(j,i,3,n) = qqk0(j,i,1,n)*qqk0(j,i,3,n) 
	qk0(j,i,4,n) = qqk0(j,i,1,n)*qqk0(j,i,4,n)
	qk0(j,i,5,n) = qqk0(j,i,5,n)/gm1+0.5*qqk0(j,i,1,n)*
     .			(qqk0(j,i,2,n)**2+qqk0(j,i,3,n)**2+qqk0(j,i,4,n)**2)  
 5002 continue
c
      do 5003 n = 1,4
      do 5003 j = 1,jdim
	do 5003 k = 1,kdim
	qi0(j,k,1,n) = qqi0(j,k,1,n)
	qi0(j,k,2,n) = qqi0(j,k,1,n)*qqi0(j,k,2,n)
	qi0(j,k,3,n) = qqi0(j,k,1,n)*qqi0(j,k,3,n) 
	qi0(j,k,4,n) = qqi0(j,k,1,n)*qqi0(j,k,4,n)
	qi0(j,k,5,n) = qqi0(j,k,5,n)/gm1+0.5*qqi0(j,k,1,n)*
     .			(qqi0(j,k,2,n)**2+qqi0(j,k,3,n)**2+qqi0(j,k,4,n)**2)  
 5003 continue
c
	do 5004 j = 1,jdim
	do 5004 k = 1,kdim
      do 5004 i = 1,idim
	q(j,k,i,1) = qq(j,k,i,1)
	q(j,k,i,2) = qq(j,k,i,1)*qq(j,k,i,2)
	q(j,k,i,3) = qq(j,k,i,1)*qq(j,k,i,3)
	q(j,k,i,4) = qq(j,k,i,1)*qq(j,k,i,4)
	q(j,k,i,5) = qq(j,k,i,5)/gm1+0.5*qq(j,k,i,1)*
     .           (qq(j,k,i,2)**2+qq(j,k,i,3)**2+qq(j,k,i,4)**2)
 5004 continue
c
c     计算谱半径
c
	do 100 j = 1,jdim1
	do 100 k = 1,kdim1
	do 100 i = 1,idim1
	c2 = gamma*qq(j,k,i,5)/qq(j,k,i,1)
	c  = sqrt(c2)
c
c     J  方向的
c
	t1 = 0.5*(sj(j,k,i,1)+sj(j+1,k,i,1))
	t2 = 0.5*(sj(j,k,i,2)+sj(j+1,k,i,2))
	t3 = 0.5*(sj(j,k,i,3)+sj(j+1,k,i,3))
	t4 = 0.5*(sj(j,k,i,4)+sj(j+1,k,i,4))
	ubar = qq(j,k,i,2)*t1+qq(j,k,i,3)*t2+qq(j,k,i,4)*t3
	ubar = ubar*t4
	r2(j,k,i) = abs(ubar)+c*t4
c
c     K 方向的
c
	t1 = 0.5*(sk(j,k,i,1)+sk(j,k+1,i,1))
	t2 = 0.5*(sk(j,k,i,2)+sk(j,k+1,i,2))
	t3 = 0.5*(sk(j,k,i,3)+sk(j,k+1,i,3))
	t4 = 0.5*(sk(j,k,i,4)+sk(j,k+1,i,4))
	ubar = qq(j,k,i,2)*t1+qq(j,k,i,3)*t2+qq(j,k,i,4)*t3
	ubar = ubar*t4
	r3(j,k,i) = abs(ubar)+c*t4 
c
c     测试程序
c
c      if(r2(j,k,i).gt.r2max)r2max=r2(j,k,i)
c	if(r2(j,k,i).lt.r2min)r2min=r2(j,k,i)
c	if(r2(j,k,i).gt.eee)icount=icount+1
  100	continue
c	write(*,*)'r2max=',r2max,'r2min=',r2min
c	ratio = 100.*icount/float(idim1*jdim1*kdim1)
c	write(*,*)'r2中大于',eee,'的所占比例为',ratio,'%'
c	stop
c
c     J = 0,jdim
c
      do 200 i = 1,idim1
	do 200 k = 1,kdim1
	c2 = gamma*qqj0(k,i,5,1)/qqj0(k,i,1,1)
	c  = sqrt(c2)

	t1 = sj(1,k,i,1)
	t2 = sj(1,k,i,2)
	t3 = sj(1,k,i,3)
	t4 = sj(1,k,i,4)
	ubar = qqj0(k,i,2,1)*t1+qqj0(k,i,3,1)*t2+qqj0(k,i,4,1)*t3
	ubar = ubar*t4
	r2(0,k,i) = abs(ubar)+c*t4
c
	c2 = gamma*qqj0(k,i,5,3)/qqj0(k,i,1,3)
	c  = sqrt(c2)

	t1 = sj(jdim,k,i,1)
	t2 = sj(jdim,k,i,2)
	t3 = sj(jdim,k,i,3)
	t4 = sj(jdim,k,i,4)
	ubar = qqj0(k,i,2,3)*t1+qqj0(k,i,3,3)*t2+qqj0(k,i,4,3)*t3
	ubar = ubar*t4
	r2(jdim,k,i) = abs(ubar)+c*t4
  200 continue
c
c     k = 0,kdim
c
      do 300 i = 1,idim1
	do 300 j = 1,jdim1
	c2 = gamma*qqk0(j,i,5,1)/qqk0(j,i,1,1)
	c  = sqrt(c2)

	t1 = sk(j,1,i,1)
	t2 = sk(j,1,i,2)
	t3 = sk(j,1,i,3)
	t4 = sk(j,1,i,4)
	ubar = qqk0(j,i,2,1)*t1+qqk0(j,i,3,1)*t2+qqk0(j,i,4,1)*t3
	ubar = ubar*t4
	r3(j,0,i) = abs(ubar)+c*t4
c
	c2 = gamma*qqk0(j,i,5,3)/qqk0(j,i,1,3)
	c  = sqrt(c2)

	t1 = sk(j,kdim,i,1)
	t2 = sk(j,kdim,i,2)
	t3 = sk(j,kdim,i,3)
	t4 = sk(j,kdim,i,4)
	ubar = qqk0(j,i,2,3)*t1+qqk0(j,i,3,3)*t2+qqk0(j,i,4,3)*t3
	ubar = ubar*t4
	r3(j,kdim,i) = abs(ubar)+c*t4
  300 continue
c
	do i = 1,idim1
	do j = 1,jdim1
	do k = 1,kdim1
	  r(j,k,i) = r3(j,k,i)/r2(j,k,i)
	  r(j,k,i) = 1.+r(j,k,i)**eta
	enddo
	enddo
	enddo
c
c     定义在cell-face
c	 j+1/2
	do 101 i = 1,idim1
	do 101 k = 1,kdim1
      r2(0,k,i) = bcj(k,i,1)*r2(0,k,i)+
     .            0.5*(1.-bcj(k,i,1))*(r2(0,k,i)+r2(1,k,i))
	do 102 j = 1,jdim2
	r2(j,k,i) = 0.5*(r(j,k,i)*r2(j,k,i)+r(j+1,k,i)*r2(j+1,k,i))
  102 continue
      r2(jdim1,k,i) = bcj(k,i,2)*r2(jdim,k,i)+
     .                0.5*(1.-bcj(k,i,2))*(r2(jdim1,k,i)+r2(jdim,k,i))
  101 continue
c
c	  k+1/2
c
	do i = 1,idim1
	do j = 1,jdim1
	do k = 1,kdim1
	  r(j,k,i) = r2(j,k,i)/r3(j,k,i)
	  r(j,k,i) = 1.+r(j,k,i)**eta
	enddo
	enddo
	enddo

      do 103 i = 1,idim1
	do 103 j = 1,jdim1
      r3(j,0,i) = bck(j,i,1)*r3(j,0,i)+
     .            0.5*(1.-bck(j,i,1))*(r3(j,0,i)+r3(j,1,i))
	do 104 k = 1,kdim2
      r3(j,k,i) = 0.5*(r(j,k,i)*r3(j,k,i)+r(j,k+1,i)*r3(j,k+1,i))
  104 continue
      r3(j,kdim1,i) = bck(j,i,2)*r3(j,kdim,i)+
     .                0.5*(1.-bck(j,i,2))*(r3(j,kdim1,i)+r3(j,kdim,i))
  103 continue
c
      if(i2d.eq.0)then
	do 401 i = 1,idim1
	do 401 k = 1,kdim1
	do 401 j = 1,jdim1
	c2 = gamma*qq(j,k,i,5)/qq(j,k,i,1)
	c  = sqrt(c2)
c
c     I 方向的
c
	t1 = 0.5*(si(j,k,i,1)+si(j,k,i+1,1))
	t2 = 0.5*(si(j,k,i,2)+si(j,k,i+1,2))
	t3 = 0.5*(si(j,k,i,3)+si(j,k,i+1,3))
	t4 = 0.5*(si(j,k,i,4)+si(j,k,i+1,4))
	ubar = qq(j,k,i,2)*t1+qq(j,k,i,3)*t2+qq(j,k,i,4)*t3
	ubar = ubar*t4
	r1(j,k,i) = abs(ubar)+c*t4
  401 continue
c
c     i = 0,idim
c
      do 400 k = 1,kdim1
	do 400 j = 1,jdim1
	c2 = gamma*qqi0(j,k,5,1)/qqi0(j,k,1,1)
	c  = sqrt(c2)

	t1 = si(j,k,1,1)
	t2 = si(j,k,1,2)
	t3 = si(j,k,1,3)
	t4 = si(j,k,1,4)
	ubar = qqi0(j,k,2,1)*t1+qqi0(j,k,3,1)*t2+qqi0(j,k,4,1)*t3
	ubar = ubar*t4
	r1(j,k,0) = abs(ubar)+c*t4
c
	c2 = gamma*qqi0(j,k,5,3)/qqi0(j,k,1,3)
	c  = sqrt(c2)

	t1 = si(j,k,idim,1)
	t2 = si(j,k,idim,2)
	t3 = si(j,k,idim,3)
	t4 = si(j,k,idim,4)
	ubar = qqi0(j,k,2,3)*t1+qqi0(j,k,3,3)*t2+qqi0(j,k,4,3)*t3
	ubar = ubar*t4
	r1(j,k,idim) = abs(ubar)+c*t4
  400 continue
c
c     i+1/2
c
      do 105 k = 1,kdim1
	do 105 j = 1,jdim1
      r1(j,k,0) = bci(j,k,1)*r1(j,k,0)+
     .            0.5*(1.-bci(j,k,1))*(r1(j,k,0)+r1(j,k,1))
	do 106 i = 1,idim2
	r1(j,k,i) = 0.5*(r1(j,k,i)+r1(j,k,i+1))
  106 continue
      r1(j,k,idim1) = bci(j,k,2)*r1(j,k,idim)+
     .                0.5*(1.-bci(j,k,2))*(r1(j,k,idim1)+r1(j,k,idim))
  105 continue
      endif
c********************************
c     J方向
c********************************      
c 
c    calculate "dz" at cell interface
c
      do 1001 m = 1,5
	do 1001 k = 1,kdim1
	do 1001 i = 1,idim1
	do 1001 j = 1,jdim2
	  dz(j,k,i,m) = q(j+1,k,i,m)-q(j,k,i,m)
 1001 continue
c
c     j = 0
c
	j = 0
	do 1002 m = 1,5
	do 1002 k = 1,kdim1
	do 1002 i = 1,idim1
	  dz(j,k,i,m) = (1.-bcj(k,i,1))*(q(1,k,i,m)-qj0(k,i,m,1))+
     .	            bcj(k,i,1)*qj0(k,i,m,2)
 1002 continue
c
c     j = -1
c
	j = -1
	do 1003 k = 1,kdim1
	do 1003 i = 1,idim1
c
	if(bcj(k,i,1).eq.0.)then
	  dz(j,k,i,1) = qj0(k,i,1,1)-qj0(k,i,1,2)
	  dz(j,k,i,2) = qj0(k,i,2,1)-qj0(k,i,2,2)
	  dz(j,k,i,3) = qj0(k,i,3,1)-qj0(k,i,3,2)
	  dz(j,k,i,4) = qj0(k,i,4,1)-qj0(k,i,4,2)
	  dz(j,k,i,5) = qj0(k,i,5,1)-qj0(k,i,5,2)
	endif
 1003 continue
c
c    j = jdim-1
c
      j = jdim1
	do 1004 m = 1,5
	do 1004 k = 1,kdim1
	do 1004 i = 1,idim1
	   dz(j,k,i,m) = (1.-bcj(k,i,2))*(qj0(k,i,m,3)-q(j,k,i,m))+
     .                  bcj(k,i,2)*qj0(k,i,m,4)
 1004 continue
c
c     j = jdim
c
	j = jdim
      do 1005 k = 1,kdim1
	do 1005 i = 1,idim1
	if(bcj(k,i,2).eq.0.)then
	  dz(j,k,i,1) = qj0(k,i,1,4)-qj0(k,i,1,3)
	  dz(j,k,i,2) = qj0(k,i,2,4)-qj0(k,i,2,3)
	  dz(j,k,i,3) = qj0(k,i,3,4)-qj0(k,i,3,3)
	  dz(j,k,i,4) = qj0(k,i,4,4)-qj0(k,i,4,3)
	  dz(j,k,i,5) = qj0(k,i,5,4)-qj0(k,i,5,3)
	endif
 1005 continue
c
c     calculate "dddz" at cell interface
c
	do 1006 m = 1,5
	do 1006 k = 1,kdim1
	do 1006 i = 1,idim1
	do 1006 j = 0,jdim1
	   d3z(j,k,i,m) = dz(j+1,k,i,m)-2.*dz(j,k,i,m)+dz(j-1,k,i,m)
 1006 continue
c
c     下面的做法详见文献"artificial dissipation and central difference schemes
c      for the Euler and Navier-Stokes equations",R.C.Swanson&E.Turkel,AIAA paper 1987,
c      87-1107,P58
c
c	do 1066 k = 1,kdim1
c	do 1066 i = 1,idim1
c	do 1066 m = 1,5
c	if(bcj(k,i,1).eq.1)then
c	   d3z(0,k,i,m) = 0.0
c	   d3z(1,k,i,m) = 0.0
c	elseif(bcj(k,i,2).eq.1)then
c	   d3z(jdim1,k,i,m) = 0.0
c	   d3z(jdim2,k,i,m) = 0.0	   
c	endif
c 1066 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     说明：psi被cancel掉 2007/10/22 by CYXU
c
c
c     计算emuj2,emuj4
c--------------------------------------------------------
c	do 1305	i = 1,idim1
c      do 1305 k = 1,kdim1
c	do 1306 j = 2,jdim2
c	t =  abs(qq(j+1,k,i,5)+2.0*qq(j,k,i,5)+qq(j-1,k,i,5)) 
c	psi(j,k,i) = abs(qq(j+1,k,i,5)-2.0*qq(j,k,i,5)+qq(j-1,k,i,5))/t    
c 1306 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     j = 1,0
c
c	if(bcj(k,i,1).eq.1.)then
c	t = abs(qq(2,k,i,5)+2.0*qq(1,k,i,5)+qq(1,k,i,5))
c	psi(1,k,i) = abs(qq(2,k,i,5)-2.0*qq(1,k,i,5)+qq(1,k,i,5))/t
c      psi(0,k,i) = -psi(1,k,i)
c	else
c	t =  abs(qq(2,k,i,5)+2.0*qq(1,k,i,5)+qqj0(k,i,5,1)) 
c	psi(1,k,i) = abs(qq(2,k,i,5)-2.0*qq(1,k,i,5)+qqj0(k,i,5,1))/t
c	t =  abs(qq(1,k,i,5)+2.0*qqj0(k,i,5,1)+qqj0(k,i,5,2)) 
c	psi(0,k,i) = abs(qq(1,k,i,5)-2.0*qqj0(k,i,5,1)+qqj0(k,i,5,2))/t
c	endif
c
c     j = jdim1,jdim
c     
c	if(bcj(k,i,2).eq.1.)then
c	t = abs(qq(jdim1,k,i,5)+2.0*qq(jdim1,k,i,5)+qq(jdim2,k,i,5))
c	psi(jdim1,k,i) = abs(qq(jdim1,k,i,5)-2.0*qq(jdim1,k,i,5)+
c     .                     qq(jdim2,k,i,5))/t
c      psi(jdim,k,i) = -psi(jdim1,k,i)
c	else
c	t =  abs(qq(jdim2,k,i,5)+2.0*qq(jdim1,k,i,5)+qqj0(k,i,5,3)) 
c	psi(jdim1,k,i) = abs(qq(jdim2,k,i,5)-2.0*qq(jdim1,k,i,5)+
c     .                     qqj0(k,i,5,3))/t
c	t =  abs(qq(jdim1,k,i,5)+2.0*qqj0(k,i,5,3)+qqj0(k,i,5,4)) 
c	psi(jdim,k,i) = abs(qq(jdim1,k,i,5)-2.0*qqj0(k,i,5,3)+
c     .                    qqj0(k,i,5,4))/t
c	endif          	
c 1305 continue
c----------------------------------------------------------
	do 1307 i = 1,idim1
	do 1307 k = 1,kdim1
	do 1307 j = 0,jdim1
c	t1 = psi(j+1,k,i)
c	t2 = psi(j,k,i)
c	t = max(t1,t2)
c	t3 = 0.0
c	t4 = eps4-eps2*t
c      t4 = eps4
c	emuj4(j,k,i) = max(t3,t4)*r2(j,k,i)
      emuj4(j,k,i) = eps4*r2(j,k,i)
c
 1307 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate "dqy"
c 
c      do 1107 m = 1,5
c	do 1107 k = 1,kdim1
c	do 1107 i = 1,idim1
c	do 1107 j = 1,jdim1
c	   t1    = emuj2(j,k,i)*dz(j,k,i,m)
c	   t2    = emuj2(j-1,k,i)*dz(j-1,k,i,m)
c	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
c	   t4    = emuj4(j-1,k,i)*d3z(j-1,k,i,m)
c	   dqy(j,k,i,m) = coef(j,k,i)*(t2-t1+t3-t4)
c 1107 continue
c     以前的
      do 1107 m = 1,5
	do 1107 k = 1,kdim1
	do 1107 i = 1,idim1
	do 1107 j = 2,jdim2
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j-1,k,i)*d3z(j-1,k,i,m)
	   dqy(j,k,i,m) = t3-t4
 1107 continue
      j = 1
      do 1108 m = 1,5
	do 1108 k = 1,kdim1
	do 1108 i = 1,idim1        
	if(bcj(k,i,1).eq.1.)then
c
c      下面4阶耗散项的添加方法详见R.C.Swanson & E.Turkel的文章P259：
c      "On Central-Difference and Upwind Schemes",JCP,1992.(下面关于
c      固壁边界条件的处理类同)
c
c	   dqy(j,k,i,m) = emuj4(j,k,i)*(dz(j+1,k,i,m)-dz(j,k,i,m))
         dqy(j,k,i,m) = emuj4(j,k,i)*d3z(j,k,i,m)-2.*emuj4(j-1,k,i)*
     .                  (dz(j,k,i,m)-dz(j-1,k,i,m))
	else	   
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j-1,k,i)*d3z(j-1,k,i,m)
	   dqy(j,k,i,m) = t3-t4
	endif
 1108 continue
      j = jdim1
      do 1110 m = 1,5
	do 1110 k = 1,kdim1
	do 1110 i = 1,idim1       
 	if(bcj(k,i,2).eq.1.)then
 	   dqy(j,k,i,m) = emuj4(j,k,i)*(dz(j,k,i,m)-dz(j-1,k,i,m))
 	else	   
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j-1,k,i)*d3z(j-1,k,i,m)
	   dqy(j,k,i,m) = t3-t4
	endif
 1110 continue
c*********************************
c     K方向
c*********************************
c 
c    calculate "dz" at cell interface
c
      do 2001 m = 1,5
	do 2001 j = 1,jdim1
	do 2001 i = 1,idim1
	do 2001 k = 1,kdim2
	  dz(j,k,i,m) = q(j,k+1,i,m)-q(j,k,i,m)
 2001 continue
c
c     k = 0
c
	k = 0
	do 2002 m = 1,5
	do 2002 j = 1,jdim1
	do 2002 i = 1,idim1
	  dz(j,k,i,m) = (1.-bck(j,i,1))*(q(j,1,i,m)-qk0(j,i,m,1))+
     .	             bck(j,i,1)*qk0(j,i,m,2)
 2002 continue
c
c     k = -1
c
	k = -1
	do 2003 j = 1,jdim1
	do 2003 i = 1,idim1
c
c     注意下面一行的做法不具有一般性，只是针对壁面B.C.,当用到其他
c     cell-face型B.C.时还需要进行修改
c
	if(bck(j,i,1).eq.0.)then
	  dz(j,k,i,1) = qk0(j,i,1,1)-qk0(j,i,1,2)
	  dz(j,k,i,2) = qk0(j,i,2,1)-qk0(j,i,2,2)
	  dz(j,k,i,3) = qk0(j,i,3,1)-qk0(j,i,3,2)
	  dz(j,k,i,4) = qk0(j,i,4,1)-qk0(j,i,4,2)
	  dz(j,k,i,5) = qk0(j,i,5,1)-qk0(j,i,5,2)
	endif
 2003 continue
c
c    k = kdim-1
c
      k = kdim1
	do 2004 m = 1,5
	do 2004 j = 1,jdim1
	do 2004 i = 1,idim1
	   dz(j,k,i,m) = (1.-bck(j,i,2))*(qk0(j,i,m,3)-q(j,k,i,m))+
     .                  bck(j,i,2)*qk0(j,i,m,4)
 2004 continue
c
c     k = kdim
c
	k = kdim
      do 2005 j = 1,jdim1
	do 2005 i = 1,idim1
	if(bck(j,i,2).eq.0.)then
	  dz(j,k,i,1) = qk0(j,i,1,4)-qk0(j,i,1,3)
	  dz(j,k,i,2) = qk0(j,i,2,4)-qk0(j,i,2,3)
	  dz(j,k,i,3) = qk0(j,i,3,4)-qk0(j,i,3,3)
	  dz(j,k,i,4) = qk0(j,i,4,4)-qk0(j,i,4,3)
	  dz(j,k,i,5) = qk0(j,i,5,4)-qk0(j,i,5,3)
	endif
 2005 continue
c
c     calculate "dddz" at cell interface
c
	do 2006 m = 1,5
	do 2006 j = 1,jdim1
	do 2006 i = 1,idim1
	do 2006 k = 0,kdim1
	   d3z(j,k,i,m) = dz(j,k+1,i,m)-2.*dz(j,k,i,m)+dz(j,k-1,i,m)
 2006 continue
c
c     下面的做法详见文献"artificial dissipation and central difference schemes
c      for the Euler and Navier-Stokes equations",R.C.Swanson&E.Turkel,AIAA paper 1987,
c      87-1107,P58
c
c	do 2066 j = 1,jdim1
c	do 2066 i = 1,idim1
c	do 2066 m = 1,5
c	if(bck(j,i,1).eq.1)then
c	   d3z(j,0,i,m) = 0.0
c	   d3z(j,1,i,m) = 0.0
c	elseif(bck(j,i,2).eq.1)then
c	   d3z(j,kdim1,i,m) = 0.0
c	   d3z(j,kdim2,i,m) = 0.0	   
c	endif
c 2066 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     计算emuj2,emuj4
c
c-------------------------------------------------------
c	do 2305	i = 1,idim1
c      do 2305 j = 1,jdim1
c	do 2306 k = 2,kdim2
c	t =  abs(qq(j,k+1,i,5)+2.0*qq(j,k,i,5)+qq(j,k-1,i,5)) 
c	psi(j,k,i) = abs(qq(j,k+1,i,5)-2.0*qq(j,k,i,5)+qq(j,k-1,i,5))/t    
c 2306 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     k = 1,0
c
c	if(bck(j,i,1).eq.1.)then
c	t = abs(qq(j,2,i,5)+2.0*qq(j,1,i,5)+qq(j,1,i,5))
c	psi(j,1,i) = abs(qq(j,2,i,5)-2.0*qq(j,1,i,5)+qq(j,1,i,5))/t
c      psi(j,0,i) = -psi(j,1,i)
c	else
c	t =  abs(qq(j,2,i,5)+2.0*qq(j,1,i,5)+qqk0(j,i,5,1)) 
c	psi(j,1,i) = abs(qq(j,2,i,5)-2.0*qq(j,1,i,5)+qqk0(j,i,5,1))/t
c	t =  abs(qq(j,1,i,5)+2.0*qqk0(j,i,5,1)+qqk0(j,i,5,2)) 
c	psi(j,0,i) = abs(qq(j,1,i,5)-2.0*qqk0(j,i,5,1)+qqk0(j,i,5,2))/t
c	endif
c
c     k = kdim1,kdim
c     
c	if(bck(j,i,2).eq.1.)then
c	t = abs(qq(j,kdim1,i,5)+2.0*qq(j,kdim1,i,5)+qq(j,kdim2,i,5))
c	psi(j,kdim1,i) = abs(qq(j,kdim1,i,5)-2.0*qq(j,kdim1,i,5)+
c     .                 qq(j,kdim2,i,5))/t
c      psi(j,kdim,i) = -psi(j,kdim1,i)
c	else
c	t =  abs(qq(j,kdim2,i,5)+2.0*qq(j,kdim1,i,5)+qqk0(j,i,5,3)) 
c	psi(j,kdim1,i) = abs(qq(j,kdim2,i,5)-2.0*qq(j,kdim1,i,5)+
c     .                 qqk0(j,i,5,3))/t
c	t =  abs(qq(j,kdim1,i,5)+2.0*qqk0(j,i,5,3)+qqk0(j,i,5,4)) 
c	psi(j,kdim,i) = abs(qq(j,kdim1,i,5)-2.0*qqk0(j,i,5,3)+
c     .                 qqk0(j,i,5,4))/t
c	endif          	
c 2305 continue
c--------------------------------------------------------
	do 2307 i = 1,idim1
	do 2307 j = 1,jdim1
	do 2307 k = 0,kdim1
c	t1 = psi(j,k+1,i)
c	t2 = psi(j,k,i)
c	t = max(t1,t2)
c	t3 = 0.0
c	t4 = eps4-eps2*t
c      t4 = eps4
c	emuj4(j,k,i) = max(t3,t4)*r3(j,k,i)
      emuj4(j,k,i) = eps4*r3(j,k,i)
 2307 continue
c
c     calculate "dqz"
c
c      do 2107 m = 1,5
c	do 2107 j = 1,jdim1
c	do 2107 i = 1,idim1
c	do 2107 k = 1,kdim1
c	   t1    = emuj2(j,k,i)*dz(j,k,i,m)
c	   t2    = emuj2(j,k-1,i)*dz(j,k-1,i,m)
c	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
c	   t4    = emuj4(j,k-1,i)*d3z(j,k-1,i,m)
c	   dqz(j,k,i,m) = coef(j,k,i)*(t2-t1+t3-t4)
c 2107 continue
c     以前的
      do 2107 m = 1,5
	do 2107 j = 1,jdim1
	do 2107 i = 1,idim1
	do 2107 k = 2,kdim2
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j,k-1,i)*d3z(j,k-1,i,m)
	   dqz(j,k,i,m) = t3-t4
 2107 continue
      k = 1
      do 2108 m = 1,5
	do 2108 j = 1,jdim1
	do 2108 i = 1,idim1
	if(bck(j,i,1).eq.1.)then
	   dqz(j,k,i,m) = emuj4(j,k,i)*(dz(j,k+1,i,m)-dz(j,k,i,m))
c	   dqz(j,k,i,m) = t2-t1+emuj4(j,k,i)*d3z(j,k,i,m)
	else
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j,k-1,i)*d3z(j,k-1,i,m)
	   dqz(j,k,i,m) = t3-t4
	endif
 2108 continue
      k = kdim1
      do 2110 m = 1,5
	do 2110 j = 1,jdim1
	do 2110 i = 1,idim1         
	if(bck(j,i,2).eq.1.)then
	   dqz(j,k,i,m) = emuj4(j,k,i)*(dz(j,k,i,m)-dz(j,k-1,i,m))
c	   dqz(j,k,i,m) = t2-t1+emuj4(j,k,i)*d3z(j,k,i,m)
	else
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j,k-1,i)*d3z(j,k-1,i,m)
	   dqz(j,k,i,m) = t3-t4
	endif
 2110 continue
c**********************************
c     I方向
c**********************************
      if(i2d.eq.0)then
c 
c    calculate "dz" at cell interface
c
      do 3001 m = 1,5
	do 3001 j = 1,jdim1
	do 3001 k = 1,kdim1
	do 3001 i = 1,idim2
	  dz(j,k,i,m) = q(j,k,i+1,m)-q(j,k,i,m)
 3001 continue
c
c     i = 0
c
	i = 0
	do 3002 m = 1,5
	do 3002 j = 1,jdim1
	do 3002 k = 1,kdim1
	  dz(j,k,i,m) = (1.-bci(j,k,1))*(q(j,k,1,m)-qi0(j,k,m,1))+
     .	             bci(j,k,1)*qi0(j,k,m,2)
 3002 continue
c
c     i = -1
c
	i = -1
	do 3003 j = 1,jdim1
	do 3003 k = 1,kdim1
c
	if(bci(j,k,1).eq.0.)then
	  dz(j,k,i,1) = qi0(j,k,1,1)-qi0(j,k,1,2)
	  dz(j,k,i,2) = qi0(j,k,2,1)-qi0(j,k,2,2)
	  dz(j,k,i,3) = qi0(j,k,3,1)-qi0(j,k,3,2)
	  dz(j,k,i,4) = qi0(j,k,4,1)-qi0(j,k,4,2)
	  dz(j,k,i,5) = qi0(j,k,5,1)-qi0(j,k,5,2)
	endif
 3003 continue
c
c    i = idim-1
c
      i = idim1
	do 3004 m = 1,5
	do 3004 j = 1,jdim1
	do 3004 k = 1,kdim1
	   dz(j,k,i,m) = (1.-bci(j,k,2))*(qi0(j,k,m,3)-q(j,k,i,m))+
     .                  bci(j,k,2)*qi0(j,k,m,4)
 3004 continue
c
c     i = idim
c
	i = idim
      do 3005 j = 1,jdim1
	do 3005 k = 1,kdim1
	if(bci(j,k,2).eq.0.)then
	  dz(j,k,i,1) = qi0(j,k,1,4)-qi0(j,k,1,3)
	  dz(j,k,i,2) = qi0(j,k,2,4)-qi0(j,k,2,3)
	  dz(j,k,i,3) = qi0(j,k,3,4)-qi0(j,k,3,3)
	  dz(j,k,i,4) = qi0(j,k,4,4)-qi0(j,k,4,3)
	  dz(j,k,i,5) = qi0(j,k,5,4)-qi0(j,k,5,3)
	endif
 3005 continue
c
c     calculate "dddz" at cell interface
c
	do 3006 m = 1,5
	do 3006 j = 1,jdim1
	do 3006 k = 1,kdim1
	do 3006 i = 0,idim1
	   d3z(j,k,i,m) = dz(j,k,i+1,m)-2.*dz(j,k,i,m)+dz(j,k,i-1,m)
 3006 continue
c
c     下面的做法详见文献"artificial dissipation and central difference schemes
c      for the Euler and Navier-Stokes equations",R.C.Swanson&E.Turkel,AIAA paper 1987,
c      87-1107,P58
c
c	do 3066 j = 1,jdim1
c	do 3066 k = 1,kdim1
c	do 3066 m = 1,5
c	if(bci(j,k,1).eq.1)then
c	   d3z(j,k,0,m) = 0.0
c	   d3z(j,k,1,m) = 0.0
c	elseif(bci(j,k,2).eq.1)then
c	   d3z(j,k,idim1,m) = 0.0
c	   d3z(j,k,idim2,m) = 0.0	   
c	endif
c 3066 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     计算emuj2,emuj4
c
c-----------------------------------------------------------
c	do 3305	j = 1,jdim1
c      do 3305 k = 1,kdim1
c	do 3306 i = 2,idim2
c	t =  abs(qq(j,k,i+1,5)+2.0*qq(j,k,i,5)+qq(j,k,i-1,5)) 
c	psi(j,k,i) = abs(qq(j,k,i+1,5)-2.0*qq(j,k,i,5)+qq(j,k,i-1,5))/t    
c 3306 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     i = 1,0
c
c	if(bci(j,k,1).eq.1.)then
c	t = abs(qq(j,k,2,5)+2.0*qq(j,k,1,5)+qq(j,k,1,5))
c	psi(j,k,1) = abs(qq(j,k,2,5)-2.0*qq(j,k,1,5)+qq(j,k,1,5))/t
c      psi(j,k,0) = -psi(j,k,1)
c	else
c	t =  abs(qq(j,k,2,5)+2.0*qq(j,k,1,5)+qqi0(j,k,5,1)) 
c	psi(j,k,1) = abs(qq(j,k,2,5)-2.0*qq(j,k,1,5)+qqi0(j,k,5,1))/t
c	t =  abs(qq(j,k,1,5)+2.0*qqi0(j,k,5,1)+qqi0(j,k,5,2)) 
c	psi(j,k,0) = abs(qq(j,k,1,5)-2.0*qqi0(j,k,5,1)+qqi0(j,k,5,2))/t
c	endif
c
c     i = idim1,idim
c     
c	if(bci(j,k,2).eq.1.)then
c	t = abs(qq(j,k,idim1,5)+2.0*qq(j,k,idim1,5)+qq(j,k,idim2,5))
c	psi(j,k,idim1) = abs(qq(j,k,idim1,5)-2.0*qq(j,k,idim1,5)+
c     .                     qq(j,k,idim2,5))/t
c      psi(j,k,idim) = -psi(j,k,idim1)
c	else
c	t =  abs(qq(j,k,idim2,5)+2.0*qq(j,k,idim1,5)+qqi0(j,k,5,3)) 
c	psi(j,k,idim1) = abs(qq(j,k,idim2,5)-2.0*qq(j,k,idim1,5)+
c     .                     qqi0(j,k,5,3))/t
c	t =  abs(qq(j,k,idim1,5)+2.0*qqi0(j,k,5,3)+qqi0(j,k,5,4)) 
c	psi(j,k,idim) = abs(qq(j,k,idim1,5)-2.0*qqi0(j,k,5,3)+
c     .                    qqi0(j,k,5,4))/t
c	endif          	
c 3305 continue
c----------------------------------------------------------------
	do 3307 j = 1,jdim1
	do 3307 k = 1,kdim1
	do 3307 i = 0,idim1
c	t1 = psi(j,k,i+1)
c	t2 = psi(j,k,i)
c	t  = max(t1,t2)
c	t3 = 0.0
c	t4 = eps4-eps2*t
c      t4 = eps4
c	emuj4(j,k,i) = max(t3,t4)*r1(j,k,i)
      emuj4(j,k,i) = eps4*r1(j,k,i)
 3307 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate "dqx"
c 
c      do 3107 m = 1,5
c	do 3107 k = 1,kdim1
c	do 3107 j = 1,jdim1
c	do 3107 i = 1,idim1
c	   t1    = emuj2(j,k,i)*dz(j,k,i,m)
c	   t2    = emuj2(j,k,i-1)*dz(j,k,i-1,m)
c	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
c	   t4    = emuj4(j,k,i-1)*d3z(j,k,i-1,m)
c	   dqx(j,k,i,m) = t2-t1+t3-t4
c 3107 continue
c	 以前的
      do 3107 m = 1,5
	do 3107 k = 1,kdim1
	do 3107 j = 1,jdim1
	do 3107 i = 2,idim2
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j,k,i-1)*d3z(j,k,i-1,m)
	   dqx(j,k,i,m) = t3-t4
 3107 continue
      i = 1
      do 3108 m = 1,5
	do 3108 k = 1,kdim1
	do 3108 j = 1,jdim1        
	if(bci(j,k,1).eq.1.)then
c
c      下面4阶耗散项的添加方法详见R.C.Swanson & E.Turkel的文章P259：
c      "On Central-Difference and Upwind Schemes",JCP,1992.(下面关于
c      固壁边界条件的处理类同)
c
	   dqx(j,k,i,m) = emuj4(j,k,i)*(dz(j,k,i+1,m)-dz(j,k,i,m))
c	   dqx(j,k,i,m) = t2-t1+emuj4(j,k,i)*d3z(j,k,i,m)
	else	   
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j,k,i-1)*d3z(j,k,i-1,m)
	   dqx(j,k,i,m) = t3-t4
	endif
 3108 continue
      i = idim1
      do 3110 m = 1,5
	do 3110 k = 1,kdim1
	do 3110 j = 1,jdim1       
	if(bci(j,k,2).eq.1.)then
	   dqx(j,k,i,m) = emuj4(j,k,i)*(dz(j,k,i,m)-dz(j,k,i-1,m))
c	   dqx(j,k,i,m) = t2-t1+emuj4(j,k,i)*d3z(j,k,i,m)
	else	   
	   t3    = emuj4(j,k,i)*d3z(j,k,i,m)
	   t4    = emuj4(j,k,i-1)*d3z(j,k,i-1,m)
	   dqx(j,k,i,m) = t3-t4
	endif
 3110 continue
c
      endif
c
      return
      end subroutine
#else
      subroutine dummy_dissipation
c**********************************************************************
c     Purpose: Provide a dummy routine to compile if HYBRID option is
c     not installed
c**********************************************************************
      return
      end subroutine
#endif


#ifdef HYBRID      
c
c    改进Jameson人工粘性中的开关函数
c	求解混合格式二进制开关函数
c    fai = 1 : 迎风格式
c    fai = 0 : 中心格式
c
	subroutine solve_fai4(jdim,kdim,idim,q,qj0,qk0,qi0,
     .                      bcj,bck,bci,ux,fai)
c
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
	dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
	dimension psi(0:jdim,0:kdim,0:idim),fai(jdim,kdim,idim)
	dimension ux(jdim-1,kdim-1,idim-1,9)
c
	common /twod/ i2d
c
      idim1 = idim-1
	jdim1 = jdim-1
	kdim1 = kdim-1
	jdim2 = jdim-2
	kdim2 = kdim-2
	idim2 = idim-2
c
c     cshock1,cshock2 :  当dq1>cshock1 & dq2>cshock2时调用迎风格式
c
      cshock1 = 0.007
      cshock2 = 1.0E-4
c      cshock2 = 0.001
	c       = 1.
	eps     = 1.0E-30
c
	do 1000 j = 1,jdim
	do 1000 k = 1,kdim
	do 1000 i = 1,idim
	   fai(j,k,i) = 0.0
 1000 continue
c
c     J 方向
c
	do 1002	i = 1,idim1
      do 1002 k = 1,kdim1
	do 1003 j = 2,jdim2
	t =  abs(q(j+1,k,i,5)+2.0*q(j,k,i,5)+q(j-1,k,i,5)) 
	psi(j,k,i) = abs(q(j+1,k,i,5)-2.0*q(j,k,i,5)+q(j-1,k,i,5))/t    
 1003 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     j = 1,0
c
	if(bcj(k,i,1).eq.1.)then
	t = abs(q(2,k,i,5)+2.0*q(1,k,i,5)+q(1,k,i,5))
	psi(1,k,i) = abs(q(2,k,i,5)-2.0*q(1,k,i,5)+q(1,k,i,5))/t
c	t = abs(qq(1,k,i,5)+2.0*qq(1,k,i,5)+qq(2,k,i,5))
c	psi(0,k,i) = abs(qq(1,k,5)-2.0*qq(1,k,5)+qq(2,k,i,5))/t
      psi(0,k,i) = psi(1,k,i)
c      psi(0,k,i) = -psi(1,k,i)
	else
	t =  abs(q(2,k,i,5)+2.0*q(1,k,i,5)+qj0(k,i,5,1)) 
	psi(1,k,i) = abs(q(2,k,i,5)-2.0*q(1,k,i,5)+qj0(k,i,5,1))/t
	t =  abs(q(1,k,i,5)+2.0*qj0(k,i,5,1)+qj0(k,i,5,2)) 
	psi(0,k,i) = abs(q(1,k,i,5)-2.0*qj0(k,i,5,1)+qj0(k,i,5,2))/t
	endif
c
c     j = jdim1,jdim
c     
	if(bcj(k,i,2).eq.1.)then
	t = abs(q(jdim1,k,i,5)+2.0*q(jdim1,k,i,5)+q(jdim2,k,i,5))
	psi(jdim1,k,i) = abs(q(jdim1,k,i,5)-2.0*q(jdim1,k,i,5)+
     .                     q(jdim2,k,i,5))/t
c	t = abs(qq(jdim2,k,i,5)+2.0*qq(jdim1,k,i,5)+qq(jdim1,k,i,5))
c	psi(jdim,k,i) = abs(qq(jdim2,k,i,5)-2.0*qq(jdim1,k,i,5)+
c     .                    qq(jdim1,k,i,5))/t
      psi(jdim,k,i) = psi(jdim1,k,i)
c      psi(jdim,k,i) = -psi(jdim1,k,i)
	else
	t =  abs(q(jdim2,k,i,5)+2.0*q(jdim1,k,i,5)+qj0(k,i,5,3)) 
	psi(jdim1,k,i) = abs(q(jdim2,k,i,5)-2.0*q(jdim1,k,i,5)+
     .                     qj0(k,i,5,3))/t
	t =  abs(q(jdim1,k,i,5)+2.0*qj0(k,i,5,3)+qj0(k,i,5,4)) 
	psi(jdim,k,i) = abs(q(jdim1,k,i,5)-2.0*qj0(k,i,5,3)+
     .                    qj0(k,i,5,4))/t
	endif          	
 1002 continue
c
c     求 fai
c
      do 1004 i = 1,idim1
	do 1004 k = 1,kdim1
	do 1004 j = 1,jdim
	t1 = psi(j-1,k,i)
	t2 = psi(j,k,i)
	dq1 = max(t1,t2)
c-------------------------------------------------------------
	if(j.eq.1)then
	du2  = (ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9))**2
	dq2  = exp(-c/(du2+eps))
	elseif(j.eq.jdim)then
	du2  = (ux(j-1,k,i,1)+ux(j-1,k,i,5)+ux(j-1,k,i,9))**2
	dq2  = exp(-c/(du2+eps))
	else
	du21 = (ux(j-1,k,i,1)+ux(j-1,k,i,5)+ux(j-1,k,i,9))**2
	du22 = (ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9))**2
	t3   = exp(-c/(du21+eps))
	t4   = exp(-c/(du22+eps))
	dq2  = max(t3,t4)
	endif
c-------------------------------------------------------------
	if(dq1.gt.cshock1.and.dq2.gt.cshock2)then
	fai(j,k,i) = 1.0
	endif
 1004 continue
c
c     K 方向
c
	do 2002	i = 1,idim1
      do 2002 j = 1,jdim1
	do 2003 k = 2,kdim2
	t =  abs(q(j,k+1,i,5)+2.0*q(j,k,i,5)+q(j,k-1,i,5)) 
	psi(j,k,i) = abs(q(j,k+1,i,5)-2.0*q(j,k,i,5)+q(j,k-1,i,5))/t    
 2003 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     k = 1,0
c
	if(bck(j,i,1).eq.1.)then
	t = abs(q(j,2,i,5)+2.0*q(j,1,i,5)+q(j,1,i,5))
	psi(j,1,i) = abs(q(j,2,i,5)-2.0*q(j,1,i,5)+q(j,1,i,5))/t
c	t = abs(qq(j,1,i,5)+2.0*qq(j,1,i,5)+qq(j,2,i,5))
c	psi(j,0,i) = abs(qq(j,1,5)-2.0*qq(j,1,5)+qq(j,2,i,5))/t
      psi(j,0,i) = psi(j,1,i)
c      psi(j,0,i) = -psi(j,1,i)
	else
	t =  abs(q(j,2,i,5)+2.0*q(j,1,i,5)+qk0(j,i,5,1)) 
	psi(j,1,i) = abs(q(j,2,i,5)-2.0*q(j,1,i,5)+qk0(j,i,5,1))/t
	t =  abs(q(j,1,i,5)+2.0*qk0(j,i,5,1)+qk0(j,i,5,2)) 
	psi(j,0,i) = abs(q(j,1,i,5)-2.0*qk0(j,i,5,1)+qk0(j,i,5,2))/t
	endif
c
c     k = kdim1,kdim
c     
	if(bck(j,i,2).eq.1.)then
	t = abs(q(j,kdim1,i,5)+2.0*q(j,kdim1,i,5)+q(j,kdim2,i,5))
	psi(j,kdim1,i) = abs(q(j,kdim1,i,5)-2.0*q(j,kdim1,i,5)+
     .                 q(j,kdim2,i,5))/t
c	t = abs(qq(j,kdim2,i,5)+2.0*qq(j,kdim1,i,5)+qq(j,kdim1,i,5))
c	psi(j,kdim,i) = abs(qq(j,kdim2,i,5)-2.0*qq(j,kdim1,i,5)+
c     .                 qq(j,kdim1,i,5))/t
      psi(j,kdim,i) = psi(j,kdim1,i)
c      psi(j,kdim,i) = -psi(j,kdim1,i)
	else
	t =  abs(q(j,kdim2,i,5)+2.0*q(j,kdim1,i,5)+qk0(j,i,5,3)) 
	psi(j,kdim1,i) = abs(q(j,kdim2,i,5)-2.0*q(j,kdim1,i,5)+
     .                 qk0(j,i,5,3))/t
	t =  abs(q(j,kdim1,i,5)+2.0*qk0(j,i,5,3)+qk0(j,i,5,4)) 
	psi(j,kdim,i) = abs(q(j,kdim1,i,5)-2.0*qk0(j,i,5,3)+
     .                 qk0(j,i,5,4))/t
	endif          	
 2002 continue
c
c     求 fai
c
      do 2004 i = 1,idim1
	do 2004 j = 1,jdim1
	do 2004 k = 1,kdim
	t1  = psi(j,k-1,i)
	t2  = psi(j,k,i)
	dq1 = max(t1,t2)
c-------------------------------------------------------------
	if(k.eq.1)then
	du2  = (ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9))**2
	dq2  = exp(-c/(du2+eps))
	elseif(k.eq.kdim)then
	du2  = (ux(j,k-1,i,1)+ux(j,k-1,i,5)+ux(j,k-1,i,9))**2
	dq2  = exp(-c/(du2+eps))
	else
	du21 = (ux(j,k-1,i,1)+ux(j,k-1,i,5)+ux(j,k-1,i,9))**2
	du22 = (ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9))**2
	t3   = exp(-c/(du21+eps))
	t4   = exp(-c/(du22+eps))
	dq2  = max(t3,t4)
	endif
c-------------------------------------------------------------
	if(dq1.gt.cshock1.and.dq2.gt.cshock2)then
	fai(j,k,i) = 1.0
	endif
 2004 continue
c
c     I 方向
c
	if(i2d.eq.0)then
	do 3002	j = 1,jdim1
      do 3002 k = 1,kdim1
	do 3003 i = 2,idim2
	t =  abs(q(j,k,i+1,5)+2.0*q(j,k,i,5)+q(j,k,i-1,5)) 
	psi(j,k,i) = abs(q(j,k,i+1,5)-2.0*q(j,k,i,5)+q(j,k,i-1,5))/t    
 3003 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     i = 1,0
c
	if(bci(j,k,1).eq.1.)then
	t = abs(q(j,k,2,5)+2.0*q(j,k,1,5)+q(j,k,1,5))
	psi(j,k,1) = abs(q(j,k,2,5)-2.0*q(j,k,1,5)+q(j,k,1,5))/t
      psi(j,k,0) = psi(j,k,1)
c      psi(j,k,0) = -psi(j,k,1)
	else
	t =  abs(q(j,k,2,5)+2.0*q(j,k,1,5)+qi0(j,k,5,1)) 
	psi(j,k,1) = abs(q(j,k,2,5)-2.0*q(j,k,1,5)+qi0(j,k,5,1))/t
	t =  abs(q(j,k,1,5)+2.0*qi0(j,k,5,1)+qi0(j,k,5,2)) 
	psi(j,k,0) = abs(q(j,k,1,5)-2.0*qi0(j,k,5,1)+qi0(j,k,5,2))/t
	endif
c
c     i = idim1,idim
c     
	if(bci(j,k,2).eq.1.)then
	t = abs(q(j,k,idim1,5)+2.0*q(j,k,idim1,5)+q(j,k,idim2,5))
	psi(j,k,idim1) = abs(q(j,k,idim1,5)-2.0*q(j,k,idim1,5)+
     .                     q(j,k,idim2,5))/t
      psi(j,k,idim) = psi(j,k,idim1)
c      psi(j,k,idim) = -psi(j,k,idim1)
	else
	t =  abs(q(j,k,idim2,5)+2.0*q(j,k,idim1,5)+qi0(j,k,5,3)) 
	psi(j,k,idim1) = abs(q(j,k,idim2,5)-2.0*q(j,k,idim1,5)+
     .                     qi0(j,k,5,3))/t
	t =  abs(q(j,k,idim1,5)+2.0*qi0(j,k,5,3)+qi0(j,k,5,4)) 
	psi(j,k,idim) = abs(q(j,k,idim1,5)-2.0*qi0(j,k,5,3)+
     .                    qi0(j,k,5,4))/t
	endif          	
 3002 continue
c
c     求 fai
c
      do 3004 j = 1,jdim1
	do 3004 k = 1,kdim1
	do 3004 i = 1,idim
	t1  = psi(j,k,i-1)
	t2  = psi(j,k,i)
	dq1 = max(t1,t2)
c-------------------------------------------------------------
	if(i.eq.1)then
	du2  = (ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9))**2
	dq2  = exp(-c/(du2+eps))
	elseif(i.eq.idim)then
	du2  = (ux(j,k,i-1,1)+ux(j,k,i-1,5)+ux(j,k,i-1,9))**2
	dq2  = exp(-c/(du2+eps))
	else
	du21 = (ux(j,k,i-1,1)+ux(j,k,i-1,5)+ux(j,k,i-1,9))**2
	du22 = (ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9))**2
	t3   = exp(-c/(du21+eps))
	t4   = exp(-c/(du22+eps))
	dq2  = max(t3,t4)
	endif
c-------------------------------------------------------------
	if(dq1.gt.cshock1.and.dq2.gt.cshock2)then
	fai(j,k,i) = 1.0
	endif
 3004 continue
      endif
c
      return
	end
c
c    L.Chern:
c    J.Fluid.Mech.(2006).vol.557.pp.29-61 by D.J.Hill et.al.
c    'LES & multiscale modelling of a R-M instability with shock'
c    fai = 1 : 迎风格式
c    fai = 0 : 中心格式
c
	subroutine solve_fai5(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,vol,
     .                      bcj,bck,bci,fai)
c
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4),
     +qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),sj(jdim,kdim,idim-1,5),
     +sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
	dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
	dimension psi1(0:jdim,0:kdim,0:idim)
	dimension psi2(0:jdim,0:kdim,0:idim)
	dimension dltx(jdim-1,kdim-1,idim-1)
	dimension fai(jdim,kdim,idim)
c
	common /twod/ i2d
c
      idim1 = idim-1
	jdim1 = jdim-1
	kdim1 = kdim-1
	jdim2 = jdim-2
	kdim2 = kdim-2
	idim2 = idim-2
c
c     cshock : 当dq大于cshock时调用迎风格式
c
	idltxc=0
	if(idltxc.ne.1)goto 901
c-------------------------------------------------------------
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
              deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
              deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
              delta = max(deltaj,deltak)
c
              if( i2d .ne. 1 ) then
                 dltx(j,k,i) = max(delta,deltai)
              end if        
            enddo
          enddo
        enddo
c------------------------------------------------------------
901   cshock = 0.002
c	epsc=2.5e3
c
	do 1000 j = 1,jdim
	do 1000 k = 1,kdim
	do 1000 i = 1,idim
	   fai(j,k,i) = 0.0
 1000 continue
c
c     J 方向
c
	do 1002	i = 1,idim1
      do 1002 k = 1,kdim1
	do 1003 j = 2,jdim2
	t =  (q(j+1,k,i,5)+2.0*q(j,k,i,5)+q(j-1,k,i,5)) 
	tt=  (q(j+1,k,i,1)+2.0*q(j,k,i,1)+q(j-1,k,i,1))
	
	psi1(j,k,i) = (q(j+1,k,i,5)-2.0*q(j,k,i,5)+q(j-1,k,i,5))/t 
	psi2(j,k,i) = (q(j+1,k,i,1)-2.0*q(j,k,i,1)+q(j-1,k,i,1))/tt   
 1003 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     j = 1,0
c
	if(bcj(k,i,1).eq.1.)then
	t = (q(2,k,i,5)+2.0*q(1,k,i,5)+q(1,k,i,5))
	tt= (q(2,k,i,1)+2.0*q(1,k,i,1)+q(1,k,i,1))
	psi1(1,k,i) = (q(2,k,i,5)-2.0*q(1,k,i,5)+q(1,k,i,5))/t
	psi2(1,k,i) = (q(2,k,i,1)-2.0*q(1,k,i,1)+q(1,k,i,1))/tt
c	t = abs(qq(1,k,i,5)+2.0*qq(1,k,i,5)+qq(2,k,i,5))
c	psi(0,k,i) = abs(qq(1,k,5)-2.0*qq(1,k,5)+qq(2,k,i,5))/t
c      psi(0,k,i) = psi(1,k,i)
      psi1(0,k,i) = psi1(1,k,i)
	psi2(0,k,i) = psi2(1,k,i)
	else
	t =  (q(2,k,i,5)+2.0*q(1,k,i,5)+qj0(k,i,5,1)) 
	tt=  (q(2,k,i,1)+2.0*q(1,k,i,1)+qj0(k,i,1,1)) 
	psi1(1,k,i) = (q(2,k,i,5)-2.0*q(1,k,i,5)+qj0(k,i,5,1))/t
	psi2(1,k,i) = (q(2,k,i,1)-2.0*q(1,k,i,1)+qj0(k,i,1,1))/tt
	t =  (q(1,k,i,5)+2.0*qj0(k,i,5,1)+qj0(k,i,5,2)) 
	tt=  (q(1,k,i,1)+2.0*qj0(k,i,1,1)+qj0(k,i,1,2))
	psi1(0,k,i) = (q(1,k,i,5)-2.0*qj0(k,i,5,1)+qj0(k,i,5,2))/t
	psi2(0,k,i) = (q(1,k,i,1)-2.0*qj0(k,i,1,1)+qj0(k,i,1,2))/tt
	endif
c
c     j = jdim1,jdim
c     
	if(bcj(k,i,2).eq.1.)then
	t = (q(jdim1,k,i,5)+2.0*q(jdim1,k,i,5)+q(jdim2,k,i,5))
	tt= (q(jdim1,k,i,1)+2.0*q(jdim1,k,i,1)+q(jdim2,k,i,1))
	psi1(jdim1,k,i) = (q(jdim1,k,i,5)-2.0*q(jdim1,k,i,5)+
     .                     q(jdim2,k,i,5))/t
	psi2(jdim1,k,i) = (q(jdim1,k,i,1)-2.0*q(jdim1,k,i,1)+
     .					 q(jdim2,k,i,1))/tt
c	t = abs(qq(jdim2,k,i,5)+2.0*qq(jdim1,k,i,5)+qq(jdim1,k,i,5))
c	psi(jdim,k,i) = abs(qq(jdim2,k,i,5)-2.0*qq(jdim1,k,i,5)+
c     .                    qq(jdim1,k,i,5))/t
c      psi(jdim,k,i) = psi(jdim1,k,i)
      psi1(jdim,k,i) = psi1(jdim1,k,i)
	psi2(jdim,k,i) = psi2(jdim1,k,i)
	else
	t =  (q(jdim2,k,i,5)+2.0*q(jdim1,k,i,5)+qj0(k,i,5,3))
	tt=  (q(jdim2,k,i,1)+2.0*q(jdim1,k,i,1)+qj0(k,i,1,3)) 
	psi1(jdim1,k,i) = (q(jdim2,k,i,5)-2.0*q(jdim1,k,i,5)+
     .                     qj0(k,i,5,3))/t
	psi2(jdim1,k,i) = (q(jdim2,k,i,1)-2.0*q(jdim1,k,i,1)+
     .                     qj0(k,i,1,3))/tt
	t =  (q(jdim1,k,i,5)+2.0*qj0(k,i,5,3)+qj0(k,i,5,4)) 
	tt=  (q(jdim1,k,i,1)+2.0*qj0(k,i,1,3)+qj0(k,i,1,4)) 
	psi1(jdim,k,i) = (q(jdim1,k,i,5)-2.0*qj0(k,i,5,3)+
     .                    qj0(k,i,5,4))/t
	psi2(jdim,k,i) = (q(jdim1,k,i,1)-2.0*qj0(k,i,1,3)+
     .                    qj0(k,i,1,4))/tt
	endif          	
 1002 continue
c
c     求 fai
c
      do 1004 i = 1,idim
	do 1004 k = 1,kdim
	do 1004 j = 1,jdim
	t1 = (psi1(j-1,k,i))
	t2 = (psi1(j,k,i))
	dq1 = max(abs(t1),abs(t2))
	aq1 = .5*(t1+t2)
	t1 = (psi2(j-1,k,i))
	t2 = (psi2(j,k,i))
	dq2 = max(abs(t1),abs(t2))
	aq2 = .5*(t1+t2)
c
	irp=0
	if((aq1*aq2).gt.0)irp=1
	if(dq1.gt.cshock.and.dq2.gt.cshock.and.irp.eq.1)then
	fai(j,k,i) = 1.0
	endif
 1004 continue
c
c     K 方向
c
	do 2002	i = 1,idim1
      do 2002 j = 1,jdim1
	do 2003 k = 2,kdim2
	t =  (q(j,k+1,i,5)+2.0*q(j,k,i,5)+q(j,k-1,i,5)) 
	tt=  (q(j,k+1,i,1)+2.0*q(j,k,i,1)+q(j,k-1,i,1))
	psi1(j,k,i) = (q(j,k+1,i,5)-2.0*q(j,k,i,5)+q(j,k-1,i,5))/t
	psi2(j,k,i) = (q(j,k+1,i,1)-2.0*q(j,k,i,1)+q(j,k-1,i,1))/tt    
 2003 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     k = 1,0
c
	if(bck(j,i,1).eq.1.)then
	t = (q(j,2,i,5)+2.0*q(j,1,i,5)+q(j,1,i,5))
	tt= (q(j,2,i,1)+2.0*q(j,1,i,1)+q(j,1,i,1))
	psi1(j,1,i) = (q(j,2,i,5)-2.0*q(j,1,i,5)+q(j,1,i,5))/t
	psi2(j,1,i) = (q(j,2,i,1)-2.0*q(j,1,i,1)+q(j,1,i,1))/tt
c	t = abs(qq(j,1,i,5)+2.0*qq(j,1,i,5)+qq(j,2,i,5))
c	psi(j,0,i) = abs(qq(j,1,5)-2.0*qq(j,1,5)+qq(j,2,i,5))/t
c      psi(j,0,i) = psi(j,1,i)
      psi1(j,0,i) = psi1(j,1,i)
	psi2(j,0,i) = psi2(j,1,i)
	else
	t =  (q(j,2,i,5)+2.0*q(j,1,i,5)+qk0(j,i,5,1)) 
	tt=  (q(j,2,i,1)+2.0*q(j,1,i,1)+qk0(j,i,1,1)) 
	psi1(j,1,i) = (q(j,2,i,5)-2.0*q(j,1,i,5)+qk0(j,i,5,1))/t
	psi2(j,1,i) = (q(j,2,i,1)-2.0*q(j,1,i,1)+qk0(j,i,1,1))/tt
	t =  (q(j,1,i,5)+2.0*qk0(j,i,5,1)+qk0(j,i,5,2)) 
	tt=  (q(j,1,i,1)+2.0*qk0(j,i,1,1)+qk0(j,i,1,2)) 
	psi1(j,0,i) = (q(j,1,i,5)-2.0*qk0(j,i,5,1)+qk0(j,i,5,2))/t
	psi2(j,0,i) = (q(j,1,i,1)-2.0*qk0(j,i,1,1)+qk0(j,i,1,2))/t
	endif
c
c     k = kdim1,kdim
c     
	if(bck(j,i,2).eq.1.)then
	t = (q(j,kdim1,i,5)+2.0*q(j,kdim1,i,5)+q(j,kdim2,i,5))
	tt= (q(j,kdim1,i,1)+2.0*q(j,kdim1,i,1)+q(j,kdim2,i,1))
	psi1(j,kdim1,i) = (q(j,kdim1,i,5)-2.0*q(j,kdim1,i,5)+
     .                 q(j,kdim2,i,5))/t
	psi2(j,kdim1,i) = (q(j,kdim1,i,1)-2.0*q(j,kdim1,i,1)+
     .                 q(j,kdim2,i,1))/tt
c	t = abs(qq(j,kdim2,i,5)+2.0*qq(j,kdim1,i,5)+qq(j,kdim1,i,5))
c	psi(j,kdim,i) = abs(qq(j,kdim2,i,5)-2.0*qq(j,kdim1,i,5)+
c     .                 qq(j,kdim1,i,5))/t
c      psi(j,kdim,i) = psi(j,kdim1,i)
      psi1(j,kdim,i) = psi1(j,kdim1,i)
	psi2(j,kdim,i) = psi2(j,kdim1,i)
	else
	t =  (q(j,kdim2,i,5)+2.0*q(j,kdim1,i,5)+qk0(j,i,5,3))
	tt=  (q(j,kdim2,i,1)+2.0*q(j,kdim1,i,1)+qk0(j,i,1,3)) 
	psi1(j,kdim1,i) = (q(j,kdim2,i,5)-2.0*q(j,kdim1,i,5)+
     .                 qk0(j,i,5,3))/t
	psi2(j,kdim1,i) = (q(j,kdim2,i,1)-2.0*q(j,kdim1,i,1)+
     .                 qk0(j,i,1,3))/tt
	t =  (q(j,kdim1,i,5)+2.0*qk0(j,i,5,3)+qk0(j,i,5,4))
	tt=  (q(j,kdim1,i,1)+2.0*qk0(j,i,1,3)+qk0(j,i,1,4))
	psi1(j,kdim,i) = (q(j,kdim1,i,5)-2.0*qk0(j,i,5,3)+
     .                 qk0(j,i,5,4))/t
	psi2(j,kdim,i) = (q(j,kdim1,i,1)-2.0*qk0(j,i,1,3)+
     .                 qk0(j,i,1,4))/tt
	endif          	
 2002 continue
c
c     求 fai
c
      do 2004 i = 1,idim1
	do 2004 j = 1,jdim1
	do 2004 k = 1,kdim
	t1 = (psi1(j,k-1,i))
	t2 = (psi1(j,k,i))
	dq1 = max(abs(t1),abs(t2))
	aq1 = .5*(t1+t2)
	t1 = (psi2(j,k-1,i))
	t2 = (psi2(j,k,i))
	dq2 = max(abs(t1),abs(t2))
	aq2 = .5*(t1+t2)
c
	irp=0
	if((aq1*aq2).gt.0)irp=1
	if(dq1.gt.cshock.and.dq2.gt.cshock.and.irp.eq.1)then
	fai(j,k,i) = 1.0
	endif
 2004 continue
c
c     I 方向
c
	if(i2d.eq.0)then
	do 3002	j = 1,jdim1
      do 3002 k = 1,kdim1
	do 3003 i = 2,idim2
	t =  (q(j,k,i+1,5)+2.0*q(j,k,i,5)+q(j,k,i-1,5)) 
	tt=  (q(j,k,i+1,1)+2.0*q(j,k,i,1)+q(j,k,i-1,1))
	psi1(j,k,i) = (q(j,k,i+1,5)-2.0*q(j,k,i,5)+q(j,k,i-1,5))/t    
	psi2(j,k,i) = (q(j,k,i+1,1)-2.0*q(j,k,i,1)+q(j,k,i-1,1))/tt   
 3003 continue
c
c     注意下面关于边界的处理方法不具有一般性，只是针对壁面B.C.,当用到
c     其他cell-face型B.C.时还需要进行修改
c
c
c     i = 1,0
c
	if(bci(j,k,1).eq.1.)then
	t = (q(j,k,2,5)+2.0*q(j,k,1,5)+q(j,k,1,5))
	tt= (q(j,k,2,1)+2.0*q(j,k,1,1)+q(j,k,1,1))
	psi1(j,k,1) = (q(j,k,2,5)-2.0*q(j,k,1,5)+q(j,k,1,5))/t
	psi2(j,k,1) = (q(j,k,2,1)-2.0*q(j,k,1,1)+q(j,k,1,1))/tt
c      psi(j,k,0) = psi(j,k,1)
      psi1(j,k,0) = psi1(j,k,1)
	psi2(j,k,0) = psi2(j,k,1)
	else
	t =  (q(j,k,2,5)+2.0*q(j,k,1,5)+qi0(j,k,5,1)) 
	tt=  (q(j,k,2,1)+2.0*q(j,k,1,1)+qi0(j,k,1,1))
	psi1(j,k,1) = (q(j,k,2,5)-2.0*q(j,k,1,5)+qi0(j,k,5,1))/t
	psi2(j,k,1) = (q(j,k,2,1)-2.0*q(j,k,1,1)+qi0(j,k,1,1))/tt
	t =  (q(j,k,1,5)+2.0*qi0(j,k,5,1)+qi0(j,k,5,2)) 
	tt=  (q(j,k,1,1)+2.0*qi0(j,k,1,1)+qi0(j,k,1,2)) 
	psi1(j,k,0) = (q(j,k,1,5)-2.0*qi0(j,k,5,1)+qi0(j,k,5,2))/t
	psi2(j,k,0) = (q(j,k,1,1)-2.0*qi0(j,k,1,1)+qi0(j,k,1,2))/tt
	endif
c
c     i = idim1,idim
c     
	if(bci(j,k,2).eq.1.)then
	t = (q(j,k,idim1,5)+2.0*q(j,k,idim1,5)+q(j,k,idim2,5))
	tt= (q(j,k,idim1,1)+2.0*q(j,k,idim1,1)+q(j,k,idim2,1))
	psi1(j,k,idim1) = (q(j,k,idim1,5)-2.0*q(j,k,idim1,5)+
     .                     q(j,k,idim2,5))/t
	psi2(j,k,idim1) = (q(j,k,idim1,1)-2.0*q(j,k,idim1,1)+
     .                     q(j,k,idim2,1))/tt
c      psi(j,k,idim) = psi(j,k,idim1)
      psi1(j,k,idim) = psi1(j,k,idim1)
	psi2(j,k,idim) = psi2(j,k,idim1)
	else
	t =  (q(j,k,idim2,5)+2.0*q(j,k,idim1,5)+qi0(j,k,5,3)) 
	tt=  (q(j,k,idim2,1)+2.0*q(j,k,idim1,1)+qi0(j,k,1,3)) 
	psi1(j,k,idim1) = (q(j,k,idim2,5)-2.0*q(j,k,idim1,5)+
     .                     qi0(j,k,5,3))/t
	psi2(j,k,idim1) = (q(j,k,idim2,1)-2.0*q(j,k,idim1,1)+
     .                     qi0(j,k,1,3))/tt
	t =  (q(j,k,idim1,5)+2.0*qi0(j,k,5,3)+qi0(j,k,5,4))
	tt=  (q(j,k,idim1,1)+2.0*qi0(j,k,1,3)+qi0(j,k,1,4))  
	psi1(j,k,idim) = (q(j,k,idim1,5)-2.0*qi0(j,k,5,3)+
     .                    qi0(j,k,5,4))/t
	psi2(j,k,idim) = (q(j,k,idim1,1)-2.0*qi0(j,k,1,3)+
     .                    qi0(j,k,1,4))/tt
	endif          	
 3002 continue
c
c     求 fai
c
      do 3004 j = 1,jdim1
	do 3004 k = 1,kdim1
	do 3004 i = 1,idim
	t1 = (psi1(j,k,i-1))
	t2 = (psi1(j,k,i))
	dq1 = max(abs(t1),abs(t2))
	aq1 = .5*(t1+t2)
	t1 = (psi2(j,k,i-1))
	t2 = (psi2(j,k,i))
	dq2 = max(abs(t1),abs(t2))
	aq2 = .5*(t1+t2)
c
	irp=0
	if((aq1*aq2).gt.0)irp=1
	if(dq1.gt.cshock.and.dq2.gt.cshock.and.irp.eq.1)then
	fai(j,k,i) = 1.0
	endif
 3004 continue
      endif
c
      return
      end
c
c     输出耗散加权函数fai
c
      subroutine faiout(jdim,kdim,idim,x,y,z,ffai,myid)
c
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension ffai(jdim,kdim,idim)
      character*2 tchar2
c
      common /twod/ i2d
c
c#ifdef HYBRID
c#ifndef DIST_MPI
      write(tchar2(1:2),'(i2.2)')myid
      ifile = 124+myid
      open(ifile,file='artif'//tchar2//'.dat',
     .status='unknown') 
c#endif
c#endif
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
      
      if(i2d.eq.1)then
      i = 1
      write(ifile,*)'variables= "x" "z" "fai"'
      write(ifile,*)'zone f= point,j=',jdim,'k=',kdim
      do 1001 k = 1,kdim
      do 1001 j = 1,jdim
      write(ifile,*)x(j,k,i),z(j,k,i),ffai(j,k,i)
 1001 continue
c
      else
c
      write(ifile,*)'variables= "x" "y" "z" "fai"'
      write(ifile,*)'zone f= point,i=',idim,'j=',jdim,'k=',kdim
      do 1002 k = 1,kdim
      do 1002 j = 1,jdim
      do 1002 i = 1,idim
      write(ifile,*)x(j,k,i),y(j,k,i),z(j,k,i),ffai(j,k,i)
 1002 continue
c
      endif
      return
      end
#else
      subroutine dummy_solve_fai
c**********************************************************************
c     Purpose: Provide a dummy routine to compile if HYBRID option is
c     not installed
c**********************************************************************
      return
      end
#endif
      end module
