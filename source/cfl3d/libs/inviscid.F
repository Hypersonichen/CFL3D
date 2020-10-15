c
c     2阶精度
c
      subroutine hybrid2(xkap,n,x1,x2,xc,npts,leq,fai)                     
c***********************************************************************
c     Purpose:  Perform monotone interpolations to the interfaces
c     of the cells.
c***********************************************************************
c
      dimension x1(n),x2(n),xc(n),fai(n) 
c
c     on input: 
c        x1 is left interface gradient del(-)
c        x2 is right interface gradient del(+) of cell
c        xc is cell center value
c     on output:  
c        x2 is gradient for right interface construction
c        x1 is gradient for left interface construction
c     for unlimited fully-upwind (kappa=-1;iswi=0)
c        x1(output)=x2(input)
c        x2(output)=x1(input)
c
      ibias = 1
      eps   = 1.0e-06
      phi   = (3.-xkap)/(1.-xkap+eps)
      eps2  = 0.5e-06
c
c     tuned k=1/3 limiter - Spekreijse - Venkat                   AIAA-90-0429
c
      delx    = 10./float(npts)
      eps2    = delx**3
c
      do 7500 izz=1,n
      if(fai(izz).eq.1.0)then
c
c     选用迎风格式
c
      t3      = x1(izz)*x1(izz)
      t4      = x2(izz)*x2(izz)
      t5      = x1(izz)*x2(izz)
      t6      = x1(izz)+x2(izz)
      term    = 0.5*(t5+eps2)/(2.*(t3 + t4) - t5 + 3.*eps2)
      x2(izz) = (x2(izz)+t6)*term 
      x1(izz) = (x1(izz)+t6)*term 
      else
c
c     选用中心格式
c
	x22 = x2(izz)
	x11	= x1(izz)
	x2(izz) = 0.125*(x11+3.0*x22)
	x1(izz) = 0.125*(3.0*x11+x22)
      endif
 7500 continue
c
c      bias gradients by average values (leq = 1 or 5)
c
      if (ibias.eq.1) then
      if (leq.eq.1 .or. leq.eq.5) then
cdir$ ivdep
      do 9000 izz=1,n
      x2(izz) = x2(izz)*xc(izz)
      x1(izz) = x1(izz)*xc(izz)
 9000 continue
      end if
      end if
      return
      end
c
      subroutine fhybrid2(ax,ay,az,area,at,f,qr,ql,n,nvtq,fai)     
c
      dimension ax(n),ay(n),az(n),area(n),at(n),
     .          qr(nvtq,5),ql(nvtq,5),f(nvtq,5)
	real,dimension(n):: fai
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,qiv(5)
     .        ,tur10,tur20
c
c      delta q across faces t(1-5)
c
      x1   = gamma/gm1
      c1   = 1.e0/gm1
      eps4 = 0.e0
c
      do 1000 i=1,n 
      t1 = qr(i,1)-ql(i,1)
      t2 = qr(i,2)-ql(i,2)
      t3 = qr(i,3)-ql(i,3)
      t4 = qr(i,4)-ql(i,4)
c
c      pressure and enthalpy
c
      t16     = 1.e0/qr(i,1)
      t5      = qr(i,5)
      qr(i,5) = x1*qr(i,5)*t16+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                              +qr(i,4)*qr(i,4)) 
c
      t15     = 1.e0/ql(i,1)
      t19     = ql(i,5) 
      ql(i,5) = x1*ql(i,5)*t15+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                              +ql(i,4)*ql(i,4)) 
c
c      unsplit contributions  f(r)+f(l)
c
      t18 = ax(i)*qr(i,2)+ay(i)*qr(i,3)+az(i)*qr(i,4)+at(i) 
      t17 = ax(i)*ql(i,2)+ay(i)*ql(i,3)+az(i)*ql(i,4)+at(i) 
      t6  = t18*qr(i,1)
      t7  = t17*ql(i,1)
      f1  = t6+t7
      f2  = t6*qr(i,2)+t7*ql(i,2)
      f3  = t6*qr(i,3)+t7*ql(i,3)
      f4  = t6*qr(i,4)+t7*ql(i,4)
      f5  = t6*qr(i,5)+t7*ql(i,5)
      t8  = t5+t19 
      f2  = f2+ax(i)*t8
      f3  = f3+ay(i)*t8
      f4  = f4+az(i)*t8
      f5  = f5-at(i)*t8
c
c      roe averaged variables
c
      t6 = qr(i,1)*t15
      t7 = sqrt(t6) 
      t6 = 1.e0/(1.e0+t7) 
      t8 = t7*t6
c
c     average density
c
      qr(i,1) = ql(i,1)*t7
c
c      u,v,w,h average
c
      t9  = ql(i,2)*t6+qr(i,2)*t8
      t10 = ql(i,3)*t6+qr(i,3)*t8 
      t11 = ql(i,4)*t6+qr(i,4)*t8 
      t12 = ql(i,5)*t6+qr(i,5)*t8 
c
c      extract sound speed
c
      t6  = (t9*t9+t10*t10+t11*t11)*0.5e0
      t7  = gm1*(t12-t6) 
      t8  = sqrt(t7) 
c
      t13 = t9*ax(i)+t10*ay(i)+t11*az(i)
c
c      the variables in t are as follows
c      1-4        delta q1-q4
c      6          q2a
c      7          c2a
c      8          ca
c      9,10,11,12 ua,va,wa,ha
c      13         ubara
c      14,15,16   alpha(1),alpha(4),alpha(5)
c      18,19,17   lambda(1),lambda(4),lambda(5)
c
c      rhoa*delta(ubar) , delta(p)/c2a
c
      ql(i,1) = qr(i,1)*(t18-t17) 
      ql(i,2) = (t5-t19)/t7 
c
      t18 = t13+at(i) 
      t18 = abs(t18)+eps4 
      t19 = t13+at(i)+t8
      t19 = abs(t19)+eps4 
      t17 = t13+at(i)-t8
      t17 = abs(t17)+eps4 
c
      t14 = t18*(t1-ql(i,2))
      t15 = .5e0*(ql(i,2)+ql(i,1)/t8) 
      t16 = (ql(i,2)-t15)*t17 
      t15 = t15*t19 
c
      qr(i,2) = t18*(t2*qr(i,1)-ax(i)*ql(i,1))
      qr(i,3) = t18*(t3*qr(i,1)-ay(i)*ql(i,1))
      qr(i,4) = t18*(t4*qr(i,1)-az(i)*ql(i,1))
      qr(i,5) = t9*qr(i,2)+t10*qr(i,3)+t11*qr(i,4)
c
      ql(i,1) = t14+t15+t16 
      ql(i,2) = t8*(t15-t16)
c
      f1 = f1+fai(i)*(-ql(i,1)) 
      f2 = f2+fai(i)*(-ql(i,1)*t9-ax(i)*ql(i,2)-qr(i,2))
      f3 = f3+fai(i)*(-ql(i,1)*t10-ay(i)*ql(i,2)-qr(i,3)) 
      f4 = f4+fai(i)*(-ql(i,1)*t11-az(i)*ql(i,2)-qr(i,4)) 
      f5 = f5+fai(i)*(-ql(i,1)*t12-t13*ql(i,2)-qr(i,5)+t7*c1*t14)
c
c      include factor one-half area
c
      t7     = .5e0*area(i) 
      f(i,1) = t7*f1
      f(i,2) = t7*f2
      f(i,3) = t7*f3
      f(i,4) = t7*f4
      f(i,5) = t7*f5
 1000 continue
      return
      end
      subroutine fhybrid3(ax,ay,az,area,at,f,qr,ql,n,nvtq,fai)                      
c***********************************************************************
c     Purpose:  Compute flux at the interface using hybrid method
c     given the left and right states at the interface.
c***********************************************************************
c
      dimension ax(n),ay(n),az(n),area(n),at(n),
     .          qr(nvtq,5),ql(nvtq,5),f(nvtq,5)
	dimension fai(n)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,qiv(5)
     .        ,tur10,tur20
c
c      delta q across faces t(1-5)
c
      x1   = gamma/gm1
      c1   = 1.e0/gm1
      eps4 = 0.e0	
	c24 = 1.0/24.0
c
	iskew = 1
c
	do 2000 i = 1,n
      if(fai(i).eq.1.0)then
c
c     迎风格式
c
c
c      delta q across faces t(1-5)
c
      t1 = qr(i,1)-ql(i,1)
      t2 = qr(i,2)-ql(i,2)
      t3 = qr(i,3)-ql(i,3)
      t4 = qr(i,4)-ql(i,4)
c
c      pressure and enthalpy
c
      t16     = 1.e0/qr(i,1)
      t5      = qr(i,5)
      qr(i,5) = x1*qr(i,5)*t16+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                              +qr(i,4)*qr(i,4)) 
c
      t15     = 1.e0/ql(i,1)
      t19     = ql(i,5) 
      ql(i,5) = x1*ql(i,5)*t15+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                              +ql(i,4)*ql(i,4)) 
c
c      unsplit contributions  f(r)+f(l)
c
      t18 = ax(i)*qr(i,2)+ay(i)*qr(i,3)+az(i)*qr(i,4)+at(i) 
      t17 = ax(i)*ql(i,2)+ay(i)*ql(i,3)+az(i)*ql(i,4)+at(i) 
      t6  = t18*qr(i,1)
      t7  = t17*ql(i,1)
      f1  = t6+t7
      f2  = t6*qr(i,2)+t7*ql(i,2)
      f3  = t6*qr(i,3)+t7*ql(i,3)
      f4  = t6*qr(i,4)+t7*ql(i,4)
      f5  = t6*qr(i,5)+t7*ql(i,5)
      t8  = t5+t19 
      f2  = f2+ax(i)*t8
      f3  = f3+ay(i)*t8
      f4  = f4+az(i)*t8
      f5  = f5-at(i)*t8
c
c      roe averaged variables
c
      t6 = qr(i,1)*t15
      t7 = sqrt(t6) 
      t6 = 1.e0/(1.e0+t7) 
      t8 = t7*t6
c
c     average density
c
      qr(i,1) = ql(i,1)*t7
c
c      u,v,w,h average
c
      t9  = ql(i,2)*t6+qr(i,2)*t8
      t10 = ql(i,3)*t6+qr(i,3)*t8 
      t11 = ql(i,4)*t6+qr(i,4)*t8 
      t12 = ql(i,5)*t6+qr(i,5)*t8 
c
c      extract sound speed
c
      t6  = (t9*t9+t10*t10+t11*t11)*0.5e0
      t7  = gm1*(t12-t6) 
      t8  = sqrt(t7) 
c
      t13 = t9*ax(i)+t10*ay(i)+t11*az(i)
c
c      the variables in t are as follows
c      1-4        delta q1-q4
c      6          q2a
c      7          c2a
c      8          ca
c      9,10,11,12 ua,va,wa,ha
c      13         ubara
c      14,15,16   alpha(1),alpha(4),alpha(5)
c      18,19,17   lambda(1),lambda(4),lambda(5)
c
c      rhoa*delta(ubar) , delta(p)/c2a
c
      ql(i,1) = qr(i,1)*(t18-t17) 
      ql(i,2) = (t5-t19)/t7 
c
      t18 = t13+at(i) 
      t18 = abs(t18)+eps4 
      t19 = t13+at(i)+t8
      t19 = abs(t19)+eps4 
      t17 = t13+at(i)-t8
      t17 = abs(t17)+eps4 
c
      t14 = t18*(t1-ql(i,2))
      t15 = .5e0*(ql(i,2)+ql(i,1)/t8) 
      t16 = (ql(i,2)-t15)*t17 
      t15 = t15*t19 
c
      qr(i,2) = t18*(t2*qr(i,1)-ax(i)*ql(i,1))
      qr(i,3) = t18*(t3*qr(i,1)-ay(i)*ql(i,1))
      qr(i,4) = t18*(t4*qr(i,1)-az(i)*ql(i,1))
      qr(i,5) = t9*qr(i,2)+t10*qr(i,3)+t11*qr(i,4)
c
      ql(i,1) = t14+t15+t16 
      ql(i,2) = t8*(t15-t16)
c
      f1 = f1-ql(i,1) 
      f2 = f2-ql(i,1)*t9-ax(i)*ql(i,2)-qr(i,2)
      f3 = f3-ql(i,1)*t10-ay(i)*ql(i,2)-qr(i,3) 
      f4 = f4-ql(i,1)*t11-az(i)*ql(i,2)-qr(i,4) 
      f5 = f5-ql(i,1)*t12-t13*ql(i,2)-qr(i,5)+t7*c1*t14
c
c      include factor one-half area
c
      t7     = .5e0*area(i) 
      f(i,1) = t7*f1
      f(i,2) = t7*f2
      f(i,3) = t7*f3
      f(i,4) = t7*f4
      f(i,5) = t7*f5
c********************************************************************
	else
c
c     中心格式
c
	if(iskew.eq.1)then
c
c      pressure and enthalpy
c
      t16     = 1.e0/qr(i,1)
      t5      = qr(i,5)
      qr(i,5) = x1*qr(i,5)*t16+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                              +qr(i,4)*qr(i,4)) 
c
      t15     = 1.e0/ql(i,1)
      t19     = ql(i,5) 
      ql(i,5) = x1*ql(i,5)*t15+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                              +ql(i,4)*ql(i,4)) 
c
c      skew-symmetric
c
      t18 = ax(i)*qr(i,2)+ay(i)*qr(i,3)+az(i)*qr(i,4)+at(i) 
      t17 = ax(i)*ql(i,2)+ay(i)*ql(i,3)+az(i)*ql(i,4)+at(i) 
      t6  = t18*qr(i,1)
      t7  = t17*ql(i,1)
      f1  = t6+t7
c      f2  = t6*qr(i,2)+t7*ql(i,2)
c      f3  = t6*qr(i,3)+t7*ql(i,3)
c      f4  = t6*qr(i,4)+t7*ql(i,4)
c      f5  = t6*qr(i,5)+t7*ql(i,5)
	f2  = 0.5*(qr(i,2)+ql(i,2))
	f3  = 0.5*(qr(i,3)+ql(i,3))
	f4  = 0.5*(qr(i,4)+ql(i,4))
      f5  = t6*qr(i,5)+t7*ql(i,5)
c     f2,f3,f4 : 对流通量项 
c	t8       : 压力
      t8  = t5+t19 
      f2  = f1*f2+ax(i)*t8
      f3  = f1*f3+ay(i)*t8
      f4  = f1*f4+az(i)*t8
      f5  = f5-at(i)*t8
c
c      include factor one-half area
c
      t7     = .5e0*area(i) 
      f(i,1) = t7*f1
      f(i,2) = t7*f2
      f(i,3) = t7*f3
      f(i,4) = t7*f4
      f(i,5) = t7*f5
      else
c
c     非斜对称形式
c
c
c      pressure and enthalpy
c
      t16     = 1.e0/qr(i,1)
      t5      = qr(i,5)
      qr(i,5) = x1*qr(i,5)*t16+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                              +qr(i,4)*qr(i,4)) 
c
      t15     = 1.e0/ql(i,1)
      t19     = ql(i,5) 
      ql(i,5) = x1*ql(i,5)*t15+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                              +ql(i,4)*ql(i,4)) 
c
c      unsplit contributions  F(qr)+F(ql)
c
      t18 = ax(i)*qr(i,2)+ay(i)*qr(i,3)+az(i)*qr(i,4)+at(i) 
      t17 = ax(i)*ql(i,2)+ay(i)*ql(i,3)+az(i)*ql(i,4)+at(i) 
      t6  = t18*qr(i,1)
      t7  = t17*ql(i,1)
      f1  = t6+t7
      f2  = t6*qr(i,2)+t7*ql(i,2)
      f3  = t6*qr(i,3)+t7*ql(i,3)
      f4  = t6*qr(i,4)+t7*ql(i,4)
      f5  = t6*qr(i,5)+t7*ql(i,5)
      t8  = t5+t19 
      f2  = f2+ax(i)*t8
      f3  = f3+ay(i)*t8
      f4  = f4+az(i)*t8
      f5  = f5-at(i)*t8
c
c      include factor one-half area
c
      t7     = .5e0*area(i) 
      f(i,1) = t7*f1
      f(i,2) = t7*f2
      f(i,3) = t7*f3
      f(i,4) = t7*f4
      f(i,5) = t7*f5
      endif
c
	endif
 2000 continue
      return
      end
