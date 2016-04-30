! The semiparametric estimator:
! Z does not appear. Adopt a discrete p(X). 
! we calculate the variance of the estimation as well.   
!     number of simulations is simus
!     number of observations is n
!     number of discrete points is m
!     number of quadrature points is len 
!     length of beta is lb
module commondata
  implicit none
  save
  integer,parameter :: simus=1000,n=1000,m=15,len=30,iseed=1
  integer,parameter :: lb=3,lw=lb*(lb*3+14)
  double precision :: sd,sig,s2
  double precision,dimension(lb) :: theta0,theta
  double precision,dimension(n) :: w
  integer,dimension(n) :: y
  double precision,dimension(m) :: x,p
  double precision,dimension(len) :: xx,ww
  double precision,dimension(len,m) :: pW
end module commondata

program woz
  use commondata
  integer :: ierr,i,j,k,iter,eflag,nbad
  double precision,dimension(lb) ::  beta,vfd,Beta2,Beta3,vfd2,beta1,vfd3
  double precision,dimension(2) :: msx,endpts
  double precision :: d,tmp,yy,tol,tmp1,tmp2,tmp3,tmp4,tmp5
  double precision,dimension(simus,lb) :: LBeta,Lvfd
  double precision,dimension(simus,lb,lb) :: Lvar
  double precision,dimension(lb,lb) :: vjac,Bn,An,Cn,Evar,An1
  double precision,dimension(m,m) :: A
  double precision,dimension(len) :: bb
  double precision,dimension(lw) :: wv
  double precision,dimension(m,lb) :: b
  external fd
  tmp=rand(iseed)
  tol=1e-4
!     true beta
!  beta=[-1.0,1.0,1.0]
  beta=0
  print*,'true beta=',beta
!     mean and standard deviation of x is msx
  msx=[-1.0,1.0]
!     standard deviation of w is s2
  s2=sqrt(0.8)
!     this is points and weights we take for quadrature	
  call gaussq(4,len,0.d0,0.d0,0,endpts,bb,xx,ww)
  !     adjust the xx and ww
  xx=sqrt(2.0)*s2*xx-1.0
  ww=sqrt(2.0)*s2*ww
  d=msx(2)*6./(m-1.0)
  tmp=0
  do i=1,m
     x(i)=msx(1)-msx(2)*3+(i-1)*d
     p(i)=exp(-(x(i)-msx(1))**2/2.0/msx(2)**2) 
     !         p(i)=1
     tmp=tmp+p(i)
  end do
  do i=1,m
     p(i)=p(i)/tmp
  end do
!      print*,'p=',p
  do i=1,len
     do j=1,m
        pW(i,j)=exp(-((xx(i)-x(j))/s2)**2/2)
     end do
  end do
!     
  iter=1
1000 if (iter.le.simus) then
     call gendata(n,lb,beta,msx,s2,w,y)
     beta1=beta
!     print*,'w=',w(1),w(2),w(n-1),w(n),y(1),y(2),y(n-1),y(n)
     call hybrd10(fd,lb,beta1,vfd,tol,eflag,wv,lw)
!     print*,'after hybrid',beta1
     if (vfd(1).eq.999) then
        eflag=-1
     end if
     if (eflag.ne.1) then
        print*,'iter=',iter,'eflag=',eflag,'failed fd','beta=',beta1
        goto 1000
     end if
     do j=1,lb
        beta2=beta1
        beta3=beta1
        beta2(j)=beta1(j)+0.0001/2
        beta3(j)=Beta1(j)-0.0001/2
        call fd(lb,Beta2,vfd2,eflag) 
        call fd(lb,Beta3,vfd3,eflag)
        do i=1,lb
           vjac(j,i)=(vfd3(i)-vfd2(i))/0.0001
        end do
     end do
     do i=1,lb
        Lvfd(iter,i)=vfd(i)
     end do
     do i=1,lb
        LBeta(iter,i)=beta1(i)
        do j=1,lb
           An(i,j)=-vjac(i,j)/n
        end do
     end do
     An1=An
     call fd1(m,lb,n,Beta1,w,y,p,x,s2,len,xx,ww,pW,Bn)
     call inv(lb,An1,An)
     call mul3(An,Bn,Evar,n) 
!            print*, 'evar=',Evar
!     The above two computes Evar=inv(An) Bn inv(An')/n
     Lvar(iter,:,:)=Evar
     iter=iter+1
     goto 1000
  end if
  open(1, file='noropt.dat')
  do iter=1,simus
     write(1,*) (LBeta(iter,j),j=1,lb)
  end do
  do iter=1,simus
     do i=1,lb
        write(1,*) (Lvar(iter,i, j), j=1,lb)
     end do
  end do
  close(1)
  return
end program woz
!
subroutine gendata(n,lb,beta,msx,s2,W,Y)
  integer,intent(in) :: n,lb
  double precision,dimension(lb),intent(in) :: beta
  double precision,dimension(2),intent(in) :: msx
  double precision,intent(in) :: s2
  double precision,dimension(n),intent(out) :: W
  integer,dimension(n),intent(out) :: Y
  double precision,dimension(n)::  X
  double precision:: tmp1,tmp
  real:: u1,u2
  integer i
!     Generate X's from Normal mean msx(1) variance msx(2)^2, They are
!     unobservable though. 
  do i=1,n
     call rnorm(u1,u2)
     X(i)=u1*msx(2)+msx(1)
!     Generate W's from Normal mean X variance sw^2. 
     W(i)=u2*s2+X(i)
!     Generate Y's from Quadratic Logistic. 
     tmp1=1./(1.+exp(beta(1)+beta(2)*X(i)+beta(3)*X(i)**2))
     tmp=rand(0)
     Y(i)=0
     if (tmp .gt. tmp1) then
        Y(i)=1
     end if
  end do
  return
end subroutine gendata
!
subroutine fd(lb,beta,vfd,eflag)
  use commondata,only : n,m,len,w,y,p,x,s2,pW,xx,ww
  integer,intent(in) :: lb
  integer,intent(out) :: eflag
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(lb),intent(out):: vfd
  integer:: i,j,k,ierr
  integer,dimension(m):: ipvt
  double precision,dimension(m,m):: A
  double precision,dimension(m):: wv,alpha01,alpha02,alpha03
  double precision,dimension(m,lb):: alpha0
  double precision,dimension(lb,m):: b
  double precision tmp,cond,zz
  double precision,dimension(2,m):: pY
  double precision,dimension(len,2):: pWY
  double precision,dimension(lb):: tmp1,tmp2,yy1,yy2
  if (sum(abs(beta)).gt.30) then
     vfd(1)=999
     return
  end if
  b=0
  pWY=0
  do j=1,m
     call logis2(0,x(j),beta,pY(1,j))
     call logis2(1,x(j),beta,pY(2,j))
  end do
  do i=1,len
     do j=1,m 
        pWY(i,1)=pWY(i,1)+p(j)*pW(i,j)*pY(1,j)
        pWY(i,2)=pWY(i,2)+p(j)*pW(i,j)*pY(2,j)
     end do
  end do
  do i=1,m
     do j=1,i-1
        call intbwy(m,lb,beta,p,s2,x,i,j,len,xx,ww,pW,pWY,pY,yy1,yy2,zz)
        A(i,j)=p(j)*zz
        A(j,i)=p(i)*zz
        do k=1,lb
           b(k,i)=b(k,i)+p(j)*yy1(k)
           b(k,j)=b(k,j)+p(i)*yy2(k)
        end do
     end do
     call intbwyl(m,lb,beta,p,s2,x,i,len,xx,ww,pW,pWY,pY,yy1,zz)
     A(i,i)=p(i)*zz
     do k=1,lb
        b(k,i)=b(k,i)+p(i)*yy1(k)
     end do
  end do
  do i=1,m
     alpha01(i)=b(1,i)
     alpha02(i)=b(2,i)
     alpha03(i)=b(3,i)
  end do
!  here, cond is condition number, the algorithm breaks down if it's too large
  call decomp(m,m,A,cond,ipvt,wv)
  !      print*, cond
  call solvels(m,m,A,alpha01,ipvt)
  call solvels(m,m,A,alpha02,ipvt)
  call solvels(m,m,A,alpha03,ipvt)
  do i=1,m
     alpha0(i,1)=alpha01(i)
     alpha0(i,2)=alpha02(i)
     alpha0(i,3)=alpha03(i)
  end do
!         print*, 'alpha0=', alpha0
  do i=1,lb
     vfd(i)=0
  end do
  do i=1,n
     call sb(lb,n,m,beta,w(i),y(i),p,x,s2,pY,tmp1)
     call alpha1(lb,n,m,beta,w(i),y(i),p,x,s2,alpha0,pY,tmp2)
!         print*, tmp2
     do j=1,lb
        vfd(j)=vfd(j)+tmp1(j)-tmp2(j)
     end do
  end do
  tmp=0
  do i=1,lb
     tmp=tmp+vfd(i)
  end do
  if (.not.((tmp.gt.0).or.(tmp.le.0))) then
     vfd(1)=999 
  end if
!  print*,'beta=',beta
!  print*,'vfd=',vfd
  return
end subroutine fd
      
subroutine sb(lb,n,m,beta,w,y,p,x,s2,pY,yy)
  integer,intent(in):: lb,n,m,y
  double precision,dimension(lb),intent(in) :: beta
  double precision,intent(in) :: w,s2
  double precision,dimension(m),intent(in) :: p,x
  double precision,dimension(2,m),intent(in) :: pY
  double precision,dimension(lb),intent(out) :: yy
  double precision,dimension(m) :: tmp
  double precision,dimension(lb,m) :: tmpp
  double precision,dimension(lb) :: tp0
  double precision :: t0,Ttmp
  integer:: i,j
  Ttmp=0
  yy=0
  do i=1,m
     tmp(i)=p(i)*pY(y+1,i)*exp(-((w-x(i))/s2)**2/2.0)
     Ttmp=Ttmp+tmp(i)
     call sbetaf(y,x(i),beta,tp0)
     do j=1,lb 
        yy(j)=yy(j)+tp0(j)*tmp(i)
     end do
  end do
  do j=1,lb
     yy(j)=yy(j)/Ttmp 
  end do
  return
end subroutine sb
!     
subroutine alpha1(lb,n,m,beta,w,y,p,x,s2,alpha0,pY,yy)
  integer,intent(in):: n,m,lb,y
  double precision,dimension(m),intent(in):: p,x
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(m,lb),intent(in):: alpha0
  double precision,dimension(2,m),intent(in):: pY
  double precision,intent(in):: w,s2
  double precision,dimension(lb),intent(out) :: yy
  double precision Ttmp
  double precision,dimension(m) :: tmp
  integer:: i,j
  Ttmp=0
  do i=1,m
     tmp(i)=p(i)*pY(y+1,i)*exp(-((w-x(i))/s2)**2/2.0)
     Ttmp=Ttmp+tmp(i)
  end do
  do i=1,lb
     yy(i)=0
     do j=1,m
        yy(i)=yy(i)+alpha0(j,i)*tmp(j)
     end do
     yy(i)=yy(i)/Ttmp
  end do
  return
end subroutine alpha1
!
subroutine sbetaf(y,x,beta,yy)
  use commondata,only : lb
  integer,intent(in):: y
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(lb),intent(out):: yy
  double precision,intent(in) :: x
  double precision tmp,a
  tmp=exp(beta(1)+beta(2)*x+beta(3)*x**2)
  if (y.eq.0) then
!         a=-tmp/(1+tmp)**2
     a=-tmp/(1+tmp)
  else
!         a=tmp/(1+tmp)**2
     a=1/(1+tmp)
  end if
  yy(1)=a
  yy(2)=a*x
  yy(3)=yy(2)*x
  return
end subroutine sbetaf
!
subroutine logis2(y,x,beta,yy)
  use commondata,only : lb
  integer,intent(in):: y
  double precision,dimension(lb),intent(in):: beta
  double precision,intent(in):: x
  double precision,intent(out):: yy
  if (y.eq.0) then
     yy=1/(1+exp(beta(1)+beta(2)*x+beta(3)*x**2))
  else
     yy=1-1/(1+exp(beta(1)+beta(2)*x+beta(3)*x**2)) 
  end if
  return
end subroutine logis2
!      
subroutine intbwy(m,lb,beta,p,s2,x,i,j,len,xx,ww,pW,pWY,pY,yy1,yy2,zz)
  integer,intent(in):: i,j,m,lb
  integer::ii
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(m),intent(in):: p,x
  double precision,intent(in):: s2
  double precision,dimension(len,2),intent(in):: pWY
  double precision,dimension(2,m),intent(in):: pY
  double precision,dimension(len,m),intent(in):: pW
  double precision,intent(out):: zz
  double precision,dimension(lb),intent(out):: yy1,yy2
  double precision,dimension(lb):: tmp0,tmp1,tmp2,tmp3
  double precision :: t1,a0,a1
  call sbetaf(0,x(j),beta,tmp0)
  call sbetaf(1,x(j),beta,tmp1)
  call sbetaf(0,x(i),beta,tmp2)
  call sbetaf(1,x(i),beta,tmp3)
  a0=pY(1,i)*pY(1,j)
  a1=pY(2,i)*pY(2,j)
  if (a0.ne.0) then
     call intw(m,lb,beta,s2,x,i,j,0,len,xx,ww,pW,pWY,t1)
     a0=t1*a0
  end if
  if (a1.ne.0) then
     call intw(m,lb,beta,s2,x,i,j,1,len,xx,ww,pW,pWY,t1)
     a1=t1*a1
  end if
  do ii=1,lb
     yy1(ii)=a0*tmp0(ii)+a1*tmp1(ii)
     yy2(ii)=a0*tmp2(ii)+a1*tmp3(ii)
  end do
  zz=a0+a1
  return
end subroutine intbwy
!      
subroutine intbwyl(m,lb,beta,p,s2,x,i,len,xx,ww,pW,pWY,pY,yy,zz)
  integer,intent(in):: i,m,lb
  integer:: ii
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(m),intent(in):: p,x
  double precision,intent(in):: s2
  double precision,dimension(len,2),intent(in):: pWY
  double precision,dimension(2,m),intent(in):: pY
  double precision,dimension(len,m),intent(in):: pW
  double precision,intent(out):: zz
  double precision,dimension(lb),intent(out):: yy
  double precision,dimension(lb):: tmp0,tmp1
  double precision :: t1,a0,a1
  call sbetaf(0,x(i),beta,tmp0)
  call sbetaf(1,x(i),beta,tmp1)
  a0=pY(1,i)*pY(1,i)
  a1=pY(2,i)*pY(2,i)
  if (a0.ne.0) then
     call intw(m,lb,beta,s2,x,i,i,0,len,xx,ww,pW,pWY,t1)
     a0=t1*a0
  end if
  if (a1.ne.0) then
     call intw(m,lb,beta,s2,x,i,i,1,len,xx,ww,pW,pWY,t1)
     a1=t1*a1
  end if
  do ii=1,lb
     yy(ii)=a0*tmp0(ii)+a1*tmp1(ii)
  end do
  zz=a0+a1
  return
end subroutine intbwyl
!
subroutine intw(m,lb,beta,s2,x,i,j,k,len,xx,ww,pW,pWY,yy)
  integer,intent(in):: i,j,k,len,m,lb
  integer:: ii
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(len),intent(in):: xx,ww
  double precision,intent(in):: s2,x
  double precision,dimension(len,2),intent(in):: pWY
  double precision,dimension(len,m),intent(in):: pW
  double precision,intent(out):: yy
  double precision :: t1,a0,a1,tmp
  yy=0
  do ii=1,len
     call integ(lb,len,m,ii,beta,s2,x,i,j,k,pW,pWY,tmp)
     yy=yy+ww(ii)*tmp*exp((xx(ii)+1)**2/s2**2/2)
  end do
  return
end subroutine intw
!  
subroutine integ(lb,len,m,ii,beta,s2,x,i,j,k,pW,pWY,yy)
  integer,intent(in):: i,j,k,len,m,lb,ii
  double precision,dimension(lb),intent(in):: beta
  double precision,intent(in):: s2,x
  double precision,dimension(len,2),intent(in):: pWY
  double precision,dimension(len,m),intent(in):: pW
  double precision,intent(out):: yy
  yy=pW(ii,i)*pW(ii,j)
  yy=yy/pWY(ii,k+1)
  return
end subroutine integ
!                
subroutine fd1(m,lb,n,beta,w,y,p,x,s2,len,xx,ww,pW,yy)
  integer,intent(in) :: m,lb,n,len
  double precision,dimension(lb),intent(in):: beta
  double precision,dimension(n),intent(in):: w
  integer,dimension(n),intent(in):: y
  double precision,dimension(m),intent(in):: p,x
  double precision,dimension(len),intent(in):: xx,ww
  double precision,intent(in):: s2
  double precision,dimension(len,m),intent(in):: pW
  double precision,dimension(lb,lb),intent(out):: yy
  integer,dimension(m):: ipvt
  integer:: i,j,k
  double precision,dimension(m,m):: A
  double precision,dimension(lb,m):: b
  double precision,dimension(m,lb):: alpha0
  double precision,dimension(lb):: tmp,tmp1,tmp2,yy1,yy2
  double precision,dimension(m):: wv,alpha01,alpha02,alpha03
  double precision temp,zz,cond
  double precision,dimension(2,m):: pY
  double precision,dimension(len,2):: pWY
  yy=0
  b=0
  do j=1,m
     call logis2(0,x(j),beta,pY(1,j))
     call logis2(1,x(j),beta,pY(2,j))
  end do
  do i=1,len
     pWY(i,1)=0
     pWY(i,2)=0
     do j=1,m
        pWY(i,1)=pWY(i,1)+p(j)*pW(i,j)*pY(1,j)
        pWY(i,2)=pWY(i,2)+p(j)*pW(i,j)*pY(2,j)
     end do
  end do
  do i=1,m
     do j=1,i-1
        call intbwy(m,lb,beta,p,s2,x,i,j,len,xx,ww,pW,pWY,pY,yy1,yy2,zz)
        A(i,j)=p(j)*zz
        A(j,i)=p(i)*zz
        do k=1,lb
           b(k,i)=b(k,i)+p(j)*yy1(k)
           b(k,j)=b(k,j)+p(i)*yy2(k) 
        end do
     end do
     call intbwyl(m,lb,beta,p,s2,x,i,len,xx,ww,pW,pWY,pY,yy1,zz)
     A(i,i)=p(i)*zz
     do k=1,lb
        b(k,i)=b(k,i)+p(i)*yy1(k)
     end do
  end do
  do i=1,m
     alpha01(i)=b(1,i)
     alpha02(i)=b(2,i)
     alpha03(i)=b(3,i)
  end do
  call decomp(m,m,A,cond,ipvt,wv)
  call solvels(m,m,A,alpha01,ipvt)
  call solvels(m,m,A,alpha02,ipvt)
  call solvels(m,m,A,alpha03,ipvt)
  do i=1,m
     alpha0(i,1)=alpha01(i)
     alpha0(i,2)=alpha02(i)
     alpha0(i,3)=alpha03(i)
  end do
  yy=0
  do i=1,n
     call sb(lb,n,m,beta,w(i),y(i),p,x,s2,pY,tmp1)
     call alpha1(lb,n,m,beta,w(i),y(i),p,x,s2,alpha0,pY,tmp2)
     do j=1,lb
        tmp(j)=tmp1(j)-tmp2(j)
     end do
     do j=1,lb
        do k=1,lb
           yy(j,k)=yy(j,k)+tmp(j)*tmp(k)/n
        end do
     end do
  end do
  return
end subroutine fd1
! invert a matrix
subroutine inv(lb,A,A1)
  integer,intent(in) :: lb
  double precision,dimension(lb,lb),intent(in) :: A
  double precision,dimension(lb,lb),intent(out) :: A1
  integer :: i
  integer,dimension(lb) :: ipvt
  double precision,dimension(lb,lb) :: B
  double precision :: cond
  double precision,dimension(lb) :: wv
  A1=0
  do i=1,lb
     A1(i,i)=1
  end do
  call decomp(lb,lb,A,cond,ipvt,wv)
  do i=1,lb
     call solvels(lb,lb,A,A1(:,i),ipvt)
  end do
  return
end subroutine inv
!     calculates T=C*B*C'/n
subroutine mul3(C,B,T,n)
  integer,intent(in):: n
  integer:: i,j,k
  double precision,dimension(3,3),intent(in):: B,C
  double precision,dimension(3,3):: A
  double precision,dimension(3,3),intent(out):: T
  do i=1,3
     do j=1,3
        A(i,j)=0
        do k=1,3
           A(i,j)=A(i,j)+C(i,k)*B(k,j)
        end do
     end do
  end do
  do i=1,3
     do j=1,3
        T(i,j)=0
        do k=1,3
           T(i,j)=T(i,j)+A(i,k)*C(j,k)
        end do
        T(i,j)=T(i,j)/n
     end do
  end do
  return
end subroutine mul3
