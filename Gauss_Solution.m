function [t,y]=Gauss_Solution(s,f,Jf,h,v,y0)
%f为函数句柄，其自变量为k维向量;Jf为f的jacobi阵;h为所需步长;v为时间终点;s为所需RK方法的长度
%本函数给出一个数值解
m=length(y0);
if s==1
    b=1;
    A=1/2;
    A0=1;
elseif s==2
    b=[1/2,1/2]';
    A=[1/4 1/4-sqrt(3)/6
        1/4+sqrt(3)/6 1/4];
    A0=[0 0
        1 0];
elseif s==3
    b=[5/18,4/9,5/18]';
    A=[5/36 2/9-sqrt(15)/15 5/36-sqrt(15)/30
        5/36+sqrt(15)/24 2/9 5/36-sqrt(15)/24
        5/36+sqrt(15)/30 2/9+sqrt(15)/15 5/36];
    A0=[0 0 0
           1/2 0 0
           -1 2 0];
end
n=round(v/h);
y=zeros(m,n);
t=h:h:v;
s=size(A,1);
%y00=y0;
F=zeros(s*m,1);
dF=zeros(s*m,s*m);
for i=1:n
    %%显式方法得初值Y0(6s*1矩阵)
    err=1;
    Y0=zeros(m*s,1);
    Y0(1:m)=y0;
    Y=f(y0);
    for k=2:s
        Y0(m*k-m+1:m*k)=y0+h*kron(A0(k,1:(k-1)),eye(m))*Y;
        Y=[Y;f(Y0(m*k-m+1:m*k))];
    end
    
   while err>=10^-13
       for k=1:s
           F(m*k-m+1:m*k)=f(Y0(m*k-m+1:m*k));
       end
       r=Y0-kron(ones(s,1),y0)-h*kron(A,eye(m))*F;
       for k=1:s
           dF(m*k-m+1:m*k,m*k-m+1:m*k)=Jf(Y0(m*k-m+1:m*k));
       end
       Y1=Y0-(eye(s*m)-h*kron(A,eye(m))*dF)\r;
       err=min(norm(Y1-Y0),norm(r));
       Y0=Y1;
   end
   for k=1:s
       F(m*k-m+1:m*k)=f(Y0(m*k-m+1:m*k));
   end
   y(:,i)=y0+h*kron(b',eye(m))*F;
   y0=y(:,i);
end