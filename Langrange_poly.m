function f=Langrange_poly(x,y)
%输入两个等‘个数’向量组，返回由其生成的拉格朗日插值多项式。
%该函数构建一个函数句柄
A=Newton_difference_quo(x,y);
g=@(t)1; %重心
f=@(t)y(1);
for i=1:length(y)-1
    g=@(t)g(t).*(t-x(i));
    f=@(t)f(t)+A(1,i+1).*g(t);
end