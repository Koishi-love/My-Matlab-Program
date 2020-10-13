function Jf=jacobi_funcdle(f,k)
%%该函数用于：以k维列向量为变量的函数句柄f，求f的雅克比矩阵句柄
X=[];
for j=1:k
    syms (['x',num2str(j)]);
    X=[X,eval(['x',num2str(j)])];
end
Jf=jacobian(f(X'),X);
Jf=matlabFunction(Jf,'vars',{X});
Jf=@(x)Jf(x');