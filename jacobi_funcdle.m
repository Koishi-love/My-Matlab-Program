function Jf=jacobi_funcdle(f,k)
%%�ú������ڣ���kά������Ϊ�����ĺ������f����f���ſ˱Ⱦ�����
X=[];
for j=1:k
    syms (['x',num2str(j)]);
    X=[X,eval(['x',num2str(j)])];
end
Jf=jacobian(f(X'),X);
Jf=matlabFunction(Jf,'vars',{X});
Jf=@(x)Jf(x');