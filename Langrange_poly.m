function f=Langrange_poly(x,y)
%���������ȡ������������飬�����������ɵ��������ղ�ֵ����ʽ��
%�ú�������һ���������
A=Newton_difference_quo(x,y);
g=@(t)1; %����
f=@(t)y(1);
for i=1:length(y)-1
    g=@(t)g(t).*(t-x(i));
    f=@(t)f(t)+A(1,i+1).*g(t);
end