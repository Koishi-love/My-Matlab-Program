function A=Newton_difference_quo(x,y)
%���������ȳ������飬����һ��ţ�ٲ�����
%����x����Ϊs����AΪs*s�������Ǿ���ÿ��Ԫ�����϶�������Ϊ������Ԫ�����϶����������ɵĲ��̣�����Ԫ����0���
s=length(x);
A=zeros(s);
A(:,1)=y;
for i=2:s
    for j=1:s-i+1  %��Ҷ�����
        A(j,i)=(A(j,i-1)-A(j+1,i-1))/(x(j)-x(j+i-1));
    end
end