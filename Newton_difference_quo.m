function A=Newton_difference_quo(x,y)
%输入两个等长向量组，返回一棵牛顿差商树
%若设x长度为s，则A为s*s左上三角矩阵，每列元素自上而下依次为其左列元素自上而下依次生成的差商，其余元素用0填充
s=length(x);
A=zeros(s);
A(:,1)=y;
for i=2:s
    for j=1:s-i+1  %层叶内深度
        A(j,i)=(A(j,i-1)-A(j+1,i-1))/(x(j)-x(j+i-1));
    end
end