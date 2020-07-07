%设n个节点的数据以数组x0（已知点的横坐标）， y0（函数值）， y1（导数值）
% 输入（注意Matlat 的数组下标从1 开始），m 个插值点以数组x 输入，输出数组y 为m
% 个插值。
% 如果对插值函数，不仅要求它在节点处与函数同值，而且要求它与函数有相同的一
% 阶、二阶甚至更高阶的导数值，这就是Hermite 插值问题

function y=hermite(x0,y0,y1,x)
n=length(x0);m=length(x);
for k=1:m
    yy=0.0;
    for i=1:n
        h=1.0;
        a=0.0;
        for j=1:n
            if j~=i
                h=h*((x(k)-x0(j))/(x0(i)-x0(j)))^2;
                a=1/(x0(i)-x0(j))+a;
            end
        end
        yy=yy+h*((x0(i)-x(k))*(2*a*y0(i)-y1(i))+y0(i));
    end
    y(k)=yy;
end