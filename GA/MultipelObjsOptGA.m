clear
clc
fitnessfcn= @my_first_multi; %适应度函数句柄
nvars= 2;%变量个数
lb=[-5,-5];%下限
ub=[5,5];%上限
A=[]; b=[];%非线性不等式约束
Aeq=[]; beq=[];%非线性等式约束

options=gaoptimset( 'ParetoFraction', 0.3, 'PopulationSize', 100,'Generations', 200,...
'StallGenLimit',200, 'TolFun' ,1e-100,'PlotFcns',@gaplotpareto) ; 
[x,fval] = gamultiobj(fitnessfcn,nvars, A,b, Aeq,beq, lb,ub,options) ;


function f= my_first_multi(x)
f(1) = x(1)^4-10* x(1)^2 +x(1) * x(2)+ x(2)^4- (x(1)^2)*(x(2)^2);
f(2) = x(2)^4- (x(1)^2)*(x(2)^2) + x(1)^4+ x(1) * x(2);
end
