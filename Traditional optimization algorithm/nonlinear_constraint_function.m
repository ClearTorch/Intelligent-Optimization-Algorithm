function [C,Ceq]=nonlinear_constraint_function(x)
C=[-x(1)^2+x(2)-x(3)^2;x(1)+x(2)^2+x(3)^3-20;]; %非线性约束不等式( <=0)
Ceq=[-x(1)-x(2)^2+2;x(2)+2*x(3)^2-3];  %非线性约束等式（ =0）