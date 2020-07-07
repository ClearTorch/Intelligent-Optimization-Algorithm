%% 无约束最优化问题
% 函数调用  [x,fval,key,c]=fminsearch(@unlimitedproblem_objectiive_fun,x0,options)
clear
x0=[0];
[x,fval,key,c]=fminsearch(@fun,x0);

%% 线性规划
%函数调用 [x,fval]linprog(F,A,B,Aeq,Beq,LB,UB,options) 
%F为目标函数系数矩阵的列向量，A为小于等于约束条件标准型的系数矩阵，B为小于等于约束条件的常数矩阵
%Aeq为等式条件系数矩阵，Beq为等式条件常数项矩阵
%LB,UB为决策变量下限和上限矩阵
%options算法选择（默认为单纯形法），单纯形、内点法、内点迭代法
%对偶理论（对偶问题时原问题的转置）和灵敏度分析(系数变化对最优解的影响程度)
clear all;
F=[-2 -3 5]';A=[-2 5 -1;1 3 1];B=[-10;12];
Aeq=[1 1 1];Beq=[7];LB=[0;0;0];UB=[inf;inf;inf];options= optimoptions('linprog','Algorithm','interior-point');
[x,fval,exitflag,output,lambda]=linprog(F,A,B,Aeq,Beq,LB,UB,options);

%% 整数线性规划
% 解都是整数的线性规划
%求解方法：求解方法分类：
% （i）分枝定界法—可求纯或混合整数线性规划。
% （ii）割平面法—可求纯或混合整数线性规划。
% （iii）隐枚举法—求解“0-1”整数规划：
% ①过滤隐枚举法；
% ②分枝隐枚举法。
% （iv）匈牙利法—解决指派问题（“0-1”规划特殊情形）。
% （v）蒙特卡洛法—求解各种类型规划。

%分枝界定法（每个决策变量为一个分支，进行线性规划求解）
clear all; %修改取整函数
F=[-40 -90]';A=[9 7;7 20];B=[56;70];
Aeq=[];Beq=[];LB0=[0;0];UB0=[inf;inf];options= optimoptions('linprog','Algorithm','interior-point');
[x0,fval0,exitflag0,output0,lambda0]=linprog(F,A,B,Aeq,Beq,LB0,UB0,options);

% 分支问题B1 x0（1）>ceil(x0(1))
LB1=[ceil(x0(1));0];UB1=[inf;inf];
[x1,fval1,exitflag1,output1,lambda1]=linprog(F,A,B,Aeq,Beq,LB1,UB1,options);

%分支问题B2 x0（1）<fix(x0(1))
LB2=[0;0];UB2=[floor(x0(1));inf];
[x2,fval2,exitflag2,output2,lambda2]=linprog(F,A,B,Aeq,Beq,LB2,UB2,options);

%分支B11 X1(2)>ceil（x1(2)）
LB11=[ceil(x0(1));ceil(x1(2))];UB11=[inf;inf];
[x11,fval11,exitflag11,output11,lambda11]=linprog(F,A,B,Aeq,Beq,LB11,UB11,options);
%分支B12 x1（2）<fix（x1（2））
LB12=[0;0];UB12=[floor(x0(1));floor(x1(2))];
[x12,fval12,exitflag12,output12,lambda12]=linprog(F,A,B,Aeq,Beq,LB12,UB12,options);
% 分支B21 X1(2)>round（x1(2)）
LB21=[0;ceil(x1(2))];UB21=[floor(x1(1));inf];
[x21,fval21,exitflag21,output21,lambda21]=linprog(F,A,B,Aeq,Beq,LB21,UB21,options);
% 分支B22 x0（2）<fix（x1（2））
LB22=[0;0];UB22=[floor(x1(1));ceil(x1(2))];
[x22,fval22,exitflag22,output22,lambda22]=linprog(F,A,B,Aeq,Beq,LB22,UB22,options);

%% 非线性规划
%非线性目标函数用.m文件编写，非线性不等式表达式用函数编写C，非线性等式条件Ceq

clear all;
A=[];B=[];Aeq=[];Beq=[];LB=[0;0;0];UB=[];
x0=[1;1;1;1];
[x,y]=fmincon('nonliner_objective_function',x0,A,B,Aeq,Beq,LB,UB,'nonlinear_constraint_function');

% 二次规划问题
% 函数调用 x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)
% 化为标准型  min  1/2*X'*H*X+f'*X
H = [1 -1; -1 2]; 
f = [-2; -6];
A = [1 1; -1 2; 2 1];
B = [2; 2; 3];
Aeq=[];
Beq=[];
LB = zeros(2,1);
UB=[];
[x,fval,exiflag,output,lambda]=quadprog(H,f,A,B,Aeq,Beq,LB,UB)

%% 插值与拟合
% 拉格朗日插值
y=lagrange(x0,y0,x);
%Hermite插值
y=lagrange(x0,y0,x);
% 分段线性插值
%'nearest' 最近项插值
% 'linear' 线性插值
% 'spline' 逐段3 次样条插值
% 'cubic' 保凹凸性3 次插值。
method='spline';
y=interp1(x0,y0,x,'method');

% 三次样条插值
y=interp1(x0,y0,x,'spline');
y=spline(x0,y0,x);
pp=csape(x0,y0,conds),y=ppval(pp,x);

%% 方差分析
% 1单因素方差分析p=anoval(x)
% 返回值p 是一个概率，当p >α 时接受0 H ，x 为m× r 的数据矩阵，x 的每一列是一个
% 水平的数据（这里各个水平上的样本容量n m i = ）。另外，还输出一个方差表和一个
% Box 图。
x=[256 254 250 248 236
242 330 277 280 252
280 290 230 305 220
298 295 302 289 252];
p=anova1(x)

% 非均衡数据
% 处理非均衡数据的用法为：
% p=anova1(x,group)
% x 为向量，从第1 组到第r 组数据依次排列；group 为与x 同长度的向量，标志x 中数
% 据的组别（在与x 第i组数据相对应的位置处输入整数i(i = 1,2,L, r)）。
x=[1620 1580 1460 1500
1670 1600 1540 1550
1700 1640 1620 1610
1750 1720 1680 1800];
x=[x(1:4),x(16),x(5:8),x(9:11),x(12:15)];
g=[ones(1,5),2*ones(1,4),3*ones(1,3),4*ones(1,4)];
p=anova1(x,g)

%2 双因素方差分析
% p=anova2(x,reps)
% 其中x 不同列的数据表示单一因素的变化情况，不同行中的数据表示另一因素的变化情
% 况。如果每种行—列对（“单元”）有不止一个的观测值，则用参数reps 来表明每个“单
% 元”多个观测值的不同标号，即reps 给出重复试验的次数t

clc,clear
x0=[58.2,52.6 56.2,41.2 65.3,60.8
    49.1,42.8 54.1,50.5 51.6,48.4
    60.1,58.3 70.9,73.2 39.2,40.7
    75.8,71.5 58.2,51.0 48.7,41.4];
x1=x0(:,1:2:5);x2=x0(:,2:2:6);
for i=1:4
    x(2*i-1,:)=x1(i,:);
    x(2*i,:)=x2(i,:);
end
[p,t,st]=anova2(x,2);

%% 回归分析
% 多元线性回归
% [b,bint,r,rint,stats]=regress(Y,X,alpha);
% ，alpha 为显著性水平（缺省时设定为0.05），b,bint 为回归系数估计值和
% 它们的置信区间，r,rint 为残差（向量）及其置信区间，stats 是用于检验回归模型的统
% 计量，有四个数值，第一个是R2，第二个是F ，第三个是与F 对应的概率p ， p <α 拒绝0 H ，回归模型成立，第四个是残差的方差s
% 残差及其置信区间可以用 rcoplot(r,rint)画图
x=0.1:0.01:0.18;
y=[42,41.5,45.0,45.5,45.0,47.5,49.0,55.0,50.0];
plot(x,y,'+')
clc,clear
x1=[0.1:0.01:0.18]';
y=[42,41.5,45.0,45.5,45.0,47.5,49.0,55.0,50.0]';
x=[ones(9,1),x1];
[b,bint,r,rint,stats]=regress(y,x);
b,bint,stats,rcoplot(r,rint)

% 即? 27.4722
% 0 β = ， ? 137.5000
% 1 β = ， 0
% ?β
% 的置信区间是[18.6851,36.2594]， 1
% ?β
% 的置信区
% 间是[75.7755,199.2245]；R2 = 0.7985，F = 27.7469， p = 0.0012，s2 = 4.0883。

% 二次曲线拟合[s,p]=polyfit(x0,y0,2)
% 上面的s是一个数据结构，用于计算函数值，如
% [y,delta]=polyconf(p,x0,s);y
% 得到y 的拟合值，及预测值y 的置信区间半径delta。
x0=17:2:29;x0=[x0,x0];
y0=[20.48 25.13 26.15 30.0 26.1 20.3 19.35...
24.35 28.11 26.3 31.4 26.92 25.7 21.3];
[p,s]=polyfit(x0,y0,2); p
[y,delta]=polyconf(p,x0,s);y
polytool(x0,y0,2)
% 用polytool(x0,y0,2)，可以得到一个如1图的交互式画面，在画面中绿色曲线为拟合
% 曲线，它两侧的红线是y 的置信区间。你可以用鼠标移动图中的十字线来改变图下方
% 的x 值，也可以在窗口内输入，左边就给出y 的预测值及其置信区间。

% 多元二项式回归
% 多元二项式回归的命令rstool，它也产生一个交互式画面，
% 并输出有关信息，用法是
% rstool(x,y,model,alpha)
% 其中输入数据x,y分别为n ×m矩阵和n维向量，alpha为显著性水平α （缺省时设定为
% 0.05），model由下列4个模型中选择1个（用字符串输入，缺省时设定为线性模型）：
% linear（线性）,purequadratic（纯二次）,interaction（交叉）,quadratic（完全二次）

x1=[120 140 190 130 155 175 125 145 180 150]';
x2=[100 110 90 150 210 150 250 270 300 250]';
y=[102 100 120 77 46 93 26 69 65 85]';
x=[x1 x2];
rstool(x,y,'purequadratic')

% 非线性回归
% 非线性回归是指因变量 y 对回归系数m β , ,β 1 L （而不是自变量）是非线性的。
% Matlab统计工具箱中的nlinfit，nlparci，nlpredci，nlintool，不仅给出拟合的回归系数，
% 而且可以给出它的置信区间，及预测值和置信区间等
% nlintool(x,y,'huaxue',beta)
% 可看到画面，并传出剩余标准差rmse

%逐步回归
% 在Matlab统计工具箱中用作逐步回归的是命令stepwise，它提供了一个交互式画
% 面，通过这个工具你可以自由地选择变量，进行统计分析，其通常用法是：
% stepwise(x,y,inmodel,alpha)
% 其中x是自变量数据，y是因变量数据，分别为n ×m和n ×1矩阵，inmodel是矩阵x的
% 列数的指标，给出初始模型中包括的子集（缺省时设定为空），alpha为显著性水平。
% Stepwise Regression 窗口，显示回归系数及其置信区间，和其它一些统计量的信
% 息。绿色表明在模型中的变量，红色表明从模型中移去的变量。在这个窗口中有Export
% 按钮，点击Export产生一个菜单，表明了要传送给Matlab工作区的参数，它们给出了统
% 计计算的一些结果。
clc,clear
x0=[1 7 26 6 60 78.5
2 1 29 15 52 74.3
3 11 56 8 20 104.3
4 11 31 8 47 87.6
5 7 52 6 33 95.9
6 11 55 9 22 109.2
7 3 71 17 6 102.7
8 1 31 22 44 72.5
9 2 54 18 22 93.1
10 21 47 4 26 115.9
11 1 40 23 34 83.8
12 11 66 9 12 113.3
13 10 68 8 12 109.4];
x=x0(:,2:5);
y=x0(:,6);
stepwise(x,y)
