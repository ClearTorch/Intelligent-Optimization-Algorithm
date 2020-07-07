%% 多种群遗传算法
%求二元函数最大值 f(x,y)=21.5+x*sin(x*pi)+y*sin(20*pi*y)
% 其中x范围[-3.0,12.1]   y的范围[4.1,5.8]

% 注意 ：只需要更改目标函数、自变量维数、译码矩阵就可以使用程序

%% 标准遗算法SGA

clc
clear all
close all
% 定义遗传算法参数
pc=0.7;
pm=0.05;
NIND= 40;   % 为个体数目
MAXGEN = 500;% 最大遗传代数
NVAR= 2;      % 变量的维数
PRECI= 20;  % 变量的二进制位数
GGAP= 0.9; % 代沟
trace= zeros(MAXGEN,1);
FieldD= [rep(PRECI,[1 ,NVAR]);[- 3,4.1;12.1,5.8];rep([1;0;1;1],[1,NVAR])]; % 建立区 城描述器
Chrom = crtbp(NIND, NVAR * PRECI);% 创建初始种群
gen= 0; % 代计数器
maxY= 0; % 最优值
ObjV=ObjectFunction(bs2rv(Chrom,FieldD)); % 计算初始种群个体的目标函数值
while gen<MAXGEN % 迭代
    FitnV= ranking(-ObjV);% 分配适应度值(assign fitness values)
    SelCh= select('sus', Chrom, FitnV, GGAP);% 选择
    SelCh= recombin('xovsp', SelCh,pc);% 重组
    SelCh= mut(SelCh,pm); % 变异
    ObjVSel = ObjectFunction(bs2rv(SelCh, FieldD));%计算子代目标函数值
    [Chrom,ObjV]= reins(Chrom,SelCh,1,1,ObjV,ObjVSel);%重插人
    gen= gen+ 1;% 代计数器增加
    if maxY< max(ObjV)% 命令行由口
        maxY = max(ObjV);
    end
    trace(gen,1)=maxY;
end
% 进化过程图
plot(1:gen,trace(:,1));
% 输出最优解
[Y,I]=max(ObjV);
X=bs2rv(Chrom,FieldD);
disp(['最优值为：',num2str(Y)])
disp(['对应的自变量取值为：',num2str(X(I,:))])

%% 多种群遗传算法
clc
clear all
close all
% 参数设置
NIND=40;                        % 个体数目
NVAR=2;                         % 变量维数
PRECI=20;                       % 变量的二进制位数
GGAP=0.9;                       %
MP= 10;	% 种群数目
FieldD = [rep(PRECI,[1,NVAR]);[- 3,4.1;12.1,5.8];rep([1;0;1;1],[1,NVAR])]; %译码矩阵
for i= 1:MP
    Chrom{i} = crtbp(NIND,NVAR*PRECI) ;	%创建初始种群
end

pc= 0.7+(0.9-0.7)*rand(MP,1);	%在[O.7,0.9]区间内随机产生交叉概率
pm=0.001+(0.05-0.001)*rand(MP,1);	%在[0.001,0. 05]区间内随机产生变异概率，
gen= 0;	%初始遗传代数
gen0= 0;	%各初始保持代数
MAXGEN= 10;	%最优个体最少保持代数
maxY= 0; %	各最优值
for i= 1:MP
    ObjV{i} = ObjectFunction(bs2rv(Chrom{i},FieldD)); %号计算各初始种群个体的目标函数值
end
Max0bjV = zeros(MP,1);	%名记录精华种群
MaxChrom = zeros(MP,PRECI*NVAR);% 各记录精华种群的编码
while gen0<= MAXGEN
    gen= gen+ 1;	%各遗传代数加1
    for i= 1:MP
        FitnV{i} = ranking(-ObjV{i});	%^各种群的适应度
        SelCh{i} = select('sus',Chrom{i}, FitnV{i} ,GGAP);%选择操作
        SelCh{i} = recombin( 'xovsp',SelCh{i}, pc(i));%交叉操作
        SelCh{i} = mut(SelCh{i},pm(i));	% 变异操作
        ObjVSel = ObjectFunction(bs2rv(SelCh{i},FieldD));%计算子代目标函数值
        [Chrom{ i},ObjV{i}] = reins(Chrom{ i},SelCh{i},1,1,ObjV{i},ObjVSel);%重插人操作
    end
    [Chrom,ObjV] = immigrant(Chrom,ObjV);	%移民操作
    [Max0bjV,MaxChrom] = EliteInduvidual(Chrom,ObjV,Max0bjV,MaxChrom); %人工选择精华种群
    YY(gen) = max(Max0bjV);	%找出精华种群中最优的个体
    
    if YY(gen)> maxY	%判断当前优化值是否与前一次优化
        %值相同
        %更新最优值
        maxY =YY(gen);
        gen0= 0;
    else
        gen0 = gen0 + 1;	%最优值保持次数加1
    end
end
%号进化过程图
plot(1:gen,YY)
xlabel('进化代数')
ylabel('最优解变化')
title('进化过程')
xlim([1,gen])
%名输出最优解
[Y,I] = max(Max0bjV);	%找出精华种群中最优的个体
X= (bs2rv(MaxChrom(I,:), FieldD));	%最优个体的解码解
disp(['最优值为: ',num2str(Y)])
disp(['对应的自变量取值;' ,num2str(X)])

%% 移民函数
function [Chrom,ObjV] = immigrant(Chrom,ObjV) %名移民算子
MP = length(Chrom);
for i=1:MP
    [MaxO,maxI] = max(ObjV{i});%找出第i种群中最优的个体
    next_i=i+1;	%目标种群(移民操作中)
    if next_i> MP;next_i = mod(next_i,MP); end
    [MinO,minI]= min(ObjV{next_i}); %找出目标种群中最劣的个体
    % 号目标种群最劣个体替换为源种群最优个体
    Chrom{next_i}(minI,:) = Chrom{ i} (maxI,:);
    ObjV{next_i}(minI)=ObjV{i}(maxI);
end
end

%% 人工选择算子
function [Max0bjV,MaxChrom] = EliteInduvidual(Chrom,ObjV,Max0bjV,MaxChrom)% 人工选择算子

MP = length( Chrom) ;%种群数
for i =1:MP
    [Max0, maxI] = max(ObjV{i});
    %找出第i种群中最优个体
    if Max0> Max0bjV(i)
        Max0bjV(i) = Max0;
        % 记录各种群的精华个体
        MaxChrom(i,:) = Chrom{ i} (maxI,:);
        %记录各种群精华个体的编码
    end
end
end

%% 目标函数
function obj=ObjectFunction(X)
col=size(X,1);
for i =1:col
    obj(i,1)=21.5+X(i,1)*sin(4*pi*X(i,1))+X(i,2)*sin(20*pi*X(i,2));
end
end
