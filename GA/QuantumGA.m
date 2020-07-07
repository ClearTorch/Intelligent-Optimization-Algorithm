%% 量子遗传算法求函数最值
clc
clear all
close all
%% 参数设置
MAXGEN = 200; % 最大遗传代数
sizepop= 40; % 种群大小
lenchrom= [20 20]; % 每个变量的二进制长度
trace = zeros(1,MAXGEN);
% 最佳个体,记录其适应度值、十进制值、二进制编码、量子比特编码
best = struct('fitness',0,'X',[],'binary',[],'chrom',[]);
% 初始化种群
chrom = InitPop( sizepop * 2, sum( lenchrom)) ;
% 对种群实施一次测量，得到二进制编码
binary = collapse( chrom);
% 求种群个体的适应度值和对应的十进制值
[fitness,X] = FitnessFunction(binary, lenchrom) ;
% 使用目标函数计算适应度
% 记录最佳个体到best
[best.fitness,bestindex] = max( fitness);
% 出最大值
best. binary = binary( bestindex,:);
best. chrom= chrom([2 * bestindex- 1:2 * bestindex],:);
best. X = X(bestindex,:);
trace(1) = best. fitness;
%fprintf('%d\n' ,1)
%% 进化
for gen = 2:MAXGEN
    %fprintf('%d\n',gen)% 提示进化代数
    % 对种群实施-次测量
    binary= collapse( chrom) ;
    % 计算适应度
    [fitness ,X] = FitnessFunction(binary,lenchrom);
    % 量子旋转门
    chrom = Qgate(chrom,fitness,best,binary);
    [newbestfitness,newbestindex] = max(fitness);  %  找到最佳值
    % 记录最佳个体到best
    if newbestfitness > best.fitness
        best.fitness = newbestfitness;
        best.binary = binary( newbestindex,:);
        best.chrom = chrom([2*newbestindex-1:2*newbestindex],:);
        best.X = X(newbestindex,:);
    end
    trace(gen) = best.fitness;
end
%% 画进化曲线
plot(1 :MAXGEN, trace);
title('进化过程');
xlabel('进化代数');
ylabel('每代的最佳适应度')
% 显示优化结果
disp(['最优解X：',num2str(best.X)]);
disp(['最大值Y：',num2str(best.fitness)]);


%% 子函数
%% 总群初始化函数 InitPop
function chrom = InitPop(M,N) % 初始化种群
% M为种群大小X2,(a和)
% N为量子比特编码长度
for i=1:M
    for j=1:N
        chrom(i,j) = 1/sqrt(2);
    end
end
end

%% 测量函数
function binary = collapse( chrom)
% 对种群实施--次测量，得到二进制编码
% 输人chrom:为量子比特编码
% 输出binary:二进制编码
[M,N]= size( chrom) ;
% 得到种群大小和编码长度
M=M/2;
% 种群大小，
binary= zeros(M,N) ;
%二进制编码大小初始化
for i=1:M
    for j=1:N
        pick = rand;
        % 产生[0,1]区间的随机数
        if pick> (chrom(2. * i- 1,j)^2)
            % 随机数大于a的平方
            binary(i,j) = 1;
        else
            binary(i,j)=0;
        end
    end
end
end

%% 量子螺旋门函数
function chrom = Qgate( chrom,fitness,best,binary) % 量子旋转门调整策略
% 输人chrom:更新前的量子比特编码
% fitness:适应度值
% best:当前种群中最优个体
% binary:二进制编码
% 各输出chrom:更新后的量子比特编码
sizepop= size(chrom,1)/2;
lenchrom= size(binary,2);
for i= 1:sizepop
    for j= 1: lenchrom
        A= chrom(2* i- 1,j);
        % a
        B= chrom(2* i,j);
        % β
        x= binary(i,j);
        b= best.binary(j);
        if ((x==0)&&(b==0))||((x== 1)&&(b== 1))
            delta= 0;% delta为旋转角的大小
            s=0; % s为旋转角的符号，即旋转方向
        elseif (x== 0)&&(b== 1)&&(fitness(i)<best.fitness)
            delta=0.01*pi;
            if A*B>0
                s=1;
            elseif A*B<0
                s=-1;
            elseif A== 0
                s=0;
            elseif B== 0
                s= sign( randn) ;
            end
        elseif (x== 0)&&(b== 1)&&(fitness(i)>= best.fitness)
            delta= 0.01*pi;
            if A* B>0
                s=-1;
            elseif A* B<0
                s=1;
            elseif A== 0
                s= sign(randn) ;
            elseif B== 0
                s=0;
            end
            
        elseif (x== 1)&&(b== 0)&&(fitness(i)< best. fitness)
            delta=0.01*pi;
            if A*B>0
                s=-1;
            elseif A*B<0
                s=1;
            elseif A== 0
                s= sign( randn) ;
            elseif B== 0
                s=0;
            end
        elseif (x== 1)&&(b== 0)&&(fitness(i)>= best. fitness)
            delta= 0.01*pi;
            if A* B>0
                s= 1;
            elseif A*B<0
                s=- 1;
            elseif A== 0
                s=0;
            elseif B== 0
                s= sign( randn) ;
            end
        end
        e= s*delta;
        %e为旋转角
        U= [cos(e) -sin(e);sin(e) cos(e)]; %量子旋转门
        y=U*[A B]';
        % y为更新后的量子位
        chrom(2* i- 1,j)= y(1);
        chrom(2* i,j) = y(2);
    end
end
end

%% 适应度函数
function [fitness,X] = FitnessFunction(binary,lenchrom)% 适应度函数
% 输人:binary二进制编码
% lenchrom各变量的二进制位数
% 输出:fitness 适应度 X十进制数(待优化参数)
sizepop = size(binary,1);
fitness = zeros(1 ,sizepop) ;

num = size( lenchrom ,2) ;
X= zeros( sizepop,num);
for i= 1:sizepop
    %使用目标函数计算适应度
    [fitness(i),X(i,:)] = Objfunction(binary(i,:),lenchrom);
end
end

%% 目标函数
function [Y,X] = Objfunction(x,lenchrom)% 目标函数
% 输人:x二进制编码  lenchrom各变量的二进制位数
% 输出:Y目标值   X十进制数
bound=[-3.0 12.1;4.1 5.8]; % 自变量的范围
% 将binary数组转换成十进制数组
X= bin2decFun(x,lenchrom,bound) ;
% 计算适应度-函数值
Y= sin(4*pi*X(1))*X(1) + sin(20*pi*X(2))*X(2);
end

%%  二进制转十进制函数
function X = bin2decFun(x,lenchrom,bound)% 二进制转化成十进制
% 输人x二进制编码
%  lenchrom各变量的二进制位数
% bound变量的范围
% 输出:X进制数
M= length( lenchrom) ;
n= 1;
X= zeros(1,M) ;
for i=1:M
    for j = lenchrom(i)-1:- 1:0
        X(i) = X(i) + x(n).*2.^j;
        n=n+ 1;
    end
end
X= bound(:,1)'+ X./(2.^lenchrom - 1).*(bound(:,2)-bound(:,1))';
end

