%% 标准遗传算法的TSP问题优化

clear all
clc
close all

X=[16.47,96.10     % 城市坐标
    16.47,94.44
    20.09,92.54
    22.39,93.37
    25.23,97.24
    22.00,96.05
    20.47,97.02
    17.20,96.29
    16.30,97.38
    14.05,98.12
    16.53,97.38
    21.52,95.59
    19.41,97.13
    20.09,92.55];
%% 参数设置
NIND=100;          % 种群大小
MAXGEN=200;        % 最大迭代数
Pc=0.9;	           % 各交叉概率
Pm=0.05;           % 变异概率
GGAP=0.9;          % 代沟(generationgap)
D=Distance(X);	   % 生成距离矩阵
N=size(D,1);	   % 34+34)
Chrom=InitPop(NIND,N);	% 初始化种群
%% 在二维图上画出所有坐标点
figure(1)
plot(X(:,1),X(:,2),'o');
% 各画出随机解的路线图
DrawPath(Chrom(1,:),X)
pause(0.0001)
% 输出随机解的路线和总距离
disp('初始种群中的一个随机值:')
OutputPath(Chrom(1,:));
Rlength=PathLength(D,Chrom(1,:));
disp(['总距离:',num2str(Rlength)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
% 优化
gen=0;
figure(2);
hold on;box on;
xlim([0,MAXGEN])
title('优化过程')
xlabel('代数')
ylabel('最优值')
ObjV=PathLength(D,Chrom);	% 计算路线长度
pre0bjV=min(ObjV);
while gen<MAXGEN
    %  计算适应度
    ObjV=PathLength(D,Chrom);	% 计算路线长度
    %fprintf('ed   者1.10f\n',gen,min(ObjV))
    line([gen - 1,gen],[pre0bjV,min(ObjV)]) ;
    pause(0.0001)
    pre0bjV = min(ObjV);
    FitnV = Fitness(ObjV) ;
    % 选择
    SelCh = Select(Chrom,FitnV,GGAP);
    % 交叉操作
    SelCh = Recombin(SelCh,Pc);
    % 变异
    SelCh = Mutate(SelCh,Pm);
    % 操作
    Se1Ch = Reverse(SelCh,D);
    % 重插人子代的新种群
    Chrom = Reins(Chrom,SelCh,ObjV);
    % 更新迭代次数
    gen=gen+1;
end
%% 画出最优解路线
ObjV=PathLength(D,Chrom);
[minObjV,minInd]=min(ObjV);
DrawPath(Chrom(minInd(1),:),X);
%输出最优解的路线和总距离
disp('最优解');
p=OutputPath(Chrom(minInd(1),:));
disp(['总距离：',num2str(ObjV(minInd(1)))]);
disp('---------------------------------------------------------------------')

%% 子函数定义
%% 种群初始化函数 InitPop
function Chrom = InitPop(NIND,N)
% 初始化种群
% 输入:
% NIND:种群大小
% N:个体染色体长度(这里为城市的个数)
% 输出:
% 初始种群

Chrom = zeros(NIND,N);% 用于存储种群
for i= 1:NIND
    Chrom(i,:) = randperm(N);% 随机生成初始种群
end
end

%% 适应度函数 Fintness
function FitnV=Fitness(len)
% 输入 len（个体的长度，即TSP距离）
% 输出 FitnV个体的适应度值
FitnV=1./len;
end

%% 进化翻转操作Reverse
function SelCh = Reverse(SelCh,D)       % 进化逆转函数
%输人:
%SelCh  被选择的个体
%D  各城市的距离矩阵
%输出:  Se1Ch进化逆转后的个体
[row,col] = size(SelCh);
ObjV= PathLength(D,SelCh); %各计算路线长度
SelCh1 = SelCh;
for i= 1:row
    r1 = randsrc(1,1,[1:col]);
    r2 = randsrc(1,1,[1:col]);
    mininverse = min([r1 r2]);
    maxinverse = max([r1 r2]);
    Se1Chl(i,mininverse:maxinverse) = SelCh1(i,maxinverse:- 1:mininverse);
end
ObjVl= PathLength(D, SelCh1); % 计算路线长度
index =ObjVl<ObjV;
SelCh( index,:) = SelCh1(index,:);
end

%% 画图函数DrawPath
function DrawPath(Chrom,X)   % 画路线图函数
%输人:
% Chrom 待画路线
% X 各城市的坐标位置
R= [Chrom(1,:) Chrom(1,1)];
% 一个随机解(个体)
figure;
hold on
plot(X(:,1),X(:,2),'o','color' ,[0.5,0.5,0.5])
plot(X(Chrom(1,1),1) ,X(Chrom(1,1) ,2), 'rv','MarkerSize',20)
for i= 1:size(X,1)
    text(X(i,1)+ 0.05,X(i,2)+0.05,num2str(i),'color',[1,0,0]);
end
A= X(R,:);
row= size(A,1);
for i= 2:row
    [arrowx,arrowy] = dsxy2figxy(gca,A(i - 1:i,1),A(i- 1:i,2));% 坐标转换
    annotation( 'textarrow' ,arrowx, arrowy, 'HeadWidth' ,8, 'color' ,[0,0,1]);
end
hold off
xlabel('横坐标')
ylabel('纵坐标')
title('轨迹图')
box on

    function varargout = dsxy2figxy(varargin)
        if length(varargin{1}) == 1 && ishandle(varargin{1}) ...
                && strcmp(get(varargin{1},'type'),'axes')
            hAx = varargin{1};
            varargin = varargin(2:end);
        else
            hAx = gca;
        end
        if length(varargin) == 1
            pos = varargin{1};
        else
            [x,y] = deal(varargin{:});
        end
        axun = get(hAx,'Units');
        set(hAx,'Units','normalized');
        axpos = get(hAx,'Position');
        axlim = axis(hAx);
        axwidth = diff(axlim(1:2));
        axheight = diff(axlim(3:4));
        if exist('x','var')
            varargout{1} = (x - axlim(1)) * axpos(3) / axwidth + axpos(1);
            varargout{2} = (y - axlim(3)) * axpos(4) / axheight + axpos(2);
        else
            pos(1) = (pos(1) - axlim(1)) / axwidth * axpos(3) + axpos(1);
            pos(2) = (pos(2) - axlim(3)) / axheight * axpos(4) + axpos(2);
            pos(3) = pos(3) * axpos(3) / axwidth;
            pos(4) = pos(4) * axpos(4 )/ axheight;
            varargout{1} = pos;
        end
        set(hAx,'Units',axun)
    end

end

%% 距离函数Distance
function D= Distance(a) % 计算两两城市之间的距离
%输入a各城市的位置坐标
% 输出D 两城市之间的距离
row= size(a,1);
D= zeros(row,row);
for i = 1:row
    for j=i+ 1:row
        D(i,j) = ((a(i,1)- a(j,1))^2+ (a(i,2)- a(j,2))^2)^0.5;
        D(j,i) = D(i,j);
    end
end
end

%% 计算各体路线长度函数 PathLength
function len=PathLength(D,Chrom)
% 输入D 两城市之间的距离  Chrom 个体的轨迹
[~,col]= size(D);
NIND= size(Chrom,1) ;
len= zeros(NIND,1);
for i = 1:NIND
    p= [Chrom(i,:) Chrom(i,1)];
    i1 = p(1:end- 1);
    i2 = p(2:end);
    len(i,1) = sum(D((i1- 1)*col + i2));
end
end

%% 输出路线函数outputParth
function p = OutputPath(R)
% 输出路线函数
% 输人R 路线
R= [R,R(1)];
N =length(R);
p = num2str(R(1));
for i= 2:N
    p=[p,'-> !',num2str(R(i)) ] ;
end
disp(p)
end

%% 选择操作函数Select
function SelCh = Select( Chrom,FitnV,GGAP) % 选择操作
% 输人:Chrom种群
% FitnV适应度值
%  GGAP选择概率
% 输出: Se1Ch被选择的个体
NIND= size(Chrom,1);
NSel = max( floor(NIND* GGAP+ .5),2);
ChrIx = Sus(FitnV,NSel);
SelCh = Chrom(ChrIx,:);
%  其中,函数Sus的代码为:
    function NewChrIx = Sus(FitnV,Nsel)
        % 输人:FitnV个体的适应度值
        % Nsel 被选择个体的数目
        % 输出:
        %NewChrIx被选择个体的索引号
        [Nind,ans]= size(FitnV);
        cumfit = cumsum(FitnV) ;
        trials= cumfit(Nind) / Nsel * (rand+ (0:Nsel- 1));
        Mf = cumfit(:, ones(1, Nsel));
        Mt= trials(:, ones(1, Nind))';
        [NewChrIx,ans] = find(Mt < Mf&[ zeros(1, Nsel); Mf(1:Nind-1, :)] <= Mt);
        [ans,shuf] = sort( rand(Nsel, 1));
        NewChrIx = NewChrIx( shuf);
    end
end
%%
function SelCh=Recombin(SelCh,Pc) % 交叉操作
% 输入:SelCh被选择的个体
% Pc交叉概率
% 输出:SelCh交叉后的个体
NSel = size( SelCh,1) ;
for i= 1:2:NSel-mod( NSel,2)
    if Pc>=rand
        % 交叉概率Pc
        [SelCh(i,:),SelCh(i+1,:)]=intercross(SelCh(i,:),SelCh(i+1,:));
    end
end
%  其中,函数intercross的代码为:
    function [a,b] = intercross(a,b)
        % 输人:a和b为两个待交叉的个体
        % 输出:a和b为交叉后得到的两个个体
        L= length(a);
        r1 = randsrc(1,1,[1:L]);
        r2 = randsrc(1,1,[1:L]);
        if r1~=r2
            a0= a;b0= b;
            s= min([r1,r2]);
            e= max([r1,r2]);
            for k= s:e
                al =a;b1= b;
                a(k)= b0(k);
                b(k) = a0(k);
                x= find(a== a(k));
                y= find(b== b(k));
                il=x(x~=k);
                i2= y(y~= k);
                if ~isempty( il)
                    a(i1) = a1(k);
                end
                if ~isempty( i2)
                    b(i2)= b1(k);
                end
            end
        end
    end
end

%% 
function SelCh = Mutate( SelCh, Pm) % 变异操作
% 输入:SelCh被选择的个体  Pm变异概率
% 输出:SelCh变异后的个体
[ NSel,L]= size(SelCh) ;
for i= 1:NSel
    if Pm>= rand
        R = randperm(L) ;
        SelCh(i,R(1:2)) = SelCh(i,R(2:-1:1));
    end
end
end

%% 
function Chrom = Reins( Chrom, SelCh,ObjV)% 重插人子代的新种群
% 输入:Chrom父代的种群  SelCh子代种群  ObjV 父代适应度
% 输出:Chrom组合父代与子代后得到的新种群
NIND = size(Chrom,1);
NSel = size(SelCh,1);
[~,index] = sort(ObjV) ;
Chrom = [Chrom( index(1 :NIND - NSel),:);SelCh];
end
