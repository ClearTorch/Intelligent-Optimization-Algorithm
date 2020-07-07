%% ��׼�Ŵ��㷨һԪ�����Ż�
% �����к�����Сֵ
% f(x)=sin(10*pi*x)/x,x��Χ[1,2] 
% ѡ������Ʊ���
% ��Ⱥ��С40  ����Ŵ�����20  ���峤��20  ����0.95  �������0.7  �������0.01

%% ��׼�Ŵ��㷨����
clc
clear all
close all
%% ��������ͼ��
figure(1);
hold on;
lb=1;ub=2;
ezplot('sin(10*pi*x)/x',[lb,ub]);
xlabel('x');
ylabel('y');
%% �����Ŵ��㷨����
NIND=20;                                 % ��Ⱥ��С
MAXGEN=20;                               % ����Ŵ�����
PRECI=20;                                % ���峤��
GGAP=0.95;                               % ����
px=0.7;                                  % �������
pm=0.01;                                 % �������
trace=zeros(2,MAXGEN);                   % Ѱ�Ž���ĳ�ʼֵ
FieldD=[PRECI;lb;ub;1;0;1;1];            % ����������
Chrom=crtbp(NIND,PRECI);                 % ����������ɢ�����Ⱥ
%% �Ż�����
gen=0;                                   % ��������
X=bs2rv(Chrom,FieldD);                   % ��ʼ��Ⱥ�����Ƶ�ʮ����ת��
ObjV=sin(10*pi*X)./X;                    % ����Ŀ��ֵ
while gen<MAXGEN
    FitnV=ranking(ObjV);                 % ������Ӧ��ֵ
    SelCh=select('sus',Chrom,FitnV,GGAP);% ѡ��
    SelCh=recombin('xovsp',SelCh,px);    % ���飨���棩
    SelCh=mut(SelCh,pm);                 % ����
    X=bs2rv(SelCh,FieldD);               % �Ӵ������ʮ����ת��
    ObjVsel=sin(10*pi*X)./X;             % �����Ӵ���Ŀ�꺯��ֵ
    [Chrom,ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVsel); % �ز����Ӵ����������õ��µ���Ⱥ
    X=bs2rv(Chrom,FieldD);               
    gen=gen+1;                           % ������������
    % ��ȡÿ�������Ž������ţ�UΪ���Ž⣬IΪ�������
    [Y,I]=min(ObjV);
    trace(1,gen)=X(I);                   % ��¼ÿ��������ֵ
    trace(2,gen)=Y;                      
end
%% ��ͼ
plot(trace(1,:),trace(2,:),'bo');        % ����ÿ�������ŵ�
grid on;
plot(X,ObjV,'b*');
% ����ͼ
figure(2);
plot(1:MAXGEN,trace(2,:));
grid on;
xlabel('�Ŵ�����');
ylabel('��ñ仯');
bestY=trace(2,end);
bestX=trace(1,end);
fprintf(['���Ž⣺\nX=',num2str(bestX),'\nY=',num2str(bestY),'\n'])


%% ��׼�Ŵ��㷨��Ԫ�����Ż�
% �������ֵ��f(x,y)=x*cos(s*pi*x)+y*sin(2*pi*x),x��Χ[-2,2],y��Χ[-2,2]
% �Ŵ��㷨�������ã�
% ��Ⱥ��С40  ����Ŵ�����50  ���峤��40��2���Ա�����ÿ����20�� ����0.95  �������0.7  �������0.01

%% �Ŵ��㷨����
clc 
clear all 
close all 
%% ��������ͼ
figure(1);
lbx=-2;ubx=2;
lby=-2;uby=2;
ezmesh('y*sin(2*pi*x)+x*cos(2*pi*y)',[lbx,ubx,lby,uby],50);%��������ͼ��
hold on;
%% �����Ŵ��㷨����
NIND=40;                                 % ��Ⱥ��С
MAXGEN=50;                               % ����Ŵ�����
PRECI=40;                                % ���峤��
GGAP=0.95;                               % ����
px=0.7;                                  % �������
pm=0.01;                                 % �������
trace=zeros(3,MAXGEN);                   % Ѱ�Ž���ĳ�ʼֵ
FieldD=[PRECI PRECI;lbx lby;ubx uby;1 1;0 0;1 1;1 1];% ����������
Chrom=crtbp(NIND,PRECI*2);               % ����������ɢ�����Ⱥ
%% �Ż�����
gen= 0;                                  % ����������
XY= bs2rv(Chrom,FieldD);                 % ����ʼ��Ⱥ��ʮ����ת��
X= XY(:,1);Y= XY(:,2);
ObjV=Y.*sin(2*pi*X) +X.*cos(2*pi* Y);    % ������Ŀ�꺯��ֵ
while gen<MAXGEN
FitnV = ranking(-ObjV);                  % ������Ӧ��ֵ
SelCh = select('sus',Chrom,FitnV,GGAP);  % ѡ��
SelCh = recombin( 'xovsp',SelCh,px);     % ����
SelCh= mut(SelCh,pm);                    % ����
XY = bs2rv( SelCh,FieldD);               % ���Ӵ������ʮ����ת��
X= XY(:,1);Y= XY(:,2);
ObjVSel= Y.*sin(2*pi*X)+X.*cos(2*pi* Y); % �������Ӵ���Ŀ�꺯��ֵ
[Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel); % �ز����Ӵ����������õ�����Ⱥ
XY = bs2rv(Chrom ,FieldD);
gen= gen+ 1;                             % ��������������
% ��ȡÿ�������Ž⼰�����,YΪ���Ž�,IΪ��������
[Y,I]= max(ObjV);
trace(1:2,gen) = XY(I,:);                % ����ÿ��������ֵ
trace(3,gen) = Y;                        % ����ÿ��������ֵ
end
%% ��ͼ
plot3(trace(1,:), trace(2,:), trace(3,:),'bo'); % ����ÿ�������ŵ�
grid on;
plot3(XY(:,1) ,XY(:,2),ObjV,'bo');      % �滭�����һ������Ⱥ
hold off
% ������ͼ
figure(2);
plot(1 :MAXGEN, trace(3,:));
grid on
xlabel('�Ŵ�����')
ylabel('��ı仯')
title('��������')
bestZ = trace(3, end) ;
bestX = trace(1 ,end) ;
bestY = trace(2,end);
fprintf(['���Ž�:\nX= ',num2str(bestX),'\nY= ',num2str(bestY),'\nZ= ',num2str(bestZ),'\n'])

