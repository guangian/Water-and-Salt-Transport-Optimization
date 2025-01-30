function [Fgb,Objx,gb,PositionValue]=CDEO(UB,LB,PVN,MH)
global MPN MIT
% [Score,Position,Convergence,PositionValue]
% %% Problem parameters
% par = Cal_par(func_num);	% Optimization problem parameters:
% PVN = par.n;					% Dimensionality
% gn = par.g;					% Number of inequality constraints
% hn = par.h;					% Number of equality constraints
% LB = par.xmin;			% Lower bound of the variables
% UB = par.xmax;			% Upper bound of the variables
% MH=@fobj1;
% %% Maximum Function Evaluations and Checkpoints
% if PVN <= 10
%     max_fes = 1e5;
% elseif PVN <= 30
%     max_fes = 2e5;
% elseif PVN <= 50
%     max_fes = 4e5;
% elseif PVN <= 150
%     max_fes = 8e5;
% else
%     max_fes = 1e6;
% end
% MPN=50; %质点数量
% MIT=max_fes/MPN; %最大迭代次数
PVP=0.976;%种群变异概率
PGWC=0.05;%水动力弥散系数权重系数
SWI=0.15;%土壤初始含水率
SWS=0.05;%土壤饱和含水率
SDR=0.0005;%饱和土壤扩散率
SHC=0.0005;%饱和土壤导水率
MDD=0.5;%机械弥散度
DSI=1.5;%弥散度指数
MPV=0.95;%平均孔隙流速
DWA=1.5;%溶质在自由水中的扩散系数
IDB=2;%溶质在自由水中的扩散系数指数
MIPG=-2.5;%最小土壤水力(水势)梯度
MAPG=2.5; %最大土壤水力(水势)梯度
UB=UB.*ones(1,PVN);
LB=LB.*ones(1,PVN);
x=rand(MPN,PVN);
WPG=rand(MPN,PVN);%初始化水动力弥散系数
Dh=rand(MPN,PVN);%初始化分子扩散系数
Dm=rand(MPN,PVN);%初始化分子扩散系数
SWC=rand(MPN,PVN);%初始化土壤含水率
gb=zeros(MIT,1);
%%%%%%%%%初始化质点个体位置和水动力弥散系数%%%%%%%%%
Gup=MAPG*ones(PVN,1);
Glp=MIPG*ones(PVN,1);
% x=initial(MPN,PVN,UB,LB);
for i=1:PVN
    for j=1:MPN
        x(j,i)=(UB(i)-LB(i))*rand+LB(i);%生成[LB UB]范围内的初始质点位置
        WPG(j,i)=WPG(j,i)*(MAPG-MIPG)+MIPG;%生成[MAPG MIPG]范围内的土水势
    end
end

%%%%%%%%%初始化个体最优位置和最优值%%%%%%%%%
gbest=x;%个体最优位置
for j=1:MPN
    Fpb(j)=MH(x(j,:));%个体最优值
end
% [Fpb]=MH(x,func_num);%个体最优值
%%%%%%%%%初始化全局最优位置和最优值%%%%%%%%%
Objx=x(1,:);
Fgb=Fpb(1);
for j=1:MPN
    if Fpb(j)<Fgb
        Objx=x(j,:);%全局最优位置
        Fgb=Fpb(j);%全局最优值
    end
end
%///////////////////////////////////////
%计算质点初始位置初始时刻土壤含水率
for i=1:PVN
    for j=1:MPN
        SWC(j,i)=SWC(j,i)*(SWS-SWI)+SWI;%生成[SWS,SWI]范围内土壤初始含水率
    end
end
%//////////////////////////////////////////////////////////////
%%%%%%%%%按照公式依次迭代直到满足精度或者迭代次数%%%%%%%%%
for jk=1:MIT
    %%%%%%%%%更新质点水动力弥散系数%%%%%%%%%
    w=PGWC*(1+log(jk/MIT)/log(MIT));
    for j=1:MPN
        %根据拉普拉斯变换后的余误差函数计算当前位置当前迭代次数下的土壤含水率
        SWC(j,:)=((SWS-SWI)/2)*(erfc((x(j,i)-SHC*jk)/(2*(MIT*SDR)^0.5))+exp(x(j,i)*SHC/SDR)*erfc((x(j,i)+SHC*jk)/(2*(MIT*SDR)^0.5)))+SWI;
        %//////////////////////////////////////////////////////////////
        R1=rand;
        R2=rand;
        Dh(j,:)=MDD*MPV^(DSI*log(jk/MIT));%dispersity(MDD),Mean pore velocity(MPV),Dispersion index(DSI)(取值2.5)
        Dm(j,:)=DWA*exp(IDB*SWC(j,:));%DMA为分子扩散系数，SWC为土壤初始含水率，IDB为分子扩散指数(取值1.25)
        WPG(j,:)=w*WPG(j,:)+Dh(j,:)*R1.*(gbest(j,:)-x(j,:))+Dm(j,:)*R2.*(Objx-x(j,:));% 更新水动力弥散系数
        %%%%%%%%%%%%%%%通过上下限Gup和Glp对压力水头(压力势)进行约束%%%%%%%%%%%%%%%
        for i=1:PVN
            if  WPG(j,i)>Gup(i)
                WPG(j,i)=Gup(i);
            end
            if  WPG(j,i)<Glp(i)
                WPG(j,i)=Glp(i);
            end
        end
        %%%%%%%%%通过边界条件约束更新质点位置%%%%%%%%%
        x(j,:)=x(j,:)+WPG(j,:);
        for i=1:PVN
            if  x(j,i)>UB(i)
                x(j,i)=UB(i);
            end
            if  x(j,i)<LB(i)
                x(j,i)=LB(i);
            end
        end
        %%%%%%%%%进行自适应变异%%%%%%%%%
        if rand>PVP
            i=ceil(PVN*rand);%将PVN*rand四舍五入到大于或等于该元素的最接近整数。
            x(j,i)=LB(i)+(UB(i)-LB(i))*rand;
        end
        %%%%%%%%%计算当前个体的适应度%%%%%%%%%
        FFb(j)=MH(x(j,:));
        %%%%%%%%%计算当前全局的适应度%%%%%%%%%
        if FFb(j)<Fpb(j)
            gbest(j,:)=x(j,:);
            Fpb(j)=FFb(j);
        end
        if Fpb(j)<Fgb
            Objx= gbest(j,:);
            Fgb=Fpb(j);
        end
    end
    gb(jk)= Fgb;%记录历代全局最优值
    FMV=min(gb);%把MIT个适应度最小值赋值给FMV
    PositionValue{1,jk}=Objx;%储存最佳参数
    PositionValue{2,jk}=Fgb;%储存最佳适应度值
    PositionValue{3,jk}=x;%储存种群参数
    PositionValue{4,jk}=Fpb;%储存种群适应度值
end
end
