%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  无迹Kalman滤波在目标跟踪中的应用
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 该地方说明：标准的卡尔曼滤波是通过直接样本X进行一步转移变化，这里通过在样本周围生成
%sigma点，sigma点的均值和X相同，sigma点的方差X相同，对sigma点进行操作，避免了又线性
%变换到非线性变换的过程。该处F虽然是线性变化，也可以用于非线性变换，该处H为非线性变换
%UKF是对非线性函数的概率密度分布进行近似，有一系列确定样本来逼近状态的后验概率密度
%而不是对非线性函数进行近似，不需要对jacobian矩阵进行求导。
%由于ukf没有把高阶项忽略掉，因此对非线性分布的高阶项有较高的计算精度。有效克服了
%扩展卡尔曼滤波估计精度低和稳定性差的问题。
function UKF
clc;clear;close all;
T=1; %雷达扫描周期
N=60/T; %总的采用次数
X=zeros(4,N); % 目标真实位置和速度
X(:,1)=[-100,2,200,20]; %目标初始位置和速度
Z=zeros(1,N); %传感器对位置的观测
delta_w=1e-3; %如果增大该参数，真实轨迹变成曲线
Q=delta_w*diag([0.5,1]) ;%过程噪声均值
G=[T^2/2,0;T,0;0,T^2/2;0,T];%过程噪声驱动矩阵
R=5; %观测噪声方差
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%状态转移矩阵
x0=200;y0=300;
Xstation=[x0,y0]; %雷达站位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UKF滤波，UT变换
w=sqrtm(R)*randn(1,N);
for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);%产生真实值
end
for t=1:N
    Z(t)=Dist(X(:,t),Xstation)+w(t);%产生观测值
end
%产生sigma点参数
L=4;
alpha=1;
kalpha=0;
belta=2;
ramda=3-L;
for j=1:2*L+1
    Wm(j)=1/(2*(L+ramda));%sigma散列点均值权重
    Wc(j)=Wm(j);%sigma散列点协方差权重
end
Wm(1)=ramda/(L+ramda);
Wc(1)=ramda/(L+ramda)+1-alpha^2+belta;%权值计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xukf=zeros(4,N);
Xukf(:,1)=X(:,1);%无迹卡尔曼滤波状态初始化，真实值第一个值
P0=eye(4);%协方差矩阵初始
for t=2:N
    xestimate= Xukf(:,t-1);
    P=P0;
    %1、获取sigma点集
    cho=(chol(P*(L+ramda)))';
    for k=1:L
        xgamaP1(:,k)=xestimate+cho(:,k);
        xgamaP2(:,k)=xestimate-cho(:,k);
    end
    Xsigma=[xestimate,xgamaP1,xgamaP2];%sigma点集
    %2、对sigma点集，进一步预测
    Xsigmapre=F*Xsigma;
    %3、利用上一步结果，计算均值和协方差
    Xpred=zeros(4,1);%均值
    for k=1:2*L+1
        Xpred=Xpred+Wm(k)*Xsigmapre(:,k);
    end
    Ppred=zeros(4,4);%协方差
    for k=1:2*L+1
        Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)';
    end
    Ppred=Ppred+G*Q*G';
    %4、根据预测值，再一次使用ut变换，得到新的sigma点集
    chor=(chol((L+ramda)*Ppred))';
    for k=1:L
        XaugsigmaP1(:,k)=Xpred+chor(:,k);%正sigma点
        XaugsigmaP2(:,k)=Xpred-chor(:,k);%负sigma点，对称变换
    end
    Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2];
    %5、观测预测
    for k=1:2*L+1
        Zsigmapre(1,k)=hfun(Xaugsigma(:,k),Xstation);
    end
    %6、计算观测预测均值和协方差
    Zpred=0; % 观测均值
    for k=1:2*L+1
        Zpred=Zpred+Wm(k)*Zsigmapre(1,k);
    end
    Pzz=0;
    for k=1:2*L+1
        Pzz=Pzz+Wc(k)*(Zsigmapre(1,k)-Zpred)*(Zsigmapre(1,k)-Zpred)';
    end
    Pzz=Pzz+R;%得到协方差Pzz

    Pxz=zeros(4,1);
    for k=1:2*L+1
        Pxz=Pxz+Wc(k)*(Xaugsigma(:,k)-Xpred)*(Zsigmapre(1,k)-Zpred)';%均值和测量值协方差矩阵
    end
    %7、计算kalman增益
    K=Pxz*inv(Pzz);
    %8、更新状态和协方差
    xestimate=Xpred+K*(Z(t)-Zpred);%更新状态
    P=Ppred-K*Pzz*K';%更新协方差
    P0=P;
    Xukf(:,t)=xestimate;
end
%误差分析
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,i),Xukf(:,i));%滤波后的误差
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');
plot(Xukf(1,:),Xukf(3,:),'-r+');%ukf轨迹
legend('真实轨迹','UKF轨迹')
figure
hold on; box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=Dist(X1,X2)
if length(X2)<=2
    d=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );
else
    d=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );
end
function [y]=hfun(x,xx) %该方程为观测方程
y=sqrt((x(1)-xx(1))^2+(x(3)-xx(2))^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
