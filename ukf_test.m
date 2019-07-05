%ukf_test
function ukf_test
clc;clear;close all;
T = 1;%��������
N = 100;
Xstation = [200,300];
X = zeros(4,N);
Z = zeros(1,N);
X(:,1) = [-100,2,200,20]';%x,vx,y,vy.
F = [1,T,0,0;
        0,1,0,0;
        0,0,1,T;
        0,0,0,1];
 delta_w = 1e-3;
 Q = delta_w*diag([0.5,1]);%������������
 G = [T^2/2,0 ;
         T,0;
         0,T^2/2;
         0,T];
 R = 5;%������������
 v = sqrtm(R)*randn(1,N);
%������ʵֵ
for i = 1:N-1
    X(:,i+1)=F*X(:,i) + G*sqrtm(Q)*randn(2,1); %�õط�ȱ����������
end

%��������ֵ
for i = 1:N
    Z(i) = Dis(Xstation,X(:,i)) + v(i);
end

%��ʼUK�������˲��ķ���
%1������UK��
L = 4;
alpha = 1;
beta = 2;
lamda = 3- L;
P0 = eye(4);
for k= 1:N
    %1������sigma��
    for i = 1:2*L+1
        Wm(i) = 1/(2*(L+lamda));
        Wc(i)=Wm(i);
    end
    Wm(1)=lamda/(lamda+L);
    Wc(1) = lamda/(lamda+L) + (1-alpha^2+beta);
    P = P0;
    cho = (chol(P*(lamda + L)))';
    Xsigma = zeros(4,2*L+1);
    Xestimate = X(:,k);
    for i = 1:L %����ֵXesitmate��Χ��Xsigma��
        XsigmaPos(:,i) = Xestimate + cho(:,i);
        XsigmaNeg(:,i) = Xestimate - cho(:,i);
    end
    Xsigma = [Xestimate,XsigmaPos,XsigmaNeg];
    %2��һ��Ԥ��
    Xpredsigma = F*Xsigma;
    %3������sigma��ľ�ֵ��Э����
    Xpred = zeros(4,1);
    for i = 1:2*L+1
        Xpred = Xpred + Wm(i)*Xpredsigma(:,i);%��Xpredsigma��ľ�ֵ
    end
    Ppred = zeros(4,4);
    for i = 1:2*L + 1
        Ppred = Ppred + Wc(i)*(Xpredsigma(:,i)-Xpred)*(Xpredsigma(:,i)-Xpred)';
    end
    Ppred = Ppred + G*Q*G';
    %4����ȡ�µ�sigma�㣬Ԥ��ֵXpred��Χ��sigma��
    chor = (chol((lamda+L)*Ppred))';
    for i = 1:L
        XaugSigmaPos(:,i) = Xpred + chor(:,i);
        XaugSigmaNeg(:,i) = Xpred - chor(:,i);
    end
    XaugSigma=[Xpred , XaugSigmaPos,XaugSigmaNeg];
    %5���۲�Ԥ��
    for i = 1:2*L+1
        ZsigmaPre(i) = Dis(Xstation,XaugSigma(:,i));
    end
    %6���۲�Ԥ��ֵ�ľ�ֵ��Э����
    Zpred = 0;
    for i = 1:2*L + 1
        Zpred = Zpred + Wm(i)*ZsigmaPre(i);
    end
    Pzz = 0;
    for i = 1:2*L+1
        Pzz = Pzz + Wc(i)*(ZsigmaPre(i)- Zpred)*(ZsigmaPre(i)- Zpred)';
    end
    Pzz = Pzz + R;
    %7������Ԥ��ֵ�Ͳ���ֵЭ����
    Pxz = zeros(4,1);
    for i = 1:2*L+1
        Pxz = Pxz + Wc(i)*(XaugSigma(:,i)- Xpred) * (ZsigmaPre(i)- Zpred)';
    end
    %8�����㿨�����˲�����
    K = Pxz * inv(Pzz);
    %9������Ԥ��ֵ��Э���
    Xestimate = Xpred + K *(Z(k)-Zpred);
    P = Ppred - K*Pzz*K';
    Xaug(:,k+1) = Xestimate;
    P0 = P;
end
%%%%
figure
hold on ;
plot(X(1,:),X(3,:),'k.');
plot(Xaug(1,:),Xaug(3,:),'r+');


function y = Dis(Xstation,X)
    y = sqrt((Xstation(1) - X(1))^2 + (Xstation(2) - X(3))^2);
end