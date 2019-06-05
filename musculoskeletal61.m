% ----------------------------------------------------------------------------
%
% file name : musculoskeletal6.m
%
% purpose : particle filter for estimation of musculoskeletal system parameter
%
% ----------------------------------------------------------------------------
 
function [] = musculoskeletal61()
 
close all;  %close all Figure already shown
clearvars;  %clear from memory
 
disp('Particle Filter (PF) estimation program start');

global dt;
dt = 0.00192; % sampling period [sec]
load angular_velocity.mat
load sourceProcess.mat
load muscle_strength_estimation.mat

time = 0;   %initial time 
endtime = 4.0051; %end time[sec]

 
nSteps = ceil((endtime - time)/dt); %mainループのfor文用. とりあえずは
                                    % 初期時刻 timeを 0と置いた時に endtime[sec] 
                                    %を dt[sec]で割ったものを nStepsとする.
                                    %ceilは正の無限大方向の丸目
                                    
%--------各時刻における数値格納用の構造体(グラフ出力に使う)-------------------
result.time=[];     %経過時刻
result.xEst=[];     %推定値
result.y=[];
result.torq=[];


%------------------実験結果-------------------------------
extension = [ex_IEMG1, ex_IEMG2, ex_Angle, ex_angular_velocity, ex_Torq];
flexion = [f_IEMG1, f_IEMG2, f_Angle, f_angular_velocity, f_Torq];
%--------------------------------------------------------------------------

%-----------------------muscle strength expectation------------------------
% [r_f, r_e, angle, angular_velocity, torq] = experiment(flexion);%flexion
[r_f, r_e, angle, angular_velocity, torq] = experiment(extension);%extension


result.obtorq=[torq];
a_f = r_f/r_f_max;
a_e = r_e/r_e_max;

u_f = a_f * f_f_o;
u_e = a_e * f_e_o;


%-----------状態の推定値--------
xEst = zeros(6,1);

xEst(1) = angle(1);
xEst(2) = angular_velocity(1);
xEst(3) = u_f(1);
xEst(4) = u_e(1);
xEst(5) = 0.2;
xEst(6) = 0.05;


NP=500; %パーティクル数(数は推定対象と条件に応じて適当に試行錯誤で決める)
NTh=NP/1.0; %resamplingを実施する有効パーティクル数

px=repmat(xEst,1,NP);%パーティクル格納変数
pw=zeros(1,NP)+1/NP;%重み変数(尤度)

va = zeros(6,1); % システム雑音の分布の標準偏差
va(1) = 0.0001;
va(2) = 0.0001;
va(3) = 1;
va(4) = 1;
va(5) = 0.000001;
va(6) = 0.000001;
vb = zeros(6,1); % システム雑音の分布平均

%観測値
  y = [angle, angular_velocity, u_f, u_e ];    % observed value y

  
tic;    %timewatch start---------------------------------------------------

% ----------------------------- Main loop --------------------------------
for i=1 : nSteps
    time = time + dt;
    
    d = 0.03;
    l = 0.0015;
    v = va.*randn(6,1)+vb;  %random number for system noise
%     v(5)=0;
%     v(6)=0;
    %a*randn+b=(標準偏差aと平均bの正規分布から得た乱数)
    
    
    %------観測雑音------------
     a = zeros(4,4); %observation noise(ob) standard deviation
    a(1,1) = 0.001;
    a(2,2) = 0.001;
%     a(3,3) = 3000; %flexion
%     a(4,4) = 100;
    a(3,3) = 100; %extension
    a(4,4) = 3000;
    b = zeros(4,1);  % ob average
    %-------------係数行列
    c = zeros(4,6);
    c(1,1) = 1;
    c(2,2) = 1;
    c(3,3) = 1;
    c(4,4) = 1;
    
%-------------------------- Particle Filter -------------------------------

    for ip=1:NP
        x=px(:,ip); % state of i th particle
        w=pw(ip);   % likelihood of i th particle
        
        vv = va.*randn(6,1)+vb; % randomnumber for system noise
        x = f(x, torq(i), d, l, vv, dt);
        
        w = Gauss4(abs(y(i,1:4)'-(c*x)), b, a);   % calculate a likelihood of a particle y=cx+w
    
        px(:,ip)=x; %save
        pw(ip)=w; %save
    end
    
    pw=Normalize(pw,NP); % Normalize
    [px,pw]=Resampling(px,pw,NTh,NP);%resampling
    xEst=px*pw'; % final estimation = expectation
    
    torqEst = -d*(xEst(1)*xEst(3)*xEst(5)+xEst(2)*xEst(3)*xEst(6)+xEst(1)*xEst(4)*xEst(5)+xEst(2)*xEst(4)*xEst(6)-xEst(3)+xEst(4));
%     torqEst = abs(torqEst);
%     k_1 = 0.2;
%     b_1 = 0.05;
%     y(i,5) = -d*(y(i,1)*y(i,3)*(k_1)+y(i,2)*y(i,3)*(b_1)+y(i,1)*y(i,4)*(k_1)+y(i,2)*y(i,4)*(b_1)-y(i,3)+y(i,4));
    %グラフ出力用の数値格納

    result.time=[result.time; time];
    result.xEst=[result.xEst;xEst'];
    result.y = [result.y;y(i,:)];
    result.torq = [result.torq;torqEst];
    
end
toc;   %stopwatch stop


DrawGraph(result);      %グラフ出力
save('musculoskeletal61.mat');   %matファイル保存(ファイル名はカッコ内で指定)
end
% ----------------------------- Main loop finish --------------------------

function [px,pw]=Resampling(px,pw,NTh,NP) %resampling function
% Low Variance Sampling algorithm
Neff=1.0/(pw*pw');
if Neff<NTh %resampling
    wcum=cumsum(pw);
    base=cumsum(pw*0+1/NP)-1/NP;%乱数を加える前のbase
    resampleID=base+rand/NP;%ルーレットを乱数分増やす
    ppx=px;%データ格納用
    ind=1;%新しいID
    for ip=1:NP
        while(resampleID(ip)>wcum(ind))
            ind=ind+1;
        end
        px(:,ip)=ppx(:,ind);%LVSで選ばれたparticleに置き換え
        pw(ip)=1/NP;%尤度の初期化
    end
end
end

function pw=Normalize(pw,NP) %likelihood normalize fuction

sumw=sum(pw);
if sumw~=0
    pw=pw/sum(pw); %正規化
else
    pw=zeros(1,NP)+1/NP;
end

end

function p=Gauss4(X,U,Sigma)
%ガウス分布(共分散)の確率密度を計算する関数
%p=1/( ((2*pi)^(1/(size(Sigma,1)))) * (sqrt(det(Sigma))) ) * (exp( (-1/2) * dot((X-repmat(U,1,size(X,2)))' * (inv(Sigma)) , (X-repmat(U,1,size(X,2)))',2)));
p=1/( ((2*pi)^((size(Sigma,1))/2)) * (sqrt(det(Sigma))) ) * (exp( (-1/2) * dot((X-repmat(U,1,size(X,2)))' * (inv(Sigma)) , (X-repmat(U,1,size(X,2)))',2)));
end

function x = f(x,T, d, i, v, dt)    %state model
    A = zeros(6, 6);
    A(1,2)=1;
    B = zeros(6, 6);
    B(2,1) = x(3)*x(5);
    B(2,2) = x(3)*x(6);
    B(2,3) = -1;
    B(2,4) = 1;
    B(2,5) = x(1)*x(4);
    B(2,6) = x(2)*x(4);
    C = zeros(6,1);
    C(2) = T;
    x = x + dt*(A*x - (d/i)*B*x-(1/i)*C) + v;
     %state equation
end


 function [a, b, c, d, e] = experiment(title)
 
      a = title(:,1);
      b = title(:,2);
      c = title(:,3);
      d = title(:,4);
      e = title(:,5);
 end
 
function []=DrawGraph(result)   % drawing graph
 
    figure(1);
    hold off;
    x=[ result.y(:,1) result.xEst(:,1)];    %観測値と推定値を一つの行列にまとめる
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;   %推定値
    plot(result.time, x(:,1),'b','linewidth', 1); hold on;  %観測値
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('Angle \theta [rad]','fontname','times new roman','fontsize', 20);
    grid on;

    figure(2);
    hold off;
    x=[ result.y(:,2) result.xEst(:,2)];    
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;   
    plot(result.time, x(:,1),'b','linewidth', 1); hold on; 
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('$${\rm Angular Velocity} \dot \theta \ {\rm [rad/s]}$$', 'interpreter', 'latex','fontname','times new roman','fontsize', 20)
    grid on;

    figure(3);
    hold off;
    x=[ result.y(:,3) result.xEst(:,3)];   
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;  
    plot(result.time, x(:,1),'b','linewidth', 1); hold on; 
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('Contractile Force {\it u_f} [N]','fontname','times new roman','fontsize', 20);
    grid on;

    figure(4);
    hold off;
    x=[ result.y(:,4) result.xEst(:,4)];  
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;   
    plot(result.time, x(:,1),'b','linewidth', 1); hold on;  
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('Contractile Force {\it u_e} [N]','fontname','times new roman','fontsize', 20);
    grid on;

    figure(5);
    hold off;
    x=[ result.xEst(:,5)];    
    plot(result.time, x(:,1),'r--o','linewidth', 1); hold on; 
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    ylim([0.18 0.22])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('Elasticity Coefficient {\it k} [N/m] ','fontname','times new roman','fontsize', 20);
    grid on;

    figure(6);
    hold off;
    x=[result.xEst(:,6)];   
    plot(result.time, x(:,1),'r--o','linewidth', 1); hold on; 
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    ylim([0.04 0.06])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('Viscosity Coefficient {\it b} [N.s/m] ','fontname','times new roman','fontsize', 20);
    grid on;
    
     figure(7);
    hold off;
    x=[result.torq result.obtorq];   
    plot(result.time, x(:,1),'r--o','linewidth', 1); hold on;    
    plot(result.time, x(:,2),'b','linewidth', 1); hold on;  
    set( gca, 'fontname','times new roman','fontsize',20 );
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 20);
    ylabel('Torque \tau [Nm]','fontname','times new roman','fontsize', 20);
    grid on;

   
end