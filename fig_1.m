function [] = fig_1()
 
close all;  %close all Figure already shownclear
load musculoskeletal6.mat

DrawGraph(result); 
end


function []=DrawGraph(result)   % drawing graph
 
    figure(1);
    hold off;
    x=[ result.y(:,1) result.xEst(:,1)];    %観測値と推定値を一つの行列にまとめる
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;   %推定値
    plot(result.time, x(:,1),'b','linewidth', 1); hold on;  %観測値
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('Angle [rad]','fontname','times new roman','fontsize', 16);
    legend('Estimate','observe');
    grid on;

    figure(2);
    hold off;
    x=[ result.y(:,2) result.xEst(:,2)];    
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;   
    plot(result.time, x(:,1),'b','linewidth', 1); hold on;  
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('Angularvalocity [rad/s]','fontname','times new roman','fontsize', 16);
    legend('Estimate','observe');
    grid on;

    figure(3);
    hold off;
    x=[ result.y(:,3) result.xEst(:,3)];   
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;  
    plot(result.time, x(:,1),'b','linewidth', 1); hold on; 
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('u_f [N]','fontname','times new roman','fontsize', 16);
    legend('Estimate','observe');
    grid on;

    figure(4);
    hold off;
    x=[ result.y(:,4) result.xEst(:,4)];  
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,2),'r--o','linewidth', 1); hold on;   
    plot(result.time, x(:,1),'b','linewidth', 1); hold on;    
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('u_e [N]','fontname','times new roman','fontsize', 16);
    legend('Estimate','observe');
    grid on;

    figure(5);
    hold off;
    x=[ result.xEst(:,5)];    
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,1),'r--o','linewidth', 1); hold on; 
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('k ','fontname','times new roman','fontsize', 16);
    legend('Estimate');
    grid on;

    figure(6);
    hold off;
    x=[result.xEst(:,6)];   
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,1),'r--o','linewidth', 1); hold on; 
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('b ','fontname','times new roman','fontsize', 16);
    legend('Estimate');
    grid on;
    
     figure(7);
    hold off;
    x=[result.torq result.obtorq];   
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.time, x(:,1),'r--o','linewidth', 1); hold on;    
    plot(result.time, x(:,2),'g','linewidth', 1); hold on;   
    xlim([0 4])
    xlabel('Time [s]','fontname','times new roman','fontsize', 16);
    ylabel('Torque [Nm]','fontname','times new roman','fontsize', 16);
    legend('Estimate', 'observe');
    grid on;

   
end