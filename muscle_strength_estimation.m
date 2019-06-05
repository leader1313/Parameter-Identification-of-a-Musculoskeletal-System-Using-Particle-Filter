clear 
close all

load experiment.mat
isometric = [iso_e_1,iso_e_2,iso_f_1,iso_f_2,time];
[is_e_1,is_e_2,is_f_1,is_f_2,cut_time] = isometricProcess(isometric);

iso_f_torq_max = 25.54;
iso_e_torq_max = -18.54;

% is_e_1 = mean(is_e_1(780:1304));
% is_e_2 = mean(is_e_2(780:1304));
% is_f_1 = mean(is_f_1(1200:1720));
% is_f_2 = mean(is_f_2(1200:1720));

% is_e_1 = max(is_e_1);
% is_e_2 = max(is_e_2);
% is_f_1 = max(is_f_1);
% is_f_2 = max(is_f_2);


r_e_max = max(is_e_2);
r_f_max = max(is_f_1);

n = find(is_f_1==r_f_max);
r_e = is_f_2(n);
m = find(is_e_2==r_e_max);
r_f = is_e_1(m);


d = 0.03 ; %moment arm
iso_f_a_f = 1;
iso_f_a_e = r_e/r_e_max;
iso_e_a_f = r_f/r_f_max;
iso_e_a_e = 1;


[f_f_o, f_e_o] = f(iso_f_torq_max,iso_e_torq_max,iso_f_a_f,iso_f_a_e,iso_e_a_f,iso_e_a_e );





save('muscle_strength_estimation.mat', 'f_f_o', 'f_e_o', 'r_e_max', 'r_f_max');


function [x, y] = f(torq_max_f,torq_max_e,a_f_f,a_e_f,a_f_e,a_e_e)
syms a b
d = 0.03;
eq_f = a_f_f*a - a_e_f*b - torq_max_f/d;
eq_e = a_f_e*a - a_e_e*b - torq_max_e/d;

[A,B] = equationsToMatrix([eq_f, eq_e], [a,b]);
z = linsolve(A,B);
temp = double(z);
x = temp(1);
y = temp(2);
end


 function [a,b,c,d,h] = isometricProcess(e)


    e_EMG_1 = e(:,1);
    e_EMG_2 = e(:,2);
    f_EMG_1 = e(:,3);
    f_EMG_2 = e(:,4);
    time = e(:,5);
 
        figure(1);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            ax1 = subplot(4,1,1); plot(time, e_EMG_1,'b','linewidth', 1); 
            ylabel(ax1, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax1, 'ech1');
            ax2 =subplot(4,1,2); plot(time, e_EMG_2,'g','linewidth', 1);   
            ylabel(ax2, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax2, 'ech2');
            ax3 =subplot(4,1,3); plot(time, f_EMG_1,'y','linewidth', 1); 
            ylabel(ax3, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax3, 'fch1');
            ax4 =subplot(4,1,4); plot(time, f_EMG_2,'m','linewidth', 1); 
            xlabel('$time(sec)$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            ylabel(ax4, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax4, 'fch2');
            grid on;
            
    e_EMG_1_rectified = abs(e_EMG_1);
    e_EMG_2_rectified = abs(e_EMG_2);
    f_EMG_1_rectified = abs(f_EMG_1);
    f_EMG_2_rectified = abs(f_EMG_2);
    
    
    
    N   = 100;        % FIR filter order
    Fp  = 1;       % 20 kHz passband-edge frequency
    Fs  = 521;       % 96 kHz sampling frequency
    Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
    Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation

 
    lowpassFilt = dsp.LowpassFilter('DesignForMinimumOrder',false, ...
            'FilterOrder',N,'PassbandFrequency',Fp,'SampleRate',Fs,...
            'PassbandRipple',0.01, 'StopbandAttenuation',80);

        
        
        
    e_EMG_1_filtered = lowpassFilt(e_EMG_1_rectified);
    e_EMG_2_filtered = lowpassFilt(e_EMG_2_rectified);
    f_EMG_1_filtered = lowpassFilt(f_EMG_1_rectified);
    f_EMG_2_filtered = lowpassFilt(f_EMG_2_rectified);
        
        
            figure(2);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            axf1 = subplot(4,1,1); plot(time, e_EMG_1_filtered ,'b','linewidth', 1); 
            ylabel(axf1, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf1, 'ech1');
            axf2 = subplot(4,1,2); plot(time, e_EMG_2_filtered,'g','linewidth', 1);  
            ylabel(axf2,'$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf2, 'ech2');
            axf3 = subplot(4,1,3); plot(time, f_EMG_1_filtered,'y','linewidth', 1);  
            ylabel(axf3, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf3, 'fch1');
            axf4 = subplot(4,1,4); plot(time, f_EMG_2_filtered,'m','linewidth', 1); 
            xlabel('$time(sec)$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            ylabel(axf4, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf4, 'fch2');
            grid on;

    cut_e_EMG_1 = [];
    cut_e_EMG_2 = [];
    cut_f_EMG_1 = [];
    cut_f_EMG_2 = [];
    cut_time = [];
    
   %cutting at 520 
   for i = 520 : 3125
   cut_e_EMG_1 = [cut_e_EMG_1;e_EMG_1_filtered(i)];
   cut_time = [cut_time;time(i)];
   cut_e_EMG_2 = [cut_e_EMG_2;e_EMG_2_filtered(i)];
   cut_f_EMG_1 = [cut_f_EMG_1;f_EMG_1_filtered(i)];
   cut_f_EMG_2 = [cut_f_EMG_2;f_EMG_2_filtered(i)];
   end
   cut_time = cut_time-0.9984;
   cut_e_EMG_1 = cut_e_EMG_1-cut_e_EMG_1(1);
   cut_e_EMG_2 = cut_e_EMG_2-cut_e_EMG_2(1);
   cut_f_EMG_1 = cut_f_EMG_1-cut_f_EMG_1(1);
   cut_f_EMG_2 = cut_f_EMG_2-cut_f_EMG_2(1);
   
   a=[];
   b=[];
   c=[];
   d=[];
   h=[];
   
   for i = 1 : 2086
   a = [a;cut_e_EMG_1(i)];
   b = [b;cut_e_EMG_2(i)];
   c = [c;cut_f_EMG_1(i)];
   d = [d;cut_f_EMG_2(i)];
   h = [h;cut_time(i)];
   end
   
   figure(3);
    hold off;
    set(gca, 'fontsize', 16, 'fontname', 'times');
    axc1 = subplot(4,1,1); plot(h, a ,'b','linewidth', 1); 
    ylabel(axc1, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc1, 'ech1');
    axc2 = subplot(4,1,2); plot(h, b,'g','linewidth', 1);  
    ylabel(axc2,'$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc2, 'ech2');
    axc3 = subplot(4,1,3); plot(h, c,'y','linewidth', 1);  
    ylabel(axc3, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc3, 'fch1');
    axc4 = subplot(4,1,4); plot(h, d,'m','linewidth', 1); 
    xlabel('$time(sec)$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    ylabel(axc4, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc4, 'fch2');
    grid on;
    
end
      
    