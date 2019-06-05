clear
close all;

load experiment.mat
extension = [e_1,e_2,e_3,e_4,time];
flexion = [f_1,f_2,f_3,f_4,time];



[ex_IEMG1,ex_IEMG2,ex_Angle,ex_Torq,time] = Process(extension);
[f_IEMG1,f_IEMG2,f_Angle,f_Torq,time] = Process(flexion);

save('sourceProcess.mat');

function [a,b,c,d,f] = Process(e)


    EMG_1 = e(:,1);
    EMG_2 = e(:,2);
    Angle = e(:,3);
    Torq = e(:,4);
    time = e(:,5);
 
        figure(1);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            ax1 = subplot(4,1,1); plot(time, EMG_1,'b','linewidth', 1); 
            ylabel(ax1, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax1, 'ch1');
            ax2 =subplot(4,1,2); plot(time, EMG_2,'g','linewidth', 1);   
            ylabel(ax2, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax2, 'ch2');
            ax3 =subplot(4,1,3); plot(time, Angle,'y','linewidth', 1); 
            ylabel(ax3, '$degree$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax3, 'ch3');
            ax4 =subplot(4,1,4); plot(time, Torq,'m','linewidth', 1); 
            xlabel('$time(sec)$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            ylabel(ax4, '$Nm$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(ax4, 'ch4');
            grid on;
            
    EMG_1_rectified = abs(EMG_1);
    EMG_2_rectified = abs(EMG_2);
    
    
    
    N   = 100;        % FIR filter order
    Fp  = 1;       % 20 kHz passband-edge frequency
    Fs  = 521;       % 96 kHz sampling frequency
    Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
    Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation

 
    lowpassFilt = dsp.LowpassFilter('DesignForMinimumOrder',false, ...
            'FilterOrder',N,'PassbandFrequency',Fp,'SampleRate',Fs,...
            'PassbandRipple',0.01, 'StopbandAttenuation',80);
    lowpassFilt2 = dsp.LowpassFilter('DesignForMinimumOrder',false, ...
            'FilterOrder',200,'PassbandFrequency',Fp,'SampleRate',Fs,...
            'PassbandRipple',0.01, 'StopbandAttenuation',80);
        
        
        
    EMG_1_filltered = lowpassFilt(EMG_1_rectified);
    EMG_2_filltered = lowpassFilt(EMG_2_rectified);
    Angle_filltered = lowpassFilt(Angle);
    Torq_filltered = lowpassFilt(Torq);
        
        
            figure(2);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            axf1 = subplot(4,1,1); plot(time, EMG_1_filltered ,'b','linewidth', 1); 
            xlim([1 5])
            ylabel(axf1, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf1, 'ch1');
            axf2 = subplot(4,1,2); plot(time, EMG_2_filltered,'g','linewidth', 1);  
            xlim([1 5])
            ylabel(axf2,'$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf2, 'ch2');
            axf3 = subplot(4,1,3); plot(time, Angle_filltered,'y','linewidth', 1);  
            xlim([1 5])
            ylabel(axf3, '$degree$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf3, 'ch3');
            axf4 = subplot(4,1,4); plot(time, Torq_filltered,'m','linewidth', 1); 
            xlim([1 5])
            xlabel('$time(sec)$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            ylabel(axf4, '$Nm$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
            title(axf4, 'ch4');
            grid on;

    cut_EMG_1 = [];
    cut_EMG_2 = [];
    cut_Angle = [];
    cut_Torq = [];
    cut_time = [];
    
   %cutting at 520 
   for i = 520 : 3125
   cut_EMG_1 = [cut_EMG_1;EMG_1_filltered(i)];
   cut_time = [cut_time;time(i)];
   cut_EMG_2 = [cut_EMG_2;EMG_2_filltered(i)];
   cut_Angle = [cut_Angle;Angle_filltered(i)];
   cut_Torq = [cut_Torq;Torq_filltered(i)];
   end
   cut_time = cut_time-0.9984;
%    cut_EMG_1 = cut_EMG_1-cut_EMG_1(1);
%    cut_EMG_2 = cut_EMG_2-cut_EMG_2(1);
   cut_Torq = -(cut_Torq-cut_Torq(1))*(10/1.2);

   cut_Angle =-(cut_Angle)*(6.28/25);
   mA = mean(cut_Angle);
    cut_Angle =cut_Angle-mA;

%     cut_Angle =-(cut_Angle-cut_Angle(1))*(6.28/25);

   a=[];
   b=[];
   c=[];
   d=[];
   f=[];
   
   for i = 1 : 2086
   a = [a;cut_EMG_1(i)];
   b = [b;cut_EMG_2(i)];
   c = [c;cut_Angle(i)];
   d = [d;cut_Torq(i)];
   f = [f;cut_time(i)];
   end
   
   figure(3);
    hold off;
    set(gca, 'fontsize', 16, 'fontname', 'times');
    axc1 = subplot(4,1,1); plot(f, a ,'b','linewidth', 1); 
    xlim([0 4])
    ylabel(axc1, '$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc1, 'EMG1');
    axc2 = subplot(4,1,2); plot(f, b,'g','linewidth', 1);  
    xlim([0 4])
    ylabel(axc2,'$V$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc2, 'EMG2');
    axc3 = subplot(4,1,3); plot(f, c,'y','linewidth', 1);  
    xlim([0 4])
    ylabel(axc3, '$degree$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc3, 'Angle');
    axc4 = subplot(4,1,4); plot(f, d,'m','linewidth', 1); 
    xlim([0 4])    
    xlabel('$time(sec)$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    ylabel(axc4, '$Nm$', 'fontsize', 16, 'fontname', 'times','Interpreter','latex');
    title(axc4, 'Torque');
    grid on;
    
end
   

    

