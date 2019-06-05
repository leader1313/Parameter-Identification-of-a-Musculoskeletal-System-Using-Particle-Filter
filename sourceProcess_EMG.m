clear
close all;

load experiment.mat
extension = [e_1,e_2,time];
flexion = [f_1,f_2,time];



[ex_IEMG1,ex_IEMG2,time] = Process(extension);
[f_IEMG1,f_IEMG2,time] = Process(flexion);

save('sourceProcess.mat');

function [a,b,c,d,f] = Process(e)


    EMG_1 = e(:,1);
    EMG_2 = e(:,2);
    time = e(:,3);
 
        figure(1);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            ax1 = subplot(2,1,1); plot(time, EMG_1,'b','linewidth', 1); 
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('CH1 [V]','fontname','times new roman','fontsize', 16);
            ax2 =subplot(2,1,2); plot(time, EMG_2,'g','linewidth', 1);   
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('CH2 [V]','fontname','times new roman','fontsize', 16);
            grid on;
            
    EMG_1_rectified = abs(EMG_1);
    EMG_2_rectified = abs(EMG_2);
    
     figure(2);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            axf1 = subplot(2,1,1); plot(time, EMG_1_rectified ,'b','linewidth', 1); 
            xlim([1 5])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('EMG2 [V]','fontname','times new roman','fontsize', 16);
            axf2 = subplot(2,1,2); plot(time, EMG_2_rectified,'g','linewidth', 1);  
            xlim([1 5])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('EMG2 [V]','fontname','times new roman','fontsize', 16);
            grid on;
    
    
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
       
    cut_EMG_1 = [];
    cut_EMG_2 = [];
    cut_time = [];
    
   %cutting at 520 
   for i = 520 : 3125
   cut_EMG_1 = [cut_EMG_1;EMG_1_filltered(i)];
   cut_time = [cut_time;time(i)];
   cut_EMG_2 = [cut_EMG_2;EMG_2_filltered(i)];
   end
   cut_time = cut_time-0.9984;

   a=[];
   b=[];
   c=[];
   d=[];
   f=[];
   
   for i = 1 : 2086
   a = [a;cut_EMG_1(i)];
   b = [b;cut_EMG_2(i)];
   f = [f;cut_time(i)];
   end
   
           figure(3);
            hold off;
            set(gca, 'fontsize', 16, 'fontname', 'times');
            axc1 = subplot(2,1,1); plot(f, a ,'b','linewidth', 1); 
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('IEMG1 [V]','fontname','times new roman','fontsize', 16);
            axc2 = subplot(2,1,2); plot(f, b,'g','linewidth', 1);  
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('IEMG2 [V]','fontname','times new roman','fontsize', 16);
            grid on;
    
end
   
