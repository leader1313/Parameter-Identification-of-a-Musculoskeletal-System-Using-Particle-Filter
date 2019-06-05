clear
close all;

load sourceProcess.mat
load angular_velocity.mat
extension = [ex_IEMG1,ex_IEMG2,ex_Angle,ex_angular_velocity,ex_Torq,time];
flexion = [f_IEMG1,f_IEMG2,f_Angle,f_angular_velocity,f_Torq,time];



%Process(extension);
Process(flexion);


function [] = Process(e)


    EMG_1 = e(:,1);
    EMG_2 = e(:,2);
    Angle = e(:,3);
    Angulr_velrocity=e(:,4);
    Torq = e(:,5);
    time = e(:,6);
 
        figure(1);
            hold off;
            plot(time, EMG_1,'b','linewidth', 1); 
            set( gca, 'fontname','times new roman','fontsize', 16 );
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('IEMG1 [V]','fontname','times new roman','fontsize', 16);
            grid on;
            
        figure(2);
            hold off;
            plot(time, EMG_2,'g','linewidth', 1); 
            set( gca, 'fontname','times new roman','fontsize',16 );
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('IEMG2 [V]','fontname','times new roman','fontsize', 16);
            grid on;
        
         figure(3);
            hold off;
            plot(time, Angle,'y','linewidth', 1); 
            set( gca, 'fontname','times new roman','fontsize',16 );
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('Angle\theta [rad]','fontname','times new roman','fontsize', 16);
            grid on;
        
         figure(4);
            hold off;
            plot(time, Angulr_velrocity,'r','linewidth', 1); 
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('$${\rm Angular Velocity} \dot \theta \ {\rm [rad/s]}$$', 'interpreter', 'latex','fontname','times new roman','fontsize', 16)
%             ylabel('AngularVelocity {\dot{\theta}} [rad/s]','fontname','times new roman','fontsize', 16);
            grid on;
            
          figure(5);
            hold off;
            plot(time, Torq,'m','linewidth', 1); 
            set( gca, 'fontname','times new roman','fontsize',16 );
            xlim([0 4])
            xlabel('Time [s]','fontname','times new roman','fontsize', 16);
            ylabel('Torque [Nm]','fontname','times new roman','fontsize', 16);
            grid on;
            
            
end
   

    
