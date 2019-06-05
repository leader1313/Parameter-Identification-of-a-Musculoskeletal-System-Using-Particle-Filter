clear

load sourceProcess.mat

dt = 0.00192;

ex_angular_velocity =[0; diff(ex_Angle)/dt];
f_angular_velocity =[0; diff(f_Angle)/dt];


save('angular_velocity.mat');

