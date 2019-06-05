clear




experiment = 'experiment(19.01.30).xlsx';
%------------Calibration-----------%

%�p�x��(360/25)�܂���(6.28/25)������
%Torq��-35.2����, (10/1.2)������

[f_1,f_2,f_3,f_4,e_1,e_2,e_3,e_4,iso_e_1,iso_e_2,iso_f_1,iso_f_2,time] = MakingMat(experiment);


save('experiment.mat');



function [a,b,c,d,e,f,g,h,i,j,k,l,t] = MakingMat(filename)

    a = xlsread(filename, 'F:F');
    b = xlsread(filename, 'G:G');
    c = xlsread(filename, 'H:H');
    d = xlsread(filename, 'I:I');
    e = xlsread(filename, 'K:K');
    f = xlsread(filename, 'L:L');
    g = xlsread(filename, 'M:M');
    h = xlsread(filename, 'N:N');
    i = xlsread(filename, 'R:R');
    j = xlsread(filename, 'S:S');
    k = xlsread(filename, 'U:U');
    l = xlsread(filename, 'V:V');
    t = xlsread(filename, 'X:X');
end


