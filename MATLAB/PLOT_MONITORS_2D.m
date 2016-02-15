%% PLOT THE RESULTS
clear all
close all
clc

path1 = '../INPUTS_DG/';
path2 = '../INPUTS_DG/';


t0 = 0;
T = 20;

Num_of_tot_mon = 2;

for i = 1 : Num_of_tot_mon
    
    if i < 10
        fileName = ['monitor0000',num2str(i),'.d'];
    elseif i < 100
        fileName = ['monitor000',num2str(i),'.d'];
    elseif i < 1000
        fileName = ['monitor00',num2str(i),'.d'];
    elseif i < 10000
        fileName = ['monitor0',num2str(i),'.d'];
    elseif i < 100000
        fileName = ['monitor',num2str(i),'.d'];
    end
    
    sol_1 = load([path1,fileName]);
    sol_2 = load([path2,fileName]);
     
    
    
    figure(i)
    subplot(211)
    hold on
    plot(sol_1(:,1),sol_1(:,2),'g-','LineWidth',2);hold on; grid on;
    plot(sol_2(:,1),sol_2(:,2),'k-','LineWidth',1);
%     plot(sol_1(:,1), -sin(sqrt(2)*pi*sol_1(:,1))*sin(pi*0.25)^2*sin(2*pi*0.5),'k-','LineWidth',1);
    xlim([t0 T]);
    %ylim([-1 1]);
    xlabel('t (s)');
    ylabel('u_x (m)');
    
    subplot(212)
    plot(sol_1(:,1),sol_1(:,3),'g-','LineWidth',2);hold on; grid on;
    plot(sol_2(:,1),sol_2(:,3),'k-','LineWidth',1);
%     plot(sol_1(:,1), sin(sqrt(2)*pi*sol_1(:,1))*sin(2*pi*0.25)*sin(pi*0.5)^2,'k-','LineWidth',1);
    xlim([t0 T]);
    %ylim([-1 1]);
    xlabel('t (s)');
    ylabel('u_y (m)');
    
    pause
    %close all
    
end