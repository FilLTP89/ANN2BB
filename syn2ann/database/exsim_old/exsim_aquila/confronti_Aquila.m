clear all
close all

st='AQK'
% load observed:
SP_AQ=load( ['D:\lavoro\DPC-INGV_progetti\progettoS2_2007-2009\Aquila\registrazioni\main\',st,'_\',st,'_PSA.DAT']);
% SP_AQ_SD=load(['..\..\..\Aquila\registrazioni\main\',st,'_\',st,'_SD.DAT']);
% TH_AQ=load( ['..\..\..\Aquila\registrazioni\main\',st,'_\',st,'_acc.DAT']);
% FS_AQ=load( ['..\..\..\Aquila\registrazioni\main\',st,'_\',st,'_FS_ACC.DAT']); 

for l=1:20
    sp_syn=load(['.\walters\',st,'\',st,'_0604__acc_s001_',num2str(l),'.sa']);

    fig3=figure(3)
    set(gca,'Fontname','Arial','Fontsize',14)
    plot(sp_syn(:,1),sp_syn(:,2)./100,'k','linewidth',2); hold on;
    plot(SP_AQ(:,1),SP_AQ(:,2:3),'linewidth',3); hold on;
    xlabel('\itT \rm[s]')
    ylabel('\itPSA \rm[m/s^2]')
    grid on
    xlim([0 3])
end

