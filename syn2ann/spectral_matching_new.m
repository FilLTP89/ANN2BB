%%%%%%%%%%%%%%%%%%%%
%Funzione per accelerogrammi spettro-compatibili da accelerogrammi reali%
%%%%%%%%%%%%%%%%%%%%
function  [out_t,out_acc,out_vel,out_dis,out_T,out_Se,out_freq,out_FAS] = ...
        spectral_matching_new(target_Se,acc_,niter,tol_u,tol_l)   

%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target_Se. acc spectrum should be in m/s2, first column should be periods
% acc_. time history should be in m/s2, first column should be time vector
% %5/12 Maria : ho passato tutto in cm/s2
% 9/12/16 AGO: please see the modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol_upp = tol_u;
tol_low = tol_l;

warning off;
 
dt = acc_(2,1)-acc_(1,1);
T_vec = target_Se(:,1)';

T_corr_ini=  T_vec(2);     % added by AGO on 9/12/16
% Comment AGO: it corrresponds to T=0.05 s and we will correct PGA later
T_corr_fin = 3;            % added by AGO 
% Comment AGO: ideally T=3 second should be governed by SPEED runs        

target = target_Se;
pga_target = target_Se(1,2); 

T_in = T_vec';T_in(1)=1e-10;
Sp_in = target_Se(:,2);
%
t = acc_(:,1);
acc = acc_(:,2).*100; %cm/s2    %Added by Maria 5/12
% acc = tukeywin(length(acc),0.10).*acc;  % added by Maria on 5/12/2016 

%
l_zeros=length(acc);
acc=[zeros(1,length(acc)) acc' zeros(1,length(acc))]';
t=[zeros(1,length(t)) t' zeros(1,length(t))]';
%

acc_in = acc;
npun = length(t);

Nfft = 2^nextpow2(length(t)); 
durata = (Nfft-1).*dt;
df = 1/durata;
freq = 0:df:1/2/dt;
acc(npun+1:Nfft) = 0;
ACC=dt*fft(acc,Nfft); 

for i=2:length(T_in)
    
Sp_acc(i)=4*pi^2*disp_spectra(acc',0.01,T_in(i),0.05)./T_in(i)^2;
freq_in(length(T_in)-i+2)=1/T_in(i);
Rapp_spe_in(length(T_in)-i+2)=Sp_acc(i)/Sp_in(i); % fattore di correzione
end

Sp_acc(1) = max(abs(acc(:)));
freq_in(length(T_in)+1)=1/2/dt; % (*)...questo funziona solo se il la fNyq è <= di 1/T_in(1)..
Rapp_spe_in(length(T_in)+1)=Sp_acc(1)/Sp_in(1);
freq_in(1)=0;
Rapp_spe_in(1)=1;

Rapp_spe = interp1(freq_in,Rapp_spe_in,freq,'linear'); % interpola il rapp alle freq dell'acc a correggere
acc_pro=acc; % si tiene l'acc non corretto!
nfreq=length(freq);


for k=1:niter  % lavora tante volte sull'accelerogramma (niter = per tutti)
	ACC_PRO=fft(acc_pro)*dt; % trasformata dell'acc da correggere
    
	for m=1:nfreq
 		if (freq(m)>=1/T_corr_fin) && (freq(m)<=1/T_corr_ini) % cambiato qui
        %if (freq(m)>=1/T_corr_fin)% cambiato qui Maria 5/12/2016
		ACC_PRO(m)=ACC_PRO(m)/Rapp_spe(m); % correzione della trasformata dell'acc
		end % cambiato qui
    end;
    
	ACC_PRO(nfreq+1)=0;
    
	for m=nfreq+2:2*nfreq
		ACC_PRO(m)=conj(ACC_PRO(2*nfreq-m+2));
    end;
    
    acc_pro=real(ifft(ACC_PRO))/dt; 
    acc_pro=detrend(acc_pro,'constant');

	for i=1:length(T_in)
    Sp_acc_pro(i,k)=4*pi^2*disp_spectra(acc_pro',0.01,T_in(i),0.05)./T_in(i)^2;
    end
	
	Sp_acc_pro(1,k)=max(abs(acc_pro));
    
        for i=1:length(T_in)
            %%% modified on 16/06/16 by AGO
            rat=Sp_acc_pro(i,k)/Sp_in(i);
%             
%             if rat>1 && rat <= 1+tol_upp || rat < 1 && rat >= 1-tol_low
%                 rat=1;
%             end
               
            Rapp_spe_pro(length(T_in)-i+2)=rat; % aggiorna il rapporto di correzione
            %%% modified on 16/06/16 by AGO
            
        end
        
        Rapp_spe_pro(1)=1;
        Rapp_spe = interp1(freq_in,Rapp_spe_pro,freq,'linear');
    
end % fine iterazioni

% added by AGO on 9/12/16
acc_pro2=acc_pro;
acc_in2=acc_in;
clear acc_pro;clear acc_in;
acc_pro_dummy=acc_pro2(l_zeros+1:length(acc_pro2));
acc_in_dummy=acc_in2(l_zeros+1:length(acc_in2));
acc_pro=acc_pro_dummy(1:length(acc_));
acc_in=acc_in_dummy(1:length(acc_));
clear acc_pro2;clear acc_in2;
clear acc_pro_dummy;clear acc_in_dummy;
clear t;
t = [0:dt:(length(acc_pro)-1)*dt]; % added by Maria on 6/12/2016 h 22.43
% end of modification by AGO on 9/12/16


% k = find(acc_in> 0.1);                                           % added by Maria on 6/12/2016 
% N1=k(1);
% w1 = zeros(N1,1);
% N2 = length(acc_pro)-N1;
% % N2 = length(acc_in)-N1;
% w2 = rectwin(N2);
% % N3 = length(acc_pro)-length(acc_in);
% % w3 = zeros(N3,1);
% acc_pro(1:N1) = acc_pro(1:N1).*w1;
% acc_pro((N1+1):(N1+N2)) = acc_pro((N1+1):(N1+N2)).*w2;
% % acc_pro((N1+N2+1):(N1+N2+N3)) = acc_pro((N1+N2+1):(N1+N2+N3)).*w3;
% % acc_pro=detrend(acc_pro,'constant');

% t_in = t; % added by Maria on 6/12/2016 
% N = min(length(acc_in),length(acc_pro));
% id = find(abs(acc_in(1:N) - acc_pro(1:N))<0.0001);



% Arias Intensity part added by Maria on 5/12/16

% Modified by AGO on 9/12/16
% Calculating the velocity time histories
vel_p=cumsum(acc_pro).*dt;
vel_in=cumsum(acc_in).*dt;
% End of modification

[t_idx,idx,Ain] = arias_intensity(acc_in,dt,0.005);% added by Maria on 6/12/2016 h.22.13 
% Comment by MI: in order to use the hybrid until the t1 for which the arias intensity is 0.01
% Comment by AGO: I found that Ia=0.5% is removing the unwanted
% initial oscillations and providing better displacement time histories


r = [(t_idx-100*dt):dt:(t_idx+100*dt)]';
a1 = acc_in((idx-100):(idx+100));
a2 = acc_pro((idx-100):(idx+100));

% modified by AGO on 9/12/16
% Comment by AGO: doing the same thing for velocities
v1 = vel_in((idx-100):(idx+100));
v2 = vel_p((idx-100):(idx+100));
% End of modification

for j = 1:length(r)
    %d(j) = abs((a1(j)-a2(j)); 
    d(j) = abs((v1(j)-v2(j))/v1(j));
    % Comment 1 by AGO: the control mechanism is velocity since it directly
    % affects the displacement time history
    % Comment 2 by AGO: rather than finding the minimum difference , it has
    % been found that relative minimum may be a better option
end
m=find(d==min(d));
t1=r(m);
s1 = length(0:dt:t1);
acc_pro(1:s1)=acc_in(1:s1);
% End of Arias Intensity part added by Maria on 5/12/16

acc_pro = tukeywin(length(acc_pro),0.10).*acc_pro;  % added by Maria on 5/12/2016

for times=1:5  % added by AGO on 18/11/2016 repeating the PGA
               % to be sure to take care of the PGA reversals

% CS on 03.05.2016                                                            
% adjust PGA in time domain 
[pga ipga] = max(abs(acc_pro)); 

if pga < (1-tol_low)*abs(pga_target) || pga > (1+tol_upp)*abs(pga_target) 

dt_cor = 0.04; % modified by AGO on 9/12/16, it was dt_cor = 0.05                                                          % added by AGO on 18/11/2016 it was dt=0.05                                   
npun_cor = round(dt_cor./(t(2)-t(1)));     
if mod(npun_cor,2)==0
	npun_cor = npun_cor+1; 
end
    
x = [ t(ipga-(npun_cor-1)/2),t(ipga),t(ipga+(npun_cor-1)/2)]; 
y = [acc_pro(ipga-(npun_cor-1)/2),pga_target.*sign(acc_pro(ipga)),...
    acc_pro(ipga+(npun_cor-1)/2)];
xi = t(ipga-npun_cor/2:ipga+npun_cor/2-1); 
yi_c = interp1(x,y,xi,'cubic'); 

yi_o=acc_pro(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2);
yi_d=yi_c'-yi_o;

acc_pro(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = yi_c;

vel_pro = cumsum(acc_pro)*dt;
% modification by AGO
% Comment by AGO: in the original version this part was not present
% we are now subtracting what we induced residual velocities to match PGA
% from the following points of the input acc-t.
acc_pro(ipga+(npun_cor+1)/2:ipga+(3*npun_cor-1)/2) = ...
acc_pro(ipga+(npun_cor+1)/2:ipga+(3*npun_cor-1)/2)-yi_d;
% end of modification by AGO
end

end
% End of PGA correction

vel_pro = cumsum(acc_pro)*dt;

% added by AGO on 9/12/16
% correcting the tail now
[t_idx2,idx2,Ain2] = arias_intensity(acc_pro,dt,0.9999)% 
% calculating the residual velocity from the last 100 samples (last 1
% seconds).
%vel_shift=mean(vel_pro(length(vel_pro)-99:length(vel_pro)));
vel_shift=mean(vel_pro(length(vel_pro)-99:length(vel_pro)));
vel_corr=-vel_shift;

% now attributing to the acceleration point at selected Arias Intensity
acc_pro(idx2+1)=acc_pro(idx2+1)+vel_corr/dt;

clear vel_pro
% end of modification


vel_pro = cumsum(acc_pro)*dt;
dis_pro = cumsum(vel_pro)*dt;


for i=2:length(T_in)
        Sp_acc_pro_final(i)=4*pi^2*disp_spectra(acc_pro,dt,T_in(i),0.05)./T_in(i)^2;
end
 	Sp_acc_pro_final(1)=max(abs(acc_pro));
    
    T_in(1) = 0;
	out_T = T_in';
	out_Se = Sp_acc_pro_final';
    out_t = [0:dt:(length(acc_pro)-1)*dt]; % added by Maria on 6/12/2016 h 22.43
	out_acc = acc_pro;
    out_vel = vel_pro;
    out_dis = dis_pro;
    out_FAS = abs(ACC_PRO(1:length(freq)));
    out_freq = freq;
end



    