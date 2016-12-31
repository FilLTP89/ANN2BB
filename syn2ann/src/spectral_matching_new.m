%%%%%%%%%%%%%%%%%%%%
%Funzione per accelerogrammi spettro-compatibili da accelerogrammi reali%
%%%%%%%%%%%%%%%%%%%%
function  [out_t,out_acc,out_vel,out_dis,out_T,out_Se,out_freq,out_FAS,n_] = ...
    spectral_matching_new(target_Se,acc_,niter,tol_u,tol_l,tol_u_pga,tol_l_pga,comp)   

%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target_Se. acc spectrum should be in m/s2, first column should be periods
% acc_. time history should be in m/s2, first column should be time vector
% 5/12/16 Maria : ho passato tutto in cm/s2
% 9/12/16 AGO   : please see the mods.
% 14/12/16 AGO  : please see the mods.
% 15/12/16 AGO  : please see the mods. Version number is added (v1.0)
% 16/12/16 AGO  : please see the mods. Version number is updated (v2.0)
% 17/12/16 AGO  : first official release (v2.1)
%                 all unused comments are removed, developers are advised
%                 to refer previous internal alpha versions.
% 18/12/16 AGO  : minor rearrangements, variability is removed with mult=1
%                 (line 102)

% Current version is v. 3.0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol_upp = tol_u;
tol_low = tol_l;

tol_upp_pga=tol_u_pga;
tol_low_pga=tol_l_pga;

T_sim=0.75; % accuracy limit period for numerical simulations
T_sc=0.6*T_sim; %represents the short periods, with an added confidence

warning off;
 
dt = acc_(2,1)-acc_(1,1);
T_vec = target_Se(:,1)';


T_corr_ini=  T_vec(2);% T_corr_ini=0.05s     

T_corr_fin = 3;       % ideally T=3s should be governed by SPEED runs        

target = target_Se;
pga_target = target_Se(1,2); 

T_in_=T_vec';
Sp_in_=target_Se(:,2);

T_in(1)=T_in_(1);T_in(2)=1/49;T_in(3:length(T_in_)+1)=T_in_(2:length(T_in_));
Sp_in(1)=Sp_in_(1);
Sp_in(2)=Sp_in_(1)+(Sp_in_(2)-Sp_in_(1))*(T_in(2)-T_in_(1))/(T_in_(2)-T_in_(1));
Sp_in(3:length(T_in_)+1)=Sp_in_(2:length(T_in_));

T_corr_ini=  T_in(2);% T_corr_ini=0.05s 

Sp_in2=Sp_in; %for variability purposes

t = acc_(:,1);
acc = acc_(:,2).*100; % passing to cm/s2

l_zeros=length(acc);
acc=[zeros(1,length(acc)) acc' zeros(1,length(acc))]';
t=[zeros(1,length(t)) t' zeros(1,length(t))]';

acc_in = acc;
npun = length(t);

Nfft = 2^nextpow2(length(t)); 
durata = (Nfft-1).*dt;
df = 1/durata;
freq = 0:df:1/2/dt;
acc(npun+1:Nfft) = 0;
ACC=dt*fft(acc,Nfft); 

mult=ones(length(Sp_in),1);

for i=2:length(T_in) 
Sp_acc(i)=newmark_sa(acc',T_in(i),0.05,0.01);

if T_in(i)<T_sc
    if i>2
        tol_upp_=tol_upp;
        tol_low_=tol_low;
    else
        tol_upp_=tol_upp_pga;
        tol_low_=tol_low_pga;
    end
if Sp_acc(i)>3*Sp_in(i)
    mult(i)=(1+0.70*tol_upp_);
    
else
    if Sp_acc(i)<0.333*Sp_in(i)
        mult(i)=(1-0.70*tol_low_);
        
    else
       logr=log10(Sp_acc(i)/Sp_in(i));
       mult(i)=((1-0.70*tol_low_)+((1+0.70*tol_upp_)-(1-0.70*tol_low_))*...
           (logr-log10(0.333))/(log10(3)-log10(0.333)));
       
    end
end

end

mult(i)=1; % no variability coeff.
Sp_in(i)=mult(i)*Sp_in(i);

freq_in(length(T_in)-i+2)=1/T_in(i);

Rapp_spe_in(length(T_in)-i+2)=Sp_acc(i)/Sp_in(i); % fattore di correzione

end

Sp_acc(1) = max(abs(acc(:)));
freq_in(length(T_in)+1)=1/2/dt; % (*)
Rapp_spe_in(length(T_in)+1)=Sp_acc(1)/Sp_in(1);
freq_in(1)=0;
Rapp_spe_in(1)=1;

Rapp_spe = interp1(freq_in,Rapp_spe_in,freq,'linear'); % interpola il rapp alle freq dell'acc a correggere
acc_pro=acc; % si tiene l'acc non corretto!
nfreq=length(freq);


for k=1:niter % iterates niter times to enhance the spectral agreement
	ACC_PRO=fft(acc_pro)*dt; % FFT of the signal to be corrected
    
	for m=1:nfreq
 		if (freq(m)>=1/T_corr_fin) && (freq(m)<=1/T_corr_ini) 
		ACC_PRO(m)=ACC_PRO(m)/Rapp_spe(m); % correzione della trasformata dell'acc
        end
    end
    
	ACC_PRO(nfreq+1)=0;
    
	for m=nfreq+2:2*nfreq
		ACC_PRO(m)=conj(ACC_PRO(2*nfreq-m+2));
    end;
    
    acc_pro=real(ifft(ACC_PRO))/dt; 
    acc_pro=detrend(acc_pro,'constant');

	for i=1:length(T_in)
    Sp_acc_pro(i,k)=newmark_sa(acc_pro',T_in(i),0.05,0.01);
    end
	
	Sp_acc_pro(1,k)=max(abs(acc_pro));
    
        for i=1:length(T_in)
            rat=Sp_acc_pro(i,k)/Sp_in(i);

            Rapp_spe_pro(length(T_in)-i+2)=rat; % aggiorna il rapporto di correzione
                  
        end
        
        Rapp_spe_pro(1)=1;
        Rapp_spe = interp1(freq_in,Rapp_spe_pro,freq,'linear');
    
end % end of iterations

Sp_in=Sp_in2; clear Sp_in2;

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
t = [0:dt:(length(acc_pro)-1)*dt];

vel_p=cumsum(acc_pro).*dt;
vel_in=cumsum(acc_in).*dt;
dis_p=cumsum(vel_p).*dt;
dis_in=cumsum(vel_in).*dt;

% correction of initial acausal oscillations caused after SPM
[t_idx,idx,Ain] = arias_intensity(acc_in,dt,0.005);% at Ia=0.5%

r = [(t_idx-100*dt):dt:(t_idx+100*dt)]';
a1 = acc_in((idx-100):(idx+100));
a2 = acc_pro((idx-100):(idx+100));

% doing the same thing for velocities
v1 = vel_in((idx-100):(idx+100));
v2 = vel_p((idx-100):(idx+100));

for j = 1:length(r) 
    d(j) = abs((v1(j)-v2(j))/v1(j));
end
m=find(d==min(d));
t1=r(m);
s1 = length(0:dt:t1);
acc_pro(1:s1)=acc_in(1:s1);

% Tapering is used to make sure that record ceases at the end of the time
% vector. Note that tukey with 20% (10% initial, 10% end) is usually a good
% compromise, and does not miss any significant motion
acc_pro = tukeywin(length(acc_pro),0.20).*acc_pro;

pga_corr = 4;
if strcmpi(comp,'ud')
    pga_corr=2;
end

for times=1:pga_corr  
% PGA correction must not be repeated many times. Default=2 times
               
[pga ipga] = max(abs(acc_pro));

dt_cor = 0.03;                                                                        
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

% taking care of potential velocity problems
flag_rev=0;

for n=ipga+(npun_cor+1)/2+1:ipga+(npun_cor+1)/2+100
if (sign(acc_pro(n))/sign(acc_pro(ipga))==-1) && ...
        (sign(acc_pro(n+1))/sign(acc_pro(ipga))==-1) && ...
        (sign(acc_pro(n-1))/sign(acc_pro(ipga))==-1) && ...
        abs(acc_pro(n)) > abs(acc_pro(n-1)) && ...
        abs(acc_pro(n)) > abs(acc_pro(n+1)) && ...
        flag_rev==0
   
    flag_rev=1;
    
    acc_pro(n-(npun_cor-1)/2:n+(npun_cor-1)/2) = ...
        acc_pro(n-(npun_cor-1)/2:n+(npun_cor-1)/2)-yi_d;
    
end
end

if flag_rev==0
    disp('Warning: PGA is only corrected, residual velocity problems may occur')
    %Comment: This warning advises user to expect residual velocity 
    %         problems, since the subtraction operation could not 
    %         be performed 
end
    


end

vel_pro = cumsum(acc_pro)*dt;

% correcting the increasing displacements caused by residual velocity
% problem after the removal of initial oscillations and PGA arrangement
% steps
[t_idx2,idx2,Ain2] = arias_intensity(acc_pro,dt,0.01);%

% calculating the residual velocity from the last 500 samples (last 5
% seconds). Note that this value may be customized.
vel_shift=mean(vel_pro(length(vel_pro)-499:length(vel_pro)));
vel_corr=-vel_shift;

% now attributing to the acceleration point at selected Arias Intensity
% (Ia=1%)
acc_corr=vel_corr/dt;
[pga ipga] = max(abs(acc_pro));

if (abs(acc_corr)> 0.2*pga)
    disp('WARNING! Potentially significant spurious acceleration is added at');
    t_idx2
    % This warning advises user to have a look at the
    % acceleration record at for points between t=[t_idx2+dt t_idx2+4dt]. 
    % Ideally even if the correction is 0.2pga, the change is distributed into
    % 4 points. Most possibly, the output will be ok, if there were not
    % significant velocity issues. User is adviced to provide an engineering
    % judgment if this message appears.
else
    disp('Info: Residual velocity is successfully corrected.');
end

acc_pro(idx2+1)=acc_pro(idx2+1)+acc_corr/4;
acc_pro(idx2+2)=acc_pro(idx2+2)+acc_corr/4;
acc_pro(idx2+3)=acc_pro(idx2+3)+acc_corr/4;
acc_pro(idx2+4)=acc_pro(idx2+4)+acc_corr/4;

clear vel_pro

vel_pro = cumsum(acc_pro)*dt;
dis_pro = cumsum(vel_pro)*dt;

n_=0;

Sp_acc_pro_final(1)=max(abs(acc_pro));

for i=1:length(T_in)
    if i>1
        Sp_acc_pro_final(i)=newmark_sa(acc_pro,T_in(i),0.05,dt);
    end
        rat_real=Sp_acc_pro_final(i)/Sp_in(i);
        if T_in(i)<=0.75
            if i>1
                if rat_real<1+tol_upp && rat_real>1-tol_low
                    if i>2
                    n_=n_+1; % convergence criterion
                    end
                end
            else
                if rat_real<1+tol_upp_pga && rat_real>1-tol_low_pga
                    n_=n_+1; % convergence criterion
                else
                    n_=n_-1; % penalty condition
                end
            end
        else
            if rat_real<1.10 && rat_real>0.90 % T>T*, a fixed 10% tolerance 
                                              % is used.
            n_=n_+1; % convergence criterion
            end
        end
end
 	  
    T_in(1) = 0;
    out_T(1)=0;
    out_T(2:length(Sp_acc_pro_final)-1)=T_in(3:length(Sp_acc_pro_final))';
	out_Se(1)=Sp_acc_pro_final(1);
    out_Se(2:length(Sp_acc_pro_final)-1) = Sp_acc_pro_final(3:length(Sp_acc_pro_final))';
    out_t = [0:dt:(length(acc_pro)-1)*dt];
	out_acc = acc_pro;
    out_vel = vel_pro;
    out_dis = dis_pro;
    ACC_upd=dt*fft(acc_pro,Nfft);
    out_FAS = abs(ACC_upd(1:length(freq)));
    out_freq = freq;
end