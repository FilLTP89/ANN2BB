%%%%%%%%%%%%%%%%%%%%
%Funzione per accelerogrammi spettro-compatibili da accelerogrammi reali%
%%%%%%%%%%%%%%%%%%%%

function  [varargout] = spectral_scaling_corr(varargin)
    dt  = varargin{1};
    tha = varargin{2};
    tar_psa = varargin{3};
    tol_upp=0.3;
    tol_low=0.1;
    
    niter=10; % number of iterations
    
    T_sd=tar_psa(:,1);
    T005=find(T_sd==0.05);
    T_vec=[];
    npt1=10; % pti tra 0 e 0.05
    passo1=(log10(0.05)-log10(0.01))/(npt1-1);
    T_vec(1:npt1)=log10(0.01):passo1:log10(0.05);
    T_vec(1:npt1)=10.^T_vec(1:npt1);
    
    if tar_psa(1,1)==0
        T_vec(1)=0.0;  %se ho T=0 sec nel target_Se non aggiungo il T=0.01 sec se no da problemi in interp1 (riga 300)
        npt2=length(T_sd)-T005;
        T_vec(npt1+1:npt1+npt2)=T_sd(T005+1:length(T_sd));
    else
        T_vec(1)=0.01;
        npt2=length(T_sd)-T005;
        T_vec(npt1+1:npt1+npt2)=T_sd(T005+1:length(T_sd));
        
        for i=length(T_vec):-1:1
            T_vec(i+1)=T_vec(i); %infittisco punti prima di T=0.04 s
        end
    end
    
    T_vec(1)=0;
    T4=find((T_vec)==4);
    T_vec_acc=T_vec(1:T4); % per T_vecc_acc mi fermo a 4 s (valori effettivi con cui correggo!!)
    
    %% vettore di periodi da usare per grafici output
    T_vec_out=tar_psa(:,1);
    T4=find(T_vec_out==4);
    T_vec_out_acc=T_vec_out(1:T4);  %T_vec_out_acc=T_vec_out(3:T4);
    %
    a=[0.01:0.005:0.04];
    b=length(a);
    for i=length(T_vec_out_acc):-1:2
        T_vec_out_acc(i+b-2)=T_vec_out_acc(i); %infittisco punti prima di T=0.04 s
    end
    T_vec_out_acc(1:b)=a;
    
    %
    Tmin=min(T_vec_out); % per grafici in logaritmo...
    if Tmin==0
        T0=find(T_vec_out==0);
        Tmin=T_vec_out(T0+1);
    end
    %% Caricamento dati%
    ipadd=1; % flag per padding iniziale
    ltg_pad=2000;
    %
    T_corr_ini=0.01;
    T_corr_fin=0.8;
    %
    target=tar_psa;
    pga_target = tar_psa(1,2); %g
    %
    Tmin=min(target(:,1));
    vmin=T_vec_acc < Tmin; a1=sum(vmin);
    Tmax=max(target(:,1));
    vmag=T_vec_acc > Tmax; a2=sum(vmag);
    
    if a2 ~= 0
        T_vec_acc([length(T_vec_acc)-a2+1:length(T_vec_acc)])=[];
    end
    
    if a1 ~= 0
        T_vec_acc([1:a1])=[];
    end
    %
    %interpolo target in periodi T_vec_acc (= T correzione)
    for tt=1:length(T_vec_acc)
        target_int(tt)=interp1(target(:,1),target(:,2),T_vec_acc(tt));
        %interpolation(target(:,1),target(:,2),T_vec_acc(tt));
    end
    
    T_in=T_vec_acc'; % periodi a cui si lavora in accelerazione!!
    Sp_in=target_int; % in m/s2...

    %
    %  fa eventualmente una decimazione del segnale
    dt1=dt;
    fac=1;
    npun1=length(tha);
    npun2=fac*length(tha);
    t1=[0:dt1:(npun1-1)*dt1];
    dt=dt1/fac;
    t=[0:dt:(npun2-1)*dt];
    acc = interp1(t1,tha,t,'spline'); % con spline calcola per ciascun valore...non mai d� NaN
    npun=length(acc);
    t=[0:dt:(npun-1)*dt];
    
    %Butterworth acausal filter parameter
    Fc=0.1; Fh=50; Nbut=4;
    % Wn=[Fc Fh];
    Wn=[Fc];
    fnyq=1/2/dt;
    ftype = 'high'; % bandpass  high
    [BuB,BuA]=butter(Nbut,Wn./fnyq,ftype);
    
    % End input parameters
    
    % baseline correction, taper and filter
    deriva=mean(acc(1:npun));
    acc=acc-deriva;
    w=tukeywin(npun,0.05);
    acc=acc.*w';
    if ipadd==1
        acc=padarray(acc,[0 ltg_pad],'both');
    end
    acc=filtfilt(BuB,BuA,acc);
    amax=max(abs(acc));
    npun=length(acc);
    t=[0:dt:(npun-1)*dt];
    
    %Calcolo della FFT dell'accelerogramma %
    ig=0;
    while (2^ig<=npun);
        ig=ig+1;
    end
    Nfft=2^ig; %numero di punti su cui viene calcolata la FFT%
    durata=(Nfft-1).*dt; %durata del segnale che include anche i punti introdotti per il calcolo della FFT%
    df=1/durata;
    freq=0:df:1/2/dt;
    
    acc(npun+1:Nfft)=0;
    vel=cumsum(acc).*dt;
    deriva=mean(vel);  % mv
    vel=vel-deriva;    % mv
    vel=filtfilt(BuB,BuA,vel);  % mv
    dis=cumsum(vel).*dt;
    npun=Nfft;
    t=[0:dt:(npun-1)*dt];
    ACC=dt*fft(acc,Nfft); %trasformata di Fourier dell'accelerogramma non corretto%
    %Fine calcolo della FFT dell'accelerogramma%
    
    ACC=fft(acc);  % trasformata di Fourier non corretta per dt..
    
    %  spettro acc non corretto calcolato ai periodi del target
    [~,~,~,Sp_acc,~]=newmark_sd(acc',dt,T_in,0.05);
    Sp_acc(1)=max(abs(acc));
    %
    for i=2:length(T_in), % parto da due per evitare di prendere T_in(1)=0 e quindi un NaN..per�..(*)
        freq_in(length(T_in)-i+2)=1/T_in(i);
        Rapp_spe_in(length(T_in)-i+2)=Sp_acc(i)/Sp_in(i); % fattore di correzione
    end
    freq_in(length(T_in)+1)=1/2/dt; % (*)...questo funziona solo se il la fNyq � <= di 1/T_in(1)..
    Rapp_spe_in(length(T_in)+1)=Sp_acc(1)/Sp_in(1);
    freq_in(1)=0;
    Rapp_spe_in(1)=1;
    
    Rapp_spe = interp1(freq_in,Rapp_spe_in,freq,'linear'); % interpola il rapp alle freq dell'acc a correggere
    
    acc_pro=acc; % si tiene l'acc non corretto!
    nfreq=length(freq);
    %
    %     for k=1:niter(jn),  % niter diversi per acc....
    for k=1:niter,  % lavora tante volte sull'accelerogramma (niter = per tutti)
        ACC_PRO=fft(acc_pro); % trasformata dell'acc da correggere
        for m=1:nfreq,
            if (freq(m)>=1/T_corr_fin) && (freq(m)<=1/T_corr_ini) % cambiato qui
                ACC_PRO(m)=ACC_PRO(m)/Rapp_spe(m); % correzione della trasformata dell'acc
            end % cambiato qui
        end;
        ACC_PRO(nfreq+1)=0;
        for m=nfreq+2:2*nfreq,
            ACC_PRO(m)=conj(ACC_PRO(2*nfreq-m+2));
        end;
        
        acc_pro=real(ifft(ACC_PRO)); % antitraforma e aggiorna l'acc da correggere
        acc_pro=detrend(acc_pro,'constant');
        
        [Sp_dis_pro(:,k),~,~,Sp_acc_pro(:,k),~]=newmark_sd(acc_pro',dt,T_in,0.05);
        Sp_acc_pro(1,k)=max(abs(acc_pro));
        
        for i=1:length(T_in),
            %%% modified on 16/06/16 by AGO
            rat=Sp_acc_pro(i,k)/Sp_in(i);
            
            if rat>1 && rat <= 1+tol_upp || rat < 1 && rat >= 1-tol_low
                rat=1;
            end
            Rapp_spe_pro(length(T_in)-i+2)=rat; % aggiorna il rapporto di correzione
        end
        
        Rapp_spe_pro(1)=1;
        Rapp_spe = interp1(freq_in,Rapp_spe_pro,freq,'linear'); % lo interpola e aggiorna il rapp_spe
        
    end % fine iterazioni
    
    % CS 03.05.2016 + 23.06.2016
    [pga,ipga] = max(abs(acc_pro));
    dt_cor = 0.05;
    npun_cor = round(dt_cor./(t(2)-t(1)));
    if mod(npun_cor,2)==0
        npun_cor = npun_cor+1;
    end
    x = [ t(ipga-(npun_cor-1)/2),t(ipga),t(ipga+(npun_cor-1)/2)];
    y = [acc_pro(ipga-(npun_cor-1)/2),pga_target.*sign(acc_pro(ipga)),...
        acc_pro(ipga+(npun_cor-1)/2)];
    xi = t(ipga-npun_cor/2:ipga+npun_cor/2-1);
    yi_c = interp1(x,y,xi,'cubic');
    acc_pro(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
        yi_c;
    if (max(abs(acc_pro))-pga_target)>1e-3
        % repeat the same operation as before
        [pga,ipga] = max(abs(acc_pro));
        x = [ t(ipga-(npun_cor-1)/2),t(ipga),t(ipga+(npun_cor-1)/2)];
        y = [acc_pro(ipga-(npun_cor-1)/2),pga_target.*sign(acc_pro(ipga)),...
            acc_pro(ipga+(npun_cor-1)/2)];
        xi = t(ipga-npun_cor/2:ipga+npun_cor/2-1);
        yi_c = interp1(x,y,xi,'cubic');
        %      yi_p = interp1(x,y,xi,'pchip');
        %      yi_l = interp1(x,y,xi,'linear');
        
        acc_pro(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
            yi_c;        
    end

    % CS 23.06.2016
    % check final PGA
    if abs((max(abs(acc_pro))-pga_target)/(pga_target))>5e-2
        disp('Check PGA adjustment!');
        stop
    end
    % Band-pass acausal Butterworth filter
    acc_pro=filtfilt(BuB,BuA,acc_pro);
    vel_pro=cumsum(acc_pro)*dt;
    vel_pro=detrend(vel_pro,'linear');
    dis_pro=cumsum(vel_pro)*dt;
    dis_pro=detrend(dis_pro,'linear');
    frac=10; %% taper
    dis_pro = taper_fun(dis_pro,frac,1,1);
    
    vel_pro=diff(dis_pro)/dt;
    vel_pro(length(t))=0;
    acc_pro=diff(vel_pro)/dt;
    acc_pro(length(t))=0;
    %% OUTPUT
    varargout{1}=dt;
    varargout{2}=acc_pro;
    varargout{3}=vel_pro;
    varargout{4}=dis_pro;
    varargout{5}=T_in;
    varargout{6}=Sp_acc_pro(:,end);
    varargout{7}=Sp_dis_pro(:,end);
    varargout{8}=freq_in;
    return
end
