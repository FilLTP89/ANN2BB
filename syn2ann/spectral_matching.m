%Programma per accelerogrammi spettro-compatibili da accelerogrammi reali%
function [varargout] = spectral_matching(varargin)
    dt    = varargin{1};
    acc1  = varargin{2};
    T_in  = varargin{3};
    Sp_in = varargin{4};
% clear all;
% close all;
% clc;
% warning off;

% Caricamento dati%

% load 'syn2ann_res_ALL';
% i=1;
% record originale (struttura dati rec.org):
        %rec.org.syn{1}.tha.e (.n,.z) : time-histories in accelerazione (tha) dell' i-esima stazione nella direzione ew (ns,z). 
        %NB: i=1 ==> MRN; i=3 ==> AQK
%         Rec_orig=rec.org.syn{i}.tha.n % : time-histories in accelerazione (tha) dell' i-esima stazione nella direzione ew (ns,z). 

%     * Numerical Simulations (SPEED) originale (struttura dati nss.org): 
%         nss.org.syn{i}.tha.e (.n,.z) : time-histories in accelerazione (tha) dell' i-esima stazione nella direzione ew (ns,z). 
%         NB: i=1 ==> MRN; i=3 ==> AQK
%          Rec_sim=nss.org.syn{i}.tha.n % : time-histories in accelerazione (tha) dell' i-esima stazione nella direzione ew (ns,z). 

%     * Hybrids con SP96 (struttura dati hbs.sps):
%         hbs.sps.syn{i}.tha.e (.n,.z) : time-histories in accelerazione (tha) dell' i-esima stazione nella direzione ew (ns,z). 
%         NB: i=1 ==> MRN; i=3 ==> AQK
%           hybrid=hbs.sps.syn{i}.tha.n % : time-histories in accelerazione (tha) dell' i-esima stazione nella direzione ew (ns,z). 


%     * Spectral Matched con ibrid basati su SP96 (struttura dati spm.sps):
%         spm.sps.e(.n,.z).syn{i}.psa.e(.n,.z) : pseudo-spectral acceleration (psa) dell' i-esima stazione nella direzione ew (ns,z). 
%         NB: i=1 ==> MRN; i=3 ==> AQK
%         NB: la direzione viene ripetuta due volte e deve coincidere : spm.sps.e.syn{i}.psa.e
    
%           Sp_in=spm.sps.n.syn{i}.psa.n %: pseudo-spectral acceleration (psa) dell' i-esima stazione nella direzione ew (ns,z). 
%     * time steps:
%         nss.org.mon.dtm(i) : delta time (dtm) dell' i-esima stazione.
%         NB: i=1 ==> MRN; i=3 ==> AQK
%         NB: chiaramente, questo vale anche per le altre strutture, come rec.org, hbs.sps
%           dt=nss.org.mon.dtm(i); 

%     * time vector:
%         nss.org.mon.vtm{i} : time vector (vtm) dell' i-esima stazione.
%         NB: i=1 ==> MRN; i=3 ==> AQK
%         NB: chiaramente, questo vale anche per le altre strutture, come rec.org, hbs.sps
%           t=nss.org.mon.vtm{i}; 

%     * natural periods
%         spm.sps.e(.n,.z).mon.vTn : natural periods vector (vTn) (lo stesso per tutte le stazioni)
%           T_in=spm.sps.n.mon.vTn;

% nomefile2='Pippo.dat'; %nome del file da salvare%

    scala=1;
    dt1=dt;
    fac=1;
    npun1=length(acc1);
    npun2=fac*length(acc1);
    t1=[0:dt1:(npun1-1)*dt1];
    dt=dt1/fac;
    t=[0:dt:(npun2-1)*dt];
    acc = interp1(t1,acc1,t,'spline');

    acc=acc*scala;
    npun=length(acc); 

    %Butterworth acausal filter parameter
    Fc=0.05; Fh=40; Nbut=2;
    % Wn=[Fc Fh];
    Wn=Fh;
    fnyq=1/2/dt;
    [BuB,BuA]=butter(Nbut,Wn./fnyq,'low');

% End input parameters

    % baseline correction 
    acc=detrend(acc, 'linear');
    acc=filtfilt(BuB,BuA,acc);
    amax=max(abs(acc));

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
    dis=cumsum(vel).*dt;
    npun=Nfft;
    t=[0:dt:(npun-1)*dt];
    ACC=dt*fft(acc,Nfft); %trasformata di Fourier dell'accelerogramma non corretto%
    %Fine calcolo della FFT dell'accelerogramma%

    % Grafici %

    % %Accelerogramma non corretto e relativo Spettro di Fourier%
    % figure(1) 
    %  
    % subplot ('position',[0.1 0.55 0.75 0.25])
    % ax1 = gca;
    % set(ax1,'FontSize',12);
    % plot(t,acc,'k','linewidth',1);
    % title ('Accelerogramma');
    % axis([0 durata -1.2*max(abs(acc)) 1.2*max(abs(acc))]);
    % xlabel('t(s)');
    % ylabel('m/s^2');
    % grid on;
    % 
    % subplot ('position',[0.1 0.15 0.75 0.25])
    % ax1 = gca;
    % set(ax1,'FontSize',12);
    % loglog(freq,abs(ACC(1:Nfft/2)),'k','linewidth',2);
    % title ('Spettro di Fourier');
    % xlabel('f(Hz)');
    % grid on;
    % hold on;
    ACC=fft(acc);

    % loglog(freq,abs(ACC(1:Nfft/2))*dt,'r','linewidth',2);

    % Compute the PSA response spectrum of the input accelerogram
    for i=2:length(T_in),
        Sp_acc(i)=4*pi^2*disp_spectra(acc,dt,T_in(i),0.05)./T_in(i)^2;
    end
    Sp_acc(1)=max(abs(acc));

    % compute the ratio of the inpute PSA vs the target PSA

    for i=2:length(T_in),
        freq_in(length(T_in)-i+2)=1/T_in(i);
        Rapp_spe_in(length(T_in)-i+2)=Sp_acc(i)/Sp_in(i);
    end
    freq_in(length(T_in)+1)=1/2/dt;
    Rapp_spe_in(length(T_in)+1)=Sp_acc(1)/Sp_in(1);
    freq_in(1)=0;
    Rapp_spe_in(1)=1;
    Rapp_spe = interp1(freq_in,Rapp_spe_in,freq,'linear');

    niter=2; % number of iterations
    acc_pro=acc;
    nfreq=length(freq);

    for k=1:niter
        ACC_PRO=fft(acc_pro);
        for m=1:nfreq
            if(k<=niter)
                ACC_PRO(m)=ACC_PRO(m)/Rapp_spe(m);
            else
                if (m*df>5)
                   ACC_PRO(m)=ACC_PRO(m)/Rapp_spe(m); 
                end
            end
        end
        ACC_PRO(nfreq+1)=0;
        for m=nfreq+2:2*nfreq
            ACC_PRO(m)=conj(ACC_PRO(2*nfreq-m+2));
        end

        acc_pro=real(ifft(ACC_PRO));
        acc_pro=detrend(acc_pro,'constant');
    
        % response spectrum of corrected waveform
        for i=2:length(T_in)
            Sp_acc_pro(i)=4*pi^2*disp_spectra(acc_pro,dt,T_in(i),0.05)./T_in(i)^2;
        end
        Sp_acc_pro(1)=max(abs(acc_pro));

        for i=1:length(T_in)
            Rapp_spe_pro(length(T_in)-i+2)=Sp_acc_pro(i)/Sp_in(i);
        end
        Rapp_spe_pro(1)=1;
        Rapp_spe = interp1(freq_in,Rapp_spe_pro,freq,'linear');
        f_amax=40;
        for i=1:length(freq)
            if (freq(i)>f_amax) 
                Rapp_spe(i)=Sp_acc_pro(1)/Sp_in(1);
            end
        end
    end

    % final correction for compatibility
    % Band-pass acausal Butterworth filter
    acc_pro=filtfilt(BuB,BuA,acc_pro);

    vel_pro=cumsum(acc_pro)*dt;
    vel_pro=detrend(vel_pro,'linear');
    dis_pro=cumsum(vel_pro)*dt;
    dis_pro=detrend(dis_pro,'linear');
    frac=5; %% taper
    dis_pro = taper_fun(dis_pro,frac,1,1);

    vel_pro=diff(dis_pro)/dt;vel_pro(length(t))=0;
    acc_pro=diff(vel_pro)/dt;acc_pro(length(t))=0;

    T_out= T_in;%[0.01, 0.03:0.01:3];
    for i=2:length(T_out)
        Sp_acc_pro(i)=4*pi^2*disp_spectra(acc_pro,dt,T_out(i),0.05)./T_out(i)^2;
    end
    Sp_acc_pro(1)=max(abs(acc_pro));

    %% *OUTPUT*
    varargout{1} = dt;
    varargout{2} = acc_pro(:);
    varargout{3} = T_out(:);
    varargout{4} = Sp_acc_pro(:);
    varargout{5} = vel_pro(:);
    varargout{6} = dis_pro(:);
    return
end
%     figure (2);
%     % plot (T_in,Sp_in,'r',T_in,Sp_acc,'k',T_out,Sp_acc_pro,'b',T_out,Sp_acc2,'g')
%     loglog (T_in,Sp_in,'r',T_in,Sp_acc,'k',T_out,Sp_acc_pro,'b'); grid on;
%     xlim([0 5])
% 
%     figure (3)
%     lRec=length(Rec_sim);
%     plot(t(1:lRec),acc(1:lRec),'k','linewidth',1);
%     hold on;
%     plot(t(1:lRec),acc_pro(1:lRec),'b','linewidth',1);
%     hold on;
%     plot(t(1:lRec),Rec_sim,'g','linewidth',1);
%     % hold on;
%     % plot(t(1:lRec),Rec_orig(1:lRec),'m','linewidth',1);
% 
%     figure (4)
%     loglog (freq, abs(ACC(1:nfreq)),'k',freq, abs(ACC_PRO(1:nfreq)),'b')
% 
%     figure (5)
%     subplot (2,1,1), plot (t,vel,'k',t,vel_pro,'b');
%     subplot (2,1,2), plot (t,dis,'k',t,dis_pro,'b');
% 
%     %Salvataggio file con accelerazione corretta%
%     fid = fopen(nomefile2,'w');
%     fprintf(fid,'%4.12f\n',acc_pro); 
%     fclose(fid); 
