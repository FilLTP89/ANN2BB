% ====================== bbsyn.m =========================================%
% the programs builds a broadband (BB) from a Low Frequency (LF) signal and 
% a High Frequency (HF) signal through a summation procedure in the
% frequency domain with suitable weight functions 

clear all
close all
clc
warning off

stat='aqk'
% ================ INPUT PARAMETERS ============================== % 
comp = 'x'; 
% frequency band 
fa = 2.5; % in Hz
fb = 3; % in Hz 

% ================ LF signal: time,comp in m/s^2 ======================== %
path_lf='..\walters\'
nomelfx = 'namlf';
nomelfx = [nomelfx,comp,'_',stat,'.dat'];
datilfx = leggi_dati(nomelfx);
dimdati=length(datilfx);
nlin_mac = 30; 
ncol_mac = 5; 
gray1=[0.75 0.75 0.75];
gray2=[0.5 0.5 0.5];

% ================ HF signal: time,comp in m/s^2 ======================== %
path_hf=['..\HF_exsim\exsim_aq\walters\',stat,'\']
nomehfx = 'namhf'; 
nomehfx = [nomehfx,comp,'_',stat,'.dat'];
datihfx = leggi_dati(nomehfx);
nlin = 12; 


nomebbx = 'nambb';
nomebbx = [nomebbx,comp,'_',stat,'.dat'];
datibbx = leggi_dati(nomebbx);

for ir = 18%1:dimdati 
     
     C=char(datihfx(ir));
     nomefile=C;
     nomefilehfx = [char(path_hf),nomefile]; 
     C=char(datilfx(ir));
     nomefile=C;
     nomefilelfx = [char(path_lf),nomefile];     
     C=char(datibbx(ir));
     nomefile=C;
     outfilebbx = [nomefile];  
     
    % loading HF time histories (acc) 
    letthfx = textread(nomefilehfx,'','headerlines',nlin); 
    thfx = letthfx(:,1); 
    dthfx = thfx(2)-thfx(1); 
    acchfx = letthfx(:,2); 
    acchfx = [acchfx;zeros(length(acchfx),1)]; 
    % Fourier spectrum 
    [SHFX,fhfx]=calcola_spettro(acchfx,dthfx); 
    
    % loading LF time histories (acc) [mac format]
    lettlfx = load(nomefilelfx);
    dtlfx = lettlfx(2,1)-lettlfx(1,1); 
    acclfx = lettlfx(:,2); 
    acclfx = [acclfx;zeros(length(acclfx),1)]; 
    tlfx = lettlfx(:,1);
    % Fourier spectrum 
    [SLFX,flfx]=calcola_spettro(acclfx,dtlfx); 
   
    
    npunlf = length(acclfx); 
    dtlf = dtlfx; 
    npunhf = length(acchfx);
    dthf = dthfx; 
    if dtlf ~= dthf
        if abs(dtlf-dthf)>1e-10
            disp('ERORR: dt_lf ~ dt_hf!!!"'); 
            disp('dt_lf =');
            disp(dtlf); 
            disp('dt_hf =');
            disp(dthf);  
            stop 
        end
    end
    dt = dtlf;
  
    if npunlf ~= npunhf
        disp('WARNING: npun_lf ~ npun_hf!!!"'); 
        disp('npun_lf =');
        disp(npunlf); 
        disp('npun_hf =');
        disp(npunhf);
    end

    npun = npunlf;
    if npunhf < npun
        acchf(1:npunhf) = acchfx(1:npunhf);
        acchf(npunhf+1:npun) = 0.0; 
    elseif npunhf > npun
        acchf = acchfx(1:npun); 
    end
    acclf = acclfx; 
   
    time = [0:dt:(npun-1)*dt]; 

  % shift of signal according to T05 of arias intensity
  I1 = 0.05; 
  [T5lf,i5lf,Ialf] = arias_intensity(acclf,dt,I1);
  [T5hf,i5hf,Iahf] = arias_intensity(acchf,dt,I1);
  
  idiff = i5hf-i5lf
  if idiff ~= 0
      if idiff > 0 
          cont = 1; 
          for i = abs(idiff):npun
              hf(cont) = acchf(i); 
              cont = cont +1; 
          end
          for i = npun+1:npun+abs(idiff)-1
              hf(cont) = 0.0; 
              cont = cont +1; 
          end
          for i = 1:npun
              lf(i) = acclf(i); 
          end
      elseif idiff < 0 
          for i = 1:abs(idiff)
              hf(i) = 0;
          end
          for i = abs(idiff)+1:npun
              hf(i) = acchf(i-abs(idiff)+1); 
          end
          for i = 1:npun
              lf(i) = acclf(i); 
          end   
      end
  else
      lf=acclf;
       hf = acchf;
  end


  y_lim = 400;
  t_0 = 0; 
  t_lim = 25; 
  figure(ir)
  subplot(311)
  ax1 = gca; 
  set(gca,'fontsize',12); 
  plot(time,lf,'Color',gray1,'linewidth',1.5);
  grid on 
  ylim([-y_lim y_lim])
%     xlim([t_0 t_lim]); 
   xlim([0 25]);
  ylabel('\it\bfacc^{LF} \rm\bf[cm/s^2]');
%   title(['Well #',lab_well]);
  subplot(312)
  ax1 = gca; 
  set(gca,'fontsize',12); 
  plot(time,hf,'Color',gray2,'linewidth',1.5);
  ylabel('\it\bfacc^{HF} \rm\bf[cm/s^2]');
  grid on 
  ylim([-y_lim y_lim])
  xlim([0 25]); 


  % ==================== Fast Fourier Transform ======================
  ig=0;
  while (2^ig<=npun); 
      ig=ig+1;
  end
  Nfft=2^ig;
  df=1/((Nfft-1).*dt);
  disp('df =');
  disp(df); 
  lf(2*npun+1:Nfft)=0; % effective duration of lf signal 
  hf(2*npun+1:Nfft)=0; % effective duration of hf signal
  LF = dt.*fft(lf,Nfft); % Fourier spectrum of LF 
  HF = dt.*fft(hf,Nfft); % Fourier spectrum of HF 
  fr = [0:df:(Nfft/2-1)*df];

  figure(100+ir)
  ax1 = gca; 
  set(gca,'fontsize',13); 
  loglog(fr,abs(LF(1:Nfft/2)),'Color',gray1,'linewidth',2.0);
  hold on 
  loglog(fr,abs(HF(1:Nfft/2)),'Color',gray2,'linewidth',2.0);
  hold on 
  xlabel('\it\bff \rm\bf[Hz]'); 
  ylabel('\it\bfFAS \rm\bf[cm/s]'); 

  xlim([0.05 20]);
  ylim([0.1 400]);  
  set(gca,'YTick',[1e-2 1e-1 1.0 10 100]); 
  set(gca,'XTick',[0.1 1 5 10 20]); 
%   title(['Well #',lab_well]);

  %======================= BROADBAND signal ================================
  nfa = round(fa/df); 
  nfb = round(fb/df); 
  fac = pi./(df*(nfb-nfa-2)); 
  
  WLF = zeros(1,Nfft/2); 
  WHF = zeros(1,Nfft/2); 

  for j = 1:Nfft/2
    if j<nfa
        WLF(j) = 1.0; % weight LF part 
    elseif ((j>= nfa)&&(j<=nfb))
        WLF(j) = 0.5+0.5*cos(fac*df*(j-nfa-1)); % weight LF part
    else
        WF(j) = 0; 
    end
    WHF(j) = 1-WLF(j);  % weight HF part
  end


  figure(1000)
  ax1 = gca; 
  set(gca,'fontsize',14); 
  plot(fr,WLF,'Color',gray1,'linewidth',3)
  hold on 
  plot(fr,WHF,'Color',gray2,'linewidth',3)
  f_l = 1.5; 
  f_h = 2.0; 
  xf_l = ones(2,1)*f_l; 
  yf_l = [0 1.05]; 
  xf_h = ones(2,1)*f_h; 
  yf_h = [0 1.05];
  xlim([0 3]); 
  ylim([0 1.1]); 
  hold on 
  plot(xf_l,yf_l,'--k'); 
  hold on 
  plot(xf_h,yf_h,'--k'); 
  grid on 
  set(ax1,'XTick',[0 0.5 1.0 1.5 2.0 2.5 3.0])
  set(ax1,'XTickLabel',{'0'; '0.5';'1.0';'1.5';'2.0';'2.5';'3.0'})
    set(ax1,'YTick',[0 0.5 1.0 ])
  set(ax1,'YTickLabel',{'0'; '0.5';'1.0'})
  
  
  for j=1:Nfft/2
      BB(j) = WLF(j).*LF(j)+ WHF(j).*HF(j);
  end
  BB(Nfft/2+1) = 0.0;
  for j=Nfft/2+2:Nfft
      BB(j) = conj(BB(Nfft-j+2));
  end

  figure(100+ir)
  hold on
  loglog(fr,abs(BB(1:Nfft/2)),'k--','linewidth',1.5); 
  legend('LF(GeoELSE)','HF(Exsim)','BB(hybrid)','Location','SouthWest'); 
%   legend('boxon');
 
  % Inverse Fourier Transform 
  bb_i = ifft(BB(:));
  for j=1:npun
      bb(j) = real(bb_i(j))./dt;
  end


  figure(ir)
  subplot(313)
  ax1 = gca; 
  set(gca,'fontsize',12); 
  plot(time,bb,'k','linewidth',1.5);
  ylabel('\it\bfacc^{BB} \rm\bf[cm/s^2]');
  grid on 
  ylim([-y_lim y_lim])
   xlim([t_0 t_lim]); 
  xlabel('\it\bft \rm\bf[s]'); 
%   fig_name =  [char(path),['figures\acc',comp,'_HF_LH_BB_',lab_well]]; 
%   nomefig = [fig_name,'.fig']; 
%   print(['-f',num2str((ir))], '-dtiff','-r300',fig_name);
%   saveas(gcf,nomefig); 
%   fig_name_eps = [char(path),['figures\acc',comp,'_HF_LH_BB_',lab_well]];
%   print(['-f',num2str((ir))], '-depsc',fig_name);
  
  
  % writing output file 
  fid = fopen(outfilebbx,'w'); 
  fprintf(fid,'%e   %e\n',[time',bb']'); 
  fclose(fid);
  
  fid = fopen(outfilebbx,'w'); 
  fprintf(fid,'%e   %e\n',[time',bb']'); 
  fclose(fid);
  
  four=[fr' abs(LF(1:Nfft/2))' abs(HF(1:Nfft/2))' abs(BB(1:Nfft/2))']
  acc_res=[time' lf' hf' bb']
  filtro=[fr' WLF' WHF']
end 


