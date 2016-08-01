%% *Spectral Scaling*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _spectral_scaling_: function to scale input response spectra so to
% matcht the target spectra
%% INPUT:
% * _inp_dtm (input signal time-step)_
% * _inp_tha (input accelerogram)_
% * _tar_psa (target psa response spectrum)_
% * _tar_vTn (target natural periods)_
%% OUTPUT:
% * _obj_dtm (matched record time-step)_
% * _obj_tha (matched accelerogram)_
% * _obj_thv (matched velocigram)_
% * _obj_thd (matched displacement record)_
% * _obj_vfr (matched frequencies)_
% * _obj_psa (matched psa response sepctrum)_
%% N.B.:
% Need for _band_pass_filter.m_
function [varargout] = spectral_scaling(varargin)
    %% *SET-UP*
    % _input/target spectra_
    inp_dtm = varargin{1};      % input time step
    inp_tha = varargin{2};      % input acceleration time history
    tar_psa = varargin{3};      % target psa
    tar_vTn = varargin{4};      % target natural periods
    tar_nT = numel(tar_vTn);   % number of target natural periods
    % _selected period_
    tar_vTn_idx = varargin{5};
    %
    % _spectral matching set up_
    %
    fac = 1;                    % resampling factor
    scl = 1;                    % scaling factor
    lfr = [];                  % corner frequency
    hfr = [];                   % cutoff frequency
    nit = 10;                   % number of iterations
    pga_target = tar_psa(1);
    %     tol_upp = 0.3; % (1+tol_upp) is the upper spec. tol.
    %     tol_low = 0.1; % (1-tol_low) is the lower spec. tol.
    %
    % _resampling and correcting_
    %
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);
    obj_tha = detrend(obj_tha, 'linear');
    obj_vtm = obj_dtm*(0:obj_ntm-1);
    %
    % _frequency response_
    %
    inp_nfr = 2^nextpow2(obj_ntm);
    inp_dfr = 1/inp_dtm/(inp_nfr-1);
    inp_vfr = (0:inp_dfr:0.5/inp_dtm)';
    obj_nfs = numel(inp_vfr);
    %
    % _refining tar_vTn vector
    % _objective frequency_
    %
    np_cor  = 10;
    idx_save_org(:,1) = find(tar_vTn>0.05);
    idx_save_new = idx_save_org+(np_cor-idx_save_org(1)+1);
    dTn_cor = (log10(0.05)-log10(2*inp_dtm))/(np_cor-1);
    tar_vTn_new(1:np_cor,1) = (log10(2*inp_dtm):dTn_cor:log10(0.05))';
    tar_vTn_new(1:np_cor,1) = 10.^tar_vTn_new(1:np_cor,1);
    tar_vTn_new(idx_save_new,1)=tar_vTn(idx_save_org,1);
    idx0 = find(tar_vTn == 0,1,'first');  % indexes
    if ~isempty(idx0)
        tar_vTn_new(1,1)=0.0;
        obj_vfr(tar_nT-idx0+1) = 0.5/obj_dtm; % Nyquist frequency
    else
        tar_vTn_new(1,1)=2*obj_dtm;
        tar_vTn_new = [0;tar_vTn_new];
    end
    tar_vTn_new = tar_vTn_new(tar_vTn_new<=5);
    idx_min = tar_vTn_new < min(tar_vTn);
    idx_max = tar_vTn_new > max(tar_vTn);
    if ~isempty(idx_min)
        tar_vTn_new(idx_min)=[];
    end
    if ~isempty(idx_max)
        tar_vTn_new(idx_max)=[];
    end
%     tar_psa_new = interp1(tar_vTn,tar_psa,tar_vTn_new,'linear');
    
    %     tar_psa = tar_psa_new;
    %     tar_vTn = tar_vTn_new;
    obj_vfr = [0;flip(1./tar_vTn(2:end));0.5/obj_dtm];
    obj_vTn = tar_vTn(:);
    vfr_cor = [1/tar_vTn(tar_vTn_idx(end));0.5/obj_dtm];
    vfr_cor = (inp_vfr>=vfr_cor(1) & inp_vfr<=vfr_cor(end));
    vfr_hfc = inp_vfr>=40;
    %% *SPECTRAL MATCHING*
    obj_psa = tar_psa;
    for i_ = 1:nit % spectral matching iterations
        %
        % _psa response spectra at target periods_
        %
        [obj_psd(tar_vTn_idx),~,~,obj_psa(tar_vTn_idx),~] = ...
            newmark_sd(obj_tha(:,1),obj_dtm,tar_vTn(tar_vTn_idx),0.05);
        %
        % _psa response spectral ratio at target periods_
        %
        obj_rra = flip(obj_psa(:,1)./tar_psa(:,1));
        %
        % _constrain 0-frequency_
        %
        obj_rra = [1;obj_rra];
        %
        % _linear interpolation_
        %
        obj_rra = interp1(obj_vfr(:,1),obj_rra(:,1),inp_vfr(:,1),'linear');
        %
        % _constrain high frequencies_
        %
        obj_rra(vfr_hfc) = obj_psa(1,1)./tar_psa(1,1);
        %
        % _fourier spectra_
        %
        obj_fsa = super_fft(obj_dtm,obj_tha(:),0,4);
        obj_fsa(vfr_cor) = obj_fsa(vfr_cor)./obj_rra(vfr_cor);
        %
        % _new accelerogram_
        %
        obj_tha = detrend(super_ifft(obj_dtm,obj_ntm,obj_fsa));
    end
    %
    % _correct PGA on time history_
    %
    try
        [pga,ipga] = max(abs(obj_tha));
        dt_cor = 0.05;
        npun_cor = round(dt_cor./(obj_vtm(2)-obj_vtm(1)));
        if mod(npun_cor,2)==0
            npun_cor = npun_cor+1;
        end
        x = [ obj_vtm(ipga-(npun_cor-1)/2),obj_vtm(ipga),obj_vtm(ipga+(npun_cor-1)/2)];
        y = [obj_tha(ipga-(npun_cor-1)/2),pga_target.*sign(obj_tha(ipga)),...
            obj_tha(ipga+(npun_cor-1)/2)];
        xi = obj_vtm(ipga-npun_cor/2:ipga+npun_cor/2-1);
        yi_c = interp1(x,y,xi,'pchip');
        obj_tha(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
            yi_c;
        if (max(abs(obj_tha))-pga_target)>1e-3
            % repeat the same operation as before
            [pga,ipga] = max(abs(obj_tha));
            x = [ obj_vtm(ipga-(npun_cor-1)/2),obj_vtm(ipga),obj_vtm(ipga+(npun_cor-1)/2)];
            y = [obj_tha(ipga-(npun_cor-1)/2),pga_target.*sign(obj_tha(ipga)),...
                obj_tha(ipga+(npun_cor-1)/2)];
            xi = obj_vtm(ipga-npun_cor/2:ipga+npun_cor/2-1);
            yi_c = interp1(x,y,xi,'pchip');
            
            obj_tha(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
                yi_c;
        end
        
        % check final PGA
        if abs((max(abs(obj_tha))-pga_target)/(pga_target))>5e-2
            disp('Check PGA adjustment!');
        end
    catch
    end
    [obj_tha,obj_thv,obj_thd] = band_pass_filter(obj_dtm,obj_tha,lfr,hfr);
    %% OUTPUT
    varargout{1} = obj_dtm;
    varargout{2} = obj_tha;
    varargout{3} = obj_thv;
    varargout{4} = obj_thd;
    varargout{5} = obj_vfr;
    varargout{6} = obj_vTn;
    varargout{7} = obj_psa;
    varargout{8} = obj_psd;
    return
end
