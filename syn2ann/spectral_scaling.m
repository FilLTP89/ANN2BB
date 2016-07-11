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
    pga_target=tar_psa(1);
    % _selected period_
    vTn_corner = varargin{5};
    %
    % _spectral matching set up_
    %
    fac = 1;                    % resampling factor
    scl = 1;                    % scaling factor
    lfr = 0.1;                   % corner frequency
    hfr = [];                   % cutoff frequency
    nit = 10;                    % number of iterations
    tol_upp=0.3; % (1+tol_upp) is the upper spec. tol.
    tol_low=0.1; % (1-tol_low) is the lower spec. tol.
    %
    % _resampling_
    %
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);
    obj_vtm = obj_dtm*(0:obj_ntm-1);
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
    tar_vTn_new = tar_vTn_new(tar_vTn_new<=4);
    idx_min = tar_vTn_new < min(tar_vTn);
    idx_max = tar_vTn_new > max(tar_vTn);
    if ~isempty(idx_min)
        tar_vTn_new(idx_min)=[];
    end
     if ~isempty(idx_max)
        tar_vTn_new(idx_max)=[];
    end
    tar_psa_new = interp1(tar_vTn,tar_psa,tar_vTn_new,'linear');

    
    if nargin>5
        sp = varargin{6};
        col = jet(nit+1);
        xpl = cell(nit+1,1);
        xpl{nit+1} = tar_vTn;
        ypl = cell(nit+1,1);
        ypl{nit+1} = tar_psa./9.81;
        set(0,'defaultaxescolororder',jet(nit+1));
        mrk = cell(nit+1,1);
        lst = cell(nit+1,1);
        mrk(nit+1) = {'p'};
        mrk(1:end-1) = {'none'};
        lst(end) = {'none'};
        lst(1:end-1) = {'-'};
        leg = cell(nit+1,1);
        leg(end) = {'TARGET'};
    end
    
    tar_psa = tar_psa_new;
    tar_vTn = tar_vTn_new;
    obj_vfr = [0;flip(1./tar_vTn(2:end));0.5/obj_dtm];
    obj_vTn = tar_vTn(:);
    vfr_cor = [1/vTn_corner,min(1/tar_vTn(2),0.5/obj_dtm)];
    %% *SPECTRAL MATCHING*
    for i_ = 1:nit % spectral matching iterations
        %
        % _psa response spectra at target periods_
        %
        [obj_psd,~,~,obj_psa,~] = ...
            newmark_sd(obj_tha(:,1),obj_dtm,tar_vTn,0.05);
        if nargin>5
            xpl{i_} = tar_vTn;
            ypl{i_} = obj_psa./9.81;
            leg(i_) = {sprintf('it-%u',i_)};
        end
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
        % _fourier spectra_
        %
        [~,~,~,obj_fsa,~,~] = super_fft(obj_dtm,obj_tha(:),0);
        obj_fsa=obj_fsa(1:obj_nfs)./obj_dtm;
        for j_ = 1:obj_nfs
            if inp_vfr(j_)>=vfr_cor(1)&&inp_vfr(j_)<=vfr_cor(2)
                obj_fsa(j_)=obj_fsa(j_)./obj_rra(j_);
            end
        end
        obj_fsa(obj_nfs+1,1)  = 0.;
        obj_fsa(obj_nfs+2:2*obj_nfs,1) = flip(conj(obj_fsa(2:obj_nfs,1)));
        %
        % _new accelerogram_
        %
        obj_tha = (real(ifft(obj_fsa(:))));
        obj_tha = obj_tha(~isnan(obj_tha(1:obj_ntm)));
    end
    if nargin>5
        fpplot('xpl',xpl,'ypl',ypl,'xlb',{'T [s]'},'ylb',{'PSA [g]'},...
            'tit',{'PSEUDO-ACCELERATION SPECTRUM (5%)'},...
            'mrk',mrk,'lst',lst,'xlm',{[0,1]},...
            'xtk',{0:.25:5},'ytk',{0:.25:2},'leg',{leg});
        saveas(gcf,sp,'epsc');
    end
    %
    % _correct PGA on time history_
    % CS 03.05.2016 + 23.06.2016
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
    yi_c = interp1(x,y,xi,'cubic');
    obj_tha(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
        yi_c;
    if (max(abs(obj_tha))-pga_target)>1e-3
        % repeat the same operation as before
        [pga,ipga] = max(abs(obj_tha));
        x = [ obj_vtm(ipga-(npun_cor-1)/2),obj_vtm(ipga),obj_vtm(ipga+(npun_cor-1)/2)];
        y = [obj_tha(ipga-(npun_cor-1)/2),pga_target.*sign(obj_tha(ipga)),...
            obj_tha(ipga+(npun_cor-1)/2)];
        xi = obj_vtm(ipga-npun_cor/2:ipga+npun_cor/2-1);
        yi_c = interp1(x,y,xi,'cubic');
        %      yi_p = interp1(x,y,xi,'pchip');
        %      yi_l = interp1(x,y,xi,'linear');
        
        obj_tha(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
            yi_c;
    end
    
    % CS 23.06.2016
    % check final PGA
    if abs((max(abs(obj_tha))-pga_target)/(pga_target))>5e-2
        disp('Check PGA adjustment!');
    end
    % _filtering
    %
    [obj_tha,obj_thv,obj_thd] = band_pass_filter(obj_dtm,obj_tha,lfr,hfr);
    obj_vfr = obj_vfr(2:end,1);
    
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