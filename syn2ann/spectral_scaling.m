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
    tar_vTn = varargin{3};      % target natural periods
    tar_psa = varargin{4};      % target psa
    
    % _selected period_
    out_idx = varargin{5};
    %
    % _spectral matching set up_
    %
    fac = 1;                    % resampling factor
    scl = 1;                    % scaling factor
    ni = 10;                     % number of iterations
    pga_target = tar_psa(1);
    tol_upp = 0.1; % (1+tol_upp) is the upper spec. tol.
    tol_low = 0.1; % (1-tol_low) is the lower spec. tol.
    %
    % _resampling and correcting_
    %
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);
    %
    % _frequency response_
    %
    ig=0;
    while (2^ig<=obj_ntm);
        ig=ig+1;
    end
    inp_nfr = 2^ig;
    inp_dfr = 1/(inp_nfr-1)/obj_dtm;
    inp_vfr = (0:inp_dfr:0.5/obj_dtm)';
    obj_tha(obj_ntm+1:inp_nfr) = 0;
    obj_ntm = numel(obj_tha);
    obj_vtm = (0:obj_ntm-1)*obj_dtm;
    inp_fsa = fft(obj_tha);
    %
    % _refining tar_vTn vector at high frequency_
    %
    vTn_cor = 0.05;
    nc  = 10;
    idxn = find(tar_vTn(tar_vTn~=0) < vTn_cor); % indexes
    
    if logical(isempty(idxn))
        idx_save_org(:,1)    = find(tar_vTn>vTn_cor);
        idx_save_new         = idx_save_org+(nc-idx_save_org(1)+1);
        
        dTn_cor              = (log10(vTn_cor)-log10(0.01))/(nc-1);
        tar_vTn_new(1:nc,1)  = (log10(0.01):dTn_cor:log10(vTn_cor))';
        tar_vTn_new(1:nc,1)  = 10.^tar_vTn_new(1:nc,1);
        tar_vTn_new(idx_save_new,1) = tar_vTn(idx_save_org,1);
        tar_vTn_new          = [0;tar_vTn_new];
        tar_psa_new          = interp1(tar_vTn,tar_psa,tar_vTn_new,'linear');
        
        out_idx_new=-999*ones(numel(out_idx),1);
        for j_=1:numel(out_idx)
            out_idx_new(j_) = find(abs(tar_vTn_new-tar_vTn(out_idx(j_)))<1e-8);
        end
        out_idx = out_idx_new;
        tar_vTn = tar_vTn_new;
        tar_psa = tar_psa_new;
    end

    obj_vfr = [0;flip(1./tar_vTn(2:end));0.5/obj_dtm];
    obj_vTn = tar_vTn(:);
    vfr_cor = [1/tar_vTn(out_idx(end)+1);1/0.01];
    vfr_cor = (inp_vfr>vfr_cor(1) & inp_vfr<=vfr_cor(end));
    keyboard
    vfr_hfc = inp_vfr>=50;
    
    %% *SPECTRAL MATCHING*
    obj_psa = SDOF_response(obj_tha(:,1),obj_dtm,tar_vTn,0.05,1);
    obj_rra = flip(obj_psa(:,1)./tar_psa(:,1));
    obj_rra = [1;obj_rra];
    obj_rra = interp1(obj_vfr(:,1),obj_rra(:,1),inp_vfr(:,1),'linear');
    
    for i_ = 1:ni % spectral matching iterations
        %
        % _fourier spectrum_
        %
        obj_fsa          = fft(obj_tha(:))*obj_dtm;
        obj_fsa(vfr_cor) = obj_fsa(vfr_cor)./obj_rra(vfr_cor);
        %
        % _inverse fourier spectrum_
        %
        obj_tha          = super_ifft(obj_dtm,obj_ntm,obj_fsa);
        obj_tha          = detrend(obj_tha);
        %
        % _new psa_
        %
        obj_psa      = SDOF_response(obj_tha(:,1),obj_dtm,tar_vTn,0.05,1);
        obj_rra      = flip(obj_psa(:,1)./tar_psa(:,1));
        obj_rra      = [1;obj_rra];
        obj_rra      = interp1(obj_vfr(:,1),obj_rra(:,1),inp_vfr(:,1),'linear');

    end
    
    %
    % _correct PGA on time history_
    %
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
        [pga,ipga] = max(abs(obj_tha));
        x = [ obj_vtm(ipga-(npun_cor-1)/2),obj_vtm(ipga),obj_vtm(ipga+(npun_cor-1)/2)];
        y = [obj_tha(ipga-(npun_cor-1)/2),pga_target.*sign(obj_tha(ipga)),...
            obj_tha(ipga+(npun_cor-1)/2)];
        xi = obj_vtm(ipga-npun_cor/2:ipga+npun_cor/2-1);
        yi_c = interp1(x,y,xi,'pchip');
        obj_tha(ipga-(npun_cor-1)/2:ipga+(npun_cor-1)/2) = ...
            yi_c;
    end
    
%     check final PGA
    if abs((max(abs(obj_tha))-pga_target)/(pga_target))>5e-2
        disp('Check PGA adjustment!');
    end
    [obj_tha,obj_thv,obj_thd] = integr_diff_avd(obj_dtm,obj_tha);
    [obj_psa,obj_psd] = SDOF_response(obj_tha(:,1),obj_dtm,tar_vTn,0.05,[1,2]);
    %% *OUTPUT*
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
