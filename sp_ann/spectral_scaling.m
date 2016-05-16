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
    
    %% SET-UP
    %%
    % _input/target spectra_
    inp_dtm = varargin{1};
    inp_tha = varargin{2};
    tar_psa = varargin{3};
    tar_vTn = varargin{4};
    tar_nTn = numel(tar_vTn);
    %%
    % _filtering properties_
    fac = 1;
    scl = 1;
    lfr = 0.2; % corner frequency
    hfr = 40;  % cutoff frequency
    frm = 40;  % max frequency considered
    %%
    % _resampling and processing_
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);
    [obj_tha,~,~] = band_pass_filter(obj_dtm,obj_tha,lfr,hfr);
    %%
    % _spectral matching set up_
    nit = 5;                    % number of iterations
    obj_nfr = tar_nTn+2;
    obj_vfr = -ones(obj_nfr,1);
    obj_vfr(end)     = 0.5/obj_dtm;
    obj_vfr(2:end-1) = flip(1./tar_vTn);
    %% SPECTRAL MATCHING
    
    for i_ = 1:nit % spectral matching iterations
        %%
        % _fourier spectra_
        [inp_vfr,obj_fsa,~,~,~,~] = super_fft(obj_dtm,obj_tha,0);
        inp_vfr = inp_vfr(find(inp_vfr<=1/2/obj_dtm));
        %%
        % _psa response spectra_
        [~,~,~,obj_psa,~] = newmark_sd(obj_tha,obj_dtm,tar_vTn,0.05);
        obj_psa(1) = max(abs(obj_tha));
        %%
        % _psa response spectral ratio_
        obj_rra          = -ones(obj_nfr,1);
        obj_rra(1)       = 1.;
        obj_rra(end)     = obj_psa(1)./tar_psa(1);
        obj_rra(2:end-1) = flip(obj_psa./tar_psa);
        obj_rra = interp1(obj_vfr,obj_rra,inp_vfr,'linear');
        obj_rra(inp_vfr>=frm) = obj_psa(1)./tar_psa(1);
        %%
        % _new (scaled) fourier spectrum_
        obj_fsa = obj_fsa./obj_rra;
        %%
        % _new accelerogram
        obj_tha = detrend(real(ifft(obj_fsa),obj_ntm,'symmetric'),'constant');
    end
    %%
    % _filtering
    [obj_tha,obj_thv,obj_thd] = band_pass_filter(obj_tha,obj_dtm,lfr,hfr);
    %% OUTPUT
    varargout{1} = obj_dtm;
    varargout{2} = obj_tha;
    varargout{3} = obj_thv;
    varargout{4} = obj_thd;
    varargout{5} = obj_vfr;
    varargout{6} = obj_psa;
    return
end
