function [varargout] = spectral_scaling_original(varargin)
    %% *SET-UP*
    inp_dtm = varargin{1};
    inp_tha = varargin{2}(:);
    tar_vTn = varargin{3}(:);
    tar_psa = varargin{4}(:);
    tar_pga = tar_psa(1);     % target pga
    %
    % _spectral matching set up_
    %
    zeta = 0.05;                 % damping coefficient
    fac  = 1;                    % resampling factor
    scl  = 1;                    % scaling factor
    ni   = 10;                   % number of iterations
    
    %% *SCALING/RESAMPLING ACCELERATION*
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);

    %% *TARGET SPECTRUM ENRICHMENT AT SP*
    vTn_cor = [2*obj_dtm,0.05];                       % upper period of refinement
    nc      = 10;                                     % number of refinement points
    idxn    = (tar_vTn(tar_vTn~=0) < vTn_cor(2)); % indexes
    
    if logical(any(idxn))
        idx_save_org(:,1)    = find(tar_vTn>vTn_cor(2));
        idx_save_new         = idx_save_org+(nc-idx_save_org(1)+1);
        %
        % _new vector of natural period_
        %
        dTn_cor              = (log10(vTn_cor(2))-log10(vTn_cor(1)))/(nc-1);
        tar_vTn_new(1:nc,1)  = (log10(vTn_cor(1)):dTn_cor:log10(vTn_cor(2)))';
        tar_vTn_new(1:nc,1)  = 10.^tar_vTn_new(1:nc,1);
        tar_vTn_new(idx_save_new,1) = tar_vTn(idx_save_org,1);
        tar_vTn_new          = [0;tar_vTn_new];
        %
        % _new vector of psa_
        %
        tar_psa_new          = interp1(tar_vTn,tar_psa,tar_vTn_new,'linear');
        tar_vTn = tar_vTn_new;
        tar_psa = tar_psa_new;
    end
    tar_nT  = numel(tar_vTn);
    
    %% *INITAL PSA AND FSA*
    inp_nfr = 2^nextpow2(obj_ntm);
    inp_dfr = 1/(inp_nfr-1)/obj_dtm;
    inp_vfr = (0:inp_dfr:0.5/obj_dtm)';
    inp_nfr = numel(inp_vfr);
    obj_vfr = [0;flip(1./tar_vTn(2:end));0.5/obj_dtm];
    %
    % _pad time-histories_
    %
    obj_tha(obj_ntm+1:inp_nfr)=0;
    %
    % _pseudo-spectral acceleration_
    %
    obj_psa = SDOF_response(obj_tha,obj_dtm,tar_vTn,zeta);
    rra_in=-999*ones(tar_nT+1,1);
    for i_ = 2:tar_nT
      rra_in(tar_nT-i_+2) = obj_psa(i_)/tar_psa(i_);
    end
    rra_in(tar_nT+1) = obj_psa(1)/tar_psa(1);
    rra_in(1)=1;
    rra = interp1(obj_vfr,rra_in,inp_vfr,'linear');
    
    
    rra_pro=-999*ones(tar_nT+1,1);
    for i_=1:ni
        obj_fsa=fft(obj_tha);
        obj_fsa(1:inp_nfr) = obj_fsa(1:inp_nfr)./rra(1:inp_nfr);
        
        for m_ = inp_nfr+2:2*inp_nfr
            obj_fsa(m_) = conj(obj_fsa(2*inp_nfr-m_+2));
        end
        obj_tha = ifft(obj_fsa(:),2*inp_nfr,1,'symmetric');
        obj_tha = detrend(obj_tha);
        obj_tha = adjust_pga(obj_dtm,obj_tha,tar_pga);
        % response spectrum of corrected waveform
        psa_pro = SDOF_response(obj_tha,obj_dtm,tar_vTn,zeta,1);
        
        rra_pro(1)=1;
        for j_=1:tar_nT
                rra_pro(tar_nT-j_+2)=psa_pro(j_)/tar_psa(j_);
        end
        
        rra = interp1(obj_vfr,rra_pro,inp_vfr,'linear');
        
        %         for j_=1:inp_nfr
        %             if (inp_vfr(j_)>f_amax)
        %                 rra(j_)=psa_pro(1)/tar_psa(1);
        %             end
        %         end
        
    end
    obj_ntm = numel(obj_tha);
%     [obj_tha,obj_thv,obj_thd,obj_vtm,pads] = ...
%         bpf_tha(obj_dtm,obj_tha,0.05,[]);
    
%     [bfb,bfa,~] = create_butter_filter(3,0.05,[],0.5/obj_dtm);
%     obj_tha = detrend(obj_tha);
%     obj_tha = filtfilt(bfb,bfa,obj_tha);
    obj_thv = cumtrapz(obj_tha)*obj_dtm;
    obj_thv = detrend(obj_thv);
    obj_thd = cumtrapz(obj_thv)*obj_dtm;
    obj_thd = detrend(obj_thd,'linear');
    frac=2.5;
    obj_thd = taper_fun(obj_thd,frac,1,0);
    
%     obj_thv=diff(obj_thd)/obj_dtm;
%     obj_thv(obj_ntm)=0;
%     obj_tha=diff(obj_thv)/obj_dtm;
%     obj_tha(obj_ntm)=0;
%         
    %% *SCALED PSA SPECTRUM*
    
    % response spectral parameters
    dTn = 0.05;
    Tmax = 5;
    Tn = (0:dTn:Tmax)';
    psa_tha_pro = SDOF_response(obj_tha',obj_dtm,Tn,zeta,1);
    
    %% *OUTPUT*
    varargout{1} = obj_dtm;
    varargout{2} = obj_tha;
    varargout{3} = Tn;
    varargout{4} = psa_tha_pro;
    varargout{5} = obj_thv;
    varargout{6} = obj_thd;
    return
end
