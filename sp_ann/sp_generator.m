%% *Generate synthetics - Sabetta & Pugliese approach*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _sp_generator_: function to produce non-stationary accelerograms
% according to the approach by Sabetta & Pugliese, 1996.
%% INPUT:
% * _mon (monitor structure)_
% * _str (string: extra metadata file)
%% OUTPUT:
% * _sps (structure of sp synthetics)_
%% N.B.:
% it has been verified that with ssc = 2 the method of sabetta
% gives as output acceleration THs with an average response spectrum
% compliant with the (median+isig/2) spectrum given by the GMPE
% Need for _sabetta.m_
%% REFERENCES
% @Article{Sabetta_Pugliese_1996,
%   author =  {Sabetta, F. and Pugliese, A.},
%   title =   {{Estimation of Response Spectra and Simulation of Nonstationary Earthquake Ground Motions}},
%   journal = {Bulletin of the Seismological Society of America},
%   year =    {1996},
%   volume =  {86},
%   number =  {2},
%   pages =   {337--352},
% }
function [varargout] = sp_generator(varargin)
    %% *SET-UP*
    % _parsing monitor metadata_
    sps.mon    = varargin{1};
    sps.mtd.fn = varargin{2};
    %% *PARSING EXTRA METADATA*
    str = importdata(sps.mtd.fn);
    % earthquake magnitude
    sps.mtd.mw  = str(1,1);
    % time-step
    sps.mtd.dtm = str(2,1);
    % site condition factor
    sps.mtd.scc = str(3,1);
    % site condition std
    sps.mtd.sst = str(4,1);
    % scale factor
    sps.mtd.scl = str(5,1);
    % output flag
    sps.mtd.ivd = str(6,1);
    [s1,s2] = size(sps.mon.cp(:));
    [~,idx] = max([s1,s2]);
    %% *SABETTA & PUGLIESE SYNTHETICS*
    sps.mon.dtm = repmat(sps.mtd.dtm,sps.mon.na,1);
    [sps.mon.vtm,tha] = arrayfun(@(x) sabetta(sps.mtd.mw,x,sps.mtd.scc,...
        sps.mtd.sst,sps.mtd.dtm,1,sps.mtd.scl),sps.mon.dep(:),'UniformOutput',0);
    sps.mon.ntm = cellfun(@(x) numel(x),sps.mon.vtm);
    bpf_out = cell(1,3);
    [bpf_out{:}] = cellfun(@(x) band_pass_filter(sps.mtd.dtm,x,sps.mon.lfr,sps.mon.hfr),tha,'UniformOutput',0);
    
    tha = arrayfun(@(x) cell2struct(repmat(x,sps.mon.nc,1),(sps.mon.cp(:)),1),...
        bpf_out{1}(:),'UniformOutput',0);
    thv = arrayfun(@(x) cell2struct(repmat(x,sps.mon.nc,1),(sps.mon.cp(:)),1),...
        bpf_out{2}(:),'UniformOutput',0);
    thd = arrayfun(@(x) cell2struct(repmat(x,sps.mon.nc,1),(sps.mon.cp(:)),1),...
        bpf_out{3}(:),'UniformOutput',0);
    
    for i_=1:sps.mon.na
        sps.syn{i_} = struct('tha',tha{i_},'thv',thv{i_},'thd',thd{i_});
        
%         [vtm_pga,pga,vtm_pgv,pgv,vtm_pgd,pgd] = cellfun(@(a,v,d) ...
%             cellfun(@(x) PGAVD_eval(sps.mtd.dtm,a.(x),v.(x),d.(x)),...
%             sps.mon.cp,'UniformOutput',0),...
%             {sps.syn{i_}.tha},{sps.syn{i_}.thv},{sps.syn{i_}.thd},'UniformOutput',0);
%         
%         sps.syn{i_}.pga = cell2struct(vtm_pga{1},sps.mon.cp(:),idx);
%         a = cell2struct(pga{1},sps.mon.cp(:),idx);
%         sps.syn{i_}.pga = cell2struct(cellfun(@(x) cat(1,sps.syn{i_}.pga.(x),a.(x)),...
%             sps.mon.cp,'UniformOutput',0),sps.mon.cp(:),idx);
%         
%         sps.syn{i_}.pgv = cell2struct(vtm_pgv{1},sps.mon.cp(:),idx);
%         v = cell2struct(pgv{1},sps.mon.cp(:),idx);
%         sps.syn{i_}.pgv = cell2struct(cellfun(@(x) cat(1,sps.syn{i_}.pgv.(x),v.(x)),...
%             sps.mon.cp,'UniformOutput',0),sps.mon.cp(:),idx);
%         
%         sps.syn{i_}.pgd = cell2struct(vtm_pgd{1},sps.mon.cp(:),idx);
%         d = cell2struct(pgd{1},sps.mon.cp(:),idx);
%         sps.syn{i_}.pgd = cell2struct(cellfun(@(x) cat(1,sps.syn{i_}.pgd.(x),d.(x)),...
%             sps.mon.cp,'UniformOutput',0),sps.mon.cp(:),idx);
    end
    %% OUTPUT
    varargout{1} = sps;
    return
end
