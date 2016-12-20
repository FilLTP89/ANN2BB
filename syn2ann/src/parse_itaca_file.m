%% *Parse ITACA record file*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _parse_itaca_file: function to parse and synchronize records from ITACA database
%% INPUT:
% * fn (file name)
%% OUTPUT:
% * dtm (time-step)
% * ntm (number of time-steps)
% * vtm (time-vector)
% * th (structure of parsed th)

function [varargout] = parse_itaca_file(varargin)
    %% *SET-UP*
    fn = varargin{1};
    nlin = 56;
    scale = 1/100; % cm/s2
    
    fid = fopen(fn);
    if fid==-1
        keyboard
        error('file does not exist');
    end
    %% *READ HEADER*
    for j=1:4
        tline = fgetl(fid);
    end
    
    eqk_time1 = tline(end-5:end);
    eqk_time = str2num(eqk_time1(1:2))*60*60+str2num(eqk_time1(3:4))*60 + ...
        str2num(eqk_time1(5:6));
    
    for j=1:22
        tline = fgetl(fid);
    end
    
    first_sample_time1 = tline(end-9:end);
    first_sample_time = str2num(first_sample_time1(1:2))*60*60+...
        str2num(first_sample_time1(3:4))*60 + ...
        str2num(first_sample_time1(5:6));
    
    time_shift = abs(first_sample_time-eqk_time);
    
    for j=1:2
        tline = fgetl(fid);
    end
    dtm = str2num(tline(end-8:end));
    
    if isempty(time_shift) == 1
        time_shift = 0;
    end
    
    nsample_shift = round(time_shift/dtm);
    %% *READ DATA*
    frewind(fid);
    th = textscan(fid,'','HeaderLines',nlin); fclose(fid);
    th = cell2mat(th).*scale;
    th = th(nsample_shift+1:end);
    ntm = numel(th);
    vtm = dtm*(0:ntm-1)';
    %% *OUTPUT*
    varargout{1} = dtm;
    varargout{2} = ntm;
    varargout{3} = vtm;
    varargout{4} = th;
    return
end