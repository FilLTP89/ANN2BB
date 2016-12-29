%% *Parse ITACA record file*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _parse_itaca_file: function to parse and synchronize records from ITACA database
%% INPUT:
% * _fn (file name)_
%% OUTPUT:
% * _dtm (time-step)_
% * _ntm (number of time-steps)_
% * vtm (time-vector)_
% * _th (structure of parsed th)_

function [varargout] = parse_itaca_file_new(varargin)
    %% *SET-UP*
    fn = varargin{1};
    tag_string(1) = {'EVENT_TIME_HHMMSS:'};
    tag_string(2) = {'DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS:'};
    tag_string(3) = {'SAMPLING_INTERVAL_S:'};
    tag_string(4) = {'UNITS:'};
    tag_string(5) = {'USER5:'};
    
    fid = fopen(fn);
    if fid==-1
        keyboard
        error('file does not exist');
    end
    
    data = textscan(fid,'%s');
    fclose all;
    for i_ = 1:numel(tag_string)
        idx(i_) = find(strcmpi(data{1},tag_string{i_})==1)+1;
    end
    first_sample_time1 = data{1}{idx(1)};
    
    first_sample_time = str2double(first_sample_time1(1:2))*60*60+...
        str2double(first_sample_time1(3:4))*60 + ...
        str2double(first_sample_time1(5:6));
    
    eqk_time1 = data{1}{idx(2)}(end-9:end-4);
    eqk_time = str2double(eqk_time1(1:2))*60*60+str2double(eqk_time1(3:4))*60 + ...
        str2double(eqk_time1(5:6));
    
    dtm = str2double(data{1}(idx(3)));
    scale = 1/100; % cm/s2
    time_shift = abs(first_sample_time-eqk_time);
    if isempty(time_shift)
        time_shift = 0;
    end
    nsample_shift = round(time_shift/dtm);
    
    %% *READ DATA*
    th = cellfun(@(x) str2num(x),data{1}(idx(5):end))*scale;
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