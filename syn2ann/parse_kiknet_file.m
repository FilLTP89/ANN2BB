%% *Parse KNET-KIKNET record filename*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% parse_kiknet_file: function to parse and synchronize records from KIKNET database
%% INPUT:
% * fn (file name)
%% OUTPUT:
% * dtm (time-step)
% * ntm (number of time-steps)
% * vtm (time-vector)
% * th (structure of parsed th)

function [varargout] = parse_kiknet_file(varargin)
    %% *SET-UP*
    fn = varargin{1};
    tag_string(1) = {'Sampling'};
    tag_string(2) = {'Scale'};
    tag_string(3) = {'Memo.'};
    tag_string(4) = {'Record'};
    tag_string(5) = {'Origin'};
    
    fid = fopen(fn);
    if fid==-1
        error('file does not exist');
    end
    
    data = textscan(fid,'%s');
    fclose all;
    for i_ = 1:numel(tag_string)
        idx(i_) = find(strcmpi(data{1},tag_string{i_})==1)+2;
    end
    
    dtm = 1/str2double(data{1}{idx(1)}(1:3));
    
    eqk_time = str2double(data{1}{idx(4)+1}(1:2))*60*60+...
        str2double(data{1}{idx(4)+1}(4:5))*60 + ...
        str2double(data{1}{idx(4)+1}(7:8));
    
    first_sample_time = str2double(data{1}{idx(5)+1}(1:2))*60*60+...
        str2double(data{1}{idx(5)+1}(4:5))*60 + ...
        str2double(data{1}{idx(5)+1}(7:8));
    
    time_shift    = abs(first_sample_time-eqk_time);
    nsample_shift = round(time_shift/dtm);
    
    str = data{1}{idx(2)};
    idx0 = findstr('(gal)',str);
    scale = str2num(str(1:idx0(1)-1))./str2num(str(idx0(end)+6:end));
    
    %% *READ DATA*
    th = cellfun(@(x) str2double(x),data{1}(idx(3)-1:end))*scale;
%     
%     th = th(nsample_shift+1:end);
%     if isempty(th)
%         keyboard
%     end
    ntm = numel(th);
    vtm = dtm*(0:ntm-1)';
    %% *OUTPUT*
    varargout{1} = dtm;
    varargout{2} = ntm;
    varargout{3} = vtm;
    varargout{4} = th;
    varargout{5} = nsample_shift;
%     dtm = cell2mat(metadata(32));dtm=str2double(dtm(1:3));
%     scale_factor=cell2mat(metadata(40));
%     idx=strfind(scale_factor,'(gal)/');
%     scale_factor = str2double(scale_factor(1:idx-1))/str2double(scale_factor(idx+6:end));
%     
%     ntm = numel(metadata)-49;    
%     for i_=50:numel(metadata)
%         tha(i_-49) = str2double(cell2mat(metadata(i_)))*scale_factor;
%     end
%     
%     vtm = (0:ntm-1)*dtm;
%     
%     varargout{1} = dtm;
%     varargout{2} = ntm;
%     varargout{3} = vtm(:);
%     varargout{4} = tha(:);
    return
end