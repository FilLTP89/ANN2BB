function [varargout]=GUI_13(pos,title,lims,sug)
    global res
    % Demonstrate how to display & change a slider's position with an edit box.
    % Slide the slider and it's position will be shown in the editbox.
    % Enter a valid number in the editbox and the slider will be moved to that
    % position.  If the number entered is outside the range of the slider,
    % the number will be reset.
    %
    %
    % Author:  Matt Fig
    % Date:  7/15/2009
    res=sug;
    S.fh = figure('units','pixels',...
        'position',pos,...
        'menubar','none',...
        'name',title,...
        'numbertitle','off',...
        'resize','off');
    S.sl = uicontrol('style','slide',...
        'unit','pix',...
        'position',[20 10 260 30],...
        'min',lims(1),'max',lims(2),'val',sug);
    S.ed = uicontrol('style','edit',...
        'unit','pix',...
        'position',[20 50 260 30],...
        'fontsize',16,...
        'string',num2str(sug));
    set([S.ed,S.sl],'call',{@ed_call,S});  % Shared Callback.
    S.pb = uicontrol('style','push',...
        'unit','pix',...
        'position',[20 90 260 30],...
        'fontsize',16,...
        'string','OK');
    set(S.pb,'call',{@pb_call,S});  % Shared Callback.
    uiwait(S.fh);
    varargout{1}=res;
    return
end

function [] = ed_call(varargin)
    global res
    % Callback for the edit box and slider.
    [h,S] = varargin{[1,3]};  % Get calling handle and structure.
    
    switch h  % Who called?
        case S.ed
            L = get(S.sl,{'min','max','value'});  % Get the slider's info.
            E = str2double(get(h,'string'));  % Numerical edit string.
            if E >= L{1} && E <= L{2}
                set(S.sl,'value',E)  % E falls within range of slider.
            else
                set(h,'string',L{3}) % User tried to set slider out of range.
            end
        case S.sl
            set(S.ed,'string',get(h,'value')) % Set edit to current slider.
        otherwise
            % Do nothing, or whatever.
    end
    
    return
end

function [] = pb_call(varargin)
    global res
    % Callback for the edit box and slider.
    [h,S] = varargin{[1,3]};  % Get calling handle and structure.
    res = str2double(get(S.ed,'string'));
    
    close(S.fh);
end