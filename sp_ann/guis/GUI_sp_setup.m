function [] = GUI_sp_setup()
    %% SYNTHETICS GENERATOR
    global mon sps Ssps
    Ssps.flag=0;
    Ssps.fh(1) = figure('units','pixels',...
        'position',[0 350 360 60],...
        'menubar','none',...
        'name','Generate Synthetics',...
        'numbertitle','off',...
        'resize','off');
    Ssps.pb(1) = uicontrol(Ssps.fh(1),'style','push',...
        'units','pix',...
        'position',[5 5 110 50],...
        'string','GENERATE');
    Ssps.pb(2) = uicontrol(Ssps.fh(1),'style','push',...
        'units','pix',...
        'position',[120 5 110 50],...
        'string','SPECTRA');
    Ssps.pb(3) = uicontrol(Ssps.fh(1),'style','push',...
        'units','pix',...
        'position',[235 5 110 50],...
        'string','CLOSE');
    
    set(Ssps.pb(1:3),'callback',{@pb_call,Ssps})  % Set callbacks.
    
    uiwait(Ssps.fh(1)) % Wait for continue or stop button.
    
    return
end

function [] = pb_call(varargin)
    global mon sps Ssps
    % Callback for pushbutton.
    press_button(varargin{1});
    
    if varargin{1} == Ssps.pb(1)  % Get the structure.
        %%
        % _define metadata_
        Ssps.fh(2) = figure('units','pixels',...
            'position',[400 300 305 110],...
            'menubar','none',...
            'name','METADATA',...
            'numbertitle','off',...
            'resize','off');
        uicontrol(Ssps.fh(2),'style','text','units','pix','position',[5 65 295 25],...
            'string','METADATA');
        Ssps.lb(1) = uicontrol(Ssps.fh(2),'style','listbox',...
            'units','pix',...
            'position',[5 35 295 25],...
            'string','Metadata file...');
        Ssps.pb(4) = uicontrol(Ssps.fh(2),'style','push',...
            'units','pix',...
            'position',[105 5 100 25],...
            'string','BROWSE');
        set(Ssps.pb(4),'callback',{@pb0_call,Ssps})  % Set callbacks.
        uiwait(Ssps.fh(2)) % Wait for continue or stop button.
        %%
        % _parse simulations_
        if Ssps.flag
            sps = sp_generator(mon,sps.mon.fn);
        else
            Ssps.flag = 0;
        end
        
    elseif varargin{1} == Ssps.pb(2)
        %%
        % _compute spectra_
        sps = sp_spectra(sps);
        disp('SPECTRA COMPUTED!');
    elseif varargin{1} == Ssps.pb(3)
        %%
        % _exit_
        sps.status=all(Ssps.flag);
        close(Ssps.fh(1))  % Found the one we are looking for.
    end
    switch varargin{1}
        case{Ssps.pb(1),Ssps.pb(2)}
            press_button(varargin{1},'OK!');
    end
    return
end

function [] = pb0_call(varargin)
    global mon sps Ssps
    % Callback for pushbutton.
    if varargin{1} == Ssps.pb(4)  % Get the structure.
        %%
        % _metadata file_
        [sps.mon.fn,pt] = uigetfile('*.dat','Select Monitor Metadata File (*.dat): ',...
            'MultiSelect','off');
        if ~strcmpi(pt(end),filesep)
            pt=strcat(pt,filesep);
        end
        set(Ssps.lb(1),'String',sps.mon.fn);
        sps.mon.fn = fullfile(pt,sps.mon.fn);
        
        close(Ssps.fh(2))  % Found the one we are looking for.
        Ssps.flag = 1;
        return
    end
end
