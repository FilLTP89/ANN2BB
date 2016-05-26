function [] = GUI_nn_setup()
    %% ANN - DATABASE
    global mon ann Sann
    Sann.flag=0;
    Sann.fh(1) = figure('units','pixels',...
        'position',[0 350 470 60],...
        'menubar','none',...
        'name','ARTIFICIAL NEURAL NETWORK',...
        'numbertitle','off',...
        'resize','off');
    Sann.pb(1) = uicontrol(Sann.fh(1),'style','push',...
        'units','pix',...
        'position',[5 5 150 50],...
        'fontsize',14,...
        'string','HORIZONTAL');
    Sann.pb(2) = uicontrol(Sann.fh(1),'style','push',...
        'units','pix',...
        'position',[160 5 150 50],...
        'fontsize',14,...
        'string','VERTICAL');
    Sann.pb(3) = uicontrol(Sann.fh(1),'style','push',...
        'units','pix',...
        'position',[315 5 150 50],...
        'fontsize',14,...
        'string','CLOSE');
    
    set(Sann.pb(1:3),'callback',{@pb_call,Sann})  % Set callbacks.
    
    uiwait(Sann.fh(1)) % Wait for continue or stop button.
    
    return
end

function [] = pb_call(varargin)
    global mon ann Sann
    % Callback for pushbutton.
    if varargin{1} == Sann.pb(1)
        % Get the structure.
        Sann.fh(2) = figure('units','pixels',...
            'position',[400 300 305 110],...
            'menubar','none',...
            'name','METADATA',...
            'numbertitle','off',...
            'resize','off');
        uicontrol(Sann.fh(2),'style','text','units','pix','position',[5 65 295 25],...
        'string','ANN FILE','fontsize',14,'fontweight','bold');
        Sann.lb(1) = uicontrol(Sann.fh(2),'style','listbox',...
            'units','pix',...
            'position',[5 35 295 25],...
            'fontsize',14,...
            'string','ANN file...');
        Sann.pb(4) = uicontrol(Sann.fh(2),'style','push',...
            'units','pix',...
            'position',[5 5 95 25],...
            'fontsize',14,...
            'string','BROWSE');
        Sann.pb(5) = uicontrol(Sann.fh(2),'style','push',...
            'units','pix',...
            'position',[105 5 95 25],...
            'fontsize',14,...
            'string','CLOSE');
        set(Sann.pb(4:5),'callback',{@pb0_call,Sann})  % Set callbacks.
        uiwait(Sann.fh(2));
        if ann.status(1) == 1;
            ann.hrz = nn_parser(ann.hrz.TnC,ann.hrz.fn);
            disp('-----> ANN(HORIZONTAL): LOADED!');
        end
        
    elseif varargin{1} == Sann.pb(2)
        % Get the structure.
        Sann.fh(3) = figure('units','pixels',...
            'position',[400 300 305 110],...
            'menubar','none',...
            'name','METADATA',...
            'numbertitle','off',...
            'resize','off');
        uicontrol(Sann.fh(3),'style','text','units','pix','position',[5 65 295 25],...
        'string','ANN FILE','fontsize',14,'fontweight','bold');
        Sann.lb(2) = uicontrol(Sann.fh(3),'style','listbox',...
            'units','pix',...
            'position',[5 35 295 25],...
            'fontsize',14,...
            'string','ANN file...');
        Sann.pb(6) = uicontrol(Sann.fh(3),'style','push',...
            'units','pix',...
            'position',[5 5 95 25],...
            'fontsize',14,...
            'string','BROWSE');
        Sann.pb(7) = uicontrol(Sann.fh(3),'style','push',...
            'units','pix',...
            'position',[105 5 95 25],...
            'fontsize',14,...
            'string','CLOSE');
        set(Sann.pb(6:7),'callback',{@pb1_call,Sann})  % Set callbacks.
        uiwait(Sann.fh(3));
        if ann.status(2) == 1;
            ann.vrt = nn_parser(ann.vrt.TnC,ann.vrt.fn);
            disp('-----> ANN(VERTICAL): LOADED!');
        end
    elseif varargin{1} == Sann.pb(3)  % Get the structure.
        %%
        % _exit_
        
        if any(ann.status)
            
            close(Sann.fh(1))  % Found the one we are looking for.
        end
    end
    return
end

function [] = pb0_call(varargin)
    global ann Sann TnC
    % Callback for pushbutton.
    if varargin{1} == Sann.pb(4)  % Get the structure.
        %%
        % _parse trained ANN network - horizontal directions_
        [ann.hrz.fn,pt] = uigetfile('*.mat','Select Monitor Metadata File (*.dat): ',...
            'MultiSelect','off');
        if ~strcmpi(pt(end),filesep)
            pt=strcat(pt,filesep);
        end
        set(Sann.lb(1),'String',ann.hrz.fn);
        ann.hrz.fn = fullfile(pt,ann.hrz.fn);
        ann.hrz.TnC  = GUI_13([500 300 305 120],'CORNER PERIOD',[.01,3],1);
        ann.status(1) = 1;
    elseif varargin{1} == Sann.pb(5)  % Get the structure.
        if ann.status(1) == 1;
            close(Sann.fh(2));
        end
    end
    return
end

function [] = pb1_call(varargin)
    global ann Sann TnC
    % Callback for pushbutton.
    if varargin{1} == Sann.pb(6)  % Get the structure.
        %%
        % _parse trained ANN network - vertical direction_
        [ann.vrt.fn,pt] = uigetfile('*.mat','Select Monitor Metadata File (*.dat): ',...
            'MultiSelect','off');
        if ~strcmpi(pt(end),filesep)
            pt=strcat(pt,filesep);
        end
        set(Sann.lb(2),'String',ann.vrt.fn);
        ann.vrt.fn = fullfile(pt,ann.vrt.fn);
        ann.vrt.TnC  = GUI_13([500 300 305 120],'CORNER PERIOD',[.01,3],1);
        ann.status(2) = 1;
    elseif varargin{1} == Sann.pb(7)  % Get the structure.
        if ann.status(2) == 1;
            close(Sann.fh(3));
        end
    end
    return
end