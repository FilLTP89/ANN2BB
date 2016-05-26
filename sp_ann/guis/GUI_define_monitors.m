function [] =  GUI_define_monitors()
    global mon Smon H
    %% MONITOR METADATA
    Smon.fh(1) = figure('units','pixels',...
        'position',[0 500 305 355],...
        'menubar','none',...
        'name','METADATA',...
        'numbertitle','off',...
        'resize','off');
    Smon.flag = zeros(3,1);
    %%
    % _WORKDIR_    
    uicontrol(Smon.fh(1),'style','text','units','pix','position',[5 325 295 25],...
        'string','WORK-DIRECTORY');
    Smon.lb(1) = uicontrol(Smon.fh(1),'style','listbox',...
        'units','pix',...
        'position',[5 295 295 25],...
        'string','Working Directory...');
    Smon.pb(1) = uicontrol('style','push',...
        'units','pix',...
        'position',[105 265 100 25],...
        'string','BROWSE');
    %%
    % _MONITOR FOLDER_
    uicontrol(Smon.fh(1),'style','text','units','pix','position',[5 235 295 25],...
        'string','MONITOR-DIRECTORY');
    Smon.lb(2) = uicontrol('style','listbox',...
        'units','pix',...
        'position',[5 205 295 25],...
        'string','Monitor Folder...','parent',Smon.fh(1));
    Smon.pb(2) = uicontrol('style','push',...
        'units','pix',...
        'position',[105 175 100 25],...
        'string','BROWSE');
    %%
    % _METADATA FILE_
    uicontrol(Smon.fh(1),'style','text','units','pix','position',[5 145 295 25],...
        'string','METADATA-FILE','fontsize',14,'fontweight','bold');
    Smon.lb(3) = uicontrol('style','listbox',...
        'units','pix',...
        'position',[5 115 295 25],...
        'string','Metadata File...','parent',Smon.fh(1));
    Smon.pb(3) = uicontrol('style','push',...
        'units','pix',...
        'position',[105 85 100 25],...
        'string','LOAD');
    
    %%
    % _OK_
    Smon.pb(4) = uicontrol(Smon.fh(1),'style','push',...
        'units','pix',...
        'position',[105 5 95 25],...
        'string','OK');
    set(Smon.pb(:),'callback',{@pb_call,Smon})  % Set callbacks.
    
    uiwait(Smon.fh(1));
    
    %% MONITOR SELECTION
    Smon.flag=zeros(3,1);
    
    Smon.fh(2) = figure('units','pixels',...
        'position',[0 5 305 350],...
        'menubar','none',...
        'name','MONITOR SELECTION',...
        'numbertitle','off',...
        'resize','off');
    %%
    % _MONITOR SELECTION_
    uicontrol(Smon.fh(2),'style','text','units','pix','position',[5 125 295 25],...
        'string','MONITOR-FILE','fontsize',14,'fontweight','bold');
    Smon.lb(4) = uicontrol(Smon.fh(2),'style','listbox',...
        'units','pix',...
        'position',[5 95 295 25],...
        'fontsize',14,...
        'string','Monitor Files...');
    Smon.pb(5) = uicontrol(Smon.fh(2),'style','push',...
        'units','pix',...
        'position',[105 65 100 25],...
        'fontsize',14,...
        'string','BROWSE',...
        'callback',{@pb_call,Smon});
    
    %%
    % _TYPE OF NUMERICAL SIMULATION_
    Smon.bg(1) = uibuttongroup(Smon.fh(2),'units','pix',...
        'pos',[5 300 295 50],'Title','Source');
    Smon.rd(1) = uicontrol(Smon.bg(1),...
        'style','rad',...
        'unit','pix',...
        'position',[85 5 110 40],...
        'string','HISADA','callback',{@tp_call,Smon});
    Smon.rd(2) = uicontrol(Smon.bg(1),...
        'style','rad',...
        'unit','pix',...
        'position',[190 5 110 40],...
        'string','SPEED','callback',{@tp_call,Smon});
    
    %%
    % _DIRECTIONS_
    Smon.bg(2) = uibuttongroup(Smon.fh(2),'units','pix',...
        'pos',[5 235 295 50],'Title','Direction');
    Smon.rd(3) = uicontrol(Smon.bg(2),...
        'style','checkbox',...
        'unit','pix',...
        'position',[15 5 90 40],...
        'string','X','callback',{@tp_call,Smon});
    Smon.rd(4) = uicontrol(Smon.bg(2),...
        'style','checkbox',...
        'unit','pix',...
        'position',[105 5 90 40],...
        'string','Y','callback',{@tp_call,Smon});
    Smon.rd(5) = uicontrol(Smon.bg(2),...
        'style','checkbox',...
        'unit','pix',...
        'position',[200 5 90 40],...
        'string','Z','callback',{@tp_call,Smon});
    
    H=Smon.bg(2).Children;
    %%
    % _confirm and close_
    Smon.pb(6) = uicontrol(Smon.fh(2),'style','push',...
        'units','pix',...
        'position',[105 5 95 25],...
        'string','OK',...
        'callback',{@pb_call,Smon});
    
    uiwait(Smon.fh(2));
    
    
    return
end

function [] = tp_call(varargin)
    global mon Smon H
    % hObject    handle to edit1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hob = varargin{1};
    hnd = varargin{3};
    switch hob
        case {Smon.rd(1),Smon.rd(2)}
            switch findobj(get(Smon.bg(1),'selectedobject'))
                case Smon.rd(1)
                    mon.tp  = 'H';
                    mon.rc  = {'v'};
                    disp('HISADA SIMULATIONS SELECTED!');
                    
                case Smon.rd(2)
                    mon.tp  = 'S';
                    mon.rc  = {'d'};
                    disp('SPEED SIMULATIONS SELECTED!');
                    
                otherwise
                    mon.tp  = 'S';
                    mon.rc  = {'d'};
                    warning('DEFAULT CHOICE: SPEED SIMULATIONS SELECTED!');
            end
            
            Smon.flag(2)=1;
            mon.nr  = numel(mon.rc);
        case {Smon.rd(3),Smon.rd(4),Smon.rd(5)}
            cp={'x';'y';'z'};
            idx = flip(arrayfun(@(x) x.get.Value,H));
            
            mon.cp = cp(logical(idx));
            mon.nc  = numel(mon.cp);
            if mon.nc>0
                Smon.flag(3)=1;
                try
                    cellfun(@(x) fprintf('%s-DIRECTION SELECTED!\n',upper(x)),mon.cp);
                catch
                    keyboard
                end
            end
    end
    
    return
end

function [] = pb_call(varargin)
    global mon Smon
    % hObject    handle to edit1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hob = varargin{1};
    hnd = varargin{3};
    
    switch hob
        case Smon.pb(1)
            %%
            % _work directory_
            mon.wd = uigetdir(userpath,'Select Working Directory: ');
            if ~strcmpi(mon.wd(end),filesep)
                strcat(mon.wd,filesep);
            end
            set(Smon.lb(1),'String',mon.wd);
            
            Smon.flag(1) = 1;
            
        case Smon.pb(2)
            %%
            % _monitor folder_
            mon.pt = uigetdir(userpath,'Select Monitor Directory: ');
            if ~strcmpi(mon.pt(end),filesep)
                mon.pt=strcat(mon.pt,filesep);
            end
            set(Smon.lb(2),'String',mon.pt);
            Smon.flag(2) = 1;
        case Smon.pb(3)
            %%
            % _metadata file_
            [mon.fn,pt] = uigetfile('*_mtd.csv','Select Monitor Metadata File (*_mtd.csv): ',...
                'MultiSelect','off');
            if ~strcmpi(pt(end),filesep)
                pt=strcat(pt,filesep);
            end
            set(Smon.lb(3),'String',mon.fn);
            mon.fn = fullfile(pt,mon.fn);
            Smon.flag(3) = 1;
        case Smon.pb(4)
            %%
            % _confirm selection_
            if all(Smon.flag)
                mon.status(1) = 1;
                close(Smon.fh(1));
            end
        case Smon.pb(5)
            %%
            % _monitor file_
            fmt = {'*.d;*_dat.*'};
            fmt=arrayfun(@(x) fullfile(mon.pt,x),fmt);
            [id,pt] = uigetfile(fmt,...
                'Select Monitor Files: ',...
                'MultiSelect','on');
            if ~strcmpi(pt(end),filesep)
                pt=strcat(pt,filesep);
            end
            if ~strcmpi(pt,mon.pt);
                warning('MONITOR PATH MISMATCH! PATH UPDATED');
                keyboard
                mon.pt=pt;
            end
            %%
            % _monitor identity_
            if isfield(mon,'tp')
                if strcmpi(mon.tp,'s')
                    if iscell(id)
                        tmp = all(cellfun(@(x) ~isempty(strfind(x,'monitor')),id));
                        idt = str2mat((cellfun(@(x) x(end-6:end-2),id,'UniformOutput',0)));
                    else
                        tmp = ~isempty(strfind(id,'monitor'));
                        idt = id(end-6:end-2);
                    end
                    
                elseif strcmpi(mon.tp,'h')
                    if iscell(id)
                        tmp = all(cellfun(@(x) ~isempty(strfind(x,'_dat')),id));
                        idt = str2mat((cellfun(@(x) x(end),id,'UniformOutput',0)));
                    else
                       tmp=~isempty(strfind(id,'_dat'));
                       idt = id(end);
                    end
                    
                end
                if tmp
                    mon.id = str2num(idt);
                    mon.na  = numel(mon.id);
                    set(Smon.lb(4),'String',id);
                    Smon.flag(1) = 1;
                else
                    disp('ERROR: FILE NOT VALID!');
                end
            else
                disp('ERROR: FIRST DEFINE TYPE OF NUMERICAL SIMULATION!');
            end
        case Smon.pb(6)
            if all(Smon.flag)
                mon.status(2) = 1;
                
                close(Smon.fh(2));
            else
                if Smon.flag(1)==0
                    disp('NO MONITOR FILE SELECTED!!');
                end
                
                if Smon.flag(2)==0
                    disp('NO NUMERICAL SIMULATION SELECTED!!');
                end
                
                if Smon.flag(3)==0
                    disp('NO DIRECTION SELECTED!!');
                end
            end
    end
    
    return
end