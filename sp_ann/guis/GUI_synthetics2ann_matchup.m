function [] = GUI_synthetics2ann_matchup()
    
    global mon nss sps hbs ann smr out S data
    
    mon.status = 0;
    S.fh = figure('units','pixels',...
        'position',[375 600 976 60],...
        'menubar','none',...
        'name','MAIN MENU',...
        'numbertitle','off',...
        'resize','off');
    
    S.pb(1) = uicontrol(S.fh,'style','push',...
        'unit','pix',...
        'position',[5 5 156 50],...
        'fontweight','bold','fontsize',14,'string','MONITOR',...
        'callback',{@pb_call});%,...
    %         'backgroundc',[0.94 .94 .94]%,...
    %         'busyaction','cancel',...% So multiple pushes don't stack.
    %         'interrupt','off');
    
    S.pb(2) = uicontrol(S.fh,'style','push',...
        'unit','pix',...
        'position',[166 5 156 50],...
        'fontweight','bold','fontsize',14,'string','SIMULATION',...
        'callback',{@pb_call},...
        'backgroundc',[0.94 .94 .94],...
        'busyaction','cancel',...% So multiple pushes don't stack.
        'interrupt','off');
    S.pb(3) = uicontrol(S.fh,'style','push',...
        'unit','pix',...
        'position',[327 5 156 50],...
        'fontweight','bold','fontsize',14,'string','S&P96',...
        'callback',{@pb_call},...
        'backgroundc',[0.94 .94 .94],...
        'busyaction','cancel',...% So multiple pushes don't stack.
        'interrupt','off');
    S.pb(4) = uicontrol(S.fh,'style','push',...
        'unit','pix',...
        'position',[488 5 156 50],...
        'fontweight','bold','fontsize',14,'string','HYBRIDIZE',...
        'callback',{@pb_call},...
        'backgroundc',[0.94 .94 .94],...
        'busyaction','cancel',...% So multiple pushes don't stack.
        'interrupt','off');
    S.pb(5) = uicontrol(S.fh,'style','push',...
        'unit','pix',...
        'position',[649 5 156 50],...
        'fontweight','bold','fontsize',14,'string','ANN',...
        'callback',{@pb_call},...
        'backgroundc',[0.94 .94 .94],...
        'busyaction','cancel',...% So multiple pushes don't stack.
        'interrupt','off');
    S.pb(6) = uicontrol(S.fh,'style','push',...
        'unit','pix',...
        'position',[810 5 156 50],...
        'fontweight','bold','fontsize',14,'string','MATCH',...
        'callback',{@pb_call},...
        'backgroundc',[0.94 .94 .94],...
        'busyaction','cancel',...% So multiple pushes don't stack.
        'interrupt','off');
    data=guihandles(S.fh);
    guidata(S.fh,data)
    uiwait(S.fh);
    return
end

function [] = pb_call(varargin)
    global mon nss sps hbs ann smr S data
    % hObject    handle to edit1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    hob = varargin{1};
    hnd = varargin{2}.Source;
    % Callback for pushbutton.
    switch hob
        case S.pb(1)  % Get the structure.
            %% DEFINE MONITOR STRUCTURE
            col = get(S.pb(1),'backg');  % Get the background color of the figure.
            set(S.pb(1),'str','RUNNING...','backg',[1 .6 .6]);
            pause(.01);
            set(S.pb(1),'str','Runnning...','backg',col);
            mon.status(1)=0;
            GUI_define_monitors;
%             count=0;
%             while ~all(mon.status)
%                 if count>0
%                     warning('    ERROR: METADATA FILES SELECTION!!');
%                 else
%                     disp('0.     SELECT METADATA:    ');
%                 end
%                 GUI_define_monitors;
%                 count=count+1;
%             end
            disp('-----> METADATA FILES : OK!');
            set(S.pb(1),'str','PARSED','backg',col)
        case S.pb(2)
            
            if all(mon.status)
                col = get(S.pb(2),'backg');  % Get the background color of the figure.
                set(S.pb(2),'str','RUNNING...','backg',[1 .6 .6]);
                pause(.01);
                set(S.pb(2),'str','Runnning...','backg',col);
                %% NUMERICAL SIMULATIONS: PARSE AND COMPUTE SPECTRA
                nss.status=0;
                GUI_ns_setup;
%                 count=0;
%                 while ~nss.status
%                     if count>0
%                         warning('    ERROR: NUMERICAL SIMULATION PROCESSING!!');
%                     else
%                         disp('1.     NUMERICAL SIMULATION:    ');
%                     end
%                     GUI_ns_setup;
%                     count=count+1;
%                 end
                disp('-----> NUMERICAL SIMULATIONS : OK!');
                set(S.pb(2),'str','LOADED','backg',col);
            else
                disp('      ERROR: FIRST DEFINE MONITORS!');
            end
        case S.pb(3)
            if all(mon.status)
                col = get(S.pb(3),'backg');  % Get the background color of the figure.
                set(S.pb(3),'str','RUNNING...','backg',[1 .6 .6]);
                pause(.01);
                set(S.pb(3),'str','Runnning...','backg',col);
                %% GENERATE SABETTA PUGLIESE SYNTHETICS
                sps.status=0;
                GUI_sp_setup;
                disp('-----> SABETTA-PUGLIESE SYNTHETICS: GENERATED!');
                set(S.pb(3),'str','GENERATED','backg',col)
            else
                disp('      ERROR: FIRST DEFINE MONITORS!');
            end
        case S.pb(4)
            
            if all(mon.status)&&nss.status&&sps.status
                col = get(S.pb(4),'backg');  % Get the background color of the figure.
                set(S.pb(4),'str','RUNNING...','backg',[1 .6 .6]);
                pause(.01);
                set(S.pb(4),'str','Runnning...','backg',col);
                %% HYBRIDIZE SYNTHETICS
                hbs.status=0;
                GUI_hb_setup;
                disp('-----> HYBRIDIZATION: DONE!');
                set(S.pb(4),'str','DONE','backg',col);
            else
                disp('      ERROR: FIRST DEFINE MONITORS!');
            end
        case S.pb(5)
            
            col = get(S.pb(5),'backg');  % Get the background color of the figure.
            set(S.pb(5),'str','RUNNING...','backg',[1 .6 .6]);
            pause(.01);
            set(S.pb(5),'str','Runnning...','backg',col);
            %% LOAD TRAINED ANN
            ann.status=zeros(2,1);
            GUI_nn_setup;
            disp('-----> ANN : LOADED!');
            set(S.pb(5),'str','LOADED','backg',col);
            
        case S.pb(6)
            if all(mon.status)&&hbs.status&&any(ann.status)
                smr.trs=[];
                col = get(S.pb(6),'backg');  % Get the background color of the figure.
                set(S.pb(6),'str','RUNNING...','backg',[1 .6 .6]);
                pause(.01);
                set(S.pb(6),'str','Runnning...','backg',col);
                %% SPECTRAL MATCHING
                % _apply trained ANN on hybrid accelerograms (horizontal)_
                if ann.status(1)
                    if ismember({'x'},hbs.mon.cp)
                        smr.trs = [smr.trs;{ann2hbs_train(hbs,ann.hrz,{'x'})}];
                        smr.status(1)=0;
                    end
                    if ismember({'y'},hbs.mon.cp)
                        smr.trs = [smr.trs;{ann2hbs_train(hbs,ann.hrz,{'y'})}];
                        smr.status(2)=0;
                    end
                    if ~ismember({'x'},hbs.mon.cp)&&~ismember({'y'},hbs.mon.cp)
                        warning('ERROR: NO HORIZONTAL COMPONENTS SELECTED!');
                    end
                    
                else
                    warning('ERROR: HORIZONTAL ANN NOT SELECTED!');
                end
                % _apply trained ANN on hybrid accelerograms (vertical)_
                if ann.status(2)
                    if ismember({'z'},hbs.mon.cp)
                        smr.trs = [smr.trs;{ann2hbs_train(hbs,ann.hrz,{'z'})}];
                        smr.status(3)=0;
                    else
                        warning('ERROR: NO VERTICAL COMPONENT SELECTED!');
                    end
                    
                else
                    warning('ERROR: VERTICAL ANN NOT SELECTED!');
                end
                
                disp('5.     SPECTRAL MATCHING: ');
                GUI_sm_setup;
                disp('-----> SPECTRA : MATCHED!');
                set(S.pb(6),'str','MATCHED','backg',col);
                if any(smr.status)
                    close(S.fh);
                end
            end
    end
    return
end
