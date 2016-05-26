function [] = GUI_ns_setup()
    %% NUMERICAL SIMULATIONS
    global mon nss Snss
    
    Snss.fh = figure('units','pixels',...
        'position',[0 350 305 70],...
        'menubar','none',...
        'name','Parse Numerical Simulations',...
        'numbertitle','off',...
        'resize','off');
    Snss.pb(1) = uicontrol('style','push',...
        'units','pix',...
        'position',[5 10 100 50],...
        'fontsize',14,...
        'string','PARSE');
    Snss.pb(2) = uicontrol('style','push',...
        'units','pix',...
        'position',[105 10 100 50],...
        'fontsize',14,...
        'string','SPECTRA');
    Snss.pb(3) = uicontrol('style','push',...
        'units','pix',...
        'position',[205 10 100 50],...
        'fontsize',14,...
        'string','CLOSE');
    
    set(Snss.pb(:),'callback',{@pb_call,Snss})  % Set callbacks.
    Snss.TF = false;  % Flag for stopping the loop.
    
    uiwait(Snss.fh) % Wait for continue or stop button.
        
    return
end

function [] = pb_call(varargin)
    global mon nss Snss 
    % Callback for pushbutton.
    if varargin{1} == Snss.pb(1)  % Get the structure.
        %%
        % _parse simulations_
        [mon,nss] = ns_parser(mon);
        Snss.flag(1) = 1;
        fprintf('PATH: %s\n',mon.fn);
        for i_=1:mon.na
            fprintf('ID: %u\n',mon.id(i_));
        end
    elseif varargin{1} == Snss.pb(2)
        %%
        % _compute spectra_
        nss = ns_spectra(nss);
%         Snss.flag(2) = 1;
        disp('SPECTRA COMPUTED!');
    elseif varargin{1} == Snss.pb(3)
        %%
        % _exit_
        nss.status=all(Snss.flag);
        close(Snss.fh)  % Found the one we are looking for.
    end
    
    return
end
