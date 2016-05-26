function [] = GUI_ns_setup()
    %% NUMERICAL SIMULATIONS
    global mon nss Snss
    Snss.flag=0;
    Snss.fh = figure('units','pixels',...
        'position',[0 350 360 60],...
        'menubar','none',...
        'name','Parse Numerical Simulations',...
        'numbertitle','off',...
        'resize','off');
    Snss.pb(1) = uicontrol('style','push',...
        'units','pix',...
        'position',[5 5 110 50],...
        'string','PARSE');
    Snss.pb(2) = uicontrol('style','push',...
        'units','pix',...
        'position',[120 5 110 50],...
        'string','SPECTRA');
    Snss.pb(3) = uicontrol('style','push',...
        'units','pix',...
        'position',[235 5 110 50],...
        'string','CLOSE');
    
    set(Snss.pb(:),'callback',{@pb_call,Snss})  % Set callbacks.
    Snss.TF = false;  % Flag for stopping the loop.
    
    uiwait(Snss.fh) % Wait for continue or stop button.
    
    return
end

function [] = pb_call(varargin)
    global mon nss Snss
    % Callback for pushbutton.
    press_button(varargin{1});
    switch varargin{1}
        case Snss.pb(1)  % Get the structure.
            %%
            % _parse simulations_
            [mon,nss] = ns_parser(mon);
            Snss.flag(1) = 1;
            fprintf('PATH: %s\n',mon.fn);
            for i_=1:mon.na
                fprintf('ID: %u\n',mon.id(i_));
            end
        case Snss.pb(2)
            %%
            % _compute spectra_
            nss = ns_spectra(nss);
            disp('SPECTRA COMPUTED!');
        case Snss.pb(3)
            %%
            % _exit_
            nss.status=all(Snss.flag);
            close(Snss.fh)  % Found the one we are looking for.
    end
    switch varargin{1}
        case {Snss.pb(1),Snss.pb(2)}
            press_button(varargin{1},'OK!');
    end
    return
end
