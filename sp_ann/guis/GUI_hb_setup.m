function [] = GUI_sp_setup()
    %% NUMERICAL SIMULATIONS
    global mon nss sps hbs Shbs
    
    Shbs.fh = figure('units','pixels',...
        'position',[0 350 305 70],...
        'menubar','none',...
        'name','LF/HF Hybridize',...
        'numbertitle','off',...
        'resize','off');
    Shbs.pb(1) = uicontrol('style','push',...
        'units','pix',...
        'position',[5 10 100 50],...
        'fontsize',14,...
        'string','MASHUP');
    Shbs.pb(2) = uicontrol('style','push',...
        'units','pix',...
        'position',[105 10 100 50],...
        'fontsize',14,...
        'string','SPECTRA');
    Shbs.pb(3) = uicontrol('style','push',...
        'units','pix',...
        'position',[205 10 100 50],...
        'fontsize',14,...
        'string','CLOSE');
    
    set(Shbs.pb(:),'callback',{@pb_call,Shbs})  % Set callbacks.
    Shbs.TF = false;  % Flag for stopping the loop.
    
    uiwait(Shbs.fh) % Wait for continue or stop button.
        
    return
end

function [] = pb_call(varargin)
    global mon nss sps hbs Shbs
    % Callback for pushbutton.
    if varargin{1} == Shbs.pb(1)  % Get the structure.
        %%
        % _parse simulations_
        hbs = lfhf_mashup(nss,sps);
        Shbs.flag(1) = 1;
    elseif varargin{1} == Shbs.pb(2)
        %%
        % _plot spectra_
%         [TODO]
%         sps = sp_spectra(sps);
        disp('SPECTRA COMPUTED!');
    elseif varargin{1} == Shbs.pb(3)
        %%
        % _exit_
        hbs.status=Shbs.flag(1);
        close(Shbs.fh)  % Found the one we are looking for.
    end
    
    return
end
