function [] = GUI_sp_setup()
    %% NUMERICAL SIMULATIONS
    global mon nss sps hbs Shbs
    Shbs.flag=0;
    Shbs.fh = figure('units','pixels',...
        'position',[0 350 235 60],...
        'menubar','none',...
        'name','LF/HF Hybridize',...
        'numbertitle','off',...
        'resize','off');
    Shbs.pb(1) = uicontrol('style','push',...
        'units','pix',...
        'position',[5 5 110 50],...
        'string','MASHUP');
    Shbs.pb(2) = uicontrol('style','push',...
        'units','pix',...
        'position',[120 5 110 50],...
        'string','CLOSE');
    
    set(Shbs.pb(:),'callback',{@pb_call,Shbs})  % Set callbacks.
    Shbs.TF = false;  % Flag for stopping the loop.
    
    uiwait(Shbs.fh) % Wait for continue or stop button.
    
    return
end

function [] = pb_call(varargin)
    global mon nss sps hbs Shbs
    % Callback for pushbutton.
    press_button(varargin{1});
    if varargin{1} == Shbs.pb(1)  % Get the structure.
        %%
        % _parse simulations_
        hbs = lfhf_mashup(nss,sps);
        Shbs.flag(1) = 1;
        press_button(varargin{1},'OK!');
    elseif varargin{1} == Shbs.pb(2)
        %%
        % _exit_
        hbs.status=Shbs.flag(1);
        close(Shbs.fh)  % Found the one we are looking for.
    end
    return
end
