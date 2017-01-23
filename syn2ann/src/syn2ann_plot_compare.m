function syn2ann_plot_compare(varargin)
    global pfg xlm xlb xtk ylm ylb ytk grd scl utd tit
    
    %% *SET-UP*
    flags = logical(varargin{1});
    xplot = varargin{2};
    identity = varargin{3};
    dlg = varargin{4};
    legplot = varargin{5};
    cpp = varargin{6};
    stdd = varargin{7};
    fnm = varargin{8};
    time_shift = varargin{9};
    mod = false;
    if nargin>10
        mod = true;
        mrkd = varargin{10};
        lstd = varargin{11};
        lwdd = varargin{12};
    end
    
    spg = numel(xplot);
    if any(flags(3:end))
        pfg.tha.c = (spg==1)*(pfg.tha.s)+(spg==2)*(pfg.tha.d)+(spg==3)*(pfg.tha.t);
        pfg.thv.c = pfg.tha.c;
        pfg.thd.c = pfg.tha.c;
    end
    spg = [spg,1];
    pax = num2cell([1:spg(1),1:spg(1)]');
    
    std(1,1) = {strcat(stdd,'-',upper(cpp))};
    for i_ = 1:spg-1
        std{1,1+i_} = {strcat(stdd,'-',upper(cpp))};
    end
    %% *PSA SPECTRUM*
    
    if flags(1)
        syn2ann_plot_psa;
        close all;
    end
    %% *FOURIER SPECTRUM*
    if flags(2)
        syn2ann_plot_fsa;
        close all;
    end
    
    %% *ACCELERATION TIME-HISTORY*
    if flags(3)
        syn2ann_plot_tha;
        close all;
    end
    %% *VELOCITY TIME-HISTORY*
    if flags(4)
        syn2ann_plot_thv;
        close all;
    end
    %% *DISPLACEMENT TIME-HISTORY*
    if flags(5)
        syn2ann_plot_thd;
        close all;
    end
    %
    return
    
end
