function press_button(varargin)
    pbt = varargin{1};
    str = 'RUNNING...';
    
    if nargin>1
        str = varargin{2};
    end
    
    bcc = get(pbt,'backg');
    rcc = [1.00 0.60 0.60];
    set(pbt,'str',str,'backg',rcc);
    pause(.01);
    set(pbt,'str',str,'backg',bcc);
    
    return
end