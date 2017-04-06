function set_axis_grid(varargin)
    if nargin>0
        hax = varargin{1};
    else
        hax = gca;
    end
    
    % The modifications start here:
    GridStyle.Color     = rgb('lightgrey');
    GridStyle.LineStyle = '-';
    GridStyle.LineWidth = 1;
    GridStyle.HitTest   = 'off';
    
    Child   = get(hax,'Children');
    XTick   = get(hax,'XTick');
    YTick   = get(hax,'YTick');
    XLimit  = get(hax,'XLim');
    YLimit  = get(hax,'YLim');
    
    if numel(XTick)>3
        idx.x = 2:numel(XTick)-1;
    else
        idx.x = 1:numel(XTick);
    end
    if numel(YTick)>3
        idx.y = 2:numel(YTick)-1;
    else
        idx.y = 1:numel(YTick);
    end
    newGrid = cat(1, ...
        line([XTick(idx.x);XTick(idx.x)],YLimit,'Parent',hax,GridStyle), ...
        line(XLimit,[YTick(idx.y);YTick(idx.y)],'Parent',hax,GridStyle));
    %% FOR MATLAB 2017A
    for i_=1:numel(newGrid)
        newGrid(i_).LineStyle='none';
        newGrid(i_).DisplayName = '';
    end
    
    % New grid on top or bottom of other objects:
    %set(hax, 'Child', [newGrid; Child(:)]);
    
    set(hax, 'Child', [Child(:); newGrid]);
    
    % Disable original dashed grid:
    set(hax, ...
        'XGrid',      'off', ...
        'YGrid',      'off', ...
        'XMinorGrid',      'off', ...
        'YMinorGrid',      'off', ...
        'YTickMode',  'manual', ...
        'TickLength', zeros(1, 2));
    return
end