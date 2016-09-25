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
    
    Child   = get(hax, 'Children');
    XTick   = get(hax, 'XTick');
    YTick   = get(hax, 'YTick');
    XLimit  = get(hax, 'XLim');
    YLimit  = get(hax, 'YLim');
    newGrid = cat(1, ...
        line([XTick; XTick], YLimit, 'Parent', hax, GridStyle), ...
        line(XLimit, [YTick; YTick], 'Parent', hax, GridStyle));
    
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