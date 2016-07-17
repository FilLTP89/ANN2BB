ccc;
cfolder=sprintf('%s%s',cd,pslash); fault_code=cfolder(end-7:end-1);
ffolder=sprintf('%sFAULT%s',cfolder,pslash);
mfolder=sprintf('%sMAP%s',cfolder,pslash);
ffault=dir(sprintf('%s*csv',ffolder));
fmap=dir(sprintf('%s*csv',mfolder));
bessel=referenceEllipsoid('Bessel 1841','meters');
mstruct = defaultm('mercator');
mstruct.geoid=[bessel.SemimajorAxis bessel.Eccentricity];
mstruct = defaultm(mstruct);
for i=1:numel(ffault)
    str=sprintf('%s%s',ffolder,ffault(i).name);
    fc=read_fault_coordinates(str);
    % draw planes
    fig=figure('name',sprintf('%s_plane_%u',fault_code,i));
    hax=axes('parent',fig,'dataaspectratio',[1 1 1]); hold(hax,'all');
    [x,y,z]=mfwdtran(mstruct,fc.fault_planes.jlat,fc.fault_planes.jlon,-fc.fault_planes.depth);
    patch(x,y,z,zeros(size(z)),'parent',hax,...
        'facealpha',.5,'facecolor','flat',...
        'linestyle','-.');
    set(hax,'xtick',[min(min(x)),max(max(x))]);
    set(hax,'ytick',[min(min(y)),max(max(y))]);
    set(hax,'ztick',[min(min(z)),max(max(z))]);
    xlim(hax,[min(min(x)),max(max(x))]);
    ylim(hax,[min(min(y)),max(max(y))]);
    zlim(hax,[min(min(z)),max(max(z))]);
    view(hax,3); format_figures(hax);
    set(hax,'fontsize',12,'fontweight','bold');
    % draw asperities
    [x,y,z]=mfwdtran(mstruct,fc.fault_asperities.jlat,fc.fault_asperities.jlon,-fc.fault_asperities.depth);
    patch(x,y,z,ones(size(z)),'parent',hax,...
        'facealpha',1,'facecolor','flat',...
        'linestyle','-.');
end