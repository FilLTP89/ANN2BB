fid=fopen('ppe2012_all_rec_pga.csv','w+');
grd=importdata('grid.xyz');
fprintf(fid,'LON [dd],LAT [dd],Rgh [m/s2],\n');
for i_=1:size(grd.data)
    fprintf(fid,'%20.10f,%20.10f,%20.10f\n',grd.data(i_,1),grd.data(i_,2),...
        grd.data(i_,3).*9.81./100);
end
fclose(fid);