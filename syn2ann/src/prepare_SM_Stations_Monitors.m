ccc;
datum=importdata('Monitors.csv');
data=-999*ones(size(datum.data,1),6);
data(:,1)=datum.data(:,1);
data(:,2)=datum.data(:,7);
data(:,3)=datum.data(:,4);
data(:,4)=datum.data(:,5);
data(:,5)=datum.data(:,2);
data(:,6)=datum.data(:,3);
fid=fopen('SM_Stations_Monitors.csv','w+');
fprintf(fid,'Station Code,Monitor ID,Repi [km],E_UTM [m],N_UTM [m],LON [dd],LAT [dd],\n');
for i_=1:size(data,1)
    fprintf(fid,'MON%u,%u,%u,%u,%u,%u,%u,\n',i_,data(i_,1),data(i_,2),data(i_,3),...
        data(i_,4),data(i_,5),data(i_,6));
end
fclose(fid);

fid=fopen('site_SP96.csv','w+');
for i_=1:size(data,1)
    fprintf(fid,'%u,\n',datum.data(i_,end));
end
fclose(fid);