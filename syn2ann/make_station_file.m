function make_station_file(varargin)
    
    file = varargin{1};
    station_list=importdata(file);
    fid=fopen('stations.txt','w+');
    
    coord = station_list.data;
    name  = station_list.textdata;
    
    for i=1:size(coord,1)
        
        if coord(i,3)==0
            fprintf(fid,'%10.4f%14.4f%14.4f\r\n',coord(i,:));
        else
            fprintf(fid,'%10.4f%14.4f%14.4f\r\n',[coord(i,1:2) -coord(i,3)]);
            
        end
        
    end
    fclose(fid);
    return
end