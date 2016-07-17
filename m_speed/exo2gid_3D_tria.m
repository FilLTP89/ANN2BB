function exo2gid_3D_tria(pathin,pathout,file_name)

nomefiler=[char(pathin),char(file_name),'.txt'];
nomefilew=[char(pathin),char(file_name),'.post.msh'];
nomefilew2=[char(pathin),char(file_name),'.monitor'];
nomefilew3=[char(pathout),char(file_name),'.out'];

fid = fopen(nomefiler,'r');
    
    disp(['Converting ', file_name,'.txt in ', file_name, '.out'])
    
    tab=7;
    ctria=0;
    grid=0;
    num_ctria=0;
    
    for i= 1:8
        tline = fgetl(fid);
    end 
    
    num_nodes=str2num(tline(20-tab:length(tline)-1));
    tline = fgetl(fid);
    num_elem=str2num(tline(19-tab:length(tline)-1));
    tline = fgetl(fid);
    num_el_blk=str2num(tline(21-tab:length(tline)-1));
    tline = fgetl(fid);
    tline = fgetl(fid);
    
    for i=1:num_el_blk
        num_el_in_blk(i)=str2num(tline(25-tab+fix(i/10):length(tline)-1));
        tline = fgetl(fid);
        num_nod_per_el(i)=str2num(tline(26-tab+fix(i/10):length(tline)-1));
        tline = fgetl(fid);
        if num_nod_per_el(i)==3
            num_ctria = num_ctria +  num_el_in_blk(i);
        end
        tline = fgetl(fid);
    end 
    
    num_ctria;
    
    l_tline=8;
    tline = fgetl(fid);
    i=0;
    
    con_ctria=zeros(num_ctria,3);
    
    grid_cord=zeros(num_nodes*3,1);
    
    ctria_tag=zeros(num_ctria,1);
    ctria_id=zeros(num_ctria,1);
    
    grid_cord=zeros(num_nodes,3);
    grid_id=zeros(num_nodes,1);
    grid_x=zeros(num_nodes,1);
    grid_y=zeros(num_nodes,1);
    grid_z=zeros(num_nodes,1);
    
    ctria_1=zeros(num_ctria,1);
    ctria_2=zeros(num_ctria,1);
    ctria_3=zeros(num_ctria,1);
       
    
    
    for i=1:num_el_blk
        
        if i==1
            tline = fgetl(fid);
            disp(['Reading ',char(num2str(i)),'st block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
        elseif i==2
            disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
            disp(['Reading ',char(num2str(i)),'nd block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
        elseif i==3
            disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
            disp(['Reading ',char(num2str(i)),'rd block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
        else
            disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
            disp(['Reading ',char(num2str(i)),'th block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
        end
        
        tini=cputime;
        
        while strcmp(tline(1:l_tline),' connect')~=1
            tline = fgetl(fid);
            if length(tline)<8
                l_tline=length(tline);
            elseif length(tline)==0
                l_tline=1;
            else
                l_tline=8;
            end
        end
        
        if num_nod_per_el(i)==3
            for j=1:num_el_in_blk(i)
                tline = fgetl(fid);
                %pos=findstr(tline,',');
                ctria=ctria+1;
                ctria_id(ctria)=ctria;
                ctria_tag(ctria)=i;
                con_ctria(ctria,:)=str2num(tline);
            end
        end 
    end
    
    disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
    
    tini=cputime;
    disp(['*** Storing data informations ***'])
    
    if ctria>0
        ctria_1(1:ctria)=con_ctria(1:ctria,1);
        ctria_2(1:ctria)=con_ctria(1:ctria,2);
        ctria_3(1:ctria)=con_ctria(1:ctria,3);
    end
    
    disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
    
    tini=cputime;
    disp(['BEGIN - Reading nodes coordinates'])
    
    %clear pos;
        
    tline = fgetl(fid);
    tline = fgetl(fid);
    vero=[];
    grid=1;
    
    while strcmp(tline(1:l_tline),' coord =')~=1
        %for i= 1:23+2*num_el_blk
            tline = fgetl(fid);
            if length(tline)<8
                l_tline=length(tline);
            elseif length(tline)==0
                l_tline=1;
            else
                l_tline=8;
            end
        end
    
    while length(vero)==0
        %clear pos;
        tline = fgetl(fid);
        vero=findstr(tline,';');
        
        if length(vero)==0
            
            %pos=findstr(tline,',');
            %pos(2:length(pos)+1)=pos(1:length(pos));
            %pos(1)=1;
            %for i=1:length(pos)-1
                %grid=grid+1;
                howmanynumb=length(str2num(tline));
                grid_cord(grid:grid+howmanynumb-1)=str2num(tline);
                grid=grid+howmanynumb;
                %grid_cord(grid)=str2num(tline(pos(i)+1:pos(i+1)-1));
                %end
            %grid=grid+1;
            %if i==[]
            %    grid_cord(grid)=str2num(tline(1:vero));
            %else
            %    grid_cord(grid)=str2num(tline(pos(i+1)+1:vero));
            %end
        else
            %pos=findstr(tline,',');
            %pos(2:length(pos)+1)=pos(1:length(pos));
            %pos(1)=1;
            %for i=1:length(pos)-1
                %grid=grid+1;
                howmanynumb=length(str2num(tline));
                grid_cord(grid:grid+howmanynumb-1)=str2num(tline);
                grid=grid+howmanynumb;
             %   grid_cord(grid)=str2num(tline(pos(i)+1:pos(i+1)-1));
             %end
        end
    end
    
    grid=num_nodes;
    for i=1:num_nodes
        grid_id(i)=i;
    end
    grid_x(1:num_nodes)=grid_cord(1:num_nodes);
    grid_y(1:num_nodes)=grid_cord(num_nodes+1:2*num_nodes);
    grid_z(1:num_nodes)=grid_cord(2*num_nodes+1:3*num_nodes);
    
    fclose(fid);

    
    disp(['END - Reading nodes coordinates in ',char(num2str(cputime-tini)),' sec.'])


    fid = fopen(nomefilew,'w');
    
    blank_space=length(num2str(grid_id(grid)));
    %fprintf(fid,' line 1: output of bending\n');
    %fprintf(fid,' line 2: output of bending\n');
    %fprintf(fid,' line 3: output of bending\n');
    %fprintf(fid,' line 4: output of bending\n');
    %fprintf(fid,' line 5: output of bending\n');
    %fprintf(fid,' n_bound_elements   n_bound_points   nelem_type\n');
    %fprintf(fid,'   %i   %i   7\n',[ctria; grid]);
    %fprintf(fid,' --- coordenadas ---\n');
    fprintf(fid,'#mesh of a table\n');
    fprintf(fid,'MESH "model" Dimension 3 ElemType Triangle Nnode 3\n');
    fprintf(fid,'Coordinates\n');
    fprintf(fid,'# node_number coordinate_x coordinate_y coordinate_z\n'); 
% fprintf(fid,'    %i      %i\n',[num_nodes; ctria]);

    
    
    tini=cputime;
    disp(['BEGIN - Writing input format'])
    
    if grid>0
            fprintf(fid,'%i  %+13.7e  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
    end

    %fprintf(fid,'  --- ielem    n1   n2\n');
    fprintf(fid,'End Coordinates\n');
    fprintf(fid,'Elements\n');
    fprintf(fid,'# element node_1 node_2 node_3 material_number\n');
    if ctria>0
            fprintf(fid,'%i    %i  %i  %i    %i\n',...
                        [ctria_id(1:ctria)';...
                        ctria_1(1:ctria)'; ctria_2(1:ctria)'; ctria_3(1:ctria)'; ctria_tag(1:ctria)']);
    end
   fprintf(fid,'End Elements\n');
    

fclose(fid);



fid = fopen(nomefilew2,'w'); 
    fprintf(fid,'MONITOR FILES\n');
    
    if grid>0
            fprintf(fid,'MONITOR  %+13.7e  %+13.7e  %+13.7e\n',[ grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
    end
    fprintf(fid,'   \n'); 
fclose(fid);

fid = fopen(nomefilew3,'w'); 
%     fprintf(fid,'%i \n',num_nodes);
    fprintf(fid,'    %i      %i\n',[num_nodes; ctria]);
    if grid>0
            fprintf(fid,'%i  %+13.7e  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
    end
    if ctria>0
        fprintf(fid,'%i    %i  %i  %i    %i\n',...
                [ctria_id(1:ctria)';...
                ctria_1(1:ctria)'; ctria_2(1:ctria)'; ctria_3(1:ctria)'; ctria_tag(1:ctria)']);
    end
fclose(fid);






