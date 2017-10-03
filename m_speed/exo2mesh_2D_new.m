function exo2mesh_2D_new(varargin)
    tstart=cputime;
    fn.mod=varargin{1};
    fn.inp=strcat(char(fn.mod),'.txt');
    fn.out=strcat(char(fn.mod),'.mesh');
    
    fid = fopen(fn.inp,'r');
    
    prop.tag  = {'num_nodes';'num_elem';'num_el_blk'};
    prop.stag = {'num_el_in_blk';'num_nod_per_el'};

    prop.nmt = numel(prop.tag);
    prop.nms = numel(prop.stag);
    
    cbar=0;
    cshell=0;
    grid=0;
    num_cshell=0;
    num_cbar=0;
    
    for i_=1:prop.nmt
        frewind(fid);
        while 1
            tline = fgetl(fid);
            idx0=strfind(tline,prop.tag{i_});
            idx1=strfind(tline,' ;');
            if ~isempty(idx0)
                prop.(prop.tag{i_})=...
                    str2num(tline(idx0+numel(prop.tag{i_})+3:idx1-1));
                break;
            end
        end
    end
    
    for b_=1:prop.(prop.tag{end})
        for i_=1:prop.nms
            tline = fgetl(fid);
            idx0=strfind(tline,prop.stag{i_});
            idx1=strfind(tline,' ;');
            prop.(prop.stag{i_})(b_,1)=...
                str2num(tline(idx0+numel(prop.stag{i_})+3:idx1-1));
        end
    end
    
    prop.num_cshell = sum(prop.num_el_in_blk(prop.num_nod_per_el==4));
    prop.num_cbar   = sum(prop.num_el_in_blk(prop.num_nod_per_el~=4));
    
    while 1
            tline = fgetl(fid);
            idx0=strfind(tline,'coord =');
            tline = fgetl(fid);
            temp=[];
            if ~isempty(idx0)
                temp=[];
                while 1
                    idx1=strfind(tline,' ;');
                    if isempty(idx1)
                        str = strread(tline, '%s','delimiter',',');
                        temp = [temp; str2double(str)];
                        tline=fgetl(fid);
                    else
                        str = strread(tline(1:end-2), '%s','delimiter',',');
                        temp = [temp; str2double(str)];
                        break;
                    end
                end
                break;
            end
    end
    
    nodes.x = temp(1:prop.num_nodes);
    nodes.y = temp(prop.num_nodes+1:2*prop.num_nodes);
    
%     
%     
%     con_cshell=zeros(num_cshell,4);
%     con_cbar=zeros(num_cbar,2);
%     
%     grid_cord=zeros(num_nodes*3,1);
%     
%     cshell_tag=zeros(num_cshell,1);
%     cshell_id=zeros(num_cshell,1);
%     
%     cbar_tag=zeros(num_cbar,1);
%     cbar_id=zeros(num_cbar,1);
%     
%     grid_cord=zeros(num_nodes,3);
%     grid_id=zeros(num_nodes,1);
%     grid_x=zeros(num_nodes,1);
%     grid_y=zeros(num_nodes,1);
%     grid_z=zeros(num_nodes,1);
%     
%     cshell_1=zeros(num_cshell,1);
%     cshell_2=zeros(num_cshell,1);
%     cshell_3=zeros(num_cshell,1);
%     cshell_4=zeros(num_cshell,1);
%     
%     cbar_1=zeros(num_cbar,1);
%     cbar_2=zeros(num_cbar,1);
%     
%     %attrib=zeros(attrib_sum,1)
%     
%     
%     
%     for i_=1:prop.num_el_blk %+ attrib_sum
%         
%         %Reading attrib
%         %if attrib(i) == 1
%         %    for j=1:num_el_in_blk(i)
%         %        tline = fgetl(fid);
%         %    end
%         %    tline = fgetl(fid);
%         %    tline = fgetl(fid);
%         %end
%         
%         if i==1
%             ['Reading ',char(num2str(i)),'st block  (',char(num2str(num_el_in_blk(i))),' elem.)']
%         elseif i==2
%             ['completed in ',char(num2str(cputime-tini)),' sec.']
%             ['Reading ',char(num2str(i)),'nd block  (',char(num2str(num_el_in_blk(i))),' elem.)']
%         elseif i==3
%             ['completed in ',char(num2str(cputime-tini)),' sec.']
%             ['Reading ',char(num2str(i)),'rd block  (',char(num2str(num_el_in_blk(i))),' elem.)']
%         else
%             ['completed in ',char(num2str(cputime-tini)),' sec.']
%             ['Reading ',char(num2str(i)),'th block  (',char(num2str(num_el_in_blk(i))),' elem.)']
%         end
%         
%         tini=cputime;
%         
%         while strcmp(tline(1:l_tline),' connect')~=1
%             tline = fgetl(fid);
%             if length(tline)<8
%                 l_tline=length(tline);
%             elseif length(tline)==0
%                 l_tline=1;
%             else
%                 l_tline=8;
%             end
%         end
%         
%         if num_nod_per_el(i)==4
%             for j=1:num_el_in_blk(i)
%                 tline = fgetl(fid);
%                 %pos=findstr(tline,',');
%                 cshell=cshell+1;
%                 cshell_id(cshell)=cshell;
%                 cshell_tag(cshell)=i;
%                 con_cshell(cshell,:)=str2num(tline);
%             end
%             
%         else
%             for j=1:num_el_in_blk(i)
%                 tline = fgetl(fid);
%                 %pos=findstr(tline,',');
%                 cbar=cbar+1;
%                 cbar_id(cbar)=cbar;
%                 cbar_tag(cbar)=i;
%                 con_cbar(cbar,:)=str2num(tline);
%             end
%         end
%     end
%     
%     ['completed in ',char(num2str(cputime-tini)),' sec.']
%     
%     tini=cputime;
%     ['*** Storing data informations ***']
%     
%     if cshell>0
%         cshell_1(1:cshell)=con_cshell(1:cshell,1);
%         cshell_2(1:cshell)=con_cshell(1:cshell,2);
%         cshell_3(1:cshell)=con_cshell(1:cshell,3);
%         cshell_4(1:cshell)=con_cshell(1:cshell,4);
%     end
%     if cbar>0
%         cbar_1(1:cbar)=con_cbar(1:cbar,1);
%         cbar_2(1:cbar)=con_cbar(1:cbar,2);
%     end
%     
%     ['completed in ',char(num2str(cputime-tini)),' sec.']
%     
%     tini=cputime;
%     ['BEGIN - Reading nodes coordinates']
%     
%     %clear pos;
%     
%     %cos� non funziona!!!
%     
%     tline = fgetl(fid);
%     tline = fgetl(fid);
%     vero=[];
%     grid=1;
%     
%     while strcmp(tline(1:l_tline),' coord =')~=1
%         %for i= 1:23+2*num_el_blk
%         tline = fgetl(fid);
%         if length(tline)<8
%             l_tline=length(tline);
%         elseif length(tline)==0
%             l_tline=1;
%         else
%             l_tline=8;
%         end
%     end
%     
%     while length(vero)==0
%         %clear pos;
%         tline = fgetl(fid);
%         vero=findstr(tline,';');
%         
%         if length(vero)==0
%             
%             %pos=findstr(tline,',');
%             %pos(2:length(pos)+1)=pos(1:length(pos));
%             %pos(1)=1;
%             %for i=1:length(pos)-1
%             %grid=grid+1;
%             howmanynumb=length(str2num(tline));
%             grid_cord(grid:grid+howmanynumb-1)=str2num(tline);
%             grid=grid+howmanynumb;
%             %grid_cord(grid)=str2num(tline(pos(i)+1:pos(i+1)-1));
%             %end
%             %grid=grid+1;
%             %if i==[]
%             %    grid_cord(grid)=str2num(tline(1:vero));
%             %else
%             %    grid_cord(grid)=str2num(tline(pos(i+1)+1:vero));
%             %end
%         else
%             %pos=findstr(tline,',');
%             %pos(2:length(pos)+1)=pos(1:length(pos));
%             %pos(1)=1;
%             %for i=1:length(pos)-1
%             %grid=grid+1;
%             howmanynumb=length(str2num(tline));
%             grid_cord(grid:grid+howmanynumb-1)=str2num(tline);
%             grid=grid+howmanynumb;
%             %   grid_cord(grid)=str2num(tline(pos(i)+1:pos(i+1)-1));
%             %end
%         end
%     end
%     
%     grid=num_nodes;
%     for i=1:num_nodes
%         grid_id(i)=i;
%     end
%     grid_x(1:num_nodes)=grid_cord(1:num_nodes);
%     grid_y(1:num_nodes)=grid_cord(num_nodes+1:2*num_nodes);
%     grid_z(1:num_nodes)=grid_cord(2*num_nodes+1:3*num_nodes);
%     
%     
%     %if test_jac==1
%     x1=grid_x(cshell_1(1:cshell));
%     x2=grid_x(cshell_2(1:cshell));
%     x3=grid_x(cshell_3(1:cshell));
%     x4=grid_x(cshell_4(1:cshell));
%     
%     y1=grid_y(cshell_1(1:cshell));
%     y2=grid_y(cshell_2(1:cshell));
%     y3=grid_y(cshell_3(1:cshell));
%     y4=grid_y(cshell_4(1:cshell));
%     
%     alfa1(1:cshell) = 0.25d0.*(-x1 +x2 +x3 -x4);
%     beta1(1:cshell) = 0.25d0.*(-x1 -x2 +x3 +x4);
%     alfa2(1:cshell) = 0.25d0.*(-y1 +y2 +y3 -y4);
%     beta2(1:cshell) = 0.25d0.*(-y1 -y2 +y3 +y4);
%     
%     jac(1:cshell) = alfa1.*beta2 - alfa2.*beta1;
% %     plot([1:cshell],jac,'b-');hold on
%     
%     %end
%     
%     for j=1:cshell
%         if jac(j)<0
%             tmp1=cshell_1(j);
%             tmp2=cshell_2(j);
%             tmp3=cshell_3(j);
%             tmp4=cshell_4(j);
%             
%             cshell_1(j)=tmp4;
%             cshell_2(j)=tmp3;
%             cshell_3(j)=tmp2;
%             cshell_4(j)=tmp1;
%         end
%     end
%     
%     x1=grid_x(cshell_1(1:cshell));
%     x2=grid_x(cshell_2(1:cshell));
%     x3=grid_x(cshell_3(1:cshell));
%     x4=grid_x(cshell_4(1:cshell));
%     
%     y1=grid_y(cshell_1(1:cshell));
%     y2=grid_y(cshell_2(1:cshell));
%     y3=grid_y(cshell_3(1:cshell));
%     y4=grid_y(cshell_4(1:cshell));
%     
%     alfa1(1:cshell) = 0.25d0.*(-x1 +x2 +x3 -x4);
%     beta1(1:cshell) = 0.25d0.*(-x1 -x2 +x3 +x4);
%     alfa2(1:cshell) = 0.25d0.*(-y1 +y2 +y3 -y4);
%     beta2(1:cshell) = 0.25d0.*(-y1 -y2 +y3 +y4);
%     
%     jac(1:cshell) = alfa1.*beta2 - alfa2.*beta1;
% %     plot([1:cshell],jac,'r--');hold on
%     
%     
%     
%     fclose(fid);
%     
%     ['END - Reading nodes coordinates in ',char(num2str(cputime-tini)),' sec.']
%     
%     fid = fopen(fn.out,'w');
%     
%     blank_space=length(num2str(grid_id(grid)));
%     %lunghezza_max+3spazi-lunghezza_corrente
%     fprintf(fid,'   %i   %i   %i   %i   %i\n',[grid; cshell+cbar; 0; 0; 0]);
%     
%     tini=cputime;
%     ['BEGIN - Writing inp format']
%     
%     if grid>0
%         fprintf(fid,'%i  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)']);
%     end
%     
%     if cbar>0
%         fprintf(fid,'%i  %i  line  %i  %i\n',...
%             [cbar_id(1:cbar)'; cbar_tag(1:cbar)';...
%             cbar_1(1:cbar)'; cbar_2(1:cbar)']);
%     end
%     
%     if cshell>0
%         fprintf(fid,'%i  %i   quad  %i  %i  %i  %i\n',...
%             [cshell_id(1:cshell)'; cshell_tag(1:cshell)';...
%             cshell_1(1:cshell)'; cshell_2(1:cshell)'; cshell_3(1:cshell)'; cshell_4(1:cshell)']);
%     end
%     
%     ['END - Writing inp format in ',char(num2str(cputime-tini)),' sec.']
%     ['TOTAL CPU time ',char(num2str(cputime-tstart)),' sec.']
%     fclose(fid);
%     
%     %get_func_value -> find out load function on Pelse3d
    
    return
end


