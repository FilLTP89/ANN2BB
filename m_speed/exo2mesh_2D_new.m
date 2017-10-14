function exo2mesh_2D_new(varargin)
    
    fn.mod=varargin{1};
    fn.inp=strcat(char(fn.mod),'.txt');
    fn.out=strcat(char(fn.mod),'.mesh');
    
    fid.in = fopen(fn.inp,'r');
    
    prop.tag  = {'num_nodes';'num_elem';'num_el_blk'};
    prop.stag = {'num_el_in_blk';'num_nod_per_el'};
    
    prop.nmt = numel(prop.tag);
    prop.nms = numel(prop.stag);
    
    for i_=1:prop.nmt
        frewind(fid.in);
        while 1
            tline = fgetl(fid.in);
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
            tline = fgetl(fid.in);
            idx0=strfind(tline,prop.stag{i_});
            idx1=strfind(tline,' ;');
            prop.(prop.stag{i_})(b_,1)=...
                str2num(tline(idx0+numel(prop.stag{i_})+3:idx1-1));
        end
    end
    
    prop.num_cshell = sum(prop.num_el_in_blk(prop.num_nod_per_el==4));
    prop.num_cbar   = sum(prop.num_el_in_blk(prop.num_nod_per_el~=4));
    
    while 1
        
        idx0=strfind(tline,'coord =');
        tline = fgetl(fid.in);
        temp=[];
        if ~isempty(idx0)
            temp=[];
            while 1
                idx1=strfind(tline,' ;');
                if isempty(idx1)
                    str = strread(tline, '%s','delimiter',',');
                    temp = [temp; str2double(str)];
                    tline=fgetl(fid.in);
                else
                    str = strread(tline(1:end-2), '%s','delimiter',',');
                    temp = [temp; str2double(str)];
                    break;
                end
            end
            break;
        end
    end
    
    nodes.coords.x = temp(1:prop.num_nodes);
    nodes.coords.y = temp(prop.num_nodes+1:2*prop.num_nodes);
    
    idxc=0;
    idxb=0;
    while 1
        tline = fgetl(fid.in);
        idx0=strfind(tline,'connect');
        if ~isempty(idx0)
            idx0=str2num(tline(idx0+7:strfind(tline,'=')-2));
            
            if prop.num_nod_per_el(idx0)==4
                idxc=idxc+1;
                elems.shell(idxc).tag = idx0;
                for e_=1:prop.num_el_in_blk(idx0)
                    tline=fgetl(fid.in);
                    idx1=strfind(tline,' ;');
                    if isempty(idx1)
                        temp = str2double(strread(tline, '%s','delimiter',','));
                    else
                        temp = str2double(strread(tline(1:end-2), '%s','delimiter',','));
                    end
                    elems.shell(idxc).con{e_} = temp;
                end
            elseif prop.num_nod_per_el(idx0)==2
                idxb=idxb+1;
                elems.bar(idxb).tag = idx0;
                for e_=1:prop.num_el_in_blk(idx0)
                    tline=fgetl(fid.in);
                    idx1=strfind(tline,' ;');
                    if isempty(idx1)
                        temp = str2double(strread(tline, '%s','delimiter',','));
                    else
                        temp = str2double(strread(tline(1:end-2), '%s','delimiter',','));
                    end
                    elems.bar(idxb).con{e_} = temp;
                end
            end
        else
            temp=strfind(tline,'qa');
            if ~isempty(temp)
                break;
            end
        end
    end
    
    %% Compute Jacobian
    for i_=1:prop.num_cshell
        
        xn=nodes.coords.x(elems.shell.con{i_});
        yn=nodes.coords.y(elems.shell.con{i_});
        elem.shell.alfa1(i_) = 0.25d0.*(-xn(1) +xn(2) +xn(3) -xn(4));
        elem.shell.beta1(i_) = 0.25d0.*(-xn(1) -xn(2) +xn(3) +xn(4));
        elem.shell.alfa2(i_) = 0.25d0.*(-yn(1) +yn(2) +yn(3) -yn(4));
        elem.shell.beta2(i_) = 0.25d0.*(-yn(1) -yn(2) +yn(3) +yn(4));
        elem.shell.jac(i_) = elem.shell.alfa1(i_).*elem.shell.beta2(i_) - ...
            elem.shell.alfa2(i_).*elem.shell.beta1(i_);
        if elem.shell.jac(i_)<0
            error('negative jacobian: element %u',i_)
        end
    end
    fclose(fid.in);
    
    %% Write speed input file
    fid.out = fopen(fn.out,'w+');
    keyboard
    fprintf(fid.out,'%15u%15u%15u%15u%15u\n',prop.num_nodes,...
        prop.num_cshell+prop.num_cbar,0,0,0);
%     for i_=1:prop.num_nodes
%         fprintf(fid.out,'%15u%15.7e%15.7e\n',
%     end
    %     if grid>0
    %         fprintf(fid.in,'%i  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)']);
    %     end
    %
    %     if cbar>0
    %         fprintf(fid.in,'%i  %i  line  %i  %i\n',...
    %             [cbar_id(1:cbar)'; cbar_tag(1:cbar)';...
    %             cbar_1(1:cbar)'; cbar_2(1:cbar)']);
    %     end
    %
    %     if cshell>0
    %         fprintf(fid.in,'%i  %i   quad  %i  %i  %i  %i\n',...
    %             [cshell_id(1:cshell)'; cshell_tag(1:cshell)';...
    %             cshell_1(1:cshell)'; cshell_2(1:cshell)'; cshell_3(1:cshell)'; cshell_4(1:cshell)']);
    %     end
    %
    %     ['END - Writing inp format in ',char(num2str(cputime-tini)),' sec.']
    %     ['TOTAL CPU time ',char(num2str(cputime-tstart)),' sec.']
    %     fclose(fid.in);
    %
    %     %get_func_value -> find out load function on Pelse3d
    
    return
end


