function exo2ucdx_3D(pathin,pathout,filename)
% *******************************************************
% *******************************************************
% ***                FILTER 3D EXODUSII -> MESH       ***
% ***                  exo2ucdx_3D_v6gamma.m          ***
% *******************************************************

nomefiler=[char(pathin),char(filename),'.txt'];
nomefilew=[char(pathout),char(filename),'.mesh'];

fid = fopen(nomefiler,'r');
disp(['Converting ', filename,'.txt in ', filename, '.mesh'])

tab=7;
cquad=0;
chexa=0;
grid=0;
num_chexa=0;
num_cquad=0;

for i= 1:8
    tline = fgetl(fid);
end

num_nodes=str2num(tline(20-tab:length(tline)-1)); %read nbnodes
tline = fgetl(fid);
num_elem=str2num(tline(19-tab:length(tline)-1)); %read nbelem
tline = fgetl(fid);
num_el_blk=str2num(tline(21-tab:length(tline)-1)); % read nbblocks
tline = fgetl(fid);
tline = fgetl(fid);

j=1;
for i=1:num_el_blk
    string = tline(2:8);
    if(strcmp(string,'num_el_'));
        num_el_in_blk(j)=str2num(tline(25-tab+fix(j/10):length(tline)-1));
        tline = fgetl(fid);
        string = tline(2:8);
    end
    if(strcmp(string,'num_nod'));
        num_nod_per_el(j)=str2num(tline(26-tab+fix(j/10):length(tline)-1));
        if num_nod_per_el(j)==8
            num_chexa = num_chexa +  num_el_in_blk(j);
        else
            num_cquad = num_cquad +  num_el_in_blk(j);
        end
        j = j + 1;
        tline = fgetl(fid);
        string = tline(2:8);
    end
    if(strcmp(string,'num_att'))
        tline = fgetl(fid);
        string = tline(2:8);
    end
    
end

num_chexa;
num_cquad;

l_tline=8;
tline = fgetl(fid);
i=0;

con_chexa=zeros(num_chexa,8);
con_cquad=zeros(num_cquad,4);

grid_cord=zeros(num_nodes*3,1);

chexa_tag=zeros(num_chexa,1);
chexa_id=zeros(num_chexa,1);

cquad_tag=zeros(num_cquad,1);
cquad_id=zeros(num_cquad,1);

grid_cord=zeros(num_nodes,3);
grid_id=zeros(num_nodes,1);
grid_x=zeros(num_nodes,1);
grid_y=zeros(num_nodes,1);
grid_z=zeros(num_nodes,1);

chexa_1=zeros(num_chexa,1);
chexa_2=zeros(num_chexa,1);
chexa_3=zeros(num_chexa,1);
chexa_4=zeros(num_chexa,1);
chexa_5=zeros(num_chexa,1);
chexa_6=zeros(num_chexa,1);
chexa_7=zeros(num_chexa,1);
chexa_8=zeros(num_chexa,1);

cquad_1=zeros(num_cquad,1);
cquad_2=zeros(num_cquad,1);
cquad_3=zeros(num_cquad,1);
cquad_4=zeros(num_cquad,1);



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
    
    if num_nod_per_el(i)==8
        for j=1:num_el_in_blk(i)
            tline = fgetl(fid);
            %pos=findstr(tline,',');
            chexa=chexa+1;
            chexa_id(chexa)=chexa;
            chexa_tag(chexa)=i;
            con_chexa(chexa,:)=str2num(tline);
        end
        
    else
        for j=1:num_el_in_blk(i)
            tline = fgetl(fid);
            %pos=findstr(tline,',');
            cquad=cquad+1;
            cquad_id(cquad)=cquad;
            cquad_tag(cquad)=i;
            con_cquad(cquad,:)=str2num(tline);
        end
    end
end

disp(['completed in ',char(num2str(cputime-tini)),' sec.'])

tini=cputime;
disp(['*** Storing data informations ***'])

if chexa>0
    chexa_1(1:chexa)=con_chexa(1:chexa,1);
    chexa_2(1:chexa)=con_chexa(1:chexa,2);
    chexa_3(1:chexa)=con_chexa(1:chexa,3);
    chexa_4(1:chexa)=con_chexa(1:chexa,4);
    chexa_5(1:chexa)=con_chexa(1:chexa,5);
    chexa_6(1:chexa)=con_chexa(1:chexa,6);
    chexa_7(1:chexa)=con_chexa(1:chexa,7);
    chexa_8(1:chexa)=con_chexa(1:chexa,8);
end
if cquad>0
    cquad_1(1:cquad)=con_cquad(1:cquad,1);
    cquad_2(1:cquad)=con_cquad(1:cquad,2);
    cquad_3(1:cquad)=con_cquad(1:cquad,3);
    cquad_4(1:cquad)=con_cquad(1:cquad,4);
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


%---------------------------------------------------------------
%if test_jac==1
%x1=grid_x(chexa_1(1:chexa));
%x2=grid_x(chexa_2(1:chexa));
%x3=grid_x(chexa_3(1:chexa));
%x4=grid_x(chexa_4(1:chexa));
%x5=grid_x(chexa_5(1:chexa));
%x6=grid_x(chexa_6(1:chexa));
%x7=grid_x(chexa_7(1:chexa));
%x8=grid_x(chexa_8(1:chexa));

%y1=grid_y(chexa_1(1:chexa));
%y2=grid_y(chexa_2(1:chexa));
%y3=grid_y(chexa_3(1:chexa));
%y4=grid_y(chexa_4(1:chexa));
%y5=grid_y(chexa_5(1:chexa));
%y6=grid_y(chexa_6(1:chexa));
%y7=grid_y(chexa_7(1:chexa));
%y8=grid_y(chexa_8(1:chexa));

%z1=grid_z(chexa_1(1:chexa));
%z2=grid_z(chexa_2(1:chexa));
%z3=grid_z(chexa_3(1:chexa));
%z4=grid_z(chexa_4(1:chexa));
%z5=grid_z(chexa_5(1:chexa));
%z6=grid_z(chexa_6(1:chexa));
%z7=grid_z(chexa_7(1:chexa));
%z8=grid_z(chexa_8(1:chexa));



x8=grid_x(chexa_1(1:chexa));
x7=grid_x(chexa_2(1:chexa));
x6=grid_x(chexa_3(1:chexa));
x5=grid_x(chexa_4(1:chexa));
x4=grid_x(chexa_5(1:chexa));
x3=grid_x(chexa_6(1:chexa));
x2=grid_x(chexa_7(1:chexa));
x1=grid_x(chexa_8(1:chexa));

y8=grid_y(chexa_1(1:chexa));
y7=grid_y(chexa_2(1:chexa));
y6=grid_y(chexa_3(1:chexa));
y5=grid_y(chexa_4(1:chexa));
y4=grid_y(chexa_5(1:chexa));
y3=grid_y(chexa_6(1:chexa));
y2=grid_y(chexa_7(1:chexa));
y1=grid_y(chexa_8(1:chexa));

z8=grid_z(chexa_1(1:chexa));
z7=grid_z(chexa_2(1:chexa));
z6=grid_z(chexa_3(1:chexa));
z5=grid_z(chexa_4(1:chexa));
z4=grid_z(chexa_5(1:chexa));
z3=grid_z(chexa_6(1:chexa));
z2=grid_z(chexa_7(1:chexa));
z1=grid_z(chexa_8(1:chexa));


%---
alfa11(1:chexa) = 0.125d0*(-x1 +x2 +x3 -x4 -x5 +x6 +x7 -x8);
alfa21(1:chexa) = 0.125d0*(-y1 +y2 +y3 -y4 -y5 +y6 +y7 -y8);
alfa31(1:chexa) = 0.125d0*(-z1 +z2 +z3 -z4 -z5 +z6 +z7 -z8);
alfa12(1:chexa) = 0.125d0*(-x1 -x2 +x3 +x4 -x5 -x6 +x7 +x8);
alfa22(1:chexa) = 0.125d0*(-y1 -y2 +y3 +y4 -y5 -y6 +y7 +y8);
alfa32(1:chexa) = 0.125d0*(-z1 -z2 +z3 +z4 -z5 -z6 +z7 +z8);
alfa13(1:chexa) = 0.125d0*(-x1 -x2 -x3 -x4 +x5 +x6 +x7 +x8);
alfa23(1:chexa) = 0.125d0*(-y1 -y2 -y3 -y4 +y5 +y6 +y7 +y8);
alfa33(1:chexa) = 0.125d0*(-z1 -z2 -z3 -z4 +z5 +z6 +z7 +z8);

beta11(1:chexa) = 0.125d0*(+x1 +x2 -x3 -x4 -x5 -x6 +x7 +x8);
beta21(1:chexa) = 0.125d0*(+y1 +y2 -y3 -y4 -y5 -y6 +y7 +y8);
beta31(1:chexa) = 0.125d0*(+z1 +z2 -z3 -z4 -z5 -z6 +z7 +z8);
beta12(1:chexa) = 0.125d0*(+x1 -x2 -x3 +x4 -x5 +x6 +x7 -x8);
beta22(1:chexa) = 0.125d0*(+y1 -y2 -y3 +y4 -y5 +y6 +y7 -y8);
beta32(1:chexa) = 0.125d0*(+z1 -z2 -z3 +z4 -z5 +z6 +z7 -z8);
beta13(1:chexa) = 0.125d0*(+x1 -x2 +x3 -x4 +x5 -x6 +x7 -x8);
beta23(1:chexa) = 0.125d0*(+y1 -y2 +y3 -y4 +y5 -y6 +y7 -y8);
beta33(1:chexa) = 0.125d0*(+z1 -z2 +z3 -z4 +z5 -z6 +z7 -z8);

%gamma1(1:chexa) = 0.125d0*(-x1 +x2 -x3 +x4 +x5 -x6 +x7 -x8);
%gamma2(1:chexa) = 0.125d0*(-y1 +y2 -y3 +y4 +y5 -y6 +y7 -y8);
%gamma3(1:chexa) = 0.125d0*(-z1 +z2 -z3 +z4 +z5 -z6 +z7 -z8);

%delta1(1:chexa) = 0.125d0*(+x1 +x2 +x3 +x4 +x5 +x6 +x7 +x8);
%delta2(1:chexa) = 0.125d0*(+y1 +y2 +y3 +y4 +y5 +y6 +y7 +y8);
%delta3(1:chexa) = 0.125d0*(+z1 +z2 +z3 +z4 +z5 +z6 +z7 +z8);


jac(1:chexa) = alfa31.*(alfa12.*alfa23 -alfa13.*alfa22)...
    - alfa32.*(alfa11.*alfa23 -alfa13.*alfa21)...
    + alfa33.*(alfa11.*alfa22 -alfa12.*alfa21);
%-----



%jac(1:cshell) = alfa1.*beta2 - alfa2.*beta1;
figure
plot([1:chexa],jac,'b-');hold on



%---------------------------------------------------------------



disp(['END - Reading nodes coordinates in ',char(num2str(cputime-tini)),' sec.'])

fid = fopen(nomefilew,'w');

blank_space=length(num2str(grid_id(grid)));
fprintf(fid,'   %i   %i   %i   %i   %i\n',[grid; chexa+cquad; 0; 0; 0]);

tini=cputime;
disp(['BEGIN - Writing inp format'])

if grid>0
    fprintf(fid,'%i  %+13.7e  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
end

if cquad>0
    fprintf(fid,'%i  %i  quad  %i  %i  %i  %i\n',...
        [cquad_id(1:cquad)'; cquad_tag(1:cquad)';...
        cquad_1(1:cquad)'; cquad_2(1:cquad)'; cquad_3(1:cquad)'; cquad_4(1:cquad)']);
end

if chexa>0
    fprintf(fid,'%i  %i   hex  %i  %i  %i  %i  %i  %i  %i  %i\n',...
        [chexa_id(1:chexa)'; chexa_tag(1:chexa)';...
        chexa_1(1:chexa)'; chexa_2(1:chexa)'; chexa_3(1:chexa)'; chexa_4(1:chexa)';...
        chexa_5(1:chexa)'; chexa_6(1:chexa)'; chexa_7(1:chexa)'; chexa_8(1:chexa)']);
end

fclose(fid);




