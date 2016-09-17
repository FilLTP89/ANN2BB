ccc;
% digitized slip
fname = ['/media/filippo/Data/Filippo/PHD_heavyweight/EMILIA_2905/',...
    'exsim_old/exsim_emilia/MRN_new/Emilia_Hisada_Herrero.srcmod'];
new_fname = ['/media/filippo/Data/Filippo/PHD_heavyweight/EMILIA_2905/',...
    'exsim_old/exsim_emilia/MRN_new/Emilia_Hisada_Herrero.txt'];
SLIP = importdata(fname);
NDIP = size(SLIP,1);
NSTR = size(SLIP,2);

fid = fopen(new_fname,'w+');

for ii_=1:NDIP
    for jj_=1:NSTR
        fprintf(fid,'%20.6f ',SLIP(ii_,jj_));
    end
    fprintf(fid,'\n');
end
fclose all;
% interpolated grid
DSTR = 22;
DDIP = 12;
DSFSTR = DSTR/NSTR;
DSFDIP = DDIP/NDIP;

Xo=DSFSTR/2:DSFSTR:DSTR;
Yo=-DSFDIP/2:-DSFDIP:-DDIP;
[Xo,Yo] = meshgrid(Xo,Yo);


colormap jet;
s=surf(Xo,Yo,SLIP,SLIP);
set(s,'EdgeColor','none');
view(2);

set(gca,'xtick',unique([0:5:DSTR,DSTR]));
set(gca,'ytick',unique([-DSTR:5:0,0]));
xlim([0,DSTR]);
ylim([-DDIP,0]);

% interpolated grid
DSFSTR = 0.5;
DSFDIP = 0.5;
NSTR = DSTR/DSFSTR;
NDIP = DDIP/DSFDIP;

Xo_new=DSFSTR/2:DSFSTR:DSTR;
Yo_new=-DSFDIP/2:-DSFDIP:-DDIP;
[Xo_new,Yo_new] = meshgrid(Xo_new,Yo_new);

SLIP_new = interp2(Xo,Yo,SLIP,Xo_new,Yo_new);
max_SLIP_new = max(max(SLIP_new));

figure
colormap jet;
s=surf(Xo_new,Yo_new,SLIP_new,SLIP_new);
set(s,'EdgeColor','none');
view(2);

set(gca,'xtick',unique([0:5:DSTR,DSTR]));
set(gca,'ytick',unique([-DSTR:5:0,0]));
xlim([0,DSTR]);
ylim([-DDIP,0]);

fid = fopen(new_fname,'w+');

for ii_=1:NDIP
    for jj_=1:NSTR
        fprintf(fid,'%12.6f ',SLIP_new(ii_,jj_)./max_SLIP_new);
    end
    fprintf(fid,'\n');
end
fclose all;