ccc;
% wd.local = '/media/filippo/Data/Filippo/PHD_heavyweight/SPEED2D_NL/INPUTS/pwn';
% sp = '/media/filippo/Data/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images';
wd.local = '/media/user/DATI/Filippo/PHD_heavyweight/SPEED2D_NL/INPUTS/pwn';
sp = '/media/user/DATI/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images';

cd(wd.local);
fn = struct2cell(dir('monitor*'));
fn = fn(1,:);

id = unique(cellfun(@(x) x(strfind(x,'r')+1:strfind(x,'.')-1),fn,'UniformOutput',0)');
nm = numel(id);

strain_label = {'\epsilon_{xx} [1]','\epsilon_{yy} [1]','\epsilon_{zz} [1]',...
    '\gamma_{xy} [1]','\gamma_{xz} [1]','\gamma_{yz} [1]'};
pstrain_label = {'\epsilon_{xx}^{pl} [1]','\epsilon_{yy}^{pl} [1]','\epsilon_{zz}^{pl} [1]',...
    '\gamma_{xy}^{pl} [1]','\gamma_{xz}^{pl} [1]','\gamma_{yz}^{pl} [1]'};
stress_label = {'\sigma_{xx} [kPa]','\sigma_{yy} [kPa]','\sigma_{zz} [kPa]',...
    '\tau_{xy} [kPa]','\tau_{xz} [kPa]','\tau_{yz} [kPa]'};

vr = {'d';'v';'a';'e';'s'};
nv = numel(vr);
time = cell(nm,nv);

for i_=1:nm
    for j_=1:nv
        try
            pt = strcat('monitor',id(i_),'.',vr{j_});
            p(i_).(vr{j_}) = importdata(pt{1});
            time{i_,j_}  = p(i_).(vr{j_})(:,1);
            p(i_).(vr{j_}) = p(i_).(vr{j_})(:,2:end);
        catch
            keyboard
        end
    end
end

idx = cellfun(@(x) min(numel(x)),time);
idx = min(min(idx));

p = arrayfun(@(x) structfun(@(y) y(1:idx,:),x,'UniformOutput',0),p);
% arrayfun(@(x) structfun(@(y) fprintf('%u\n',size(y,1)),x,'UniformOutput',0),p);
time = time(1);
sc = 1:5;
ns = numel(sc);

xpl = cell(ns,1);
ypl = cell(ns,1);

for i_ = 1:ns
    xpl(i_) = {p(sc(i_)).e(:,3)};
    ypl(i_) = {p(sc(i_)).s(:,4)};
end

fpplot('xpl',xpl,'ypl',ypl,'xlb',{'\gamma_{xy} [1]'},'ylb',{'\tau_{xy} [kPa]'},...
    'vfg','on');