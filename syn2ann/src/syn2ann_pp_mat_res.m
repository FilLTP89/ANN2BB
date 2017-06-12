ccc;
pfx = 'syn2ann_res_thess_ALL_';
fld = '/tmp1/gattif/syn2ann_maps_thess';
cpn = {'ew';'ns';'ud'};
tmp = dir(fullfile(fld,strcat(pfx,'*.mat')));
nfm = numel(tmp);
[fnm{1:nfm,1}] = deal(tmp.name);
[~,id.prc] = ...
    sort(str2double(cellfun(@(x) x((regexpi(x,pfx,'end'))+1:...
    (regexpi(x,'.mat'))-1),fnm,'UniformOutput',false)));
pbs = cell(2985,1);
hbs = cell(2985,1);
spm = cell(2985,1);
mon = cell(nfm,1);
bhr = cell(nfm,1);
keyboard
cnt = 1;
for i_=1:nfm
    temp = load(fullfile(fld,fnm{id.prc(i_)}));
    [mon{i_,1}] = deal(temp.mon);
    [bhr{i_,1}] = deal(temp.bhr);
    for j_=1:temp.mon.na
        for k_=1:numel(cpn)
            pbs{cnt,1}.psa.(cpn{k_}) = temp.pbs.bst.syn{j_}.psa.(cpn{k_});
            hbs{cnt,1}.psa.(cpn{k_}) = temp.hbs.bst.syn{j_}.psa.(cpn{k_});
            spm{cnt,1}.psa.(cpn{k_}) = temp.spm.sps.(cpn{k_}).syn{j_}.psa.(cpn{k_});
        end
        cnt=cnt+1;
    end
end

save('pbs_thessaloniki.mat','pbs');
save('hbs_thessaloniki.mat','hbs');
save('spm_thessaloniki.mat','spm');
save('mon_thessaloniki.mat','mon');
save('bhr_thessaloniki.mat','bhr');
