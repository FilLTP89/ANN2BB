ccc;
xpl = cell(1,1);
ypl = cell(1,1);
for i_=1
    temp = textread(sprintf('MRN_exsim_acc_s001.%u.1',i_),'',...
        'headerlines',12);
    xpl{i_,1} = temp(:,1);
    ypl{i_,1} = temp(:,2);
end

fpplot('xpl',xpl,'ypl',ypl,'xlm',{[0,650]},'vfg','on');