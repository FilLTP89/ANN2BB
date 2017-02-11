Np = 3;
for k_ = 1:tst.mtd.nr
    
%     [inp.vTn,~,~,~] = trann_define_inout(tst.mtd.TnC(k_));
    idx = find(rec.org.mon.vTn==tst.mtd.Tno(k_));
    %     [~,idx,~] = intersect(rec.org.mon.vTn,inp.vTn);
    for j_=1:rec.org.mon.nc
        cpp = rec.org.mon.cp{j_};
        for i_ = 1:rec.org.mon.na
            rgr{k_}.psa.(cpp)(i_,1) = abs(rec.org.syn{i_}.psa.(cpp)(idx));
            rgr{k_}.pga.(cpp)(i_,1) = abs(rec.org.syn{i_}.pga.(cpp)(2));
        end
        rgr{k_}.cfc.(cpp) = polyfit(rgr{k_}.psa.(cpp),rgr{k_}.pga.(cpp),Np);
        
        fpplot('xpl',{rgr{k_}.psa.(cpp),sort(rgr{k_}.psa.(cpp))},...
            'ypl',{rgr{k_}.pga.(cpp),polyval(rgr{k_}.cfc.(cpp),sort(rgr{k_}.psa.(cpp)))},...
            'vfg','on','pfg',[0,0,10,10],'mrk',{'o','none'},'lst',{'none','-'});
        keyboard
    end
end