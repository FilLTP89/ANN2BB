%% *SPECTRAL SCALING*
fprintf('---------------------\n5. SPECTRAL MATCHING\n---------------------\n');
%% *MATCHING*
fprintf('--> Matching\n');
spm = syn2ann_sm(hbs,trs);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
for j_ = 1:hbs.mon.nc
    spm.(hbs.mon.cp{j_}) = syn2ann_thp(spm.(hbs.mon.cp{j_}));
end
%% *PGA-PGV-PGD & ARIAS INTENSITY & SPECTRA*
fprintf('--> Peak Values and Arias\n');
for j_ = 1:hbs.mon.nc
    spm.(hbs.mon.cp{j_}) = syn2ann_thp(spm.(hbs.mon.cp{j_}));
%     spm.(hbs.mon.cp{j_}) = syn2ann_spp(spm.(hbs.mon.cp{j_}));
    spm.(hbs.mon.cp{j_}).mon.vfr = cell(spm.(hbs.mon.cp{j_}).mon.na,1);
    %% *FOURIER SPECTRUM*
    for i_ = 1:spm.(hbs.mon.cp{j_}).mon.na
        % Nyquist frequency
        fr_max      = 1/2/spm.(hbs.mon.cp{j_}).mon.dtm(i_);
        % fft points
        spm.(hbs.mon.cp{j_}).mon.nfr(i_) = 2^nextpow2(spm.(hbs.mon.cp{j_}).mon.ntm(i_));
        % frequency period step
        spm.(hbs.mon.cp{j_}).mon.dfr(i_) = 1/spm.(hbs.mon.cp{j_}).mon.dtm(i_)/(spm.(hbs.mon.cp{j_}).mon.nfr(i_)-1);
        % Nyquist frequency index
        spm.(hbs.mon.cp{j_}).mon.nNy(i_) = floor(fr_max/spm.(hbs.mon.cp{j_}).mon.dfr(i_))+1;
        % frequency vector
        
        spm.(hbs.mon.cp{j_}).mon.vfr{i_} = spm.(hbs.mon.cp{j_}).mon.dfr(i_)*(0:spm.(hbs.mon.cp{j_}).mon.nfr(i_)-1)';

        spm.(hbs.mon.cp{j_}).syn{i_}.fsa.(spm.(hbs.mon.cp{j_}).mon.cp{1}) = ...
            spm.(hbs.mon.cp{j_}).mon.dtm(i_)*fft(spm.(hbs.mon.cp{j_}).syn{i_}.tha.(spm.(hbs.mon.cp{j_}).mon.cp{1}),...
            spm.(hbs.mon.cp{j_}).mon.nfr(i_));
    end

end