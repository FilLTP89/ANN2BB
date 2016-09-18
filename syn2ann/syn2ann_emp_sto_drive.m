fprintf('---------------------\n3. EMPIRICAL\STOCHASTIC\n---------------------\n');

switch lower(hybrid_type)
    case 'sp96'
        % _SABETTA & PUGLIESE 1996_
        syn2ann_sp96_drive;
    case 'exsim'
        % _EXSIM_
        syn2ann_exsim_drive;
    case 'both'
        % _SABETTA & PUGLIESE 1996_
        syn2ann_sp96_drive;
        % _EXSIM_
        syn2ann_exsim_drive;
end