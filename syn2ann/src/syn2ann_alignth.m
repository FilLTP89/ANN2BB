% dominant direction

% FILTERED RECORD
PGV = 0;
for j_ = 1:numel(cpp)
    if abs(rec.fil.syn{mm_(i_)}.pgv.(cpp{j_})(2))>PGV
        PGV = abs(rec.fil.syn{mm_(i_)}.pgv.(cpp{j_})(2));
        PGHI = j_;
    end
end
vtm_shift(1) = hbs.syn{mm_(i_)}.pgv.(cpp{PGHI})(1) -...
    rec.fil.syn{mm_(i_)}.pgv.(cpp{PGHI})(1);

% ORIGINAL RECORDS
PGV = 0;
for j_ = 1:numel(cpp)
    if abs(rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(2))>PGV
        PGV = abs(rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(2));
        PGHI = j_;
    end
end
vtm_shift(2) = hbs_org.syn{mm_(i_)}.pgv.(cpp{PGHI})(1) -...
    rec.org.syn{mm_(i_)}.pgv.(cpp{PGHI})(1);
