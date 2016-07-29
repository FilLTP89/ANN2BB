vtm_shift = zeros(numel(mm_),1);

% dominant direction
for i_ = 1:numel(mm_)
    PGV = 0;
    for j_ = 1:numel(cpp)
        if abs(rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(2))>PGV
            PGV = abs(rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(2));
            PGHI = j_;
        end
    end
    vtm_shift(mm_(i_)) = hbs.syn{mm_(i_)}.pgv.(cpp{PGHI})(1) -...
            rec.org.syn{mm_(i_)}.pgv.(cpp{PGHI})(1);
end