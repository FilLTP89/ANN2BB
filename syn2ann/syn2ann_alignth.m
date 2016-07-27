vtm_shift = zeros(numel(mm_),1);

% mean among directions
% for i_ = 1:numel(mm_)
%     for j_ = 1:numel(cpp)
%         vtm_shift(mm_(i_)) = vtm_shift(mm_(i_)) + hbs.syn{mm_(i_)}.pgv.(cpp{j_})(1) -...
%             rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(1);
%     end
%     vtm_shift(mm_(i_)) = vtm_shift(mm_(i_))/numel(cpp);
% end

% dominant direction
for i_ = 1:numel(mm_)
    pgv = 0;
    for j_ = 1:numel(cpp)
        pgv_temp = max([abs(hbs.syn{mm_(i_)}.pgv.(cpp{j_})(2)),...
            abs(rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(2))]);
        if pgv_temp>pgv
            pgv = pgv_temp;
            vtm_shift(mm_(i_)) = hbs.syn{mm_(i_)}.pgv.(cpp{j_})(1) -...
                rec.org.syn{mm_(i_)}.pgv.(cpp{j_})(1);
        end
    end
end