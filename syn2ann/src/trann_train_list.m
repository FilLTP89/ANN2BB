%% *ANN METADATA ann*
ann.trn.nr = 4;
% % gh <---> AB
% ann.trn.mtd(1).TnC = 0.75;
% ann.trn.mtd(1).cp = 'gh';
% ann.trn.mtd(1).scl = 'AB';
% % ud <---> AB
% ann.trn.mtd(2).TnC = 0.75;
% ann.trn.mtd(2).cp = 'ud';
% ann.trn.mtd(2).scl = 'AB';
% % gh <---> CD
% ann.trn.mtd(3).TnC = 0.75;
% ann.trn.mtd(3).cp = 'gh';
% ann.trn.mtd(3).scl = 'CD';
% % ud <---> CD
% ann.trn.mtd(4).TnC = 0.75;
% ann.trn.mtd(4).cp = 'ud';
% ann.trn.mtd(4).scl = 'CD';
% gh <---> ALL
ann.trn.mtd(1).TnC = 0.5;
ann.trn.mtd(1).cp = 'gh';
ann.trn.mtd(1).scl = 'ALL';
ann.trn.mtd(1).train_strategy = 'classic';
% ud <---> ALL
ann.trn.mtd(2).TnC = 0.5;
ann.trn.mtd(2).cp = 'ud';
ann.trn.mtd(2).scl = 'ALL';
ann.trn.mtd(2).train_strategy = 'classic';

% gh <---> ALL
ann.trn.mtd(3).TnC = 0.5;
ann.trn.mtd(3).cp = 'gh';
ann.trn.mtd(3).scl = 'ALL';
ann.trn.mtd(3).train_strategy = 'bootstrap';

% ud <---> ALL
ann.trn.mtd(4).TnC = 0.5;
ann.trn.mtd(4).cp = 'ud';
ann.trn.mtd(4).scl = 'ALL';
ann.trn.mtd(4).train_strategy = 'bootstrap';

% _database_
for i_ = 1:ann.trn.nr
    ann.trn.mtd(i_).dbn = dbn;
end
