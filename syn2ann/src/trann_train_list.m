%% *ANN METADATA ann*
ann.trn.nr = 18;
% gh <---> ALL
ann.trn.mtd(1).TnC = 0.50;
ann.trn.mtd(1).cp = 'gh';
ann.trn.mtd(1).scl = 'ALL';
ann.trn.mtd(1).train_strategy = 'classic';
ann.trn.mtd(1).nhn = 30;

% ud <---> ALL
ann.trn.mtd(2).TnC = 0.75;
ann.trn.mtd(2).cp = 'gh';
ann.trn.mtd(2).scl = 'ALL';
ann.trn.mtd(2).train_strategy = 'classic';
ann.trn.mtd(2).nhn = 30;

% gh <---> ALL
ann.trn.mtd(3).TnC = 1.00;
ann.trn.mtd(3).cp = 'gh';
ann.trn.mtd(3).scl = 'ALL';
ann.trn.mtd(3).train_strategy = 'classic';
ann.trn.mtd(3).nhn = 30;

% ud <---> ALL
ann.trn.mtd(4).TnC = 0.50;
ann.trn.mtd(4).cp = 'ud';
ann.trn.mtd(4).scl = 'ALL';
ann.trn.mtd(4).train_strategy = 'classic';
ann.trn.mtd(4).nhn = 30;

% ud <---> ALL
ann.trn.mtd(5).TnC = 0.75;
ann.trn.mtd(5).cp = 'ud';
ann.trn.mtd(5).scl = 'ALL';
ann.trn.mtd(5).train_strategy = 'classic';
ann.trn.mtd(5).nhn = 30;

% gh <---> ALL
ann.trn.mtd(6).TnC = 1.00;
ann.trn.mtd(6).cp = 'ud';
ann.trn.mtd(6).scl = 'ALL';
ann.trn.mtd(6).train_strategy = 'classic';
ann.trn.mtd(6).nhn = 30;

% gh <---> AB
ann.trn.mtd(1).TnC = 0.50;
ann.trn.mtd(1).cp = 'gh';
ann.trn.mtd(1).scl = 'AB';
ann.trn.mtd(1).train_strategy = 'classic';
ann.trn.mtd(1).nhn = 30;

% ud <---> AB
ann.trn.mtd(2).TnC = 0.75;
ann.trn.mtd(2).cp = 'gh';
ann.trn.mtd(2).scl = 'AB';
ann.trn.mtd(2).train_strategy = 'classic';
ann.trn.mtd(2).nhn = 30;

% gh <---> AB
ann.trn.mtd(3).TnC = 1.00;
ann.trn.mtd(3).cp = 'gh';
ann.trn.mtd(3).scl = 'AB';
ann.trn.mtd(3).train_strategy = 'classic';
ann.trn.mtd(3).nhn = 30;

% ud <---> AB
ann.trn.mtd(4).TnC = 0.50;
ann.trn.mtd(4).cp = 'ud';
ann.trn.mtd(4).scl = 'AB';
ann.trn.mtd(4).train_strategy = 'classic';
ann.trn.mtd(4).nhn = 30;

% ud <---> AB
ann.trn.mtd(5).TnC = 0.75;
ann.trn.mtd(5).cp = 'ud';
ann.trn.mtd(5).scl = 'AB';
ann.trn.mtd(5).train_strategy = 'classic';
ann.trn.mtd(5).nhn = 30;

% gh <---> AB
ann.trn.mtd(6).TnC = 1.00;
ann.trn.mtd(6).cp = 'ud';
ann.trn.mtd(6).scl = 'AB';
ann.trn.mtd(6).train_strategy = 'classic';
ann.trn.mtd(6).nhn = 30;

% gh <---> CD
ann.trn.mtd(1).TnC = 0.50;
ann.trn.mtd(1).cp = 'gh';
ann.trn.mtd(1).scl = 'CD';
ann.trn.mtd(1).train_strategy = 'classic';
ann.trn.mtd(1).nhn = 30;

% ud <---> CD
ann.trn.mtd(2).TnC = 0.75;
ann.trn.mtd(2).cp = 'gh';
ann.trn.mtd(2).scl = 'CD';
ann.trn.mtd(2).train_strategy = 'classic';
ann.trn.mtd(2).nhn = 30;

% gh <---> CD
ann.trn.mtd(3).TnC = 1.00;
ann.trn.mtd(3).cp = 'gh';
ann.trn.mtd(3).scl = 'CD';
ann.trn.mtd(3).train_strategy = 'classic';
ann.trn.mtd(3).nhn = 30;

% ud <---> UD
ann.trn.mtd(4).TnC = 0.50;
ann.trn.mtd(4).cp = 'ud';
ann.trn.mtd(4).scl = 'CD';
ann.trn.mtd(4).train_strategy = 'classic';
ann.trn.mtd(4).nhn = 30;

% ud <---> UD
ann.trn.mtd(5).TnC = 0.75;
ann.trn.mtd(5).cp = 'ud';
ann.trn.mtd(5).scl = 'CD';
ann.trn.mtd(5).train_strategy = 'classic';
ann.trn.mtd(5).nhn = 30;

% gh <---> UD
ann.trn.mtd(6).TnC = 1.00;
ann.trn.mtd(6).cp = 'ud';
ann.trn.mtd(6).scl = 'CD';
ann.trn.mtd(6).train_strategy = 'classic';
ann.trn.mtd(6).nhn = 30;

% _database_
for i_ = 1:ann.trn.nr
    ann.trn.mtd(i_).dbn = dbn;
end
