%% *ANN METADATA ann*
% define series of ANN to train - alternatively define each ANN apart
% _number of ann_
ann.trn.nr = 2; 
% _corner periods_
tmp.TnC = [0.5,0.5];
% _direction (gh=geometrical mean horizontal;ud=vertical)
tmp.cp  = {'gh','ud'};
% _site class (ALL,AB,CD)_
tmp.scl = {'ALL','ALL'};
% _train strategy (classic)_
tmp.str = {'classic','classic'};
% _number of neurons_
tmp.nnr = [10,10];

for i_=1:ann.trn.nr
    ann.trn.mtd(i_).TnC = tmp.TnC(i_);
    ann.trn.mtd(i_).cp  = tmp.cp{i_};
    ann.trn.mtd(i_).scl = tmp.scl{i_};
    ann.trn.mtd(i_).train_strategy = tmp.str{i_};
    ann.trn.mtd(i_).nhn = tmp.nnr(i_);
end
clear tmp;

% % ALL - GH - TnC = 0.50 s
% ann.trn.mtd(1).TnC = 0.50;
% ann.trn.mtd(1).cp  = 'gh';
% ann.trn.mtd(1).scl = 'ALL';
% ann.trn.mtd(1).train_strategy = 'classic';
% ann.trn.mtd(1).nhn = 30;
% 
% % ALL - GH - TnC = 0.75 s
% ann.trn.mtd(2).TnC = 0.75;
% ann.trn.mtd(2).cp  = 'gh';
% ann.trn.mtd(2).scl = 'ALL';
% ann.trn.mtd(2).train_strategy = 'classic';
% ann.trn.mtd(2).nhn = 30;
% 
% % ALL - GH - TnC = 1.00 s
% ann.trn.mtd(3).TnC = 1.00;
% ann.trn.mtd(3).cp  = 'gh';
% ann.trn.mtd(3).scl = 'ALL';
% ann.trn.mtd(3).train_strategy = 'classic';
% ann.trn.mtd(3).nhn = 30;
% 
% % ALL - UD - TnC = 0.50 s
% ann.trn.mtd(4).TnC = 0.50;
% ann.trn.mtd(4).cp  = 'ud';
% ann.trn.mtd(4).scl = 'ALL';
% ann.trn.mtd(4).train_strategy = 'classic';
% ann.trn.mtd(4).nhn = 30;
% 
% % ALL - UD - TnC = 0.75 s
% ann.trn.mtd(5).TnC = 0.75;
% ann.trn.mtd(5).cp  = 'ud';
% ann.trn.mtd(5).scl = 'ALL';
% ann.trn.mtd(5).train_strategy = 'classic';
% ann.trn.mtd(5).nhn = 30;
% 
% % ALL - UD - TnC = 1.00 s
% ann.trn.mtd(6).TnC = 1.00;
% ann.trn.mtd(6).cp  = 'ud';
% ann.trn.mtd(6).scl = 'ALL';
% ann.trn.mtd(6).train_strategy = 'classic';
% ann.trn.mtd(6).nhn = 30;
% 
% % AB - GH - TnC = 0.50 s
% ann.trn.mtd(7).TnC = 0.50;
% ann.trn.mtd(7).cp  = 'gh';
% ann.trn.mtd(7).scl = 'AB';
% ann.trn.mtd(7).train_strategy = 'classic';
% ann.trn.mtd(7).nhn = 30;
% 
% % AB - GH - TnC = 0.75 s
% ann.trn.mtd(8).TnC = 0.75;
% ann.trn.mtd(8).cp  = 'gh';
% ann.trn.mtd(8).scl = 'AB';
% ann.trn.mtd(8).train_strategy = 'classic';
% ann.trn.mtd(8).nhn = 30;
% 
% % AB - GH - TnC = 1.00 s
% ann.trn.mtd(9).TnC = 1.00;
% ann.trn.mtd(9).cp  = 'gh';
% ann.trn.mtd(9).scl = 'AB';
% ann.trn.mtd(9).train_strategy = 'classic';
% ann.trn.mtd(9).nhn = 30;
% 
% % AB - UD - TnC = 0.50 s
% ann.trn.mtd(10).TnC = 0.50;
% ann.trn.mtd(10).cp  = 'ud';
% ann.trn.mtd(10).scl = 'AB';
% ann.trn.mtd(10).train_strategy = 'classic';
% ann.trn.mtd(10).nhn = 30;
% 
% % AB - UD - TnC = 0.75 s
% ann.trn.mtd(11).TnC = 0.75;
% ann.trn.mtd(11).cp  = 'ud';
% ann.trn.mtd(11).scl = 'AB';
% ann.trn.mtd(11).train_strategy = 'classic';
% ann.trn.mtd(11).nhn = 30;
% 
% % AB - UD - TnC = 1.00 s
% ann.trn.mtd(12).TnC = 1.00;
% ann.trn.mtd(12).cp  = 'ud';
% ann.trn.mtd(12).scl = 'AB';
% ann.trn.mtd(12).train_strategy = 'classic';
% ann.trn.mtd(12).nhn = 30;
% 
% % CD - GH - TnC = 0.50 s
% ann.trn.mtd(13).TnC = 0.50;
% ann.trn.mtd(13).cp  = 'gh';
% ann.trn.mtd(13).scl = 'CD';
% ann.trn.mtd(13).train_strategy = 'classic';
% ann.trn.mtd(13).nhn = 30;
% 
% % CD - GH - TnC = 0.75 s
% ann.trn.mtd(14).TnC = 0.75;
% ann.trn.mtd(14).cp  = 'gh';
% ann.trn.mtd(14).scl = 'CD';
% ann.trn.mtd(14).train_strategy = 'classic';
% ann.trn.mtd(14).nhn = 30;
% 
% % CD - GH - TnC = 1.00 s
% ann.trn.mtd(15).TnC = 1.00;
% ann.trn.mtd(15).cp  = 'gh';
% ann.trn.mtd(15).scl = 'CD';
% ann.trn.mtd(15).train_strategy = 'classic';
% ann.trn.mtd(15).nhn = 30;
% 
% % CD - UD - TnC = 0.50 s
% ann.trn.mtd(16).TnC = 0.50;
% ann.trn.mtd(16).cp  = 'ud';
% ann.trn.mtd(16).scl = 'CD';
% ann.trn.mtd(16).train_strategy = 'classic';
% ann.trn.mtd(16).nhn = 30;
% 
% % CD - UD - TnC = 0.75 s
% ann.trn.mtd(17).TnC = 0.75;
% ann.trn.mtd(17).cp  = 'ud';
% ann.trn.mtd(17).scl = 'CD';
% ann.trn.mtd(17).train_strategy = 'classic';
% ann.trn.mtd(17).nhn = 30;
% 
% % CD - UD - TnC = 1.00 s
% ann.trn.mtd(18).TnC = 1.00;
% ann.trn.mtd(18).cp  = 'ud';
% ann.trn.mtd(18).scl = 'CD';
% ann.trn.mtd(18).train_strategy = 'classic';
% ann.trn.mtd(18).nhn = 30;

% _database_
for i_ = 1:ann.trn.nr
    ann.trn.mtd(i_).dbn = dbn;
end
