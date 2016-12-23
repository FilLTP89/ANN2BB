function [varargout] = train_ann_valid(varargin)
    %% *SET-UP*
    ann = varargin{1};
    
    %% *ANN VALIDATION*
    fprintf('VALIDATING...\n');
    ann.out_trn.all = ann.net(ann.inp.trn);
    ann.out_vld.all = ann.net(ann.inp.vld);
    
%     %% *DEFINE INPUTS/TARGETS SUBSETS FROM TRAIN PROCEDURE*
%     % _TRANING/VALIDATION/TEST SUBSET VALUES_
%     ann.out_trn.trn = ann.out_trn.all(:,ann.trs.trainInd);
%     ann.out_trn.vld = ann.out_trn.all(:,ann.trs.valInd);
%     ann.out_trn.tst = ann.out_trn.all(:,ann.trs.testInd);
%     
%     % _TRANING/VALIDATION/TEST SUBSET VALUES_
%     ann.out_tar.trn = ann.tar.trn(:,ann.trs.trainInd);
%     ann.out_tar.vld = ann.tar.trn(:,ann.trs.valInd);
%     ann.out_tar.tst = ann.tar.trn(:,ann.trs.testInd);

    %% *COMPUTE PERFORMANCE*
    fprintf('COMPUTING PERFORMANCE...\n');
    ann.prf.trn = perform(ann.net,ann.tar.trn,ann.out_trn.all);
    ann.prf.vld = perform(ann.net,ann.tar.vld,ann.out_vld.all);
    
    %% *OUTPUT*
    varargout{1} = ann;
    return
end