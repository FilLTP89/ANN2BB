%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_run_maps_: function to perform classical LF/HF hybridization between
% SPEED results (PBS) and empirical/stochastic synthetics (SPS/EXS) and spectral
% match them upon ANN predictions. Fast version of the code: suitable for a few 
% demonstrative cases and to plot spectra/time-histories.
%% *N.B.*
% Working-Flow:
% 1). SPLIT IN CHUNKS: divide monitors on several jobs 
% 2). RUN SYN2ANN: syn2ann_run on each chunk of monitor
% 3). WRITE MAPS: write shake maps per each chunk of monitor
% Need for:

%% *1). SPLIT IN CHUNKS
case_chunks.nb = floor(numel(monn.id)/job_nb);

case_chunks.id = cell(job_nb,1);
for NJB = 1:job_nb-1
    case_chunks.id{NJB,1} = monn.id(case_chunks.nb*(NJB-1)+1:case_chunks.nb*NJB);
end     
case_chunks.id{end,1} = monn.id(case_chunks.nb*(NJB)+1:end);

for NJB = 1:job_nb
    selected_case = case_chunks.id{NJB};
    syn2ann_setup_maps_rewind;
end

keyboard

for NJB = 1:job_nb

    job{NJB} = batch('syn2ann_write_maps', 1, {NJB});
end

