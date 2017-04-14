
%% *WORKDIR*
%wd = fullfile(filesep,'Users','maria','Documents','PHD','Work','Generation_BB_Groundmotion_with_ANN');
    %% **EMILIA_120529**
    %wd_data =  fullfile(wd,'database_earthqks','EMILIA_120529');
    %wd_results = fullfile(wd,'database_earthqks','EMILIA_120529','PPT_SPEED_results');
    %% **AQUILA_090406**
    % wd_data =  fullfile(wd,'database_earthqks','AQUILA_090406');
    % wd_results = fullfile(wd,'database_earthqks','AQUILA_090406','PPT_SPEED_results');
    %% **SALONICCO_20061978**
    % wd_data =  fullfile(wd,'database_earthqks','SALONICCO_20061978');
    % wd_results = fullfile(wd,'database_earthqks','SALONICCO_20061978','PPT_SPEED_results');
    %% **ISTANBUL**
    % wd_data =  fullfile(wd,'database_earthqks','ISTANBUL');
    % wd_results = fullfile(wd,'database_earthqks','ISTANBUL','PPT_SPEED_results');
     %% **PECHINO**
    % wd_data =  fullfile(wd,'database_earthqks','PECHINO');
    % wd_results = fullfile(wd,'database_earthqks','PECHINO','PPT_SPEED_results');


%% *DEFINE INPUT PARAMETERS*
    %% **EMILIA_120529**
    % MRN0
    % Mw = 6;
    % lon = 11.062;
    % lat = 44.878;
    % R_epi = 4.1; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    % MIR3
    % Mw = 6;
    % lon = 11.105;
    % lat = 44.938;
    % R_epi = 9.82; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    % MIR5
    % Mw = 6;
    % lon = 11.107;
    % lat = 44.981;
    % R_epi = 14.52; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    % SAN0
    % Mw = 6;
    % lon = 11.143;
    % lat = 44.838;
    % R_epi = 4.7; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    % FIN0
    % Mw = 6;
    % lon = 11.287;
    % lat = 44.83;
    % R_epi = 16; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    % MOG0
    % Mw = 6;
    % lon = 10.912;
    % lat = 44.932;
    % R_epi = 16.4; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    % RAV0
    % Mw = 6;
    % lon = 11.143;
    % lat = 44.716;
    % R_epi = 15.7; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    % MIR8
    % Mw = 6;
    % lon = 11.09;
    % lat = 44.917;
    % R_epi = 7.33; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    % MIR1
    % Mw = 6;
    % lon = 11.071;
    % lat = 44.844;
    % R_epi = 1.41; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    % MIR2
    % Mw = 6;
    % lon = 11.073;
    % lat = 44.887;
    % R_epi = 4.12; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    % SMS0
    % Mw = 6;
    % lon = 11.235;
    % lat = 44.934;
    % R_epi = 15; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %% **AQUILA_090406**
    % %AQK-mon.15320
    % Mw = 6.2;
    % lon = ;
    % lat = ;
    % R_epi = 1.787; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %AQU-mon.15321
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 1.949; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %Paganica-mon.14811
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 7.584; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %AQV-mon.15324
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 4.216; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %Onna-mon.15062
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 7.304; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %Castelnuovo-mon.15392
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 21.446; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.15121
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 13.797; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %GSA-mon.15397
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 13.879; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %CLN-mon.15398
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 31.595; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %MTR-mon.15401
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 22.343; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.13002
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 5.201; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.8365
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 12.380; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.11305
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 19.017; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.11650
    % Mw = 6.2;
    % lon = 
    % lat = 
    % R_epi = 26.626; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================	
    %% **ISTANBUL**
    % % Airport Istanbul-mon.4560
    % lon = 28.81058;
    % lat = 40.98558;
    % R_epi = 33.157; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Ayasofya-mon.5751
    % lon = 28.98158;
    % lat = 41,00809;
    % R_epi = 47.312; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Prince Island Burgazadasi-mon.6568
    % lon = 29.06708;
    % lat = 40.87758;
    % R_epi = 51.268; % epicentral distance [Km] 
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Bosforo Bridge-mon.6597
    % lon = 29.06708;
    % lat = 41.08908;
    % R_epi = 57.593; % epicentral distance [Km] 
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Stadium-mon.6207
    % lon = 29.03558;
    % lat = 40.98558;
    % R_epi = 50.772; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Istanbul Technical University-mon.6118
    % lon = 29.02658;
    % lat = 41.10709;
    % R_epi = 55.586; % epicentral distance [Km] 
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Sport Arena-mon.6944
    % lon = 29.10308;
    % lat = 40.99458;
    % R_epi = 56.500; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Fatih University-mon.4400
    % lon = 28.78808;
    % lat = 41.03058;
    % R_epi = 34.151; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % Prince Island Adalar-mon.3459
    % lon = 28.62158;
    % lat = 41.08908;
    % R_epi = 29.854; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % mon.7183
    % lon = 29.13008;
    % lat = 40.87308;
    % R_epi = 56.547; % epicentral distance [Km] 
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % mon.8885
    % lon = 29.40910;
    % lat = 40.76510;
    % R_epi = 80.626; % epicentral distance [Km] 
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % mon.6874
    % lon = 29.09860;
    % lat = 40.87310;
    % R_epi = 53.894; % epicentral distance [Km] 
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % mon.7661
    % lon = 29.19760;
    % lat = 40.88660;
    % R_epi = 62.307; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % mon.8528
    % lon = 29.35960;
    % lat = 40.83260;
    % R_epi = 75.885; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % mon.4039
    % lon = 28.72060;
    % lat = 40.97660;
    % R_epi = 26.079; % epicentral distance [Km] 
    % scc = 1; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %% **SALONICCO_20061978**
    % % THE - mon.16456
    % Mw = 6.50;
    % lon = 22.933;
    % lat = 40.633;
    % R_epi =  28.67; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.15462
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 15.91; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.6039
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 17.16; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % %mon.14302
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 21.88; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %mon.15271
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 4.64; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %mon.11832
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 26.65; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %mon.13469
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 10.51; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %mon.15829
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 9.58; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %mon.15314
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 27.80; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %mon.12252
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 10.49; % epicentral distance [Km]
    % scc = 0; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % % AMP - mon.16458
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 29.67; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %LAB - mon.16460
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 27.30; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
    %======================================================================================================
    % %POL-mon.16464
    % Mw = 6.50;
    % lon = ;
    % lat = ;
    % R_epi = 28.95; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %% **PECHINO**
    % % mon.19622
    % lon = 116.52;
    % lat = 39.62;
    % R_epi =  15.058; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.20536
    % lon = 116.40;
    % lat = 39.60;
    % R_epi =  6.180; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.20574
    % lon = 116.43;
    % lat = 39.60;
    % R_epi =  7.186; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.21465
    % lon = 116.31;
    % lat = 39.57;
    % R_epi =  6.322; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.21422
    % lon = 116.40;
    % lat = 39.56;
    % R_epi =  1.854; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
     % % mon.21028
    % lon = 116.30;
    % lat = 39.54;
    % R_epi =  7.439; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.19997
    % lon = 116.31;
    % lat = 39.50;
    % R_epi =  11.071; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.21454
    % lon = 116.36;
    % lat = 39.56;
    % R_epi =  2.425; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.20631
    % lon = 116.33;
    % lat = 39.49;
    % R_epi =  10.473; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.19842
    % lon = 116.30;
    % lat = 39.49;
    % R_epi =  11.978; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.21689
    % lon = 116.36;
    % lat = 39.53;
    % R_epi =  4.017; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.21578
    % lon = 116.38;
    % lat = 39.55;
    % R_epi =  0.323; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.20689
    % lon = 116.49;
    % lat = 39.59;
    % R_epi =  10.984; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.20962
    % lon = 116.44;
    % lat = 39.59;
    % R_epi =  6.982; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.20855
    % lon = 116.49;
    % lat = 39.59;
    % R_epi = 10.647; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.19113
    % lon = 116.55;
    % lat = 39.62;
    % R_epi = 17.164; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
     % % mon.19614
    % lon = 116.52;
    % lat = 39.60;
    % R_epi = 13.156; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
     % % mon.20369
    % lon = 116.28;
    % lat = 39.53;
    % R_epi = 9.581; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
     % % mon.19272
    % lon = 116.47;
    % lat = 39.63;
    % R_epi = 12.648; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.21521
    % lon = 116.45;
    % lat = 39.54;
    % R_epi = 5.832; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.16265
    % lon = 116.21;
    % lat = 39.49;
    % R_epi = 17.506; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.16446
    % lon = 116.21;
    % lat = 39.48;
    % R_epi = 18.127; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.16023
    % lon = 116.26;
    % lat = 39.46;
    % R_epi = 18.391; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
    % % mon.16455
    % lon = 116.21;
    % lat = 39.51;
    % R_epi = 16.062; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
     % % mon.17985
    % lon = 116.24;
    % lat = 39.46;
    % R_epi = 18.593; % epicentral distance [Km]
    % scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
    %======================================================================================================
%% *NUMERICAL SIMULATION*
% simulationID
    %% **EMILIA_120529**
    % % MRN0
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16928';
    %======================================================================================================
    % % MIR3
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '13249';
    %======================================================================================================
    % MIR5
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '10265';
    %======================================================================================================
    % % SAN0
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '17670';
    %======================================================================================================
    % %RAV0
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '12306';
    %======================================================================================================
    % FIN0
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '17869';
    %======================================================================================================
    % MOG0
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '14120';
    %======================================================================================================
    % MIR8
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15045';
    %======================================================================================================
    % MIR1
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '18446';
    %======================================================================================================
    % MIR2
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '18437';
    %======================================================================================================
    %SMS0
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '13497';
    %% **AQUILA_090406**
    % %AQK-mon.15320
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15320';
    %======================================================================================================
    % %AQU-mon.15321
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15321';
    %======================================================================================================
    %AQV-mon.15324
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15324';
    %======================================================================================================
    %Paganica-mon.14811
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '14811';
    %======================================================================================================
    % %Onna-mon.15062
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15062';
    %======================================================================================================
    % %Castelnuovo-mon.15062
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15392';
    %======================================================================================================
    % % mon.15121
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15121';
    %======================================================================================================
    % % GSA-mon.15397
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15397';
    %======================================================================================================
    % % CLN-mon.15398
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15398';
    %======================================================================================================
    % % MTR-mon.15401
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15401';
    %======================================================================================================
    % % mon.13002
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '13002';
    %======================================================================================================
    % % mon.8365
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '8365';
    %======================================================================================================
    % % mon.11305
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '11305';
    %======================================================================================================
    % % mon.11650
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '11650';
    %% **ISTANBUL**
    % % Airport Istanbul-mon.4560
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '4560';
    % Mw = 7;
    %======================================================================================================
    % % Ayasofya-mon.5751
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '5751';
    % Mw = 7;
    %======================================================================================================
    % % Prince Island Burgazadasi-mon.6568
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '6568';
    % Mw = 7;
    %======================================================================================================
    % % Bosforo Bridge-mon.6597
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '6597';
    % Mw = 7;
    %======================================================================================================
    % % Stadium-mon.6207
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '6207';
    % Mw = 7;
    %======================================================================================================
    % % Istanbul Technical University-mon.6118
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '6118';
    % Mw = 7;
    %======================================================================================================
    % % Sport Arena-mon.6944
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '6944';
    % Mw = 7;
    %======================================================================================================
    % % Fatih University-mon.4400
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '4400';
    % Mw = 7;);
    %======================================================================================================
    % % Prince Island Adalar-mon.3459
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '3459';
    % Mw = 7;
    %======================================================================================================
    % % mon.7183
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '7183';
    % Mw = 7;
    %======================================================================================================
    % % mon.8885
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '8885';
    % Mw = 7;
    %======================================================================================================
    % % mon.6874
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '6874';
    % Mw = 7;
    %======================================================================================================
    % % mon.7661
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '7661';
    % Mw = 7;
    %======================================================================================================
    % % mon.8528
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '8528';
    % Mw = 7;
    %======================================================================================================
    % % mon.4039
    % num_sim.simID = 'E00508'; 
    % num_sim.monID = '4039';
    % Mw = 7;
    %% **SALONICCO_20061978**
    % % THE - mon.16456
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16456';
    %======================================================================================================
    % %mon.15462
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15462';
    %======================================================================================================
    % %mon.6039
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '6039'; 
    %======================================================================================================
    % %mon.14302
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '14302'; 
    %======================================================================================================
    % %mon.15271
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15271';
    %======================================================================================================
    % %mon.11832
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '11832';
    %======================================================================================================
    % %mon.13469
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '13469';
    %======================================================================================================
    % %mon.15829
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15829';
    %======================================================================================================
    % %mon.15314
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '15314';
    %======================================================================================================
    % %mon.12252
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '12252';
    %======================================================================================================
    % % AMP - mon.16458
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16458';
    %======================================================================================================
    % %LAB - mon.16460
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16460';
    %======================================================================================================
    % %POL-mon.16464
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16464';
    %======================================================================================================
    %% **PECHINO**
    % % mon.19622
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '19622';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20536
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20536';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20574
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20574';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.21465
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21465';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.21422
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21422';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.21028
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21028';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.19997
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '19997';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.21454
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21454';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20631
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20631';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.19842
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '19842';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.21689
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21689';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.21578
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21578';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20689
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20689';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20962
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20962';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20855
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20855';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.19113
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '19113';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.19614
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '19614';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.20369
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '20369';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.19272
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '19272';
    % Mw = 6.5;
    %======================================================================================================
     % % mon.21521
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '21521';
    % Mw = 6.5;
    %======================================================================================================
     % % mon.16265
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16265';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.16446
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16446';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.16023
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16023';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.16455
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '16455';
    % Mw = 6.5;
    %======================================================================================================
    % % mon.17985
    % num_sim.simID = 'E00001'; 
    % num_sim.monID = '17985';
    % Mw = 6.5;
    %======================================================================================================