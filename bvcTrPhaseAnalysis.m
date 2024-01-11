function [RM] = bvcTrPhaseAnalysis( RM )
% Analysis of trace cell spike phase (of local theta), in baseline and trace fields.
% NOTE: requires raw data in SCAN open on base workspace, and SCAN data needs to have been
% processed by bvcTrGetBestThetaEEGForTet first.

% Define main field.
% prms.mapNormMode  = 'Z';
% prms.BslFieldThr  = 0.5; 
% EEG
prms.minSNR        = 2;
prms.EEGChooseMode = 'DIST';  % 'AZ2';  %   'Tet';  %  ' HS';  % 'AZ'
% Dir maps (directionality in field).
prms.dirBin       = 6;
% AC
prms.runAC        = 0;
prms.ACBin        = 0.001;
prms.ACWin        = 0.020;
prms.useGPU       = 1;
% Pre-allocate results fields. A bit long-winded, as we need all combinations of 1) probe,
% 2) field type, 3) actual analysis score. All have the same format, though, 
fieldNames  = {'AllSpk',     'BslF',     'BarrF',     'PostF',     'BkGr'    ,  'BarrPostF'};
scoreNames  = {'Ph', 'RV', 'NSpk', 'MR', 'pR', 'pHA', 'DirRV', 'DirSpk', 'DirPos', 'DirRte'};
scoreIsCell = [ 0     0      0      0      0     0       0         1        1          1   ];
for itPr = 1:3
    for itFN = 1:length(fieldNames)
        for itSN = 1:length(scoreNames)
            if scoreIsCell(itSN)
                RM.( [scoreNames{itSN} fieldNames{itFN} '_Pr' num2str(itPr)]) = cell( height(RM), 3 );
            else
                RM.( [scoreNames{itSN} fieldNames{itFN} '_Pr' num2str(itPr)]) = nan( height(RM), 3 );
            end
        end
        if prms.runAC
            RM.( ['AC' fieldNames{itFN} '_Pr' num2str(itPr)])  = cell( height(RM), 3 );
            RM.( ['ISI' fieldNames{itFN} '_Pr' num2str(itPr)]) = cell( height(RM), 3 );
        end
    end
end
%
RM.trMinSNR     = nan( height(RM), 1 );
RM.trMaxPerRail = nan( height(RM), 1 );
RM.meanThFreq   = nan( height(RM), 1 );
RM.EEGTet       = nan( height(RM), 1 );
RM.EEGAnatZone  = nan( height(RM), 1 );
RM.EEGIsMatched = true( height(RM), 1 );


RM.Properties.UserData.Phase = prms;

hWait = waitbar(0,'Please wait .. ');
for itCl = 1:height( RM )
    waitbar(itCl/height(RM),hWait);
    
    % Get raw data from workspace.
    scanDSName = ['R' num2str(RM.Rat(itCl)) '_' RM.Date{itCl}];
%     if strcmp( scanDSName, 'R460_131115' );   continue;   end   % Specific problem with missing EEG files, skip unitl resolved.
    data = evalin( 'base', scanDSName );   % Get raw data in scan from base workspace.
    
    % Get the correct EEG (best SNR for either Tet or hemisphere).
    if strcmp( prms.EEGChooseMode, 'Tet' )       % Match tet
        eegInd = find( cat(1, data.trials(1).eeg.tet) == RM.Tet(itCl)    &   cat(1, data.trials(1).eeg.bestForTet)  );
    elseif strcmp( prms.EEGChooseMode, 'HS' )    % Match hemisphere
        eegInd = find( cat(1, data.trials(1).eeg.HS) == RM.HS(itCl)      &   cat(1, data.trials(1).eeg.bestForHS)  );
    elseif strcmp( prms.EEGChooseMode, 'AZ' )    % Match anat zones 1 - 3 (P-M-D)
        eegInd = find( cat(1, data.trials(1).eeg.HS) == RM.HS(itCl)      &   cat(1, data.trials(1).eeg.anatZone) == RM.anatZone(itCl)     &    cat(1, data.trials(1).eeg.bestForZone)  );
    elseif strcmp( prms.EEGChooseMode, 'AZ2' )   % Match anat zone 1 vs 2/3 (P - M/D)
        if RM.anatZone(itCl) == 1
            eegInd = find( cat(1, data.trials(1).eeg.HS) == RM.HS(itCl)      &   cat(1, data.trials(1).eeg.anatZone) == 1     &    cat(1, data.trials(1).eeg.bestForZone)  );
        else
            eegInd = find( cat(1, data.trials(1).eeg.HS) == RM.HS(itCl)      &   cat(1, data.trials(1).eeg.bestForDist)      );
        end
    elseif strcmp( prms.EEGChooseMode, 'DIST' ) % Best of zones 2/3 (MD) in hemisphere
        eegInd = find( cat(1, data.trials(1).eeg.HS) == RM.HS(itCl)      &   cat(1, data.trials(1).eeg.bestForDist)      ); 
    end
    if isempty(eegInd) || data.trials(1).eeg(eegInd).trMinSNR < prms.minSNR
        eegInd                = find( cat(1, data.trials(1).eeg.HS) == RM.HS(itCl)    &    cat(1, data.trials(1).eeg.bestForHS)  );
        RM.EEGIsMatched(itCl) = false;
    end
    % If we have managed to find an EEG, record the quality characterisitcs here.
    if ~isempty( eegInd )
        RM.trMinSNR(itCl)     = data.trials(1).eeg(eegInd).trMinSNR;
        RM.trMaxPerRail(itCl) = data.trials(1).eeg(eegInd).trMaxPerRail;
        RM.meanThFreq(itCl)   = data.trials(1).eeg(eegInd).meanThFreq;
        RM.EEGTet(itCl)       = data.trials(1).eeg(eegInd).tet;
        RM.EEGAnatZone(itCl)  = data.trials(1).eeg(eegInd).anatZone;
    else
        continue   % If we don't even have a hemisphere EEG, skip phase analysis for this cell.
    end
    
    for itPr = 1:3
        if isempty( RM.(['RMProbe' num2str(itPr)]){itCl,1} );   continue;   end
        
        trInds = RM.(['IndScnPr' num2str(itPr)])(itCl,:);
        for itTr = 1:3  % itTr = iterator across pre-, probe, post-
            
            % Need to check if this cell is one for which only pre-probe baseline exists.
            if isempty( RM.(['RMProbe' num2str(itPr)]){itCl,2} )
                if itTr==1
                    % If so, if this *is* the pre-probe trial, create null field masks for all expect BslF, and continue with loop iteration.
                    RM.BFMask{itCl,itPr}    = zeros( size(RM.BslFMask{itCl,itPr}) );
                    RM.PoPrFMask{itCl,itPr} = zeros( size(RM.BslFMask{itCl,itPr}) );
                else
                    % Otherwise, skip this loop iteration.
                    continue;
                end
            end
            
            % Get LFP, get phases.
            EEG                 = data.trials( trInds(itTr) ).eeg(eegInd);
            EEGRailInd          = EEG.eeg==127 | EEG.eeg==-128;
            [~,eegPh,~]         = eeg_instfrequency( EEG, EEG.thetaFreq + [-3 3] );
            eegPh( EEGRailInd ) = nan;
            % Filter for speed (set phase to NaN when speed is slow).
            speed                 = double( data.trials( trInds(itTr) ).speed );
            pos2eegInd            = reshape( (1:length(EEG.eeg))', 5, [] );
            immIndForEEG          = reshape( pos2eegInd( :, speed<5 ), [], 1 );
            eegPh( immIndForEEG ) = nan;
            % Get phase for spike
            spkTAsEEGInd        = ceil( data.trials( trInds(itTr) ).cells(RM.IndClScn(itCl)).st .* EEG.sample_rate );
            spkPh               = eegPh( spkTAsEEGInd );
                        
            % Convert raw XY to a (linear) bin index.
            binX         = ceil( data.trials( trInds(itTr) ).x ./ RM.Properties.UserData.binSize );
            binY         = ceil( data.trials( trInds(itTr) ).y ./ RM.Properties.UserData.binSize );
            binXYLin     = sub2ind( size( RM.(['RMProbe' num2str(itPr)]){itCl,itTr}), binY, binX );
            % For each spike, get its corresponding position (as a linear index into rate map).
            spkTAsPosInd = ceil( data.trials( trInds(itTr) ).cells(RM.IndClScn(itCl)).st .* data.trials( trInds(itTr) ).sample_rate );
            spkPosBin    = binXYLin( spkTAsPosInd );
            % .. and the same thing for head direction (for making dir maps).
            dir           = data.trials( trInds(itTr) ).dir;
            dir( dir<=0 ) = 1;
            dir( dir>360) = 360;
            binDir        = ceil( dir ./ prms.dirBin );
            spkDirBin     = binDir( spkTAsPosInd );
            

            % Get defintions of baseline field, barrier field and post-probe field.
            linIndBslF   = find( RM.BslFMask{itCl,itPr}  == 1 );
            linIndBarrF  = find( RM.BFMask{itCl,itPr}    == 1 );
            linIndPostF  = find( RM.PoPrFMask{itCl,itPr} == 1 );
            % 'linIndBkGr' needs to have a different definition depending on itTr. Otherwise, we are excluding
            %  spikes in union of the barrier and post-barr fields, in ALL trials in the probe set (which is wrong, I think).
%             if itTr==1
%                 linIndBkGr = find( ~mainFMask  );
%             elseif itTr==2
                linIndBkGr = find( ~(  RM.BslFMask{itCl,itPr}==1 | RM.BFMask{itCl,itPr}==1  ) );
%             elseif itTr==3
%                 linIndBkGr = find( ~mainFMask & RM.BFMask{itCl,itPr}~=1 & RM.PoPrFMask{itCl,itPr}~=1 );
%             end
            % Get also 'union of barr and post' as a field type, as well (increase N Spk for trace analysis).
            linIndBarrAndPost = union( linIndBarrF, linIndPostF );
            % Get an index that allows you to also select all spikes
            linIndVisEnv = find(  ~isnan( RM.(['RMProbe' num2str(itPr)]){itCl,1} )   );

            
            % Get phase by position: get overall mean phase for trial, phase in baseline field, phase in trace field.
            % IMPORTANT: make sure order of 'fieldInds' defined below matches with that defined at the preallocation stage
            %            at the top of the function. Always add new field defintions to the end of the list.
            fieldInds  = {linIndVisEnv, linIndBslF, linIndBarrF, linIndPostF, linIndBkGr, linIndBarrAndPost};  
%           fieldNames = {'AllSpk',     'BslF',     'BarrF',     'PostF',     'BkGr'    ,  'BarrPostF'};
            for itFT = 1:length(fieldInds)
                spkInFieldInd                                                    = ismember( spkPosBin, fieldInds{itFT} );
                posInFieldInd                                                    = ismember( binXYLin, fieldInds{itFT} );
                % Measures of phase and phase concentration.
                RM.(['Ph'   fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = circ_mean( spkPh( spkInFieldInd ), [], 1, '+ve' );
                RM.(['RV'   fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = circ_r( spkPh( spkInFieldInd ) );
                RM.(['pR'   fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = circ_rtest( spkPh( spkInFieldInd ) );
                RM.(['pHA'  fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = circ_otest( spkPh( spkInFieldInd ) );
                % Controls for mean rate, speed, etc
                RM.(['NSpk' fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = sum( ~isnan( spkPh( spkInFieldInd ) ) );  % Number of spikes with an actual contributing phase.
                RM.(['MR'   fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = (sum( spkInFieldInd ) / sum( posInFieldInd )) * 50;
                RM.(['Spd'  fieldNames{itFT} '_Pr'   num2str(itPr)])(itCl,itTr)  = nanmean( data.trials( trInds(itTr) ).speed( posInFieldInd ) );
                % Dir maps: sampling of pos and spks is the same principle as above, but now instead 
                % of getting one number per field, we need to get a histogram of dwell/spks in each bin.
                RM.(['DirPos' fieldNames{itFT} '_Pr'  num2str(itPr)]){itCl,itTr} = accumarray( binDir(posInFieldInd), 1, [360/prms.dirBin, 1] );
                RM.(['DirSpk' fieldNames{itFT} '_Pr'  num2str(itPr)]){itCl,itTr} = accumarray( spkDirBin(spkInFieldInd), 1, [360/prms.dirBin, 1] );
                RM.(['DirRte' fieldNames{itFT} '_Pr'  num2str(itPr)]){itCl,itTr} = RM.(['DirSpk' fieldNames{itFT} '_Pr'   num2str(itPr)]){itCl,itTr}  ./   RM.(['DirPos' fieldNames{itFT} '_Pr'   num2str(itPr)]){itCl,itTr};
                RM.(['DirRV'  fieldNames{itFT} '_Pr'  num2str(itPr)])(itCl,itTr) = dir_rayleighvector( RM.(['DirRte' fieldNames{itFT} '_Pr'  num2str(itPr)]){itCl,itTr} );
                % Burstiness measures.
                if prms.runAC && any( strcmp( fieldNames{itFT}, {'AllSpk','BslF'} ) )
                    spkTr       = data.trials( trInds(itTr) ).cells(RM.IndClScn(itCl)).st( spkInFieldInd );
                    % AC.
                    [AC, lags]  = spk_crosscorr(  spkTr,  'AC',  prms.ACBin, prms.ACWin, data.trials( trInds(itTr) ).dur );
                    RM.(['AC'   fieldNames{itFT} '_Pr'   num2str(itPr)]){itCl,itTr}  = AC( lags>0 );
                    % ISI.
                    diffs       = diff( spkTr );
                    diffs       = diffs( diffs<10 );
                    RM.(['ISI' fieldNames{itFT} '_Pr' num2str(itPr)]){itCl,itTr}  = histcounts( diffs, logspace(-3,1,50), 'Normalization', 'probability' );
                end
            end

        end
        
    end
end
delete(hWait);

