function [rtn] = bvcTrGetBestThetaEEGForTet(RM)
% Run through scan data EEGs, get that with the best theta signal on each tet (best SNR, mean across trials), save
% that one and delete the rest. Mark the preserved EEGs with '.tet' for convenience.

SD     = gss;

% Need to get a list of which tetrodes in which anatomical zone.
CL      = readtable(  'VTC_Metafile_8_020720.xlsx'  );   %  'DOUBLING_BVCsONLY_TraceCellMetaFile_TW_col_names.xlsx'  'Vector Trace Cells Metafile 21-02-19 1808 TW SLP.xlsx'
ratList = unique( CL.Rat );
for itRt = 1:length(ratList)
   CLRat = CL( CL.Rat==ratList(itRt), : );
   tetList = unique( CLRat.Tet );
   for itTt = 1:length( tetList )
       tetAZ{ itRt }{1}( itTt ) = tetList(itTt);
       tetAZ{ itRt }{2}{ itTt } = CLRat.SUBIC_ZONE{ find(CLRat.Tet==tetList(itTt), 1 ) }; 
   end 
end


[snrAll, perRailAll] = deal( nan( 0, 64 ) );
[trialAll, DSIndAll] = deal( nan( 0, 1 ) );
for itDS=1:length( SD.selData )
    
    data = evalin( 'base', SD.selData{itDS} );
    disp( SD.selData{itDS} );
    
    % Get a list of the actual trials used in probe analysis.
    DSInRMInd    = RM.Rat==str2double(SD.selData{itDS}(2:4))  &  strcmp( RM.Date, SD.selData{itDS}(6:end) );
    totTrialList = [RM.IndScnPr1( DSInRMInd, : );   RM.IndScnPr2( DSInRMInd, : );   RM.IndScnPr3( DSInRMInd, : )];
    totTrialList = unique( totTrialList( ~isnan(totTrialList) ) );
    
    % Run an FFT on every EEG, on every trial, to get the theta SNR.
    [snr, perRail, thFreq] = deal( nan( length(data.trials), 64 ) );
    if ~isfield( data.trials(  totTrialList(1)  ).eeg(1), 'thetaFreq' )
        runFFT = 1;    
    else
        runFFT = 0;    
    end
    for itTr = 1:length(data.trials)
        for itEeg = 1:length( data.trials(itTr).eeg )
            if isempty( data.trials(itTr).eeg(itEeg).eeg );   continue;   end
            if ~any( itTr == totTrialList );                  continue;   end
            
            %%%%% Metadata for recording tet. This all goes in trial 1. %%%%
            % Define tet and hemisphere for signal source of this EEG.
            data.trials(1).eeg(itEeg).tet      = ceil(data.trials(1).eeg(itEeg).channel / 4);
            if data.trials(1).eeg(itEeg).tet > 16
                data.trials(1).eeg(itEeg).tet = data.trials(1).eeg(itEeg).tet - 16;  % Design of 128-DACQ was that ch 1-64 could be 'mirrored' onto 65-128. This is how EEGs were recorded in some rats.
            end
            % Define hemisphere.
            if any( str2double(SD.selData{itDS}(2:4)) == [462 463] )
                data.trials(1).eeg(itEeg).HS   = ceil( data.trials(1).eeg(itEeg).tet  /  4 );
            else
                data.trials(1).eeg(itEeg).HS   = ceil( data.trials(1).eeg(itEeg).tet  /  8 );
            end
            % Mark anatomical zone for this tetrode.
            ratAZ                                  = tetAZ{ str2double(SD.selData{itDS}(2:4))==ratList };
            tetIndInAZKey                          = data.trials(1).eeg(itEeg).tet==ratAZ{1};
            if ~any(tetIndInAZKey)
                data.trials(1).eeg(itEeg).anatZone = 4;  % This is for the case of EEGs that exist in scan, but there are no cells on that tet in the master sheet. Don't know their anatomy, therefore.
            else
                data.trials(1).eeg(itEeg).anatZone = find( strcmp( ratAZ{2}{tetIndInAZKey}, {'PROX','MID','DIST','NONE'} ) );  % Convert to numerical code
            end
            
            %%%% SNR etc for each trial each tetrode.
            % Get theta SNR (run FFT if necessary, ifnot run before).
            if runFFT
                [thFreq(itTr,itEeg), ~, snr(itTr,itEeg)] = eeg_powerspec(data.trials(itTr).eeg(itEeg), 'hfCutOff', 25, 'freqRes', 0.02, 'fig', 0);
                perRail(itTr,itEeg)                      = sum( data.trials(itTr).eeg(itEeg).eeg==127  |  data.trials(itTr).eeg(itEeg).eeg==-128 ) / length( data.trials(itTr).eeg(itEeg).eeg ) * 100;
                data.trials(itTr).eeg(itEeg).thetaSNR    = snr(itTr,itEeg);
                data.trials(itTr).eeg(itEeg).perRail     = perRail(itTr,itEeg);
                data.trials(itTr).eeg(itEeg).thetaFreq   = thFreq(itTr,itEeg);
            else
                snr(itTr,itEeg)                       = data.trials(itTr).eeg(itEeg).thetaSNR;
                perRail(itTr,itEeg)                   = data.trials(itTr).eeg(itEeg).perRail;
                thFreq(itTr,itEeg)                    = data.trials(itTr).eeg(itEeg).thetaFreq;
            end
            
        end
    end

    
    % Find the minimum SNR and maximum %-on-rails over all trials - this will be used to pick an EEG channel. Alao mean theta frequency, for checking.
    trMinSNR     = nanmin( snr, [], 1 );
    trMaxPerRail = nanmax( perRail, [], 1 );
    meanThFreq   = nanmean( thFreq, 1 );
    for itEeg = 1:length(data.trials(1).eeg)
        data.trials(1).eeg(itEeg).trMinSNR     = trMinSNR(itEeg);
        data.trials(1).eeg(itEeg).trMaxPerRail = trMaxPerRail(itEeg);
        data.trials(1).eeg(itEeg).meanThFreq   = meanThFreq(itEeg);
    end
    
    
    % Put all the data together to make overall characterisation plots.
    snrAll     = cat( 1, snrAll, snr );
    perRailAll = cat( 1, perRailAll, perRail );
    trialAll   = cat( 1, trialAll, (1:length(data.trials))' );
    DSIndAll   = cat( 1, DSIndAll, ones(length(data.trials),1).*itDS );
    
    % Get the list of tets for each EEG.
    tetAll    = cat( 1, data.trials(1).eeg.tet);
    HSAll     = cat( 1, data.trials(1).eeg.HS);
    azAll     = cat( 1, data.trials(1).eeg.anatZone);

    
    % For each tet, find the index of the EEG channel with the best SNR, both on each tetrode, and within each hemisphere.
    eegToKeep = [];
    for itEeg = 1:length( data.trials(1).eeg )
        data.trials(1).eeg(itEeg).bestForTet  = 0;
        data.trials(1).eeg(itEeg).bestForHS   = 0;
        data.trials(1).eeg(itEeg).bestForZone = 0;
        data.trials(1).eeg(itEeg).bestForDist = 0;
        data.trials(1).eeg(itEeg).bestForSess = 0;
    end
    % Best for tet ..
    for itTet = 1:16
        tetInd =  tetAll == itTet;
        if any( tetInd )
            tetTmpSNR                              = trMinSNR;
            tetTmpSNR( ~tetInd )                   = nan;
            tetTmpSNR( trMaxPerRail>1 )            = nan;
            [~, bestEEG]                           = nanmax( tetTmpSNR );
            data.trials(1).eeg(bestEEG).bestForTet = 1;
            eegToKeep                              = [eegToKeep, bestEEG]; %#ok<AGROW>
        end
    end
    % Best for anat zone (matched by HS) ..
    for itHS = 1:2
        for itAZ = 1:3
            azInd =  azAll==itAZ & HSAll==itHS;
            if any( azInd )
                azTmpSNR                                = trMinSNR;
                azTmpSNR( ~azInd )                      = nan;
                azTmpSNR( trMaxPerRail>1 )              = nan;
                [~, bestEEG]                            = nanmax( azTmpSNR );
                data.trials(1).eeg(bestEEG).bestForZone = 1;
            end
        end
    end
    % Best for distal+mid zone (for each HS) ..
    for itHS = 1:2
        azInd =  (azAll==2 | azAll==3)  & HSAll==itHS;
        if any( azInd )
            azTmpSNR                                = trMinSNR;
            azTmpSNR( ~azInd )                      = nan;
            azTmpSNR( trMaxPerRail>1 )              = nan;
            [~, bestEEG]                            = nanmax( azTmpSNR );
            data.trials(1).eeg(bestEEG).bestForDist = 1;
        end
    end
    % Best for HS ..
    [SNRByHS, bestEEGByHS] = deal( nan(1,2) );
    for itHS = 1:2
        HSInd =  HSAll == itHS;
        if any( HSInd )
            HSTmpSNR                              = trMinSNR;
            HSTmpSNR( ~HSInd )                    = nan;
            HSTmpSNR( trMaxPerRail>1.5 )          = nan;
            [bestSNR, bestEEG]                    = nanmax( HSTmpSNR );
            data.trials(1).eeg(bestEEG).bestForHS = 1;
            bestEEGByHS(itHS)                     = bestEEG;
            SNRByHS(itHS)                         = bestSNR;
        end
    end
    % Best for expt ..
    [~,bestHSForSess] = nanmax( SNRByHS );
    data.trials(1).eeg(   bestEEGByHS(bestHSForSess)   ).bestForSess = 1;
    
    
    % For every trial in dataset, select only EEGs to be kept, and then reassign data to workspace.
%     for itTr = 1:length(data.trials)
%         data.trials(itTr).eeg = data.trials(itTr).eeg( eegToKeep );
%         for itEeg = 1:length(data.trials(itTr).eeg)
%             data.trials(itTr).eeg(itEeg).eeg = int8( data.trials(itTr).eeg(itEeg).eeg );
%         end
%     end
    
    
    assignin( 'base', SD.selData{itDS}, data );
    
end

rtn.snr     = snrAll;
rtn.perRail = perRailAll;
rtn.trial   = trialAll;
rtn.DSInd   = DSIndAll;
rtn.DSName  = SD.selData;