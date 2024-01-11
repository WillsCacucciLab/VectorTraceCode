function [ RM ] = bvcTrPopPlots( RM, varargin )
% Population plots of scores for BVC trace.


prms.barrRespScore  = 'BFSum';
prms.barrRespThr    = 70;
prms.traceScoreType = 'Pro';   %  '';  %
prms.traceMin       = 0.2;
prms.traceMax       = 2;
prms.overlapMin     = 0.4;
prms.overlapMax     = 1;
prms.traceUsePerc   = [90];
prms.overlapUsePerc = [90];
prms.withinCellThr  = 0;
prms.usePr1Always   = 0;
prms.useNewCell     = 1;
prms.useOldCell     = 1;
prms.incBslOnly     = 0;
prms.incNonWallVect = 0;  % Includes *all* cells, including stuff not classified as wall- or vector-responsive by hand. Should be 0 except for prox vs dist phase.
% Phase filtering.
prms.minNSpk       = 50;
prms.maxPValR      = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %


% Define trace cells.
RM = bvcTrDefineTrace( RM, prms );

% Run distance to wall analysis.
% RM = bvcTrDistToBarrAndWall( RM, 'useBFMaskInPostTr', 0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEAL WITH PHASE SEPARATELY: filter by overall phase modulation and nSpks 
% in each field, and then select phase for 'best' trial.
if prms.usePr1Always
    prToUse = ones(height(RM),1);
else
    prToUse = RM.BstPrNum;
end
% fieldNames = {'AllSpk',     'BslF',     'BarrF',     'PostF',     'BkGr'}; % ,  'BarrPostF'   };
% scoreNames = {'Ph', 'RV'};
% for itFT = 1:length(fieldNames)
%     for itSc = 1:length(scoreNames)   
%         RM.( [ scoreNames{itSc} fieldNames{itFT} ] ) = nan( height(RM), 3 );
%         for itCl = 1:height(RM)
%             RM.( [ scoreNames{itSc} fieldNames{itFT} ] )(  itCl,      :      )     = RM.( [  scoreNames{itSc} fieldNames{itFT}  '_Pr' num2str(prToUse(itCl))] )(  itCl,  :  );
%             
%             cellFilterInd = RM.([  'NSpk'   fieldNames{itFT}  '_Pr' num2str(prToUse(itCl))])(  itCl,  :  )   <   prms.minNSpk      |    ...
%                             RM.([  'pR'     'AllSpk'          '_Pr' num2str(prToUse(itCl))])(  itCl,  :  )   >   prms.maxPValR     ;
%                         
%             RM.( [ scoreNames{itSc} fieldNames{itFT} ] )(  itCl,  cellFilterInd  ) = nan;
%         end
%     end
% end

% For convenience, get some 'best probe' values pre-calculated.
% The following scores are all 'one per probe set', and hence are stored in 
% a (nCell,3) array, where column = probe. Getting the best probe is a trivial 
% job, just reference into this array with (pre-calculated) 'RM.BstPrInd'.
scoreNameInd = strncmp( RM.Properties.VariableNames, 'BF', 2 )        |  ...
               strncmp( RM.Properties.VariableNames, 'TrF', 3 )       |  ...
               strncmp( RM.Properties.VariableNames, 'TrF', 3 )       |  ...
               strncmp( RM.Properties.VariableNames, 'Dist', 4 )      |  ... 
               strncmp( RM.Properties.VariableNames, 'Pk_',  3 )      |  ...
               strncmp( RM.Properties.VariableNames, 'Cen_', 4 )     |  ...
               strncmp( RM.Properties.VariableNames, 'TrMR_', 5 )        ;
scoreNameInd( strcmp( RM.Properties.VariableNames, 'BFMask' ) ) = false;   % Shoudl have really called this 'MaskBF'
scoreNames   = RM.Properties.VariableNames( scoreNameInd );
for itSc = 1:length(scoreNames)
    if prms.usePr1Always
        RM.( scoreNames{itSc} )               = RM.( scoreNames{itSc} )( :, 1 );
    else
        RM.( scoreNames{itSc} )               = RM.( scoreNames{itSc} )( RM.BstPrLin );
    end
    % Remove results where there is no barrier response.
    if ~iscell( RM.( scoreNames{itSc} ) )
        RM.( scoreNames{itSc} )( ~RM.isBarrResp ) = nan;
    end
end


% These scores all have 'one score per trial', and are hence stored in 3 sets
% of (nCell,3) arrays, the relevant probe referenced by '_PrN' in the field name.
% Getting best probe more involved, dynamically reference field name in 1:nCell loop.
scoreList = {'MR_v2BarrF','MR_v2BslF','MR_v2AllSpk','MR_v2BkGr', ... 
                'ISIAllSpk', 'ACAllSpk','ISIBslF', 'ACBslF','ISIBarrF', 'ACBarrF', ...
                'ObjNewName', 'ObjType', ...
                'VMAng','VMDist','VMAngExt','VMDistExt', ...
                'IntraStab', 'IntraStabBFMask', 'IntraStabBslFMask'};  %  ,  'EnvProbe',
for itSc=1:length(scoreList)
    if ~any(strcmp( RM.Properties.VariableNames, [scoreList{itSc} '_Pr1'] ) );   continue;   end
    if iscell( RM.([scoreList{itSc} '_Pr1']) )
        RM.(scoreList{itSc}) = cell( size( RM.([scoreList{itSc} '_Pr1']) ) );
    else
        RM.(scoreList{itSc}) = nan( size( RM.([scoreList{itSc} '_Pr1']) ) );
    end
end
if prms.usePr1Always
    prToUse = ones(height(RM),1);
else
    prToUse = RM.BstPrNum;
end
for itCl = 1:height(RM)
    for itSc = 1:length(scoreList)
        if ~any(strcmp( RM.Properties.VariableNames, [scoreList{itSc} '_Pr1'] ) );   continue;   end
        RM.( scoreList{itSc} )(itCl,:) = RM.( [scoreList{itSc} '_Pr' num2str(prToUse(itCl))] )( itCl, : );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots: different plots are in separate subfunctions.
if 0;   RM = vectorFieldScores(RM);                 end
if 0;   RM = pre2PostTraceFieldCorrs(RM);           end
if 1;   RM = intraTrialStabs(RM);                   end
if 0;   RM = compare_BAR_BARCUEBRK(RM);             end
if 0;   RM = traceScoreHistograms(RM);              end
if 0;   RM = meanTraceByAnatZone(RM);               end
if 0;   RM = overlapVersusTraceIncShufStats(RM);    end
if 0;   RM = manualRatingVersusTrace(RM);           end
if 0;   RM = rawRateChangeVersusTrace(RM);          end
if 0;   RM = vectDistTuningHistograms(RM);          end
if 0;   RM = traceFieldShiftToWall(RM);             end
if 0;   RM = percTraceByObjType(RM);                end
if 0;   RM = exportToExcel(RM);                     end
if 0;   RM = ACBurstingPlots(RM);                   end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = vectorFieldScores(RM)
% Plots of vector map scores.
if 1
    cueTypeInd = strcmp(RM.ObjType, '1SqObj' ) | strcmp(RM.ObjNewName, 'Bottle1' );
else
    cueTypeInd = strcmp(RM.ObjType, 'indObjs' ) & ~strcmp(RM.ObjNewName, 'Bottle1' );
end

anatInd = true(height(RM),1); %   strcmp( RM.SUBIC_ZONE, 'DIST' );   %   

CTInd    = { RM.isTrace & cueTypeInd & anatInd, ~RM.isTrace & RM.isBarrResp & cueTypeInd & anatInd};
% CTInd{3} = RM.isBarrResp & cueTypeInd & anatInd;

hFig                   = gra_multiplot( length(CTInd)*3, 5, 'figborder', [1.5 0.5 0.5 0.5] );   axArr = getappdata(hFig, 'axesHandles');
[D_abs, D_bnd, A_abs ] = deal( cell(2,3) );
[D_diff, A_diff ]      = deal( cell(2,2) );
anatKey                = cell(2,1);
barColours             = {[112, 48, 160]./255, [255, 192, 0]./255};
histBinsD              = 0:2.5:60;
for itCT = 1:length(CTInd)
    
    % Plot basic properties.
    for itTr=1:3
        
        % 'Trace' properties of non-trace cells in POST marked in a different colour.
        if itCT==2 && itTr==3
            histFace = [1 1 1].*0.7;
        else
            histFace = barColours{itCT};
        end
                   
        % Distance.
        D_abs{itCT,itTr} = RM.VMDist(CTInd{itCT},itTr);        
        ax               = axArr(itCT,itTr);
        H                = histogram( ax, D_abs{itCT,itTr}, histBinsD );
        H.FaceColor      = histFace;
        D_bnd{itCT,itTr} = H.BinCounts;   % Save for later use - e.g. histogram with cell types interleaved bars.
        
        % Angle.
        A_abs{itCT,itTr} = mod(  (( RM.VMAng(CTInd{itCT},itTr) )./pi).*180,   360   );
        ax               = axArr(itCT+length(CTInd),itTr);
        A                = mod( RM.VMAng(CTInd{itCT},itTr), pi*2 );
        H                = histogram( ax, (A./pi).*180, linspace(0,360,60) );
        ax.XLim          = [0 360];
        H.FaceColor      = histFace;
        
        % Angle & Distance scatter.
        ax  = axArr(itCT+(length(CTInd)*2),itTr);
        axes(ax);
        hPP = polarplot(RM.VMAng(CTInd{itCT},itTr), D_abs{itCT,itTr}, 'bo' );
        set(hPP, 'markersize', 2, 'markeredgecolor', histFace);
        ax  = hPP.Parent;
        set(ax,'RLim',[0 60],'RTickLabel', cat(2, cell(1,length(ax.RTick)-1), {60} ),  ...
              'ThetaTickLabel', {'N','','','E','','','S','','','W','',''}, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');

        % Generate key to anatomy for XL sheet.
        isDistAnat    = strcmp( RM.SUBIC_ZONE, 'DIST' );
        anatKey{itCT} = isDistAnat( CTInd{itCT} );
        
    end
    
    % Plot changes in tuning between trials.
    for itTr = 1:2
        
        % 'Trace' properties of non-trace cells in POST marked in a different colour.
        if itCT==2 && itTr==2
            histFace = [1 1 1].*0.7;
        else
            histFace = barColours{itCT};
        end
        
        % Distance.
        diffDist       = RM.VMDist(:,itTr+1) - RM.VMDist(:,itTr);
        if itCT==1 && itTr==2
            bigDistDiffInd = (find(diffDist<-30 & CTInd{itCT}))'
        end
        D_diff{itCT,itTr} = diffDist( CTInd{itCT} );
        ax                = axArr(itCT,itTr+3);
        if 0
            H                 = histogram( ax, D_diff{itCT,itTr}, -20:2.5:40 );
            H.FaceColor       = histFace;
        else
            X  = RM.VMDist(CTInd{itCT},itTr+1);
            Y  = RM.VMDist(CTInd{itCT},itTr);
            hL = plot( ax, X,  Y, 'bo' );
            set(hL, 'MarkerEdgeColor', histFace, 'markersize', 3 );            
            hold(ax,'on');
            plot(ax,[0 60], [0 60], 'k:');
            set(ax,'xtick',0:20:60, 'ytick',0:20:60 );
        end

        % Angle.
        diffAng  = circ_dist( RM.VMAng(:,itTr), RM.VMAng(:,itTr+1) ); 
        diffAng  = (diffAng ./ pi) .* 180;
        if itCT==1 && itTr==2
            bigAngShiftInd = (find(diffAng>90 & CTInd{itCT}))';  % uncomment this to get the list of cells with big angular shifts
        end
        diffAng = diffAng( CTInd{itCT} );
        ax = axArr(itCT+length(CTInd),itTr+3);
        if 0
            H  = histogram( ax, abs(diffAng), linspace(0,180,30) );
            H.FaceColor = histFace;
            ax.XLim = [0 180];
        else
            X              = (RM.VMAng(CTInd{itCT},itTr+1)./pi).*180;
            Y              = (RM.VMAng(CTInd{itCT},itTr)./pi).*180;
            X( X>360 )     = X( X>360 ) - 360;
            Y( Y>360 )     = Y( Y>360 ) - 360;
            wrI            = abs(X-Y)>180  &  X>Y;
            X( wrI )       = X( wrI ) - 360;
            wrI            = abs(X-Y)>180  &  Y>X;
            Y( wrI )       = Y( wrI ) - 360;
            hL = plot( ax, X, Y, 'bo' );
            set(hL, 'MarkerEdgeColor', histFace, 'markersize', 3 );
            hold(ax,'on');
            plot(ax,[-60 360], [-60 360], 'k:');
            set(ax,'xtick',[0 120 240 360], 'ytick',[0 120 240 360] );
        end
        A_diff{itCT,itTr} = diffAng;   %  abs(diffAng);
        
        % Angle and distance.
        ax       = axArr(itCT+length(CTInd)*2,itTr+3);
        axes(ax);
        clNumInd = find( CTInd{itCT} );
        for itCl=1:length(clNumInd)
            dd = RM.VMDist(clNumInd(itCl),[0 1]+itTr);
            da = RM.VMAng(clNumInd(itCl),[0 1]+itTr);
            if ~any( isnan( [dd da] ) )
                hPP = polarplot(da,dd,'b');
                set(hPP,'linewidth',0.1,'color',histFace);
                hold on;
                hPP = polarplot(da(2),dd(2),'bo');
                set(hPP, 'markersize', 2, 'markeredgecolor',histFace);
            end
        end
        ax  = hPP.Parent;
        set(ax,'RLim',[0 60],'RTickLabel', cat(2, cell(1,length(ax.RTick)-1), {60} ),  ...
              'ThetaTickLabel', {'N','','','E','','','S','','','W','',''}, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    end
end
gra_multilabel( hFig, 'row', {{'Distance','VTC'},{'Distance','Non-trace'},{'Angle','VTC'},{'Angle','Non-trace'}} );
gra_multilabel( hFig, 'col', {'PRE','CUE','POST',{'DIFF','PRE->CUE'},{'DIFF','CUE->POST'}} );

% Summary Mean+-SEM figure: absolute distance, mean angle and distance changes CUE->POST.
hFig = gra_multiplot( 1, 3, 'figborder', [1.5 0.5 0.5 0.5] );   axArr = getappdata(hFig, 'axesHandles');
% a) overall mean distance.
D_m  = cellfun(@nanmean, D_abs);
D_e  = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), D_abs);
hB   = gra_groupedbars( D_m', D_e', axArr(1) );
hB(1).FaceColor = barColours{1};
hB(2).FaceColor = barColours{2};
set(axArr(1), 'xticklabel', {'PRE','CUE','POST'});
title(axArr(1), 'Vector cell distance tuning');
% b) absolute and signed distance shift CUE->POST.
D_diff_abs = cellfun(@abs,D_diff,'UniformOutput',0);
D_m  = [cellfun(@nanmean, D_diff); cellfun(@nanmean, D_diff_abs)];
D_e  = [cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), D_diff); cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), D_diff_abs)];
hB   = gra_groupedbars( D_m, D_e, axArr(2) );
hB(1).FaceColor = barColours{1};
hB(2).FaceColor = barColours{2};
set(axArr(2), 'xticklabel', {'Pr>C','C>Po'});
title(axArr(2), {'Vector Distance Shift', 'Signed Diffs     Abs Diffs'});
% c) absolute angle change CUE->POST.
% D_m  = cellfun(@nanmean, A_diff);
% D_e  = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), A_diff);
D_m  = cellfun(@(x)     (circ_mean(circ_ang2rad(x))./pi).*180, A_diff);
D_e  = cellfun(@(x) (circ_confmean(circ_ang2rad(x))./pi).*180, A_diff);
hB   = gra_groupedbars( D_m', D_e', axArr(3) );
hB(1).FaceColor = barColours{1};
hB(2).FaceColor = barColours{2};
set(axArr(3), 'xticklabel', {'PRE>CUE','CUE>POST'});
title(axArr(3), 'Trace Field Abs Angle Shift');

% Export data to Excel sheet
if 1
    T          = table( cell2mat(D_abs(:,1)),cell2mat(D_abs(:,2)),cell2mat(D_abs(:,3)), 'VariableNames', {'Dist_Tr1','Dist_Tr2','Dist_Tr3'} );
    T.Ang_1    = cell2mat(A_abs(:,1));
    T.Ang_2    = cell2mat(A_abs(:,2));
    T.Ang_3    = cell2mat(A_abs(:,3));
    T.AngDiff_1 = cell2mat(A_diff(:,1));
    T.AngDiff_2 = cell2mat(A_diff(:,2));
    T.AngDiff_3 = circ_rad2ang( circ_dist( circ_ang2rad(cell2mat(A_abs(:,1))), circ_ang2rad(cell2mat(A_abs(:,3))) ) );
    T.CellType = [ones( length(D_abs{1,1}), 1); ones( length(D_abs{2,1}), 1).*2];
    T.AnatKey  = double( cell2mat( anatKey ) );
    writetable( T, 'vector_output.xlsx');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = pre2PostTraceFieldCorrs(RM)
% Plots which compare firing in cue-field in pre-cue versus post-cue trial. 
% Trying to answer Ref1 NN rev1 - are trace fields potentiations of pre-existing firing?
% Z-score the rate maps.
% for itPr = 1:3
%     RM.( ['RMProbe' num2str(itPr)] ) = cellfun( @(x) (x-nanmean(x(:)))./nanstd(x(:)), RM.( ['RMProbe' num2str(itPr)] ), 'UniformOutput', 0 ); 
% end

cellInds = {RM.isTrace, ~RM.isTrace & RM.isBarrResp};
trComps  = 3;
hFig     = gra_multiplot( length(cellInds), length(trComps)*3, 'plotsize', [4 4], 'figborder', [2 1 1 1] );   axArr = getappdata(hFig, 'axesHandles');
for itCI = 1:length(cellInds)
    for itTr = 1:length(trComps)  % Compare pre- trial to cue or post-

        % Get rates in cue and post-cue field, in pre- and post-     
        CFRates   = cell( height(RM), 2, 2 );
        fieldList = {'BFMask','PoPrFMask'};
        for itCl = 1:height(RM)
            for itFd = 1:2
                mask                 = RM.( fieldList{itFd} ){ itCl, RM.BstPrNum(itCl) }==1;
                CFRates{itCl,1,itFd} = RM.(  ['RMProbe' num2str(RM.BstPrNum(itCl))]  ){itCl,1}(  mask  );
                CFRates{itCl,2,itFd} = RM.(  ['RMProbe' num2str(RM.BstPrNum(itCl))]  ){itCl,trComps(itTr)}(  mask  );   
            end
        end
        CFRates = CFRates( cellInds{itCI}, :, : );
        
        % Test Cue field versus Post-cue field corrs.
        F_corr{1} = cellfun( @map_spatialcorr, CFRates(:,1,1), CFRates(:,2,1) );
        F_corr{2} = cellfun( @map_spatialcorr, CFRates(:,1,2), CFRates(:,2,2) );
        for itCr=1:numel(F_corr)
            ind          = F_corr{itCr}>0.999 | F_corr{itCr}<-0.99 | isnan(F_corr{itCr});
            F_corr{itCr} = F_corr{itCr}( ~ind );
        end
        [~,p,~,stats] = ttest2( atanh(F_corr{1}), atanh(F_corr{2}) );
        fprintf( 1, 'Test CF vs PoF, celltype %d, df=%d t=%3.2f p=%4.3f \n', itCI, stats.df, stats.tstat, p );

        % Now go back to data for one field type only (as code above was last-minute hack).
        CFRates = CFRates(:,:,2);   % .. or as a hack insert '2' for analysis of PoF.
        
        % Plots.
        % Plot 1 - correlation between overall means of firing in cue field.
        CF_MR      = cellfun( @nanmean, CFRates );
        nanInd     = any( isnan(CF_MR), 2 );
        [rho,pVal] = corr(CF_MR( ~nanInd, : ));
        ax         = axArr(itCI, ((itTr-1)*3)+1);
        plot(ax, CF_MR(:,1), CF_MR(:,2), 'k.');
        title(ax, sprintf( 'r^2=%3.2f, p=%4.3f' ,rho(2)^2, pVal(2)) );
        hold(ax,'on');
        plot(ax,[0 20], [0 20], 'k:' );
%         set(ax,'ylim',[-1.5 2.5], 'xlim', [-1.5 0]);
        
        % Plot 2 **OLD** - Comparison of MR-in-field correlations (black dashed line; same r-value as in plot 1) 
%         %          with MR-in-field correlations derived from pos-shuffled cue fields.
%         shMR1      = cell2mat( RM.shufCF_MR(cellInds{itCI},1) );
%         shMR2      = cell2mat( RM.shufCF_MR(cellInds{itCI},trComps(itTr)) );
%         shufMRCorr = nan(size(shMR1,1),1);
%         for itSh = 1:size(shMR1,2)
%             nanInd           = isnan( shMR1(:,itSh) ) | isnan( shMR2(:,itSh) );
%             shufMRCorr(itSh) = corr2( shMR1(~nanInd,itSh), shMR2(~nanInd,itSh) );
%         end
%         % Plot shuffled corrs first.
%         ax         = axArr(itCI, ((itTr-1)*3)+2);
%         histogram( ax, shufMRCorr, -0.5:0.05:1, 'Normalization', 'probability' );
%         hold(ax,'on');
%         % Plot real r-value.
%         hold(ax, 'on');
%         plot(ax, [1 1].*rho(2), ax.YLim, 'k:');

        % Plot 2 - r-value of bin-by-bin correlations (across trial, within cue field).
        %          Real data and r-values from pos-shuffled cue fields both plotted.
        CF_corr     = cellfun( @map_spatialcorr, CFRates(:,1), CFRates(:,2) );
        CF_corrShuf = cell2mat( RM.shufPoF_Corr( cellInds{itCI}, trComps(itTr)-1 ) );
        ax          = axArr(itCI, ((itTr-1)*3)+2);
        % Plot shuffled corrs first, looks better visually.
        histogram( ax, CF_corrShuf, -0.5:0.05:1, 'Normalization', 'probability', 'DisplayStyle', 'Stairs' );
        hold(ax,'on');
        % Plot actual data
        histogram( ax, CF_corr, -0.5:0.05:1, 'Normalization', 'probability' );
        hold(ax, 'on');
        % Stats - parametric difference of CF correlations from zero.
        CF_corr( CF_corr>0.999 | CF_corr<-0.999 ) = nan;
        CF_corr                                   = CF_corr( ~isnan(CF_corr) );
        [~,p,~,stats] = ttest( atanh(CF_corr), 0 );
        fprintf( 1, 'Diff to 0: CellType %d, df=%d t=%3.2f p=%4.3f \n', itCI, stats.df, stats.tstat, p );
        % Stats - parametric difference of CF correlations from shuffled mean.
        CF_corrShuf( CF_corrShuf>0.999 | CF_corrShuf<-0.999 ) = nan;
        z_shuf        = atanh( CF_corrShuf(:) );
        [~,p,~,stats] = ttest( atanh(CF_corr), nanmean(z_shuf) );
        fprintf( 1, 'Diff to Shuf: CellType %d, df=%d t=%3.2f p=%4.3f \n', itCI, stats.df, stats.tstat, p );
        
%         [p,~,stats] = signrank( CF_corr, 0 );
%         fprintf( 1, 'CellType %d, Z=%3.2f p=%4.3f \n', itCI,  stats.zval, p ); 
        
        % Plot 3 - boxplot summary of histograms plotted in Plot 2, plus shuffled median.
        ax = axArr(1, ((itTr-1)*3)+3);
        boxplot( ax, CF_corr, 'Notch', 'off', 'positions', itCI, 'symbol', 'o', 'outliersize', 2 );
        hold(ax,'on');
        plot(ax, [-0.3 0.3]+itCI, [1 1].*nanmedian(CF_corrShuf(:)), 'k:');
        set(ax,'xlim',[0 3], 'xtick', [1 2], 'xticklabel', {'VTC','Non-trace'}, 'ylim', [-0.75 1]);
        
        
%         m          = nanmedian(CF_corr);
% %         e          = [-1 1] .* (nanstd(CF_corr)./sqrt(sum(~isnan(CF_corr))));
%         e          = bootci( 100, {@nanmedian, CF_corr} );
%         plot(ax, m, ax.YLim(2)*0.95, 'ks' );
%         plot(ax, e, [1 1].*ax.YLim(2)*0.95, 'k-' );
%         plot(ax, [m m], ax.YLim, 'k:');
%         plot(ax, [1 1].*nanmedian(CF_corrShuf(:)), ax.YLim, 'b:');
%         text( ax, m, ax.YLim(2)*1.05, sprintf( '%3.2f%s%4.3f', m, 177, e(2) ) );
%         
    
    end
end
plotTypeLabs = {'Corr Cue Field Mean Rates','Corr CF MRs: comp to shuf','Within Cue Field Corrs'};
trCompLabs   = {'PRE vs CUE', 'PRE vs POST'};
colCnt       = 1;
for itTC = 1:length(trCompLabs)
    for itPT = 1:length(plotTypeLabs)
        colLabs{colCnt} = {plotTypeLabs{itPT} trCompLabs{itTC}};
        colCnt          = colCnt + 1;
    end
end
gra_multilabel(hFig,'row',{'Trace','Non-Trace'});
gra_multilabel(hFig,'col',colLabs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = intraTrialStabs(RM)
% Plot intra-trial stabilities.
anatInd = true(height(RM),1); % strcmp( RM.SUBIC_ZONE, 'DIST' );
CTInd   = { RM.isTrace & anatInd, ~RM.isTrace & RM.isBarrResp & anatInd};
[M,E]   = deal( nan(length(CTInd), 3 ) );
R       = cell( 1, 2 );
for itCT=1:length(CTInd)
    d           = RM.IntraStabBFMask(CTInd{itCT},:);    % IntraStabBFMask
    M(itCT,1:3) = nanmean(d,1);
    E(itCT,1:3) = nanstd(d,1) ./ sqrt(sum(~isnan(d),1));
    R{itCT}     = d;
end
hB = gra_groupedbars( M, E );
ax = hB.Parent;
ax.YLim = [0 1];
figure;
subplot( 1, 2, 1 );     boxplot( R{1},'notch', 'on' );
subplot( 1, 2, 2 );     boxplot( R{2},'notch', 'on'  );
rawForSPSS = [RM.isTrace( RM.isBarrResp ), RM.IntraStab( RM.isBarrResp , : ) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = compare_BAR_BARCUEBRK(RM)
% Compare cue and trace responses elicited by BAR (as a single trial) and BAR when part of BAR-CUE-BRK.
CTInd{1} = RM.isBarrResp & (RM.ProbeSetType==1 | RM.ProbeSetType==1) & strncmp(RM.ObjNewName,'Barrier',7);
CTInd{2} = RM.isBarrResp & RM.ProbeSetType==2;
for itCT=1:2
    fprintf( 1, '%d/%d. \t', sum( RM.isTrace(CTInd{itCT}) ), sum(CTInd{itCT} ) );
    tr = RM.TrFMeanPro(CTInd{itCT} & RM.isTrace); 
    fprintf( 1, 'Trace=%3.2f \x00B1%3.2f \n', mean(tr), std(tr)/sqrt(length(tr)) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = traceScoreHistograms(RM)
% Histograms of actual and shuffled trace scores.
sc2Plt   = {'TrFMean',  'TrFMeanPro', 'BFPoOverlap' };
histBins = { -0.5:0.05:2,  -0.25:0.05:1.2,    0:0.05:1     };
shufSc   = {'shufTrSc',    'shufTrScPro',   'shufOvLp'    };
hFig     = gra_multiplot( 2, length(sc2Plt) );   axArr = getappdata(hFig, 'axesHandles');
for itSc = 1:length(sc2Plt)  
    % Plain histogram.
    ax = axArr( 1, itSc );
    histogram( ax,  RM.( sc2Plt{itSc} )(RM.isBarrResp), histBins{itSc} );        
    title( ax, sc2Plt{itSc} );
    
    % Inc shuffled scores, both with prob normalisation.
    ax      = axArr( 2, itSc );
    histogram( ax,  RM.( sc2Plt{itSc} ), histBins{itSc}, 'normalization', 'probability' );
    hold( ax, 'on' );
    randAll = cell2mat( RM.( shufSc{itSc} )( RM.isBarrResp )  );  % Note that we are only using the cells for which there is a barrier response.
    p95     = prctile( randAll(:), 95 );
    histogram( ax, randAll(:), histBins{itSc}, 'normalization', 'probability', 'DisplayStyle', 'Stairs' );   
    title( ax, {sc2Plt{itSc}, ['Orange=shuffle, 95th %ile ' num2str(p95,'%3.3f')] } );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = meanTraceByAnatZone(RM)
% Overall mean trace score (including in anat zones).
sc2Plt   = {'TrFMeanPro', 'BFPoOverlap' };
anatInd  = {RM.anatZone==1, RM.anatZone==2 | RM.anatZone==3};  % true(size(RM.anatZone)), 
CTInd    = {RM.isTrace, ~RM.isTrace & RM.isBarrResp};
hFig     = gra_multiplot( 1, length(sc2Plt) );   axArr = getappdata(hFig, 'axesHandles');
for itSc=1:length(sc2Plt)
    [M,E] = deal( nan(length(CTInd), length(anatInd) ) );
    for itCT = 1:length(CTInd)
        for itAn = 1:length(anatInd)
            D = RM.( sc2Plt{itSc} )( CTInd{itCT} & anatInd{itAn} & RM.isBarrResp );
            
            M(itCT,itAn) = nanmean( D, 1 );
            E(itCT,itAn) = nanstd(D,[],1) ./ sqrt( sum(~isnan(D),1) );
            N(itCT,itAn) = sum(~isnan(D),1);
        end
    end
    % Plot
    M  = M';  E = E';  N = N';
    ax = axArr( itSc );
    hB = bar(ax, M);
    X  = bsxfun( @plus, cell2mat(get(hB,'XData')).',  [hB.XOffset] );
    hold(ax, 'on');
    hL = errorbar( ax, X, M, E, 'k-' );
    [hL.LineStyle] = deal( 'none' );   
end
% Pie charts of % cells in each zone is trace
hFig = gra_multiplot( 1, length(anatInd) );   axArr = getappdata(hFig, 'axesHandles');
for itZn = 1:length(anatInd)
    zoneCnt = nan(1,length(CTInd));
    for itCT = 1:length(CTInd)
        zoneCnt(itCT) = sum( anatInd{itZn} & CTInd{itCT} );
    end
    hP = pie( axArr(itZn), zoneCnt );
    for itLb = 2:2:length(hP)
        hP(itLb).String = [hP(itLb).String ' (' num2str(zoneCnt(itLb/2)) ')'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = overlapVersusTraceIncShufStats(RM)
% Overlap versus trace score.
TrType = {'Pro'};
hFig   = gra_multiplot( length(TrType), 2 );   axArr = getappdata(hFig, 'axesHandles');
for itTT = 1:length(TrType)
    trInd  = RM.isTrace;
    nTrInd = ~RM.isTrace & RM.isBarrResp;
    hP = plot( axArr(itTT,1), RM.(['TrFMean' TrType{itTT}])(trInd), RM.BFPoOverlap(trInd), 'or' );
    hP.MarkerEdgeColor = [112, 48, 160]./255;
    hold( axArr(itTT,1), 'on' );
    Y = RM.BFPoOverlap(nTrInd);
    Y( Y==0 ) = Y( Y==0 ) + (rand(sum(Y==0),1)./100 - 0.01);  % Add jitter to points with overlap == 0
    hP = plot( axArr(itTT,1), RM.(['TrFMean' TrType{itTT}])(nTrInd), Y, '.k' );
    hP.MarkerEdgeColor = [255, 192, 0]./255;
    hP.MarkerSize      = 10;
    set( axArr(itTT,1), 'xgrid', 'on', 'ygrid', 'on', 'ylim', [-0.1 1], 'xlim', [-0.2 1.2] );
    xlabel( axArr(itTT,1), {'Trace score', '(Proportion of barrier response)'} );
    ylabel( axArr(itTT,1), {'Overlap score', '(Mean of [field A bins in B] and [field B bins in A])'} );
    title( axArr(itTT,1), { sprintf( 'Mean trace of trace = %3.2f', nanmean( RM.(['TrFMean' TrType{itTT}])(trInd) ) ), sprintf( 'Mean trace of non-trace = %3.2f', nanmean( RM.(['TrFMean' TrType{itTT}])(nTrInd) ) ) } );
    % Compare real data to control, at different overlap and trace thresholds.
    % (Z-test of real N trace cells versus proportion in shuffled (latter treated as pop mean, not sample).
    DTr      = RM.(['TrFMean' TrType{itTT}]);
    DOL      = RM.BFPoOverlap;
    DNanInd  = isnan(DTr) | isnan(DOL);
    DTr      = DTr(  ~DNanInd );
    DOL      = DOL( ~DNanInd );
    ShTr     = cell2mat( RM.( ['shufTrSc' TrType{itTT}] )( ~isnan( RM.TrFMean ) )  );    % Note that we are only using the cells for which there is a barrier response.
    ShOL     = cell2mat( RM.shufOvLp( ~isnan( RM.TrFMean ) )  );                      % Note that we are only using the cells for which there is a barrier response.
    ShNanInd = isnan(ShTr) | isnan(ShOL);
    ShTr     = ShTr( ~ShNanInd );
    ShOL     = ShOL( ~ShNanInd );
    binsOL   = [0.205 0.206]; %(floor( nanmin(DOL).*10 )/10) : 0.1 : (ceil( nanmax(DOL).*10 )/10);
    binsTr   = [0.38  0.381];  % ((floor( nanmin(DTr).*10 )) : 1 : (ceil( nanmax(DTr).*10 ))) ./ 10;    
    [R2ShRatio,ZScore] = deal( nan( length(binsOL)-1, length(binsTr)-1 ) );
    for itTr = 1:(length( binsTr )-1)
        for itOL = 1:(length( binsOL )-1)
            ShPro                = sum( ShTr>=binsTr(itTr) &  ShOL>=binsOL(itOL) )  /   length( ShTr );
            DPro                 = sum( DTr>=binsTr(itTr)  &   DOL>=binsOL(itOL) )  /   length( DTr );
            if ShPro~=0
                R2ShRatio(itOL,itTr) = DPro / ShPro;
                ZScore(itOL,itTr)    = (DPro-ShPro)  /  sqrt(   (ShPro*(1-ShPro))/length(DOL)   );
            end
            if 1 % binsTr(itTr)==0.2 && binsOL(itOL)==0.4
               disp( ZScore(itOL,itTr) )
            end
        end
    end
    % Plot the two-variable real-shuffle comparison generated above
    ax = axArr(itTT,2);
    im = ZScore;
    imagesc( im, 'parent', ax );
    set(ax,'YDir','normal','xtick',1:length(binsTr),'xticklabel',binsTr,'ytick',1:length(binsOL),'yticklabel',binsOL);
    v = version('-release');
    if num2str(v(1:4)) > 2017
    	xtickangle(ax,45);
    end
    title( ax, 'Z stat for real vs shuffled % trace'  );
    colorbar(ax);
    xlabel( ax, 'Threshold Trace score' );
    ylabel( ax, 'Threshold Overlap score' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = manualRatingVersusTrace(RM)
% Manual versus trace score comparison.
hFig              = gra_multiplot( 1, 1 );   axArr = getappdata(hFig, 'axesHandles');
RM.Ratings_Joint  = RM.Ratings_Colin+RM.Ratings_Steve;
RM.BarrTraceJoint = double( ~RM.isTrace );
RM.BarrTraceJoint(~RM.isBarrResp) = 2;
histogram2( axArr(1),  RM.BarrTraceJoint, RM.Ratings_Joint );
set(axArr(1), 'xtick', 0:2, 'xticklabel', {'tr','non-tr','non-barrR'}, 'zlim',[0 80],'CameraPosition',[11.27,17.26,250],'CameraTarget',[1,2,35],'CameraViewAngle',11.4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = rawRateChangeVersusTrace(RM)
% Comparison of 'trace' score versus using 'raw' rate in barrier field to measure trace effect.
% First get a 'raw' (plain mean rate) equivalent of the trace score.
% Need to collapse the 'MR' (or 'MR_2') score from by field*probe 
% analysis to one score per probe set.
fNames = {'TrMR_Post2BarrRatio', 'TrMR_Post2BarrSubPre', 'TrMR_Post2BarrBkGrCorr'};
for itFd = 1:length(fNames);   RM.( fNames{itFd} ) = nan( height(RM), 3 );   end
for itPr = 1:3
    MRBFTr1   = RM.( ['MR_v2BarrF_Pr' num2str(itPr)] )( :, 1 );
    MRBFTr2   = RM.( ['MR_v2BarrF_Pr' num2str(itPr)] )( :, 2 );
    MRBFTr3   = RM.( ['MR_v2BarrF_Pr' num2str(itPr)] )( :, 3 );
    MRBkGrTr1 = RM.( ['MR_v2BkGr_Pr' num2str(itPr)] )( :, 1 );
    MRBkGrTr3 = RM.( ['MR_v2BkGr_Pr' num2str(itPr)] )( :, 3 );
    
    RM.TrMR_Post2BarrRatio(:,itPr)    = MRBFTr3 ./ MRBFTr2;
    RM.TrMR_Post2BarrSubPre(:,itPr)   = (MRBFTr3-MRBFTr1) ./ (MRBFTr2-MRBFTr1);
    RM.TrMR_Post2BarrBkGrCorr(:,itPr) = ((MRBFTr3-MRBFTr1) ./ (MRBFTr2-MRBFTr1))    ./ ( MRBkGrTr3./MRBkGrTr1 ) ;
end
% Plot trace versus mean rate scores.
fNames = {'TrMR_Post2BarrRatio', 'TrMR_Post2BarrSubPre', 'TrMR_Post2BarrBkGrCorr'};
hFig   = gra_multiplot( 1, length(fNames) );   axArr = getappdata(hFig, 'axesHandles');
for itSc = 1:length(fNames)
    
    ax = axArr( 1, itSc );
    plot( ax, RM.TrFMean, RM.( fNames{itSc} ), 'bo' );
    set(ax,'xlim',[-1 2],'ylim',[-1 2])
    hold(ax, 'on' );
    plot( ax, [-1 2], [-1 2], 'k:' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = vectDistTuningHistograms(RM)
% Histogram of vectors, for trace and not trace (and other variables).
dist     = RM.DistBFPk2Barr./4;
cInd1    = { RM.anatZone==1 | RM.anatZone==2 | RM.anatZone==3 };  %  { RM.anatZone==1, RM.anatZone==2 | RM.anatZone==3 };
cInd2    = { ~RM.isTrace & RM.isBarrResp, RM.isTrace};
barColours = {[255, 192, 0]./255, [112, 48, 160]./255};
hFig     = gra_multiplot( 1, length(cInd1)+1 );   axArr = getappdata(hFig, 'axesHandles');
[M,E] = deal( nan(  length(cInd1), length(cInd2) ) );
for itC1 = 1:length( cInd1 )
    for itC2 = 1:length( cInd2 )
        D            = dist( cInd1{itC1} & cInd2{itC2} );
        ax = axArr( itC1 );
        hH = histogram(ax, D, 0:2.5:70, 'Normalization', 'probability');
        hH.FaceColor = barColours{itC2}; 
        hold(ax,'on');
        M(itC1,itC2) = nanmean(D);
        E(itC1,itC2) = nanstd(D)/sqrt(sum(~isnan(D)));
    end
    xlabel(ax, 'Barrier Field Peak to Barrier Distance (cm)' );
    set( ax, 'xlim', [0 70] );
end
% Plot Mean and SEM.
M = fliplr(M);  E = fliplr(E);
ax = axArr(end);
hB = bar( ax, M );
if length(hB)>1
    X  = bsxfun( @plus, cell2mat(get(hB,'XData')).',  [hB.XOffset] );
else
    X  = bsxfun( @plus, get(hB,'XData').',  [hB.XOffset] );
end
hold(ax,'on');
hL = errorbar(ax, X(:), M(:), E(:), 'k-');   hL.LineStyle = 'none';
% [~,pval,~,STATS] = ttest2( D{1}, D{2} );
% title( ax, sprintf( 't(%d)=%3.3f, p=%5.4f', STATS.df, STATS.tstat, pval ) );
% set(ax, 'xticklabel', {'Trace', 'Non-Trace'});
ylabel( ax, 'Dist to barrier (cm)' );
% gra_multilabel( hFig, 'col', {'ANAT=1', 'ANAT=2/3' , ''} );   %  , 'ANAT=3',

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = traceFieldShiftToWall(RM)
%%% Does the post-field shift towards the walls, in trace cells? %%%
% For shuffled data, need to run the analysis through one time first, to get actual shift distances
% between barr and post trials. Then run through a second time, getting shuffled dist to wall/barr,
% where shuffled shifts are drawn from the same distribution as hte real ones.
shiftDist = { RM.DistFShift_Pk, RM.DistFShift_Cen};
RM        = bvcTrDistToBarrAndWall( RM, 'useBFMaskInPostTr', 0, 'BF2BarrFShiftDist', shiftDist );
% Make plots.
fieldPos = {'Pk', 'Cen'};
envObj   = {'Wall', 'Barr', 'Ratio'};
trNames  = {'BF', 'PoPrF', 'ShufMean'};
cSpec    = {'',   'ro',     'rx' };
% Get mean shuffled distance from nShuf raw shuffles. 
for itFP = 1:length(fieldPos)
    for itEO = 1:2
        RM.( ['DistShufMean' fieldPos{itFP} '2' envObj{itEO}] ) = cellfun( @nanmean, RM.( ['DistShufRaw' fieldPos{itFP} '2' envObj{itEO}] ) );
    end
end
% Calculate dist2wall/dist2barr ratio.
for itFP = 1:length(fieldPos)
    for itTr = 1:length(trNames)
        A                                                     = RM.( ['Dist' trNames{itTr} fieldPos{itFP} '2Barr'] );
        B                                                     = RM.( ['Dist' trNames{itTr} fieldPos{itFP} '2Wall'] );
        RM.( ['Dist' trNames{itTr} fieldPos{itFP} '2Ratio'] ) = B  ./ (A + B);
        if any(itTr==[2 3])
        	RM.( ['Dist' trNames{itTr} fieldPos{itFP} '2RatioDiff'] ) = RM.( ['Dist' 'BF' fieldPos{itFP} '2Ratio'] ) - RM.( ['Dist' trNames{itTr} fieldPos{itFP} '2Ratio'] );
        end
    end
end
% Plot: scatter pos in barr trial versus post trial.
hFig     = gra_multiplot( length(fieldPos), length(envObj), 'plotsize', [3 3], 'axesborder', [2 2 2 2] );   axArr = getappdata(hFig, 'axesHandles');
filtInd  = RM.isTrace; % & RM.(['Dist' 'BF' 'Pk' '2Barr'] ) > 40;
for itFP = 1:length(fieldPos)    
    for itEO = 1:length(envObj)
        ax   = axArr( itFP, itEO );
        for itPt = 1:2  % itPt is iterator for plotting BF vs PoPrF, and then BF vs Shuf.
            Y      = RM.( ['Dist' trNames{itPt} fieldPos{itFP} '2' envObj{itEO}] )(filtInd) + 5;  % + half a bin (otherwise bin next to wall counted as 0).
            if ~strcmp( envObj{itEO}, 'RatioDiff' )
                X  = RM.( ['DistBF'    fieldPos{itFP} '2' envObj{itEO}] )(filtInd) + 5;
            else
                X  = (rand(size(Y)).*0.1) + itPt - 1;
            end
            % For 'fieldPos' = 'Pk', need to add jitter so we can see overlapping points better.
            if strcmp( fieldPos{itFP}, 'Pk' )
                X = round( X.*10 )./10;
                Y = round( Y.*10 )./10;
                jitter = [-1 -1; -1 1; 1 -1; 0 1; 0 -1; 1 0; -1 0; 1 1] .* 2;
                uniqXY = unique( [X,Y], 'rows' );
                for itUV = 1:size(uniqXY,1)
                    ind = find( X==uniqXY(itUV,1) & Y==uniqXY(itUV,2) );
                    for itDP = 2:length(ind)
                        X( ind(itDP) ) = X( ind(itDP) ) + jitter(itDP,1);
                        Y( ind(itDP) ) = Y( ind(itDP) ) + jitter(itDP,2);
                    end
                end
            end
            % Convert pixels to cm. 
            if ~strncmp( envObj{itEO}, 'Ratio', 5 )
               X = X./4;   Y = Y./4;   
            end
            % Convert pixels to cm. 
            hL = plot( ax, X, Y, cSpec{itPt} );
            hL.MarkerSize = 4;
            hold( ax, 'on' );
            % Plot mean and SEM, only for RatioDiff.
            if strcmp( envObj{itEO}, 'RatioDiff' )
                M = nanmean(Y);
                E = nanstd(Y) ./ sqrt( sum(~isnan(Y)) );
                errorbar( ax, itPt-0.5, M, E, cSpec{itPt} );
            end
            % Paired T-Test (on real data only)
            if itPt==2 && ~strcmp( envObj{itEO}, 'RatioDiff' )
                [~,pVal,~,S] = ttest(X,Y);
                title( ax, sprintf( 'df=%d, t=%2.3f, p=%4.3f', S.df, S.tstat, pVal ) );
            else
                RatioDiffTTestData{itPt} = Y;
            end
        end
        if strcmp( envObj{itEO}, 'Ratio' )
            plot( ax, [0 1], [0 1], 'k:' );
        elseif ~strcmp( envObj{itEO}, 'RatioDiff' )
            plot( ax, [0 65], [0 65], 'k:' );
        end
        xlabel( ax, {'Distance in','Barrier Trial (cm)'} );   ylabel( ax, {'Distance in','Post Trial (cm)'} );
        set(ax,'xlim',[0 60], 'ylim', [0 60] );
        % Paired T-Test.
        if strcmp( envObj{itEO}, 'RatioDiff' )
            [~,pVal,~,S] = ttest(RatioDiffTTestData{2},RatioDiffTTestData{3});
            title( ax, sprintf( 'df=%d, t=%2.3f, p=%4.3f', S.df, S.tstat, pVal ) );
        end
    end
end
gra_multilabel( hFig, 'col', {'Dist to Wall','Dist to Barrier', {'Dist2Barr', '------------------', '(Dist2Barr + Dist2Wall)'} } );
gra_multilabel( hFig, 'row', fieldPos );

% Plot 2: dedicated plot of mean shift to wall/shift from barrier.
hFig = gra_multiplot( length(fieldPos), 4, 'plotsize', [3 3], 'axesborder', [2 2 2 2] );   axArr = getappdata(hFig, 'axesHandles');
for itFP = 1:length(fieldPos)
    for itEO = 1:2  % Just for 'Wall' and 'Barr', no 'Ratio' type measure.
        rawForStats = cell(1,3);
        for itTr = 2:3  % Just for 'PoPrF' and 'ShufMean'
            
            D                 = RM.( ['Dist' 'BF' fieldPos{itFP} '2'  envObj{itEO}] ) - RM.( ['Dist' trNames{itTr} fieldPos{itFP} '2' envObj{itEO}] );
            D                 = D( filtInd );
            D                 = D ./ 4;
            rawForStats{itTr} = D;
            
            % Plot this raw data in histogram.
            ax = axArr(itFP,((itEO-1)*2)+1);
            histogram( ax, D, -30:2:30 );
            hold(ax,'on');
            plot( ax, [0 0], ax.YLim, 'k-' );
            
            % Plot mean + SEM.
            M  = nanmean(D);
            E  = nanstd(D) ./sqrt( sum(~isnan(D)) );
            ax = axArr(itFP,((itEO-1)*2)+2);
            bar(ax,itTr-1,M,cSpec{itTr}(1));
            hold(ax,'on');
            hL = errorbar( ax, itTr-1, M, E, 'k-' );

        end
        % Label histogram.
        ax = axArr(itFP,((itEO-1)*2)+1);
        xlabel(ax, ['Distance to ' envObj{itEO} ' (Barr-Post; cm)']);
        ylabel(ax,'Count');
        set(ax,'xlim',[-30 30]);
        % Label bar chart.
        ax = axArr(itFP,((itEO-1)*2)+2);
        ylabel(ax, ['Distance to ' envObj{itEO} ' (Barr-Post; cm)']);
        set(ax,'xtick',1:2, 'xticklabel', {'Data','Shuffle'});
        % T-test on bar chart data.
        [~,pVal,~,S] = ttest(rawForStats{2},rawForStats{3});
        title( ax, sprintf( 'df=%d, t=%2.3f, p=%4.3f', S.df, S.tstat, pVal ) );
        
    end
end
gra_multilabel(hFig, 'row', fieldPos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = percTraceByObjType(RM)
% Quantification of env types in which trace cells found.
% Rationalise names a bit: collapse both sizes of barrier to one, and also all independent bottles to one group.
RM.ObjNewName( strncmp(RM.ObjNewName,'Barrier',7) )                    = {'Barrier'};
RM.ObjNewName( RM.ProbeSetType==2 & strcmp(RM.ObjNewName,'Barrier')  ) = {'Barr-BCB'};  %   Separate Barr single trials from those in BAR-CAR-BRK sequences.
RM.ObjNewName( strncmp(RM.ObjNewName,'Bottle',6)  &  ~strcmp(RM.ObjNewName,'BottleBarrier')  ) = {'Bottle'};   % &  ~strcmp(RM.ObjNewName,'Bottle1') 
envList = unique( RM.ObjNewName( RM.isBarrResp ) );
% envList = envList( [5 1 2 4 3] );  % Manual re-sort by %trace resp.
for itEn = 1:length(envList)
    totCell(itEn) = sum( strcmp( RM.ObjNewName( RM.isBarrResp ), envList{itEn} ) );
    traCell(itEn) = sum( strcmp( RM.ObjNewName( RM.isTrace ), envList{itEn} ) );
    txtLab{itEn}  = sprintf( '%d/%d', traCell(itEn), totCell(itEn) );
end
proTraInEnv         = (traCell ./ totCell) .* 100;
[proTrSort,sortInd] = sort(proTraInEnv,'descend');
txtLab              = txtLab(sortInd);
figure;
bar( proTrSort );
text( 1:length(envList), proTrSort+5, txtLab, 'HorizontalAlignment', 'center' );
set( gca, 'XTickLabel', envList(sortInd), 'ylim', [0 80] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = exportToExcel(RM)
% Write results table to excel.
% First any special manipulation of data. 1) Conversion of Pk_ and Cen_ values from 1x2 vectors to separate X and Y entries.
fPosType = {'Pk_','Cen_'};   trialType = {'BF', 'PoPrF'};
if 0
    for itFP = 1:2
        for itTT = 1:2
            temp = cellfun(@transpose, RM.( [fPosType{itFP} trialType{itTT}] ), 'UniformOutput', 0 );
            temp = cell2mat(temp);
            RM.( [fPosType{itFP} trialType{itTT} '_X'] ) = temp(:,1);
            RM.( [fPosType{itFP} trialType{itTT} '_Y'] ) = temp(:,2);
        end
    end
end
% Quick hack: only report cells run in extended return-2-baseline
if 0
    RM.ER2BExists = cellfun( @(x) ~isempty(x), RM.RMExtR2B(:,1) );
    RM.isTrER2B   = RM.isTrace & RM.ER2BTrFPro_maskBF(:,1) >=0.2 & RM.ER2BOvLp_maskBF(:,1)>0.4;
    RM            = RM( RM.isTrER2B & RM.ER2BExists,  :  );
end
% Now define fields of interest, 
field2Keep = {'Rat', 'Date', 'Tet', 'Cell_Num', 'Ratings_Colin', 'Ratings_Steve','SUBIC_ZONE','anatZone', ...
              'BFSum','BFMean',...
              'TrFMean','TrFMeanPro','BFPoOverlap' ...
              'isBarrResp','isTrace','bslDataOnly','Wall_Resp_Vect',...
              'ER2BTrF_maskBF'...
              'DistBFPk2Barr','DistBFPk2Wall','DistBFCen2Barr','DistBFCen2Wall'...
              'Pk_BF_X','Pk_BF_Y','Cen_BF_X','Cen_BF_Y','Pk_PoPrF_X','Pk_PoPrF_Y','Cen_PoPrF_X','Cen_PoPrF_Y', ...
              'PhBarrF','PhBslF','PhAllSpk','PhBkGr',...
              'RVAllSpk', 'pRAllSpk' ...
              'MR_v2BarrF','MR_v2BslF','MR_v2AllSpk','MR_v2BkGr',...
             };
field2Rem  = setdiff( RM.Properties.VariableNames, field2Keep );
% Convert radians to degrees for phase results.
for itFd = 1:length(field2Keep)
    if strncmp(field2Keep{itFd},'Ph',2)
        d                       = (RM.( field2Keep{itFd} ) ./(2*pi)) .* 360;
        d(d<0)                  = d(d<0) + 360;
        RM.( field2Keep{itFd} ) = d;
    end
end
% Remove fields not wanted, then use 'writetable' to autowrite whole table to xl.
for itFd = 1:length(field2Rem)
    RM.( field2Rem{itFd} ) = [];
end
writetable( RM, 'bvcTrOutput.xlsx' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RM] = ACBurstingPlots(RM)
% Plots of burst score
RM.ACAllSpk = cellfun( @transpose, RM.ACAllSpk, 'UniformOutput', 0 );
RM.ACBslF   = cellfun( @transpose, RM.ACBslF, 'UniformOutput', 0 );
RM.ACBarrF  = cellfun( @transpose, RM.ACBarrF, 'UniformOutput', 0 );
CTInds      = {RM.isTrace, ~RM.isTrace&RM.isBarrResp};
cSpec       = {[0    0.4470    0.7410], [0.8500    0.3250    0.0980]};
scNms       = {'AC', 'ISI'};
scNms       = scNms( [1 2] );
tr2Plt      = [1 2 3];
hFig        = gra_multiplot( length(CTInds), length(tr2Plt)*length(scNms), 'figborder', [1 0.5 0.5 0.5] );   axArr = getappdata(hFig, 'axesHandles');
for itTr = 1:length(tr2Plt)
    for itCT = 1:length( CTInds )
        for itSc = 1:2
            
        	A = cell2mat( RM.( [scNms{itSc} 'BslF'] )( CTInds{itCT}, tr2Plt(itTr) ) );
            
            A = bsxfun( @rdivide, A, nansum(A,2) );
            
            A(A>0.2) = nan;
            
            if itSc==1
                [~,pkRow]  = nanmax( A, [], 2 );
                [~,srtInd] = sort(pkRow);
%                 [~,srtInd] = sort( A(:,5) );
            end
            A          = A(srtInd,:);
            
            ax = axArr( itCT, ((itTr-1)*length(scNms))+itSc );
            imagesc(ax, A );
            
            if strcmp( scNms{itSc}, 'ISI' )
                bins = logspace( -3, 1, 50 );
                bins = (round(bins.*100))./100;
            else
                bins = (1:20) ./ 1000;
            end
            ax.XTickLabel = bins( ax.XTick );
            
            if itCT==2
                xlabel(ax,scNms{itSc});
            end
            if itTr==1 && itSc==1
                ylabel( ax, 'Cell #' );
            end
            
%             % Histograms of burst score.
%             AC RM.BstBrst( CTInds{itCT}, itTr )
% 
%             ax = axArr(itTr,1);
%             histogram( ax, RM.BstBrst( CTInds{itCT}, itTr ), 0:0.1:3, 'normalization', 'probability' );
%             hold(ax, 'on');
%             set(ax,'xlim',[0 3]);
%             % Mean AC plots.
%             ax    = axArr(itTr,2);
%             acArr = cell2mat( RM.BstAC( CTInds{itCT}, itTr ) );
%             if 1
%                 acArr = bsxfun( @rdivide, acArr, nanmean(acArr,2) );
%             end
%             mAC   = nanmean( acArr );
%             stdAC = nanstd( acArr );
%             hL = plot( ax, 1:size(acArr,2), mAC, 'b-' );    hL.Color = cSpec{itCT};
%             hold( ax, 'on' );
%             hL = plot( ax, 1:size(acArr,2), mAC+stdAC, 'b:' );    hL.Color = cSpec{itCT};
%             hL = plot( ax, 1:size(acArr,2), mAC-stdAC, 'b:' );    hL.Color = cSpec{itCT};
        
        end
    end
end
gra_multilabel( hFig, 'row', {'Trace', 'Non-trace'});
gra_multilabel( hFig, 'col', {'Pre', '', 'Cue', '', 'Post'} );





