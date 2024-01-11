function [] = bvcTrPhPlots( RM, varargin )
% Plots for phase analysis.

% Define barrier-induced doubling.
prms.barrRespScore = 'BFSum';
prms.barrRespThr   = 70;
% Define trace cell.
prms.traceScoreType = 'Pro';   %  '';  %
prms.traceMin       = 0.2;
prms.traceMax       = 2;
prms.traceUsePerc   = [90];
prms.overlapMin     = 0.4;
prms.overlapMax     = 1;
prms.overlapUsePerc = [90];
prms.useNewCell     = 1;
prms.useOldCell     = 1;
% Phase
prms.minNSpk       = 50;
prms.maxPValR      = 0.01;
prms.maxPValHA     = 1;
prms.minRV         = 0;
prms.usePr1Always  = 0;
prms.AZ2use        = [2 3];
prms.useOnlyTetEEG = 0;
prms.phBins        = linspace( 0, 2*pi, 50 );
prms.phBinsInDeg   = round( prms.phBins./pi.*180 );
prms.RVBins        = 0:0.05:0.7;
prms.compToPlot    = [ 7:10 ];
prms.statTypeForDiffs = 'l';  %  'l';   % 'c' ot 'l' for circular or linear.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a list of comparisons to run - 
%     Index, FType1,    CellInd1, Trial1, FType2,   CellInd2,   Trial2; ..
CL  =  { 1, 'BslF',   'trace' ,      1, 'BslF',   'non-trace',    1;  ... % Compare trace and no trace in main field
         2, 'AllSpk', 'trace',       1, 'AllSpk', 'non-trace',    1;  ... % Compare trace and no trace for all spikes fired
         3, 'BslF',   'trace',       2, 'BslF',    'non-trace',   2;  ...    % Do the same as in 1, but now in cue trial.
         4, 'BarrF',   'trace',       2, 'BarrF',   'non-trace',    2;  ...    % .. for this AllSpk not so meaningful, need to split BslF and BarrF
         
         5, 'BslF',   'trace',       2, 'BarrF',    'trace',      2;  ...    % Main field and barrier field firing in barrier trial for trace cells 
         6, 'BslF',  'non-trace',    2, 'BarrF',   'non-trace',   2;  ...    %  .. as above, for non-trace cells.
         
         7, 'BarrF',  'trace',       2, 'BarrF',     'trace',     1;  ...    %  Barrier field, barrier versus pre-trial.
         8, 'BarrF',  'non-trace',   2,  'BarrF',  'non-trace',   1;  ...    %  .. as above, for non-trace cells.
         9, 'BarrF',  'trace',       2, 'BarrF',     'trace',     3;  ...    % As the two above, but now comparing cue with post-cue.
         10, 'BarrF',  'non-trace',  2, 'BarrF',  'non-trace',    3;  ...    %  
         
         11, 'BslF',   'non-trace',   3, 'PostF',  'non-trace',    3;  ...  % Bsl versus Post-Barr, in post-probe trial, non-trace cells only.
         12, 'BslF',  'non-trace',   3, 'BkGr',  'non-trace',    3;  ...      % Bsl versus Background firing, in post-probe trial, non-trace cells only.
         13, 'BarrF',  'trace',       2, 'BarrF',     'trace',     1;  ...    %  Barrier field, barrier versus pre-trial.
         14, 'BarrF',  'non-trace',   2,  'BarrF',  'non-trace',   1;  ...    %  .. as above, for non-trace cells.
         15, 'BarrF',  'trace',       2, 'BarrF',   'non-trace',   2;  ...    %  Barrier field, trace versus non trace.
         
         
        };
CL         =  CL( prms.compToPlot, : );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define trace cells.
RM = bvcTrDefineTrace( RM, prms );

if prms.useOnlyTetEEG
    RM = RM(  RM.EEGIsMatched,  :  );
end

% Get the phase results for the best overall probe. This can either be probe 1 in day
% (to look at running order novelty effects), or probe with best barrier response.
if prms.usePr1Always
    prToUse = ones(height(RM),1);
else
    prToUse = RM.BstPrNum;
end
fieldNames = {'AllSpk',     'BslF',     'BarrF',     'PostF',     'BkGr'}; % ,  'BarrPostF'   };
scoreNames = {'Ph', 'RV'};
for itFT = 1:length(fieldNames)
    for itSc = 1:length(scoreNames)
        
         RM.( [ scoreNames{itSc} fieldNames{itFT} ] ) = nan( height(RM), 3 );

        for itCl = 1:height(RM)


            RM.( [ scoreNames{itSc} fieldNames{itFT} ] )(  itCl,      :      )     = RM.( [  scoreNames{itSc} fieldNames{itFT}  '_Pr' num2str(prToUse(itCl))] )(  itCl,  :  );
            
            cellFilterInd = RM.([  'NSpk'   fieldNames{itFT}  '_Pr' num2str(prToUse(itCl))])(  itCl,  :  )   <   prms.minNSpk      |    ...
                            RM.([  'RV'     'AllSpk'          '_Pr' num2str(prToUse(itCl))])(  itCl,  :  )   <   prms.minRV        |    ...
                            RM.([  'pR'     'AllSpk'          '_Pr' num2str(prToUse(itCl))])(  itCl,  :  )   >   prms.maxPValR     |    ...
                            RM.([  'pHA'    'AllSpk'          '_Pr' num2str(prToUse(itCl))])(  itCl,  :  )   >   prms.maxPValHA;
                        
            RM.( [ scoreNames{itSc} fieldNames{itFT} ] )(  itCl,  cellFilterInd  ) = nan;
            
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots:

% Filter by anatomy.
RM = RM(  any( bsxfun(@eq, RM.anatZone, prms.AZ2use), 2),  :  );



% Define cell selection inds (these needs to be cross-checked against the listing in 'CL', above.
cellInds   = { 'trace',       RM.isTrace;                 ...
               'non-trace',  ~RM.isTrace & RM.isBarrResp; ...
              };

% Pre-allocate a bunch of pairs of samples for paired comparison.
for itCp=1:size(CL,1)
    D(itCp).Ph   = cell(2,1);
    D(itCp).RV   = cell(2,1);
    D(itCp).Diff = cell(1,1);
end
% Now get comparison data.
for itCp = 1:size(CL,1)
    cellInd1      = cellInds{  strcmp( CL{itCp,3}, cellInds(:,1) ), 2 };
    cellInd2      = cellInds{  strcmp( CL{itCp,6}, cellInds(:,1) ), 2 };
    D(itCp).Ph{1} = RM.( ['Ph' CL{itCp,2}] )(  cellInd1,  CL{itCp,4}  );
    D(itCp).Ph{2} = RM.( ['Ph' CL{itCp,5}] )(  cellInd2,  CL{itCp,7}  );
    D(itCp).RV{1} = RM.( ['RV' CL{itCp,2}] )(  cellInd1,  CL{itCp,4}  );
    D(itCp).RV{2} = RM.( ['RV' CL{itCp,5}] )(  cellInd2,  CL{itCp,7}  );
end

% Plot Comparisons %
hFig           = gra_multiplot( length(D), 8, 'figborder', [1.75 0.5 0.5 0.5] );    axArray = getappdata( hFig, 'axesHandles' );
labStr         = cell(1,size(CL,1));
wtnCellPhDiffs = cell( length(D) );   % Within cell phase diffs (calculated *within* each comparison pair) are saved, to be later compared *across* comparison pairs.
for itCp = 1:length(D)

    % Draw phase x RV 2D histograms (one per group)
    phaseXRVHist_SF( D(itCp), axArray(itCp,[1 2]), prms );
    
    % Draw overlapping phase histograms.
    phaseHist_SF( D(itCp), axArray( itCp, 3 ), prms );
    
    % Draw overlapping RV histograms.
    RVHist_SF( D(itCp), axArray( itCp, 4 ), prms );
    
    % Mean vector arrow plots.
    arrowPlot1stOrder_SF( D(itCp), axArray( itCp, 5 ), prms );
    arrowPlot2ndOrder_SF( D(itCp), axArray( itCp, 6 ), prms );

    % Stats
    allStats_SF( D(itCp), axArray( itCp, 7 ), prms );
    
    % Subtrative difference (closest distance) plot.
    wtnCellPhDiffs{itCp} = phaseDiff_SF( D(itCp), axArray( itCp, 8 ), prms );
    
    % Labels for each comparision.
    labStr{itCp} = { ['Gr1: ' CL{itCp,3} 'Tr' num2str(CL{itCp,4}) ' ' CL{itCp,2}  ], 'vs', ...
                     ['Gr2: ' CL{itCp,6} 'Tr' num2str(CL{itCp,7}) ' ' CL{itCp,5}  ] };  
end
gra_multilabel( hFig, 'row', labStr );
titleF = {'maxPValR','usePr1Always','useOnlyTetEEG'};
for itF = 1:length( titleF )
    tStr{itF} = [titleF{itF} '=' num2str( prms.(titleF{itF}) )];
end
gra_multilabel( hFig, 'col',  {'Ph_X_RV group 1', 'Ph_X_RV group 2', 'Ph: Gr1=blue, Gr2=orange', 'RV: Gr1=blue, Gr2=orange',  {'Mean of angles','(1st Order)'}, {'Mean of Vectors','(2nd Order)'}, 'Stats'} );
gra_multilabel( hFig, 'title', tStr );

% Some further stats: compare within pair phase diffs *across* comparison pairs.
phDiffPVal =  nan( length(D) );
for itCp1=1:length(D)
    for itCp2=(itCp1+1):length(D)
        if ~isempty(wtnCellPhDiffs{itCp1}) && ~isempty(wtnCellPhDiffs{itCp2})
            if strcmp( prms.statTypeForDiffs,'l');
                [~,phDiffPVal(itCp1,itCp2),~,~] = ttest2( wtnCellPhDiffs{itCp1}, wtnCellPhDiffs{itCp2} );
            elseif strcmp( prms.statTypeForDiffs,'c');
                [phDiffPVal(itCp1,itCp2),~]   = circ_wwtest( wtnCellPhDiffs{itCp1}, wtnCellPhDiffs{itCp2} );
            end
        end
    end
end
disp( phDiffPVal );


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 2: within tetrode (or session) plot of phase differences between trace and non-trace cells.
RM.DateAsNum = cellfun( @str2double, RM.Date, 'UniformOutput', 1 );
IDParams     = [RM.Rat, RM.DateAsNum, RM.HS];
uniqIDList  = unique( IDParams, 'rows' );
CLwt         = { 1, 'BslF',   'trace' ,      1, 'BslF',   'non-trace',    1 };
Dwt.Ph       = repmat( {nan(size(uniqIDList,1),1)}, 1, 2 );
for itCp = 1:size(CLwt,1)
    cellInd1      = cellInds{  strcmp( CL{itCp,3}, cellInds(:,1) ), 2 };
    cellInd2      = cellInds{  strcmp( CL{itCp,6}, cellInds(:,1) ), 2 };
    for itUT = 1:size(uniqIDList,1)
        indForTet = all(    bsxfun( @eq, IDParams, uniqIDList(itUT,:) ),    2    );
        if any( cellInd1(indForTet)  )    &&    any( cellInd2(indForTet)  )
            Dwt(itCp).Ph{1}(itUT) = circ_mean2ndOrder(    RM.( ['Ph' CL{itCp,2}] )(  cellInd1&indForTet,  CL{itCp,4}  ), RM.( ['RV' CL{itCp,2}] )(  cellInd1&indForTet,  CL{itCp,4}  )   );
            Dwt(itCp).Ph{2}(itUT) = circ_mean2ndOrder(    RM.( ['Ph' CL{itCp,5}] )(  cellInd2&indForTet,  CL{itCp,7}  ), RM.( ['Ph' CL{itCp,5}] )(  cellInd2&indForTet,  CL{itCp,7}  )   );          
%         	  Dwt(itCp).Ph{1}(itUT) = circ_mean(    RM.( ['Ph' CL{itCp,2}] )(  cellInd1&indForTet,  CL{itCp,4}  ), [], 1, '+ve'   );
%             Dwt(itCp).Ph{2}(itUT) = circ_mean(    RM.( ['Ph' CL{itCp,5}] )(  cellInd2&indForTet,  CL{itCp,7}  ), [], 1, '+ve'   );
        end 
    end
end
figure;  plot( Dwt(itCp).Ph{1}, Dwt(itCp).Ph{2}, 'kx' );
hold on;
plot( [2 5], [2 5], 'k:' );
[p,F]      = circ_pairedhotellingtest( Dwt.Ph{1}, Dwt.Ph{2} )

disp('dum');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = phaseHist_SF( D, ax, prms, varargin )
% Plot just phase, two groups overlaid.
if ~isempty(varargin) && strcmp(varargin{1},'diff')
    histBins = linspace( -pi, pi, 50 );
else
    histBins = prms.phBins;
end
for itGr=1:2
    histogram(ax, D.Ph{itGr}, histBins );
    hold( ax, 'on' );
end
plot( ax, [1 1].*(max(histBins)/2), ax.YLim, 'k:' );
set(  ax, 'xlim', [min(histBins) max(histBins)]);
set(  ax, 'XTickLabel', round( (ax.XTick ./ pi) * 180 )  );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = phaseXRVHist_SF( D, axArr, prms)
% Plot just RV, two groups overlaid.
for itGr=1:2
    ax  = axArr(itGr);
    HC2 = histcounts2( D.Ph{itGr}, D.RV{itGr}, prms.phBins, prms.RVBins );
    imagesc( ax, HC2' );
    set( ax, 'XTickLabel', prms.phBinsInDeg(ax.XTick), 'YTickLabel', prms.RVBins(ax.YTick)  );
    set( ax, 'YDir', 'normal' );
    % Plot 180 degrees.
    hold(ax,'on');
    plot(ax,[1 1].*(ax.XLim(2)/2), ax.YLim, 'w:' );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = RVHist_SF( D, ax, prms)
% Plot just RV, two groups overlaid.
for itGr=1:2
    h = histogram(ax, D.RV{itGr}, prms.RVBins );
    hold( ax, 'on' );
end
set(ax, 'xlim', [0 max(prms.RVBins)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = arrowPlot2ndOrder_SF( D, ax, prms )
% Calculate mean RV and angle, using mean-of-means procedure, and plot for each group.
% cSpec = {[100 100 255]./255, [255 165 0]./255};
cSpec = {[0    0.4470    0.7410], [0.8500    0.3250    0.0980]};
for itGr=1:2
    [ma,mr]  = circ_mean2ndOrder( D.Ph{itGr}, D.RV{itGr} );
    [V,U]    = pol2cart( ma, mr );
    hL       = compass( ax, V, U );
    hL.Color = cSpec{itGr};
    hold(ax,'on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = arrowPlot1stOrder_SF( D, ax, prms )
% Calculate mean RV and angle, treating cell means angles as data, and plot for each group.
% cSpec = {'b', 'm'};
cSpec = {[0    0.4470    0.7410], [0.8500    0.3250    0.0980]};
cflm  = {[],[],};
for itGr=1:2
    [ma,cflm{1},cflm{2}] = circ_mean( D.Ph{itGr}, [], 1, '+ve' );
    [V,U]      = pol2cart( ma, 1 );
    compass( ax, V, U, [cSpec{itGr} '-'] );
    hold(ax,'on');
    for itCf = 1:2
        [V,U]    = pol2cart( cflm{itCf}, 1 );
        hL       = compass( ax, V, U, 'k:' );
        hL.Color = cSpec{itGr};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d] = phaseDiff_SF( D, ax, prms )
% Calculate within cell phase diff.
if length(D.Ph{1})~=length( D.Ph{2} )
    d = [];
    return
end
if strcmp( prms.statTypeForDiffs,'l');
    d     = circ_dist( D.Ph{1}, D.Ph{2} );
    d     = (d./pi) .* 180;
    M     = nanmean( d );
    E     = nanstd( d ) / sqrt( sum(~isnan(d)) );
    errorbar(ax, 1, M, E, 'ks');
    set(ax,'ylim', [-45 -10], 'ygrid', 'on');
    d     = d(~isnan(d));
    fprintf( 1, '%4.3f err=%4.3f\n', M, E );
elseif strcmp( prms.statTypeForDiffs,'c');
    d     = circ_dist( D.Ph{1}, D.Ph{2} );
    M     = circ_mean( d );
    E     = circ_confmean(d, 0.05);
    errorbar(ax, 1, M, E, 'ks');
    set(ax,'ylim', [-0.8 -0.2], 'ygrid', 'on');
    d     = d(~isnan(d));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = allStats_SF( D, ax, prms )
% Run all potentially useful stats: 1st versus 2nd order on angles, paired and/or not paired.
%%%
str{  1  } = 'Ind samp 1st Order (Watson-Williams):';
[p,T]      = circ_wwtest( D.Ph{1}(~isnan(D.Ph{1})), D.Ph{2}(~isnan(D.Ph{2})) );
str{end+1} = ['N=' num2str(sum(~isnan(D.Ph{1}))) ',' num2str(sum(~isnan(D.Ph{2}))) ', F = ' num2str(T{2,5},'%4.2f') ', p = ' num2str(p,'%4.3f')];
% str{end+1} = ' ';
%%%
% str{end+1} = 'Ind samp 2nd Order (Batschelet-Hotelling):';
% [p,F]      = circ_2samphotellingtest( D.Ph{1}, D.Ph{2},  D.RV{1}, D.RV{2} );
% str{end+1} = ['F = ' num2str(F,'%4.2f') ', p = ' num2str(p,'%4.3f')];
% str{end+1} = ' ';
%%%
if length( D.Ph{1} ) == length( D.Ph{2} )
    %%%
    str{end+1} = 'Paired Hotelling 1st Order:';
    [p,F]      = circ_pairedhotellingtest( D.Ph{1}, D.Ph{2} );
    str{end+1} = ['N=' num2str(sum(~isnan(D.Ph{1})&~isnan(D.Ph{1}))) ', F = ' num2str(F,'%4.2f') ', p = ' num2str(p,'%4.3f')];
%     str{end+1} = ' ';
    %%%
%     str{end+1} = 'paired Hotelling 2nd Order:'; 
%     [p,F]      = circ_pairedhotellingtest( D.Ph{1}, D.Ph{2}, D.RV{1}, D.RV{2} );
%     str{end+1} = ['F = ' num2str(F,'%4.2f') ', p = ' num2str(p,'%4.3f')];
    %%%
end
text(ax,0,1,str,'verticalalignment','top');
axis(ax,'off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = phaseMeanCI_SF( D, ax, varargin )
% Calculate mean and CI of a bunch of phases, plot as bars in groups.
diffMode = 0;
if ~isempty(varargin) && strcmp(varargin{1},'diff')
    diffMode = 1;
end
D = D';
[M,CI] = deal( nan( size(D) ) );
for itGr=1:numel(D)
    if isempty(D{itGr});   continue;   end
    if diffMode
        M(itGr)  = circ_mean( D{itGr} );        
    else
        M(itGr)  = circ_mean( D{itGr}, [], 1, '+ve'  );
    end
    CI(itGr) = circ_confmean(D{itGr}, 0.05, ones(size(D{itGr})), [], 1);
end
hB = bar( ax, M );
hold(ax, 'on');
if diffMode
    X  = bsxfun( @plus, get(hB,'XData')',  [hB.XOffset] );
else
    X  = bsxfun( @plus, cell2mat(get(hB,'XData')).',  [hB.XOffset] );
end
hL = errorbar(ax, X, M, CI, 'k-');  
[hL.LineStyle] = deal( 'none' );
% Format Plot
if diffMode
    set(ax,'xlim',[0 size(D,1)+1],'ylim',[nanmin(M(:)-CI(:)) nanmax(M(:)+CI(:))].*[0.8 1.2] );
else
    set(ax,'xlim',[0 size(D,1)+1],'ylim',[nanmin(M(:)-CI(:)) nanmax(M(:)+CI(:))].*[0.8 1.2] );
end



