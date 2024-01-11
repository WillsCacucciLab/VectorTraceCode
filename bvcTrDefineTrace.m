function [RM] = bvcTrDefineTrace( RM, varargin )
% Define trace cells.
% Currently just does population thresholding, after filtering for whether there is a barrier field or not.

prms.barrRespScore  = 'BFSum';
prms.barrRespThr    = 70;
prms.traceScoreType = 'Pro';   %  '';  %
prms.traceMin       = 0.2;
prms.traceMax       = 100;
prms.traceUsePerc   = [];
prms.overlapMin     = 0.4;
prms.overlapMax     = 100;
prms.overlapUsePerc = [];
prms.withinCellThr  = 0;
prms.write2Excel    = 0;
prms.useNewCell     = 1;
prms.useOldCell     = 1;
prms.incBslOnly     = 0;
prms.incNonWallVect = 0;  % Includes *all* cells, including stuff not classified as wall- or vector-responsive by hand. Should be 0 except for prox vs dist phase.
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
% Define barrier responsive cells.
bstBarrScore  = nanmax( RM.( prms.barrRespScore ), [], 2 );
barrRespInd   = bstBarrScore > prms.barrRespThr;
% We can exclude 'new' cells (those added in Feb/Mar 2020 for Nat Neurosci revision 1), so that
% we can sanity check if we still get the original results, using later pipeline. Or even just plot those cells.
% Practically, we do this by setting 'isBarrResp' to 0 where relevant, so cells remain in the strucuture in order.
if ~prms.useNewCell
    barrRespInd   = barrRespInd & ~RM.NEW_CELL;
elseif ~prms.useOldCell
    barrRespInd   = barrRespInd & RM.NEW_CELL;
end
% If requested, can mark pre-probe bsl only data is barrier responsive. This means that the pre-probe
% data for this is included in analysis. Probe and post-probe scores should always be NaN.
RM.bslDataOnly = cellfun( @(x) ~isempty(x), RM.RMProbe1(:,2) );
if prms.incBslOnly
    barrRespInd   = barrRespInd | RM.bslDataOnly;
end
% Set non-wall vector cells to 'non-barrier responsive', so they are not included ...
% unless option to include is active, e.g. for prox vs dist baseline phase.
if ~prms.incNonWallVect && any(strcmp( RM.Properties.VariableNames, 'Wall_Resp_Vect' ))
    barrRespInd( ~RM.Wall_Resp_Vect ) = false;
end
% Set 'isBarrResp' variable in output table.
RM.isBarrResp = barrRespInd;


% Get trace and overlap scores in convenient format.
traceScore     = RM.( ['TrFMean' prms.traceScoreType] )( RM.BstPrLin );
overlapScore   = RM.BFPoOverlap( RM.BstPrLin );

if ~isempty( prms.traceUsePerc )
    if prms.withinCellThr
        prms.traceMin = cellfun( @(x) prctile(x, prms.traceUsePerc), RM.( ['shufTrSc' prms.traceScoreType] ) );
    else
        prms.traceMin = prctile(  reshape(cell2mat(RM.( ['shufTrSc' prms.traceScoreType] )(barrRespInd)), 1, []), prms.traceUsePerc );
        fprintf(1, 'Trace thr (%dth%%-ile) = %4.3f, \t', prms.traceUsePerc, prms.traceMin );
    end
    prms.traceMax = repmat( 100, size( prms.traceMin ) );
end
if ~isempty( prms.overlapUsePerc )
    if prms.withinCellThr
        prms.overlapMin = cellfun( @(x) prctile(x, prms.overlapUsePerc), RM.shufOvLp );
    else
        prms.overlapMin = prctile(  reshape(cell2mat(RM.shufOvLp(barrRespInd)), 1, []), prms.overlapUsePerc );
        fprintf(1, 'Ovelap thr (%dth%%-ile) = %4.3f, \t', prms.overlapUsePerc, prms.overlapMin );
    end
    prms.overlapMax = repmat( 100, size( prms.overlapMin ) );
end

RM.isTrace               = traceScore>=prms.traceMin  &  traceScore<=prms.traceMax  &  overlapScore>=prms.overlapMin  &  overlapScore<=prms.overlapMax;
RM.isTrace(~barrRespInd) = false;

fprintf( 1, 'N barrier responsive cells = %d, N trace cells = %d (%d%%) \n', sum(RM.isBarrResp), sum(RM.isTrace), round(sum(RM.isTrace)/sum(RM.isBarrResp)*100) );

% Save trace thresholds in .Properties (so calling functions have access to percentile thresholds).
RM.Properties.UserData.CurrentThrTrace   = prms.traceMin;
RM.Properties.UserData.CurrentThrOverlap = prms.overlapMin;

% Write to excel file.
if prms.write2Excel
    RM2Wr = RM;
    f2Rem = setdiff( RM.Properties.VariableNames, {'Rat','Date','Tet','Cell_Num','Ratings_Colin','Ratings_Steve','isTrace','isBarrResp'} );
    for ii=1:length(f2Rem)
        RM2Wr.( f2Rem{ii} ) = [];
    end
    writetable( RM2Wr, 'traceCellDefs.xlsx', 'FileType', 'spreadsheet' );
end


% RM.isTracePro(~barrRespInd) = false;
    
% RM.isTrace     = RM.TrFMean( RM.BstPrLin )    > prctile(  reshape(cell2mat(RM.shufTrSc(barrRespInd)), 1, []), 95 );
% RM.isTracePro  = RM.TrFMeanPro( RM.BstPrLin ) > prctile(  reshape(cell2mat(RM.shufTrScPro(barrRespInd)), 1, []), 95 );
