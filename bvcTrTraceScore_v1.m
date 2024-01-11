function [BFMean,BFSum,TrFMean,TrFSum,BFMask,PoPrFMask,BslFMask] = bvcTrTraceScore_v1( M, varargin )
% Version 1 of trace score, find the barrier field simply as a large-enough
% contiguous block of firing.
%
% Input is a cell array of 3 maps, {pre-bsl, probe, post-bsl}

% Analysis parameters.
prms.bstProbeDef   = 'BFSum';
prms.mapNormMode   = 'norm2Pk';   %   'Z';
prms.BslFieldThr   = 0.5;      % 'Pk/2';
prms.BarrFieldThr  = 0.5;      % 'Pk/2';
prms.PostFieldThr  = 0.5;
prms.subBslFiring  = 0;
prms.BFInvalidIfBslFOverlap = 1;        % If 1, barrF regions which overlap with main field are excluded outright. If 0, main field pixels are still removed from barr F regions, but then these can still count as barr fields.
prms.BFSortParam            = 'Area';   %  'MeanIntensity', 'Area'  % Select 'the' barrier field from many by size or summed rate? (Param names reflect regionprops arguments). 
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


[BFMean,BFSum,TrFMean,TrFSum] = deal(nan(1,1));
[BFMask,PoPrFMask,BslFMask]   = deal( [] );

% If the 'M' input (maps) is empty, it means input is coming from a probe not run. Return the values defined above. 
if isempty( M{1} );   return;    end

% Define 'baseline' field (main field pre-probe.
BslFMask   = bvcTrDefineMainField( M{1}, prms );  % Note raw rate map passed, normalised in that function depending on mode.

% If M{2} is empty, this is a cell without a probe - process BslFMask for output, and return (all other outputs are NaN).
if isempty( M{2} );   
    BslFMask              = double( BslFMask );  % For consisteny, make BslFMask the same format os the other masks,
    BslFMask(isnan(M{1})) = 3;                   % i.e. double, where 1 is field, 0 no field and 3 unvisited.
    return   
end

% Make Z-scored, pre-probe subtracted rate maps. 
if strcmp( prms.mapNormMode, 'Z' )
    normM  = cellfun(    @(x) (x-nanmean(x(:))) ./ nanstd(x(:)),     M, 'UniformOutput', 0 );  % Z-score the maps
elseif strcmp( prms.mapNormMode, 'norm2Pk' )
    normM  = cellfun(    @(x) x./nanmax(x(:)),     M, 'UniformOutput', 0 );  % Normalise to peak rate.
end
Ms = { normM{2}-normM{1},  normM{3}-normM{1} };                                                   % Subtract pre-bsl from probe and post.

% Define 'new fields' (new as compared to pre-probe) in cue and post-probe trials.
bfMaskFin = cell(1,2);
for itTr=1:2
    
    % !!!This is the definition of the barrier and post-barrier fields!!!
    %%% TODO - this is assuming firing subtracted, this has to be the case, to remove cells with pre-barr fields, but then what happens to thr?
    if itTr==1
        fMask = Ms{itTr} >= prms.BarrFieldThr;
    else
        fMask = Ms{itTr} >= prms.PostFieldThr;
    end
    
    
    % Remove 'bsl field' bins from barr field mask (as long as we are not just going to exclude them completely, which will happen below).
    if ~prms.BFInvalidIfBslFOverlap
        fMask( BslFMask ) = false;
    end

    % Label (so as to facilitate bsl F overlap check).
    fLabels = bwlabel( fMask, 4 );
        
    % If requested, exclude as invalid any contiguous area that overlaps the bsl field.
    if prms.BFInvalidIfBslFOverlap
        labelList = setdiff( unique(fLabels), 0 );
        for ii=1:length(labelList)
            tempMask = fLabels==labelList(ii);
            if any( tempMask(:) & BslFMask(:) )
                fLabels( tempMask ) = 0;
            end
        end
    end
    
    % If no barrier field, no trace defined. Return an 'all-false' field mask in this case.
    if all( fLabels(:)==0 );  bfMaskFin{itTr}=zeros( size(fMask) );  continue;   end
    
    %%% TODO: sort by summed firing, not area?
    RP          = regionprops( fLabels, Ms{itTr}, 'Area', 'MeanIntensity' );
    [~,sortInd] = sort( cat(1, RP(1:end).( prms.BFSortParam )), 'descend' );
    bfLabelFin  = sortInd(1);
    
    bfMask                                     = zeros( size(fMask) );  % Default = 0
    bfMask( fLabels==bfLabelFin )              = 1;                     % The actual used mask = 1
    bfMask( fLabels~=bfLabelFin & fLabels~=0 ) = 2;                     % Other fields, but not that chosen as mask are 2.
    bfMask( isnan(normM{3}) )                  = 3;                     % Non-visited = 3.
    
    bfMaskFin{itTr}                            = bfMask;
end
BFMask                    = bfMaskFin{1};
PoPrFMask                 = bfMaskFin{2};
BslFMask                  = double( BslFMask );  % For consisteny, make BslFMask the same format os the other two,
BslFMask(isnan(normM{3})) = 3;                   % i.e. double, where 1 is field, 0 no field and 3 unvisited.

if ~any( BFMask(:) );   return;   end  % Don't calculate trace if there is no barrier field.

% Barrier field and trace scores are diff in mean Z between pre and post, within BF(probe) mask. Can be mean or summed firing rate.
BFMean  = mean( normM{2}( BFMask==1 ) );
BFSum   = sum( normM{2}(  BFMask==1 ) );
TrFMean = mean( normM{3}(  BFMask==1 ) );
TrFSum  = sum( normM{3}(  BFMask==1 ) );
if prms.subBslFiring
    BFMean  = BFMean  - mean( normM{1}( BFMask==1 ) );
    BFSum   = BFSum   -  sum( normM{1}( BFMask==1 ) );
    TrFMean = TrFMean - mean( normM{1}(  BFMask==1 ) );
    TrFSum  = TrFSum  -  sum( normM{1}(  BFMask==1 ) );
end