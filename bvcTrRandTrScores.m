function [randMeanZ, randOvLp, randMeanR, randCFCorr] = bvcTrRandTrScores(M,BFField,PoPrField,varargin)
% Generate a 'shuffled' BVC trace score, by random sub-sampling of a matching number of
% bins from the return to baseline trial.
% 
%   [randMeanZ, randOvLp] = bvcTrRandTrScores( {RM_pre, RM_PoPr}, BFField, PoPrField, prms );
% 
% Rate Maps should be raw, not Z-scored yet. 'BFField' and 'PoPrField' should be the 'label' (numerical) masks,
% so the logical mask for the field is found via BFField==1, for example.
%
% Outputs are vectors of size (1,prms.nShuf).

% Analysis parameters.
prms.nShuf        = 100;
prms.shufType     = 'shiftRotBF';  % 'randBins';  %
prms.subBslFiring = 0;
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

%%%% IMPORTANT - hacked here to test the post-probe, instead of cue, field.
if 1
    fieldMaskToTest = BFField;
else
    fieldMaskToTest = PoPrField;
end

if strcmp( prms.mapNormMode, 'Z' )
    normM   = cellfun( @(x) (x-nanmean(x(:))) ./ nanstd(x(:)),   M, 'UniformOutput', 0 );  % Z-score the maps
elseif strcmp( prms.mapNormMode, 'norm2Pk' )
    normM   = cellfun(    @(x) x./nanmax(x(:)),   M, 'UniformOutput',  0  );  % Z-score the maps
end
normDiffMap = normM{3}-normM{1};    % Subtract pre-bsl from post.
normMap     = normM{3};             %  .. or not. The version we use is selected below, depending on setting of prms.subBslFiring

% mainFMask = M{1}>(nanmax(M{1}(:))*0.5);
mainFMask = bvcTrDefineMainField( M{1}, prms );
mapVisNF  = normDiffMap( ~isnan(M{3}) & ~mainFMask );
nBinInF   = sum( fieldMaskToTest(:)==1 );

if strcmp( prms.shufType, 'randBins' )
    randInd   = nan(  nBinInF, prms.nShuf  );
    for ii=1:prms.nShuf
       randInd(:,ii) = randperm(  length(mapVisNF), nBinInF  )';
    end
    randZVals = mapVisNF( randInd );
    randMeanZ = mean( randZVals );
% elseif 0
%     bfMask       = double( maskLab==1 ) ./ sum( maskLab(:)==1 );
%     invalPosMask = isnan(zDiffMap) | mainF;
%     zDiffForConv = zDiffMap;
%     zDiffForConv(isnan(zDiffForConv)) = 0;
%     bfConvZDiff  = conv2( bfMask, zDiffForConv, 'full' );
%     bfConvInValP = conv2( bfMask, invalPosMask );
%     bfConvZDiff( bfConvInValP>0 ) = nan;
%     randMeanZ    = bfConvZDiff
elseif strcmp( prms.shufType, 'shiftRotBF' )
    
    % Initialise.
    randRots  = [0:5:355];
    % Find the BF and crop matrix close
    bfMask    = fieldMaskToTest==1;
    [bfR,bfC] = ind2sub( size(bfMask), find(bfMask) );
    bfCropped = bfMask( min(bfR):max(bfR), min(bfC):max(bfC) );
    randMeanZ = nan( 1, prms.nShuf );
    randOvLp  = nan( 1, prms.nShuf );
    randMeanR  = repmat( {nan(1,prms.nShuf)}, 1, 3 );  % These two have different format, output is a {1,3} or {1,2}
    randCFCorr = repmat( {nan(1,prms.nShuf)}, 1, 2 );  % cell array, where the columns show trials, or trial comparisons.
    
    itSh      = 1;  
    while itSh < prms.nShuf
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find a shuffled position 
        % Apply random rotation, and crop again.
        bfCrRot   = imrotate( bfCropped, randRots(randperm(length(randRots),1)) );
        [bfR,bfC] = ind2sub( size(bfCrRot), find(bfCrRot) );
        bfCrRot   = bfCrRot( min(bfR):max(bfR), min(bfC):max(bfC) );
        % Get the bounding box of the visited env.
        visMask       = ~isnan(M{3});
        [visR,visC]   = ind2sub( size(visMask), find(visMask) );
        % Mark out the square where it is OK to place the rot BF
        valForBFPlace    = false( size(visMask) );
        valForBFPlace( min(visR):(max(visR)-size(bfCrRot,1)), min(visC):(max(visC)-size(bfCrRot,2)) ) = true;
        valForBFPlaceInd = find( valForBFPlace );
        % Choose a random bin in square for TL corner, and make rand pos/rot BF mask.
        [randBFPlaceR,randBFPlaceC] = ind2sub( size(valForBFPlace), valForBFPlaceInd( randperm(length(valForBFPlaceInd),1) ) );
        randBFMask                  = false( size( visMask) );
        randBFMask( (1:size(bfCrRot,1))+randBFPlaceR-1, (1:size(bfCrRot,2))+randBFPlaceC-1 ) = bfCrRot;
        % If this BF placement is >50% overlapping with main field, try again. If not, reassign overlapping
        % bins to the surrounding bins contiguous to the BF (not in main field or outside vis env).
        mfBfOverN      = sum(sum( randBFMask & (mainFMask | ~visMask) ));
        if mfBfOverN > sum( randBFMask(:) )*0.5;   continue;   end
        while mfBfOverN > 1
            nHoodMask = imdilate( randBFMask, [0 1 0; 1 1 1; 0 1 0] )  &  ~randBFMask;
            nHoodInd  = find( nHoodMask & ~mainFMask & visMask );
            if isempty(nHoodInd);   break;   end
            randBFMask( nHoodInd( 1:min( [length(nHoodInd) mfBfOverN] ) ) ) = true;
            randBFMask(  (mainFMask | ~visMask)    )                        = false;
            mfBfOverN = mfBfOverN - length(nHoodInd);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Having finalised the position of the random BF mask, get the diff-Z score for this shuffle.
        if prms.subBslFiring
            randMeanZ(itSh)  = nanmean( normDiffMap( randBFMask ) );
        else
            randMeanZ(itSh)  = nanmean( normMap( randBFMask ) );            
        end
        itSh             = itSh + 1;
        % .. and also get the overlap score, between the actual post-probe field, and the random-shifted barrier field.
        randOvLp(itSh)   = bvcTrTraceOverlapScore( randBFMask, PoPrField, M{3} );  % Note that the rate-weighting of overlap is based on RAW post-probe field, not Z-scored, not post-pre subtracted map.
        
        % Get some scores decribing the similarity of firing across PRE-POST, PRE-CUE trials, *in the shuffled cue-field*.
        % This is a shuffled control for an analysis for NN rev 1: does POST trace firing reflect existing PRE firing.
        % 1) Get the overall mean (z-scored) firing rate in the shuffled barrier field, in all three trials.
        for itTr=1:3
            randMeanR{1,itTr}(itSh)  = nanmean( normM{itTr}( randBFMask ) );
        end
        % Get the bin-by-bin correlation *within cue field* between PRE-CUE and PRE_POST.
        for itTr = 1:2
            corrMask                 = randBFMask & ~isnan( normM{1} ) & ~isnan( normM{itTr+1} );
            randCFCorr{1,itTr}(itSh) = corr2( normM{1}( corrMask ), normM{itTr+1}( corrMask ) );
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save random BF for bebug plot.
%         randBFMask      = double( randBFMask );
%         randBFMask( (mainFMask | ~visMask) ) = 2;
%         rBFAll{itSh}    = randBFMask;
    end
end
if 0
    hFig = gra_multiplot( 10, 10 );   axArr = getappdata(hFig, 'axesHandles');
    for ii=1:100
       imagesc(  rBFAll{ii}, 'parent', axArr(ii) );
    end
end

