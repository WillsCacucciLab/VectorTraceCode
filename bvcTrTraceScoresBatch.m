function [RM_out] = bvcTrTraceScoresBatch( RM, varargin )
% Batch run of scores across all recorded cells.

% Analysis parameters.
% Trace field definition.
prms.bstProbeDef   = 'BFSum';
prms.mapNormMode   =  'Z';     % 'norm2Pk';   %  
prms.BslFieldThr   = 0.5;      % 'Pk/2';
prms.BarrFieldThr  = 1;        % 'Pk/2';
prms.PostFieldThr  = 1;
prms.subBslFiring  = 1;
prms.BFInvalidIfBslFOverlap = 0;        % If 1, barrF regions which overlap with main field are excluded outright. If 0, main field pixels are still removed from barr F regions, but then these can still count as barr fields.
prms.BFSortParam            = 'Area';   %  'MeanIntensity', 'Area'  % Select 'the' barrier field from many by size or summed rate? (Param names reflect regionprops arguments). 
% Shuffling
prms.forceRunShuf = 1;     % If 1, will run shuffling even if shuffled results exist already in table 'RM'.
prms.nShuf        = 500;
prms.shufType     = 'shiftRotBF';   % 'randBins';  %
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
RM.Properties.UserData.Trace = prms;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-assign score fields.
[RM.TrFMean,RM.TrFSum,RM.BFMean,RM.BFSum,RM.BFPoOverlap] = deal( nan(height(RM),3) );      % Mean and sum firing rates in Trace or Barrier fields. 
[RM.BFMask, RM.PoPrFMask, RM.BslFMask]    = deal(cell( height(RM), 3 ));    % Masks of barrier field, and post-probe field. Note that the latter is not necessarily the 'trace' (which is generally defined using BFMask).
RM.BstPrNum                               = nan( height(RM), 1 );           % 'Best' probe when multiple are present: keep this as both 
RM.BstPrLog                               = false( height(RM), 1 );         % logical and numerical array for indexing convenience.
% Trace scores for Ext R2B trials come in two flavours, using either the BF or PoPrF Masks to define trace area.
maskNames = {'BF','PoPrF'};
for itMT = 1:2
    RM.( ['ER2BTrF_mask' maskNames{itMT}] )    = nan(height(RM),size( RM.RMExtR2B, 2 ));
    RM.( ['ER2BTrFPro_mask' maskNames{itMT}] ) = nan(height(RM),size( RM.RMExtR2B, 2 ));
    RM.( ['ER2BOvLp_mask' maskNames{itMT}] )   = nan(height(RM),size( RM.RMExtR2B, 2 ));
end


for itCl = 1:height(RM)
    
%     if itCl==140
%         dbstop in bvcTrTraceScore_v1 at 32
%     end
    
    % For each probe ..
    for itPr = 1:3
        
        PrMaps = RM.(['RMProbe' num2str(itPr)])( itCl, : );
        
        % Find 'new' firing fields (compared to pre-probe), get mean and summed firing in those fields, for both barrier and post-probe trials.
        [RM.BFMean(itCl,itPr),RM.BFSum(itCl,itPr),RM.TrFMean(itCl,itPr),RM.TrFSum(itCl,itPr),RM.BFMask{itCl,itPr},RM.PoPrFMask{itCl,itPr},RM.BslFMask{itCl,itPr}] = ...
            bvcTrTraceScore_v1( PrMaps, prms );
        
        % Barrier field overlap with post-probe field (i.e. are non-bsl fields in post-probe actually trace?)
        RM.BFPoOverlap(itCl,itPr) = bvcTrTraceOverlapScore( RM.BFMask{itCl,itPr}, RM.PoPrFMask{itCl,itPr}, PrMaps{3} );
        
        
        % Extended return to baseline sequences. Use the same functions, now replacing PoPr with ExtR2B trial N.
        % We also calculate two versions of both scores: 1) using the barr trial to define the trace mask, 2) using the PoPr trial instead.
        % Note that we are including a catch that this block only runs for the probe set which the ER2B follows.
        if ~isempty( RM.RMExtR2B{itCl,1} )  &&   RM.ExtR2BPPr(itCl)==itPr
            for itEB=1:size( RM.RMExtR2B, 2 )
                if isempty( RM.RMExtR2B{itCl,itEB} );   continue;   end
                for itMT = 1:2   % MT = 'Mask Trial', i.e. is mask quantifying ER2B coming from BF or PoPrF.

                    [~,~,RM.( ['ER2BTrF_mask' maskNames{itMT}] )(itCl,itEB),~,~,EBMask,~] = bvcTrTraceScore_v1( { PrMaps{1}, PrMaps{itMT+1}, RM.RMExtR2B{itCl,itEB} }, prms );
                    RM.( ['ER2BOvLp_mask' maskNames{itMT}] )(itCl,itEB)                   = bvcTrTraceOverlapScore( RM.( [maskNames{itMT} 'Mask'] ){itCl,itPr}, EBMask, RM.RMExtR2B{itCl,itEB} );

                end
            end
        end
                
    end
    
    % Define 'best' probe.
    if all( isnan( RM.( prms.bstProbeDef )(itCl,:) ) )
        % This condition is necessary to catch cells without any probes run. Define best probe as 1: 
        % needs to be something or breaks linear indexing across whole array, I think.
        RM.BstPrNum(itCl)     = 1;
    else
        % This is the normal case - get the highest score across all probes run.
        [~,RM.BstPrNum(itCl)] = nanmax(  RM.( prms.bstProbeDef )(itCl,:) );
    end
    RM.BstPrLog(itCl, RM.BstPrNum(itCl)) = true;
    
end
% Make 'proportion' trace scores = trace field as a proportion of barrier field. 
RM.TrFMeanPro = RM.TrFMean ./ RM.BFMean;
RM.TrFSumPro  = RM.TrFSum ./ RM.BFSum;
% Also for extended return to baselines:
erbTrBF                    = RM.BFMean(:,1);
erbTrBF( RM.ExtR2BPPr==2 ) = RM.BFMean(RM.ExtR2BPPr==2,2);
for itMT = 1:2
    RM.( ['ER2BTrFPro_mask' maskNames{itMT}] )  = bsxfun( @rdivide, RM.( ['ER2BTrF_mask' maskNames{itMT}] ), erbTrBF );
end
% Another 'Best' probe index - linear index that get get correct entries from nCellx3 arrays in one line.
RM.BstPrLin   = sub2ind( [height(RM), 3], (1:height(RM))', RM.BstPrNum );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin-shuffled control trace scores.
if prms.forceRunShuf  ||  ~any(strcmp(RM.Properties.VariableNames,'shufTrSc'))
    
    [RM.shuf95, RM.shuf99]                   = deal( nan( height(RM), 1 ) );
    [RM.shufTrSc,RM.shufTrScPro,RM.shufOvLp] = deal(cell( height(RM), 1 ));
    RM.shufER2B                              = cell( height(RM), 4 );
    RM.shufCF_MR                             = cell( height(RM), 3 );
    RM.shufCF_Corr                           = cell( height(RM), 2 );
  
    hWait = waitbar(0, 'Running Shuffle');        
    for itCl = 1:height(RM)
        if ~any( RM.BFMask{itCl,RM.BstPrNum(itCl)}(:)==1 );   continue;    end
        % Run main shuffle - spatial shuffle of BF position in POST trial.
        [RM.shufTrSc{itCl}, RM.shufOvLp{itCl}, RM.shufCF_MR(itCl,:), RM.shufCF_Corr(itCl,:)] = ...
            bvcTrRandTrScores( RM.(['RMProbe' num2str(RM.BstPrNum(itCl))])( itCl, : ), RM.BFMask{itCl,RM.BstPrNum(itCl)},  RM.PoPrFMask{itCl,RM.BstPrNum(itCl)}, prms );
        % Run shuffle for extended return to baseline - substitute ER2B trial for normal POST.
        for itTr = 1:4
            if isempty( RM.RMExtR2B{itCl,itTr} );  continue;   end
            RM.shufER2B{itCl,itTr}   = bvcTrRandTrScores( cat(2, RM.(['RMProbe' num2str(RM.BstPrNum(itCl))])( itCl, 1:2 ),  RM.RMExtR2B(itCl,itTr) ), ...
                                           RM.BFMask{itCl,RM.BstPrNum(itCl)},  RM.PoPrFMask{itCl,RM.BstPrNum(itCl)}, prms );
            RM.shufER2B{itCl,itTr}   = RM.shufER2B{itCl,itTr} ./ RM.BFMean( RM.BstPrNum(itCl) );
        end
        RM.shuf95(itCl)              = prctile( RM.shufTrSc{itCl}, 95 );
        RM.shuf99(itCl)              = prctile( RM.shufTrSc{itCl}, 99 );
        RM.shufTrScPro{itCl}         = RM.shufTrSc{itCl} ./ RM.BFMean( itCl, RM.BstPrNum(itCl) );
        waitbar( itCl/height(RM), hWait );
    end
    delete(hWait);
end


RM_out = RM;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
