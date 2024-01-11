function [RM, stForAn] = bvcTrExtRtn2Bsl_v2( RM )
% Analysis of extended return to baseline.

% Define trace cells.
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
% Params specific to 
prms.useRealTimeRunOrder = 1;
prms.plotCellLines  = 0;
prms.plotTr3PlusGroup = 1;

nTr2Plot = 3;


% Define trace cells.
RM = bvcTrDefineTrace( RM, prms );

% Get shuffle 95% thr.
sh    = RM.shufTrScPro;
sh    = sh( cellfun(@(x) ~isempty(x), sh ) );
shPop = cell2mat( cellfun( @transpose, sh, 'UniformOutput', 0 ) );


% Select the subset of data where there is a ER2B run on trace cells ..
RM.ER2BExists = cellfun( @(x) ~isempty(x), RM.RMExtR2B(:,1) );
RM.isTrER2B   = RM.ER2BTrFPro_maskBF(:,1) >=RM.Properties.UserData.CurrentThrTrace & RM.ER2BOvLp_maskBF(:,1)>RM.Properties.UserData.CurrentThrOverlap;   % IMPORTANT: defining here cells that are trace or probe preceding ER2B. 
RM.isBarrResp( ~RM.isTrER2B & RM.ER2BTrFPro_maskBF(:,1) >=RM.Properties.UserData.CurrentThrTrace) = false;  % % Exclude non-trace cells with high trace low overlap.


% Calculate the time elapsed between probe end and ER2B start.
RM.TSExtR2B_hr                = cellfun( @str2HrFun, RM.TSExtR2B, 'UniformOutput', 1 );
RM.TSProbe                    = cellfun( @str2HrFun, RM.TSProbe1 );
RM.TSProbe2_hr                = cellfun( @str2HrFun, RM.TSProbe2 );
RM.TSProbe( RM.ExtR2BPPr==2 ) = RM.TSProbe2_hr( RM.ExtR2BPPr==2 );
RM.TimeElapsed                = bsxfun( @minus, RM.TSExtR2B_hr, RM.TSProbe) - 0.3;   %  -0.3 as we are taking the interval from the end of the (20 min) probe trial to the start of the baseline. 

% Get the 'real time' run order (i.e. take into account if the ER2B sequence skips a trial).
if prms.useRealTimeRunOrder
    RM.IndExtR2B_RT = bsxfun( @minus, RM.IndExtR2B,  RM.IndExtR2B(:,1) ) + 1;
else
    RM.IndExtR2B_RT = repmat( 1:4, height(RM), 1 );
    RM.IndExtR2B_RT( isnan( RM.IndExtR2B ) )  = nan;
end

% Plot mean trace.
scNames   = {'TrFPro'}; %  ,'TrFPro','OvLp'
maskNames = {'BF'};  % ,'PoPrF'
grInds    = { RM.isTrER2B & RM.isBarrResp, ~RM.isTrER2B & RM.isBarrResp  };    % 
hFig = gra_multiplot( length(scNames), length(maskNames) );   axArr = getappdata(hFig, 'axesHandles');
lineCol = {'b','r'};
for itSc = 1:length(scNames)
    for itMT = 1:length(maskNames)
        [M,E,N,T,thr95Mn,thr95Med] = deal( nan(length(grInds),nTr2Plot) ); 
        stForAn = nan( 0, 3 );
        for itGr = 1:length(grInds)
            D         = RM.( ['ER2B' scNames{itSc} '_mask' maskNames{itMT}] );
            
            if strcmp( scNames{itSc}, 'OvLp' )
                D( isnan(D)  & ~isnan(RM.ER2BTrFPro_maskBF ) ) = 0;       %   
            end
            
            for itTr = 1:nTr2Plot
                ind_grTr     = bsxfun( @and, RM.IndExtR2B_RT==itTr, grInds{itGr} );
                if itTr==3 && prms.plotTr3PlusGroup
                    ind_grTr     = ind_grTr | bsxfun( @and, RM.IndExtR2B_RT>3, grInds{itGr} );
                end
                D_grTr         = D( ind_grTr);
                M(itGr,itTr)   = nanmean( D_grTr );
                E(itGr,itTr)   = nanstd(D_grTr) ./ sqrt( sum(~isnan(D_grTr)) );
                N(itGr,itTr)   = sum(~isnan(D_grTr));
                T(itGr,itTr)   = nanmean(  RM.TimeElapsed(  ind_grTr )   );
                Tma(itGr,itTr) = nanmax(  RM.TimeElapsed(  ind_grTr )   );
                Tmi(itGr,itTr) = nanmin(  RM.TimeElapsed(  ind_grTr )   );
                
                            
                % 95% population values for trace score
                rInd                = randi( numel(shPop), N(itGr,itTr), 1000 );
                rTr                 = shPop( rInd );
                rMns                = nanmean( rTr, 1 );
                thr95Mn(itGr,itTr)  = prctile( rMns, 95 );
                rMdns               = nanmedian( rTr, 1 );
                thr95Med(itGr,itTr) = prctile( rMdns, 95 );
                
            end

            
            % Record raw values for anova:
            nthTrInd  = RM.IndExtR2B_RT;
            if prms.plotTr3PlusGroup
                nthTrInd( nthTrInd>3 ) = 3;
            end
            isDataInd = ~isnan(D)   &   nthTrInd<=nTr2Plot   &   repmat( grInds{itGr}, 1, size(D,2) );
            stForGr   = [D(isDataInd) nthTrInd(isDataInd) ones(sum(isDataInd(:)),1).*itGr];
            stForAn   = [stForAn; stForGr];
        end
        % Plot
        M  = M';  E = E';  N = N';  thr95Mn = thr95Mn';  thr95Med = thr95Med';
        ax   = axArr( itSc, itMT );
        hB   = bar(ax, M);
        barX = bsxfun( @plus, cell2mat(get(hB,'XData')).',  [hB.XOffset] );
        hold(ax, 'on');
        hL = errorbar( ax, barX, M, E, 'k-' );
        [hL.LineStyle]  = deal( 'none' );
        hB(1).FaceColor = [112, 48, 160] ./ 255;
        hB(2).FaceColor = [255, 192, 0] ./ 255;
%         Plot indivdual cell lines.
        if prms.plotCellLines
            for itGr = 1:length(grInds)
                Y = RM.( ['ER2B' scNames{itSc} '_mask' maskNames{itMT}] )(  grInds{itGr}, : );
                X = RM.IndExtR2B_RT(  grInds{itGr}, : );
                X(X>nTr2Plot) = nan;
                for itCl = 1:size(X,1)
                    plot( X(itCl,:)+hB(itGr).XOffset, Y(itCl,:), ['x' lineCol{itGr} '-'] );
                end
            end
        end
        % Pop mean 95% threshold line. 
        for itBr=1:numel(barX)
            plot(ax, [-0.15 0.15]+barX(itBr), [1 1].*thr95Mn(itBr), 'k:' );
        end
        % Add annotaion:
        if 1
            % Add N
            text(ax, barX(:), M(:)+0.1, cellstr(num2str( N(:) )) );
        end
        % Add mean time delay (in tick labels).
        mnT  = nanmean(T,1);
        maxT = nanmax(Tma,[],1);
        minT = nanmin(Tmi,[],1);
        for itTr=1:nTr2Plot
            xLab{itTr} = sprintf('Bsl %d: %2.1f Hr (%2.1f-%2.1f)', itTr, mnT(itTr), minT(itTr), maxT(itTr) ); %  (ax, X(:), M(:)+0.1, cellstr(num2str( T(:), '%2.1f' )) );
        end
        ax.XTickLabel = xLab;
        if ~verLessThan( 'matlab', '2018a' )
            xtickangle(ax, 45);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Boxplot figure.
        figure;
        boxplot( stForAn(:,1), { stForAn(:,2), stForAn(:,3) }, 'symbol', 'o', 'OutlierSize' , 3 );


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANOVA.
        % Run ANOVA.
        [~,~,STATS,~]=anovan( stForAn(:,1), { stForAn(:,2), stForAn(:,3) }, 'varnames', {'trial','cellType'}, 'model', 'interaction' );
        figure;
         C = multcompare( STATS, 'Dimension', 1:2 );
%          disp( C );
        % t-test trial N trace versus non-trace.
        for itCp=1:nTr2Plot
            [~,P,~,STATS] = ttest2( stForAn( stForAn(:,2)==itCp & stForAn(:,3)==1, 1 ), ...
                                    stForAn( stForAn(:,2)==itCp & stForAn(:,3)==2, 1 )  );
            fprintf(1, 'p=%4.3f\n', P);
        end
%         [~,P,~,STATS] = ttest( stForAn( stForAn(:,2)==1 & stForAn(:,3)==1, : ),  0.11  );

    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB_FUNCTION.
function [t] = str2HrFun(str)
if isempty(str) || all(isnan(str))
    t = nan;
else
    t = str2double(str(1:2)) + (str2double(str(4:5))/60);
end
