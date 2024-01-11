function [arg_out] = bvcTrScoreVectMaps( arg_in, varargin )
% Score and characterise vector maps.
% Can use this to either score one map, or batch score all maps in table.
% Will also make circularised vector response maps, for easier visualisation.


if istable( arg_in )
    
    
    RM = arg_in;
    
    % Preallocate.
    for itPr = 1:3
        RM.(['VMAng_Pr' num2str(itPr)]) = nan( height(RM), 3 );
        RM.(['VMDist_Pr' num2str(itPr)]) = nan( height(RM), 3 );
        RM.(['VMAngExt_Pr' num2str(itPr)]) = nan( height(RM), 3 );
        RM.(['VMDistExt_Pr' num2str(itPr)]) = nan( height(RM), 3 );
    end
    
    
    for itCl = 1:height(RM)
        for itPr = 1:3
            for itTr = 1:3
                
                % Crop the distance axis of the vector map to half the width of the box:
                vMap      = RM.(['VMProbe' num2str(itPr)]){itCl,itTr}( 2:end, : );  % First row is always only nan.
                
                % Input check - RM.VMProbeN is empty where there is no map, all nan where there is no cue field.
                if isempty( vMap )  ||  all( isnan( vMap(:)) );    
                    continue;    
                end
                
                % Vector dwell map - crop same as rate map.
%                 vMapDwell = RM.(['VDwellProbe' num2str(itPr)]){itCl,itTr}( 2:end, : );
%                 if ~strcmp( RM.(['ObjNewName_Pr' num2str(itPr)]){itCl,1}, 'Barrier100' )
%                     vMap      = vMap( 1:20, : );
%                     vMapDwell = vMapDwell( 1:20, : );
%                 end
                
                % Do scoring and assign.
                RP               = scoreVMapSF( vMap );  % , vMapDwell
                pkSubs           = fliplr( RP.WeightedCentroid );
                pkDist           = (pkSubs(1)*2.5) - 1.25;                        % Subtract 1/2 bin to set distance to bin centre not edge.
                pkAng            = ( (pkSubs(2)*(6*pi/180))-(3*pi/180) ) + pi/2;  % hard-code ang bin size 6 deg. Orientation of circular axis is E-S-W-N-E
                pkAng(pkAng>360) = pkAng(pkAng>360) - 360;
                
                RM.(['VMAng_Pr' num2str(itPr)])(itCl,itTr)   = pkAng;
                RM.(['VMDist_Pr' num2str(itPr)])(itCl,itTr)  = pkDist;
                
                RM.(['VMAngExt_Pr' num2str(itPr)])(itCl,itTr)   = RP.BoundingBox(1)*(6*pi/180);
                RM.(['VMDistExt_Pr' num2str(itPr)])(itCl,itTr)  = RP.BoundingBox(2)*2.5;
            end
        end
    end
 
    % Output
    arg_out = RM;

    
elseif ismatrix( arg_in )
    
    arg_out = scoreVMapSF( arg_in, varargin );
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RP] = scoreVMapSF( vMap, varargin )
% Get correct distance range

% Exclude vector bins by dwell.
if ~isempty( varargin ) && ~isempty( varargin{1} )
    vMapDwell          = varargin{1};
    vMap(vMapDwell<0) = nan;
end


% Shift and wrap the map circularly so that the peak bin is in the centre (deals with issue of field wrapping around map ends).
[pk, pkInd] = nanmax(vMap(:));
[pkR, pkC]  = ind2sub(size(vMap), pkInd);
wrapBin     = pkC+(size(vMap,2)/2);   % wrapBin is the bin that will come to the edge of the wrapped vMap
wrapInd     = mod( (1:size(vMap,2)) + wrapBin, size(vMap,2) );
wrapInd(wrapInd==0) = size(vMap,2);
vMWrap      = vMap(  :,  wrapInd );

% Threshold the map to get main field.
if 1
    mn          = nanmean(vMap(:));
    std         = nanstd(vMap(:));
%     thr         = mn + ((pk - mn)/2);
    thr         = mn + std;
else
    thr         = pk - (pk-nanmin(vMap(:)))*0.2;
end

% Get the main field and associated scores.
RP          = regionprops( vMWrap>=thr, vMWrap,  'Convexhull','MaxIntensity', 'WeightedCentroid','BoundingBox' );
[~,pkFdID]  = max(  cat(1,RP.MaxIntensity)  );
RP          = RP( pkFdID );


% 'Unwrap' any scores based on column indices.
uwFunc                   = @(x,y) mod( x+y, size(vMap,2) );
RP.ConvexHull(:,1)       = feval(uwFunc, RP.ConvexHull(:,1), wrapBin );
RP.WeightedCentroid(:,1) = feval(uwFunc, RP.WeightedCentroid(:,1), wrapBin );



