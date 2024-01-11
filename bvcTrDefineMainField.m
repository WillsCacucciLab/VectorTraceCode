function [out] = bvcTrDefineMainField( map, varargin )
% Define main field in pre-barrier trial: this is a couple of lines, but needs 
% to be consistent across several different calling functions.
%
%   [mask] = bvcTrDefineMainField( map, prms )
%
% map is the bsl trial rate map.
% mask is a logical mask for the main field.

prms.BslFieldThr  = 0.5;
prms.mapNormMode  = 'Z';
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

if strcmp( prms.mapNormMode, 'Z' )
    map = (map-nanmean(map(:))) ./ nanstd(map(:));
elseif strcmp( prms.mapNormMode, 'norm2Pk' )
    map = map ./ nanmax(map(:));
end

out = map >= prms.BslFieldThr;
