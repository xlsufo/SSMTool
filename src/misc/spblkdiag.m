function S= spblkdiag(varargin)
%SPBLKDIAG make sparse block diagonal matrix
%    S = spblkdiag(X1,X2,X3,...,Xn);
%    S = spblkdiag(X3d);
%    will make a block diagonal matrix out of 2D slices 

%  Copyright 2000-2020 The MathWorks, Inc. and Ford Global Technologies, Inc.

if nargin == 1 && ndims(varargin{1})==3
    % convert 3D matrix to a list of 2D matrices
    varargin= num2cell(varargin{1},1:2);
end
%call internal function
S= matlab.internal.math.blkdiag(varargin{:});