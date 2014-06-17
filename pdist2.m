% This function belongs to Piotr Dollar's Toolbox
% http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html
% Please refer to the above web page for definitions and clarifications
%
% Calculates the distance between sets of vectors.
%
% Let X be an m-by-p matrix representing m points in p-dimensional space
% and Y be an n-by-p matrix representing another set of points in the same
% space. This function computes the m-by-n distance matrix D where D(i,j)
% is the distance between X(i,:) and Y(j,:).  This function has been
% optimized where possible, with most of the distance computations
% requiring few or no loops.
%
% The metric can be one of the following:
%
% 'euclidean' / 'sqeuclidean':
%   Euclidean / SQUARED Euclidean distance.  Note that 'sqeuclidean'
%   is significantly faster.
%
% 'chisq'
%   The chi-squared distance between two vectors is defined as:
%    d(x,y) = sum( (xi-yi)^2 / (xi+yi) ) / 2;
%   The chi-squared distance is useful when comparing histograms.
%
% 'cosine'
%   Distance is defined as the cosine of the angle between two vectors.
%
% 'emd'
%   Earth Mover's Distance (EMD) between positive vectors (histograms).
%   Note for 1D, with all histograms having equal weight, there is a simple
%   closed form for the calculation of the EMD.  The EMD between histograms
%   x and y is given by the sum(abs(cdf(x)-cdf(y))), where cdf is the
%   cumulative distribution function (computed simply by cumsum).
%
% 'L1'
%   The L1 distance between two vectors is defined as:  sum(abs(x-y));
%
%
% USAGE
%  D = pdist2( X, Y, [metric] )
%
% INPUTS
%  X        - [m x p] matrix of m p-dimensional vectors
%  Y        - [n x p] matrix of n p-dimensional vectors
%  metric   - ['sqeuclidean'], 'chisq', 'cosine', 'emd', 'euclidean', 'L1'
%
% OUTPUTS
%  D        - [m x n] distance matrix
%
% EXAMPLE
%  [X,IDX] = demoGenData(100,0,5,4,10,2,0);
%  D = pdist2( X, X, 'sqeuclidean' );
%  distMatrixShow( D, IDX );
%
% See also PDIST, DISTMATRIXSHOW

% Piotr's Image&Video Toolbox      Version 2.0
% Copyright (C) 2007 Piotr Dollar.  [pdollar-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

function D = pdist2( X, Y, metric )
    
    if( nargin<3 || isempty(metric) ); metric=0; end;
    
    switch metric
        case {0,'sqeuclidean'}
            D = distEucSq( X, Y );
        case 'euclidean'
            D = sqrt(distEucSq( X, Y ));
        case 'L1'
            D = distL1( X, Y );
        case 'cosine'
            D = distCosine( X, Y );
        case 'emd'
            D = distEmd( X, Y );
        case 'chisq'
            D = distChiSq( X, Y );
        case 'correlation'
            D = distCorr( X, Y );
        otherwise
            error(['pdist2 - unknown metric: ' metric]);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distL1( X, Y )
    
    m = size(X,1);  n = size(Y,1);
    mOnes = ones(1,m); D = zeros(m,n);
    for i=1:n
        yi = Y(i,:);  yi = yi( mOnes, : );
        D(:,i) = sum( abs( X-yi),2 );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distCorr( X, Y )
    D = 1-corr(X,Y);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distCosine( X, Y )
    
    if( ~isa(X,'double') || ~isa(Y,'double'))
        error( 'Inputs must be of type double'); end;
    
    p=size(X,2);
    XX = sqrt(sum(X.*X,2)); X = X ./ XX(:,ones(1,p));
    YY = sqrt(sum(Y.*Y,2)); Y = Y ./ YY(:,ones(1,p));
    D = 1 - X*Y';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distEmd( X, Y )
    
    Xcdf = cumsum(X,2);
    Ycdf = cumsum(Y,2);
    
    m = size(X,1);  n = size(Y,1);
    mOnes = ones(1,m); D = zeros(m,n);
    for i=1:n
        ycdf = Ycdf(i,:);
        ycdfRep = ycdf( mOnes, : );
        D(:,i) = sum(abs(Xcdf - ycdfRep),2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distChiSq( X, Y )
    
    m = size(X,1);  n = size(Y,1);
    mOnes = ones(1,m); D = zeros(m,n);
    for i=1:n
        yi = Y(i,:);  yiRep = yi( mOnes, : );
        s = yiRep + X;    d = yiRep - X;
        D(:,i) = sum( d.^2 ./ (s+eps), 2 );
    end
    D = D/2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = distEucSq( X, Y )
    
    m = size(X,1); n = size(Y,1);
    XX = sum(X.*X,2);
    YY = sum(Y'.*Y',1);
    D = XX(:,ones(1,n)) + YY(ones(1,m),:) - 2*X*Y';