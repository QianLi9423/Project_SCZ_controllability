function cohiMatrix = wtcMatrix(XMat,TR,bandPass,varargin)
%% Wavelet coherence
%
% USAGE: [Rsq,period,scale,coi,sig95]=wtcMatrix(X[,settings])
%        X: nxT matrix, where n is the number of regions and T is the
%        number of time points.
%
% Settings: Pad: pad the time series with zeros?
% .         Dj: Octaves per scale (default: '1/12')
% .         S0: Minimum scale
% .         J1: Total number of scales
% .         Mother: Mother wavelet (default 'morlet')
% .         MaxScale: An easier way of specifying J1
% .         MakeFigure: Make a figure or simply return the output.
% .         BlackandWhite: Create black and white figures
% .         AR1: the ar1 coefficients of the series
% .              (default='auto' using a naive ar1 estimator. See ar1nv.m)
% .         MonteCarloCount: Number of surrogate data sets in the significance calculation. (default=300)
% .         ArrowDensity (default: [30 30])
% .         ArrowSize (default: 1)
% .         ArrowHeadSize (default: 1)
%
% Settings can also be specified using abbreviations. e.g. ms=MaxScale.
% For detailed help on some parameters type help wavelet.
%
% Example:
%    X =randn(100,50);
%    wtc(sin(t),sin(t.*cos(t*.01)),'ms',16)
%
% Please acknowledge the use of this software in any publications:
%   "Crosswavelet and wavelet coherence software were provided by
%   A. Grinsted."
%
% (C) Aslak Grinsted 2002-2014
%
% http://www.glaciology.net/wavelet-coherence


% -------------------------------------------------------------------------
%The MIT License (MIT)
%
%Copyright (c) 2014 Aslak Grinsted
%
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%
%The above copyright notice and this permission notice shall be included in
%all copies or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
%THE SOFTWARE.
%---------------------------------------------------------------------------



% ------validate and reformat timeseries.
[nRegion,nTPoint] = size(XMat);
% [x,dt]=formatts(x);
% [y,dty]=formatts(y);
% if dt~=dty
%     error('timestep must be equal between time series')
% end
for i = 1:nRegion
    [x{i},dt] = formatts(XMat(i,:));
end



t=x{1}(:,1); %common time period
if length(t)<4
    error('The two time series must overlap.')
end



n=length(t);

%----------default arguments for the wavelet transform-----------
Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
            'Dj',1/12, ...    % this will do 12 sub-octaves per octave
            'S0',2*dt,...    % this says start at a scale of 2 years
            'J1',[],...
            'Mother','Morlet', ...
            'MaxScale',[],...   %a more simple way to specify J1
            'MakeFigure',(nargout==0),...
            'MonteCarloCount',300,...
            'AR1','auto',...
            'ArrowDensity',[30 30],...
            'ArrowSize',1,...
            'ArrowHeadSize',1);
Args=parseArgs(varargin,Args,{'BlackandWhite'});
if isempty(Args.J1)
    if isempty(Args.MaxScale)
        Args.MaxScale=(n*.17)*2*dt; %auto maxscale
    end
    Args.J1=round(log2(Args.MaxScale/Args.S0)/Args.Dj);
end

ad=mean(Args.ArrowDensity);
Args.ArrowSize=Args.ArrowSize*30*.03/ad;
%Args.ArrowHeadSize=Args.ArrowHeadSize*Args.ArrowSize*220;
Args.ArrowHeadSize=Args.ArrowHeadSize*120/ad;

if ~strcmpi(Args.Mother,'morlet')
    warning('WTC:InappropriateSmoothingOperator','Smoothing operator is designed for morlet wavelet.')
end

nx=size(x{1},1);



%-----------:::::::::::::--------- ANALYZE ----------::::::::::::------------

for i = 1:nRegion
    [X{i},period,scale,~] = wavelet(x{i}(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
end
%Smooth X and Y before truncating!  (minimize coi)
sinv=1./(scale');

for i = 1:nRegion
    sX{i}=smoothwavelet(sinv(:,ones(1,nx)).*(abs(X{i}).^2),dt,period,Args.Dj,scale);
end


cohiMatrix = zeros(nRegion,nRegion);
% -------- Cross wavelet -------
for i = 1:(nRegion-1)
    for j = (i+1):nRegion
        Wxy=X{i}.*conj(X{j});
        sWxy=smoothwavelet(sinv(:,ones(1,n)).*Wxy,dt,period,Args.Dj,scale);
        Rsq =abs(sWxy).^2./(sX{i}.*sX{j});
        freq = 1./(period*TR);
                    % timeWave = mean(mean(Rsq(find(freq<bandPass(2)&freq>bandPass(1)),:)));
        timewave = mean(mean(Rsq(freq<bandPass(2)&freq>bandPass(1),:)));
        if(isnan(timewave))
             error('Too few elements, choose larger window.');
        end
        cohiMatrix(i,j) = timewave;
    end
end
cohiMatrix = cohiMatrix + cohiMatrix';
% ----------------------- Wavelet coherence ---------------------------------

