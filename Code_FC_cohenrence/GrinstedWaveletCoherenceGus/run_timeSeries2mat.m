function Aij = run_timeSeries2mat(timeSeriesInput,inputTR,passBand,windowSizeInput,windowOverlapInput,analysisTypeInput,outputDir,subjectNameInput,taskNameInput,timeLag)
% PURPOSE: Converts time series data to multilayer network
%
% INPUT:
%          timeSeriesInput: Time series data
%                           * Rows | node/voxel/region/sensor
%                           * Cols | sample value
%
%                  inputTR: Repetition time, duration per sample (seconds)
%                           (Default: 1)
%
%          windowSizeInput:	Number of samples per window/block
%                           (Default: Length of time series)
%
%       windowOverlapInput:	Number of samples windows overlap
%                           (Default: 0)
%
%                outputDir: Output directory for adjacency matrix MAT file
%                           (Default: Current directory)
%
%             analysisType: Choice of analysis
%                           (Default: Pearson's correlation)
%                           * Pearson's correlation ('Pearson','corr')
%                           * Spearman's rank correlation ('Spearman')
%                           * Kendall's tau correlation ('Kendall')
%                           * Wavelet coherence ('wavelet','wtc')
%                           * Cross correlation ('crosscorr','xcorr')
%
%         subjectNameInput: Subject name (Default: Leave empty)
%
%            taskNameInput: Task name (Default: Leave empty)
%
%                  timeLag: Time lag for cross-correlation
%
% OUTPUT:
%	Aij:	Multilayer network
%	Network saved to output directory in the form:
%	Aij_<taskName>_t<windowSize>_o<windowOverlap>.mat
%                       - OR -
%	Aij_<taskName>_<analysisType>_t<windowSize>_o<windowOverlap>.mat (If analysisFlag = 1: Line 100)
%
%--------------------------------------------------------------------------
%
%  Author: Qawi Telesford
%    Date: 2015-03-19
% Version: 1.0
%
%--------------------------------------------------------------------------
%% Error Checking
if(nargin < 1 || isempty(timeSeriesInput))
    error('Missing time series data, please input time series data.');
end

if(nargin < 2 || isempty(inputTR))
    disp('Input TR not provided, using default of 1s.'); inputTR = 1;
end

if(nargin < 3 || isempty(passBand))
    disp('Input passBand not provided, using default of [0.01,0.08].');
end

if(nargin < 4 || isempty(windowSizeInput))
    disp(['Window size not provided, using length of time series: ',num2str(size(timeSeriesInput,2))]);
end

if(nargin < 5 || isempty(windowOverlapInput))
    windowOverlapInput = 0;
end

if(nargin < 6 || isempty(analysisTypeInput))
    disp('Analysis type not selected, using default: Pearson''s');
    analysisTypeInput = 'correlation';
end

if(nargin < 7 || isempty(outputDir))
    outputDir = pwd;
    % disp(['No Output directory, using current directory: ',outputDir]);
end

if(strcmp(analysisTypeInput,'corr') || strcmp(analysisTypeInput,'pearson'))
    analysisTypeInput = 'Pearson';
elseif(strcmp(analysisTypeInput,'spearman'))
    analysisTypeInput = 'Spearman';
elseif(strcmp(analysisTypeInput,'kendall'))
    analysisTypeInput = 'Kendall';
elseif(strcmp(analysisTypeInput,'wtc'))
    analysisTypeInput = 'wavelet';
elseif(strcmp(analysisTypeInput,'xcorr'))
    analysisTypeInput = 'crosscorr';
end

if(nargin < 8 || isempty(subjectNameInput))
    subjectNameInput = [];
end

if(nargin < 9 || isempty(taskNameInput))
    taskNameInput = [];
end

if(nargin < 10 || isempty(timeLag))
    timeLag = 1;
end

analysisOutputFlag = 0; % Flag to include analysis type output filename

%% Initialize Output Variables
Aij            = {};%#ok<*NASGU>         % multilayer adjacency matrix
timeSeriesData = timeSeriesInput;	     % time series data
bandPass       = [];                     % band-pass filter frequency range
dataRange      = [];                     % sample range for each window
dataLength     = size(timeSeriesData,2); % length of time series
windowSize     = windowSizeInput;        % window size
windowOverlap  = windowOverlapInput;	 % window overlap
samplingRate   = 1/inputTR;              % sampling rate 1/TR
TR             = inputTR;                % repetition time
N              = size(timeSeriesData,1); % number of nodes
analysisType   = analysisTypeInput;      % analysis type
subjectName    = subjectNameInput;       % subject name
taskName       = taskNameInput;          % task name

% Clear input variables
clear('timeSeriesInput','inputTR',....
    'windowSizeInput','windowOverlapInput',...
    'analysisTypeInput','subjectNameInput','taskNameInput');

% Band-pass filter for wavelet coherence analysis, assuming data is fMRI
if(strcmp(analysisType,'wavelet'))
    % tasks
    bandPass =  passBand; % [0.06, 0.125];
    
    % resting
    % bandPass = [0.01 0.08];
end

%% Data range matrix, specifying index for time series windows
stepSize = round(windowSize*(1-(windowOverlap/windowSize)));

if(windowSize ~= dataLength)
    timeWinStart = round(1:stepSize:dataLength);
    timeWinEnd   = round(timeWinStart + windowSize - 1);
else
    timeWinStart = 1;
    timeWinEnd   = dataLength;
end

% Truncate time windows larger than time series
timeWinTrunc = find(timeWinEnd>dataLength);
timeWinStart(timeWinTrunc) = [];
timeWinEnd(timeWinTrunc) = [];

dataRange = [timeWinStart; timeWinEnd]'; %#ok<NASGU>

%% Generate multilayer network
if(~strcmp(analysisType,'crosscorr'))
    Aij = cell(length(timeWinStart),1);
else
    clear('Aij');
    Aij = struct();
end

for kk = 1:length(timeWinStart)
    
    d = timeSeriesData(:,timeWinStart(kk):timeWinEnd(kk))'; % Each column is a channel.
    mca = bsxfun(@minus,d,mean(d));
    
    if(strcmp(analysisType,'correlation') || strcmp(analysisType,'corr') || strcmp(analysisType,'Pearson') || strcmp(analysisType,'pearson'))
        % Pearson's correlation
        Aij{kk} = corr(mca,mca,'type','Pearson'); %#ok<*AGROW>
        Aij{kk}(logical(eye(N))) = 0;
        
    elseif(strcmp(analysisType,'Spearman') || strcmp(analysisType,'spearman'))
        % Spearman's rank correlation
        Aij{kk} = corr(mca,mca,'type','Spearman');
        Aij{kk}(logical(eye(N))) = 0;
        
    elseif(strcmp(analysisType,'Kendall') || strcmp(analysisType,'kendall'))
        % Spearman's rank correlation
        Aij{kk} = corr(mca,mca,'type','Kendall');
        Aij{kk}(logical(eye(N))) = 0;
        
    elseif(strcmp(analysisType,'crosscorr') || strcmp(analysisType,'xcorr'))
        % Cross correlation
        [Aij_crossCorr,Aij_crossCorrSign,...
            Aij_crossCorrPos,Aij_crossCorrNeg,...
            Aij_noZeroLag,...
            Aij_crossCorr_noZeroLag,Aij_crossCorrSign_noZeroLag,...
            Aij_crossCorrPos_noZeroLag,Aij_crossCorrNeg_noZeroLag] = ...
            deal(cell(numel(timeWinStart),1));
        
        Aij_noZeroLag{kk,1} = ones(N);
        
        crossCorrSignal = xcorr(mca,timeLag,'coef');
        [maxCC, maxCC_index] = max(abs(crossCorrSignal));
        Aij_crossCorr{kk,1} = reshape(maxCC,N,N);
        
        % cross corr sign is the signed value of max abs
        ccs_index = sub2ind(size(crossCorrSignal),maxCC_index,1:size(crossCorrSignal,2));
        Aij_crossCorrSign{kk,1} = reshape(crossCorrSignal(ccs_index),N,N);
        
        % positive and negative adjacency matrices
        Aij_crossCorrPos{kk,1} = Aij_crossCorrSign{kk,1}.*double(Aij_crossCorrSign{kk,1}>0);
        Aij_crossCorrNeg{kk,1} = Aij_crossCorrSign{kk,1}.*double(Aij_crossCorrSign{kk,1}<0);
        
        % Aij_noZeroLag is 0 if the max crossCorrSignal is at timelag 0 or is on
        % the diagonal
        is_zero = maxCC_index==timeLag+1;
        Aij_noZeroLag{kk,1}(is_zero) = 0;
        Aij_noZeroLag{kk,1}(logical(eye(N))) = 0;
        Aij_crossCorr{kk,1}(logical(eye(N))) = 0;
        Aij_crossCorrSign{kk,1}(logical(eye(N))) = 0;
        Aij_crossCorrPos{kk,1}(logical(eye(N))) = 0;
        Aij_crossCorrNeg{kk,1}(logical(eye(N))) = 0;
        
        Aij_crossCorr_noZeroLag{kk,1}     = Aij_noZeroLag{kk,1}.*Aij_crossCorr{kk,1};
        Aij_crossCorrSign_noZeroLag{kk,1}	= Aij_noZeroLag{kk,1}.*Aij_crossCorrSign{kk,1};
        Aij_crossCorrPos_noZeroLag{kk,1}	= Aij_noZeroLag{kk,1}.*double(Aij_crossCorrPos{kk,1}>0);
        Aij_crossCorrNeg_noZeroLag{kk,1}	= Aij_noZeroLag{kk,1}.*double(Aij_crossCorrNeg{kk,1}<0);
        
        Aij.crossCorr{kk,1}                 = Aij_crossCorr{kk,1};
        Aij.crossCorrSign{kk,1}             = Aij_crossCorrSign{kk,1};
        Aij.crossCorrPos{kk,1}              = Aij_crossCorrPos{kk,1};
        Aij.crossCorrNeg{kk,1}              = Aij_crossCorrNeg{kk,1};
        Aij.noZeroLag{kk,1}                 = Aij_noZeroLag{kk,1};
        Aij.crossCorr_noZeroLag{kk,1}       = Aij_crossCorr_noZeroLag{kk,1};
        Aij.crossCorrSign_noZeroLag{kk,1}	= Aij_crossCorrSign_noZeroLag{kk,1};
        Aij.crossCorrPos_noZeroLag{kk,1}	= Aij_crossCorrPos_noZeroLag{kk,1};
        Aij.crossCorrNeg_noZeroLag{kk,1}	= Aij_crossCorrNeg_noZeroLag{kk,1};
        
    elseif(strcmp(analysisType,'wavelet') || strcmp(analysisType,'wtc'))
        Aij{kk} = zeros(N);
        for ii = 1:N
            for jj = 1:N
                timeData1 = timeSeriesData(ii,timeWinStart(kk):timeWinEnd(kk));
                timeData2 = timeSeriesData(jj,timeWinStart(kk):timeWinEnd(kk));
                
                if(ii ~= jj && ii < jj)
                    % Wavelet coherence
                    if any(isnan(timeData1)) || any(isnan(timeData2)) || std(timeData1) == 0 || std(timeData2) == 0
                        Aij{kk,1}(ii,jj)=nan;
                        continue;
                    end
                    [Rsq,period,~] = wtc(timeData1,timeData2,'mcc',0);
                    freq = 1./(period*TR);
                    % timeWave = mean(mean(Rsq(find(freq<bandPass(2)&freq>bandPass(1)),:)));
                    timeWave = mean(mean(Rsq(freq<bandPass(2)&freq>bandPass(1),:)));
                    
                    if(isnan(timeWave))
                        error('Too few elements, choose larger window.');
                    end
                    Aij{kk,1}(ii,jj) = timeWave;
                end

            end
        end
        Aij{kk} = Aij{kk} + Aij{kk}';
    end
end