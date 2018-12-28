%% An Initialisation
clear all;
filedate = lower(datestr(floor(now-1)));
global HEIGHT PIXELSIZE DETECTOR_NORTH;

% Detector design:
DETECTOR_NORTH = 141; % used to correct for position of detector against real North 
                       %TODO verify real north direction
HEIGHT = 0.7;
PIXELSIZE = 0.25;

try
    filepath = ['D:\Data\2018\' filedate '.hdf5'];
    info = h5info(filepath);
    groups = info.Groups;
    lG = length(groups);   
catch
    disp('File not found!')
    quit
    return
end
%%
for v = 1:lG
    runname = string(groups(v).Name);
    data = h5read(filepath,[char(runname) '/coincidences']);
    coinc = length(data.Pixel);
    CoincWindow = 5;                                                       % coincidence window in ns
    cw = ceil(CoincWindow / 0.25);                                         % number of samples in coincidence window
    RunData = struct();
    LowToHiRes = int64(zeros(1,coinc));
    Te = texp;                                                             % Creates a matrix containing all possible travel times
    
    % In this loop the high and low resolution timestamps are combined and the
    % .time value is the high resolution (0.25 ns intervals) since the start of the run.
    
    parfor n = 1:coinc
        LowToHiRes(n) = int64(data.LowResHitTime(n)-data.LowResHitTime(1))*4e+9;
        RunData(n).time = int64(data.HiResHitTime(n)) + LowToHiRes(n);
        RunData(n).pixel = data.Pixel(n);
        RunData(n).realtimeLow = int32(data.LowResHitTime(n));
        RunData(n).realtimeHi = int64(data.HiResHitTime(n));
    end
    
    
    % Splitting the hit list to top and bottom detector hits.
    UpData = RunData([RunData.pixel] > 15);
    DownData = RunData([RunData.pixel] < 16);
    
    
    %%
    uL = length([UpData.pixel]);
    timePairs = int64(zeros(uL,2));
    pixPairs = single(zeros(uL,2));
    realtimeLowPairs = uint32(zeros(uL,2));
    realtimeHiPairs = uint64(zeros(uL,2));
    P=0;
    
    
    Closest = nearestpoint([UpData.time],[DownData.time]);
    dL = length([DownData.pixel]);
    hitIndex = zeros(uL,1);
    parfor u = 1:uL
        
        a = max(Closest(u)-10, 1);
        b = min(Closest(u)+10, dL);
        dt = [DownData(a:b).time];
        ut = UpData(u).time;
        
        Pindex = find((dt>ut-4) .* (dt<ut+cw)) + a-1;
        
        if ~isempty(Pindex)
            dCheckHit = [DownData(Pindex).pixel]+1;
            uCheckHit = UpData(u).pixel-15;
            timeCheckHit = Te(uCheckHit,dCheckHit);
            timeRealHit = abs(double([DownData(Pindex).time] - UpData(u).time)*0.25);
            [~,x] = min(abs(timeCheckHit - timeRealHit));
            hitIndex(u) = Pindex(x);
        end
    end
    
    matchIndex = hitIndex~=0;
    hitIndex = hitIndex(matchIndex);
    timePairs = [[UpData(matchIndex).time]', [DownData(hitIndex).time]'];
    pixPairs = [[UpData(matchIndex).pixel]'-15, [DownData(hitIndex).pixel]'+1];
    realtimeLowPairs = [[UpData(matchIndex).realtimeLow]',[DownData(hitIndex).realtimeLow]'];
    realtimeHiPairs = [[UpData(matchIndex).realtimeHi]',[DownData(hitIndex).realtimeHi]'];
    
    %% Confidence Calculation
    
    treal = (double(timePairs(:,1)-timePairs(:,2))*0.25);
    texpected = Te(sub2ind(size(Te),pixPairs(:,1),pixPairs(:,2)));
    td = abs(treal)-texpected;
    
    
    %histogram(td)
    sigma = std(td,'omitnan');
    
    %% Direction Detection
    hitDirection = getHitDirection; % Load hit direction and inclination map
    
    %% Muon Candidates
    lP = length(realtimeHiPairs);
    muonCand = struct();
    muonCand.timeLow = uint32(zeros(lP,1));
    muonCand.timeHi = uint64(zeros(lP,1));
    muonCand.confidence = double(zeros(lP,1));
    muonCand.direction = single(zeros(lP,1));
    muonCand.inclination = single(zeros(lP,1));
    
    parfor i = 1:lP
        if realtimeLowPairs(i,1)== realtimeLowPairs(i,2)
            muonCand(i).timeLow = realtimeLowPairs(i,1);
            muonCand(i).timeHi = mean([realtimeHiPairs(i,1),realtimeHiPairs(i,2)],'native');
        else
            muonCand(i).timeLow = min(realtimeLowPairs(i,1),realtimeLowPairs(i,2));
            muonCand(i).timeHi = mean([realtimeHiPairs(i,1),realtimeHiPairs(i,2)+2e+9],'native');
        end
        z = td(i)/sigma;
        muonCand(i).confidence = 1-normcdf(z)
        
        muonCand(i).direction = hitDirection.direction(sub2ind(size(hitDirection.direction),pixPairs(i,1),pixPairs(i,2)));
        muonCand(i).inclination = hitDirection.inclination(sub2ind(size(hitDirection.inclination),pixPairs(i,1),pixPairs(i,2)));
    end
    %% Hit Rate
    hitRate = diff(find([1, diff([muonCand.timeLow]), 1]));
    hitLabels = unique([muonCand.timeLow]);
    
    % Fill no-hit times
    hitLabelsMissing = setdiff([min(hitLabels):max(hitLabels)], hitLabels);
    hitLabels = [hitLabels hitLabelsMissing];
    hitRate = [hitRate zeros(1,length(hitLabelsMissing))];
    
    % Timewise sorting 
    [hitLabels, sortIndex] = sort(hitLabels);
    hitRate = hitRate(sortIndex);
    
    HitRate = struct('time', num2cell(uint32(hitLabels)), 'rate', num2cell(uint16(hitRate)));
    
    
    save(['D:\Data\EEE Analysis\' 'BERG-01-' datestr((now-1),'YYYY-mm-dd') '-RUN' num2str(v) '_hitrates.bin'],'HitRate')
    save(['D:\Data\EEE Analysis\' 'BERG-01-' datestr((now-1),'YYYY-mm-dd') '-RUN' num2str(v) '.bin'],'muonCand')
end
%% End Cleanup
clear all;
close all;
