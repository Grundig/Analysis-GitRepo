%% Initialisation
tic;
filepath = 'C:\Users\laptop\Desktop\2018\2018-08-07.hdf5';
run = '120';
data = h5read(filepath, '/RUN 120/coincidences');
coinc = length(data.Pixel);
CoincWindow = 5;                                                           % coincidence window in ns
cw = ceil(CoincWindow / 0.256);                                             % number of samples in coincidence window
RunData = struct();
LowToHiRes = int64(zeros(1,coinc));
texp;
load('texp.mat');                                                           % Loads a matrix containing all possible travel times
%parpool;

% In this loop the high and low resolution timestamps are combined and the
% .time value is the high resolution (0.256 ns intervals) since the start of the run.

parfor n = 1:coinc
    LowToHiRes(n) = int64(data.LowResHitTime(n)-data.LowResHitTime(1))*3.90625e+10;
    RunData(n).time = int64(data.HiResHitTime(n)) + LowToHiRes(n);
    RunData(n).pixel = data.Pixel(n);
end



%UpIndex = find([RunData.Pixel] > 15);
%DownIndex = find([RunData.Pixel] < 16);

% Splitting the hit list to top and bottom detector hits.
UpData = RunData([RunData.pixel] > 15);
DownData = RunData([RunData.pixel] < 16);


%%
uL = length([UpData.pixel]);
timePairs = int64(zeros(uL,2));
pixPairs = single(zeros(uL,2));
P=0;


Closest = nearestpoint([UpData.time],[DownData.time]);
dL = length([DownData.pixel]);
hitIndex = zeros(uL,1);
parfor u = 1:uL

    a = max(Closest(u)-10, 1);
    b = min(Closest(u)+10, dL);
    dt = [DownData(a:b).time];
    ut = UpData(u).time;
    
    Pindex = find((dt>ut-cw) .* (dt<ut+cw)) + a-1;
    
   if ~isempty(Pindex)
        dCheckHit = [DownData(Pindex).pixel]+1; 
        uCheckHit = UpData(u).pixel-15;
        timeCheckHit = Te(uCheckHit,dCheckHit);
        timeRealHit = abs(double([DownData(Pindex).time] - UpData(u).time)*0.256);
        [~,x] = min(abs(timeCheckHit - timeRealHit));
        hitIndex(u) = Pindex(x);
    
        
    end
end

matchIndex = hitIndex~=0;
hitIndex = hitIndex(matchIndex);
timePairs = [[UpData(matchIndex).time]', [DownData(hitIndex).time]'];
pixPairs = [[UpData(matchIndex).pixel]'-15, [DownData(hitIndex).pixel]'+1];

RunTimeFind = toc
%% Confidence Calculation



%pixPairs = pixPairs(matchIndex,:);
%timePairs = timePairs(matchIndex,:);
%treal = abs(double(timePairs(:,1)-timePairs(:,2))*0.256);
% parfor i = [1:length(pixPairs)]
%     try
%     pix1 = pixPairs(i,1);
%     pix2 = pixPairs(i,2);
%     texpect = Te(pix1,pix2);
%     treal = abs(double(timePairs(i,1)-timePairs(i,2))*0.256);
%     td(i) = treal-texpect;
%     catch
%         td(i) = NaN;
%     end
% end
treal = (double(timePairs(:,1)-timePairs(:,2))*0.256);
texpected = Te(sub2ind(size(Te),pixPairs(:,1),pixPairs(:,2)));
td = abs(treal)-texpected;
histogram(td)
sigma = std(td,'omitnan')
plot(treal,texpected,'.')
RunTimeTot = toc