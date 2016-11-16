function dataOut=pulse_id_thrAreaNormMyo(dataIn,dataOut,thr)

for i=1:size(dataIn,2)
    tI=dataOut(i).endPts(1);
    tF=dataOut(i).endPts(2);
    tFF=tF-tI+1;
    
    dataOut(i).areaNormMyo=dataOut(i).normMyo./dataIn(i).smooArea;
    dataOut(i).rateMyo=rate(dataOut(i).areaNormMyo,dataIn(i).timeRes);
    dataOut(i).smooRateMyo=smooth2a(dataOut(i).rateMyo(tI:tF,:),1,0);
    dataOut(i).smooRateMyo(:,any(isnan(dataOut(i).smooRateMyo),1))=NaN;
    pulseDim=[tFF,dataIn(i).numCells,size(thr,2)];
    % Holds pulse sequence,x has value of 1 where rate is greater than thr
    dataOut(i).timeMyo=NaN(pulseDim);
    % Holds location and amp of pulse
    dataOut(i).maxMyo=NaN(pulseDim);
    % Holds number of pulses per cell
    dataOut(i).numMyo=NaN(size(thr,2),dataIn(i).numCells);
    dataOut(i).numMyo(isnan(dataOut(i).smooRateMyo(1:size(thr,2),:))==0)=0;
end

% Using smoothed rate, find where rate is greater than the threshold
for k=1:size(dataIn,2)
    tI=dataOut(k).endPts(1);
    tF=dataOut(k).endPts(2);
    tFF=tF-tI+1;
%     Establish cutoff for pulses using mean and thr
    meanRateMyo=nanmean(dataOut(k).smooRateMyo);
    stdRateMyo=nanstd(dataOut(k).smooRateMyo);
    dataOut(k).cutoffMyo=vertcat(meanRateMyo+thr(1)*stdRateMyo,...
        meanRateMyo+thr(2)*stdRateMyo,...
        meanRateMyo+thr(3)*stdRateMyo,...
        meanRateMyo+thr(4)*stdRateMyo);
    for i=1:size(thr,2)
        repThr=repmat(dataOut(k).cutoffMyo(i,:),tFF,1);
        dataOut(k).timeMyo(1:tFF,:,i)=dataOut(k).smooRateMyo>repThr;
    end
end
%  Get count for number of pulses per cell, also find amplitude and 
% location of pulse
for k=1:size(dataIn,2)
    tI=dataOut(k).endPts(1);
    tF=dataOut(k).endPts(2);
    tFF=tF-tI+1;
    for t=1:4
        for c=1:dataIn(k).numCells
            startPulse=1;
            for i=2:tFF-1
                if dataOut(k).timeMyo(i,c,t)==1 && dataOut(k).timeMyo(i-1,c,t)==0
                    startPulse=i;
                    if isnan(dataOut(k).numMyo(t,c))==1
                        dataOut(k).numMyo(t,c)=0;
                    end
                    dataOut(k).numMyo(t,c)=dataOut(k).numMyo(t,c)+1;
                elseif dataOut(k).timeMyo(i,c,t)==1 && dataOut(k).timeMyo(i+1,c,t)==0
                    endPulse=i;
                    [pulseAmp, pulseIndex]=max(dataOut(k).smooRateMyo(startPulse:endPulse,c));
                    dataOut(k).maxMyo(startPulse+pulseIndex-1,c,t)=pulseAmp;
                end
            end
        end
    end
end

end