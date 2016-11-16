function dataOut=smoo_rate_xcorr(dataIn)
% Smooths data using +/-1 frame, normalizes myosin to area, calcualtes rate
% using rate function and then uses nanxcorr for cross correlation. 
% Note: Uses window of 20 for xcorrelation

for i=1:size(dataIn,2) 
%     smooths data
    dataIn(i).smooArea=smooth2a(dataIn(i).area,1,0);
    dataIn(i).smooMyo=smooth2a(dataIn(i).myo,1,0);
%     normalizes myosin to area
    dataIn(i).normMyo=dataIn(i).smooMyo./dataIn(i).smooArea;
%     calculates constriction rate
    dataIn(i).rateArea=-rate(dataIn(i).smooArea,dataIn(i).timeRes);
    dataIn(i).rateMyo=rate(dataIn(i).smooMyo,dataIn(i).timeRes);
    dataIn(i).rateNormMyo=rate(dataIn(i).normMyo,dataIn(i).timeRes);
%     calculates cross correlation
    dataIn(i).xcorr=nanxcorr(dataIn(i).rateArea(dataIn(i).corrTime,:),...
        dataIn(i).rateNormMyo(dataIn(i).corrTime,:),...
        20);
    dataIn(i).xcorr=dataIn(i).xcorr(0==sum(isnan(dataIn(i).xcorr),2),:);
    dataIn(i).meanxcorr=nanmean(dataIn(i).xcorr);
end
dataOut=dataIn;
end