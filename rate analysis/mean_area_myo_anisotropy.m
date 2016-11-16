function dataOut=mean_area_myo_anisotropy(dataIn)
% Calculates mean and SD for myosin, area, anisotropy and anistropy XY.
for i=1:size(dataIn,2)
    dataIn(i).meanMyo=nanmean(dataIn(i).myo,2);
    dataIn(i).stdMyo=nanstd(dataIn(i).myo,0,2);
    dataIn(i).meanArea=nanmean(dataIn(i).area,2);
    dataIn(i).stdArea=nanstd(dataIn(i).area,0,2);
    dataIn(i).meanAnisotropyXY=nanmean(dataIn(i).anisotropyXY,2);
    dataIn(i).stdAnisotropyXY=nanstd(dataIn(i).anisotropyXY,0,2);
    dataIn(i).stdAnisotropy=nanstd(dataIn(i).anisotropy,0,2);  
end
dataOut=dataIn;
end