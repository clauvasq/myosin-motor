function [dataOut]=myo_min_max_norm(dataIn,dataOut)

for i=1:size(dataIn,2)
    tI=dataOut(i).endPts(1);
    tF=dataOut(i).endPts(2);
    tFF=tF-tI+1;
    dataIn(i).myoMin=min(dataIn(i).meanMyo(tI:tF));
    dataIn(i).myoMax=max(dataIn(i).meanMyo(tI:tF));
    dataOut(i).normMyo=(dataIn(i).smooMyo-dataIn(i).myoMin)./(dataIn(i).myoMax-dataIn(i).myoMin);    
end

end