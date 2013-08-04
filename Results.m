clear all
load('SurvillanceVideo1_RawData_V01')
k = k-1;
MSE = zeros(k,1);
NonZeroMeas = zeros(k,1);
spatMeas = zeros(k,1);
temporalMeas = 0;
temporalComp = zeros(k,1);
for i = 1:k
    NonZeroMeas(i,1) = sum(sum(MeasMatrix(i).cdata));
    spatMeas(i,1) = .35 * NonZeroMeas(i,1);
    temporalComp(i,1) = NonZeroMeas(i,1)./4096; 
    temporalMeas = NonZeroMeas(i,1) + temporalMeas;
    MSE(i,1) = mean(mean((255.*movdstore(i+1).cdata - 255.*ReconsIm(i).cdata).^2));
end

temporalMeas = temporalMeas./k;
spatMeasAvg = sum(spatMeas)./k;
TempCompAvg = sum(temporalMeas)./k;
MSEAvg = sum(MSE)./k;

MaxIm = max(max(ReconsIm(1).cdata));
MSE = mean(mean((255.*movdr(1).cdata - 255.*ReconsIm(1).cdata).^2));
PSNRim = 10*log10((255*MaxIm).^2/MSE);
