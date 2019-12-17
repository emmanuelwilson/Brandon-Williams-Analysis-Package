function [rate] = FiringRateSort(shuf,itterations)
% %{
rate = zeros(length(shuf.mrall),1);

for i = 1 : length(shuf.mrall)/itterations
    rate(i:length(shuf.mrall)/itterations:length(shuf.mrall),1) = sum(shuf.firing(:,i))/length(shuf.firing(:,1))*30; %length(find(shuf.firing(:,i)))/length(shuf.firing(:,1));
end
%}
%{
rate = zeros(length(shuf.deconvolvedSig(1,:)),1);

for i = 1 : length(shuf.deconvolvedSig(1,:))/itterations
    rate(i:length(shuf.deconvolvedSig(1,:))/itterations:length(shuf.deconvolvedSig(1,:)),1) = sum(shuf.deconvolvedSig(:,i))/length(shuf.deconvolvedSig(:,1))*30; 
end
%}
end

%{
MRL2 = cat(1,out21.mrall,out22.mrall,out23.mrall,out24.mrall,out25.mrall,out26.mrall,out28.mrall,out29.mrall,out210.mrall,out211.mrall,out212.mrall,out214.mrall);
MRLW = cat(1,outW1.mrall,outW2.mrall,outW3.mrall,outW4.mrall,outW5.mrall,outW6.mrall,outW7.mrall,outW8.mrall,outW9.mrall,outW14.mrall);
MRLA1 = cat(1,outA11.mrall,outA12.mrall,outA13.mrall,outA14.mrall,outA15.mrall,outA16.mrall,outA17.mrall,outA18.mrall,outA19.mrall,outA110.mrall,outA111.mrall,outA112.mrall,outA113.mrall,outA114.mrall);
MRLA2 = cat(1,outA21.mrall,outA24.mrall,outA25.mrall,outA28.mrall,outA29.mrall,outA210.mrall,outA211.mrall,outA212.mrall,outA213.mrall,outA214.mrall);
MRL6 = cat(1,out61.mrall,out62.mrall,out63.mrall,out64.mrall,out65.mrall,out69.mrall,out610.mrall,out611.mrall,out612.mrall,out613.mrall,out614.mrall);
MRL10 = cat(1,out101.mrall,out102.mrall,out103.mrall,out104.mrall,out106.mrall,out108.mrall,out1012.mrall,out1014.mrall);
MRL8 = cat(1,out805.mrall,out806.mrall,out807.mrall);
MRL = cat(1,MRL2,MRLW,MRLA1,MRL6,MRL10,MRLA2,MRL8);

FrateA2 = cat(1, FiringRateSort(outA21,100),FiringRateSort(outA24,100),FiringRateSort(outA28,100),FiringRateSort(outA29,100),FiringRateSort(outA210,100),FiringRateSort(outA211,100),FiringRateSort(outA212,100),FiringRateSort(outA213,100),FiringRateSort(outA214,100));
Frate6 = cat(1,FiringRateSort(out61,100),FiringRateSort(out62,100),FiringRateSort(out63,100),FiringRateSort(out64,100),FiringRateSort(out65,100),FiringRateSort(out612,100),FiringRateSort(out613,100),FiringRateSort(out614,100));
Frate2 = cat(1,FiringRateSort(out21,100),FiringRateSort(out22,100),FiringRateSort(out23,100),FiringRateSort(out24,100),FiringRateSort(out25,100),FiringRateSort(out26,100),FiringRateSort(out28,100),FiringRateSort(out29,100),FiringRateSort(out210,100),FiringRateSort(out211,100),FiringRateSort(out212,100),FiringRateSort(out214,100));
FrateW = cat(1,FiringRateSort(outW1,100),FiringRateSort(outW2,100),FiringRateSort(outW3,100),FiringRateSort(outW4,100),FiringRateSort(outW5,100),FiringRateSort(outW6,100),FiringRateSort(outW7,100),FiringRateSort(outW8,100),FiringRateSort(outW9,100),FiringRateSort(outW14,100));
FrateA1 = cat(1,FiringRateSort(outA11,100),FiringRateSort(outA12,100),FiringRateSort(outA13,100),FiringRateSort(outA14,100),FiringRateSort(outA15,100),FiringRateSort(outA16,100),FiringRateSort(outA17,100),FiringRateSort(outA18,100),FiringRateSort(outA19,100),FiringRateSort(outA110,100),FiringRateSort(outA111,100),FiringRateSort(outA112,100),FiringRateSort(outA113,100),FiringRateSort(outA114,100));
Frate10 = cat(1, FiringRateSort(out101,100),FiringRateSort(out102,100),FiringRateSort(out103,100),FiringRateSort(out104,100),FiringRateSort(out106,100),FiringRateSort(out108,100),FiringRateSort(out1012,100),FiringRateSort(out1014,100));
Frate8 = cat(1,FiringRateSort(out805,100),FiringRateSort(out806,100),FiringRateSort(out807,100));

Cor2  = cat(1,out21.correlationEO,out22.correlationEO,out23.correlationEO,out24.correlationEO,out25.correlationEO,out26.correlationEO,out28.correlationEO,out29.correlationEO,out210.correlationEO,out211.correlationEO,out212.correlationEO,out214.correlationEO);
CorW = cat(1,outW1.correlationEO,outW2.correlationEO,outW3.correlationEO,outW4.correlationEO,outW5.correlationEO,outW6.correlationEO,outW7.correlationEO,outW8.correlationEO,outW9.correlationEO,outW14.correlationEO);
CorA1 = cat(1,outA11.correlationEO,outA12.correlationEO,outA13.correlationEO,outA14.correlationEO,outA15.correlationEO,outA16.correlationEO,outA17.correlationEO,outA18.correlationEO,outA19.correlationEO,outA110.correlationEO,outA111.correlationEO,outA112.correlationEO,outA113.correlationEO,outA114.correlationEO);
CorA2 = cat(1,outA21.correlationEO,outA24.correlationEO,outA28.correlationEO,outA29.correlationEO,outA210.correlationEO,outA211.correlationEO,outA212.correlationEO,outA213.correlationEO,outA214.correlationEO);
Cor6 = cat(1,out61.correlationEO,out62.correlationEO,out63.correlationEO,out64.correlationEO,out65.correlationEO,out612.correlationEO,out613.correlationEO,out614.correlationEO);
Cor10 = cat(1,out101.correlationEO,out102.correlationEO,out103.correlationEO,out104.correlationEO,out106.correlationEO,out108.correlationEO,out1012.correlationEO,out1014.correlationEO);
Cor8 = cat(1,out805.correlationEO,out806.correlationEO,out807.correlationEO);
%}