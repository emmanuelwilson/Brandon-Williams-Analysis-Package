%%Splithalf comparison of delays 
function [d1, d2] = delaySplitHalf(ms,dctrace, ditrace,dcctrace,dicctrace,dc,di,dcc,dic, delay, delay2) 

dctrace1 = [];
ditrace1 = [];
dcctrace1 = [];
dicctrace1 = [];

dctrace2 = [];
ditrace2 = [];
dcctrace2 = [];
dicctrace2 = [];

%Find the halfway seperation point, making sure it doesn't split a delay period
half = round(length(ms.FiltTraces(1,:))/2);
mind = 0;
% for i = 1: length(delay)
%     if i >= length(delay2)
%         break
%     elseif half>= delay(i) && half<delay2(i)
%         [~,mind] = min(abs(delay-half));
%         half = delay(mind)-1;
%     end
% end
if mind == 0 
    [~,mind] = min((delay-half)>0);
    half = delay(mind)-1;
end
    
for i =1: length(delay)
    if i < mind
        if dc == delay(i)
            dctrace1(i,:,dccount) = dctrace(:,:,find(dc == delay(i)));
            dccount1 = dccount1 + 1;
        elseif di == delay(i)
            ditrace1(i,:,dicount) = ditrace(:,:,find(di == delay(i)));
            dicount1 = dicount1 + 1;
        elseif dcc == delay(i)
            dcctrace1(i,:,dcccount) = dcctrace(:,:,find(dcc == delay(i)));
            dcccount1 = dcccount1 + 1;
        elseif di == delay(i)
            dicctrace1(i,:,dicccount) = dicctrace(:,:,find(dic == delay(i)));
            dicccount1 = dicccount1 + 1;
        end
    else
        if dc == delay(i)
            dctrace2(i,:,dccount) = dctrace(:,:,find(dc == delay(i)));
            dccount2 = dccount2 + 1;
        elseif di == delay(i)
            ditrace2(i,:,dicount) = ditrace(:,:,find(di == delay(i)));
            dicount2 = dicount2 + 1;
        elseif dcc == delay(i)
            dcctrace2(i,:,dcccount) = dcctrace(:,:,find(dcc == delay(i)));
            dcccount2 = dcccount2 + 1;
        elseif di == delay(i)
            dicctrace2(i,:,dicccount) = dicctrace(:,:,find(dic == delay(i)));
            dicccount2 = dicccount2 + 1;
        end
    end
end
d1.dctrace = dctrace1;
d1.ditrace = ditrace1;
d1.dcctrace = dcctrace1;
d1.dicctrace = dicctrace1;

d2.dctrace = dctrace2;
d2.ditrace = ditrace2;
d2.dcctrace = dcctrace2;
d2.dicctrace = dicctrace2;

end