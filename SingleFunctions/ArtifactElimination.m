%%Point and click light artifact remover. 

pixY = length(msCam1(1).cdata(:,1));
pixX = length(msCam1(1).cdata(1,:));
pixelNum = pixX*pixY;
corr = zeros(size(msCam1(1).cdata));
mask = ones(size(corr));

pixelCorr = zeros(pixY,pixX,250);
a = [];
b = zeros(size(corr));
for j = 1 : pixelNum
    for i = 1 : 250
        
        a(i) = msCam1(i).cdata(j);
        
    end
    b(j) =  mean(a);
    corr(j) = mean(diff(a));
end
% b = -b;
bsmooth = [];
bdiff = [];
bdiff1 = [];
bzdiff = bdiff;
ypeaks = zeros(1,500);
yloc = zeros(1,500);
ywidth = zeros(1,500);
yprom = zeros(1,500);
for i  = 1 : pixY
    bsmooth(i,:) = smooth(b(i,:));
    bsmooth2(i,:) = smoothdata(bsmooth(i,:),'SmoothingFactor',0.01);
end

for i = 1 : pixX
    bsmooth2(:,i) = smoothdata(bsmooth2(:,i),'SmoothingFactor',0.01);
%     [high, low] = envelope(bsmooth2(:,i), 10,'peak');
%     bsmooth3(:,i) = (high+low)/2;
    bdiff1(:,i) = diff(bsmooth2(:,i));
    bdiff2 = bdiff1;
    bdiff2(end +1,i) = 0;
    bdiff(:,i) = diff(bdiff2(:,i));
    bzdiff(:,i) = zscore(bdiff(:,i));
end


figure
imagesc(bsmooth2)
colormap(gray)

prompt = 'To Keep it [ENTER]/Delete [d]/Finish(point not saved) [f]:';
more = true;
wrongUI = true;
pts = [];
count = 0;

while more == true
    count = count +1;
    h = impoint;
    ui = input(prompt,'s');
    wrongUI = true;
    while wrongUI == true
        if strcmp(ui,'')
            pts(count,:) = getPosition(h) ;
            wrongUI = false;
        elseif strcmp(ui,'d')
            count = count - 1;
            delete(h);
            wrongUI = false;
        elseif strcmp(ui,'f')
            delete(h);
            wrongUI = false;
            more = false;
        else            
            ui = input(prompt,'s');
        end
    end
end

data = -bsmooth2;

for i = 1 : length(pts(:,1))
    objR = true;
    objL = true;
    [pks,locs,w,p] = findpeaks(data(round(pts(i,2)),:));
    xloc1 = abs(locs - pts(i,1));
    [~,I] = min(xloc1);    
    xpeak = abs(round(locs(I)));    
    ypoint = pts(i,2);
    width = w(I)/2;
    mask(round(ypoint-width):round(ypoint+width),xpeak) = 0;
    xpeak = xpeak-1;
    xpeak1 = xpeak;
    first = true;
    while objR
        xpeak= xpeak+1;
        [pks,locs,w,p] = findpeaks(data(:,xpeak));
        ploc = abs(locs - ypoint);
        [closest,I] = min(ploc);        
        if first
            first = false;
            ypoint = abs(round(locs(I)));
            width = (w(I)/2);
            width1 = width;
            mask(round(ypoint-width):round(ypoint+width),xpeak) = 0;
        elseif closest < width       
        ypoint = abs(round(locs(I)));
        width = (w(I)/2);
        mask(round(ypoint-width):round(ypoint+width),xpeak) = 0;   
        else
            objR = false;
        end
    end
    ypoint = round(pts(i,2));
    first = true;
    xpeak = xpeak1+1;
    while objL
        xpeak= xpeak-1;
        [pks,locs,w,p] = findpeaks(data(:,xpeak));
        ploc = abs(locs - ypoint);
        [closest,I] = min(ploc);
        if first && closest < width1
            first = false;
            ypoint = abs(round(locs(I)));
            width = (w(I)/2);
            mask(round(ypoint-width):round(ypoint+width),xpeak) = 0;
        elseif closest < width
            ypoint = abs(round(locs(I)));
            width = (w(I)/2);
            mask(round(ypoint-width):round(ypoint+width),xpeak) = 0;
        else
            objL = false;
        end
    end
end

b = -bsmooth;
xpeaks = [];
xloc = [];
xwidth = [];
xprominance = [];
for i  = 1 : pixX
    [xpeaks(i,:), xloc(i,:), xwidth(i,:), xprom(i,:)] = findpeaks(b(i,:));
end

Correlation = mean(pixelCorr,3);
