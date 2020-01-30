%% Generates figure with EBC ratemap across contexts and saves in designated folder

function [] = EBCratemapCompAcrossDaysGrid_Object(fullpass1,fullpass2,fullpass3,fullpass4,ebc1,ebc2,ebc3,ebc4,ebc5,ebc6,ebc7,ebc8,ebc9,ebc10,regmap,name)
mkdir(name);
if ebc1.dimX > ebc1.dimY
    FOVsizeA = round(ebc1.dimY/2);
else
    FOVsizeA = round(ebc1.dimX/2);
end

if ebc3.dimX > ebc3.dimY
    FOVsizeB = round(ebc3.dimY/2);
else
    FOVsizeB = round(ebc3.dimX/2);
end

thetaBins3d = deg2rad(linspace(-180,180,round(360)/6));
distanceBinsA = 0:1.2:FOVsizeA;
distanceBinsB = 0:1.2:FOVsizeB;
c = 1;
figure
pause
for i = 1 : length(regmap)  
    isreg = find(regmap(i,:));    
    if find(isreg == 1)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc1.rm(:,:,regmap(i,1))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc1.rm(:,:,regmap(i,1))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('Context A')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,1))];
        text(20,coordlims(4),info);
        if ~isempty(find(fullpass1 == regmap(i,1)))
            text(20,coordlims(3),'EBC');
        end
    end
    c = c+1;
    if find(isreg == 2)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc2.rm(:,:,regmap(i,2))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc2.rm(:,:,regmap(i,2))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('Context A2')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,2))];
        text(20,coordlims(4),info);
        if ~isempty(find(fullpass2 == regmap(i,2)))
            text(20,coordlims(3),'EBC');
        end
    end
    c = c+1;
    if find(isreg == 3)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsB(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc3.rm(:,:,regmap(i,3))); shading interp        
        xlim([-36 36])
        ylim([-36 36])
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc3.rm(:,:,regmap(i,3))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('Context B')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,3))];
        text(32,coordlims(4),info);
        if ~isempty(find(fullpass3 == regmap(i,3)))
            text(32,coordlims(3),'EBC');
        end
    end
    c = c+1;
    if find(isreg == 4)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc4.rm(:,:,regmap(i,4))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc4.rm(:,:,regmap(i,4))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('Context A3')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,4))];
        text(20,coordlims(4),info);
        if ~isempty(find(fullpass4 == regmap(i,4)))
            text(20,coordlims(3),'EBC');
        end
    end
    c = c+1;
    if find(isreg == 5)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc5.ratemapO(:,:,regmap(i,5))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc5.ratemapO(:,:,regmap(i,5))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('No Milk')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,5))];
        text(20,coordlims(4),info);
    end
    c = c+1;
    if find(isreg == 6)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc6.ratemapO(:,:,regmap(i,6))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc6.ratemapO(:,:,regmap(i,6))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('NW Milk')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,6))];
        text(20,coordlims(4),info);
    end
    c = c+1;
    if find(isreg == 7)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc7.ratemapO(:,:,regmap(i,7))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc7.ratemapO(:,:,regmap(i,7))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('SE Milk')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,7))];
        text(20,coordlims(4),info);
    end
    c = c+1;
    if find(isreg == 8)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc8.ratemapO(:,:,regmap(i,8))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc8.ratemapO(:,:,regmap(i,8))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('NW Milk')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,8))];
        text(20,coordlims(4),info);
    end
    c = c+1;
    if find(isreg == 9)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc9.ratemapO(:,:,regmap(i,9))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc9.ratemapO(:,:,regmap(i,9))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('SE Milk')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,9))];
        text(20,coordlims(4),info);
    end
    c = c+1;
    if find(isreg == 10)
        subplot(5,10,c)
        % the +pi/2 brings "forwards" to "up"
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsA(1:end-1));
        [x, y] = pol2cart(t2,r2);
        h=surface(x,y, ebc10.ratemapO(:,:,regmap(i,10))); shading interp        
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis square
        axis off
        colormap(jet)
        caxis([0 max(max(ebc10.ratemapO(:,:,regmap(i,10))))])
        set(gca,'YDir','Normal')
        freezeColors
        title('NW Milk')
        coordlims = axis;
        info = ['Cell:' num2str(regmap(i,10))];
        text(20,coordlims(4),info);
    end
    c = c+1;
    if mod(i,5) == 0 || i == length(regmap)
        saveas(gcf,[name,'/',num2str(i),'EBCdays.jpg']);
        pause(0.01)
        clf
        c = 1;
    end
end
end


