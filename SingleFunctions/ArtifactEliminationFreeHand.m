figure
imagesc(bsmooth2)
colormap(gray)

prompt = 'To Keep it [ENTER]/Delete [d]/Finish(point not saved) [f]:';
more = true;
wrongUI = true;
pts = zeros(500,2);
count = 0;

while more == true
    count = count +1;
    h = imfreehand;
    ui = input(prompt,'s');
    wrongUI = true;
    while wrongUI == true
        if strcmp(ui,'')
            pts(1:length(getPosition(h)),:,count) = getPosition(h) ;
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

pt = reshape(permute(pts,[2 1 3]), size(pts,2),[]);
pt1 = pt(1,:);
pt2 = pt(2,:);
pt1 = pt1(pt1>0);
pt2 = pt2(pt2>0);
pt = cat(1,pt1,pt2);

data = -bsmooth2;