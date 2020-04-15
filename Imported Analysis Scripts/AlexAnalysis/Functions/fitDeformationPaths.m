function [dp1x dp2x dp1y dp2y] = fitDeformationPaths(p1,p2)
    tp1 = p1;
    tp2 = p2;
    
    tp1f = floor(tp1);
    tp2f = floor(tp2);
    
    boundsMapX = [];
    boundsMapY = [];
    for i = 0:nanmax([tp1f(1,:) tp2f(1,:)])
        Ay = [nanmin(tp1f(2,tp1f(1,:)==i)) nanmax(tp1f(2,tp1f(1,:)==i))];
        By = [nanmin(tp2f(2,tp2f(1,:)==i)) nanmax(tp2f(2,tp2f(1,:)==i))];
        boundsMapX = [boundsMapX; i Ay By];
    end
    for j = 0:nanmax([tp1f(2,:) tp2f(2,:)])
        Ax = [nanmin(tp1f(1,tp1f(2,:)==j)) nanmax(tp1f(1,tp1f(2,:)==j))];
        Bx = [nanmin(tp2f(1,tp2f(2,:)==j)) nanmax(tp2f(1,tp2f(2,:)==j))];
        boundsMapY = [boundsMapY; j Ax Bx];
    end
    
    %%% Compress on Y
    
    tmp = [diff(boundsMapX(:,[2 3]),[],2) diff(boundsMapX(:,[4 5]),[],2)];
    rescale = bsxfun(@rdivide,nanmin(tmp,[],2),tmp);
    shift = nanmax(boundsMapX(:,[2 4]),[],2);
    
    dp1x = p1;
    [a b] = ismember(tp1f(1,:),boundsMapX(:,1));
    dp1x(2,:) = dp1x(2,:).*rescale(b,1)' + [shift(b)-boundsMapX(b,2)]';
    
    dp2x = p2;
    [a b] = ismember(tp2f(1,:),boundsMapX(:,1));
    dp2x(2,:) = dp2x(2,:).*rescale(b,2)' + [shift(b)-boundsMapX(b,4)]';
    
    %%% Compress on X
    
    tmp = [diff(boundsMapY(:,[2 3]),[],2) diff(boundsMapY(:,[4 5]),[],2)];
    rescale = bsxfun(@rdivide,nanmin(tmp,[],2),tmp);
    shift = nanmax(boundsMapY(:,[2 4]),[],2);
    
    dp1y = p1;
    [a b] = ismember(tp1f(2,:),boundsMapY(:,1));
    dp1y(1,:) = dp1y(1,:).*rescale(b,1)' + [shift(b)-boundsMapY(b,2)]';
    
    dp2y = p2;
    [a b] = ismember(tp2f(2,:),boundsMapY(:,1));
    dp2y(1,:) = dp2y(1,:).*rescale(b,2)' + [shift(b)-boundsMapY(b,4)]';
end




















