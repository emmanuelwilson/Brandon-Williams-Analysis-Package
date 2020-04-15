function isMostRecent = nthMostRecent(isIn,n)
    isIn = [false(length(isIn(:,1)),1) diff(isIn,[],2)>0];
    isMostRecent = isIn;
    list = false(length(isIn(:,1)),n);
    for i = 2:length(isIn(1,:))
        if all(isIn(:,i)==0)
            isMostRecent(:,i) = list(:,n);
        else
            for q = 2:n
                list(:,q) = list(:,q-1);
            end
            list(:,1) = isMostRecent(:,i);
        end
    end
end