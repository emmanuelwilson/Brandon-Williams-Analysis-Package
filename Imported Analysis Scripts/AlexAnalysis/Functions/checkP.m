function checkP(p)
    for i = find(ismember(p,'\/'))
        if ~isdir(p(1:i-1))
            mkdir(p(1:i-1))
        end
    end
end