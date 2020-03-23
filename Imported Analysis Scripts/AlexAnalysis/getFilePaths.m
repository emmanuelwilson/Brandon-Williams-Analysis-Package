function p = getFilePaths(root,type)
    p = [];
    stack = dir(root);
    stack = cat(1,{stack(:).name});
    stack(ismember(stack,[{'.'} {'..'}])) = [];
    for i = 1:length(stack)
        stack{i} = [root '/' stack{i}];
    end    
    while ~isempty(stack)
        if isdir(stack{1})
            add = dir(stack{1});
            add = cat(1,{add(:).name});
            add(ismember(add,[{'.'} {'..'}])) = [];
            for i = 1:length(add)
                add{i} = [stack{1} '/' add{i}];
            end
            stack = [stack add];
        elseif ismember({stack{1}(end-(length(type)-1):end)},{type})
            p = [p; stack(1)];
        end
        stack(1) = [];
    end
end