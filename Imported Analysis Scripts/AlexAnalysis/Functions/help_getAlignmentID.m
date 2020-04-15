function id = help_getAlignmentID(alignments,nFold,paths)
    id = nan;
    isFold = cat(1,alignments(:).nFold)==nFold;
    isPaths = false(length(alignments),1);
    for i = 1:length(paths)
        paths{i} = paths{i}(find(ismember(paths{i},'/'),1,'last')+1:end);
    end
    for i = 1:length(alignments)
        if length(alignments(i).sessions)~=length(paths)
            continue
        end
        tp = alignments(i).sessions;
        for j = 1:length(tp)
            tp{j} = tp{j}(find(ismember(tp{j},'/'),1,'last')+1:end);
        end
        if all(cellfun(@strcmp,tp,paths))
            isPaths(i) = true;
        end
    end
    if any(isFold&isPaths)
        id = find(isFold&isPaths,1,'first');
    end
end