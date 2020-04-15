function gm = groupMat(m,label,collapse)
    if nargin < 3 || isempty(collapse)
        collapse = true;
    end
    
    m(logical(eye(size(m)))) = {[]};

    gm = repmat({[]},[nanmax(label)]);
    for i = 1:nanmax(label)
        for j = 1:nanmax(label)
            if collapse
                gm{i,j} = cat(1,m{label==i,label==j},m{label==j,label==i});
            else
                a = cat(1,m{label==i,label==j});
                b = cat(1,m{label==j,label==i});
                if ~isempty(b)
                    b = b(:,[3 4 1 2]);
                end
                gm{i,j} = cat(1,a,b);
            end
        end
    end
end