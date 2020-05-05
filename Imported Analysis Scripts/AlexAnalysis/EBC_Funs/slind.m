function ns = slind(s,num)
    tmp = find(ismember(s,'/'));
    if num(1)==0
        ns = s(1:tmp(num(2))-1);
    elseif num(2)==0
        ns = s(tmp(num(1))+1:end);
    else
        ns = s(tmp(num(1))+1:tmp(num(2))-1);
    end
end