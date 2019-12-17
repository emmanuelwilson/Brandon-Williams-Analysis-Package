%% Probablity of registration vs MRL difference

function [out] = RegProbVsMRLplot(Singlemap,alignmap,probmap,fullpass1,fullpass2,fullpass3,fullpass4, EBC1, EBC2, EBC3, EBC4)

[~, ind1, ~] = intersect(Singlemap(:,1),fullpass1);
[~, ind2, ~] = intersect(Singlemap(:,2),fullpass2);
[~, ind3, ~] = intersect(Singlemap(:,3),fullpass3);
[~, ind4, ~] = intersect(Singlemap(:,4),fullpass4);

pass12 = intersect(ind1,ind2);
pass13 = intersect(ind1,ind3);
pass14 = intersect(ind1,ind4);
pass23 = intersect(ind2,ind3);
pass24 = intersect(ind2,ind4);
pass34 = intersect(ind3,ind4);

cell12 = Singlemap(pass12,1);
cell21 = Singlemap(pass12,2);

cell13 = Singlemap(pass13,1);
cell31 = Singlemap(pass13,3);

cell14 = Singlemap(pass14,1);
cell41 = Singlemap(pass14,4);

cell23 = Singlemap(pass23,2);
cell32 = Singlemap(pass23,3);

cell24 = Singlemap(pass24,2);
cell42 = Singlemap(pass24,4);

cell34 = Singlemap(pass34,3);
cell43 = Singlemap(pass34,4);

out.prob12 = probfind(alignmap, probmap,cell12,1,2);
out.prob13 = probfind(alignmap, probmap,cell13,1,3);
out.prob14 = probfind(alignmap, probmap,cell14,1,4);
out.prob23 = probfind(alignmap, probmap,cell23,2,3);
out.prob24 = probfind(alignmap, probmap,cell24,2,4);
out.prob34 = probfind(alignmap, probmap,cell34,3,4);

out.mrldiff12 = abs(EBC1.mrall(cell12) - EBC2.mrall(cell21));
out.mrldiff13 = abs(EBC1.mrall(cell13) - EBC3.mrall(cell31));
out.mrldiff14 = abs(EBC1.mrall(cell14) - EBC4.mrall(cell41));
out.mrldiff23 = abs(EBC2.mrall(cell23) - EBC3.mrall(cell32));

out.mrldiff24 = abs(EBC2.mrall(cell24) - EBC4.mrall(cell42));
out.mrldiff34 = abs(EBC3.mrall(cell34) - EBC4.mrall(cell43));

end

function prob = probfind(alignmap, probmap,cells1,s1,s2)
prob = zeros(length(cells1),1);
probmap = probmap{s1,s2};
for i = 1 : length(cells1)
    cellind = find(alignmap{s1,s2}(:,1) == cells1(i));
    t = probmap(cellind);
    prob(i) = t{1}(2,1);
    if isnan(prob(i))
        prob(i) = t{1}(1,1);
    end
end

end