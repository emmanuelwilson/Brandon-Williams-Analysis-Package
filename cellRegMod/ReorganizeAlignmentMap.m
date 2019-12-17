%% Organizes Alignment Map into a single matrix

%INPUT: 
%   - map: Cell registration map in cell matrix form for each session
%OUTPUT:
%   - alignmentMap: n by m matrix where n is the number of cells pairs and
%   m the session number. Each number represents the cell index for in that
%   session.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignmentMap = ReorganizeAlignmentMap(map)
nmap = [];

%look through all cells
for i = 1 : length(map)-1
    for j = i+1 : length(map)
        if j == i +1 && i == 1
            nmap = map{i,j};
        else
            %add cell combinations
            for k = 1 : length(map{i,j}(:,1))
                if all(map{i,j}(k,:)~=0)
                    nmap(find(nmap(:,i)== map{i,j}(k,1)),j) = map{i,j}(k,2);
                elseif map{i,j}(k,1) == 0
                    addon = zeros(1,length(nmap(1,:)));
                    addon(j) = map{i,j}(k,2);
                    temp = cat(1,nmap,addon);
                    nmap = temp;
                else
                    nmap(find(nmap(:,i)== map{i,j}(k,1)),j) = 0;                    
                end
            end
        end
    end
end
%reorganize results and eliminate repeated rows
nmap = unique(nmap,'first','rows');
nmap = flip(nmap,1);

%Elimates rows that contain parts of a more complete sequence
for i = 1 : length(nmap(1,:))
    for j = 1 : length(nmap(:,1))
        if (nmap(j,i) > 0) && (length(find(nmap(:,i) == j)) > 1)
            ind = find(nmap(:,i) == j);
            elim = [];
            count = 1; 
            for a = 1 : length(ind)-1                
                for b = a+1 : length(ind)
                    comp = nmap(ind(a),i:end)-nmap(ind(b),i:end);
                    if sum(comp) == 0 && sum(nmap(ind(b), 1: i-1)) == 0
                        elim(count) = ind(b);
                        count = count +1;
                    elseif sum(comp) == 0 && sum(nmap(ind(a), 1: i-1)) == 0
                        elim(count) = ind(a);
                        count = count +1;
                    end
                end
            end
            nmap(elim,:) = 0;
        end
    end
end



clean = sum(nmap,2);
elim = find(clean==0);
nmap(elim,:) = [];



alignmentMap = nmap;

end