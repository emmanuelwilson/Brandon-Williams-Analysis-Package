function h1 = mkSimplePie(d)
%     if nargin < 2 || isempty(labels)
%         labels = length(d{1}(1,:));
%     end
    if ~iscell(d)
        d = {d};
    end
    
    degStep = 0.001;
    radii = [0.5 1 1.5];
    figure
    set(gcf,'position',[50 50 250.*size(d')],'color','w')
    colors = [0.9 0.3 0.3; 0.9 0.6 0.6; 0.1 0.1 0.9 ;0.5 0.3 0.9; 0.7 0.3 0.9];
    
    for i = 1:length(d(1,:))
        for j = 1:length(d(:,1))
            subplot(length(d(:,1)),length(d(1,:)),(i-1).*length(d(1,:))+j)
            hold on
            for k = 1:nanmax(d{i,j})
                cts = nansum(d{i,j}==k)./nansum(~isnan(d{i,j}));
                for i = 1:length(cts)
                    h1(i) = patch('Faces',[1:length([start:degStep:start+cts(i)]).*2],...
                        'Vertices',[[radii(1).*cos(2.*pi.*[start:degStep:start+cts(i)]) ...
                        radii(2).*cos(2.*pi.*[start+cts(i):-degStep:start])]',...
                        [radii(1).*sin(2.*pi.*[start:degStep:start+cts(i)]) ...
                        radii(2).*sin(2.*pi.*[start+cts(i):-degStep:start])]'],'FaceColor',colors(i,:),'edgecolor','w');

                end
            end
            axis equal
            axis off
        end
    end
end