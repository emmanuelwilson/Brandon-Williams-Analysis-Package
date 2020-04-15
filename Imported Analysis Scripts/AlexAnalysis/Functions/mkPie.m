function h1 = mkPie(d)
%     if nargin < 2 || isempty(labels)
%         labels = length(d{1}(1,:));
%     end

    degStep = 0.001;
    radii = [0.5 1 1.5];
%     figure
%     set(gcf,'position',[50 50 250.*size(d')],'color','w')
    colors = [{[0.9 0.4 0.4; 0.9 0.6 0.6]} {[0.2 0.2 0.2]} {[0.4 0.4 0.9; 0.6 0.6 0.9]}];
    colors = [{[0.9 0.7 0.7; 0.2 0.2 0.2; 0.7 0.7 0.9]} {[0.9 0.4 0.4; 0.9 0.6 0.6; ...
        0.2 0.2 0.2; 0.4 0.4 0.9; 0.6 0.6 0.9]}];
    
    for ai = 1:length(d(1,:))
        for coni = 1:length(d(:,1))
%             subplot(length(d(:,1)),length(d(1,:)),(coni-1).*length(d(1,:))+ai)
            
            curD = d{coni,ai};
%             for i = 2:length(curD(:,1))
%                 for j = 1:length(curD(1,:))
%                     curD{1,j} = [curD{1,j}; curD{i,j}];
%                 end
%             end
%             curD = cellfun(@nansum,curD(1,:),'uniformoutput',false);
            
            innerCount = 0;
            cts = cellfun(@nansum,curD);
            cts = cts./nansum(cts);
            for i = 1:length(cts)
                start = nansum(cts(1:i-1));
                if isempty(start)
                    start = 0;
                end
                h1(i) = patch('Faces',[1:length([start:degStep:start+cts(i)]).*2],...
                    'Vertices',[[radii(1).*cos(2.*pi.*[start:degStep:start+cts(i)]) ...
                    radii(2).*cos(2.*pi.*[start+cts(i):-degStep:start])]',...
                    [radii(1).*sin(2.*pi.*[start:degStep:start+cts(i)]) ...
                    radii(2).*sin(2.*pi.*[start+cts(i):-degStep:start])]'],'FaceColor',colors{1}(i,:),'edgecolor','w');
                innerCount = innerCount+1;
                if length(curD{i})>1
                    cts2 = curD{i}./nansum(curD{i});
                    for j = 1:length(cts2)
                        start2 = nansum(cts2(1:j-1));
                        if isempty(start2)
                            start2 = 0;
                        end
                        finish = start+cts(i).*start2;
                        h2(j) = patch('Faces',[1:length([start+cts(i).*start2:degStep:start+cts(i).*nansum(cts2(1:j))]).*2],...
                            'Vertices',[[radii(2).*cos(2.*pi.*[start+cts(i).*start2:degStep:start+cts(i).*nansum(cts2(1:j))]) ...
                            radii(3).*cos(2.*pi.*[start+cts(i).*nansum(cts2(1:j)):-degStep:start+cts(i).*start2])]',...
                            [radii(2).*sin(2.*pi.*[start+cts(i).*start2:degStep:start+cts(i).*nansum(cts2(1:j))]) ...
                            radii(3).*sin(2.*pi.*[start+cts(i).*nansum(cts2(1:j)):-degStep:start+cts(i).*start2])]'],'FaceColor',colors{2}(innerCount,:),...
                            'edgecolor','w');
                        if j <length(curD{i})
                            innerCount = innerCount+1;
                        end
                    end
                else
                    
                end
            end
            axis equal
            axis off
%             legend([h1 h2],[labels{1} labels{2}],'box','off','location','northwest');
        end
    end
end