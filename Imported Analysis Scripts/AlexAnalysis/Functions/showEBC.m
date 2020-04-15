function showEBC(m)
    
    tstep = 2.*pi./length(m(:,1));
    dstep = 1;
    for i = 1:length(m(:,1))
        for j = 1:length(m(1,:))
            patch([(j-1).*dstep.*sin(((i-1).*tstep)) (j-1).*dstep.*sin(((i).*tstep)) ...
                (j).*dstep.*sin(((i).*tstep)) (j).*dstep.*sin(((i-1).*tstep))], ...
                [(j-1).*dstep.*cos(((i-1).*tstep)) (j-1).*dstep.*cos(((i).*tstep)) ...
                (j).*dstep.*cos(((i).*tstep)) (j).*dstep.*cos(((i-1).*tstep))],m(i,j),'edgecolor','none');
        end
    end
    axis equal
end