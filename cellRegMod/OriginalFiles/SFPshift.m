%Apply FOV shift on "ms" SFP
%INPUT:
%   -ms: minscope structure, needs SFPs field
%   -shift: 1 by 2 vector, the first value is the horizontal shift and the
%   second value is the vertical shift.
%OUTPUT:
%   -ms: same structure as previous but with a shifted SFPs field.
function ms = SFPshift(ms,shift)
    sfp = circshift(ms.SFPs,[round(shift(1)) round(shift(2))]);             %Apply shifts
    %eliminate "spill over" footprints
    %horizontal
    if shift(1)>0
        sfp(1:round(shift(1)),:,:) = 0;
    else
        sfp(end+round(shift(1)):end,:,:) = 0;
    end
    %vertical
    if shift(2)>0
        sfp(:,1:round(shift(2)),:) = 0;
    else
        sfp(:,end+round(shift(2)):end,:) = 0;
    end
        
    ms.SFPs = sfp;                                                          %replace previous footprint
end