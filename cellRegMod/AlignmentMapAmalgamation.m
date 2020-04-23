%% Creates a single registration map acorss days using pairwise registration values
function [Singlemap] = AlignmentMapAmalgamation(alignment)

Singlemap = ReorganizeAlignmentMap(alignment.alignmentMap);
[Singlemap, avg_psame] = ElimConflict2020(Singlemap,alignment.alignmentMap,alignment.scoreMap);
Singlemap = FindMissingCells2020(Singlemap,[],alignment);
end