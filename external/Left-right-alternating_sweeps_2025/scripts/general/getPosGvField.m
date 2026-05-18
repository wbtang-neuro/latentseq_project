function [fdFine, fdCoarse, arenaName] = getPosGvField(sessionType)
% retrieve the fieldnames for the "extended" position binning grids, given
% a specified session type

 % pos_of_dark and pos_of_novel are equivalent to pos_of, but pos_4m is
 % different
strs = split(sessionType, ["_dark", "_novel"]);

arenaName = strs(1);

fdFine = join(["pos", arenaName, "fine", "extended"], "_");
fdCoarse = join(["pos", arenaName, "coarse", "extended"], "_");
end