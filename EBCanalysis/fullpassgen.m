%% fullpass generation, takes the cells which passed the MRL crit and applies the 45deg pref firing crit across splithalves. 

EBC = EBCA21;
p = passedA21';
splitangle = rad2deg(EBC.headangle1 - EBC.headangle2);
fullpassA21 = p((splitangle(p) > -45 & splitangle(p) < 45) |  (splitangle(p) > -315 & splitangle(p) < -405) | (splitangle(p) > 315 & splitangle(p) < 405));
EBC = EBCA22;
p = passedA22';
splitangle = rad2deg(EBC.headangle1 - EBC.headangle2);
fullpassA22 = p((splitangle(p) > -45 & splitangle(p) < 45) |  (splitangle(p) > -315 & splitangle(p) < -405) | (splitangle(p) > 315 & splitangle(p) < 405));
EBC = EBCA23;
p = passedA23';
splitangle = rad2deg(EBC.headangle1 - EBC.headangle2);
fullpassA23 = p((splitangle(p) > -45 & splitangle(p) < 45) |  (splitangle(p) > -315 & splitangle(p) < -405) | (splitangle(p) > 315 & splitangle(p) < 405));
EBC = EBCA24;
p = passedA24';
splitangle = rad2deg(EBC.headangle1 - EBC.headangle2);
fullpassA24 = p((splitangle(p) > -45 & splitangle(p) < 45) |  (splitangle(p) > -315 & splitangle(p) < -405) | (splitangle(p) > 315 & splitangle(p) < 405));