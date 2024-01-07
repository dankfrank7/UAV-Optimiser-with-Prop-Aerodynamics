% clc
clear
warning off

filename = 'inputs/X57_Cruise.vap';
%filename = 'inputs/vap_test.vap'
%filename = 'inputs/TMotor.vap''
filename = 'inputs/Goland_Wing.vap';
%filename = 'Tom_WINGTIProot.vap';
%filename = 'inputs/Misc/StandardCirrusTail2.vap'

VAP_IN = [];
OUTP = fcnVAP_MAIN(filename, VAP_IN);

