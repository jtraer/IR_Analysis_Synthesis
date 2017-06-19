function [R,C,Mt,Amnd]=Input_IR_Survey_2;

% Path to Recordings and Golay codes
rcnt=0;
rcnt=rcnt+1; R(rcnt).Pth='RecordedAudio/Reverb_Survey_2/4022'; R(rcnt).Gpth='RawGolay/golay_44kHz_N17_2min_24bits'; 
rcnt=rcnt+1; R(rcnt).Pth='RecordedAudio/Reverb_Survey_2/4078'; R(rcnt).Gpth='RawGolay/golay_44kHz_N17_2min_24bits'; 

%==> Calibration recordings
ccnt=0;
ccnt=ccnt+1; C(ccnt).Pth='RecordedAudio/Reverb_Survey_2/CAL_TscmZpp'; C(ccnt).Gpth='RawGolay/golay_44kHz_N17_2min_24bits'; 

% Specify MetaData to be saved with recording
mcnt=0;
%==> Recorder info
mcnt=mcnt+1;Mt{mcnt}='Meta.App.Mic';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.Recorder';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.Gain';
%==> Speaker info
mcnt=mcnt+1;Mt{mcnt}='Meta.App.Speaker';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.Volume';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.PolarAngle_fromTop';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.AzimuthalAngle_fromFront';
%==> Recording info
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Class';
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.SpaceName';
%===> Room reverb
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Size';
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.WallProximity';
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Distance';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.NoPeople';
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Door';
%===> Object reverb
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Size';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Material';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Location';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Damping';

% amend the criteria for plots
acnt=0;
acnt=acnt+1; Amnd(acnt).Prm='Meta.App.Mic';
vcnt=0;
vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.App.PolarAngle_fromTop'; Amnd(acnt).Val(vcnt).Exp='90';
vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.App.AzimuthalAngle_fromFront'; Amnd(acnt).Val(vcnt).Exp='0';
