function [R,C,Mt,Amnd,html_tmp]=Input_IRSurvey_NatStats;

                                % Path to html template
  html_tmp='RvrbStrct_Summary.html';

% Path to Recordings and Golay codes
rcnt=0;
rcnt=rcnt+1; R(rcnt).Pth='RecordedAudio/FullSurvey/'; R(rcnt).Gpth='RawGolay/golay_44kHz_N19_3min_24bits'; 

%==> Calibration recordings
ccnt=0;
ccnt=ccnt+1; C(ccnt).Pth='RecordedAudio/FullSurvey/CAL_IonTscm/'; C(ccnt).Gpth='RawGolay/golay_44kHz_N19_3min_24bits'; 
%ccnt=ccnt+1; C(ccnt).Pth='RecordedAudio/Reverb_Survey_2/CAL_Tscm2Zpp'; C(ccnt).Gpth='RawGolay/golay_44kHz_N17_2min_24bits'; 

% Specify MetaData to be saved with recording
mcnt=0;
%==> Recorder info
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Mic';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Recorder';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Gain';
%==> Speaker info
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Speaker';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Volume';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.PolarAngle_fromTop';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.AzimuthalAngle_fromFront';
%==> Recording info
mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Class';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.SpaceName';
%===> Room reverb
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Size';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.WallProximity';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Distance';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.NoPeople';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Door';
%===> Object reverb
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Size';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Material';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Location';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Damping';

% amend the criteria for plots
Amnd=[];
%acnt=0;
%acnt=acnt+1; Amnd(acnt).Prm='Meta.App.Mic';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.App.PolarAngle_fromTop'; Amnd(acnt).Val(vcnt).Exp='90';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.App.AzimuthalAngle_fromFront'; Amnd(acnt).Val(vcnt).Exp='0';
%
%acnt=acnt+1; Amnd(acnt).Prm='Meta.App.AzimuthalAngle_fromFront';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Class'; Amnd(acnt).Val(vcnt).Exp='CAL';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.App.PolarAngle_fromTop'; Amnd(acnt).Val(vcnt).Exp='90';
%
%acnt=acnt+1; Amnd(acnt).Prm='Meta.App.PolarAngle_fromTop';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Class'; Amnd(acnt).Val(vcnt).Exp='CAL';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.App.AzimuthalAngle_fromFront'; Amnd(acnt).Val(vcnt).Exp='0';
%
%acnt=acnt+1; Amnd(acnt).Prm='Meta.Env.Class';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Distance'; Amnd(acnt).Val(vcnt).Exp='1.6';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.WallProximity'; Amnd(acnt).Val(vcnt).Exp='NA';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Door'; Amnd(acnt).Val(vcnt).Exp='NA';
%
%acnt=acnt+1; Amnd(acnt).Prm='Meta.Env.Distance';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.WallProximity'; Amnd(acnt).Val(vcnt).Exp='NA';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Door'; Amnd(acnt).Val(vcnt).Exp='NA';
%
%acnt=acnt+1; Amnd(acnt).Prm='Meta.Env.Size';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Distance'; Amnd(acnt).Val(vcnt).Exp='1.6';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.WallProximity'; Amnd(acnt).Val(vcnt).Exp='NA';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Door'; Amnd(acnt).Val(vcnt).Exp='NA';
%
%acnt=acnt+1; Amnd(acnt).Prm='Meta.Env.SpaceName';
%vcnt=0;
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Distance'; Amnd(acnt).Val(vcnt).Exp='1.6';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.WallProximity'; Amnd(acnt).Val(vcnt).Exp='NA';
%vcnt=vcnt+1; Amnd(acnt).Var(vcnt).Exp='Meta.Env.Door'; Amnd(acnt).Val(vcnt).Exp='NA';
