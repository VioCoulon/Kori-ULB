figure;
load ASEfor_m2_80_toto
SLC0=SLC;
time0=time;
plot(time,gradient(SLC0)/ctr.dt); hold on; % Frank.
%plot(time,SLC0); hold on;
load ASEfor_m2_80a_toto
SLC(SLC==0)=NaN;
SLC(1)=0;
plot(time+time0(end),gradient(SLC)/ctr.dt); % Frank.
%plot(time+time0(end),SLC);
%plot(time,gradient(SLC)/ctr.dt);