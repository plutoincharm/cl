clc;
clearvars -except GRFx GRFy;
%%%
% global  A B Nx Nu pert MI L m nx ny tx ty g  r lam  vars misc alp alpval indic kc lamall xdata lamx lamy af fx fy
%%%
global  L r m MI misc g 
%Left  leg
L2 = 0.4418;
L3 = 0.4033;
r2 = 0.1872;
r3 = 0.1738;
m2 = 7.4;
m3 = 3.411;
m4 = 1.073;
MI2 = 0.1038;
MI3 = 0.05916;
MI4 = 0.01;

%Right leg
L5 = 0.4418;
L6 = 0.4033;
r5 = 0.1872;
r6 = 0.1738;
m5 = 7.4;
m6 = 3.411;
m7 = 1.073;
MI5 = 0.1038;
MI6 = 0.05916;
MI7 = 0.01; 

% HAT
m1 = 50.172;
MI1 = 1.95;


%%%
g = 9.81;
L = [L2 L3 L5 L6];
r = [r2 r3 r5 r6];
m = [m1 m2 m3 m4 m5 m6 m7];
MI = [MI1 MI2 MI3 MI4 MI5 MI6 MI7];
%%%
coord = readmatrix("JointCoords3D_edit.xlsx");
GRdata = readmatrix("GRdata.xlsx");

gtc = 2:126;  % gait cycle time is from 0 to 1.24 seconds
time = coord(:,1);
% frame starts from 2 , ends at 126

% total 125 frames

% tfit is the time taken 
tfit = time(gtc);
dt = 0.01; % dt should be set to 0.01 only
tt = 0:dt:1.24;

% GRdata has values from 1 to 216
fx = GRdata(gtc,8);
fy = GRdata(gtc,9);

%
lcop_x = 2-0.001*coord(:,38);
lcop_y = 0.001*coord(:,39)+0.02;
rcop_x = 2-0.001*coord(:,35);
rcop_y = 0.001*coord(:,36)+0.02;
lpost_x = 2-0.001*coord(:,14);
lmeta_x = 2-0.001*coord(:,17);
%{
dcopx = diff(rcop_x(77:159));
lcop_x(145) = lpost_x(145);
lcop_x(91) = lmeta_x(91);
for ii = 146:159
    lcop_x(ii) = lcop_x(ii-1) + dcopx(-145 + ii);
end
for ii = 90:-1:77
    lcop_x(ii) = lcop_x(ii+1) - dcopx(end + ii - 90);
end
%}
xtrunk = 2-0.001*coord(gtc,2);
ytrunk = 0.001*coord(gtc,3);
%%% left side
lxhip = 2-0.001*coord(gtc,5);
lyhip = 0.001*coord(gtc,6);
lxknee = 2-0.001*coord(gtc,8);
lyknee = 0.001*coord(gtc,9);
lxankle = 2-0.001*coord(gtc,11);
lyankle = 0.001*coord(gtc,12);
lxpost = 2-0.001*coord(gtc,14);
lypost = 0.001*coord(gtc,15);
lxmeta = 2-0.001*coord(gtc,17);
lymeta = 0.001*coord(gtc,18);
lxcop = lcop_x(gtc);
lycop = lcop_y(gtc);
%%% right side
rxhip = 2-0.001*coord(gtc,20);
ryhip = 0.001*coord(gtc,21);
rxknee = 2-0.001*coord(gtc,23);
ryknee = 0.001*coord(gtc,24);
rxankle = 2-0.001*coord(gtc,26);
ryankle = 0.001*coord(gtc,27);
rxpost = 2-0.001*coord(gtc,29);
rypost = 0.001*coord(gtc,30);
rxmeta = 2-0.001*coord(gtc,32);
rymeta = 0.001*coord(gtc,33);
rxcop = rcop_x(gtc);
rycop = rcop_y(gtc);
%%% mid-point of pelvis
x_p = (lxhip+rxhip)/2;
y_p = (lyhip+ryhip)/2;
rxfoot = 1/3*(rxankle + rxpost + rxmeta);
ryfoot = 1/3*(ryankle + rypost + rymeta);
lxfoot = 1/3*(lxankle + lxpost + lxmeta);
lyfoot = 1/3*(lyankle + lypost + lymeta);

a4 = mean(sqrt((rxankle-rxfoot).^2 + (ryankle-ryfoot).^2));
b4 = mean(sqrt((rxpost-rxfoot).^2 + (rypost-ryfoot).^2));
c4 = mean(sqrt((rxmeta-rxfoot).^2 + (rymeta-ryfoot).^2));
d4 = mean(sqrt((rxpost-rxankle).^2 + (rypost-ryankle).^2));
e4 = mean(sqrt((rxpost-rxmeta).^2 + (rypost-rymeta).^2));
u4 = mean(sqrt((rxmeta-rxankle).^2 + (rymeta-ryankle).^2));
f4 = mean(sqrt((rxmeta-rxankle).^2 + (rymeta-ryankle).^2));   % angle between line from foot to COP

a7 = mean(sqrt((lxankle-lxfoot).^2 + (lyankle-lyfoot).^2));
b7 = mean(sqrt((lxpost-lxfoot).^2 + (lypost-lyfoot).^2));
c7 = mean(sqrt((lxmeta-lxfoot).^2 + (lymeta-lyfoot).^2));
d7 = mean(sqrt((lxpost-lxankle).^2 + (lypost-lyankle).^2));
e7 = mean(sqrt((lxpost-lxmeta).^2 + (lypost-lymeta).^2));
f7 = mean(sqrt((rxmeta-rxankle).^2 + (rymeta-ryankle).^2)); % angle between line from foot to COP

gamma61 = acos((a4^2 + b4^2 - d4^2)/(2*a4*b4));
gamma62 = gamma61 + acos((b4^2+c4^2-e4^2)/(2*b4*c4));
gamma71 = acos((a7^2+b7^2-d7^2)/(2*a7*b7));
gamma72 = gamma71 + acos((b7^2+c7^2-e7^2)/(2*b7*c7));
%%% Sequence: SS2 -> DS1 -> SS1 -> 
%%% SSP -2
r4SS2 = a7;
r7SS2 = a4;
r4t = c7;         
%r4h = b7;
r4h = b7;
gamma61 = gamma71;
gamma62 = gamma72;

tht7SS2 =  atan2(ryankle(55:76)-ryfoot(55:76),rxankle(55:76)-rxfoot(55:76));
tht4SS2 =  atan2(lyankle(55:76)-lyfoot(55:76),lxankle(55:76)-lxfoot(55:76));



lt = mean(sqrt((x_p-xtrunk).^2 + (y_p-ytrunk).^2));
d1 = mean(sqrt((lxmeta-lxankle).^2));
d3 = mean(sqrt((lxpost-lxankle).^2 ));
h1 = mean(sqrt((lymeta-lyankle).^2));
h3 = mean(sqrt((lypost-lyankle).^2 ));


lycop(20:76) = lypost(20:76);                                
%lycop(20:76) = zeros(57,1);
r4 = a7;
r4c = mean(sqrt((lxfoot(20:76)-lxcop(20:76)).^2 + (lyfoot(20:76)-lycop(20:76)).^2));
r4d = mean(sqrt((lxankle(20:76)-lxcop(20:76)).^2 + (lyankle(20:76)-lycop(20:76)).^2));
gamma43 = 2*pi -acos((r4^2 + r4c^2 - r4d^2)/(2*r4*r4c));
varsSS2 = [r4SS2,r7SS2,r4t,r4h,r4c,gamma43];





%%% angles for inverse dynamics code
x1 = xtrunk;
y1 = ytrunk;
tht1 = atan2(ytrunk-y_p,xtrunk-x_p);
tht5 = atan2(ryhip-ryknee,rxhip-rxknee);
tht2 = atan2(lyhip-lyknee,lxhip-lxknee);
tht6 = atan2(ryknee-ryankle,rxknee-rxankle);
tht3 = atan2(lyknee-lyankle,lxknee-lxankle);
%tht4 =  atan2(lyankle-lypost ,lxankle-lxpost );
tht4 = atan2(lyankle-lyfoot,lxankle-lxfoot);
tht7 =  atan2(ryankle-ryfoot ,rxankle-rxfoot );
%tht4 = tht4SS2; %[tht4SS2 ; tht6DS1 ; tht6SS1 ; tht6DS2];
%tht7 = [tht7SS2 ; tht7DS1 ; tht7SS1 ; tht7DS2];
%%% fitting fourier series for theta
xk1 = lxknee;
yk1 = lyknee;
f1 = fit(tfit,tht1,'fourier8');
tht1 = f1(tt);
f2 = fit(tfit,tht2,'fourier8');
tht2 = f2(tt);
f3 = fit(tfit,tht3,'fourier8');
tht3 = f3(tt);
f4 = fit(tfit,tht4,'fourier8');
tht4 = f4(tt);
f5 = fit(tfit,tht5,'fourier8');
tht5 = f5(tt);
f6 = fit(tfit,tht6,'fourier8');
tht6 = f6(tt);
f7 = fit(tfit,tht7,'fourier8');
tht7 = f7(tt);
fx1 = fit(tfit,x1,'fourier8');
x1 = fx1(tt);
fy1 = fit(tfit,y1,'fourier8');
y1 = fy1(tt);

fkx1 = fit(tfit,xk1,'fourier8');
xk1 = fkx1(tt);
fky1 = fit(tfit,yk1,'fourier8');
yk1 = fky1(tt);

%%% misc
% del2 = mean(tht1 - atan2(y1-ryhip,x1-rxhip));
% del3 = mean(pi - (pi+tht1) + atan2(y1-lyhip,x1-lxhip));
R3 = mean(sqrt((rxhip-xtrunk).^2 + (ryhip-ytrunk).^2));
R2 = mean(sqrt((lxhip-xtrunk).^2 + (lyhip-ytrunk).^2));
del2 = abs(mean(pi-tht1 -(pi- atan2(y1-ryhip,x1-rxhip))));
del3 = abs(mean(tht1 - atan2(y1-ryhip,x1-rxhip)));

midlx  = ((lxpost+lxmeta)/2);
midly = ((lypost+lymeta)/2);
gammamid =  mean(atan2(lyfoot-midly, lxfoot-midlx));    
misc = [R2 R3 del2 del3 gamma61 gamma71 gammamid];

%%% accelerations and velocity
[vx1,ax1] = differentiate(fx1,tt);
[vy1,ay1] = differentiate(fy1,tt);
[omg1,alp1] = differentiate(f1,tt);
[omg2,alp2] = differentiate(f2,tt);
[omg3,alp3] = differentiate(f3,tt);
[omg4,alp4] = differentiate(f4,tt);
[omg5,alp5] = differentiate(f5,tt);
[omg6,alp6] = differentiate(f6,tt);
[omg7,alp7] = differentiate(f7,tt);

[vxk1,axk1] = differentiate(fkx1,tt);
[vyk1,ayk1] = differentiate(fky1,tt);

omg = [vx1';vy1';omg1';omg2';omg3';omg4';omg5';omg6';omg7'];
bc =  [x1';y1';tht1';tht2';tht3';tht4';tht5';tht6';tht7';lxcop';lycop'];
alpval =  [ax1';ay1';alp1';alp2';alp3';alp4';alp5';alp6';alp7'];


%%% printing the position at each time stamp
%{
c = 1;
for k = 1:length(tt)
%for k = 1:214
    
     axis([0 2 -0 2])
   % base =line([-1000 2000],[20 20],'LineWidth',1,'Color','black');
       % pelvic =line([lhx(i) rhip_x(i)],[lhip_z(i) rhip_z(i)],'LineWidth',1,'Color','black');
        T=line([x1(k) lxhip(k)],[y1(k) lyhip(k) ],'LineWidth',1,'Color','black');
         u=line([lxhip(k) lxknee(k)],[lyhip(k) lyknee(k)],'LineWidth',1,'Color','red');
        v=line([lxknee(k) lxankle(k)],[lyknee(k) lyankle(k)],'LineWidth',1,'Color','red');
        w=line([lxankle(k)  lxmeta(k)],[lyankle(k)  lymeta(k)],'LineWidth',1,'Color','red');
        x=line([lxankle(k) lxpost(k)],[lyankle(k) lypost(k)],'LineWidth',1,'Color','red');
        y=line([lxpost(k)  lxmeta(k)],[lypost(k)  lymeta(k)],'LineWidth',1,'Color','red');
        %plot(lxcop(k),lycop(k),'g*')
       % legend('Left leg')
       hold on
        plot(lxcop(k),lycop(k),'g*')

        
       
        Te=line([x1(k) rxhip(k)],[y1(k) ryhip(k) ],'LineWidth',1,'Color','black');
         ue=line([rxhip(k) rxknee(k)],[ryhip(k) ryknee(k)],'LineWidth',1,'Color','blue');
        ve=line([rxknee(k) rxankle(k)],[ryknee(k) ryankle(k)],'LineWidth',1,'Color','blue');
        we=line([rxankle(k)  rxmeta(k)],[ryankle(k)  rymeta(k)],'LineWidth',1,'Color','blue');
        xe=line([rxankle(k) rxpost(k)],[ryankle(k) rypost(k)],'LineWidth',1,'Color','blue');
        ye=line([rxpost(k)  rxmeta(k)],[rypost(k)  rymeta(k)],'LineWidth',1,'Color','blue');
       
       
pause()
end
%}
%hold on
%plot(lxcop,lycop,'g*')


%%% plotting
%{
GRFxr = zeros(136,1);
GRFxr(54:68)=sol(54:68,1);
GRFxr(69:121)=sol(69:121,1)+sol(69:121,2);
GRFxr(122:136)=sol(122:136,1);
figure;
plot(tt,GRFxr,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,GRFx,'r-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_x (N) \rightarrow');
legend('2D Lagrangian','Experimental');
%
GRFyr = zeros(136,1);
GRFyr(54:68)=sol(54:68,3);
GRFyr(69:121)=sol(69:121,3)+sol(69:121,4);
GRFyr(122:136)=sol(122:136,3);
figure;
plot(tt,GRFyr,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,GRFy,'r-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_y (N) \rightarrow');
legend('2D Lagrangian','Experimental');

%}



%%% plotting angles 
%{
figure; hold on; grid on;
plot(tfit,tht1,'-o');
title('Theta left thigh ');
ylabel('Angle in degrees ');
xlabel('time frame(s)');

figure; hold on; grid on;
plot(tfit,tht2,'-o');
title('Theta left leg');
ylabel('Angle in degrees ');
xlabel('time frame(s)');

figure; hold on; grid on;
plot(tfit,tht3,'-o');
title('Theta right thigh');
ylabel('Angle in degrees ');
xlabel('time frame(s)');

figure; hold on; grid on;
plot(tfit,tht4,'-o');
title('Theta right leg');
ylabel('Angle in degrees ');
xlabel('time frame(s)');
%time=cord(1:213,1);

%}
writematrix(bc,'BC.xlsx');
writematrix(omg,'omg.xlsx');
%writematrix(alp,'alp.xlsx');



vars = varsSS2;
%end
