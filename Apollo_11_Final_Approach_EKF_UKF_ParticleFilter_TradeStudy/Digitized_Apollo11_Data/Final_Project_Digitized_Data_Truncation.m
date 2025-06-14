%% read in digitized measurement data and truncate it%%
clear all;close all;clc;
%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('Range_vs_Time.mat')
Time_R = Range_vs_Time(:,1); %Time sec
Range = Range_vs_Time(:,2); %Range ft
clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_R);
    indx2 = find(Time_R >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_R_trunc(ii,:) = Time_R(indx,:);
        Range_trunc(ii,:) = Range(indx,:);
    end
end
iter = 0;
for ii = 1: length(Time_R_trunc)
    if Time_R_trunc(ii) == 0
        iter = iter + 1;
    end
end
for ii = 1: length(Time_R_trunc)-iter
    Time_R_trunc(ii) = Time_R_trunc(iter+ii);
    Range_trunc(ii) = Range_trunc(iter+ii);
end
figure; plot(Time_R_trunc,Range_trunc,'r-'); ...hold on; 
    figure; plot(Time_R,Range,'b-')

%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('Altitude_vs_Time.mat')
Time_A = Altitude_vs_Time(:,1); %Time sec
Alt = Altitude_vs_Time(:,2); %Altitude ft
clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_A);
    indx2 = find(Time_A >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_A_trunc(ii,:) = Time_A(indx,:);
        Alt_trunc(ii,:) = Alt(indx,:);
    end
end
iter = 0;
for ii = 1: length(Time_A_trunc)
    if Time_A_trunc(ii) == 0
        iter = iter + 1;
    end
end
for ii = 1: length(Time_A_trunc)-iter
    Time_A_trunc(ii) = Time_A_trunc(iter+ii);
    Alt_trunc(ii) = Alt_trunc(iter+ii);
end
figure; plot(Time_A_trunc,Alt_trunc,'r-'); ...hold on; 
    figure; plot(Time_A,Alt,'b-')
%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('Pitch_vs_Time.mat')
Time_P = Pitch_vs_Time(:,1); %Time sec
Pitch = Pitch_vs_Time(:,2); %Pitch rad
clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_P);
    indx2 = find(Time_P >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_P_trunc(ii,:) = Time_P(indx,:);
        Pitch_trunc(ii,:) = Pitch(indx,:);
    end
end
iter = 0;
for ii = 1: length(Time_P_trunc)
    if Time_P_trunc(ii) == 0
        iter = iter + 1;
    end
end
for ii = 1: length(Time_P_trunc)-iter
    Time_P_trunc(ii) = Time_P_trunc(iter+ii);
    Pitch_trunc(ii) = Pitch_trunc(iter+ii);
end

figure; plot(Time_P_trunc,Pitch_trunc,'r-');...hold on; 
figure; plot(Time_P,Pitch,'b-')

%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('Thrust_vs_TimefromIgnition.mat')
Time_T = Thrust_vs_TimefromIgnition(:,1)*60; %Time sec
Thrust = Thrust_vs_TimefromIgnition(:,2); %Thrust lbf

indxT = find(Time_T > 550);

for ii = 1:(length(Time_T)-indxT(1))
    Time_T_trunc1(ii,:) = Time_T(indxT(1)+ii);
    Thrust_trunc1(ii,:) = Thrust(indxT(1)+ii);
end

Time_T_trunc1 = Time_T_trunc1-Time_T_trunc1(1);

clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_T_trunc1);
    indx2 = find(Time_T_trunc1 >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_T_trunc(ii,:) = Time_T_trunc1(indx,:);
        Thrust_trunc(ii,:) = Thrust_trunc1(indx,:);
    end
end

figure; plot(Time_T_trunc,Thrust_trunc,'r-o');...hold on; 
figure; plot(Time_T,Thrust,'b-o')


%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('AltRate_vs_Time.mat')
Time_AR = AltRate_vs_Time(:,1); %Time sec
AltRate = AltRate_vs_Time(:,2); %AltRate ft/sec
clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_AR);
    indx2 = find(Time_AR >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_AR_trunc(ii,:) = Time_AR(indx,:);
        AltRate_trunc(ii,:) = AltRate(indx,:);
    end
end
iter = 0;
for ii = 1: length(Time_AR_trunc)
    if Time_AR_trunc(ii) == 0
        iter = iter + 1;
    end
end
for ii = 1: length(Time_AR_trunc)-iter
    Time_AR_trunc(ii) = Time_AR_trunc(iter+ii);
    AltRate_trunc(ii) = AltRate_trunc(iter+ii);
end
figure; plot(Time_AR_trunc,AltRate_trunc,'r-'); ...hold on; 
    figure; plot(Time_AR,AltRate,'b-')

%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('AltAccel_vs_Time.mat')
Time_AA = AltAccel_vs_Time(:,1); %Time sec
AltAccel = AltAccel_vs_Time(:,2); %AltAccel ft/sec
clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_AA);
    indx2 = find(Time_AA >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_AA_trunc(ii,:) = Time_AA(indx,:);
        AltAccel_trunc(ii,:) = AltAccel(indx,:);
    end
end
iter = 0;
for ii = 1: length(Time_AA_trunc)
    if Time_AA_trunc(ii) == 0
        iter = iter + 1;
    end
end
for ii = 1: length(Time_AA_trunc)-iter
    Time_AA_trunc(ii) = Time_AA_trunc(iter+ii);
    AltAccel_trunc(ii) = AltAccel_trunc(iter+ii);
end
figure; plot(Time_AA_trunc,AltAccel_trunc,'r-'); ...hold on; 
    figure; plot(Time_AA,AltAccel,'b-')

%%
tArray = [];
dt = 0.2;
tf = 150;
iter = 0;
for t = 0 : dt : tf
    iter = iter + 1;
    tArray = [tArray t];
end
tArray = tArray';
load('RangeAccel_vs_Time.mat')
Time_RA = RangeAccel_vs_Time(:,1); %Time sec
RangeAccel = RangeAccel_vs_Time(:,2); %RangeAccel ft/sec
clear indx
for ii = 1: length(tArray)
    indx1 = find(tArray(ii) >= Time_RA);
    indx2 = find(Time_RA >= tArray(ii));
    try
        indx = find(indx1 == indx2);
    end
    try
        indx = indx1(end);
        Time_RA_trunc(ii,:) = Time_RA(indx,:);
        RangeAccel_trunc(ii,:) = RangeAccel(indx,:);
    end
end
iter = 0;
for ii = 1: length(Time_RA_trunc)
    if Time_RA_trunc(ii) == 0
        iter = iter + 1;
    end
end
for ii = 1: length(Time_RA_trunc)-iter
    Time_RA_trunc(ii) = Time_RA_trunc(iter+ii);
    RangeAccel_trunc(ii) = RangeAccel_trunc(iter+ii);
end
figure; plot(Time_RA_trunc,RangeAccel_trunc,'r-'); ...hold on; 
    figure; plot(Time_RA,RangeAccel,'b-')

%%

save('Time_T_trunc.mat', 'Time_T_trunc');
save('Thrust_trunc.mat', 'Thrust_trunc');

save('Time_P_trunc.mat', 'Time_P_trunc');
save('Pitch_trunc.mat', 'Pitch_trunc');

save('Time_A_trunc.mat', 'Time_A_trunc');
save('Alt_trunc.mat', 'Alt_trunc');

save('Time_AR_trunc.mat', 'Time_AR_trunc');
save('AltRate_trunc.mat', 'AltRate_trunc');

save('Time_AA_trunc.mat', 'Time_AA_trunc');
save('AltAccel_trunc.mat', 'AltAccel_trunc');

save('Time_R_trunc.mat', 'Time_R_trunc');
save('Range_trunc.mat', 'Range_trunc');

save('Time_RA_trunc.mat', 'Time_RA_trunc');
save('RangeAccel_trunc.mat', 'RangeAccel_trunc');

