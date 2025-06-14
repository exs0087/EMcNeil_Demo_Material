close all;

t_out = 0.1:0.1:length(Attitude)*0.1;

p1 = figure(1);
plot(t_out,UKF_Range,'c-o');
hold on
plot(t_out,Sim_Range,'k-');
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Range (ft)');xlabel('Time (sec)');
legend('UKF Range','Sim Range','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);
grid on;

p2 = figure(2);
plot(t_out,y_meas_2ws,'mo'); 
hold on;
plot(t_out,UKF_Alt,'c-o');
plot(t_out,Sim_Alt,'k-');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');xlabel('Time (sec)');
legend('Meas Alt','UKF Alt','Sim Alt','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);

figure(3); 
plot(t_out,y_meas_2ws,'mo'); 
hold on;
plot(t_out,UKF_Alt,'co');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');xlabel('Time (sec)');
legend('Meas Alt','UKF Alt','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);

figure(4);
plot(t_out,UKF_Alt,'c-o');
hold on;
plot(t_out,Sim_Alt,'ko');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');xlabel('Time (sec)');
legend('UKF Alt','Sim Alt','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);

p5 = figure(5); 
subplot(3,1,1)
plot(t_out,(LatLongAlt_tows(:,1).*99486),'r-');
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Latitude (ft)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

subplot(3,1,2);
plot(t_out,LatLongAlt_tows(:,2)*99486,'b-');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Longitude (ft)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

subplot(3,1,3);
plot(t_out,LatLongAlt_tows(:,3),'k-');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude (ft)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

p6 = figure(6);
plot(t_out,FXYZ);
legend('Fx','Fy','Fz');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('FX,FY,FZ (lbf)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

p7 = figure(7); 
plot(t_out,Attitude(:,2)*(180/pi),'m');
grid on;
title('Unscented Kalman Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Pitch (deg)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

saveas(p1,'UKF_Range_state_estimates.png');
saveas(p2,'UKF_Altitude_state_estimates.png');
saveas(p5,'UKF_LatLonAlt.png');
saveas(p6,'UKF_Fxyz.png');
saveas(p7,'UKF_Pitch_Angle.png');