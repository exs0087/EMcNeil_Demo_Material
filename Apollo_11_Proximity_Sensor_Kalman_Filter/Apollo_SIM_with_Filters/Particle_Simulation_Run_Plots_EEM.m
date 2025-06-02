close all;

t_out = 0.1:0.1:length(Attitude)*0.1;

p1 = figure(1);
plot(t_out,Particle_Range,'c-o');
hold on
plot(t_out,Sim_Range,'k-');
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Range (ft)');xlabel('Time (sec)');
legend('Particle Range','Sim Range','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);
grid on;

p2 = figure(2);
plot(t_out,y_meas_2ws,'mo'); 
hold on;
plot(t_out,Particle_Alt,'c-o');
plot(t_out,Sim_Alt,'k-');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');xlabel('Time (sec)');
legend('Meas Alt','Particle Alt','Sim Alt','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);

figure(3); 
plot(t_out,y_meas_2ws,'mo'); 
hold on;
plot(t_out,Particle_Alt,'co');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');xlabel('Time (sec)');
legend('Meas Alt','Particle Alt','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);

figure(4);
plot(t_out,Particle_Alt,'c-o');
hold on;
plot(t_out,Sim_Alt,'ko');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Est Altitude (ft)');xlabel('Time (sec)');
legend('Particle Alt','Sim Alt','Location','SouthEast');
xlim([0,length(Attitude)*0.1]);

p5 = figure(5); 
subplot(3,1,1)
plot(t_out,(LatLongAlt_tows(:,1).*99486),'r-');
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Latitude (ft)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

subplot(3,1,2);
plot(t_out,LatLongAlt_tows(:,2)*99486,'b-');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Longitude (ft)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

subplot(3,1,3);
plot(t_out,LatLongAlt_tows(:,3),'k-');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Altitude (ft)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

p6 = figure(6);
plot(t_out,FXYZ);
legend('Fx','Fy','Fz');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('FX,FY,FZ (lbf)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

p7 = figure(7); 
plot(t_out,Attitude(:,2)*(180/pi),'m');
grid on;
title('Particle Filter');
set(gca,'FontSize',12); set(gcf,'Color','White');
ylabel('Pitch (deg)');xlabel('Time (sec)');
xlim([0,length(Attitude)*0.1]);

saveas(p1,'Particle_Range_state_estimates.png');
saveas(p2,'Particle_Altitude_state_estimates.png');
saveas(p5,'Particle_LatLonAlt.png');
saveas(p6,'Particle_Fxyz.png');
saveas(p7,'Particle_Pitch_Angle.png');