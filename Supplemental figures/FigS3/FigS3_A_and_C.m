filename = '1X_540min.mat';
plot_particles(filename);

filename = '1X_1000min.mat';
plot_particles(filename);

filename = '1X_5000min.mat';
plot_particles(filename);

filename = '2X_400min.mat';
plot_particles(filename);

filename = '2X_540min.mat';
plot_particles(filename);

filename = '2X_1000min.mat';
plot_particles(filename);

function plot_particles(filename)
load(filename);
SPEN_size = 10;
figure('units','pixels','position',[0 0 600 600])

% Plot distribution of bound Xist with multiple SPEN

% Xist with 1 SPEN
plot3(XS_x,XS_y,XS_z,'g.','MarkerSize',SPEN_size); hold on;
% Xist with 2 SPEN
plot3(XS2_x,XS2_y,XS2_z,'g.','MarkerSize',SPEN_size*2);
% Xist with 3 SPEN
plot3(XS3_x,XS3_y,XS3_z,'g.','MarkerSize',SPEN_size*3);
% Xist with 4 SPEN ...
plot3(XS4_x,XS4_y,XS4_z,'g.','MarkerSize',SPEN_size*4);

plot3(XS5_x,XS5_y,XS5_z,'g.','MarkerSize',SPEN_size*5); 
plot3(XS6_x,XS6_y,XS6_z,'g.','MarkerSize',SPEN_size*6); 
plot3(XS7_x,XS7_y,XS7_z,'g.','MarkerSize',SPEN_size*7); 
plot3(XS8_x,XS8_y,XS8_z,'g.','MarkerSize',SPEN_size*8); 
plot3(XS9_x,XS9_y,XS9_z,'g.','MarkerSize',SPEN_size*9); 

% Xist with 10 SPEN
plot3(XS10_x,XS10_y,XS10_z,'g.','MarkerSize',SPEN_size*10); 

plot3(x_total,y_total,z_total,'k.','MarkerSize',5);
plot3(Xistb_x,Xistb_y,Xistb_z,'r.','MarkerSize',SPEN_size);
ylim([-8,8]); xlim([-8 8]); zlim([-8 8]); xticks([]);yticks([]);zticks([]);
xticklabels([]);yticklabels([]);zticklabels([]);
view(0,0);

end

