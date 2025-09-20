close all; clear;

filename = '40min.mat';
plot_particles(filename);

% filename = '100min.mat';
% plot_particles(filename);

filename = '190min.mat';
plot_particles(filename);

% filename = '350min.mat';
% plot_particles(filename);

filename = '1000min.mat';
plot_particles(filename);

function plot_particles(filename)
load(filename);
SPEN_size = 30;
figure('units','pixels','position',[0 0 600 600])

% Plot distribution of bound Xist with multiple SPEN
RGB = orderedcolors("glow");
green = RGB(5,:);

hold on
% Xist with 1–10 SPEN
scatter3(XS_x,   XS_y,   XS_z,   SPEN_size,  green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS2_x,  XS2_y,  XS2_z,  SPEN_size*2,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS3_x,  XS3_y,  XS3_z,  SPEN_size*3,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS4_x,  XS4_y,  XS4_z,  SPEN_size*4,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS5_x,  XS5_y,  XS5_z,  SPEN_size*5,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS6_x,  XS6_y,  XS6_z,  SPEN_size*6,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS7_x,  XS7_y,  XS7_z,  SPEN_size*7,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS8_x,  XS8_y,  XS8_z,  SPEN_size*8,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS9_x,  XS9_y,  XS9_z,  SPEN_size*9,green, 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(XS10_x, XS10_y, XS10_z, SPEN_size*10,   green, 'filled', 'MarkerFaceAlpha', 0.8);

% Total and bound Xist
scatter3(x_total, y_total, z_total, 5, 'k', 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(Xistb_x, Xistb_y, Xistb_z, SPEN_size, 'r', 'filled', 'MarkerFaceAlpha', 0.8);

% Plot settings
ylim([-8, 8]); xlim([-8 8]); zlim([-8 8]);
xticks([]); yticks([]); zticks([]);
xticklabels([]); yticklabels([]); zticklabels([]);
view(0, 0);

end

