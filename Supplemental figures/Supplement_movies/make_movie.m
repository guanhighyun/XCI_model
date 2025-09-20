function make_movie(filename,nframes,moviedir,moviename)
x_total = []; y_total = []; z_total = [];
SPEN_size = 15;
for j = 1:3 % Change to 1:2 if two chromosomes, change to 1 if one chromosome
    matname = sprintf('%s.cfg_chr%d.mat.mat',filename,j);
    load(matname);
    x_total = [x,x_total];
    y_total = [y,y_total];
    z_total = [z,z_total];
end
filename = sprintf('%s.txt',filename);
[t,positions]=read_molPos3(filename,nframes);
[~,positions]=read_molPos3_SPEN(sprintf('XS_%s',filename),nframes,"XS",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS2_%s',filename),nframes,"XS2",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS3_%s',filename),nframes,"XS3",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS4_%s',filename),nframes,"XS4",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS5_%s',filename),nframes,"XS5",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS6_%s',filename),nframes,"XS6",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS7_%s',filename),nframes,"XS7",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS8_%s',filename),nframes,"XS8",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS9_%s',filename),nframes,"XS9",positions);
[~,positions]=read_molPos3_SPEN(sprintf('XS10_%s',filename),nframes,"XS10",positions);
%[~,positions]=read_molPos3_SPEN(sprintf('Activator_%s',filename),nframes,"Activator",positions);

movObj = VideoWriter(sprintf('%s/%s.mov',moviedir,moviename),'MPEG-4');
movObj.FrameRate = 50;
movObj.Quality = 100;
    open(movObj)

figure('units','pixels','position',[0 0 400 400])
t(1) = 0;

nXistb1 = nan(1,nframes);
nXistb2 = nan(1,nframes);
nXistb3 = nan(1,nframes);

nSPENb1 = nan(1,nframes);
nSPENb2 = nan(1,nframes);
nSPENb3 = nan(1,nframes);

nActivator = nan(1,nframes);

RGB = orderedcolors("glow");
green = RGB(5,:);

for i= 1:10:nframes;
    if isnan(t(i))
       break;
    end
    cla;  
    [Xistb_x,Xistb_y,Xistb_z] = read_coordinates('Xistb',i,positions);
    [XS_x,XS_y,XS_z] = read_coordinates('XS',i,positions);
    [XS2_x,XS2_y,XS2_z] = read_coordinates('XS2',i,positions);
    [XS3_x,XS3_y,XS3_z] = read_coordinates('XS3',i,positions);
    [XS4_x,XS4_y,XS4_z] = read_coordinates('XS4',i,positions);
    [XS5_x,XS5_y,XS5_z] = read_coordinates('XS5',i,positions);
    [XS6_x,XS6_y,XS6_z] = read_coordinates('XS6',i,positions);
    [XS7_x,XS7_y,XS7_z] = read_coordinates('XS7',i,positions);
    [XS8_x,XS8_y,XS8_z] = read_coordinates('XS8',i,positions);
    [XS9_x,XS9_y,XS9_z] = read_coordinates('XS9',i,positions);
    [XS10_x,XS10_y,XS10_z] = read_coordinates('XS10',i,positions);

    scatter3(XS_x,   XS_y,   XS_z,   SPEN_size,  green, 'filled', 'MarkerFaceAlpha', 0.8); hold on;
    scatter3(XS2_x,  XS2_y,  XS2_z,  SPEN_size*2,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS3_x,  XS3_y,  XS3_z,  SPEN_size*3,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS4_x,  XS4_y,  XS4_z,  SPEN_size*4,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS5_x,  XS5_y,  XS5_z,  SPEN_size*5,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS6_x,  XS6_y,  XS6_z,  SPEN_size*6,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS7_x,  XS7_y,  XS7_z,  SPEN_size*7,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS8_x,  XS8_y,  XS8_z,  SPEN_size*8,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS9_x,  XS9_y,  XS9_z,  SPEN_size*9,green, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(XS10_x, XS10_y, XS10_z, SPEN_size*10,   green, 'filled', 'MarkerFaceAlpha', 0.8);

    scatter3(x_total, y_total, z_total, 5, 'k', 'filled', 'MarkerFaceAlpha', 0.8);
    scatter3(Xistb_x, Xistb_y, Xistb_z, SPEN_size, 'r', 'filled', 'MarkerFaceAlpha', 0.8);

    ylim([-10,10]); xlim([-10 10]); zlim([-10 10]); xticks([]);yticks([]);zticks([]);xticklabels([]);yticklabels([]);zticklabels([]);
    time = round(t(i)/60);
    title(sprintf('%i hours',time),'fontsize',25); view(0,0);
    view(0,0)
    drawnow;
    writeVideo(movObj,getframe(gcf));  
end
close(movObj);

end

function [x,y,z] = read_coordinates(name,i,positions)
if ~isempty(positions.(name){i}) 
    x = positions.(name){i}(:,1);
    y = positions.(name){i}(:,2);
    z = positions.(name){i}(:,3);
else
    x = nan;
    y = nan;    
    z = nan; 
end
end

function plot_species(t,positions,Name)
    for i = 1:find(~isnan(t),1,'last')
        if ~isempty(positions.(Name))
            n(i) = size(positions.(Name){i},1);
        end
    end; figure;plot(n)
    title(Name)
end