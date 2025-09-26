function [x,y,z] = random_walk_gene(filename,center)
% The Xist transcription site is at 60% of the chromosome length.
% Centered on the Xist transcription site, we applied random walk algorithm
% to 
N = 60; [x,y,z] = generate_coordinates(N,center); % Binding site 1-60
N = 101-60; [x2,y2,z2] = generate_coordinates(N,center); % Binding site 61-100

x2(1) = []; y2(1) = []; z2(1) = [];
x = [flip(x),x2]; y = [flip(y),y2]; z = [flip(z),z2];
save([filename '.mat'],"x","y","z")
end

% Random walk process
function [x,y,z] = generate_coordinates(N,center)
s = center;
x = nan(1,N); y = x; z = x; 

x(1) = s(1);
y(1) = s(2);
z(1) = s(3);
step = 1;0.5;
i = 2;
while i <= N
    theta = 2*pi*rand();
    phi = pi*rand();
    xi = x(i-1) + step*sin(phi)*cos(theta);
    yi = y(i-1) + step*sin(phi)*sin(theta);
    zi = z(i-1) + step*cos(phi);
    if vecnorm([xi,yi,zi]-center) <= 2.5
        x(i) = xi;
        y(i) = yi;
        z(i) = zi;
        v(i) = vecnorm([x(i),y(i),z(i)]-[x(i-1),y(i-1),z(i-1)]);
        i = i+1;
    end
end
end