
s = [0,0,0]; N = 100;
x = nan(1,N); y = x; z = x;
x(1) = s(1);
y(1) = s(2);
z(1) = s(3);
step = 1;
i = 2;
while i <= N
    theta = 2*pi*rand();
    phi = pi*rand();
    xi = x(i-1) + step*sin(phi)*cos(theta);
    yi = y(i-1) + step*sin(phi)*sin(theta);
    zi = z(i-1) + step*cos(phi);
    if vecnorm([xi,yi,zi]-s) <= 2.5 - 0.01
        x(i) = xi;
        y(i) = yi;
        z(i) = zi;
        i = i+1;
    end
end
plot3(x,y,z,'.--','color',[0.7,0.7,0.7],'MarkerSize',15,'markerfacecolor','k','markeredgecolor','k','LineWidth',0.3);
axis square; set(gca,'fontsize',35)
xticklabels('x (\mum)'); yticklabels('y (\mum)'); zticklabels('z (\mum)'); grid on;
