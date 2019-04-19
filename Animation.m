function Animation(Y_ref)
figure
Obstacle1_x = [1,3,3,1,1,3,1,3];
Obstacle1_y = [-1,-1,1,1,-1,1,1,-1] * 2.25;
Obstacle2_x = [1,3,3,1,1,3,1,3] * -1;
vehicle_shape = [4.2, 1.6]; % length x width
for t = 1:10:length(Y_ref)
    plot(Obstacle1_x, Obstacle1_y,'k');hold on;
    plot(Obstacle2_x, Obstacle1_y,'k');hold on;
    x = Y_ref(1,t);
    y = Y_ref(2,t);
    phi = Y_ref(3,t);
    l = vehicle_shape(1) / 2;
    w = vehicle_shape(2) / 2;
    vehiclex = [x + l*cos(phi) + w * sin(phi), x + l * cos(phi) - w * sin(phi),  x - l * cos(phi) - w * sin(phi),  x - l * cos(phi) + w * sin(phi),  x + l * cos(phi) + w * sin(phi)];
    vehicley = [y + l*sin(phi) - w *cos(phi),y + l*sin(phi) + w *cos(phi),y - l*sin(phi) + w *cos(phi),y - l*sin(phi) - w *cos(phi),y + l*sin(phi) - w *cos(phi)];
    plot(vehiclex, vehicley,'b');
    drawnow 
end

