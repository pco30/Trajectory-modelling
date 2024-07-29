function EarthProjection(apogee,v,dt)
% Earth projection uses the data from ShootingMethod to plot 2D animated
% representations of the data, and 3D representations of projectile motion
% on earth as seen from space, accounting for Earths' rotation.

r=6.3878*10^6; % Radius of Earth

[theta,z,t] = ShootingMethod(apogee,v,dt); % Obtaining data

% Setting up figure limits, colours, labels and Earth plot
figure(1)
hold on
set(gca,'Color','b')
set(gcf, 'Position', get(0, 'Screensize'));
plot(r*cos(0:2*pi/1000:2*pi),r*sin(0:2*pi/1000:2*pi)-r,'k')
xlim([z(1,1) z(1,end-1)+10^4])
ylim([z(3,end-1)-10^4 max(z(3,:))+10^4])
xlabel('X-Displacement (m)')
ylabel('Y-Displacement (m)')
fill(r*cos(0:2*pi/1000:2*pi),r*sin(0:2*pi/1000:2*pi)-r,'g')


% Animated plot to represent velocity of the projectile (reccomend
% selecting stepsize as 0.1 for this, but 0.04 for accuracy)
h = animatedline(Color='r',LineWidth=1.5);
for k = 1:length(z)
    addpoints(h,z(1,k),z(3,k));
    drawnow
end

% Locating and labelling apogee
[a,b]=max(z(5,:));
text(z(1,b)+10^4,z(3,b)+3000,'Apogee','Color','w')
plot(z(1,b),z(3,b),'p','color','w')

% Locating and labelling parachute release point
d = z(5,:)<10000;
e = diff(d(1,:));
[f,g]=max(e);
text(z(1,g)-5*10^4,z(3,g),'Parachute Release','Color','w')
plot(z(1,g),z(3,g),'p','color','w')

% Locating and labelling impact point
text(z(1,end)-2*10^4,z(3,end)+5*10^3,'Impact','Color','w')
plot(z(1,end),z(3,end),'p','color','w')

% Requested variables displayed in textbox
annotation('textbox', [0.3,0.1,0.1,0.1], 'BackgroundColor','w','String',{"Angle of elevation: "+theta(end)+" °", "Time of flight: "+t(end)+" s", "Distance travelled over ground : "+z(1,end)/1000+" Km","Impact velocity: "+sqrt(z(2,end)^2+z(4,end)^2)+" m/s"})

legend('','Earth Surface','Projectile Path',color='w')
hold off

%% 3D-projectile motion with earth rotation taken into account

omega = 7.2921150*10^-5; % Angular velocity of earth 
v = omega * z(1,:); % Velocity of earth at a given radius from the z-axis

%Euler method to find distance earth rotates
z(6,1) = 0;
for n=1:length(v)-1
z(6,n+1) = z(6,n) + v(n)*0.1; % Displacement of projectile due to rotation of earth
end

%Computing points for surface of a sphere for Earth
[x,y,zz] = sphere(120);
x = x*r;
y = y*r;
zz = zz*r;

%Setting up figures
figure('color','k');
set(gcf, 'Position', get(0, 'Screensize'));
tiledlayout(1,2)

nexttile

% Plotting Earth with projectile trajectory
globe=surf(x,y,zz-r); % Plotting sphere
title('Earth View','color','w')
hold on
axis equal
image = 'https://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg'; % See link for original photo
earthtexture = imread(image);
set(globe, 'FaceColor', 'texturemap', 'CData', earthtexture,'EdgeColor', 'none') % Applying image to sphere and removing lines
plot3(z(1,:),z(6,:),z(3,:),'r',LineWidth=1) % Projectile path in 3D
set(gca,'Visible','off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
hold off


nexttile
% Same plot above but zoomed in
globe1=surf(x,y,zz-r);
hold on
title('Close up view','color','w')
set(globe1, 'FaceColor', 'texturemap', 'CData', earthtexture,'EdgeColor', 'none')
plot3(z(1,:),z(6,:),z(3,:),'r',LineWidth=1)
set(gca, 'Visible','off');
set(findall(gca, 'type', 'text'), 'visible', 'on')
axis equal
zlim([-0.1*r 0.1*r])
xlim([-0.2*r 0.2*r])
ylim([-0.2*r 0.2*r])

hold off