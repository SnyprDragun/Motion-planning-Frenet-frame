%% MOTION PLANNING IN FRENET FRAME WITH TRAILER KINEMATICS
clear; close all; clc;

dt = 0.05; % time step
T = 50;    % simulation time

%% 1) Figure-of-eight waypoints in Cartesian frame
t = linspace(0,2*pi,500);
a = 10; b = 5; % scaling
x_ref = a * sin(t);
y_ref = b * sin(t).*cos(t);  % Lissajous curve
waypoints = [x_ref' y_ref'];

figure; plot(x_ref,y_ref,'b'); axis equal; grid on;
title('Figure-of-Eight in Cartesian Frame');

%% Transform to Frenet frame
[s, d] = cartesian2frenet(waypoints);
figure; plot(s,d,'r'); grid on;
title('Figure-of-Eight in Frenet Frame');
xlabel('s (along road)'); ylabel('d (lateral offset)');

%% 2) Track figure-of-eight in Frenet frame
% Robot parameters
x = [0; -1; 0]; % [x; y; theta]
v = 2;          % forward velocity
L = 2;          % wheelbase
k = 0.5;        % pure pursuit gain
traj = [];

for i=1:length(s)
    % Target point
    target = waypoints(i,:)';
    alpha = atan2(target(2)-x(2), target(1)-x(1)) - x(3);
    omega = k * alpha;
    
    % Unicycle model
    x(1) = x(1) + dt*v*cos(x(3));
    x(2) = x(2) + dt*v*sin(x(3));
    x(3) = x(3) + dt*omega;
    
    traj = [traj; x'];
end

figure; plot(waypoints(:,1),waypoints(:,2),'k--'); hold on;
plot(traj(:,1),traj(:,2),'b'); axis equal; grid on;
title('Tracking Figure-of-Eight in Frenet Frame');

%% 3a) Introduce obstacle and avoid
obs = [2 0]; r_obs = 1; % obstacle location and radius
% Simple obstacle avoidance by adjusting target d when close
traj_obs = [];
x = [0; -1; 0];

for i=1:length(s)
    target = waypoints(i,:)';
    dist_obs = norm(x(1:2)-obs');
    if dist_obs < 3
        target(2) = target(2) + 2; % shift lateral offset
    end
    
    alpha = atan2(target(2)-x(2), target(1)-x(1)) - x(3);
    omega = k * alpha;
    
    x(1) = x(1) + dt*v*cos(x(3));
    x(2) = x(2) + dt*v*sin(x(3));
    x(3) = x(3) + dt*omega;
    
    traj_obs = [traj_obs; x'];
end

figure; plot(waypoints(:,1),waypoints(:,2),'k--'); hold on;
viscircles(obs,r_obs); plot(traj_obs(:,1),traj_obs(:,2),'r'); axis equal; grid on;
title('Obstacle Avoidance in Frenet Frame');

%% 3b) Tracking with single trailer
Lh = 1; D = 2; % hitch and trailer distances
x = [0; -1; 0]; phi = 0; % trailer hitch angle
traj_trailer = [];

for i=1:length(s)
    target = waypoints(i,:)';
    alpha = atan2(target(2)-x(2), target(1)-x(1)) - x(3);
    omega = k * alpha;
    
    % Robot dynamics
    x(1) = x(1) + dt*v*cos(x(3));
    x(2) = x(2) + dt*v*sin(x(3));
    x(3) = x(3) + dt*omega;
    
    % Trailer dynamics
    vt = v*cos(x(3)-phi) + Lh*omega*sin(x(3)-phi);
    phi = phi + dt*((v*sin(x(3)-phi) - Lh*omega*cos(x(3)-phi))/D);
    
    traj_trailer = [traj_trailer; x' phi];
end

figure; plot(waypoints(:,1),waypoints(:,2),'k--'); hold on;
plot(traj_trailer(:,1),traj_trailer(:,2),'g'); axis equal; grid on;
title('Figure-of-Eight Tracking with Trailer');

%% 4) Reverse tracking with trailer (MPC placeholder)
% Here we assume you will implement MPC or LQR separately.
% For demonstration, we simulate reverse driving by negating v.
v = -1;
x = [0; -1; 0]; phi = 0;
traj_rev = [];
for i=1:length(s)
    target = waypoints(i,:)';
    alpha = atan2(target(2)-x(2), target(1)-x(1)) - x(3);
    omega = k * alpha;
    
    % Robot reverse dynamics
    x(1) = x(1) + dt*v*cos(x(3));
    x(2) = x(2) + dt*v*sin(x(3));
    x(3) = x(3) + dt*omega;
    
    % Trailer dynamics
    vt = v*cos(x(3)-phi) + Lh*omega*sin(x(3)-phi);
    phi = phi + dt*((v*sin(x(3)-phi) - Lh*omega*cos(x(3)-phi))/D);
    
    traj_rev = [traj_rev; x' phi];
end

figure; plot(waypoints(:,1),waypoints(:,2),'k--'); hold on;
plot(traj_rev(:,1),traj_rev(:,2),'m'); axis equal; grid on;
title('Reverse Tracking with Trailer');

%% --- Supporting Function ---
function [s,d] = cartesian2frenet(waypoints)
    s = [0];
    d = [0];
    for i=2:size(waypoints,1)
        ds = norm(waypoints(i,:)-waypoints(i-1,:));
        s(end+1) = s(end)+ds;
        d(end+1) = 0; % assuming no lateral offset from path centerline
    end
end
