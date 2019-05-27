%Initializing the agents to random positions with barrier certificates 
%and data plotting.  This script shows how to initialize robots to a
%particular point
%Paul Glotfelter 
%3/24/2016
clear all, close all, clc

N_cycle = 5;
N_patrol = 3;
N_enemy = 1;
N = N_cycle + N_patrol + N_enemy;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);

videoFLag = 0;

% Initialize x so that we don't run into problems later.  This isn't always
% necessary
x = r.get_poses();
r.step();

if videoFLag 
    vid = VideoWriter('RoboBrigade_Defense.mp4', 'MPEG-4');
    vid.Quality = 100;
    vid.FrameRate = 72;
    open(vid);
    writeVideo(vid, getframe(gcf));
end

% Create a barrier certificate so that the robots don't collide
si_barrier_certificate = create_si_barrier_certificate('SafetyRadius', 0.15);
si_to_uni_dynamics = create_si_to_uni_mapping2('LinearVelocityGain', 0.25, 'AngularVelocityLimit', 7.5);
        
% Variables Needed
circle_radius = 0.2;
circle_center = [-1.1;0];

circular_ic = [circle_radius*cos(0:2*pi/N_cycle:2*pi*(1-1/N_cycle)) + circle_center(1)...
    * ones(1,N_cycle) ; circle_radius*sin(0:2*pi/N_cycle:2*pi*(1- 1/N_cycle)) + ...
    circle_center(2)*ones(1,N_cycle)];
patrol_ic = [-0.4,-0.1;0,-0.1;-0.2,0]';
enemy_ic = [1.4,0]';
initial_conditions = [circular_ic patrol_ic enemy_ic];
initial_conditions = [initial_conditions ; zeros(1,N)];

%% Initial Position Set-up
args = {'PositionError', 0.1, 'RotationError', 50};
init_checker = create_is_initialized(args{:});
controller = create_si_position_controller();

while(~init_checker(x, initial_conditions))
    x = r.get_poses();
    dxi = controller(x(1:2, :), initial_conditions(1:2, :));
    for i = 1 : N
        if (norm(dxi(:,i))>r.max_linear_velocity)
            dxi(:,i) = r.max_linear_velocity*(dxi(:,i)./norm(dxi(:,i)));
        end
    end
    
    dxi = si_barrier_certificate(dxi, x(1:2, :));
    dxu = si_to_uni_dynamics(dxi, x);
    
    r.set_velocities(1:N, dxu);
    r.step();
end

%% Grab tools we need to convert from single-integrator to unicycle dynamics

% Single-integrator -> unicycle dynamics mapping
si_to_uni_dyn = create_si_to_uni_mapping2('LinearVelocityGain', 0.7, 'AngularVelocityLimit', 0.32*r.max_angular_velocity);
% Single-integrator barrier certificates
si_barrier_cert = create_si_barrier_certificate('SafetyRadius', 0.15);
% Single-integrator position controller
si_pos_controller = create_si_position_controller();

%% Parameters for Iteration & Algorithms

iterations = 3000;
dxi = zeros(2, N);
formation_control_gain = 15;

% waypoints for enemy movement and initial position of the patrol team
waypoints_enemy = [1.35,0.8;0.5,0.8;0.5,-0.8;1.35,-0.8]';
close_enough = 0.05;
enemy_state = 1;

%% Graph topology for patrol formation
A = [0,1,1;...
    1,0,1;...
    1,1,0];
delta = [2,0,0;
    0,2,0;
    0,0,2];
L = delta-A;
weights = [0,0.32,0.16;
    0.32,0,0.16;
    0.16,0.16,0];

%% Algorithm

for t = 1:iterations
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();

    % dynamics for circlic pursuit around the flag for robot 1-5
    for i = 1:N_cycle
        dxi(:, i) = [0 ; 0];
        error_radius = circle_radius - norm(x(1:2,i) - circle_center);
        theta = pi/N_cycle + error_radius;
        R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
        
        if i == N_cycle
            dxi(:, i) = 0.2*R*(x(1:2,1) - x(1:2,N_cycle));
        else
            dxi(:, i) = 0.2*R*(x(1:2,i+1) - x(1:2,i));
        end 
        
        if (norm(dxi(:,i))>r.max_linear_velocity)
            dxi(:,i) = r.max_linear_velocity*(dxi(:,i)./norm(dxi(:,i)));
        end
    end
    
    % dynamic for the patrol team
    for i = N_cycle+1 : N_cycle + N_patrol
        dxi(:, i) = [0 ; 0];
        neighbors = topological_neighbors(L, i-N_cycle);
        if (norm(x(1:2, N_cycle + N_patrol) - x(1:2, i))^2 < weights(i-N_cycle,N_patrol))
            dxi(:,i) = 0;
        else
            for j = neighbors
                if (i == N_cycle + N_patrol)
                    dxi(:, i) = dxi(:, i) + 0.25*formation_control_gain*...
                    (norm(x(1:2, j+N_cycle) - x(1:2, i))^2 - weights(i-N_cycle,j)^2)*...
                    (x(1:2, j+N_cycle) - x(1:2, i)) + ...
                    2*((circle_center - x(1:2,i)) + 1.3*(x(1:2,N) - x(1:2,i)));
                else
                    dxi(:, i) = dxi(:, i) + ...
                    formation_control_gain*2*(norm(x(1:2, j+N_cycle) - ...
                    x(1:2, i))^2 - weights(i-N_cycle,j)^2)*(x(1:2, j+N_cycle) - ...
                    x(1:2, i));
                end
            end  
        end
        
        if (norm(dxi(:,i))>r.max_linear_velocity)
            dxi(:,i) = r.max_linear_velocity*(dxi(:,i)./norm(dxi(:,i)));
        end
    end
    
    % dynamics for the enemy (robot 9)
    waypoint_enemy = waypoints_enemy(:,enemy_state);
    enemMove = si_pos_controller(x(1:2, N), waypoint_enemy);
    dxi(:,N) = 0.2*enemMove/norm(enemMove);
    if(norm(x(1:2, N) - waypoint_enemy) < close_enough)
        enemy_state = enemy_state + 1;
        if (enemy_state == 5)
            enemy_state = 1;
        end
    end
    if (norm(dxi(:,N))>r.max_linear_velocity)
            dxi(:,N) = r.max_linear_velocity*(dxi(:,N)./norm(dxi(:,N)));
    end
    
    %% Use barrier certificate and convert to unicycle dynamics
    dxi = si_barrier_cert(dxi, x);
    dxu = si_to_uni_dyn(dxi, x);
    
    %% Send velocities to agents
    
    %Set velocities
    r.set_velocities(1:N, dxu);
    
    %Iterate experiment
    r.step();
    if videoFLag && mod(t,10)                               % Record a video frame every 10 iterations
            writeVideo(vid, getframe(gcf)); 
    end
end

if videoFLag; close(vid); end

% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();