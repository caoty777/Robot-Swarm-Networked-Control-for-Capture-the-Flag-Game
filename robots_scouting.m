%Initializing the agents to random positions with barrier certificates
%and data plotting.  This script shows how to initialize robots to a
%particular point
%Paul Glotfelter
%3/24/2016

N = 3;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);

videoFLag = 0;

% Initialize x so that we don't run into problems later.  This isn't always
% necessary

x = r.get_poses();
r.step();

if videoFLag 
    vid = VideoWriter('RoboBrigade_Recon.mp4', 'MPEG-4');
    vid.Quality = 100;
    vid.FrameRate = 72;
    open(vid);
    writeVideo(vid, getframe(gcf));
end

% Create a barrier certificate so that the robots don't collide
si_barrier_certificate = create_si_barrier_certificate('SafetyRadius', 0.22);
si_to_uni_dynamics = create_si_to_uni_mapping2();

%Get randomized initial conditions in the robotarium arena
initial_conditions = generate_initial_conditions(N, 'Width', 1, 'Height',1, 'Spacing', 0.2);


% We'll make the rotation error huge so that the initialization checker
% doesn't care about it
args = {'PositionError', 0.1, 'RotationError', 50};
init_checker = create_is_initialized(args{:});
controller = create_si_position_controller();

% while(~init_checker(x, initial_conditions))
%     
%     x = r.get_poses();
%     dxi = controller(x(1:2, :), initial_conditions(1:2, :));
%     dxi = si_barrier_certificate(dxi, x(1:2, :));
%     dxu = si_to_uni_dynamics(dxi, x);
%     
%     r.set_velocities(1:N, dxu);
%     r.step();
% end\

%% Create the desired Laplacian
iterations = 6600;
l = .53;
p = 1.5;
q = .9183681179;
%Graph laplacian

L = zeros(N, N);
A = [0,1,1;...
    1,0,1;...
    1,1,0];
delta = [2,0,0;
    0,2,0;
    0,0,2];
L = delta-A;
weights = [0,0.3,0.16;
    0.3,0,0.16;
    0.16,0.16,0];


%Initialize velocity vector
dxi = zeros(2, N);

%State for leader
state = 1;

% These are gains for our formation control algorithm
formation_control_gain = 15;



%% Grab tools we need to convert from single-integrator to unicycle dynamics

% Single-integrator -> unicycle dynamics mapping
si_to_uni_dyn = create_si_to_uni_mapping2('LinearVelocityGain', 0.7, 'AngularVelocityLimit', 0.2*r.max_angular_velocity);
% Single-integrator barrier certificates
si_barrier_cert = create_si_barrier_certificate('SafetyRadius', 1.5*r.robot_diameter);
% Single-integrator position controller
si_pos_controller = create_si_position_controller();

% waypoints = [1 0.5; -1 0.5; -1 -0.5; 1 -0.5]';
waypoints = [0.7 0.6; -.7 .6; -1.3 0; -.7 -.6; .7 -.6; 1.3 0]';
close_enough = 0.08;
state = 1;
for t = 1:iterations
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();
    
    %% Algorithm
    waypoint = waypoints(:, state);
    for i = 1:N
        %Zero velocity and get the topological neighbors of agent i
        dxi(:, i) = [0 ; 0];
        
        neighbors = topological_neighbors(L, i);
        
        for j = neighbors
            difsq = norm(x(1:2, j) - ...
                x(1:2, i))^2 - weights(i,j)^2;
            if (norm(x(1:2, 3) - x(1:2, i))^2 < weights(i,3))
                dxi(:,i) = 0;
            elseif (i==3)
                wptVec = si_pos_controller(x(1:2, 3), waypoint);
                wptDir = wptVec/norm(wptVec);
                dxi(:, i) = dxi(:, i) + ...
                    0.5*formation_control_gain*...
                    norm(wptVec)/(N-1)*...
                    (difsq)*(x(1:2, j) - ...
                    x(1:2, i)) + ...
                    0.4*wptDir/(N-1);
            else
                dxi(:, i) = dxi(:, i) + ...
                    formation_control_gain*(difsq)*(x(1:2, j) - ...
                    x(1:2, i));
            end
            
        end
        if (norm(dxi(:,i))>r.max_linear_velocity)
            dxi(:,i) = r.max_linear_velocity*(dxi(:,i)./norm(dxi(:,i)));
        end
       
    end
    if(norm(x(1:2, 3) - waypoint) < close_enough)
        state = state+1;
        if (state == 7)
            state=1;
        end
    end
    %% Make the leader travel between waypoints
    
%    waypoint = waypoints(:, state);
%     
%     switch state        
%         case 1
%             dxi(:, 3) = si_pos_controller(x(1:2, 3), waypoint);
%             if(norm(x(1:2, 3) - waypoint) < close_enough)
%                 state = 2;
%             end
%         case 2
%             dxi(:, 3) = si_pos_controller(x(1:2, 3), waypoint);
%             if(norm(x(1:2, 3) - waypoint) < close_enough)
% %                 state = 3;
%                   state = 1;
%             end
% %         case 3
% %             dxi(:, 3) = si_pos_controller(x(1:2, 3), waypoint);
% %             if(norm(x(1:2, 3) - waypoint) < close_enough)
% %                 state = 4;
% %             end
% %         case 4
% %             dxi(:, 3) = si_pos_controller(x(1:2, 3), waypoint);
% %             if(norm(x(1:2, 3) - waypoint) < close_enough)
% %                 state = 1;
% %             end
%     end
    
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