classdef Joint < matlab.mixin.SetGet
    %JOINT An object of a body attached to a proximal body via a hinge
    %joint
    properties
        %Static values
        joint_name;       % Name of the joint
        file_name;        % Name of the associated simulation file        
        m;                % Mass of the distal segment
        
        p_joint;          % Location of joint, in distal frame
        euler_ang_joint;  % Orientation of the joint, in proximal frame
        p_dis_frame;      % Location of the distal frame, in proximal frame
        euler_ang_dis;    % Orientation of the distal segment, in proximal frame
        
        axis1;            % Length of x axis of distal segment
        axis2;            % Length of y axis of distal segment
        axis3;            % Length of z axis of distal segment
        
        theta;            % Rotation of the joint
        velocity;         % Rate of rotation of the joint
        stiffness;        % Stiffness of the joint
        max_torque;       % Maximum torque expected on the joint
        num_pts;          % Number of positions (also equal to the number of stiffnesses) to train on
        to_plot;          % Boolean of whether to plot the residuals or not
        
        dt_motion;        % Time 
        theta_motion;     % Joint anlge
        theta_dot_motion; % Joint velocity
        
        joint_motion;     % Position of Joint during walking
        body_motion;      % Position of distal body during walking
        
    end
    
    methods
        function obj = Joint(p_joint,euler_ang_joint,p_body,euler_ang_body,theta_range,theta_offset,velocity_range,max_torque,joint_name,file_name,to_plot)
            obj.joint_name = joint_name;
            obj.file_name = file_name;
            obj.to_plot = to_plot;
            
            obj.p_joint = p_joint;
            obj.euler_ang_joint = euler_ang_joint;
            obj.p_dis_frame = p_body;
            obj.euler_ang_dis = euler_ang_body;
            
            obj.num_pts = 10;
            obj.theta = linspace(theta_range(1),theta_range(2),obj.num_pts)+theta_offset;
            obj.velocity = linspace(velocity_range(1),velocity_range(2),obj.num_pts+1);
            
            %Stiffness is the desired torque divided by the angular
            %displacement requires to obtain that torque. If we assume we
            %can generate the user input torque with a deflection of 10
            %degrees, that equals about 0.15 radian.
            if isempty(max_torque)
                max_torque = 1;
            end
            obj.max_torque = abs(max_torque);
            %For a reasonable deflection of 10 degrees
            obj.stiffness = linspace(0,max_torque/.175,obj.num_pts);
            
        end
%% Calculate joint axis rotation matrix based on angles
        function result = axis_angle_rotation(obj,i)
            c = cos(obj.theta_motion(i));
            s = sin(obj.theta_motion(i));
            
            a1 = obj.axis1;
            a2 = obj.axis2;
            a3 = obj.axis3;
            
            result = [c+a1^2*(1-c), a1*a2*(1-c)-a3*s, a1*a3*(1-c)+a2*s;...
                 a1*a2*(1-c)+a3*s, c+a2^2*(1-c), a2*a3*(1-c)-a1*s;...
                 a1*a3*(1-c)-a2*s, a2*a3*(1-c)+a1*s, c+a3^2*(1-c)];
        end
    end
    
    
end