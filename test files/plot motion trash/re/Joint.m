classdef Joint < matlab.mixin.SetGet
    %JOINT An object of a body attached to a proximal body via a hinge
    %joint
    properties
        %Static values
        joint_name; %Name of the joint
        file_name; %Name of the associated simulation file
        p_joint_dis; %location of joint, in distal frame
        euler_ang_joint; %orientation of the joint, in proximal frame
        
        m; %mass of the distal segment
        axis1; %length of x axis of distal segment
        axis2; %length of y axis of distal segment
        axis3; %length of z axis of distal segment
        p_dis_frame; %location of the distal frame, in proximal frame
        euler_ang_dis; %orientation of the distal segment, in proximal frame
        
        theta; %rotation of the joint
        velocity; %rate of rotation of the joint
        stiffness; %stiffness of the joint
        max_torque; %maximum torque expected on the joint
        num_pts; %number of positions (also equal to the number of stiffnesses) to train on
        to_plot; %boolean of whether to plot the residuals or not
        
        theta_motion;
        theta_dot_motion;
        dt_motion;
        
        joint_motion; % position of Joint during walking
        body_motion; % position of distal body during walking
        
    end
    
    methods
        function obj = Joint(p_joint_dis,euler_ang_joint,p_dis_frame,euler_ang_dis_frame,theta_range,theta_offset,velocity_range,max_torque,joint_name,file_name,to_plot)
            obj.joint_name = joint_name;
            obj.file_name = file_name;
            obj.to_plot = to_plot;
            
            obj.p_joint_dis = p_joint_dis;
            obj.euler_ang_joint = euler_ang_joint;
            obj.p_dis_frame = p_dis_frame;
            obj.euler_ang_dis = euler_ang_dis_frame;
            
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
        function result = axis_angle_rotation(obj,i,joint_axis)
            c = cos(obj.theta_motion(i));
            s = sin(obj.theta_motion(i));
            
            a1 = joint_axis(1);
            a2 = joint_axis(2);
            a3 = joint_axis(3);
            
            result = [c+a1^2*(1-c), a1*a2*(1-c)-a3*s, a1*a3*(1-c)+a2*s;...
                 a1*a2*(1-c)+a3*s, c+a2^2*(1-c), a2*a3*(1-c)-a1*s;...
                 a1*a3*(1-c)-a2*s, a2*a3*(1-c)+a1*s, c+a3^2*(1-c)];
        end
    end
    
    
end