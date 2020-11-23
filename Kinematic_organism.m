classdef Kinematic_organism < matlab.mixin.SetGet
    properties
        proj_file;
        organism_name;
        body_weight;
        leg_obj = {};
        mirrored;
        num_legs;
        joint_limits;
    end
    
    methods
        function obj = Kinematic_organism(proj_file, name, bodies, joints, joint_limits, min_step_period, body_weight, theta_offset, torque_range, mirrored, to_plot )
            constructor = tic;
            
            obj.proj_file = proj_file;
            obj.organism_name = name;
            obj.body_weight = body_weight;
            
            obj.mirrored = mirrored;
            
            obj.num_legs = size(bodies,2);
            obj.joint_limits = joint_limits;
            
            %First, make sure all "bodies" chains start with the same body,
            %which will be the main body segment.
            all_body_chains_acceptable = 1;
            i = 1;
            while i < size(bodies,2) && all_body_chains_acceptable
                if ~strcmp(bodies{1,i},bodies{1,i+1})
                    all_body_chains_acceptable = 0;
                    disp('All body segment chains must start with the same segment.')
                end
                i = i + 1;
            end
            
            if all_body_chains_acceptable
                for i=1:obj.num_legs
                    %Each column represents another leg.
                    %num_segments = length(cellfun(@isempty,bodies(:,i))==0)-1;
                    try temp_joints = joints(:,i);
                    catch err
                        err_str = strcat('The joint list has thrown the error ',err.identifier);
                        disp(err_str);
                    end
                    
                    try temp_joint_limits = joint_limits(:,2*i-1:2*i);
                    catch err
                        keyboard
                        err_str = strcat('The joint position limits have thrown the error ',err.indentifier);
                        disp(err_str);
                        disp('The organism will be instantiated, but muscle properties and control networks cannot be designed until joint limits are provided.');
                    end
                    
                    try temp_theta_offset = theta_offset(:,i);
                    catch err
                        err_str = strcat('The joint offsets have thrown the error ',err.indentifier);
                        disp(err_str);
                    end
                    
                    %Velocity limits should be scaled to the joint's range
                    %of motion. Assume that each joint is capable of making
                    %a complete "round trip" over the period specified by
                    %the user (which encodes the maximum walking speed).
                    velocity_limits = zeros(size(joint_limits));
                    for j=1:size(velocity_limits,1)
                        for k=1:size(velocity_limits,2)/2
                            max_vel = (joint_limits(j,2*k)-joint_limits(j,2*k-1))/(1/2*min_step_period);
                            velocity_limits(j,2*k-1) = -max_vel;
                            velocity_limits(j,2*k) = max_vel;
                        end
                    end
                    
                    try temp_vel_limits = velocity_limits(:,2*i-1:2*i);
                    catch err
                        err_str = strcat('The joint velocity limits have thrown the error ',err.indentifier);
                        disp(err_str);
                    end

                    obj.leg_obj{i} = Leg(proj_file,name,bodies(:,i),temp_joints,temp_joint_limits,temp_theta_offset,temp_vel_limits,[],to_plot);
                    disp(['Leg ',num2str(i),' constructed.']);
                end
            end
            
            t_construction = toc(constructor)/60;
            fprintf('It took %f minutes to construct the organism object.\n',t_construction);

        end
        %%   Load the kinematc and dynamic data  
        function success = load_and_process_kinematic_and_dynamic_data(obj,Data_file,StrideTime,to_plot)
            load (Data_file)
            load('JointTorque.mat');
            
            current_size = size(Theta,1);
            num_of_samples = 100; %number of data points we wish to have per step
%             num_of_samples = 590;
            
            if exist('AllTime','var')
                total_length = AllTime(end); %Length of time that the data represents, in seconds
            else
                total_length = StrideTime;
            end            
            dt = 1*total_length/(num_of_samples-1); %time steps 
                      
            time = linspace(0,1,current_size); %make a normalized time vector (how the data is loaded)
            newtime = 0:dt:total_length; %make the desired time vector (how the data will be interpolated)
            
            Theta = Theta-pi/2; % change coordinations
                       
            Theta = interp1(linspace(0,1,current_size),Theta,linspace(0,1,num_of_samples));
            Theta_dot = interp1(linspace(0,1,current_size),Theta_dot,linspace(0,1,num_of_samples));
            Theta_doubledot = interp1(linspace(0,1,current_size),Theta_doubledot,linspace(0,1,num_of_samples));
            Torque = interp1(linspace(0,1,current_size),JointTorque,linspace(0,1,num_of_samples));
            
            for i=1:obj.num_legs
                obj.leg_obj{i}.dt_motion = newtime';
                obj.leg_obj{i}.theta_motion = Theta;
                obj.leg_obj{i}.theta_dot_motion = Theta_dot;
                obj.leg_obj{i}.theta_doubledot_motion = Theta_doubledot;
                obj.leg_obj{i}.torque_motion = Torque;
            end
            
            if to_plot == 1
                     
                figure(1)
                plot(Theta,'Linewidth',2)
                title('Joint Angles (rads)')
                legend('Hip','Knee','Ankle')
                set(gcf,'Position',[5 150 500 400])
                
                figure(2)
                plot(Theta_dot,'Linewidth',2)
                title('Joint Velocities (rads/sec)')
                legend({'Hip','Knee','Ankle'},'Location','southeast')
                set(gcf,'Position',[518 150 500 400])
                
                figure(3)
                plot(Theta_doubledot,'Linewidth',2)
                title('Joint Accelarations (rads/sec^2)')
                legend({'Hip','Knee','Ankle'},'Location','southeast')
                set(gcf,'Position',[1031 150 500 400])
                
                %Decide if needs to manually export Figure
                choice = input ('Manually export figure?(Y/N): ','s');
                switch choice
                    case {'Y','y'}
                        keyboard 
                        close all
                    case {'N','n'}
                        close all
                    otherwise
                        disp('Invalid input!')
                        close all
                end
            end
            
            if to_plot == 2
                T = {'Hip','Knee','Ankle'};
                for i = 1:3
                subplot(1,3,i)
                plot(Theta(:,i),'Linewidth',2)
                set(gca,'XTick',0:50:100);
                xlabel('stance                         swing')
                title(T{i})
                end
                suptitle('Joint Angles (radians)')
            end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.set_muscle_properties();
%             obj.compute_all_jacobians([]);
%             obj.compute_all_EOM([]);
%             obj.compute_all_MN();
            obj.design_neuron_controller();
            keyboard
            success = 1;
        end
        %%   Set the muscle properties
        function  set_muscle_properties(obj)
            to_plot = 0;
            for i=1:obj.num_legs
                obj.leg_obj{i}.store_jointbody_position(to_plot)
                obj.leg_obj{i}.store_muscle_lengthvelocity(to_plot);
                obj.leg_obj{i}.store_muscle_parameters(to_plot);
                
                % Call Plot Functions
                obj.leg_obj{i}.plot_walking_motion(to_plot);
                obj.leg_obj{i}.plot_phase_passive_tension(to_plot);
            end
        end
        %%   Calculate All Jacobian
        function compute_all_jacobians(obj,configs) 
            if isempty(configs)
                configs = zeros(size(obj.joint_limits,1),obj.num_legs);
            elseif size(configs,1) == size(obj.joint_limits,1) && size(configs,2) == 1
                configs = repmat(configs,1,obj.num_legs);
            elseif size(configs,1) == size(obj.joint_limits,1) && size(configs,2) == obj.num_legs
                %configs is the proper size
            end
            
            for i=1:obj.num_legs
                obj.leg_obj{i}.compute_jac(configs);                        
            end
        end
        %%   Calculate Inertia Forces and Gravitational Forces (EOM)
        function compute_all_EOM(obj,configs)
            if isempty(configs)
                configs = zeros(size(obj.joint_limits,1),obj.num_legs);
            elseif size(configs,1) == size(obj.joint_limits,1) && size(configs,2) == 1
                configs = repmat(configs,1,obj.num_legs);
            elseif size(configs,1) == size(obj.joint_limits,1) && size(configs,2) == obj.num_legs
                %configs is the proper size
            end
            to_plot = 0;
            
            for i=1:obj.num_legs
                obj.leg_obj{i}.compute_EOM1(configs);
                obj.leg_obj{i}.compute_EOM2(configs);
                obj.leg_obj{i}.compute_ground_EOM(configs,0);
                
                % Call Plot functions
                obj.leg_obj{i}.plot_ground_walking_motion(to_plot);
                obj.leg_obj{i}.plot_EOM(to_plot);
            end
        end
        %%   Calculate Motoneuron activation during walking
        function compute_all_MN(obj)
            to_plot = 1;                   
            for i=1:obj.num_legs
%                 obj.leg_obj{i}.compute_joint_moment_arm(to_plot);
                obj.leg_obj{i}.compute_MN_act_for_motion(to_plot);
            end
        end
        %%  Design Synthetic Nervous System
        function design_neuron_controller(obj)
           
             for i=1:obj.num_legs
                obj.leg_obj{i}.load_neurons();
                obj.leg_obj{i}.design_synapse();
            end
        end
    end
end