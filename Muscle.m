classdef Muscle < matlab.mixin.SetGet
    properties
      muscle_name              % Muscle name
      muscle_index             % Muscle index in .aproj
      pos_attachments          % Muscle attachment position in local frame
      
      % Muscle parameters
      x_off                     % X offset in V (-.04=-40mV)
      max_force                 % Stimulus-Tension Amplitude in Newtons
      steepness                 % Steepness
      y_off                     % Y offset in Newtons
      damping                   % Damping constant in Ns/m
      Kse                       % Serial stiffness in N/m
      Kpe                       % Parallel stiffness in N/m
      l_min                     % Minimum muscle length
      l_max                     % Maximum muscle length
      RestingLength             % Resting length in meters-
      l_width                   % Muscle width
      
      % Muscle data from literature or experiment
      mass                           % Muscle mass in [mg]
      lf_lm                          % Lf/Lm, fiber length to muscle length ratio from Charles 2013
      opt_fiber_length               % Optimal muscle fiber length in [mm]
      opt_muscle_length              % Optimal muscle length calculated based on lf/lm and muscle length [mm]
      pennation_angle                % Pinnation angle in [deg]
      Po                             % Maximum isometric force in [g]
      tendon_sl                      % Tendon slack length in [mm]
      vmax_fiber                     % Maximum fiber shortening velocity in [mm/s]
      
      
      
      % Muscle motion profile
      muscle_length_profile          % Muscle length during walking 
      muscle_velocity_profile        % Muscle Velocity during walking 
    end
    methods
        function obj = Muscle()
            
        end
        %% Muscle Experimental data load function
        function load_par(obj,name,mass,opt_fiber_length,pennation_angle,Po,vmax_fiber,tendon_sl,opt_muscle_length,lf_lm)
            % see if the name matches
            if contains(obj.muscle_name,name)
                obj.mass = mass; 
                obj.opt_fiber_length = opt_fiber_length;
                obj.pennation_angle = pennation_angle;
                obj.Po = Po;
                obj.vmax_fiber = vmax_fiber;
                obj.tendon_sl = tendon_sl;
                obj.opt_muscle_length = opt_muscle_length;
                obj.lf_lm = lf_lm;
            else
                disp('Muscle name not match,please check data file.')
            end
        end
        %% Muscle Length and vector function
        function musc_len = muscle_length(obj,step)
            musc_len = 0;
            for i = 2:size(obj.pos_attachments,1)
                musc_len = musc_len + norm(obj.pos_attachments{i-1,4}(step,:) - obj.pos_attachments{i,4}(step,:));
            end
        end 
        function musc_vec = muscle_unit_vector(obj,step)
            for i = 2:size(obj.pos_attachments,1)
                temp_vec = obj.pos_attachments{i-1,4}(step,:) - obj.pos_attachments{i,4}(step,:);
                musc_vec = temp_vec/norm(temp_vec);
            end
        end
        %% Muscle Max/Min length function
        function compute_min_max_muscle_lengths(obj)
            obj.l_min = min(obj.muscle_length_profile);
            obj.l_max = max(obj.muscle_length_profile);         
        end
        %% Muscle reseting length function
        function compute_reseting_length(obj,Cabs_bodies,pos_bodies)
            pos_attach = cell(size(obj.pos_attachments,1),1);
            musc_len = 0;
            for i = 1:size(obj.pos_attachments,1)
                body = obj.pos_attachments{i,3};
                if body == 1
                    pos_attach{i,:} = Cabs_bodies(:,:,1)*obj.pos_attachments{i,1};
                elseif body == 2
                    pos_attach{i,:} = Cabs_bodies(:,:,1)*pos_bodies(:,2) + Cabs_bodies(:,:,2)*obj.pos_attachments{i,1};
                elseif body == 3
                    pos_attach{i,:} = Cabs_bodies(:,:,1)*pos_bodies(:,2) + Cabs_bodies(:,:,2)*pos_bodies(:,3) + ...
                        Cabs_bodies(:,:,3)*obj.pos_attachments{i,1};
                elseif body == 4
                    pos_attach{i,:} = Cabs_bodies(:,:,1)*pos_bodies(:,2) + Cabs_bodies(:,:,2)*pos_bodies(:,3) + ...
                        Cabs_bodies(:,:,3)*pos_bodies(:,4) + Cabs_bodies(:,:,4)*obj.pos_attachments{i,1};
                end
                
                if i > 1
                    musc_len = musc_len + norm(pos_attach{i,:}-pos_attach{i-1,:});
                end
            end
            obj.RestingLength = musc_len;
        end
        function setup_reseting_length(obj)
            %Sets the resting length at the maximum length so it always
            %gets stronger as the muscle gets longer
            obj.RestingLength = obj.l_max;
        end
        %% Muscle Lwidth length function
        function compute_l_width(obj)
            %Sets the l-width propery such that at the smallest used
            %length, the strength is at the desired percent
            
            %Huristic gotten from estimating strength of Ekeberg\Pearson
            %cat muscles at their weakest
            %Changed to smaller 15/1/13 to make sure muscles don't pull too
            %strongly when they are short (especially helpful for ankle
            %extension), immediately changed back
            percent_at_edge = .7;
            
            musc_range = obj.l_max-obj.l_min;
            obj.l_width = abs(musc_range)/sqrt(1-percent_at_edge); %m
        end
        %% Muscle damping function
        function compute_damping(obj)
            % Linear damping coefficient af = a/Fmax = b/vmax
            obj.damping = obj.Po/(obj.vmax_fiber)*9.81;
        end
        %% Muscle Series Elastic function
        function compute_Kse(obj,max_musc_deflection)
            %KSE is calculated based on the idea that under max load, the
            %muscle deflects a certain percent: 1/C*rest_length
            factor = (1-obj.lf_lm)*max_musc_deflection; % legnth change ratio caused by series elasticity.
            obj.Kse = (obj.Po/1000*9.81)/(obj.RestingLength*factor);
        end
        %% Muscle Parallel Elastic function
        function compute_Kpe(obj,max_musc_deflection)
            %KPE is calculated based on the idea that under max load, the
            %SEC and PEC absorb all the force
            
            %factor = 1;
            %obj.Kpe = factor*obj.Kse*(obj.Po/1000*9.81)/((obj.Po/1000*9.81)+obj.Kse*obj.l_max*max_musc_deflection);
            
            obj.Kpe = obj.Kse*(obj.Po/1000*9.81)/(obj.Kse*obj.l_max*max_musc_deflection-obj.Po/1000*9.81);
        end
        %% Muscle Passive Tension function
        function [Tension]= compute_passive_tension(obj,dt)
            x = [obj.muscle_length_profile;obj.muscle_length_profile;obj.muscle_length_profile];
            xdot = [obj.muscle_velocity_profile;obj.muscle_velocity_profile;obj.muscle_velocity_profile];
            T = NaN(length(x),1);
            deltax = max(0,x - obj.RestingLength);
            
%             for i = 1:length(x)
%                 dT = obj.Kse*obj.Kpe/obj.damping*x(i) + obj.Kse*xdot(i) - T(i)*(obj.Kse + obj.Kpe)/obj.damping;
%                 T(i+1) = T(i)+ dT*dt;
%             end

            True = 1;
            a = 0;
            while True ==1
                Ten = Tension_approximation(deltax,xdot,dt,obj.Kpe,obj.Kse,obj.damping,T);
                T(1) = Ten(end);
                a = a+1;
                if  Ten(1) == Ten(length(x)/3+1)|| a > 10000
                    True = 0;
                end
            end
            Tension = Ten(length(x)/3+1:2*length(x)/3);
            % resampling
            num_of_samples = 100; %number of data points we wish to have per step
            current_size = size(Tension,1);
            Tension = interp1(linspace(0,1,current_size),Tension,linspace(0,1,num_of_samples));
%             keyboard
        end
    end
end