classdef Leg < matlab.mixin.SetGet
    %JACOBIAN Create a Jacobian manipulator object from an AnimatLab
    %simulation file.
    %   Detailed explanation goes here
    %% Display output (Ctrl+G+line or F2 to jump)
    % Line 247-241    Display bodies and joints
    % Line 306-308    Display muscle names
    % Line 358-363    Display attachment associate with muscle
    % Line 452-460    Display attachment associate with body and it's position
    
    properties
        proj_file; %project file
        original_text; %content of the .asim file
        organism_name; %name of the organism in the simulation file that we are interested in
        
        num_bodies; %number of bodies in the leg, as given by the user's list 
        bodies; %names of the segments in the leg
        joints; %names of the joints in the leg
        joint_types; %type of joints in the leg, as read from the .asim file
        joint_obj; %cell array of joint objects (for each joint in the leg)
        musc_inds;
        musc_obj;%cell array of muscle objects (for each muscle in the leg)
        to_plot; %boolean of whether or not to plot graphs of the leg.
        vec_len; %length of the axis vector drawn in the kinematic maps.
        
        Cabs_bodies; %Absolute rotation of the bodies
        Cabs_joints; %Absolute orientation of the joints
        C_joints; %Rotation of the joints as per theta
        CR_bodies; %Relative rotation of the bodies with respect to the proximal one
        CR_joints; %Relative rotation of the joints with respect to the local frame
        cum_C_joints; %Cumulative rotation (product) of proximal joints
        pos_bodies; %position of the distal body in the proximal frame
        pos_joints; %location of the distal joint in the proximal frame
        pos_attachments; %location of muscle attachments in their own frame   
        euler_angs_joints;
        euler_angs_bodies;
        theta_range; %min and max rotation for each joint
        velocity_range;

        p_st; %Position of the body relative to the foot on the ground.
        leg_attach_pt;
        foot_pos_residual;
        control_current; %Stimuli to joint position controllers necessary for each "config."
        body_twists; %Body frame transformations desired from a standing still posture. These correspond to "config."
        body_thetas;
        toe_pos_known = 0; %Boolean about whether the toe position is known or not.
        foot_vec; %vector that describes the toe wrt the most distal joint. Important for calculating forces from the foot.
        ankle_factor; %multiply distal max torque by 3 if it operates like an ankle.
    
        
        dt_motion;
        theta_motion;
        theta_dot_motion;
        torque_motion;
        
        
        gsyn_hip_knee_out;
        elo_hip_knee_out;
        ehi_hip_knee_out;
    end

    methods
        function obj = Leg(proj_file,organism_name,bodies,joints,theta_range,theta_offset,velocity_range,torque_range,to_plot)
        %%  Open the simulation files 
            obj.joint_types = cell(1,length(joints));
            obj.organism_name = organism_name;
            obj.bodies = bodies;
            obj.joints = joints;
            obj.to_plot = to_plot;
            
            obj.num_bodies = 0;
            i = 1;
            while i <= length(bodies) && ~isempty(bodies{i})
                obj.num_bodies = obj.num_bodies + 1;
                i = i + 1;
            end
            obj.theta_range = theta_range;
            obj.velocity_range = velocity_range;

            %Read the directory provided by the user.
            dir_file_division = max(strfind(proj_file,'\'));
            proj_folder = proj_file(1:dir_file_division-1);
            sim_files = dir(proj_folder);
            obj.proj_file = proj_file;

            %How many files are in that folder?
            len = length(sim_files);
            sim_file_str = [];

            %Pull out the simulation file
            for i=1:len
                name_str = sim_files(i).name;        
                if length(name_str) > 14 && strcmp(name_str(end-14:end),'Standalone.asim')
                    %We've found the file we need
                    sim_file_str = [proj_folder,'\',name_str];
                end
            end

            %Load the text of the simulation file
            try obj.original_text = importdata(sim_file_str);
            catch
                disp('No simulation file exists. Open AnimatLab and export a standalone simulation.')
                keyboard
            end
            
            %Find where the list of organisms begins
            organisms_found = contains(obj.original_text,'<Organisms>');
            organisms_lower_limit = find(organisms_found,1,'first');

            %Find where the organism we want begins
            organism_found = contains(obj.original_text(organisms_lower_limit:end),['<Name>',obj.organism_name,'</Name>']);
            organism_lower_limit = find(organism_found,1,'first') + organisms_lower_limit - 1;

         %% Initialize Body and Joint Rotations for Inertial and Local Frame
            %Orientation of the joint axis, in the inertial frame
            obj.Cabs_joints = zeros(3,3,obj.num_bodies);
            obj.Cabs_joints(:,:,1) = eye(3);

            %Rotation of the joint axis, in the local frame
            obj.C_joints = zeros(3,3,obj.num_bodies);
            obj.C_joints(:,:,1) = eye(3);

            %Rotation of each body relative to the previous
            obj.CR_bodies = zeros(3,3,obj.num_bodies); 
            obj.CR_joints = zeros(3,3,obj.num_bodies);

            %Rotation of each body in the inertial frame
            obj.Cabs_bodies = zeros(3,3,obj.num_bodies);

            %Cumulative body rotation
            cum_body_rot = eye(3);

         %% Initialize Body and Joint Positions in the Local Reference Frames
            %Position of the joints and bodies, in their local frames.
            obj.pos_joints = zeros(3,obj.num_bodies);
            obj.pos_bodies = zeros(3,obj.num_bodies);
            obj.pos_attachments = cell(obj.num_bodies,1); %one for each joint
            euler_angs_bodies = zeros(3,obj.num_bodies);
            obj.euler_angs_joints = zeros(3,obj.num_bodies);
            
        %% Fill in body and joint details          
            i = 1;
            while i <= obj.num_bodies
                % Find the body of interest
                if ~isempty(obj.bodies{i})
                    body_found = contains(obj.original_text,['<Name>',obj.bodies{i},'</Name>']);
                    next_body_ind = find(body_found,1,'first');
                    chain_lower_limit = next_body_ind;

                    if isempty(chain_lower_limit)
                        disp('An improper body or joint was entered, or perhaps out of order.')
                    end
                    
                    %Find the string that lists its rotation
                    body_rot_found = contains(obj.original_text(chain_lower_limit:end),'<Rotation');
                    body_rot_ind = find(body_rot_found,1,'first') + chain_lower_limit - 1;
                    body_rot_str = obj.original_text(body_rot_ind);

                    %Read that body's rotation angles and construct a rotation matrix.
                    quote_locs = cell2mat(strfind(body_rot_str,'"'));
                    for j=1:length(quote_locs)/2
                        euler_angs_bodies(j,i) = str2double(body_rot_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
                    end
                    obj.euler_angs_bodies(:,i) = euler_angs_bodies(:,i);
                    
                    %Compute the rotation matrix for each axis.
                    rotx_b = [1,0,0;0,cos(euler_angs_bodies(1,i)),-sin(euler_angs_bodies(1,i));0,sin(euler_angs_bodies(1,i)),cos(euler_angs_bodies(1,i))];
                    roty_b = [cos(euler_angs_bodies(2,i)),0,sin(euler_angs_bodies(2,i));0,1,0;-sin(euler_angs_bodies(2,i)),0,cos(euler_angs_bodies(2,i))];
                    rotz_b = [cos(euler_angs_bodies(3,i)),-sin(euler_angs_bodies(3,i)),0;sin(euler_angs_bodies(3,i)),cos(euler_angs_bodies(3,i)),0;0,0,1];

                    %Record the relative rotation of this frame with regard to the
                    %previous, as well as the absolute rotation with respect to the ground.
                    obj.CR_bodies(:,:,i) = rotx_b*roty_b*rotz_b;
                    cum_body_rot = cum_body_rot * obj.CR_bodies(:,:,i);
                    obj.Cabs_bodies(:,:,i) = cum_body_rot;

                    %Find the position of the body
                    body_pos_found = contains(obj.original_text(chain_lower_limit:end),'<Position');
                    body_pos_ind = find(body_pos_found,1,'first') + chain_lower_limit - 1;
                    body_pos_str = obj.original_text(body_pos_ind);

                    quote_locs = cell2mat(strfind(body_pos_str,'"'));
                    for j=1:length(quote_locs)/2
                        obj.pos_bodies(j,i) = str2double(body_pos_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
                    end
       
                    
                    %%% Now for the joint
                    if i > 1
                        %Move up the minimum line value now that that body is done.
                        joint_found = contains(obj.original_text,['<Name>',obj.joints{i},'</Name>']);
                        next_joint_ind = find(joint_found,1,'first');
                        chain_lower_limit = next_joint_ind;

                        %Find the orientation of that joint
                        joint_rot_found = contains(obj.original_text(chain_lower_limit:end),'<Rotation');
                        joint_rot_ind = find(joint_rot_found,1,'first') + chain_lower_limit - 1;
                        joint_rot_str = obj.original_text(joint_rot_ind);

                        quote_locs = cell2mat(strfind(joint_rot_str,'"'));
                        for j=1:length(quote_locs)/2
                            obj.euler_angs_joints(j,i) = str2double(joint_rot_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
                        end

                        %Compute the rotation of the joint's axis within the local frame.
                        rotx_j = [1,0,0;0,cos(obj.euler_angs_joints(1,i)),-sin(obj.euler_angs_joints(1,i));0,sin(obj.euler_angs_joints(1,i)),cos(obj.euler_angs_joints(1,i))];
                        roty_j = [cos(obj.euler_angs_joints(2,i)),0,sin(obj.euler_angs_joints(2,i));0,1,0;-sin(obj.euler_angs_joints(2,i)),0,cos(obj.euler_angs_joints(2,i))];
                        rotz_j = [cos(obj.euler_angs_joints(3,i)),-sin(obj.euler_angs_joints(3,i)),0;sin(obj.euler_angs_joints(3,i)),cos(obj.euler_angs_joints(3,i)),0;0,0,1];

                        %Quick side track: find what kind of joint this is. If it is
                        %prismatic, it does not rotate anything. If it is a hinge, then we
                        %need to account for that rotation.
                        joint_type_found = contains(obj.original_text(chain_lower_limit:end),'<PartType');
                        joint_type_ind = find(joint_type_found,1,'first') + chain_lower_limit - 1;
                        joint_type_str = obj.original_text(joint_type_ind);

                        %pull out the joint type
                        type_begin = strfind(joint_type_str,'.');
                        type_end = strfind(joint_type_str,'<');
                        obj.joint_types{i} = joint_type_str{:}(type_begin{:}(end)+1:type_end{:}(end)-1);

                        %Calculate the orientation of the axis in space, Cabs_joints, and
                        %the rotation of the joint, C_joints
                        obj.Cabs_joints(:,:,i) = obj.Cabs_bodies(:,:,i)*rotx_j*roty_j*rotz_j;
                        obj.CR_joints(:,:,i) = rotx_j*roty_j*rotz_j;
                        obj.C_joints(:,:,i) = eye(3);

                        %Find the location of the joint within the frame
                        joint_pos_found = contains(obj.original_text(chain_lower_limit:end),'<Position');
                        joint_pos_ind = find(joint_pos_found,1,'first') + chain_lower_limit - 1;
                        joint_pos_str = obj.original_text(joint_pos_ind);

                        quote_locs = cell2mat(strfind(joint_pos_str,'"'));
                        for j=1:length(quote_locs)/2
                            obj.pos_joints(j,i) = str2double(joint_pos_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
                        end

                        %Move up the minimum line value now that that body is done.
                        try
                            joint_found = contains(obj.original_text,['<Name>',obj.bodies{i+1},'</Name>']);
                            next_joint_ind = find(joint_found,1,'first');
                            chain_lower_limit = next_joint_ind;
                        catch
                            disp('End of joint chain reached.')
                        end
                    end
                end
                i = i + 1;
            end
            
%             %Display body and joint created
%             disp('Following leg body was created:')
%             disp(obj.bodies)
%             disp('Following leg joint was created:')
%             disp(obj.joints)
            
        %% Find muscles, their attachments, and associate them with the proper bodies and joints.
            tstart = tic;
            musc_found = contains(obj.original_text,'<Type>LinearHillMuscle</Type>');
            musc_inds = find(musc_found);
            musc_name_inds = musc_inds-2;
                 
            name_true = 1;
            
            for i=2:obj.num_bodies-1
                if ~strcmp(obj.joints{i}(1:2),obj.joints{i+1}(1:2))
                    name_true = 0;
                end
            end
            
            if name_true
                leg_name = obj.joints{2}(1:2);
            else
                leg_name = input('What is the prefix of objects belonging to this leg?\nA cell array of prefixes may be entered. ');
            end
            
            if ischar(leg_name)
                musc_for_this_leg = strfind(obj.original_text(musc_name_inds),['<Name>',leg_name]);
            elseif iscell(leg_name)
                musc_for_this_leg = cell(length(musc_name_inds),1);
                for i=1:length(leg_name)
                    temp_musc_for_this_leg = strfind(obj.original_text(musc_name_inds),['<Name>',leg_name{i}]);
                    for j=1:length(temp_musc_for_this_leg)
                        if ~isempty(temp_musc_for_this_leg{j})
                            musc_for_this_leg{j} = 1;
                        end
                    end
                end
            end
            %Logically pick all the muscle from this leg
            musc_for_this_leg_inds = cellfun(@isempty,musc_for_this_leg) == 0;
            
            %These are the indices of the names of muscles.
            obj.musc_inds = musc_name_inds(musc_for_this_leg_inds);
            attach_to_find = cell(length(obj.musc_inds),1);
            
            for i = 1:length(obj.musc_inds)
                obj.musc_obj{i,1} = Muscle();
            end
            
            %Find names of muscles.
            musc_names = obj.original_text(obj.musc_inds);
            for i=1:length(musc_names)
                musc_names{i} = strrep(musc_names{i},'<Name>','');
                musc_names{i} = strrep(musc_names{i},'</Name>','');
                obj.musc_obj{i}.muscle_name = musc_names{i};
                obj.musc_obj{i}.muscle_index = obj.musc_inds(i);
            end   
           
%             %display muscle created
%             disp('Following muscle was created:')
%             disp(musc_names)
            
            %Now that we know where the muscles are saved, we need to
            %extract all of the attachment IDs, find their names (to figure
            %out which body they belong to), and save their locations to
            %create joint objects.
            for i=1:length(obj.musc_inds)
                attachment_start = contains(obj.original_text(obj.musc_inds(i):end),'<Attachments>');
                attachment_start_ind = find(attachment_start,1,'first') + obj.musc_inds(i) - 1;
                attachment_end = contains(obj.original_text(obj.musc_inds(i):end),'</Attachments>');
                attachment_end_ind = find(attachment_end,1,'first') + obj.musc_inds(i) - 1;

                attach_to_find{i} = obj.original_text(attachment_start_ind + 1:attachment_end_ind - 1);
                %remove the word "attach" from the strings
                for j=1:length(attach_to_find{i})
                    attach_to_find{i}{j} = strrep(attach_to_find{i}{j},'Attach','');
                end                
            end
            
            % now we got the ID for all those attachments
            attach_names = cell(size(attach_to_find));
            attach_locs = cell(size(attach_to_find));
            attach_names_str = cell(size(attach_to_find));
            
            %find indices of the names of these attachments.
            for i=1:length(obj.musc_inds)
                for j=1:length(attach_to_find{i})
                    %save the names
                    id_loc = contains(obj.original_text,attach_to_find{i}{j});
                    id_ind = find(id_loc,1,'first');
                    attach_names{i}{j} = id_ind - 1;  
                    
                    AN = obj.original_text(attach_names{i}{j});
                    AN = strrep(AN,'<Name>','');
                    AN = strrep(AN,'</Name>','');
                    attach_names_str{i}{j} = AN;
                    
                    %Find the position of that attachment. This mapping is
                    %identical to the names.
                    attach_pos_found = contains(obj.original_text(id_ind:end),'<Position');
                    attach_pos_ind = find(attach_pos_found,1,'first') + id_ind - 1;
                    attach_pos_str = obj.original_text(attach_pos_ind);

                    quote_locs = cell2mat(strfind(attach_pos_str,'"'));
                    for k=1:length(quote_locs)/2
                        attach_locs{i}{j}(k,1) = str2double(attach_pos_str{1}(quote_locs(2*k-1)+1:quote_locs(2*k)-1));
                    end
                end
            end
            
%             % Display attachment associate with muscles.
%             disp('Now we created attachments, below is which muscle it asscoiate with:')
%             for i = 1:length(obj.musc_inds)
%                 musc_attach = [musc_names(i);string(attach_names_str{i}')];
%                 disp(musc_attach)
%             end
            
            %We will need to remove redundant attachments. If multiple muscles use the same
            %attachment, it will show up in our list twice. Therefore we
            %will start a list used_indices, and every time we add an
            %attachment's position to our record, we add the index of its
            %name to this list. Attachments whose indices already appear on
            %the list will be ignored.
            used_indices = [];
            
            %Tally up the total number of attachments on each body so that
            %we can set the size for obj.pos_attachments
            %Cells don't allow you to increment the size
            pelvis_count = 0;
            femur_count = 0;
            tibia_count = 0;
            foot_count = 0;
            
            for j=1:length(attach_names)
                for k=1:length(attach_names{j})
                    % check is the attachment is already been used
                    if isempty(find(used_indices == attach_names{j}{k},1,'first'))
                        temp_name_str = char(obj.original_text(attach_names{j}{k}));
                        if contains(temp_name_str,' Pelvis ')
                            pelvis_count = pelvis_count + 1;
                        elseif contains(temp_name_str,' Femur ')
                            femur_count = femur_count + 1;
                        elseif contains(temp_name_str,' Tibia ')
                            %spaces around Tibia important bc muscle
                            %'Tibialis Anterior' contains Tibia
                            tibia_count = tibia_count + 1;
                        elseif contains(temp_name_str,' Foot ')
                            foot_count = foot_count + 1;
                        else
                            %attach_names{j}{k}
                            keyboard
                            error('Attachment not labeled with a body part')
                        end
                        used_indices = [used_indices,attach_names{j}{k}];
                    end
                end
            end
            
            used_indices = used_indices';
            
            %Set the size of obj.pos_attachments for each body (Position,Name,Body)
            obj.pos_attachments{1} = cell(pelvis_count,3);
            obj.pos_attachments{2} = cell(femur_count,3);
            obj.pos_attachments{3} = cell(tibia_count,3);
            obj.pos_attachments{4} = cell(foot_count,3);
            
            %Store information about attachment points for each body in
            %obj.pos_attachments. Resultant cell array holds position
            %information for XYZ and name strings of every unique
            %attachment point on each body
            for j=1:length(attach_names)
                for k=1:length(attach_names{j})
                    temp_name_str = attach_names_str{j}{k};
                    if contains(temp_name_str,' Pelvis ')
                            body = 1;
                        elseif contains(temp_name_str,' Femur ')
                            body = 2;
                        elseif contains(temp_name_str,' Tibia ')
                            body = 3;
                        elseif contains(temp_name_str,' Foot ')
                            body = 4;
                        else
                            keyboard
                            error('Attachment not labeled with a body part')
                    end
                    
                    % Saving information for muscle properties
                    obj.musc_obj{j}.pos_attachments{k,1} = attach_locs{j}{k};
                    obj.musc_obj{j}.pos_attachments{k,2} = attach_names_str{j}{k};
                    obj.musc_obj{j}.pos_attachments{k,3} = body;
                    obj.musc_obj{j}.pos_attachments{k,4} = [];
                    % Saving informations for attachment, follows
                    % position,name.body attached to 
                    if ~isempty(find(used_indices == attach_names{j}{k},1,'first'))
                        used_indices(find(used_indices == attach_names{j}{k},1,'first')) = 0;
                        i = sum(~cellfun(@isempty,obj.pos_attachments{body}(:,1)))+1;
                        obj.pos_attachments{body}{i,1} = attach_locs{j}{k};
                        obj.pos_attachments{body}{i,2} = attach_names_str{j}{k};
                        obj.pos_attachments{body}{i,3} = obj.bodies(body); 
                    else
                        %skip this reused attachment
                    end               
                end
            end
            
%             % Display attachment position and which body part it associate with.
%             disp('Position of attachements and which body part it attached with:')
%             for i = 1:length(obj.pos_attachments)
%                 disp(string(obj.pos_attachments{i}(1,3)))
%                 for j = 1:length(obj.pos_attachments{i})
%                 attach_pos = [string(obj.pos_attachments{i}(j,2)),cell2mat(obj.pos_attachments{i}(j,1))'];
%                 disp(attach_pos)
%                 end
%             end
            
            telapsed = toc(tstart);
            disp(['Attachments Stored.',' (',num2str(telapsed),'s)'])
            
        %% Store joint properties   
            %Now, we have all the information we need to construct a cell
            %of joint objects, which make up the leg.
            obj.joint_obj = cell(obj.num_bodies,1);
            for i = 2:obj.num_bodies
                if isempty(torque_range)
                    obj.joint_obj{i} = Joint(obj.pos_joints(:,i),obj.euler_angs_joints(:,i),obj.pos_bodies(:,i),euler_angs_bodies(:,i),theta_range(i-1,:),theta_offset(i-1),velocity_range(i-1,:),[],joints{i},proj_file,0);
                else
                    obj.joint_obj{i} = Joint(obj.pos_joints(:,i),obj.euler_angs_joints(:,i),obj.pos_bodies(:,i),euler_angs_bodies(:,i),joints{i},theta_range(i-1,:),theta_offset(i-1),velocity_range(i-1,:),torque_range(i-1,:),proj_file,0);
                end       
            end
        
        end 
        %% Calculate joints and bodies position during walking
        function store_jointbody_position(obj,to_plot)
            tstart = tic;
            % Store joints motion in joint objects 
            for i = 2:length(obj.joints)
               obj.joint_obj{i}.dt_motion = obj.dt_motion;
               obj.joint_obj{i}.theta_motion = obj.theta_motion(:,i-1);
               obj.joint_obj{i}.theta_dot_motion = obj.theta_dot_motion(:,i-1);
            end
            
            %Store joint axis during walking in joit objects
            for i = 2:length(obj.joints)
                uu_joint = obj.CR_bodies(:,:,i)*obj.CR_joints(:,:,i)*[-1;0;0];
                obj.joint_obj{i}.axis1 = uu_joint(1);
                obj.joint_obj{i}.axis2 = uu_joint(2);
                obj.joint_obj{i}.axis3 = uu_joint(3);
            end
            
            %Initialize Body and Joint Positions in joint objects            
            sizer = zeros(length(obj.dt_motion),3);
            pos_femur = sizer;
            pos_tibia = sizer;
            pos_foot = sizer;
            pos_knee = sizer;
            pos_ankle = sizer;
            
            % Position of bodies and joints during walking!
            pos_hip_in = obj.Cabs_bodies(:,:,1)*obj.pos_bodies(:,2)+ obj.Cabs_bodies(:,:,2)*obj.pos_joints(:,2);
            pos_hip = sizer + pos_hip_in' ;% Hip joint is located in Pelvis, so the position is constant.
            for i = 1:length(obj.dt_motion)
                a = obj.joint_obj{2}.axis_angle_rotation(i); % axis rotation of hip
                b = obj.joint_obj{3}.axis_angle_rotation(i); % axis rotation of knee
                c = obj.joint_obj{4}.axis_angle_rotation(i); % axis rotation of ankle
                
                hip_rel = obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2);
                pos_femur(i,:) = pos_hip(i,:)' - hip_rel*obj.pos_joints(:,2);
                
                knee_rel = hip_rel*b*obj.CR_bodies(:,:,3);
                pos_knee(i,:) = pos_femur(i,:)' + hip_rel*obj.pos_bodies(:,3)+ hip_rel*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
                pos_tibia(i,:) = pos_knee(i,:)' - knee_rel*obj.pos_joints(:,3);
                
                ankle_rel = knee_rel*c*obj.CR_bodies(:,:,4);
                pos_ankle(i,:) = pos_tibia(i,:)' + knee_rel*obj.pos_bodies(:,4) + knee_rel*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
                pos_foot(i,:) = pos_ankle(i,:)' - ankle_rel*obj.pos_joints(:,4);                
            end
                       
            %Store joints position during walking in obj.joint
            obj.joint_obj{2}.joint_motion = pos_hip;
            obj.joint_obj{3}.joint_motion = pos_knee;
            obj.joint_obj{4}.joint_motion = pos_ankle;
            %Store bodies position during walking in obj.joint
            obj.joint_obj{2}.body_motion = pos_femur;
            obj.joint_obj{3}.body_motion = pos_tibia;
            obj.joint_obj{4}.body_motion = pos_foot;
            
            telapsed = toc(tstart);
            disp(['Walking motion captured.',' (',num2str(telapsed),'s)'])
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
            % Save data and plot walking motion
            if to_plot == 1
                dt = obj.dt_motion;
                save ('motion.mat','dt','pos_hip','pos_femur','pos_knee','pos_tibia','pos_ankle','pos_foot')
                Plot_motion()
                %Decide if needs to manually export Figure
                keyboard 
                close all
            end
        end
        %% Calculate muscle length and velocity during walking
        function store_muscle_lengthvelocity(obj,to_plot)
            %Uses the joint rotation matrices at each timestep to calculate the muscle length and velocity. Stores this as a profile in the muscle object.
            tstart = tic;
            
            for i = 1:length(obj.dt_motion)
                % Calculate joint ratation matrices
                a = obj.joint_obj{2}.axis_angle_rotation(i); % axis rotation of hip
                b = obj.joint_obj{3}.axis_angle_rotation(i); % axis rotation of knee
                c = obj.joint_obj{4}.axis_angle_rotation(i); % axis rotation of ankle
                
                for j = 1:length(obj.musc_obj)
                    % Calculate attachment world position during walking
                    for k = 1:size(obj.musc_obj{j}.pos_attachments,1)
                        body = obj.musc_obj{j}.pos_attachments{k,3};
                        if body == 1
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.CR_bodies(:,:,1)*obj.musc_obj{j}.pos_attachments{k,1};
                        elseif body == 2
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.joint_obj{2}.body_motion(i,:)'+ ...
                                obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*obj.musc_obj{j}.pos_attachments{k,1};
                        elseif body == 3
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.joint_obj{3}.body_motion(i,:)'+ ...
                                obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*b*obj.CR_bodies(:,:,3)*...
                                obj.musc_obj{j}.pos_attachments{k,1};
                        elseif body == 4    
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.joint_obj{4}.body_motion(i,:)'+ ...
                                obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*b*obj.CR_bodies(:,:,3)*...
                                c*obj.CR_bodies(:,:,4)*obj.musc_obj{j}.pos_attachments{k,1};
                        end                       
                    end
                    % Calculate muscle length during walking
                    obj.musc_obj{j}.muscle_length_profile(i,:) = obj.musc_obj{j}.muscle_length(i);
                end                
            end

            % Calculate mucle velocity
            dt = obj.dt_motion(2)- obj.dt_motion(1);
            for i = 1:length(obj.musc_obj)
                length_all = [obj.musc_obj{i}.muscle_length_profile;obj.musc_obj{i}.muscle_length_profile;obj.musc_obj{i}.muscle_length_profile];
                velocity_all = diff(length_all)/dt;
              
%                 %smoothing
%                 avgblocks = 3;
%                 coeffblocks = ones(1,avgblocks)/avgblocks;
%                
%                 for k=1:10
%                     velocity_all(:,1) = filtfilt(coeffblocks,1,velocity_all(:,1));
%                 end
                t = length(length_all)/3;
                velocity = velocity_all(t+1:2*t);
                obj.musc_obj{i}.muscle_velocity_profile = velocity;
%                 obj.musc_obj{i}.muscle_velocity_profile = diff(obj.musc_obj{i}.muscle_length_profile)/dt;
            end
            
            telapsed = toc(tstart);
            disp(['Muscle walking profile stored.',' (',num2str(telapsed),'s)'])
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
            % Save data and plot muscle length&velocity during walking
            if to_plot == 1
                dt = obj.dt_motion;
                muscle_name = cell(length(obj.musc_obj),1);
                muscle_length = muscle_name;
                muscle_velocity = muscle_name;
                for i = 1:length(obj.musc_obj)
                    muscle_name{i} = obj.musc_obj{i}.muscle_name;
                    muscle_length{i} = obj.musc_obj{i}.muscle_length_profile;
                    muscle_velocity{i} = obj.musc_obj{i}.muscle_velocity_profile;
                end
                save ('muscle.mat','muscle_name','muscle_length','muscle_velocity')
                Plot_muscle()
                %Decide if needs to manually export Figure
                keyboard 
                close all
            end
        end
        %% Calculate muscle properties(Kpe,Kse,B,Lwith,Restinglength......)
        function store_muscle_parameters(obj,to_plot)
            tstart = tic;
            load('BiarMuscleParameters.mat'); 
            
            dt = obj.dt_motion(2)- obj.dt_motion(1);
            Tension = cell(length(obj.musc_obj),1);
            muscle_name = Tension;
            
            %Huristic from Meijer 1998
            max_musc_deflection = .04;
            
            for i = 1:length(obj.musc_obj)
                % Load muscle experimental parameters(For a compare and
                % guildance)
                obj.musc_obj{i}.load_par(Musc_par{i+1,:});
                
                % Calculate muscle Max/Min length
                obj.musc_obj{i}.compute_min_max_muscle_lengths();
                % Setup muscle reseting length
%                 obj.musc_obj{i}.compute_reseting_length(obj.Cabs_bodies,obj.pos_bodies);
                obj.musc_obj{i}.setup_reseting_length();
                % Calculate muscle Lwidth length
                obj.musc_obj{i}.compute_l_width();
                
                % Calculate muscle Damping B
                obj.musc_obj{i}.compute_damping();
                % Calculate muscle Series Elastic Kse
                obj.musc_obj{i}.compute_Kse(max_musc_deflection);
                % Calculate muscle Parallel Elastic Kpe
                obj.musc_obj{i}.compute_Kpe(max_musc_deflection);
                % Calculate muscle Activation Curve
                obj.musc_obj{i}.calc_activation_curve();
                
                % Calculate muscle passive tension
                Tension{i} = obj.musc_obj{i}.compute_passive_tension(dt);
                muscle_name{i} = obj.musc_obj{i}.muscle_name;
            end
            
            telapsed = toc(tstart);
            disp(['Muscle properties stored.',' (',num2str(telapsed),'s)'])
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Save data and plot muscle passive-tension
            if to_plot == 1
                save ('Passive_tension.mat','muscle_name','Tension')
                Plot_Passive_tension()
                %Decide if needs to manually export Figure                
                keyboard 
                close all
            end
        end
        %% Calculate Jacobian matrix(Spatial Manipulator)
        function [Jac,p_foot] = compute_jac(obj,theta,new_figure,rotate_plot)
            %To clarify, this Jacobian is not the classical numerical jacobian
            %(differential of a function for each generalized coordinate). This
            %is the spatial manipulator Jacobian (see "A Mathematical Introduction to Robotic Manipulation" Murray, Li, and Sastry 1994)
            if length(theta) == obj.num_bodies - 1
                theta = [0;theta];
            elseif length(theta) == obj.num_bodies
                %Do Nothing!
            elseif isempty(theta)
                theta = zeros(obj.num_bodies,1);
            else
                disp('Orientation vector is not proper length.Should be a nx1 vector where n=num_bodies (first element is time)')
                p_foot = -1;
                Jac = -1;
                return
            end
            
            % Test sets
            theta = [0; obj.theta_motion(1,:)'];
            Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(1,:)', obj.joint_obj{3}.body_motion(1,:)', obj.joint_obj{4}.body_motion(1,:)'];
            Cj_N  = [obj.joint_obj{2}.joint_motion(1,:)', obj.joint_obj{3}.joint_motion(1,:)', obj.joint_obj{4}.joint_motion(1,:)'];
            
            % Initialization
            r_N = zeros(3,obj.num_bodies);         % Bodies world position
            j_N = zeros(3,obj.num_bodies-1);       % Joint world position
            omega = zeros(3,obj.num_bodies-1);     % Axes of joint rotation
            Jac = zeros(6,obj.num_bodies-1);       % Jacobian matrix
            
            
            for i=2:obj.num_bodies
                %This matrix describes how the limb rotates the distal segments.
                %Compute the absolute axis of rotation for the joint, in
                %the 0 position.
                u = obj.Cabs_bodies(:,:,i)*obj.CR_joints(:,:,i)*[-1;0;0];
                %Calculate the orientation of the axis in space, Cabs_joints, and
                %the rotation of the joint (joint angle), C_joints
                if strcmp(obj.joint_types{i},'Hinge')
                    %Now, find the rotation 
                    obj.Cabs_joints(:,:,i) = obj.Cabs_bodies(:,:,i)*obj.CR_joints(:,:,i)*[1,0,0;0,cos(-theta(i)),-sin(-theta(i));0,sin(-theta(i)),cos(-theta(i))];
                    
                    obj.C_joints(:,:,i) = [cos(theta(i))+u(1)^2*(1-cos(theta(i))),u(1)*(2)*(1-cos(theta(i)))-u(3)*sin(theta(i)),u(1)*u(3)*(1-cos(theta(i)))+u(2)*sin(theta(i));u(2)*u(1)*(1-cos(theta(i)))+u(3)*sin(theta(i)),cos(theta(i))+u(2)^2*(1-cos(theta(i))),u(2)*u(3)*(1-cos(theta(i)))-u(1)*sin(theta(i));u(3)*u(1)*(1-cos(theta(i)))-u(2)*sin(theta(i)),u(3)*u(2)*(1-cos(theta(i)))+u(1)*sin(theta(i)),cos(theta(i))+u(3)^2*(1-cos(theta(i)))];
                elseif strcmp(obj.joint_types{i},'Prismatic')
                    % Do Nothing!
                end 
            end
            
            
            for i=2:obj.num_bodies
                %pos_bodies should be counted from 2 to end. Each column is that body's
                %position relative to the first body (which is 0s for the first).
                %Second body's position is given in the frame of the first.
                
                %r_N is the "World Position" in Animatlab of each body
                r_N(:,i) = r_N(:,i-1) + obj.Cabs_bodies(:,:,i-1)*obj.pos_bodies(:,i);
                %Similarly, j_N is the "World Position" in Animatlab of each joint. 
                j_N(:,i-1) = r_N(:,i) + obj.Cabs_bodies(:,:,i)*obj.pos_joints(:,i);
            end
            keyboard
            
            %In order to compound the effects of rotating a proximal joint, we need to
            %compute the cumulative matrix product of the rotation of each of the
            %joints. z element of cum_C_joints holds the rotation to apply to each
            %joint-COM or joint-joint relative position vector.
            obj.cum_C_joints = zeros(size(obj.C_joints));
            obj.cum_C_joints(:,:,1) = eye(3);
            
            %vectors that point joint to body_cm and joint to joint. These are
            %generated with the points above.
            cm_vec = zeros(3,obj.num_bodies-1);
            j_vec = zeros(3,obj.num_bodies-2);
            
            %The cumulative sum of the above vectors will be stored in cm_cum_vec
            %and j_cum_vec.
            cm_cum_vec = zeros(3,obj.num_bodies-1);
            j_cum_vec = zeros(3,obj.num_bodies-1);
            
            for i=1:obj.num_bodies-1
                %Now that we have the absolute positions of each of the joints and CMs,
                %we need to make vectors that describe the position of joints and CMs
                %with respect to the proximal joint. In this way we can then rotate the
                %members into their actual orientation.

                %Since the rotation of a proximal joint affects the location and
                %orientation of distal pieces, we need to compute the cumulative
                %rotation.
                obj.cum_C_joints(:,:,i+1) = obj.cum_C_joints(:,:,i) * obj.C_joints(:,:,i+1);

                if i == 1
                    cm_vec(:,i) = r_N(:,i+1) - r_N(:,i);
                    cm_cum_vec(:,i) = cm_vec(:,i);
                    j_cum_vec(:,i) = cm_vec(:,i);
                else
                    cm_vec(:,i) = r_N(:,i+1) - j_N(:,i);
                    j_vec(:,i-1) = j_N(:,i) - j_N(:,i-1);
                    %Calculate the "cumulative vectors" of positions of the joints and
                    %centers of mass. This is important for building our forward
                    %kinematic map. We can subtract this matrix from any column within
                    %it to yield the Jacobian at that point.
                    j_cum_vec(:,i) = j_cum_vec(:,i-1) + obj.cum_C_joints(:,:,i)*j_vec(:,i-1);
                    cm_cum_vec(:,i) = j_cum_vec(:,i) + obj.cum_C_joints(:,:,i+1)*cm_vec(:,i);
                end
            end
            
            if rotate_plot
                j_cum_vec = [1,0,0;0,0,-1;0,1,0]*j_cum_vec; %%%%%%%%% To rotate plot
                cm_cum_vec = [1,0,0;0,0,-1;0,1,0]*cm_cum_vec; %%%%%%%%% To rotate plot
            end
            
            Jac_fake = Jac;
            
             %Assuming that the COM of the last segment is halfway
            %between the distal joint and the end, add that segment to
            %it for drawing.
%             keyboard
            if ~obj.toe_pos_known
                disp(obj.bodies{obj.num_bodies})

%                 toe_pos = input('Where is the most distal point \non the foot? If it is opposite \nthe center of the distal frame \nfrom the most distal joint, simply hit [Enter].\nOtherwise, please input a vector that \ndescribes the toe location with respect \nto the ankle. This should be in \nthe global frame. ');
                toe_pos = [-17.899*10^-3;-85.102*10^-3;21.338*10^-3]-[-21.997*10^-3;-63.961*10^-3;21.34*10^-3];
%                 toe_pos = [-65.09*10^-3;-40.694*10^-2;28.41*10^-2]-[-65.09*10^-3;-28.329*10^-2;25.41*10^-2];
                
                if isempty(toe_pos)
%                     j_N(:,end+1) = j_N(:,end) - 2*obj.Cabs_bodies(:,:,end)*obj.pos_joints(:,end);
                    obj.foot_vec = -2*obj.pos_joints(:,end);
                else
%                     j_N(:,end+1) = j_N(:,end) + toe_pos;
                    obj.foot_vec = obj.Cabs_bodies(:,:,end)'*toe_pos;
                end
                obj.toe_pos_known = 1;
            end            
            
            j_N(:,end+1) = j_N(:,end) + obj.Cabs_bodies(:,:,end)*obj.foot_vec;
                
            obj.p_st = j_N - repmat(j_N(:,1),1,size(j_N,2));
            obj.leg_attach_pt = j_N(:,1);

            if obj.to_plot
                figure
                clf
                hold on
                
                if rotate_plot
                    j_N = [1,0,0;0,0,-1;0,1,0]*j_N; %%%%%%%% To rotate plot
                    r_N = [1,0,0;0,0,-1;0,1,0]*r_N; %%%%%%%% To rotate plot
                end
                
                plot3(j_N(1,:),j_N(2,:),j_N(3,:),'bo','Linewidth',2)
                plot3(j_N(1,:),j_N(2,:),j_N(3,:),'b','Linewidth',2)
%                 plot3(r_N(1,:),r_N(2,:),r_N(3,:),'bo','Linewidth',2)
            end
            
            obj.vec_len = .002;
            
            for i=1:obj.num_bodies-1
                if rotate_plot
                    axis_vec = obj.vec_len*[1,0,0;0,0,-1;0,1,0]*obj.Cabs_joints(:,:,i+1)*[-1;0;0]; %%%%%%%  To rotate plot
                else
                    axis_vec = obj.vec_len*obj.Cabs_joints(:,:,i+1)*[-1;0;0];
                end
                axis_pts = [j_N(:,i),j_N(:,i) + axis_vec];

                %The axes of rotation have been rotated by proximal joint movements, so
                %account for that.
                omega(:,i) = obj.cum_C_joints(:,:,i+1)*obj.Cabs_joints(:,:,i+1)*[-1;0;0];
                
                %Normalize the rotation direction
                omega(:,i) = omega(:,i)/norm(omega(:,i));                
                axis_vec_rotated = obj.vec_len*omega(:,i);
                axis_pts_rotated = [j_cum_vec(:,i),j_cum_vec(:,i) + axis_vec_rotated];

                if obj.to_plot
                    plot3(axis_pts(1,:),axis_pts(2,:),axis_pts(3,:),'r')
                    plot3(axis_pts_rotated(1,:),axis_pts_rotated(2,:),axis_pts_rotated(3,:),'cyan')
                end
            end

            if obj.to_plot
                % axis equal
                xlabel('x')
                ylabel('y')
                zlabel('z')
                if rotate_plot
                    view(-42,36)
                else
                    view(47,26)
                end
                view(120,30)
                grid on
                hold off
            end

            %Now we need to add the tarsus location to and remove the body location
            %from j_cum_vec. Since there is not actually a joint at the center of mass
            %of the body, we do not have any torque to calculate about that point, so
            %that lever arm is useless to us. Instead, we'll assume that the center of
            %mass of the most proximal segment is in its geometical center, meaning the
            %final radius is 2*cum_C_joints(:,:,num_bodies)*cm_vec(:,num_bodies-1).
            
            if rotate_plot
%                 j_cum_vec(:,end+1) = j_cum_vec(:,obj.num_bodies-1) + 2*[1,0,0;0,0,-1;0,1,0]*obj.cum_C_joints(:,:,obj.num_bodies)*cm_vec(:,obj.num_bodies-1); %%%%%%%%
                j_cum_vec(:,end+1) = j_cum_vec(:,obj.num_bodies-1) + [1,0,0;0,0,-1;0,1,0]*obj.cum_C_joints(:,:,obj.num_bodies)*obj.Cabs_bodies(:,:,end)*obj.foot_vec; %%%%%%%%
            else
%                 j_cum_vec(:,end+1) = j_cum_vec(:,obj.num_bodies-1) + 2*obj.cum_C_joints(:,:,obj.num_bodies)*cm_vec(:,obj.num_bodies-1);
                j_cum_vec(:,end+1) = j_cum_vec(:,obj.num_bodies-1) + obj.cum_C_joints(:,:,obj.num_bodies)*obj.Cabs_bodies(:,:,end)*obj.foot_vec;
            
            end
            
            for i=1:obj.num_bodies-1
                if strcmp(obj.joint_types{i+1},'Hinge')
                    Jac(1:3,i) = -cross(omega(:,i),j_cum_vec(:,i)); 
            %         Jac(1:3,i) = -cross(omega(:,i),j_cum_vec(:,i) - j_cum_vec(:,end)); %This is the Jacobian for a force on the tarsus

                    Jac_fake(1:3,i) = -cross(omega(:,i),j_cum_vec(:,i) - j_cum_vec(:,end));
                    Jac(4:6,i) = omega(:,i);
                    Jac_fake(4:6,i) = omega(:,i);
                elseif strcmp(obj.joint_types{i+1},'Prismatic')
                    Jac(1:3,i) = omega(:,i);
                    Jac_fake(1:3,i) = omega(:,i);
                else
                    disp('Unidentified joint type.')
                end
            end

            if obj.to_plot
                figure(gcf)
                hold on
                plot3(j_cum_vec(1,:),j_cum_vec(2,:),j_cum_vec(3,:),'magenta--','Linewidth',2)
                plot3(j_cum_vec(1,:),j_cum_vec(2,:),j_cum_vec(3,:),'magentao','Linewidth',2)
                
                pminx = 2*min([min(j_N(1,:)),min(r_N(1,:)),min(j_cum_vec(1,:)),min(cm_cum_vec(1,:))]);
                pminy = 2*min([min(j_N(2,:)),min(r_N(2,:)),min(j_cum_vec(2,:)),min(cm_cum_vec(2,:))]);
                pminz = 2*min([min(j_N(3,:)),min(r_N(3,:)),min(j_cum_vec(3,:)),min(cm_cum_vec(3,:))]);
                pmaxx = 2*max([max(j_N(1,:)),max(r_N(1,:)),max(j_cum_vec(1,:)),max(cm_cum_vec(1,:))]);
                pmaxy = 2*max([max(j_N(2,:)),max(r_N(2,:)),max(j_cum_vec(2,:)),max(cm_cum_vec(2,:))]);
                pmaxz = 2*max([max(j_N(3,:)),max(r_N(3,:)),max(j_cum_vec(3,:)),max(cm_cum_vec(3,:))]);
                
%                 axis([pminx,pmaxx,pminy,pmaxy,pminz,pmaxz])

%                 plot3(cm_cum_vec(1,:),cm_cum_vec(2,:),cm_cum_vec(3,:),'magenta*')
                axis equal
                hold off
            end
            
            p_foot = j_cum_vec(:,end);
        end
        
        %% Function!
        % *Below is modular function that can call by user if necessary*
        %% Function: Read Muscle Parameters from Animatlab
        function read_muscle_params(obj)
            num_musc = length(obj.musc_inds);
            load('AnimatlabProperties.mat');
            
%             property_names = repmat(Properties{4,1}(:,1),num_musc,1);
            property_names = Properties{4,1}(:,1);
            numProps = length(property_names);
            property_values = zeros(num_musc*numProps,1);
            property_inds = zeros(num_musc*numProps,3);
            
            %for each property
            for i=1:num_musc

                %Define its line as the "lower limit." All properties will be found
                %below it.
                lower_limit = obj.musc_inds(i);
                for j=1:numProps
                    if strcmp(property_names{j},'damping')
                        %We should find muscle damping, which is actually called "B,"
                        %just like the muscle activation's curve. Therefore we need to
                        %be a little more specific. We want the second B that shows up
                        %for a particular muscle.
                        prop_to_find = '<B>';
                        which_to_find = 2;
                    else
                        %Now find the desired property, formatted like the file.
                        prop_to_find = ['<',property_names{j},'>'];
                        which_to_find = 1;
                    end

                    %Find that string in the file, the first time after the lower limit
                    %(where the named object is found). 
                    
                    %prop_found = strfind(obj.original_text(lower_limit:end),prop_to_find);
                    prop_found = contains(obj.original_text(lower_limit:end),prop_to_find);

                    %Find the index at which this occurs, and save this for quick
                    %reference later. Remember that we were only examining
                    %original_text after the lower_limit, so we need to add that back
                    %on. -1 makes the indexing work properly.
                    
                    %temp = find(~cellfun(@isempty,prop_found)) + lower_limit - 1;
                    temp = find(prop_found) + lower_limit - 1;
                    property_inds(numProps*(i-1)+j,1) = temp(which_to_find);

                    %Find the final index of the row to keep before the number begins
                    %Number we're looking for is formatted like
                    %'<A>-0.04<A>'
                    %Index of > before the number we want
                    
                    property_inds(numProps*(i-1)+j,2) = length(prop_to_find);

                    %Find the first index of the row to keep after the number begins
                    %Index of < after number we want
                    property_inds(numProps*(i-1)+j,3) = cell2mat(strfind(obj.original_text(property_inds(numProps*(i-1)+j,1)),'</'));
                end
            end
            
            %For each property that we're changing, read in the value from
            %the current simulation
            for i=1:size(property_inds,1) %Each property
                property_values(i) = str2double(obj.original_text{property_inds(i,1)}(property_inds(i,2)+1:property_inds(i,3)-1));
            end
            
            %Make a case-insensitive list of all the muscles on this leg
            musc_names = obj.original_text(obj.musc_inds);
            for i=1:num_musc
                obj.musc_obj{i,1}.muscle_name = lower(musc_names{i}(7:end-7));
                obj.musc_obj{i,1}.muscle_index = obj.musc_inds(i);
                %Store the muscle properties in individual muscle objects
                obj.musc_obj{i}.x_off = property_values(numProps*i-9);
                obj.musc_obj{i}.ST_max = property_values(numProps*i-8);
                obj.musc_obj{i}.steepness = property_values(numProps*i-7);
                obj.musc_obj{i}.y_off = property_values(numProps*i-6);
                obj.musc_obj{i}.RestingLength = property_values(numProps*i-5);
                obj.musc_obj{i}.l_width = property_values(numProps*i-4);
                obj.musc_obj{i}.Kse = property_values(numProps*i-3);
                obj.musc_obj{i}.Kpe = property_values(numProps*i-2);
                obj.musc_obj{i}.damping = property_values(numProps*i-1);
                obj.musc_obj{i}.max_force = property_values(numProps*i);
            end
           
            disp('end')
        end
        %% Function: Write Muscle Parameters to Animatlab
        function [parameters] = write_parameters_to_animatlab(obj)
            %%% Search for the muscle index of interested!
            num_muscles = size(obj.musc_obj,1);
            parameters = cell(num_muscles,1);
            project_file = importdata(obj.proj_file);
            muscle_addresses = contains(project_file,'<Type>LinearHillMuscle</Type>');
            muscle_indices = find(muscle_addresses)-2;
            
            % Pick Leg
            name_true = 1;
            for i=2:obj.num_bodies-1
                if ~strcmp(obj.joints{i}(1:2),obj.joints{i+1}(1:2))
                    name_true = 0;
                end
            end
            
            if name_true
                leg_name = obj.joints{2}(1:2);
            else
                leg_name = input('What is the prefix of objects belonging to this leg?\nA cell array of prefixes may be entered. ');
            end
            
            if ischar(leg_name)
                musc_for_this_leg = contains(project_file(muscle_indices),['<Name>',leg_name]);
            elseif iscell(leg_name)
                musc_for_this_leg = cell(length(muscle_indices),1);
                for i=1:length(leg_name)
                    temp_musc_for_this_leg = contains(project_file(muscle_indices),['<Name>',leg_name{i}]);
                    for j=1:length(temp_musc_for_this_leg)
                        if ~isempty(temp_musc_for_this_leg{j})
                            musc_for_this_leg{j} = 1;
                        end
                    end
                end
            end            
            %Logically pick all the muscle from this leg
            muscle_indices = muscle_indices(musc_for_this_leg);
            
            %%% Parameters setting
            %Some paramters terms have placeholders because we want to write that parameter to our Matlab objects but not overwrite them in the simulation file
            parameter_terms = {'NamePlaceholder';...
                               '<B Value';...
                               '<Lwidth Value';...
                               'VmaxPlaceHolder';...
                               '<Kse Value';...
                               '<Kpe Value';...
                               '<B Value';...
                               '<D Value';...
                               '<LowerLimitScale Value';...
                               '<UpperLimitScale Value';...
                               '<RestingLength'};
            for i=1:num_muscles
                parameters{1,1} = 'Muscle name';
                parameters{i+1,1} = obj.musc_obj{i}.muscle_name;
                parameters{1,2} = 'Maximum force';
                parameters{i+1,2} = obj.musc_obj{i}.max_force;
                parameters{1,3} = 'L_width (mm)';
                parameters{i+1,3} = obj.musc_obj{i}.l_width*100;
                parameters{1,4} = 'V_max Muscle';
                parameters{i+1,4} = obj.musc_obj{i}.vmax_fiber;
                parameters{1,5} = 'Kse';
                parameters{i+1,5} = obj.musc_obj{i}.Kse;
                parameters{1,6} = 'Kpe';
                parameters{i+1,6} = obj.musc_obj{i}.Kpe;
                parameters{1,7} = 'B';
                parameters{i+1,7} = obj.musc_obj{i}.damping;
                parameters{1,8} = 'Yoffset (mN)';
                parameters{i+1,8} =obj.musc_obj{i}.y_off;
                parameters{1,9} = 'l_min (mm)';
                parameters{i+1,9} = obj.musc_obj{i}.l_min*1000;
                parameters{1,10} = 'l_max (mm)';
                parameters{i+1,10} = obj.musc_obj{i}.l_max*1000;
                parameters{1,11} = 'l_rest (mm)';
                parameters{i+1,11} = obj.musc_obj{i}.RestingLength*1000;
                for j=1:size(parameters,2)
                    lower_limit = muscle_indices(i);
                    scale = 1;
                    if j ~= 1 && j ~= 4 
                        if j == 7
                            %For some dumb reason, the creator fo Animatlab has two parameters named 'B'. In order to put damping in the right place, we have to
                            %skip over the first B.
                            prop_addresses = contains(project_file(lower_limit:end),parameter_terms{j});
                            lower_limit = find(prop_addresses,1,'first')+lower_limit;
                        end
                        if j == 9 || j == 10
                            %We need to do a similar skip over w Lmin and Lmax since they're used for two different things.
                            prop_addresses = contains(project_file(lower_limit:end),parameter_terms{j});
                            lower_limit = find(prop_addresses,1,'first')+lower_limit;
                        end
                    prop_addresses = contains(project_file(lower_limit:end),parameter_terms{j});
                    prop_index = find(prop_addresses,1,'first')+lower_limit-1;
                    %Find the line with the parameter
                    line_of_interest = project_file{prop_index};
                    if contains(line_of_interest,'None') == 1
                    else
                        if contains(line_of_interest,'milli') == 1
                            scale = 1000;
                        elseif contains(line_of_interest,'centi') == 1
                            scale = 100;
                        end
                    end
                    quotelocs = strfind(line_of_interest,'"');
                    modified_line = strcat(line_of_interest(1:quotelocs(1)),num2str(parameters{i+1,j}),line_of_interest(quotelocs(2):quotelocs(end-1)),num2str(parameters{i+1,j}/scale),line_of_interest(quotelocs(end):end));
                    %Replace the line with the modified parameter
                    project_file{prop_index} = modified_line;
                    end
                end
            end
            keyboard
            carry_on = input('You are about to overwrite the Animatlab project file you''re using with new parameters.\n This could permanently ruin your project file if there are errors.\n If this is what you want to do, type YES. Otherwise, the program will not overwrite.\n','s');
            if strcmp(carry_on,'YES')
                filename = 'C:\Users\kaiyu\Desktop\Rat Kinematics and dynamics\CalculatedPropertiesForAnimatlab.xlsx';
                xlswrite(filename,parameters);
                %file_path = strcat(obj.proj_file(1:end-6),'_2',obj.proj_file(end-5:end));
                file_path = strcat(obj.proj_file);
                fileID = fopen(file_path,'w');
                formatSpec = '%s\n';
                nrows = size(project_file);
                for row = 1:nrows
                    fprintf(fileID,formatSpec,project_file{row,:});
                end
                fclose(fileID);
            end
        end
    end
end