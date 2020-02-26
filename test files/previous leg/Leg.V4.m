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
        original_text; %content of the .asim file
        joint_types; %type of joints in the leg, as read from the .asim file
        num_bodies; %number of bodies in the leg, as given by the user's list
        organism_name; %name of the organism in the simulation file that we are interested in
        bodies; %names of the segments in the leg
        joints; %names of the joints in the leg
        musc_inds;
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
        joint_obj; %cell array of joint objects (for each joint in the leg)
        euler_angs_joints;
        euler_angs_bodies;
        theta_range; %min and max rotation for each joint
        velocity_range;
        proj_file; %project file
        p_st; %Position of the body relative to the foot on the ground.
        leg_attach_pt;
        foot_pos_residual;
        control_current; %Stimuli to joint position controllers necessary for each "config."
        body_twists; %Body frame transformations desired from a standing still posture. These correspond to "config."
        body_thetas;
        toe_pos_known = 0; %Boolean about whether the toe position is known or not.
        foot_vec; %vector that describes the toe wrt the most distal joint. Important for calculating forces from the foot.
        ankle_factor; %multiply distal max torque by 3 if it operates like an ankle.
    
        musc_obj;%cell array of muscle objects (for each muscle in the leg)
        
        torque_motion;
        theta_motion;
        theta_dot_motion;
        dt_motion;
        
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
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.CR_bodies(:,:,1)*(obj.pos_bodies(:,2) + ...
                                a*obj.CR_bodies(:,:,2)*obj.musc_obj{j}.pos_attachments{k,1});
                        elseif body == 3
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.CR_bodies(:,:,1)*(obj.pos_bodies(:,2) + ...
                                a*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+b*obj.CR_bodies(:,:,3)*...
                                obj.musc_obj{j}.pos_attachments{k,1}));
                        elseif body == 4    
                            obj.musc_obj{j}.pos_attachments{k,4}(i,:) = obj.CR_bodies(:,:,1)*(obj.pos_bodies(:,2) + ...
                                a*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+b*obj.CR_bodies(:,:,3)*((obj.pos_bodies(:,4) + ...
                                c*obj.CR_bodies(:,:,4)*obj.musc_obj{j}.pos_attachments{k,1}))));
                        end                       
                    end
                    % Calculate muscle length during walking
                    obj.musc_obj{j}.muscle_length_profile(i,:) = obj.musc_obj{j}.muscle_length(i);
                end                
            end
            % Calculate mucle velocity
            dt = obj.dt_motion(2)-obj.dt_motion(1);
            for i = 1:length(obj.musc_obj)
                obj.musc_obj{i}.muscle_velocity_profile = diff(obj.musc_obj{i}.muscle_length_profile)/dt;
            end
            
            telapsed = toc(tstart);
            disp(['Muscle walking profile stored.',' (',num2str(telapsed),'s)'])
            
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
    end
end