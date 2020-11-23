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
        theta_doubledot_motion;
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
            if to_plot 
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
            if to_plot 
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
                
%                 SC{i+1,1} = obj.musc_obj{i}.muscle_name;
%                 SC{i+1,2} = obj.musc_obj{i}.damping;
%                 SC{i+1,3} = obj.musc_obj{i}.Kse;
%                 SC{i+1,4} = obj.musc_obj{i}.Kpe;
%                 SC{i+1,5} = obj.musc_obj{i}.Kse/obj.musc_obj{i}.Kpe;
%                 SC{i+1,6} = obj.musc_obj{i}.damping/(obj.musc_obj{i}.Kse+obj.musc_obj{i}.Kpe);
                
                % Calculate muscle passive tension
                Tension{i} = obj.musc_obj{i}.compute_passive_tension(dt);
                muscle_name{i} = obj.musc_obj{i}.muscle_name;
            end
            
            telapsed = toc(tstart);
            disp(['Muscle properties stored.',' (',num2str(telapsed),'s)'])
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Save data and plot muscle passive-tension
            if to_plot 
                save ('Passive_tension.mat','muscle_name','Tension')
                Plot_Passive_tension()
                %Decide if needs to manually export Figure                
                keyboard 
                close all
            end
        end
        %% Calculate Jacobian matrix(Spatial Manipulator)
        function [Jac,p_foot] = compute_jac(obj,theta)
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
%             theta = [0; obj.theta_motion(1,:)'];   % Begin of Stance phase 
%             Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(1,:)', obj.joint_obj{3}.body_motion(1,:)', obj.joint_obj{4}.body_motion(1,:)'];
%             Cj_N  = [obj.joint_obj{2}.joint_motion(1,:)', obj.joint_obj{3}.joint_motion(1,:)', obj.joint_obj{4}.joint_motion(1,:)'];

%             theta = [0; obj.theta_motion(46,:)'];  % Begin of Swing phase 
%             Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(46,:)', obj.joint_obj{3}.body_motion(46,:)', obj.joint_obj{4}.body_motion(46,:)'];
%             Cj_N  = [obj.joint_obj{2}.joint_motion(46,:)', obj.joint_obj{3}.joint_motion(46,:)', obj.joint_obj{4}.joint_motion(46,:)'];
            
            % Initialization
            r_N = zeros(3,obj.num_bodies);         % Bodies world position
            j_N = zeros(3,obj.num_bodies-1);       % Joint world position
            omega = zeros(3,obj.num_bodies-1);     % Axes of joint rotation
            Jac = zeros(6,obj.num_bodies-1);       % Jacobian matrix
            
            %In order to compound the effects of rotating a proximal joint, we need to
            %compute the cumulative matrix product of the rotation of each of the
            %joints. z element of cum_C_joints holds the rotation to apply to each
            %joint-COM or joint-joint relative position vector.
            obj.cum_C_joints = zeros(size(obj.C_joints));
            obj.cum_C_joints(:,:,1) = eye(3)*obj.CR_bodies(:,:,1);
            
            
            for i=2:obj.num_bodies
                %This matrix describes how the limb rotates the distal segments.
                %Compute the absolute axis of rotation for the joint, in
                %the 0 position.
                u = obj.CR_bodies(:,:,i)*obj.CR_joints(:,:,i)*[-1;0;0];
                %Calculate the orientation of the axis in space, Cabs_joints, and
                %the rotation of the joint (joint angle), C_joints
                if strcmp(obj.joint_types{i},'Hinge')
                    %Now, find the rotation 
                    % C_joint here like a,b,c in previous world position
                    obj.C_joints(:,:,i) = [cos(theta(i))+u(1)^2*(1-cos(theta(i))),u(1)*u(2)*(1-cos(theta(i)))-u(3)*sin(theta(i)),u(1)*u(3)*(1-cos(theta(i)))+u(2)*sin(theta(i));u(2)*u(1)*(1-cos(theta(i)))+u(3)*sin(theta(i)),cos(theta(i))+u(2)^2*(1-cos(theta(i))),u(2)*u(3)*(1-cos(theta(i)))-u(1)*sin(theta(i));u(3)*u(1)*(1-cos(theta(i)))-u(2)*sin(theta(i)),u(3)*u(2)*(1-cos(theta(i)))+u(1)*sin(theta(i)),cos(theta(i))+u(3)^2*(1-cos(theta(i)))];
                    obj.cum_C_joints(:,:,i) = obj.cum_C_joints(:,:,i-1)*obj.C_joints(:,:,i)*obj.CR_bodies(:,:,i);
                    obj.Cabs_joints(:,:,i) = obj.cum_C_joints(:,:,i)*obj.CR_joints(:,:,i)*[1,0,0;0,cos(-theta(i)),-sin(-theta(i));0,sin(-theta(i)),cos(-theta(i))];
                elseif strcmp(obj.joint_types{i},'Prismatic')
                    % Do Nothing!
                end 
            end
            
            
            for i=2:obj.num_bodies
                %pos_bodies should be counted from 2 to end. Each column is that body's
                %position relative to the first body (which is 0s for the first).
                %Second body's position is given in the frame of the first.
                
               
                % j_N is the "World Position" in Animatlab of each joint. 
                j_N(:,i-1) = r_N(:,i-1)+obj.cum_C_joints(:,:,i-1)*obj.pos_bodies(:,i)+obj.cum_C_joints(:,:,i-1)*obj.CR_bodies(:,:,i)*obj.pos_joints(:,i);
                % r_N is the "World Position" in Animatlab of each body
                r_N(:,i) = j_N(:,i-1)-obj.cum_C_joints(:,:,i)*obj.pos_joints(:,i);
            end
            
        
            %Assuming that the COM of the last segment is halfway
            %between the distal joint and the end, add that segment to
            %it for drawing.

            if ~obj.toe_pos_known
               % Toe = ['Toe position:  ',obj.bodies{obj.num_bodies}];
               % disp(Toe)

%                 toe_pos = input('Where is the most distal point \non the foot? If it is opposite \nthe center of the distal frame \nfrom the most distal joint, simply hit [Enter].\nOtherwise, please input a vector that \ndescribes the toe location with respect \nto the ankle. This should be in \nthe global frame. ');
%                 toe_pos = [-17.899*10^-3;-85.102*10^-3;21.338*10^-3]-[-21.997*10^-3;-63.961*10^-3;21.34*10^-3];
                toe_pos = [-17.708;-73.113;-21.412]/1000 - [-22.027;-63.924;-21.424]/1000; %Local position of the toe.
                
                if isempty(toe_pos)
%                     j_N(:,end+1) = j_N(:,end) - 2*obj.Cabs_bodies(:,:,end)*obj.pos_joints(:,end);
                    obj.foot_vec = -2*obj.pos_joints(:,end);
                else
%                     j_N(:,end+1) = j_N(:,end) + toe_pos;
                    obj.foot_vec = obj.Cabs_bodies(:,:,end)'*toe_pos;
                end
                obj.toe_pos_known = 1;
            end            
            
            j_N(:,end+1) = j_N(:,end) + obj.cum_C_joints(:,:,end)*obj.foot_vec;
            
            % Stores the world positions of the joint minus the hip joint
            obj.p_st = j_N - repmat(j_N(:,1),1,size(j_N,2));
            %Hip joint
            obj.leg_attach_pt = j_N(:,1);                            
            obj.vec_len = .005;
            
             for i=1:obj.num_bodies-1
                %omega(:,i) = obj.Cabs_bodies(:,:,i)*(obj.CR_bodies(:,:,i+1)*obj.CR_joints(:,:,i+1)*[-1;0;0]);
                omega(:,i) = obj.Cabs_joints(:,:,i+1)*[-1;0;0];
                
                if strcmp(obj.joint_types{i+1},'Hinge')
                    %This is the joint twist. Crossing the axis of rotation
                    %with the member being rotated provides the twist
                    %Joint twists make up the top half of the spatial
                    %manipulator jacobian                    
                    Jac(1:3,i) = -cross(omega(:,i),j_N(:,i));
                    Jac(4:6,i) = omega(:,i);
                elseif strcmp(obj.joint_types{i+1},'Prismatic')
                    Jac(1:3,i) = omega(:,i);
                else
                    disp('Unidentified joint type.')
                end
             end
            
            p_foot = j_N(:,end);
            
        end
        %% Calculate EOM of Leg during Walking(1)
        function [Inertia,Grav] = compute_EOM1(obj,theta)
            % This section used to calculate the intertia forces and
            % gravitational forces of leg during walking.
            % The output Inertia is the inertia matrix, should times
            % acceleration be become forces.
            
            if length(theta) == obj.num_bodies - 1
                theta = [0;theta];
            elseif length(theta) == obj.num_bodies
                %Do Nothing!
            elseif isempty(theta)
                theta = zeros(obj.num_bodies,1);
            else
                disp('Orientation vector is not proper length.Should be a nx1 vector where n=num_bodies (first element is time)')
                Inertia = -1;
                Grav = -1;
                return
            end
            
            % Test sets
%             theta = [0; obj.theta_motion(1,:)'];   % Begin of Stance phase (Inertia*obj.theta_doubledot_motion(1,:)'+ Grav)
%             Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(1,:)', obj.joint_obj{3}.body_motion(1,:)', obj.joint_obj{4}.body_motion(1,:)'];
%             Cj_N  = [obj.joint_obj{2}.joint_motion(1,:)', obj.joint_obj{3}.joint_motion(1,:)', obj.joint_obj{4}.joint_motion(1,:)'];

%             theta = [0; obj.theta_motion(46,:)'];  % Begin of Swing phase (Inertia*obj.theta_doubledot_motion(46,:)'+ Grav)
%             Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(46,:)', obj.joint_obj{3}.body_motion(46,:)', obj.joint_obj{4}.body_motion(46,:)'];
%             Cj_N  = [obj.joint_obj{2}.joint_motion(46,:)', obj.joint_obj{3}.joint_motion(46,:)', obj.joint_obj{4}.joint_motion(46,:)'];

            
            % First let's get the mass of each segment in Kg
            M = zeros((obj.num_bodies-1),1); % Mass matrix
            for i = 2:obj.num_bodies
                % Find the right body
                body_found = contains(obj.original_text,['<Name>',obj.bodies{i},'</Name>']);
                next_body_ind = find(body_found,1,'first');
                chain_lower_limit = next_body_ind;
                % Searching for the mass
                mass_found = contains(obj.original_text(chain_lower_limit:end),'<Mass>');
                mass_ind = find(mass_found,1,'first') + chain_lower_limit - 1;
                mass = obj.original_text(mass_ind);
                % Extract the value 
                mass = strrep(mass,'<Mass>','');
                mass = strrep(mass,'</Mass>','');
                obj.joint_obj{i}.m = str2double(mass)/1000;
                M(i-1) = str2double(mass)/1000;
            end
            
            %%% Now calculate EOM
            r_N = zeros(3,obj.num_bodies);         % Bodies world position
            j_N = zeros(3,obj.num_bodies-1);       % Joint world position
            obj.cum_C_joints = zeros(size(obj.C_joints));
            C_def = obj.cum_C_joints;
            obj.cum_C_joints(:,:,1) = eye(3)*obj.CR_bodies(:,:,1);
            
            for i=2:obj.num_bodies
                % Rotate Matrix based on theta
                u = obj.CR_bodies(:,:,i)*obj.CR_joints(:,:,i)*[-1;0;0];
                obj.C_joints(:,:,i) = [cos(theta(i))+u(1)^2*(1-cos(theta(i))),u(1)*u(2)*(1-cos(theta(i)))-u(3)*sin(theta(i)),u(1)*u(3)*(1-cos(theta(i)))+u(2)*sin(theta(i));u(2)*u(1)*(1-cos(theta(i)))+u(3)*sin(theta(i)),cos(theta(i))+u(2)^2*(1-cos(theta(i))),u(2)*u(3)*(1-cos(theta(i)))-u(1)*sin(theta(i));u(3)*u(1)*(1-cos(theta(i)))-u(2)*sin(theta(i)),u(3)*u(2)*(1-cos(theta(i)))+u(1)*sin(theta(i)),cos(theta(i))+u(3)^2*(1-cos(theta(i)))];
                obj.cum_C_joints(:,:,i) = obj.cum_C_joints(:,:,i-1)*obj.C_joints(:,:,i)*obj.CR_bodies(:,:,i);
                C_def(:,:,i) = [sin(theta(i))*(u(1)^2-1),u(1)*u(2)*sin(theta(i))-u(3)*cos(theta(i)),u(1)*u(3)*sin(theta(i))+u(2)*cos(theta(i));u(1)*u(2)*sin(theta(i))+u(3)*cos(theta(i)),sin(theta(i))*(u(2)^2-1),u(2)*u(3)*sin(theta(i))-u(1)*cos(theta(i));u(3)*u(1)*sin(theta(i))-u(2)*cos(theta(i)),u(3)*u(2)*sin(theta(i))+u(1)*cos(theta(i)),sin(theta(i))*(u(3)^2-1)];
               
                % j_N is the "World Position" in Animatlab of each joint. 
                j_N(:,i-1) = r_N(:,i-1)+obj.cum_C_joints(:,:,i-1)*obj.pos_bodies(:,i)+obj.cum_C_joints(:,:,i-1)*obj.CR_bodies(:,:,i)*obj.pos_joints(:,i);
                % r_N is the "World Position" in Animatlab of each body
                r_N(:,i) = j_N(:,i-1)-obj.cum_C_joints(:,:,i)*obj.pos_joints(:,i);
            end
            
            % Calculate the length of the sigment
            R1 = norm(r_N(:,2)-j_N(:,1)); % Femur COM to Hip
            L1 = norm(j_N(:,2)-j_N(:,1)); % Length of Femur (Hip to knee distance)
            R2 = norm(r_N(:,3)-j_N(:,2)); % Tibia COM to Knee
            L2 = norm(j_N(:,3)-j_N(:,2)); % Length of Tibia (knee to ankle distance)
            R3 = norm(r_N(:,4)-j_N(:,3)); % foot COM to Ankle
            
            % The Inertia Matrix 
            Inertia = [ M(1)*R1^2 + M(2)*L1^2 + M(3)*L1^2, M(2)*L1*R2 + M(3)*L1*L2,  M(3)*L1*R3;...
                        M(2)*L1*R2 + M(3)*L1*L2, M(2)*R2^2 + M(3)*L2^2, M(3)*L2*R3;...
                        M(3)*L1*R3, M(3)*L2*R3, M(3)*R3^2];
                    
            Grav1 = -M(1)*9.81*obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2)...
                -M(2)*9.81*(obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2)-obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+obj.CR_bodies(:,:,3)*obj.pos_joints(:,3))+obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3))...
                -M(3)*9.81*(obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2)-obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+obj.CR_bodies(:,:,3)*obj.pos_joints(:,3))+obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)-...
                obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+obj.CR_bodies(:,:,4)*obj.pos_joints(:,4))+ obj.CR_bodies(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            Grav2 =-M(2)*9.81*obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)...
                -M(3)*9.81*(obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)-obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+obj.CR_bodies(:,:,4)*obj.pos_joints(:,4))+obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            Grav3 =-M(3)*9.81*obj.cum_C_joints(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4); 
            Grav = [Grav1(2);Grav2(2);Grav3(2)];
%             keyboard
        end
        %% Calculate EOM of Leg during Walking(2)
        function [Inertia,Int,Moment,Grav] = compute_EOM2(obj,theta)
            % This section used to calculate the intertia forces and
            % gravitational forces of leg during walking.
            % The output Inertia is the inertia matrix, should times
            % acceleration be become forces.
            
            if length(theta) == obj.num_bodies - 1
                theta = [0;theta];
            elseif length(theta) == obj.num_bodies
                %Do Nothing!
            elseif isempty(theta)
                theta = zeros(obj.num_bodies,1);
            else
                disp('Orientation vector is not proper length.Should be a nx1 vector where n=num_bodies (first element is time)')
                Inertia = -1;
                Int = -1;
                Moment = -1;
                Grav = -1;
                return
            end
            
            % Test sets
%             theta = [0; obj.theta_motion(1,:)'];   % Begin of Stance phase 
%             Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(1,:)', obj.joint_obj{3}.body_motion(1,:)', obj.joint_obj{4}.body_motion(1,:)'];

%             theta = [0; obj.theta_motion(46,:)'];  % Begin of Swing phase 
%             Cr_N  = [[0;0;0],obj.joint_obj{2}.body_motion(46,:)', obj.joint_obj{3}.body_motion(46,:)', obj.joint_obj{4}.body_motion(46,:)'];

            
            % First let's get the mass of each segment in Kg
            M = zeros((obj.num_bodies-1),1); % Mass matrix
            for i = 2:obj.num_bodies
                % Find the right body
                body_found = contains(obj.original_text,['<Name>',obj.bodies{i},'</Name>']);
                next_body_ind = find(body_found,1,'first');
                chain_lower_limit = next_body_ind;
                % Searching for the mass
                mass_found = contains(obj.original_text(chain_lower_limit:end),'<Mass>');
                mass_ind = find(mass_found,1,'first') + chain_lower_limit - 1;
                mass = obj.original_text(mass_ind);
                % Extract the value 
                mass = strrep(mass,'<Mass>','');
                mass = strrep(mass,'</Mass>','');
                obj.joint_obj{i}.m = str2double(mass)/1000;
                M(i-1) = str2double(mass)/1000;
            end
            
            %%% Now calculate EOM
            r_N = zeros(3,obj.num_bodies);         % Bodies world position
            obj.cum_C_joints = zeros(size(obj.C_joints));
            C_def = obj.cum_C_joints;
            C_dot =  C_def;
            obj.cum_C_joints(:,:,1) = eye(3)*obj.CR_bodies(:,:,1);
            
            for i=2:obj.num_bodies
                % Rotate Matrix based on theta
                u = obj.CR_bodies(:,:,i)*obj.CR_joints(:,:,i)*[-1;0;0];
                obj.C_joints(:,:,i) = [cos(theta(i))+u(1)^2*(1-cos(theta(i))),u(1)*u(2)*(1-cos(theta(i)))-u(3)*sin(theta(i)),u(1)*u(3)*(1-cos(theta(i)))+u(2)*sin(theta(i));u(2)*u(1)*(1-cos(theta(i)))+u(3)*sin(theta(i)),cos(theta(i))+u(2)^2*(1-cos(theta(i))),u(2)*u(3)*(1-cos(theta(i)))-u(1)*sin(theta(i));u(3)*u(1)*(1-cos(theta(i)))-u(2)*sin(theta(i)),u(3)*u(2)*(1-cos(theta(i)))+u(1)*sin(theta(i)),cos(theta(i))+u(3)^2*(1-cos(theta(i)))];
                obj.cum_C_joints(:,:,i) = obj.cum_C_joints(:,:,i-1)*obj.C_joints(:,:,i)*obj.CR_bodies(:,:,i);
                C_def(:,:,i) = [sin(theta(i))*(u(1)^2-1),u(1)*u(2)*sin(theta(i))-u(3)*cos(theta(i)),u(1)*u(3)*sin(theta(i))+u(2)*cos(theta(i));u(1)*u(2)*sin(theta(i))+u(3)*cos(theta(i)),sin(theta(i))*(u(2)^2-1),u(2)*u(3)*sin(theta(i))-u(1)*cos(theta(i));u(3)*u(1)*sin(theta(i))-u(2)*cos(theta(i)),u(3)*u(2)*sin(theta(i))+u(1)*cos(theta(i)),sin(theta(i))*(u(3)^2-1)];
                C_dot(:,:,i) =[cos(theta(i))*(u(1)^2-1),u(1)*u(2)*cos(theta(i))+u(3)*sin(theta(i)),u(1)*u(3)*cos(theta(i))-u(2)*sin(theta(i));u(1)*u(2)*cos(theta(i))-u(3)*sin(theta(i)),cos(theta(i))*(u(2)^2-1),u(2)*u(3)*cos(theta(i))+u(1)*sin(theta(i));u(1)*u(3)*cos(theta(i))+u(2)*sin(theta(i)),u(2)*u(3)*cos(theta(i))-u(1)*sin(theta(i)),cos(theta(i))*(u(3)^2-1)]; 
                % r_N is the "World Position" in Animatlab of each body
                r_N(:,i) = r_N(:,i-1)+ obj.cum_C_joints(:,:,i-1)*obj.pos_bodies(:,i)+obj.cum_C_joints(:,:,i-1)*(eye(3)-obj.C_joints(:,:,i))*obj.CR_bodies(:,:,i)*obj.pos_joints(:,i);       
            end
            
            %%% Velocity factor of each body (for compare with EOM1)
            R1 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2);
            L21 = obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(-obj.pos_joints(:,2)+obj.pos_bodies(:,3)+(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3));
            R2 = -obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
            L31 = L21+obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            L32 = obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(-obj.pos_joints(:,3)+obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            R3 = -obj.cum_C_joints(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
            
            %%%Double dot factor
            % dev Theta1 factorof R1
            R11 = -obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2); 
            % dev Theta1 factor of L21
            L211 = obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*(-obj.pos_joints(:,2)+obj.pos_bodies(:,3)+(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3));
            % dev Theta2 factor of L21
            L212 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
            % dev Theta1 factor of R2
            R21 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
            % dev Theta2 factor of R2
            R22 = -obj.cum_C_joints(:,:,2)*C_dot(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3); 
            % dev Theta1 factor of L31
            L311 = L211+obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            % dev Theta2 factor of L31
            L312 = L212+obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            % dev Theta3 factor of L31
            L313 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
            % dev Theta1 factor of L32
            L321 = obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(-obj.pos_joints(:,3)+obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            % dev Theta2 factor of L32
            L322 = obj.cum_C_joints(:,:,2   )*C_dot(:,:,3)*obj.CR_bodies(:,:,3)*(-obj.pos_joints(:,3)+obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4));
            % dev Theta3 factor of L32
            L323= -obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
            % dev Theta1 factor of R3
            R31 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
            % dev Theta1 factor of R3
            R32 = -obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
            % dev Theta1 factor of R3
            R33 = -obj.cum_C_joints(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
            
            
            % The Inertia Matrix
            Inertia = [M(1)*(R1'*R1) + M(2)*(L21'*L21) + M(3)*(L31'*L31), M(2)*L21'*R2 + M(3)*L31'*L32,  M(3)*L31'*R3;...
                       M(2)*L21'*R2 + M(3)*L31'*L32, M(2)*(R2'*R2) + M(3)*(L32'*L32), M(3)*L32'*R3;...
                       M(3)*L31'*R3, M(3)*L32'*R3, M(3)*(R3'*R3)];
            % Nonhomo Inertia Matrix
            Int = [2*M(1)*R1'*R11+2*M(2)*L21'*L211+2*M(3)*L31'*L311, 2*M(2)*L21'*L212+M(2)*L21'*R21+2*M(3)*L31'*L312+M(2)*R2'*L211+M(3)*L31'*L321+M(3)*L32'*L311, 2*M(3)*L31'*L313+M(3)*L31'*R31+M(3)*R3'*L311, M(2)*L21'*R22+M(2)*R2'*L212+M(3)*L31'*L322+M(3)*L32'*L312, M(3)*L31'*L323+M(3)*L32'*L313+M(3)*L31'*R32+M(3)*R3'*L312, M(3)*L31'*R33+M(3)*R3'*L313;...
                   M(2)*L21'*R21+M(2)*R2'*L211+M(3)*L31'*L321+M(3)*L32'*L311, 2*M(2)*R2'*R21+M(2)*L21'*R22+M(2)*R2'*L212+M(3)*L31'*L322+M(3)*L32'*L312+2*M(3)*L32'*L321, M(3)*L31'*L323+M(3)*L32'*L313+M(3)*R3'*L321, 2*M(2)*R2'*R22+2*M(2)*L32'*L322, 2*M(3)*L32'*L323+M(3)*L32'*R32+M(3)*R3'*L322, M(3)*L32'*R33+M(3)*R3'*L323;...
                   M(3)*L31'*R31+M(3)*R3'*L311, M(3)*L31'*R32+M(3)*R2'*L312+M(3)*L32'*R31+M(3)*R3'*L321, M(3)*L31'*R33+M(3)*R3'*L313+2*M(3)*R3'*R31, M(3)*L32'*R32+M(3)*R3'*L322, M(3)*L32'*R33+M(3)*R3'*L323+2*M(3)*R3'*R32, 2*M(3)*R3'*R33];
               
            % The Momentum Matrix
            Moment = [M(1)*R1'*R11+M(2)*L21'*L211+M(3)*L31'*L311, M(2)*R2'*L211+M(2)*L21'*R21+M(3)*L31'*L321+M(3)*L32'*L311, M(3)*R3'*L311+M(3)*L31'*R31, M(2)*R2'*R21+M(3)*L32'*L321, M(3)*R3'*L321+M(3)*L32'*R31, M(3)*R3'*R31;...
                      M(2)*L21'*L212+M(3)*L31'*L312, M(2)*R2'*L212+M(2)*L21'*R22+M(3)*L31'*L322+M(3)*L32'*L312, M(3)*R3'*L312+M(3)*L31'*R32, M(2)*R2'*R22+M(3)*L32'*L322, M(3)*R3'*L322+M(3)*L32'*R32, M(3)*R3'*R32;...
                      M(3)*L31'*L313, M(3)*L31'*L323+M(3)*L32'*L313, M(3)*R3'*L313+M(3)*L31'*R33, M(3)*L32'*L323, M(3)*R3'*L323+M(3)*L32'*R33, M(3)*R3'*R33] ;
            
            % The Gravitational Matrix          
            Grav = 9.81*[[0 1 0]*(-(M(1)+M(2)+M(3))*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2)+(M(2)+M(3))*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3))+M(3)*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)));...
                         [0 1 0]*(-(M(2)+M(3))*obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)+M(3)*obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+(eye(3)-obj.C_joints(:,:,4))*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)));...
                         [0 1 0]*(-M(3)*obj.cum_C_joints(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4))];
           
%             keyboard
        end
        %% Calculate Gravitational force of ground walking
        function [Inertia,Int,Moment,Grav] = compute_ground_EOM(obj,theta,phase)
            % Calculate EOM for ground walking is little bit different from
            % previous, as it needs to be split into two phase:
            % Stance and Swing. As the Stance phase containes pelvis into
            % the Lagrange Equation of motion and Swing phase don't
            
            if length(theta) == obj.num_bodies - 1
                theta = [0;theta];
            elseif length(theta) == obj.num_bodies
                %Do Nothing!
            elseif isempty(theta)
                theta = zeros(obj.num_bodies,1);
            else
                disp('Orientation vector is not proper length.Should be a nx1 vector where n=num_bodies (first element is time)')
                Inertia = -1;
                Int = -1;
                Moment = -1;
                Grav = -1;
                return
            end
            
            % First let's get the mass of each segment in Kg
            % 1: Pelvis; 2: Femur; 3: Tiabi; 4: Foot
            M = zeros((obj.num_bodies),1); % Mass matrix
            for i = 1:obj.num_bodies
                % Find the right body
                body_found = contains(obj.original_text,['<Name>',obj.bodies{i},'</Name>']);
                next_body_ind = find(body_found,1,'first');
                chain_lower_limit = next_body_ind;
                % Searching for the mass
                mass_found = contains(obj.original_text(chain_lower_limit:end),'<Mass>');
                mass_ind = find(mass_found,1,'first') + chain_lower_limit - 1;
                mass = obj.original_text(mass_ind);
                % Extract the value 
                mass = strrep(mass,'<Mass>','');
                mass = strrep(mass,'</Mass>','');
                M(i) = str2double(mass)/1000;
            end
            
             %%% Now calculate EOM
            r_N = zeros(3,obj.num_bodies);         % Bodies world position
            obj.cum_C_joints = zeros(size(obj.C_joints));
            C_def = obj.cum_C_joints;
            C_dot =  C_def;
            obj.cum_C_joints(:,:,1) = eye(3)*obj.CR_bodies(:,:,1);
            
            for i=2:obj.num_bodies
                % Rotate Matrix based on theta
                u = obj.CR_bodies(:,:,i)*obj.CR_joints(:,:,i)*[-1;0;0];
                obj.C_joints(:,:,i) = [cos(theta(i))+u(1)^2*(1-cos(theta(i))),u(1)*u(2)*(1-cos(theta(i)))-u(3)*sin(theta(i)),u(1)*u(3)*(1-cos(theta(i)))+u(2)*sin(theta(i));u(2)*u(1)*(1-cos(theta(i)))+u(3)*sin(theta(i)),cos(theta(i))+u(2)^2*(1-cos(theta(i))),u(2)*u(3)*(1-cos(theta(i)))-u(1)*sin(theta(i));u(3)*u(1)*(1-cos(theta(i)))-u(2)*sin(theta(i)),u(3)*u(2)*(1-cos(theta(i)))+u(1)*sin(theta(i)),cos(theta(i))+u(3)^2*(1-cos(theta(i)))];
                obj.cum_C_joints(:,:,i) = obj.cum_C_joints(:,:,i-1)*obj.C_joints(:,:,i)*obj.CR_bodies(:,:,i);
                C_def(:,:,i) = [sin(theta(i))*(u(1)^2-1),u(1)*u(2)*sin(theta(i))-u(3)*cos(theta(i)),u(1)*u(3)*sin(theta(i))+u(2)*cos(theta(i));u(1)*u(2)*sin(theta(i))+u(3)*cos(theta(i)),sin(theta(i))*(u(2)^2-1),u(2)*u(3)*sin(theta(i))-u(1)*cos(theta(i));u(3)*u(1)*sin(theta(i))-u(2)*cos(theta(i)),u(3)*u(2)*sin(theta(i))+u(1)*cos(theta(i)),sin(theta(i))*(u(3)^2-1)];
                C_dot(:,:,i) = [cos(theta(i))*(u(1)^2-1),u(1)*u(2)*cos(theta(i))+u(3)*sin(theta(i)),u(1)*u(3)*cos(theta(i))-u(2)*sin(theta(i));u(1)*u(2)*cos(theta(i))-u(3)*sin(theta(i)),cos(theta(i))*(u(2)^2-1),u(2)*u(3)*cos(theta(i))+u(1)*sin(theta(i));u(1)*u(3)*cos(theta(i))+u(2)*sin(theta(i)),u(2)*u(3)*cos(theta(i))-u(1)*sin(theta(i)),cos(theta(i))*(u(3)^2-1)]; 
                % r_N is the "World Position" in Animatlab of each body
                r_N(:,i) = r_N(:,i-1)+ obj.cum_C_joints(:,:,i-1)*obj.pos_bodies(:,i)+obj.cum_C_joints(:,:,i-1)*(eye(3)-obj.C_joints(:,:,i))*obj.CR_bodies(:,:,i)*obj.pos_joints(:,i);       
            end      
            
            if ~obj.toe_pos_known
                toe_pos = [-17.708;-73.113;-21.412]/1000 - [-22.027;-63.924;-21.424]/1000; %Local position of the toe.
                obj.foot_vec = obj.Cabs_bodies(:,:,end)'*toe_pos;
            end
            obj.toe_pos_known = 1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
            if phase == 1
                %%% EOM for Stance phase
                %%% Velocity factor of each body
                R1 = obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_joints(:,2)-obj.pos_bodies(:,3)-(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)-obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_bodies(:,4)-obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)-obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec);
                L21 = R1-obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2);
                R2 = obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_joints(:,3)-obj.pos_bodies(:,4)-obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)-obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec);
                L31 = L21+obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3));
                L32 = R2-obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
                R3 = -obj.cum_C_joints(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                L41 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                L42 = -obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                
                %%%Double dot factor
                R11 = obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_joints(:,2)-obj.pos_bodies(:,3)-(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)-obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_bodies(:,4)-obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)-obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec);
                R12 = obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_joints(:,3)-obj.pos_bodies(:,4)-obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)-obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec);
                R13 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                R22 = obj.cum_C_joints(:,:,2)*C_dot(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_joints(:,3)-obj.pos_bodies(:,4)-obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)-obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec);
                R23 = -obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                R33 = -obj.cum_C_joints(:,:,3)*C_dot(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                
                L211 = R11-obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2);  
                L311 = L211+obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3));
                L312 = R12-obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
                L322 = R22 - obj.cum_C_joints(:,:,2)*C_dot(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
                L411 = -obj.cum_C_joints(:,:,1)*C_dot(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                L412 = -obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                L422 = -obj.cum_C_joints(:,:,2)*C_dot(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec;
                  
                % The Inertia Matrix
                Inertia = [M(1)*(R1'*R1)+M(2)*(L21'*L21)+M(3)*(L31'*L31)+M(4)*(L41'*L41), M(1)*R1'*R2+M(2)*L21'*R2+M(3)*L31'*L32+M(4)*L41'*L42, M(1)*R1'*R3+M(2)*L21'*R3+M(3)*L31'*R3+M(4)*L41'*R3;...
                           M(1)*R1'*R2+M(2)*L21'*R2+M(3)*L31'*L32+M(4)*L41'*L42, M(1)*(R2'*R2)+M(2)*(R2'*R2)+M(3)*(L32'*L32)+M(4)*(L42'*L42), M(1)*R2'*R3+M(2)*R2'*R3+M(3)*L32'*R3+M(4)*L42'*R3;...
                           M(1)*R1'*R3+M(2)*L21'*R3+M(3)*L31'*R3+M(4)*L41'*R3, M(1)*R2'*R3+M(2)*R2'*R3+M(3)*L32'*R3+M(4)*L42'*R3, (M(1)+M(2)+M(3)+M(4))*(R3'*R3)];
                
                % Nonhomo Inertia Matrix
                Int = [2*M(1)*R1'*R11+2*M(2)*L21'*L211+2*M(3)*L31'*L311+2*M(4)*L41'*L411, M(1)*(3*R1'*R12+R11'*R2)+M(2)*(3*L21'*R12+L211'*R2)+M(3)*(3*L31'*L312+L311'*L32)+M(4)*(3*L41'*L412+L411'*L42), M(1)*(3*R1'*R13+R11'*R3)+M(2)*(3*L21'*R13+L211'*R3)+M(3)*(3*L31'*R13+L311'*R3)+M(4)*(3*L41'*R13+L411'*R3), (M(1)+M(2))*R12'*R2+M(1)*R1'*R22+M(2)*L21'*R22+M(3)*(L312'*L32+L322'*L31)+M(4)*(L412'*L42+L422'*L41), (M(1)+M(2))*(R13'*R2+R12'*R3)+2*M(1)*R1'*R23+2*M(2)*L21'*R23+M(3)*(R13'*L32+2*L31'*R23+L312'*R3)+M(4)*(R13'*L42+2*R23'*L41+L412'*R3), M(1)*(R13'*R3+R1'*R33)+M(2)*(R13'*R3+L21'*R33)+M(3)*(R13'*R3+L31'*R33)+M(4)*(R13'*R3+L41'*R33);...
                       M(1)*(R1'*R12+R2'*R11)+M(2)*(L21'*R12+L211'*R2)+M(3)*(L31'*L312+L32'*L311)+M(4)*(L41'*L412+L42'*L411), M(1)*(3*R2'*R12+R1'*R22)+M(2)*(3*R2'*R12+L21'*R22)+M(3)*(3*L32'*L312+L31'*L322)+M(4)*(3*L42'*L412+L41'*L422), M(1)*(R1'*R23+2*R2'*R13+R3'*R12)+M(2)*(L21'*R23+2*R2'*R13+R3'*R12)+M(3)*(L31'*R23+2*L32'*R13+R3'*L312)+M(4)*(L41'*R23+2*L42'*R13+R3'*L412), 2*(M(1)+M(2))*R2'*R22+2*M(3)*L32'*L322+2*M(4)*L42'*L422, (M(1)+M(2))*(3*R2'*R23+R3'*R22)+M(3)*(3*L32'*R23+R3'*L322)+M(4)*(3*L42'*R23+R3'*L422), (M(1)+M(2)+M(3)+M(4))*R3'*R23+(M(1)+M(2))*R2'*R33+M(3)*L32'*R33+M(4)*L42'*R33;...
                       M(1)*(R11'*R3+R1'*R13)+M(2)*(L211'*R3+L21'*R13)+M(3)*(L311'*R3+L31'*R13)+M(4)*(L411'*R3+L41'*R13), M(1)*(2*R12'*R3+R2'*R13+R1'*R23)+M(2)*(2*R12'*R3+R2'*R13+L21'*R23)+M(3)*(2*L312'*R3+L31'*R23+L32'*R13)+M(4)*(2*L412'*R3+L41'*R23+L42'*R13), 3*(M(1)+M(2)+M(3)+M(4))*R3'*R13+M(1)*R1'*R33+M(2)*L21'*R33+M(3)*L31'*R33+M(4)*L41'*R33, (M(1)+M(2))*(R22'*R3+R23'*R2)+M(3)*(L322'*R3+L32'*R23)+M(4)*(L422'*R3+L42'*R23), 3*(M(1)+M(2)+M(3)+M(4))*R3'*R23+(M(1)+M(2))*R2'*R33+M(3)*L32'*R33+M(4)*L42'*R33, 2*(M(1)+M(2)+M(3)+M(4))*R3'*R33];
                
                % The Momentum Matrix
                Moment = [M(1)*R1'*R11+M(2)*L21'*L211+M(3)*L31'*L311+M(4)*L41'*L411, M(1)*(R1'*R12+R11'*R2)+M(2)*(L21'*R12+L211'*R2)+M(3)*(L31'*L312+L311'*L32)+M(4)*(L41'*L412+L411'*L42), M(1)*(R3'*R12+R13'*R2)+M(2)*(L21'*R13+L211'*R3)+M(3)*(L31'*R13+L311'*R3)+M(4)*(L41'*R13+L411'*R3), (M(1)+M(2))*R12'*R2+M(3)*L312'*L32+M(4)*L412'*L42, (M(1)+M(2))*(R13'*R2+R12'*R3)+M(3)*(R13'*L32+L312'*R3)+M(4)*(R13'*L42+L412'*R3), (M(1)+M(2)+M(3)+M(4))*R3'*R13;...
                          M(1)*R1'*R12+M(2)*L21'*R12+M(3)*L31'*L312+M(4)*L41'*L412, M(1)*(R2'*R12+R1'*R22)+M(2)*(R2'*R12+L21'*R22)+M(3)*(L32'*L312+L31'*L322)+M(4)*(L42'*L412+L41'*L422), M(1)*(R2'*R23+R3'*R22)+M(2)*(L21'*R23+R3'*R12)+M(3)*(L31'*R23+R3'*L312)+M(4)*(L41'*R23+R3'*L412), (M(1)+M(2))*R2'*R22+M(3)*L32'*L322+M(4)*L42'*L422, (M(1)+M(2))*(R2'*R23+R3'*R22)+M(3)*(L32'*R23+R3'*L322)+M(4)*(L42'*R23+R3'*L422), (M(1)+M(2)+M(3)+M(4))*R3'*R23;...
                          M(1)*R1'*R13+M(2)*L21'*R13+M(3)*L31'*R13+M(4)*L41'*R13, M(1)*(R2'*R13+R1'*R23)+M(2)*(R2'*R13+L21'*R23)+M(3)*(L32'*R13+L31'*R23)+M(4)*(L42'*R12+L41'*R23), M(1)*(R2'*R33+R3'*R23)+M(2)*(L21'*R33+R3'*R13)+M(3)*(L31'*R33+R3'*R13)+M(4)*(L41'*R33+R3'*R13), (M(1)+M(2))*R2'*R23+M(3)*L32'*R23+M(4)*L42'*R23, (M(1)+M(2))*(R2'*R33+R3'*R23)+M(3)*(L32'*R33+R3'*R23)+M(4)*(L42'*R33+R3'*R23), (M(1)+M(2)+M(3)+M(4))*R3'*R33];
                
                % The Gravitational Matrix   
                Grav = 9.81*[[0 1 0]*(M(1)*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2)-(M(1)+M(2))*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*(obj.pos_bodies(:,3)+(eye(3)-obj.C_joints(:,:,3))*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3))-(M(1)+M(2)+M(3))*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)+obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec)-M(4)*obj.cum_C_joints(:,:,1)*C_def(:,:,2)*obj.CR_bodies(:,:,2)*obj.C_joints(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*(obj.pos_joints(:,4)+obj.foot_vec));...
                    [0 1 0]*((M(1)+M(2))*obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3)-(M(1)+M(2)+M(3))*obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*(obj.pos_bodies(:,4)+obj.CR_bodies(:,:,4)*obj.pos_joints(:,4)+obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*obj.foot_vec)-M(4)*obj.cum_C_joints(:,:,2)*C_def(:,:,3)*obj.CR_bodies(:,:,3)*obj.C_joints(:,:,4)*obj.CR_bodies(:,:,4)*(obj.pos_joints(:,4)+obj.foot_vec));...
                    [0 1 0]*(-(M(1)+M(2)+M(3)+M(4))*obj.cum_C_joints(:,:,3)*C_def(:,:,4)*obj.CR_bodies(:,:,4)*(obj.pos_joints(:,4)+obj.foot_vec))];
            end
            
            if phase == 2
                % EOM for Swing phase
                
            end
%             keyboard
        end
        %% Calculate JTF
        function [foot_force_joint_torques] = compute_JTF(obj)
            %Using a list of GRF at the foot, compute the maximum torques
            %that the joints should be able to apply.
            num_steps = 3;  %assume the data is for 3 steps
            phase_length = floor(size(obj.dt_motion,1)/(num_steps*4));
            
            %Midstance, ToeOff, MidSwing, ToeContact
            load('ToeForces.mat')
            for i=1:4 %four phases of a step
                switch i
                    case 1
                        foot_forces = ToeForce_MidStance;
                    case 2
                        foot_forces = ToeForce_ToeOff;
                    case 3
                        foot_forces = ToeForce_MidSwing;
                    case 4
                        foot_forces = ToeForce_ToeContact;
                end
                foot_forces = foot_forces/9.81;
                num_foot_forces = size(foot_forces,1);
                foot_force_joint_torques = zeros(obj.num_bodies-1,num_foot_forces);
                n = 1;
                
                current_config = obj.theta_motion((i)*phase_length,:)';
                
                %Input the joint angles that we just computed, and
                %insert them into the configuration of the leg.
                [J,foot_position(:,i)] = obj.compute_jac(current_config);

                for k=1:num_foot_forces
                    
                    force_vec(:,k) = -[foot_forces(k,:)';0];
                    
                    %Calculate J'*F to find the joint torques. J is in
                    %the body frame, so F also needs to be the wrench
                    %acting on the most proximal joint.

                    foot_wrench = [force_vec(:,k);cross(foot_position(:,i),force_vec(:,k))];
                    current_torques = J'*foot_wrench;
                    foot_force_joint_torques(:,n) = current_torques;
                    n = n + 1;
                end
            end
            plotting = 1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if plotting 
                load('Kinematics and Dynamics.mat')
                current_size1 = size(foot_force_joint_torques,2);
                current_size2 = size(TorqueAll,1);
                num_of_samples = 100;
                
                foot_force_joint_torques = interp1(linspace(0,1,current_size1),foot_force_joint_torques',linspace(0,1,num_of_samples));
                TorqueAll =  interp1(linspace(0,1,current_size2),TorqueAll,linspace(0,1,num_of_samples));
                
                subplot(2,3,1)
                plot(foot_force_joint_torques(:,1),'b-','Linewidth',1.5)
                title('Hip JTF','FontSize',9)
                subplot(2,3,2)
                plot(foot_force_joint_torques(:,2),'r-','Linewidth',1.5)
                title('Knee JTF','FontSize',9)
                subplot(2,3,3)
                plot(foot_force_joint_torques(:,3),'m-','Linewidth',1.5)
                title('Ankle JTF','FontSize',9)
                subplot(2,3,4)
                plot(TorqueAll(:,1),'b-','Linewidth',1.5)
                title('Hip Torque','FontSize',9)
                subplot(2,3,5)
                plot(TorqueAll(:,2),'r-','Linewidth',1.5)
                title('Knee Torque','FontSize',9)
                subplot(2,3,6)
                plot(TorqueAll(:,3),'m-','Linewidth',1.5)
                title('Ankle Torque','FontSize',9)

                keyboard
            end
        end
        %% Calculate MN activation during walking
        function compute_MN_act_for_motion(obj,to_plot)
            keyboard
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
        %% Function: Plot rotate axis during motion
        function plot_Omega(obj)
            % Plot bodies and joints motion during walking
            load('motion.mat')
            num_steps = length(obj.dt_motion);
            
            
            % Position of the joints
            xh = pos_hip(:,1);
            yh = pos_hip(:,2);
            zh = pos_hip(:,3);
            xk = pos_knee(:,1);
            yk = pos_knee(:,2);
            zk = pos_knee(:,3);
            xa = pos_ankle(:,1);
            ya = pos_ankle(:,2);
            za = pos_ankle(:,3);
            
            % Posion of the bodies (COM)
            x1 = pos_femur(:,1);
            y1 = pos_femur(:,2);
            z1 = pos_femur(:,3);
            x2 = pos_tibia(:,1);
            y2 = pos_tibia(:,2);
            z2 = pos_tibia(:,3);
            x3 = pos_foot(:,1);
            y3 = pos_foot(:,2);
            z3 = pos_foot(:,3);
            
            % Initailizing the Jac(omega) and toe position
            toe_position = zeros(num_steps,3);
            
            
            for i = 1:num_steps
                % Calculate Jac and toe position
                [J,toe] = obj.compute_jac(obj.theta_motion(i,:)');
                toe_position(i,:) = toe';
                xt = toe(1);
                yt = toe(2);
                zt = toe(3);
                
                for j=1:obj.num_bodies-1
                    omega(:,j) = J(4:6,j);
                    W(:,j) = obj.joint_obj{j+1}.joint_motion(i,:)'+omega(:,j)*obj.vec_len;
                end

                
                
                % Plot 3D motion
                plot3([zh(i),zk(i)],[xh(i),xk(i)],[yh(i),yk(i)],'k-','Linewidth',2);
                hold on
                plot3([zk(i),za(i)],[xk(i),xa(i)],[yk(i),ya(i)],'k-','Linewidth',2);
                hold on
                plot3([za(i),zt],[xa(i),xt],[ya(i),yt],'k-','Linewidth',2);
                hold on
                plot3([zh(i),W(3,1)],[xh(i),W(1,1)],[yh(i),W(2,1)],'c--','Linewidth',1);
                hold on
                plot3([zk(i),W(3,2)],[xk(i),W(1,2)],[yk(i),W(2,2)],'c--','Linewidth',1);
                hold on
                plot3([za(i),W(3,3)],[xa(i),W(1,3)],[ya(i),W(2,3)],'c--','Linewidth',1);
                hold off
                
                axis equal
                axis([-0.03 0 -0.05 0.05 -0.07 0])
                pbaspect([1 1 1])
%                 view([-120 30])
                view([-90 0])
                set(gcf,'Position',[400 300 700 600])  
                grid on
                
                % COM
                text(z1(i),x1(i),y1(i),'o','color','r');
                text(z2(i),x2(i),y2(i),'o','color','r');
                text(z3(i),x3(i),y3(i),'o','color','r');
                
                % Omega
                text(W(3,1),W(1,1),W(2,1),'x','color','g');
                text(W(3,2),W(1,2),W(2,2),'x','color','g');
                text(W(3,3),W(1,3),W(2,3),'x','color','g')

                L(i) = getframe;               
            end
            movie(L,5,24)
            keyboard
        end
        %% Function: Plot EOM
        function plot_EOM(obj,to_plot)
            % THis is for compare of two different method, EOM1 and EOME2
            % and plot the difference
            
            to_plot = 2;
            
            % Initialization 
            num_steps = length(obj.dt_motion);
            Sizer = zeros(num_steps,3);
            IF = Sizer;
            GF = Sizer;
            INF = Sizer;
            MF = Sizer;
            
            for i = 1:num_steps
                
                NT = [obj.theta_dot_motion(i,1)^2;obj.theta_dot_motion(i,1)*obj.theta_dot_motion(i,2);obj.theta_dot_motion(i,1)*obj.theta_dot_motion(i,3);obj.theta_dot_motion(i,2)^2;obj.theta_dot_motion(i,2)*obj.theta_dot_motion(i,3);obj.theta_dot_motion(i,3)^2];
                
%                 %For EOM1
%                 [Inertia,Grav] = obj.compute_EOM1(obj.theta_motion(i,:)');
%                 IF1(i,:) = Inertia*obj.theta_doubledot_motion(i,:)';
%                 GF1(i,:) = Grav;
% 
%                 %For EOM2
%                 [Inertia,Int,Moment,Grav] = obj.compute_EOM2(obj.theta_motion(i,:)');
%                 IF2(i,:) = Inertia*obj.theta_doubledot_motion(i,:)';
%                 INF(i,:) = Int*NT;
%                 MF(i,:) = Moment*NT;
%                 GF2(i,:) = Grav;
                
                %For Ground walking
                if i <= num_steps/2
                    [Inertia,Int,Moment,Grav] = obj.compute_ground_EOM(obj.theta_motion(i,:)',1);
                else
                    [Inertia,Int,Moment,Grav] = obj.compute_ground_EOM(obj.theta_motion(i,:)',1);
                end
                    IF(i,:) = Inertia*obj.theta_doubledot_motion(i,:)';
                    INF(i,:) = Int*NT;
                    MF(i,:) = Moment*NT;
                    GF(i,:) = Grav;
            end
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Save data and plot muscle passive-tension
            T = {'Hip';'Knee';'Ankle'};
            
            if to_plot ==1
                figure(1)
                for i = 1:obj.num_bodies-1
                    subplot(1,3,i)
                    plot(IF(:,i),'b','Linewidth',1.5)
                    title(T{i});
                end
                suptitle('Intertia Force interms of Theta double dot')
                
                figure(2)
                for i = 1:obj.num_bodies-1
                    subplot(1,3,i)
                    plot(INF(:,i),'b','Linewidth',1.5)
                    title(T{i});
                end
                suptitle('Intertia Force interms of Theta double multiple ')
                
                figure(3)
                for i = 1:obj.num_bodies-1
                    subplot(1,3,i)
                    plot(MF(:,i),'b','Linewidth',1.5)
                    title(T{i});
                end
                suptitle('Momentum force')
                
                figure(4)
                for i = 1:obj.num_bodies-1
                    subplot(1,3,i)
                    plot(GF(:,i),'b','Linewidth',1.5)
                    title(T{i});
                end
                suptitle('Gravitational Force ')
            end
            
            if to_plot ==2
                for i = 1:obj.num_bodies-1
                    subplot(1,3,i)
                    plot(IF(:,i)+INF(:,i)+GF(:,i)-MF(:,i),'b','Linewidth',1.5)
                    hold off
                    title(T{i});
                end
                suptitle('EOM force during walking')
            end
            
%             keyboard
%             EOM = IF2+INF+GF2-MF;
%             save ('EOM.mat','EOM')
        end
        %% Function: Plot walking motion for two Leg walking
        function plot_walking_motion(obj,to_plot)  
%             to_plot = 1;
            %%% Initialization
            sizer = cell(3,1);
            LH_joints = sizer;
            LH_bodies = sizer;
            RH_joints = sizer;
            RH_bodies = sizer;
            step_num = length(obj.dt_motion);
            Stance2Swing = step_num/2;
            
            %%% Record motion data
            
            % Joints motion
            LH_motion = obj.theta_motion;
            RH_motion = circshift(LH_motion,Stance2Swing,1);
            
            % Joints and bodies position during walking
            for i = 1: obj.num_bodies-1
                LH_joints{i} = obj.joint_obj{i+1}.joint_motion;
                LH_bodies{i} = obj.joint_obj{i+1}.body_motion;
                RH_joints{i} = circshift(LH_joints{i},Stance2Swing,1);
                RH_bodies{i} = circshift(LH_bodies{i},Stance2Swing,1);
            end
            
            if to_plot ==1
                
                % Plot joint motion 
                figure (1)
                T = {'Hip';'Knee';'Ankle'};
                for i = 1: obj.num_bodies-1
                    subplot(1,3,i)
                    plot(LH_motion(:,i),'b','Linewidth',2)
                    hold on 
                    plot(RH_motion(:,i),'r','Linewidth',2)
                    hold off
                    title(T{i})
                end
                set(gcf,'Position',[200 300 1200 400])
                suptitle('joint motions for hind limb')

                
                % Plot Walking motion
                figure(2)
                for i = 1:step_num
                    % Calculate Jac and toe position
                    [~,LH_toe] = obj.compute_jac(LH_motion(i,:)');
                    [~,RH_toe] = obj.compute_jac(RH_motion(i,:)');
                    
                    % Plot walking motion 
                    plot([LH_joints{1}(i,1),LH_joints{2}(i,1)],[LH_joints{1}(i,2),LH_joints{2}(i,2)],'k-','Linewidth',2);
                    hold on
                    plot([RH_joints{1}(i,1),RH_joints{2}(i,1)],[RH_joints{1}(i,2),RH_joints{2}(i,2)],'r-','Linewidth',2);
                    hold on
                    plot([LH_joints{2}(i,1),LH_joints{3}(i,1)],[LH_joints{2}(i,2),LH_joints{3}(i,2)],'k-','Linewidth',2);
                    hold on
                    plot([RH_joints{2}(i,1),RH_joints{3}(i,1)],[RH_joints{2}(i,2),RH_joints{3}(i,2)],'r-','Linewidth',2);
                    hold on
                    plot([LH_joints{3}(i,1),LH_toe(1)],[LH_joints{3}(i,2),LH_toe(2)],'k-','Linewidth',2);
                    hold on
                    plot([RH_joints{3}(i,1),RH_toe(1)],[RH_joints{3}(i,2),RH_toe(2)],'r-','Linewidth',2);
                    hold off
                    
                    title('Waling motion(Side View)')
                    
                    axis([-0.1 0.1 -0.08 0])
                    axis equal
                    
                    hline = refline([0 -0.0618]);
                    hline.Color = 'k';
                    
                    set(gca,'XDir','reverse');
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                    set(gcf,'Position',[500 300 800 600])
                    

                    L(i) = getframe;              
                end
%                 movie(L,5,30)
                close
                
                % Plot for walking cartoon
                figure (3)
                Phase = 1:step_num/4:step_num;
                for i = 1:4
                    j = Phase(i);
                    
                    % Calculate Jac and toe position
                    [~,LH_toe] = obj.compute_jac(LH_motion(j,:)');
                    [~,RH_toe] = obj.compute_jac(RH_motion(j,:)');
                    
                    % store the new joint positon
                    basline = -0.0618;
                    Lhx = LH_joints{1}(j,1) ;
                    Lhy = LH_joints{1}(j,2) - basline;
                    Lkx = LH_joints{2}(j,1);
                    Lky = LH_joints{2}(j,2) - basline;
                    Lax = LH_joints{3}(j,1);
                    Lay = LH_joints{3}(j,2) - basline;
                    Ltx =  LH_toe(1);
                    Lty = LH_toe(2) - basline;
                    
                    Rhx = RH_joints{1}(j,1) ;
                    Rhy = RH_joints{1}(j,2) - basline;
                    Rkx = RH_joints{2}(j,1);
                    Rky = RH_joints{2}(j,2) - basline;
                    Rax = RH_joints{3}(j,1);
                    Ray = RH_joints{3}(j,2) - basline;
                    Rtx =  RH_toe(1);
                    Rty = RH_toe(2) - basline;
                    
                    %plot ground motion
                    subplot(1,4,i)
                    plot([Lhx,Lkx],[Lhy,Lky],'k-','Linewidth',2);
                    hold on
                    plot([Rhx,Rkx],[Rhy,Rky],'r-','Linewidth',2);
                    hold on
                    plot([Lkx,Lax],[Lky,Lay],'k-','Linewidth',2);
                    hold on
                    plot([Rkx,Rax],[Rky,Ray],'r-','Linewidth',2);
                    hold on
                    plot([Lax,Ltx],[Lay,Lty],'k-','Linewidth',2);
                    hold on
                    plot([Rax,Rtx],[Ray,Rty],'r-','Linewidth',2);
                    hold off
                    
%                     axis equal
                    axis([-0.05 0.05 0 0.08])
                    pbaspect([1 1 1])
                    
                    set(gca,'XDir','reverse');
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    set(gca,'YColor','none')
                    box off
                end
                keyboard
            end  
        end
        %% Function: Plot passive tension around walking phases
        function plot_phase_passive_tension(obj,to_plot)
%             to_plot = 1;
            
            step_num = length(obj.dt_motion);
            musc_num = length(obj.musc_inds);
            
            PT = zeros(step_num,musc_num);
            holder = zeros(20,8);
            
            for i = 1:musc_num
                PT(:,i) = obj.musc_obj{i}.passive_tension_profile;
            end
            
            PT = [PT;PT];
            Phase = [100,25,50,75];
            
            if to_plot ==1
                MN = {'BFA','BFP','VA','GA','TA','SO','IP','RF'};
                PN = {'Touch down','Mid stance','Toe off', 'Mid swing'};
                for i = 1:4
%                     subplot (1,4,i)
                    st = Phase(i);
                    t=0:1:20;
                    
                    for j = 1:musc_num
                        
                        if i ==1
                            PT(st-10:st+10,j) = smooth(PT(st-10:st+10,j));
                        end
                        
                        subplot (8,4,i+4*(j-1))
                        plot(t,PT(st-10:st+10,j));
  
                        xticks([0 10 20])
                        set(gca,'xticklabel',[]);
                        xlim auto
                        set(gca,'YTick',[]);
%                         yMin = round(min(PT(st-10:st+10,j)),2);
%                         yMax = round(max(PT(st-10:st+10,j)),2);
%                         yticklabels([yMin yMax])
%                         ylim auto
                        if i ==1 
                            ylabel(MN{j},'fontsize', 11)
                        end
                    end
                    set(gca, 'xticklabel', {'',PN{i},''}, 'fontsize', 12);
                end

            end
           
        end
        %% Function: Plot walking motion for ground walking
        function plot_ground_walking_motion(obj,to_plot)
            % For the ground walking, we need split into two phase to
            % calculate the walking motion and Gravitational Force.I.E.the
            % stance and swing phase as two groups. the stance phase
            % calculate the jacobian of foot position as base line. For the
            % swing phase we use end of stance phase pelvis postion to
            % calculate the swing motion.
                     
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step_num = length(obj.dt_motion);
            Stance2Swing = step_num/2;
            toe_position = zeros(step_num,3);
            
            %%% Calculate toe_position during walking
            for i = 1:step_num
                [~,toe_position(i,:)] = obj.compute_jac(obj.theta_motion(i,:)');
            end
            
            %%% For Stance phase, the toe is on the ground, so we set
            %%% toe_position as basepoint.
            
            
            % Calculate bodies and joints postion during Stance phase
            ST_toe = toe_position(1:Stance2Swing,:);
            ST_joint_motion = cell(obj.num_bodies,1);
            ST_joint_motion{4} = zeros(Stance2Swing,3);
            ST_body_motion = cell(obj.num_bodies,1);
            ST_body_motion{1} = zeros(Stance2Swing,3) - ST_toe;
            
            for i = 1:obj.num_bodies-1
                ST_joint_motion{i} = obj.joint_obj{i+1}.joint_motion(1:Stance2Swing,:) - ST_toe;
                ST_body_motion{i+1} = obj.joint_obj{i+1}.body_motion(1:Stance2Swing,:) - ST_toe;
            end
            
            %%% For Swing phase, the pelvis is fixed at air, so we set the
            %%% pelvis as the basepoint.
            
            % Record the pelvis position for last phase
%             Pelvis_position = ST_body_motion{1}(end,:);
            Pelvis_position(:,1) = ST_body_motion{1}(:,1) - ST_body_motion{1}(1,1) + ST_body_motion{1}(end,1);
            Pelvis_position(:,2) = ST_body_motion{1}(:,2);
            Pelvis_position(:,3) = flipud(ST_body_motion{1}(:,3));
           
            
            % Calculate bodies and joints postion during Swing phase
            SW_toe = toe_position(Stance2Swing+1:step_num,:);
            SW_joint_motion = cell(obj.num_bodies,1);
            SW_joint_motion{4} = zeros(Stance2Swing,3) + Pelvis_position + SW_toe;
            SW_body_motion = cell(obj.num_bodies,1);
            SW_body_motion{1} = zeros(Stance2Swing,3) + Pelvis_position;
            
            for i = 1:obj.num_bodies-1
                SW_joint_motion{i} = obj.joint_obj{i+1}.joint_motion(Stance2Swing+1:step_num,:) + Pelvis_position;
                SW_body_motion{i+1} = obj.joint_obj{i+1}.body_motion(Stance2Swing+1:step_num,:) + Pelvis_position;
            end
            
            %%% Combine Stance and Swing Phase
            for i = 1:obj.num_bodies
                body_motion{i} = [ST_body_motion{i};SW_body_motion{i}];
                joint_motion{i} = [ST_joint_motion{i};SW_joint_motion{i}];                
            end
            
            % Calculate RH_motion
            for i = 1:obj.num_bodies -1
                RH_motion{i}(1:Stance2Swing,:) = obj.joint_obj{i+1}.joint_motion(Stance2Swing+1:end,:) - ST_toe;
                RH_motion{i}(Stance2Swing+1:step_num,:) = obj.joint_obj{i+1}.joint_motion(1:Stance2Swing,:) + Pelvis_position;
            end
            RH_motion{4}(1:Stance2Swing,:) = SW_toe- ST_toe;
            RH_motion{4}(Stance2Swing+1:step_num,:) = ST_toe  + Pelvis_position;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if to_plot == 1
                % Plot decesion: 0.Nothing; 1.All phase; 2.Stance Phase; 3.Swing phase

                Chose_plot = 1;

                if Chose_plot ==1                    
                    for i = 1:step_num
                        plot([joint_motion{1}(i,1),joint_motion{2}(i,1)],[joint_motion{1}(i,2),joint_motion{2}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([RH_motion{1}(i,1),RH_motion{2}(i,1)],[RH_motion{1}(i,2),RH_motion{2}(i,2)],'r-','Linewidth',2);
                        hold on
                        plot([joint_motion{2}(i,1),joint_motion{3}(i,1)],[joint_motion{2}(i,2),joint_motion{3}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([RH_motion{2}(i,1),RH_motion{3}(i,1)],[RH_motion{2}(i,2),RH_motion{3}(i,2)],'r-','Linewidth',2);
                        hold on
                        plot([joint_motion{3}(i,1),joint_motion{4}(i,1)],[joint_motion{3}(i,2),joint_motion{4}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([RH_motion{3}(i,1),RH_motion{4}(i,1)],[RH_motion{3}(i,2),RH_motion{4}(i,2)],'r-','Linewidth',2);
                        hold off
                        
                        axis equal
                        axis([-0.05 0.16 0 0.08])
                        
                        set(gca,'XDir','reverse');
                        set(gcf,'Position',[500 300 800 400])
                        
                        text(joint_motion{1}(i,1),joint_motion{1}(i,2),'x','color','r');
                        text(body_motion{1}(i,1),body_motion{1}(i,2),'o','color','b')
                        text(joint_motion{2}(i,1),joint_motion{2}(i,2),'x','color','r');
                        text(body_motion{2}(i,1),body_motion{2}(i,2),'o','color','b');
                        text(joint_motion{3}(i,1),joint_motion{3}(i,2),'x','color','r');
                        text(body_motion{3}(i,1),body_motion{3}(i,2),'o','color','b');
                        text(body_motion{4}(i,1),body_motion{4}(i,2),'o','color','b');
                        
                        W(i) = getframe;
                    end
                    movie(W,5,30)
                end
                
                if Chose_plot == 2
                    for i = 1:Stance2Swing
                        plot([ST_joint_motion{1}(i,1),ST_joint_motion{2}(i,1)],[ST_joint_motion{1}(i,2),ST_joint_motion{2}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([ST_joint_motion{2}(i,1),ST_joint_motion{3}(i,1)],[ST_joint_motion{2}(i,2),ST_joint_motion{3}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([ST_joint_motion{3}(i,1),ST_joint_motion{4}(i,1)],[ST_joint_motion{3}(i,2),ST_joint_motion{4}(i,2)],'k-','Linewidth',2);
                        hold off
                        
                        axis equal
                        axis([-0.06 0.06 0 0.09])
                        
                        set(gca,'XDir','reverse');
                        
                        ST(i) = getframe;
                    end
                    movie(ST,5,30)
                end
                
                if Chose_plot == 3
                    for i = 1:Stance2Swing
                        plot([SW_joint_motion{1}(i,1),SW_joint_motion{2}(i,1)],[SW_joint_motion{1}(i,2),SW_joint_motion{2}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([SW_joint_motion{2}(i,1),SW_joint_motion{3}(i,1)],[SW_joint_motion{2}(i,2),SW_joint_motion{3}(i,2)],'k-','Linewidth',2);
                        hold on
                        plot([SW_joint_motion{3}(i,1),SW_joint_motion{4}(i,1)],[SW_joint_motion{3}(i,2),SW_joint_motion{4}(i,2)],'k-','Linewidth',2);
                        hold off
                        
                        axis equal
                        axis([-0.01 0.16 0 0.09])
                        
                        set(gca,'XDir','reverse');
                        
                        
                        SW(i) = getframe;
                    end
                    movie(SW,5,30)
                end
                keyboard
            end
        end
    end
end