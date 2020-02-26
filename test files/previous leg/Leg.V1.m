classdef Leg < matlab.mixin.SetGet
    %JACOBIAN Create a Jacobian manipulator object from an AnimatLab
    %simulation file.
    %   Detailed explanation goes here
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
        Cabs_joints; %Absolute orientation of the joints
        C_joints; %Rotation of the joints as per theta
        CR_bodies; %Relative rotation of the bodies with respect to the proximal one
        CR_joints; %Relative rotation of the joints with respect to the local frame
        cum_C_joints; %Cumulative rotation (product) of proximal joints
        Cabs_bodies; %Absolute rotation of the bodies
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
        %%  Initialize some variables. 
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
                            disp('End of chain reached. Leg body and joint constucted.')
                        end
                    end
                end
                i = i + 1;
            end
            
        %% Find muscles, their attachments, and associate them with the proper bodies and joints.
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
            
            %find indices of the names of these attachments.
            for i=1:length(obj.musc_inds)
                for j=1:length(attach_to_find{i})
                    %save the names
                    id_loc = contains(obj.original_text,attach_to_find{i}{j});
                    id_ind = find(id_loc,1,'first');
                    attach_names{i}{j} = id_ind - 1;  
                     
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
                    % check if the attachment already been used.
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
            
            %Set the size of obj.pos_attachments for each body
            obj.pos_attachments{2} = cell(pelvis_count,3);
            obj.pos_attachments{3} = cell(femur_count,3);
            obj.pos_attachments{4} = cell(tibia_count,3);
            obj.pos_attachments{5} = cell(foot_count,3);
            
            %Store information about attachment points for each body in
            %obj.pos_attachments. Resultant cell array holds position
            %information for XYZ and name strings of every unique
            %attachment point on each body
            
            %Store information about attachment points for each body in
            %obj.pos_attachments. Resultant cell array holds position
            %information for XYZ and name strings of every unique
            %attachment point on each body
            for j=1:length(attach_names)
                obj.musc_obj{j}.muscle_index = musc_name_inds(j);
                for k=1:length(attach_names{j})
                    if isequal(size(attach_locs{j}{k}),[3,1])
                        temp_name_str = char(obj.original_text(attach_names{j}{k}));
                        temp_name_str = strrep(temp_name_str,'<Name>','');
                        temp_name_str = strrep(temp_name_str,'</Name>','');

                        if contains(temp_name_str,' Pelvis ')
                            body = 2;
                        elseif contains(temp_name_str,' Femur ')
                            body = 3;
                        elseif contains(temp_name_str,' Tibia ')
                            body = 4;
                        elseif contains(temp_name_str,' Foot ')
                            body = 5;
                        else
                            keyboard
                            error('Attachment not labeled with a body part')
                        end
                        %Following if statement necessary to save
                        %information only for non-redundant attachment
                        %points
                        obj.musc_obj{j,1}.pos_attachments{k,1} = attach_locs{j}{k};
                        num_attach = strcat(num2str(size(obj.musc_obj{j,1}.pos_attachments{k,1},2)),') ');
                        obj.musc_obj{j,1}.pos_attachments{k,2} = [num_attach,obj.original_text{attach_names{j}{k}}(7:end-7),'\n'];
                        obj.musc_obj{j,1}.pos_attachments{k,3} = body-1;
                        obj.musc_obj{j,1}.pos_attachments{k,4} = [];
                        if ~isempty(find(used_indices == attach_names{j}{k},1,'first'))
                            used_indices(find(used_indices == attach_names{j}{k},1,'first')) = 0;
                            jj = sum(~cellfun(@isempty,obj.pos_attachments{body,1}(:,1)))+1;
                            obj.pos_attachments{body,1}{jj,1} = [obj.pos_attachments{body,1}{jj,1},attach_locs{j}{k}];
                            num_attach = strcat(num2str(size(obj.pos_attachments{body,1}{jj,1},2)),') ');
                            obj.pos_attachments{body,1}{jj,2} = [obj.pos_attachments{body,1}{jj,2},num_attach,obj.original_text{attach_names{j}{k}}(7:end-7),'\n'];
                        else
                            %skip this reused attachment
                        end
                    end  
                end
            end
            
            keyboard
        end
       
    end
end