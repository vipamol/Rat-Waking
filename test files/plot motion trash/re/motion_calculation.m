function motion_calculation(obj)

% Calculate uu_joint, the rotatation matrix of joint axis
for i = 1:3
uu_joint{i} = obj.CR_bodies(:,:,i+1)*obj.CR_joints(:,:,i+1)*[-1;0;0];
end

% Calculate joints and bodies location during walking.
for i = 1:length(obj.dt_motion)
    a = obj.joint_obj{2}.axis_angle_rotation(i,uu_joint{1}); % axis rotation of hip
    b = obj.joint_obj{3}.axis_angle_rotation(i,uu_joint{2}); % axis rotation of knee
    c = obj.joint_obj{4}.axis_angle_rotation(i,uu_joint{3}); % axis rotation of ankle
 
    pos_hip(i,:) = obj.Cabs_bodies(:,:,1)*obj.pos_bodies(:,2)+ obj.Cabs_bodies(:,:,2)*obj.pos_joints(:,2);
    pos_femur(i,:) = pos_hip(i,:)' - obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*obj.pos_joints(:,2);
    
    pos_knee(i,:) = pos_femur(i,:)' + obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*obj.pos_bodies(:,3)+obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
    pos_tibia(i,:) = pos_knee(i,:)' - obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*b*obj.CR_bodies(:,:,3)*obj.pos_joints(:,3);
    
    pos_ankle(i,:) = pos_tibia(i,:)' + obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*b*obj.CR_bodies(:,:,3)*obj.pos_bodies(:,4)+obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*b*obj.CR_bodies(:,:,3)*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
    pos_foot(i,:) = pos_ankle(i,:)' - obj.CR_bodies(:,:,1)*a*obj.CR_bodies(:,:,2)*b*obj.CR_bodies(:,:,3)*c*obj.CR_bodies(:,:,4)*obj.pos_joints(:,4);
   
end

dt = obj.dt_motion;
save ('motion.mat','dt','pos_hip','pos_femur','pos_knee','pos_tibia','pos_ankle','pos_foot')
end

