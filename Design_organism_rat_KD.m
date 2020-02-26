clear workspace
close all
clc

%Decide whether to load a previously made animal or load a new one.
%0- Do Nothing
%1- Clear out all animal information and make a new one
%2- Load in an animal object file based on AniLoc described below
%3- Change the location of the project file within the animal object
%(useful for reading in and using new synapse parameters without making an
%entire new animal)
makerat = 1;

%Text for describing the animal object that is made and trained this time.
%This text is added to the save file name.
Text_append = '_InverseKinematics';
Text_to_append = '_InverseKinematics';
vers = '2020';
RatDate = '29-Mar';
RatJoint = 'LH_Leg';
RatLoc = ['C:\Users\kaiyu\Desktop\Rat Kinematics and dynamics\',vers,'_Total_Optimizer\OutputFigures-',RatDate,'-2020-Rat\',RatJoint,Text_append,'.mat'];

%What joints and legs to train
jnts = (3);
legz = (1);

%The project file and data of which your Animatlab Animal is being trained.
proj_file = 'C:\Users\Kaiyu\Desktop\Rat Kinematics and dynamics\Biarticular.aproj';

Data_file = 'C:\Users\kaiyu\Desktop\Rat Kinematics and dynamics\JointKinematics.mat';

if makerat == 1
    
    organism_name = 'rat';
    joints = {[],'LH_HipZ','LH_Knee','LH_AnkleZ'}';
    bodies = {'Pelvis','LH_Femur','LH_Tibia','LH_Foot'}';
    joint_limits = [-60,50;-80,60;-80,25]*pi/180;
    body_weight = 0.30 * 9.81; %N
    theta_offset = zeros(3,1);

    %Joint angle indices to use for inverse kinematics. 0 should be used to pad
    %the matrix to be of proper form.
    primary_angles = [];

    still_standing_pose = [];
    body_twists = [];
    mirrored = [];

    tune = tic;

    min_step_period = 1/300; 
    
    org = Kinematic_organism(proj_file,organism_name,bodies,joints,joint_limits,min_step_period,body_weight,theta_offset,[],mirrored,0);
    org.load_and_process_kinematic_and_dynamic_data(Data_file,0.28,0)
    
elseif makerat == 2
    loadedrat = load(RatLoc);
    org=loadedrat.obj;

elseif makerat == 3
    org.proj_file = proj_file;
end
