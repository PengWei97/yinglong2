clear all
close all
clc

%% Step 1: Extract FeatureID_NodeID data
% File path
file_path = 'd.stelset';

% Read the file
file_data = fileread(file_path);

% Split the file data into lines
lines = strsplit(file_data, '\n');

% Initialize the map
FeatureID_NodeID = containers.Map('KeyType', 'double', 'ValueType', 'double');

% Process each line of data
for i = 1:length(lines)
    line = strtrim(lines{i});
    
    % Skip empty lines
    if isempty(line)
        continue;
    end
    
    % Split the line into tokens
    tokens = strsplit(line);
    
    % Get the Feature_ID and node_ID
    Feature_ID = str2double(tokens{1});
    node_IDs = str2double(tokens(2:end));
    
    % Map each node_ID to Feature_ID
    for j = 1:length(node_IDs)
        FeatureID_NodeID(node_IDs(j)) = Feature_ID;
    end
end

%% Step 2: Extract NodeID_xyz data
% Read the file content
file_path = 'd.stnode';
file_data = fileread(file_path);

% Split the data lines
data_lines = split(file_data, newline);

% Create the NodeID_xyz table
NodeID_xyz = table('Size', [length(data_lines)-1, 4], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames', {'NodeID', 'x', 'y', 'z'});

% Parse each line of data
for i = 1:length(data_lines)-1
    line_data = str2double(split(data_lines{i}, ' '));
    NodeID_xyz(i, :) = array2table(line_data', 'VariableNames', {'NodeID', 'x', 'y', 'z'});
end

%% Step 3: Extract FeatureID_EulerAngle data
% File path
file_path = 'd.msh';

% Open the file
fileID = fopen(file_path, 'r');

% Find the line with "$ElsetOrientations"
line = fgetl(fileID);
while ischar(line)
    if contains(line, '$ElsetOrientations')
        break;
    end
    line = fgetl(fileID);
end

% Skip a line
line = fgetl(fileID);

% Read the data and store it in a cell array
data = {};
line = fgetl(fileID);
while ischar(line)
    if contains(line, '$EndElsetOrientations')
        break;
    end
    data = [data; strsplit(line)];
    line = fgetl(fileID);
end

% Close the file
fclose(fileID);
FeatureID_EulerAngle = table('Size', [length(data), 4], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames', {'FeatureID', 'phi1', 'PHI', 'phi2'});

% Parse each line of data
FeatureID_EulerAngle.FeatureID = str2double(data(:, 1));
FeatureID_EulerAngle.phi1 = str2double(data(:, 2));
FeatureID_EulerAngle.PHI = str2double(data(:, 3));
FeatureID_EulerAngle.phi2 = str2double(data(:, 4));

%% Step 4: Create the INL_data table
INL_data = table('Size', [length(NodeID_xyz.x), 9], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'int32', 'int32', 'int32'}, 'VariableNames', {'phi1', 'PHI', 'phi2', 'x', 'y', 'z', 'FeatureID', 'PhaseId', 'Symmetry'});

for i = 1:length(NodeID_xyz.x)
    index_Gr = FeatureID_NodeID(i);

    INL_data.phi1(i) = FeatureID_EulerAngle.phi1(index_Gr) * (pi/180);
    INL_data.PHI(i) = FeatureID_EulerAngle.PHI(index_Gr) * (pi/180);
    INL_data.phi2(i) = FeatureID_EulerAngle.phi2(index_Gr) * (pi/180);

    INL_data.x(i) = NodeID_xyz.x(i);
    INL_data.y(i) = NodeID_xyz.y(i);
    INL_data.z(i) = NodeID_xyz.z(i);    

    INL_data.FeatureID(i) = index_Gr;      
end

INL_data.PhaseId = ones(length(NodeID_xyz.x),1) * 0;
INL_data.Symmetry = ones(length(NodeID_xyz.x),1) * 62;

%% Write the INL_data table to a file (demo_msh2inl.txt)
file_path = 'case_msh2inl.inl';
fileID = fopen(file_path, 'w');

% Write the header information
fprintf(fileID, '# File written from OrientationAnalysis Version 6.5.163.1998a502a\n');
fprintf(fileID, ['# DateTime: ', datestr(now), '\n']);
fprintf(fileID, ['# X_STEP: %.6f\n'], max(INL_data.x)/(sqrt(length(INL_data.x))-1));
fprintf(fileID, ['# Y_STEP: %.6f\n'], max(INL_data.y)/(sqrt(length(INL_data.y))-1));
fprintf(fileID, ['# Z_STEP: %.6f\n'], max(INL_data.z)/(sqrt(length(INL_data.z))-1));
fprintf(fileID, '#\n');
fprintf(fileID, ['# X_MIN: %.6f\n'], min(INL_data.x)- 0.1e-2);
fprintf(fileID, ['# Y_MIN: %.6f\n'], min(INL_data.y)- 0.1e-2);
fprintf(fileID, ['# Z_MIN: %.6f\n'], min(INL_data.z));
fprintf(fileID, '#\n');
fprintf(fileID, ['# X_MAX: %.6f\n'], max(INL_data.x));
fprintf(fileID, ['# Y_MAX: %.6f\n'], max(INL_data.y));
fprintf(fileID, ['# Z_MAX: %.6f\n'], max(INL_data.z));
fprintf(fileID, '#\n');
fprintf(fileID, ['# X_DIM: %d\n'], sqrt(length(INL_data.x)));
fprintf(fileID, ['# Y_DIM: %d\n'], sqrt(length(INL_data.y)));
fprintf(fileID, ['# Z_DIM: 0\n']);
fprintf(fileID, '#\n');
fprintf(fileID, '# Phase_1:\n');
fprintf(fileID, ['# Symmetry_1: %d\n'], 62);
fprintf(fileID, ['# Features_1: %d\n'], max(INL_data.FeatureID));
fprintf(fileID, '#\n');
fprintf(fileID, ['# Num_Features: %d\n'], max(INL_data.FeatureID));
fprintf(fileID, '#\n');
fprintf(fileID, '# phi1 PHI phi2 x y z FeatureId PhaseId Symmetry\n');

% Write the INL_data table
for i = 1:height(INL_data)
    fprintf(fileID, '%.6g %.6g %.6g %.6g %.6g %.6g %d %d %d\n', ...
        INL_data.phi1(i), INL_data.PHI(i), INL_data.phi2(i), ...
        INL_data.x(i), INL_data.y(i), INL_data.z(i), ...
        INL_data.FeatureID(i), INL_data.PhaseId(i), INL_data.Symmetry(i));
end

fclose(fileID);

disp(['File "', file_path, '" created successfully.']);