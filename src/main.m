%% Finite Element Analysis for Cantilever Beam Bending
% Euler Beam
%--------------------------------------------------------------------------
close all
clc
clear
%--------------------------------------------------------------------------
%%  1. Problem Parameters (SI)

% Beam Geometry (Solid Semicircle Cross Section)
shape_id = 8;
L = input('What is the length of the Beam in meters? ');                       % Length
a = input([ ...
    'What is the vertex angle of the isosceles ', ...                          % Angle - (See the Cross Section image)
    'triangle for the solid semicircle in degrees? ' ...
    ]);                                                             

% Material Properties
E = 71e9;                       % Young Modulus
v = 0.33;                       % Poisson Ratio
G = E/(2*(1+v));                % Shear Modulus

% Load
C = 1.6;                        % Load Coefficient
Pz = -input(['What is the Magnitude of the applied', ...
    ' force on the Tip of the Beam (positive is down, negative is up) ']);     % Load 

%% Area and Area Moment of Inertia for Solid Semicircle

h = input(['What is the length of each of ',...
    'the isosceles triangle side in meters? ' ]);         

[A_beam_nominal, I_beam_nominal] = Area_and_Moment_of_Inertia(h,a);

fprintf('--- Nominal Beam Parameters ---\n');
fprintf('Nominal Height parameter h_nominal: %.4f m\n', h);
fprintf('Nominal Moment of Inertia Izz: %.3e m^4\n', I_beam_nominal);
fprintf('Load Pz: %.2f N\n', Pz);

% Theoretical Euler tip displacement for reference
Tip_disp_Euler_theoretical = Pz * L^3 / (3 * E * I_beam_nominal);
fprintf('Theoretical Euler tip displacement: %.6e m\n\n', Tip_disp_Euler_theoretical);
%--------------------------------------------------------------------------
%%  2. Discritization, Mesh and Stifness Matrix

Element_No = input('How many Finite Elements do you want to use? ');
Node_DoF = 2;
Le = L/Element_No;
Node_No = Element_No + 1;
Total_DoFs = Node_No * Node_DoF;

%% Define Load Vector
Active_DoF_Count = Total_DoFs - Node_DoF;

F = zeros(Active_DoF_Count, 1);
F(Active_DoF_Count - 1) = Pz;                                                  % Apply load to transverse DoF of the last node

%% Methods Used and Solution for user specified values
method_vector = ["Euler", "Shear Reduced", "Shear Full"];

T = cell(1, length(method_vector));                                            % To store full displacement vectors
R = cell(1, length(method_vector)); % To store full rotation vectors


for i = 1:length(method_vector)
    fprintf('Processing method: %s\n', method_vector(i));
    K_e = Local_Stifness(E, I_beam_nominal, Le, G, A_beam_nominal, method_vector(i));
    K_Global_System = Global_Stifness(Node_DoF, K_e, Element_No);

    %% 3. Boundary Conditions & Solve
    K_e_Global_Active = K_Global_System;
    K_e_Global_Active(1:Node_DoF,:) = [];
    K_e_Global_Active(:,1:Node_DoF) = [];
    
    % Solve for active (free) DoFs
    Global_Displacements_Active = K_e_Global_Active \ F;
    
    % Check if the solver failed and returned NaN or Inf
    if any(isnan(Global_Displacements_Active)) || any(isinf(Global_Displacements_Active))
        fprintf('WARNING: Solver produced NaN/Inf for method %s. Results will not be plotted.\n', method_vector(i));
        num_active_nodes = Node_No - 1;
        active_displacements = nan(num_active_nodes, 1);
        active_rotations = nan(num_active_nodes, 1);
    else
        active_displacements = Global_Displacements_Active(1:Node_DoF:end);
        active_rotations = Global_Displacements_Active(2:Node_DoF:end);
    end

    % Build the complete result vectors, including the fixed node (w=0, theta=0)
    T{i} = [0; active_displacements]; % Now the vector has span = Node_No
    R{i} = [0; active_rotations];     % Now the vector has span = Node_No
    fprintf('---------------------------------\n');
end

%% Plotting

% Define x-coordinates for ALL nodes
x_values_for_nodes = linspace(0, L, Node_No);

% Tip transverse displacement for each method
fprintf('Euler tip displacement: %.6e m\n\n', T{1}(end));
fprintf('Shear Reduced tip displacement: %.6e m\n\n', T{2}(end));
fprintf('Shear Full tip displacement: %.6e m\n\n', T{3}(end));

% Plotting Transverse Displacements with Dual Y-Axes
figure;
color_order = get(gca, 'colororder'); 
yyaxis left;
h1 = plot(x_values_for_nodes, T{2}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(2)), 'Color', color_order(2,:));
hold on;
h2 = plot(x_values_for_nodes, T{3}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(3)), 'Color', color_order(3,:));
ylabel('Displacement w (m) - Timoshenko');
xlabel('x (m)');
ax = gca; ax.YColor = 'k';

yyaxis right;
h3 = plot(x_values_for_nodes, T{1}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(1)), 'Color', color_order(1,:));
ylabel('Displacement w (m) - Euler');
ax = gca; ax.YColor = color_order(1,:);

hold off; title('Cantilever Beam Deflection (Dual Axis)'); grid on;
legend([h1, h2, h3], 'Location', 'best');

% Plotting Rotations with Dual Y-Axes
figure;
yyaxis left;
h1_rot = plot(x_values_for_nodes, R{2}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(2)), 'Color', color_order(2,:));
hold on;
h2_rot = plot(x_values_for_nodes, R{3}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(3)), 'Color', color_order(3,:));
ylabel('Rotation \theta (rad) - Timoshenko');
xlabel('x (m)');
ax = gca; ax.YColor = 'k';

yyaxis right;
h3_rot = plot(x_values_for_nodes, R{1}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(1)), 'Color', color_order(1,:));
ylabel('Rotation \theta (rad) - Euler');
ax = gca; ax.YColor = color_order(1,:);

hold off; title('Cantilever Beam Rotation (Dual Axis)'); grid on;
legend([h1_rot, h2_rot, h3_rot], 'Location', 'best');



%% 5. Solution with L/h = 100 (h = L*100) and Nelm = [8 24 144]

Nelm = [8 24 144];
h_100 = L/100;
[A_beam_100, I_beam_100] = Area_and_Moment_of_Inertia(h_100,a);

T_100 = cell(length(Nelm), length(method_vector));
R_100 = cell(length(Nelm), length(method_vector));

for i = 1:3
    fprintf("Number of Elements used for calculation= %d\n\n", Nelm(i));
    for j = 1:length(method_vector)
        
        Le_100 = L/Nelm(i);
        Node_No_100 = Nelm(i) + 1;
        Total_DoFs_100 = Node_No_100 * Node_DoF;

        fprintf('Processing method: %s\n', method_vector(j));
        K_e_100 = Local_Stifness(E, I_beam_100, Le_100, G, A_beam_100, method_vector(j));
        K_e_Global_System_100 = Global_Stifness(Node_DoF, K_e_100, Nelm(i));
    
        %% Boundary Conditions & Solve
        K_e_Global_Active = K_e_Global_System_100;
        K_e_Global_Active(1:Node_DoF,:) = [];
        K_e_Global_Active(:,1:Node_DoF) = [];
        
        % Solve for active (free) DoFs   
        Active_DoF_Count_100 = Total_DoFs_100 - Node_DoF;
        
        F_100 = zeros(Active_DoF_Count_100, 1);
        F_100(Active_DoF_Count_100 - 1) = Pz; % Apply load to transverse DoF of the last node

        Global_Displacements_Active = K_e_Global_Active \ F_100;
        
        % *** ROBUSTNESS FIX ***
        % Check if the solver failed and returned NaN or Inf
        if any(isnan(Global_Displacements_Active)) || any(isinf(Global_Displacements_Active))
            fprintf('WARNING: Solver produced NaN/Inf for method %s. Results will not be plotted.\n', method_vector(j));
            num_active_nodes = Node_No_100 - 1;
            active_displacements = nan(num_active_nodes, 1);
            active_rotations = nan(num_active_nodes, 1);
        else
            active_displacements = Global_Displacements_Active(1:Node_DoF:end);
            active_rotations = Global_Displacements_Active(2:Node_DoF:end);
        end
    
        % Build the complete result vectors, including the fixed node (w=0, theta=0)
        T_100{i,j} = [0; active_displacements]; % Now the vector has span = Node_No
        R_100{i,j} = [0; active_rotations];     % Now the vector has span = Node_No
        fprintf('---------------------------------\n');
    end
    %% Plotting for each Element No
    
    % Define x-coordinates for ALL nodes
    x_values_for_nodes_100 = linspace(0, L, Node_No_100);
    
    % Tip transverse displacement for each method
    fprintf('Euler tip displacement for L/h = 100: %.6e m\n\n', T_100{i,1}(end));
    fprintf('Shear Reduced tip displacement for L/h = 100: %.6e m\n\n', T_100{i,2}(end));
    fprintf('Shear Full tip displacement for L/h = 100: %.6e m\n\n', T_100{i,3}(end));
    
    % Plotting Transverse Displacements with Dual Y-Axes
    figure;
    color_order = get(gca, 'colororder'); 
    yyaxis left;
    h1 = plot(x_values_for_nodes_100, T_100{i,2}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(2)), 'Color', color_order(2,:));
    hold on;
    h2 = plot(x_values_for_nodes_100, T_100{i,3}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(3)), 'Color', color_order(3,:));
    ylabel('Displacement w (m) - Timoshenko');
    xlabel('x (m)');
    ax = gca; ax.YColor = 'k';
    
    yyaxis right;
    h3 = plot(x_values_for_nodes_100, T_100{i,1}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(1)), 'Color', color_order(1,:));
    ylabel('Displacement w (m) - Euler');
    ax = gca; ax.YColor = color_order(1,:);
    
    hold off; title(sprintf('Beam Transverse Displacement with %d Elements', Nelm(i))); grid on;
    legend([h1, h2, h3], 'Location', 'best');
    
    % Plotting Rotations with Dual Y-Axes
    figure;
    yyaxis left;
    h1_rot = plot(x_values_for_nodes_100, R_100{i,2}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(2)), 'Color', color_order(2,:));
    hold on;
    h2_rot = plot(x_values_for_nodes_100, R_100{i,3}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(3)), 'Color', color_order(3,:));
    ylabel('Rotation \theta (rad) - Timoshenko');
    xlabel('x (m)');
    ax = gca; ax.YColor = 'k';
    
    yyaxis right;
    h3_rot = plot(x_values_for_nodes_100, R_100{i,1}, 'LineWidth', 1.5, 'DisplayName', char(method_vector(1)), 'Color', color_order(1,:));
    ylabel('Rotation \theta (rad) - Euler');
    ax = gca; ax.YColor = color_order(1,:);
    
    hold off; title(sprintf('Beam Rotation with %d Elements', Nelm(i))); grid on;
    %title(['Iteration ' num2str(k) ', value = ' num2str(value)]);
    legend([h1_rot, h2_rot, h3_rot], 'Location', 'best');
    

end    


%% 5. Solution for L/h = (2:150) and user specified Element No

L_h = 2:1:150;                                                                   % L/h from 2 to 150
h_loop = L ./ L_h;                                                               % thickness for each L/h

w_tip = zeros(length(method_vector), numel(L_h));

for j = 1:length(method_vector)
    for k = 1:length(h_loop)
        [A_beam_loop, I_beam_loop] = Area_and_Moment_of_Inertia(h_loop(k),a);
        K_e_loop = Local_Stifness(E, I_beam_loop, Le, G, A_beam_loop, method_vector(j));
        K_Global_System = Global_Stifness(Node_DoF, K_e_loop, Element_No);
    
        %% 3. Boundary Conditions & Solve
        K_e_Global_Active = K_Global_System;
        K_e_Global_Active(1:Node_DoF,:) = [];
        K_e_Global_Active(:,1:Node_DoF) = [];
        
        % Solve for active (free) DoFs
        U_Active = K_e_Global_Active \ F;
        
        % Check if the solver failed and returned NaN or Inf
        if any(isnan(U_Active)) || any(isinf(U_Active))
            fprintf('WARNING: Solver produced NaN/Inf for method %s. Results will not be plotted.\n', method_vector(j));
            num_active_nodes = Node_No - 1;
            w_vec = nan(num_active_nodes, 1);
            w_tip(j,k) = w_vec(end);                                            % tip deflection

        else
            w_vec = U_Active(1:Node_DoF:end);
            w_tip(j,k) = w_vec(end);                                            % tip deflection

        end
    end
end

% Normalize by Euler (method_vector(1))
w_norm = w_tip ./ w_tip(1,:); 

% Plot deflection
figure;
hold on;                                    
color_order = get(gca,'ColorOrder');

h_sr = plot(L_h, w_tip(2,:), '-', 'LineWidth',1.5, 'Color',color_order(2,:), 'DisplayName','Shear Reduced');
h_sf = plot(L_h, w_tip(3,:), '-', 'LineWidth',1.5, 'Color',color_order(3,:), 'DisplayName','Shear Full');
h_e  = plot(L_h, w_tip(1,:), '--', 'DisplayName','Euler', 'Color',color_order(1,:),'LineWidth',2);

xlabel('L/h');
ylabel('Tip deflection w (m)');
title('Tip deflection vs L/h');
grid on;
legend([h_e,h_sr,h_sf],'Location','best');  
hold off;

% Plot normalized deflection
figure;
hold on;                                    
color_order = get(gca,'ColorOrder');

h_en  = plot(L_h, w_norm(1,:), '--', 'LineWidth',1.5, 'Color',color_order(1,:), 'DisplayName','Euler w/w_e');
h_srn = plot(L_h, w_norm(2,:), '--', 'LineWidth',1.5, 'Color',color_order(2,:), 'DisplayName','Shear Reduced w/w_e');
h_sfn = plot(L_h, w_norm(3,:), '--', 'LineWidth',1.5, 'Color',color_order(3,:), 'DisplayName','Shear Full w/w_e');

xlabel('L/h');
ylabel('Normalized deflection w/w_{Euler}');
title('Normalized tip deflection vs L/h');
grid on;
legend([h_en,h_srn,h_sfn],'Location','best');  
hold off;

disp('Script finished.');
