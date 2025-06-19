% Global_Stifness.m
function[K_Global] = Global_Stifness(Node_DoF, K_e, Element_No)
    Node_No = Element_No + 1;
    Total_System_DoF = Node_No * Node_DoF;
% I create an empty Matrix sized as the Global Stifness Matrix
    K_Global = zeros(Total_System_DoF);

% For each i element
    for i = 1:Element_No
    
    % Add the Local K_e to the (2*i - 1) position
        j = 2*i - 1;
        K_Global( j:(j+Node_DoF*2 - 1), j:(j+Node_DoF*2 - 1) ) = ...
            K_Global( j:(j+Node_DoF*2 - 1), j:(j+Node_DoF*2 - 1) ) + K_e;

    end
