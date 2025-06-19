% Local_Stifness.m
function[K_e] = Local_Stifness(E, I, Le, G, A, Method)
    GA = G*A; 
    EI = E*I; 
        
    % Euler Beam Method Stiffness Matrix
    if strcmp(Method, "Euler")
        K_e = (E*I/Le)*[ ...
             12/Le^2,   6/Le,   -12/Le^2,   6/Le; 
              6/Le,       4,    -6/Le,      2; 
            -12/Le^2,  -6/Le,    12/Le^2,  -6/Le; 
              6/Le,       2,    -6/Le,      4 ...
        ];

    % Shear (Timoshenko) Full Integration Method Stiffness Matrix
    elseif strcmp(Method, "Shear Full") 
        % Matrix from PDF page 6 (Shear element, full integration)
        K_e = [ ...
            GA/Le,           -GA/2,        -GA/Le,         -GA/2;
           -GA/2,     (GA*Le)/3 + EI/Le,    GA/2,     (GA*Le)/6 - EI/Le;
           -GA/Le,            GA/2,         GA/Le,          GA/2;
           -GA/2,     (GA*Le)/6 - EI/Le,    GA/2,     EI/Le + (GA*Le)/3  ...
        ];

    % Shear (Timoshenko) Reduced Integration Method Stiffness Matrix
    elseif strcmp(Method, "Shear Reduced") 
        K_e = [ ...
            GA/Le,          -GA/2,         -GA/Le,         -GA/2; 
           -GA/2,    (GA*Le)/4 + EI/Le,     GA/2,    (GA*Le)/4 - EI/Le; 
           -GA/Le,           GA/2,          GA/Le,          GA/2; 
           -GA/2,    (GA*Le)/4 - EI/Le,     GA/2,    (GA*Le)/4 + EI/Le  ...
        ];
    else
        error('Unknown method specified for local stiffness matrix: %s', Method);
    end
end
