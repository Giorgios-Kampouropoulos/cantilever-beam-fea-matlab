function[A, I_x] = Area_and_Moment_of_Inertia(h,a)
    
    r = h*sin(a/2);                 % Radius
    A = (1/2)*pi*(r^2);             % Area
    y_c = (4*r)/(3*pi);             % Centroid of the Semicircle
    
    Ix_c = pi*(r^4)/4;              % Area Moment of Inertia for Full Circle
    
    % Area Moment of Inertia for Semicircle with respect to the base (diameter)
    Ix_sc_base = Ix_c/2;             
    
    % Area Moment of Inertia for Semicircle with respect to it's centroid
    I_x= Ix_sc_base - A*(y_c^2); 
