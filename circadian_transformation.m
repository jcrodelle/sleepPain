function [Cp] = circadian_transformation(C)
%this function transforms the sleep circadian process to the pain circadian
%process.

% new amplitude and vertical shift
new_alpha = 0.3411;
new_theta = 0.03;

% 12-hour shift by reflecting points
C_at_0 = C ;
C_refl = -1*C_at_0;
C_up = C_refl ;
Cp = new_alpha*C_up + new_theta;

end