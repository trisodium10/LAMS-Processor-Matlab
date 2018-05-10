function vect = RPY_Matrix_Vector(x,vL)
% Estimation function written for EstimateAverageRPY_Error_SDN_Honewell.m
v = RPYmat(x(1),x(2),x(3))*vL;
vect = v(:);