function Y = RPYmat(Roll,Pitch,Yaw)
% Y = RPYmat(Roll,Pitch,Yaw)
% Computes Roll, Pitch Yah Transformation matrix using Lenschow's
% coordinate frame where:
% x is directed toward the aircraft nose
% y is directed toward starboard (right) wing
% z is directed to the aircraft bottom

Y = [cos(Yaw) -sin(Yaw) 0; sin(Yaw) cos(Yaw) 0; 0 0 1]*...
    [cos(Pitch) 0 sin(Pitch); 0 1 0; -sin(Pitch) 0 cos(Pitch)]*...
    [1 0 0; 0 cos(Roll) -sin(Roll); 0 sin(Roll) cos(Roll)];