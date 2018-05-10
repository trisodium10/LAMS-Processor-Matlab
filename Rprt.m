function rP = Rprt(theta,v)
% rP = Rprt(theta,v)
% Returns 3D PRT rotation matrix for rotation angle theta about the 3D
% vector v.  
% Rotation angle is defined looking into v.  i.e.  If v points
% in the positive x direction, positive rotation couples y into z and z
% into negative y.
v = v/norm(v);
rP = zeros(3);
for i = 1:3
    for j = 1:3
        sumLC = 0;
        for k = 1:3
            sumLC = sumLC+LeviCivita(i,j,k)*v(k);
        end
        rP(i,j) = (i == j)*cos(theta)+v(i)*v(j)*(1-cos(theta))-sin(theta)*sumLC;
    end
end