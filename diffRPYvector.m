function diffVector2 = diffRPYvector(r1,p1,h1,r2,p2,h2)
% diffVector = diffRPYvector(r1,p1,h1,r2,p2,h2)
%  Computes the transformation matrix from first roll, pitch, yaw
%  coordinates to the second set.
% Compuation is done in vector format to accomidate vectorized compuation.

diffVector = [cos(h1).*cos(h2).*cos(p1).*cos(p2) + cos(p1).*cos(p2).*sin(h1).*sin(h2) + ...
   sin(p1).*sin(p2); ...
     (-cos(p1)).*sin(p2).*sin(r1) + ...
   cos(h2).*cos(p2).*((-cos(r1)).*sin(h1) + ...
            cos(h1).*sin(p1).*sin(r1)) + ...
   cos(p2).*sin(h2).*(cos(h1).*cos(r1) + ...
            sin(h1).*sin(p1).*sin(r1)); (-cos(p1)).*cos(r1).*sin(p2) + ...
       cos(p2).*sin(h1).*(cos(r1).*sin(h2).*sin(p1) + cos(h2).*sin(r1)) + ...
       cos(h1).*cos(p2).*(cos(h2).*cos(r1).*sin(p1) - sin(h2).*sin(r1)); ...
   (-cos(p2)).*sin(p1).*sin(r2) + cos(h1).*cos(p1).*((-cos(r2)).*sin(h2) + ...
            cos(h2).*sin(p2).*sin(r2)) + ...
   cos(p1).*sin(h1).*(cos(h2).*cos(r2) + ...
            sin(h2).*sin(p2).*sin(r2)); ...
  cos(p1).*cos(p2).*sin(r1).*sin(r2) + ...
       (cos(r1).*sin(h1) - ...
      cos(h1).*sin(p1).*sin(r1)).*(cos(r2).*sin(h2) - ...
            cos(h2).*sin(p2).*sin(r2)) + (cos(h1).*cos(r1) + ...
      sin(h1).*sin(p1).*sin(r1)).* ...
         (cos(h2).*cos(r2) + sin(h2).*sin(p2).*sin(r2)); ...
  cos(p1).*cos(p2).*cos(r1).*sin(r2) + ...
       (cos(h1).*cos(r1).*sin(p1) + ...
      sin(h1).*sin(r1)).*((-cos(r2)).*sin(h2) + ...
            cos(h2).*sin(p2).*sin(r2)) + (cos(r1).*sin(h1).*sin(p1) - ...
      cos(h1).*sin(r1)).* ...
         (cos(h2).*cos(r2) + sin(h2).*sin(p2).*sin(r2)); ...
   (-cos(p2)).*cos(r2).*sin(p1) + ...
       cos(p1).*(sin(h1).*(cos(r2).*sin(h2).*sin(p2) - cos(h2).*sin(r2)) + ...
            cos(h1).*(cos(h2).*cos(r2).*sin(p2) + sin(h2).*sin(r2))); ...
     cos(p1).*cos(p2).*cos(r2).* ...
    sin(r1) + (cos(h1).*cos(r1) + sin(h1).*sin(p1).*sin(r1)).* ...
         (cos(r2).*sin(h2).*sin(p2) - cos(h2).*sin(r2)) + ...
       ((-cos(r1)).*sin(h1) + ...
      cos(h1).*sin(p1).*sin(r1)).*(cos(h2).*cos(r2).*sin(p2) + ...
            sin(h2).*sin(r2)); cos(p1).*cos(p2).*cos(r1).*cos(r2) + ...
       (cos(r1).*sin(h1).*sin(p1) - ...
      cos(h1).*sin(r1)).*(cos(r2).*sin(h2).*sin(p2) - ...
            cos(h2).*sin(r2)) + (cos(h1).*cos(r1).*sin(p1) + ...
      sin(h1).*sin(r1)).* ...
         (cos(h2).*cos(r2).*sin(p2) + sin(h2).*sin(r2))];
     
     
diffVector2 = diffVector;
diffVector2(2,:) = diffVector(4,:);
diffVector2(3,:) = diffVector(7,:);
diffVector2(4,:) = diffVector(2,:);
diffVector2(6,:) = diffVector(8,:);
diffVector2(7,:) = diffVector(3,:);
diffVector2(8,:) = diffVector(6,:);
