function Rm = AvgTransformationMatrix(T)
% Rm = AvgTransformationMatrix(T)
% Computes the average transformation of T assuming T is purely rotational
% by averaging over all Quarternions.
% The matrix T is expected in vector format and the averager rotation
% matrix Rm is returned in vector format.
% To format a matrix T, use T(:).
% to obtain Rm in matrix format use reshape(Rm,3,3)



% T = Tl2aVec(:,iUseData_Tl2a);

qr = sqrt(1+T(1,:)+T(5,:)+T(9,:))/2;
qx = (T(6,:)-T(8,:))./(4*qr);
qy = (T(7,:)-T(3,:))./(4*qr);
qz = (T(2,:)-T(4,:))./(4*qr);

% theta = 2*acos(qr);
% qnorm = zeros(3,length(qr));
% qnorm(1,:) = qx./sin(0.5*theta);
% qnorm(2,:) = qy./sin(0.5*theta);
% qnorm(3,:) = qz./sin(0.5*theta);
% 
% qMeanExp = nanmean(0.5*theta.*qnorm,2);
% 
% thetaMean = 2*asin(norm(qMeanExp));

% qrM = cos(thetaMean/2);
% qxM = sin(thetaMean/2)*qMeanExp(1);
% qyM = sin(thetaMean/2)*qMeanExp(2);
% qzM = sin(thetaMean/2)*qMeanExp(3);

qrM = mean(qr);
qxM = mean(qx);
qyM = mean(qy);
qzM = mean(qz);


Rm(1,:) = 1-2*qyM.^2-2*qzM.^2;
Rm(2,:) = 2*(qxM.*qyM+qzM.*qrM);
Rm(3,:) = 2*(qxM.*qzM-qyM.*qrM);
Rm(4,:) = 2*(qxM.*qyM-qzM.*qrM);
Rm(5,:) = 1-2*qxM.^2-2*qzM.^2;
Rm(6,:) = 2*(qyM.*qzM+qxM.*qrM);
Rm(7,:) = 2*(qxM.*qzM+qyM.*qrM);
Rm(8,:) = 2*(qyM.*qzM-qxM.*qrM);
Rm(9,:) = 1-2*qxM.^2-2*qyM.^2;