%Additional Code for Cavity Project for AOE 5434
%By William Brown
%Last updated: 11/23/22

%% Determine coefficients for left sweep artificial viscosity pressure term
clear, clc
A = [1, 1, 1, 1, 1;-1, 0, 1, 2, 3;1, 0, 1, 4, 9;-1, 0, 1, 8, 27;1, 0, 1, 16, 81];
b = [0;0;0;0;24];
Aright = [1, 1, 1, 1, 1;-3, -2, -1, 0, 1;9, 4, 1, 0, 1;-27, -8, -1, 0, 1;81, 16, 1, 0, 1];
inv(Aright)*b;

%% Testing Logic for 4th order pressure term
clear, clc
imax = 33;
jmax = 33;
dx = 1;
dy = 1;
u = zeros(imax,jmax,3);
four = 4.0;
six = 6.0;

for i = 2:imax-1
    for j = 2:jmax-1
        %If statement for special FDE conditions on the fourth derivative
        %near boundaries
        if i == 2 && j == 2 %Bottom Left Corner
            d4pdx4(i,j) = (u(i-1,j,1)-four*u(i,j,1)+six*u(i+1,j,1)-four*u(i+2,j,1)+u(i+3,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-1,1)-four*u(i,j,1)+six*u(i,j+1,1)-four*u(i,j+2,1)+u(i,j+3,1))/dy;
        elseif i == 2 && j == jmax-1 %Top Left Corner
            d4pdx4(i,j) = (u(i-1,j,1)-four*u(i,j,1)+six*u(i+1,j,1)-four*u(i+2,j,1)+u(i+3,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-3,1)-four*u(i,j-2,1)+six*u(i,j-1,1)-four*u(i,j,1)+u(i,j+1,1))/dy;
        elseif i == imax-1 && j == 2 %Bottom Right Corner
            d4pdx4(i,j) = (u(i-3,j,1)-four*u(i-2,j,1)+six*u(i-1,j,1)-four*u(i,j,1)+u(i+1,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-1,1)-four*u(i,j,1)+six*u(i,j+1,1)-four*u(i,j+2,1)+u(i,j+3,1))/dy;
        elseif i == imax-1 && j == jmax-1 %Top Right Corner
            d4pdx4(i,j) = (u(i-3,j,1)-four*u(i-2,j,1)+six*u(i-1,j,1)-four*u(i,j,1)+u(i+1,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-3,1)-four*u(i,j-2,1)+six*u(i,j-1,1)-four*u(i,j,1)+u(i,j+1,1))/dy;
        elseif i == 2 %at left boundary (special condition for d4pdx4 not d4pdy4)
            d4pdx4(i,j) = (u(i-1,j,1)-four*u(i,j,1)+six*u(i+1,j,1)-four*u(i+2,j,1)+u(i+3,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-2,1)-four*u(i,j-1,1)+six*u(i,j,1)-four*u(i,j+1,1)+u(i,j+2,1))/dy;
        elseif i == imax-1 %at right boundary (special condition for d4pdx4 not d4pdy4)
            d4pdx4(i,j) = (u(i-3,j,1)-four*u(i-2,j,1)+six*u(i-1,j,1)-four*u(i,j,1)+u(i+1,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-2,1)-four*u(i,j-1,1)+six*u(i,j,1)-four*u(i,j+1,1)+u(i,j+2,1))/dy;
        elseif j == 2 %at bottom boundary (special condition for d4pdy4, not d4pdx4)
            d4pdx4(i,j) = (u(i-2,j,1)-four*u(i-1,j,1)+six*u(i,j,1)-four*u(i+1,j,1)+u(i+2,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-1,1)-four*u(i,j,1)+six*u(i,j+1,1)-four*u(i,j+2,1)+u(i,j+3,1))/dy;
        elseif j == jmax-1 %at top boundary (speical condition for d4pdy4 not d4pdx4)
            d4pdx4(i,j) = (u(i-2,j,1)-four*u(i-1,j,1)+six*u(i,j,1)-four*u(i+1,j,1)+u(i+2,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-3,1)-four*u(i,j-2,1)+six*u(i,j-1,1)-four*u(i,j,1)+u(i,j+1,1))/dy;
        else %Not at any boundary (no special condition)
            d4pdx4(i,j) = (u(i-2,j,1)-four*u(i-1,j,1)+six*u(i,j,1)-four*u(i+1,j,1)+u(i+2,j,1))/dx;
            d4pdy4(i,j) = (u(i,j-2,1)-four*u(i,j-1,1)+six*u(i,j,1)-four*u(i,j+1,1)+u(i,j+2,1))/dy;            
        end
        beta2(i,j) = max(u(i,j,2)*u(i,j,2)+u(i,j,3)*u(i,j,3),rkappa*ulid*ulid);     %calculate beta^2        
        lambda_x(i,j) = half*(abs(u(i,j,2))+sqrt(u(i,j,2)*u(i,j,2)+four*beta2));    %find x eigenvalue
        lambda_y(i,j) = half*(abs(u(i,j,3))+sqrt(u(i,j,3)*u(i,j,3)+four*beta2));    %find y eigenvalue

        artviscx(i,j) = -((abs(lambda_x(i,j))*Cx)/beta2(i,j))*d4pdx4(i,j);
        artviscy(i,j) = -((abs(lambda_y(i,j))*Cy)/beta2(i,j))*d4pdy4(i,j);
    end
end

%% Array subtraction
clear, clc
imax = 33;
jmax = 33;
u = zeros(imax,jmax,3);
ummsArray = zeros(imax,jmax,3);

u = u+1;
ummsArray(:,:,1) = ummsArray(:,:,1) + 4;

rLinf(1) = max(max(abs(u(:,:,1)-ummsArray(:,:,1))));
rL2norm(1) = sqrt(sum(abs(u(:,:,1)-ummsArray(:,:,1)).^2,'all')/(imax*jmax));

%5 Array division
clear, clc

res = [1, 1, 1];
resinit = [2,2,2];

res./resinit;