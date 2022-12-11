%AOE 5434 .dat Compiler Code
%By William Brown
%Last updated: 12/8/22

clear, clc

%% Load DE information and calculate norms
grids = ["33x33","65x65","129x129","257x257"];
schemes = ["Point_Jacobi","SGS"];
filenames = ["DE-p.csv","DE-u.csv","DE-v.csv"];

L1norm = zeros(length(grids),length(schemes),length(filenames));
L2norm = zeros(length(grids),length(schemes),length(filenames));
Linfnorm = zeros(length(grids),length(schemes),length(filenames));

for a = 1:length(grids)
    cd(grids(a))
    cd("imms=1")
    for b = 1:length(schemes)
        cd(schemes(b))
            for c = 1:length(filenames)
                A = readmatrix(filenames(c)); %Temp variable to hold onto data
                L1norm(a,b,c) = sqrt((sum(abs(A(:,:)),'all'))./(size(A,1)*size(A,2)));
                L2norm(a,b,c) = sqrt((sum(abs(A(:,:)).^2,'all'))./(size(A,1)*size(A,2)));
                Linfnorm(a,b,c) = max(max(abs(A(:,:))));
            end
        cd ..
    end
    cd ..
    cd ..
end

%% Plot DE Norms versus grid resolution
h = [257/33, 257/65, 257/129, 257/257];
figure(1) %Point Jacobi
%Pressure
loglog(h,L1norm(:,1,1),'^-r','Linewidth',2) %L1 norm
hold on
loglog(h,L2norm(:,1,1),'^--r','Linewidth',2) %L2 norm
hold on
loglog(h,Linfnorm(:,1,1),'^-.r','Linewidth',2) %Linf norm
hold on
%U-velocity
loglog(h,L1norm(:,1,2),'sq-.b','Linewidth',2) %L1 norm
hold on
loglog(h,L2norm(:,1,2),'sq-.b','Linewidth',2) %L2 norm
hold on
loglog(h,Linfnorm(:,1,2),'sq-.b','Linewidth',2) %Linf norm
%V-velocity
loglog(h,L1norm(:,1,3),'o-.g','Linewidth',2) %L1 norm
hold on
loglog(h,L2norm(:,1,3),'o-.g','Linewidth',2) %L2 norm
hold on
loglog(h,Linfnorm(:,1,3),'o-.g','Linewidth',2) %Linf norm
hold off
title("DE norms vs Mesh Resoluion for Point Jacobi")
xlabel("h",fontsize=20)
ylabel("DE Norms",fontsize=20)
grid on
legend("p: L1","p: L2","p: Linf","u: L1","u: L2","u: Linf","v: L1","v: L2","v: Linf")
set(gcf,'position',[50 50 800 600])
figure(2) %Symmetric Gauss-Seidel
%Pressure
loglog(h,L1norm(:,2,1),'^-r','Linewidth',2) %L1 norm
hold on
loglog(h,L2norm(:,2,1),'^--r','Linewidth',2) %L2 norm
hold on
loglog(h,Linfnorm(:,2,1),'^-.r','Linewidth',2) %Linf norm
hold on
%U-velocity
loglog(h,L1norm(:,2,2),'sq-.b','Linewidth',2) %L1 norm
hold on
loglog(h,L2norm(:,2,2),'sq-.b','Linewidth',2) %L2 norm
hold on
loglog(h,Linfnorm(:,2,2),'sq-.b','Linewidth',2) %Linf norm
%V-velocity
loglog(h,L1norm(:,2,3),'o-.g','Linewidth',2) %L1 norm
hold on
loglog(h,L2norm(:,2,3),'o-.g','Linewidth',2) %L2 norm
hold on
loglog(h,Linfnorm(:,2,3),'o-.g','Linewidth',2) %Linf norm
hold off
title("DE norms vs Mesh Resoluion for Symmetric Gauss-Seidel")
xlabel("h",fontsize=20)
ylabel("DE Norms",fontsize=20)
grid on
legend("p: L1","p: L2","p: Linf","u: L1","u: L2","u: Linf","v: L1","v: L2","v: Linf")
set(gcf,'position',[50 50 800 600])
%% Calculate Observed Order of Accuracy
for a = 1:length(grids)-1
    for b = 1:length(schemes)
        for c = 1:length(filenames)
            pL1(a,b,c) = log(L1norm(a+1,b,c)/L1norm(a,b,c))/log(h(a+1)/h(a));
            pL2(a,b,c) = log(L2norm(a+1,b,c)/L2norm(a,b,c))/log(h(a+1)/h(a));
            pLinf(a,b,c) = log(Linfnorm(a+1,b,c)/Linfnorm(a,b,c))/log(h(a+1)/h(a));
        end
    end
end
%% Plot Observed Order of Accuracy
figure(3) %Point Jacobi Order of Accuracy
%Pressure
plot(h(2:4),pL1(:,1,1),'^-r','Linewidth',2)
hold on
plot(h(2:4),pL2(:,1,1),'^--r','Linewidth',2)
hold on
plot(h(2:4),pLinf(:,1,1),'^-.r','Linewidth',2)
hold on
%U-velocity
plot(h(2:4),pL1(:,1,2),'sq-b','Linewidth',2) %L1 norm
hold on
plot(h(2:4),pL2(:,1,2),'sq--b','Linewidth',2) %L2 norm
hold on
plot(h(2:4),pLinf(:,1,2),'sq-.b','Linewidth',2) %Linf norm
%V-velocity
plot(h(2:4),pL1(:,1,3),'o-g','Linewidth',2) %L1 norm
hold on
plot(h(2:4),pL2(:,1,3),'o--g','Linewidth',2) %L2 norm
hold on
plot(h(2:4),pLinf(:,1,3),'o-.g','Linewidth',2) %Linf norm
hold off
title("Observed Order of Accuracy vs Mesh Refinement for Point Jacobi")
xlabel("h",fontsize=20)
ylabel("$$\hat{p}$$",'Interpreter','Latex',fontsize=20)
grid on
legend("p: L1","p: L2","p: Linf","u: L1","u: L2","u: Linf","v: L1","v: L2","v: Linf")
set(gcf,'position',[50 50 800 600])
figure(4) %Symmetric Gauss-Seidel Order of Accuracy
%Pressure
plot(h(2:4),pL1(:,2,1),'^-r','Linewidth',2)
hold on
plot(h(2:4),pL2(:,2,1),'^--r','Linewidth',2)
hold on
plot(h(2:4),pLinf(:,2,1),'^-.r','Linewidth',2)
hold on
%U-velocity
plot(h(2:4),pL1(:,2,2),'sq-b','Linewidth',2) %L1 norm
hold on
plot(h(2:4),pL2(:,2,2),'sq--b','Linewidth',2) %L2 norm
hold on
plot(h(2:4),pLinf(:,2,2),'sq-.b','Linewidth',2) %Linf norm
%V-velocity
plot(h(2:4),pL1(:,2,3),'o-g','Linewidth',2) %L1 norm
hold on
plot(h(2:4),pL2(:,2,3),'o--g','Linewidth',2) %L2 norm
hold on
plot(h(2:4),pLinf(:,2,3),'o-.g','Linewidth',2) %Linf norm
hold off
title("Observed Order of Accuracy vs Mesh Refinement for Symmetric Gauss-Seidel")
xlabel("h",fontsize=20)
ylabel("$$\hat{p}$$",'Interpreter','Latex',fontsize=20)
grid on
legend("p: L1","p: L2","p: Linf","u: L1","u: L2","u: Linf","v: L1","v: L2","v: Linf")
set(gcf,'position',[50 50 800 600])