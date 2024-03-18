%% Compute overlap coefficient between a flat gaussian beam and a waveguide mode

clear variables
close all
um = 1e-6;

data_rib = load("qte_s_rib_Triplex.txt"); 
% 1um*1um rib waveguide with 1.98 Si3N4 core and 1.45 SiO2 substrate

data_rid = load("qte_s_ridge_InP.txt");
% 2um * (1+1)um ridge waveguide with InGaAsP core and InP substrate

data_bur = load("qte_s_bur_InGaAsP.txt");
% 0.9 um *0.4 um buried waveguide in InGaAsP core and InP substrate
delta = 0.01*um;

x_rib = -1:0.01:1.99;
xd_rib = numel(x_rib);
y_rib = -2.5:0.01:2.5;
yd_rib = numel(y_rib);
[X_rib,Y_rib] = meshgrid(x_rib,y_rib);

x_rid = -2:0.02:3.98;
xd_rid = numel(x_rid);
y_rid = -5:0.02:5;
yd_rid = numel(y_rid);
[X_rid,Y_rid] = meshgrid(x_rid,y_rid);

x_bur = -1:0.01:0.99;
xd_bur = numel(x_bur);
y_bur = -2.25:0.01:2.25;
yd_bur = numel(y_bur);
[X_bur,Y_bur] = meshgrid(x_bur,y_bur);



%gaussian waist (is the half of the spot)
w0 = 10*um/2;
E0 = 1;

%return the value of a gaussian with waist w0 dislocated by posx, posy at (x,y)
flat_gauss = @(x,y,w,posx,posy) E0 *exp(-((x-posx)^2 + (y-posy)^2)/(w^2));

%reshape the input array into a yd by xd matrix
mat_data_rib = reshape(data_rib(:,3), yd_rib, xd_rib);
mat_data_rid = reshape(data_rid(:,3), yd_rid, xd_rib);
mat_data_bur = reshape(data_bur(:,3), yd_bur, xd_bur);

%calculate the integral of the mode squared
mode_int_rib = trapz(delta, trapz(delta,abs(mat_data_rib).^2,2));
mode_int_rid = trapz(delta*2, trapz(delta*2,abs(mat_data_rid).^2,2));
mode_int_bur = trapz(delta, trapz(delta,abs(mat_data_bur).^2,2));


%find the peak of the waveguide mode
i = find(mat_data_rib == max(mat_data_rib, [], "all"));
max_x_rib = X_rib(i);
max_y_rib = Y_rib(i);

i = find(mat_data_rid == max(mat_data_rid, [], "all"));
max_x_rid = X_rid(i);
max_y_rid = Y_rid(i);

i = find(mat_data_bur == max(mat_data_bur, [], "all"));
max_x_bur = X_bur(i);
max_y_bur = Y_bur(i);



%plot the mode
figure(1);
s = surf(Y_rib,X_rib,mat_data_rib./max(mat_data_rib, [], "all"));
s.EdgeColor ="flat";
title("Si3N4/Si02 Rib Waveguide")
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
xlim([-2.5 2.5])
ylim([-1 1.99])
ax = gca;
ax.FontSize = 13;

figure(2);
s = surf(Y_rid,X_rid,mat_data_rid./max(mat_data_rid, [], "all"));
s.EdgeColor ="flat";
title("InGaAsP/InP Ridge Waveguide")
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
xlim([-5 5])
ylim([-2 3.98])
ax = gca;
ax.FontSize = 13;

figure(3);
s = surf(Y_bur,X_bur,mat_data_bur./max(mat_data_bur, [], "all"));
s.EdgeColor ="flat";
title("InGaAsP/InP Buried Waveguide")
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
xlim([-2.25 2.25])
ylim([-1 0.99])
ax = gca;
ax.FontSize = 13;

%return the gaussian's matrix
gauss_mat_rib = zeros(yd_rib,xd_rib);
for i = 1:1:numel(X_rib)
    gauss_mat_rib(i) = flat_gauss(X_rib(i)*um, Y_rib(i)*um, w0, max_x_rib*um, max_y_rib*um);
end

gauss_mat_rid = zeros(yd_rid,xd_rid);
for i = 1:1:numel(X_rid)
    gauss_mat_rid(i) = flat_gauss(X_rid(i)*um, Y_rid(i)*um, w0, max_x_rid*um, max_y_rid*um);
end

gauss_mat_bur = zeros(yd_bur,xd_bur);
for i = 1:1:numel(X_bur)
    gauss_mat_bur(i) = flat_gauss(X_bur(i)*um, Y_bur(i)*um, w0, max_x_bur*um, max_y_bur*um);
end
k = -5:0.01:5;
[X,Y] = meshgrid(k,k);
gauss_mat_def = zeros(numel(k),numel(k));
for i = 1:1:numel(X)
    gauss_mat_def(i) = flat_gauss(X(i)*um, Y(i)*um, 5*um, 0, 0);
end


%plot the gaussian
figure(4)
s = surf(X,Y, gauss_mat_def);
s.EdgeColor ="flat";
title("10\mum-spot gaussian beam");
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
ax = gca;
ax.FontSize = 13;


%calculate the integral of the gaussian squared 
gauss_int_rib = trapz(delta, trapz(delta,abs(gauss_mat_rib).^2,2));
gauss_int_rid = trapz(delta*2, trapz(delta*2,abs(gauss_mat_rid).^2,2));
gauss_int_bur = trapz(delta, trapz(delta,abs(gauss_mat_bur).^2,2));


%calculate the overlapping integral
mu_rib = (abs(trapz(delta, trapz(delta,gauss_mat_rib.*mat_data_rib,2))).^2) / (gauss_int_rib*mode_int_rib)
mu_rid = (abs(trapz(delta*2, trapz(delta*2,gauss_mat_rid.*mat_data_rid,2))).^2) / (gauss_int_rid*mode_int_rid)
mu_bur = (abs(trapz(delta, trapz(delta,gauss_mat_bur.*mat_data_bur,2))).^2) / (gauss_int_bur*mode_int_bur)


