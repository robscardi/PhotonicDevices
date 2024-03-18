
%% Compute the optimal gaussian spot for the given waveguide 
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

E0 = 1;
%return the value of a gaussian with waist w0 dislocated by posx, posy at (x,y)
flat_gauss = @(x,y,w0,posx,posy) E0 *exp(-((x-posx)^2 + (y-posy)^2)/(w0^2));

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
max_x_ind_rib = find(x_rib == max_x_rib);
max_y_ind_rib = find(y_rib == max_y_rib);

i = find(mat_data_rid == max(mat_data_rid, [], "all"));
max_x_rid = X_rid(i);
max_y_rid = Y_rid(i);
max_x_ind_rid = find(x_rid == max_x_rid);
max_y_ind_rid = find(y_rid == max_y_rid);

i = find(mat_data_bur == max(mat_data_bur, [], "all"));
max_x_bur = X_bur(i);
max_y_bur = Y_bur(i);
max_x_ind_bur = find(x_rib == max_x_bur);
max_y_ind_bur = find(y_rib == max_y_bur);

%plot the mode
% figure(1);
% s = surf(Y_rib,X_rib,mat_data_rib./max(mat_data_rib, [], "all"));
% s.EdgeColor ="flat";
% title("Si3N4/Si02 Rib Waveguide")
% xlabel("y[\mum]", Interpreter="tex");
% ylabel("x[\mum]", Interpreter="tex");
% ax = gca;
% ax.FontSize = 13;
% 
% figure(2);
% s = surf(Y_rid,X_rid,mat_data_rid./max(mat_data_rid, [], "all"));
% s.EdgeColor ="flat";
% title("InGaAsP/InP Ridge Waveguide")
% xlabel("y[\mum]", Interpreter="tex");
% ylabel("x[\mum]", Interpreter="tex");
% ax = gca;
% ax.FontSize = 13;
% 
% figure(3);
% s = surf(Y_bur,X_bur,mat_data_bur./max(mat_data_bur, [], "all"));
% s.EdgeColor ="flat";
% title("InGaAsP/InP Buried Waveguide")
% xlabel("y[\mum]", Interpreter="tex");
% ylabel("x[\mum]", Interpreter="tex");
% ax = gca;
% ax.FontSize = 13;

%range of w0
w0 = 0.1:0.001:5;
%overlap coefficient vector
mu_rib = zeros(1,numel(w0));
mu_rid = zeros(1,numel(w0));
mu_bur = zeros(1,numel(w0));

gauss_mat_rib = zeros(yd_rib,xd_rib);
gauss_mat_rid = zeros(yd_rid,xd_rid);
gauss_mat_bur = zeros(yd_bur,xd_bur);

%iterate over w0 and calculate the coefficient vector
for j = 1:1:numel(w0)
    w = w0(j)*um;

    for i = 1:1:numel(X_rib)
        gauss_mat_rib(i) = flat_gauss(X_rib(i)*um, Y_rib(i)*um, w, max_x_rib*um, max_y_rib*um);
    end
    for i = 1:1:numel(X_rid)
        gauss_mat_rid(i) = flat_gauss(X_rid(i)*um, Y_rid(i)*um, w, max_x_rid*um, max_y_rid*um);
    end
    for i = 1:1:numel(X_bur)
        gauss_mat_bur(i) = flat_gauss(X_bur(i)*um, Y_bur(i)*um, w, max_x_bur*um, max_y_bur*um);
    end

    gauss_int_rib = trapz(delta, trapz(delta,abs(gauss_mat_rib).^2,2));
    gauss_int_rid = trapz(delta*2, trapz(delta*2,abs(gauss_mat_rid).^2,2));
    gauss_int_bur = trapz(delta, trapz(delta,abs(gauss_mat_bur).^2,2));

    mu_rib(j) = (abs(trapz(delta, trapz(delta,conj(gauss_mat_rib).*mat_data_rib,2))).^2) / (gauss_int_rib*mode_int_rib);
    mu_rid(j) = (abs(trapz(delta*2, trapz(delta*2,gauss_mat_rid.*mat_data_rid,2))).^2) / (gauss_int_rid*mode_int_rid);
    mu_bur(j) = (abs(trapz(delta, trapz(delta,gauss_mat_bur.*mat_data_bur,2))).^2) / (gauss_int_bur*mode_int_bur);

end

%find maximum mu index
opt_mu_ind_rib = find(mu_rib == max(mu_rib, [], "all"));
opt_mu_ind_rid = find(mu_rid == max(mu_rid, [], "all"));
opt_mu_ind_bur = find(mu_bur == max(mu_bur, [], "all"));
%find optimal waist
opt_w_rib = w0(opt_mu_ind_rib);
opt_w_rid = w0(opt_mu_ind_rid);
opt_w_bur = w0(opt_mu_ind_bur);

%return the vector of the gaussian slice parallel to the y axis in 
%correspondece of the peak of the mode
gauss_sil_rib = zeros(1,yd_rib);
gauss_sil_rid = zeros(1,yd_rid);
gauss_sil_bur = zeros(1,yd_bur);

for i = 1:1:yd_rib
    gauss_sil_rib(i) = flat_gauss(max_x_rib*um, y_rib(i)*um, opt_w_rib*um, max_x_rib*um, max_y_rib*um);
end

for i = 1:1:yd_rid
    gauss_sil_rid(i) = flat_gauss(max_x_rid*um, y_rid(i)*um, opt_w_rid*um, max_x_rid*um, max_y_rid*um);
end

for i = 1:1:yd_bur
    gauss_sil_bur(i) = flat_gauss(max_x_bur*um, y_bur(i)*um, opt_w_bur*um, max_x_bur*um, max_y_bur*um);
end

%mode slice
data_sil_rib = mat_data_rib(:,max_x_ind_rib)/max(mat_data_rib, [],"all");
data_sil_rid = mat_data_rid(:,max_x_ind_rid)/max(mat_data_rid, [],"all");
data_sil_bur = mat_data_bur(:,max_x_ind_bur)/max(mat_data_bur, [],"all");

%maximum difference
max_perc_difference_rib = max(abs(gauss_sil_rib-data_sil_rib'),[], "all")
max_perc_difference_bur = max(abs(gauss_sil_bur-data_sil_bur'),[], "all")
max_perc_difference_rid = max(abs(gauss_sil_rid-data_sil_rid'),[], "all")

%plot gauss and mode 
figure(4)
hold on
plot(y_rib, gauss_sil_rib)
plot(y_rib, data_sil_rib)
title("Rib WG Mode VS Gaussian")
xlabel("y[\mum]", Interpreter="tex");
ylabel("E/$E_{max}$", Interpreter="tex");
ax = gca;
ax.FontSize = 13;

figure(5)
hold on
plot(y_rid, gauss_sil_rid)
plot(y_rid, data_sil_rid)
title("Ridge WG Mode VS Gaussian")
xlabel("y[\mum]", Interpreter="tex");
ylabel("E/$E_{max}$", Interpreter="tex");
ax = gca;
ax.FontSize = 13;


figure(6)
hold on
plot(y_bur, gauss_sil_bur)
plot(y_bur, data_sil_bur)
title("Buried WG Mode VS Gaussian")
xlabel("y[\mum]", Interpreter="tex");
ylabel("E/$E_{max}$", Interpreter="tex");
ax = gca;
ax.FontSize = 13;

%mu over w0*2
figure(7)
p = plot(w0*2,mu_rib);
title("\mu rib waveguide over gaussian spot", Interpreter="tex")
xlabel("y[\mum]", Interpreter="tex");
ylabel("\mu", Interpreter="tex");
xlim([0.2 10])
datatip(p, opt_w_rib*2, max(mu_rib, [], "all"));
ax = gca;
ax.FontSize = 13;

%mu over w0*2
figure(8)
p = plot(w0*2,mu_rid);
title("\mu ridge waveguide over gaussian spot", Interpreter="tex")
xlabel("y[\mum]", Interpreter="tex");
ylabel("\mu", Interpreter="tex");
xlim([0.2 10])
datatip(p, opt_w_rid*2, max(mu_rid, [], "all"));
ax = gca;
ax.FontSize = 13;

%mu over w0*2
figure(9)
p = plot(w0*2,mu_bur);
title("\mu buried waveguide over gaussian spot", Interpreter="tex")
xlabel("y[\mum]", Interpreter="tex");
ylabel("\mu", Interpreter="tex");
xlim([0.2 10])
datatip(p, opt_w_bur*2, max(mu_bur, [], "all"));
ax = gca;
ax.FontSize = 13;

