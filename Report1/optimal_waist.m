
%% Compute the optimal gaussian spot for the given waveguide 
clear variables
close all
um = 1e-6;
data_rib = load("qte_s_rib_Triplex.txt");
% 1um*1um rib waveguide with 1.98 Si3N4 core and 1.45 SiO2 substrate

%data_channel = load

%data_ridge = 

x = -1:0.01:1.99;
xd = numel(x);
y = -2.5:0.01:2.5;
yd = numel(y);
delta = 0.01*um;

x_um = x*um;
y_um = y*um;
[X,Y] = meshgrid(x,y);

E0 = 1;
%return the value of a gaussian with waist w0 dislocated by posx, posy at (x,y)
flat_gauss = @(x,y,w0,posx,posy) E0 *exp(-((x-posx)^2 + (y-posy)^2)/(w0^2));

%reshape the input array into a yd by xd matrix
mat_data_rib = reshape(data_rib(:,3), yd, xd);

%find the peak of the waveguide mode
i = find(mat_data_rib == max(mat_data_rib, [], "all"));
max_x = X(i);
max_y = Y(i);
max_x_ind = find(x == max_x);
max_y_ind = find(y == max_y);

%calculate the integral of the mode squared
mode_int = trapz(delta, trapz(delta,abs(mat_data_rib).^2,2));

%plot the mode
figure(1)
s = surf(Y,X,mat_data_rib);
s.EdgeColor ="flat";
title("Si3N4 Rib Waveguide")
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
ax = gca;
ax.FontSize = 13;


%range of w0
w0 = 0.1:0.01:10;
%overlap coefficient vector
mu = zeros(1,numel(w0));

%iterate over w0 and calculate the coefficient vector
for j = 1:1:numel(w0)
    gauss_mat = zeros(yd,xd);
    w = w0(j)*um;

    for i = 1:1:numel(X)
        gauss_mat(i) = flat_gauss(X(i)*um, Y(i)*um, w, max_x*um, max_y*um);
    end

    gauss_int = trapz(delta, trapz(delta,abs(gauss_mat).^2,2));

    mu(j) = (abs(trapz(delta, trapz(delta,conj(gauss_mat).*mat_data_rib,2))).^2) / (gauss_int*mode_int);
end

%find maximum mu index
opt_mu_ind = find(mu == max(mu, [], "all"));
%find optimal waist
opt_w = w0(opt_mu_ind);

%return the vector of the gaussian slice parallel to the y axis in 
%correspondece of the peak of the mode
gauss_sil = zeros(1,numel(y));
for i = 1:1:numel(y)
    gauss_sil(i) = flat_gauss(max_x*um, y(i)*um, opt_w*um, max_x*um, max_y*um);
end

%mode slice
data_sil = mat_data_rib(:,max_x_ind)/max(mat_data_rib, [],"all");

%maximum difference
max_perc_difference = max(abs(gauss_sil-data_sil'),[], "all")

%plot gauss and mode 
figure(2)
hold on
plot(y, gauss_sil)
plot(y, data_sil)
title("Mode VS Gaussian")
xlabel("y[\mum]", Interpreter="tex");
ylabel("E/$E_{max}$", Interpreter="tex");
ax = gca;
ax.FontSize = 13;

%mu over w0*2
figure(3)
p = plot(w0*2,mu);
title("\mu over gaussian spot", Interpreter="tex")
xlabel("y[\mum]", Interpreter="tex");
ylabel("\mu", Interpreter="tex");
datatip(p, opt_w*2, max(mu, [], "all"));
ax = gca;
ax.FontSize = 13;