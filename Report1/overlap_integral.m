%% Compute overlap coefficient between a flat gaussian beam and a waveguide mode

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
[X,Y] = meshgrid(x,y);

%gaussian waist (is the half of the spot)
w0 = 10*um/2;
E0 = 1;

%return the value of a gaussian with waist w0 dislocated by posx, posy at (x,y)
flat_gauss = @(x,y,w,posx,posy) E0 *exp(-((x-posx)^2 + (y-posy)^2)/(w^2));

%reshape the input array into a yd by xd matrix
mat_data_rib = reshape(data_rib(:,3), yd, xd);

%calculate the integral of the mode squared
mode_int_rib = trapz(delta, trapz(delta,abs(mat_data_rib).^2,2));

%find the peak of the waveguide mode
i = find(mat_data_rib == max(mat_data_rib, [], "all"));
max_x = X(i);
max_y = Y(i);

%plot the mode
figure(1);
s = surf(X,Y,mat_data_rib./max(mat_data_rib, [], "all"));
s.EdgeColor ="flat";
title("Si3N4 Rib Waveguide")
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
ax = gca;
ax.FontSize = 13;

%return the gaussian's matrix
gauss_mat = zeros(yd,xd);
for i = 1:1:numel(X)
    gauss_mat(i) = flat_gauss(X(i)*um, Y(i)*um, w0, max_x*um, max_y*um);
end

%plot the gaussian
figure(2)
s = surf(X,Y, gauss_mat./max(gauss_mat, [], "all"));
s.EdgeColor ="flat";
title("" + w0*2/um + "\mum-spot gaussian beam");
xlabel("y[\mum]", Interpreter="tex");
ylabel("x[\mum]", Interpreter="tex");
ax = gca;
ax.FontSize = 13;


%calculate the integral of the gaussian squared 
gauss_int = trapz(delta, trapz(delta,abs(gauss_mat).^2,2));

%calculate the overlapping integral
mu_rib = (abs(trapz(delta, trapz(delta,gauss_mat.*mat_data_rib,2))).^2) / (gauss_int*mode_int_rib)



