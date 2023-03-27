clc
clear all
close all

%% Datos

Nx = 64;
Ny = 256;
Lx = 1;
Ly = 2;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

X = linspace(0, Lx, Nx);
Y = linspace(0, Ly, Ny);

[x, y] = meshgrid(X, Y);

%% Matrices de diferenciacion

ex = ones(Nx,1);
Dxx = (1/dx^2)*spdiags([ex -2*ex ex], [-1 0 1], Nx, Nx);

ey = ones(Ny,1);
Dyy = (1/dy^2)*spdiags([ey -2*ey ey], [-1 0 1], Ny, Ny);
 
L = kron(Dxx, speye(Ny)) + kron(speye(Nx), Dyy);

%% Condiciones de contorno tipo Neumann

bcx = sparse(zeros(Nx)); bcx(1,2) = 1/dx^2; bcx(end,end-1) = 1/dx^2; 
bcy = sparse(zeros(Ny)); bcy(1,2) = 1/dy^2; bcy(end,end-1) = 1/dy^2; 
bc = kron(bcx, speye(Ny)) + kron(speye(Nx), bcy);

%% Solucion

% Pinta solucion inicial
% figure (1),
% subplot(211), pcolor(x, y, u0)
% colormap jet
% shading interp;
% caxis([-2*pi 2*pi])
% colorbar

% Solucion incial
F0 = F(x,y);
F0 = reshape(F0, Ny*Nx, 1);

% Solucion numerica
u = (L + bc)\F0;

% Redistribuye las matrices
u = reshape(u, Ny, Nx);

% Solucion teorica
uTeo = fTeo(x,y);

% Comparacion de soluciones
figure,
subplot(121), pcolor(x, y, u)
title('Num.')  
shading interp
colormap jet
colorbar
set(gca,'layer','top')
box on
%axis square

subplot(122), pcolor(x, y, uTeo)
title('Theo.')
shading interp
colormap jet
colorbar
set(gca,'layer','top')
box on
%axis square

% Error
figure,
subplot(121), pcolor(x, y, abs(uTeo - u)./(1 + abs(uTeo)))
title('Error')
shading interp
colormap jet
colorbar
set(gca,'layer','top')
box on
%axis square

