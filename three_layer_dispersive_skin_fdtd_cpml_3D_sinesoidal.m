% FDTD 3D in three layer dispersive skin with CPML
% Author:shayan dodge
% Email address:dodgeshayan@gmail.com
%% initialize the matlab workspace
close all;
clear all;
clc;
%profile on
%v=parpool(32)
%% some constants
mu_0 = 1.2566370614359173e-06;
eps_0= 8.8541878176203892e-12;
c=299792458.0;% speed of light
%% wave definition
waveforms.sinusoidal.frequency=1E10;
T=1/waveforms.sinusoidal.frequency;
lambda=c*T;
omega=2*pi*waveforms.sinusoidal.frequency;
%% FDTD variables
number_of_cells_per_wavelength=29;
dx=lambda/number_of_cells_per_wavelength;
dy=lambda/number_of_cells_per_wavelength;
dz=lambda/number_of_cells_per_wavelength;
totalTime=1000*T;
courant_factor=0.9;
dt=1/(c*sqrt((1/dx^2)+(1/dy^2)+(1/dz^2)));
dt=courant_factor*dt;
totalTimeStep=floor(totalTime/dt);
%% boundary conditions
boundary.type_xn='cpml';
boundary.air_buffer_number_of_cells_xn=0; 
boundary.cpml_number_of_cells_xn=-10;

boundary.type_xp = 'cpml';
boundary.air_buffer_number_of_cells_xp=0;
boundary.cpml_number_of_cells_xp=-10;

boundary.type_yn = 'cpml';
boundary.air_buffer_number_of_cells_yn=0;
boundary.cpml_number_of_cells_yn=-10;

boundary.type_yp = 'cpml';
boundary.air_buffer_number_of_cells_yp=0;
boundary.cpml_number_of_cells_yp=-10;

boundary.type_zn = 'cpml';
boundary.air_buffer_number_of_cells_zn=0;
boundary.cpml_number_of_cells_zn=-10;

boundary.type_zp='cpml';
boundary.air_buffer_number_of_cells_zp=0;
boundary.cpml_number_of_cells_zp=-10;

boundary.cpml_order = 4;
boundary.cpml_sigma_max = 1;
boundary.cpml_kappa_max = 15;
boundary.cpml_alpha_order = 1; 
boundary.cpml_alpha_max = 0.24;
boundary.cpml_eps_R= 1;
%% materialtype
%here define and initialize the arrays of material types
%air
material_type(1).eps_r=1;
material_type(1).mu_r=1;
material_type(1).sigma_e=0;
material_type(1).sigma_m=1e-20;
material_type(1).color=[1 1 1];
material_type(1).debye_kesi=0;


% ************************************
% ************************************
% Contact me to see the complete code.
% Shayan Dodge 
% dodgeshayan@gmail.com
% ************************************
% ************************************

tic
%% run_fdtd_time_marching_loop
for number=1:10^3%totalTimeStep
%% update_magnetic_fields_CPML
Psi_hyx_xn(:,:,:) = cpml_b_mx_xn.*Psi_hyx_xn+cpml_a_mx_xn.*...
    ( Ez (2:n_cpml_xn+1,:,:) - Ez(1:n_cpml_xn,:,:));
Psi_hzx_xn(:,:,:) = cpml_b_mx_xn.* Psi_hzx_xn+cpml_a_mx_xn.*...
    ( Ey (2:n_cpml_xn+1,:,:) - Ey(1:n_cpml_xn ,:,:));
Psi_hyx_xp (:,:,:) = cpml_b_mx_xp.* Psi_hyx_xp+...
    cpml_a_mx_xp.*(Ez (nx-n_cpml_xp+2:nx+1,:,:)-Ez(nx-n_cpml_xp+1:nx,:,:));
Psi_hzx_xp (:,:,:) = cpml_b_mx_xp.* Psi_hzx_xp+...
    cpml_a_mx_xp.*(Ey (nx-n_cpml_xp+2:nx+1,:,:)-Ey(1+nx-n_cpml_xp:nx,:,:));

%% update_magnetic_fields

Hx = ((Chxh.* Hx)+(Chxey .*(Ey (1:nxp1,1:ny,2:nzp1)- Ey (1:nxp1,1:ny,1:nz)))+...
    (Chxez.*(Ez(1:nxp1,2:nyp1,1:nz)- Ez (1:nxp1,1:ny,1:nz))));

%% update_electric_fields_for_CPML
% 
Psi_eyx_xn(:,:,:)= cpml_b_ex_xn.* Psi_eyx_xn(:,:,:)+ cpml_a_ex_xn.*...
    (Hz(2:n_cpml_xn+1,:,:)-Hz(1:n_cpml_xn,:,:));
Psi_ezx_xn(:,:,:)= cpml_b_ex_xn.* Psi_ezx_xn(:,:,:)+ cpml_a_ex_xn.*...
    (Hy(2:n_cpml_xn+1,:,:)-Hy(1:n_cpml_xn,:,:));

Ey (2:n_cpml_xn+1,:,:) = Ey (2:n_cpml_xn+1,:,:)+CPsi_eyx_xn.*Psi_eyx_xn;
Ez (2:n_cpml_xn+1,:,:) = Ez (2:n_cpml_xn+1,:,:)+CPsi_ezx_xn.*Psi_ezx_xn; 


%% update_electric_fields

 Ex(1:nx,2:ny,2:nz)=Cexe(1:nx,2:ny,2:nz).*Ex(1:nx,2:ny,2:nz)+...
     Cexhz(1:nx,2:ny,2:nz).*(Hz(1:nx,2:ny,2:nz)-Hz(1:nx,1:ny-1,2:nz))+...
     Cexhy(1:nx,2:ny,2:nz).*(Hy(1:nx,2:ny,2:nz)-Hy(1:nx,2:ny,1:nz-1));

%% update the additional debye terms
 end
toc
%delete(v);
%profile report
