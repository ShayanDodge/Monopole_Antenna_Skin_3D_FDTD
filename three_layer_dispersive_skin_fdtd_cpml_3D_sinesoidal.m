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

%a dielectric_1
material_type(2).debye_dispersive=true;
material_type(2).eps_r=2.165;
material_type(2).mu_r=1;
material_type(2).sigma_e=7.6;
material_type(2).sigma_m=1e-20;
material_type(2).color=[255/255 20/255 147/255];
material_type(2).debye_eps_s=4.520;
material_type(2).debye_A=1;
material_type(2).debye_tau=7.386;
material_type(2).debye_kesi=material_type(2).debye_A.*eps_0.*(material_type(2).debye_eps_s-material_type(2).eps_r);

%a dielectric_2
material_type(3).debye_dispersive=true;
material_type(3).eps_r=5.12;
material_type(3).mu_r=1;
material_type(3).sigma_e=2.68;
material_type(3).sigma_m=1e-20;
material_type(3).color=[1 0.87 0.88];
material_type(3).debye_eps_s=37.60;
material_type(3).debye_A=1.074;
material_type(3).debye_tau=6.20;
material_type(3).debye_kesi=material_type(3).debye_A.*eps_0.*(material_type(3).debye_eps_s-material_type(3).eps_r);

%a dielectric_3
material_type(4).debye_dispersive=true;
material_type(4).eps_r=3.24;
material_type(4).mu_r=1;
material_type(4).sigma_e=0.7;
material_type(4).sigma_m=1e-20;
material_type(4).color=[1 1 0];
material_type(4).debye_eps_s=4.96;
material_type(4).debye_A=1.074;
material_type(4).debye_tau=2.14;
material_type(4).debye_kesi=material_type(4).debye_A.*eps_0.*(material_type(4).debye_eps_s-material_type(4).eps_r);

%pec
material_type(5).eps_r=1;
material_type(5).mu_r=1;
material_type(5).sigma_e=10^10;
material_type(5).sigma_m=1e-20;
material_type(5).color=[1 1 0];
material_type(5).debye_kesi=0;

%a dielectric_4
material_type(6).eps_r=4;
material_type(6).mu_r=3;
material_type(6).sigma_e=2;
material_type(6).sigma_m=1;
material_type(6).color=[1 0 0];
material_type(6).debye_kesi=0;

%indices of material types defining air, pec, and pmc
material_type_index_air=1;
material_type_index_dielectric_1=2;
material_type_index_dielectric_2=3;
material_type_index_dielectric_3=4;
material_type_index_pec=5;
material_type_index_dielectric_4=6;
number_of_material_type=6;
%% define_geometry
%define a brick

brick(1).min_x=-31e-3;
brick(1).min_y=-31e-3;
brick(1).min_z= 0;
brick(1).max_x= 30e-3;
brick(1).max_y= 30e-3;
brick(1).max_z= 2e-3;
brick(1).material_type=2;

brick(2).min_x=-31e-3;
brick(2).min_y=-31e-3;
brick(2).min_z= 2e-3;
brick(2).max_x= 30e-3;
brick(2).max_y= 30e-3;
brick(2).max_z= 10e-3;
brick(2).material_type=3;

brick(3).min_x=-31e-3;
brick(3).min_y=-31e-3;
brick(3).min_z= 10e-3;
brick(3).max_x= 30e-3;
brick(3).max_y= 30e-3;
brick(3).max_z= 60e-3;
brick(3).material_type=4;

%define a cylinder
cylinders(1).radius=2*10^-3;
cylinders(1).center_x=0;
cylinders(1).center_y=0;
cylinders(1).min_z=-15.*10^-3;
cylinders(1).min_x=cylinders(1).center_x-cylinders(1).radius;
cylinders(1).min_y=cylinders(1).center_y-cylinders(1).radius;
cylinders(1).max_z=-10.*10^-3;
cylinders(1).max_x=cylinders(1).center_x+cylinders(1).radius;
cylinders(1).max_y=cylinders(1).center_y+cylinders(1).radius;
cylinders(1).material_type=5;

%define a cylinder
cylinders(2).radius=3*10^-3;
cylinders(2).center_x=0;
cylinders(2).center_y=0;
cylinders(2).min_z=-15.*10^-3;
cylinders(2).min_x=cylinders(1).center_x-cylinders(1).radius;
cylinders(2).min_y=cylinders(1).center_y-cylinders(1).radius;
cylinders(2).max_z=-10.*10^-3;
cylinders(2).max_x=cylinders(1).center_x+cylinders(1).radius;
cylinders(2).max_y=cylinders(1).center_y+cylinders(1).radius;
cylinders(2).material_type=6;
%% calculating the domain size
number_of_brick=size(brick,2);
number_of_cylinder=size(cylinders,2);
%find the minimum and maximum coordinates of a box encapsulating the object
number_of_objects=1;
for i=1:number_of_cylinder
    min_x(number_of_objects)=cylinders(i).min_x;
    min_y(number_of_objects)=cylinders(i).min_y;
    min_z(number_of_objects)=cylinders(i).min_z;
    max_x(number_of_objects)=cylinders(i).max_x;
    max_y(number_of_objects)=cylinders(i).max_y; 
    max_z(number_of_objects)=cylinders(i).max_z;
    number_of_objects=number_of_objects+1;
end
for i=1:number_of_brick
    min_x(number_of_objects)=brick(i).min_x;
    min_y(number_of_objects)=brick(i).min_y;
    min_z(number_of_objects)=brick(i).min_z;
    max_x(number_of_objects)=brick(i).max_x;
    max_y(number_of_objects)=brick(i).max_y; 
    max_z(number_of_objects)=brick(i).max_z;
    number_of_objects=number_of_objects+1;
end
fdtd_domain.min_x=min(min_x);
fdtd_domain.min_y=min(min_y);    
fdtd_domain.min_z=min(min_z);    
fdtd_domain.max_x=max(max_x);    
fdtd_domain.max_y=max(max_y);    
fdtd_domain.max_z=max(max_z);     
% Determine the problem space boundaries including air buffers
fdtd_domain.min_x = fdtd_domain.min_x-dx *...
    boundary.air_buffer_number_of_cells_xn;
fdtd_domain.min_y = fdtd_domain.min_y-dy *...
    boundary.air_buffer_number_of_cells_yn ;
fdtd_domain.min_z = fdtd_domain.min_z-dz*...
    boundary.air_buffer_number_of_cells_zn ; 
fdtd_domain.max_x = fdtd_domain.max_x+dx *...
    boundary.air_buffer_number_of_cells_xp; 
fdtd_domain.max_y = fdtd_domain. max_y+dy *...
    boundary.air_buffer_number_of_cells_yp ;
fdtd_domain. max_z= fdtd_domain. max_z+dz *...
    boundary.air_buffer_number_of_cells_zp;
% Determine the problem space boundaries including cpml layers
if strcmp (boundary.type_xn,'cpml') &&(boundary.cpml_number_of_cells_xn>0)
    fdtd_domain.min_x = fdtd_domain.min_x -dx *...
        boundary.cpml_number_of_cells_xn;
end
if strcmp (boundary.type_xp, 'cpml') &&(boundary.cpml_number_of_cells_xp >0)
    fdtd_domain.max_x = fdtd_domain. max_x + dx *...
        boundary.cpml_number_of_cells_xp;
end
if strcmp( boundary.type_yn, 'cpml') &&(boundary.cpml_number_of_cells_yn >0)
    fdtd_domain.min_y = fdtd_domain.min_y - dy *...
        boundary.cpml_number_of_cells_yn;
end
if strcmp (boundary.type_yp, 'cpml') &&(boundary.cpml_number_of_cells_yp >0)
   fdtd_domain. max_y = fdtd_domain. max_y+ dy *...
       boundary.cpml_number_of_cells_yp ;
end
if strcmp( boundary.type_zn, 'cpml') &&(boundary.cpml_number_of_cells_zn >0)
    fdtd_domain. min_z = fdtd_domain.min_z - dz *...
        boundary.cpml_number_of_cells_zn ;
end
if strcmp (boundary.type_zp, 'cpml') && (boundary.cpml_number_of_cells_zp>0) 
    fdtd_domain. max_z = fdtd_domain. max_z + dz*...
        boundary.cpml_number_of_cells_zp;
end
%detemining the problem space size
fdtd_domain.size_x=fdtd_domain.max_x-fdtd_domain.min_x;
fdtd_domain.size_y=fdtd_domain.max_y-fdtd_domain.min_y;
fdtd_domain.size_z=fdtd_domain.max_z-fdtd_domain.min_z;
%number of cells in x, y, and z directions
nx=(round(fdtd_domain.size_x/dx));
ny=round(fdtd_domain.size_y/dy);
nz=round(fdtd_domain.size_z/dz);
%adjust domain size by snapping to cells
fdtd_domain.size_x=nx*dx;
fdtd_domain.size_y=ny*dy;
fdtd_domain.size_z=nz*dz;
fdtd_domain.max_x=fdtd_domain.min_x+fdtd_domain.size_x;
fdtd_domain.max_y=fdtd_domain.min_y+fdtd_domain.size_y;
fdtd_domain.max_z=fdtd_domain.min_z+fdtd_domain.size_z;
%some frequently used auxiliary parametrs
nxp1=nx+1;  nyp1=ny+1;  nzp1=nz+1;
nxm1=nx-1;  nym1=ny-1;  nzm1=nz-1;
nxm2=nx-2;  nym2=ny-2;  nzm2=nz-2;
%create arrays storing the center coordinates of the cells in
fdtd_domain.cell_center_coordinates_x=zeros(nx,ny,nz,'gpuArray');
fdtd_domain.cell_center_coordinates_y=zeros(nx,ny,nz,'gpuArray');
fdtd_domain.cell_center_coordinates_z=zeros(nx,ny,nz,'gpuArray');
for ind =1:nx
    fdtd_domain.cell_center_coordinates_x(ind,:,:)=...
       ((ind-0.5)*dx+fdtd_domain.min_x);
end
for ind =1:ny
    fdtd_domain.cell_center_coordinates_y(:,ind,:)=...
        (ind-0.5)*dy+fdtd_domain.min_y;
end
for ind =1:nz
    fdtd_domain.cell_center_coordinates_z(:,:,ind)=...
        (ind-0.5)*dz+fdtd_domain.min_z;
end
cx=fdtd_domain.cell_center_coordinates_x;
cy=fdtd_domain.cell_center_coordinates_y;
cz=fdtd_domain.cell_center_coordinates_z;
xcoor=linspace (fdtd_domain. min_x, fdtd_domain.max_x, nxp1);
ycoor=linspace (fdtd_domain. min_y, fdtd_domain.max_y, nyp1);
zcoor=linspace (fdtd_domain. min_z, fdtd_domain.max_z, nzp1);

%% material_3d_space
material_3d_space=ones(nx,ny,nz);
%% creating_brick
%creat the 3d object in problem space by
for ind=1:number_of_brick
%%convert brick end coordinates to node indices
blx=round((brick(ind).min_x-fdtd_domain.min_x)/dx)+1;
bly=round((brick(ind).min_y-fdtd_domain.min_y)/dy)+1;
blz=round((brick(ind).min_z-fdtd_domain.min_z)/dz)+1;
bux=round((brick(ind).max_x-fdtd_domain.min_x)/dx)+1;
buy=round((brick(ind).max_y-fdtd_domain.min_y)/dy)+1;
buz=round((brick(ind).max_z-fdtd_domain.min_z)/dz)+1;
%%assign material type of brick to the cells
material_3d_space(blx:bux-1,bly:buy-1,blz:buz-1)=brick(ind).material_type;
end

for in=1:nz 
distance(:,:,in)=sqrt((cylinders(1).center_y-cy(:,:,in)).^2+...
    (cylinders(1).center_x-cx(:,:,in)).^2);
end
l=find( distance>=cylinders(1).radius&distance<=1.5.*cylinders(1).radius);
material_3d_space(l)=cylinders(1).material_type;

for in=1:nz 
distance(:,:,in)=sqrt((cylinders(2).center_y-cy(:,:,in)).^2+...
    (cylinders(2).center_x-cx(:,:,in)).^2);
end
l=find(distance<=3.*cylinders(2).radius& distance>=1.5.*cylinders(2).radius);
 material_3d_space(l)=cylinders(2).material_type;

kesi=(ones(nxp1,nyp1,nzp1,'gpuArray'));
eps_r_x=(ones(nx,nyp1,nzp1,'gpuArray'));
eps_r_y=(ones(nxp1,ny,nzp1,'gpuArray'));
eps_r_z=(ones(nxp1,nyp1,nz,'gpuArray'));
mu_r_x=(ones(nxp1,ny,nz,'gpuArray'));
mu_r_y=(ones(nx,nyp1,nz,'gpuArray'));
mu_r_z=(ones(nx,ny,nzp1,'gpuArray'));
sigma_e_x=(zeros(nx,nyp1,nzp1,'gpuArray'));
sigma_e_y=(zeros(nxp1,ny,nzp1,'gpuArray'));
sigma_e_z=(zeros(nxp1,nyp1,nz,'gpuArray'));
sigma_m_x=(zeros(nxp1,ny,nz,'gpuArray'));
sigma_m_y=(zeros(nx,nyp1,nz,'gpuArray'));
sigma_m_z=(zeros(nx,ny,nzp1,'gpuArray'));
%% calculate_material_component_values
%calculate material component values by averaging
for ind=1:size(material_type,2)
    t_eps_r(ind)=material_type(ind).eps_r;
    t_mu_r(ind)=material_type(ind).mu_r;
    t_sigma_e(ind)=material_type(ind).sigma_e;
    t_sigma_m(ind)=material_type(ind).sigma_m;
    t_kesi(ind)=material_type(ind).debye_kesi;
end
%assign negligibly small values to t_mu_r and t_sigma_m where they are zero
%in order to prevent division by zerro error
t_mu_r(find(t_mu_r==0))=1e-20;
t_sigma_m(find(t_sigma_m==0.0000))=1e-20;
disp('calculating_kesi')
kesi(1:nx,1:ny,1:nz)=t_kesi(material_3d_space(1:nx,1:ny,1:nz));
disp('calculating_eps_r_x');
%eps_r_x(i,j,k)is average of four cells(i,j,k)(i,j-1,k)(i,j,k-1)(i,j-1,k-1)
eps_r_x(1:nx,2:ny,2:nz)=0.25*(t_eps_r(material_3d_space(1:nx,2:ny,2:nz))+...
    t_eps_r(material_3d_space(1:nx,1:ny-1,2:nz))+...
    t_eps_r(material_3d_space(1:nx,2:ny,1:nz-1))+...
    t_eps_r(material_3d_space(1:nx,1:ny-1,1:nz-1)));
disp('calculating_eps_r_y');
%%eps_r_y(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j,k-1)(i-1,j,k-1)
eps_r_y(2:nx,1:ny,2:nz)=0.25*(t_eps_r(material_3d_space(2:nx,1:ny,2:nz))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny,2:nz))+...
    t_eps_r(material_3d_space(2:nx,1:ny,1:nz-1))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny,1:nz-1)));
disp('calculating_eps_r_z');
%eps_r_z(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j-1,k)(i-1,j-1,k)
eps_r_z(2:nx,2:ny,1:nz)=0.25*(t_eps_r(material_3d_space(2:nx,1:ny-1,1:nz))+...
    t_eps_r(material_3d_space(1:nx-1,2:ny,1:nz))+...
    t_eps_r(material_3d_space(2:nx,2:ny,1:nz))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny-1,1:nz)));
disp('calculating_sigma_e_x');
%%sigma_e_x(i,j,k)is average of four cells(i,j,k)(i,j-1,k)(i,j,k-1)(i,j-1,k-1)
sigma_e_x(1:nx,2:ny,2:nz)=0.25*(t_sigma_e(material_3d_space(1:nx,2:ny,2:nz))+...
    t_sigma_e(material_3d_space(1:nx,1:ny-1,2:nz))+...
    t_sigma_e(material_3d_space(1:nx,2:ny,1:nz-1))+...
    t_sigma_e(material_3d_space(1:nx,1:ny-1,1:nz-1)));
disp('calculating_sigma_e_y');
%%sigma_e_y(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j,k-1)(i-1,j,k-1)
sigma_e_y(2:nx,1:ny,2:nz)=0.25*(t_sigma_e(material_3d_space(2:nx,1:ny,2:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny,2:nz))+...
    t_sigma_e(material_3d_space(2:nx,1:ny,1:nz-1))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny,1:nz-1)));
disp('calculating_sigma_e_z');
%sigma_e_z(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j-1,k)(i-1,j-1,k)
sigma_e_z(2:nx,2:ny,1:nz)=0.25*(t_sigma_e(material_3d_space(2:nx,1:ny-1,1:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,2:ny,1:nz))+...
    t_sigma_e(material_3d_space(2:nx,2:ny,1:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny-1,1:nz)));
disp('calculating_sigma_m_x');
%%sigma_m_x(i,j,k)is average of two cells(i,j,k)(i-1,j,k)
sigma_m_x(2:nx,1:ny,1:nz)=0.5*(t_sigma_m(material_3d_space(2:nx,1:ny,1:nz))+...
    t_sigma_m(material_3d_space(1:nx-1,1:ny,1:nz)));
disp('calculating_sigma_m_y');
%sigma_e_y(i,j,k)is average of two cells(i,j,k)(i,j-1,k)
sigma_m_y(1:nx,2:ny,1:nz)=0.5*(t_sigma_m(material_3d_space(1:nx,2:ny,1:nz))+...
    t_sigma_m(material_3d_space(1:nx,1:ny-1,1:nz)));
disp('calculating_sigma_m_z');
%sigma_e_z(i,j,k)is average of two cells(i,j,k)(i,j,k-1)
sigma_m_z(1:nx,1:ny,2:nz)=0.5*(t_sigma_m(material_3d_space(1:nx,1:ny,2:nz))+...
    t_sigma_m(material_3d_space(1:nx,1:ny,1:nz-1)));
disp('calculating_mu_r_x');
%mu_r_x(i,j,k)is average of two cells(i,j,k)(i-1,j,k)
mu_r_x(2:nx,1:ny,1:nz)=0.5*(t_mu_r(material_3d_space(2:nx,1:ny,1:nz))+...
    t_mu_r(material_3d_space(1:nx-1,1:ny,1:nz)));
disp('calculating_mu_r_y');
%mu_r_y(i,j,k)is average of two cells(i,j,k)(i,j-1,k)
mu_r_y(1:nx,2:ny,1:nz)=0.5*(t_mu_r(material_3d_space(1:nx,2:ny,1:nz))+...
    t_mu_r(material_3d_space(1:nx,1:ny-1,1:nz)));
disp('calculating_mu_r_z');
%mu_r_z(i,j,k)is average of two cells(i,j,k)(i,j,k-1)
mu_r_z(1:nx,1:ny,2:nz)=0.5*(t_mu_r(material_3d_space(1:nx,1:ny,2:nz))+...
    t_mu_r(material_3d_space(1:nx,1:ny,1:nz-1)));
% Arrays to store material indices distribution of dispersive objects 
dispersive_Ex = gpuArray(zeros(nx,nyp1,nzp1)); 
dispersive_Ey = gpuArray(zeros(nxp1,ny,nzp1)); 
dispersive_Ez = gpuArray(zeros(nxp1,nyp1,nz));

for ind = 1:number_of_brick
mat_ind = brick(ind).material_type;
material = material_type(mat_ind); 
% convert brick end coordinates to node indices
is =round((brick(ind).min_x-fdtd_domain.min_x)/dx)+1;
js =round((brick(ind).min_y-fdtd_domain.min_y)/dy)+1;
ks=round((brick(ind).min_z-fdtd_domain.min_z)/dz)+1;
ie =round((brick(ind).max_x-fdtd_domain.min_x)/dx)+1;
je =round((brick(ind).max_y-fdtd_domain.min_y)/dy)+1; 
ke=round((brick(ind).max_z-fdtd_domain.min_z)/dz)+1;
dispersive_Ex (is:ie-1, js:je, ks:ke) = mat_ind;
dispersive_Ey (is:ie, js:je-1, ks:ke) = mat_ind;
dispersive_Ez (is:ie, js:je, ks: ke-1) = mat_ind; 
end 
% create a new structure array to store indices of field components 
% that need to be updated using Debye formulation
number_of_debye = 0;
for ind = 1:number_of_material_type
material = material_type(ind);
if material.debye_dispersive==1
number_of_debye = number_of_debye + 1;
debye(number_of_debye).material = ind;
I= find (dispersive_Ex == ind);
debye(number_of_debye).Ex_indices = I; 
eps_r_x(I) = material.eps_r; 
sigma_e_x(I) = material.sigma_e;
I = find (dispersive_Ey == ind);
debye(number_of_debye). Ey_indices = I;
eps_r_y(I) = material.eps_r;
sigma_e_y(I) = material.sigma_e;
I = find (dispersive_Ez == ind); 
debye(number_of_debye).Ez_indices = I;
eps_r_z(I) = material.eps_r;
sigma_e_z(I) = material.sigma_e;
end
end
clear dispersive_Ex dispersive_Ey dispersive_Ez

%% creating_field_array
%create and initialize field and current arrays
Hx=(zeros(nxp1,ny,nz,'gpuArray'));
Hy=(zeros(nx,nyp1,nz,'gpuArray'));
Hz=(zeros(nx,ny,nzp1,'gpuArray'));
Ex=(zeros(nx,nyp1,nzp1,'gpuArray'));
Ey=(zeros(nxp1,ny,nzp1,'gpuArray'));
Ez=(zeros(nxp1,nyp1,nz,'gpuArray'));
Exp=(zeros(nx,nyp1,nzp1,'gpuArray'));
Eyp=(zeros(nxp1,ny,nzp1,'gpuArray'));
Ezp=(zeros(nxp1,nyp1,nz,'gpuArray'));
%% initialize_updating_coefficients
% Coeffiecients updating Ex 
Cexe=(2*eps_r_x*eps_0-dt*sigma_e_x)./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexhz=(2*dt/dy)./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexhy=-(2*dt/dz)./(2*eps_r_x*eps_0+dt*sigma_e_x);
% Coeffiecients updating Ey 
Ceye=(2*eps_r_y*eps_0-dt*sigma_e_y)./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyhx=(2*dt/dz)./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyhz=-(2*dt/dx)./(2*eps_r_y*eps_0+dt*sigma_e_y);
% Coeffiecients updating Ez
Ceze =(2*eps_r_z*eps_0-dt*sigma_e_z)./(2*eps_r_z.*eps_0+dt*sigma_e_z);
Cezhy=(2*dt/dx)./(2*eps_r_z*eps_0+dt*sigma_e_z);
Cezhx=-(2*dt/dy)./(2*eps_r_z.*eps_0+dt*sigma_e_z);
%general magnetic field updating coefficients
% Coeffiecients updating Hx
Chxh =((2*mu_r_x*mu_0-dt*sigma_m_x)./(2*mu_r_x*mu_0+dt*sigma_m_x));
Chxez=(-(2*dt/dy)./(2*mu_r_x*mu_0+dt*sigma_m_x));
Chxey=((2*dt/dz)./(2*mu_r_x*mu_0+dt*sigma_m_x));
% Coeffiecients updating Hy
Chyh=(2*mu_r_y*mu_0-dt*sigma_m_y)./(2*mu_r_y.*mu_0+dt*sigma_m_y);
Chyex=(-2*dt/dz)./(2*mu_r_y*mu_0+dt*sigma_m_y);
Chyez=(2*dt/dx)./(2*mu_r_y.*mu_0+dt*sigma_m_y);
% Coeffiecients updating Hz
Chzh =(2*mu_r_z*mu_0-dt*sigma_m_z)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzey=-(2*dt/dx)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzex=(2*dt/dy)./(2*mu_r_z*mu_0+dt*sigma_m_z);
%%initialize_dispersive_coefficients
% initialize Debye updating coefficients
for ind = 1:number_of_debye
material = material_type(debye(ind).material);
n_elements = size (debye(ind).Ex_indices,1); 
debye(ind).J_x = zeros(n_elements,1);
n_elements = size(debye(ind). Ey_indices,1); 
debye(ind). J_y = zeros(n_elements,1); 
n_elements = size(debye(ind). Ez_indices,1); 
debye(ind).J_z = zeros(n_elements,1);
kesi = eps_0 *material.debye_A* (material.debye_eps_s -material.eps_r);
debye(ind).C_jj = (2.*material_type(debye(ind).material).debye_tau - dt)./(2.*material_type(debye(ind).material).debye_tau +dt); 
debye(ind).C_je = 2.*kesi./(2.*material_type(debye(ind).material).debye_tau+dt); 
end
%% initialize_boundary_conditions_3d
% define logical parameters for the conditions that will be used often
n_cpml_xn = abs( boundary.cpml_number_of_cells_xn);
n_cpml_xp = abs( boundary.cpml_number_of_cells_xp);
n_cpml_yn = abs( boundary.cpml_number_of_cells_yn);
n_cpml_yp = abs( boundary.cpml_number_of_cells_yp);
n_cpml_zn= abs ( boundary.cpml_number_of_cells_zn); 
n_cpml_zp = abs( boundary.cpml_number_of_cells_zp);
% Call CPML initialization routine if any side is CPML
% Initialize CPML boundary condition
pml_order = boundary.cpml_order;% order of the polynomial distribution 
sigma_max = boundary.cpml_sigma_max;
kappa_max = boundary.cpml_kappa_max;
alpha_order = boundary.cpml_alpha_order;
alpha_max = boundary.cpml_alpha_max;
eps_R = boundary.cpml_eps_R;
% Initialize cpml for xn region
sigma_opt=sigma_max*(n_cpml_xn+1)/(sqrt(eps_R)*150*pi*dx);
rho_e=((n_cpml_xn:-1:1)-0.75)/n_cpml_xn;
rho_m=((n_cpml_xn:-1:1)-0.25)/n_cpml_xn;
sigma_ex_xn=sigma_opt*abs(rho_e).^pml_order;
sigma_mx_xn=sigma_opt*abs(rho_m).^pml_order;
kappa_ex_xn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mx_xn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ex_xn=alpha_max*abs(rho_e).^pml_order;
alpha_mx_xn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ex_xn=exp((-dt/eps_0)...
    *((sigma_ex_xn./kappa_ex_xn)+alpha_ex_xn));
cpml_a_ex_xn=(1/dx)*(cpml_b_ex_xn-1).*sigma_ex_xn ...
    ./(kappa_ex_xn.*(sigma_ex_xn+kappa_ex_xn.*alpha_ex_xn));
cpml_b_mx_xn=exp((-dt/eps_0)...
    *((sigma_mx_xn./kappa_mx_xn)+alpha_mx_xn));
cpml_a_mx_xn=(1/dx)*(cpml_b_mx_xn-1).*sigma_mx_xn ...
    ./(kappa_mx_xn.*(sigma_mx_xn+kappa_mx_xn.*alpha_mx_xn));

Psi_eyx_xn = gpuArray(zeros (n_cpml_xn , ny, nzp1));
Psi_ezx_xn = gpuArray(zeros (n_cpml_xn , nyp1, nz));
Psi_hyx_xn = gpuArray(zeros (n_cpml_xn , nyp1, nz));
Psi_hzx_xn = gpuArray(zeros (n_cpml_xn , ny, nzp1 ));

CPsi_eyx_xn = gpuArray(Ceyhz (2:n_cpml_xn +1,:,:)* dx);
CPsi_ezx_xn = gpuArray(Cezhy (2:n_cpml_xn +1,:,:)* dx);
CPsi_hyx_xn = gpuArray(Chyez (1:n_cpml_xn , : , :)*dx);
CPsi_hzx_xn = gpuArray(Chzey (1:n_cpml_xn,:,:) * dx) ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_xn
Ceyhz (i+1,:,:) = Ceyhz (i+1,:,:)/ kappa_ex_xn (i); 
Cezhy (i+1,:,:) = Cezhy(i+1,:,:)/kappa_ex_xn (i);
Chyez(i,:,:) = Chyez(i,:,:)/ kappa_mx_xn(i);
Chzey(i,:,:) = Chzey(i, :,:)/ kappa_mx_xn(i);
end
% Initialize cpml for xp region
sigma_opt=sigma_max*(n_cpml_xp+1)/(sqrt(eps_R)*150*pi*dx);
rho_e=((1:1:n_cpml_xp)-0.75)/n_cpml_xp;
rho_m=((1:1:n_cpml_xp)-0.25)/n_cpml_xp;
sigma_ex_xp=sigma_opt*abs(rho_e).^pml_order;
sigma_mx_xp=sigma_opt*abs(rho_m).^pml_order;
kappa_ex_xp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mx_xp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ex_xp=alpha_max*abs(rho_e).^pml_order;
alpha_mx_xp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ex_xp=exp((-dt/eps_0)...
    *((sigma_ex_xp./kappa_ex_xp)+alpha_ex_xp));
cpml_a_ex_xp=(1/dx)*(cpml_b_ex_xp-1).*sigma_ex_xp ...
    ./(kappa_ex_xp.*(sigma_ex_xp+kappa_ex_xp.*alpha_ex_xp));
cpml_b_mx_xp=exp((-dt/eps_0)...
    *((sigma_mx_xp./kappa_mx_xp)+alpha_mx_xp));
cpml_a_mx_xp=1/dx*(cpml_b_mx_xp-1).*sigma_mx_xp ...
    ./(kappa_mx_xp.*(sigma_mx_xp+kappa_mx_xp.*alpha_mx_xp));

Psi_eyx_xp = zeros (n_cpml_xp , ny, nzp1,'gpuArray');
Psi_ezx_xp = zeros (n_cpml_xp , nyp1, nz,'gpuArray');
Psi_hyx_xp = zeros (n_cpml_xp , nyp1, nz,'gpuArray');
Psi_hzx_xp = zeros (n_cpml_xp , ny, nzp1 ,'gpuArray');

CPsi_eyx_xp = gpuArray(Ceyhz (nxp1-n_cpml_xp:nx,:,:)* dx);
CPsi_ezx_xp = gpuArray(Cezhy (nxp1-n_cpml_xp:nx,:,:)* dx);
CPsi_hyx_xp = gpuArray(Chyez (nxp1-n_cpml_xp:nx,:,:)* dx);
CPsi_hzx_xp = gpuArray(Chzey (nxp1-n_cpml_xp:nx,:,:)* dx) ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1:n_cpml_xp
Ceyhz (nx-n_cpml_xp+i,:,:) = Ceyhz(nx-n_cpml_xp+i,:,:)/ kappa_ex_xp (i);
Cezhy (nx-n_cpml_xp+i,:,:) = Cezhy(nx-n_cpml_xp+i,:,:)/ kappa_ex_xp (i);
Chyez (nx-n_cpml_xp+i,:,:) = Chyez(nx-n_cpml_xp+i,:,:)/ kappa_mx_xp (i);
Chzey (nx-n_cpml_xp+i,:,:) = Chzey(nx-n_cpml_xp+i,:,:)/ kappa_mx_xp (i);
end

% Initialize cpml for yn region
sigma_opt=sigma_max*(n_cpml_yn+1)/(sqrt(eps_R)*150*pi*dy);
rho_e=((n_cpml_yn:-1:1)-0.75)/n_cpml_yn;
rho_m=((n_cpml_yn:-1:1)-0.25)/n_cpml_yn;
sigma_ey_yn=sigma_opt*abs(rho_e).^pml_order;
sigma_my_yn=sigma_opt*abs(rho_m).^pml_order;
kappa_ey_yn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_my_yn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ey_yn=alpha_max*abs(rho_e).^pml_order;
alpha_my_yn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ey_yn=exp((-dt/eps_0)...
    *((sigma_ey_yn./kappa_ey_yn)+alpha_ey_yn));
cpml_a_ey_yn=1/dy*(cpml_b_ey_yn-1).*sigma_ey_yn ...
    ./(kappa_ey_yn.*(sigma_ey_yn+kappa_ey_yn.*alpha_ey_yn));
cpml_b_my_yn=exp((-dt/eps_0)...
    *((sigma_my_yn./kappa_my_yn)+alpha_my_yn));
cpml_a_my_yn=1/dy*(cpml_b_my_yn-1).*sigma_my_yn ...
    ./(kappa_my_yn.*(sigma_my_yn+kappa_my_yn.*alpha_my_yn));


Psi_exy_yn = zeros (nx,n_cpml_yn,nzp1,'gpuArray'); 
Psi_ezy_yn = zeros (nxp1,n_cpml_yn,nz,'gpuArray');
Psi_hxy_yn = zeros (nxp1,n_cpml_yn,nz,'gpuArray');
Psi_hzy_yn = zeros(nx,n_cpml_yn, nzp1 ,'gpuArray');
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_exy_yn = Cexhz (:,2:n_cpml_yn+1,:)*dy;
CPsi_ezy_yn = Cezhx (:,2:n_cpml_yn+1,:)*dy;
CPsi_hxy_yn = Chxez (:,1:n_cpml_yn  ,:)*dy;
CPsi_hzy_yn = Chzex (:,1:n_cpml_yn  ,:)*dy;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for j = 1: n_cpml_yn
Cexhz (:,j+1,:) = Cexhz (:,j+1,:)/ kappa_ey_yn (j);
Cezhx (:,j+1,:) = Cezhx (:,j+1,:)/ kappa_ey_yn (j);
Chxez (:,j  ,:) = Chxez (:,j  ,:)/ kappa_my_yn (j);
Chzex (:,j  ,:) = Chzex (:,j  ,:)/ kappa_my_yn (j);
end
% Initialize cpml for yp region
sigma_opt=sigma_max*(n_cpml_yp+1)/(sqrt(eps_R)*150*pi*dy);
rho_e=((1:1:n_cpml_yp)-0.75)/n_cpml_yp;
rho_m=((1:1:n_cpml_yp)-0.25)/n_cpml_yp;
sigma_ey_yp=sigma_opt*abs(rho_e).^pml_order;
sigma_my_yp=sigma_opt*abs(rho_m).^pml_order;
kappa_ey_yp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_my_yp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ey_yp=alpha_max*abs(rho_e).^pml_order;
alpha_my_yp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ey_yp=exp((-dt/eps_0)...
    *((sigma_ey_yp./kappa_ey_yp)+alpha_ey_yp));
cpml_a_ey_yp=1/dx*(cpml_b_ey_yp-1).*sigma_ey_yp ...
    ./(kappa_ey_yp.*(sigma_ey_yp+kappa_ey_yp.*alpha_ey_yp));
cpml_b_my_yp=exp((-dt/eps_0)...
    *((sigma_my_yp./kappa_my_yp)+alpha_my_yp));
cpml_a_my_yp=1/dy*(cpml_b_my_yp-1).*sigma_my_yp ...
    ./(kappa_my_yp.*(sigma_my_yp+kappa_my_yp.*alpha_my_yp));



Psi_exy_yp = zeros (nx,n_cpml_yp,nzp1,'gpuArray'); 
Psi_ezy_yp = zeros (nxp1,n_cpml_yp,nz,'gpuArray');
Psi_hxy_yp = zeros (nxp1,n_cpml_yp,nz,'gpuArray');
Psi_hzy_yp = zeros (nx,n_cpml_yp,nzp1,'gpuArray');
% Create and initialize 2D cpml convolution coefficients
% Notice that Ey (nxp1,:,:) and Ez (nxp1 ,:,:) are not updated by cmpl 
CPsi_exy_yp = Cexhz (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_ezy_yp = Cezhx (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_hxy_yp = Chxez (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_hzy_yp = Chzex (:,nyp1-n_cpml_yp:ny,:)* dy;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (nxp1,:,:) and Ez (nxp1 ,:,:) are not updated by cmpl 
for j = 1:n_cpml_yp
Cexhz (:,ny-n_cpml_yp+j,:) = Cexhz (:,ny-n_cpml_yp+j,:)/ kappa_ey_yp (j);
Cezhx (:,ny-n_cpml_yp+j,:) = Cezhx (:,ny-n_cpml_yp+j,:)/ kappa_ey_yp (j);
Chxez (:,ny-n_cpml_yp+j,:) = Chxez (:,ny-n_cpml_yp+j,:)/ kappa_my_yp (j);
Chzex (:,ny-n_cpml_yp+j,:) = Chzex (:,ny-n_cpml_yp+j,:)/ kappa_my_yp (j);
end

% Initialize cpml for zn region
sigma_opt=sigma_max*(n_cpml_zn+1)/(sqrt(eps_R)*150*pi*dz);
rho_e=((n_cpml_zn:-1:1)-0.75)/n_cpml_zn;
rho_m=((n_cpml_zn:-1:1)-0.25)/n_cpml_zn;
sigma_ez_zn=sigma_opt*abs(rho_e).^pml_order;
sigma_mz_zn=sigma_opt*abs(rho_m).^pml_order;
kappa_ez_zn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mz_zn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ez_zn=alpha_max*abs(rho_e).^pml_order;
alpha_mz_zn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ez_zn=exp((-dt/eps_0)...
    *((sigma_ez_zn./kappa_ez_zn)+alpha_ez_zn));
cpml_a_ez_zn=1/dz*(cpml_b_ez_zn-1).*sigma_ez_zn ...
    ./(kappa_ez_zn.*(sigma_ez_zn+kappa_ez_zn.*alpha_ez_zn));
cpml_b_mz_zn=exp((-dt/eps_0)...
    *((sigma_mz_zn./kappa_mz_zn)+alpha_mz_zn));
cpml_a_mz_zn=1/dz*(cpml_b_mz_zn-1).*sigma_mz_zn ...
    ./(kappa_mz_zn.*(sigma_mz_zn+kappa_mz_zn.*alpha_mz_zn));

Psi_eyz_zn = zeros (nxp1,ny,n_cpml_zn,'gpuArray');
Psi_exz_zn = zeros (nx,nyp1,n_cpml_zn,'gpuArray');
Psi_hyz_zn = zeros (nx,nyp1,n_cpml_zn,'gpuArray');
Psi_hxz_zn = zeros (nxp1,ny,n_cpml_zn ,'gpuArray');
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_eyz_zn = Ceyhx (:,:,2: n_cpml_zn+1)* dz;
CPsi_exz_zn = Cexhy (:,:,2: n_cpml_zn+1)* dz;
CPsi_hyz_zn = Chyex (:,:,1: n_cpml_zn)* dz;
CPsi_hxz_zn = Chxey (:,:,1: n_cpml_zn)* dz ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for j = 1: n_cpml_zn
Cexhy (:,:,j+1) = Cexhy (:,:,j+1)/ kappa_ez_zn (j); 
Ceyhx (:,:,j+1) = Ceyhx(:,:,j+1)/kappa_ez_zn (j);
Chxey (:,:,j) = gpuArray(Chxey(:,:,j)/ kappa_mz_zn(j));
Chyex (:,:,j) = Chyex(:,:,j)/ kappa_mz_zn(j);
end
% Initialize cpml for zp region
sigma_opt=sigma_max*(n_cpml_zp+1)/(sqrt(eps_R)*150*pi*dz);
rho_e=((1:1:n_cpml_zp)-0.75)/n_cpml_zp;
rho_m=((1:1:n_cpml_zp)-0.25)/n_cpml_zp;
sigma_ez_zp=sigma_opt*abs(rho_e).^pml_order;
sigma_mz_zp=sigma_opt*abs(rho_m).^pml_order;
kappa_ez_zp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mz_zp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ez_zp=alpha_max*abs(rho_e).^pml_order;
alpha_mz_zp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ez_zp=exp((-dt/eps_0)...
    *((sigma_ez_zp./kappa_ez_zp)+alpha_ez_zp));
cpml_a_ez_zp=1/dz*(cpml_b_ez_zp-1).*sigma_ez_zp ...
    ./(kappa_ez_zp.*(sigma_ez_zp+kappa_ez_zp.*alpha_ez_zp));
cpml_b_mz_zp=exp((-dt/eps_0)...
    *((sigma_mz_zp./kappa_mz_zp)+alpha_mz_zp));
cpml_a_mz_zp=1/dz*(cpml_b_mz_zp-1).*sigma_mz_zp ...
    ./(kappa_mz_zp.*(sigma_mz_zp+kappa_mz_zp.*alpha_mz_zp));
Psi_eyz_zp = zeros (nxp1,ny,n_cpml_zp,'gpuArray');
Psi_exz_zp = zeros (nx,nyp1,n_cpml_zp,'gpuArray');
Psi_hyz_zp = zeros (nx,nyp1,n_cpml_zp,'gpuArray');
Psi_hxz_zp = zeros (nxp1,ny,n_cpml_zp,'gpuArray');
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_eyz_zp = Ceyhx (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_exz_zp = Cexhy (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_hyz_zp = Chyex (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_hxz_zp = Chxey (:,:,nzp1-n_cpml_zp:nz)* dz;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_zp
Cexhy(:,:,nz-n_cpml_zp+i) = Cexhy(:,:,nz-n_cpml_zp+i)/kappa_ez_zp(i); 
Ceyhx(:,:,nz-n_cpml_zp+i) = Ceyhx(:,:,nz-n_cpml_zp+i)/kappa_ez_zp(i);
Chxey(:,:,nz-n_cpml_zp+i) = gpuArray(Chxey(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i));
Chyex(:,:,nz-n_cpml_zp+i) = Chyex(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i);
end
%% source
r=1.*10^-3;R=2.*10^-3;
amptidute=6000;
source.center_x=0;
source.center_y=0;
source.center_z=dz.*(abs(boundary.cpml_number_of_cells_zn)+2);
source.center_z_in=abs(boundary.cpml_number_of_cells_zn)+2;
source.Radius=10^-3;
L_Ey=0.5.*log(dy./source.Radius);
L_Ex=0.5.*log(dx./source.Radius);
S_Hy=0.5.*log(dx./source.Radius);
S_Hx=0.5.*log(dy./source.Radius);
a_L_Ex=amptidute./L_Ex;
a_L_Ey=amptidute./L_Ey;
% Coeffiecients updating Ez
Ceze((nx+1)/2,(ny+1)/2,n_cpml_zn+2)=0;
Cezhy((nx+1)/2,(ny+1)/2,n_cpml_zn+2)=0;
Cezhx((nx+1)/2,(ny+1)/2,n_cpml_zn+2)=0;
% Coeffiecients updating Hx
Chxh((nx+1)/2,(ny+1)/2:((ny+1)/2)-1,n_cpml_zn+2)=gpuArray(1);
Chxez((nx+1)/2,(ny+1)/2:((ny+1)/2)-1,n_cpml_zn+2)=gpuArray(dt./(mu_0.*mu_r_x((nx+1)/2,(ny+1)/2:((ny+1)/2)-1,n_cpml_zn+2).*dy.*S_Hx));
Chxey((nx+1)/2,(ny+1)/2:((ny+1)/2)-1,n_cpml_zn+2)=gpuArray(dt.*L_Ey./(mu_0.*mu_r_x((nx+1)/2,(ny+1)/2:((ny+1)/2)-1,n_cpml_zn+2).*dz.*S_Hx));
% Coeffiecients updating Hy
Chyh((nx+1)./2:((nx+1)./2)-1,(ny+1)/2,n_cpml_zn+2)=1;
Chyex((nx+1)./2:((nx+1)./2)-1,(ny+1)/2,n_cpml_zn+2)=dt.*L_Ex./(mu_0.*mu_r_y((nx+1)/2:((nx+1)./2)-1,(ny+1)/2,n_cpml_zn+2).*dz.*S_Hy);
Chyez((nx+1)./2:((nx+1)./2)-1,(ny+1)/2,n_cpml_zn+2)=dt./(mu_0.*mu_r_y((nx+1)/2:((nx+1)/2)-1,(ny+1)/2,n_cpml_zn+2).*dx.*S_Hy);
%% initialize_output_parameters_3d
Hx_sample=zeros(nxp1,ny,nz)  ;
Hy_sample=zeros(nx,nyp1,nz)  ;
Hz_sample=zeros(nx,ny,nzp1)  ;
Ex_sample=zeros(nx,nyp1,nzp1);
Ey_sample=zeros(nxp1,ny,nzp1);
Ez_sample=zeros(nxp1,nyp1,nz);
Ex_sample_2D=zeros(1,nzp1);
Ez_sample_2D=zeros(1,nz);

 h_1=figure(1)
for ind=1:number_of_brick
hold on
patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x ...
    brick(ind).max_x],[brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y brick(ind).min_y],[brick(ind).min_z...
    brick(ind).min_z brick(ind).min_z brick(ind).min_z],...
    material_type(ind+1).color)
patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x], [brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y brick(ind).min_y], [brick(ind).max_z...
    brick(ind).max_z brick(ind).max_z brick(ind).max_z],...
    material_type(ind+1).color)
patch([brick(ind).min_x brick(ind).min_x brick(ind).min_x...
    brick(ind).min_x], [brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y brick(ind).min_y], [brick(ind).min_z...
    brick(ind).min_z brick(ind).max_z brick(ind).max_z],...
    material_type(ind+1).color)
patch([brick(ind).max_x brick(ind).max_x brick(ind).max_x...
    brick(ind).max_x], [brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y brick(ind).min_y], [brick(ind).min_z...
    brick(ind).min_z brick(ind).max_z brick(ind).max_z],...
    material_type(ind+1).color)
patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x], [brick(ind).min_y brick(ind).min_y...
    brick(ind).min_y brick(ind).min_y], [brick(ind).min_z...
    brick(ind).max_z brick(ind).max_z brick(ind).min_z],...
    material_type(ind+1).color)
patch([brick(ind).min_x brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x], [brick(ind).max_y brick(ind).max_y...
    brick(ind).max_y brick(ind).max_y], [brick(ind).min_z...
    brick(ind).max_z brick(ind).max_z brick(ind).min_z],...
    material_type(ind+1).color)
view([60 9])
alpha(.1)
end
[X,Y,Z]=cylinder2P(cylinders(1).radius,50,[cylinders(1).center_y cylinders(1).center_x...
cylinders(1).min_z],[cylinders(1).center_y cylinders(1).center_x cylinders(1).max_z]);
surf(X,Y,Z,'edgecolor', 'none','FaceColor',[0,0,0]);
alpha(.1)
caxis([-10000 10000])
load MyColormap
colormap(mymap)
colorbar('parent',h_1)
drawnow
v = VideoWriter('3D.avi');
open(v);

% h_3=figure(2)
% h_4=plot(zcoor,Ex_sample_2D(:));
% xlabel('z','fontsize', 12);ylabel('Ex','fontsize', 12);
% %ylim([-5.*10^2 5.*10^2]);
% grid on
% 
% h_5=figure(3)
% h_6=plot(zcoor(1:nz),Ez_sample_2D(:));
% xlabel('z','fontsize', 12);ylabel('Ez','fontsize', 12);
% %ylim([-5.*10^2 5.*10^2]);
% grid on
cpml_b_mx_xn=cpml_b_mx_xn';
cpml_a_mx_xn=cpml_a_mx_xn';
cpml_b_mx_xp=cpml_b_mx_xp';
cpml_a_mx_xp=cpml_a_mx_xp';
cpml_b_ex_xn=cpml_b_ex_xn';
cpml_a_ex_xn=cpml_a_ex_xn';
cpml_b_ex_xp=cpml_b_ex_xp';
cpml_a_ex_xp=cpml_a_ex_xp';

cpml_b_mz_zn=permute(cpml_b_mz_zn,[1 3 2]);
cpml_a_mz_zn=permute(cpml_a_mz_zn,[1 3 2]);
cpml_b_mz_zp=permute(cpml_b_mz_zp,[1 3 2]);
cpml_a_mz_zp=permute(cpml_a_mz_zp,[1 3 2]);

cpml_b_ez_zn=permute(cpml_b_ez_zn,[1 3 2]);
cpml_a_ez_zn=permute(cpml_a_ez_zn,[1 3 2]);
cpml_b_ez_zp=permute(cpml_b_ez_zp,[1 3 2]);
cpml_a_ez_zp=permute(cpml_a_ez_zp,[1 3 2]);
% spmd
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

Psi_hxy_yn(:,:,:)=cpml_b_my_yn.*Psi_hxy_yn+cpml_a_my_yn.*...
    (Ez(:,2:n_cpml_yn+1,:) - Ez(:,1:n_cpml_yn,:));
Psi_hzy_yn(:,:,:)=cpml_b_my_yn.*Psi_hzy_yn+cpml_a_my_yn.*...
    (Ex(:,2:n_cpml_yn+1,:) - Ex(:,1:n_cpml_yn,:));
Psi_hxy_yp(:,:,:)=cpml_b_my_yp.*Psi_hxy_yp(:,:,:)+cpml_a_my_yp.*...
    (Ez(:,ny-n_cpml_yp+2:ny+1,:)-Ez(:,ny-n_cpml_yp+1:ny,:));
Psi_hzy_yp(:,:,:)=cpml_b_my_yp.*Psi_hzy_yp(:,:,:)+cpml_a_my_yp.*...
    (Ex(:,ny-n_cpml_yp+2:ny+1,:)-Ex(:,1+ny-n_cpml_yp:ny,:));

Psi_hxz_zn(:,:,:)=(cpml_b_mz_zn).*Psi_hxz_zn(:,:,:)+(cpml_a_mz_zn).*...
    (Ey(:,:,2:n_cpml_zn+1) - Ey(:,:,1:n_cpml_zn));
Psi_hyz_zn(:,:,:)=(cpml_b_mz_zn).*Psi_hyz_zn(:,:,:)+(cpml_a_mz_zn).*...
    (Ex(:,:,2:n_cpml_zn+1) - Ex(:,:,1:n_cpml_zn));
Psi_hyz_zp (:,:,:) = (cpml_b_mz_zp).*Psi_hyz_zp(:,:,:)+(cpml_a_mz_zp).*...
    (Ex(:,:,nz-n_cpml_zp+2:nz+1)-Ex(:,:,nz-n_cpml_zp+1:nz));
Psi_hxz_zp (:,:,:) = (cpml_b_mz_zp).*Psi_hxz_zp(:,:,:)+(cpml_a_mz_zp).*...
    (Ey(:,:,nz-n_cpml_zp+2:nz+1)-Ey(:,:,nz-n_cpml_zp+1:nz));

Hy(1:n_cpml_xn,:,:)=Hy(1:n_cpml_xn,:,:)+CPsi_hyx_xn(:,:,:).*Psi_hyx_xn(:,:,:);
Hz(1:n_cpml_xn,:,:)=Hz(1:n_cpml_xn,:,:)+CPsi_hzx_xn(:,:,:).*Psi_hzx_xn(:,:,:);

Hy(nx-n_cpml_xp+1:nx,:,:)=Hy(nx-n_cpml_xp+1:nx,:,:)+CPsi_hyx_xp(:,:,:).*Psi_hyx_xp(:,:,:);
Hz(nx-n_cpml_xp+1:nx,:,:)=Hz(nx-n_cpml_xp+1:nx,:,:)+CPsi_hzx_xp(:,:,:).*Psi_hzx_xp(:,:,:);

Hx(:,1:n_cpml_yn,:)=Hx(:,1:n_cpml_yn,:)+CPsi_hxy_yn(:,:,:).*Psi_hxy_yn(:,:,:);
Hz(:,1:n_cpml_yn,:)=Hz(:,1:n_cpml_yn,:)+CPsi_hzy_yn(:,:,:).*Psi_hzy_yn(:,:,:);

Hx(:,ny-n_cpml_yp+1:ny,:)=Hx(:,ny-n_cpml_yp+1:ny,:)+CPsi_hxy_yp(:,:,:).*Psi_hxy_yp(:,:,:);
Hz(:,ny-n_cpml_yp+1:ny,:)=Hz(:,ny-n_cpml_yp+1:ny,:)+CPsi_hzy_yp(:,:,:).*Psi_hzy_yp(:,:,:);

Hx(:,:,1:n_cpml_zn)=Hx(:,:,1:n_cpml_zn)+CPsi_hxz_zn(:,:,:).*Psi_hxz_zn(:,:,:);
Hy(:,:,1:n_cpml_zn)=Hy(:,:,1:n_cpml_zn)+CPsi_hyz_zn(:,:,:).*Psi_hyz_zn(:,:,:);

Hx(:,:,nz-n_cpml_zp+1:nz)=Hx(:,:,nz-n_cpml_zp+1:nz)+CPsi_hxz_zp(:,:,:).* Psi_hxz_zp(:,:,:);
Hy(:,:,nz-n_cpml_zp+1:nz)=Hy(:,:,nz-n_cpml_zp+1:nz)+CPsi_hyz_zp(:,:,:).* Psi_hyz_zp(:,:,:);
%% update_magnetic_fields

Hx = ((Chxh.* Hx)+(Chxey .*(Ey (1:nxp1,1:ny,2:nzp1)- Ey (1:nxp1,1:ny,1:nz)))+...
    (Chxez.*(Ez(1:nxp1,2:nyp1,1:nz)- Ez (1:nxp1,1:ny,1:nz))));
Hy = (Chyh.* Hy)+(Chyez .*(Ez (2:nxp1,1:nyp1,1:nz)- Ez (1:nx,1:nyp1,1:nz)))+...
      (Chyex.*(Ex(1:nx,1:nyp1,2:nzp1)- Ex (1:nx,1:nyp1,1:nz)));
Hz = Chzh.* Hz+Chzex .*(Ex (1:nx,2:nyp1,1:nzp1)- Ex (1:nx,1:ny,1:nzp1))+...
     Chzey.*(Ey(2:nxp1,1:ny,1:nzp1)- Ey (1:nx,1:ny,1:nzp1));
%% update_electric_fields_for_CPML
% 
Psi_eyx_xn(:,:,:)= cpml_b_ex_xn.* Psi_eyx_xn(:,:,:)+ cpml_a_ex_xn.*...
    (Hz(2:n_cpml_xn+1,:,:)-Hz(1:n_cpml_xn,:,:));
Psi_ezx_xn(:,:,:)= cpml_b_ex_xn.* Psi_ezx_xn(:,:,:)+ cpml_a_ex_xn.*...
    (Hy(2:n_cpml_xn+1,:,:)-Hy(1:n_cpml_xn,:,:));

Ey (2:n_cpml_xn+1,:,:) = Ey (2:n_cpml_xn+1,:,:)+CPsi_eyx_xn.*Psi_eyx_xn;
Ez (2:n_cpml_xn+1,:,:) = Ez (2:n_cpml_xn+1,:,:)+CPsi_ezx_xn.*Psi_ezx_xn; 

Psi_eyx_xp(:,:,:) = cpml_b_ex_xp.*Psi_eyx_xp(:,:,:)+cpml_a_ex_xp.*...
    (Hz(1+nx-n_cpml_xp:nx ,:,:) - Hz(nx-n_cpml_xp:nx-1,:,:));
Psi_ezx_xp(:,:,:) = cpml_b_ex_xp.*Psi_ezx_xp(:,:,:)+cpml_a_ex_xp.*...
    (Hy(1+nx-n_cpml_xp:nx,:,:) - Hy(nx-n_cpml_xp:nx-1,:,:));

Ey(nx-n_cpml_xp+1:nx,:,:)=Ey(nx-n_cpml_xp+1:nx,:,:)+CPsi_eyx_xp.*Psi_eyx_xp;
Ez(nx-n_cpml_xp+1:nx,:,:)=Ez(nx-n_cpml_xp+1:nx,:,:)+CPsi_ezx_xp.*Psi_ezx_xp;

Psi_exy_yn(:,:,:) = cpml_b_ey_yn.* Psi_exy_yn(:,:,:)+ cpml_a_ey_yn.*...
    (Hz(:,2:n_cpml_yn+1,:) - Hz(:,1:n_cpml_yn,:));
Psi_ezy_yn(:,:,:) = cpml_b_ey_yn.* Psi_ezy_yn(:,:,:)+ cpml_a_ey_yn.*...
    (Hx(:,2:n_cpml_yn+1,:) - Hx(:,1:n_cpml_yn,:));

Ex (:,2:n_cpml_yn+1,:)=Ex(:,2:n_cpml_yn +1,:)+ CPsi_exy_yn .* Psi_exy_yn;
Ez (:,2:n_cpml_yn+1,:)=Ez(:,2:n_cpml_yn +1,:)+ CPsi_ezy_yn .* Psi_ezy_yn; 

Psi_exy_yp(:,:,:) = cpml_b_ey_yp.*Psi_exy_yp(:,:,:)+cpml_a_ey_yp.*...
    (Hz(:,1+ny-n_cpml_yp:ny,:) - Hz(:,ny-n_cpml_yp:ny-1,:));
Psi_ezy_yp(:,:,:) = cpml_b_ey_yp.*Psi_ezy_yp(:,:,:)+cpml_a_ey_yp.*...
    (Hx(:,1+ny-n_cpml_yp:ny,:) - Hx(:,ny-n_cpml_yp:ny-1,:));

Ex(:,ny-n_cpml_yp+1:ny,:)=Ex(:,ny-n_cpml_yp+1:ny,:)+CPsi_exy_yp.*Psi_exy_yp;
Ez(:,ny-n_cpml_yp+1:ny,:)=Ez(:,ny-n_cpml_yp+1:ny,:)+CPsi_ezy_yp.*Psi_ezy_yp;

Psi_exz_zn(:,:,:) = cpml_b_ez_zn.* Psi_exz_zn(:,:,:)+ cpml_a_ez_zn.*...
    (Hy(:,:,2:n_cpml_zn+1) - Hy(:,:,1:n_cpml_zn));
Psi_eyz_zn(:,:,:) = cpml_b_ez_zn.* Psi_eyz_zn(:,:,:)+ cpml_a_ez_zn.*...
    (Hx(:,:,2:n_cpml_zn+1) - Hx(:,:,1:n_cpml_zn));

Ex (:,:,2:n_cpml_zn+1) = Ex(:,:,2:n_cpml_zn+1)+ CPsi_exz_zn .* Psi_exz_zn;
Ey (:,:,2:n_cpml_zn +1)= Ey(:,:,2:n_cpml_zn+1)+ CPsi_eyz_zn .* Psi_eyz_zn; 

Psi_exz_zp(:,:,:) = cpml_b_ez_zp .*Psi_exz_zp(:,:,:)+cpml_a_ez_zp.*...
    (Hy(:,:,1+nz-n_cpml_zp:nz) - Hy(:,:,nz-n_cpml_zp:nz-1));
Psi_eyz_zp(:,:,:) = cpml_b_ez_zp.*Psi_eyz_zp(:,:,:)+cpml_a_ez_zp.*...
    (Hx(:,:,1+nz-n_cpml_zp:nz) - Hx(:,:,nz-n_cpml_zp:nz-1));

Ex(:,:,nz-n_cpml_zp+1:nz)=Ex(:,:,nz-n_cpml_zp+1:nz)+CPsi_exz_zp.*Psi_exz_zp;
Ey(:,:,nz-n_cpml_zp+1:nz)=Ey(:,:,nz-n_cpml_zp+1:nz)+CPsi_eyz_zp.*Psi_eyz_zp;
%% update_electric_fields

 Ex(1:nx,2:ny,2:nz)=Cexe(1:nx,2:ny,2:nz).*Ex(1:nx,2:ny,2:nz)+...
     Cexhz(1:nx,2:ny,2:nz).*(Hz(1:nx,2:ny,2:nz)-Hz(1:nx,1:ny-1,2:nz))+...
     Cexhy(1:nx,2:ny,2:nz).*(Hy(1:nx,2:ny,2:nz)-Hy(1:nx,2:ny,1:nz-1));
 Ey(2:nx,1:ny,2:nz)=Ceye(2:nx,1:ny,2:nz).*Ey(2:nx,1:ny,2:nz)+...
     Ceyhx(2:nx,1:ny,2:nz).*(Hx(2:nx,1:ny,2:nz)-Hx(2:nx,1:ny,1:nz-1))+...
     Ceyhz(2:nx,1:ny,2:nz).*(Hz(2:nx,1:ny,2:nz)-Hz(1:nx-1,1:ny,2:nz));
 Ez (2:nx,2:ny,1:nz)=Ceze(2:nx,2:ny,1:nz).*Ez(2:nx,2:ny,1:nz)+...
     Cezhy(2:nx,2:ny,1:nz).*(Hy(2:nx,2:ny,1:nz)-Hy(1:nx-1,2:ny,1:nz))+...
     Cezhx(2:nx,2:ny,1:nz).*(Hx(2:nx,2:ny,1:nz)-Hx(2:nx,1:ny-1,1:nz));
%% update the additional debye terms
for ind =1:number_of_debye
C_jj= debye(ind).C_jj;
C_je= debye(ind).C_je; 
J_x = debye(ind).J_x;
debye(ind). J_x=C_je.*(Ex(debye(ind).Ex_indices)-Exp(debye(ind).Ex_indices));
J_y = debye(ind).J_y;
debye(ind). J_y=C_jj .* debye(ind).J_y+C_je.*(Ey(debye(ind).Ey_indices)-Eyp(debye(ind).Ey_indices));
J_z = debye(ind).J_z;
debye(ind). J_z=C_jj .* debye(ind).J_z+C_je.*(Ez(debye(ind).Ez_indices)-Ezp(debye(ind).Ez_indices));
end
%add additional debye terms to electric field components
for ind=1:number_of_debye
    indices=debye(ind).Ex_indices;
    Ex(indices)=Ex(indices)-(Cexhz(indices).*dy).*(2*material_type(debye(ind).material).debye_tau./(2*material_type(debye(ind).material).debye_tau+dt)).*debye(ind). J_x;
    indices=debye(ind).Ey_indices;
    Ey(indices)=Ey(indices)-(Ceyhx(indices).*dz).*(2*material_type(debye(ind).material).debye_tau./(2*material_type(debye(ind).material).debye_tau+dt)).*debye(ind). J_y;
    indices=debye(ind).Ez_indices;
    Ez(indices)=Ez(indices)-(Cezhy(indices).*dx).*(2*material_type(debye(ind).material).debye_tau./(2*material_type(debye(ind).material).debye_tau+dt)).*debye(ind). J_z;
end
Exp=Ex;
Eyp=Ey;
Ezp=Ez;
%% update source
  Ex((nx+1)/2,(ny+1)/2,n_cpml_zn+2)=Ex((nx+1)/2,(ny+1)/2,n_cpml_zn+2)+((a_L_Ex.*sin(number*dt*omega))); % Sine source 
  Ex(((nx+1)/2)-1,(ny+1)/2,n_cpml_zn+2)=Ex(((nx+1)/2)-1,(ny+1)/2,n_cpml_zn+2)-(a_L_Ex.*sin(number*dt*omega)); % Sine source 
  Ey((nx+1)/2,(ny+1)/2,n_cpml_zn+2)=Ey((nx+1)/2,(ny+1)/2,n_cpml_zn+2)+(a_L_Ey.*sin(number*dt*omega)); % Sine source 
  Ey(((nx+1)/2),((ny+1)/2)-1,n_cpml_zn+2)=Ey(((nx+1)/2),((ny+1)/2)-1,n_cpml_zn+2)-(a_L_Ey.*sin(number*dt*omega)); % Sine source 
%end
 %% capture_and_display_sampled_fields
if mod(number,10)==0||number==totalTimeStep
Hx_sample=Hx;
Hy_sample=Hy;
Hz_sample=Hz;
Ex_sample=Ex;
Ey_sample=Ey;
Ez_sample=Ez;

[x,y,z]=meshgrid(ycoor,xcoor(1:nx),zcoor);
h_2=slice(x,y,z,Ex_sample,0,0,[],'cubic','parent',h_1);
xlabel('x','fontsize', 12);ylabel('y','fontsize', 12);zlabel('z','fontsize', 12);
title({'Ez'}); 
set(h_2, 'edgecolor', 'none')
load MyColormap
colormap(mymap)
colorbar('parent',h_1)
caxis([-900 900])
drawnow
% shading interp;
frame = getframe(h_1);
writeVideo(v,frame);

% for i=1:nzp1
% Ex_sample_2D(i)=Ex((nx+1)/2,(ny+1)/2,i) ;
% end
% set(h_4,'YData',Ex_sample_2D(:))
% 
% for i=1:nz
% Ez_sample_2D(i)=Ez((nx+1)/2,(ny+1)/2,i) ;
% end
% set(h_6,'YData',Ez_sample_2D(:))

% disp(number)
% savefig(h_1,'three_layer_dispersive_skin_fdtd_cpml_sinesoidal_3D_3D.fig')
% savefig(h_3,'three_layer_dispersive_skin_fdtd_cpml_sinesoidal_3D_2D_Ex.fig')
% savefig(h_5,'three_layer_dispersive_skin_fdtd_cpml_sinesoidal_3D_2D_Ez.fig')

     end
    
 end
toc
%delete(v);
%profile report