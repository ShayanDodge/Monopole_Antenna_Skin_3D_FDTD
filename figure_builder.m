clc
clear all
close all
uiopen('C:\Users\Shayan-pc\Desktop\3D\three_layer_dispersive_skin_fdtd_cpml_sinesoidal_3D_3D.fig',1)
for i=0:26
    saveas(gcf,['three_layer_dispersive_skin_fdtd_cpml_sinesoidal_3D_3D_' num2str(27-i) '.tiff'])
    close(gcf)
end