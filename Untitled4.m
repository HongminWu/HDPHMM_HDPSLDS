clear all
clc
clf

%% main
fig_hei=0.4;
fig_wei=0.8;
lef_cor_x=0.15;
lef_cor_y=0.1;

subplot(2,1,1,'position', [lef_cor_x lef_cor_y+fig_hei fig_wei fig_hei]);
plot(rand(5));
set(gca, 'XTick', []);

subplot(2,1,2,'position', [lef_cor_x lef_cor_y         fig_wei fig_hei]);
plot(rand(5));
set(gca, 'XTick', []);