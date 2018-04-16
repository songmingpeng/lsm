% colormap_plot 
% author: taiping.z(email:taiping.z@outlook.com)
% date: Web Apr 26 2017
% computational neuroscience lab at SIA
% colormap_plot 

clc;
clear all;
close all;

max_color_value = 10;

% jet_color = colormap(hsv(max_color_value));
% jet_color = colormap(cool(max_color_value));
 jet_color = colormap(hot(max_color_value));
% jet_color = colormap(pink(max_color_value));
% jet_color = colormap(gray(max_color_value));
% jet_color = colormap(pink(max_color_value));
% jet_color = colormap(bone(max_color_value));
% jet_color = colormap(jet(max_color_value));
% jet_color = colormap(copper(max_color_value));
% jet_color = colormap(prim(max_color_value));
% jet_color = colormap(flag(max_color_value));

x = linspace(0,6*pi,1000);

y = sin(x);

color_index = ceil(y*5 + 5);

for i = 1:1:length(y),

    selected_color = jet_color(color_index(i),:);

    plot(x(i), y(i), 'o','color',selected_color);

    hold on;
end