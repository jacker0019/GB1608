clear all
close all
clc

r = 0.05;
L = 0.4;
a_1 = 0.2;
a_2 = 0.2;
a_3 = 0.57;
x = 0.2112;
y = 0.2197;
z = 0.0000 - 0.3404i;
d_1 = ((128*3^(1/2))/5 - 271/25)^(1/2)/40 + 3/8; %0.5197;
d_2 = ((32*3^(1/2))/5 + 3569/25)^(1/2)/40 + 3/8; %0.6851;
d_3 = 831^(1/2)/100 - 3/25; %0.1683;

G_1 = (x - d_1 - r/2)^2 + (y + a_1 - r*sqrt(3)/2)^2 + (z)^2 - L^2;
G_2 = (x - d_2 - r/2)^2 + (y - a_2 + r*sqrt(3)/2)^2 + (z)^2 - L^2;
G_3 = (x + r - a_3 - d_3)^2 + (y)^2 + (z)^2 - L^2;
