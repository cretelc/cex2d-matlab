close all; clear all; clc;

Pr_Pf = linspace(0, 0.99);
vswr = (1 + sqrt(Pr_Pf)) / (1 - sqrt(Pr_Pf));

plot(Pr_Pf, vswr)