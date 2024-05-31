close all; clear all; clc;

Pr_Pf = linspace(0, 0.25);
vswr = (1 + sqrt(Pr_Pf)) ./ (1 - sqrt(Pr_Pf));

semilogy(Pr_Pf, vswr)
grid("on")