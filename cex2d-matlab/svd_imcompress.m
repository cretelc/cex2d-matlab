% Steve Brunton - Singular Value Decomposition 

close all; clear all; clc;


A = imread('thruster.jpeg');
X = double(rgb2gray(A));
nx = size(X,1); ny = size(X,2);

subplot(2,2,1)
imagesc(A), axis off
title('Original');

set(gcf, "Position", [1400 100 1200 1600])

%% SCD Approximations
[U, S, V] = svd(X, 'econ');


plotind = 2; 
for r = [5 20 100]
    Xapprox = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';
    subplot(2, 2, plotind), plotind = plotind+1;
    imagesc(Xapprox), axis off;
    %title ('r=', num2str(r, '%d'));
end 


%% Singular Values 
figure, subplot(1,2,1)
semilogy(diag(S), 'k', 'Linewidth',2), grid on
xlabel('r')
ylabel('Singular value, \sigma_r')
xlim([-50, 500])
set(gca, 'FontSize', 14)
subplot(1, 2, 2)
plot(cumsum(diag(S))/sum(diag(S)), 'k', 'LineWidth', 2), grid on
xlabel('r')
ylabel('Cumulative Energy')
xlim([-50 500]);ylim([0 1.1])
set(gca, 'FontSize', 14)
set(gcf, 'Position', [1400 100 3*550 3*240])
%% Correlation Values
% The U and V columsn are eigenvectors 

XXt = X*X';
XtX = X'*X;
subplot(1,2,1)
imagesc(XtX);
axis off;
subplot(1,2,2)
imagesc(XXt)
axis off;