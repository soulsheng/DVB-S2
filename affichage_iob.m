% Affichage
clc
clear all
close all

echelle = 0:0.25:3;
% QPSK
load('ibo_16APSK_0.mat');
plot(echelle,TEB_LDPC)
hold on 

% 8-PSK
load('ibo_16APSK_1.mat');
plot(echelle,TEB_LDPC,'r')
hold on

% 16-APSK
load('ibo_16APSK_2.mat');
plot(echelle,TEB_LDPC,'g')

grid on
hold off

xlabel('Es/No');
ylabel('Capacite');
legend('ibo = 0','ibo = -1','ibo = -2');
title('TEB en fonction de Es/No');