% housekeeping
clear all
close all
clc

% Importing data test

a = 'Sim1a.txt';
aID = fopen(a,'r');
testval = fscanf(aID,'%s %s %s',[1,3]);
A = fscanf(aID,'%f %f %f',[3,inf]);

g = 'Sim1g.txt';
gID = fopen(g,'r');
testval = fscanf(gID,'%s %s',[1,2]);
G = fscanf(gID,'%f %f',[2,inf]);
TA = G(2,:)'