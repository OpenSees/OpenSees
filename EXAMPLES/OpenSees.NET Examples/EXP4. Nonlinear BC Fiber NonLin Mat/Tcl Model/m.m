clc; clear all
load data\node2disp.out
load data\node1reac.out
load sapfd.txt
plot(abs(node2disp(:,3)),abs(node1reac(:,3))); hold on
plot(abs(sapfd(:,2)),abs(sapfd(:,3)),':r');