function plotlist = addplotlist(plotlist, name, x, y, s)

% function plotlist = addplotlist(plotlist, name, x, y, s)
%
% voeg een grafiek toe
% invoer
% x - x-coordinaten
% y - y-coordinaten
% s - tekenstijl
% 
% zie ook : plotlist, newplotlist, doplotlist, plot

if nargin < 5, s = ''; end
item.name = name;
item.x = x;
item.y = y;
item.s = s;
n = length(plotlist);
plotlist{end+1} = item;
