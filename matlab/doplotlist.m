function doplotlist(plotlist, command, pos)

% function doplotlist(plotlist, command, pos)
%
% maak een grafiek met daarop alle items in plotlist
% 
% zie ook : plotlist, newplotlist, addplotlist, plot, legend

if nargin == 2
  pos = [];
end

if length(plotlist) < 1, return; end
string = '';
sep = '';
names = {};
for i = 1:length(plotlist)
  string = sprintf('%s%splotlist{%d}.x, plotlist{%d}.y, plotlist{%d}.s', ...
		   string, sep, i, i, i);
  sep = ', ';
  names{i} = plotlist{i}.name;
end
%fprintf('%s\n', string);
string = sprintf('%s(%s)', command, string);
%fprintf('%s\n', string);
eval(string);

if length(pos) ~= 0
  %pos = 4; % lower right-hand corner
  legend(names, pos);
end

