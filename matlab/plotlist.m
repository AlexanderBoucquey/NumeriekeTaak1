% plotlist
%
% enkele handige functies om meerdere grafieken samen te voegen
%
% newplotlist - maak een nieuwe lijst van grafieken
% addplotlist - voeg een grafiek toe
% doplotlist - geef de grafieken weer
% testplotlist - een voorbeeld
%
% voorbeeld
% pl = newplotlist;
% x = 1:5;
% y = x.^2;
% pl = addplotlist(pl, 'een', x, y, 'b');
% x = 1:10;
% y = x + 3;
% pl = addplotlist(pl, 'twee', x, y, 'g');
% doplotlist(pl, 'plot');
%
% zie ook : plot, semilogy, semilogx, loglog