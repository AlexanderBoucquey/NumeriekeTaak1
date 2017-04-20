pl = newplotlist;
x = 1:5;
y = x.^2;
pl = addplotlist(pl, 'een', x, y, 'b');
x = 1:10;
y = x + 3;
pl = addplotlist(pl, 'twee', x, y, 'g');
figure(1)
doplotlist(pl, 'plot');
figure(2)
doplotlist(pl, 'plot', 4);
