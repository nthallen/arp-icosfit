%%
% Map useful values for R1, R2 vs L based on equations
Rmax = 3;
L = 1;
x1 = [L L Rmax Rmax L];
y1 = [0 -Rmax -Rmax L-Rmax 0];
x2 = [L L 0 L];
y2 = [0 L L 0];
x3 = [L Rmax Rmax L L];
y3 = [L L Rmax Rmax L];
x4 = y1;
y4 = x1;

clf;
fill(x1,y1,'b');
hold on;
fill(x2,y2,'b');
fill(x3,y3,'b');
fill(x4,y4,'b');
hold off;
set(gca,'xtick',[0 L Rmax],'XTickLabel', {0, 'L', 'R_{max}'});
set(gca,'ytick',[-Rmax 0 L Rmax], 'YTickLabel', {'-R_{max}', 0, 'L', 'R_{max}'});
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('R_1');
ylabel('R_2');
grid;