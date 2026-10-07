clear; clc; close all;

% --- Dados (colunas: A=y0, B=h, F=t_medio, H=sigma_tot_t) ---
M  = readmatrix('dados.xlsx','Range','A2:I6');
h  = M(:,2);   t = M(:,6);   st = M(:,8);

sh = 0.005 * ones(size(h));       % <-- incerteza de h (m): AJUSTE AQUI

% --- Linearização ---
x = h;   sx = sh;
y = t.^2;  sy = 2*t.*st;

% --- Ajuste de York ---
wx = 1./sx.^2;  wy = 1./sy.^2;
p = polyfit(x,y,1);  A = p(1);          % chute inicial
for it = 1:1000
    W    = wx.*wy ./ (wx + A^2*wy); 
    xbar = sum(W.*x)/sum(W);  ybar = sum(W.*y)/sum(W);
    U = x - xbar;  V = y - ybar;
    beta = W .* (U./wy + A*V./wx);       % (corrigido)
    Anew = sum(W.*beta.*V) / sum(W.*beta.*U);
    if abs(Anew-A) < 1e-15, A = Anew; break; end
    A = Anew;
end
B = ybar - A*xbar;

xi = xbar + beta;  xa = sum(W.*xi)/sum(W);  u = xi - xa;
sA = sqrt(1/sum(W.*u.^2));
sB = sqrt(1/sum(W) + xa^2*sA^2);
chi2red = sum(W.*(y - A*x - B).^2)/(numel(x)-2);

fprintf('A = %.6f +/- %.6f\n', A, sA);
fprintf('B = %.6f +/- %.6f\n', B, sB);
fprintf('chi2 reduzido = %.4f\n', chi2red);

% --- Gráfico ---
figure; hold on; box on; grid on;
errorbar(x, y, sy, sy, sx, sx, 'ko', 'MarkerFaceColor','k');
xf = linspace(min(x)-0.1, max(x)+0.1, 200);
plot(xf, A*xf + B, 'r-', 'LineWidth',1.5);
xlabel('h (m)'); ylabel('t^2 (s^2)');
legend('Dados','Ajuste de York','Location','northwest');