clear all
close all 
clc
%%
silver = 1; % 1 for silver and 0 for gold
%%
h = 6.626070040e-34;	% Planck's constant [J s]
joule_to_evs =  6.242e+18;  % Joule = Watt * sec
h_bar = h/(2*pi) * joule_to_evs; % h_bar in [eV s];

if silver == 0 %'gold'
      w_p_ev = 8.89;
      Gamma_ev = 0.07088;
      w_p = w_p_ev/h_bar;% * ev_to_angular_frequency;
      Gamma = Gamma_ev/h_bar;% * ev_to_angular_frequency;
      filename = 'gold';
elseif silver == 1 %'silver'
      w_p_ev = 9.04;
      Gamma_ev = 0.02125;
      w_p = w_p_ev/h_bar;% * ev_to_angular_frequency;
      Gamma = Gamma_ev/h_bar;% * ev_to_angular_frequency;
      filename = 'silver';
end
fTHz = linspace(300,1000,3000);
c = physconst('LightSpeed');
mum = c./fTHz * 1e-6;
nm = mum * 1e3;

fHz = fTHz * 1e12;

omega = 2*pi * fHz ;

eps_rel = 1 - w_p.^2 ./ (omega.^2 + 1i * Gamma .* omega);
f = figure;
f.Position = [1012 589 541 185];
plot(fTHz, real(eps_rel),'--k','LineWidth',3)
ax = gca;

for kk = 1:length(nm)
    if silver == 0%'gold'
        eps_c(kk) = geteps_gold(nm(kk)*1e-3);
    elseif silver == 1%'silver'
        eps_c(kk) = geteps(nm(kk)*1e-3);
    end
end
hold on
plot(fTHz, real(eps_c),'-b','LineWidth',3)

% 
ax.YLim = [-60,10];
ax.XLabel.Interpreter = 'Latex';
ax.XLabel.String = 'frequency in $\mathrm{THz}$';
ax.YLabel.Interpreter = 'Latex';
ax.YLabel.String = '$\mathrm{Re}(\varepsilon_r)$';

ax.FontSize = 14;
print(gcf, '-depsc', strcat(filename,'_real'));

g = figure;
g.Position = [1012 308 541 185];

plot(fTHz, imag(eps_rel),'--k','LineWidth',3)
hold on
plot(fTHz, imag(eps_c),'-b','LineWidth',3)
ax = gca;
ax.YLim = [0,6];
ax.XLabel.Interpreter = 'Latex';
ax.XLabel.String = 'frequency in $\mathrm{THz}$';
ax.YLabel.Interpreter = 'Latex';
ax.YLabel.String = '$\mathrm{Im}(\varepsilon_r)$';
legend('Drude model','Experimental data','Location','NorthWest')
ax.FontSize = 14;
% print(gcf, '-depsc', strcat(filename,'_imag'));