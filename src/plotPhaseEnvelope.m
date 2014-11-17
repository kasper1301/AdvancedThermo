clear all
close all
dewPoints    = importdata('../data/dewPoints.csv');
bubblePoints = importdata('../data/bubblePoints.csv');

% Convert to SI-units
dewP = dewPoints(:,1)*6895; % [Pa]
dewT = (dewPoints(:,2)-32)/1.8 + 273.15; % [K]
dewRho = dewPoints(:,3)*1e3; % [kg/m^3]

bubbleP = bubblePoints(:,1)*6895; % [Pa]
bubbleT = (bubblePoints(:,2)-32)/1.8 + 273.15; % [K]
bubbleRho = bubblePoints(:,3)*1e3; % [kg/m^3]

n_tot   = [1.60; 94.50; 2.60; 0.81; 0.52]/100.03; % [mol]
mM      = [28.0134; 16.0425; 30.069; 44.0956; 58.1222]*1e-3; % [kg/mol]
m_tot   = mM'*n_tot;

dewV = m_tot./dewRho;
bubbleV = m_tot./bubbleRho;

semilogx(dewV,dewT,'bo')
hold on
semilogx(bubbleV,bubbleT,'rx')
legend({'Experimental dew points','Experimental bubble points'})

xlabel('$V\quad\left[\SI{}{\cubic\meter}\right]$')
ylabel('$T\quad\left[\SI{}{\kelvin}\right]$')
matlab2tikz('../fig/TVexperimental.tex','parseStrings',false, 'width', '\textwidth')

load('results.mat')

clear dewV bubbleV
dewV = [];
dewT = [];
bubbleV = [];
bubbleT = [];

for i = 1:length(V)
    if isempty(V{i})
        continue
    end
    [Vmax, idx] = max(V{i});
    if X{1}{i}(idx)/Vmax > .99
        % Dew point
        dewV = [dewV; Vmax];
        dewT = [dewT; T{i}(idx)];
    elseif X{1}{i}(idx)/Vmax < .01
        % Bubble point
        bubbleV = [bubbleV; Vmax];
        bubbleT = [bubbleT; T{i}(idx)];
    end
    [Vmin, idx] = min(V{i});
    if X{1}{i}(idx)/Vmin > .99
        % Dew point
        dewV = [dewV; Vmin];
        dewT = [dewT; T{i}(idx)];
    elseif X{1}{i}(idx)/Vmin < .01
        % Bubble point
        bubbleV = [bubbleV; Vmin];
        bubbleT = [bubbleT; T{i}(idx)];
    end
end
[dewV, idx] = sort(dewV);
dewT = dewT(idx);
[bubbleV, idx] = sort(bubbleV);
bubbleT = bubbleT(idx);

semilogx(dewV, dewT, 'b', 'LineWidth', 2)
semilogx(bubbleV, bubbleT, 'r', 'LineWidth', 2)

legend({'Experimental dew points','Experimental bubble points',...
        'Calculated dew curve','Calculated bubble curve'},...
        'Location', 'SouthEast')
    
matlab2tikz('../fig/TV.tex','parseStrings',false, 'width', '\textwidth')
%%
close all
Tvec = zeros(length(T),1);
for i = 1:length(T)
    Tvec(i) = T{i}(1);
end
[Tvec,Tidx] = sort(Tvec);

color = jet(length(Tidx));

legendstring = {};

for j = 1:floor(length(Tidx)/10):length(Tidx)
    i = Tidx(j);
    [Vi, idx] = sort(V{i});
    Si = S{i}(idx);
    idx = 1:floor(length(Vi)/100):length(Vi);
    Vi = Vi(idx);
    Si = Si(idx);
    semilogx(Vi,Si,'Color',color(j,:),'LineWidth', 2)
    legendstring = {legendstring{:}, ['$T = \SI{', num2str(Tvec(j)), '}{\kelvin}$']};
    hold on
end
legend(legendstring, 'Location', 'EastOutside')
xlabel('$V\quad\left[\SI{}{\cubic\meter}\right]$')
ylabel('$S\eos\quad\left[\SI{}{\joule\per\kelvin\per\mole}\right]$')
axis tight
matlab2tikz('../fig/SV.tex','parseStrings', 'width', '\textwidth')

