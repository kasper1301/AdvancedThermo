path(genpath('~/.matlab'),path)

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
legend({'Dew points','Bubble points'})

xlabel('$V\quad\left[\SI{}{\cubic\meter}\right]$')
ylabel('$T\quad\left[\SI{}{\kelvin}\right]$')
matlab2tikz('../fig/TVexperimental.tex','parseStrings',false)

exit