%% define constants for plant
samplingRate = 40;
dt = 1/samplingRate; % sensor/paramEst sampling rate

ni = 3;
no = 3;
u_max = [80; 80; 80];    % mbar
u_o = [0; 0; 0]; % mbar
y_max = 500*ones(2*no, 1);  % px
% y_max = [628; 455; 433];  % px
y_o = zeros(2*no, 1);
% y_o = y_max/2;    % px (for sysID)
% y_o = [y_max(1)/2; y_max(2) + y_max(1)/4; y_max(3) + y_max(1)/4];    % px (for control)

densityScale = 15e6;

ut = symunit;
px = newUnit('px', ut.um * 33/220);

%% oil channel constants
d_oil = 100 * ut.um;    % channel diameter
A_oil = pi * d_oil^2/4;     % cross-sectional area
l_oil = 10 * ut.mm;
nu_oil = 5 * ut.cSt;
rho_oil = 0.913 * ut.g/ut.ml * densityScale;
mu_oil = nu_oil * rho_oil;
mu_oil = unitConvert(mu_oil, ut.Pa * ut.s);

R_oil = 32 * mu_oil * l_oil / d_oil^2;
% R_oil = unitConvert(R_oil, ut.mbar / ut.um / ut.s); % units of R should be mbar/(um/s)
R_oil = unitConvert(R_oil, ut.mbar / ut.px / ut.s); % convert to px for convenience
R_oil = double(separateUnits(R_oil));

L_oil = rho_oil * l_oil;
% L_oil = unitConvert(L_oil, ut.mbar / ut.um / ut.s^2); % units of L should be mbar/(um/s^2)
L_oil = unitConvert(L_oil, ut.mbar / ut.px / ut.s^2); % convert to px for convenience
L_oil = double(separateUnits(L_oil));

E = 2.97 * ut.MPa;     % Young's modulus for PDMS (usually 1.32 to 2.97 MPa)
k_oil = E*A_oil/l_oil;
K_oil = rho_oil * (1500 * ut.m/ut.s)^2; % bulk modulus estimate
C_oil = A_oil / k_oil + l_oil / K_oil;
% C_oil = unitConvert(C_oil, ut.um / ut.mbar); % units of C should be um/mbar
C_oil = unitConvert(C_oil, ut.px / ut.mbar); % convert to px for convenience
C_oil = double(separateUnits(C_oil));

%% aqueous channel constants
d_aq = 100 * ut.um;    % channel diameter
A_aq = pi * d_oil^2/4;     % cross-sectional area
l_aq = 10 * ut.mm;
nu_aq = 1 * ut.cSt;
rho_aq = 0.997 * ut.g/ut.ml * densityScale;
mu_aq = nu_aq * rho_aq;
mu_aq = unitConvert(mu_aq, ut.Pa * ut.s);

R_aq = 32 * mu_aq * l_aq / d_aq^2;
R_aq = unitConvert(R_aq, ut.mbar / ut.px / ut.s); % convert to px for convenience
R_aq = double(separateUnits(R_aq));

L_aq = rho_aq * l_aq;
L_aq = unitConvert(L_aq, ut.mbar / ut.px / ut.s^2); % convert to px for convenience
L_aq = double(separateUnits(L_aq));

C0 = 1e-5;
Rc0 = 0;
L0 = L_aq;
Rl0 = 0;
Rr0 = R_aq;
C1 = C_oil;
Rc1 = 0;
L1 = L_oil;
Rl1 = 0;
Rr1 = R_oil;
C2 = C_oil;
Rc2 = 0;
L2 = L_oil;
Rl2 = 0;
Rr2 = R_oil;
