function total_ldg_calculation()
clc
% Parameters regressed from Freeman's viscosity data
par = [26.1, 0.02532, 6.7524, -8.212, 3.755];  % Regressed by Di (2015 Q1 Report)
%par = [26.1592, 0.0265, 7.6949, -7.8028, 3.3663]; % Regress by Gao
writefile = 'UT Austin Data 3-26-18 master.xlsx';
readfile = writefile;
%readfile= 'UT Austin Data 3-25-18.xlsx';
sheet = 1;

lean_Tset = xlsread(readfile,'Sheet1', 'FX5:FX1443');
lean_rhoset = xlsread(readfile, 'JM5:JM1443');

mid_Tset = xlsread(readfile,'Sheet1', 'CI5:CI1443');
mid_rhoset = xlsread(readfile,'Sheet1', 'JR5:FJR1443');

rich_Tset = xlsread(readfile,'Sheet1', 'JQ5:JQ1443');
rich_rhoset= xlsread(readfile,'Sheet1', 'CS5:CS1443');
miuset =xlsread(readfile,'Sheet1', 'CT5:CT1443');

% Unit conversion
lean_Tset = (lean_Tset+ 459.67) .* 5./9;
mid_Tset = (mid_Tset+ 459.67) .* 5./9;
rich_Tset = (rich_Tset+ 459.67) .* 5./9;
lean_rhoset = lean_rhoset ./0.062427960576145;
mid_rhoset = mid_rhoset ./0.062427960576145;
rich_rhoset = rich_rhoset./0.062427960576145;

Cpzset = zeros(length(lean_Tset), 1);
Cco2set = zeros(length(lean_Tset), 1);
lldgset = zeros(length(lean_Tset), 1);
mldgset = zeros(length(lean_Tset),1);
rldgset = zeros(length(lean_Tset), 1);
mset = zeros(length(lean_Tset), 1);
exitflags = zeros(length(lean_Tset), 1);

for i =1:length(lean_Tset)
    T_lean = lean_Tset(i);
    T_mid = mid_Tset(i);
    T_rich = rich_Tset(i);
    rho_lean = lean_rhoset(i);
    rho_mid = mid_rhoset(i);
    rho_rich = rich_rhoset(i);
    miu = miuset(i);
    [Cpz, Cco2, rldg, exitflag] = ldg_solver(T_rich, rho_rich, miu);
    Cpzset(i) = Cpz;
    Cco2set(i) = Cco2;
    rldgset(i) = rldg;
    m = 1000/(1000/Cpz - 86.136 - 88*rldg);
    mset(i) = m;   
    lldg = ldg_calc(T_lean, m, rho_lean);
    lldgset(i) = lldg;
    mldg = ldg_calc(T_mid, m, rho_mid);
    mldgset(i) = mldg;
    exitflags(i) =exitflag;
    
end

x = [mset,lldgset, rldgset, mldgset];
header = {'calc molarity', 'calc lldg', 'calc rldg', 'calc mldg'};
xlswrite(writefile, header, sheet, 'JX2')
xlswrite(writefile, x, sheet, 'JX5')

    function [CPZ, CCO2, rldg, exitflag] = ldg_solver(Tk, rho_kg, miu)
% viscosity and density at given temperature
v = @(x)(viscosity(x, Tk)-miu);
d = @(x)(density(x, Tk)-rho_kg);

guess = [3.21, 1.78];

fun1 = @(x) [v(x), d(x)];
options = optimoptions('fsolve','Display','off');
[x, fval, exitflag, output] = fsolve(fun1, guess,options);
CPZ = x(1);      % PZ concentration in mol/kg
CCO2 = x(2);    % CO2 conc in mol/kg
rldg = CCO2/(2*CPZ);

    end

    function vis = viscosity(C, T) 
    %function to calculate solvent viscosity on C&T
        CPZ = C(1);
        CCO2 = C(2);
        v_water = 2.414e-5*10.^(247.8./(T-140))*1000;  % water viscosity
        A = par(1);
        B = par(2);
        C = par(3);
        D = par(4);
        E = par(5);
        vis = v_water * exp((A/T - B)*(C*CPZ + D*CCO2 + E*CPZ.*CCO2));
    end

    function den = density (C, T)
    % function to calculate solvent density on C&T
        CPZ = C(1);
        CCO2 = C(2);
        Tc = T - 273.15;
        rho_water = 1000*(1-(Tc+288.9414)/(508929.2*(Tc+68.12963))*(Tc-3.9863)^2);
        rho_solvent = rho_water * (0.0408 * CCO2 + 0.008 * CPZ + 0.991);
        den = rho_solvent;
    end

    function y = lldg_density(x, T, M)
    % function to calculate solvent density at given loading, temperature,
    % and molarity
    pzwt = (M*86.136)/((44*x *2*M)+1000+M*86.136);
    CPZ = pzwt*1000/86.136;
    CCO2 = 2*CPZ*x;

    Tc = T - 273.15;
    rho_water = 1000*(1-(Tc+288.9414)/(508929.2*(Tc+68.12963))*(Tc-3.9863)^2);
    rho_solvent = rho_water * (0.0408 * CCO2 + 0.008 * CPZ + 0.991);
    y = rho_solvent;

    end

    function ldg = ldg_calc( T, M, rho)
    % function to solve solvent loading on given temperature, molarity, and
    % density.
        guess = 0.21;
        fun2 = @(x) lldg_density(x, T, M)- rho;
        options = optimoptions('fsolve','Display','off');
        ldg = fsolve(fun2, guess,options);
    end
end
