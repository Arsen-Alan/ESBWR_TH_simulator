
% Reactor parameters

R = 5.88/2; %m                                       (OK)
H = 3.048; %m                                        (OK)

Nfa       = 1132;  % Total number of FA              (OK)
%Nfr       = 92;    % Total number of FR per FA       (OK)
Nfr       = 78 + 14*2/3;% Total number of FR per FA  (OK) (considering part-lenght rods
dout      = 1.026; % cm outer diameter of FR         (OK)
cladThick = 0.071; % cm                              (OK) (est)            
gapThick  = 0.01;  % cm                              (OK) (est) 
din       = dout - 2*cladThick; % cm                 (OK) (est)
dpellet   = din  - 2*gapThick;  % cm                 (OK) (est)
pitch     = 1.681; % cm pitch between fuel rods      (OK) (according to datasheet 0.1inch more than standard BWR)
% Dh = 4*0.01*(pitch^2-pi*dout^2/4)/(4*pitch+pi*dout); %Hydraulic diameter m
Dh = dout*( 4/pi() * (pitch/dout)^2 - 1)*0.01;
phi = 1;
g = 9.81;         % gravity m/s^2
 
wt   = 9570;  % kg/s primary coolant flow rate       (OK)
p    = 71.7;  % bar  reactor operating pressure      (OK)
Tin  = 276.2; % ºC   core coolant inlet temperature  (OK)
Tout = 287.7; % ºC   core coolant outlet temperature (OK)
safety = 0;
Pth = 4500; % Mw   total thermal power             (OK)
Tfl = 2867; % ºC Max temperature fuel
Tcl = 1200; % ºC Max temp clad


% Other parameters

dToExtd = 5/6; % ratio between dimesion and extrapolated dimension

lambdaf = 2.5; %w/mk (OK) (est)
lambdag = 0.6; %w/mk (OK)
lambdac = 11;  %w/mk (OK)

% Simulation parameters

Nbin = 50; % Number of axial bins
deltaH = H/Nbin;
Hbin = [0]
%Hbin = [deltaH/2];

for i = 2:Nbin
    
    Hbin(i) = Hbin(i-1) + deltaH;
    
end

%
%
%

% Main 

%
%
%

% Part one: Power distribution 
r = r1(R); 
Rp = powerRatio(r, R, H, dToExtd);
q2av = aveFlux(Pth*Rp(1), Pth*Rp(2), Nfa, Nfr, dout, H); % kw/cm2
flag1 = 1
Rq = heatRatio(Nbin , H, dToExtd); 
flag2 = 1
Pf = peakingF(r, R, dToExtd);

% Pth = 4940 max power
DeltaPth = 5; %Mw 
n= 1; 
Pth = 4500

    MCPRVec = [];
    maxTc1  = [];
    maxTc2  = [];
    maxTf1  = [];
    maxTf2  = [];
    PthVec  = [];

while and(safety == 0 , n <= 100)
    
    q2av = aveFlux(Pth*Rp(1), Pth*Rp(2), Nfa, Nfr, dout, H); % kw/cm2
 
% Part two: Flow profiling (careful with the iteration parameters)

[w,CPR] = coolFlowDet2(wt, p, Nfa, Nfr, q2av, Tin, H, dout, pitch)


% Part three

%lambdaF = thermalCondF(1400); % w/m K

    [MCPR, Tf, Tc, Tw, xAct] = solverSafety(wt, p, Nfa, Nfr, q2av, Rq, Pf, Tin, Nbin, H, pitch, dout, din, dpellet, lambdaf, lambdag, lambdac);
    [CPR, Tfav, Tcav, Twav, xActav,w,G] = solverAv(wt, p, Nfa, Nfr, q2av, Rq, Tin, Nbin, H, pitch, dout, din, dpellet, lambdaf, lambdag, lambdac);
    Lb = [boilL(p, Tin, G(1), q2av(1), Rq, Nbin, H, dout, pitch),boilL(p, Tin, G(2), q2av(2), Rq, Nbin, H, dout, pitch)];
%menosDeltaP = [presdrop(Twav(1,:), xActav(1,Nbin), Nbin, Dh, H, phi,g,w(1),Lb(1),pitch, dout,Nfa, Nfr,p),presdrop(Twav(2,:), xActav(2,Nbin), Nbin, Dh, H, phi,g,w(2),Lb(2),pitch, dout, Nfa, Nfr, p)];
    
    MCPRVec(n,1) = MCPR(1);
    MCPRVec(n,2) = MCPR(2);
    maxTc1(n)  = max(Tc(1,:));
    maxTc2(n)  = max(Tc(2,:));
    maxTf1(n)  = max(Tf(1,:));
    maxTf2(n)  = max(Tf(2,:));
    PthVec(n)  = Pth;
    
    if min(MCPR)<=1.4
        safety = 1;
        fprintf("CPR less than minimum")
    end
    if or(max(Tc(1,:))>=Tcl*0.9, max(Tc(2,:))>=Tcl*0.9)
        safety = 1;
        fprintf("Temperature cladding more than maximum")
    end
    if or(max(Tf(1,:))>=Tfl*0.9, max(Tf(2,:)>=Tfl*0.9))
        fprintf("Temperature fuel more than maximum")
        safety = 1;
    end
    
    Pth = Pth + DeltaPth
    
    n = n + 1
    
end

[menosDeltaP1,menosDeltaPBoil,menosDeltaPsinglep,U] = presdrop(Twav(1,:), xActav(1,:), Nbin, Dh, H, phi,g,w(1),Lb(1),pitch, dout,Nfa, Nfr,p);
[menosDeltaP2,menosDeltaPBoil2,menosDeltaPsinglep2,U2] = presdrop(Twav(2,:), xActav(2,:), Nbin, Dh, H, phi,g,w(2),Lb(2),pitch, dout,Nfa, Nfr,p);

Coeforif = orifice(w,menosDeltaP1,menosDeltaP2,Tin,p,pitch,dout,Nfa,Nfr);
    
figure(1);
plot(Hbin, xAct(1,:));
hold on;
plot(Hbin, xAct(2,:));
hold on
grid on;
grid minor;

figure(2);
plot(Hbin, Tf(1,:));
hold on;
plot(Hbin, Tf(2,:));
hold on
grid on;
grid minor;

figure(3);
plot(Hbin, Tc(1,:));
hold on;
plot(Hbin, Tc(2,:));
hold on
grid on;
grid minor;

figure(4);
plot(Hbin, Tw(1,:));
hold on;
plot(Hbin, Tw(2,:));
hold on
grid on;
grid minor;

MCPR

%h1 = hSingleFlow(0.0713, Tout, p, pitch, dout)

%
%

% Functions


function r = r1(R) % to get equal volumes in both radial zones
    r = R/sqrt(2); 
end

function Rp = powerRatio(r1, R, H, dToExtd) % Ratio of power in radial zones / total power
    
    syms r theta z real
    Pi = sym('pi');
    f = besselj(0,2.405*r*dToExtd/R)*cos(Pi*z*dToExtd/H);        % integral f*r*dr*dtheta*dz

    I11 = int(f*r,z,-H/2,H/2);         % z integral
    I12 = int(I11,r,0,r1);             % r-integral
    I1 = int(I12,theta,0,2*Pi);        % theta-integral
    I1d = double(I1);                  % Int zone 1
    
    I21 = int(f*r,z,-H/2,H/2);         % z integral
    I22 = int(I21,r,r1,R);             % r-integral
    I2 = int(I22,theta,0,2*Pi);        % theta-integral
    I2d = double(I2);                  % Int zone 2
    
    It1 = int(f*r,z,-H/2,H/2);         % z integral
    It2 = int(It1,r,0,R);              % r-integral
    It = int(It2,theta,0,2*Pi);        % theta-integral
    Itd = double(It);                  % Int total
    
    Rp = [I1d/Itd,I2d/Itd];

end

function q2 = aveFlux(P1, P2, Nfa, Nfr, dout, H) % Av heat flux in both radial zones assuming equal volumes 

q2 = 1000*[2*P1 / (Nfa*Nfr*(pi()*dout*H*100)), 2*P2 / (Nfa*Nfr*(pi()*dout*H*100))]; %kw/cm2

end 

function Rq = heatRatio(Nbin , H, dToExtd) % Ratio q"i(axial)/q"(average)
    
    deltaH = H/Nbin;
    hi = -H/2; 

    syms z real
    Pi = sym('pi');
    f = cos(Pi*z*dToExtd/H);   
    
    RqAux = [];
    
    for i = 1:Nbin
    
        RqAux(i) = (1/deltaH)*int(f,hi, hi+deltaH) / ((1/H)*int(f,-H/2,H/2));
        
        hi = hi+deltaH;
        
    
    end 
    
    Rq = RqAux;
    
end

function Pf = peakingF(r1, R, dToExtd) % Peaking factor

    syms r theta z real
    Pi = sym('pi');
    f = besselj(0,2.405*r*dToExtd/R);  
    
    Pf1 = besselj(0,2.405*0*dToExtd/R) / (1/(Pi*r1^2)* int(f*2*Pi*r,0,r1));
 
    Pf2 = besselj(0,2.405*r1*dToExtd/R) / (1/(Pi*R^2-Pi*r1^2)* int(f*2*Pi*r,r1,R));

    Pf = [double(Pf1), double(Pf2)];
end 

function Xdry = dryout(p, G, Pf, Lb )  % dry out quality GE-CISE correlation
    
    G  = G*1e4;   % from kg/cm2s to kg/m2s
    Gr = G/1356.23;
    
    p  = p*1e5;   % from bar to Pa
    pr = p/6894.757;
    
    Lb     = Lb*1e-2; % from cm to m
    Lbstar = Lb/0.0254; % L boiling corrected 
    
    A = 1.055 - 0.013*((pr-600)/400)^2 - 1.233*Gr +0.907*Gr^2 -0.285*Gr^3;
    B = 17.98 + 78.873*Gr - 35.464*Gr^2;
    B = B / 1.43; %Correction for lattice of 10x10 instead of 7x7
    
    Xdry = (A*Lbstar*(1.24/Pf)/(B+Lbstar));
    
end 

function Xdry2 = dryout2(p , G, dout, pitch)  % dry out quality Levitan and Lantsman correlation
    
    G  = G *1e4; %kg/cm2s to kg/m2s
    
    Dh = dout*(4/pi() * (pitch/dout)^2-1);
    Dh = Dh*10;  %cm      to mm
    
    xdry8  = (0.39 + 1.57*p/98 - 2.04*(p/98)^2+0.68*(p/98)^3)*(G/1000)^(-0.5);
    Xdry2  = xdry8*(8/Dh)^(0.15); 

end

function qCr = qcrit(p, Tin, G, Xdry, dout, pitch)

    i_in = XSteam("h_pT", p, Tin); % kJ/kg
    i_f  = XSteam("hL_p", p);      % kJ/kg
    i_g  = XSteam("hV_p", p);      % kJ/kg

    i_cr = Xdry*(i_g-i_f) + i_f;   % kJ/kg

    A_cool = pitch^2 - pi()*(dout/2)^2; % cm2

    qCr = (i_cr - i_in)*G*A_cool; % kw

end % Critical power

function Lb = boilL(p, Tin, G, q2av, Rq, Nbin, H, dout, pitch)

    i_in = XSteam("h_pT", p, Tin); % kJ/kg
    i_f  = XSteam("hL_p", p);      % kJ/kg
    i_g  = XSteam("hV_p", p);      % kJ/kg
    
    deltaH = H/Nbin * 1e2; % cm
    A_cool = pitch^2 - pi()*(dout/2)^2; % cm2
    A_rod  = deltaH * pi()*dout; % cm2
    
    i      = i_in;
    i_aux  = i_in;
    N      = 1;
    
    while and( i < i_f , N <= Nbin)
        
        i = i_aux;
        i_aux = i_aux + (q2av*Rq(N)*A_rod) / (G*A_cool);
        
        N = N + 1;
        
    end
    
    N = N - 1;
    
    if and(i ~= i_f, N < Nbin)
    
        l = (i_f - i) * (G*A_cool) / (q2av*Rq(N)*pi()*dout);
        Lb = (deltaH * N + l)*0.01; %m

    else
        fprintf("\n The coolant doesn't boil \n\n");
        Lb = -1;
    end
    
    
end % Boiling length

function LbCrit = boilLcrit(p, Tin, G1, G2, Rq, Nbin, H, dout, pitch, Pf) % To get Lb by iterations

    Lb_aux   = [0,0];     %cm
    Lb_it    = [100,100]; %cm
    marginLb = 1;    %cm
    G        = [G1,G2]; % kg/cm2s
    
    n = [1,1];
    nIter = 100;

    for i = 1:2

        while and(abs(Lb_it(i) - Lb_aux(i)) > marginLb, n(i)<nIter)
    
            Xdry = dryout(p, G(i), Pf(i), Lb_it(i))

            qCr = qcrit(p, Tin, G(i), Xdry, dout, pitch) % kw

            q2Cr = qCr / (H*1e2*pi()*dout) %kw/cm2
        
            Lb_aux(i) = Lb_it(i)
            Lb_it(i)  = boilL(p, Tin, G(i), q2Cr, Rq, Nbin, H, dout, pitch);
            
            n(i) = n(i) + 1
        end
    
    end
     
    LbCrit = Lb_it;
end

function [w,CPR] = coolFlowDet(wt, p, Nfa, Nfr, q2av, Rq, Pf, Tin, Nbin, H, dout, pitch)
    
    
    G_in   = wt/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr); % Kg/(cm^2 s) Mass flux 
    
    G      = [G_in, G_in]; % Initial mass flux vector 
    
    w1     = wt/2; % Initial mass flow radial zone 1
    w2     = wt/2; % Initial mass flow radial zone 2

    LbCrit = boilLcrit(p, Tin, G(1), G(2), Rq, Nbin, H, dout, pitch, Pf);
    Xdry   = [dryout(p, G(1), Pf(1), LbCrit(1)), dryout(p, G(2), Pf(2), LbCrit(2))];
    
    qCr    = [qcrit(p, Tin, G(1), Xdry(1), dout, pitch) , qcrit(p, Tin, G(2), Xdry(2), dout, pitch)];
    
    
    
    CPR1   = qCr(1)/(q2av(1)*pi()*dout*H*1e2);
    CPR2   = qCr(2)/(q2av(2)*pi()*dout*H*1e2);
    
    marginCPR = 0.0002; % Maximum difference accepted in CPR
    wStep     = 0.2; % kg/s (change of water flow in the iteration)
    
    i     = 0;
    nIter = 15000;
    
    while and(abs(CPR1 - CPR2) > marginCPR, i<nIter)
        
        if CPR1 < CPR2
            
            w1 = w1 + wStep;
            w2 = w2 - wStep;
            
        else
            
            w1 = w1 - wStep;
            w2 = w2 + wStep;
        
        end
        
        G = [2*w1/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr), 2*w2/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr)];
        
        LbCrit = boilLcrit(p, Tin, G(1), G(2), Rq, Nbin, H, dout, pitch, Pf);
        Xdry   = [dryout(p, G(1), Pf(1), LbCrit(1)), dryout(p, G(2), Pf(2), LbCrit(2))];
        qCr    = [qcrit(p, Tin, G(1), Xdry(1), dout, pitch),qcrit(p, Tin, G(2), Xdry(2), dout, pitch)];
        
        CPR1   = qCr(1)/(q2av(1)*pi()*dout*H*1e2)
        CPR2   = qCr(2)/(q2av(2)*pi()*dout*H*1e2)
        
        i = i + 1;
        
    end
    
    if and(i == nIter, abs(CPR1 - CPR2) > marginCPR)
        
        fprintf("Warning: The maximum number of iterations was reached, w did not converge")
    end
   
    w   = [w1,w2];
    CPR = [CPR1,CPR2];

end

function [w,CPR] = coolFlowDet2(wt, p, Nfa, Nfr, q2av, Tin, H, dout, pitch)
    
    
    G_in   = wt/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr); % Kg/(cm^2 s) Mass flux 
    
    G      = [G_in, G_in]; % Initial mass flux vector 
    
    w1     = wt/2; % Initial mass flow radial zone 1
    w2     = wt/2; % Initial mass flow radial zone 2

    Xdry   = [dryout2(p , G(1), dout, pitch), dryout2(p , G(2), dout, pitch)];
    qCr    = [qcrit(p, Tin, G(1), Xdry(1), dout, pitch) , qcrit(p, Tin, G(2), Xdry(2), dout, pitch)];
    
    CPR1   = qCr(1)/(q2av(1)*pi()*dout*H*1e2);
    CPR2   = qCr(2)/(q2av(2)*pi()*dout*H*1e2);
    
    marginCPR = 0.0002; % Maximum difference accepted in CPR
    wStep     = 0.2; % kg/s (change of water flow in the iteration)
    
    i     = 0;
    nIter = 15000;
    
    while and(abs(CPR1 - CPR2) > marginCPR, i<nIter)
        
        if CPR1 < CPR2
            
            w1 = w1 + wStep;
            w2 = w2 - wStep;
            
        else
            
            w1 = w1 - wStep;
            w2 = w2 + wStep;
        
        end
        
        G = [2*w1/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr), 2*w2/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr)];
        
        Xdry   = [dryout2(p , G(1), dout, pitch), dryout2(p , G(2), dout, pitch)];
        
        for j = 1:2
            
            if Xdry(j)>1
                Xdry(j) = 1;
            end
            
        end
        
        qCr    = [qcrit(p, Tin, G(1), Xdry(1), dout, pitch),qcrit(p, Tin, G(2), Xdry(2), dout, pitch)];
        
        CPR1   = qCr(1)/(q2av(1)*pi()*dout*H*1e2);
        CPR2   = qCr(2)/(q2av(2)*pi()*dout*H*1e2);
        
        i = i + 1;
        
    end
    
    if and(i == nIter, abs(CPR1 - CPR2) > marginCPR)
        
        fprintf("Warning: The maximum number of iterations was reached, w did not converge")
    end
   
    w   = [w1,w2];
    CPR = [CPR1,CPR2];

end

function [Tg, Tfc] = tempinc(lambdaf, lambdag, lambdac, rfo,rgo, rco, h, qhot, Tw) %qhot: linear power density, Tw:temperature water
    
    DeltaTf = qhot/(4*pi()*lambdaf);
    DeltaTg = qhot*log(rgo/rfo)/(2*lambdag);
    DeltaTc = qhot*log(rco/rfo)/(2*lambdac);
    DeltaTl = qhot/(2*pi()*rco*h);
    
    Tc = DeltaTl + Tw;
    Tg = DeltaTc + Tc; % Temperature max cladding
    Tf = DeltaTg + Tg; 
    Tfc = DeltaTf + Tf; % Temperature max fuel
    
end 

function lambdaF = thermalCondF(Tf) % thermal conductivity UO2 (density 95%)
    
    Tf = Tf + 273.15; % ºC to K
    
    tau = Tf/1000;
    
    lambdaF = 100/(7.5408+17.692*tau+3.614*tau^2) + (6400/tau^(5/2)) * exp(-16.35/tau); % W/(m K)

end 

function h1Phase = hSingleFlow(G, Tw, p, pitch, dout) % heat transfer coeficient single phase (Markoczy correlation)

    Dh = dout*( 4/pi() * (pitch/dout)^2 - 1); % cm
    Dh = Dh * 1e-2; % cm to m
    
    G  = G * 1e4; % kg/cm2s to kg/m2s

    
    mu       = XSteam("my_pT", p, Tw); % Pa s
    Cp       = XSteam("Cp_pT", p, Tw); % kJ/(kg °C)
    Cp       = Cp * 1000;              % J/(kg ºC)
    lambda   = XSteam("tc_pT", p, Tw); % W/(m °C)
    
    Re = G*Dh/mu;
    Pr = Cp*mu/lambda;
    
    NuDB = 0.023 * Re^0.8*Pr^0.4;
    
    B =  4/pi() * (pitch/dout)^2 - 1;
    
    NuBundle = (1 + 0.91*Re^(-0.1)*Pr^0.4*(1-2*exp(-B)))*NuDB;
    
    h1Phase = NuBundle * lambda / Dh;

end 

function h2Phase = hTwoPhaseFlow(G, p, q2, pitch, dout, x)


    Dh = dout*( 4/pi() * (pitch/dout)^2 - 1);
    Dh = Dh * 1e-2; % cm to m
    
    G  = G * 1e4;        % kg/cm2s to kg/m2s
    q2 = q2 * 1e3 * 1e4; % kw/cm2  to w/m2
    
    Tsat     = XSteam("Tsat_p", p);          % saturation temperature
    rho_f    = 1/XSteam("vL_p", p);          % kg/ m3 the density function doesnt work
    rho_g    = 1/XSteam("vV_p", p);          % kg/ m3 the density function doesnt work
    mu_f     = XSteam("my_pT", p, Tsat-0.1); % Pa s
    mu_g     = XSteam("my_pT", p, Tsat+0.1); % Pa s
    Cp_f     = XSteam("CpL_p", p);           % kJ/(kg °C)
    Cp_f     = Cp_f * 1000;                  % J/(kg ºC)
    lambda_f = XSteam("tcL_p", p);           % W/(m °C)
    i_f      = XSteam("hL_p", p);            % kJ/kg
    i_g      = XSteam("hV_p", p);            % kJ/kg
    i_fg     = i_g - i_f;                    % KJ/kg
    i_fg     = i_fg * 1e3;                   % J/kg
    sigma    = XSteam("st_p", p);            % N/m
    
    Rel = G*(1-x)*Dh/mu_f;    % liquid Re
    Prl = Cp_f*mu_f/lambda_f; % liquid Pr
    
    Xtt = ((1-x)/x)^0.9*(rho_g/rho_f)^0.5*(mu_f/mu_g)^0.1;
    
    if 1/Xtt > 0.1
    
        F = 2.35*(0.213+1/Xtt)^0.736;
    
    else
        
        F = 1;
        
    end
    
    S = (1 + 2.56*1e-6*F^1.463*Rel^1.17)^(-1);
 
    
    h_mac = 0.023*lambda_f*Rel^0.8*Prl^0.4*F/Dh;
   
    deltaTaux    = 0;
    deltaT       = 10; 
    marginDeltaT = 0.1;
    nIter        = 200;
    n            = 1;
    
    h_mic_aux = 0.00122*S*(lambda_f^0.79*Cp_f^0.45*rho_f^0.49);
    h_mic_aux = h_mic_aux/(sigma^0.5*mu_f^0.29*i_fg^0.24*rho_g^0.24);
    
    while and(abs(deltaT - deltaTaux) > marginDeltaT, n<nIter)
        
        ps = XSteam("psat_T", Tsat+deltaT);
        
        deltaP = ps - p;     %bar
        deltaP = deltaP*1e5;  %Pa
        
        h_mic = h_mic_aux*deltaT^0.24*deltaP^0.75;
        
        h = h_mac + h_mic;
        
        deltaTaux = deltaT;
        deltaT    = q2/h;
        n = n + 1;
        
    end
    
    if abs(deltaT - deltaTaux) < marginDeltaT
        
        h2Phase = h;
    
    else
        
        fprintf("Warning; h2Phase didn't converge")
        h2Phase = -1;
        
    end
    
end

function [MCPR, Tf, Tc, Tw, xAct,G] = solverSafety(wt, p, Nfa, Nfr, q2av, Rq, Pf, Tin, Nbin, H, pitch, dout, din, dpellet, lambdaf, lambdag, lambdac)

    i_f      = XSteam("hL_p", p);            % kJ/kg
    i_g      = XSteam("hV_p", p);            % kJ/kg
    i_fg     = i_g - i_f;                    % KJ/kg
   
    
    [w,CPR] = coolFlowDet2(wt, p, Nfa, Nfr, q2av, Tin, H, dout, pitch);
    
    G = [2*w(1)/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr), 2*w(2)/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr)];
    
    i_cell = []; % matrix of av entalpy per cell
    x_act  = []; % matrix of av quality per cell
    T_w    = []; % matrix of av coolant temp per cell
    T_c    = []; % matrix of max cladding temp per cell
    T_f    = []; % matrix of max fuel temp per cell
    
    deltaH = H/Nbin * 1e2; % cm
    A_cool = pitch^2 - pi()*(dout/2)^2; % cm2
    A_rod  = deltaH * pi()*dout; % cm2

    rfo = (dpellet/2)*1e-2; % m
    rgo = (din/2)*1e-2;     % m
    rco = (dout/2)*1e-2;    % m

    for zone = 1:2
        
        i_aux  = 0;
        i_in   = XSteam("h_pT", p, Tin); % kJ/kg
        
        flag = 0;
        for j = 1:Nbin
            
            q2jHot = q2av(zone)*Rq(j)*Pf(zone);
            
            i_aux  = i_in;
            i_in   = i_in + (q2jHot*A_rod) / (G(zone)*A_cool);
            i_cell(zone,j) = (i_in+i_aux)/2;
            
            x_act (zone,j) = (i_cell(zone,j)-i_f)/i_fg;
            
            T_w   (zone,j) = XSteam("T_ph",p,i_cell(zone,j));
            
            qhot = q2jHot * pi() * dout; % kw/cm
            qhot = qhot * 1e3 * 1e2; % w/m
            
            if i_cell(zone,j) < i_f
                
                h = hSingleFlow(G(zone), T_w(zone,j), p, pitch, dout); % W/(m2K)
                
            else
                if flag == 0
                  
                    h1 = hSingleFlow(G(zone), T_w(zone,j)-0.1, p, pitch, dout); % W/(m2K)
                    h2 = hTwoPhaseFlow(G(zone), p, q2jHot, pitch, dout, x_act(zone,j)); % W/(m2K)
                    
                    h  = h1 * abs((i_f-i_aux)/(i_in-i_aux)) + h2 * abs((i_f-i_in)/(i_in-i_aux));
                    
                    flag = 1;
                    
                else
                    
                    h = hTwoPhaseFlow(G(zone), p, q2jHot, pitch, dout, x_act(zone,j)); % W/(m2K)
                                     
                end
                
            end
            h;
            [T_c(zone,j),T_f(zone,j)] = tempinc(lambdaf, lambdag, lambdac, rfo,rgo, rco, h, qhot, T_w(zone,j));
            
            
        end
        
    end 

     
    MCPR = [CPR(1)/Pf(1),CPR(2)/Pf(2)]; 
    Tf   = T_f;
    Tc   = T_c;
    Tw   = T_w;
    xAct = x_act;
    
end

function [CPR, Tfav, Tcav, Twav, xActav,w,G] = solverAv(wt, p, Nfa, Nfr, q2av, Rq, Tin, Nbin, H, pitch, dout, din, dpellet, lambdaf, lambdag, lambdac)

    i_f      = XSteam("hL_p", p);            % kJ/kg
    i_g      = XSteam("hV_p", p);            % kJ/kg
    i_fg     = i_g - i_f;                    % KJ/kg
   
    
    [w,CPR] = coolFlowDet2(wt, p, Nfa, Nfr, q2av, Tin, H, dout, pitch);
   
    G = [2*w(1)/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr), 2*w(2)/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr)];
    
    i_cell = []; % matrix of av entalpy per cell
    x_act  = []; % matrix of av quality per cell
    T_w    = []; % matrix of av coolant temp per cell
    T_c    = []; % matrix of max cladding temp per cell
    T_f   = []; % matrix of max fuel temp per cell
    
    deltaH = H/Nbin * 1e2; % cm
    A_cool = pitch^2 - pi()*(dout/2)^2; % cm2
    A_rod  = deltaH * pi()*dout; % cm2

    rfo = (dpellet/2)*1e-2; % m
    rgo = (din/2)*1e-2;     % m
    rco = (dout/2)*1e-2;    % m

    for zone = 1:2
        
        i_aux  = 0;
        i_in   = XSteam("h_pT", p, Tin); % kJ/kg
        
        flag = 0;
        for j = 1:Nbin
            
            q2jHot = q2av(zone)*Rq(j);
            
            i_aux  = i_in;
            i_in   = i_in + (q2jHot*A_rod) / (G(zone)*A_cool);
            i_cell(zone,j) = (i_in+i_aux)/2;
            
            x_act (zone,j) = (i_cell(zone,j)-i_f)/i_fg;
            
            T_w   (zone,j) = XSteam("T_ph",p,i_cell(zone,j));
            
            qhot = q2jHot * pi() * dout; % kw/cm
            qhot = qhot * 1e3 * 1e2; % w/m
            
            if i_cell(zone,j) < i_f
                
                h = hSingleFlow(G(zone), T_w(zone,j), p, pitch, dout); % W/(m2K)
                
            else
                if flag == 0
                  
                    h1 = hSingleFlow(G(zone), T_w(zone,j)-0.1, p, pitch, dout); % W/(m2K)
                    h2 = hTwoPhaseFlow(G(zone), p, q2jHot, pitch, dout, x_act(zone,j)); % W/(m2K)
                    
                    h  = h1 * abs((i_f-i_aux)/(i_in-i_aux)) + h2 * abs((i_f-i_in)/(i_in-i_aux));
                    
                    flag = 1;
                    
                else
                    
                    h = hTwoPhaseFlow(G(zone), p, q2jHot, pitch, dout, x_act(zone,j)); % W/(m2K)
                                     
                end
                
            end
           
            [T_c(zone,j),T_f(zone,j)] = tempinc(lambdaf, lambdag, lambdac, rfo,rgo, rco, h, qhot, T_w(zone,j));
            
            
        end
        
    end 

     
    CPR = [CPR(1),CPR(2)]; 
    Tfav   = T_f;
    Tcav   = T_c;
    Twav   = T_w;
    xActav = x_act;
    
end

function [menosDeltaP,menosDeltaPBoil,menosDeltaPsinglep,U] = presdrop(Tw, Xout, Nbin, Dh, H, phi,g,w,Lb, pitch, dout, Nfa, Nfr,p )

deltaH = H/Nbin; 
menosDeltaPsinglep = 0;   
menosDeltaPBoil = 0;
Haux= 0;
Tsat     = XSteam("Tsat_p", p);
G = w/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr); % igual el G hay que ponerle el ./ en las operaciones al ser un vector



for i = 1:Nbin
    
        if Tw(i) < Tsat
            rhoL = 1/XSteam("v_pT",p, Tw(i));
            Tw(i);
            muL = XSteam("my_pT",p, Tw(i));
            U = w/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr*rhoL);
            Cf=32*muL/(rhoL*Dh*U);
            deltaPsingle= Cf*(4*(deltaH)*G^2)/(Dh*2*rhoL) + (deltaH)*rhoL*g*sin(phi);
            menosDeltaPsinglep = menosDeltaPsinglep + deltaPsingle;
            Haux = Haux + deltaH;
       
        else
            
            rhoL    = 1/XSteam("vL_p", p);          % kg/ m3 the density function doesnt work
            rhoV    = 1/XSteam("vV_p", p);          % kg/ m3 the density function doesnt work
            muL     = XSteam("my_pT", p, Tsat-0.1); % Pa s
            muV     = XSteam("my_pT", p, Tsat+0.1); % Pa s
            
            U = w/((pitch^2-(dout/2)^2*pi())*Nfa*Nfr*rhoL);
            Cf=32*muL/(rhoL*Dh*U);
            a = Xout(i);
            r2 = (rhoL/rhoV-1)*(Xout(i)-Xout(i-1));
            syms eta 
            f=(1+Xout(i)*(rhoL/rhoV-1)*eta)/(1+Xout(i)*(muL/muV-1)*eta)^0.25;
            r3s = int(f,eta,0,1);
            r3 = double (r3s);
            r4 = rhoV*log((rhoV+(rhoL- rhoV)*Xout)/rhoV)/((rhoL-rhoV)*Xout);
            %deltaPtwo = r2*G^2/rhoL + r3*Cf*(4*deltaH*G^2)/(Dh*2*rhoL) + r4*deltaH*rhoL*g*sin(phi)
            deltaPtwo = r2*G^2/rhoL;  %phi = 1 para flujo ascendente en este caso
            deltaPtwo = deltaPtwo + r3*Cf*(4*deltaH*G^2)/(Dh*2*rhoL);
            deltaPtwo = deltaPtwo + r4*deltaH*rhoL*g*sin(phi);
            
            menosDeltaPBoil = menosDeltaPBoil + deltaPtwo;
            Haux = Haux + deltaH;
        end
end
    
    menosDeltaP = menosDeltaPsinglep + menosDeltaPBoil;


end
function Coeforif = orifice(w,menosDeltaP1, menosDeltaP2,Tin,p,pitch,dout,Nfa,Nfr)
 rhoL = 1/XSteam("v_pT",p, Tin);
 G = w./((pitch^2-(dout/2)^2*pi())*Nfa*Nfr);
 Pdiff = menosDeltaP1-menosDeltaP2;
Coeforif = Pdiff*2*rhoL/(G(1)*10000)^2 ; 
end