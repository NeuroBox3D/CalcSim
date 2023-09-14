function [v,o] = jvdcc(cc,os,t0,dt,V)
% cc is the current cytosolic calcium 
% os is a vector of the current states
% t0 is the current time
% dt is the time step size
% V this is the vector of voltages

num_steps=25;      % number of time steps for ode solve
k = dt/num_steps;   % time step size
co=1.0;               % concentration of extracellular calcium
rhovdcc = 1.0e12;     % density of vdcc channels

% compute the voltge using GHK equation, using
% voltage V, the current cytosolic calcium and extracellular calcium

F = GHK(V,cc,co);

% get the current states, assign to variable
s1 = os(:,1); s2 = os(:,2);
ts1_0 = 1.7e-3; ts2_0 = 70.0e-3;

% I am using FE for vdcc ode equations
% fprintf('\n###################################\n')
% fprintf("\t Solving VDCC ODE equations...\n")
% fprintf("\t Start t0 = %0.7f [s]\n",t0)
% fprintf("\t\t Iteration..")
[a1,b1]=rates_s1(V);
[a2,b2]=rates_s2(V);
for i=1:num_steps
%     fprintf("%i..",i)    
%     if ((mod(i,20)==0) && (i<num_steps))
%         fprintf("\n\t\t Iteration..")
%     end
    F1 = f(s1,a1,b1,ts1_0); F2 = f(s2,a2,b2,ts2_0);
    s1 = s1 + k.*F1; s2 = s2 + k.*F2;
end
% fprintf("\n\t End tf = %0.7f [s]\n",(num_steps*k+t0))

G = s1.*(s2).^2;
o = [s1, s2];

v = rhovdcc.*G.*F;
end

function F = GHK(V,cc,co)
pca = 3.8e-20;    % permeability of calcium in vdcc channels N-type
Rgas = 8.314;     % universal gas constant Joules/(K.mol)
Far = 96485.0;    % Faraday constant Coulombs/mol
T = 310.0;          % Temperataure in Kelvin
z = 2.0;            % valence

c1 = z.*Far.*V./(Rgas*T);
F =-1.0.*pca.*c1.*(cc-co.*exp(-c1))./(1.0-exp(-c1));

ind = ((abs(V)<=1e-8) | (V==0));
F(ind) = pca.* ((co-cc(ind))-(Far./(Rgas.*T)).*(co+cc(ind)).*V(ind));

end

%% right hand side of ode equations
function v = f(s,a,b,tau0)

xinf = a./(a+b);
tau_s = 1./(a+b)+tau0;

v = (xinf-s)./tau_s;
end

%% rate function on first state
function [a,b]=rates_s1(V)
Rgas = 8.314;   % universal gas constant Joules/(K.mol)
Far = 96485.0;    % Faraday constant Coulombs/mol
T = 310.0;        % Temperataure in Kelvin

K = 1.7e-3; gm = 0; Vs = -21.0e-3; z=2.0;

a = K.*exp(z.*gm.*(V-Vs).*Far./(Rgas.*T));
b = K.*exp(-z.*(1-gm).*(V-Vs).*Far./(Rgas.*T));
end

%% rate function on second state
function [a,b]=rates_s2(V)
Rgas = 8.314;   % universal gas constant Joules/(K.mol)
Far = 96485.0;    % Faraday constant Coulombs/mol
T = 310.0;        % Temperataure in Kelvin

K = 70.0e-3; gm = 0.0; Vs = -40.0e-3; z=1.0;

a = K.*exp(z.*gm.*(V-Vs).*Far./(Rgas.*T));
b = K.*exp(-z.*(1-gm).*(V-Vs).*Far./(Rgas.*T));
end