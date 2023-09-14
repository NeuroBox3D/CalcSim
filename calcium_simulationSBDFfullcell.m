function calcium_simulationSBDFfullcell(filename,outfolder,datafolder,load_syn,vtk)
%% This is the default input morphology
addpath(genpath([pwd, filesep, 'codes' ]));
if (nargin == 0)
    filename = 'morphos/NMO_150605/ref2smallsoma.swc';
    outfolder = 'vtk0_fullcellcalb40_5Hz';
    datafolder = 'results/fullcellcalb40_5Hz';
    load_syn = 0;
    vtk = 0;
end

%make output folders
if ((~isfolder(outfolder)) && (vtk ~=0))
    %fprintf(sprintf('Making %s output folder\n',outfolder))
    mkdir(outfolder)
end

if ~isfolder(datafolder)
    %fprintf(sprintf('Making %s output folder\n',datafolder))
    mkdir(datafolder)
end

% First get the graph from the .swc file
G = get_graph_from_swc(filename);
x = G.Nodes.Coord(:,1);
y = G.Nodes.Coord(:,2);
z = G.Nodes.Coord(:,3);

[~,~,~,brchLst,~,~, ~,~,~,~]=getGraphStructure(filename,0,0);
brchLst = cell2mat(brchLst)';
%G.Nodes.Coord(brchLst,:)
%length(brchLst)
%pause
global wNaK wVDCC leakyon

% for turning on or off sodium and potassium and vdcc
wNaK = 1;
wVDCC = 1;
leakyon = 0;

% this will assign the radii size to the geometry
global dendR ermR
dendR = abs(G.Nodes.Size);
%dendR =dendR.*0 + 0.125;%0.4; 
dendR = dendR.*1e-6;     % dendritic radius [m]
ermR = dendR.*(15.0/40.0);                      % er radius [m]
n = length(ermR);

% simulation parameters
dt = 10.0e-6;                                   % delta t in seconds;
endTime = 4.0;                                  % end time in seconds;
nT = ceil(endTime/dt);                          % number of time steps

global save_every
save_every = 10;                                % how often should an output file be made
%% calcium simulation constants
% Diffusion constants
Diff_c = 220e-12;                               % cytosolic calcium [m2/s]
Diff_b = 20e-12;                                % calbindin         [m2/s]
Diff_ce = 10e-12;                               % ER calcium        [m2/s]
Diff_be = 27e-12;                               % Calreticulin      [m2/s]

% rate constants
global kbpos kbneg kbpose kbnege co btot cceq ceeq betot
kbpos = 27e3;                                   % 1/((mol/m3).s)
kbneg = 19.0;                                   % 1/s
kbpose = 100;                                   % 1/((mol/m3).s)
kbnege = 200;                                   % 1/s

% initialize solution vector
btot = 40e-3;                                   % mol/m3 total Calbindin28k
cceq = 50e-6;                                   % mol/m3 Equilibrium Cytosolic Ca conc.
ceeq = 250e-3;                                  % mol/m3 Equilibrium ER Ca conc.
betot = 3600e-3;                                % mol/m3 total ER buffer
b = kbneg*btot/(kbneg + kbpos*cceq);
be = kbnege*betot/(kbnege + kbpose*ceeq);
co = 1.0;                                       % mol/m3, this is the extracellular calcium

% now make basic diffusion stencil using edge lengths this is for the calcium
% diffusion matrices
M = make_diffusion_stencil(G);

% now make diffusion matrices using coefficients
DC = Diff_c.*M; DE = Diff_ce.*M;                        % the diffusion matrices are the same for CYT and ER
DB = Diff_b.*M; DBE = Diff_be.*M;                       % the diffusion matrix for calbindin and ER buffer

%% Specify Synaptic Input Indices, load .mat file
load(sprintf('synapses%i.mat',load_syn),'synapses');

% synapses.start_time = 0.020;
% synapses.end_time = 0.030;
% synapses.at_node = 1;
% synapses.amp = 1e-6;

%% Here assemble the system matrix for all diffusion terms
% System matrix is a large diagonal block matrix
%  [DC       ] for Cytosolic calcium
%  [  DE     ] for ER calcium
%  [    DB   ] for Calbindin
%  [      DBE] for ER buffer
S = blkdiag(DC,DE,DB,DBE);                      % this function makes a giant block diagonal matrix
[nr,~] = size(S);                               % what is the number of rows of this matrix

% this is based on SBDF2 method
LHS = eye(nr) -(2.0/3.0).*dt.*S;
%RHS = eye(nr);

% using decomposition(..) allows MatLab to efficiently solve systems
dLHS = decomposition(LHS);     

%% Initialize the solution vector
% Make these global variables for quickly accessing the indices of the
% solution for different variables.
global cci cei bbi bbei
cci = 1:n;                          % indices of cyt solution
cei = n + cci;                      % indices of er calcium solution
bbi = n + cei;                      % indices of calbin solution
bbei = n + bbi;                     % indices of ER buffer solution

% initialize the solution arrays with initial concentrations
u = zeros(nr,1);                    % this initializes the solution vector
u(cci) = cceq; u(cei) = 25e-3;       % set initial concentrations 
u(bbi) = b;  u(bbei) = be;          % set initial concentrations

%for SBDF2 need prior state
u0 = u;

% these are the initial states for the Ryanodine receptors, scaled
% accordingly
RyRProb = getInitProb(cceq);
o1 = ones(n,1).*RyRProb(1); o2 =ones(n,1).*RyRProb(2);
c1 = ones(n,1).*RyRProb(3); c2 = 1-o1-o2-c1;

ryr_states = [o1, o2, c1, c2];
ryr_states0 = ryr_states;           %for SBDF2 need prior state
cnvtoL = 1e-3;                      % convert back to mol/L

% Leakage due to Amyloid Beta pores
global a m kbeta
a = 5e-6;               % nmol/m3, Conc. of Amyloid Beta
m = 4;                  % Cooperative factor
kbeta = 1;              % 1/s, Rate constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------- Setup Voltage Problem --------------------- %%
% voltage simulation constants
global R C gk ek gna ena gl el

R=750*1e-2;                % resistance                   [ohm.m]
C=1e-2;                    % capacitance                  [F/m2]
gk=5*1e1;                  % potassium ion conductance    [S/m2]
ek=-90*1e-3;               % potassium reversal potential [V]
gna=50*1e1;                % sodium ion conductance       [S/m2]
ena=50*1e-3;               % sodium reversal potential    [V]
gl= 0.05;                  % leak conductance             [S/m2]
el=-60*1e-3;               % leak potential               [V]

vstart = -60.0 *1e-3;      % initial cell voltage [V]
ni=0.00654;                % these are ion gating variables dimensionless
mi=0.00654;                % these are ion gating variables dimensionless
hi=0.9997;                 % these are ion gating variables dimensionless
qi=0.975;                  % these are ion gating variables dimensionless (for Fast Calcium Current)

% Makes diffusion matrix for second derivative in HH model equation
% Makes diffusion matrix for the second derivative of the voltage
DV = stencilMaker(n,R,dendR,C,filename);

% Setup system matrices for Voltage problem
VLHS = eye(n) -(2.0/3.0).*dt.*DV;
%VRHS = eye(n);
dVLHS = decomposition(VLHS);

% initialize Voltage and HH states and corresponding states for SBDF2
v =  u(cci).*0+vstart; v0 = v;
nn = u(cci).*0+ni; nn0 = nn;
mm = u(cci).*0+mi; mm0 = mm;
hh = u(cci).*0+hi; hh0 = hh;
qq = u(cci).*0+qi; qq0 = qq;

vdcc_states = initialize_vdcc_states(v);
vdcc_states0 = vdcc_states;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocities
global vle vlp
u_tmp = zeros(nr,1);              
u_tmp(cci) = cceq; u_tmp(cei) = ceeq;    
u_tmp(bbi) = b;  u_tmp(bbei) = be;  

jp = jP(u_tmp);                                 % PMCA Calcium Flux               
jn = jN(u_tmp);                                 % NCX Calcum Flux
[jr,~] = jR(u_tmp,ryr_states,dt);  % RyR Calcium Flux
js = jS(u_tmp);                                 % SERCA Pump Flux
[jvdc,~] = jvdcc(u(cci),vdcc_states,1*dt,dt,v); % VDCC flux

vle = (js(1) - jr(1))/(ceeq - cceq);
vlp = (jp(1) + jn(1) - jvdc(1))/(co - cceq);

vec1 = zeros(length(brchLst),nT/save_every);
vec2 = zeros(length(brchLst),nT/save_every);
vec3 = zeros(length(brchLst),nT/save_every);
vec4 = zeros(length(brchLst),nT/save_every);
vec5 = zeros(length(brchLst),nT/save_every);
vec6 = zeros(length(brchLst),nT/save_every);
vec7 = zeros(length(brchLst),nT/save_every);

rec_ind = brchLst; % measuring at the branch points
for i=1:nT
    if (mod(i,1)==0)
        fprintf('time index = %i ',i)
        fprintf('Solving at t = %0.6f [ms]..\n',(i-1)*dt*1e3)
    end
    
    if (mod(i,save_every*10)==0)
        if(vtk==1)
            W = [ u(cci).*cnvtoL, u(cei).*cnvtoL, u(bbi).*cnvtoL, u(bbei).*cnvtoL , v*1e3];
            vtkwrite(x,y,z,W,i,outfolder);
        end
    end
    
    if (mod(i,save_every)==0) 
          vec1(:,i/save_every) = u(cci(rec_ind)).*cnvtoL;
          vec2(:,i/save_every) = u(bbi(rec_ind)).*cnvtoL;
          vec3(:,i/save_every) = u(cei(rec_ind)).*cnvtoL;
          vec4(:,i/save_every) = u(bbei(rec_ind)).*cnvtoL;
          vec5(:,i/save_every) = v(rec_ind).*1e3;
          vec6(:,i/save_every) = ryr_states(rec_ind,1);
          vec7(:,i/save_every) = ryr_states(rec_ind,2);
    end
    
    % get r.h.s. component for SBDF2 solve of calcium dynamics
    [Fu0,~,~]= system_reaction(u0,v,i-1,dt,ryr_states0,vdcc_states0,synapses);
    ryr_states0 = ryr_states; vdcc_states0 = vdcc_states; % update the prior states for RyRs and VDCCs
    [Fu,ryr_states,vdcc_states] = system_reaction(u,v,i,dt,ryr_states,vdcc_states,synapses);
    
    % assemble r.h.s. for calcium dynamics
    b = (4.0/3.0).*u+(4.0/3.0)*dt.*Fu-(1.0/3.0).*u0-(2.0/3.0)*dt.*Fu0;
    
    % Get the r.h.s components for SBDF2 solve step for Voltage problem
    % this must be done before updating the concentrations since we use FFC
    % and Ica
    Fv0 = cable_reaction(v0,nn0,mm0,hh0,qq0,u0(cci),i);
    Fv = cable_reaction(v,nn,mm,hh,qq,u(cci),i);
    
    % update calcium concentration
    u0 = u;
    u = dLHS\b;   
    
    % update HH states
    tempstate = nn;
    nn = (4.0/3.0).*nn + (4.0/3.0)*dt.*stateF(nn,an(v),bn(v))-(1.0/3.0).*nn0-(2.0/3.0)*dt.*stateF(nn0,an(v0),bn(v0));
    nn0 = tempstate;
    
    tempstate = mm;
    mm = (4.0/3.0).*mm + (4.0/3.0)*dt.*stateF(mm,am(v),bm(v))-(1.0/3.0).*mm0-(2.0/3.0)*dt.*stateF(mm0,am(v0),bm(v0));
    mm0 = tempstate;
    
    tempstate = hh;
    hh = (4.0/3.0).*hh + (4.0/3.0)*dt.*stateF(hh,ah(v),bh(v))-(1.0/3.0).*hh0-(2.0/3.0)*dt.*stateF(hh0,ah(v0),bh(v0));
    hh0 = tempstate;
    
    % update FFC state
    tempstate = qq;
    qq = (4.0/3.0).*qq + (4.0/3.0)*dt.*FFC(v,qq)-(1.0/3.0).*qq0-(2.0/3.0)*dt.*FFC(v0,qq0);
    qq0 = tempstate;
    
    % assemble r.h.s. for SBDF2 step of the voltage problem
    b = (4.0/3.0).*v+(4.0/3.0)*dt.*Fv-(1.0/3.0).*v0-(2.0/3.0)*dt.*Fv0;
    v0 = v;
    v = dVLHS\b;
end
save(sprintf('%s/cc_%i.mat',datafolder,load_syn),'vec1');
save(sprintf('%s/bb_%i.mat',datafolder,load_syn),'vec2');
save(sprintf('%s/ce_%i.mat',datafolder,load_syn),'vec3');
save(sprintf('%s/bbe_%i.mat',datafolder,load_syn),'vec4');
save(sprintf('%s/vv_%i.mat',datafolder,load_syn),'vec5');
save(sprintf('%s/o1_%i.mat',datafolder,load_syn),'vec6');
save(sprintf('%s/o2_%i.mat',datafolder,load_syn),'vec7');
end

%% System Reaction term
function [yout,ryr_states,vdcc_states] = system_reaction(u,v,i,dt,ryr_states,vdcc_states,syn)
global cci wVDCC
global co dendR ermR 

n = length(cci);
% assemble the reaction terms in the update
jp = jP(u);                             % PMCA Calcium Flux               
jn = jN(u);                             % NCX Calcum Flux
jsyn = u(cci).*0;                   
if ~isempty(syn)
    jsyn = jSyn(n,i*dt,syn);         % Synaptic Calcium Flux
end
jlp = jLP(u,co,syn,dt*i);                        % Plasma Mem. Leak Flux

% plot jsyn at index = 1
% plot(1:length(jsyn),jsyn)
% ylim([0 7e-6])
% drawnow

[jvdc,vdcc_states] = jvdcc(u(cci),vdcc_states,i*dt,dt,v);
if (wVDCC == 0)
    jvdc = 0.0;
end

% combine all flux for the Plasma Membrane +jsyn
JPM = (2.*dendR./(dendR.^2-ermR.^2)).*(-jp-jn+jlp+jsyn+jvdc);

% Flux for Ryanodine receptors
[jr,ryr_states] = jR(u,ryr_states,dt);
jle = jLe(u);                           % ER Leak Flux
js = jS(u);                             % SERCA Pump Flux
JERM = jr+jle-js;                       % combine all  ER fluxes

jsoc = jSOC(u); 

yout = [cyt_reaction(u)+JPM+(2.*ermR./(dendR.^2-ermR.^2)).*JERM; ...
    cyt_reaction_er(u)-1.*(2./(ermR)).*JERM + 1.*(2./(ermR)).*jsoc; ...
    cyt_reaction(u); ...
    cyt_reaction_er(u);
    ];
end

%% calcium simulation functions
function v = cyt_reaction(u)
global cci bbi
global kbneg btot kbpos

v = kbneg.*(btot - u(bbi))-kbpos.*u(bbi).*u(cci);
end

function v = cyt_reaction_er(u)
global cei bbei
global kbnege betot kbpose

v = kbnege.*(betot - u(bbei))-kbpose.*u(bbei).*u(cei);
end

function v = jSyn(n,t,syn)
v = zeros(n,1);

% twindow=[0.02 0.025];
% wend = twindow(2)+0.01;
% wstart =wend + 0.2;
wstart = 0.2; % this corresponds to the frequency, i.e. 10 --> 1 pulse every 10 seconds
               % i.e. 0.2 --> 5 pules every 1 second
               
for i = 1:length(syn)
    tt = mod(t,wstart);
    
    if (( tt>= syn(i).start_time) && (tt <= syn(i).end_time))
        m = -syn(i).amp / ( syn(i).end_time - syn(i).start_time  );
        v(syn(i).at_node) = v(syn(i).at_node) + m*(tt-syn(i).start_time) + syn(i).amp;
        
        %if (t>= 1.0)
        %    v(syn(i).at_node) = 0.0;
        %end
    end    
end

end

function v = jS(u)
global cci cei
cc = u(cci); ce = u(cei);

rhos = 2390e12;             % 1/m2
Is = 6.5e-24;               % mol2/(m3.s) 
Ks = 180e-6;                % mol/m3

v = rhos.*Is.*cc./((Ks+cc).*ce);
end

function v = jLe(u)
global cci cei vle a m kbeta
cc = u(cci); ce = u(cei);

flag = 0;

v = vle.*(ce-cc)+ kbeta*a^m*flag;
end

function [v,o] = jR(u,os,dt)
global cci cei
cc = u(cci); ce = u(cei);

% parameters for current calculation
rhor = 3e12;             %1/m2
Iref = 3.5e-18;            % mol/s
ceref = 250e-3;            % mol/m3 
Ir =Iref.*(ce-cc)./ceref;

% parameters for kinetics
kaneg = 28.8;           % 1/s
kapos = 1500e12;        % m12/(mol4.s)
kbneg = 385.9;          % 1/s
kbpos = 1500e9;         % m9/(mol3.s)
kcneg = 0.1;            % 1/s
kcpos = 1.75;           % 1/s

% these are the incoming state values
o1 = os(:,1); o2 = os(:,2); c1 = os(:,3); c2 = os(:,4); 

num_steps = 1000;
k = dt/num_steps; % this is the time step size for solving RyR ode.

% using BE, that is backward euler
for i=1:num_steps
    c1 = (k*kaneg.*o1 + c1)./(1+k*kapos.*cc.^4);
    o2 = (k*kbpos.*(cc.^3).*(o1)+o2)./(1+k*kbneg);
    c2 = (k*kcpos.*o1+c2)./(1+k*kcneg);    
    o1 = 1-c1-o2-c2; % update o1 last, because you already have it at the beginning of the step.
end
o = [o1 o2 c1 c2];
pr0 = o1 + o2;
v = rhor.*pr0.*Ir;
end

function v = jP(u)
global cci
rhop = 500e12;      % m-2 
Ip =1.7e-23;        % mol/2
Kp = 60e-6;         % mol/m3         

cc = u(cci);

v=rhop.*Ip.*cc.^2./(Kp*Kp+cc.^2);
end

function v = jN(u)
global cci
rhon = 15e12;       % m-2  
In =2.5e-21;        % mol/s
Kn = 1.8e-3;        % mol/m3
cc = u(cci);

v = rhon.*In.*cc./(Kn+cc);
end

function v=jLP(u,co,syn,t)
global cci vlp a m kbeta leakyon
cc =u(cci);
flag = cci'*0;

if t>=0.000
    for i=1:length(syn)
        if syn(i).amp==0
            flag(syn(i).at_node)=1;
        end
    end
end

leakyI = kbeta*a^m*flag;
%plot(leakyI);
%drawnow

v = vlp.*(co-cc) + leakyI.*leakyon;
end

function v = jSOC(u)
global cci cei ceeq co dendR
Isoc = 2.1*10^-15; % Ampere 
z = 2; % valency of Ca
Ao = 0.25*10^-18; % m^2 Area of ORAI channel
F = 96485; % C/mol Faraday's constant       

flag = zeros(length(cci),1);

for i = 1:length(cci)
    if 247.602e-3 < u(i+length(cci))
        flag(i) = exp(ceeq*1000-u(i+length(cci))*1000)-1;
    else
        flag(i) = 10;
    end
end

rhosoc = 10^-8.*flag;   % m-2, active SOC
rhosoc = rhosoc.*dendR/(0.4*1e-6); % scale for radius 40 mum -> rhosoc, R mum -> ?
jSOC = rhosoc*Isoc.*log(10*co./u(cei))./(F*z*Ao);

% v = jSOC.*flag;
v = jSOC;
end

function RyRProb = getInitProb(cc)
% This returns intital/equilibirum probabilites for RyR ODEs
kaneg = 28.8;           % 1/s
kapos = 1500e12;        % m12/(mol4.s)
kbneg = 385.9;          % 1/s
kbpos = 1500e9;         % m9/(mol3.s)
kcneg = 0.1;            % 1/s
kcpos = 1.75;           % 1/s

A = [1 1 1 1; kaneg 0 -kapos*cc^4 0; kbpos*cc^3 -kbneg 0 0; kcpos 0 0 -kcneg];
B = [1; 0; 0; 0];
RyRProb = linsolve(A,B);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Voltage Problem Functions %%%%%%%%%%%%%%%%%%%
function v = initialize_vdcc_states(v)
Rgas = 8.314;   % universal gas constant Joules/(K.mol)
Far = 96485;    % Faraday constant Coulombs/mol
T = 310;        % Temperataure in Kelvin

Va = -21*1e-3; Vb = -40*1e-3;
za = 2 ; zb =1;

s1 = 1.0./(1.0 + exp(-za.*(v-Va).*Far./(Rgas*T))); 
s2 = 1.0./(1.0 + exp(-zb.*(v-Vb).*Far./(Rgas*T))); 

v = [s1,s2];
end

%% voltage simulation functions
function yout = cable_reaction(V,n,m,h,q,cc,i)
global C gk ek gna ena gl el
global wNaK
cyt_thres = 5e-3;

yout = (gl/C).*(V-el);  

if wNaK == 1
    yout = yout + (gk/C).*n.^4.*(V-ek);
    yout = yout + (gna/C).*m.^3.*h.*(V-ena);
end

 Ic = Ica(V,cc,q);
% plot(Ic);
% drawnow

yout = -1.*yout+ Ic.*sigmoid(cc,cyt_thres);%(cc>cyt_thres);
end

function v = sigmoid(cc,cyt_thres)

a =cyt_thres;
c = -1e6;

v= 1./(1+exp(c.*(cc-a)));
end

%Fast Calcium Current for HH
function v = Ica(V,c,q)
global C co

Eca = 12.5.*log(co./c).*1e-3;
gca = 5.16*6; % [S/m2]
K = 0.01;   % [mol/m3]
h = K./(K+c);

v =(-1.*gca/C).*q.*h.*(V-Eca); 
end

function v = FFC(V,q)
minf = (1+exp(-1.*((V.*1e3)+6)./8)).^(-1);
taum = (7.8*1e-3)./( exp(((V.*1e3)+6)./16)+exp(-1.*((V.*1e3)+6)./16));

v = (minf - q)./taum;
end

function yout = stateF(s,a,b)
    yout = a.*(1-s)-b.*s;
end
% These functions are for the gating variables in the Hodgkin-Huxle formulism
% carefully notice that for this simulation I am using MKS, these gating
% function use [mV] and output [ms]^-1 so they need to be properly scaled!!
function out=an(vin)
vin = vin.*1e3;
out=(-0.032).*(vin-15)./(exp(-1.*(vin-15)./5)-1);
out = out*1e3;
end

function out=bn(vin)
vin = vin.*1e3;
out=(0.5).*exp(-1.*(vin-10)./40);
out = out*1e3;
end

function out=am(vin)
vin = vin.*1e3;
out=(-0.32).*(vin-13)./(exp(-1.*(vin-13)./4)-1);
out = out*1e3;
end

function out=bm(vin)
vin = vin.*1e3;
out=(0.28).*(vin-40)./(exp((vin-40)./5)-1);
out = out*1e3;
end

function out=ah(vin)
vin = vin.*1e3;
out=(0.128).*exp(-1.*(vin-17)./18);
out = out*1e3;
end

function out=bh(vin)
vin = vin.*1e3;
out=4./(exp((40-vin)./5)+1);
out = out*1e3;
end
