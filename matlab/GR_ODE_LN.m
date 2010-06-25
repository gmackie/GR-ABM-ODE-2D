function yODE=GR_ODE_LN(inputs) % inputs: vector of 43 variables

%% inputs is a string
inputs=str2num(inputs);

%% Set time span of the simulation
t_init=0; 
t_end=1/144; % length of the simulation=10 minutes
tspan=[t_init t_end];

%% RATES and MAX RATES (k)  
%GR_ODE_LN_parameter_list

f=@GR_ODE_LN_Model; %
options=odeset('MaxStep',t_end-t_init);%,'OutputFcn','odeplot','OutputSel',[39:43]);% 14=MDC
[t,y]=ode113(@(t,y)f(t,y),tspan,inputs,options);

%% Tgamma: y(20)+y(22), Ctl: y(23), Tregs: 10% of Tgamma
[a b]=size(y);
%yODE=y(a,:);

% yODE is the empty string, as I want to print the output in my own way
yODE='';

for i=1:43
    fprintf(1, '%d\n',y(a,i));
end
