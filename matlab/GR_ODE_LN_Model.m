function dydt=GR_ODE_LN_Model(t,y)

%% PARAMETERS %%
GR_ODE_LN_parameter_list

%% Equations

%%  DENDRITIC CELLS

%% MDC LN [y(14)] - Mature/Licensed Dendritic Cell (in the Lymph Node)
%% IDCmigration/scaling + IDCLN_maturation1 + IDCLN_maturation2 - MDCdeath
MDCdeath=muMDC_LN*y(14);

MDC_LN = - MDCdeath;

%% CD4+ and CD8+ T cells in the LN
%% Naive CD4+ T cells in the Lymph Node  N4 - [y(15)] --> N4_recruit - N4_death - N4_diff_T0
N4_recruit = sn4 + k13*y(15)*(y(14)/(y(14)+hs13)); 
N4_death = muN4*y(15);
N4_diff_T0 = k14*y(15)*y(14);

N4_LN = N4_recruit - N4_death - N4_diff_T0; % Naive CD4 LN (N4) - [15]

%% Precursor Th1 LN  [y(16)] ----> N4_diff_T0 + T0LN_prolif - T0LN_diff_T1LN_MDC - T0LN_diff_T1LN_Mac - T0LN_migr 
T0LN_prolif=k15*y(16)*(1-y(16)/rho2)*(1-(y(31)/(y(31)+hsI10_T0LN))); % 8 hrs doubling time, they double only 4 times
T0LN_migr=csi1*y(16);
T0LN_diff_T1LN_MDC = k20a*y(16)*(y(14)/(y(14)+f2*y(31)+hs20a)); 
T0LN_diff_T1LN_Mac = 0;%k29a*y(30)*y(16)*(mac_CC_LN/(mac_CC_LN+hs29a)); 

TH0_LN = N4_diff_T0 + T0LN_prolif - T0LN_migr - T0LN_diff_T1LN_Mac -...
    T0LN_diff_T1LN_MDC;% precursor TH1 LN - [16]  ...

%% Th1 LN  [y(32)]
T1LN_migr=csi1a*y(32);
%T1LN_TNF_apop = 0;%k22a*y(32)*(y(28)/(y(28)+hs22a));  
%T1LN_death_restim = muT1*y(32)*(1-k21*(mac_CC_LN/(mac_CC_LN+hs21)));  

TH1_LN = T0LN_diff_T1LN_MDC + T0LN_diff_T1LN_Mac - ...
    T1LN_migr;  %%Th1 LN LN [32]- T1LN_death_restim T1LN_TNF_apop -

%% Naive CD8+ T cells in the Lymph Node  N8 - [y(17)]  --> N8_recruit - N8_death - N8_diff_T80
N8_recruit = sn8 + k16*y(17)*(y(14)/(y(14)+hs16)); ... model 2COMP MODEL
N8_death = muN8*y(17);
wT80=5e-1;
N8_diff_T80 = k17*y(14)*y(17)*((y(32)+wT80*y(16))/(y(32)+wT80*y(16)+hs17));

N8_LN = N8_recruit - N8_death - N8_diff_T80; %Naive CD8 LN (N8) - [17]

%% CTL precursor cells in the Lymph Node  T80LN  [y(18)] --> N8_diff_T80 + T80LN_prolif - T0LN_migr
T80LN_prolif = k18*y(18)*(1-y(18)/rho3)*(1-(y(31)/(y(31)+hsI10_T80LN))); % 8 hrs doubling time, they double only 4 times
T80LN_migr = csi2*y(18);
T80LN_diff_CTLLN_MDC=k24a*y(18)*(y(14)/(y(14)+f2*y(31)+hs24a));
T80LN_diff_CTLLN_Mac=0;%k30a*y(30)*y(18)*(mac_CC_LN/(mac_CC_LN+hs30a));

T80_LN = N8_diff_T80 + T80LN_prolif - T80LN_migr - T80LN_diff_CTLLN_MDC - ...
    T80LN_diff_CTLLN_Mac;  ... T80LN - [18 ]

%% T8 [y(33)] cells in the Lymph Node,  T8LN
T8LN_migr = csi2a*y(33);
T8LN_TNF_apop=  0;%k26a*y(33)*(y(28))/(y(28)+hs26);
%T8LN_death_restim=  muT8*y(33)*(1-k25*(mac_CC_LN/(mac_CC_LN+hs25)));

T8_LN = m*(T80LN_diff_CTLLN_MDC + T80LN_diff_CTLLN_Mac) - ...
    T8LN_TNF_apop - T8LN_migr;  %%T8 LN [33] - T8LN_death_restim

%% CTL [y(34)] cells in the Lymph Node,  TCLN
CTLLN_migr = csi2b*y(34);
CTLLN_TNF_apop = 0;% k28a*y(34)*(y(28))/(y(28)+hs28a);
%CTLLN_death_restim=  muTC*y(34)*(1-k27*(mac_CC_LN/(mac_CC_LN+hs27)));

CTL_LN = m*(T80LN_diff_CTLLN_MDC + T80LN_diff_CTLLN_Mac) - CTLLN_TNF_apop - CTLLN_migr; % - CTLLN_death_restim


dydt = ... 
[0;%MR_LUNG; ... MR LUNG [1]
0;%MI_LUNG; ... MI LUNG [2]
0;%MA_LUNG; ... MA LUNG [3]
0;%MR_LN; ... MR LN [4]
0;%MI_LN; ... MI LN [5]
0;%MA_LN; ... MA LN [6]
0;%BE_LUNG;  ... BE LUNG [7] , BE_migration
0;%BI_M_LUNG; ... BI_M_LUNG [8]
0;%BE_LN; ... BE LN [9] 
0;%BI_M_LN; ... BI-M LN [10]
0;%IDC_LUNG; ... IDC LUNG - [11] 
0;%MDC_LUNG; ... MDCL, MDC LUNG [12]
0;%IDC_LN; ... IDC LN - [13] 
MDC_LN; ... MDC - LN [14]
N4_LN; ... Naive CD4 LN (N4) - [15]
TH0_LN; ... precursor TH1 LN - [16] 
N8_LN; ... ,  Naive CD8 LN (N8) - [17]
T80_LN;%;  ... T80LN - [18 ]
0;%TH0_LUNG;% Th0 Lung [19]
0;%TH1_LUNG;% Th1 Lung [20]
0;%T80_LUNG;% T80 Lung [21]
0;%T8_LUNG;% T8 Lung [22]
0;%CTL_LUNG;% TC Lung [23]
0;%TNF_LUNG; ... TNF LUNG [24]
0;%IFN_LUNG; %%... IFN-g LUNG [25]
0;%I12_LUNG; %%... IL12 LUNG [26]
0;%I10_LUNG; %%... IL10 LUNG [27]
0;%TNF_LN; ... TNF LN [28]
0;%IFN_LN; %%... IFN-g LN [29]
0;%I12_LN; %% IL12 LN [30]
0;%I10_LN;  %%IL10 LN [31]
TH1_LN;  %%Th1 LN LN [32]
T8_LN;  %%T8 LN [33]
CTL_LN; %%CTL LN [34]
0;%M2_LUNG;  % M2 in the Lung [35]
0;%M2_LN; % M2 in the LN [36]
0;%delta; %[37]
0;%deltaLN];  % [38]
scaling*T0LN_migr;% -csi1*() precursor TH1 migration from the LN [39]
scaling*T1LN_migr;% -csi1a*() TH1 migration from the LN [40]
scaling*T80LN_migr;% -csi2*() precursor T8 migration from the LN [41]
scaling*T8LN_migr;% -csi2a*() T8 migration from the LN [42]
scaling*CTLLN_migr];%CTLLN_migr;% -csi2b*() CTL migration from the LN[43]
%{
TNF_dep_recr_M; % TNF_recr_MR [35]
TNF_dep_recr_M_LN; % TNF_recr_MRLN [36]
TNF_dep_recr_IDC; % TNF_recr_IDC [37]
TNF_dep_recr_IDC_LN; % TNF_recr_IDCLN [38]
scaling*T0LN_migr*(TNF_dep_recr_T0); % TNF_recr_T0 [39]
0;%scaling*T1LN_migr*(TNF_dep_recr_T1); % TNF_recr_T1 [40]
0;%scaling*T80LN_migr*(TNF_dep_recr_T80); % TNF_recr_T80 [41]
0;%scaling*T8LN_migr*(TNF_dep_recr_T8); % TNF_recr_T8 [42]
0;%scaling*CTLLN_migr*(TNF_dep_recr_CTL); % TNF_recr_CTL [43]
0;%TNF_indep_recr_M; % Recr_MR [44]
0;%TNF_indep_recr_M_LN; % Recr_MRLN [45]
0;%TNF_indep_recr_IDC; % Recr_IDC [46]
0;%TNF_indep_recr_IDC_LN; % Recr_IDCLN [47]
0;%N4_recruit - sn4; % Recr_N4 [48]
0;%N8_recruit - sn8; % ,  Recr_N8 [49]
0;%scaling*T0LN_migr*(TNF_indep_recr_T0); % Recr_T0 [50]
0;%scaling*T1LN_migr*(TNF_indep_recr_T1); % Recr_T1 [51]
0;%scaling*T80LN_migr*(TNF_indep_recr_T80); % Recr_T80 [52]
0;%scaling*T8LN_migr*(TNF_indep_recr_T8); % Recr_T8 [53]
0;%scaling*CTLLN_migr*(TNF_indep_recr_CTL); % Recr_CTL [54]
0;%deltaTNF*N*MI_TNF_apop; % TNF_apop_MI [55]
0;%deltaTNF*N*MI_TNF_apopLN; % TNF_apop_MILN [56]
0;%T1LN_TNF_apop; % TNF_apop_T1LN [57]
0;%T1_TNF_apop; % TNF_apop_T1 [58]
0;%T8LN_TNF_apop; % TNF_apop_T8LN [59]
0;%CTLLN_TNF_apop; % TNF_apop_CTLLN [60]
0;%T8_TNF_apop; % TNF_apop_T8 [61]
0;%CTL_TNF_apop; % TNF_apop_CTL [62]
0;%M2_infection; % Infection_MR [63]
0;%M2_infectionLN; % Infection_MRLN [64]
0;%M1_activation; % Activation_MR [65]
0;%M1_activationLN; % Activation_MRLN [66]
0;%IDC_uptake; % Infection_IDC [67]
0;%N*MI_burst; % Bursting_MI [68]
0;%N*MI_burstLN; % Bursting_MILN [69]
0;%deltaFAS*N*MI_FAS_apop; % Fas_MI [70]
0;%deltaFAS*N*MI_FAS_apopLN; % FAS_MILN [71]
0;%N*MI_CTL_killing; % CTLkilling_MI [72]
0;%N*MI_CTL_killingLN; % CTLkilling_MILN [73]
0;%muMI*N*y(2); % Nat_Death_MI [74]
0;%muMI*N*y(5); % Nat_Death_MILN [75]
0;%muMI*N*y(12); % Nat_Death_MDCL [76]
0;%muMI*N*y(14); % Nat_Death_MDC [77]
0;%IDCLN_maturation1 + IDCLN_maturation2; % Maturation_IDCLN [78]
0;%fi*IDCmigration/scaling; % Migration_IDC [79]
0;%IDCLN_maturation1; % IDC uptake in the LN [80]
0;%IDCLN_maturation2; % Exosome in the LN [81]
0;%BI_DC_LUNG; % BI_DC_LUNG [82]
0;%BI_DC_LN; % BI_DC_LN [83] 
0;%M2_LUNG;  % M2 in the Lung [84]
M2_LN; % M2 in the LN [85]
delta; %[86]
deltaLN];  % [87]
%}