function [out]  = hp(parametre)

if nargin == 0
  close all
  clc
  clear all
parametre.glide_cd=14;%14%glide of temperature of condenser [K]
parametre.v_cmp_swept=120*2/100^3;%2*120/100^3;%512/100^3;%30/100^3*6.4*8/3; %swept volume of the compressor [m3]
parametre.t_cd_htf_su=62; % secondary fluid supply temperature of the condenser [°C]
parametre.t_ev_htf_su=62; % secondary fluid supply temperature of the evaporator [°C]
parametre.fluid='R1233zd(E)';%R1234yf';%'R1233zd(E)';%working fluid
parametre.rpm_cmp=6000;%3000;%compressor speed [RPM]
parametre.pinch_ev=2;%pinch evaporator [K]
parametre.pinch_cd=2;% pinch condenser [K]
parametre.glide_ev=5;%glide of temperature of the evaporator [K]
parametre.rv_optim=0;%rv_optim=0 means the volume ratio of the compressor is considered leading to under and overexpansion losses.
parametre.rv=1.7;%volume ratio of the compressor (only considered if rv_optim=0)
parametre.eff_cmp=0.7;%nominal maximum efficieny of the compressor
parametre.ts=1;%ts=1 plot the ts diagram
parametre.sc=5;%subcooling [K]
parametre.heat=1;%heat=1 is the hot configuration. heat=0 is the cold configuration.
parametre.epsilon=0.01;%convergence parameter. keep as low as possible.
parametre.dp_ev=0.1;%pressure drop in the evaporator line [bar]
parametre.dp_cd=0.1;%pressure drop in the condenser line [bar]
end

flag_t_crit=0;
t_cd_htf_su=parametre.t_cd_htf_su;
t_ev_htf_su=parametre.t_ev_htf_su;
v_cmp_swept=parametre.v_cmp_swept;
fluid=parametre.fluid;
rpm_cmp=parametre.rpm_cmp;
pinch_ev__=parametre.pinch_ev;
pinch_cd__=parametre.pinch_cd;
glide_ev=parametre.glide_ev;
glide_cd=parametre.glide_cd;
rv=parametre.rv;
eff_cmp=parametre.eff_cmp;
FF=1; %17.426*rpm_cmp^(-0.36);
sc=parametre.sc;
    dp_ev=parametre.dp_ev;
    dp_cd=parametre.dp_cd;
%inputs
marge=10;
epsilon=parametre.epsilon;%%%%%%%%%
critere_flag=0.2;

oh=5;
cp_sf=4200;

v_exp_swept=v_cmp_swept/rv;


%parameters


% this part of code ensure the code to converge in case the temperature is above critical temperature
p_crit_fluid=CoolProp.PropsSI('Pcrit','T',t_cd_htf_su+273.15,'Q',0,fluid)/1e5;
t_crit_fluid=CoolProp.PropsSI('Tcrit','P',10e5,'Q',0,fluid)-273.15;
if t_crit_fluid<t_cd_htf_su+marge;
    flag_out=1;
cop=0;
w_dot_cmp_el=0;
q_dot_cd=0;
eff_cmp_cor=0;
q_dot_cd=0;
q_dot_ev=0;
Re_cmp_su=0;
q_dot_ev=0;
m_dot_wf=0;
w_dot_pp_aux=0;
m_dot_cd_htf=0;
FF=0;
p_cd_guess=0;
p_ev_guess=0;
m_dot_htf_ev=0;
rho_cmp_ex=0;
rho_ev_ex=0;
q_ev_su=0;
else
    flag_out=0;
    p_cd_guess=CoolProp.PropsSI('P','T',t_cd_htf_su+273.15,'Q',0,fluid)/1e5; %minimum pas atteignable

end


pinch_cd=0;%%%%%%%%%%%%


while pinch_cd<pinch_cd__ & flag_out==0 

%condenser
t_cd_sat=CoolProp.PropsSI('T','P',p_cd_guess*1e5,'Q',0,fluid)-273.15;
t_cd_ex=t_cd_sat-sc;
h_cd_ex=CoolProp.PropsSI('H','T',t_cd_ex+273.15,'P',p_cd_guess*1e5,fluid);
s_cd_ex=CoolProp.PropsSI('S','T',t_cd_ex+273.15,'P',p_cd_guess*1e5,fluid);
s_cd_1=CoolProp.PropsSI('S','Q',1,'P',p_cd_guess*1e5,fluid);
s_cd_0=CoolProp.PropsSI('S','Q',0,'P',p_cd_guess*1e5,fluid);

%evaporator
p_ev_guess=CoolProp.PropsSI('P','T',t_ev_htf_su+273.15,'Q',0,fluid)/1e5; %minimum pas atteignable


pinch_ev=-10;%%%%%%%%%%%%%%%
while pinch_ev<pinch_ev__
t_ev_sat=CoolProp.PropsSI('T','P',p_ev_guess*1e5,'Q',0,fluid)-273.15;
t_ev_ex=t_ev_sat+oh;
rho_ev_ex=CoolProp.PropsSI('D','T',t_ev_ex+273.15,'P',p_ev_guess*1e5,fluid);
mu_ev_ex=CoolProp.PropsSI('V','P',p_ev_guess*1e5,'T',t_ev_ex+273.15,fluid);
v_dot_cmp=rpm_cmp/60*v_cmp_swept;
v_dot_cmp_m3h=v_dot_cmp*3600;
m_dot_wf=rho_ev_ex*v_dot_cmp;
h_ev_su=h_cd_ex;
q_ev_su=CoolProp.PropsSI('Q','H',h_ev_su,'P',p_ev_guess*1e5,fluid);
h_ev_ex=CoolProp.PropsSI('H','T',t_ev_ex+273.15,'P',p_ev_guess*1e5,fluid);
s_ev_ex=CoolProp.PropsSI('S','T',t_ev_ex+273.15,'P',p_ev_guess*1e5,fluid);
s_ev_su=CoolProp.PropsSI('S','H',h_ev_su,'P',p_ev_guess*1e5,fluid);
s_ev_1=CoolProp.PropsSI('S','Q',1,'P',p_ev_guess*1e5,fluid);
q_dot_ev=m_dot_wf*(h_ev_ex-h_ev_su);
t_ev_htf_ex=t_ev_htf_su-glide_ev;
m_dot_htf_ev=q_dot_ev/(glide_ev*cp_sf);
%pinch_ev=min(t_ev_ex-t_ev_htf_ex,t_ev_sat-t_ev_htf_su);
pinch_ev=min(t_ev_htf_su-t_ev_ex,t_ev_htf_ex-t_ev_sat);
p_ev_guess=(1-epsilon)*p_ev_guess;
end

%compressor
p_cmp_su=p_ev_guess-dp_ev;
p_cmp_ex=p_cd_guess+dp_cd;


h_cmp_su=CoolProp.PropsSI('H','T',t_ev_ex+273.15,'P',p_cmp_su*1e5,fluid);
s_cmp_su=CoolProp.PropsSI('S','T',t_ev_ex+273.15,'P',p_cmp_su*1e5,fluid);
rho_cmp_su=CoolProp.PropsSI('D','T',t_ev_ex+273.15,'P',p_cmp_su*1e5,fluid);
try
h_cmp_ex_s=CoolProp.PropsSI('H','S',s_cmp_su,'P',p_cmp_ex*1e5,fluid);
catch
h_cmp_ex_s=CoolProp.PropsSI('H','S',s_cmp_su*0.9,'P',p_cmp_ex*1e5,fluid);
flag_out=1;
flag_t_crit=1;
end
w_dot_cmp_el_=m_dot_wf*(h_cmp_ex_s-h_cmp_su)/eff_cmp*FF;
%h_cmp_ex=h_cmp_su+w_dot_cmp_el_/eff_cmp/m_dot_wf;
h_cmp_ex=h_cmp_su+w_dot_cmp_el_/m_dot_wf;
t_cmp_ex=CoolProp.PropsSI('T','H',h_cmp_ex,'P',p_cmp_ex*1e5,fluid)-273.15;
s_cmp_ex=CoolProp.PropsSI('S','H',h_cmp_ex,'P',p_cmp_ex*1e5,fluid);
v_cmp_su=1/rho_cmp_su; %CoolProp.PropsSI('V','P',p_ev_guess*1e5,'T',t_ev_ex+273.15,fluid);
v_cmp_in=v_cmp_su/rv;
rho_cmp_in=1/v_cmp_in;
h_cmp_in=CoolProp.PropsSI('H','D',rho_cmp_in,'S',s_cmp_su,fluid);
p_cmp_in=CoolProp.PropsSI('P','D',rho_cmp_in,'S',s_cmp_su,fluid)/1e5;
w_1=h_cmp_in-h_cmp_su;
w_2=v_cmp_in*(p_cmp_ex-p_cmp_in)*1e5;
w_dot_1=w_1*m_dot_wf;
w_dot_2=w_2*m_dot_wf;

eff_cmp_cor=eff_cmp*(m_dot_wf*(h_cmp_ex_s-h_cmp_su))/(w_dot_1+w_dot_2);

if parametre.rv_optim==1
w_dot_cmp_el_cor=w_dot_cmp_el_;
else
w_dot_cmp_el_cor=w_dot_cmp_el_/eff_cmp_cor;%(w_dot_1+w_dot_2)/eff_cmp*FF;
end
h_cmp_ex_=h_cmp_su+w_dot_cmp_el_cor/eff_cmp/m_dot_wf;


Re_cmp_su=rho_ev_ex*v_dot_cmp/mu_ev_ex;
rho_cmp_ex=CoolProp.PropsSI('D','T',t_cmp_ex+273.15,'P',p_cd_guess*1e5,fluid);


%condenser
h_cd_1=CoolProp.PropsSI('H','Q',1,'P',p_cd_guess*1e5,fluid);
q_dot_cd_=m_dot_wf*(h_cmp_ex_-h_cd_ex);
q_dot_cd_tpl=m_dot_wf*(h_cd_1-h_cd_ex);
m_dot_cd_htf=q_dot_cd_/(cp_sf*glide_cd);
t_pinch_cd_1=t_cd_htf_su+q_dot_cd_tpl/(m_dot_cd_htf*cp_sf);
pinch_cd=min(t_cd_sat-t_pinch_cd_1,t_cd_ex-t_cd_htf_su);
t_cd_htf_ex=t_cd_htf_su+glide_cd;

rho_htf=1000;
v_dot_htf=m_dot_htf_ev/rho_htf;
w_dot_pp_aux=v_dot_htf*0.4*1e5/0.7;

w_dot_cmp_el_cor=w_dot_cmp_el_cor+w_dot_pp_aux;

%cycle
if parametre.heat==1
cop_=q_dot_cd_/w_dot_cmp_el_cor;
else
cop_=q_dot_ev/w_dot_cmp_el_cor;
end

p_cd_guess=p_cd_guess*(1+epsilon);

 if p_cd_guess>p_crit_fluid
     flag_out=1;
     cop=0;
end

%[pinch_ev p_ev_guess pinch_cd p_cd_guess]
if p_ev_guess<p_cd_guess
cop=cop_;
w_dot_cmp_el=w_dot_cmp_el_cor;
q_dot_cd=q_dot_cd_;
else
    cop=0;
    w_dot_cmp_el=0;
    q_dot_cd=0;
end
end



if parametre.ts==1
    t=[t_ev_sat t_ev_sat t_ev_ex t_cmp_ex t_cd_sat t_cd_sat t_cd_ex t_ev_sat];
    s=[s_ev_su s_ev_1 s_ev_ex s_cmp_ex s_cd_1 s_cd_0 s_cd_ex s_ev_su];
    t_h=[ t_cd_htf_ex t_cd_htf_su];
    s_h=[s_cmp_ex s_cd_ex];
    t_c=[t_ev_htf_ex t_ev_htf_su];
    s_c=[s_ev_su s_ev_ex];

for tt=1:100
    t_sat(tt)=-20+(tt-1)*(90-(-20))/(100-1);
    t_sat_iter=t_sat(tt);
    s_sat_0(tt)=CoolProp.PropsSI('S','Q',0,'T',t_sat_iter+273.15,fluid);
    s_sat_1(tt)=CoolProp.PropsSI('S','Q',1,'T',t_sat_iter+273.15,fluid);
end

figure
plot(s,t,'g')
hold on
plot(s_h,t_h,'r')
hold on
plot(s_c,t_c,'b')
hold on
plot(s_sat_0,t_sat,'k')
hold on
plot(s_sat_1,t_sat,'k')
grid on

 xlabel('Entropy [J/(K.kg]','FontSize',12) % x-axis label
 ylabel('Temperature [°C]','FontSize',12) % y-axis label
 set(gcf,'color','w')
 ylim([0 95])
out.s=s;
out.t=t;
out.t_c=t_c;
out.t_h=t_h;
out.s_h=s_h;
out.s_c=s_c;
out.s_sat_0=s_sat_0;
out.t_sat=t_sat;
out.s_sat_1=s_sat_1;
end


if flag_t_crit==0;

out.eff_cmp=eff_cmp_cor;
out.cop=cop;
out.w_dot_cmp_el=w_dot_cmp_el;
out.q_dot_cd=q_dot_cd;
out.q_dot_ev=q_dot_ev;
out.Re_cmp_su=Re_cmp_su;
out.q_dot_ev=q_dot_ev;
out.m_dot_wf=m_dot_wf;
out.w_dot_pp_aux=w_dot_pp_aux;
out.m_dot_cd_htf=m_dot_cd_htf;
%out.overcharge=w_dot_cmp_el/parametre.p_el_scroll_nom;
out.flag_out=flag_out;
out.FF=FF;
out.w_dot_pp_aux=w_dot_pp_aux;
out.p_cd=p_cd_guess;
out.p_ev=p_ev_guess;
out.m_dot_ev_sf=m_dot_htf_ev;
out.rv_cycle=rho_cmp_ex/rho_ev_ex;
out.q_ev_su=q_ev_su;
out.rho_cmp_su=rho_cmp_su;
out.rho_cmp_ex=rho_cmp_ex;
out.v_dot_cmp=v_dot_cmp;

else
out.eff_cmp=0;
out.w_dot_cmp_el=0;
out.q_dot_cd=0;
out.q_dot_ev=0;
out.Re_cmp_su=0;
out.q_dot_ev=0;
out.m_dot_wf=0.1;
out.w_dot_pp_aux=0;
out.m_dot_cd_htf=0;
out.overcharge=0;
out.flag_out=flag_out;
out.FF=0;
out.w_dot_pp_aux=0;
out.p_cd=0;
out.p_ev=0;
out.m_dot_ev_sf=0;
out.rv_cycle=0;
out.q_ev_su=0;
out.v_dot_cmp_m3h=v_dot_cmp_m3h;
out.cop=0;
end