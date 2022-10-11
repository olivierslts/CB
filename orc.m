function [out] = orc(parametre)

if nargin == 0
    close all
    clear all
    clc
parametre.v_cmp_swept=3*120/100^3;%2*250/2.45/100^3; 9.1667e-04/2.5;%1050/100^3;%3*120/100^3;%3*120/100^3;%174*4/3/100^3;%expander swept volume [m3]
parametre.glide_ev=8;%condenser glide [K]
parametre.fluid='R1233zd(E)';%'R1233zd(E)';%Working fluid
parametre.pinch_ev=2;% evaporator pinch point [K]
parametre.epsilon=0.01;%convergence criteria [-]
parametre.pinch_cd=2;% condenser pinch point [K]
parametre.glide_cd=8;%evaporator glide [K]
parametre.rv=1.7;%2.45;% expander volume ratio [-]
parametre.eff_cmp=0.7; %expander nominal efficiency [-]
parametre.t_cd_htf_su=15; %condenser secondary fluid supply temperature [°C]
parametre.t_ev_htf_su=70; %evaporator secondary fluid supply temperature [°C]
parametre.ts=1;%plot ts diagram
parametre.rpm_cmp=6000;%3000;%expander speed[RPM]
parametre.re_egal_un=0;%1=the flow is constrained in order to obtainn the same Reynolds ratio, 0 if calculating the flow based on the expander caracteristics
parametre.m_dot_wf_from_hp=99999945*1.3;%69;%mass flow rate of the heat pump in order to fix the same reynolds ratio
parametre.dp_ev=0.1;%pressure drop in the evaporator line [bar]
parametre.dp_cd=0.1;%pressure drop in the condenser line [bar]
parametre.rv_optim=0;%1 = the volumeratio of the expander is optimized, 0=the model takes into account under and over expansion losses
end

    dp_ev=parametre.dp_ev;
    dp_cd=parametre.dp_cd;
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
eff_exp=parametre.eff_cmp;


%inputs
epsilon=parametre.epsilon;
oh=5;
sc=5;
eff_pp=0.5;
cp_sf=4200;
FF=1;%17.426*rpm_cmp^(-0.36);
v_exp_swept=v_cmp_swept/rv;


p_crit_fluid=CoolProp.PropsSI('Pcrit','T',t_ev_htf_su+273.15+20,'Q',0,fluid)/1e5;
t_crit_fluid=CoolProp.PropsSI('Tcrit','P',10e5,'Q',0,fluid)-273.15;

if t_crit_fluid<t_ev_htf_su
    flag_out=1;
   % w_dot_net=0;
    eta=0;
    pinch_cd=-10;
    p_ev_guess=-10;
    p_cd_guess=-10;
    eff_exp_cor=0;
w_dot_exp_el_cor=0;
w_dot_net=0;
q_dot_cd=0;
q_dot_ev=0;
m_dot_wf=0;
FF=0;
Re_exp_ex=0;
m_dot_wf=0;
w_dot_pp_aux=0;
w_dot_pp_el=0;
m_dot_cd_sf=0;
m_dot_ev_sf=0;
rho_exp_su=1;
rho_exp_ex=1;
else
    flag_out=0;
    p_cd_guess=CoolProp.PropsSI('P','T',t_cd_htf_su+273.15,'Q',0,fluid)/1e5; %minimum pas atteignable
pinch_cd=-10;
end

%guess


while pinch_cd<pinch_cd__ & flag_out==0
t_cd_sat=CoolProp.PropsSI('T','P',p_cd_guess*1e5,'Q',0,fluid)-273.15;
p_ev_guess=CoolProp.PropsSI('P','T',t_ev_htf_su+273.15,'Q',0,fluid)/1e5; %maximum pas atteignable


pinch_ev=-10;%guess

while pinch_ev<pinch_ev__
    %evaporator
t_ev_sat=CoolProp.PropsSI('T','P',p_ev_guess*1e5,'Q',0,fluid)-273.15;
t_ev_ex=t_ev_sat+oh;
t_ev_su=t_cd_sat-sc;
h_ev_su=CoolProp.PropsSI('H','P',p_ev_guess*1e5,'T',t_ev_su+273.15,fluid);
rho_ev_su=CoolProp.PropsSI('D','P',p_ev_guess*1e5,'T',t_ev_su+273.15,fluid);
%visc_ev_su=CoolProp.PropsSI('V','P',p_ev_guess*1e5,'T',t_ev_su+273.15,fluid);

h_ev_0=CoolProp.PropsSI('H','P',p_ev_guess*1e5,'Q',0,fluid);
h_ev_1=CoolProp.PropsSI('H','P',p_ev_guess*1e5,'Q',1,fluid);
h_ev_ex=CoolProp.PropsSI('H','P',p_ev_guess*1e5,'T',t_ev_ex+273.15,fluid);
s_ev_su=CoolProp.PropsSI('S','P',p_ev_guess*1e5,'T',t_ev_su+273.15,fluid);
s_ev_0=CoolProp.PropsSI('S','P',p_ev_guess*1e5,'Q',0,fluid);
s_ev_1=CoolProp.PropsSI('S','P',p_ev_guess*1e5,'Q',1,fluid);
s_ev_ex=CoolProp.PropsSI('S','P',p_ev_guess*1e5,'T',t_ev_ex+273.15,fluid);


rho_ev_ex=CoolProp.PropsSI('D','P',p_ev_guess*1e5,'T',t_ev_ex+273.15,fluid);
v_dot_exp_su=rpm_cmp/60*v_cmp_swept/rv;
v_dot_exp_su_m3h=v_dot_exp_su*3600;

if parametre.re_egal_un==0
m_dot_wf=v_dot_exp_su*rho_ev_ex;
rpm_exp=rpm_cmp;
else
   m_dot_wf=parametre.m_dot_wf_from_hp;
   v_dot_exp_su=m_dot_wf/rho_ev_ex;
   rpm_exp=v_dot_exp_su*60/v_cmp_swept/rv;
end


q_dot_ev=m_dot_wf*(h_ev_ex-h_ev_su);
m_dot_ev_sf=q_dot_ev/(cp_sf*glide_ev);
q_dot_ev_tpv=m_dot_wf*(h_ev_ex-h_ev_0);
t_pinch_ev_0=t_ev_htf_su-q_dot_ev_tpv/(m_dot_ev_sf*cp_sf);
t_ev_htf_ex=t_ev_htf_su-glide_ev;

pinch_ev=min((t_ev_htf_su-t_ev_ex),t_pinch_ev_0-t_ev_sat);

p_ev_guess=p_ev_guess*(1-epsilon);

end

%expander
p_exp_su=p_ev_guess-dp_ev;
p_exp_ex=p_cd_guess+dp_cd;

%if parametre.rv_optim==1
  %  rv_test=(p_exp_su/p_exp_ex)^(1/1.3);%%%%%%%%%%%%%%%%%%%
%else
%    rv=rv;
%end
h_exp_su=CoolProp.PropsSI('H','P',p_exp_su*1e5,'T',t_ev_ex+273.15,fluid);
s_exp_su=CoolProp.PropsSI('S','P',p_exp_su*1e5,'T',t_ev_ex+273.15,fluid);
rho_exp_su=CoolProp.PropsSI('D','P',p_exp_su*1e5,'T',t_ev_ex+273.15,fluid);
h_exp_ex_is=CoolProp.PropsSI('H','P',p_exp_ex*1e5,'S',s_exp_su,fluid);
w_dot_exp_el=eff_exp*m_dot_wf*(h_exp_su-h_exp_ex_is)/FF;
h_exp_ex=h_exp_su-w_dot_exp_el*FF/m_dot_wf;
t_exp_ex=CoolProp.PropsSI('T','H',h_exp_ex,'P',p_exp_ex*1e5,fluid)-273.15;
rho_exp_ex_=CoolProp.PropsSI('D','H',h_exp_ex,'P',p_exp_ex*1e5,fluid);
v_dot_exp_ex_m3h=m_dot_wf/rho_exp_ex_*3600;
rv_opt=rho_exp_su/rho_exp_ex_;
v_exp_su=1/rho_exp_su;
v_exp_in=v_exp_su*rv;
rho_exp_in=1/v_exp_in;
h_exp_in=CoolProp.PropsSI('H','D',rho_exp_in,'S',s_exp_su,fluid);
p_exp_in=CoolProp.PropsSI('P','D',rho_exp_in,'S',s_exp_su,fluid)/1e5;
 w_1=h_exp_su-h_exp_in;
 w_2=v_exp_in*(p_exp_in-p_exp_ex)*1e5;
 w_dot_1=w_1*m_dot_wf;
 w_dot_2=w_2*m_dot_wf;
if parametre.rv_optim==0
        w_dot_exp_el_cor=(w_dot_1+w_dot_2)*eff_exp/FF;
else
    w_dot_exp_el_cor=eff_exp*m_dot_wf*(h_exp_su-h_exp_ex_is)/FF;

end


 eff_exp_cor=(w_dot_1+w_dot_2)/(m_dot_wf*(h_exp_su-h_exp_ex_is))*eff_exp;

h_exp_ex_=h_exp_su-w_dot_exp_el_cor/m_dot_wf;
t_exp_ex_=CoolProp.PropsSI('T','H',h_exp_ex_,'P',p_exp_ex*1e5,fluid)-273.15;
rho_exp_ex=CoolProp.PropsSI('D','P',p_exp_ex*1e5,'T',t_exp_ex+273.15,fluid);
mu_exp_ex=CoolProp.PropsSI('V','P',p_exp_ex*1e5,'T',t_exp_ex+273.15,fluid);
v_dot_exp_ex=m_dot_wf/rho_exp_ex;
Re_exp_ex=rho_exp_ex*v_dot_exp_ex/mu_exp_ex;



%condenser
h_cd_ex=CoolProp.PropsSI('H','P',p_cd_guess*1e5,'T',t_ev_su+273.15,fluid);
h_cd_1=CoolProp.PropsSI('H','P',p_cd_guess*1e5,'Q',1,fluid);
s_cd_su=CoolProp.PropsSI('S','P',p_cd_guess*1e5,'T',t_exp_ex+273.15,fluid);
s_cd_0=CoolProp.PropsSI('S','P',p_cd_guess*1e5,'Q',0,fluid);
s_cd_1=CoolProp.PropsSI('S','P',p_cd_guess*1e5,'Q',1,fluid);
s_cd_ex=CoolProp.PropsSI('S','P',p_cd_guess*1e5,'T',t_ev_su+273.15,fluid);

q_dot_cd=m_dot_wf*(h_exp_ex_-h_cd_ex);
q_dot_cd_tpl=m_dot_wf*(h_cd_1-h_cd_ex);
m_dot_cd_sf=q_dot_cd_tpl/(cp_sf*glide_cd);
t_cd_pinch_1=t_cd_htf_su+q_dot_cd_tpl/(m_dot_cd_sf*cp_sf);
pinch_cd=min(t_cd_sat-t_cd_pinch_1,t_ev_su-t_cd_htf_su);
t_cd_htf_ex=t_cd_htf_su+glide_cd;

rho_htf=1000;
v_dot_htf=m_dot_cd_sf/rho_htf;
w_dot_pp_aux=v_dot_htf*0.4*1e5/0.7;


p_cd_guess=p_cd_guess*(1+epsilon);


[p_cd_guess p_ev_guess pinch_ev pinch_cd];

%pump
rho_pp=CoolProp.PropsSI('D','P',p_cd_guess*1e5,'T',t_ev_su+273.15,fluid);
w_dot_pp_el=m_dot_wf/rho_pp*(p_ev_guess-p_cd_guess)*1e5/eff_pp;

%cycle
w_dot_net=w_dot_exp_el_cor-w_dot_pp_el-w_dot_pp_aux;
eta_=w_dot_net/q_dot_ev;

end

if p_ev_guess>p_cd_guess
eta=max(0,eta_);
else
    eta=0;
end





if parametre.ts==1
    t=[t_ev_su t_ev_sat t_ev_sat t_ev_ex t_exp_ex t_cd_sat t_cd_sat t_ev_su t_ev_su];
    s=[s_ev_su s_ev_0 s_ev_1 s_ev_ex s_cd_su s_cd_1 s_cd_0 s_cd_ex s_ev_su];
     t_h=[ t_ev_htf_ex t_ev_htf_su];
     s_h=[s_ev_su s_ev_ex];
     t_c=[t_cd_htf_ex t_cd_htf_su];
     s_c=[s_cd_su s_cd_ex];
% 
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

out.eff_exp=eff_exp_cor;
out.w_dot_exp_el_cor=w_dot_exp_el_cor;
out.eta=eta;
out.w_dot_net=w_dot_net;
out.q_dot_cd=q_dot_cd;
out.q_dot_ev=q_dot_ev;
out.m_dot_wf=m_dot_wf;
out.FF=FF;
out.Re_exp_ex=Re_exp_ex;
out.m_dot_wf=m_dot_wf;
out.w_dot_pp_aux=w_dot_pp_aux;
out.w_dot_pp_el=w_dot_pp_el;
out.p_cd=p_cd_guess;
out.p_ev=p_ev_guess;
out.m_dot_cd_sf=m_dot_cd_sf;
out.m_dot_ev_sf=m_dot_ev_sf;
out.rv_cycle=rho_exp_su/rho_exp_ex;
out.rpm_exp=rpm_exp;
out.v_exp_swept=v_exp_swept;
out.rho_exp_su=rho_ev_ex;
out.rho_exp_ex=rho_exp_ex;
out.v_dot_exp_su=v_dot_exp_su;
out.w_dot_pp_el=w_dot_pp_el;
out.w_dot_exp_el_cor=w_dot_exp_el_cor;
out.eff_exp_cor=eff_exp_cor;
out.v_dot_exp_ex_m3h=v_dot_exp_ex_m3h;