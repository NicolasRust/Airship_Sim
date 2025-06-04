%--- parameters
Ts = 0.1;           % match your Simulink sample time
t_end = 70;         % total sim time

%--- simulate Simulink model (adjust model name & data logging)
t     = out.u_ref.Time;
u_ref = out.u_ref.Data;
u_act = out.u.Data;

%--- build experiments
exp1 = iddata( u_act(t>=0 & t<20), u_ref(t>=0 & t<20), Ts );
exp2 = iddata( u_act(t>=20 & t<40), u_ref(t>=20 & t<40), Ts );
exp3 = iddata( u_act(t>=40 & t<60), u_ref(t>=40 & t<60), Ts );
data = merge(exp2, exp3);

%--- estimate TF
sys_id1 = tfest(exp1, 1, 0)
sys_id2 = tfest(data, 1, 0)

%--- display results
disp(sys_id)
figure; bode(sys_id); grid on;