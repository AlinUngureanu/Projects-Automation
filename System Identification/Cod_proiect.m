t = ungureanu.X.Data;
u = double(ungureanu.Y(1,3).Data'); % intrare
w = double(ungureanu.Y(1,2).Data'); % viteza
y = double(ungureanu.Y(1,1).Data'); % pozitie

sgtitle('Semnalele initiale vazute separat');
subplot(311); plot(t,u); title('Factor de umplere'); xlabel('t(sec)'); ylabel('u(%)');
subplot(312); plot(t,w); title('Viteza unghiulara'); xlabel('t(sec)'); ylabel('w(rad/s)');
subplot(313); plot(t,y); title('Pozitia unghiulara'); xlabel('t(sec)'); ylabel('y(impulsuri)');

figure
plot(t,[u*200,w,y/2]); title('Semnalele initiale'); xlabel('t(sec)'); 
ylabel('u(%), w(rad/s), y(impulsuri)'); grid
legend('u(%)', 'w(rad/s)', 'y(impulsuri)')

tm = t(2769)-t(2717);
wst = mean(w(3081:3246));
ust = mean(u(3081:3246));
K = wst/ust
w63 = 0.63*(wst);

figure
plot(t,[u*200,w,y/2]); title('Semnalele initiale'); xlabel('t(sec)'); 
ylabel('u(%), w(rad/s), y(impulsuri)'); grid
hold on
plot(t,w63*ones(1,length(t)),'--c')
T = t(2802)-t(2717);
N = round(tm/(t(2)-t(1)));
uN = [u(1)*ones(N,1); u(1:end-N)];

A = -1/T;
B = K/T;
C = 1;
D = 0;

ysim = lsim(A,B,C,D,uN,t,w(1));
hold on
plot(t,ysim,'g')
legend('u(%)', 'w(rad/s)', 'y(impulsuri)','w63','w dupa identificare')
emin = norm(w-ysim)/norm(w-mean(w))*100

%% cu interpolare
tm = t(2769)-t(2717);
wst = mean(w(3081:3246));
ust = mean(u(3081:3246));
K = wst/ust
w63 = 0.63*(wst);

figure
plot(t,[u*200,w,y/2]); title('Semnalele initiale'); xlabel('t(sec)'); 
ylabel('u(%), w(rad/s), y(impulsuri)'); grid
hold on
plot(t,w63*ones(1,length(t)),'--c')

t63 = (w63-w(2801))*(t(2802)-t(2801))/(w(2802)-w(2801))+t(2801);

T = t63-t(2717);
N = round(tm/(t(2)-t(1)));
uN = [u(1)*ones(N,1); u(1:end-N)];

A = -1/T;
B = K/T;
C = 1;
D = 0;

ysim = lsim(A,B,C,D,uN,t,w(1));
hold on
plot(t,ysim,'g')
legend('u(%)', 'w(rad/s)', 'y(impulsuri)','w63','w dupa identificare')
emin = norm(w-ysim)/norm(w-mean(w))*100

%% autocorelatia la treapta identificata prin metoda neparametrica
% % for i = 1:length(t)
% %     e(i) = y(i)-ysim(i); % y_real - y_prezis
% % end
% % 
% % s=0;
% % for i = 1:length(t)
% %     s = s+e(i)*e(i);
% % end
% % 
% % R0 = 1/length(t)*s;
% % 
% % % R0 = 1/length(t)*sum(e.*e)
% % % RN = [1,zeros(N,)]
% % 
% % R = [R0; zeros(length(t),1)]';
% % 
% % for i = 1:length(t)-1
% %     s = 0;
% %     for j = i+1:length(t)
% %             s = s+e(j)*e(j-i);
% %     end
% %     R(i) = s;
% % end
% % 
% % RN = [1;zeros(length(t),1)]';
% % for i = 1:length(t)
% %     RN(i) = R(i)/R0;
% % end
% % n=0;
% % for i=1:length(t)
% %     if abs(RN(i))<=2.17/sqrt(length(t))
% %         n=i;
% %         break;
% %     end
% % end

%% met neparametrica cu regresia liniara
figure
plot(t,[u*200,w,y/2]); title('Semnale'); xlabel('t(sec)'); 
ylabel('u(%), w(rad/s), y(impulsuri)'); grid
tm = t(2769)-t(2717);
wst = mean(w(3081:3246));
ust = mean(u(3081:3246));
% w0 = mean(w(2482:2649));
% u0 = mean(u(2482:2649));
K = wst/ust
i5 = 2778;
i6 = 2918;
tk = t(i5:i6);
xk = log(abs(wst-w(i5:i6)));

Areg = [sum(tk.^2) sum(tk);
        sum(tk) length(tk)]
Breg = [sum(tk*xk) ; sum(xk)];
sol = Areg \ Breg;
T = -1/sol(1)
N = round(tm/(t(2)-t(1)));
uN = [u(1)*ones(N,1); u(1:end-N)];

hold on

A = -1/T;
B = K/T;
C = 1;
D = 0;

wsim = lsim(A,B,C,D,uN,t,w(1));
plot(t,wsim,'g')
legend('u(%)', 'w(rad/s)', 'y(impulsuri)','w dupa identificare')

figure
plot(tk,xk), title('Dreapta de regresie'), xlabel('t_k'); ylabel('x_k'); grid

emin = norm(w-wsim)/norm(w-mean(w))*100

%% Autocorelatia erorilor pentru id. pe baza raspunsului la treapta
j1 = 2482;
j2 = 2649;
j3 = 3081;
j4 = 3286;
e = w(j1:j4)-wsim(j1:j4); %erorile dintre viteza masurata si cea simulata, doar pe portiunea cu treapta
N = length(w(j1:j4));  % = 805
R0 = sum(e.^2)/N; % eroarea medie patratica
n = 0;
R1 = sum([e(1:(N-1))]'* e(2:N))/N;
RN1 = R1/R0;   % 0.9641
R2 = sum([e(1:(N-2))]'* e(3:N))/N;
RN2 = R2/R0     % 0.9298
%fac un for ca sa vad cate valori nu indeplinesc conditia
for i = 1:N
    R(i) = sum([e(1:N-i)]' * e((i+1):N))/N;
    RN(i) = R(i)/R0;

     if abs(RN(i)) <= 2.17/sqrt(N)  % 2.17/sqrt(N) = 0.0765
         n = i; % n = 249
         break;
     end
end
% suma de la t=0 la N din e(t)*e(t-1)
%verificare
R248 = sum([e(1:(N-248))]' *e(249:N))/N;
RN248 = R248/R0; % = 0.0770 > 0.0765
R249 = sum([e(1:(N-249))]' *e(250:N))/N;
RN249= R249/R0; % = 0.0753 < 0.0765
R250 = sum([e(1:(N-250))]' *e(251:N))/N;
RN250 = R250/R0 % = 0.0736 < 0.0765
%concluzie: primele 248 de masuratori de pe treapta nu se afla in banda de
%incredere a autocorelatiei


%% Indici
% indici pentru identificare
i1 = 3458;
i2 = 5186;

% indici pentru validare
i3 = 6764;
i4 = 8313;

idx_id = (i1:i2);
idx_vd = (i3:i4);

figure
sgtitle('Datele pentru identificare u->w');
subplot(211); plot(t(idx_id),u(idx_id)); title('Factor de umplere'); xlabel('t(sec)'); ylabel('u(%)');
subplot(212); plot(t(idx_id),w(idx_id)); title('Viteza unghiulara'); xlabel('t(sec)'); ylabel('w(rad/s)');

figure
sgtitle('Datele pentru identificare w->y');
subplot(211); plot(t(idx_id),w(idx_id)); title('Viteza unghiulara'); xlabel('t(sec)'); ylabel('w(rad/s)');
subplot(212); plot(t(idx_id),y(idx_id)); title('Pozitia unghiulara'); xlabel('t(sec)'); ylabel('y(impulsuri)');

figure
sgtitle('Datele pentru validare u->w');
subplot(211); plot(t(idx_vd),u(idx_vd)); title('Factor de umplere'); xlabel('t(sec)'); ylabel('u(%)');
subplot(212); plot(t(idx_vd),w(idx_vd)); title('Viteza unghiulara'); xlabel('t(sec)'); ylabel('w(rad/s)');

figure
sgtitle('Datele pentru validare w->y');
subplot(211); plot(t(idx_vd),w(idx_vd)); title('Viteza unghiulara'); xlabel('t(sec)'); ylabel('w(rad/s)');
subplot(212); plot(t(idx_vd),y(idx_vd)); title('Pozitia unghiulara'); xlabel('t(sec)'); ylabel('y(impulsuri)');

%% viteza cu verificare prin autocorelatie cu ARX

dt = t(2);

dw_id = iddata(w(idx_id,1),u(idx_id,1),dt);
dw_vd = iddata(w(idx_vd,1),u(idx_vd,1),dt);
% dw_gen = iddata(w(1:end,1),u(1:end,1),dt);
figure

Hw_arx = arx(dw_id,[1,1,9])
compare(dw_vd,Hw_arx), shg;
figure
resid(dw_vd,Hw_arx),shg;

Hd1 = tf(Hw_arx.B,Hw_arx.A,dt,'variable','z^-1')
Hc1 = d2c(Hw_arx,'zoh')
Mc1 = tf(Hc1)

sys_w = ss(Hc1);
figure
plot(t,u);
hold on
w_sim = lsim(sys_w.A,sys_w.B,sys_w.C,sys_w.D,u,t,w(1)./sys_w.C);
plot(t,[w,w_sim]),grid
e_Mpn = norm(w-w_sim)/norm(w-mean(w))*100

%% viteza cu verificare prin autocorelatie cu ARMAX

dt = t(2);

dw_id = iddata(w(idx_id,1),u(idx_id,1),dt);
dw_vd = iddata(w(idx_vd,1),u(idx_vd,1),dt);
% dw_gen = iddata(w(1:end,1),u(1:end,1),dt);
figure

Hw_armax = armax(dw_id,[1,1,1,8])
compare(dw_vd,Hw_armax), shg;
figure
resid(dw_vd,Hw_armax),shg;

Hd2 = tf(Hw_armax.B,Hw_armax.A,dt,'variable','z^-1')
Hc2 = d2c(Hw_armax,'zoh')
Mc2 = tf(Hc2)

sys_w = ss(Hc2);
figure
% plot(t,u);
hold on
w_sim = lsim(sys_w.A,sys_w.B,sys_w.C,sys_w.D,u,t,w(1)/sys_w.C);
plot(t,[w,w_sim]),grid
title('Viteza unghiulara')
xlabel('t(sec)'); ylabel('w(rad/s)')
legend('\omega', '\omega_s_i_m');
e_Mpn = norm(w-w_sim)/norm(w-mean(w))*100

%% viteza cu verificare prin intercorelatie cu IV

dt = t(2);

dw_id = iddata(w(idx_id,1),u(idx_id,1),dt);
dw_vd = iddata(w(idx_vd,1),u(idx_vd,1),dt);
% dw_gen = iddata(w(1:end,1),u(1:end,1),dt);

figure

Hw_iv = iv4(dw_id,[1,1,11])
compare(dw_vd,Hw_iv), shg;
figure
resid(dw_vd,Hw_iv),shg;

Hd3 = tf(Hw_iv.B,Hw_iv.A,dt,'variable','z^-1')
Hc3 = d2c(Hw_iv,'zoh')
Mc3 = tf(Hc3)


sys_w = ss(Hc3);
figure
% plot(t,u);
hold on
w_sim = lsim(sys_w.A,sys_w.B,sys_w.C,sys_w.D,u,t,w(1)/sys_w.C);
plot(t,[w,w_sim]),grid
title('Viteza unghiulara')
xlabel('t(sec)'); ylabel('w(rad/s)')
legend('\omega', '\omega_s_i_m');
e_Mpn = norm(w-w_sim)/norm(w-mean(w))*100

%% viteza cu verificare prin intercorelatie cu OE

dt = t(2);

dw_id = iddata(w(idx_id,1),u(idx_id,1),dt);
dw_vd = iddata(w(idx_vd,1),u(idx_vd,1),dt);
% dw_gen = iddata(w(1:end,1),u(1:end,1),dt);

figure

Hw_oe = oe(dw_id,[1,1,1])
compare(dw_vd,Hw_oe), shg;
figure
resid(dw_vd,Hw_oe),shg;

Hd4 = tf(Hw_oe.B,Hw_oe.F,dt,'variable','z^-1')
Hc4 = d2c(Hd4,'zoh')

sys_w = ss(Hc4);
figure
plot(t,u);
hold on
w_sim = lsim(sys_w.A,sys_w.B,sys_w.C,sys_w.D,u,t,w(1)/sys_w.C)
plot(t,[w,w_sim]),grid
e_Mpn = norm(w-w_sim)/norm(w-mean(w))*100

%% pozitie prin autocorelatie unde nici ARX nici cu ARMAX nu trece testul

dy_id = iddata(y(idx_id,1),w(idx_id,1),dt);
dy_vd = iddata(y(idx_vd,1),w(idx_vd,1),dt);
% dy_gen = iddata(y(1:end,1),w(1:end,1),dt);

figure

Hy_arx = arx(dy_id,[1,1,0]);
compare(dy_vd,Hy_arx), shg;
figure
resid(dy_vd,Hy_arx),shg;

figure
Hy_armax = armax(dy_id,[1,1,1,0])
compare(dy_vd,Hy_armax), shg;
figure
resid(dy_vd,Hy_armax),shg;

Hd10 = tf(Hy_armax.B,Hy_armax.A,dt,'variable','z^-1')
Hc10 = tf(Hy_armax.B/dt,[1 0])

%% decimare date
N = 12;
u1 = u(i1:N:i2);
w1 = w(i1:N:i2);
y1 = y(i1:N:i2);
u1_v = u(i3:N:i4);
w1_v = w(i3:N:i4);
y1_v = y(i3:N:i4);
% u1_gen = u(1:N:end);
% w1_gen = w(1:N:end);
% y1_gen = y(1:N:end);

%% pozitie prin autocorelatie cu ARX decimat
dy1_id = iddata(y1,w1,N*dt);
dy1_vd = iddata(y1_v,w1_v,N*dt);
% dy1_gen = iddata(y1_gen,w1_gen,N*dt);
figure

Hy1_arx = arx(dy1_id,[1,1,0])
compare(dy1_vd,Hy1_arx), shg;
figure
resid(dy1_vd,Hy1_arx),shg;

Hd5 = tf(Hy1_arx.B,Hy1_arx.A,N*dt,'variable','z^-1')
Hc5_1 = d2c(Hd5,'zoh')

Hc5 = tf(Hy1_arx.B(end)/(N*dt),[1 0])
sys_y = ss(Hc5);
figure
plot(t,w);
hold on
y_sim = lsim(sys_y.A,sys_y.B,sys_y.C,sys_y.D,w,t,y(1)/sys_y.C)
plot(t,[y,y_sim]),grid
e_Mpn = norm(y-y_sim)/norm(y-mean(y))*100
%% pozitie prin autocorelatie cu ARMAX decimat
dy1_id = iddata(y1,w1,N*dt);
dy1_vd = iddata(y1_v,w1_v,N*dt);
% dy1_gen = iddata(y1_gen,w1_gen,N*dt);
figure

Hy1_armax = armax(dy1_id,[1,1,1,0])
compare(dy1_vd,Hy1_armax), shg;
figure
resid(dy1_vd,Hy1_armax),shg;

Hd6 = tf(Hy1_armax.B,Hy1_armax.A,N*dt,'variable','z^-1')
Hc6_1 = d2c(Hd6,'zoh')

Hc6 = tf(Hy1_armax.B/N/dt,[1 0])
sys_y = ss(Hc6);
figure
% plot(t,w);
hold on
y_sim = lsim(sys_y.A,sys_y.B,sys_y.C,sys_y.D,w,t,y(1)/sys_y.C);
plot(t,[y,y_sim]),grid
title('Pozitia unghiulara')
xlabel('t(sec)'); ylabel('y(impulsuri)')
legend('\theta', '\theta_s_i_m');
e_Mpn = norm(y-y_sim)/norm(y-mean(y))*100
%% pozitie prin intercorelatie cu OE
dy_id = iddata(y(idx_id,1),w(idx_id,1),dt);
dy_vd = iddata(y(idx_vd,1),w(idx_vd,1),dt);
% dy_gen = iddata(y(1:end,1),w(1:end,1),dt);
figure

Hy1_oe = oe(dy_id,[1,1,0])
compare(dy_vd,Hy1_oe), shg;
figure
resid(dy_vd,Hy1_oe),shg;

Hd7 = tf(Hy1_oe.B,Hy1_oe.F,dt,'variable','z^-1')
Hc7_1 = d2c(Hd7,'zoh')

Hc7 = tf(Hy1_oe.B/dt,[1 0])
sys_y = ss(Hc7);
figure
% plot(t,w);
hold on
y_sim = lsim(sys_y.A,sys_y.B,sys_y.C,sys_y.D,w,t,y(1)/sys_y.C);
plot(t,[y,y_sim]),grid
title('Pozitia unghiulara')
xlabel('t(sec)'); ylabel('y(impulsuri)')
legend('\theta', '\theta_s_i_m');
e_Mpn = norm(y-y_sim)/norm(y-mean(y))*100

%% pozitie prin intercorelatie cu IV
dy_id = iddata(y(idx_id,1),w(idx_id,1),dt);
dy_vd = iddata(y(idx_vd,1),w(idx_vd,1),dt);
% dy_gen = iddata(y(1:end,1),w(1:end,1),dt);
figure

Hy1_iv = iv4(dy_id,[1,1,0])
compare(dy_vd,Hy1_iv), shg;
figure
resid(dy_vd,Hy1_iv),shg;

Hd8 = tf(Hy1_iv.B,Hy1_iv.A,dt,'variable','z^-1')
Hc8_1 = d2c(Hd8,'zoh')
Hc8 = tf(Hy1_iv.B(end)/dt,[1 0])
sys_y = ss(Hc8);
figure
plot(t,w);
hold on
y_sim = lsim(sys_y.A,sys_y.B,sys_y.C,sys_y.D,w,t,y(1)/sys_y.C)
plot(t,[y,y_sim]),grid
e_Mpn = norm(y-y_sim)/norm(y-mean(y))*100
%% pozitie prin intercorelatie cu IV decimat
dy1_id = iddata(y1,w1,N*dt);
dy1_vd = iddata(y1_v,w1_v,N*dt);
% dy1_gen = iddata(y1_gen,w1_gen,N*dt);
figure

Hy1_iv4 = iv4(dy1_id,[1,1,0])
compare(dy1_vd,Hy1_iv4),shg;
figure
resid(dy1_vd,Hy1_iv4),shg;

Hd9 = tf(Hy1_iv4.B,Hy1_iv4.A,N*dt,'variable','z^-1')
Hc9_1 = d2c(Hd9,'zoh')
Hc9 = tf(Hy1_iv4.B(end)/N/dt,[1 0])
sys_y = ss(Hc9);
figure
plot(t,w);
hold on
y_sim = lsim(sys_y.A,sys_y.B,sys_y.C,sys_y.D,w,t,y(1)/sys_y.C)
plot(t,[y,y_sim]),grid
e_Mpn = norm(y-y_sim)/norm(y-mean(y))*100