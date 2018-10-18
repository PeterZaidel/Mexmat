clear all
close all
clc

trj_data=load('Data/trj.dat');
imu_data=load('Data/imu.dat');

trj_time = trj_data(:,1);
trj_geo = trj_data(:, 2:4);

trj_geo(:,1) = deg2rad(trj_geo(:,1));
trj_geo(:,2) = deg2rad(trj_geo(:,2));

trj_v = trj_data(:, 5:7);
trj_eulers = deg2rad(trj_data(:, 8:10));


imu_time = imu_data(:, 1);
imu_omega_z = deg2rad(imu_data(:, 2:4));
imu_f_z = imu_data(:, 5:7);

dt = imu_time(2) - imu_time(1);

geo_0 = trj_geo(2,:);

u_z=mean(imu_omega_z(1:100/dt,:));
f_zm=mean(imu_f_z(1:100/dt,:));
gx=-[0; 0; norm(f_zm)];

L=zeros(3);
L(:,3)=f_zm/norm(f_zm);
L(:,3) = L(:, 3)/ norm(L(:, 3));
L(:,2)=(u_z'-L(:,3)*norm(u_z)*sin( geo_0(1) ))/(norm(u_z)*cos(geo_0(1)));
L(:,2) = L(:, 2)/ norm(L(:, 2));
 
% ортогонализация грамма-шмидта
l32 = L(:,3)' * L(:,2) /(L(:,3)' * L(:,3));
L(:,2) = L(:,2) - l32 * L(:,3);
L(:,2) = L(:, 2)/ norm(L(:, 2));

L(:,1)=cross(L(:,2),L(:,3));


a=6378137;
b=6356752;
R=1/((1/a^2+3/b^2)^(1/2));
e2=6.6943799901413*10^(-3);
g_e = 9.78030;
uzemli = 2*pi/86164.090530833;

Az=L;
Ax=eye(3);
V= trj_v(1, :)';

geo = zeros(89701, 3);
V_d = zeros(89701, 3);
eulers = zeros(89701, 3);

geo(1,:) = geo_0;
V_d(1, :) = V;

eulers(1, :) = trj_eulers(1, :);

for i=1:size(imu_time)-1
    %i
    dt = imu_time(i+1) - imu_time(i);
    trj_idx = find(trj_time == floor( imu_time(i) ));
    
    phi_i = geo(i, 1);
    lam_i = geo(i, 2);
    h_i = trj_geo(trj_idx, 3);
    
    % тут была ошибка - Козлов сказал, что нужно добавть 2e-3
    mod_g = g_e * (1 - 2 * h_i/a + 3/4 *e2 * (sin(phi_i))^2 ) + 2e-3;
    gx = [0;0;-mod_g];
    
    Re = a/sqrt(1-e2*(sin(phi_i))^2);
    Rn = a*(1-e2)/(sqrt(1-e2*sin(phi_i)*sin(phi_i)))^3;

    Az = integratePoisson(Az, imu_omega_z(i,:), dt); 
    omega_x = [ -V(2)/(Rn+h_i), V(1)/(Re+h_i), V(1)*tan(phi_i)/(Re+h_i)];
    ux = [ 0, norm(u_z)*cos(phi_i), norm(u_z)*sin(phi_i) ];
    Ax= integratePoisson(Ax,omega_x + ux,dt);
    L=Az*Ax';
    %L'*L
    
    % нормализация матрицы ориентации - чисто для дополнительной проверки
    L(:,3) = L(:, 3)/ norm(L(:, 3));
    L(:,2) = L(:, 2)/ norm(L(:, 2));
    % ортогонализация грамма-шмидта
    l32 = L(:,3)' * L(:,2) /(L(:,3)' * L(:,3));
    L(:,2) = L(:,2) - l32 * L(:,3);
    L(:,1)=cross(L(:,2),L(:,3));
    
    lam_real = trj_geo(trj_idx, 2);
    
    phi_next = phi_i + ( V(2) / ( Rn + h_i ) )*dt;
    lam_next = lam_i + ( V(1) /(  (Re+ h_i)*cos(phi_i) ) )*dt;
    h_next = h_i + V(3)*dt; 
    
    geo_next = [phi_next, lam_next, h_next];
    geo(i+1,:) = geo_next;

    kos=koso(omega_x+2*ux);
    
    prev_V = V;
    
    trj_V = trj_v(trj_idx, :);
    V=V+(kos*V + gx + L'*imu_f_z(i,:)')*dt;
    V_d(i+1, :) = V;
    
    %gx + L'*imu_f_z(i,:)'
    
    %[trj_V; V']
    
    cur_eulers = L2eulers(L);
    
    real_eulers = trj_eulers(trj_idx, :);
    eulers(i+1, :) = L2eulers(L);
end

% отрисовка координат: широта, долгота, высота
% plot geo
figure()
plot(imu_time,geo(:,1),'r', 'DisplayName' ,'my phi')
hold
plot(trj_time,trj_geo(:,1),'g', 'DisplayName' ,'original phi')
legend
title('geo phi')

figure()
plot(imu_time,geo(:,2),'r', 'DisplayName' ,'my lam')
hold
plot(trj_time,trj_geo(:,2),'g', 'DisplayName' ,'original lam')
legend
title('geo lam')

figure()
plot(imu_time,geo(:,3),'r', 'DisplayName' ,'my h')
hold
plot(trj_time,trj_geo(:,3),'g', 'DisplayName' ,'original h')
legend
title('geo h')


% отрисовка скоростей V1, V2, V3
% plot V
figure()
plot(imu_time,V_d(:,1),'r', 'DisplayName' ,'my v1')
hold
plot(trj_time,trj_v(:,1),'g', 'DisplayName' ,'original v1')
legend
title('V1')

figure()
plot(imu_time,V_d(:,2),'r', 'DisplayName' ,'my v2')
hold
plot(trj_time,trj_v(:,2),'g', 'DisplayName' ,'original v2')
legend
title('V2')

figure()
plot(imu_time,V_d(:,3),'r', 'DisplayName' ,'my v3')
hold
plot(trj_time,trj_v(:,3),'g', 'DisplayName' ,'original v3')
legend
title('V3')

% отрисовка углов
% plot eulers
figure()
plot(imu_time,eulers(:,1),'r', 'DisplayName' ,'my csi')
hold
plot(trj_time,trj_eulers(:,1),'g', 'DisplayName' ,'original csi')
legend
title('eulers csi')

figure()
plot(imu_time,eulers(:,2),'r', 'DisplayName' ,'my gam')
hold
plot(trj_time,trj_eulers(:,2),'g', 'DisplayName' ,'original gam')
legend
title('eulers gam')

figure()
plot(imu_time,eulers(:,3),'r', 'DisplayName' ,'my th')
hold
plot(trj_time,trj_eulers(:,3),'g', 'DisplayName' ,'original th')
legend
title('eulers th')


% ошибки в точках измерения реальных данных

delta_geo = zeros(1, 3);
delta_v = zeros(1, 3);
delta_eulers = zeros(1, 3);

% примерный радиус земли для перевода ошибок по широте и долготе в метры
R_zemli = 6371 * 10^3;
for i=1:size(trj_time)
    imu_idx = find(imu_time == trj_time(i));
    if isempty(imu_idx)
        continue
    end
    err = trj_geo(i, :) - geo(imu_idx, :);
    err(1) = err(1) * R_zemli;
    err(2) = err(2) * R_zemli * cos(trj_geo(i, 1));
    delta_geo(i, :) = err;
    
    delta_v(i, :) = trj_v(i,:) - V_d(imu_idx, :);
    delta_eulers(i, :) = trj_eulers(i, :) - eulers(imu_idx, :);
end

%plot geo errors
figure()
subplot(3,1,1)
plot(delta_geo(:,1),'r', 'DisplayName' ,'phi error')
subplot(3,1,2)
plot(delta_geo(:,2),'g', 'DisplayName' ,'lam error')
subplot(3,1,3)
plot(delta_geo(:,3),'g', 'DisplayName' ,'h error')
legend
title('geo errors')


%plot V errors
figure()
subplot(3,1,1)
plot(delta_v(:,1),'r', 'DisplayName' ,'v1 error')
subplot(3,1,2)
plot(delta_v(:,2),'g', 'DisplayName' ,'v2 error')
subplot(3,1,3)
plot(delta_v(:,3),'g', 'DisplayName' ,'v3 error')
legend
title('V errors')

%plot eulers errors
figure()
subplot(3,1,1)
plot(delta_eulers(:,1),'r', 'DisplayName' ,'delta_eulers 1')
subplot(3,1,2)
plot(delta_eulers(:,2),'g', 'DisplayName' ,'delta_eulers 2')
subplot(3,1,3)
plot(delta_eulers(:,3),'g', 'DisplayName' ,'delta_eulers 3')
legend
title('eulers errors')
