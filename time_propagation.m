% 우주궤도역학 term project#1 (i=0, RAAN=0, w=0)

%변수 선언(지구 중력 상수, 반장축, 이심률)
mu = 398600;  %지구중력상수 뮤(km^3/s^2)
a = 20000;    % 반장축 a(km)
e = 0.8;      % 이심률
nu = 0;       % 진근점각(True Anomaly)


p = a * (1 - e^2);                                  % p 계산
r = p / (1 + e * cos(nu));                          % r계산
r_vec = [r * cos(nu); r * sin(nu); 0];              % r벡터 변환
v_vec = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];    % v벡터 변환

state0 = [r_vec; v_vec];       % 초기 상태 벡터
T = 2 * pi * sqrt(a^3 / mu);   % 주기 계산
n=4;                           % n번
time_sampling=n*700;           % time_sampling
tspan = [0, n*T];        % 주기를 n 번 도는 궤도 propagation 범위
t_eval = linspace(tspan(1), tspan(2), time_sampling); % 시간 분할(0~주기, 각 700개 구간으로 나눠서)
% 주기가 늘어날 수록 time sampling도 커져야함

% two_body 함수 정의
two_body = @(t, y) [y(4:6); -mu * y(1:3) / norm(y(1:3))^3];

% 오차 보정
opts=odeset('RelTol',1e-9,'AbsTol',1e-9);
%적분
[t, state] = ode45(two_body, t_eval, state0, opts);

%그래프용 변수 재정의
x = state(:,1);
y = state(:,2);
vx = state(:,4);
vy = state(:,5);
vz = state(:,6);
v_mag = sqrt(vx.^2 + vy.^2 + vz.^2);
time_min = t / 60;

% 2*2 tiledLayout setting
figure;
tiledlayout(2,2);

% 1번째 tile(propagation of X)
nexttile;
plot(time_min, x, 'b');
title('X position by time');
xlabel('Time (min)');
ylabel('X [km]');
grid on;

% 2번째 tile(propagation of Y)
nexttile;
plot(time_min, y, 'g');
title('Y position by time');
xlabel('Time (min)');
ylabel('Y (km)');
grid on;

% 3번째 tile(propagation of velocity)
nexttile;
plot(time_min, v_mag, 'm');
title('Velocity by time');
xlabel('Time (min)');
ylabel('Velocity (km/s)');
grid on;

% 4번째 tile(orbit)
nexttile;
plot(x, y, 'k');
title('Orbit shape');
xlabel('X (km)');
ylabel('Y (km)');
axis equal;
grid on