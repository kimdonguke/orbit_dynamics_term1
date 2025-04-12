% 우주궤도역학 term project#1 (i=0, RAAN=0, w=0)

%변수 선언(지구 중력 상수, 반장축, 이심률)
mu = 398600; %지구중력상수 뮤(km^3/s^2)
a = 20000;    % 반장축 a(km)
e = 0.3;      % 이심률
nu = 0;       % 진근점각(True Anomaly)


p = a * (1 - e^2);                                  % p 계산
r = p / (1 + e * cos(nu));                          % r계산
r_vec = [r * cos(nu); r * sin(nu); 0];              % r벡터 변환
v_vec = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];    % v벡터 변환

state0 = [r_vec; v_vec];       % 초기 상태 벡터
T = 2 * pi * sqrt(a^3 / mu);  % 주기 계산
tspan = [0, 2*T];        % 주기를 두 번 도는 궤도 propagation 범위
t_eval = linspace(tspan(1), tspan(2), 1000); % 시간 분할(0~궤도 2바퀴, 1000개 구간으로 나눠서)

% two_body 함수 정의
two_body = @(t, y) [y(4:6); -mu * y(1:3) / norm(y(1:3))^3];

% 오차 보정
opts=odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, state] = ode45(two_body, t_eval, state0,opts);


figure;
plot3(state(:,1), state(:,2), state(:,3));
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbit Propagation in ECI Frame');
grid on; axis equal;