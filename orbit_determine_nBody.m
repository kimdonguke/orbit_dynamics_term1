% 우주궤도역학 term project#1

% 상수 선언
G = 6.67430e-20;
labels = {'Earth', 'Moon', 'Satellite', 'Sun'};
mu=0;

% masses
m_earth = 5.972e24;
m_moon  = 7.342e22;
m_sat   = 1000;       % small satellite(임의질량)
m_sun   = 1.989e30;

% Initial positions
r_earth = [0; 0; 0];
r_moon  = [384400; 0; 0];
r_sat   = [10000; 0; 0];
r_sun   = [1.496e8; 0; 0];

% Initial velocities
v_earth = [0; 29.78; 0];
v_moon  = v_earth + [0; 1.022; 0];
v_sat   = v_earth + [0; 7; 0];
v_sun   = [0; 0; 0]; % helio centric inertial frame

% 위치,속도,질량 배열
state0 = [r_earth; r_moon; r_sat; r_sun; v_earth; v_moon; v_sat; v_sun];
masses = [m_earth, m_moon, m_sat, m_sun]; %지구 달 위성 태양

%% orbit determine 시작 %%
% 1. Time Sampling 
n=10; % n 일
T = n*86160; 
time_sampling = n*50; % 궤도 변화를 보는 것이기 때문에 곱해지는 인자를 크게 할 필요가 없음
t_eval=linspace(0, T, time_sampling);

% 2. time propagation
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, y] = ode45(@(t,y) computeNBody(t,y,masses,G),t_eval, state0, opts);

% 3. centric frame이 바뀌는지 관찰, 바뀐다면 출력
frame_switch(y,masses,G,labels);


%궤도 결정 함수
function orbit_type = orbit_determine(r, v, mu)
    h = cross(r, v);
    r_norm = norm(r);
    e_vec = (1/mu) * (cross(v, h) - mu * (r / r_norm));
    e = norm(e_vec);

    if abs(e) < 1e-6
        orbit_type = 'Circular';
    elseif e < 1
        orbit_type = 'Elliptic';
    elseif abs(e - 1) < 1e-6
        orbit_type = 'Parabolic';
    else
        orbit_type = 'Hyperbolic';
    end
end

% time propagation 함수 선언
function dydt = computeNBody(t, y, masses, G)
    N = length(masses); %천체 갯수
    r = reshape(y(1:3*N), 3, N); 
    v = reshape(y(3*N+1:end), 3, N);
    a = zeros(3, N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                % F=ma, F=GMm*r/|r|^3
                diff = r(:,j) - r(:,i); %거리차
                dist3 = norm(diff)^3 + 1e-9; % 0방지(syntax error)
                a(:,i) = a(:,i) + G * masses(j) * diff / dist3; %
                %fprintf("%d",i); %debug
            end
        end
    end
    dydt = [v(:); a(:)]; %dydt 반환[Vn, An]
end

% centric frame 결정 함수
% function [r_rel, v_rel, primary_idx] = select_frame(r_sat, v_sat, r_state, v_state, masses, G)
%     accels = G .* masses ./ vecnorm(r_state - r_sat, 2, 1).^2; %상대거리 슬라이싱 후 a=GM/r^2 계산
%     [~, primary_idx] = max(accels); % 가속도가 가장 큰 원소의 인덱스 추출
%     r_rel = r_sat - r_state(:, primary_idx);
%     v_rel = v_sat - v_state(:, primary_idx);
% end
% 폐기;;

% centric frame이 바뀌는지 안 바뀌는지 검사하는 함수
function frame_switch(y, masses, G, labels)
    % masses: 각 천체의 질량 배열
    % G: 중력상수
    % labels: 각 천체 이름 문자열 배열

    N = length(masses);      % 천체 수
    steps = size(y, 1);      % 시간 스텝 수

    % 위성 인덱스
    idx_sat = 3;

    % 이전 중심체 인덱스 초기화
    previous_idx = -1;

    r_full = reshape(y(:, 1:3*N)', 3, N, steps);
    v_full = reshape(y(:, 3*N+1:end)', 3, N, steps);

    for k = 1:steps
        % 위치/속도 재구성
        r_all = r_full(:, :, k); %위치 행렬 참조
        v_all = v_full(:, :, k); %속도 행렬 참조

        r_sat = r_all(:, idx_sat);
        v_sat = v_all(:, idx_sat);

        % 모든 천체와의 중력 가속도 크기 계산
        accels = G .* masses ./ vecnorm(r_all - r_sat, 2, 1).^2;
        accels(idx_sat) = -inf; %자기 자신을 모든 행에서 지움

        [~, primary_idx] = max(accels); % 값은 무시하고 인덱스만

        % 중심체가 바뀌는 순간에만 출력
        if primary_idx ~= previous_idx
            % 상대 위치/속도 계산
            r_rel = r_sat - r_all(:, primary_idx);
            v_rel = v_sat - v_all(:, primary_idx);
            mu = G * masses(primary_idx);

            % 궤도 형태 판별
            orbit = orbit_determine(r_rel, v_rel, mu);
            fprintf('Step %d | Frame: %s | Orbit: %s\n', k, labels{primary_idx}, orbit);

            previous_idx = primary_idx;
        end
    end
end