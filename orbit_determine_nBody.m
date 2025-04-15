% 우주궤도역학 term project#1

% 변수 선언
G = 6.67430e-20;

% masses
m_earth = 5.972e24;
m_moon  = 7.342e22;
m_sat   = 1000;       % small satellite(임의질량)
m_sun   = 1.989e30;

% Initial positions
r_earth = [0; 0; 0];
r_moon  = [384400; 0; 0];
r_sat   = [50000; 0; 0];
r_sun   = [1.496e8; 0; 0];

% Initial velocities
v_earth = [0; 29.78; 0];
v_moon  = v_earth + [0; 1.022; 0];
v_sat   = v_earth + [0; 10; 0];
v_sun   = [0; 0; 0]; % helio centric inertial frame


masses = [m_earth, m_moon, m_sat, m_sun]; %지구 달 위성 태양


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

% centric frame 결정 함수
function [r_rel, v_rel, primary_idx] = select_frame(r_sat, v_sat, r_state, v_state, masses, G)
    accels = G .* masses ./ vecnorm(r_state - r_sat, 2, 1).^2; %상대거리 슬라이싱 후 a=GM/r^2 계산
    [~, primary_idx] = max(accels); % 가속도가 가장 큰 원소의 인덱스 추출
    r_rel = r_sat - r_state(:, primary_idx);
    v_rel = v_sat - v_state(:, primary_idx);
end