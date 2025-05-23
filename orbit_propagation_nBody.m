% % 우주궤도역학 term project#1 n-body propagation

clear;
clc;

% Constants
G = 6.67430e-20; % 중력상수
mu = 398600;     % 지구중력상수 
earth_radius = 6371; %지구 반지름

% 변수선언
% masses
m_earth = 5.972e24;
m_moon  = 7.342e22;
m_sat   = 1000;       % small satellite(임의질량)
m_sun   = 1.989e30;

% Initial positions
r_earth = [0; 0; 0];
r_moon  = [384400; 0; 0];
%r_sat   = [310000; 0; 0]; %swingby distance
r_sat   = [100000; 0; 0];
r_sun   = [1.496e8; 0; 0];

% Initial velocities
v_earth = [0; 29.78; 0];
v_moon  = v_earth + [0; 1.022; 0];
v_sat   = v_earth + [0; 10; 0];
% cv_sat   =sqrt(mu/norm(r_sat)); %circula orbit으로 계산하고 싶을 때 사용 할 수 있는 식
% v_sat   = v_earth + [0; cv_sat; 0];
v_sun   = [0; 0; 0]; % helio centric inertial frame

% 4-body initial state vector
state0 = [r_earth; r_moon; r_sat; r_sun; v_earth; v_moon; v_sat; v_sun];
masses = [m_earth, m_moon, m_sat, m_sun];

% Time span, sampling
n=180; % n 일
T = n*86160; 
time_sampling = n*200; %하루를 200번으로 쪼개서 샘플링, 일수가 작다면 잘게 쪼개도 좋음
t_eval=linspace(0, T, time_sampling);

% 적분
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, y] = ode45(@(t,y) computeNBody(t,y,masses,G),t_eval, state0, opts);
[t_rel,y_rel] = ode45(@(t_rel,y_rel) rel_computeNBody(t_rel,y_rel,masses,G),t_eval, state0, opts);

heliocentric_plot(y); % 태양 중심 좌표계
%eci_plot(y); % 지구중심좌표계

% 함수 선언
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
            end
        end
    end
    dydt = [v(:); a(:)]; %dydt 반환[Vn, An]
end

function dydt = rel_computeNBody(t, y, masses, G)
    N = length(masses);              % 천체 수
    r = reshape(y(1:3*N), 3, N);     % 위치 행렬
    v = reshape(y(3*N+1:end), 3, N); % 속도 행렬
    a = zeros(3, N);                 % 가속도 행렬 초기화
    
    c = 299792.458; % 빛의 속도 [km/s]

    for i = 1:N
        for j = 1:N
            if i ~= j
                diff = r(:,j) - r(:,i);   % 거리 벡터
                dist = norm(diff);       % 거리 크기
                dist3 = dist^3 + 1e-9;    % 0 방지용

                % 뉴턴 역학 중력 가속도
                a(:,i) = a(:,i) + G * masses(j) * diff / dist3;

                % 일반 상대성 이론 보정항 (근사)
                mu_j = G * masses(j);       % 중심천체 j의 중력상수
                v_i = v(:,i);               % 관측천체 i의 속도
                gr_term1 = (4 * mu_j / dist - dot(v_i, v_i)) * diff;
                gr_term2 = 4 * dot(diff, v_i) * v_i;
                a_gr = (mu_j / (c^2 * dist^3)) * (gr_term1 + gr_term2);

                % GR 보정 추가
                a(:,i) = a(:,i) + a_gr;
            end
        end
    end
    dydt = [v(:); a(:)]; % 상태 벡터 도함수 반환
end

function heliocentric_plot(y)
    % 상수
    N = 4;  % 지구, 달, 위성, 태양(n-body 갯수)
    % 위치 벡터 추출 및 재배열
    r = y(:, 1:3*N);               % 위치 벡터만 사용
    r = reshape(r', 3, N, []);     % [3 x N x length(t)]
    
    % 천체별 위치 분리
    
    r_earth = squeeze(r(:,1,:));   % 지구
    r_moon  = squeeze(r(:,2,:));   % 달
    r_sat   = squeeze(r(:,3,:));   % 위성
    r_sun   = squeeze(r(:,4,:));   % 태양
    
    %Helio-centric 그래프
    figure;
    nexttile;
    plot(r_earth(1,:), r_earth(2,:), 'b', 'DisplayName', 'Earth'); hold on;
    plot(r_moon(1,:), r_moon(2,:), 'g', 'DisplayName', 'Moon');
    plot(r_sat(1,:), r_sat(2,:), 'r', 'DisplayName', 'Satellite');
    scatter(r_sun(1,1), r_sun(2,1), 80, 'filled');
    xlabel('X [km]'); ylabel('Y [km]');
    title('Heliocentric Frame : n-body time propagation(n=4)');
    legend;
    grid on;
    axis equal;
end

function eci_plot(y)
    % tips : ECI frame을 관찰 할 때에는 지구-위성-달의 관계를 보기 때문에
    % sat_r을 크게 설정하는 것이 좋습니다(r_moon = 384400)
    N = 3;

    r = y(:, 1:3*N);
    r = reshape(r', 3, N, []);
    
    r_earth = squeeze(r(:,1,:));
    r_moon  = squeeze(r(:,2,:));
    r_sat   = squeeze(r(:,3,:));

    % 상대 거리 선언
    r_moon_rel = r_moon - r_earth;   % 달 - 지구
    r_sat_rel  = r_sat  - r_earth;   % 위성 - 지구
    
    plot(r_moon_rel(1,:), r_moon_rel(2,:), 'g', 'DisplayName', 'Moon'); hold on;
    plot(r_sat_rel(1,:), r_sat_rel(2,:), 'r', 'DisplayName', 'Satllite');
    scatter(0, 0, 20, 'yellow', 'filled');  % 지구 질점
    xlabel('X [km]'); ylabel('Y [km]');
    title('ECI Frame : n-body time propagation(n=4)');
    axis equal;
    grid on; 
    legend;
    hold off;
end


