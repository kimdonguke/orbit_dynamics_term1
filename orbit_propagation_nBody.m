% % 우주궤도역학 term project#1 n-body propagation

% Constants
G = 6.67430e-20; % 중력상수
mu = 398600;     % 지구중력상수 
earth_radius = 6371; %지구 반지름

% Masses
m_earth = 5.972e24;
m_moon  = 7.342e22;
m_sat   = 1000;       % small satellite
m_sun   = 1.989e30;

% Initial positions in ECI
r_earth = [0; 0; 0];
r_moon  = [384400; 0; 0];
r_sat   = [7000; 0; 0];
r_sun   = [1.496e8; 0; 0]; % 1 AU from Earth

% Initial velocities
v_earth = [0; 29.78; 0];
v_moon  = v_earth + [0; 1.022; 0];
v_sat   = v_earth + [0; 10; 0]; %7.8
v_sun   = [0; 0; 0]; % assume stationary in short sim

% 4-body initial state vector
state0 = [r_earth; r_moon; r_sat; r_sun; v_earth; v_moon; v_sat; v_sun];
masses = [m_earth, m_moon, m_sat, m_sun];

% Time span
n=1200; % n 일
T = n*86160; 

% 적분
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, y] = ode45(@(t,y) computeNBody(t,y,masses,G), [0 T], state0, opts);

heliocentric_plot(t,y);



% 함수 선언
function dydt = computeNBody(t, y, masses, G)
    N = length(masses);
    r = reshape(y(1:3*N), 3, N);
    v = reshape(y(3*N+1:end), 3, N);
    a = zeros(3, N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                diff = r(:,j) - r(:,i);
                dist3 = norm(diff)^3 + 1e-9;
                a(:,i) = a(:,i) + G * masses(j) * diff / dist3;
            end
        end
    end
    dydt = [v(:); a(:)];
end

function heliocentric_plot(t, y)
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
    title('Heliocentric Frame: 4-body Trajectories'); legend; grid on; axis equal;
end

function eci_plot(t,y)
    earth_radius = 6371;
    %Earth-Centered Inertial (ECI) 시점 시각화
    r_moon_rel = r_moon - r_earth;   % 달 - 지구
    r_sat_rel  = r_sat  - r_earth;   % 위성 - 지구
    
    % 지구 반지름 표시용 원
    theta = linspace(0, 2*pi, 300);
    earth_circle_x = earth_radius * cos(theta);
    earth_circle_y = earth_radius * sin(theta);
    
    plot(r_moon_rel(1,:), r_moon_rel(2,:), 'g', 'DisplayName', 'Moon (rel. Earth)'); hold on;
    plot(r_sat_rel(1,:)*25, r_sat_rel(2,:)*25, 'r', 'DisplayName', 'Satellite (x25)');
    plot(earth_circle_x, earth_circle_y, 'b', 'LineWidth', 1.5, 'DisplayName', 'Earth Radius');
    scatter(0, 0, 2, 'yellow', 'filled');  % 지구 중심점
    xlabel('X [km]'); ylabel('Y [km]');
    title('ECI Frame: Relative Trajectories around Earth');
    axis equal; grid on; legend; hold off;
end


