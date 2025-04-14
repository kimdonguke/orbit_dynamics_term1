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
v_sat   = v_earth + [0; 7.8; 0];
v_sun   = [0; 0; 0]; % assume stationary in short sim

% 4-body initial state vector
state0 = [r_earth; r_moon; r_sat; r_sun; v_earth; v_moon; v_sat; v_sun];
masses = [m_earth, m_moon, m_sat, m_sun];

% N-body dynamics function
n_body = @(t, y) computeNBody(t, y, masses, G);

% Time span
T = 120*24*3600; % simulate for 1 day

% Integrate n-body system
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, y] = ode45(@(t,y) computeNBody(t,y,masses,G), [0 T], state0, opts);

% --- Function for computing n-body derivatives ---
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
