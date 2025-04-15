% 우주궤도역학 term project#1

% 변수 선언(r, v vector, 지구 중력 상수 μ)
r = [7000; 300; 70];
v = [1; 7.5; 2];
mu = 398600; 
r_norm = norm(r);
v_norm = norm(v);

%angular moment h vector(const)
h = cross(r, v);
h_norm = norm(h);

%이심률(eccentricity)
e_vec = (1/mu) * (cross(v, h) - mu * (r / r_norm));
e = norm(e_vec);


% Inclination
i = acos(h(3) / h_norm);

% Accending Node N, RAAN(Ω)
k_hat = [0; 0; 1];
N = cross(k_hat, h);
N_norm = norm(N);
if N_norm ~= 0
    RAAN = acos(N(1) / N_norm);
    if N(2) < 0
        RAAN = 2*pi - RAAN;
    end
else
    RAAN = 0;
end

% Argument of perigee
if N_norm ~= 0 && e > 1e-6
    w = acos(dot(N, e_vec) / (N_norm * e));
    if e_vec(3) < 0
        w = 2*pi - w;
    end
else
    w = 0;
end

% True Anomaly
if e > 1e-6
    nu = acos(dot(e_vec, r) / (e * r_norm));
    if dot(r, v) < 0
        nu = 2*pi - nu;
    end
else
    nu = 0;
end

% Orbit determine
if abs(e) < 1e-6
    type = 'Circular orbit';
elseif e < 1
    type = 'Elliptical orbit';
elseif abs(e - 1) < 1e-6
    type = 'Parabolic orbit';
else
    type = 'Hyperbolic orbit';
end

fprintf('Orbit type: %s\n', type);
fprintf('Inclination (i): %.6f deg\n', rad2deg(i));
fprintf('RAAN: %.6f deg\n', rad2deg(RAAN));
fprintf('Argument of Perigee (w): %.6f deg\n', rad2deg(w));
fprintf('True Anomaly (nu): %.6f deg\n', rad2deg(nu));