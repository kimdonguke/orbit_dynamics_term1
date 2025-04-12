% 우주궤도역학 term project#1 (i=0, RAAN=0, w=0)
r = [7000; 0; 0];
v = [0; 7.5; 0];
mu = 398600;


r_norm = norm(r);
v_norm = norm(v);

%angular moment vector
h = cross(r, v);
h_norm = norm(h);

%이심률(eccentricity)
e_vec = (1/mu) * (cross(v, h) - mu * (r / r_norm));
e = norm(e_vec);

%specific mechanical energy
epsilon = v_norm^2 / 2 - mu / r_norm;

%orbit determine
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