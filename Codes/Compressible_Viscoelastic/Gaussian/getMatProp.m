function Gvec = getMatProp(mesh)

x = mesh.P(:,1);
y = mesh.P(:,2);

mu = 10 + 5 * exp(-((x-60).^2/(20^2) + (y-60).^2/(20^2)));      % in kPa
eta = 1 + 0.5 * exp(-((x-60).^2/(20^2) + (y-60).^2/(20^2)));    % in kPa

Gvec = mu + 1i*eta;

Gvec = 0.001*Gvec;

end
