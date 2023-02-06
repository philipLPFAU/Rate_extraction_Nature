% Minimal example for the rate extraction procedure from a wavefunction propagated by numerically solving the time-dependent Schrödinger equation as used in
% "Tracing attosecond electron emission from a nanometric needle tip" 
% by Philip Dienstbier, Lennart Seiffert, Timo Paschen, Andreas Liehl, Alfred Leitenstorfer, Thomas Fennel and Peter Hommelhoff. 
%
% The corresponding discussions on the evaluation of the resonances, their out-projection from a propagated wave function and the reate extraction
% are detailed in the methods section and supplementary information of the manuscript.
%
% The underlying theory model for the nanotip is adapted from [L. Seiffert et al., J. Phys. B: At. Mol. Opt. Phys. 51, 134001 (2018)].
%
% Note: The grid parameters (spatial & temporal) have been reduced w.r.t. the calculations performed for the published work
% to minimize the calculation time for this minimal example
%
% Tested with MATLAB version R2017b and above

%% Clear workspace and command window
clear;
clc;


%% Define physical constants and conversion factors
c_SI = 3.0e8;                   % speed of light in [m/s] 
epsilon_0_SI = 8.85418782e-12;  % vacuum permittivity in [As/Vm]
E_au_to_ev = 27.211386245988;   % energy conversion [a.u.] -> [eV]
x_au_to_ang = 0.529177249;      % position conversion [a.u.] -> [angstroms]
t_au_to_fs =  0.024188843;      % time conversion [a.u.] -> [fs]
F_au_to_SI = 5.1422e11;         % electric field strength conversion [a.u.] -> [V/m]


%% Define spatial and temporal grids
x = linspace(-50, 250, 5001) / x_au_to_ang;  % spatial grid (1xN vector) given in [angstroem] and converted to [a.u.]
t = linspace(-25, 15, 5000) / t_au_to_fs;   % temporal grid (1xM vector) given in [fs] and converted to [a.u.]
    

%% Calculate ground state and corresponding binding potential for a given pair of work function and Fermi energy
% The binding potential is determined by adjusting the potential width via bisection 
% such that the ground state (determined via imaginary time propagation) reflects the desired work function and fermi energy
E_W = 6.6 / E_au_to_ev; % work function given in [eV] and converted to [a.u.]
E_F = 7.0 / E_au_to_ev; % Fermi energy given in [eV] and converted to [a.u.]

% Find ground state PSI_0 (1xN vector) with energy expectation value E_PSI_0 in [a.u.] 
% together with respective binding potential V_0 (1xN vector) in [a.u.] with width w in [a.u.]
[PSI_0, E_PSI_0, V_0, w] = Find_ground_state_and_potential(x, E_W, E_F);

% Plot binding potential
figure('Units', 'normalized', 'Position', [0 0 1 1], 'color', 'w');
subplot(2,4,[1,5]);
Plot_potential(x*x_au_to_ang, V_0*E_au_to_ev)


%% Calculate resonances for set of (negative) field strengths and truncated perturbed potential
E_0_SI_static = (0:-0.1:-0.3)*1e9; % Set of field strengths in V/m (1xL vector)
E_0_SI_static_truncation = -1e9; % Truncation field strength in V/m

% Calculate and plot resonances PSI_RES (LxN matrix) as described in the methods section of the paper
PSI_res = Calculate_resonances(E_0_SI_static, x, V_0, E_PSI_0, E_0_SI_static_truncation, E_W);


%% Laser pulse parameters
I_0_SI = 7.5e16;    % intensity in [W/m^2] 
lambda_SI = 1560;   % wavelength in [nm]
tau_SI = 10;        % pulse duration in [fs]
phi_cep = pi;       % carrier-envelope phase in [rad]

% Calculate field strength and cycle duration and convert to atomic units
E_0 = sqrt(2*I_0_SI/(epsilon_0_SI*c_SI)) / F_au_to_SI;  % peak field strength in [a.u.]
T = (lambda_SI*1e-9)/c_SI*1e15/t_au_to_fs;              % cycle duration in [a.u.]
tau = tau_SI/t_au_to_fs;                                % pulse duration in [a.u.]
omega = 2*pi/T;                                         % angular frequency in [a.u.]

% Calculate time-dependent waveform
E = E_0*exp(-2*log(2)*(t/(tau)).^2).*cos(omega*t+phi_cep); % optical waveform in [a.u.] (1xM vector)


%% Calculate spatial potential profile
% Homogeneous field vanishing for x < 0 (that is inside the metal)
V_spatial = x; V_spatial(x<0) = 0; % sptial profile (1xN vector)

% Introduce constant potential after 20 atomic units to reduce backscattering
x_const_index = find(x>20,1,'first');
V_spatial(x_const_index:end) = V_spatial(x_const_index); 


%% Propagate wavefunction 
% Time propagation provides wave function PSI_reduced ~(1xM/reduction_factor vector) on reduced time grid ~(1xM/reduction_factor vector)
% Note: a small time step is required for the propagation, but not for the rate extraction 
subplot(2,4,[2,6]);
handle = Plot_wavefunction(x, t, PSI_0*0, 'Time-dependent probability density');
cb = colorbar('Location', 'west'); 
cb_pos = get(cb, 'Position');
set(cb, 'Position', [cb_pos(1), cb_pos(2) + cb_pos(4)/2, cb_pos(3), cb_pos(4)/2]);
zlab = get(cb,'ylabel'); 
set(zlab,'String',{'Log. prob. density [arb. units]'}); 
plot(t*t_au_to_fs, 50 + E*F_au_to_SI/1e9, 'r', 'Linewidth', 2, 'DisplayName', 'Waveform E(t)');
legend show;
drawnow;

PSI_propagated = Propagate_wave_function(x, t, PSI_0, V_0, V_spatial, E, handle);


%% Remove resonances via out-projection
% Initialize wave function for resonance removal
PSI_after_projection = PSI_propagated;

for l = 1:size(PSI_res,2)
    % Remove resonance number (l) for each time step
    for k = 1:length(t)
        % Projection step
        PSI_after_projection(:,k) = PSI_after_projection(:,k)-projection(x,PSI_after_projection(:,k),PSI_res(:,l));
     
        % Calculate yield = norm of wavefunction after projection step restricted to vacuum side (x>0)
        yield(l,k) = trapz(x(x>0), abs(PSI_after_projection(x>0,k)).^2);
    end
    
    % Calculate rate as time derivative of the yield
    rate(l,:) = gradient(yield(l,:)) ./ gradient(t);
    
    % Show wavefunction after each projection step depending on position x and time t
    subplot(2,4, l + 2 + 2*(l>2));
    Plot_wavefunction(x, t, PSI_after_projection,['Density after projecting out \Psi_{res}^{(',num2str(l-1),')}']);
    drawnow;
end

% plot the yield and rate for the last wave function
plot(t*t_au_to_fs, 100 * yield(l,:) / max(yield(l,:)), 'k', 'Linewidth', 2, 'DisplayName', 'Yield Y(t) [arb. units]');
plot(t*t_au_to_fs, 100 * rate(l,:) / max(rate(l,:)), 'r', 'Linewidth', 2, 'DisplayName', 'Rate R(t) [arb. units]');
legend('Location', 'northwest');











%% function definitions
function [PSI_propagated] = Propagate_wave_function(x, t, PSI_0, V_0_input, V_spatial, E_t, handle)
    % Description: The function propagates an initial wave function considering its binding potential and an interaction potential from the laser field.
    % The implementation is based on the Crank-Nicolson scheme and closely follows Chapter I of "Dieter Bauer - Computational Strong-field quantum dynamics (De Gruyter, 2017)".
    % Parameters:   x - spatial grid in au (1xN vector)
    %               t - time grid in au (1xM vector)
    %               PSI_0 - initial wave function (1xN vector)
    %               V_0_input - binding potential in au (1xN vector)
    %               V_spatial - spatial potential profile in au (1xN vector)
    %               E_t - optical waveform in au (1xM vector)
    %               reduction_factor - reduces number of saved time steps by given factor (scalar)
    % Output:       t_reduced - reduced time grid ~(1xM/reduction_factor)
    %               PSI_propagated_reduced - propagated wavefunction saved on reduced time grid ~(1xM/reduction_factor)
     
    % Spatial and temporal step width in au
    Delta_x = x(2)-x(1); 
    Delta_t = t(2)-t(1);
    
    % Spatial and temporal grid sizes
    M = length(x);
    N = length(t);
    
    % Reduced time axis for saving the wave function
    PSI_propagated = zeros(M, length(t));
    
    % Absorbing boundary conditions to avoid reflections at spatial boundaries
    V_OI = -10*(1i); alpha = 7.5;
    V_r = V_OI ./ (cosh((x(end)-x)/alpha)).^2;
    V_l = V_OI ./ (cosh((x-x(1))/alpha)).^2;
    
    % Second spatial derivative in matrix form (d^2/dx^2)
    D = 1/(Delta_x^2)*(sparse(1:M,1:M,-2,M,M)+sparse(2:M,1:M-1,1,M,M)+sparse(1:M-1,2:M,1,M,M)); 
    
    % Potential part of Hamiltonian
    V_0 = sparse(1:M,1:M,V_0_input+V_r+V_l,M,M); 
    
    % Precalculated unit matrix
    unity = sparse(1:M,1:M,1,M,M);
    
    % Initial wave function
    PSI = PSI_0;

    % Time propagation
    for n=1:N
        % Interaction potential
        V_int = sparse(1:M, 1:M, V_spatial*E_t(n), M, M);
        
        % Effective potential V(x,t) = V_0(x) + V_int(x,t)
        V = V_0 + V_int;
    
        % Propagation matrices
        A_plus  = unity - (1i*Delta_t/2)*(-(1/2)*D+V);
        A_minus = unity + (1i*Delta_t/2)*(-(1/2)*D+V);
    
        % Solve equation A_minus PSI_(n+1) = A_plus PSI_n = b
        b = A_plus*PSI; 
        PSI = A_minus\b;
    
        % Only save wave function on reduced time grid
        PSI_propagated(:,n) = PSI;

        if mod(n,10)==0
            t_au_to_fs =  0.024188843; % time conversion au -> fs
            if t(n)*t_au_to_fs>=-15.5 && t(n)*t_au_to_fs<=15.5
                x_reduction = 1:10:length(x);
                set(handle, 'CData', log10(abs(PSI_propagated(x_reduction,:)).^2));
                drawnow;
            end
        end
    
        % Status message
        msg = sprintf('Propagating wave function: %3.1f / 100.0 ... ', n/N*100);
        if(n==1)
            fprintf(['\n',msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        else
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
    fprintf('done\n');
end

function [PSI_min,E_exp_min,V_0,w] = Find_ground_state_and_potential(x,E_W,E_f)
    % Description: The function finds the ground state and corresponding potential for a given work function and Fermi-energy.
    % Parameters:   x - spatial grid in au (1xN vector)
    %               E_W - work function in au (scalar)
    %               E_F - Fermi-level in au (scalar)
    % Output:       PSI_min - ground state fullfilling E_W and E_F (1xN)
    %               E_exp_min - energy expectation value of determined ground state in au (scalar)
    %               V_0 - binding potential in au (1xN)
    %               w - width of binding potential in au (scalar)

    E_au_to_ev = 27.211386245988; % energy conversion au -> eV
    
    % Determined ground state energy needs to fullfil this error boundary
    max_error_SI = 0.001; 
    max_error = max_error_SI/E_au_to_ev; % error boundary in au
    
    % Maximum amount of bisection trials
    max_trials = 50; 
    
    % Depth of box potential in au
    U_0 = E_W+E_f;
    
    % Initial search interval encapsuling the correct potential width
    % w_initial = pi/sqrt(2*E_f);   % Initial guess of potential width in au using an infinite box potential
    % w_int = [0, w_initial];       % for arbitrary combinations of work function and fermi energy
    w_int = [2.1, 2.2]; % optimized initial guessfor the parameters in the paper
    
    fprintf('Calculating ground state and corresponding potential ... ');
    
    % In each iteration the middle of the search interval is used. The lower or
    % upper half of the interval is chosen in the next iteration depending on the 
    % estimated ground state energy (=bisection).
    for n = 1:max_trials
        
        % Middle of interval
        w = mean(w_int);
        
        % A box potential with smooth edges (Fermi-box) is chosen for the binding potential.
        s = 0.01;
        V_0 = U_0*(1-1./(exp((x)/s)+1)-1./(exp((-(x+w))/s)+1));

        % Imaginary time propagation provides the ground state wavefunction
        % PSI and its energy expectation value E_exp.
        [PSI,E_exp] = Imaginary_time_propagation_ground_state(x, V_0, -E_W);
        
        % Compute energy difference from target value.
        E_diff = (-E_W)-E_exp;
    
        % If current expectation value is above the target value, use upper
        % half of search interval (=large potential widths). Otherwise, use
        % lower half of search interval (=smaller potential widths).
        if(E_exp>-(E_W)  || isnan(E_exp))
            w_int = [w,w_int(2)];
        else
            w_int = [w_int(1),w];         
        end
    
        % Keep track of solution with smallest deviation from target value
         if(n==1)
             E_diff_min = abs(E_diff);
             PSI_min = PSI;
             E_exp_min = E_exp;
         end
         
         if(isnan(E_diff_min)&&~isnan(E_diff))
             E_diff_min = abs(E_diff);
             PSI_min = PSI;
             E_exp_min = E_exp;
         end
         
         if(E_diff_min>abs(E_diff))
             E_diff_min = min(E_diff_min,abs(E_diff));
             PSI_min = PSI;
             E_exp_min = E_exp;
         end
        
        % Stop iteration if error boundary is fullfilled
        if(abs(E_diff)<max_error)
            break
        end
    end

    fprintf('done\n\n');

    % Show if maximum number of trials has been reached
    if(n == max_trials)
        disp(['Maximum of trials reached, returning best fit!']);
    end
end

function [PSI,E_exp] = Imaginary_time_propagation_ground_state(x, V_0, E_guess)
    % Description: The function uses imaginary time propagation to
    % find the ground state for a given binding potential. The implementation is based
    % on the Crank-Nicolson scheme and closely follows Chapter I of " Dieter Bauer - Computational
    % Strong-field quantum dynamics (De Gruyter, 2017)".
    % Parameters:   x - spatial grid in au (1xN vector)
    %               V_0 - binding potential in au (1xN vector)
    %               E_guess - guess of the eigen energy (to set the imaginary time step for fast convergence)
    % Output:       PSI - ground state wavefunction (1xN)
    %               E_exp - energy expectation value of ground state in au (scalar)

    % Length and step width of spatial grid
    M = length(x);
    Delta_x = x(2)-x(1);
    
    % Imaginary time step and number of time steps
    Delta_t = 2i / E_guess;
    N_t_max = 3000;
    
    % Second spatial derivative in matrix form (d^2/dx^2)
    D = 1/(Delta_x^2)*(sparse(1:M,1:M,-2,M,M)+sparse(2:M,1:M-1,1,M,M)+sparse(1:M-1,2:M,1,M,M));
    
    % Potential part of Hamiltonian
    V_0 = sparse(1:M,1:M,V_0,M,M); 
    
    % Complete Hamiltonian with kinetic term -1/2D and potential V_0 
    H_0 = (-(1/2)*D+V_0);
    
    % Propagation matrices
    unity = sparse(1:M,1:M,1,M,M);
    A_plus = unity-(1i*Delta_t/2)*H_0;
    A_minus = unity+(1i*Delta_t/2)*H_0;
    
    % Start with random PSI 
    PSI = rand(M,1);
    E_exp = nan;
    
    % Imaginary time propagation 
    for n=1:N_t_max
            % Solve equation A_minus PSI_(n+1) = A_plus PSI_n = b
            b = A_plus*PSI; 
            PSI = A_minus\b;
    
            % Normalize (required in imaginary time propagation)
            PSI = PSI / sqrt(trapz(x, abs(PSI).^2));
            
            % Energy expectation value <psi|H_0|psi>
            E_exp = trapz(x, conj(PSI).*(H_0*PSI));
            
            % Calculate residual <psi_res|psi_res> using |psi_res> = (H_0-E)|psi>
            PSI_res = (H_0-E_exp*unity)*PSI;
            R = trapz(x,(PSI_res).^2);
            
            % Check if wavefunction properly decays quickly outside the binding potential using the median
            R_base_level = median(abs(PSI).^2);
    
            % If the residual is small and the wavefunction decays properly we found the ground state accurately
            if(R<1e-25 && R_base_level<1e-75) 
                break
            end
    end
end

function [PSI_res] = Calculate_resonances(E_0_SI_static,x,V_0,E_PSI_0,E_0_SI_static_truncation, E_W)
    % Description: The function calculates a set of resonances later used to extract the instantaneous rate. The resonances form an orthonormal set
    % of functions, which describe any static perturbation of the ground state and are a good approximation for the dynamical case. For details see
    % methods section of the manuscript.
    % Parameters:   E_0_SI_static - set of static field strengths in au (1xL vector)
    %               x - spatial grid in au (1xN vector)
    %               V_0 - unperturbed binding potential in au (1xN vector)
    %               E_PSI_0 - binding energy of unperturbed ground state in au (scalar)
    %               E_0_SI_static_truncation - truncation field strength in au (scalar)
    % Output:       PSI_res - resonances (MxL matrix)
    
    F_au_to_SI = 5.1422e11; % electric field strength conversion au -> V/m
    x_au_to_ang = 0.529177249; % position conversion au -> angstroms
    E_au_to_ev = 27.211386245988; % energy conversion au -> eV

    % Convert static field strengths and truncation field strength to au
    E_0_static = E_0_SI_static/F_au_to_SI;
    E_0_static_truncation = E_0_SI_static_truncation/F_au_to_SI;

    % Compute perturbed ground states for the set of field strengths via imaginary time propagation 
    for i = 1:length(E_0_SI_static)
        fprintf('Calculating resonance Psi_res^(%d) (%d / %d) ... ', i-1, i, length(E_0_SI_static));

        % Calculate perturbed potentials
        V_E(:,i) = Compute_perturbed_potential(x,V_0,E_0_static(i),E_PSI_0,E_0_static_truncation);

        % Apply imaginary time propagation to find perturbed ground states
        [PSI_res(:,i),~] = Imaginary_time_propagation_ground_state(x,V_E(:,i), -E_W);
        drawnow;

        % Remove all previous resonances via modified Gram-Schmidt orthogonalization 
        for j = 1:(i-1)
            PSI_res(:,i) = PSI_res(:,i)-projection(x,PSI_res(:,i),PSI_res(:,j));
        end

        % Normalize wave function
        PSI_res(:,i) = PSI_res(:,i)./sqrt(trapz(x,abs(PSI_res(:,i)).^2));

        % Calculate energy expecation value
        M = length(x);
        Delta_x = x(2)-x(1);
        D = 1/(Delta_x^2)*(sparse(1:M,1:M,-2,M,M)+sparse(2:M,1:M-1,1,M,M)+sparse(1:M-1,2:M,1,M,M));
        V0 = sparse(1:M,1:M,V_0,M,M); 
        H0 = (-(1/2)*D+V0);
        E_res(i) = trapz(x, conj(PSI_res(:,i)) .* (H0*PSI_res(:,i)));

        
        % For visualization we chose that each resonance has positive amplitude coming from positive position values
        [~,x_index] = min(abs((x*x_au_to_ang)-5));
        if sum(PSI_res(x_index,i))<0
            PSI_res(:,i) = -PSI_res(:,i);
        end

        % Plot resonance
        plot(x*x_au_to_ang, E_res(i)*E_au_to_ev + PSI_res(:,i), 'Linewidth', 2, 'DisplayName', ['\Psi_{res}^{(',num2str(i-1),')}']);
        drawnow;

        fprintf('done\n');
    end
end

function [] = Plot_potential(x, V)
    % Description: The function plots the binding potential
    % Parameters:   x - spatial grid in [angstroem]
    %               V - binding potential in [eV]
    plot(x, V, 'k', 'Linewidth', 2, 'DisplayName', 'Potential');
    hold on;
    xlim([-5,10]);
    xlabel('Position (angstrom)');
    ylabel('Energy (eV) & Wavefunction \Psi_{res}^{(n)}');
    legend('Location', 'southeast');
    drawnow;
end

function [handle] = Plot_wavefunction(x,t,PSI,title_string)
    % Description: The function plots the propagated wave function depending on position and time.
    % Parameters:   x - spatial grid in au (1xN vector)
    %               t - (reduced) time grid in au (1xM vector)
    %               PSI - propagated wave function (NxM)
    %               title_string - title of plot (char array/string)

    % Conversion constants (from atomic units to SI units)
    x_au_to_ang = 0.529177249; % position conversion au -> angstroms
    t_au_to_fs =  0.024188843; % time conversion au -> fs
    
    % Reduce image size along position axis
    x_reduction = 1:10:length(x);
    
    % Plot wavefunction depending on position x and time t using SI units
    handle = imagesc(t*t_au_to_fs,x*x_au_to_ang,log10(abs(PSI(x_reduction,:)).^2));
    hold on;
    colormap(IBcolor(128));
    set(gca,'YDir','normal');
    ylabel('Position (angstrom)'); 
    xlabel('Time (fs)');
    xlim([-15,15]);
    ylim([-10,150]);
    caxis([-12,-2]);
    title(title_string);
end

function [V] = Compute_perturbed_potential(x,V_0,E_0,E_PSI_0,E_0_boundary)
    % Description: The function adds the potential from a static field to the
    % binding potential. The resulting potential is truncated by the truncation
    % field strength E_0_boundary. The value for the truncation field strength 
    % needs to be chosen such that the ground state remains within the perturbed potential.
    % Parameters:   x - spatial grid in au (1xN vector)
    %               V_0 - bingin potential in au (1xN vector)
    %               E_0 - static field strength in au (scalar)
    %               E_PSI_0 - ground state energy in au (scalar)
    %               E_0_boundary - truncation field strength in au (scalar)
    % Output:       V - perturbed and truncated static potential (1xN vector)

    % Calculate perturbing potential for a homogeneous field in au
    V_E = x*E_0;V_E(x<0) = 0;
    
    % Add perturbation to binding potential
    V = V_0+V_E;
    
    % Truncate perturbed potential to keep the perturbed ground state within 
    % the initial binding potential.
    V_E_boundary = x*E_0_boundary;
    x_extension = find(V_E_boundary<E_PSI_0,1,'first')-1;
         
    if(~isempty(x_extension))
        V(x>x(x_extension)) = 0;
    end
end

function [p] = projection(x, y1, y2)
    % Description: The function calculates the projection of the wave function
    % (or generally two complex vectors) y1 onto y2. 
    % Parameters:   x - spatial grid in au (1xN vector)
    %               y1 - first wave function (1xN vector)
    %               y2 - second wave function (1xN vector)
    % Output:       p - projection of y1 onto y2 (1xN)

    p = trapz(x,conj(y2).*y1).*y2./(trapz(x,abs(y2).^2));
end

function y = IBcolor(colorN)
    % Description: This function creates colorbar function in IB style
    % Parameters:  colorN - number of color steps (scalar)
    % Output:      y - colormap (colorNx3 matrix)
    
    % List of color points
    Rx = [0    ,18   ,28   ,40   ,50    ,55    ,60    ,70    ,85    ,92 ,100 ,105 ,118]' / 118;
    Ry = [0.75 ,0.92 ,0.975 ,0.6 ,0.135 ,0.125 ,0.125 ,0.125 ,0.125 ,0.12 ,0.2 ,0.55 ,1]';
    Gx = [0    ,5     ,8     ,10  ,15    ,20   ,28   ,32  ,40  ,48 ,60   ,65   ,70    ,80    ,90    ,95    ,98, 105 ,118]' / 118;
    Gy = [0.09 ,0.092 ,0.125 ,0.2 ,0.575 ,0.71 ,0.91 ,0.9 ,0.8 ,0.625 ,0.63 ,0.65 ,0.675 ,0.675 ,0.425 ,0.175 ,0.15 ,0.55 ,1]';
    Bx = [0    ,20   ,30   ,40   ,50    ,55   ,58    ,59  ,61 ,65 ,70   ,74    ,81  ,92  ,100 ,105 ,118]' / 118;
    By = [0.12 ,0.12 ,0.12 ,0.12 ,0.135 ,0.15 ,0.17 ,0.2 ,0.3 ,0.5 ,0.75 ,0.925 ,0.9 ,0.5 ,0.5 ,0.7 ,1]';
    
    % Interpolation
    Xq = (0:colorN-1)'/(colorN-1);
    R = interp1(Rx,Ry,Xq,'linear');
    G = interp1(Gx,Gy,Xq,'linear');
    B = interp1(Bx,By,Xq,'linear');
    
    y = flipud([R G B]);
end