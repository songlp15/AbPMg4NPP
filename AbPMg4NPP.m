function PP = AbPMg4NPP(PAR0, aph443, a490, bb490)
%% ========================================================================
%  Absorption-based Primary Productivity Model (AbPM)
%  Spectrally integrated implementation 
%
%  Created by Luping Song and Jinghui Wu
%  Apr 6, 2026
%
%  Affiliation:
%    - Xiamen University (XMU)
%    - Zhejiang Ocean University (ZJOU)
%    - Lamont-Doherty Earth Observatory (LDEO), Columbia University
%
%  DESCRIPTION
%  This function estimates depth-integrated primary production (PP)
%  using the absorption-based primary productivity model
%  (AbPM), in which phytoplankton absorption is spectrally averaged
%  over the photosynthetically active radiation (PAR) range.
%
%  In this formulation, implicitly assuming a depth-invariant 
% light–absorption relationship.
%  ------------------------------------------------------------------------
%  INPUTS
%
%   PAR0  : Surface Photosynthetically available radiation (PAR) above
%   sea surface (Einstein m⁻² d⁻¹)；PAR = PAR0*0.96
%
%   aph443  : Phytoplankton absorption at 443 nm (m⁻¹)
%
%   a490    : Total absorption coefficient at 490 nm (m⁻¹)
%
%   bb490   : Total backscattering coefficient at 490 nm (m⁻¹)
%  ------------------------------------------------------------------------
%  OUTPUT
%
%    PP    : Depth-integrated primary production (mg C m⁻² d⁻¹)
%
%  ------------------------------------------------------------------------
%   other key parameters:
%
%      Phi_z               Energy Absorption
%      fm           maximum quantum yield of photosynthesis                   mol C (mol photons)-1
%      Kf           the irradiance when Phi corresponds to a half of fm       mol m-2 d-1
%      a0,a1 : Spectral coefficients used to reconstruct aph(λ) from Lee et
%               al., 1998
%       wl    : Wavelength (nm), e.g., 400:1:700
%  ------------------------------------------------------------------------
%  References and Credits
%
% Lee, Z., et al. (2011), An assessment of optical properties and primary 
% production derived from remote sensing in the Southern Ocean (SO GasEx), 
% J. Geophys. Res., 118, 4241?4255, doi:10.1002/jgrc.20308.
%   
% ========================================================================
%% ---------------- Spectral grid (400–700 nm) ----------------
wl = (400:1:700)';     % nm

%% ---------------- HOPE / Hydrolight spectral parameters ----------------
% These parameters are treated as fixed spectral shapes

load('a0_a1.mat','a0','a1');   % a0/a1 from HOPE
a0 = a0(:); a1 = a1(:); %301*1
%% ---------------- Model parameters ----------------
fm  = 0.06;     % default maximum quantum yield (mol C mol photons⁻¹)
kf  = 10;       % half-saturation irradiance (Ein m⁻² d⁻¹)
v   = 0.01;     % photoinhibition parameter (Ein m⁻² d⁻¹)⁻¹
phi_sun = 45;   % solar zenith angle (degrees)

%% ---------------- Step 1: IOPs ----------------
if ~isfinite(aph443) || aph443 <= 0 || ...
   ~isfinite(a490)   || a490   <= 0 || ...
   ~isfinite(bb490)  || bb490  <= 0

    PP = nan;
    return
end
%% ----------------Step 2: Phytoplankton absorption ----------------
% Phytoplankton absorption (HOPE model) (Lee et al., 1999)
aph = aph443 .* (a0 + a1 .* log(aph443));      % aph(λ)

% old scheme (spectrally  averaged aph)
aph_mean = trapz(wl, aph) ./ (700 - 400);

%% ---------------- Step 3: Vertical Attenuation ----------------
% Lee et al. (2007)

PAR = PAR0 * 0.96;% PAR just below surface

chi0   = -0.057;
chi1   = 0.48245;
chi2   = 4.2213;

theta0 = 0.18319;
theta1 = 0.70238;
theta2 = -2.5673;

alpha0 = 0.09;
alpha1 = 1.46545;
alpha2 = -0.6666;

sza_rad = phi_sun * pi / 180;

K1 = (chi0 + chi1 * sqrt(a490) + chi2 * bb490) ...
    * (1 + alpha0 * sin(sza_rad));

K2 = (theta0 + theta1 * a490 + theta2 * bb490) ...
    * (alpha1 + alpha2 * cos(sza_rad));

%% ---------------- Step 4-6: Depth-integrated PP ----------------
z_eu = 200;%Fixed depth integration (0–200 m) 

NPP_z = zeros(z_eu,1);

for z = 1:z_eu

    % Spectral irradiance at depth
     E_z = PAR .* exp(-(K1 + K2 ./ sqrt(1 + z)) .* z); % PAR(z)
     % E_z = trapz(wl, E_lambda);  

    % Quantum yield
    Phi_z = (kf ./ (kf + E_z)) .* fm;

    % Conventional AbPM: aph_mean × PAR(z)
    NPP_z(z) = Phi_z .* aph_mean .* E_z .* exp(-v .* E_z);

end
%% ---------------- Step 7:  Depth integration----------------

PP = trapz(NPP_z) * 12.0107 * 1000;   % mg C m⁻² d⁻¹

end