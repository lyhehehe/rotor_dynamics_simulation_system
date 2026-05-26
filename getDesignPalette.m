%% getDesignPalette - Return project color palette as a nested structure
%
% This function provides a centralised set of RGB color values for consistent
% plot styling across all visualization functions in the project.
%
%% Syntax
%  c = getDesignPalette()
%
%% Description
% |getDesignPalette| returns a nested struct |c| organized by category,
% color family, and shade level. Each leaf field holds a 1×3 RGB row vector.
% Access pattern: |c.Category.Color.sShade| (e.g., |c.Accents.Blue.s4|).
% The palette contains:
% * Background: Stone (warm grey), Grey (cool grey) — each with 6 shades (s1–s6)
% * Accents: Red, Blue, Yellow — primary highlight colors
% * Extended: Olive, Green, Teal, Purple, Orange, Skin — supplementary colors
%   (also aliases Red, Blue, Yellow from Accents for uniform access)
%
%% Output Arguments
% * |c| - Nested color palette structure [struct]
%   Access: |c.Category.Color.sShade| where shade index s1 (lightest) to s6 (darkest)
%
% Copyright (c) 2021-2026 Haopeng Zhang, Northwestern Polytechnical University, Politecnico di Milano
% This code is licensed under the MIT License. See the LICENSE file in the project root for the full text of the license.
%

function c = getDesignPalette()
    % Use: c.Category.Color.sShade (e.g., c.Accents.Blue.s4)

    % Helper function to convert a 6-row matrix into a shade structure
    makeShades = @(m) struct('s1', m(1,:), 's2', m(2,:), 's3', m(3,:), ...
                             's4', m(4,:), 's5', m(5,:), 's6', m(6,:));

    %% 1. MAIN BACKGROUND
    c.Background.Stone = makeShades([0.97,0.96,0.93; 0.91,0.89,0.84; 0.82,0.80,0.74; 0.69,0.67,0.60; 0.53,0.52,0.44; 0.38,0.37,0.31]);
    c.Background.Grey  = makeShades([0.94,0.95,0.97; 0.85,0.87,0.91; 0.72,0.76,0.82; 0.55,0.59,0.67; 0.38,0.42,0.49; 0.20,0.23,0.28]);

    %% 2. MAIN ACCENTS
    c.Accents.Red    = makeShades([0.98,0.85,0.85; 0.94,0.65,0.65; 0.89,0.40,0.40; 0.80,0.20,0.20; 0.65,0.15,0.15; 0.48,0.08,0.08]);
    c.Accents.Blue   = makeShades([0.82,0.92,0.99; 0.61,0.82,0.94; 0.35,0.62,0.84; 0.00,0.44,0.74; 0.00,0.32,0.58; 0.00,0.20,0.40]);
    c.Accents.Yellow = makeShades([1.00,0.95,0.80; 1.00,0.88,0.58; 0.95,0.78,0.30; 0.82,0.65,0.15; 0.64,0.50,0.08; 0.45,0.35,0.04]);

    %% 3. EXTENDED PALETTE
    c.Extended.Olive  = makeShades([0.96,0.97,0.73; 0.89,0.92,0.42; 0.76,0.80,0.00; 0.58,0.64,0.00; 0.41,0.48,0.06; 0.27,0.32,0.08]);
    c.Extended.Green  = makeShades([0.89,0.94,0.83; 0.69,0.84,0.59; 0.47,0.70,0.36; 0.23,0.53,0.23; 0.13,0.40,0.16; 0.07,0.26,0.11]);
    c.Extended.Teal   = makeShades([0.88,0.95,0.96; 0.65,0.86,0.88; 0.35,0.73,0.78; 0.00,0.60,0.68; 0.00,0.45,0.52; 0.00,0.30,0.36]);
    c.Extended.Purple = makeShades([0.93,0.86,0.94; 0.82,0.66,0.84; 0.69,0.46,0.71; 0.55,0.25,0.58; 0.43,0.15,0.46; 0.29,0.09,0.32]);
    c.Extended.Orange = makeShades([1.00,0.90,0.80; 1.00,0.78,0.55; 0.98,0.58,0.25; 0.93,0.42,0.00; 0.75,0.31,0.00; 0.53,0.20,0.00]);
    c.Extended.Skin   = makeShades([1.00,0.93,0.88; 0.91,0.80,0.71; 0.77,0.64,0.54; 0.62,0.50,0.42; 0.50,0.39,0.33; 0.35,0.27,0.22]);
    
    % Reusing Accents in Extended for completeness
    c.Extended.Red    = c.Accents.Red;
    c.Extended.Blue   = c.Accents.Blue;
    c.Extended.Yellow = c.Accents.Yellow;
end