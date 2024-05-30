%% Compute rupture sizes (length, width, slip) given an earthquake magnitude
%% Based on:
%% 1) Well & Coppersmith, 1994: New Empirical Relationships among Magnitude, Rupture Length,Rupture Width, Rupture Area, and Surface Displacement.
%% and
%% 2) Blaser et al., 2010: Scaling Relations of Earthquake Source Parameter Estimates with Special Focus on Subduction Environment.
%%
%% Rino Salman, June 2020
%% EOS-RS Lab, NTU, Singapore
clear all

%% Magnitude and type of your earthquake
Mw = 7.3;
earthquake_type = 'strike-slip'; %options: 'normal', 'thrust', 'strike-slip'
laws = 'Blaser_et_all'; %options: 'Wells_and_Coppersmith', 'Blaser_et_all'



%% FROM THIS LINE, READ ONLY, DO NOT CHANGE ! %%

%% Compute the rupture sizes
[L, W, slip] = rupture_size_slip(earthquake_type, Mw, laws);


%% Display to screen
fprintf('Scaling laws chosen: %s \n',laws)
fprintf('Earthquake type: %s earthquake\n',earthquake_type)
fprintf('Magnitude size: %.1f\n',Mw)
fprintf('\n')
fprintf('Estimated length: %.0f m\n',L)
fprintf('Estimated width: %.0f m\n',W)
fprintf('Estimated slip: %.2f m\n',slip)


function [L, W, slip] = rupture_size_slip(earthquake_type, Mw, laws)

    %% Compute rupture length (Table 1 of Blaser et al., 2010)
    if strcmp(laws, 'Wells_and_Coppersmith')
        
        switch earthquake_type
            
            case 'thrust'
                La = -2.42;
                Lb = 0.58;
                
            case 'normal'
                La = -1.88;
                Lb = 0.5;
                
            case 'strike-slip'
                La = -2.57;
                Lb = 0.62;
                
        end
        
    else % Blaser et al., 2010
        
        switch earthquake_type
            
            case 'thrust'
                La = -2.37;
                Lb = 0.57;
                
            case 'normal'
                La = -1.91;
                Lb = 0.52;
                
            case 'strike-slip'
                La = -2.69;
                Lb = 0.64;
                
        end
        
    end

    % Compute rupture length based on an equation written in Table 1 (unit: km)
    % log10L = a+bxMw;
    L = 10^(La+Lb*Mw);

    % Convert to meter
    L = L * 1000;


    %% Compute rupture width (Table 2 of Blaser et al., 2010)
    if strcmp(laws, 'Wells_and_Coppersmith')
        
        switch earthquake_type
            
            case 'thrust'
                Wa = -1.61;
                Wb = 0.41;
                
            case 'normal'
                Wa = -1.14;
                Wb = 0.35;
                
            case 'strike-slip'
                Wa = -0.76;
                Wb = 0.27;
                
        end
        
    else % Blaser et al., 2010
        
        switch earthquake_type
            
            case 'thrust'
                Wa = -1.86;
                Wb = 0.46;
                
            case 'normal'
                Wa = -1.20;
                Wb = 0.36;
                
            case 'strike-slip'
                Wa = -1.12;
                Wb = 0.33;
                
        end
        
    end

    % Compute rupture width based on an equation written in Table 2 (unit: km)
    % log10W = a+bxMw;
    W = 10^(Wa+Wb*Mw);

    % Convert to meter
    W = W * 1000;


    %% Compute moment magnitude (unit: dyne-cm)
    % Mw = 2/3log10Mo - 10.7
    Mo = 10^((Mw+10.7)*1.5);

    % Convert to newton-meter (pascal)
    Mo = Mo * 1e-7;

    %% Compute slip
    % Mo = length * width * slip * miu
    % miu assumed 30 GPa
    miu = 3e10;
    slip = Mo/(L*W*miu);

end
