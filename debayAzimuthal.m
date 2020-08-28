classdef debayAzimuthal < RichWolf
    %DEBAY Summary of this class goes here
    %   Detailed explanation goes here
    % continuos polarization and SPP
    
    properties
        zCross = 0;
        rCross = 0;
        
%         aCircE = 1;     % linear by default
%         bCircE = 0;
%         aCircH = 0;
%         bCircH = 1;
        
        mVortex = 0;    % without vortex by default
        mAzimBeam = 1;
    end

    methods(Static)
        function [r, psi] = calcrpsi(x,y)
            
            r = (x.^2 + y.^2).^0.5;
            
            if y >= 0 && x >= 0
                psi = atan(y./x);
            elseif y > 0 && x <= 0
                psi = atan(y./x) + pi;
            elseif y < 0 && x < 0
                psi = atan(y./x) + pi;
            elseif y <= 0 && x >= 0
                psi = atan(y./x) + 2*pi;
            end
        end
        
        function m = mMatrixEx(theta, phi, a, b)
            p11 = 1 + cos(phi).*cos(phi) .* (cos(theta) - 1);
            p12 = sin(phi) .* cos(phi) .* (cos(theta)-1);
            m = p11 .* a + p12 .* b;
        end
        function m = mMatrixEy(theta, phi, a, b)
            p21 = sin(phi) .* cos(phi) .* (cos(theta)-1);
            p22 = 1 + sin(phi).*sin(phi) .* (cos(theta)-1);
            m = p21 .* a + p22 .* b;
        end
        function m = mMatrixEz(theta, phi, a, b)
            p31 = -1 * sin(theta) .* cos(phi);
            p32 = -1 * sin(theta) .* sin(phi);
            m = p31 .* a + p32 .* b;
        end
        
        function m = mMatrixHx(theta, phi, a, b)
            p11 = 1 + cos(phi).*cos(phi) .* (cos(theta) - 1);
            p12 = sin(phi) .* cos(phi) .* (cos(theta)-1);
            m = p11 .* a + p12 .* b;
        end
        function m = mMatrixHy(theta, phi, a, b)
            p21 = sin(phi) .* cos(phi) .* (cos(theta)-1);
            p22 = 1 + sin(phi).*sin(phi) .* (cos(theta)-1);
            m = p21 .* a + p22 .* b;
        end
        function m = mMatrixHz(theta, phi, a, b)
            p31 = -1 * sin(theta) .* cos(phi);
            p32 = -1 * sin(theta) .* sin(phi);
            m = p31 .* a + p32 .* b;
        end
    end
    
    methods
        % calcFWHMs
        function obj = calcFWHMx(obj)
            obj = obj.calcIntensityAlongX();
            obj.FWHMx = obj.calcFWHM(obj.rCoord'/obj.lambda, obj.intensityAlongX');
        end
        function obj = calcFWHMy(obj)
            obj = obj.calcIntensityAlongY();
            obj.FWHMy = obj.calcFWHM(obj.rCoord'/obj.lambda, obj.intensityAlongY');
        end
        function obj = calcFWHMz(obj)
            obj = obj.calcIntensityAlongZ();
            obj.FWHMz = obj.calcFWHM(obj.zCoord'/obj.lambda, obj.intensityAlongZ');
        end
        % <- calcFWHMs

        % calc amplitude E in focal point
        function obj = calcAmpExFocus(obj)
            ampEx = zeros(obj.rNum, obj.rNum);
            for k = 1:obj.rNum
                for kk = 1:obj.rNum
                    [r, psi] = obj.calcrpsi(obj.rCoord(kk), obj.rCoord(k));
                    ampEx(k,kk) = ampEx(k,kk) + integral2(@(theta, phi)underintEx(obj, theta, r, obj.zCross, psi, phi, -1*sin(obj.mAzimBeam*phi), cos(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampExFocus = ampEx;
        end
        function obj = calcAmpEyFocus(obj)
            ampEy = zeros(obj.rNum, obj.rNum);
            for k = 1:obj.rNum
                for kk = 1:obj.rNum
                    [r, psi] = obj.calcrpsi(obj.rCoord(kk), obj.rCoord(k));
                    ampEy(k,kk) = ampEy(k,kk) + integral2(@(theta, phi)underintEy(obj, theta, r, obj.zCross, psi, phi, -1*sin(obj.mAzimBeam*phi), cos(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampEyFocus = ampEy;
        end
        function obj = calcAmpEzFocus(obj)
            ampEz = zeros(obj.rNum, obj.rNum);
            for k = 1:obj.rNum
                for kk = 1:obj.rNum
                    [r, psi] = obj.calcrpsi(obj.rCoord(kk), obj.rCoord(k));
                    ampEz(k,kk) = ampEz(k,kk) + integral2(@(theta, phi)underintEz(obj, theta, r, obj.zCross, psi, phi, -1*sin(obj.mAzimBeam*phi), cos(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampEzFocus = ampEz;
        end           
        % <- calc amplitude E in focal point        
        
        % calc H in focal spot
        function obj = calcAmpHxFocus(obj)
            ampHx = zeros(obj.rNum, obj.rNum);
            for k = 1:obj.rNum
                for kk = 1:obj.rNum
                    [r, psi] = obj.calcrpsi(obj.rCoord(kk), obj.rCoord(k));
                    ampHx(k,kk) = ampHx(k,kk) + integral2(@(theta, phi)underintHx(obj, theta, r, obj.zCross, psi, phi, -1*cos(obj.mAzimBeam*phi), -1*sin(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampHxFocus = ampHx;
        end
        function obj = calcAmpHyFocus(obj)
            ampHy = zeros(obj.rNum, obj.rNum);
            for k = 1:obj.rNum
                for kk = 1:obj.rNum
                    [r, psi] = obj.calcrpsi(obj.rCoord(kk), obj.rCoord(k));
                    ampHy(k,kk) = ampHy(k,kk) + integral2(@(theta, phi)underintHy(obj, theta, r, obj.zCross, psi, phi, -1*cos(obj.mAzimBeam*phi), -1*sin(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampHyFocus = ampHy;
        end
        function obj = calcAmpHzFocus(obj)
            ampHz = zeros(obj.rNum, obj.rNum);
            for k = 1:obj.rNum
                for kk = 1:obj.rNum
                    [r, psi] = obj.calcrpsi(obj.rCoord(kk), obj.rCoord(k));
                    ampHz(k,kk) = ampHz(k,kk) + integral2(@(theta, phi)underintHz(obj, theta, r, obj.zCross, psi, phi, -1*cos(obj.mAzimBeam*phi), -1*sin(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampHzFocus = ampHz;
        end                   
        % <- calc H in focal spot
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% !!!!!!!!!!!!!!!!!!!!
        
        % calc amplitude E in focal point
        function obj = calcAmpExAlongZ(obj)
            ampEx = zeros(obj.rNum, obj.zNum);
            for k = 1:obj.rNum
                for kk = 1:obj.zNum
                    [r, psi] = obj.calcrpsi(0, obj.rCoord(k));
                    ampEx(k,kk) = ampEx(k,kk) + integral2(@(theta, phi)underintEx(obj, theta, r, obj.zCoord(kk), psi, phi, -1*sin(obj.mAzimBeam*phi), cos(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampExAlongZ = ampEx;
        end
        function obj = calcAmpEyAlongZ(obj)
            ampEy = zeros(obj.rNum, obj.zNum);
            for k = 1:obj.rNum
                for kk = 1:obj.zNum
                    [r, psi] = obj.calcrpsi(0, obj.rCoord(k));
                    ampEy(k,kk) = ampEy(k,kk) + integral2(@(theta, phi)underintEy(obj, theta, r, obj.zCoord(kk), psi, phi, -1*sin(obj.mAzimBeam*phi), cos(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampEyAlongZ = ampEy;
        end
        function obj = calcAmpEzAlongZ(obj)
            ampEz = zeros(obj.rNum, obj.zNum);
            for k = 1:obj.rNum
                for kk = 1:obj.zNum
                    [r, psi] = obj.calcrpsi(0, obj.rCoord(k));
                    ampEz(k,kk) = ampEz(k,kk) + integral2(@(theta, phi)underintEz(obj, theta, r, obj.zCoord(kk), psi, phi, -1*sin(obj.mAzimBeam*phi), cos(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampEzAlongZ = ampEz;
        end           
        % <- calc amplitude E in focal point        
        
        % calc H in focal spot
        function obj = calcAmpHxAlongZ(obj)
            ampHx = zeros(obj.rNum, obj.zNum);
            for k = 1:obj.rNum
                for kk = 1:obj.zNum
                    [r, psi] = obj.calcrpsi(0, obj.rCoord(k));
                    ampHx(k,kk) = ampHx(k,kk) + integral2(@(theta, phi)underintHx(obj, theta, r, obj.zCoord(kk), psi, phi, -1*cos(obj.mAzimBeam*phi), -1*sin(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampHxAlongZ = ampHx;
        end
        function obj = calcAmpHyAlongZ(obj)
            ampHy = zeros(obj.rNum, obj.zNum);
            for k = 1:obj.rNum
                for kk = 1:obj.zNum
                    [r, psi] = obj.calcrpsi(0, obj.rCoord(k));
                    ampHy(k,kk) = ampHy(k,kk) + integral2(@(theta, phi)underintHy(obj, theta, r, obj.zCoord(kk), psi, phi, -1*cos(obj.mAzimBeam*phi), -1*sin(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampHyAlongZ = ampHy;
        end
        function obj = calcAmpHzAlongZ(obj)
            ampHz = zeros(obj.rNum, obj.zNum);
            for k = 1:obj.rNum
                for kk = 1:obj.zNum
                    [r, psi] = obj.calcrpsi(0, obj.rCoord(k));
                    ampHz(k,kk) = ampHz(k,kk) + integral2(@(theta, phi)underintHz(obj, theta, r, obj.zCoord(kk), psi, phi, -1*cos(obj.mAzimBeam*phi), -1*sin(obj.mAzimBeam*phi)), 0, obj.thetaMax, 0, 2*pi);
                end
            end
            obj.ampHzAlongZ = ampHz;
        end                   
        % <- calc H in focal spot        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% <-------- !!!!!!!!!!!!!!!!!!!!
        
        
        % Ex calc
        function Ex = underintEx(obj, theta, r, z, psi, phi, a, b)
            Ex = underinteq(obj, theta, phi, z, r, psi);
            Ex = Ex .* obj.mMatrixEx(theta, phi, a, b);    % expression under integral
            Ex = -1i * Ex;      % not exact constant
        end
        % <- Ex calc
        
        % Ey calc
        function Ey = underintEy(obj, theta, r, z, psi, phi, a, b)
            Ey = underinteq(obj, theta, phi, z, r, psi);
            Ey = Ey .* obj.mMatrixEy(theta, phi, a, b);    % expression under integral
            Ey = -1i * Ey;      % not exact constant
        end
        % <- Ey calc
        
        % Ez calc
        function Ez = underintEz(obj, theta, r, z, psi, phi, a, b)
            Ez = underinteq(obj, theta, phi, z, r, psi);
            Ez = Ez .* obj.mMatrixEz(theta, phi, a, b);    % expression under integral
            Ez = -1i * Ez;      % not exact constant
        end        
        % <- Ez calc
        
        % Hx calc
        function Hx = underintHx(obj, theta, r, z, psi, phi, a, b)
            Hx = underinteq(obj, theta, phi, z, r, psi);
            Hx = Hx .* obj.mMatrixHx(theta, phi, a, b);    % expression under integral
            Hx = -1i * Hx;      % not exact constant
        end
        % <- Hx calc
        
        % Hy calc
        function Hy = underintHy(obj, theta, r, z, psi, phi, a, b)
            Hy = underinteq(obj, theta, phi, z, r, psi);
            Hy = Hy .* obj.mMatrixHy(theta, phi, a, b);    % expression under integral
            Hy = -1i * Hy;      % not exact constant
        end
        % <- Hy calc        
        
        % Hz calc
        function Hz = underintHz(obj, theta, r, z, psi, phi, a, b)
            Hz = underinteq(obj, theta, phi, z, r, psi);
            Hz = Hz .* obj.mMatrixHz(theta, phi, a, b);    % expression under integral
            Hz = -1i * Hz;      % not exact constant
        end        
        % <- Hz calc      

        function undeq = underinteq(obj, theta, phi, z, r, psi)     % the equal part under the integral
            T = obj.apodization(theta);
            L = obj.source(theta, phi);         % выкинуть ???
            L = exp(1i*phi*obj.mVortex);   % use in the case of SPP
%             undeq = sin(theta) .* exp(1i*obj.kWave*ones(size(theta)) .* (z*cos(theta) + r * sin(theta) .* cos(phi-psi)) ) .* L .* T;
%             undeq = sin(theta) .* exp(1i*obj.kWave*ones(size(theta)) .* (z*cos(theta) + r * sin(theta) .* cos(phi-psi)) ) .* T; % uniform source
            undeq = sin(theta) .* exp(1i*obj.kWave*ones(size(theta)) .* (z*cos(theta) + r * sin(theta) .* cos(phi-psi)) ) .* T .* L; % uniform source
        end
    end
end


