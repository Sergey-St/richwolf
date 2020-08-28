classdef radial < debayRadial
    %ABGMODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods

        function obj = radial(lambda, thetaMin, NA)
            obj.lambda = lambda;
            obj.kWave = 2*pi/lambda;
            obj.NA = NA;
            obj.thetaMax = asin(NA/obj.refrMedia);
            obj.thetaMin = thetaMin;
        end

        function L = source(obj, theta, phi)
            L = 1;
        end

        function T = apodization(obj, theta)
            T = (cos(theta)).^(-1.5);      % apodization function for flat difractive lens 
        end

        function plotsource(obj)
%             R = 25;
%             alpha = obj.bettacoef / R;
%             n = obj.ncoef;
%             c = obj.ccoef;
% 
%             rinMin = -2*R;
%             rinMax = 2*R;
%             rinNum = 200;
%             rin = linspace(rinMin, rinMax, rinNum);
%             ein = zeros(size(rin,2), size(rin,2));
%             for kk = 1:rinNum
%                 for k = 1:rinNum
%                     [r, phi] = obj.calcrpsi(rin(k),rin(kk));
%                     
%                     ein(kk,k) = (alpha*r * (alpha*r - 2*c*exp(1i*phi)) ) ^ 0.5;
%                     ein(kk,k) = besselj(n, ein(kk,k));
%                     ein(kk,k) = ein(kk,k) .* (alpha*r ./ (alpha*r - 2*c*exp(1i*phi))) .^ (n/2);
%                     ein(kk,k) = ein(kk,k) .* exp(1i*n*phi);
%                 end
%             end
%             figure;
%                 imagesc(rin, rin, ein.*conj(ein));
%                 axis xy;
%                 colormap('hot');
%                 colorbar;
%                 L1 = xlabel('x, \mum');
%                 L2 = ylabel('y, \mum');
%                 modiffigures(gcf, gca);
%                 modiflabels([L1 L2]);
%             xb = linspace(-1*R, R, 100);
%             yb = zeros(1,size(xb,2));
%             for xc = 1:100
%                 yb(xc) = (R^2 - xb(xc)^2)^0.5;
%             end
%                 hold on; plot(xb, yb, '--w', 'LineWidth', 2.0);
%                 hold on; plot(xb, -1*yb, '--w', 'LineWidth', 2.0);
%             figure;
%                 plot(rin, ein(:,rinNum/2).*conj(ein(:,rinNum/2)));
%                 L1 = xlabel('y, \mum');
%                 L2 = ylabel('Intensity, a.u.');
%                 modiffigures(gcf, gca);
%                 modiflabels([L1 L2]);
%                 border = get(gca,'YLim');
%                 hold on; plot([-1*R -1*R], border, '--k');
%                 hold on; plot([R R], border, '--k');
%             figure;
%                 plot(rin, ein(rinNum/2,:).*conj(ein(rinNum/2,:)));
%                 L1 = xlabel('x, \mum');
%                 L2 = ylabel('Intensity, a.u.');
%                 modiffigures(gcf, gca);
%                 modiflabels([L1 L2]);
%                 border = get(gca,'YLim');
%                 hold on; plot([-1*R -1*R], border, '--k');
%                 hold on; plot([R R], border, '--k');
        end
    end
end

