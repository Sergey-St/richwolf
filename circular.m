classdef circular < debayCircular
    %ABGMODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods

        function obj = circular(lambda, thetaMin, NA)
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

        end
    end
end

