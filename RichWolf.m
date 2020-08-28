classdef RichWolf
    %RICHWOLF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda;         % wavelength of the focused light (in mum)
        kWave;          % wavenumber
        NA;                 % numerical aperture of the lens
        refrMedia = 1.0;    % refraction index of media (air by default)
        thetaMax;           % maximal angle
        thetaMin = 0.0;     % minimal angle (zero by default)

        rMin = -2;
        rMax = 2;
        rNum = 200;
        rCoord;
        zMin = -20;
        zMax = 20;
        zNum = 200;
        zCoord;
        outFolder = strcat(['..\..\data\', date, '-', num2str(round(rand*1000)), '-']);
       
        FWHMx;
        FWHMy;
        FWHMz;
        HMA;
        focalVolume;
        
        % %%%%%%%% distribution in focus %%%%%%%%
        intensityFocus;
        intensityHFocus;
        
        intExFocus
        intEyFocus
        intEzFocus
        
        intHxFocus
        intHyFocus
        intHzFocus
        
        ampExFocus
        ampEyFocus
        ampEzFocus

        ampHxFocus
        ampHyFocus
        ampHzFocus        
        
        sxFocus
        syFocus
        szFocus
        % %%%%%%%% <- distribution in focus %%%%%%%%
        
        % %%%%%%%% distribution along Z (xz plane) %%%%%%%%
        intensityAlongZ;
        intensityHAlongZ;
        
        intExAlongZ
        intEyAlongZ
        intEzAlongZ
        
        intHxAlongZ
        intHyAlongZ
        intHzAlongZ
        
        ampExAlongZ
        ampEyAlongZ
        ampEzAlongZ

        ampHxAlongZ
        ampHyAlongZ
        ampHzAlongZ       
        
        sxAlongZ
        syAlongZ
        szAlongZ
        % %%%%%%%% <- distribution along Z (xz plane) %%%%%%%%
        
        storeFlag = false;
    end
    
    methods(Static)
        function out = calcSquare(fwhmx, fwhmy)
            out = pi * fwhmx * fwhmy /4;       % /4 because radiuses not diameters
        end
        function out = calcVolume(fwhmx, fwhmy, fwhmz)
            out = 4/3 * pi * fwhmx * fwhmy * fwhmz /8;
        end
        function fwhm = calcFWHM(arrayX, arrayY)
            fwhm = computefwhm(arrayX,arrayY,false);
        end
    end
    
    methods
        
        % constructor
        function obj = RichWolf()
            obj.rCoord = linspace(obj.rMin, obj.rMax, obj.rNum);    % in mum
            obj.zCoord = linspace(obj.zMin, obj.zMax, obj.zNum);    % in mum
        end
        % <- constructor
        
        % %%%%%%%%%%%%%%%%% setters %%%%%%%%%%%%%%%%%
        function obj = setThetaMin(obj, thetaMin)
            obj.thetaMin = thetaMin;
        end
        function obj = setNA(obj, NA, refrMedia)
            obj.NA = NA;
            obj.refrMedia = refrMedia;
            obj.thetaMax = asin(NA/refrMedia);
        end
        function obj = setThetaMax(obj, thetaMax, refrMedia)
            obj.thetaMax = thetaMax;
            obj.NA = refrMedia * sin(thetaMax);
            obj.refrMedia = refrMedia;
        end
        function obj = setRBorders(obj, rMin, rMax, rNum)
            obj.rMin = rMin;
            obj.rMax = rMax;
            obj.rNum = rNum;
            obj.rCoord = linspace(rMin, rMax, rNum);
        end
        function obj = setZBorders(obj, zMin, zMax, zNum)
            obj.zMin = zMin;
            obj.zMax = zMax;
            obj.zNum = zNum;
            obj.zCoord = linspace(zMin, zMax, zNum);
        end
        % <-- %%%%%%%%%%%%%%%%% setters %%%%%%%%%%%%%%%%%

        function obj = calcAmpExFocus(obj)
        end
        function obj = calcAmpEyFocus(obj)
        end
        function obj = calcAmpEzFocus(obj)
        end
        function obj = calcAmpHxFocus(obj)
        end
        function obj = calcAmpHyFocus(obj)
        end
        function obj = calcAmpHzFocus(obj)
        end
        
        function obj = calcAmpExAlongZ(obj)
        end
        function obj = calcAmpEyAlongZ(obj)
        end
        function obj = calcAmpEzAlongZ(obj)
        end
        function obj = calcAmpHxAlongZ(obj)
        end
        function obj = calcAmpHyAlongZ(obj)
        end
        function obj = calcAmpHzAlongZ(obj)
        end
        
        % calc intensity in focal point
        function obj = calcIntExFocus(obj)
            obj.intExFocus = obj.ampExFocus.*conj(obj.ampExFocus);
        end
        function obj = calcIntEyFocus(obj)
            obj.intEyFocus = obj.ampEyFocus.*conj(obj.ampEyFocus);
        end
        function obj = calcIntEzFocus(obj)
            obj.intEzFocus = obj.ampEzFocus.*conj(obj.ampEzFocus);
        end
        % <- calc intensity in focal point
        
        % calc intensity H in focal point
        function obj = calcIntHxFocus(obj)
            obj.intHxFocus = obj.ampHxFocus.*conj(obj.ampHxFocus);
        end
        function obj = calcIntHyFocus(obj)
            obj.intHyFocus = obj.ampHyFocus.*conj(obj.ampHyFocus);
        end
        function obj = calcIntHzFocus(obj)
            obj.intHzFocus = obj.ampHzFocus.*conj(obj.ampHzFocus);
        end        
        % <- calc intensity H in focal point
        
        % calc intensity along Z
        function obj = calcIntExAlongZ(obj)
            obj.intExAlongZ = obj.ampExAlongZ.*conj(obj.ampExAlongZ);
        end
        function obj = calcIntEyAlongZ(obj)
            obj.intEyAlongZ = obj.ampEyAlongZ.*conj(obj.ampEyAlongZ);
        end
        function obj = calcIntEzAlongZ(obj)
            obj.intEzAlongZ = obj.ampEzAlongZ.*conj(obj.ampEzAlongZ);
        end
        % <- calc intensity along Z
        
        % calc intensity H along Z
        function obj = calcIntHxAlongZ(obj)
            obj.intHxAlongZ = obj.ampHxAlongZ.*conj(obj.ampHxAlongZ);
        end
        function obj = calcIntHyAlongZ(obj)
            obj.intHyAlongZ = obj.ampHyAlongZ.*conj(obj.ampHyAlongZ);
        end
        function obj = calcIntHzAlongZ(obj)
            obj.intHzAlongZ = obj.ampHzAlongZ.*conj(obj.ampHzAlongZ);
        end        
        % <- calc intensity H along Z
        
        function obj = calcAllFocus(obj)
            fileID = fopen(strcat([obj.outFolder,'log.log']), 'a');
%             fprintf(fileID,'r.aCircE = %d + %d *i\r\n',  real(obj.aCircE), imag(obj.aCircE));
%             fprintf(fileID,'r.bCircE = %d + %d *i\r\n',  real(obj.bCircE), imag(obj.bCircE));
%             fprintf(fileID,'r.aCircH = %d + %d *i\r\n',  real(obj.aCircH), imag(obj.aCircH));
%             fprintf(fileID,'r.bCircH = %d + %d *i\r\n',  real(obj.bCircH), imag(obj.bCircH));
            fprintf(fileID,'r.mVortex = %d\r\n', obj.mVortex);
            fclose(fileID);
            
            obj = obj.calcAmpExFocus();
            obj = obj.calcAmpEyFocus();
            obj = obj.calcAmpEzFocus();
            
            obj = obj.calcAmpHxFocus();
            obj = obj.calcAmpHyFocus();
            obj = obj.calcAmpHzFocus();
            
            obj = obj.calcIntExFocus();
            obj = obj.calcIntEyFocus();
            obj = obj.calcIntEzFocus();

            obj = obj.calcIntHxFocus();
            obj = obj.calcIntHyFocus();
            obj = obj.calcIntHzFocus();
            
            obj.intensityFocus  = obj.intExFocus + obj.intEyFocus + obj.intEzFocus;
            obj.intensityHFocus = obj.intHxFocus + obj.intHyFocus + obj.intHzFocus;
            
            obj.sxFocus = real(obj.ampEyFocus.*conj(obj.ampHzFocus) - obj.ampEzFocus.*conj(obj.ampHyFocus));
            obj.syFocus = real(obj.ampEzFocus.*conj(obj.ampHxFocus) - obj.ampExFocus.*conj(obj.ampHzFocus));
            obj.szFocus = real(obj.ampExFocus.*conj(obj.ampHyFocus) - obj.ampEyFocus.*conj(obj.ampHxFocus));
        end
        
        function obj = calcAllAlongZ(obj)
            fileID = fopen(strcat([obj.outFolder,'log.log']), 'a');
%             fprintf(fileID,'r.aCircE = %d + %d *i\r\n',  real(obj.aCircE), imag(obj.aCircE));
%             fprintf(fileID,'r.bCircE = %d + %d *i\r\n',  real(obj.bCircE), imag(obj.bCircE));
%             fprintf(fileID,'r.aCircH = %d + %d *i\r\n',  real(obj.aCircH), imag(obj.aCircH));
%             fprintf(fileID,'r.bCircH = %d + %d *i\r\n',  real(obj.bCircH), imag(obj.bCircH));
            fprintf(fileID,'r.mVortex = %d\r\n', obj.mVortex);
            fclose(fileID);
            
            obj = obj.calcAmpExAlongZ();
            obj = obj.calcAmpEyAlongZ();
            obj = obj.calcAmpEzAlongZ();
            
            obj = obj.calcAmpHxAlongZ();
            obj = obj.calcAmpHyAlongZ();
            obj = obj.calcAmpHzAlongZ();
            
            obj = obj.calcIntExAlongZ();
            obj = obj.calcIntEyAlongZ();
            obj = obj.calcIntEzAlongZ();

            obj = obj.calcIntHxAlongZ();
            obj = obj.calcIntHyAlongZ();
            obj = obj.calcIntHzAlongZ();
            
            obj.intensityAlongZ  = obj.intExAlongZ + obj.intEyAlongZ + obj.intEzAlongZ;
            obj.intensityHAlongZ = obj.intHxAlongZ + obj.intHyAlongZ + obj.intHzAlongZ;
            
            obj.sxAlongZ = real(obj.ampEyAlongZ.*conj(obj.ampHzAlongZ) - obj.ampEzAlongZ.*conj(obj.ampHyAlongZ));
            obj.syAlongZ = real(obj.ampEzAlongZ.*conj(obj.ampHxAlongZ) - obj.ampExAlongZ.*conj(obj.ampHzAlongZ));
            obj.szAlongZ = real(obj.ampExAlongZ.*conj(obj.ampHyAlongZ) - obj.ampEyAlongZ.*conj(obj.ampHxAlongZ));
        end
        
        % plotters
        % plot in focus
        function obj = plotIntensityFocus(obj, componentsFlag)
            figure;
                imagesc(obj.rCoord, obj.rCoord, obj.intensityFocus)
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('x, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);            
                if obj.storeFlag
                    temp = [0 obj.rCoord; obj.rCoord' obj.intensityFocus];
                    save(strcat([obj.outFolder,'intensityFocus.dat']),'temp','-ascii')
                end
                
            if componentsFlag
                figure;
                    imagesc(obj.rCoord, obj.rCoord, obj.intExFocus)
                    colorbar;
                    colormap('hot');
                    axis xy;
                    L1 = xlabel('x, \mum');
                    L2 = ylabel('y, \mum');
                    modiffigures(gcf, gca);
                    modiflabels([L1 L2]);
                    if obj.storeFlag
                        temp = [0 obj.rCoord; obj.rCoord' obj.intExFocus + obj.intEyFocus];
                        save(strcat([obj.outFolder,'erFocus.dat']),'temp','-ascii')
                    end
                figure;
                    imagesc(obj.rCoord, obj.rCoord, obj.intEyFocus)
                    colorbar;
                    colormap('hot');
                    axis xy;
                    L1 = xlabel('x, \mum');
                    L2 = ylabel('y, \mum');
                    modiffigures(gcf, gca);
                    modiflabels([L1 L2]);
                    if obj.storeFlag
                        temp = [0 obj.rCoord; obj.rCoord' obj.intExFocus + obj.intEyFocus];
                        save(strcat([obj.outFolder,'erFocus.dat']),'temp','-ascii')
                    end
                figure;
                    imagesc(obj.rCoord, obj.rCoord, obj.intEzFocus)
                    colorbar;
                    colormap('hot');
                    axis xy;
                    L1 = xlabel('x, \mum');
                    L2 = ylabel('y, \mum');
                    modiffigures(gcf, gca);
                    modiflabels([L1 L2]);
                    if obj.storeFlag   
                        temp = [0 obj.rCoord; obj.rCoord' obj.intEzFocus];
                        save(strcat([obj.outFolder,'ezFocus.dat']),'temp','-ascii')
                    end
            end
        end
        function obj = plotIntensityHFocus(obj, componentsFlag)
            figure;
                imagesc(obj.rCoord, obj.rCoord, obj.intensityHFocus)
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('x, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);            
                if obj.storeFlag
                    temp = [0 obj.rCoord; obj.rCoord' obj.intensityHFocus];
                    save(strcat([obj.outFolder,'intensityHFocus.dat']),'temp','-ascii')
                end
                
            if componentsFlag
                figure;
                    imagesc(obj.rCoord, obj.rCoord, obj.intHxFocus)
                    colorbar;
                    colormap('hot');
                    axis xy;
                    L1 = xlabel('x, \mum');
                    L2 = ylabel('y, \mum');
                    modiffigures(gcf, gca);
                    modiflabels([L1 L2]);
                    if obj.storeFlag
                        temp = [0 obj.rCoord; obj.rCoord' obj.intHxFocus + obj.intHyFocus];
                        save(strcat([obj.outFolder,'hrFocus.dat']),'temp','-ascii')
                    end
                figure;
                    imagesc(obj.rCoord, obj.rCoord, obj.intHyFocus)
                    colorbar;
                    colormap('hot');
                    axis xy;
                    L1 = xlabel('x, \mum');
                    L2 = ylabel('y, \mum');
                    modiffigures(gcf, gca);
                    modiflabels([L1 L2]);
                    if obj.storeFlag
                        temp = [0 obj.rCoord; obj.rCoord' obj.intHxFocus + obj.intHyFocus];
                        save(strcat([obj.outFolder,'hrFocus.dat']),'temp','-ascii')
                    end   
                figure;
                    imagesc(obj.rCoord, obj.rCoord, obj.intHzFocus)
                    colorbar;
                    colormap('hot');
                    axis xy;
                    L1 = xlabel('x, \mum');
                    L2 = ylabel('y, \mum');
                    modiffigures(gcf, gca);
                    modiflabels([L1 L2]);
                    if obj.storeFlag   
                        temp = [0 obj.rCoord; obj.rCoord' obj.intHzFocus];
                        save(strcat([obj.outFolder,'hzFocus.dat']),'temp','-ascii')
                    end
            end
                
        end
        function obj = plotPoyntingFocus(obj)
            figure;
                imagesc(obj.rCoord, obj.rCoord, real(obj.szFocus))
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('x, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);
                if obj.storeFlag
                    temp = [0 obj.rCoord; obj.rCoord' real(obj.szFocus)];
                    save(strcat([obj.outFolder,'szFocus.dat']),'temp','-ascii')    
                end
            figure;
                imagesc(obj.rCoord, obj.rCoord, real(obj.sxFocus))
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('x, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);
                if obj.storeFlag
                    temp = [0 obj.rCoord; obj.rCoord' real(obj.sxFocus)];
                    save(strcat([obj.outFolder,'sxFocus.dat']),'temp','-ascii')    
                end   
            figure;
                imagesc(obj.rCoord, obj.rCoord, real(obj.syFocus))
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('x, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);
                if obj.storeFlag
                    temp = [0 obj.rCoord; obj.rCoord' real(obj.syFocus)];
                    save(strcat([obj.outFolder,'syFocus.dat']),'temp','-ascii')    
                end              
        end

        function obj = plotAllFocus(obj)
            
            obj = obj.calcAllFocus();
            
            obj = obj.plotIntensityFocus(true);
            obj = obj.plotIntensityHFocus(true);
            obj = obj.plotPoyntingFocus();
 
        end  
        % <-- plot in focus
        
        % plot along Z
        function obj = plotAllAlongZ(obj)
            
            obj = obj.calcAllAlongZ();

            obj = obj.plotIntensityAlongZ(true);
            obj = obj.plotIntensityHAlongZ(true);
            obj = obj.plotPoyntingAlongZ();

        end  
        
        function obj = plotIntensityAlongZ(obj, componentsFlag)
            
            figure;
                imagesc(obj.zCoord, obj.rCoord, obj.intensityAlongZ)
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('z, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);            
                if obj.storeFlag
%                     temp = [0 obj.rCoord; obj.zCoord' obj.intensityAlongZ];
%                     save(strcat([obj.outFolder,'intensityAlongZ.dat']),'temp','-ascii')
                end   
                
                if componentsFlag
                
                    figure;
                        imagesc(obj.zCoord, obj.rCoord, obj.intExAlongZ + obj.intEyAlongZ)
                        colorbar;
                        colormap('hot');
                        axis xy;
                        L1 = xlabel('z, \mum');
                        L2 = ylabel('y, \mum');
                        modiffigures(gcf, gca);
                        modiflabels([L1 L2]);
                        if obj.storeFlag
        %                     temp = [0 obj.rCoord; obj.zCoord' obj.intExAlongZ + obj.intEyAlongZ];
        %                     save(strcat([obj.outFolder,'erAlongZ.dat']),'temp','-ascii')
                        end
                    figure;
                        imagesc(obj.zCoord, obj.rCoord, obj.intEzAlongZ)
                        colorbar;
                        colormap('hot');
                        axis xy;
                        L1 = xlabel('z, \mum');
                        L2 = ylabel('y, \mum');
                        modiffigures(gcf, gca);
                        modiflabels([L1 L2]);
                        if obj.storeFlag   
        %                     temp = [0 obj.rCoord; obj.zCoord' obj.intEzAlongZ];
        %                     save(strcat([obj.outFolder,'ezAlongZ.dat']),'temp','-ascii')
                        end
                end
                
        end
        function obj = plotIntensityHAlongZ(obj, componentsFlag)
            figure;
                imagesc(obj.zCoord, obj.rCoord, obj.intensityHAlongZ)
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('z, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);            
                if obj.storeFlag
%                     temp = [0 obj.rCoord; obj.zCoord' obj.intensityHAlongZ];
%                     save(strcat([obj.outFolder,'intensityHAlongZ.dat']),'temp','-ascii')
                end            
                
                if componentsFlag
                    
                    figure;
                        imagesc(obj.zCoord, obj.rCoord, obj.intHxAlongZ+obj.intHyAlongZ)
                        colorbar;
                        colormap('hot');
                        axis xy;
                        L1 = xlabel('z, \mum');
                        L2 = ylabel('y, \mum');
                        modiffigures(gcf, gca);
                        modiflabels([L1 L2]);            
                        if obj.storeFlag
        %                     temp = [0 obj.rCoord; obj.zCoord' obj.intensityHAlongZ];
        %                     save(strcat([obj.outFolder,'intensityHAlongZ.dat']),'temp','-ascii')
                        end 
                                    figure;
                        imagesc(obj.zCoord, obj.rCoord, obj.intHzAlongZ)
                        colorbar;
                        colormap('hot');
                        axis xy;
                        L1 = xlabel('z, \mum');
                        L2 = ylabel('y, \mum');
                        modiffigures(gcf, gca);
                        modiflabels([L1 L2]);            
                        if obj.storeFlag
        %                     temp = [0 obj.rCoord; obj.zCoord' obj.intensityHAlongZ];
        %                     save(strcat([obj.outFolder,'intensityHAlongZ.dat']),'temp','-ascii')
                        end 
                    
                    
                end
        end
        function obj = plotPoyntingAlongZ(obj)
            figure;
                imagesc(obj.zCoord, obj.rCoord, real(obj.szAlongZ))
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('z, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);
                if obj.storeFlag
%                     temp = [0 obj.rCoord; obj.zCoord' real(obj.szAlongZ)];
%                     save(strcat([obj.outFolder,'szAlongZ.dat']),'temp','-ascii')    
                end      
            figure;
                imagesc(obj.zCoord, obj.rCoord, real(obj.sxAlongZ))
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('z, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);
                if obj.storeFlag
%                     temp = [0 obj.rCoord; obj.zCoord' real(obj.szAlongZ)];
%                     save(strcat([obj.outFolder,'szAlongZ.dat']),'temp','-ascii')    
                end 
            figure;
                imagesc(obj.zCoord, obj.rCoord, real(obj.syAlongZ))
                colorbar;
                colormap('hot');
                axis xy;
                L1 = xlabel('z, \mum');
                L2 = ylabel('y, \mum');
                modiffigures(gcf, gca);
                modiflabels([L1 L2]);
                if obj.storeFlag
%                     temp = [0 obj.rCoord; obj.zCoord' real(obj.szAlongZ)];
%                     save(strcat([obj.outFolder,'szAlongZ.dat']),'temp','-ascii')    
                end                 
                
                
        end        
        
        % <-- plot along Z
        % <-- plotters 
    end
    
end