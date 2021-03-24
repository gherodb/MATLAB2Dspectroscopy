%% This code is written by Soroush Khosravi and Daniel Busa to simulate the Hetrodyne detection method
%% and also examine our 2D analysis code performance in retrieving the signal.

%% Read Files
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
clear;
lambda=dlmread('.\filessim\lambda02fs.txt'); %read wavelenght axis (nm)
tau=dlmread('.\filessim\tau02fs.txt'); %read coherence time axis (nm)
cntlambdatau=dlmread('.\filessim\spectra02fs.txt'); %read raw experimental data from detector
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Constants 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
    c= 2.99792458*(10.^8); % m./s.
    f0=0.5765;             % PHz.
    w0=2*pi*f0;            % from f0
    T=35;                  % fs.
    npix=length(lambda);
    ntau=length(tau);
    tau0indx=ceil(ntau/2);
    ln1psqrt2=log(1+sqrt(2))/log(exp(1)); % factor
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
%% Axis and Conversions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
    cnttaulambda=transpose(cntlambdatau);
    cnttau0lambda=cnttaulambda(tau0indx,:); % spectral interferogram at tau=0
    cnttau0freq=fliplr(cnttau0lambda); % to go from wavelength to frequency
%    
    fc(1:npix)=0; % create frequency matrix (1*1600)
    for kt=1:npix
    fc(kt)=(c*10^-6)/(lambda(npix+1-kt)); % set frequency axis from small to large PHz
    end
%
    dfc=(fc(npix)-fc(1))/((npix-1));
    dfsc=1;
    df=(dfsc)*dfc; % PHz.
%   
    Nsc=16.4267;
%     Nsc=10; % scale number of pixels by factor of 10
    N1=ceil(((fc(npix))*(npix))/(fc(npix)-fc(1))); 
    N=Nsc*N1; 
%    
    dt=1/(N*df); % fs. (NOTE: 'dt' is scaled inversely by 'Nsc' factor)
%
    f=0:df:(N-1)*df; % PHz. Full 'f' (frequency) "axis" created
    w=2*pi*f; % 'w' (angular frequency) axis created from 'f' for pulse equations
%
    t=(((1-N)/2)*dt):dt:(((N-1)/2)*dt); % fs. Full time domain "axis" created
%    
    pix=1:1:npix; % Create "axis" based on camera pixels
    pixN=1:1:N; % Create "axis" based on "theoretical camera"
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Phase
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
    phisig=0;             %temporal phase sig (NOT USED)
    philo=0;              %temporal phase lo  (NOT USED)
%    
%   phase at sample position
    sphsigS0=0; sphsigS1=0; sphsigS2=0; sphsigS3=0; sphsigS4=0;
    sphloS0=0; sphloS1=600; sphloS2=0; sphloS3=0; sphloS4=0;
%
%   phase due to path
    sphP0=0; sphP1=0; sphP2=0; sphP3=0; sphP4=0;
%
%   total phase of sig
    sphsigT0=sphsigS0+sphP0;
    sphsigT1=sphsigS1+sphP1;
    sphsigT2=sphsigS2+sphP2;
    sphsigT3=sphsigS3+sphP3;
    sphsigT4=sphsigS4+sphP4;
%
%   total phase of lo
    sphloT0=sphloS0+sphP0;
    sphloT1=sphloS1+sphP1;
    sphloT2=sphloS2+sphP2;
    sphloT3=sphloS3+sphP3;
    sphloT4=sphloS4+sphP4;
%    
%   full ph(w) for sig and lo at sample position
    sphsigS=(sphsigS0)+(sphsigS1.*(w-w0))+(sphsigS2.*(((w-w0).^2)./2))+(sphsigS3.*(((w-w0).^3)./6))+(sphsigS4.*(((w-w0).^4)./24));
    sphloS=(sphloS0)+(sphloS1.*(w-w0))+(sphloS2.*(((w-w0).^2)./2))+(sphloS3.*(((w-w0).^3)./6))+(sphloS4.*(((w-w0).^4)./24));
%
%   full ph(w) for sig and lo at camera position
    sphsigT=(sphsigT0)+(sphsigT1.*(w-w0))+(sphsigT2.*(((w-w0).^2)./2))+(sphsigT3.*(((w-w0).^3)./6))+(sphsigT4.*(((w-w0).^4)./24));
    sphloT=(sphloT0)+(sphloT1.*(w-w0))+(sphloT2.*(((w-w0).^2)./2))+(sphloT3.*(((w-w0).^3)./6))+(sphloT4.*(((w-w0).^4)./24));
%  
%   correction factor
    sphloCOR=(sphloT0)+((sphloT1-sphloS1).*(w-w0))+(sphloT2.*(((w-w0).^2)./2))+(sphloT3.*(((w-w0).^3)./6))+(sphloT4.*(((w-w0).^4)./24));
%   
    sphcorr=sphloCOR-(sphloT-sphloS);
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Spectral Interferogram Simulation A: Building Pulse Shapes 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   Envelope equations for simulated laser pulses
    AGf=sqrt(pi/(2*0.69314718))*T*exp(-1.*((T^2)*power((w-w0),2))/(8*0.69314718)); % Gaussian in frequency domain
    AGt=exp((-2.*0.69314718).*((power(t,2)))/((T)^2)); % Gaussian in time domain
    AHt=sech(2*ln1psqrt2*(t/T)); % sech^2 in time domain
%
    GSIGf0=AGf.*exp(-1i*w*0); %GSIGf0=AGf.*exp(-1i*w*sphsigS1);
%
%   Signal equations with oscillations
    GSIGt0=AGt.*exp(1i*w0*(t));
    HSIGt0=AHt.*exp(1i*w0*(t));  
    GLOt=AGt.*exp(1i*w0*t);
%
%   Find the indices marker for the midpoint of the Gaussian pulse
    [maxAGt,ind]=max(abs(AGt));
    midP=(ind);
%
%   Create three-part piecewise pulse equation (matrix) on 'N'
    ESIGtP(1:N)=0;
%   Cubic (growth)
    for ip=midP-1000:midP-100
        ESIGtP(ip)=((ip-(midP-1000))/(900))^3;
    end
%   Constant    
    for ipp=midP-100:midP+100
        ESIGtP(ipp)=1;
    end
%   Exponential (decay)    
    for ippp=midP+100:N
        ESIGtP(ippp)=(exp(-(ippp-(midP+100))/150));
    end
%
%   Create a half Gaussian, half sech^2 pulse equation (matrix) on 'N'
    ESIGtW(1:N)=0;      
%   First half (to midpoint index)
    for ih=1:midP
        ESIGtW(ih)=GSIGt0(ih);
    end
%   Second half (from midpoint index)
    for ihh=midP:N
        ESIGtW(ihh)=HSIGt0(ihh);
    end
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Choosing Pulse Shape
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
    choice=4;
%
%   Loop to decide pulse signal equation
    if choice==1;
        ESIGt=0;
        disp('Using frequency-domain equations');
    elseif choice==2;
        ESIGt=ESIGtP.*exp(1i*w0*t);
        disp('Using piecewise built time-domain function...')
    elseif choice==3;
        ESIGt=ESIGtW;    
        disp('Using 1/2 Gaussian 1/2 sech^2 time-domain function...')
    elseif choice==4;
        ESIGt=GSIGt0;
        disp('Strictly Gaussian Pulse time-domain equation...')    
    else
        disp('INVALID CHOICE, TRY AGAIN...')
    end
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Spectral Interferogram Simulation B: Applying Phase
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   Frequency-domain equations
%
%   Applying phase
    GSIGfS=GSIGf0.*exp(-1i.*sphsigS); % signal pulse with only sample phase
    GSIGfC=GSIGf0.*exp(-1i.*sphsigT); % signal pulse with total phase
    GLOfC=GSIGf0.*exp(-1i.*sphloT); % LO pulse with total phase
%
%
%   Time-domain equations
%
%   Transform time-domain equations to frequency domain
    FTGSIGt=fft(fftshift(ESIGt)); % Signal
    FTGLOt=fft(fftshift(GLOt)); % LO
%
%   Applying phase
    FTGSIGtS=FTGSIGt.*exp(-1i.*sphsigS); % signal pulse with only sample phase
    FTGSIGtC=FTGSIGt.*exp(-1i.*sphsigT); % signal pulse with total phase
    FTGLOtC=FTGLOt.*exp(-1i.*sphloT); % LO pulse with total phase
%
%   Transform phased signals from frequency-domain to time-domain
%
%   (Equations originally in frequency-domain)
    FTGSIGf0=ifftshift(ifft(GSIGf0)); % signal pulse with no phase
    FTGSIGfS=ifftshift(ifft(GSIGfS)); % signal pulse with phase @sample
    FTGSIGfC=ifftshift(ifft(GSIGfC)); % signal pulse with phase @camera
    FTGLOfC=ifftshift(ifft(GLOfC)); % LO pulse with phase @camera
%    
%   (Equations originally in time-domain)
    ESIGt0=ESIGt.*exp(1i*w0*t); % signal pulse with no phase
    iFTFTGSIGtS=ifftshift(ifft(FTGSIGtS)); % signal pulse with phase @sample
    iFTFTGSIGtC=ifftshift(ifft(FTGSIGtC)); % signal pulse with phase @camera
    iFTFTGLOtC=ifftshift(ifft(FTGLOtC)); % LO pulse with phase @camera
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Spectral Interferogram Simulation C: Building Spectral Interferogram
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   Loop to determine signal addition, display message indicating choice
    if choice==1
        Esumt=FTGSIGfC+FTGLOfC; % using frequency-domain equations
    else
        Esumt=iFTFTGSIGtC+iFTFTGLOtC; % using time-domain equations
    end
%
    FTEsumt=fft(Esumt); % Fourier transform simulates grating
%    
    IFTEsumt=FTEsumt.*conj(FTEsumt); % Multiplying by conjugate simulates camera reading
%
%   Files written from equations (frequency-domain equations only)
    dlmwrite('Signal_1_freq.txt',GSIGfS);
    dlmwrite('Signal_1_time.txt',FTGSIGfS);
    dlmwrite('LO_freq.txt',GLOfC);
    dlmwrite('LO_time.txt',FTGLOfC);
    dlmwrite('Spectral_Interferogram_1_freq.txt',IFTEsumt);     
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Amplitude of LO 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%    
%   Change zero values to 1 to avoid infinity when dividing the term in 2D analysis
    if choice==1
    ALOf=sqrt(GLOfC.*conj(GLOfC)); % Amplitude of LO in frequency-domain (from frequency-domain equations)
    else
    ALOf=sqrt(FTGLOtC.*conj(FTGLOtC)); % Amplitude of LO in frequency-domain (from time-domain equations)
    end
%    
%   Finds the first and last values greater than 1 in LO amplitude function
%   Lowcut and highcut indices are used to rewrite the LO amplitude matrix 
    AGtFIND=abs(ALOf);
    [AGto1,ind]=find(AGtFIND>1,1,'first');
    Alowcut=pixN(ind);
    [AGto2,ind]=find(AGtFIND>1,1,'last');
    Ahighcut=pixN(ind);
%
    for kii=1:Alowcut
        ALOf(kii)=1;
    end
    for kiii=Ahighcut:N
        ALOf(kiii)=1;
    end
%    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% 2D Code Begins
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Plotting [1]: Pulse Equations and Spectral Interferogram
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   Plot LO amplitude with adjusted end values
    figure
    plot(pixN,ALOf,'g');
    title('Infinity Avoidance by Changing Small Values beyond the edges to 1');
    xlim([Alowcut-1000 Ahighcut+1000]);
    xlabel('pixN')
    line([Alowcut,Alowcut],[0,200]);
    line([Ahighcut,Ahighcut],[0,200]);
%    
%   Plot real parts of signal and LO pulses on time-domain and frequency-domain
    if choice==1 % using frequency-domain equations
        figure
        subplot(2,1,1)
        plot(t,real(FTGSIGfC),t,real(FTGLOfC),'g');
        legend('Signal [real] from equation','LO [real] from equation');
        title('Signal and LO electric field in time domain');
        xlim([-100 1100]);
        xlabel('Time in fs.')
        subplot(2,1,2)
        plot(f,real(GSIGfC),f,real(GLOfC),'g');
        title('Signal and LO electric field in frequency domain');
        legend('Signal [real] from equation','LO [real] from equation');
        xlim([f0-0.03 f0+0.03])
        xlabel('Frequency in PHz.')    
    else % using time-domain equations
        figure
        subplot(2,1,1)
        plot(t,real(iFTFTGSIGtC),t,real(iFTFTGLOtC),'g');
        legend('Signal [real] from equation','LO [real] from equation');
        title('Signal and LO electric field in time domain');
        xlim([-100 1100]);
        xlabel('Time in fs.')
        subplot(2,1,2)
        plot(f,real(FTGSIGtC),f,real(FTGLOtC),'g');
        title('Signal and LO electric field in frequency domain');
        legend('Signal [real] from equation','LO [real] from equation');
        xlim([f0-0.03 f0+0.03])
        xlabel('Frequency in PHz.')    
    end
%    
%   Plot simulated interferogram versus experimental interferogram
    figure
    plot(f,IFTEsumt,'b',fc,(1)*cnttau0freq,'k');
    legend('Simulated Spectral Interferogram (Intensity)','Rhodamine Experimental Spectral Interferogram');
    title('Compare Experimental data and Simulation');
    xlim([f0-0.03 f0+0.03])
    xlabel('Frequency in PHz.')    
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% 2D code 1st step (ifft and shift)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%    shiftifftSI=ifftshift(ifft(IFTEsumt)); % Fourier transform to time domain
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% 2D code 2nd step (Window in time domain to get rid of mirror and zero frequency component)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
    winshiftifftSI=shiftifftSI;
%    
%   Finds endpoints and midpoint indices on 'pixN' and time-domain
    winshiftFIND=abs(winshiftifftSI);
    [maxwinshift,ind]=max(winshiftFIND);
    piwinMAX=pixN(ind);
    twinMAX=t(ind);
    [wino1,ind]=find(winshiftFIND>1,1,'first');
    piwinFIND1=pixN(ind);
    twinFIND1=t(ind);
    [wino2,ind]=find(winshiftFIND>1,1,'last');
    piwinFIND2=pixN(ind);
    twinFIND2=t(ind);
%
%   Define high and low cuts based on endpoints and maximum    
    WINlowcut00=ceil(piwinFIND2-piwinMAX)/2+piwinMAX;
    WINhighcut00=ceil((1.005)*piwinFIND2);
    WINlowcut=ceil((0.99)*piwinFIND1);
    WINhighcut=piwinMAX-ceil((piwinMAX-piwinFIND1)/2.2);
%   
%   Overwrite extraneous matrix values to '0'
    for kwin=1:WINlowcut
    winshiftifftSI(kwin)=0;
    end
    for kwin2=WINhighcut:N
    winshiftifftSI(kwin2)=0;
    end
%
%   Determine midpoint index of isolated pulse component
    AFTERwinshiftFIND=abs(winshiftifftSI);   
    [maxAFTERwinshift,ind]=max(AFTERwinshiftFIND);
    piwinMAX2=pixN(ind);
    twinMAX2=t(ind);
    timediff=twinMAX-twinMAX2;
%   
%   Create piecewise equation of line between maxima
    wLINE(1:N)=0;
    for iw=piwinMAX2:piwinMAX
        wLINE(iw)=maxwinshift;
    end
    
%%  Plotting [2]
%   Plot spectral interferogram with adjustments on time-domain
    figure
    subplot(3,1,1)
    plot(real(shiftifftSI));
    title('Spectral Interferogram SI(t) [real] in time-domain');
    line([WINlowcut,WINlowcut],[0,maxwinshift*1.25]);
    line([WINhighcut,WINhighcut],[0,maxwinshift*1.25]);
    xlim([WINlowcut-1000 WINhighcut00+1000])
    xlabel('pixN')
    subplot(3,1,2)
    plot(t,wLINE,'r',t,real(shiftifftSI));
    title(['width = ',num2str(timediff),' fs']);
    legend('Width line between maxima','SI(t) [real] in time-domain')
    line([twinMAX2,twinMAX2],[0,maxwinshift*1.25]);
    line([twinMAX,twinMAX],[0,maxwinshift*1.25]);
    xlim([twinMAX2-500 twinMAX+500])
    subplot(3,1,3)
    plot(t,real(winshiftifftSI),'b');
    title({'SI(t)after isolation adjustment','Esig(t)*Elo(t)*exp(i*w*DELAY)'});
    xlim([twinFIND1-500 twinFIND2+500]);
    xlabel('Time fs.');
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% 2D code 3rd step (fft and fftshift back to frequency domain)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
    fftshiftwinshiftifftSI=fft(fftshift(winshiftifftSI)); % Transform newly adjusted SI(t) to frequency domain
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% 2D code last step (Phasing and Removing Delay)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   Define 'delay' equation
    delay(1:N)=0;
    for kff=1:npix 
    delay(kff+N1-npix)=(exp(-1i*(((2*pi)*f(kff+N1-npix))-w0).*sphloS1)); %2*pi
    end
%    
%   Apply delay, LO amplitude, and phase corrections
    delayfftshiftwinshiftifftSI=(fftshiftwinshiftifftSI).*(delay); % removes DELAY factor
    ampdelayfftshiftwinshiftifftSI=delayfftshiftwinshiftifftSI./ALOf; % removes LO amplitude factor
    phiampdelayfftshiftwinshiftifftSI=ampdelayfftshiftwinshiftifftSI.*exp(-1i*sphcorr); % removes phasing factor
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Plotting [3]
%   Plot adjusted SI(w)
    figure
    subplot(3,1,1)
    plot(f,real(fftshiftwinshiftifftSI),'r');
    title({'Simulated SI(w) [real] after isolation-adjustment','Esig(w)*Elo(w)*exp(i*w*DELAY)*exp(i*phi(w))'});
    xlim([f0-0.03 f0+0.03])   
    subplot(3,1,2)
    plot(f,real(delay),'--m');
    title('delay [real]');
    xlim([f0-0.03 f0+0.03])
    subplot(3,1,3)
    plot(f,real(delayfftshiftwinshiftifftSI),'--r');
    title({'Simulated SI(w) [real] after multiplying by delay factor','Esig(w)*Elo(w)*exp(i*phi(w))'});
    xlim([f0-0.03 f0+0.03])
%     
    figure
    plot(f,(0.3)*real(IFTEsumt),'k',f,real(delayfftshiftwinshiftifftSI),'-.g',f,(300)*real(ampdelayfftshiftwinshiftifftSI),'--m',f,(250)*real(phiampdelayfftshiftwinshiftifftSI),'r');
    legend('(0.3)*Simulated SI(w)','SI(w)*delay [real]','(300)*SI(w).*delay./LOamplitude [real]','(250)*SI(w).*delay./LOamplitude*phasing [real]');
    title('Comparing adjustments');
    xlim([f0-0.03 f0+0.03])
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
%% Plotting [4]
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   Plot and compare retrieved signal against original simulated signal in frequency-domain
    if choice==1 % using frequency-domain equations
    figure
    subplot(1,3,1)
    plot(f,(1)*real(GSIGf0),'k',f,(1)*real(GSIGfS),'m',f,(1).*real(phiampdelayfftshiftwinshiftifftSI),'--r');
    legend('Simulated pulse without phasing','Simulated pulse with phasing','Retrieved signal from 2D Code');
    title('[real]')
    xlim([f0-0.03 f0+0.03])
    subplot(1,3,2)
    plot(f,(1)*imag(GSIGf0),'k',f,(1)*imag(GSIGfS),'m',f,(1).*imag(phiampdelayfftshiftwinshiftifftSI),'--r');
    title('[imaginary]')
    xlim([f0-0.03 f0+0.03])
    subplot(1,3,3)
    plot(f,(1)*abs(GSIGf0),'k',f,(1)*abs(GSIGfS),'m',f,(1).*abs(phiampdelayfftshiftwinshiftifftSI),'--r');
    title('[abs]')
    xlim([f0-0.03 f0+0.03])
    xlabel('Frequency PHz.')
    else % using time-domain equations
    figure
    subplot(1,3,1)
    plot(f,(1)*real(FTGSIGt),'k',f,(1)*real(FTGSIGtS),'m',f,(1).*real(phiampdelayfftshiftwinshiftifftSI),'--r');
    legend('Simulated pulse without phasing','Simulated pulse with phasing','Retrieved signal from 2D Code');
    title('[real]')
    xlim([f0-0.03 f0+0.03])
    subplot(1,3,2)
    plot(f,(1)*imag(FTGSIGt),'k',f,(1)*imag(FTGSIGtS),'m',f,(1).*imag(phiampdelayfftshiftwinshiftifftSI),'--r');
    title('[imaginary]')
    xlim([f0-0.03 f0+0.03])
    subplot(1,3,3)
    plot(f,(1)*abs(FTGSIGt),'k',f,(1)*abs(FTGSIGtS),'m',f,(1).*abs(phiampdelayfftshiftwinshiftifftSI),'--r');
    title('[abs]')
    xlim([f0-0.03 f0+0.03])
    xlabel('Frequency PHz.')
    end
%    
%   Transform retrieved signal to time-domain
    shiftifftampphidelayfftshiftwinshiftifftSI2=ifftshift(ifft(phiampdelayfftshiftwinshiftifftSI));
%    
%   Plot and compare retrieved signal against original simulated signal in time-domain
    if choice==1 % using frequency-domain equations
    figure
    plot(t,real(FTGSIGf0),'k',t,real(FTGSIGfS),'c',t,real(shiftifftampphidelayfftshiftwinshiftifftSI2),'--b');
    legend('phi=1 signal with specphi=0','Simulated signal at sample position','signal from 2D code in time domain');
    xlim([-200 200]);
    xlabel('Time fs.')
    else % using time-domain equations
    figure
    plot(t,real(ESIGt0),'k',t,real(iFTFTGSIGtS),'c',t,(1)*real(shiftifftampphidelayfftshiftwinshiftifftSI2),'--b');
    legend('signal with specphi=0','Simulated signal at sample position','signal from 2D code in time domain');
    xlim([-200 200]);
    xlabel('Time fs.')
    end
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %