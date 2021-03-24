
%% Read files I
clear;

lambda=dlmread('.\filessim\lambda02fs.txt'); %read wavelenght axis (nm)
%     cntlambdatau=dlmread('.\filessim\spectra03fs.txt'); %read raw data from detector
%     cnttaulambda=transpose(cntlambdatau);
%     cnttaulambda=fliplr(cnttaulambda);
%     cnttau0lambda=cnttaulambda(601,:);


%% Constants 
    c= 2.99792458*(10.^8); %m/s
    f0=0.576;
    w0=2*pi*f0;
    T=35;
    Tj=600;
    phi=0;
    D1=200;
    D2=D1+Tj;
    npix=length(lambda);
    
    
%% Axis
    fc(1:npix)=0;
    for kt=1:npix
    fc(kt)=(c*10^-6)/(lambda(npix+1-kt)); %set frequency axis from small to large PHz
    end

    Nsc=1;
    N=ceil((Nsc)*(fc(npix)*(npix))/(fc(npix)-fc(1))); 
    h_N=ceil((N/Nsc));
  
    dfc=(fc(npix)-fc(1))/((npix-1));
    dfsc=1;
    df=(dfsc)*dfc; %PHz
    dtsc=1;
    dt=1/(dtsc*h_N*df); %fs
%
% We can now define our frequency and time domains:
%
    f=0:df:(h_N-1)*df; %PHz
    w=2*pi*f;
%
    t=(((1-h_N)/2)*dt):dt:(((h_N-1)/2)*dt); %fs
    
    pix=1:1:npix;
    pixN=1:1:N;

    %% Spectral Interferogram Simulation
    
    GSIGt=exp(-1.*(power((t-D1),2))/(2*(T)^2)).*exp(1i*w0*(t-D1)).*exp(-1i*phi);
    GLOt=exp(-1.*(power((t-D2),2))/(2*(T)^2)).*exp(1i*w0*(t-D2)).*exp(-1i*phi);
    GSIGtTj=exp(-1.*(power((t-Tj),2))/(2*(T)^2)).*exp(1i*w0*(t-Tj)).*exp(-1i*phi);

    
    EsumGt=GSIGt+GLOt;
    FTEsumGt=fft((EsumGt)); %Grating
%     FTEsumGt=fft(fftshift(EsumGt)); %Grating
    FTGSIGtTj=fft(fftshift(GSIGtTj));
    
    IFTEsumGt=FTEsumGt.*conj(FTEsumGt); % Detector
    
    HSIGt=sech(0.881374*((t-D1)/T).^2).*exp(1i*w0*(t-D1)).*exp(-1i*phi);
    HLOt=sech(0.881374*((t-D2)/T).^2).*exp(1i*w0*(t-D2)).*exp(-1i*phi);
    HSIGtTj=sech(0.881374*((t-Tj)/T).^2).*exp(1i*w0*(t-Tj)).*exp(-1i*phi);
 
    EsumHt=HSIGt+HLOt;
    FTEsumHt=fft((EsumHt)); %Grating
%     FTEsumHt=fft(fftshift(EsumHt)); %Grating
    FTHSIGtTj=fft(fftshift(HSIGtTj));
    
    IFTEsumHt=FTEsumHt.*conj(FTEsumHt); % Detector
    
 

    dlmwrite('.\filessim\spectra03fs.txt',IFTEsumGt);
    dlmwrite('.\filessim\spectra04fs.txt',IFTEsumHt);
    
%%  Read Files II
    cntlambdatau=dlmread('.\filessim\spectra03fs.txt'); %read raw data from detector
    cnttaulambda=(cntlambdatau);
%     cnttaulambda=transpose(cntlambdatau);
%     cnttaulambda=fliplr(cnttaulambda);
%     cnttau0lambda=cnttaulambda(601,:);
   
    figure
    plot(f,IFTEsumGt,'k',f,(0.8)*cnttaulambda,'g');
    legend('Simulated Spectral Interferogram (Gaussian)','Rhodamine Experimental Spectral Interferogram');
    title('Compare Experimental data and Simulation');
    xlim([f0-0.03 f0+0.03])
    
    figure
    plot(f,IFTEsumHt,'k',f,(0.8)*cnttaulambda,'c');
    legend('Simulated Spectral Interferogram (sech^2)','Rhodamine Experimental Spectral Interferogram');
    title('Compare Experimental data and Simulation');
    xlim([f0-0.03 f0+0.03])



    
    %% Signal and LO in frequency domain 
    
    FTGSIGt=fft(fftshift(GSIGt));
    FTGLOt=fft(fftshift(GLOt));
    FTGLOtcon=conj(FTGLOt);
    IFTGLOt=sqrt(FTGLOt.*FTGLOtcon);
%     FTGSIGt2=(FTGSIGt).*(FTGLOtcon);
%     FTGSIGt2=(FTGSIGt);
    FTGSIGt2=(FTGSIGtTj);
    
    FTHSIGt=fft(fftshift(HSIGt));
    FTHLOt=fft(fftshift(HLOt));
    FTHLOtcon=conj(FTGLOt);
    IFTHLOt=sqrt(FTGLOt.*FTHLOtcon);
%     FTHSIGt2=(FTHSIGt).*(FTHLOtcon);
%     FTHSIGt2=(FTHSIGt);
    FTHSIGt2=(FTHSIGtTj);
    
    figure
    subplot(2,1,1)
    plot(real(GSIGt));
    title('Gaussian')
    subplot(2,1,2)
    plot(real(FTGSIGt));
    %xlim([700 1400]);
    
    figure
    subplot(2,1,1)
    plot(real(HSIGt));
    title('sech^2');
    subplot(2,1,2)
    plot(real(FTHSIGt));
    %xlim([700 1400]);
    
%% 2D code first step (ifft and shift)

    shiftifftSI=ifftshift(ifft(IFTEsumGt));

    shiftifftSI2=ifftshift(ifft(IFTEsumHt));


%% Window in time domain to get rid of mirror and zero frequency components

    winshiftifftSI=shiftifftSI;
    winshiftifftSI2=shiftifftSI2;

    lowcut=5300;
    highcut=6000;

    for kwin=1:lowcut
    winshiftifftSI(kwin)=0;
    winshiftifftSI2(kwin)=0;
    end

    for kwin2=highcut:N
    winshiftifftSI(kwin2)=0;
    winshiftifftSI2(kwin2)=0;    
    end

    figure;
    subplot(2,1,1)
    plot(t,real(shiftifftSI));
    title('real(shiftifftcnttaulambda) Gaussian');
    subplot(2,1,2)
    plot(real(winshiftifftSI));
    title('real(shiftifftcnttaulambda) Gaussian after window');
    xlabel('1/wavelength (nm-1)');
    ylabel('counts');

    figure;
    subplot(2,1,1)
    plot(t,real(shiftifftSI2));
    title('real(shiftifftcnttaulambda) sech^2');
    subplot(2,1,2)
    plot(real(winshiftifftSI2));
    title('real(shiftifftcnttaulambda) sech^2 after window');
    xlabel('1/wavelength (nm-1)');
    ylabel('counts');

%%  fft and fftshift back to frequency domain

    fftshiftwinshiftifftcnttaulambda=fft(fftshift(winshiftifftSI));
    
    fftshiftwinshiftifftcnttaulambdaSECH=fft(fftshift(winshiftifftSI2));


%% Phasing

    lo(1:N)=0;

    for kff=1:npix 
    lo(N-npix+kff)=exp(-1i*(((2*pi)*c*(10^-6))/(lambda(npix)))*(600)); %2*pi
    end
    
    lofftshiftwinshiftifftcnttaulambda=(fftshiftwinshiftifftcnttaulambda).*(lo);
    lofftshiftwinshiftifftcnttaulambdaSECH=(fftshiftwinshiftifftcnttaulambdaSECH).*(lo);

    phlofftshiftwinshiftifftcnttaulambda=lofftshiftwinshiftifftcnttaulambda.*exp(1i*phi);
    phlofftshiftwinshiftifftcnttaulambda=phlofftshiftwinshiftifftcnttaulambda./IFTGLOt;
 
    phlofftshiftwinshiftifftcnttaulambdaSECH=lofftshiftwinshiftifftcnttaulambdaSECH.*exp(1i*phi);
    phlofftshiftwinshiftifftcnttaulambdaSECH=phlofftshiftwinshiftifftcnttaulambdaSECH./IFTHLOt;
    
    figure;
    plot(pixN,real(IFTEsumGt),'k',pixN,real(fftshiftwinshiftifftcnttaulambda),'r',pixN,real(phlofftshiftwinshiftifftcnttaulambda),'g');
    legend('raw data','real after fft(fftshift) back to frequency domain','real of pure electric field of signal after phasing');
    title('compare Gaussian');
    xlim([N-npix N]);
    ylim([-4000 12000]);

    figure;
    plot(pixN,real(IFTEsumHt),'k',pixN,real(fftshiftwinshiftifftcnttaulambdaSECH),'m',pixN,real(phlofftshiftwinshiftifftcnttaulambdaSECH),'c');
    legend('raw data','real after fft(fftshift) back to frequency domain','real of pure electric field of signal after phasing');
    title('compare sech^2');
    xlim([N-npix N]);
    ylim([-4000 12000]);


%% comparing the original signal and the retrieved signal usng 2D code

    figure
    plot(pixN,(1)*real(FTGSIGt2),'k',pixN,real(phlofftshiftwinshiftifftcnttaulambda),'g');
    legend('original signal','signal from 2D code');
    title('Gaussian');
    xlim([N-npix N]);
    ylim([-60 60]);

    figure
    plot(pixN,(1)*real(FTHSIGt2),'k',pixN,real(phlofftshiftwinshiftifftcnttaulambdaSECH),'c');
    legend('original signal','signal from 2D code');
    title('sech^2');
    xlim([N-npix N]);
    ylim([-60 60]);
    
    
