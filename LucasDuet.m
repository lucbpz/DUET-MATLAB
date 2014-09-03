clear
close all
%% Importar mezcla

[x,FS,NBITS]=wavread('roomsim_D004.wav');

x1=x(:,1);
x2=x(:,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables iniciales

windowsize = 1024;
window = hamming(windowsize);
window2=window(1:1023);

overlap=0.25*windowsize; %overlap del 75% 
overlap2=0.75*windowsize; %y este es el overlap que usaba para antitransformar
nfft=1024;

fvec=[linspace(1/(nfft/2),pi,nfft/2-1) linspace(-pi,0,nfft/2)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transformar a STFT.

[tf1,F,T]=spectrogram(x1,window,overlap,fvec);
tf2 = spectrogram(x2,window,overlap,fvec);
[sizex,sizey]=size(tf1);

fmat=repmat(fvec,1,length(T)); %fmat=[fmat;-flipud(fmat)];


%% Cálculo de datos

R21=(tf2)./(tf1);
att=abs(R21);
attsim = att - 1./att; %atenuacionsimetrica
delta = -imag(log(R21))./fmat ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histograma ponderado

p = 1; q = 0; % powers used to weight histogram 	
tfweight = (abs(tf1).*abs(tf2)).^p.*abs(fmat).^q ; %weights 	
maxa = 0.7 ; maxd = 3.6 ; % hist boundaries 	
abins = 35; dbins = 50; % number o f h i s t bins 	
% 
% %only consider time?freq  points yielding estimates in bounds 	
amask=(abs(attsim)<maxa)&(abs(delta)<maxd) ; 	
alpha_vec = attsim(amask); 	
delta_vec = delta(amask); 	
tfweight = tfweight(amask); 	
% 
% % determine histogram indices 	
alpha_ind = round(1+(abins-1)*(alpha_vec+maxa)/(2*maxa)); 	
delta_ind = round(1+(dbins-1)*(delta_vec+maxd)/(2*maxd));
% 
% % FULL?SPARSE TRICK TO CREATE 2D WEIGHTED HISTOGRAM 	
 A = full(sparse(alpha_ind , delta_ind , tfweight , abins , dbins )); 	
% 
% % smooth the histogram ?local average 3?by?3 neighboring bins 	
% A = twoDsmooth(A,3); 	
% 
% %plot 2?D histogram 
 mesh(linspace(-maxd,maxd,dbins), linspace(-maxa,maxa,abins),A);	
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determinar picos de histograma manualmente
numfuentes = 3;
delta_p = [0 -1 1];
alpha_p = [0 0 0];
% convert alpha to a 	
a_p = (alpha_p+sqrt(alpha_p.^2+4))/2; 	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construccion de máscaras
for i = 1:numfuentes
    jft(i,:,:) = (abs(((a_p(i)*exp(-sqrt(-1)*fmat*delta_p(i))).*tf1)-tf2).^2)/(1+a_p(i)^2); 
end

[minimos,iminimos]=min(jft,[],1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Antitransformar

mask=[];
overlap2=0.75*windowsize;
row=zeros(1,length(x1)+2*windowsize);

for i = 1:numfuentes
    mask=squeeze(iminimos==i);
    
    fuenteest=(tf1+a_p(i)*exp(sqrt(-1)*fmat*delta_p(i)).*tf2)./(1+a_p(i)^2);
    fuenteest=fuenteest.*squeeze(mask);
    
%     prueba1=[fuenteest; conj(fuenteest(end-1:-1:2))];
%         fuenteinv=real(ifft(prueba1));

    fuenteinv=real(ifft(fuenteest));
    
    for j = 1:sizey % overlap , window , and add 	
        ventanaoverlap = ((j-1)*overlap2+1):((j-1)*overlap2+windowsize-1);
        row(ventanaoverlap)=row(ventanaoverlap)+((fuenteinv(:,j)).*window2)';
    end
end


%%%%%%%%%%%%%%%%probar
wavwrite(row,FS,'lucasprueba.wav');

