% Implementación del Algoritmo DUET para separación ciega de fuentes
% Lucas Bernalte Pérez 
% Proyecto Fin de Carrera. Ingeniería de Telecomunicación.
% Universidad de Sevilla.
% 25/10/2014
%
% Notación:
%
% s1, s2, ... -> señales fuente.
% stf1, stf2, ... -> transformada de fourier de las señales fuente.
% SNR -> Relación señal a ruido de las fuentes (dB).
% FS -> frecuencia de muestreo.
% windowsize -> tamaño de la ventana para análisis STFT.
% window -> tipo de ventana usada (en nuestro caso, hamming).
% nfft -> resolución en frecuencia.
% overlap -> número de muestras que se solaparán.
% f -> vector de frecuencias
% fmat -> matriz que contiene el vector de frecuencias repetido en
% columnas, normalizada.
% a1, a2, ... -> atenuaciones que tendrán las fuentes en la mezcla.
% d1, d2, ... -> retrasos que tendrán las fuentes en la mezcla.
% tf1, tf2 -> observaciones (señales recogidas por los micrófonos) creadas
% mediante la mezcla instantánea.
% att, attsim, delta -> atenuación, atenuación simétrica y retrasos
% estimados.


%% PREPARACIÓN DE FUENTES

% SNR = 37.5; 
[s1,FS,NBITS]=wavread('dev1_male3_src_1.wav');
% s1 = awgn(s1,SNR);
[s2,FS,NBITS]=wavread('dev1_male3_src_2.wav');
% s2 = awgn(s2,SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREAR LA MEZCLA INSTANTÁNEA

windowsize = 1024;
window = hamming(windowsize);
nfft=1024;
overlap=0.25*windowsize; %Tiene este valor para que la antitransformada sea directa.

stf1 = stft(s1, window, overlap, nfft, FS);
stf2 = stft(s2, window, overlap, nfft, FS);

f = (0:size(stf1,1)-1)*FS/nfft;
fmat=repmat(f',1,size(stf1,2));
fmat=fmat/8000*pi;

%parámetros en la mezcla. Según DUET, ponemos a1=1 y d1=1.
a1=1;
d1=0.02;
a2=1;
d2=-0.02; %si tenemos más fuentes podemos seguir, a3, d3, a4, d4, etc.


n1 = 0*randn(size(s1,1),1);
n2 = stft(n1, window, overlap, nfft, FS);
SNR=20*log10(mean((s1+n1).^2)/mean(n1.^2))



tf1=stf1+stf2+n2; %Si tenemos más fuentes: +stf3+..., etc
tf2=a1*exp(-1i.*fmat*d1).*stf1+a2*exp(-1i.*fmat*d2).*stf2+n2; %Aquí haríamos lo mismo: +a3*exp(-1i.*fmat*d3).*stf3+...;



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

figure(1)
plot(alpha_ind,delta_ind,'.')
xi=linspace(min(alpha_ind),max(alpha_ind),abins);
yi=linspace(min(delta_ind),max(delta_ind),dbins);
xr=interp1(xi,1:numel(xi),alpha_ind,'nearest')';
yr=interp1(yi,1:numel(yi),delta_ind,'nearest')';
Z=accumarray([alpha_ind delta_ind],tfweight,[abins dbins]);

xaxis=linspace(-maxa,maxa,35);
yaxis=linspace(-maxd,maxd,50);
figure(2)
surf(yaxis,xaxis,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Introducir número de fuentes y su posición en el histograma manualmente

numfuentes = 2;
delta_p = [0 2];
alpha_p = [0 -0.45];
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

for i = 1:numfuentes
    mask=squeeze(iminimos==i);
    
    fuenteest=(tf1+a_p(i)*exp(sqrt(-1)*fmat*delta_p(i)).*tf2)./(1+a_p(i)^2);
    fuenteest=fuenteest.*squeeze(mask);
    
    y = zeros(1, length(s1));


    for j = 0:overlap:(overlap*(size(fuenteest, 2)-1))

        aux = fuenteest(:, 1+j/overlap);
        
        
        aux = [aux; conj(aux(end-1:-1:2))];
        aux = real(ifft(aux));
        
        % overlap, window and add
        y((j+1):(j+nfft)) = y((j+1):(j+nfft)) + (aux.*window)';
    end
    if i==1
        salida1=istft(fuenteest,overlap,window,nfft,FS);
        [SDR1]=bss_eval_sources(salida1,s1') %poner ; al final si no se desea mostrar la calidad de la separación de esta fuente.
    end
    if i==2
        salida2=istft(fuenteest,overlap,window,nfft,FS);
        [SDR2]=bss_eval_sources(salida2,s2')
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Si se van a pasar las fuentes estimadas a un archivo wav:

wavwrite(salida1,FS,'fuenteestimada1.wav');
wavwrite(salida2,FS,'fuenteestimada2.wav');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


