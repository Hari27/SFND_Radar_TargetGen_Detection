Read me: 
Step 1:
As a first step first define the specification of the radar with range max and range resolution
and define the target's initial position and target's speed.

%% Radar Specifications 

range_max = 200;
range_res = 1;
vmax = 70;
c = 3e8;
% define the target's initial position and velocity. Note : Velocity remains contant
initial_Position = 100; % Initial distance of the target should be within 200m.
v = 30; % speed of the target should be within -70 to 70.

Step 2: Second step is to generate the FMCW waveform chirp, we are calculating the BW (B),
chirp time (Tchirp) and Slope (slope) of the FMCW. Define the carrier frequency, the number
of chirps (Nd) and the number of samples on each chirp.

% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

B = c/2/range_res; % Bandwidth formula
Tchirp = 5.5*(2*range_max)/c; % formula for chirp taken from lecture
slp = B/Tchirp; % slope formula also taken from lecture

%Operating carrier frequency of Radar 
fc= 77e9; %carrier freq
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


Step 3: Prepare an initial vector where you can store the beat signal index wise.

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


Step 4: Loop for simulation of the Tx, Rx and Mix signals just multiply element wise 
between Tx and Rx.

for i=1:length(t)          

    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = initial_Position + (v*t(i));
    td(i) = 2*r_t(i)/c;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i)   = cos(2*pi*(fc*t(i) + 0.5*slp*t(i)^2));
    Rx(i)  = cos(2*pi*(fc*(t(i) - td(i)) + 0.5*slp*(t(i) - td(i))^2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i); % this is the beat signal.
    
end

Step 5: Range FFT for calculating the distance of the target.

reshaped_Mix = reshape(Mix,Nr*Nd,1);

% 1 D FFT
% *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
sig_range_fft = fft(reshaped_Mix(1:Nr)); 
sig_range_fft = sig_range_fft/max(sig_range_fft);


% *%TODO* :
% Take the absolute value of FFT output
sig_range_fft = abs(sig_range_fft);

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_range_fft_plt = sig_range_fft(1:Nr/2+1);
%fs = 1/(t(1000) - t(999));

%plotting the range
figure ('Name','Range from First FFT');
plot(sig_range_fft_plt);
title('Range from fft 1');
ylabel('Amplitude normalized');
xlabel('Range');


Step 6: 2D FFT for creating the Range Doppler Map

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

Step 7: CFAR implementation 
Tr is the number of training cells. Td is the number of training bands. Gr is the number of guard cells. Gd is the number of guard bands. 
Course provided parameters were used.

Tr = 12;
Td = 8;
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;


% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.2;


% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1); %zeros(Nr/2-2*(Td + Gd), Nd-2*(Tr + Gr));
gridSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);  % Follow the lecture @ 2D CFAR
trainingCells = gridSize - (2*Gr+1)*(2*Gd+1); % Follow the lecture @ 2D CFAR

Step 8: 2D CFAR is now ready, we can proceed to make the RDM (Range Doppler Mapping) 

RDM = RDM/max(max(RDM));
tot_den = 2*(Tr+Gr+1)*2*(Td+Gd+1) - (Gr*Gd-1);
for range_idx=Tr+Gr+1:(Nr/2)-(Gr+Tr)
    for doppler_idx=Td+Gd+1:Nd-(Gd+Td)
        count = 0;
        % Slide the CUT across the complete matrix. 
         for p = range_idx-(Tr+Gr):range_idx+Tr+Gr
            for q = doppler_idx-(Td+Gd):doppler_idx+(Td+Gd)
                if (abs(range_idx-p)>Gr || abs(doppler_idx-q)>Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                    count = count + 1;
                end   
            end
         end
        
% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

          threshold = pow2db(noise_level/tot_den);
          threshold = threshold + offset;
          Cell_Under_Test = RDM(range_idx, doppler_idx);
          if(Cell_Under_Test < threshold)
              RDM(range_idx,doppler_idx) = 0;
          else
              RDM(range_idx,doppler_idx) = 1; 
          end
        noise_level = 0; 
    end 
end

Last step: We need to overcome the map size problem by making the edges value to 0.

% Edges: initial part of the range axis and the final part of the range axis to be
% assigned to zero.
RDM(1:Tr+Gr,:) = 0;
RDM((Nr-(Tr+Gr))+1:end,:) = 0;
% Edges: initial part of the doppler axis and the final part of the doppler axis to be
% assigned to zero.
RDM(:,1:Td+Gd+1) = 0;
RDM(RDM~= 0 & RDM~= 1) = 0;


% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
title('threshold');
figure,surf(doppler_axis,range_axis,RDM);
colorbar;