%  Example of use of the new Field II program running under Matlab
%
%  This example shows how a linear array B-mode system scans an image
%  when doing color flow mapping
%
%  This script assumes that the field_init procedure has been called
%  Here the field simulation is performed and the data is stored
%  in rf-files; one for each rf-line done. The data must then
%  subsequently be processed to yield the image. The data for the
%  scatteres are read from the file pht_data.mat, so that the procedure
%  can be started again or run for a number of workstations.
%
%  Example by Joergen Arendt Jensen and Peter Munk, March 14, 1997.
%  Version 2.1 by Joergen Arendt Jensen, August 14, 1998 for Matlab 5.

%  Generate the transducer apertures for send and receive

f0=5e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=0.05/1000;          %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]
N_elements=196;          %  Number of physical elements
rec_N_active=64;         %  Number of active elements in receive
xmit_N_active=64;        %  Number of active elements in transmit

% Use triangles

set_field('use_triangles',0);

%  Set the sampling frequency

set_sampling(fs);

%  Generate aperture for emission

xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);

%  Set the impulse response and excitation of the xmit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (xmit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (xmit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%  Do for the number of CFM lines

Ncfm=10;
for k=1:Ncfm

%   Load the computer phantom

  cmd=['load sim_flow/scat_',num2str(k),'.mat']
  eval(cmd)

  %   Do linear array imaging

  no_lines=20;                    %  Number of lines in image
  image_width=40/1000;           %  Size of image sector
  d_x=image_width/(no_lines-1);   %  Increment for image

  %  Set the different focal zones for reception

  rec_zone_start=30/1000;
  rec_zone_stop=100/1000;
  rec_zone_size=10/1000;

  focal_zones_center=[rec_zone_start:rec_zone_size:rec_zone_stop]';
  focal_zones=focal_zones_center-0.5*rec_zone_size;
  Nf=max(size(focal_zones));
  focus_times=focal_zones/1540;

  %  Set a Hanning apodization on the receive aperture
  %  Dynamic opening aperture is used.

  Fnumber=2.0;
  rec_N_active_dyn=round(focal_zones_center./(Fnumber*(width+kerf)));

  for ii=1:Nf
    if rec_N_active_dyn(ii)>rec_N_active 
      rec_N_active_dyn(ii)=rec_N_active; 
      end
    rec_N_pre_dyn(ii) = ceil(rec_N_active/2  - rec_N_active_dyn(ii)/2);
    rec_N_post_dyn(ii) = rec_N_active - rec_N_pre_dyn(ii) - rec_N_active_dyn(ii);
    rec_apo=(ones(1,rec_N_active_dyn(ii)));
    rec_apo_matrix_sub(ii,:)=[zeros(1,rec_N_pre_dyn(ii)) rec_apo zeros(1,rec_N_post_dyn(ii))];
    end

  %  Transmit focus

  z_focus=40/1000;   

  %   Set a Hanning apodization on the xmit aperture

  xmit_apo=hanning(xmit_N_active)';
  %xmit_apo=ones(1,xmit_N_active);

  % Do imaging line by line

  i_start=1;
  x= -image_width/2 +(i_start-1)*d_x;

  for i=i_start:no_lines
  i
    %   Set the focus for this direction

    xdc_center_focus (emit_aperture, [x 0 0]);
    xdc_focus (xmit_aperture, 0, [x 0 z_focus]);
    xdc_center_focus (receive_aperture, [x 0 0]);
    xdc_focus (receive_aperture, focus_times, [x*ones(Nf,1), zeros(Nf,1), focal_zones]);

    %  Calculate the apodization 
   
    xmit_N_pre  = round(x/(width+kerf) + N_elements/2 - xmit_N_active/2);
    xmit_N_post = N_elements - xmit_N_pre - xmit_N_active;
    xmit_apo_vector=[zeros(1,xmit_N_pre) xmit_apo zeros(1,xmit_N_post)];

    rec_N_pre(i) = round(x/(width+kerf) + N_elements/2 - rec_N_active/2);
    rec_N_post(i) = N_elements - rec_N_pre(i) - rec_N_active;
 
    rec_apo_matrix=[zeros(size(focus_times,1),rec_N_pre(i)) rec_apo_matrix_sub zeros(size(focus_times,1),rec_N_post(i))];

    xdc_apodization (xmit_aperture, 0, xmit_apo_vector);
    xdc_apodization (receive_aperture, focus_times , rec_apo_matrix);
   
    %   Calculate the received response

    [rf_data, tstart]=calc_scat(xmit_aperture, receive_aperture, positions, amp);

    %  Store the result

    cmd=['save sim_flow/rft',num2str(k),'l',num2str(i),'.mat rf_data tstart']
    eval(cmd)

    %  Steer in another direction

    x = x + d_x;
    
    end  %  Loop for lines

  end  %  CFM loop
  
%   Free space for apertures

xdc_free (xmit_aperture)
xdc_free (receive_aperture)

