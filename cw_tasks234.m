beta=1;
phi=0.663225; %38 degrees = 0.663225 radians
id=fopen('input1.txt','r');
formatSpec='%c';
A=fscanf(id,formatSpec);
ascivals=double(A);
[~ , cols]=size(ascivals);
fclose(id);

%construct ansi code - 8 bits per character
encoder=strings(1,cols);
st1='0';
for i=1:1:cols
    temp=ascivals(i);
    temp=dec2bin(temp);
    if(length(temp)~=8)
        for j=1:(8-length(temp))
            temp=append(st1,temp);
        end
    end
    encoder(1,i)=temp;
end
%create input bit stream
x=zeros(1,cols*8);
for i=1:cols
  temp=encoder(1,i);
  b=str2double(regexp(num2str(temp),'\d','match'));
  for j=1:1:length(b)
      x(1,(i-1)*8+j)=b(j);
  end
  clear b;
end
%map the bitstream to symbols
r=sqrt(2);
s1=r*exp(1i*phi);
s2=r*exp(1i*(phi+(pi/2)));
s3=r*exp(1i*(phi+pi));
s4=r*exp(1i*(phi+(3*pi/2)));
symbs=[];
for i=1:2:length(x)
    if(x(i)==0 && x(i+1)==0)
       symbs=[symbs s1];
    elseif (x(i)==0 && x(i+1)==1)
        symbs=[symbs s2];
    elseif(x(i)==1 && x(i+1)==1)
         symbs=[symbs s3];
    else
        symbs=[symbs s4];
    end
end
%Adding noise-adjustable based on SNR
Ps= dot(symbs,symbs)/(length(symbs));
SNRin=input('Enter SNR value in dBs');
SNRin=10^(SNRin/10);
Pn=Ps/SNRin; %we keep the constellation points fixed and change the noise power accordingly
noise=sqrt(Pn)* (randn(size(symbs)) + 1i*randn(size(symbs)))/sqrt(2);
r_x=symbs+noise;

%constellation diagram
sz=30;
scatter(real(r_x),imag(r_x),sz);
hold on;
scatter(real(symbs),imag(symbs),sz,'r*');
title ('Constellation Diagram with Noise');
subtitle ('Decision region:adjacent lines and locus bounded by two adjacent symbols');
hold on;
line1x=[0 real(s1)];
line1y=[0 imag(s1)];
plot(line1x,line1y);
hold on;
line2x=[0 real(s2)];
line2y=[0 imag(s2)];
plot(line2x,line2y);
hold on;
line3x=[0 real(s3)];
line3y=[0 imag(s3)];
plot(line3x,line3y);
hold on;
line4x=[0 real(s4)];
line4y=[0 imag(s4)];
plot(line4x,line4y,'color','g');
legend ('noise-corrupted symbol','reference symbols','00','01','11','10', 'FontSize',7);
xlabel('Real Part'); ylabel('Imaginary Part');

%demodulation- %Decision function- choose s(i) based on each element's
%distance angle-wise from closest symbols (map it to the closest s(i))
ref_angles=[angle(s1) angle(s2) angle(s3) angle(s4)];
ref_angles=wrapTo2Pi(ref_angles);
angles= angle(r_x);
angles=wrapTo2Pi(angles);
dm=[];

for ind=1:length(angles)
    if(ref_angles(1)<=angles(ind))&&(angles(ind)<=ref_angles(2))
        if(abs(angles(ind)-ref_angles(1))<= abs(ref_angles(2)-angles(ind)))
            dm=[dm 0,0];
        else
            dm=[dm 0,1];
        end
    elseif (ref_angles(2)<=angles(ind))&&(angles(ind)<=ref_angles(3))
        if(abs(angles(ind)-ref_angles(2))<= abs(ref_angles(3)-angles(ind)))
            dm=[dm 0,1];
        else
            dm=[dm 1,1];
        end
    elseif (ref_angles(3)<=angles(ind))&&(angles(ind)<=ref_angles(4))
        if(abs(angles(ind)-ref_angles(3))<= abs(ref_angles(4)-angles(ind)))
            dm=[dm 1,1];
        else
            dm=[dm 1,0];
        end
    elseif(angles(ind)>ref_angles(4))
        dm=[dm 1,0];
    else
        dm=[dm 0,0];
    end
end

%decoding 
decoder=strings(1,length(dm)/8);
temparr=[];
counter=1;
for k=1:length(dm)
    if(k==counter*8+1)
        temp_num=num2str(temparr);
        decoder(1,counter)=temp_num;
        clear temparr; 
        temparr=[]; %must redeclare the array once cleared
        counter=counter+1;
    end
       no=dm(k);
       temparr=[temparr no];
end
decoder = regexprep(decoder, '\s+', '') ; %we remove the spaces
ansi_out = bin2dec(decoder); %output ansi code
recovered_message=char(ansi_out);   %convert ansi code to original message
dlmwrite('recovered_message.txt',recovered_message,'delimiter','');

%number of bits in error
bits_in_error=0;
for i=1:length(x)
    if(x(1,i)~=dm(1,i))
        bits_in_error=bits_in_error+1;
    end
end
error_prob=bits_in_error/length(x); %practical ber 
%theoretical ber
pe_theoretical= qfunc(sqrt(SNRin));



