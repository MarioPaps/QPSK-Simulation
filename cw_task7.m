beta=1;
phi=0.663225; %38 degrees = 0.663225 radians
id=fopen('input1.txt','r');
formatSpec='%c';
A=fscanf(id,formatSpec);
ascivals=double(A);
[~ , cols]=size(ascivals);

id2=fopen('jammer.txt','r');
formatSpec='%c';
B=fscanf(id2,formatSpec);
ascivalsjam=double(B);

%construct ansi code -  8 bits per character
encoder=strings(1,cols);
jamencoder=strings(1,cols);
st1='0';
for i=1:1:cols
    temp=ascivals(i);
    temp=dec2bin(temp);
    jamtemp=ascivalsjam(i);
    jamtemp=dec2bin(jamtemp);
    if(length(temp)~=8)
        for j=1:(8-length(temp))
            temp=append(st1,temp);
        end
    end
    encoder(1,i)=temp;
    if(length(jamtemp)~=8)
        for k=1:(8-length(jamtemp))
            jamtemp=append(st1,jamtemp);
        end
    end
    jamencoder(1,i)=jamtemp;
end

%create input bit stream and jammer bit stream
x=zeros(1,cols*8);
xjam=zeros(1,cols*8);
for m=1:cols
  temp=encoder(1,m);
  b=str2double(regexp(num2str(temp),'\d','match'));
  jamtemp=jamencoder(1,m);
  c=str2double(regexp(num2str(jamtemp),'\d','match'));
  for j=1:1:length(b)
      x(1,(m-1)*8+j)=b(j);
      xjam(1,(m-1)*8+j)=c(j);
  end
  clear b; clear c; 
end
fclose('all');

%map the bitstream to symbols
r=sqrt(2);
s1=r*exp(1i*phi);
s2=r*exp(1i*(phi+(pi/2)));
s3=r*exp(1i*(phi+pi));
s4=r*exp(1i*(phi+(3*pi/2)));
symbs=[];
jamsymbs=[];
r2= sqrt(20);
ref_jam=[r2*exp(1i*phi)  r2*exp(1i*(phi+(pi/2))) r2*exp(1i*(phi+pi)) r2*exp(1i*(phi+(3*pi/2)))];
for m=1:2:length(x)
    %message signal symbols
    if(x(m)==0 && x(m+1)==0)
        symbs=[symbs s1];
    elseif (x(m)==0 && x(m+1)==1)
        symbs=[symbs s2];
    elseif(x(m)==1 && x(m+1)==1)
         symbs=[symbs s3];
    else
        symbs=[symbs s4];
    end
    %jammer symbols
     if(xjam(m)==0 && xjam(m+1)==0)
        jamsymbs=[jamsymbs ref_jam(1)];
    elseif (xjam(m)==0 && xjam(m+1)==1)
        jamsymbs=[jamsymbs ref_jam(2)];
    elseif(xjam(m)==1 && xjam(m+1)==1)
         jamsymbs=[jamsymbs ref_jam(3)];
    else
        jamsymbs=[jamsymbs ref_jam(4)];
    end
end

%message signal and jammer signal gold code generator
seq1=[1 1 1 1 1]; seq2=[1 1 1 1 1];
seq1_next=zeros(1,5); seq2_next=zeros(1,5);
Nc= 2^length(seq1)-1;
des1=[seq1(1,5)]; des2=[seq2(1,5)];
jam1=[1 1 1 1 1]; jam2=[1 1 1 1 1];
jam1_next=zeros(1,5); jam2_next=zeros(1,5);
jam1_out=[jam1(1,5)]; jam2_out=[jam2(1,5)];
for k=1:(Nc-1)
     %desired message: D^5+D^2+1
     seq1_next(1,1)= xor(seq1(k,3),seq1(k,5));
     seq1_next(1,2:5)= seq1(k, 1:4);
     seq1=[seq1; seq1_next];
     des1(1,k+1)= seq1(k+1,5);
     clear seq1_next;
     %desired message: % D^5 + D^3 + D^2 + D + 1 
     seq2_next(1,1)=xor(xor(seq2(k,2),seq2(k,3)), xor(seq2(k,4),seq2(k,5)));
     seq2_next(1,2:5)=seq2(k,1:4);
     seq2=[seq2; seq2_next];
     des2(1,k+1)=seq2(k+1,5);
     clear seq2_next;
     %jammer message: D^5+D^3+D^2+D+1
     jam1_next(1,1)= xor(jam1(k,2),jam1(k,5));
     jam1_next(1,2:5)= jam1(k, 1:4);
     jam1=[jam1;jam1_next];
     jam1_out(1,k+1)=jam1(k+1,5);
     clear jam1_next;
     %jammer message: D^5+D^4+D^2+D+1
     jam2_next(1,1)= xor(xor(jam2(k,3),jam2(k,1)), xor(jam2(k,4),jam2(k,5)));
     jam2_next(1,2:5)=jam2(k,1:4);
     jam2=[jam2; jam2_next];
     jam2_out(1,k+1)=jam2(k+1,5);
     clear jam2_next;
end

% message signal gold code -> delay that satisfies both conditions is delay=27
des2=circshift(des2,27);
mes_gold=bitxor(des1,des2); %equivalent to xor function->preferable to obtain a double array

%jammer signal gold code ->delay that satisfies both conditions is delay=28
jam2_out= circshift(jam2_out,28);
jam_gold=bitxor(jam1_out,jam2_out); %equivalent to xor function->preferable to obtain a double array


%spread spectrum ; we map binary 0s to +1s and  binary 1s to -1s
mes_gold= -2*mes_gold+1;
symbs_ext= repelem(symbs,Nc);
mes_gold_ext=repmat(mes_gold,[1 length(symbs_ext)/Nc]);
symbs_pn= symbs_ext.*mes_gold_ext;

jam_gold= -2*jam_gold+1;
jamsymbs_ext= repelem(jamsymbs,Nc);
jam_gold_ext=repmat(jam_gold,[1 length(jamsymbs_ext)/Nc]);
jam_pn= jamsymbs_ext.*jam_gold_ext;

%channel: message signal+ jammer interference + noise
SNRin=1000;
Ps=1;
Pn= Ps/SNRin; 
Pj=10;
r_x= symbs_pn+jam_pn + sqrt(Pn)*((randn(size(symbs_pn)))+1i*randn(size(symbs_pn)))/ sqrt(2); 


%constellation diagram
sz=30; %constellation diagram
scatter(real(r_x),imag(r_x),sz);
hold on;
scatter(real(symbs),imag(symbs),sz,'r*');
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

title('Constellation Diagram for SS QPSK');
subtitle ('Decision region:adjacent lines and locus bounded by two adjacent symbols');
xlabel('Real Part'); ylabel('Imaginary Part');
legend ('noise-corrupted symbol','reference symbols','00','01','11','10','FontSize',7);


%despreader code
mult=r_x.*mes_gold_ext; %multiply by pn sequence
despread_rx=zeros(1,length(mult)/Nc);
for m=1:(length(mult)/Nc)
    sum=0;
    for p=1:Nc
        sum=sum+mult((m-1)*Nc+p);
    end
    average=sum/Nc;
    despread_rx(m)=average;
end

%demodulation- %Decision function- choose s(i) based on each element's
%distance angle-wise from closest symbols (map it to the closest s(i))
ref_angles=[angle(s1) angle(s2) angle(s3) angle(s4)];
ref_angles=wrapTo2Pi(ref_angles);
angles= angle(despread_rx);
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
recovered_message=char(ansi_out);   %converted ansi code to original message
dlmwrite('recovered_message.txt',recovered_message,'delimiter','');

%number of bits in error
bits_in_error=0;
for m=1:length(x)
    if(x(1,m)~=dm(1,m))
        bits_in_error=bits_in_error+1;
    end
end
error_prob=bits_in_error/length(x); %practical ber prob








    