close all;clear all;clc;

x=randn(500,1);
M=length(x);
mean = mean(x);
variance=var(x);
a=0.1;
[y1,fs,n] = wavread('s2.wav'); % i/p audio file
v = y1(1:length(y1)-mod(length(y1),20),1);
D = dct(v);
D_abs = abs(D);
seg = buffer(D_abs,length(v)/20);
en = sum(seg(:,1:end).^2);

% y = seg(:,find(en==max(en)));
ind=find(en==max(en));
y = seg(:,ind);
[Y_sort,in]=sort(y);
N = length(Y_sort);
z_in = in(N:-1:N-M+1);
% z = Y_sort(N:-1:N-M+1);
% z1 = y(z_in);
y_new=y;
z = y(z_in);
y_new(z_in)=y(z_in).*(1+a.*x);
D_new=[];
for i=1:20
   if (i==ind) 
      D_new=[D_new;y_new] ;
   else
       D_new=[D_new;seg(:,i)];
   end
end
%convert absolute back to real values
D_ind = find(D<0);
D_real = D_new;
for i=1:length(D_ind)
   D_real(D_ind(i)) = -D_new(D_ind(i));
end

% IDCT
v_re = idct(D_real);
figure(1),
subplot 311,plot(v(:,1));title('i/p audio file');axis([0 length(v) -(max(abs(v))+0.2) (max(abs(v))+0.2)]);
subplot 312,plot(v_re);title('watermarked audio file');axis([0 length(v_re) -(max(abs(v_re))+0.2) (max(abs(v_re))+0.2)]);
subplot 313,plot(v_re-v);title('difference signal');axis([0 length(v) -0.5 0.5])

wavwrite(v_re,fs,n,'s2_reconstructed');

% Gaussian attack
v_at = awgn(v_re,10,'measured','linear');
% v_at = v_re;
D_at = abs(dct(v_at));
seg_at = buffer(D_at,length(v_at)/20);
en_at = sum(seg_at(:,1:end).^2);
% subplot 313,plot(v_at);
ind_at=find(en_at==max(en_at));
y_at = seg_at(:,ind_at);
[Y_at,in_at]=sort(y_at);
N = length(Y_at);
% z_in_at = in_at(N:-1:N-M+1);
z_at = y_at(z_in);
X = ((z_at./z)-1)./a;

figure(2),
subplot 311,plot(v_at);title('AWGN attack on watermarked signal');axis([0 length(v_at) -(max(abs(v_at))+0.2) (max(abs(v_at))+0.2)]);
subplot 312,plot(x);title('original watermarking sequence');axis([0 length(x) -(max(abs(x))+0.2) (max(abs(x))+0.2)]);
subplot 313,plot(X);title('reconstructed watermark sequence');axis([0 length(x) -(max(abs(X))+0.2) (max(abs(X))+0.2)]);


%SNR

o=sum(v.^2);
O=sum((v-v_re).^2);
SNR = 10*log10(o/O)

figure(3)
plot(y);title('highest energy segment');
