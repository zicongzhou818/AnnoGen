clc
clear;
close all;
format long
%% MonaLisa Test
N=128;h=2;
m=N;n=N;
x=1:n;
y=1:m;
[X,Y]=ndgrid(x,y);
White_bg = zeros(m,n)+255;
%% 

T_original=double(rgb2gray(imread('holefilled1DrA.jpg')));
R_original=double(rgb2gray(imread('holefilled1DrB.jpg')));
S_original=double(rgb2gray(imread('holefilled1DrC.jpg')));

T=imresize(T_original,[m,n]);
R=imresize(R_original,[m,n]);
S=imresize(S_original,[m,n]);

%% 1st avg
tic
[postr_x,postr_y,Tn_tr,U1tr,U2tr]=Img_Reg2D_NewDev(T,R,N);
toc
tic
[posts_x,posts_y,Tn_ts,U1ts,U2ts]=Img_Reg2D_NewDev(T,S,N);
toc

figure(1)
subplot(2,4,1), imshow(R,[]), title('R');
subplot(2,4,2), imshow(Tn_tr,[]), title('T(\phi_{tr})');
subplot(2,4,3), imshow(abs(Tn_tr-R),[]), title('|T(\phi_{tr})-R|');
subplot(2,4,4), imshow(White_bg,[]), title('\phi_{tr}'), hold on 
for i = 1:h:N
    plot(postr_y(i,1:+h:end),postr_x(i,1:+h:end),'r-'), hold on
    plot(postr_y(1:+h:end,i),postr_x(1:+h:end,i),'r-'), hold on
end
subplot(2,4,5), imshow(S,[]), title('S');
subplot(2,4,6), imshow(Tn_ts,[]), title('T(\phi_{ts})');
subplot(2,4,7), imshow(abs(Tn_ts-S),[]), title('|T(\phi_{ts})-S|');
subplot(2,4,8), imshow(White_bg,[]), title('\phi_{ts}'), hold on 
for i = 1:h:N
    plot(posts_y(i,1:+h:end),posts_x(i,1:+h:end),'r-'), hold on
    plot(posts_y(1:+h:end,i),posts_x(1:+h:end,i),'r-'), hold on
end

[ftr, gtr] = compute_JD_and_Curl(postr_x,postr_y,1);
[fts, gts] = compute_JD_and_Curl(posts_x,posts_y,1);
jdnt = nthroot(ftr.*fts,3);
cvnt = (gtr+gts)/3;
tic
[avgt_x,avgt_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnt, cvnt,N,X,Y);
toc
avgt = interpn(X, Y, T, avgt_x,avgt_y,'makima'); 


% 2nd avg
tic
[posavgtt_x,posavgtt_y,Tn_avgtt,U1avgtt,U2avgtt]=Img_Reg2D_NewDev(avgt,T,N);
toc
tic
[posavgtr_x,posavgtr_y,Tn_avgtr,U1avgtr,U2avgtr]=Img_Reg2D_NewDev(avgt,R,N);
toc
tic
[posavgts_x,posavgts_y,Tn_avgts,U1avgts,U2avgts]=Img_Reg2D_NewDev(avgt,S,N);
toc

[favgtt, gavgtt] = compute_JD_and_Curl(posavgtt_x,posavgtt_y,1);
[favgtr, gavgtr] = compute_JD_and_Curl(posavgtr_x,posavgtr_y,1);
[favgts, gavgts] = compute_JD_and_Curl(posavgts_x,posavgts_y,1);
jdn1 = nthroot(favgtt.*favgtr.*favgts,3);
cvn1 = (gavgtt+gavgtr+gavgts)/3;

tic
[avgt2_x,avgt2_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdn1, cvn1,N,X,Y);
toc
avgt2 = interpn(X, Y, avgt, avgt2_x,avgt2_y,'makima'); 

figure(2)
subplot(3,4,1), imshow(T,[]);
subplot(3,4,2), imshow(Tn_avgtt,[]);
subplot(3,4,3), imshow(abs(Tn_avgtt-T),[]);
subplot(3,4,4), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgtt_y(i,1:+h:end),posavgtt_x(i,1:+h:end),'r-'), hold on
    plot(posavgtt_y(1:+h:end,i),posavgtt_x(1:+h:end,i),'r-'), hold on
end
subplot(3,4,5), imshow(R,[]);
subplot(3,4,6), imshow(Tn_avgtr,[]);
subplot(3,4,7), imshow(abs(Tn_avgtr-R),[]);
subplot(3,4,8), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgtr_y(i,1:+h:end),posavgtr_x(i,1:+h:end),'r-'), hold on
    plot(posavgtr_y(1:+h:end,i),posavgtr_x(1:+h:end,i),'r-'), hold on
end
subplot(3,4,9), imshow(S,[]);
subplot(3,4,10), imshow(Tn_avgts,[]);
subplot(3,4,11), imshow(abs(Tn_avgts-S),[]);
subplot(3,4,12), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgts_y(i,1:+h:end),posavgts_x(i,1:+h:end),'r-'), hold on
    plot(posavgts_y(1:+h:end,i),posavgts_x(1:+h:end,i),'r-'), hold on
end


tic
[posavgt2t_x,posavgt2t_y,Tn_avgt2t,U1avgt2t,U2avgt2t]=Img_Reg2D_NewDev(avgt2,T,N);
toc
tic
[posavgt2r_x,posavgt2r_y,Tn_avgt2r,U1avgt2r,U2avgt2r]=Img_Reg2D_NewDev(avgt2,R,N);
toc
tic
[posavgt2s_x,posavgt2s_y,Tn_avgt2s,U1avgt2s,U2avgt2s]=Img_Reg2D_NewDev(avgt2,S,N);
toc


figure(3)
subplot(3,4,1), imshow(T,[]);
subplot(3,4,2), imshow(Tn_avgt2t,[]);
subplot(3,4,3), imshow(abs(Tn_avgt2t-T),[]);
subplot(3,4,4), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgt2t_y(i,1:+h:end),posavgt2t_x(i,1:+h:end),'r-'), hold on
    plot(posavgt2t_y(1:+h:end,i),posavgt2t_x(1:+h:end,i),'r-'), hold on
end
subplot(3,4,5), imshow(R,[]);
subplot(3,4,6), imshow(Tn_avgt2r,[]);
subplot(3,4,7), imshow(abs(Tn_avgt2r-R),[]);
subplot(3,4,8), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgt2r_y(i,1:+h:end),posavgt2r_x(i,1:+h:end),'r-'), hold on
    plot(posavgt2r_y(1:+h:end,i),posavgt2r_x(1:+h:end,i),'r-'), hold on
end
subplot(3,4,9), imshow(S,[]);
subplot(3,4,10), imshow(Tn_avgt2s,[]);
subplot(3,4,11), imshow(abs(Tn_avgt2s-S),[]);
subplot(3,4,12), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgt2s_y(i,1:+h:end),posavgt2s_x(i,1:+h:end),'r-'), hold on
    plot(posavgt2s_y(1:+h:end,i),posavgt2s_x(1:+h:end,i),'r-'), hold on
end

%form distribution
[favgt2t, gavgt2t] = compute_JD_and_Curl(posavgtt_x,posavgtt_y,1);
[favgt2r, gavgt2r] = compute_JD_and_Curl(posavgtr_x,posavgtr_y,1);
[favgt2s, gavgt2s] = compute_JD_and_Curl(posavgts_x,posavgts_y,1);
jdnt2 = nthroot(favgt2t.*favgt2r.*favgt2s,3);
cvnt2 = (gavgt2t+gavgt2r+gavgt2s)/3;
tic
[avgt3_x,avgt3_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnt2, cvnt2,N,X,Y);
toc
avgt3 = interpn(X, Y, avgt2, avgt3_x, avgt3_y,'makima'); 

tic
[posavgt3t_x,posavgt3t_y,Tn_avgt3t,U1avgt3t,U2avgt3t]=Img_Reg2D_NewDev(avgt3,T,N);
toc
tic
[posavgt3r_x,posavgt3r_y,Tn_avgt3r,U1avgt3r,U2avgt3r]=Img_Reg2D_NewDev(avgt3,R,N);
toc
tic
[posavgt3s_x,posavgt3s_y,Tn_avgt3s,U1avgt3s,U2avgt3s]=Img_Reg2D_NewDev(avgt3,S,N);
toc
[favgt3t, gavgt3t] = compute_JD_and_Curl(posavgt3t_x,posavgt3t_y,1);
[favgt3r, gavgt3r] = compute_JD_and_Curl(posavgt3r_x,posavgt3r_y,1);
[favgt3s, gavgt3s] = compute_JD_and_Curl(posavgt3s_x,posavgt3s_y,1);
jdnt3 = nthroot(favgt3t.*favgt3r.*favgt3s,3);
cvnt3 = (gavgt3t+gavgt3r+gavgt3s)/3;
aa=log(favgt3t./jdnt3);
bb=log(favgt3t./jdnt3);
cc=log(favgt3s./jdnt3);
stdev=nthroot((aa.^2+bb.^2+cc.^2)./3,2);
sdjdnt = exp(stdev);
sdcvnt = sqrt((cvnt3-gavgt3t).^2+(cvnt3-gavgt3r).^2+(cvnt3-gavgt3s).^2)/2;


tic
[avgt4_x,avgt4_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnt3, cvnt3,N,X,Y);
toc
avgt4 = interpn(X, Y, avgt3, avgt4_x, avgt4_y,'makima'); 

figure(4)
subplot(1,4,1), imshow(avgt,[]);
subplot(1,4,2), imshow(avgt2,[]);
subplot(1,4,3), imshow(avgt3,[]);
subplot(1,4,4), imshow(avgt4,[]);


tic
[avgnt1sd_x,avgnt1sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(sdjdnt, sdcvnt,N,X,Y);
toc
avgnt1sd = interpn(X, Y, avgt3, avgnt1sd_x,avgnt1sd_y,'makima'); 
tic
[avgnt2sd_x,avgnt2sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdnt.^2), 2*sdcvnt,N,X,Y);
toc
avgnt2sd = interpn(X, Y, avgt3, avgnt2sd_x,avgnt2sd_y,'makima'); 
tic
[avgntb1sd_x,avgntb1sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdnt.^(-1)), -sdcvnt,N,X,Y);
toc
avgntb1sd = interpn(X, Y, avgt3, avgntb1sd_x,avgntb1sd_y,'makima'); 
tic
[avgntb2sd_x,avgntb2sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdnt.^(-2)), -2*sdcvnt,N,X,Y);
toc
avgntb2sd = interpn(X, Y, avgt3, avgntb2sd_x,avgntb2sd_y,'makima'); 

figure(5)
h=3;
subplot(2,5,5), imshow(White_bg,[]), title('2\sigma'), hold on 
for i = 1:h:N
    plot(avgnt2sd_y(i,1:+h:end),avgnt2sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnt2sd_y(1:+h:end,i),avgnt2sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,4), imshow(White_bg,[]), title('\sigma'), hold on 
for i = 1:h:N
    plot(avgnt1sd_y(i,1:+h:end),avgnt1sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnt1sd_y(1:+h:end,i),avgnt1sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,3), imshow(White_bg,[]), title('Id'), hold on 
for i = 1:h:N
    plot(Y(i,1:+h:end),X(i,1:+h:end),'r-'), hold on
    plot(Y(1:+h:end,i),X(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,2), imshow(White_bg,[]), title('-\sigma'), hold on
for i = 1:h:N
    plot(avgntb1sd_y(i,1:+h:end),avgntb1sd_x(i,1:+h:end),'r-'), hold on
    plot(avgntb1sd_y(1:+h:end,i),avgntb1sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,1), imshow(White_bg,[]), title('-2\sigma'), hold on
for i = 1:h:N
    plot(avgntb2sd_y(i,1:+h:end),avgntb2sd_x(i,1:+h:end),'r-'), hold on
    plot(avgntb2sd_y(1:+h:end,i),avgntb2sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

subplot(2,5,10), imshow(avgnt2sd,[]), title('I_{2\sigma}');
subplot(2,5,9), imshow(avgnt1sd,[]), title('I_{\sigma}');
subplot(2,5,8), imshow(avgt3,[]), title('I_{avg1}');
subplot(2,5,7), imshow(avgntb1sd,[]), title('I_{-\sigma}');
subplot(2,5,6), imshow(avgntb2sd,[]), title('I_{-2\sigma}');

%% 2nd avg
tic 
[posrt_x,posrt_y,Tn_rt,U1rt,U2rt]=Img_Reg2D_NewDev(R,T,N);
toc
tic
[posrs_x,posrs_y,Tn_rs,U1rs,U2rs]=Img_Reg2D_NewDev(R,S,N);
toc

figure(6)
subplot(2,4,1), imshow(T,[]), title('T');
subplot(2,4,2), imshow(Tn_rt,[]), title('T(\phi_{rt})');
subplot(2,4,3), imshow(abs(Tn_rt-T),[]), title('|T(\phi_{rt})-T|');
subplot(2,4,4), imshow(White_bg,[]), title('\phi_{tr}'), hold on 
for i = 1:h:N
    plot(posrt_y(i,1:+h:end),posrt_x(i,1:+h:end),'r-'), hold on
    plot(posrt_y(1:+h:end,i),posrt_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,4,5), imshow(S,[]), title('S');
subplot(2,4,6), imshow(Tn_rs,[]), title('T(\phi_{rs})');
subplot(2,4,7), imshow(abs(Tn_rs-S),[]), title('|T(\phi_{rs})-S|');
subplot(2,4,8), imshow(White_bg,[]), title('\phi_{rs}'), hold on 
for i = 1:h:N
    plot(posrs_y(i,1:+h:end),posrs_x(i,1:+h:end),'r-'), hold on
    plot(posrs_y(1:+h:end,i),posrs_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

[frt, grt] = compute_JD_and_Curl(posrt_x,posrt_y,1);
[frs, grs] = compute_JD_and_Curl(posrs_x,posrs_y,1);
jdnr = nthroot(frt.*frs,3);
cvnr = (grt+grs)/3;
tic
[avgr_x,avgr_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnr, cvnr,N,X,Y);
toc
avgr = interpn(X, Y, R, avgr_x,avgr_y,'makima'); 

% 2nd avg
tic
[posavgrt_x,posavgrt_y,Tn_avgrt,U1avgrt,U2avgrt]=Img_Reg2D_NewDev(avgr,T,N);
toc
tic
[posavgrr_x,posavgrr_y,Tn_avgrr,U1avgrr,U2avgrr]=Img_Reg2D_NewDev(avgr,R,N);
toc
tic
[posavgrs_x,posavgrs_y,Tn_avgrs,U1avgrs,U2avgrs]=Img_Reg2D_NewDev(avgr,S,N);
toc

[favgrt, gavgrt] = compute_JD_and_Curl(posavgrt_x,posavgrt_y,1);
[favgrr, gavgrr] = compute_JD_and_Curl(posavgrr_x,posavgrr_y,1);
[favgrs, gavgrs] = compute_JD_and_Curl(posavgrs_x,posavgrs_y,1);
jdnr1 = nthroot(favgrt.*favgrr.*favgrs,3);
cvnr1 = (gavgrt+gavgrr+gavgrs)/3;

tic
[avgr2_x,avgr2_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnr1, cvnr1,N,X,Y);
toc
avgr2 = interpn(X, Y, avgr, avgr2_x,avgr2_y,'makima'); 

figure(7)
subplot(3,4,1), imshow(T,[]);
subplot(3,4,2), imshow(Tn_avgrt,[]);
subplot(3,4,3), imshow(abs(Tn_avgrt-T),[]);
subplot(3,4,4), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgrt_y(i,1:+h:end),posavgrt_x(i,1:+h:end),'r-'), hold on
    plot(posavgrt_y(1:+h:end,i),posavgrt_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,5), imshow(R,[]);
subplot(3,4,6), imshow(Tn_avgrr,[]);
subplot(3,4,7), imshow(abs(Tn_avgrr-R),[]);
subplot(3,4,8), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgrt_y(i,1:+h:end),posavgrt_x(i,1:+h:end),'r-'), hold on
    plot(posavgrt_y(1:+h:end,i),posavgrt_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,9), imshow(S,[]);
subplot(3,4,10), imshow(Tn_avgrs,[]);
subplot(3,4,11), imshow(abs(Tn_avgrs-S),[]);
subplot(3,4,12), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgrs_y(i,1:+h:end),posavgrs_x(i,1:+h:end),'r-'), hold on
    plot(posavgrs_y(1:+h:end,i),posavgrs_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

tic
[posavgr2t_x,posavgr2t_y,Tn_avgr2t,U1avgr2t,U2avgr2t]=Img_Reg2D_NewDev(avgr2,T,N);
toc
tic
[posavgr2r_x,posavgr2r_y,Tn_avgr2r,U1avgr2r,U2avgr2r]=Img_Reg2D_NewDev(avgr2,R,N);
toc
tic
[posavgr2s_x,posavgr2s_y,Tn_avgr2s,U1avgr2s,U2avgr2s]=Img_Reg2D_NewDev(avgr2,S,N);
toc

figure(8)
subplot(3,4,1), imshow(T,[]);
subplot(3,4,2), imshow(Tn_avgr2t,[]);
subplot(3,4,3), imshow(abs(Tn_avgr2t-T),[]);
subplot(3,4,4), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgr2t_y(i,1:+h:end),posavgr2t_x(i,1:+h:end),'r-'), hold on
    plot(posavgr2t_y(1:+h:end,i),posavgr2t_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,5), imshow(R,[]);
subplot(3,4,6), imshow(Tn_avgr2r,[]);
subplot(3,4,7), imshow(abs(Tn_avgr2r-R),[]);
subplot(3,4,8), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgr2r_y(i,1:+h:end),posavgr2r_x(i,1:+h:end),'r-'), hold on
    plot(posavgr2r_y(1:+h:end,i),posavgr2r_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,9), imshow(S,[]);
subplot(3,4,10), imshow(Tn_avgr2s,[]);
subplot(3,4,11), imshow(abs(Tn_avgr2s-S),[]);
subplot(3,4,12), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgr2s_y(i,1:+h:end),posavgr2s_x(i,1:+h:end),'r-'), hold on
    plot(posavgr2s_y(1:+h:end,i),posavgr2s_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

%form distribution
[favgr2t, gavgr2t] = compute_JD_and_Curl(posavgrt_x,posavgrt_y,1);
[favgr2r, gavgr2r] = compute_JD_and_Curl(posavgrr_x,posavgrr_y,1);
[favgr2s, gavgr2s] = compute_JD_and_Curl(posavgrs_x,posavgrs_y,1);
jdnr2 = nthroot(favgr2t.*favgr2r.*favgr2s,3);
cvnr2 = (gavgr2t+gavgr2r+gavgr2s)/3;
tic
[avgr3_x,avgr3_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnr2, cvnr2,N,X,Y);
toc
avgr3 = interpn(X, Y, avgr2, avgr3_x, avgr3_y,'makima'); 


tic
[posavgr3t_x,posavgr3t_y,Tn_avgr3t,U1avgr3t,U2avgr3t]=Img_Reg2D_NewDev(avgr3,T,N);
toc
tic
[posavgr3r_x,posavgr3r_y,Tn_avgr3r,U1avgr3r,U2avgr3r]=Img_Reg2D_NewDev(avgr3,R,N);
toc
tic
[posavgr3s_x,posavgr3s_y,Tn_avgr3s,U1avgr3s,U2avgr3s]=Img_Reg2D_NewDev(avgr3,S,N);
toc
[favgr3t, gavgr3t] = compute_JD_and_Curl(posavgr3t_x,posavgr3t_y,1);
[favgr3r, gavgr3r] = compute_JD_and_Curl(posavgr3r_x,posavgr3r_y,1);
[favgr3s, gavgr3s] = compute_JD_and_Curl(posavgr3s_x,posavgr3s_y,1);

jdnr3 = nthroot(favgr3t.*favgr3r.*favgr3s,3);
cvnr3 = (gavgr3t+gavgr3r+gavgr3s)/3;
aa=log(favgr3t./jdnr3);
bb=log(favgr3t./jdnr3);
cc=log(favgr3s./jdnr3);
stdev=nthroot((aa.^2+bb.^2+cc.^2)./3,2);
sdjdnr = exp(stdev);
sdcvnr = sqrt((cvnr3-gavgr3t).^2+(cvnr3-gavgr3r).^2+(cvnr3-gavgr3s).^2)/2;

tic
[avgr4_x,avgr4_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdnr3, cvnr3,N,X,Y);
toc
avgr4 = interpn(X, Y, avgr3, avgr4_x, avgr4_y,'makima'); 

figure(9)
subplot(1,4,1), imshow(avgr,[]);
subplot(1,4,2), imshow(avgr2,[]);
subplot(1,4,3), imshow(avgr3,[]);
subplot(1,4,4), imshow(avgr4,[]);

tic
[avgnr1sd_x,avgnr1sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(sdjdnr, sdcvnr,N,X,Y);
toc
avgnr1sd = interpn(X, Y, avgr3, avgnr1sd_x,avgnr1sd_y,'makima'); 
tic
[avgnr2sd_x,avgnr2sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdnr.^2), 2*sdcvnr,N,X,Y);
toc
avgnr2sd = interpn(X, Y, avgr3, avgnr2sd_x,avgnr2sd_y,'makima'); 
tic
[avgnrb1sd_x,avgnrb1sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdnr.^(-1)), -sdcvnr,N,X,Y);
toc
avgnrb1sd = interpn(X, Y, avgr3, avgnrb1sd_x,avgnrb1sd_y,'makima'); 
tic
[avgnrb2sd_x,avgnrb2sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdnr.^(-2)), -2*sdcvnr,N,X,Y);
toc
avgnrb2sd = interpn(X, Y, avgr3, avgnrb2sd_x,avgnrb2sd_y,'makima'); 


figure(10)
h=3;
subplot(2,5,5), imshow(White_bg,[]), title('2\sigma'), hold on 
for i = 1:h:N
    plot(avgnr2sd_y(i,1:+h:end),avgnr2sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnr2sd_y(1:+h:end,i),avgnr2sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,4), imshow(White_bg,[]), title('\sigma'), hold on 
for i = 1:h:N
    plot(avgnr1sd_y(i,1:+h:end),avgnr1sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnr1sd_y(1:+h:end,i),avgnr1sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,3), imshow(White_bg,[]), title('Id'), hold on 
for i = 1:h:N
    plot(Y(i,1:+h:end),X(i,1:+h:end),'r-'), hold on
    plot(Y(1:+h:end,i),X(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,2), imshow(White_bg,[]), title('-\sigma'), hold on
for i = 1:h:N
    plot(avgnrb1sd_y(i,1:+h:end),avgnrb1sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnrb1sd_y(1:+h:end,i),avgnrb1sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,1), imshow(White_bg,[]), title('-2\sigma'), hold on
for i = 1:h:N
    plot(avgnrb2sd_y(i,1:+h:end),avgnrb2sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnrb2sd_y(1:+h:end,i),avgnrb2sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

subplot(2,5,10), imshow(avgnr2sd,[]), title('I_{2\sigma}');
subplot(2,5,9), imshow(avgnr1sd,[]), title('I_{\sigma}');
subplot(2,5,8), imshow(avgr3,[]), title('I_{avg2}');
subplot(2,5,7), imshow(avgnrb1sd,[]), title('I_{-\sigma}');
subplot(2,5,6), imshow(avgnrb2sd,[]), title('I_{-2\sigma}');



%% 3rd avg
tic
[possr_x,possr_y,Tn_sr,U1sr,U2sr]=Img_Reg2D_NewDev(S,R,N);
toc
tic
[posst_x,posst_y,Tn_st,U1st,U2st]=Img_Reg2D_NewDev(S,T,N);
toc

figure(11)
subplot(2,4,1), imshow(R,[]), title('R');
subplot(2,4,2), imshow(Tn_sr,[]), title('T(\phi_{sr})');
subplot(2,4,3), imshow(abs(Tn_sr-R),[]), title('|T(\phi_{sr})-R|');
subplot(2,4,4), imshow(White_bg,[]), title('\phi_{sr}'), hold on 
for i = 1:h:N
    plot(possr_y(i,1:+h:end),possr_x(i,1:+h:end),'r-'), hold on
    plot(possr_y(1:+h:end,i),possr_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,4,5), imshow(T,[]), title('T');
subplot(2,4,6), imshow(Tn_st,[]), title('T(\phi_{st})');
subplot(2,4,7), imshow(abs(Tn_st-T),[]), title('|T(\phi_{st})-T|');
subplot(2,4,8), imshow(White_bg,[]), title('\phi_{st}'), hold on 
for i = 1:h:N
    plot(posst_y(i,1:+h:end),posst_x(i,1:+h:end),'r-'), hold on
    plot(posst_y(1:+h:end,i),posst_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])


[fsr, gsr] = compute_JD_and_Curl(possr_x,possr_y,1);
[fst, gst] = compute_JD_and_Curl(posst_x,posst_y,1);
jdns = nthroot(fsr.*fst,3);
cvns = (gsr+gst)/3;
tic
[avgs_x,avgs_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdns, cvns,N,X,Y);
toc
avgs = interpn(X, Y, S, avgs_x,avgs_y,'makima'); 

% 2nd avg
tic
[posavgst_x,posavgst_y,Tn_avgst,U1avgst,U2avgst]=Img_Reg2D_NewDev(avgs,T,N);
toc
tic
[posavgsr_x,posavgsr_y,Tn_avgsr,U1avgsr,U2avgsr]=Img_Reg2D_NewDev(avgs,R,N);
toc
tic
[posavgss_x,posavgss_y,Tn_avgss,U1avgss,U2avgss]=Img_Reg2D_NewDev(avgs,S,N);
toc

[favgst, gavgst] = compute_JD_and_Curl(posavgst_x,posavgst_y,1);
[favgsr, gavgsr] = compute_JD_and_Curl(posavgsr_x,posavgsr_y,1);
[favgss, gavgss] = compute_JD_and_Curl(posavgss_x,posavgss_y,1);
jdns1 = nthroot(favgst.*favgsr.*favgss,3);
cvns1 = (gavgst+gavgsr+gavgss)/3;

tic
[avgs2_x,avgs2_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdns1, cvns1,N,X,Y);
toc
avgs2 = interpn(X, Y, avgs, avgs2_x,avgs2_y,'makima'); 

figure(12)
subplot(3,4,1), imshow(T,[]);
subplot(3,4,2), imshow(Tn_avgst,[]);
subplot(3,4,3), imshow(abs(Tn_avgst-T),[]);
subplot(3,4,4), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgst_y(i,1:+h:end),posavgst_x(i,1:+h:end),'r-'), hold on
    plot(posavgst_y(1:+h:end,i),posavgst_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,5), imshow(R,[]);
subplot(3,4,6), imshow(Tn_avgsr,[]);
subplot(3,4,7), imshow(abs(Tn_avgsr-R),[]);
subplot(3,4,8), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgsr_y(i,1:+h:end),posavgsr_x(i,1:+h:end),'r-'), hold on
    plot(posavgsr_y(1:+h:end,i),posavgsr_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,9), imshow(S,[]);
subplot(3,4,10), imshow(Tn_avgss,[]);
subplot(3,4,11), imshow(abs(Tn_avgss-S),[]);
subplot(3,4,12), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgss_y(i,1:+h:end),posavgss_x(i,1:+h:end),'r-'), hold on
    plot(posavgss_y(1:+h:end,i),posavgss_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])


tic
[posavgs2t_x,posavgs2t_y,Tn_avgs2t,U1avgs2t,U2avgs2t]=Img_Reg2D_NewDev(avgs2,T,N);
toc
tic
[posavgs2r_x,posavgs2r_y,Tn_avgs2r,U1avgs2r,U2avgs2r]=Img_Reg2D_NewDev(avgs2,R,N);
toc
tic
[posavgs2s_x,posavgs2s_y,Tn_avgs2s,U1avgs2s,U2avgs2s]=Img_Reg2D_NewDev(avgs2,S,N);
toc


figure(13)
subplot(3,4,1), imshow(T,[]);
subplot(3,4,2), imshow(Tn_avgs2t,[]);
subplot(3,4,3), imshow(abs(Tn_avgs2t-T),[]);
subplot(3,4,4), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgs2t_y(i,1:+h:end),posavgs2t_x(i,1:+h:end),'r-'), hold on
    plot(posavgs2t_y(1:+h:end,i),posavgs2t_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,5), imshow(R,[]);
subplot(3,4,6), imshow(Tn_avgs2r,[]);
subplot(3,4,7), imshow(abs(Tn_avgs2r-R),[]);
subplot(3,4,8), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgs2r_y(i,1:+h:end),posavgs2r_x(i,1:+h:end),'r-'), hold on
    plot(posavgs2r_y(1:+h:end,i),posavgs2r_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(3,4,9), imshow(S,[]);
subplot(3,4,10), imshow(Tn_avgs2s,[]);
subplot(3,4,11), imshow(abs(Tn_avgs2s-S),[]);
subplot(3,4,12), imshow(White_bg,[]), hold on 
for i = 1:h:N
    plot(posavgs2s_y(i,1:+h:end),posavgs2s_x(i,1:+h:end),'r-'), hold on
    plot(posavgs2s_y(1:+h:end,i),posavgs2s_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

%form distribution
[favgs2t, gavgs2t] = compute_JD_and_Curl(posavgst_x,posavgst_y,1);
[favgs2r, gavgs2r] = compute_JD_and_Curl(posavgsr_x,posavgsr_y,1);
[favgs2s, gavgs2s] = compute_JD_and_Curl(posavgss_x,posavgss_y,1);

jdns2 = nthroot(favgs2t.*favgs2r.*favgs2s,3);
cvns2 = (gavgs2t+gavgs2r+gavgs2s)/3;
tic
[avgs3_x,avgs3_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdns2, cvns2,N,X,Y);
toc
avgs3 = interpn(X, Y, avgs2, avgs3_x, avgs3_y,'makima'); 


tic
[posavgs3t_x,posavgs3t_y,Tn_avgs3t,U1avgs3t,U2avgs3t]=Img_Reg2D_NewDev(avgs3,T,N);
toc
tic
[posavgs3r_x,posavgs3r_y,Tn_avgs3r,U1avgs3r,U2avgs3r]=Img_Reg2D_NewDev(avgs3,R,N);
toc
tic
[posavgs3s_x,posavgs3s_y,Tn_avgs3s,U1avgs3s,U2avgs3s]=Img_Reg2D_NewDev(avgs3,S,N);
toc
[favgs3t, gavgs3t] = compute_JD_and_Curl(posavgs3t_x,posavgs3t_y,1);
[favgs3r, gavgs3r] = compute_JD_and_Curl(posavgs3r_x,posavgs3r_y,1);
[favgs3s, gavgs3s] = compute_JD_and_Curl(posavgs3s_x,posavgs3s_y,1);

jdns3 = nthroot(favgs3t.*favgs3r.*favgs3s,3);
cvns3 = (gavgs3t+gavgs3r+gavgs3s)/3;
aa=log(favgs3t./jdns3);
bb=log(favgs3t./jdns3);
cc=log(favgs3s./jdns3);
stdev=nthroot((aa.^2+bb.^2+cc.^2)./3,2);
sdjdns = exp(stdev);
sdcvns = sqrt((cvns3-gavgs3t).^2+(cvns3-gavgs3r).^2+(cvns3-gavgs3s).^2)/2;

tic
[avgs4_x,avgs4_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(jdns3, cvns3,N,X,Y);
toc
avgs4 = interpn(X, Y, avgs3, avgs4_x, avgs4_y,'makima'); 

figure(14)
subplot(1,4,1), imshow(avgs,[]);
subplot(1,4,2), imshow(avgs2,[]);
subplot(1,4,3), imshow(avgs3,[]);
subplot(1,4,4), imshow(avgs4,[]);


tic
[avgns1sd_x,avgns1sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(sdjdns, sdcvns,N,X,Y);
toc
avgns1sd = interpn(X, Y, avgs3, avgns1sd_x,avgns1sd_y,'makima'); 
tic
[avgns2sd_x,avgns2sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdns.^2), 2*sdcvns,N,X,Y);
toc
avgns2sd = interpn(X, Y, avgs3, avgns2sd_x,avgns2sd_y,'makima'); 
tic
[avgnsb1sd_x,avgnsb1sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdns.^(-1)), -sdcvns,N,X,Y);
toc
avgnsb1sd = interpn(X, Y, avgs3, avgnsb1sd_x,avgnsb1sd_y,'makima'); 
tic
[avgnsb2sd_x,avgnsb2sd_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast((sdjdns.^(-2)), -2*sdcvns,N,X,Y);
toc
avgnsb2sd = interpn(X, Y, avgs3, avgnsb2sd_x,avgnsb2sd_y,'makima'); 

figure(15)
h=3;
subplot(2,5,5), imshow(White_bg,[]), title('2\sigma'), hold on 
for i = 1:h:N
    plot(avgns2sd_y(i,1:+h:end),avgns2sd_x(i,1:+h:end),'r-'), hold on
    plot(avgns2sd_y(1:+h:end,i),avgns2sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,4), imshow(White_bg,[]), title('\sigma'), hold on 
for i = 1:h:N
    plot(avgns1sd_y(i,1:+h:end),avgns1sd_x(i,1:+h:end),'r-'), hold on
    plot(avgns1sd_y(1:+h:end,i),avgns1sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,3), imshow(White_bg,[]), title('Id'), hold on 
for i = 1:h:N
    plot(Y(i,1:+h:end),X(i,1:+h:end),'r-'), hold on
    plot(Y(1:+h:end,i),X(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,2), imshow(White_bg,[]), title('-\sigma'), hold on
for i = 1:h:N
    plot(avgnsb1sd_y(i,1:+h:end),avgnsb1sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnsb1sd_y(1:+h:end,i),avgnsb1sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])
subplot(2,5,1), imshow(White_bg,[]), title('-2\sigma'), hold on
for i = 1:h:N
    plot(avgnsb2sd_y(i,1:+h:end),avgnsb2sd_x(i,1:+h:end),'r-'), hold on
    plot(avgnsb2sd_y(1:+h:end,i),avgnsb2sd_x(1:+h:end,i),'r-'), hold on
end
axis([1,N,1,N])

subplot(2,5,10), imshow(avgns2sd,[],'border','tight'), title('I_{2\sigma}');
subplot(2,5,9), imshow(avgns1sd,[]), title('I_{\sigma}');
subplot(2,5,8), imshow(avgs3,[]), title('I_{avg3}');
subplot(2,5,7), imshow(avgnsb1sd,[]), title('I_{-\sigma}');
subplot(2,5,6), imshow(avgnsb2sd,[]), title('I_{-2\sigma}');

figure(2021)
imshow(imresize(avgns2sd, [255 255]),[ ]);


%% check KL divegence
threshold=0.5;
[KL_divTR3,countTR]=globalKLdiv(jdnt3,cvnt3,sdjdnt,sdcvnt,jdnr3,cvnr3,sdjdnr,sdcvnr,threshold,N);
[KL_divRS3,countRS]=globalKLdiv(jdnr3,cvnr3,sdjdnr,sdcvnr,jdns3,cvns3,sdjdns,sdcvns,threshold,N);
[KL_divST3,countST]=globalKLdiv(jdns3,cvns3,sdjdns,sdcvns,jdnt3,cvnt3,sdjdnt,sdcvnt,threshold,N);

[KL_divRT3,countRT]=globalKLdiv(jdnr3,cvnr3,sdjdnr,sdcvnr,jdnt3,cvnt3,sdjdnt,sdcvnt,threshold,N);
[KL_divSR3,countSR]=globalKLdiv(jdns3,cvns3,sdjdns,sdcvns,jdnr3,cvnr3,sdjdnr,sdcvnr,threshold,N);
[KL_divTS3,countTS]=globalKLdiv(jdnt3,cvnt3,sdjdnt,sdcvnt,jdns3,cvns3,sdjdns,sdcvns,threshold,N);

maxdiff_KL_TR=max(max(KL_divTR3))
countTR
maxdiff_KL_RS=max(max(KL_divRS3))
countRS
maxdiff_KL_ST=max(max(KL_divST3))
countST

maxdiff_KL_RT=max(max(KL_divRT3))
countRT
maxdiff_KL_SR=max(max(KL_divSR3))
countSR
maxdiff_KL_TS=max(max(KL_divTS3))
countTS

% maxdiff_KL_TR =
%    1.418907605433187
% countTR =
%     11
% maxdiff_KL_RS =
%    0.674951604855725
% countRS =
%      1
% maxdiff_KL_ST =
%    0.541437480368126
% countST =
%      1
% maxdiff_KL_RT =
%    0.719203541024865
% countRT =
%      1
% maxdiff_KL_SR =
%    1.479030643632085
% countSR =
%     10
% maxdiff_KL_TS =
%    0.794047121648590
% countTS =
%      2

avg_KL_TR=sum(sum(KL_divTR3))/(N^2)
avg_KL_RT=sum(sum(KL_divRT3))/(N^2)
% avg_KL_TR =
%      7.615212123200882e-04
% avg_KL_RT =
%      1.014296206889639e-04

avg_KL_RS=sum(sum(KL_divRS3))/(N^2)
avg_KL_SR=sum(sum(KL_divSR3))/(N^2)
% avg_KL_RS =
%      1.389882578966381e-04
% avg_KL_SR =
%      8.477706289912203e-04

avg_KL_ST=sum(sum(KL_divST3))/(N^2)
avg_KL_TS=sum(sum(KL_divTS3))/(N^2)
% avg_KL_ST =
%      1.016346751186368e-04
% avg_KL_TS =
%      1.771604121667297e-04

figure(16)
imshow(imresize(KL_divTR3, [255 255]),[ ])
figure(17)
imshow(imresize(KL_divRS3, [255 255]),[ ])
figure(18)
imshow(imresize(KL_divST3, [255 255]),[ ])

figure(162)
imshow(imresize(KL_divRT3, [255 255]),[ ])
figure(172)
imshow(imresize(KL_divSR3, [255 255]),[ ])
figure(182)
imshow(imresize(KL_divTS3, [255 255]),[ ])

figure(19)
imshow(imresize(abs(avgt4-avgr4), [255 255]),[ ]);
figure(20)
imshow(imresize(abs(avgr4-avgs4), [255 255]),[ ]);
figure(21)
imshow(imresize(abs(avgs4-avgt4), [255 255]),[ ]);



I_raw=imresize(double(rgb2gray(imread('rawsample1.jpg'))), [N N]);
I_seg1=I_raw.*avgt3;
I_seg2=I_raw.*avgr3; 
I_seg3=I_raw.*avgs3;

figure(22)
imshow(imresize(I_seg1, [255 255]),[ ]);
figure(23)
imshow(imresize(I_seg2, [255 255]),[ ]);
figure(24)
imshow(imresize(I_seg3, [255 255]),[ ]);

figure(25)
imshow(imresize(avgt3, [255 255]),[ ]);
figure(26)
imshow(imresize(avgr3, [255 255]),[ ]);
figure(27)
imshow(imresize(avgs3, [255 255]),[ ]);
sum(sum(jdns3))
%% draw samples
rdm=randn(1,7);
deno=max(abs(rdm));
rdm=(rdm/deno)*0.72;
I_raw=imresize(double(rgb2gray(imread('rawsample1.jpg'))), [N N]);

porp=rdm(1)
tic
[avgnsrand_x,avgnsrand_y,~,~,~, ~,~,~]=PJDC_on_given_mesh2fast(sdjdns.^(0.75+porp), (0.75+porp)*sdcvns,N,X,Y);
toc
avgnsrand = interpn(X, Y, avgs3, avgnsrand_x,avgnsrand_y,'makima'); 
figure(100+i)
imshow(imresize(avgnsrand, [255 255]),[ ]);
figure(200+i)
imshow(imresize(avgnsrand.*I_raw, [255 255]),[ ]);
