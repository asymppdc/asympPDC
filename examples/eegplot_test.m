figure
nPoints = 2000;
nDiscard = 1;
P_TIME0 = 1;
V_TIME0 = P_TIME0;
V_TIME1 = nPoints

clear ID_CHANNEL

Ch= {'x1'; 'x2'; 'x3'; 'x4';'x5'};

Yout = fbaccala2001a_ex3( nPoints, nDiscard );

%[a b] =  size(Yout);

nCh =  min(size(Yout));

V_FA = 1;

T_obs='Testing';
T=1/V_FA;

N=length(Yout(:,1));
t = P_TIME0:T:(P_TIME0+(N-1)*T); % x-axis scale

for Kch=1:nCh,
   subplot(nCh+1,1,Kch);
   plot(t,Yout(:,Kch))
   H_CHPLOT(Kch)=gca;
   axis([V_TIME0 V_TIME1 min(Yout(:,Kch)) max(Yout(:,Kch))]);
   set(gca,'Fontsize', [ 8 ]);
   YMIN(Kch+1)=min(Yout(:,Kch));
   YMAX(Kch+1)=max(Yout(:,Kch));
   disp('YMIN and YMAX updated 2.')
   eegplot2;
end;

Kch=Kch+1;
subplot(2*Kch,1,2*Kch-1)
plot([V_TIME0 V_TIME1],[0 0],'color',[1 1 1])
%axis([V_TIME0 V_TIME1 -1 1]);
axis([0 V_TIME1 -1 1]);
eegplot2

shg
