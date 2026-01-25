* Decadimenti del B in due pseudoscalari
* In questa versione c`e' un editing finale che permette una rapido
* uso dell'output secondo la parametrizzazione di Neubert e 
* per la fattorizzazione.
*#-
Off stat;

S LL,LR,SP;
S i,j,k,o,[Sqr(2)],[Sqr(3)],[Sqr(6)],A,lam,sig,[_*Exp[I*del]*_];
S Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb;
S Vdu,Vsu,Vbu,Vdc,Vsc,Vbc,Vdt,Vst,Vbt;
S C1,C2,C3,C4,C5,C6,C7,C8,C9,C10;
S O1,O2,O3,O4,O5,O6,O7,O8,O9,O10,O1C,O2C;
S O1P,O2P,O1CP,O2CP,O3P,O4P,O5P,O6P,O7P,O8P,O9P,O10P,Obd,Obs;
S P1,P2,P1C,P2C,P1P,P2P,P1CP,P2CP;
S I00,I20,x,y,z;
S P,B,K,Kb,D,Db,R,KS,KbS,DS,DbS,P8,V8;

CF S,SB,Conj,GDP,GCP,GDPE,GCPE,GDPA,GCPA,SS,fs,fq;
CF t1,aa1;
F aq,aqb,q,qb,f,fp,In,Out;
F u,d,c,s,b;
F ub,db,cb,sb,bb;
F D0,D0b,DPl,DM,DsP,DsM;
F D0S,D0bS,DPS,DMS,DsPS,DsMS;
F Bd,Bs,BP;
F PP,P0,PM,KP,KM,K0,K0b,ET1,ET8,ETc;
F RP,R0,RM,KPS,KMS,K0S,K0bS,PH,OM,JPSI,OM1,OM8;
CF CP,DP,CPA,DPA,CPAB,DPAB,CPE,DPE,DEA,CEA,DCA,DDA;
CF CE,DE,CA,DA,DELL,DELR,DESP;
F e1,e2,e1t,e2t;
F a1,a2,ea1,ea2,a1t,a2t,ea1,ea2t;
F p1,p2,p3,p4,dp1ew,dp2ew,dp3ew,dp4ew,p1t,p2t,p3t,p4t;
F gim1,gim2,gim3,gim4,g1t,g2t,g3t,g4t;

CF Blob;

.global

* Scrivo l'H_eff come nei miei appunti ( escluso -G_F/sqrt(2) )

L H=Vbu*Vud*(C1*O1 + C2*O2)
   +Vbc*Vcd*(C1*O1C + C2*O2C)
   -Vbt*Vtd*(C3*O3 + C4*O4 + C5*O5 + C6*O6 + C7*O7 + C8*O8 + C9*O9 + C10*O10 )
   +Vbu*Vus*(C1*O1P + C2*O2P)
   +Vbc*Vcs*(C1*O1CP + C2*O2CP)
   -Vbt*Vts*(C3*O3P + C4*O4P + C5*O5P + C6*O6P + C7*O7P + C8*O8P + 
	C9*O9P + C10*O10P )
   +Vbu*Vcd*(C1*P1 + C2*P2)
   +Vbc*Vud*(C1*P1C + C2*P2C)
   +Vbu*Vcs*(C1*P1P + C2*P2P)
   +Vbc*Vus*(C1*P1CP + C2*P2CP)
   -Vbt*Vtd*Obd-Vbt*Vts*Obs;


Id,O1=bb(0)*d(0)*ub(o)*u(o)*LL;
Id,O1C=bb(0)*d(0)*cb(o)*c(o)*LL;
id O1P=bb(0)*s(0)*ub(o)*u(o)*LL;
Id,O1CP=bb(0)*s(0)*cb(o)*c(o)*LL;
Id,O2=bb(0)*u(0)*ub(o)*d(o)*LL;
Id,O2C=bb(0)*c(0)*cb(o)*d(o)*LL;
id O2P=bb(0)*u(0)*ub(o)*s(o)*LL;
Id,O2CP=bb(0)*c(0)*cb(o)*s(o)*LL;
Id,O3=bb(0)*d(0)*(ub(o)*u(o)+db(o)*d(o)+cb(o)*c(o)+sb(o)*s(o))*LL;
Id,O4=bb(0)*(u(0)*ub(o)+d(0)*db(o)+c(0)*cb(o)+s(0)*sb(o))*d(o)*LL;
Id,O5=bb(0)*d(0)*(ub(o)*u(o)+db(o)*d(o)+cb(o)*c(o)+sb(o)*s(o))*LR;
Id,O6=-2*( bb(0)*(u(0)*ub(o)+d(0)*db(o)+c(0)*cb(o)+s(0)*sb(o))*d(o) )*SP;
Id,O3P=bb(0)*s(0)*(ub(o)*u(o)+db(o)*d(o)+cb(o)*c(o)+sb(o)*s(o))*LL;
Id,O4P=bb(0)*(u(0)*ub(o)+d(0)*db(o)+c(0)*cb(o)+s(0)*sb(o))*s(o)*LL;
Id,O5P=bb(0)*s(0)*(ub(o)*u(o)+db(o)*d(o)+cb(o)*c(o)+sb(o)*s(o))*LR;
Id,O6P=-2*( bb(0)*(u(0)*ub(o)+d(0)*db(o)+c(0)*cb(o)+s(0)*sb(o))*s(o) )*SP;
Id,P1=bb(0)*d(0)*cb(o)*u(o)*LL;
Id,P2=bb(0)*u(0)*cb(o)*d(o)*LL;
Id,P1C=bb(0)*d(0)*ub(o)*c(o)*LL;
Id,P2C=bb(0)*c(0)*ub(o)*d(o)*LL;
Id,P1P=bb(0)*s(0)*cb(o)*u(o)*LL;
Id,P2P=bb(0)*u(0)*cb(o)*s(o)*LL;
Id,P1CP=bb(0)*s(0)*ub(o)*c(o)*LL;
Id,P2CP=bb(0)*c(0)*ub(o)*s(o)*LL;
*Id,Obd=(C5+C6)*bb(0)*d(0);
*Id,Obs=(C5+C6)*bb(0)*s(0);
Id,Obd=0;
Id,Obs=0;


* Elettropinguini
Id,O7=3/2*bb(0)*d(0)*(2/3*ub(o)*u(o)-1/3*db(o)*d(o)+
     2/3*cb(o)*c(o)-1/3*sb(o)*s(o))*LR;
Id,O8=-2*3/2*bb(0)*(2/3*u(0)*ub(o)-1/3*d(0)*db(o)+
     2/3*c(0)*cb(o)-1/3*s(0)*sb(o))*d(o)*SP;
Id,O9=3/2*bb(0)*d(0)*(2/3*ub(o)*u(o)-1/3*db(o)*d(o)+
     2/3*cb(o)*c(o)-1/3*sb(o)*s(o))*LL;
Id,O10=3/2*bb(0)*(2/3*u(0)*ub(o)-1/3*d(0)*db(o)+
     2/3*c(0)*cb(o)-1/3*s(0)*sb(o))*d(o)*LL;
Id,O7P=3/2*bb(0)*s(0)*(2/3*ub(o)*u(o)-1/3*db(o)*d(o)+
     2/3*cb(o)*c(o)-1/3*sb(o)*s(o))*LR;
Id,O8P=-2*3/2*bb(0)*(2/3*u(0)*ub(o)-1/3*d(0)*db(o)+
     2/3*c(0)*cb(o)-1/3*s(0)*sb(o))*s(o)*SP;
Id,O9P=3/2*bb(0)*s(0)*(2/3*ub(o)*u(o)-1/3*db(o)*d(o)+
     2/3*cb(o)*c(o)-1/3*sb(o)*s(o))*LL;
Id,O10P=3/2*bb(0)*(2/3*u(0)*ub(o)-1/3*d(0)*db(o)+
     2/3*c(0)*cb(o)-1/3*s(0)*sb(o))*s(o)*LL;

* Sostituisco gli angoli CKM nella parametrizzazione di wolfenstein
* tenendo i termini fino a lam^3 inclusi

* Unitarieta'
Id,Vcs*Vbc=-Vts*Vbt-Vus*Vbu;
Id,Vcd*Vbc=-Vtd*Vbt-Vud*Vbu;

*Id,Vub=lam^3*Vub;
*Id,Vbu=lam^3*Vbu;
*Id,Vcb=lam^2*Vcb;
*Id,Vbc=lam^2*Vbc;
*Id,Vcd=lam*Vcd;
*Id,Vdc=lam*Vdc;
*Id,Vtd=lam^3*Vtd;
*Id,Vdt=lam^3*Vdt;
*Id,Vts=lam^2*Vts;
*Id,Vst=lam^2*Vst;
*Id,Vus=lam*Vus;
*Id,Vsu=lam*Vsu;

.sort

* Ora scrivo le ampiezze per i vari decadimenti

* + -> + 0

#do i={BP}
#do j={PP|KP|DPl|DsP|RP|KPS|DPS|DsPS}
*#do j={PP|KP|DPl|DsP}
#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|OM1|OM8|PH|OM|D0S|D0bS|JPSI}
*#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|PH|OM|D0S|D0bS|JPSI}
*#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc}


G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*Out;

#enddo
#enddo
#enddo

.sort

* Bd -> + -

* #do i={Bs}
#do i={Bd}
#do j={PP|KP|DPl|DsP|RP|KPS|DPS|DsPS}
*#do j={PP|KP|DPl|DsP}
#do k={PM|KM|DM|DsM|RM|KMS|DMS|DsMS}
*#do k={PM|KM|DM|DsM}

G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*Out;

#enddo
#enddo
#enddo

.sort


* Bd -> 0 0

*#do i={Bs}
#do i={Bd}
#Define jj "0"
#do j={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|OM1|OM8|PH|OM|D0S|D0bS|JPSI}
*#do j={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|PH|OM|D0S|D0bS|JPSI}
*#do j={P0|K0|K0b|D0|D0b|PH|JPSI}
*#do j={P0|K0|K0b|ET1|ET8|D0|D0b|ETc}
#Define jj1 "{`jj'+1}"
#Define jj "`jj1'"
#Define kk "0"
#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|OM1|OM8|PH|OM|D0S|D0bS|JPSI}
*#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|PH|OM|D0S|D0bS|JPSI}
*#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc}
*#do k={P0|K0|K0b|D0|D0b|PH|JPSI}
#Define kk1 "{`kk'+1}"
#Define kk "`kk1'"
#If `kk' > `jj' 

G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*Out;

#endif

#If `kk' = `jj' 

G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*
		Out/[Sqr(2)];

#endif
#enddo
#enddo
#enddo

.sort

* Bs -> + -

#do i={Bs}
*#do i={Bd}
#do j={PP|KP|DPl|DsP|RP|KPS|DPS|DsPS}
*#do j={PP|KP|DPl|DsP}
#do k={PM|KM|DM|DsM|RM|KMS|DMS|DsMS}
*#do k={PM|KM|DM|DsM}

G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*Out;

#enddo
#enddo
#enddo

.sort

drop H;

* Bs -> 0 0

#do i={Bs}
*#do i={Bd}
#Define jj "0"
#do j={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|OM1|OM8|PH|OM|D0S|D0bS|JPSI}
*#do j={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|PH|OM|D0S|D0bS|JPSI}
*#do j={P0|K0|K0b|ET1|ET8|D0|D0b|ETc}
*#do j={P0|K0|K0b|D0|D0b|PH|JPSI}
#Define jj1 "{`jj'+1}"
#Define jj "`jj1'"
#Define kk "0"
#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|OM1|OM8|PH|OM|D0S|D0bS|JPSI}
*#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc|R0|K0S|K0bS|PH|OM|D0S|D0bS|JPSI}
*#do k={P0|K0|K0b|ET1|ET8|D0|D0b|ETc}
*#do k={P0|K0|K0b|D0|D0b|PH|JPSI}
#Define kk1 "{`kk'+1}"
#Define kk "`kk1'"
#If `kk' > `jj' 

G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*Out;

#endif
#If `kk' = `jj' 

G Am`i'`j'`k'(`i',`j',`k')=In*`i'(1)*H*`j'(2)*`k'(3)*fs(`i',`j',`k')*
	Out/[Sqr(2)];

#endif
#enddo
#enddo
#enddo

* Adesso esprimo i vari mesoni in termini di quarks

Id,BP(i?)=ub(i)*b(i);
Id,Bd(i?)=db(i)*b(i);
Id,Bs(i?)=sb(i)*b(i);
Id,KP(i?)=u(i)*sb(i);
Id KPS(i?)=u(i)*sb(i);
Id,KM(i?)=-s(i)*ub(i);
Id,KMS(i?)=-s(i)*ub(i);
Id,K0(i?)=d(i)*sb(i);
Id,K0S(i?)=d(i)*sb(i);
Id,K0b(i?)=s(i)*db(i);
Id,K0bS(i?)=s(i)*db(i);
Id,PP(i?)=u(i)*db(i);
Id,RP(i?)=u(i)*db(i);
Id,PM(i?)=-d(i)*ub(i);
Id,RM(i?)=-d(i)*ub(i);
Id,P0(i?)=(d(i)*db(i)-u(i)*ub(i))/[Sqr(2)];
Id,R0(i?)=(d(i)*db(i)-u(i)*ub(i))/[Sqr(2)];
Id,ET1(i?)=(s(i)*sb(i)+u(i)*ub(i)+d(i)*db(i))/[Sqr(3)];
Id,OM1(i?)=(s(i)*sb(i)+u(i)*ub(i)+d(i)*db(i))/[Sqr(3)];
Id,PH(i?)=-s(i)*sb(i);
Id,ET8(i?)=-(2*s(i)*sb(i)-u(i)*ub(i)-d(i)*db(i))/[Sqr(6)];
Id,OM8(i?)=-(2*s(i)*sb(i)-u(i)*ub(i)-d(i)*db(i))/[Sqr(6)];
Id,OM(i?)=(u(i)*ub(i)+d(i)*db(i))/[Sqr(2)];
Id DsP(i?)=c(i)*sb(i);
Id DsPS(i?)=c(i)*sb(i);
Id DsM(i?)=s(i)*cb(i);
Id DsMS(i?)=s(i)*cb(i);
Id DPl(i?)=c(i)*db(i);
Id DPS(i?)=c(i)*db(i);
Id DM(i?)=d(i)*cb(i);
Id DMS(i?)=d(i)*cb(i);
Id D0(i?)=-c(i)*ub(i);
Id D0S(i?)=-c(i)*ub(i);
Id,D0b(i?)=u(i)*cb(i);
Id,D0bS(i?)=u(i)*cb(i);
Id ETc(i?)=c(i)*cb(i);
Id JPSI(i?)=c(i)*cb(i);

* Per Prima cosa rinomino i campi di quark

#do f={u|d|c|s|b}
  Id,`f'(i?)=q(`f',i);
  Id,`f'b(i?)=qb(`f',i);
#enddo

*print;
.sort

repeat;
  Id,In*q(f?,i?)=In*aq(f,i);
  repeat;
    Id,aq(f?,i?)*qb(f?,j?)=S(f,i,j)-qb(f,j)*aq(f,i);
  endrepeat;
  Id,aq(f?,i?)*qb(fp?,j?)=-qb(fp,j)*aq(f,i);
  Id,aq(f?,i?)*q(fp?,j?)=-q(fp,j)*aq(f,i);
  Id,aq(f?,i?)*Out=0;

  Id,In*qb(f?,i?)=In*aqb(f,i);
  repeat;
    Id,aqb(f?,i?)*q(f?,j?)=-S(f,j,i)-q(f,j)*aqb(f,i);
  endrepeat;
  Id,aqb(f?,i?)*q(fp?,j?)=-q(fp,j)*aqb(f,i);
  Id,aqb(f?,i?)*qb(fp?,j?)=-qb(fp,j)*aqb(f,i);
  Id,aqb(f?,i?)*Out=0;
endrepeat;

Id,In*Out=1;

*print;
.sort

Id,S(b,i?,j?)=SB(i,j);
Id,S(f?,i?,i?)=Blob(f,i,i);
Id,S(f?,0,o)=Blob(f,0,o);
Id,S(f?,o,0)=Blob(f,o,0);
*Id,S(c,i?,j?)=SC(i,j);
*Id,S(f?,i?,j?)=S(i,j);
*Id,Blob(d,i?,j?)=Blob(u,i,j);
*Id,Blob(s,i?,j?)=Blob(u,i,j);
Id,[Sqr(2)]^-2=1/2;
Id,[Sqr(6)]^-2=1/6;
Id,[Sqr(2)]^-1*[Sqr(3)]^-1=[Sqr(6)]^-1;
Id,[Sqr(2)]^-1*[Sqr(6)]^-1=[Sqr(3)]^-1/2;
Id,[Sqr(3)]^-1*[Sqr(6)]^-1=[Sqr(2)]^-1/3;

*print;
.sort

* Identifico le differenti topologie delle contrazioni

Id,SB(1,0)*S(u?,0,i?)*S(d?,i?,o)*S(s?,o,j?)*S(c?,j?,1)=-CE(u,d,s,c,i);
Id,SB(1,0)*S(u?,0,j?)*S(d?,j?,1)*S(s?,o,i?)*S(c?,i?,o)=DE(s,c,u,d,i);
Id,SB(1,0)*S(u?,0,i?)*S(d?,i?,j?)*S(s?,j?,o)*S(c?,o,1)=-CA(u,d,s,c,i);
Id,SB(1,0)*S(u?,0,1)*S(d?,o,i?)*S(s?,i?,j?)*S(c?,j?,o)=DA(d,s,c,u,i);
Id,SB(1,0)*Blob(f?,0,o)*S(u?,o,i?)*S(d?,i?,j?)*S(s?,j?,1)=-CP(f,u,d,s,i);
Id,SB(1,0)*S(u?,0,i?)*S(d?,i?,j?)*S(s?,j?,1)*Blob(f?,o,o)=DP(f,u,d,s,i);

id DE(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DE(u,d,s,c)*fs(z,x,y);
id DE(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DE(u,d,s,c)*fs(z,y,x);
id CE(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CE(u,d,s,c)*fs(z,x,y);
id CE(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CE(u,d,s,c)*fs(z,y,x);

id DP(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DP(u,d,s,c)*fs(z,x,y);
id DP(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DP(u,d,s,c)*fs(z,y,x);
id CP(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CP(u,d,s,c)*fs(z,x,y);
id CP(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CP(u,d,s,c)*fs(z,y,x);

id DA(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DA(u,d,s,c)*fs(z,x,y);
id DA(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DA(u,d,s,c)*fs(z,y,x);
id CA(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CA(u,d,s,c)*fs(z,x,y);
id CA(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CA(u,d,s,c)*fs(z,y,x);

*  cose incalcolabili

Id,SB(1,0)*Blob(f?,0,o)*S(u?,o,1)*S(d?,i?,j?)*S(s?,j?,i?)=CPA(f,d,s,u,i);
Id,SB(1,0)*Blob(f?,o,o)*S(u?,0,1)*S(d?,i?,j?)*S(s?,j?,i?)=-DPA(f,d,s,u,i);
Id,SB(1,0)*Blob(f?,0,o)*S(u?,o,i?)*Blob(fp?,j?,j?)*S(d?,i?,1)=CPE(f,fp,u,d,j);
Id,SB(1,0)*S(u?,0,i?)*Blob(fp?,j?,j?)*S(d?,i?,1)*Blob(f?,o,o)=-DPE(f,fp,u,d,j);
Id,SB(1,0)*S(u?,0,1)*S(d?,o,i?)*Blob(f?,j?,j?)*S(s?,i?,o)=-DEA(d,s,f,u,i);
Id,SB(1,0)*S(u?,0,i?)*Blob(f?,j?,j?)*S(d?,i?,o)*S(s?,o,1)=CEA(u,d,f,s,i);
Id,SB(1,0)*S(u?,o,1)*Blob(c?,0,o)*Blob(d?,i?,i?)*Blob(s?,j?,j?)=-CPAB(c,d,s,u,i);
Id,SB(1,0)*S(u?,0,1)*Blob(c?,o,o)*Blob(d?,i?,i?)*Blob(s?,j?,j?)=DPAB(c,d,s,u,i);

*Id,SB(1,0)*S(u?,0,i?)*S(c?,i?,j?)*S(d?,j?,1)*fs(z?,x?,y?)=
* -eta*DE(u,d,u,z,x,y);
*Id,SB(1,0)*S(d?,0,1)*S(u?,i?,j?)*S(s?,j?,i?)=0;
*Id,SB(1,0)*S(d?,0,i?)*S(s?,i?,1)*Blob(u?,j?,j?)=0;
*Id,SB(1,0)*S(d?,0,1)*Blob(s?,i?,i?)*Blob(u?,j?,j?)=0;

id DPE(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DPE(u,d,s,c)*fs(z,x,y);
id DPE(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DPE(u,d,s,c)*fs(z,y,x);
id CPE(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CPE(u,d,s,c)*fs(z,x,y);
id CPE(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CPE(u,d,s,c)*fs(z,y,x);

id DEA(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DEA(u,d,s,c)*fs(z,x,y);
id DEA(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DEA(u,d,s,c)*fs(z,y,x);
id CEA(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CEA(u,d,s,c)*fs(z,x,y);
id CEA(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CEA(u,d,s,c)*fs(z,y,x);

id DPA(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DPA(u,d,s,c)*fs(z,x,y);
id DPA(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DPA(u,d,s,c)*fs(z,y,x);
id CPA(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CPA(u,d,s,c)*fs(z,x,y);
id CPA(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CPA(u,d,s,c)*fs(z,y,x);

id DPAB(u?,d?,s?,c?,2)*fs(z?,x?,y?)=DPAB(u,d,s,c)*fs(z,x,y);
id DPAB(u?,d?,s?,c?,3)*fs(z?,x?,y?)=DPAB(u,d,s,c)*fs(z,y,x);
id CPAB(u?,d?,s?,c?,2)*fs(z?,x?,y?)=CPAB(u,d,s,c)*fs(z,x,y);
id CPAB(u?,d?,s?,c?,3)*fs(z?,x?,y?)=CPAB(u,d,s,c)*fs(z,y,x);

* B C1,C2,C3,C4,C5,C6,C9,C10;
* B DELL,DELR,DESP,DA;
B lam,Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb,
  Vdu,Vsu,Vbu,Vdc,Vsc,Vbc,Vdt,Vst,Vbt,LL,LR,SP,[Sqr(2)],[Sqr(3)],[Sqr(6)],fs;

*print;
.sort

* Esprimo tutto in termini dei parametri fattorizzati

id DE(c?,u?,s?,d?)*C2=e1(c,u,s,d)-CE(c,u,s,d)*C1;
id DE(c?,u?,s?,d?)*C1=e2(c,u,s,d)-CE(c,u,s,d)*C2;

id DA(c?,u?,s?,d?)*C2=a1(c,u,s,d)-CA(c,u,s,d)*C1;
id DA(c?,u?,s?,d?)*C1=a2(c,u,s,d)-CA(c,u,s,d)*C2;

id DEA(c?,u?,s?,d?)*C2=ea1(c,u,s,d)-CEA(c,u,s,d)*C1;
id DEA(c?,u?,s?,d?)*C1=ea2(c,u,s,d)-CEA(c,u,s,d)*C2;

id CP(c,f?,b?,fp?)*C2*Vbt= (p1(f,b,fp)
                  -DP(c,f,b,fp)*C1-CE(f,b,b,fp)*C3-CE(f,b,b,fp)*C5*LR/LL-
		  DE(f,b,b,fp)*C4-DE(f,b,b,fp)*C6*(-2*SP)/LL+fq(b)*(
		  -CE(f,b,b,fp)*C9-CE(f,b,b,fp)*C7*LR/LL-
		  DE(f,b,b,fp)*C10-DE(f,b,b,fp)*C8*(-2*SP)/LL)-
		  ((DP(u,f,b,fp)+DP(d,f,b,fp)+CP(f,f,b,fp)+DP(s,f,b,fp)+DP(c,f,b,fp))*C3+
		  (DP(u,f,b,fp)+DP(d,f,b,fp)+CP(f,f,b,fp)+DP(s,f,b,fp)+DP(c,f,b,fp))*C5*LR/LL+
		  (CP(u,f,b,fp)+CP(d,f,b,fp)+DP(f,f,b,fp)+CP(s,f,b,fp)+CP(c,f,b,fp))*C4+
		  (CP(u,f,b,fp)+CP(d,f,b,fp)+DP(f,f,b,fp)+CP(s,f,b,fp)+CP(c,f,b,fp))*C6*
			(-2*SP)/LL)-
                  ((DP(u,f,b,fp)-1/2*DP(d,f,b,fp)-1/2*CP(f,f,b,fp)-1/2*DP(s,f,b,fp)+
			DP(c,f,b,fp))*C9+
		  (DP(u,f,b,fp)-1/2*DP(d,f,b,fp)-1/2*CP(f,f,b,fp)-1/2*DP(s,f,b,fp)+
			DP(c,f,b,fp))*C7*LR/LL+
		  (CP(u,f,b,fp)-1/2*CP(d,f,b,fp)-1/2*DP(f,f,b,fp)-1/2*CP(s,f,b,fp)+
			CP(c,f,b,fp))*C10+
		  (CP(u,f,b,fp)-1/2*CP(d,f,b,fp)-1/2*DP(f,f,b,fp)-1/2*CP(s,f,b,fp)+
			CP(c,f,b,fp))*C8*
			(-2*SP)/LL)+fq(fp)*(-
		  CA(f,b,fp,fp)*C9-DA(f,b,fp,fp)*C10-CA(f,b,fp,fp)*C7*LR/LL-
                  DA(f,b,fp,fp)*C8*(-2*SP)/LL) -
		  CA(f,b,fp,fp)*C3-DA(f,b,fp,fp)*C4-CA(f,b,fp,fp)*C5*LR/LL-
                  DA(f,b,fp,fp)*C6*(-2*SP)/LL)*Vbt;

id DP(u,f?,b?,fp?)*C1=-gim1(f,b,fp)+DP(c,f,b,fp)*C1 
		  - (CP(u,f,b,fp)-CP(c,f,b,fp))*C2;

id CPE(c,f?,b?,fp?)*fs(Bd?,K0?,KP?)*C2*Vbt= ((p2(b,f,fp)
                  -DPE(c,f,b,fp)*C1
		  -DE(f,f,b,fp)*C3-DE(f,f,b,fp)*C5*LR/LL-
		  CE(f,f,b,fp)*C4-CE(f,f,b,fp)*C6*(-2*SP)/LL+fq(f)*(
                  -DE(f,f,b,fp)*C9-DE(f,f,b,fp)*C7*LR/LL-
		  CE(f,f,b,fp)*C10-CE(f,f,b,fp)*C8*(-2*SP)/LL))*fs(Bd,K0,KP)+ 
		  (-CEA(b,fp,f,fp)*C3-CEA(b,fp,f,fp)*C5*LR/LL-
		  DEA(b,fp,f,fp)*C4-DEA(b,fp,f,fp)*C6*(-2*SP)/LL+fq(fp)*(
                  -CEA(b,fp,f,fp)*C9-CEA(b,fp,f,fp)*C7*LR/LL-
		  DEA(b,fp,f,fp)*C10-DEA(b,fp,f,fp)*C8*(-2*SP)/LL))*fs(Bd,KP,K0) -
		  (((DPE(u,f,b,fp)+DPE(d,f,b,fp)+CPE(b,f,b,fp)+DPE(s,f,b,fp)+DPE(c,f,b,fp))*
			C3+
		  (DPE(u,f,b,fp)+DPE(d,f,b,fp)+CPE(b,f,b,fp)+DPE(s,f,b,fp)+DPE(c,f,b,fp))*
			C5*LR/LL+
		  (CPE(u,f,b,fp)+CPE(d,f,b,fp)+DPE(b,f,b,fp)+CPE(s,f,b,fp)+CPE(c,f,b,fp))*C4+
		  (CPE(u,f,b,fp)+CPE(d,f,b,fp)+DPE(b,f,b,fp)+CPE(s,f,b,fp)+CPE(c,f,b,fp))*C6*
		  (-2*SP)/LL+
                  (DPE(u,f,b,fp)-1/2*DPE(d,f,b,fp)-1/2*CPE(b,f,b,fp)-1/2*DPE(s,f,b,fp)+
			DPE(c,f,b,fp))*C9+
		  (DPE(u,f,b,fp)-1/2*DPE(d,f,b,fp)-1/2*CPE(b,f,b,fp)-1/2*DPE(s,f,b,fp)+
			DPE(c,f,b,fp))*C7*LR/LL+
		  (CPE(u,f,b,fp)-1/2*CPE(d,f,b,fp)-1/2*DPE(b,f,b,fp)-1/2*CPE(s,f,b,fp)+
			CPE(c,f,b,fp))*C10+
		  (CPE(u,f,b,fp)-1/2*CPE(d,f,b,fp)-1/2*DPE(b,f,b,fp)-1/2*CPE(s,f,b,fp)+
			CPE(c,f,b,fp))*C8*
			(-2*SP)/LL))*fs(Bd,K0,KP))*Vbt;

id DPE(u,f?,b?,fp?)*C1=-gim2(b,f,fp)+DPE(c,f,b,fp)*C1 
		  - (CPE(u,f,b,fp)-CPE(c,f,b,fp))*C2;

id fs(Bd?,KP?,KM?)*CPA(c,f?,b?,fp?)*C2*Vbt= (p3(f,b,fp)*fs(Bd,KP,KM)-
		  DPA(c,f,b,fp)*C1*fs(Bd,KP,KM)-
		  (DA(f,b,f,fp)*C3 + DA(f,b,f,fp)*C5*LR/LL +
		  CA(f,b,f,fp)*C4 +CA(f,b,f,fp)*C6*(-2*SP)/LL +fq(f)*(
                  DA(f,b,f,fp)*C9 + DA(f,b,f,fp)*C7*LR/LL +
		  CA(f,b,f,fp)*C10 +CA(f,b,f,fp)*C8*(-2*SP)/LL))*fs(Bd,KM,KP)-
		  (DA(b,f,b,fp)*C3 + DA(b,f,b,fp)*C5*LR/LL +
		  CA(b,f,b,fp)*C4 +CA(b,f,b,fp)*C6*(-2*SP)/LL +fq(b)*(
                  DA(b,f,b,fp)*C9 + DA(b,f,b,fp)*C7*LR/LL +
		  CA(b,f,b,fp)*C10 +CA(b,f,b,fp)*C8*(-2*SP)/LL))*fs(Bd,KP,KM)-
		  ((DPA(u,f,b,fp)+DPA(d,f,b,fp)+CPA(fp,f,b,fp)+DPA(s,f,b,fp)+DPA(c,f,b,fp))*C3+
		  (DPA(u,f,b,fp)+DPA(d,f,b,fp)+CPA(fp,f,b,fp)+DPA(s,f,b,fp)+DPA(c,f,b,fp))*
			C5*LR/LL+
		  (CPA(u,f,b,fp)+CPA(d,f,b,fp)+DPA(fp,f,b,fp)+CPA(s,f,b,fp)+CPA(c,f,b,fp))*C4+
		  (CPA(u,f,b,fp)+CPA(d,f,b,fp)+DPA(fp,f,b,fp)+CPA(s,f,b,fp)+CPA(c,f,b,fp))*C6*
		  (-2*SP)/LL +
                  (DPA(u,f,b,fp)-1/2*DPA(d,f,b,fp)-1/2*CPA(fp,f,b,fp)-1/2*DPA(s,f,b,fp)+
			DPA(c,f,b,fp))*C9+
		  (DPA(u,f,b,fp)-1/2*DPA(d,f,b,fp)-1/2*CPA(fp,f,b,fp)-1/2*DPA(s,f,b,fp)+
			DPA(c,f,b,fp))*C7*LR/LL+
		  (CPA(u,f,b,fp)-1/2*CPA(d,f,b,fp)-1/2*DPA(fp,f,b,fp)-1/2*CPA(s,f,b,fp)+
			CPA(c,f,b,fp))*C10+
		  (CPA(u,f,b,fp)-1/2*CPA(d,f,b,fp)-1/2*DPA(fp,f,b,fp)-1/2*CPA(s,f,b,fp)+
			CPA(c,f,b,fp))*C8*
			(-2*SP)/LL)*fs(Bd,KP,KM))*Vbt;

id DPA(u,f?,b?,fp?)*C1=-gim3(f,b,fp)+DPA(c,f,b,fp)*C1 
		  - (CPA(u,f,b,fp)-CPA(c,f,b,fp))*C2;

id fs(Bd?,KP?,KM?)*CPAB(c,f?,b?,fp?)*C2*Vbt= (p4(f,b,fp)*fs(Bd,KP,KM)-
		  DPAB(c,f,b,fp)*C1*fs(Bd,KP,KM)-
		  (DEA(f,f,b,fp)*C3 + DEA(f,f,b,fp)*C5*LR/LL +
		  CEA(f,f,b,fp)*C4 +CEA(f,f,b,fp)*C6*(-2*SP)/LL +fq(f)*(
                  DEA(f,f,b,fp)*C9 + DEA(f,f,b,fp)*C7*LR/LL +
		  CEA(f,f,b,fp)*C10 +CEA(f,f,b,fp)*C8*(-2*SP)/LL))*fs(Bd,KP,KM)-
		  (DEA(b,b,f,fp)*C3 + DEA(b,b,f,fp)*C5*LR/LL +
		  CEA(b,b,f,fp)*C4 +CEA(b,b,f,fp)*C6*(-2*SP)/LL +fq(b)*(
                  DEA(b,b,f,fp)*C9 + DEA(b,b,f,fp)*C7*LR/LL +
		  CEA(b,b,f,fp)*C10 +CEA(b,b,f,fp)*C8*(-2*SP)/LL))*fs(Bd,KM,KP)-
		  ((DPAB(u,f,b,fp)+DPAB(d,f,b,fp)+CPAB(fp,f,b,fp)+DPAB(s,f,b,fp)+
			DPAB(c,f,b,fp))*C3+
		  (DPAB(u,f,b,fp)+DPAB(d,f,b,fp)+CPAB(fp,f,b,fp)+DPAB(s,f,b,fp)+
			DPAB(c,f,b,fp))*
			C5*LR/LL+
		  (CPAB(u,f,b,fp)+CPAB(d,f,b,fp)+DPAB(fp,f,b,fp)+CPAB(s,f,b,fp)+
			CPAB(c,f,b,fp))*C4+
		  (CPAB(u,f,b,fp)+CPAB(d,f,b,fp)+DPAB(fp,f,b,fp)+CPAB(s,f,b,fp)+
			CPAB(c,f,b,fp))*C6*
		  (-2*SP)/LL +
                  (DPAB(u,f,b,fp)-1/2*DPAB(d,f,b,fp)-1/2*CPAB(fp,f,b,fp)-1/2*DPAB(s,f,b,fp)+
			DPAB(c,f,b,fp))*C9+
		  (DPAB(u,f,b,fp)-1/2*DPAB(d,f,b,fp)-1/2*CPAB(fp,f,b,fp)-1/2*DPAB(s,f,b,fp)+
			DPAB(c,f,b,fp))*C7*LR/LL+
		  (CPAB(u,f,b,fp)-1/2*CPAB(d,f,b,fp)-1/2*DPAB(fp,f,b,fp)-1/2*CPAB(s,f,b,fp)+
			CPAB(c,f,b,fp))*C10+
		  (CPAB(u,f,b,fp)-1/2*CPAB(d,f,b,fp)-1/2*DPAB(fp,f,b,fp)-1/2*CPAB(s,f,b,fp)+
			CPAB(c,f,b,fp))*C8*
			(-2*SP)/LL)*fs(Bd,KP,KM))*Vbt;

id DPAB(u,f?,b?,fp?)*C1=-gim4(f,b,fp)+DPAB(c,f,b,fp)*C1 
		  - (CPAB(u,f,b,fp)-CPAB(c,f,b,fp))*C2;

id fq(u)=1;
id fq(d)=-1/2;
id fq(s)=-1/2;
id fq(c)=1;

id LL=1;
id LR=1;
id SP=-1/2;

*print;
.sort

id,e1(u?,d?,s?,c?)*fs(Bd?,KP?,KM?)=e1(u,d,s,c,Bd,KP,KM);
id,e2(u?,d?,s?,c?)*fs(Bd?,KP?,KM?)=e2(u,d,s,c,Bd,KP,KM);
id,a1(u?,d?,s?,c?)*fs(Bd?,KP?,KM?)=a1(u,d,s,c,Bd,KP,KM);
id,a2(u?,d?,s?,c?)*fs(Bd?,KP?,KM?)=a2(u,d,s,c,Bd,KP,KM);
id,ea1(u?,d?,s?,c?)*fs(Bd?,KP?,KM?)=ea1(u,d,s,c,Bd,KP,KM);
id,ea2(u?,d?,s?,c?)*fs(Bd?,KP?,KM?)=ea2(u,d,s,c,Bd,KP,KM);
id,p1(u?,d?,s?)*fs(Bd?,KP?,KM?)=p1(u,d,s,Bd,KP,KM);
id,p2(u?,d?,s?)*fs(Bd?,KP?,KM?)=p2(u,d,s,Bd,KP,KM);
id,p3(u?,d?,s?)*fs(Bd?,KP?,KM?)=p3(u,d,s,Bd,KP,KM);
id,p4(u?,d?,s?)*fs(Bd?,KP?,KM?)=p4(u,d,s,Bd,KP,KM);
id,gim1(u?,d?,s?)*fs(Bd?,KP?,KM?)=gim1(u,d,s,Bd,KP,KM);
id,gim2(u?,d?,s?)*fs(Bd?,KP?,KM?)=gim2(u,d,s,Bd,KP,KM);
id,gim3(u?,d?,s?)*fs(Bd?,KP?,KM?)=gim3(u,d,s,Bd,KP,KM);
id,gim4(u?,d?,s?)*fs(Bd?,KP?,KM?)=gim4(u,d,s,Bd,KP,KM);

B lam,Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb,
  Vdu,Vsu,Vbu,Vdc,Vsc,Vbc,Vdt,Vst,Vbt,LL,LR,SP,[Sqr(2)],[Sqr(3)],[Sqr(6)],fs;
print;
.sort

* SU(2) Symmetry

* Special treatment of EW penguins

Id,p1(f?,s?,u,Bd?,KP?,KM?)=p1(f,s,d,Bd,KP,KM)+dp1ew(f,s,u,Bd,KP,KM);
Id,p1(f?,u,s?,Bd?,KP?,KM?)=p1(f,d,s,Bd,KP,KM)+dp1ew(f,u,s,Bd,KP,KM);
Id,p2(f?,s?,u,Bd?,KP?,KM?)=p2(f,s,d,Bd,KP,KM)+dp2ew(f,s,u,Bd,KP,KM);
Id,p2(f?,u,s?,Bd?,KP?,KM?)=p2(f,d,s,Bd,KP,KM)+dp2ew(f,u,s,Bd,KP,KM);
Id,p3(f?,s?,u,Bd?,KP?,KM?)=p3(f,s,d,Bd,KP,KM)+dp3ew(f,s,u,Bd,KP,KM);
Id,p3(f?,u,s?,Bd?,KP?,KM?)=p3(f,d,s,Bd,KP,KM)+dp3ew(f,u,s,Bd,KP,KM);
Id,p3(u,f?,s?,Bd?,KP?,KM?)=p3(d,f,s,Bd,KP,KM)+dp3ew(u,f,s,Bd,KP,KM);
Id,p4(f?,s?,u,Bd?,KP?,KM?)=p4(f,s,d,Bd,KP,KM)+dp4ew(f,s,u,Bd,KP,KM);
Id,p4(f?,u,s?,Bd?,KP?,KM?)=p4(f,d,s,Bd,KP,KM)+dp4ew(f,u,s,Bd,KP,KM);
Id,p4(u,f?,s?,Bd?,KP?,KM?)=p4(d,f,s,Bd,KP,KM)+dp4ew(u,f,s,Bd,KP,KM);

.sort

Set Allowed: e1,e2,a1,a2,ea1,ea2,p1,p2,p3,p4,gim1,gim2,gim3,gim4;

Argument Allowed;
Id,u=d;
EndArgument;
Argument Allowed;
Id,Bd=B;
Id,BP=B;
Id,P0=P;
Id,PM=P;
Id,PP=P;
Id,R0=R;
Id,RM=R;
Id,RP=R;
Id,K0=K;
Id,KP=K;
Id,K0b=Kb;
Id,KM=Kb;
Id,K0S=KS;
Id,KPS=KS;
Id,K0bS=KbS;
Id,KMS=KbS;
Id,D0=D;
Id,DPl=D;
Id,D0b=Db;
Id,DM=Db;
Id,D0S=DS;
Id,DPS=DS;
Id,D0bS=DbS;
Id,DMS=DbS;
EndArgument;

* Unitarieta'
Id,Vts*Vbt=-Vcs*Vbc-Vus*Vbu;
Id,Vtd*Vbt=-Vcd*Vbc-Vud*Vbu;

B lam,Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb,
  Vdu,Vsu,Vbu,Vdc,Vsc,Vbc,Vdt,Vst,Vbt,LL,LR,SP,[Sqr(2)],[Sqr(3)],[Sqr(6)],fs;
print;
.sort

* Special substitutions for the study of golden modes

Id,e1(s?,c,c,f?,B?,D?,Db?)=e1t(s,c,c,f,B,D,Db)-p1(s,c,f,B,D,Db);
Id,gim1(s?,c,f?,B?,D?,Db?)=g1t(s,c,f,B,D,Db)+p1(s,c,f,B,D,Db);
Id,gim3(s?,c,f?,B?,D?,Db?)=g3t(s,c,f,B,D,Db)+p3(s,c,f,B,D,Db);
Id,e2(c,c,u?,f?,B?,JPSI,P?)=e2t(c,c,u,f,B,JPSI,P)-p2(u,c,f,B,JPSI,P);
Id,gim2(u?,c,f?,B?,JPSI,P?)=g2t(u,c,f,B,JPSI,P)+p2(u,c,f,B,JPSI,P);
Id,ea2(c,c,s?,f?,B?,JPSI,ET1?)=ea2t(c,c,s,f,B,JPSI,ET1) - p4(c,s,f,B,JPSI,ET1);
Id,a2(c,s?,c,f?,B?,D?,Db?)=a2t(c,s,c,f,B,D,Db)-p3(s,c,f,B,D,Db);

.sort

* p3, p4, gim3 and gim4 are symmetric under the interchange of q1<->q2 and M1<->M2
Id,p3(u?,c,s?,B?,ET1?,ET8?)=p3(c,u,s,B,ET8,ET1);
Id,gim3(u?,c,s?,B?,ET1?,ET8?)=gim3(c,u,s,B,ET8,ET1);
Id,g3t(u?,c,s?,B?,ET1?,ET8?)=g3t(c,u,s,B,ET8,ET1);
Id,p4(u?,c,f?,B?,ET1?,JPSI?)=p4(c,u,f,B,JPSI,ET1);
Id,gim4(u?,c,f?,B?,ET1?,JPSI?)=gim4(c,u,f,B,JPSI,ET1);

.sort 

Id,gim3(c,s?,u?,B?,D?,Db?)=g3t(c,s,u,B,D,Db)+p3(c,s,u,B,D,Db);
Id,gim4(c,u?,f?,B?,JPSI?,ET1?)=g4t(c,u,f,B,JPSI,ET1)+p4(c,u,f,B,JPSI,ET1);

B lam,Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb,
  Vdu,Vsu,Vbu,Vdc,Vsc,Vbc,Vdt,Vst,Vbt,LL,LR,SP,[Sqr(2)],[Sqr(3)],[Sqr(6)],fs;
print;
.sort

* SU(3) Symmetry
Argument;
Id,s=d;
EndArgument;
Argument;
Id,P=P8;
Id,K=P8;
Id,Kb=P8;
Id,ET8=P8;
Id,R=V8;
Id,KS=V8;
Id,KbS=V8;
Id,OM8=V8;
Id,DsP=D;
Id,DsM=Db;
Id,DsPS=DS;
Id,DsMS=DbS;
EndArgument;

B lam,Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb,
  Vdu,Vsu,Vbu,Vdc,Vsc,Vbc,Vdt,Vst,Vbt,LL,LR,SP,[Sqr(2)],[Sqr(3)],[Sqr(6)],fs;
print;
.end


* 2->3 converted
