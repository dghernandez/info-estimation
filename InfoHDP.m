(* ::Package:: *)

BeginPackage["InfoHDP`"];

Unprotect["InfoHDP`*"];
ClearAll["InfoHDP`*", "InfoHDP`Private`*"];

(******************************************************************************)

genPriorPij::usage = "genPriorPij[alfa, beta, ndist:1, Ns:10^3] generates probability distributions such that q_i ~ DP(alfa_Ns) and q_j|i ~ Beta(beta,beta). It gives ndist probability vectors pij={p1_1, p1_0, p2_1, p2_0, ..., pNs_1, pNs_0}.";
genNastyPij::usage = "genNastyPij[alfa, psure, ndist_:1, Ns:10^4] generates probability distributions such that q_i ~ DP(alfa_Ns) and q_j|i = {psure(Prob=0.25), 0.5(Prob=0.50), 1.-psure(Prob=0.25)}.";
genNastyPij2::usage = "genNastyPij2[psure_, type_:2, ndist_:1, Ns_:10^4] generates probability distributions such that q_i ~ 1/Ns and q_j|i= {psure(Prob=0.33), 0.5(Prob=0.33), 1.-psure(Prob=0.33)}.";
genSamplesPrior::usage = "genSamplesPrior[qij, M, Ns:10^4] generates M samples from prior qij=pij[k].";
Strue::usage = "Strue[p] gives the joint entropy for probability p.";
Sxtrue::usage = "Strue[p] gives the X entropy for probability p.";
Sytrue::usage = "Strue[p] gives the Y entropy for probability p.";
Itrue::usage = "Itrue[p] gives the MI for probability p, assuming pij={p1_1, p1_0, ..., pNs_1, pNs_0}.";
DKm2::usage = "DKm2[samples] -> list {m, Delta_Km}.";
Snaive::usage = "Snaive[N, DKm] gives the maximum likelihood estimator for entropy.";
Smaxlik::usage = "Smaxlik[samples] gives ML entropy estimator from samples.";
Inaive::usage = "Inaive[samples] gives ML mutual information estimator from samples.";

logLa::usage = "logLa[alfa, N, k] posterior log-probability for alfa (NSB). No prior over alfa (like -log(alfa)), but integration in exp(alfa).";
asol::usage = "asol[N, k] gives the alfa solution for nsb, assuming N_states -> Infinity. Also, no prior on alfa";
Spost::usage = "Spost[alfa, N, dkm] gives the posterior entropy at alfa (NSB).";
Smap::usage = "Smap[samples] gives the maximum a posteriori for entropy (NSB).";
D2expalogL::usage = "D2expalogL[Exp[alfa], N, k] second derivative in exp(alfa), to get intervals of integration.";
intEa::usage = "intEa[alfa, N, k, nsig:3.] intervals for integration in exp(alfa).";
Sint::usage = "Sint[samples]={sint, std_sint} gives NSB entropy estimator with the integration from samples. Neglecting error coming from variance for fixed alfa.";
Insb::usage = "Insb[samples] gives NSB estimate for I=Sx+Sy-Sxy.";
InsbCon::usage = "Insb[samples] gives NSB estimate for I=Sx-Sx/y.";

n10sam::usage = "n10sam[samples] gives a vector of size Kx with components {ni_1, ni_0}.";
logLb::usage = "logLb[beta, Kx, n10] gives the posterior log-probability for hyperparameter beta.";
bsol::usage = "bsol[Kx, n10] gives the solution for the hyperparameter beta.";
bsolE::usage = "bsolE[Kx, n10] gives the solution for the hyperparameter beta (looking in Log[beta]).";
SYconX::usage = "SYconX[alfa, beta, N, n10] gives the posterior for the conditional entropy Sy/x.";
IhdpMAP::usage = "IhdpMAP[samples, onlyb_:0.] gives the map estimate of mutual info using InfoHDP. If onlyb==1, then it just using alfa=0. (no pseudocounts)";
D2expblogL::usage = "D2expblogL[Exp[b], Kx, n10] gives the second derivative of logL(b) with respect to Log[b]=eta.";
intEb::usage = "intEb[bx_,kx_,n10_,nsig_:3.] gives interval for integration in log(b).";
IhdpIntb::usage = "IhdpIntb[samples] gives the InfoHDP estimator integrating over peak of posterior (only in b) -> {I_est, dI_est, Sy|x_est}. Variance only coming from uncertainty in beta, and Sy calculated as ML.";
varSYx::usage = "varSYx[beta, n, n'] gives the variance of S(Y|x) for fixed beta and a specific state x with counts n_x=n+n'.";
varSYconX::usage = "varSYconX[beta, nn, n10] gives the variance S(Y|X) for fixed beta.";
SYx2::usage = "SYx2[beta, n, n'] gives the second moment of S(Y|x) for fixed beta and a specific state x with counts n_x=n+n'.";
SYconX2::usage = "SYconX2[beta , n, n'] gives the second moment of S(Y|X) for fixed beta.";

genPriorPijT::usage = "genPriorPijT[alfa_,beta_,qy_,Ns_:10^4] generates probabilities {qx, qy|x, qxy} with prior and marginal qy.";
genSamplesPriorT::usage = "genSamplesPriorT[pi_,pjdadoi_,M_,Ns_:10^4] generates samples {state x, state y} given the probabilities.";
nxysam::usage = "nxysam[sam_,Ny_] counts for each x that occurs, the number of y samples (result in matrix Kx Ny).";
logLbT::usage = "logLbT[beta, qy, nxy] gives the marginal log-likelihood for beta, given a marginal (estimated) qy.";
bsolT::usage = "bsolT[qy_,nxy_] gives the beta that maximazes the marginal log-likelihood, given a estimated qy and counts.";
SYconXT::usage = "SYconXT[bb_,nn_,qy_,nxy_] gives the posterior for the conditional entropy Sy/x.";
IhdpMAPT::usage = "IhdpMAPT[sam_,ML_:0] gives the map estimate of mutual info using InfoHDP. Samples as {Nsam, 2} matrix and each row is an element {xi,yi} in integers.";
InaiveT::usage = "InaiveT[sam_] gets the ML estimate for Ixy, assuming samples as {Nsam, 2} matrix and each row is an element {xi,yi} in integers.";
InsbT::usage = "InsbT[samples] gives NSB estimate for I=Sx+Sy-Sxy, assuming samples as {Nsam, 2} matrix and each row is an element {xi,yi} in integers.";

addline::usage = "addline[plot, lineValue, Mmax:10^5] adds a constant dashed line at lineValue to a plot.";

Begin["`Private`"];

(******************************************************************************)

genPriorPij[alfa_, beta_, ndist_:1, Ns_:10^4] := Block[{pij, alist, bes, pi},
	alist=ConstantArray[alfa/Ns,Ns];
	bes=RandomVariate[BetaDistribution[beta,beta],Ns];
	pi=RandomVariate[DirichletDistribution[alist],ndist];
	pi=Join[pi,{1-(pi//Transpose//Total)}//Transpose,2];
	pij=Table[
		Table[{pi[[k,i]]bes[[i]],pi[[k,i]](1.-bes[[i]])},{i,Ns}]//Flatten
		,{k,ndist}];
	pij];
genSamplesPrior[qij_, M_, Ns_:10^4] := Block[{sam, 
	nameStates={Range[Ns],-Range[Ns]}//Transpose//Flatten},
	sam=RandomChoice[(qij//Abs) -> nameStates,M];(*numerical fluctuactions may make some <0*)
	sam];
genNastyPij[alfa_, psure_, type_:1, ndist_:1, Ns_:10^4]:=Block[{pij, alist, bes, pi,prdel},
	alist=ConstantArray[alfa/Ns,Ns];
	prdel=Which[type==1,{0.25,0.5,0.25},type==2,{1./3,1./3,1./3}];
	bes=RandomChoice[prdel->{psure,0.5,1.-psure},Ns];
	pi=RandomVariate[DirichletDistribution[alist],ndist];
	pi=Join[pi,{1-(pi//Transpose//Total)}//Transpose,2];
	pij=Table[
		Table[{pi[[k,i]]bes[[i]],pi[[k,i]](1.-bes[[i]])},{i,Ns}]//Flatten
		,{k,ndist}];
	pij];
genNastyPij2[psure_, type_:2, ndist_:1, Ns_:10^4]:=Block[{pij, alist, pi, prdel, bes},
	prdel=Which[type==1,{0.25,0.5,0.25},type==2,{1./3,1./3,1./3}];
	bes=RandomChoice[prdel->{psure,0.5,1.-psure},Ns];
	pi=ConstantArray[1./Ns,{ndist,Ns-1}];
	pi=Join[pi,{1-(pi//Transpose//Total)}//Transpose,2];
	pij=Table[
		Table[{pi[[k,i]]bes[[i]],pi[[k,i]](1.-bes[[i]])},{i,Ns}]//Flatten
		,{k,ndist}];
	pij];

(******************************************************************************)

Strue[p_] := (If[#>0.,-# Log[#],0.]&/@p)//Total;
Sxtrue[p_] := (p//ArrayReshape[#,{Length@p/2,2}]&//Transpose//Total//Strue);
Sytrue[p_] := (p//ArrayReshape[#,{Length@p/2,2}]&//Total//Strue);
Itrue[p_] := -Strue[p]+Sytrue[p]+Sxtrue[p];
DKm2[sam_] := (sam//Tally)[[All,2]]//Tally//Sort[#,#1[[1]]<#2[[1]]&]&;
Snaive[nn_,dkm2_] := 1.Sum[-dkm2[[i,2]](dkm2[[i,1]]/nn)Log[dkm2[[i,1]]/nn],{i,Length@dkm2}];
Smaxlik[sam_] := Snaive[sam//Length,DKm2[sam]];
Inaive[sam_] := Block[{Iml=0.,nn=Length@sam,samxz=Abs[#]&/@sam,samyz=Sign[#]&/@sam,
	dkmzX,dkmz,dkmzY},
	dkmz=DKm2[sam];dkmzX=DKm2[samxz];dkmzY=DKm2[samyz];
	Iml=Snaive[nn,dkmzX]+Snaive[nn,dkmzY]-Snaive[nn,dkmz];
	Iml];

(******************************************************************************)

(*NSB functions*)
(*we need to coordinate functions, such that asol takes logL, ans same for D2expalogL, etc*)
logLa[x_,n_,k_] := (k-1.)Log[x]+LogGamma[1 +x]-LogGamma[n+x](*-1.Log[x]*);
asol[nn_, k_] := Block[{solx,x1=1.},
	x1=1. (nn (k/nn)^(3/2)/Sqrt[2 (1-(k/nn))]);(*smart inicial point\[Rule]see NSB_seudocounts_sol.nb*)
	solx=FindRoot[(k-1)/x+PolyGamma[1.+x]-PolyGamma[nn+x]==0,{x,x1,0.001,100. nn}];
	x/.solx];
Spost[x_,nn_,dkm_] := PolyGamma[nn+x+1.]-(1/( x+nn ))(1.Sum[dkm[[i,2]]dkm[[i,1]]PolyGamma[dkm[[i,1]]+1],{i,Length@dkm}]);
Smap[sam_] := Block[{smap=0.,nn=Length@sam,dkmz,kz,ifz,az},
	dkmz=DKm2[sam];
	kz=sam//DeleteDuplicates//Length;
	az=asol[nn,kz];
	(*now let's calculate the maximum a posteriori S_MAP*)
	smap=Spost[az,nn,dkmz];
	smap];
D2expalogL[ex_,n_,k_]:=1.E^ex (PolyGamma[0,1+E^ex]-PolyGamma[0,E^ex+n]+E^ex (PolyGamma[1,1+E^ex]-PolyGamma[1,E^ex+n]));
intEa[xx_,nz_,kz_,nsig_:3.]:=Block[{sigea=0.,ead=0.,eau=0.},(*Intervals for integration in exp(alfa)*)
	sigea=(-D2expalogL[Log[xx],nz,kz])^-0.5;
	ead=Log[xx]-nsig sigea;
	eau=Log[xx]+nsig sigea;
	{ead,eau}];
Sint[sam_]:=Block[{sint=0.,s2int=0.,dsint=0.,
	nn=Length@sam,dkmz,kz,az,ead,eau,npoints=25.,logLaz,listEa,listLogL},(*Integrating Entropy*)
	dkmz=DKm2[sam];
	kz=sam//DeleteDuplicates//Length;
	az=asol[nn,kz];
	logLaz=logLa[az,nn,kz];
	{ead,eau}=intEa[az,nn,kz];
	listEa=Range[ead,eau,(eau-ead)/(npoints-1.)];
	listLogL=(logLa[Exp[#],nn,kz]-logLaz)& /@ listEa;
	listLogL=listLogL//Exp//Normalize[#,Total]&;
	sint=(Spost[Exp[#],nn,dkmz]& /@ listEa).listLogL;
	s2int=((Spost[Exp[#],nn,dkmz]& /@ listEa)^2).listLogL;
	dsint=Sqrt[s2int-sint^2];(*variance only coming from uncertainty in alfa*)
	{sint,dsint}];
Insb[sam_]:=Block[{insb=0.,i2nsb=0.,dinsb=0.,sx,sy,sxy,
	nn=Length@sam,samxz=Abs[#]&/@sam,samyz=Sign[#]&/@sam},
	(*Mutual information with NSB, I=Sx+Sy-Sxy*)
	sx=Sint[samxz][[1]];sy=Sint[samyz][[1]];sxy=Sint[sam][[1]];
	insb=sx+sy-sxy;
	{insb,sx,sy,sxy}];
InsbCon[sam_]:=Block[{insb=0.,i2nsb=0.,dinsb=0.,sx,sx1,sx0,
	nn=Length@sam,nn1,nn0,
	samxz=Abs[#]&/@sam,samx1z=sam//Select[#,#>0&]&,samx0z=sam//Select[#,#<0&]&},
	(*Mutual information with NSB, I=Sx-Sx|y*)
	nn1=Length@samx1z;nn0=Length@samx0z;
	sx=Sint[samxz][[1]];sx1=Sint[samx1z][[1]];sx0=Sint[samx0z][[1]];
	insb=sx-(sx1 nn1+sx0 nn0)/nn;
	{insb,sx,sx1,sx0}];

(******************************************************************************)

(*InfoHDP specific functions*)
n10sam[sam_] := Block[{n10, samx, tsamx},
	samx=Abs[#]&/@sam;
	tsamx=samx//Tally;
	n10=Table[Block[{ni1=0,ni0},
		ni1=Count[sam,tsamx[[tt,1]]];
		ni0=tsamx[[tt,2]]-ni1;{ni1,ni0}],{tt,Length@tsamx}];
	n10];
logLb[b_,kx_,n10_,noprior_:0.] := kx(LogGamma[2.b]-2.LogGamma[b])+Sum[1.LogGamma[1.b+n10[[i,1]]]+LogGamma[1.b+n10[[i,2]]]-LogGamma[2.b+n10[[i,1]]+n10[[i,2]]],{i,Length@n10}]+(1.-noprior)(Log[1.b]+Log[2.PolyGamma[1,2.b+1]-PolyGamma[1,b+1]]);
bsol[kx_,n10_,noprior_:0.]:=(*NArgMax[logLb[Exp[ebb],kx,n10,noprior],ebb]//Exp*)NArgMax[{logLb[Exp[ebb],kx,n10,noprior]&&ebb>Log[10^-3]&&ebb<Log[10^3]},ebb,MaxIterations->1000]//Exp;
bsolE[kx_,n10_]:=Block[{solb},
	solb=FindRoot[0.==1.+Exp[x](2.kx(PolyGamma[2.Exp[x]]-PolyGamma[Exp[x]])+(4.PolyGamma[2,2.Exp[x]+1]-PolyGamma[2,Exp[x]+1])/(2.PolyGamma[1,2.Exp[x]+1]-PolyGamma[1,Exp[x]+1])+Sum[1.PolyGamma[Exp[x]+n10[[i,1]]]+PolyGamma[Exp[x]+n10[[i,2]]]-2.PolyGamma[2Exp[x]+n10[[i,1]]+n10[[i,2]]],{i,Length@n10}]),{x,1.}];
	(x/.solb)//Exp];
SYconX[x_,bb_,nn_,n10_]:=1.(x/(x+nn))(PolyGamma[2.bb+1]-PolyGamma[bb+1.])+(1./(x+nn))Sum[(n10[[i,1]]+n10[[i,2]])(PolyGamma[(n10[[i,1]]+n10[[i,2]])+2.bb+1.]-(n10[[i,1]]+bb)/((n10[[i,1]]+n10[[i,2]])+2.bb) PolyGamma[n10[[i,1]]+bb+1.]-(n10[[i,2]]+bb)/((n10[[i,1]]+n10[[i,2]])+2.bb) PolyGamma[n10[[i,2]]+bb+1.]),{i,Length@n10}];
IhdpMAP[sam_,onlyb_:0,noprior_:0.] := Block[{ihdp=0., a1=0., b1=0., nn=Length@sam, kk, kx, n10, samx,
	sycx=0., sy=0.},
	a1=0.;
	If[onlyb!=1,
		kk=sam//DeleteDuplicates//Length;
		a1=asol[nn,kk];];
	samx=Abs[#]&/@sam;
	kx=samx//DeleteDuplicates//Length;
	n10=n10sam[sam];
	b1=(*NArgMax[logLb[Exp[ebb],kx,n10,noprior],ebb]//Exp*)bsol[kx,n10,noprior];
	sy=(Sign[#]&/@sam)//Smaxlik;(*--------\[Rule]mejorar luego*)
	sycx=SYconX[a1,b1,nn,n10];
	ihdp=sy-sycx;
	ihdp];
D2expblogL[eb_,kx_,n10_,noprior_:0.]:=kx (2 E^eb PolyGamma[0,2 E^eb]-2 (E^eb PolyGamma[0,E^eb]+E^(2 eb) PolyGamma[1,E^eb])+4 E^(2 eb) PolyGamma[1,2 E^eb])+Sum[-2 E^eb PolyGamma[0,2 E^eb+n10[[i,1]]+n10[[i,2]]]+E^eb PolyGamma[0,E^eb+n10[[i,2]]]+E^eb PolyGamma[0,E^eb+n10[[i,1]]]-4 E^(2 eb) PolyGamma[1,2 E^eb+n10[[i,1]]+n10[[i,2]]]+E^(2 eb) PolyGamma[1,E^eb+n10[[i,2]]]+E^(2 eb) PolyGamma[1,E^eb+n10[[i,1]]],{i,Length@n10}]+(1.-noprior)(-(-E^eb PolyGamma[2,1+E^eb]+4 E^eb PolyGamma[2,1+2 E^eb])^2/(-PolyGamma[1,1+E^eb]+2 PolyGamma[1,1+2 E^eb])^2+(-E^eb PolyGamma[2,1+E^eb]+4 E^eb PolyGamma[2,1+2 E^eb]-E^(2 eb) PolyGamma[3,1+E^eb]+8 E^(2 eb) PolyGamma[3,1+2 E^eb])/(-PolyGamma[1,1+E^eb]+2 PolyGamma[1,1+2 E^eb]));
intEb[bx_,kx_,n10_,nsig_:3.,noprior_:0.]:=Block[{sigeb=0.,ebd=0.,ebu=0.},(*Intervals for integration in log[b]*)
	sigeb=(-D2expblogL[Log[bx],kx,n10,noprior])^-0.5;
	ebd=Log[bx]-nsig sigeb;
	ebu=Log[bx]+nsig sigeb;
	{ebd,ebu}];
(*integrating over beta (no alfa here)*)
IhdpIntb[sam_,onlyb_:0,noprior_:0.]:=Block[{ihdp=0., bz=0.,az=0., nn=Length@sam, kk, kx, n10, samx,
	sycx=0., sy=0.,
	sint=0.,s2int=0.,dsint=0.,varsint=0,
	ebd,ebu,npoints=25.,logLbz,listEb,listLogL},(*Integrating Info*)
	If[onlyb!=1,
		kk=sam//DeleteDuplicates//Length;
		az=asol[nn,kk];];	
	samx=Abs[#]&/@sam;
	kx=samx//DeleteDuplicates//Length;
	n10=n10sam[sam];
	bz=(*NArgMax[logLb[Exp[ebb],kx,n10,noprior],ebb]//Exp*)bsol[kx,n10,noprior];
	sy=(Sign[#]&/@sam)//Smaxlik;(*--------\[Rule]mejorar luego*)
	
	logLbz=logLb[bz,kx,n10,noprior];
	{ebd,ebu}=intEb[bz,kx,n10,3.,noprior];
	listEb=Range[ebd,ebu,(ebu-ebd)/(npoints-1.)];
	listLogL=(logLb[Exp[#],kx,n10,noprior]-logLbz)& /@ listEb;
	listLogL=listLogL//Exp//Normalize[#,Total]&;
	
	sint=(SYconX[(*0.*)az,Exp[#],nn,n10]& /@ listEb).listLogL;
	s2int=(SYconX2[Exp[#],nn,n10]& /@ listEb).listLogL;
	(*varsint=(varSYconX[Exp[#],nn,n10]& /@ listEb).listLogL;*)(*this is average variance with fixed beta ~20% of total*)
	dsint=Sqrt[s2int-sint^2];
	(*dsint=Sqrt[varsint];*)
	ihdp=sy-sint;
	{ihdp,dsint,sint}];
varSYx[b_,n0_,n1_]:=(2(b+n1)(b+n0))/((2b+n1+n0)(2b+n1+n0+1)) ((PolyGamma[b+n1+1]-PolyGamma[2b+n1+n0+2])(PolyGamma[b+n0+1]-PolyGamma[2b+n1+n0+2])-PolyGamma[1,2b+n1+n0+2])+((b+n1)(b+n1+1))/((2b+n1+n0)(2b+n1+n0+1)) ((PolyGamma[b+n1+2]-PolyGamma[2b+n1+n0+2])^2+PolyGamma[1,b+n1+2]-PolyGamma[1,2b+n1+n0+2])+((b+n0)(b+n0+1))/((2b+n1+n0)(2b+n1+n0+1)) ((PolyGamma[b+n0+2]-PolyGamma[2b+n1+n0+2])^2+PolyGamma[1,b+n0+2]-PolyGamma[1,2b+n1+n0+2])-(PolyGamma[2b+n1+n0+1]-(b+n1)/(2b+n1+n0) PolyGamma[b+n1+1]-(b+n0)/(2b+n1+n0) PolyGamma[b+n0+1])^2;
SYx2[b_,n0_,n1_]:=(2(b+n1)(b+n0))/((2b+n1+n0)(2b+n1+n0+1)) ((PolyGamma[b+n1+1]-PolyGamma[2b+n1+n0+2])(PolyGamma[b+n0+1]-PolyGamma[2b+n1+n0+2])-PolyGamma[1,2b+n1+n0+2])+((b+n1)(b+n1+1))/((2b+n1+n0)(2b+n1+n0+1)) ((PolyGamma[b+n1+2]-PolyGamma[2b+n1+n0+2])^2+PolyGamma[1,b+n1+2]-PolyGamma[1,2b+n1+n0+2])+((b+n0)(b+n0+1))/((2b+n1+n0)(2b+n1+n0+1)) ((PolyGamma[b+n0+2]-PolyGamma[2b+n1+n0+2])^2+PolyGamma[1,b+n0+2]-PolyGamma[1,2b+n1+n0+2]);
varSYconX[bb_,nn_,n10_]:=(1./(nn^2))Sum[(n10[[i,1]]+n10[[i,2]])^2 varSYx[bb,n10[[i,1]],n10[[i,2]]],{i,Length@n10}];
SYconX2[bb_,nn_,n10_]:=varSYconX[bb,nn,n10]+(SYconX[0.,bb,nn,n10])^2;

(******************************************************************************)
(******************************************************************************)
(*Beyond Binary Labels*)
genPriorPijT[alfa_,beta_,qy_,Ns_:10^4]:=Block[{pij,alist,pjdadoi,pi},alist=ConstantArray[alfa/Ns,Ns];
	pjdadoi=RandomVariate[DirichletDistribution[beta qy],Ns];
	pjdadoi=Join[pjdadoi,{1-(pjdadoi//Transpose//Total)}//Transpose,2];
	pi=RandomVariate[DirichletDistribution[alist]];
	pi=pi~Join~{1.-(pi//Total)};
	pij=Table[pi[[i]]pjdadoi[[i,j]],{j,Length@qy},{i,Ns}];
	{pi,pjdadoi,pij}];
genSamplesPriorT[pi_,pjdadoi_,M_,Ns_:10^4]:=Block[{sam,samx,nameXStates=Range[Ns],Ny=(Dimensions@pjdadoi)[[2]],samy},
	samx=RandomChoice[(pi//Abs)->nameXStates,M];
	sam=Table[{samx[[m]],RandomChoice[(pjdadoi[[samx[[m]] ]]//Abs)-> Range[Ny]]},{m,M}];(*numerical fluctuactions may make some<0*)sam];

nxysam[sam_,Ny_]:=Block[{nxy,samx,tsamx},samx=sam[[All,1]];
	tsamx=samx//Tally;
	nxy=Table[Count[sam,{tsamx[[tt,1]],yy}],{yy,Ny},{tt,Length@tsamx}];
	nxy//Transpose];
logLbT[b_,qy_,nxy_]:=Block[{kx=Length@nxy,Ny=Length@qy,ll=0.},
	ll=kx(LogGamma[b]-Sum[LogGamma[b qy[[j]]],{j,Ny}])+Sum[1.Sum[LogGamma[1.b qy[[j]]+nxy[[i,j]]],{j,Ny}]-LogGamma[b+(nxy[[i]]//Total) ],{i,kx}];
	ll];
bsolT[qy_,nxy_]:=NArgMax[logLbT[Exp[ebb],qy,nxy],ebb,MaxIterations->1000]//Exp;
SYconXT[bb_,nn_,qy_,nxy_]:=Block[{kx=Length@nxy,Ny=Length@qy,ss=0.},
	ss=(1./nn)Sum[(nxy[[i]]//Total)(PolyGamma[(nxy[[i]]//Total)+bb+1.]-Sum[(bb qy[[j]]+nxy[[i,j]])PolyGamma[1.bb qy[[j]]+nxy[[i,j]]+1],{j,Ny}]/((nxy[[i]]//Total)+bb)),{i,Length@nxy}];
	ss];
IhdpMAPT[sam_,ML_:0]:= Block[{ihdp=0., b1=0., nn=Length@sam, ny=sam[[All,2]]//Max, nxy, qye, kx, n10, samx,
	sycx=0., sy=0.},
	nxy=nxysam[sam,ny];
	qye=1.((nxy//Total)+1./ny)/((nxy//Flatten//Total)+1.);(*we have add SG correction to prob*)
	If[ML==1,qye=1.(nxy//Total)/(nxy//Flatten//Total)];
	b1=bsolT[qye,nxy];
	sy=(qye//Strue);(*pa mejorar, con PolyGamma*)
	sycx=SYconXT[b1,nn,qye,nxy];
	ihdp=sy-sycx;
	ihdp];	
InaiveT[sam_] := Block[{Iml=0.,nn=Length@sam,samxz=sam[[All,1]],samyz=sam[[All,2]],
	dkmzX,dkmz,dkmzY},
	dkmz=DKm2[sam];dkmzX=DKm2[samxz];dkmzY=DKm2[samyz];
	Iml=Snaive[nn,dkmzX]+Snaive[nn,dkmzY]-Snaive[nn,dkmz];
	Iml];
InsbT[sam_]:=Block[{insb=0.,i2nsb=0.,dinsb=0.,sx,sy,sxy,
	nn=Length@sam,samxz=sam[[All,1]],samyz=sam[[All,2]]},
	(*Mutual information with NSB, I=Sx+Sy-Sxy*)
	sx=Sint[samxz][[1]];sy=Sint[samyz][[1]];sxy=Sint[sam][[1]];
	insb=sx+sy-sxy;
	{insb,sx,sy,sxy}];

(******************************************************************************)

(*For plots*)
addline[plt_,true_,Mmax_:10^5]:=Show[plt,Plot[true,{m,0,Mmax},PlotStyle->{Dashed,Gray}]];

(******************************************************************************)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["InfoHDP`*"], Head[#] === Symbol &]];

End[];
EndPackage[];
