(* ::Package:: *)

(* ::Section:: *)
(*All*)


(* ::Subsection:: *)
(*Constants*)


tmax=100;
steps=500;
times=Range[0,tmax,tmax/(steps-1)];


runs=1000;


length=12;


bsites=Range[2,length,2];


j=-1;


fieldsW10={-6.810118375580723`,5.399321183414791`,7.047572832814368`,-4.414079687154078`,-8.41626615919514`,7.9393892160016355`,-5.806999450409989`,-8.53697581757193`,0.6676914555529989`,-2.5431706107186014`,6.844895559224515`,3.891957559769441`};


fieldsW5={-3.4050591877903615`,2.6996605917073957`,3.523786416407184`,-2.207039843577039`,-4.20813307959757`,3.9696946080008177`,-2.9034997252049943`,-4.268487908785965`,0.33384572777649946`,-1.2715853053593007`,3.4224477796122574`,1.9459787798847206`};


(* ::Subsection:: *)
(*Dyn funcs*)


add3[num_]:=Mod[num-1,3]+1


addl[num_]:=Mod[num-1,length]+1


bosdot[su_,s_,a_,ham_]:=-I D[ham,boscc[su][s][a][t]]/.boscc[su][ss_][xx_][t]->bos[su][ss][xx][t]\[Conjugate]


(*bos2dot[s_,a_,ham_]:=-\[ImaginaryI] D[ham,sbos2[s][a][t]]/.sbos2[ss_][xx_][t]->bos2[ss][xx][t]\[Conjugate]*)


(* ::Subsection:: *)
(*Matrices*)


smat[s_]:=PauliMatrix[s]/2


bmat[s1_,s2_]:=KroneckerProduct[PauliMatrix[s1],PauliMatrix[s2]]/If[s1==0||s2==0,2,4]


vecboscc[su_][ss_][t_]:=boscc[su][ss][#][t]&/@Range[su]


vecbos[su_][ss_][t_]:=bos[su][ss][#][t]&/@Range[su]


(*vecbos4nt[ss_]:=bos4[ss][#]&/@Range[4]*)


(*vecbos2nt[ss_]:=bos2[ss][#]&/@Range[2]*)


cS[ss_][xx_][t_]:=Which[MemberQ[bsites,addl[ss]],vecboscc[4][ss][t].bmat[xx,0].vecbos[4][ss][t],MemberQ[addl/@(bsites+1),addl[ss]],vecboscc[4][addl[ss-1]][t].bmat[0,xx].vecbos[4][addl[ss-1]][t],True,vecboscc[2][ss][t].smat[xx].vecbos[2][ss][t]]


(*cS[ss_][xx_]:=Which[MemberQ[bsites,addl[ss]],vecbos4nt[ss]\[Conjugate].bmat[xx,0].vecbos4nt[ss],MemberQ[addl/@(bsites+1),addl[ss]],vecbos4nt[addl[ss-1]]\[Conjugate].bmat[0,xx].vecbos4nt[addl[ss-1]],True,vecbos2nt[ss]\[Conjugate].smat[xx].vecbos2nt[ss]]*)


cB[ss_][xx_,yy_][t_]:=vecboscc[4][ss][t].bmat[xx,yy].vecbos[4][ss][t]


(* ::Subsection:: *)
(*eqns*)


(* ::Subsubsection:: *)
(*SU(4)*)


hamself[ss_]:=fieldsW10[[addl[ss]]]cS[addl[ss]][3][t]


sscoup[ss_]:=j(cS[addl[ss]][1][t]cS[addl[ss+1]][1][t]+cS[addl[ss]][2][t]cS[addl[ss+1]][2][t]+cS[addl[ss]][3][t]cS[addl[ss+1]][3][t])


bcoup[ss_]:=j(cB[addl[ss]][1,1][t]+cB[addl[ss]][2,2][t]+cB[addl[ss]][3,3][t])


hamcoupsu4[ss_]:=If[MemberQ[bsites,addl[ss]],bcoup[addl[ss]],sscoup[addl[ss]]]


(*hamcoupsu2[ss_]:=sscoup[addl[ss]]*)


(*eqss4=Table[bos[2][addl[ss]][aa]'[t]==bos2dot[addl[ss],aa,hamself[addl[ss]]+hamcoupsu4[addl[ss]]+hamcoupsu4[addl[ss-1]]],{ss,Complement[Range[length],bsites,addl/@(bsites+1)]},{aa,2}];*)


eqsb4=Table[bos[4][addl[ss]][aa]'[t]==bosdot[4,addl[ss],aa,hamself[addl[ss]]+hamself[addl[ss+1]]+hamcoupsu4[addl[ss]]+hamcoupsu4[addl[ss-1]]+hamcoupsu4[addl[ss+1]]],{ss,bsites},{aa,4}];


eqall4=Flatten[{(*eqss4,*)eqsb4}];


(*eqss2=Table[bos2[addl[ss]][aa]'[t]==bos2dot[addl[ss],aa,hamself[addl[ss]]+hamcoupsu4[addl[ss]]+hamcoupsu4[addl[ss-1]]],{ss,length},{aa,2}];*)


(*eqall2=Flatten[{eqss2}];*)


(* ::Subsection:: *)
(*Inits*)


initspin={2,2,2,2,2,1,2,1,1,1,1,1};


(* ::Input:: *)
(*vecud[s1_,s2_]:=Flatten[Normal[KroneckerProduct[SparseArray[{initspin[[addl[s1]]]}->1,{2}],SparseArray[{initspin[[addl[s2]]]}->1,{2}]]]]*)


(* ::Input:: *)
(*(*inits4:=Table[Thread[bos[4][addl[ss]][#][0]&/@Range[4]\[Equal]((RandomVariate[NormalDistribution[#,1/2]]+\[ImaginaryI] RandomVariate[NormalDistribution[0,1/2]])&/@vecud[addl[ss],addl[ss+1]])],{ss,bsites}]*)*)


(* ::Input:: *)
(*inits4:=Table[Thread[bos[4][addl[ss]][#][0]&/@Range[4]==(If[#==0,(RandomVariate[NormalDistribution[0,1/2]]+I RandomVariate[NormalDistribution[0,1/2]]),Exp[I RandomReal[{0,2\[Pi]}]]]&/@vecud[addl[ss],addl[ss+1]])],{ss,bsites}]*)


(* ::Subsection:: *)
(*run TWA*)


(* ::Input:: *)
(*through2[funclist_,value_]:=Map[#[value]&,funclist,{2}]*)


(* ::Input:: *)
(*eachTWA4=Table[*)
(*solv=NDSolveValue[Flatten[{eqall4,inits4}],Table[bos[4][addl[ss]][bn],{ss,bsites},{bn,4}],{t,0,tmax}];*)
(*bosvals=Map[through2[solv,#]&,times];*)
(*szvals=Flatten[Map[{#\[Conjugate].bmat[3,0].#,#\[Conjugate].bmat[0,3].#}&,bosvals,{2}],{2,3}];*)
(*szcor=szvals szvals[[All,1]];*)
(*Chop[{szvals,Total[szcor]/length}]*)
(*,{rr,runs}];*)
(*fullTWA4=Total[eachTWA4]/runs;*)
(*sqTWA4=Total[eachTWA4^2]/runs;*)


(*eachTWA2=Table[solv=NDSolveValue[Flatten[{eqall2,initsS}],Flatten[Table[cS[addl[ss]][sp],{ss,length},{sp,3,3}]],{t,0,tmax}];{(Through[solv[#]]&/@times)\[Transpose],Total[(Through[solv[#]]Through[solv[0]]&/@times)\[Transpose]]/length},{rr,runs}];
fullTWA2=Total[eachTWA2]/runs;
sqTWA2=Total[eachTWA2^2]/runs;*)


(*mmu=MaxMemoryUsed[]/10.^6;*)


(*Save["12site.dat",{mmu,fullTWA2,sqTWA2,fullTWA4,sqTWA4}];*)
