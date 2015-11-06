(* ::Package:: *)

(* ::Section:: *)
(*All*)


(* ::Subsection:: *)
(*Constants*)


tmax=100;
steps=500;
times=Range[0,tmax,tmax/(steps-1)];


runs=1;


length=12;


bsites=Range[2,length,2];


j=-1;


fieldsW10={-6.810118375580723`,5.399321183414791`,7.047572832814368`,-4.414079687154078`,-8.41626615919514`,7.9393892160016355`,-5.806999450409989`,-8.53697581757193`,0.6676914555529989`,-2.5431706107186014`,6.844895559224515`,3.891957559769441`};


fieldsW5={-3.4050591877903615`,2.6996605917073957`,3.523786416407184`,-2.207039843577039`,-4.20813307959757`,3.9696946080008177`,-2.9034997252049943`,-4.268487908785965`,0.33384572777649946`,-1.2715853053593007`,3.4224477796122574`,1.9459787798847206`};


(* ::Subsection:: *)
(*Dyn funcs*)


add3[num_]:=Mod[num-1,3]+1


addl[num_]:=Mod[num-1,length]+1


bos4dot[s_,a_,ham_]:=-I D[ham,sbos4[s][a][t]]/.sbos4[ss_][xx_][t]->bos4[ss][xx][t]\[Conjugate]


bos2dot[s_,a_,ham_]:=-I D[ham,sbos2[s][a][t]]/.sbos2[ss_][xx_][t]->bos2[ss][xx][t]\[Conjugate]


(* ::Subsection:: *)
(*Matrices*)


smat[s_]:=PauliMatrix[s]/2


bmat[s1_,s2_]:=KroneckerProduct[PauliMatrix[s1],PauliMatrix[s2]]/If[s1==0||s2==0,2,4]


vecsbos2[ss_]:=sbos4[ss][#][t]&/@Range[2]


vecbos2[ss_]:=bos4[ss][#][t]&/@Range[2]


vecsbos4[ss_]:=sbos4[ss][#][t]&/@Range[4]


vecbos4[ss_]:=bos4[ss][#][t]&/@Range[4]


cS[ss_][xx_][t]:=Which[MemberQ[bsites,addl[ss]],vecsbos4[ss].bmat[xx,0].vecbos4[ss],MemberQ[addl/@(bsites+1),addl[ss]],vecsbos4[ss].bmat[0,xx].vecbos4[ss],True,vecsbos2[ss].smat[xx].vecbos2[ss]]


cB[ss_][xx_,yy_][t]:=vecsbos4[ss].bmat[xx,yy].vecbos4[ss]


(* ::Subsection:: *)
(*eqns*)


(* ::Subsubsection:: *)
(*SU(4)*)


hamself[ss_]:=fieldsW10[[addl[ss]]]cS[addl[ss]][3][t]


sscoup[ss_]:=j(cS[addl[ss]][1][t]cS[addl[ss+1]][1][t]+cS[addl[ss]][2][t]cS[addl[ss+1]][2][t]+cS[addl[ss]][3][t]cS[addl[ss+1]][3][t])


bcoup[ss_]:=j(cB[addl[ss]][1,1][t]+cB[addl[ss]][2,2][t]+cB[addl[ss]][3,3][t])


hamcoupsu4[ss_]:=If[MemberQ[bsites,addl[ss]],bcoup[addl[ss]],sscoup[addl[ss]]]


(*hamcoupsu2[ss_]:=sscoup[addl[ss]]*)


eqss4=Table[bos2[addl[ss]][aa]'[t]==bos2dot[addl[ss],aa,hamself[addl[ss]]+hamcoupsu4[addl[ss]]+hamcoupsu4[addl[ss-1]]],{ss,Complement[Range[length],bsites,addl/@(bsites+1)]},{aa,2}];


eqsb4=Table[bos4[addl[ss]][aa]'[t]==bos4dot[addl[ss],aa,hamself[addl[ss]]+hamself[addl[ss+1]]+hamcoupsu4[addl[ss]]+hamcoupsu4[addl[ss-1]]+hamcoupsu4[addl[ss+1]]],{ss,bsites},{aa,4}];


eqall4=Flatten[{eqss4,eqsb4}];


(*eqss2=Table[bos2[addl[ss]][aa]'[t]==bos2dot[addl[ss],aa,hamself[addl[ss]]+hamcoupsu4[addl[ss]]+hamcoupsu4[addl[ss-1]]],{ss,length},{aa,2}];*)


(*eqall2=Flatten[{eqss2}];*)


(* ::Subsection:: *)
(*Inits*)


initspin={2,2,2,2,2,1,2,1,1,1,1,1};


random[mean_,var_]:=If[var!=0,RandomVariate[NormalDistribution[mean,Sqrt[var]]],mean]


(* ::Subsubsection:: *)
(*spins*)


(*mean[op_,vec_]:=vec\[Conjugate].op.vec*)


(*cov[op1_,op2_,vec_]:=vec\[Conjugate].(op1.op2+op2.op1).vec/2-vec\[Conjugate].op1.vec vec\[Conjugate].op2.vec*)


(*spinmean=Table[mean[PauliMatrix[sp]/2,#],{sp,3}]&/@{{1,0},{0,1}};*)


(*spincov=Table[cov[PauliMatrix[sp1]/2,PauliMatrix[sp2]/2,#],{sp1,3},{sp2,3}]&/@{{1,0},{0,1}};*)


(*spincov=Table[0,{2},{3},{3}];*)


(*initsS:=Table[cS[addl[ss]][sp][0]==random[spinmean[[initspin[[addl[ss]]],sp]],spincov[[initspin[[addl[ss]]],sp,sp]]],{ss,length},{sp,3}]*)


(* ::Subsubsection:: *)
(*bispins*)


(*bindlist=Tuples[Range[3],2];*)


(*cor[op1_,op2_,vec_]:=vec\[Conjugate].(op1.op2).vec*)


(*spincor=Table[cor[PauliMatrix[sp1]/2,PauliMatrix[sp2]/2,#],{sp1,3},{sp2,3}]&/@{{1,0},{0,1}};*)


(*bimean[ss_,sp1_,sp2_]:=spinmean[[initspin[[addl[ss]]],sp1]]spinmean[[initspin[[addl[ss+1]]],sp2]]*)


(*bicov[ss_,sp1_,sp2_,sp3_,sp4_]:=(spincor[[initspin[[addl[ss]]],sp1,sp3]]spincor[[initspin[[addl[ss+1]]],sp2,sp4]]+spincor[[initspin[[addl[ss]]],sp3,sp1]]spincor[[initspin[[addl[ss+1]]],sp4,sp2]])/2-
spinmean[[initspin[[addl[ss]]],sp1]]spinmean[[initspin[[addl[ss+1]]],sp2]]spinmean[[initspin[[addl[ss]]],sp3]]spinmean[[initspin[[addl[ss+1]]],sp4]]*)


(*covmat[ss_]:=Outer[bicov@@Join[{addl[ss]},#1,#2]&,bindlist,bindlist,1]*)


(*rotmat=Normalize/@Eigenvectors[covmat[addl[#]]]&/@Range[length];*)


(*rotcov=Eigenvalues[covmat[addl[#]]]&/@Range[length];*)


(*rotcov=Table[0,{12},{9}];*)


(*rotmean=Table[rotmat[[ss]].(bimean@@Join[{addl[ss]},#]&/@bindlist),{ss,length}];*)


(*initsB:=(
initsrot=Table[random[rotmean[[ss,bv]],rotcov[[ss,bv]]],{ss,length},{bv,9}];
initsorigbasis=MapThread[Dot,{Transpose[rotmat,{1,3,2}],initsrot},1];
Table[((cB[addl[ss]])@@bindlist[[bi]])[0]\[Equal]initsorigbasis[[ss,bi]],{ss,bsites},{bi,9}]
)*)


inits


(* ::Subsection:: *)
(*run TWA*)


eachTWA4=Table[solv=NDSolveValue[Flatten[{eqall4,initsS,initsB}],Flatten[Table[cS[addl[ss]][sp],{ss,length},{sp,3,3}]],{t,0,tmax}];{(Through[solv[#]]&/@times)\[Transpose],Total[(Through[solv[#]]Through[solv[0]]&/@times)\[Transpose]]/length},{rr,runs}];
fullTWA4=Total[eachTWA4]/runs;
sqTWA4=Total[eachTWA4^2]/runs;


eachTWA2=Table[solv=NDSolveValue[Flatten[{eqall2,initsS}],Flatten[Table[cS[addl[ss]][sp],{ss,length},{sp,3,3}]],{t,0,tmax}];{(Through[solv[#]]&/@times)\[Transpose],Total[(Through[solv[#]]Through[solv[0]]&/@times)\[Transpose]]/length},{rr,runs}];
fullTWA2=Total[eachTWA2]/runs;
sqTWA2=Total[eachTWA2^2]/runs;


mmu=MaxMemoryUsed[]/10.^6;


Save["12site.dat",{mmu,fullTWA2,sqTWA2,fullTWA4,sqTWA4}];
