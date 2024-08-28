(* ::Package:: *)

BeginPackage["NumericRay`"];


NRay::usage = "{a, \[Lambda], \[Eta], \!\(\*SuperscriptBox[SubscriptBox[\(\[Nu]\), \(r\)], \(s\)]\), \!\(\*SuperscriptBox[SubscriptBox[\(\[Nu]\), \(\[Theta]\)], \(s\)]\), \!\(\*SubscriptBox[\(r\), \(s\)]\), \!\(\*SubscriptBox[\(\[Theta]\), \(s\)]\), \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\), \!\(\*SubscriptBox[\(t\), \(s\)]\), pt0->1, precisiongoal->15, accuracygoal->15, precision->20, smax->20000, toprint->True, toplot->True, splot->200, imagesize->250}";
Demo::usage = "{a, \[Lambda], \[Eta], \!\(\*SuperscriptBox[SubscriptBox[\(\[Nu]\), \(r\)], \(s\)]\), \!\(\*SuperscriptBox[SubscriptBox[\(\[Nu]\), \(\[Theta]\)], \(s\)]\), \!\(\*SubscriptBox[\(r\), \(s\)]\), \!\(\*SubscriptBox[\(\[Theta]\), \(s\)]\), \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\), smax}";


Begin["`Private`"]


(* ::Section:: *)
(*NRay*)


Options[NRay]={"pt0"->1,"precisiongoal"->15,"accuracygoal"->15,"precision"->20,"smax"->20000,"toprint"->True,"toplot"->True,"splot"->200,"imagesize"->250};
NRay[a_,\[Lambda]_,\[Eta]_,\[Nu]r_,\[Nu]\[Theta]_,rs_,\[Theta]s_,\[Phi]s_,ts_,OptionsPattern[]]:=NG[a,\[Lambda],\[Eta],\[Nu]r,\[Nu]\[Theta],rs,\[Theta]s,\[Phi]s,ts,
	OptionValue["pt0"],OptionValue["precisiongoal"],OptionValue["accuracygoal"],OptionValue["precision"],
	OptionValue["smax"],OptionValue["toprint"],OptionValue["toplot"],OptionValue["splot"],OptionValue["imagesize"]];


NG[a_,\[Lambda]_,\[Eta]_,\[Nu]r_,\[Nu]\[Theta]_,rs_,\[Theta]s_,\[Phi]s_,ts_,pt0input_,precisiongoal_,accuracygoal_,precision_,smax_,toprint_,toplot_,splot_,imagesize_]:=Block[
{chris,eqt,eqr,eq\[Theta],eq\[Phi],momentum,eqpt,eqpr,eqp\[Theta],eqp\[Phi],
gtt,gt\[Phi],g\[Phi]\[Phi],g\[Theta]\[Theta],grr,pt0,pr0,p\[Theta]0,p\[Phi]0,pdes,sol,coordinates},

(*christoffel symbols*)
chris={
	{
		{0,-((2 (a^2+\[ScriptR]^2) (a^2-2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]]))/((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)),-((4 a^2 \[ScriptR] Sin[2 \[CurlyTheta]])/(a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),0},
		{-((2 (a^2+\[ScriptR]^2) (a^2-2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]]))/((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)),0,0,(2 a (a^4-3 a^2 \[ScriptR]^2-6 \[ScriptR]^4+a^2 (a^2-\[ScriptR]^2) Cos[2 \[CurlyTheta]]) Sin[\[CurlyTheta]]^2)/((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)},
		{-((4 a^2 \[ScriptR] Sin[2 \[CurlyTheta]])/(a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),0,0,(8 a^3 \[ScriptR] Cos[\[CurlyTheta]] Sin[\[CurlyTheta]]^3)/(a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2},
		{0,(2 a (a^4-3 a^2 \[ScriptR]^2-6 \[ScriptR]^4+a^2 (a^2-\[ScriptR]^2) Cos[2 \[CurlyTheta]]) Sin[\[CurlyTheta]]^2)/((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),(8 a^3 \[ScriptR] Cos[\[CurlyTheta]] Sin[\[CurlyTheta]]^3)/(a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2,0}
	},
	{
		{-(((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2-2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]]))/(2 (\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3)),0,0,(a (a^2+(-2+\[ScriptR]) \[ScriptR]) (-\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2) Sin[\[CurlyTheta]]^2)/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3},
		{0,((a^2-\[ScriptR]) \[ScriptR]+a^2 (1-\[ScriptR]) Cos[\[CurlyTheta]]^2)/((a^2+(-2+\[ScriptR]) \[ScriptR]) (\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)),-((a^2 Cos[\[CurlyTheta]] Sin[\[CurlyTheta]])/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)),0},
		{0,-((a^2 Cos[\[CurlyTheta]] Sin[\[CurlyTheta]])/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)),-((\[ScriptR] (a^2+(-2+\[ScriptR]) \[ScriptR]))/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)),0},
		{(a (a^2+(-2+\[ScriptR]) \[ScriptR]) (-\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2) Sin[\[CurlyTheta]]^2)/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3,0,0,-(((a^2+(-2+\[ScriptR]) \[ScriptR]) Sin[\[CurlyTheta]]^2 (Cos[\[CurlyTheta]]^2 (2 a^2 \[ScriptR] (a^2+\[ScriptR]^2)+a^4 (1-\[ScriptR]) Sin[\[CurlyTheta]]^2)+\[ScriptR] (-a^4+\[ScriptR]^4+a^2 (a^2-\[ScriptR]) Sin[\[CurlyTheta]]^2)))/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3)}
	},
	{
		{-((a^2 \[ScriptR] Sin[2 \[CurlyTheta]])/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3),0,0,(a \[ScriptR] (a^2+\[ScriptR]^2) Sin[2 \[CurlyTheta]])/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3},
		{0,(a^2 Cos[\[CurlyTheta]] Sin[\[CurlyTheta]])/((a^2+(-2+\[ScriptR]) \[ScriptR]) (\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)),\[ScriptR]/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2),0},
		{0,\[ScriptR]/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2),-((a^2 Cos[\[CurlyTheta]] Sin[\[CurlyTheta]])/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)),0},
		{(a \[ScriptR] (a^2+\[ScriptR]^2) Sin[2 \[CurlyTheta]])/(\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3,0,0,-1/(16 (\[ScriptR]^2+a^2 Cos[\[CurlyTheta]]^2)^3) (3 a^6+10 a^4 \[ScriptR]+11 a^4 \[ScriptR]^2+16 a^2 \[ScriptR]^3+16 a^2 \[ScriptR]^4+8 \[ScriptR]^6+4 a^2 (a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2) Cos[2 \[CurlyTheta]]+a^4 (a^2-2 \[ScriptR]+\[ScriptR]^2) Cos[4 \[CurlyTheta]]) Sin[2 \[CurlyTheta]]}
	},
	{
		{0,-((2 a (a^2-2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]]))/((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)),-((8 a \[ScriptR] Cot[\[CurlyTheta]])/(a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),0},
		{-((2 a (a^2-2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]]))/((a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)),0,0,(a^4+3 a^4 \[ScriptR]-12 a^2 \[ScriptR]^2+8 a^2 \[ScriptR]^3-16 \[ScriptR]^4+8 \[ScriptR]^5+4 a^2 \[ScriptR] (a^2+\[ScriptR] (-1+2 \[ScriptR])) Cos[2 \[CurlyTheta]]-a^4 (1-\[ScriptR]) Cos[4 \[CurlyTheta]])/(2 (a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)},
		{-((8 a \[ScriptR] Cot[\[CurlyTheta]])/(a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),0,0,((3 a^4+8 a^2 \[ScriptR]+8 a^2 \[ScriptR]^2+8 \[ScriptR]^4+4 a^2 (a^2+2 (-1+\[ScriptR]) \[ScriptR]) Cos[2 \[CurlyTheta]]+a^4 Cos[4 \[CurlyTheta]]) Cot[\[CurlyTheta]])/(2 (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2)},
		{0,(a^4+3 a^4 \[ScriptR]-12 a^2 \[ScriptR]^2+8 a^2 \[ScriptR]^3-16 \[ScriptR]^4+8 \[ScriptR]^5+4 a^2 \[ScriptR] (a^2+\[ScriptR] (-1+2 \[ScriptR])) Cos[2 \[CurlyTheta]]-a^4 (1-\[ScriptR]) Cos[4 \[CurlyTheta]])/(2 (a^2+(-2+\[ScriptR]) \[ScriptR]) (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),((3 a^4+8 a^2 \[ScriptR]+8 a^2 \[ScriptR]^2+8 \[ScriptR]^4+4 a^2 (a^2+2 (-1+\[ScriptR]) \[ScriptR]) Cos[2 \[CurlyTheta]]+a^4 Cos[4 \[CurlyTheta]]) Cot[\[CurlyTheta]])/(2 (a^2+2 \[ScriptR]^2+a^2 Cos[2 \[CurlyTheta]])^2),0}
	}
};


(*Equations of motion*)
eqt=(t'[s]==pt[s]);
eqr=(r'[s]==pr[s]);
eq\[Theta]=(\[Theta]'[s]==p\[Theta][s]);
eq\[Phi]=(\[Phi]'[s]==p\[Phi][s]);
momentum[i_]:={pt[s],pr[s],p\[Theta][s],p\[Phi][s]}[[i]];
eqpt=(pt'[s]==-Sum[chris[[1,ii,jj]]momentum[ii]momentum[jj],{ii,1,4},{jj,1,4}]/.{\[ScriptT]->t[s],\[ScriptR]->r[s],\[CurlyTheta]->\[Theta][s],\[CurlyPhi]->\[Phi][s]});
eqpr=(pr'[s]==-Sum[chris[[2,ii,jj]]momentum[ii]momentum[jj],{ii,1,4},{jj,1,4}]/.{\[ScriptT]->t[s],\[ScriptR]->r[s],\[CurlyTheta]->\[Theta][s],\[CurlyPhi]->\[Phi][s]});
eqp\[Theta]=(p\[Theta]'[s]==-Sum[chris[[3,ii,jj]]momentum[ii]momentum[jj],{ii,1,4},{jj,1,4}]/.{\[ScriptT]->t[s],\[ScriptR]->r[s],\[CurlyTheta]->\[Theta][s],\[CurlyPhi]->\[Phi][s]});
eqp\[Phi]=(p\[Phi]'[s]==-Sum[chris[[4,ii,jj]]momentum[ii]momentum[jj],{ii,1,4},{jj,1,4}]/.{\[ScriptT]->t[s],\[ScriptR]->r[s],\[CurlyTheta]->\[Theta][s],\[CurlyPhi]->\[Phi][s]});


(*Covariant metrics*)
gtt=-((a^2-2 M r+r^2-a^2 Sin[\[CurlyTheta]]^2)/(r^2+a^2 Cos[\[CurlyTheta]]^2))/.{M->1};
gt\[Phi]=-((2 a M r Sin[\[CurlyTheta]]^2)/(r^2+a^2 Cos[\[CurlyTheta]]^2))/.{M->1};
g\[Phi]\[Phi]=(Sin[\[CurlyTheta]]^2 ((a^2+r^2)^2-a^2 (a^2+r (-2 M+r)) Sin[\[CurlyTheta]]^2))/(r^2+a^2 Cos[\[CurlyTheta]]^2)/.{M->1};
g\[Theta]\[Theta]=r^2+a^2 Cos[\[CurlyTheta]]^2;
grr=(r^2+a^2 Cos[\[CurlyTheta]]^2)/(a^2-2 M r+r^2)/.{M->1};


(*Contravariant momenta*)
pt0=pt0input;
p\[Phi]0=-((gt\[Phi]+\[Lambda] gtt)/(g\[Phi]\[Phi]+\[Lambda] gt\[Phi])) pt0/.{r->rs,\[CurlyTheta]->\[Theta]s};
p\[Theta]0=\[Nu]\[Theta] Sqrt[\[Eta]+a^2 Cos[\[CurlyTheta]]^2-\[Lambda]^2 Cot[\[CurlyTheta]]^2] (-(gtt pt0+gt\[Phi] p\[Phi]0))/g\[Theta]\[Theta]/.{r->rs,\[CurlyTheta]->\[Theta]s};
pr0=\[Nu]r Sqrt[-(gtt pt0^2+g\[Theta]\[Theta] p\[Theta]0^2+g\[Phi]\[Phi] p\[Phi]0^2+2 gt\[Phi] pt0 p\[Phi]0)/grr]/.{r->rs,\[CurlyTheta]->\[Theta]s};
pdes={eqt,eqr,eq\[Theta],eq\[Phi],eqpt,eqpr,eqp\[Theta],eqp\[Phi],t[0]==ts,r[0]==rs,\[Theta][0]==\[Theta]s,\[Phi][0]==\[Phi]s,pt[0]==pt0,pr[0]==pr0,p\[Theta][0]==p\[Theta]0,p\[Phi][0]==p\[Phi]0};


(*Numeric integration*)
sol=NDSolve[SetPrecision[{eqt,eqr,eq\[Theta],eq\[Phi],eqpt,eqpr,eqp\[Theta],eqp\[Phi],t[0]==ts,r[0]==rs,\[Theta][0]==\[Theta]s,\[Phi][0]==\[Phi]s,pt[0]==pt0,pr[0]==pr0,p\[Theta][0]==p\[Theta]0,p\[Phi][0]==p\[Phi]0},precision],{t,r,\[Theta],\[Phi],pt,pr,p\[Theta],p\[Phi]},{s,0,smax},
PrecisionGoal->precisiongoal,
AccuracyGoal->accuracygoal,
WorkingPrecision->precision
];

coordinates[i_]:={t,r,\[Theta],\[Phi]}[[i]];

(*If[toprint,

Print["Initial momenta:\n{\!\(\*FractionBox[SuperscriptBox[\(p\), \(r\)], \(\(|\)\*SuperscriptBox[\(p\), \(r\)]\(|\)\)]\),\!\(\*FractionBox[\(\(\\\ \)\(r\\\ \*SuperscriptBox[\(p\), \(\[Theta]\)]\)\), \(\(|\)\*SuperscriptBox[\(p\), \(r\)]\(|\)\)]\),\!\(\*FractionBox[\(r\\\ Sin[\[Theta]] \*SuperscriptBox[\(p\), \(\[Phi]\)]\), \(\(|\)\*SuperscriptBox[\(p\), \(r\)]\(|\)\)]\)} = ",{pr0/Abs[pr0],(rs p\[Theta]0)/Abs[pr0],(rs Sin[\[Theta]s]p\[Phi]0)/Abs[pr0]}];

Print["\!\(\*SubscriptBox[\(r\), \(f\)]\) = ",First[r[smax]/.sol]];
Print["\!\(\*SubscriptBox[\(\[Theta]\), \(\[Infinity]\)]\) = ",First[\[Theta][smax]/.sol]];
Print["\!\(\*SubscriptBox[\(\[Phi]\), \(\[Infinity]\)]\) = ",First[\[Phi][smax]/.sol]];
Print["Mod[\!\(\*SubscriptBox[\(\[Phi]\), \(\[Infinity]\)]\), 2\[Pi]] = ",Mod[First[\[Phi][smax]/.sol],2\[Pi]]];
];*)

If[toplot,
Echo@Table[Plot[coordinates[ii][s]/.sol,{s,0,splot},PlotRange->All,Frame->True,FrameStyle->Directive[14,Black],FrameLabel->{"s",{"t","r","\[Theta]","\[Phi]"}[[ii]]},ImageSize->imagesize],{ii,2,4}]];


sol[[1,;;,2]]
]


(* ::Section:: *)
(*3D Demo*)


Options[Demo] = {"length$spin" -> 7, "arrowheads$spin" -> 0.02, "arrowthickness$spin" -> 0.005,
   "size$source" -> 0.7, "size$bh" -> "default", "plotstyle" -> Automatic, "boxed" -> True,
   "plotrange" -> All, "label" -> "", "labelpos" -> {0.5, 0.5, 1.2}, "regionfunc" -> (True &), "coords?" -> False, "coords$length" -> 10,
   "viewpoint" -> {1.3, -2.4, 2.}, "viewprojection" -> "Perspective"};
Demo[a_, \[Lambda]_, \[Eta]_, \[Nu]rs_, \[Nu]\[Theta]s_, rs_, \[Theta]s_, \[Phi]s_, smax_, OptionsPattern[]] := Module[
  {size$source, size$bh, length$spin, arrowheads$spin, arrowthickness$spin, plotstyle, boxed, plotrange, regionfunc,
   bg, coordinates, ray, xyz, plt},
  size$source = OptionValue["size$source"];
  size$bh = OptionValue["size$bh"]; If[size$bh == "default", size$bh = 1 + Sqrt[1 - a^2]];
  length$spin = OptionValue["length$spin"];
  arrowheads$spin = OptionValue["arrowheads$spin"];
  arrowthickness$spin = OptionValue["arrowthickness$spin"];
  plotstyle = OptionValue["plotstyle"];
  boxed = OptionValue["boxed"];
  plotrange = OptionValue["plotrange"];
  regionfunc = OptionValue["regionfunc"];
  
  (*Source point, BH, coordinate system*)
  bg := Graphics3D[{
     Red, Sphere[{rs  Sin[\[Theta]s] Cos[\[Phi]s], rs  Sin[\[Theta]s] Sin[\[Phi]s], rs  Cos[\[Theta]s]}, size$source],
     Black, Sphere[{0, 0, 0}, size$bh],
     Black, Arrowheads[arrowheads$spin], Thickness[arrowthickness$spin], Arrow[{{0, 0, 0}, {0, 0, length$spin}}],
     Text[Style[OptionValue["label"], Directive[Black, 15]], Scaled[OptionValue["labelpos"]]]
     }
    , Boxed -> boxed, Axes->True, AxesLabel->{"x", "y", "z"}];
  coordinates = Graphics3D[{Black, Arrowheads[0.04], Thickness[0.005],
     	Arrow[OptionValue["coords$length"] {{0, 0, 0}, {1, 0, 0}}],
     	Arrow[OptionValue["coords$length"] {{0, 0, 0}, {0, 1, 0}}],
     	Arrow[OptionValue["coords$length"] {{0, 0, 0}, {0, 0, 1}}],
     	Text[Style[TraditionalForm@"x", Directive[Black, Italic, 14]], OptionValue["coords$length"] {1.2, 0, 0}],
     	Text[Style[TraditionalForm@"y", Directive[Black, Italic, 14]], OptionValue["coords$length"] {0, 1.35, 0}],
     	Text[Style[TraditionalForm@"z", Directive[Black, Italic, 14]], OptionValue["coords$length"] {0, 0, 1.2}]
     	}];
  
  (*Geodesics*)
  ray = NRay[a, \[Lambda], \[Eta], \[Nu]rs, \[Nu]\[Theta]s, rs, \[Theta]s, \[Phi]s, 0, "toprint" -> False, "toplot" -> False];
  xyz[ray_, s_] := Module[{r, \[Theta], \[Phi]},
    r = ray[[2]][s];
    \[Theta] = ray[[3]][s];
    \[Phi] = ray[[4]][s];
    {r  Sin[\[Theta]] Cos[\[Phi]], r  Sin[\[Theta]] Sin[\[Phi]], r  Cos[\[Theta]]}];
  
  
  (*Plot*)
  plt = {bg};
  If[OptionValue["coords?"] == True, AppendTo[plt, coordinates]];
  plt = Join[
    	plt,
    	{
     	ParametricPlot3D[xyz[ray, s], {s, 0, smax}, PlotStyle -> plotstyle, RegionFunction -> regionfunc],
     	PlotRange -> plotrange,
     	ViewPoint -> OptionValue["viewpoint"],
     	ViewProjection -> OptionValue["viewprojection"]
     	}];
  	
  (*Print the end point*)
  Print["End point: {r, \[Theta], \[Phi]} = ", {ray[[2]][smax],ray[[3]][smax],ray[[4]][smax]}, " for s = ",smax];
  
  Apply[Show, plt]	
  ]


(* ::Section:: *)
(*End*)


End[];


EndPackage[];
