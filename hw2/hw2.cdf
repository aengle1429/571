(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 9.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1063,         20]
NotebookDataLength[      3368,        125]
NotebookOptionsPosition[      3804,        115]
NotebookOutlinePosition[      4264,        136]
CellTagsIndexPosition[      4221,        133]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{
  RowBox[{"f", "[", 
   RowBox[{"b0_", ",", "b1_", ",", "b2_", ",", "b3_"}], "]"}], ":=", 
  RowBox[{"{", 
   RowBox[{"b0", ",", 
    RowBox[{"b1", "/", "b0"}], ",", 
    RowBox[{"Sqrt", "[", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"b1", "+", "b2"}], ")"}], "/", "b1"}], "]"}], ",", 
    RowBox[{"Sqrt", "[", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"b1", "+", "b3"}], ")"}], "/", "b1"}], "]"}]}], 
   "}"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"FullSimplify", "[", 
   RowBox[{"D", "[", 
    RowBox[{
     RowBox[{"f", "[", 
      RowBox[{"b0", ",", "b1", ",", "b2", ",", "b3"}], "]"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"b0", ",", "b1", ",", "b2", ",", "b3"}], "}"}], ",", "1"}], 
      "}"}]}], "]"}], "]"}], "//", "MatrixForm"}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1", "0", "0", "0"},
     {
      RowBox[{"-", 
       FractionBox["b1", 
        SuperscriptBox["b0", "2"]]}], 
      FractionBox["1", "b0"], "0", "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["b2", 
        RowBox[{"2", " ", 
         SuperscriptBox["b1", "2"], " ", 
         SqrtBox[
          FractionBox[
           RowBox[{"b1", "+", "b2"}], "b1"]]}]]}], 
      FractionBox["1", 
       RowBox[{"2", " ", "b1", " ", 
        SqrtBox[
         FractionBox[
          RowBox[{"b1", "+", "b2"}], "b1"]]}]], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["b3", 
        RowBox[{"2", " ", 
         SuperscriptBox["b1", "2"], " ", 
         SqrtBox[
          FractionBox[
           RowBox[{"b1", "+", "b3"}], "b1"]]}]]}], "0", 
      FractionBox["1", 
       RowBox[{"2", " ", "b1", " ", 
        SqrtBox[
         FractionBox[
          RowBox[{"b1", "+", "b3"}], "b1"]]}]]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output"]
}, Open  ]]
},
WindowSize->{740, 651},
Visible->True,
ScrollingOptions->{"VerticalScrollRange"->Fit},
ShowCellBracket->False,
Deployed->True,
CellContext->Notebook,
TrackCellChangeTimes->False,
FrontEndVersion->"9.0 for Mac OS X x86 (32-bit, 64-bit Kernel) (January 25, \
2013)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[1463, 33, 463, 15, 28, "Input"],
Cell[CellGroupData[{
Cell[1951, 52, 375, 11, 28, "Input"],
Cell[2329, 65, 1459, 47, 154, "Output"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature kwD0zQlsDwGa7AgYUgwzSwMh *)
