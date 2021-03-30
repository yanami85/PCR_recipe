# PCR_recipe
PCRで必要な各種試薬の量やTm値から求めるthermal cyclerの温度設定を自動的に出力するソフト。(required python > 3.7)

## 必要なもの
+ dataclasses
+ Biopython
+ pandas
+ numpy
+ pySimpleGUI

いずれも2020-03-29現在で最新版を使用している。pip or anacondaで適宜入れて下さい。

## 使用方法
1. PCR GUI.pyとpcr_recipe.pyを同じディレクトリに格納する。
2. PCR GUI.pyを走らせる
3. 以下の画面が出てくるので、primer、プライマー濃度、増幅領域、テンプレート濃度 etc. を入力
   (使うメーカーは、2020-03-29現在で**KOD Plus, KOD One, PrimeSTAR**に対応している。KOD FXとか細かいものは適宜増やしていく予定。)
4. 「次へ」をクリックすれば、Thermal cyclerの設定と必要試薬量が一覧で見られる。
5. 「HTMLに出力」をクリックすると1. のディレクトリにHTMLファイルが保存されるので、ブラウザで印刷することも可能。
